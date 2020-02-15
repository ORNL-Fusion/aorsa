!module METS2AORSA

!  private
!  public:: WMATPRECALC_AORSA,GETNONMAXSWMAT_AORSA,GETNONMAXSWMAT_AORSA_NEW

!contains

!
!*************************************************************************
!


    subroutine getnonmaxsigma_aorsa_newi_1(w,zspec,aspec,dens,bmag, &
       & k1,xi1,jnxi1,k2,xi2,jnxi2,nbessj,enorm,uparmin,uparmax, &
       & nupar,nuper,uper,upar,dfduper,dfdupar,wspec,ifail, &
	 & l_first, l_interp, kperp_max, nkperp, xkprl0)

!   ---------------------------------------------------------
!   lee's version: interpolates for l_interp = .true
!           does full integrals for l_interp = .false.
!   ---------------------------------------------------------
	
    implicit none


    logical, optional :: l_interp
    logical, optional :: l_first
    integer, optional :: nkperp
    real, optional :: kperp_max
    real :: kperp_max_l
    real :: dk, kperp, p
    integer :: nk, nkperp_l
    integer :: ll, mm, nn
    !nupar, -2 to nkper =2, 0 to nharm=1
    real, save, allocatable, dimension(:,:,:) :: fpint_11, fpint_33
    real, save, allocatable, dimension(:,:) :: fpint_31   ! no nharm
    real, save, allocatable, dimension(:,:) :: pint_11, pint_33
    real, save, allocatable, dimension(:) :: pint_31
    real, intent(in) :: w,zspec,aspec,dens,bmag
    real, dimension(3) :: k1,k2
    integer, intent(in) :: nupar,nuper,nbessj
    real, dimension(nuper) :: xi1,xi2
    real, dimension(nuper, nbessj), intent(in) :: jnxi1,jnxi2
    real, intent(in) :: enorm,uparmin,uparmax
    real, dimension(nuper), intent(in) :: uper
    real, dimension(nupar), intent(in):: upar
    real, dimension(nuper,nupar), intent(in) :: dfduper,dfdupar
    real, save, allocatable, dimension(:,:) ::  dfdthp
    complex, dimension(3,3), intent(inout) :: wspec
    complex :: sig_11, sig_21, sig_22, sig_31, sig_32, sig_33
    complex :: sig_fact
    integer, intent(inout) :: ifail
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07

    real, parameter:: wp2fact=1.745915
    real, parameter:: mpc2=938271998.38
    real, parameter:: c=2.99792458e8
    real, parameter:: pi=3.141592653597932384
    real :: ga_33(nupar)
    real :: w2,wp2,wpfact,wcw,rrp,rrm,rr,irr, du_3, g_33j, xkprl0
    real, save :: g_33, n_save, b_save
    real :: mut0,sqmut0,kpara1, kpara0, npara1, npara0, verp_norm
    real :: isq2,cupar,cuper,nwcw,dfactpar,dfactper, dfactper0, lf0,lnf0uper
    real :: dfactperi, a1, a2, a3, a4, a5,beta1,beta2
    real :: upar0, factor, kperp2, vperp_norm2
    real :: ssgg_11, ssgg_31, ssgg_33, du, temp, kperp_f
    integer :: nharm, iharm, j, k, m, n, dn
    complex :: wref(3,3)
    real :: rerror, ierror, i_fact, r_fact, vperp_norm
    real :: thtrp_11,thtrp_21,thtrp_22,thtrp_31,thtrp_32,thtrp_33
    real :: o_thtrp_11,o_thtrp_21,o_thtrp_22,o_thtrp_31,o_thtrp_32,o_thtrp_33
    real :: e_thtrp_11,e_thtrp_21,e_thtrp_22,e_thtrp_31,e_thtrp_32,e_thtrp_33
    real :: thtip_11,thtip_21,thtip_22,thtip_31,thtip_32,thtip_33
    real :: thtrm_11,thtrm_21,thtrm_22,thtrm_31,thtrm_32,thtrm_33
    real :: o_thtrm_11,o_thtrm_21,o_thtrm_22,o_thtrm_31,o_thtrm_32,o_thtrm_33
    real :: e_thtrm_11,e_thtrm_21,e_thtrm_22,e_thtrm_31,e_thtrm_32,e_thtrm_33
    real :: thtim_11,thtim_21,thtim_22,thtim_31,thtim_32,thtim_33
    real, dimension(nupar) :: temp_11,temp_21,temp_22,temp_31,temp_32,temp_33
    real, dimension(nupar) :: sg_11,sg_21,sg_22,sg_31,sg_32,sg_33




    logical, parameter :: use_ppart6 = .true.
    logical :: is_uniform 
    integer, parameter :: nfxmax = 9
    real*8, dimension(nfxmax) ::  vint
    real*8, dimension(nupar,nfxmax) :: fx
    integer :: nfx
    real*8 :: dx,dh,tol
    integer :: i

!  check to make sure that the v_par mesh is uniform
!  if we need to check par, why not perp as well?

    is_uniform = .true.
    tol = 1.0d-7
    dh = upar(2)-upar(1)
    do i=1,nupar-1
       dx = upar(i+1)-upar(i)
       is_uniform = abs(dx-dh).le. tol*dx
       if (.not.is_uniform) exit
    enddo


!  create a local value of nkperp
    nkperp_l = nkperp !don't change the calling arguement.

    !are the storage arrays allocated?
    !nupar; -2 to nkperp + 2, 0 to nharm + 1 (0 to nbessj -1)
    if(.not. allocated(fpint_11)) then
       if(nkperp_l .lt. 0) then
          nkperp_l = 0  !need to have at least one slot in the nk dimension
       end if !end if set nkperp_1
    
!  allcate storate arrays for interpolation/integration
!  fpint group stores values for interpolation
!  three are needed:  one each for sig_11, sig_f33, and sig_31
!  pint stores the values of the perpendicular integral for the actual k_perp*vperp/omega_c
!  for fpint, the range is:  the parallel mesh; the perpendicular mesh +/- 2 for interpolation;
!  and the harmonic index  (no offset)
!  dfdthp is a common factor for all of sigma that is closely related to dfdtheta

       allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_31(1:nupar, -2:nkperp_l+2),  &
            &   pint_11(1:nupar, 0:nbessj-1),	&
            &   pint_33(1:nupar, 0:nbessj-1),	&
            &   pint_31(1:nupar),               &
            &   dfdthp(1:nuper,1:nupar))
    else  !check for size
       if((size(fpint_11,2) .ne. nkperp_l + 5)  &
          &   .or. (size(fpint_11,1) .ne. nupar)  &
          &   .or. (size(fpint_11,3) .ne. nbessj)) then
          deallocate(fpint_11,fpint_33,fpint_31,pint_11,pint_33,pint_31)
          print*, 'arrays changed size, reallocated'
          allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_31(1:nupar, -2:nkperp_l+2),  &
         &   pint_11(1:nupar, 0:nbessj-1),	&
         &   pint_33(1:nupar, 0:nbessj-1),	&
         &   pint_31(1:nupar),			&
         &   dfdthp(1:nuper,1:nupar))
       end if !what size
    end if  !not allocated
!have b or n changed?  if yes, then l_first must be true
    if(((b_save - bmag)**2/bmag**2 + (n_save - dens)**2/dens**2)  &
       &  .gt. 10d-6) then
       l_first = .true.
    end if  !end if b or n changed
    n_save = dens
    b_save = bmag
    

    ifail=0

    w2=w*w    !omega rf squared
    wp2=dens*zspec**2*wp2fact/aspec !plasma frequency squared
    wpfact=wp2/w2  !multiplying factor for sigma
    wcw=bmag*zspec*EOVERMH/aspec/w  !cyclotron frequency/rf frequency
    beta1=atan2(k1(2),k1(1))  !angle of perp k1
    beta2=atan2(k2(2),k2(1))  !angle of perp k2
    
    kpara1 = k1(3)  ! parallel wave vector
    kpara0 = xkprl0 ! no upshift version
    
    npara1 = kpara1 * c/w  ! parallel index of refraction
    npara0 = kpara0 * c/w  ! no upshift version
    
    mut0=0.5*mpc2*aspec/enorm !a factor for converting
    sqmut0=sqrt(mut0)
    
    dfactper  = npara1 / sqmut0   !puts in the velocity nomalization
    dfactper0 = npara0 / sqmut0   !puts in the velocity nomalization
    
    dfactperi = 1.0 / dfactper
    isq2=sqrt(0.5)
    nharm=nbessj-2  !bessel routine starts with array index 1 == j0 
       !and one order higher than nharm is needed
    kperp2 = k1(1)*k1(1)+ k1(2)*k1(2)
    kperp = sqrt(kperp2)
    r_fact = -dfactperi
    i_fact = -pi * abs(dfactperi)

    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  !used in bessel idendities for pependicular velocity integrals.
    vperp_norm = sqrt(vperp_norm2)  ! need to fix kperp = 0 bess argument normalization need to check.

    sig_fact = 2.0 * pi * wpfact
    kperp_max_l = kperp_max  !keep a local copy that can be changed

!      -------------------------------------------------------
!      do extra integral (stix) to convert all elements to "u"
!      do perpendicular integral by inline trapezoidal rule
!      and parallel integral by simpson's rule
!      -------------------------------------------------------
       

!  independent of interpolation or not, we need the extra factor per
!  stix
    if(l_first .or. .not. l_interp) then
       du = (uper(nuper) - uper(1))/(nuper - 1)
          do j = 1, nupar  !par loop for g_33
             g_33j = 0.0
             do k = 2, nuper - 1  !perp loop for g_33
                g_33j = g_33j +            &
             &     upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
             end do !nperp loop for 33 integral

          k = nuper  !do the last point of trapazoidal integratin
          ga_33(j) = g_33j +                       &
   &          0.5 * upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
       end do  ! end vpar loop all values are filled
       call eqsimpson1d_2(nupar, upar, ga_33, g_33)
       g_33 = g_33 * du
       
    end if !l_first or not l_interp)       

!  the derivative combination in the conductivity can be precomputed
       do j = 1, nupar
          do k = 1 , nuper
             dfdthp(k,j) = upar(j) * dfduper(k,j) - uper(k) * dfdupar(k,j)
	     
	     
!            --------------------------------------------
!            Approximate 2nd term with no upshift, xkprl0
!           ---------------------------------------------	     	     
             dfdthp(k,j) = dfduper(k,j) - dfactper0 * dfdthp(k,j)
	     
	     
!t  his isn't really dfdth now, but a combination that is used
          end do  ! vperp loop
       end do !vpar loop


!      -----------------------------------------------
!      loop over harmonics:  this loop fills integrals
!      -----------------------------------------------

!  logic for dealing with interp/no_interp--
!  for no interp, set kperp to the value we want, and indecies so we execute only once
    if(l_interp) then  !going to interpolate--set up loop
       nk = 0 !starting point for fill loop
       dn = 1 !index increment for fill loop
       dk = kperp_max_l/nkperp_l  !k_perp increment for fill loop
    else
       nk = 1  !start at nk = 1
       dn = 3  !big jump to exit after first loop
       nkperp_l = 1  !one time through loop
       dk = kperp  !k_perp for "exact" calculation
    end if
    if(l_first .or. .not. l_interp) then  !for these cases, we need to do perp integrals
!the loop starts with 0 and increments to nkperp+2 for l_interp, but with 1 for
!one loop for .not. l_interp.
       do  !this loop will either fill the kperp tables, or calculate one value
          if (nk .gt. nkperp_l + 2) exit  !need to fill two past nkperp
             kperp_f = nk*dk
             if(kperp_f .lt. 1.0e-10) then
                kperp_f = 1.0e-8  ! don't let kperp too near zero
             end if
   
   !over all logic comments:
   ! sigma has nine elements, and has a sum over harmonic numbers plus to minus
   ! a max.  when sigma is written for e+ and e-, and has kperp in the x direction,
   ! as is the case here becuase other directions are handled via rotations,
   ! the elements involve the following bessel functions for a harmonic number l
   !              1                        2                        3
   !  1   j(l+1)**2*vperp**2       j(l+1)j(l-1)*vperp**2    j(l+1)j(l)*vperp*vpar
   !  2   j(l+1)j(l-1)*vperp**2    j(l-1)j(l-1)*vperp**2    j(l-1)j(l)*vperp*vpar
   !  1   j(l+1)*j(l)*vperp*vpar   j(l-1)j(l)*vperp*vpar    j(l)j(l)*vpar**2
   !  two types of symmetries and the use of bessel function identities allow us
   !  to reduce the 2*nharm*9 terms to nharm*2
   !  first we note the symmetry in sigma when using +/- coordinates
   !  this reduces the nine elements to six
   !  we can use bessel function identities to get 12, and 13 form 11 and 33
   !  except for nharm = 1
   ! for 12, observe that 2 * l * j(l,z) / z = j(l+1,z) + j(l-1,)
   ! if we square this term (33) and then substract the correct entries for j(l-1), j(l+1)
   ! we can get the cross term.  an extra factor normalization factor must be used to correct
   ! the vperp**2 for the enorm and for kperp that is in the z.
   !  for 13, we do the same thing except use
   ! for harmonics, we note that on changing the sign of l,  l+1=> -(l-1); l=>-l; l-1=>-(l+1)
   ! this allows the following idendifications:
   !  11<=>22 with no sign change--(pairs of bessel functions with the same odd/even
   !  index parity do not change sign on index sign changes
   !  12=21=>12=21 , 33(no changes)
   !  -31=>-32; -32=>32
   !..all symmetries preserved.
   
             iharm = 0  !need 31 for iharm = 0
   !
   ! -- build array with integrand -- !
   ! -- bessel functions are now calculated inside
   
             call wmatprecalc_aorsa(zspec, aspec, enorm, bmag, kperp_f, uper,  &
                 &  nuper, nbessj, 2*nbessj+8,xi1, jnxi1, ifail)
             if(ifail .ne. 0) then
                stop 'bessel function failure in interp intialize'
             end if  !bessel function failure
   
             do j=1,nupar   !start vapr loop for filling iharm = 0
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,1)
                   ssgg_33 = ssgg_33 + lf0 * jnxi1(k,1)
                   lf0 = uper(k) * lf0
                   ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                   ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
                end do !vperp
   
                k = nuper  ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(k,j) * jnxi1(k,1)
                ssgg_33 = ssgg_33 + lf0 * jnxi1(k, 1)
                lf0 = uper(k) * lf0
                ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
   
                ssgg_11 = ssgg_11 * 0.5
                ssgg_31 = ssgg_31 * isq2
   
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
                fpint_31(j,nk) = ssgg_31 * du
             end do !vpar
   !  the 11 and 33 terms are sufficient for all but iharm = 0
   
             do iharm = 1 , nharm + 1
   !do preliminaries for harmonic loop
   
   
   ! -- build array with integrand for the rest of nharms--need one more than nharm!
   
             do j=1,nupar
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
   
   
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,iharm+1)**2
                   ssgg_33 = ssgg_33 + lf0
                   ssgg_11 = ssgg_11 + uper(k) * uper(k) * lf0
                end do !vperp
   
   ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(nuper,j) * jnxi1(nuper,iharm+1)**2
   
                ssgg_33 = ssgg_33 + lf0
   
                ssgg_11 = ssgg_11 + uper(nuper) * uper(nuper) * lf0
                ssgg_11 = ssgg_11 * 0.5
   ! now copy into fipint arrays
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
             end do !vpar
          end do !iharm for generating integrals
   
   
          nk = nk + dn  !i'm doing the do loop indexing myself to control 1st loop
       end do !fill do loop over kperp (one kperp for exact)
    if(l_interp .and. l_first) then  !fill fpint with extra points for interp
!  use symmetries about origin to add extra terms for interpolation
       fpint_11(1:nupar,-2,0:nharm+1) =   fpint_11(1:nupar,2,0:nharm+1)
       fpint_11(1:nupar,-1,0:nharm+1) =   fpint_11(1:nupar,1,0:nharm+1)
       fpint_33(1:nupar,-2,0:nharm) =   fpint_33(1:nupar,2,0:nharm)
       fpint_33(1:nupar,-1,0:nharm) =   fpint_33(1:nupar,1,0:nharm)
       fpint_31(1:nupar,-2) = - fpint_31(1:nupar,2)
       fpint_31(1:nupar,-1) = - fpint_31(1:nupar,1)
    end if !put in end points at kperp = 0
    end if !if we do perp integrals
    sig_11 = cmplx(0.,0.)
    sig_21 = cmplx(0.,0.)
    sig_22 = cmplx(0.,0.)
    sig_31 = cmplx(0.,0.)
    sig_32 = cmplx(0.,0.)
    sig_33 = cmplx(0.,0.)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!  the remaining code either uses the exact integrals, not l_interp, or interpolates in kperp.	
    du_3 = (upar(nupar)-upar(1))/real(nupar-1)/3.
    if(l_interp) then
       dk = kperp_max_l/nkperp_l
       nk = kperp/dk  + 0.5
       p = kperp/dk - nk
! -- build array with interpolation -- !
       if(.true.) then  !fourth order
          a1 = (p*p-1.)*p*(p-2.)/24.
          a2 = -(p-1.)*p*(p*p-4.)/6.
          a3 = (p*p-1.)*(p*p-4.)/4.
          a4 = -(p+1.)*p*(p*p-4.)/6.
          a5 = (p*p-1.)*p*(p+2.)/24.
          pint_11(1:nupar,0:nharm+1) =   &
            &    fpint_11(1:nupar,nk-2,0:nharm+1) * a1 &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2 &
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3 &
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4 &
            & +  fpint_11(1:nupar,nk+2,0:nharm+1) * a5

          pint_33(1:nupar,0:nharm) =   &
            &    fpint_33(1:nupar,nk-2,0:nharm) * a1	&
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2	&
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3     &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4	&
            & +  fpint_33(1:nupar,nk+2,0:nharm) * a5

          pint_31(1:nupar) =   &
            &    fpint_31(1:nupar, nk-2) * a1   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk  ) * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4   &
            & +  fpint_31(1:nupar, nk+2) * a5

       else	! second order
          a2 = p*(p-1.)/2.
          a3 = (1.-p*p)
          a4 = p*(p+1.)/2.
          pint_11(1:nupar,0:iharm+1) =   &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2	&
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3	&
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4

          pint_33(1:nupar,0:iharm+1) =   &
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2   &
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3   &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4

          pint_31(1:nupar) =   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk)   * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4
       end if

    else !integrals are already exactly done for exact--just copy
       pint_11(1:nupar,0:nharm+1) = fpint_11(1:nupar,1,0:nharm+1)
       pint_33(1:nupar,0:nharm) = fpint_33(1:nupar,1,0:nharm)
       pint_31(1:nupar) = fpint_31(1:nupar,1)
    end if

!      the perp integrals have been calculated/interpolated now do
!      parallel, but reuse them for negative harmonics

    do iharm = 0, nharm
!do preliminaries for harmonic loop
       nwcw=real(iharm)*wcw  !l*omega_cyc
       sg_11(1:nupar) = pint_11(1:nupar,iharm + 1)

! do the 31, 32, 21, 22 compression
	
       select case(iharm)
          case(0) !perp blocks are all j1*j1
             sg_22(1:nupar) =   sg_11(1:nupar)  ! no -2 for 0 but they are the same
             sg_21(1:nupar) = - sg_22(1:nupar)     ! 12--22 are anti symmetric
             sg_31(1:nupar) =   pint_31(1:nupar) !this time it has the vpar in it
             sg_32(1:nupar) = - pint_31(1:nupar)
          case default
             sg_22(1:nupar) =   pint_11(1:nupar,iharm - 1)!  look two back for 22
             sg_21(1:nupar) = -0.5*(sg_11(1:nupar)  &
             & +  sg_22(1:nupar))                &
             & +  vperp_norm2*iharm*iharm*pint_33(1:nupar,iharm)  ! use bessel function identity for 21
             sg_32(1:nupar) = sg_31(1:nupar)  !  this is from the previous harmonic
             sg_31(1:nupar) = (sg_11(1:nupar) +  sg_21(1:nupar))*upar(1:nupar)  &
             &  *isq2/iharm/vperp_norm
       end select

       sg_33(1:nupar) = pint_33(1:nupar,iharm) * upar(1:nupar)**2

! now have all the pieces, do the parallel integrals
! -- resonance relation -- !

       rrp=1. - nwcw - dfactper * uparmax
       rrm=1. - nwcw - dfactper * uparmin
	
       if (rrp*rrm.gt.0) then
! -- no resonance here -- !
! -- form the itegrands for the non-singular parallel integral
	
          do j=1,nupar
             rr = 1.- nwcw - dfactper * upar(j)
             irr = 1. / rr
             temp_11(j) = sg_11(j) * irr
             temp_21(j) = sg_21(j) * irr
             temp_22(j) = sg_22(j) * irr
             temp_31(j) = sg_31(j) * irr
             temp_32(j) = sg_32(j) * irr
             temp_33(j) = sg_33(j) * irr
          end do

!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

          e_thtrp_11 = sum(temp_11(2:nupar-1:2))
          e_thtrp_21 = sum(temp_21(2:nupar-1:2))
          e_thtrp_31 = sum(temp_31(2:nupar-1:2))
          e_thtrp_22 = sum(temp_22(2:nupar-1:2))
          e_thtrp_32 = sum(temp_32(2:nupar-1:2))
          e_thtrp_33 = sum(temp_33(2:nupar-1:2))

          o_thtrp_11 = sum(temp_11(3:nupar-1:2))
          o_thtrp_21 = sum(temp_21(3:nupar-1:2))
          o_thtrp_31 = sum(temp_31(3:nupar-1:2))
          o_thtrp_22 = sum(temp_22(3:nupar-1:2))
          o_thtrp_32 = sum(temp_32(3:nupar-1:2))
          o_thtrp_33 = sum(temp_33(3:nupar-1:2))

          thtrp_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrp_11 + 2.*e_thtrp_11))*du_3
          thtrp_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrp_21 + 2.*e_thtrp_21))*du_3
          thtrp_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrp_31 + 2.*e_thtrp_31))*du_3
          thtrp_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrp_22 + 2.*e_thtrp_22))*du_3
          thtrp_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrp_32 + 2.*e_thtrp_32))*du_3
          thtrp_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrp_33 + 2.*e_thtrp_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ --!
	
          sig_11 = sig_11 +  cmplx(thtrp_11,0.0)
          sig_21 = sig_21 +  cmplx(thtrp_21,0.0)
          sig_22 = sig_22 +  cmplx(thtrp_22,0.0)
          sig_31 = sig_31 +  cmplx(thtrp_31,0.0)
          sig_32 = sig_32 +  cmplx(thtrp_32,0.0)
          sig_33 = sig_33 +  cmplx(thtrp_33,0.0)
	
       else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
          upar0 = (1. - nwcw) * dfactperi
	
          if (use_ppart6) then
             fx(1:nupar,1) = sg_11(1:nupar)
             fx(1:nupar,2) = sg_21(1:nupar)
             fx(1:nupar,3) = sg_31(1:nupar)
             fx(1:nupar,4) = sg_22(1:nupar)
             fx(1:nupar,5) = sg_32(1:nupar)
             fx(1:nupar,6) = sg_33(1:nupar)
             nfx = 6
             call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

             thtrp_11 = vint(1)
             thtrp_21 = vint(2)
             thtrp_31 = vint(3)
             thtrp_22 = vint(4)
             thtrp_32 = vint(5)
             thtrp_33 = vint(6)
          else
             call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrp_11)
             call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrp_21)
             call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrp_31)
             call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrp_22)
             call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrp_32)
             call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrp_33)

          endif

! -- 2. anti-hermitian part -- !

          call winterp1d_2(upar,nupar,sg_11,upar0,thtip_11)
          call winterp1d_2(upar,nupar,sg_21,upar0,thtip_21)
          call winterp1d_2(upar,nupar,sg_31,upar0,thtip_31)
          call winterp1d_2(upar,nupar,sg_22,upar0,thtip_22)
          call winterp1d_2(upar,nupar,sg_32,upar0,thtip_32)
          call winterp1d_2(upar,nupar,sg_33,upar0,thtip_33)

          sig_11 = sig_11 +  cmplx(r_fact*thtrp_11,i_fact*thtip_11)
          sig_21 = sig_21 +  cmplx(r_fact*thtrp_21,i_fact*thtip_21)
          sig_22 = sig_22 +  cmplx(r_fact*thtrp_22,i_fact*thtip_22)
          sig_31 = sig_31 +  cmplx(r_fact*thtrp_31,i_fact*thtip_31)
          sig_32 = sig_32 +  cmplx(r_fact*thtrp_32,i_fact*thtip_32)
          sig_33 = sig_33 +  cmplx(r_fact*thtrp_33,i_fact*thtip_33)

       end if

       if(iharm .eq. 0) cycle  !  skip the negative for zero

          nwcw = - nwcw
 ! -- resonance relation -- !

          rrp = 1. - nwcw - dfactper * uparmax
          rrm = 1. - nwcw - dfactper * uparmin
	
          if (rrp*rrm.gt.0) then
! -- no resonance here -- !
	
             do j=1,nupar
                rr = 1.- nwcw - dfactper * upar(j)
                irr = 1. / rr
                temp_11(j) = sg_11(j) * irr
                temp_21(j) = sg_21(j) * irr
                temp_22(j) = sg_22(j) * irr
                temp_31(j) = sg_31(j) * irr
                temp_32(j) = sg_32(j) * irr
                temp_33(j) = sg_33(j) * irr
             end do

	!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

             e_thtrm_11 = sum(temp_11(2:nupar-1:2))
             e_thtrm_21 = sum(temp_21(2:nupar-1:2))
             e_thtrm_31 = sum(temp_31(2:nupar-1:2))
             e_thtrm_22 = sum(temp_22(2:nupar-1:2))
             e_thtrm_32 = sum(temp_32(2:nupar-1:2))
             e_thtrm_33 = sum(temp_33(2:nupar-1:2))
	
             o_thtrm_11 = sum(temp_11(3:nupar-1:2))
             o_thtrm_21 = sum(temp_21(3:nupar-1:2))
             o_thtrm_31 = sum(temp_31(3:nupar-1:2))
             o_thtrm_22 = sum(temp_22(3:nupar-1:2))
             o_thtrm_32 = sum(temp_32(3:nupar-1:2))
             o_thtrm_33 = sum(temp_33(3:nupar-1:2))

             thtrm_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrm_11 + 2.*e_thtrm_11))*du_3
             thtrm_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrm_21 + 2.*e_thtrm_21))*du_3
             thtrm_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrm_31 + 2.*e_thtrm_31))*du_3
             thtrm_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrm_22 + 2.*e_thtrm_22))*du_3
             thtrm_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrm_32 + 2.*e_thtrm_32))*du_3
             thtrm_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrm_33 + 2.*e_thtrm_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ is zero!
!  signs and idex variations in the following reflect the negative harmonic
             sig_11 = sig_11 +  cmplx(thtrm_22,0.0)
             sig_21 = sig_21 +  cmplx(thtrm_21,0.0)
             sig_22 = sig_22 +  cmplx(thtrm_11,0.0)
             sig_31 = sig_31 -  cmplx(thtrm_32,0.0)
             sig_32 = sig_32 -  cmplx(thtrm_31,0.0)
             sig_33 = sig_33 +  cmplx(thtrm_33,0.0)


          else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
	
             upar0 = dfactperi * (1. - nwcw)


             if (use_ppart6) then

                fx(1:nupar,1) = sg_11(1:nupar)
                fx(1:nupar,2) = sg_21(1:nupar)
                fx(1:nupar,3) = sg_31(1:nupar)
                fx(1:nupar,4) = sg_22(1:nupar)
                fx(1:nupar,5) = sg_32(1:nupar)
                fx(1:nupar,6) = sg_33(1:nupar)

                nfx = 6
                call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

                thtrm_11 = vint(1)
                thtrm_21 = vint(2)
                thtrm_31 = vint(3)
                thtrm_22 = vint(4)
                thtrm_32 = vint(5)
                thtrm_33 = vint(6)


             else

                call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrm_11)
                call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrm_21)
                call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrm_31)
                call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrm_22)
                call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrm_32)
                call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrm_33)

             endif

! -- 2. anti-hermitian part -- !
                call winterp1d_2(upar,nupar,sg_11,upar0,thtim_11)
                call winterp1d_2(upar,nupar,sg_21,upar0,thtim_21)
                call winterp1d_2(upar,nupar,sg_31,upar0,thtim_31)
                call winterp1d_2(upar,nupar,sg_22,upar0,thtim_22)
                call winterp1d_2(upar,nupar,sg_32,upar0,thtim_32)
                call winterp1d_2(upar,nupar,sg_33,upar0,thtim_33)


!  signs and idex variations in the following reflect the negative harmonic 
                sig_11 = sig_11 +  cmplx(r_fact*thtrm_22,i_fact*thtim_22)
                sig_21 = sig_21 +  cmplx(r_fact*thtrm_21,i_fact*thtim_21)
                sig_22 = sig_22 +  cmplx(r_fact*thtrm_11,i_fact*thtim_11)
                sig_31 = sig_31 -  cmplx(r_fact*thtrm_32,i_fact*thtim_32)
                sig_32 = sig_32 -  cmplx(r_fact*thtrm_31,i_fact*thtim_31)
                sig_33 = sig_33 +  cmplx(r_fact*thtrm_33,i_fact*thtim_33)

             end if  !resonance exists

         end do !harmonic sum

   ! add in extra term in sig33 from the w to u conversion

         wspec(1,1) = sig_fact *sig_11
         wspec(2,1) = sig_fact *sig_21
         wspec(2,2) = sig_fact *sig_22
         wspec(3,1) = sig_fact *sig_31
         wspec(3,2) = sig_fact *sig_32
         wspec(3,3) = sig_fact * (sig_33 + cmplx(g_33,0.0))
         wspec(1,2) = sig_fact *sig_21
         wspec(1,3) = sig_fact *sig_31
         wspec(2,3) = sig_fact *sig_32

!	 if(l_first) then
!	  print*, 'normalized errors for exact = ', (.not. l_interp), 'calculations'
!	  print*, ' i  j   real    imaginary'
!	  do mm = 1,3
!	     do nn = 1,3
!		  rerror = (real(wspec(mm,nn))- real(wref(mm,nn)))  &
!		   &  /real(wref(mm,nn))
!		   ierror = (imag(wspec(mm,nn))- imag(wref(mm,nn)))  &
!		   &  /imag(wref(mm,nn))
!		   print '(i3, i3, 2e15.3)' , mm, nn, rerror ,ierror
!	     end do
!	  end do
!	 end if

         call wrotate_aorsa(beta1,beta2,wspec)
    return  !finished with orgiginal non iteration version
    end subroutine getnonmaxsigma_aorsa_newi_1


!
!*************************************************************************
!


    subroutine getnonmaxsigma_aorsa_newi_2(w,zspec,aspec,dens,bmag, &
       & k1,xi1,jnxi1,k2,xi2,jnxi2,nbessj,enorm,uparmin,uparmax, &
       & nupar,nuper,uper,upar,dfduper,dfdupar,wspec,ifail, &
	 & l_first, l_interp, kperp_max, nkperp, xkprl0)

!   ---------------------------------------------------------
!   lee's version: interpolates for l_interp = .true
!           does full integrals for l_interp = .false.
!   ---------------------------------------------------------
	
    implicit none


    logical, optional :: l_interp
    logical, optional :: l_first
    integer, optional :: nkperp
    real, optional :: kperp_max
    real :: kperp_max_l
    real :: dk, kperp, p
    integer :: nk, nkperp_l
    integer :: ll, mm, nn
    !nupar, -2 to nkper =2, 0 to nharm=1
    real, save, allocatable, dimension(:,:,:) :: fpint_11, fpint_33
    real, save, allocatable, dimension(:,:) :: fpint_31   ! no nharm
    real, save, allocatable, dimension(:,:) :: pint_11, pint_33
    real, save, allocatable, dimension(:) :: pint_31
    real, intent(in) :: w,zspec,aspec,dens,bmag
    real, dimension(3) :: k1,k2
    integer, intent(in) :: nupar,nuper,nbessj
    real, dimension(nuper) :: xi1,xi2
    real, dimension(nuper, nbessj), intent(in) :: jnxi1,jnxi2
    real, intent(in) :: enorm,uparmin,uparmax
    real, dimension(nuper), intent(in) :: uper
    real, dimension(nupar), intent(in):: upar
    real, dimension(nuper,nupar), intent(in) :: dfduper,dfdupar
    real, save, allocatable, dimension(:,:) ::  dfdthp
    complex, dimension(3,3), intent(inout) :: wspec
    complex :: sig_11, sig_21, sig_22, sig_31, sig_32, sig_33
    complex :: sig_fact
    integer, intent(inout) :: ifail
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: wp2fact=1.745915
    real, parameter:: mpc2=938271998.38
    real, parameter:: c=2.99792458e8
    real, parameter:: pi=3.141592653597932384
    real :: ga_33(nupar)
    real :: w2,wp2,wpfact,wcw,rrp,rrm,rr,irr, du_3, g_33j, xkprl0
    real, save :: g_33, n_save, b_save
    real :: mut0,sqmut0,kpara1, kpara0, npara1, npara0, verp_norm
    real :: isq2,cupar,cuper,nwcw,dfactpar,dfactper, dfactper0, lf0,lnf0uper
    real :: dfactperi, a1, a2, a3, a4, a5,beta1,beta2
    real :: upar0, factor, kperp2, vperp_norm2
    real :: ssgg_11, ssgg_31, ssgg_33, du, temp, kperp_f
    integer :: nharm, iharm, j, k, m, n, dn
    complex :: wref(3,3)
    real :: rerror, ierror, i_fact, r_fact, vperp_norm
    real :: thtrp_11,thtrp_21,thtrp_22,thtrp_31,thtrp_32,thtrp_33
    real :: o_thtrp_11,o_thtrp_21,o_thtrp_22,o_thtrp_31,o_thtrp_32,o_thtrp_33
    real :: e_thtrp_11,e_thtrp_21,e_thtrp_22,e_thtrp_31,e_thtrp_32,e_thtrp_33
    real :: thtip_11,thtip_21,thtip_22,thtip_31,thtip_32,thtip_33
    real :: thtrm_11,thtrm_21,thtrm_22,thtrm_31,thtrm_32,thtrm_33
    real :: o_thtrm_11,o_thtrm_21,o_thtrm_22,o_thtrm_31,o_thtrm_32,o_thtrm_33
    real :: e_thtrm_11,e_thtrm_21,e_thtrm_22,e_thtrm_31,e_thtrm_32,e_thtrm_33
    real :: thtim_11,thtim_21,thtim_22,thtim_31,thtim_32,thtim_33
    real, dimension(nupar) :: temp_11,temp_21,temp_22,temp_31,temp_32,temp_33
    real, dimension(nupar) :: sg_11,sg_21,sg_22,sg_31,sg_32,sg_33




    logical, parameter :: use_ppart6 = .true.
    logical :: is_uniform 
    integer, parameter :: nfxmax = 9
    real*8, dimension(nfxmax) ::  vint
    real*8, dimension(nupar,nfxmax) :: fx
    integer :: nfx
    real*8 :: dx,dh,tol
    integer :: i

!  check to make sure that the v_par mesh is uniform
!  if we need to check par, why not perp as well?

    is_uniform = .true.
    tol = 1.0d-7
    dh = upar(2)-upar(1)
    do i=1,nupar-1
       dx = upar(i+1)-upar(i)
       is_uniform = abs(dx-dh).le. tol*dx
       if (.not.is_uniform) exit
    enddo


!  create a local value of nkperp
    nkperp_l = nkperp !don't change the calling arguement.

    !are the storage arrays allocated?
    !nupar; -2 to nkperp + 2, 0 to nharm + 1 (0 to nbessj -1)
    if(.not. allocated(fpint_11)) then
       if(nkperp_l .lt. 0) then
          nkperp_l = 0  !need to have at least one slot in the nk dimension
       end if !end if set nkperp_1
    
!  allcate storate arrays for interpolation/integration
!  fpint group stores values for interpolation
!  three are needed:  one each for sig_11, sig_f33, and sig_31
!  pint stores the values of the perpendicular integral for the actual k_perp*vperp/omega_c
!  for fpint, the range is:  the parallel mesh; the perpendicular mesh +/- 2 for interpolation;
!  and the harmonic index  (no offset)
!  dfdthp is a common factor for all of sigma that is closely related to dfdtheta

       allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_31(1:nupar, -2:nkperp_l+2),  &
            &   pint_11(1:nupar, 0:nbessj-1),	&
            &   pint_33(1:nupar, 0:nbessj-1),	&
            &   pint_31(1:nupar),               &
            &   dfdthp(1:nuper,1:nupar))
    else  !check for size
       if((size(fpint_11,2) .ne. nkperp_l + 5)  &
          &   .or. (size(fpint_11,1) .ne. nupar)  &
          &   .or. (size(fpint_11,3) .ne. nbessj)) then
          deallocate(fpint_11,fpint_33,fpint_31,pint_11,pint_33,pint_31)
          print*, 'arrays changed size, reallocated'
          allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_31(1:nupar, -2:nkperp_l+2),  &
         &   pint_11(1:nupar, 0:nbessj-1),	&
         &   pint_33(1:nupar, 0:nbessj-1),	&
         &   pint_31(1:nupar),			&
         &   dfdthp(1:nuper,1:nupar))
       end if !what size
    end if  !not allocated
!have b or n changed?  if yes, then l_first must be true
    if(((b_save - bmag)**2/bmag**2 + (n_save - dens)**2/dens**2)  &
       &  .gt. 10d-6) then
       l_first = .true.
    end if  !end if b or n changed
    n_save = dens
    b_save = bmag
    

    ifail=0

    w2=w*w    !omega rf squared
    wp2=dens*zspec**2*wp2fact/aspec !plasma frequency squared
    wpfact=wp2/w2  !multiplying factor for sigma
    wcw=bmag*zspec*EOVERMH/aspec/w  !cyclotron frequency/rf frequency
    beta1=atan2(k1(2),k1(1))  !angle of perp k1
    beta2=atan2(k2(2),k2(1))  !angle of perp k2
    
    kpara1 = k1(3)  ! parallel wave vector
    kpara0 = xkprl0 ! no upshift version
    
    npara1 = kpara1 * c/w  ! parallel index of refraction
    npara0 = kpara0 * c/w  ! no upshift version
    
    mut0=0.5*mpc2*aspec/enorm !a factor for converting
    sqmut0=sqrt(mut0)
    
    dfactper  = npara1 / sqmut0   !puts in the velocity nomalization
    dfactper0 = npara0 / sqmut0   !puts in the velocity nomalization
    
    dfactperi = 1.0 / dfactper
    isq2=sqrt(0.5)
    nharm=nbessj-2  !bessel routine starts with array index 1 == j0 
       !and one order higher than nharm is needed
    kperp2 = k1(1)*k1(1)+ k1(2)*k1(2)
    kperp = sqrt(kperp2)
    r_fact = -dfactperi
    i_fact = -pi * abs(dfactperi)

    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  !used in bessel idendities for pependicular velocity integrals.
    vperp_norm = sqrt(vperp_norm2)  ! need to fix kperp = 0 bess argument normalization need to check.

    sig_fact = 2.0 * pi * wpfact
    kperp_max_l = kperp_max  !keep a local copy that can be changed

!      -------------------------------------------------------
!      do extra integral (stix) to convert all elements to "u"
!      do perpendicular integral by inline trapezoidal rule
!      and parallel integral by simpson's rule
!      -------------------------------------------------------
       

!  independent of interpolation or not, we need the extra factor per
!  stix
    if(l_first .or. .not. l_interp) then
       du = (uper(nuper) - uper(1))/(nuper - 1)
          do j = 1, nupar  !par loop for g_33
             g_33j = 0.0
             do k = 2, nuper - 1  !perp loop for g_33
                g_33j = g_33j +            &
             &     upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
             end do !nperp loop for 33 integral

          k = nuper  !do the last point of trapazoidal integratin
          ga_33(j) = g_33j +                       &
   &          0.5 * upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
       end do  ! end vpar loop all values are filled
       call eqsimpson1d_2(nupar, upar, ga_33, g_33)
       g_33 = g_33 * du
       
    end if !l_first or not l_interp)       

!  the derivative combination in the conductivity can be precomputed
       do j = 1, nupar
          do k = 1 , nuper
             dfdthp(k,j) = upar(j) * dfduper(k,j) - uper(k) * dfdupar(k,j)
	     
	     
!            --------------------------------------------
!            Approximate 2nd term with no upshift, xkprl0
!           ---------------------------------------------	     	     
             dfdthp(k,j) = dfduper(k,j) - dfactper0 * dfdthp(k,j)
	     
	     
!t  his isn't really dfdth now, but a combination that is used
          end do  ! vperp loop
       end do !vpar loop


!      -----------------------------------------------
!      loop over harmonics:  this loop fills integrals
!      -----------------------------------------------

!  logic for dealing with interp/no_interp--
!  for no interp, set kperp to the value we want, and indecies so we execute only once
    if(l_interp) then  !going to interpolate--set up loop
       nk = 0 !starting point for fill loop
       dn = 1 !index increment for fill loop
       dk = kperp_max_l/nkperp_l  !k_perp increment for fill loop
    else
       nk = 1  !start at nk = 1
       dn = 3  !big jump to exit after first loop
       nkperp_l = 1  !one time through loop
       dk = kperp  !k_perp for "exact" calculation
    end if
    if(l_first .or. .not. l_interp) then  !for these cases, we need to do perp integrals
!the loop starts with 0 and increments to nkperp+2 for l_interp, but with 1 for
!one loop for .not. l_interp.
       do  !this loop will either fill the kperp tables, or calculate one value
          if (nk .gt. nkperp_l + 2) exit  !need to fill two past nkperp
             kperp_f = nk*dk
             if(kperp_f .lt. 1.0e-10) then
                kperp_f = 1.0e-8  ! don't let kperp too near zero
             end if
   
   !over all logic comments:
   ! sigma has nine elements, and has a sum over harmonic numbers plus to minus
   ! a max.  when sigma is written for e+ and e-, and has kperp in the x direction,
   ! as is the case here becuase other directions are handled via rotations,
   ! the elements involve the following bessel functions for a harmonic number l
   !              1                        2                        3
   !  1   j(l+1)**2*vperp**2       j(l+1)j(l-1)*vperp**2    j(l+1)j(l)*vperp*vpar
   !  2   j(l+1)j(l-1)*vperp**2    j(l-1)j(l-1)*vperp**2    j(l-1)j(l)*vperp*vpar
   !  1   j(l+1)*j(l)*vperp*vpar   j(l-1)j(l)*vperp*vpar    j(l)j(l)*vpar**2
   !  two types of symmetries and the use of bessel function identities allow us
   !  to reduce the 2*nharm*9 terms to nharm*2
   !  first we note the symmetry in sigma when using +/- coordinates
   !  this reduces the nine elements to six
   !  we can use bessel function identities to get 12, and 13 form 11 and 33
   !  except for nharm = 1
   ! for 12, observe that 2 * l * j(l,z) / z = j(l+1,z) + j(l-1,)
   ! if we square this term (33) and then substract the correct entries for j(l-1), j(l+1)
   ! we can get the cross term.  an extra factor normalization factor must be used to correct
   ! the vperp**2 for the enorm and for kperp that is in the z.
   !  for 13, we do the same thing except use
   ! for harmonics, we note that on changing the sign of l,  l+1=> -(l-1); l=>-l; l-1=>-(l+1)
   ! this allows the following idendifications:
   !  11<=>22 with no sign change--(pairs of bessel functions with the same odd/even
   !  index parity do not change sign on index sign changes
   !  12=21=>12=21 , 33(no changes)
   !  -31=>-32; -32=>32
   !..all symmetries preserved.
   
             iharm = 0  !need 31 for iharm = 0
   !
   ! -- build array with integrand -- !
   ! -- bessel functions are now calculated inside
   
             call wmatprecalc_aorsa(zspec, aspec, enorm, bmag, kperp_f, uper,  &
                 &  nuper, nbessj, 2*nbessj+8,xi1, jnxi1, ifail)
             if(ifail .ne. 0) then
                stop 'bessel function failure in interp intialize'
             end if  !bessel function failure
   
             do j=1,nupar   !start vapr loop for filling iharm = 0
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,1)
                   ssgg_33 = ssgg_33 + lf0 * jnxi1(k,1)
                   lf0 = uper(k) * lf0
                   ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                   ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
                end do !vperp
   
                k = nuper  ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(k,j) * jnxi1(k,1)
                ssgg_33 = ssgg_33 + lf0 * jnxi1(k, 1)
                lf0 = uper(k) * lf0
                ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
   
                ssgg_11 = ssgg_11 * 0.5
                ssgg_31 = ssgg_31 * isq2
   
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
                fpint_31(j,nk) = ssgg_31 * du
             end do !vpar
   !  the 11 and 33 terms are sufficient for all but iharm = 0
   
             do iharm = 1 , nharm + 1
   !do preliminaries for harmonic loop
   
   
   ! -- build array with integrand for the rest of nharms--need one more than nharm!
   
             do j=1,nupar
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
   
   
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,iharm+1)**2
                   ssgg_33 = ssgg_33 + lf0
                   ssgg_11 = ssgg_11 + uper(k) * uper(k) * lf0
                end do !vperp
   
   ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(nuper,j) * jnxi1(nuper,iharm+1)**2
   
                ssgg_33 = ssgg_33 + lf0
   
                ssgg_11 = ssgg_11 + uper(nuper) * uper(nuper) * lf0
                ssgg_11 = ssgg_11 * 0.5
   ! now copy into fipint arrays
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
             end do !vpar
          end do !iharm for generating integrals
   
   
          nk = nk + dn  !i'm doing the do loop indexing myself to control 1st loop
       end do !fill do loop over kperp (one kperp for exact)
    if(l_interp .and. l_first) then  !fill fpint with extra points for interp
!  use symmetries about origin to add extra terms for interpolation
       fpint_11(1:nupar,-2,0:nharm+1) =   fpint_11(1:nupar,2,0:nharm+1)
       fpint_11(1:nupar,-1,0:nharm+1) =   fpint_11(1:nupar,1,0:nharm+1)
       fpint_33(1:nupar,-2,0:nharm) =   fpint_33(1:nupar,2,0:nharm)
       fpint_33(1:nupar,-1,0:nharm) =   fpint_33(1:nupar,1,0:nharm)
       fpint_31(1:nupar,-2) = - fpint_31(1:nupar,2)
       fpint_31(1:nupar,-1) = - fpint_31(1:nupar,1)
    end if !put in end points at kperp = 0
    end if !if we do perp integrals
    sig_11 = cmplx(0.,0.)
    sig_21 = cmplx(0.,0.)
    sig_22 = cmplx(0.,0.)
    sig_31 = cmplx(0.,0.)
    sig_32 = cmplx(0.,0.)
    sig_33 = cmplx(0.,0.)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!  the remaining code either uses the exact integrals, not l_interp, or interpolates in kperp.	
    du_3 = (upar(nupar)-upar(1))/real(nupar-1)/3.
    if(l_interp) then
       dk = kperp_max_l/nkperp_l
       nk = kperp/dk  + 0.5
       p = kperp/dk - nk
! -- build array with interpolation -- !
       if(.true.) then  !fourth order
          a1 = (p*p-1.)*p*(p-2.)/24.
          a2 = -(p-1.)*p*(p*p-4.)/6.
          a3 = (p*p-1.)*(p*p-4.)/4.
          a4 = -(p+1.)*p*(p*p-4.)/6.
          a5 = (p*p-1.)*p*(p+2.)/24.
          pint_11(1:nupar,0:nharm+1) =   &
            &    fpint_11(1:nupar,nk-2,0:nharm+1) * a1 &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2 &
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3 &
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4 &
            & +  fpint_11(1:nupar,nk+2,0:nharm+1) * a5

          pint_33(1:nupar,0:nharm) =   &
            &    fpint_33(1:nupar,nk-2,0:nharm) * a1	&
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2	&
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3     &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4	&
            & +  fpint_33(1:nupar,nk+2,0:nharm) * a5

          pint_31(1:nupar) =   &
            &    fpint_31(1:nupar, nk-2) * a1   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk  ) * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4   &
            & +  fpint_31(1:nupar, nk+2) * a5

       else	! second order
          a2 = p*(p-1.)/2.
          a3 = (1.-p*p)
          a4 = p*(p+1.)/2.
          pint_11(1:nupar,0:iharm+1) =   &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2	&
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3	&
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4

          pint_33(1:nupar,0:iharm+1) =   &
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2   &
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3   &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4

          pint_31(1:nupar) =   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk)   * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4
       end if

    else !integrals are already exactly done for exact--just copy
       pint_11(1:nupar,0:nharm+1) = fpint_11(1:nupar,1,0:nharm+1)
       pint_33(1:nupar,0:nharm) = fpint_33(1:nupar,1,0:nharm)
       pint_31(1:nupar) = fpint_31(1:nupar,1)
    end if

!      the perp integrals have been calculated/interpolated now do
!      parallel, but reuse them for negative harmonics

    do iharm = 0, nharm
!do preliminaries for harmonic loop
       nwcw=real(iharm)*wcw  !l*omega_cyc
       sg_11(1:nupar) = pint_11(1:nupar,iharm + 1)

! do the 31, 32, 21, 22 compression
	
       select case(iharm)
          case(0) !perp blocks are all j1*j1
             sg_22(1:nupar) =   sg_11(1:nupar)  ! no -2 for 0 but they are the same
             sg_21(1:nupar) = - sg_22(1:nupar)     ! 12--22 are anti symmetric
             sg_31(1:nupar) =   pint_31(1:nupar) !this time it has the vpar in it
             sg_32(1:nupar) = - pint_31(1:nupar)
          case default
             sg_22(1:nupar) =   pint_11(1:nupar,iharm - 1)!  look two back for 22
             sg_21(1:nupar) = -0.5*(sg_11(1:nupar)  &
             & +  sg_22(1:nupar))                &
             & +  vperp_norm2*iharm*iharm*pint_33(1:nupar,iharm)  ! use bessel function identity for 21
             sg_32(1:nupar) = sg_31(1:nupar)  !  this is from the previous harmonic
             sg_31(1:nupar) = (sg_11(1:nupar) +  sg_21(1:nupar))*upar(1:nupar)  &
             &  *isq2/iharm/vperp_norm
       end select

       sg_33(1:nupar) = pint_33(1:nupar,iharm) * upar(1:nupar)**2

! now have all the pieces, do the parallel integrals
! -- resonance relation -- !

       rrp=1. - nwcw - dfactper * uparmax
       rrm=1. - nwcw - dfactper * uparmin
	
       if (rrp*rrm.gt.0) then
! -- no resonance here -- !
! -- form the itegrands for the non-singular parallel integral
	
          do j=1,nupar
             rr = 1.- nwcw - dfactper * upar(j)
             irr = 1. / rr
             temp_11(j) = sg_11(j) * irr
             temp_21(j) = sg_21(j) * irr
             temp_22(j) = sg_22(j) * irr
             temp_31(j) = sg_31(j) * irr
             temp_32(j) = sg_32(j) * irr
             temp_33(j) = sg_33(j) * irr
          end do

!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

          e_thtrp_11 = sum(temp_11(2:nupar-1:2))
          e_thtrp_21 = sum(temp_21(2:nupar-1:2))
          e_thtrp_31 = sum(temp_31(2:nupar-1:2))
          e_thtrp_22 = sum(temp_22(2:nupar-1:2))
          e_thtrp_32 = sum(temp_32(2:nupar-1:2))
          e_thtrp_33 = sum(temp_33(2:nupar-1:2))

          o_thtrp_11 = sum(temp_11(3:nupar-1:2))
          o_thtrp_21 = sum(temp_21(3:nupar-1:2))
          o_thtrp_31 = sum(temp_31(3:nupar-1:2))
          o_thtrp_22 = sum(temp_22(3:nupar-1:2))
          o_thtrp_32 = sum(temp_32(3:nupar-1:2))
          o_thtrp_33 = sum(temp_33(3:nupar-1:2))

          thtrp_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrp_11 + 2.*e_thtrp_11))*du_3
          thtrp_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrp_21 + 2.*e_thtrp_21))*du_3
          thtrp_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrp_31 + 2.*e_thtrp_31))*du_3
          thtrp_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrp_22 + 2.*e_thtrp_22))*du_3
          thtrp_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrp_32 + 2.*e_thtrp_32))*du_3
          thtrp_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrp_33 + 2.*e_thtrp_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ --!
	
          sig_11 = sig_11 +  cmplx(thtrp_11,0.0)
          sig_21 = sig_21 +  cmplx(thtrp_21,0.0)
          sig_22 = sig_22 +  cmplx(thtrp_22,0.0)
          sig_31 = sig_31 +  cmplx(thtrp_31,0.0)
          sig_32 = sig_32 +  cmplx(thtrp_32,0.0)
          sig_33 = sig_33 +  cmplx(thtrp_33,0.0)
	
       else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
          upar0 = (1. - nwcw) * dfactperi
	
          if (use_ppart6) then
             fx(1:nupar,1) = sg_11(1:nupar)
             fx(1:nupar,2) = sg_21(1:nupar)
             fx(1:nupar,3) = sg_31(1:nupar)
             fx(1:nupar,4) = sg_22(1:nupar)
             fx(1:nupar,5) = sg_32(1:nupar)
             fx(1:nupar,6) = sg_33(1:nupar)
             nfx = 6
             call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

             thtrp_11 = vint(1)
             thtrp_21 = vint(2)
             thtrp_31 = vint(3)
             thtrp_22 = vint(4)
             thtrp_32 = vint(5)
             thtrp_33 = vint(6)
          else
             call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrp_11)
             call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrp_21)
             call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrp_31)
             call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrp_22)
             call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrp_32)
             call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrp_33)

          endif

! -- 2. anti-hermitian part -- !

          call winterp1d_2(upar,nupar,sg_11,upar0,thtip_11)
          call winterp1d_2(upar,nupar,sg_21,upar0,thtip_21)
          call winterp1d_2(upar,nupar,sg_31,upar0,thtip_31)
          call winterp1d_2(upar,nupar,sg_22,upar0,thtip_22)
          call winterp1d_2(upar,nupar,sg_32,upar0,thtip_32)
          call winterp1d_2(upar,nupar,sg_33,upar0,thtip_33)

          sig_11 = sig_11 +  cmplx(r_fact*thtrp_11,i_fact*thtip_11)
          sig_21 = sig_21 +  cmplx(r_fact*thtrp_21,i_fact*thtip_21)
          sig_22 = sig_22 +  cmplx(r_fact*thtrp_22,i_fact*thtip_22)
          sig_31 = sig_31 +  cmplx(r_fact*thtrp_31,i_fact*thtip_31)
          sig_32 = sig_32 +  cmplx(r_fact*thtrp_32,i_fact*thtip_32)
          sig_33 = sig_33 +  cmplx(r_fact*thtrp_33,i_fact*thtip_33)

       end if

       if(iharm .eq. 0) cycle  !  skip the negative for zero

          nwcw = - nwcw
 ! -- resonance relation -- !

          rrp = 1. - nwcw - dfactper * uparmax
          rrm = 1. - nwcw - dfactper * uparmin
	
          if (rrp*rrm.gt.0) then
! -- no resonance here -- !
	
             do j=1,nupar
                rr = 1.- nwcw - dfactper * upar(j)
                irr = 1. / rr
                temp_11(j) = sg_11(j) * irr
                temp_21(j) = sg_21(j) * irr
                temp_22(j) = sg_22(j) * irr
                temp_31(j) = sg_31(j) * irr
                temp_32(j) = sg_32(j) * irr
                temp_33(j) = sg_33(j) * irr
             end do

	!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

             e_thtrm_11 = sum(temp_11(2:nupar-1:2))
             e_thtrm_21 = sum(temp_21(2:nupar-1:2))
             e_thtrm_31 = sum(temp_31(2:nupar-1:2))
             e_thtrm_22 = sum(temp_22(2:nupar-1:2))
             e_thtrm_32 = sum(temp_32(2:nupar-1:2))
             e_thtrm_33 = sum(temp_33(2:nupar-1:2))
	
             o_thtrm_11 = sum(temp_11(3:nupar-1:2))
             o_thtrm_21 = sum(temp_21(3:nupar-1:2))
             o_thtrm_31 = sum(temp_31(3:nupar-1:2))
             o_thtrm_22 = sum(temp_22(3:nupar-1:2))
             o_thtrm_32 = sum(temp_32(3:nupar-1:2))
             o_thtrm_33 = sum(temp_33(3:nupar-1:2))

             thtrm_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrm_11 + 2.*e_thtrm_11))*du_3
             thtrm_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrm_21 + 2.*e_thtrm_21))*du_3
             thtrm_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrm_31 + 2.*e_thtrm_31))*du_3
             thtrm_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrm_22 + 2.*e_thtrm_22))*du_3
             thtrm_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrm_32 + 2.*e_thtrm_32))*du_3
             thtrm_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrm_33 + 2.*e_thtrm_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ is zero!
!  signs and idex variations in the following reflect the negative harmonic
             sig_11 = sig_11 +  cmplx(thtrm_22,0.0)
             sig_21 = sig_21 +  cmplx(thtrm_21,0.0)
             sig_22 = sig_22 +  cmplx(thtrm_11,0.0)
             sig_31 = sig_31 -  cmplx(thtrm_32,0.0)
             sig_32 = sig_32 -  cmplx(thtrm_31,0.0)
             sig_33 = sig_33 +  cmplx(thtrm_33,0.0)


          else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
	
             upar0 = dfactperi * (1. - nwcw)


             if (use_ppart6) then

                fx(1:nupar,1) = sg_11(1:nupar)
                fx(1:nupar,2) = sg_21(1:nupar)
                fx(1:nupar,3) = sg_31(1:nupar)
                fx(1:nupar,4) = sg_22(1:nupar)
                fx(1:nupar,5) = sg_32(1:nupar)
                fx(1:nupar,6) = sg_33(1:nupar)

                nfx = 6
                call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

                thtrm_11 = vint(1)
                thtrm_21 = vint(2)
                thtrm_31 = vint(3)
                thtrm_22 = vint(4)
                thtrm_32 = vint(5)
                thtrm_33 = vint(6)


             else

                call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrm_11)
                call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrm_21)
                call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrm_31)
                call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrm_22)
                call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrm_32)
                call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrm_33)

             endif

! -- 2. anti-hermitian part -- !
                call winterp1d_2(upar,nupar,sg_11,upar0,thtim_11)
                call winterp1d_2(upar,nupar,sg_21,upar0,thtim_21)
                call winterp1d_2(upar,nupar,sg_31,upar0,thtim_31)
                call winterp1d_2(upar,nupar,sg_22,upar0,thtim_22)
                call winterp1d_2(upar,nupar,sg_32,upar0,thtim_32)
                call winterp1d_2(upar,nupar,sg_33,upar0,thtim_33)


!  signs and idex variations in the following reflect the negative harmonic 
                sig_11 = sig_11 +  cmplx(r_fact*thtrm_22,i_fact*thtim_22)
                sig_21 = sig_21 +  cmplx(r_fact*thtrm_21,i_fact*thtim_21)
                sig_22 = sig_22 +  cmplx(r_fact*thtrm_11,i_fact*thtim_11)
                sig_31 = sig_31 -  cmplx(r_fact*thtrm_32,i_fact*thtim_32)
                sig_32 = sig_32 -  cmplx(r_fact*thtrm_31,i_fact*thtim_31)
                sig_33 = sig_33 +  cmplx(r_fact*thtrm_33,i_fact*thtim_33)

             end if  !resonance exists

         end do !harmonic sum

   ! add in extra term in sig33 from the w to u conversion

         wspec(1,1) = sig_fact *sig_11
         wspec(2,1) = sig_fact *sig_21
         wspec(2,2) = sig_fact *sig_22
         wspec(3,1) = sig_fact *sig_31
         wspec(3,2) = sig_fact *sig_32
         wspec(3,3) = sig_fact * (sig_33 + cmplx(g_33,0.0))
         wspec(1,2) = sig_fact *sig_21
         wspec(1,3) = sig_fact *sig_31
         wspec(2,3) = sig_fact *sig_32

!	 if(l_first) then
!	  print*, 'normalized errors for exact = ', (.not. l_interp), 'calculations'
!	  print*, ' i  j   real    imaginary'
!	  do mm = 1,3
!	     do nn = 1,3
!		  rerror = (real(wspec(mm,nn))- real(wref(mm,nn)))  &
!		   &  /real(wref(mm,nn))
!		   ierror = (imag(wspec(mm,nn))- imag(wref(mm,nn)))  &
!		   &  /imag(wref(mm,nn))
!		   print '(i3, i3, 2e15.3)' , mm, nn, rerror ,ierror
!	     end do
!	  end do
!	 end if

         call wrotate_aorsa(beta1,beta2,wspec)
    return  !finished with orgiginal non iteration version
    end subroutine getnonmaxsigma_aorsa_newi_2


!
!*************************************************************************
!

    subroutine GETNONMAX_SIGMA_AORSA_NEW(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL)
    implicit none
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    real  DFDTH(NUPER,NUPAR)
    complex, dimension(3,3), intent(inout):: WSPEC
    complex zi
    integer, intent(inout):: IFAIL

    complex:: BETAFACT
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPER,3,3):: SSGG
    real, dimension(NUPAR,3,3):: SGG

    real ga_33(nupar)
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_31
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_11

    real, dimension(NUPAR) :: sgg_31_last, sgg_21p
    real, dimension(3,3):: THETARE,THETAIM
    real:: W2,WP2,WPFACT,WCW,RRP,RRM,RR,IRR, du_3, g_33, g_33j
    real:: MUT0,SQMUT0,BETA1,BETA2,KPARA1,NPARA1
    real:: ISQ2,CUPAR,CUPER,NWCW,DFACTPAR,DFACTPER,LF0,LNF0UPER
    real:: DFACTPERI, reson
    real:: UPAR0,SFACT0,SFACTP1,SFACTM1, factor, kperp2, vperp_norm2
    real:: ssgg_11,ssgg_22,ssgg_31,ssgg_33,ssgg_32,ssgg_21,du, temp
    integer:: NHARM,IHARM,NJ,NJP1,NJM1,J,K,M,N

    real, dimension(nupar) :: irr_j, rr_j

    integer :: kdim,  jdim
    real, dimension(nuper,nupar) :: LF0kj
    real :: ssgg_11a,ssgg_33a, ssgg_31a
    real :: maxabserr, maxrelerr
    real :: abserr_11, abserr_33, abserr_31
    real :: relerr_11, relerr_33, relerr_31

    real, dimension(:,:), pointer :: JN0N0,JN0N0U,JN0N1U
    real, dimension(:,:), pointer :: LF0JN0N0, LF0JN0N0U, LF0JN0N1U
    real, dimension(nupar,3*nbessj-1), target :: Jmat
    real, dimension(nuper,3*nbessj-1), target :: Kmat

    integer :: mm,nn,kk,ld1,ld2,ld3
    real :: alpha,beta

    integer :: mn
    integer, dimension(6) :: mlist,nlist
    real, dimension(3,3) :: sumsgg_even, sumsgg_odd
    real :: sum_even, sum_odd, JNP1XI1_K, JN0XI1_K

    logical :: is_uniform 

    integer, parameter :: nfxmax = 9
    integer :: nfx 
    real*8, dimension(nfxmax) ::  vint
    real*8, dimension(nupar,nfxmax) :: fx
    real*8 :: dh,dx,tol
    integer :: i

    is_uniform = .true.
    tol = 1.0d-7
    dh = upar(2)-upar(1)
    do i=1,nupar-1
       dx = upar(i+1)-upar(i)
       is_uniform = abs(dx-dh).le. tol*dx
       if (.not.is_uniform) exit
    enddo

    jdim = nupar
    kdim = nuper

    mn = 1
    do m=1,3
       do n=1,m
          mlist(mn) = m
          nlist(mn) = n
          mn = mn+1
       enddo
    enddo

    zi = cmplx(0.0, 1.0)

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    DFACTPER = NPARA1 / SQMUT0
    DFACTPERI = 1.0 / DFACTPER
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2
    kperp2 = k1(1)*k1(1)+ k1(2)*k1(2)
    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  ! need to fix kperp = 0

!   -------------------------------------------------------
!   do extra integral (STIX) to convert all elements to "U"
!   do perpendicular integral by inline trapezoidal rule
!   and parallel integral by simpson's rule
!   -------------------------------------------------------

    do j = 1, nupar
       do k = 1 , nuper
	  DFDTH(k,j) = upar(j) * DFDUPER(k,j) - uper(k) * DFDUPAR(k,j)
       end do
    end do

    du = (uper(nuper) - uper(1))/(nuper - 1)
    
    do J = 1, NUPAR
       g_33j =  0.0
       
       do K = 2, NUPER - 1
          g_33j = g_33j + DFDTH(k,j)
       end do
       
       k = 1
       g_33j = g_33j + 0.5d0*DFDTH(k,j)
       k = nuper
       g_33j = g_33j + 0.5d0*DFDTH(k,j)

       g_33j = upar(j)*(-g_33j)

       ga_33(j) = g_33j  * du

    enddo


    g_33  = 0.0

!   the uperp integral is done; now do parallel integral by simpson's rule

    call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)

    JN0N0 => Kmat(1:kdim,1:nbessj)
    JN0N0U => Kmat(1:kdim,(nbessj+1):(2*nbessj) )
    JN0N1U => Kmat(1:kdim,(2*nbessj+1):(2*nbessj+nbessj-1) )

    LF0JN0N0 => Jmat(1:jdim,1:nbessj)
    LF0JN0N0U => Jmat(1:jdim,(nbessj+1):(2*nbessj) )
    LF0JN0N1U => Jmat(1:jdim,(2*nbessj+1):(2*nbessj+nbessj-1) )

    do j=1,jdim
       do k=1,kdim
          LF0kj(k,j) = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
       enddo
    enddo

!   -------------------------------
!   precompute LF0kj(k,j) * vectors
!   -------------------------------
    do nj=1,nbessj
       do k=1,kdim
          JN0N0(k,nj) = JNXI1(k,nj)**2
          JN0N0U(k,nj) = (JNXI1(k,nj) * UPER(k))**2
       enddo
    enddo
    
    do nj=1,nbessj-1
       do k=1,kdim
          JN0N1U(k,nj) = JNXI1(k,nj)*JNXI1(k,nj+1)*UPER(k)
       enddo
    enddo

!   -----------------------------
!   one big happy matrix multiply
!   -----------------------------

    mm = jdim
    nn = 3*nbessj-1
    kk = (kdim-1) - 2 + 1
    ld1 = size(LF0kj,1)
    ld2 = size(Kmat,1)
    ld3 = size(Jmat,1)
    alpha = 1.0
    beta = 0.0
    call dgemm( 'T', 'N', mm,nn,kk,                 &
        alpha, LF0kj(2,1), ld1,   Kmat(2,1),ld2,    &
        beta,  Jmat(1,1), ld3 )

    maxabserr = 0.0
    maxrelerr = 0.0

!   -------------------
!   Loop over harmonics
!   -------------------

    do IHARM = -1, NHARM

       NWCW=real(IHARM)*WCW  !L*omega_cyc
       
       reson = 1.0 - nwcw

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)

       ! -- Build array with integrand -- !
	
       do J=1,NUPAR
	
          ssgg_11 = 0.0
          ssgg_31 = 0.0
          ssgg_33 = 0.0

!         --------------------------------------------
!         sums already computed in the matrix multiply
!         --------------------------------------------

          ssgg_33a = LF0JN0N0(j,nj) * (sfact0**2)
          ssgg_11a = LF0JN0N0U(j,njp1) * (sfactp1**2)

!         -------------------------------------
!         note it is possible for njp1 .lt. nj
!         -------------------------------------
          ssgg_31a = LF0JN0N1U(j,min(nj,njp1)) * (sfact0*sfactp1) * upar(j)

          ssgg_33 = ssgg_33a
          ssgg_11 = ssgg_11a
          ssgg_31 = ssgg_31a
	
          ssgg_11 = ssgg_11 * 0.5
	  ssgg_31 = ssgg_31 * ISQ2		

!         ---------------------------------------------------
!         no need to form vectors for JN0XI1(:) or JNP1XI1(:)
!         ---------------------------------------------------

          k = 1

          LF0 = LF0kj(k,j)
          JN0XI1_K = (sfact0*JNXI1(k,nj))
          ssgg_33 = ssgg_33 + 0.5 * JN0XI1_K**2 *LF0

          JNP1XI1_k = (sfactp1*JNXI1(k,njp1))
          LF0 = UPER(K) * JNP1XI1_K * LF0
          ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1_K * LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1_K  * LF0

          k = nuper

          LF0 = LF0kj(k,j)
          JN0XI1_K = (sfact0*JNXI1(k,nj))
          ssgg_33 = ssgg_33 + 0.5 * LF0 * JN0XI1_K**2

          JNP1XI1_K = (sfactp1 * JNXI1(k,njp1))
          LF0 = UPER(K) * JNP1XI1_K  * LF0
          ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1_K * LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1_K  * LF0

!         end of the kperp loop
	
	  du = (uper(nuper) - uper(1))/(nuper - 1)
	  sgg_11(j,iharm) =  ssgg_11 * du
	  sgg_21p(j) =  ssgg_33 * du  ! this is a piece of 21
	  sgg_31(j,iharm) =	ssgg_31 * du
	  SGG(j,3,3) =	UPAR(J)**2 * sgg_21p(j)

       end do
!      -----------------------------------------------
!      The perp integrals have been calculated; now do
!      parallel, but reuse them for negative harmonics
!      ----------------------------------------------
       if(iharm .eq. -1) cycle
	
          du_3 = (UPAR(NUPAR)-UPAR(1))/real(nupar-1)/3.	
	  SGG(1:nupar,1,1) = sgg_11(1:nupar,iharm)
	
!         do the 31--32, 21, 22 compression
	
	  select case(iharm)
	  case(0)
             SGG(1:nupar,2,2) = sgg_11(1:nupar,iharm)  ! no -2 for 0 but they are the same
	     SGG(1:nupar,2,1) = - SGG(1:nupar,2,2)     ! 12--21 are anti symmetric
	  case default
	     SGG(1:nupar,2,2) = sgg_11(1:nupar,iharm-2)!  look two back for 22
             SGG(1:nupar,2,1) = -0.5*(SGG(1:nupar,1,1)+ SGG(1:nupar,2,2))&
	       &+ vperp_norm2*iharm*iharm*(sgg_21p(1:nupar))  ! use bessel function identity for 21
	  end select

	  SGG(1:nupar,3,2) = sgg_31(1:nupar,iharm-1)
	  SGG(1:nupar,3,1) = sgg_31(1:nupar,iharm)
	
!         ------------------------	 	
!         -- Resonance relation -- 
!         ------------------------
          RRP=1. - NWCW - DFACTPER * UPARMAX
          RRM=1. - NWCW - DFACTPER * UPARMIN
		
          if (RRP*RRM.GT.0) then	  
!            -----------------------	  
!            -- No resonance here -- 
!            -----------------------	

             do J=1,NUPAR
                RR = 1.- NWCW - DFACTPER * UPAR(J)
                IRR = 1. / RR
                rr_j(j) = rr
                irr_j(j) = irr
             enddo

             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)

               do j=1,nupar
                  IRR = irr_j(j)
	          SGG(J,M,N) = SGG(J,M,N) * IRR
               end do	    
             end do
	
             !-- Hermitian part of $\Theta/(2\pi)$ --!
	
             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)

                sum_even = sum(sgg(2:(nupar-1):2,m,n))
                sum_odd  = sum(sgg(3:(nupar-1):2,m,n))

                thetare(m,n) = (sgg(1,m,n) + sgg(nupar,m,n) +    &
                   2.0d0*(sum_odd + sum_even) + 2.0d0*sum_even ) * du_3
             end do
	
             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)

                do j=1,nupar
                   SGG(J,M,N) = SGG(J,M,N) * rr_j(j)
                enddo
             end do
	
 	     THETARE(1,2) = THETARE(2,1)
             THETARE(1,3) = THETARE(3,1)
             THETARE(2,3) = THETARE(3,2)

             !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
             THETAIM(1:3,1:3)=0.

          else
!            -------------------------------	  
!            There is a resonance all right 
!            -------------------------------
             ! -- 1. Hermitian part -- !
             UPAR0 = (1. - NWCW) * DFACTPERI	

             fx(1:nupar,1) = SGG(1:nupar,1,1)
             fx(1:nupar,2) = SGG(1:nupar,2,1)
             fx(1:nupar,3) = SGG(1:nupar,3,1)
             fx(1:nupar,4) = SGG(1:nupar,2,2)
             fx(1:nupar,5) = SGG(1:nupar,3,2)
             fx(1:nupar,6) = SGG(1:nupar,3,3)

             nfx = 6
             call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

             thetare(1,1) = vint(1)
             thetare(2,1) = vint(2)
             thetare(3,1) = vint(3)
             thetare(2,2) = vint(4)
             thetare(3,2) = vint(5)
             thetare(3,3) = vint(6)

             THETARE(1,2) = THETARE(2,1)
             THETARE(1,3) = THETARE(3,1)
             THETARE(2,3) = THETARE(3,2)

             ! -- 2. Anti-hermitian part -- !
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
             call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

             THETAIM(1,2) = THETAIM(2,1)
             THETAIM(1,3) = THETAIM(3,1)
             THETAIM(2,3) = THETAIM(3,2)
	

             do M=1,3
                do N=1,3
                   THETARE(M,N) = -THETARE(M,N) * DFACTPERI
                   THETAIM(M,N) = -PI * THETAIM(M,N) * abs(DFACTPERI) 
                end do
             end do

          end if
	  

          do N=1,3
             do M=1,3
                WSPEC(M,N) = WSPEC(M,N) + 2. * PI * WPFACT * BETAFACT * &
!	           & cmplx(THETARE(M,N), THETAIM(M,N))
		   & (THETARE(M,N) + zi * THETAIM(M,N))	
             end do
          end do

	  if(iharm .ne. 0) then  !  skip the negative for zero

	     nwcw = - nwcw

      	     ! -- Resonance relation -- !		
             RRP = 1. - NWCW - DFACTPER * UPARMAX
             RRM = 1. - NWCW - DFACTPER * UPARMIN
	
             if (RRP*RRM.GT.0) then
!            -----------------------	  
!            -- No resonance here -- 
!            -----------------------
	
             do J=1,NUPAR
                RR = 1.- NWCW - DFACTPER * UPAR(J)
                IRR = 1. / RR
                rr_j(j) = rr
                irr_j(j) = irr
             enddo

             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)
                do j=1,nupar
                   IRR = irr_j(j)
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                enddo
             end do
		
!            Hermitian part of $\Theta/(2\pi)$ --

             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)

                sum_even = sum(sgg(2:(nupar-1):2,m,n))
                sum_odd  = sum(sgg(3:(nupar-1):2,m,n))

                thetare(m,n) = (sgg(1,m,n)+sgg(nupar,m,n) +     &
                   2.0d0*(sum_even + sum_odd) + 2.0d0*sum_even) * du_3
             enddo

             do mn=1,6
                m = mlist(mn)
                n = nlist(mn)

                do j=1,nupar
                   SGG(J,M,N) = SGG(J,M,N)  * rr_j(j)
                enddo

             end do

	     temp = THETARE(3,2)
	     THETARE(3,2) = -THETARE(3,1)
	     THETARE(3,1) = -temp
		
	     temp = THETARE(1,1)
	     THETARE(1,1) = THETARE(2,2)
	     THETARE(2,2) = temp

 	     THETARE(1,2) = THETARE(2,1)
      	     THETARE(1,3) = THETARE(3,1)
      	     THETARE(2,3) = THETARE(3,2)

!            Anti-hermitian part of $\Theta/(2\pi)$ 
      	     THETAIM(1:3,1:3)=0.

          else
!            -------------------------------	  
!            There is a resonance all right 
!            -------------------------------
      	     ! -- 1. Hermitian part -- !	
	
      	     UPAR0 = DFACTPERI * (1. - NWCW)

             fx(1:nupar,1) = SGG(1:nupar,1,1)
             fx(1:nupar,2) = SGG(1:nupar,2,1)
             fx(1:nupar,3) = SGG(1:nupar,3,1)
             fx(1:nupar,4) = SGG(1:nupar,2,2)
             fx(1:nupar,5) = SGG(1:nupar,3,2)
             fx(1:nupar,6) = SGG(1:nupar,3,3)

             nfx = 6
             call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

             thetare(1,1) = vint(1)
             thetare(2,1) = vint(2)
             thetare(3,1) = vint(3)
             thetare(2,2) = vint(4)
             thetare(3,2) = vint(5)
             thetare(3,3) = vint(6)

             temp = THETARE(3,2)
	     THETARE(3,2) = -THETARE(3,1)
	     THETARE(3,1) = -temp
	     temp = THETARE(1,1)
	     THETARE(1,1) = THETARE(2,2)
	     THETARE(2,2) = temp
		
      	     THETARE(1,2) = THETARE(2,1)
      	     THETARE(1,3) = THETARE(3,1)
      	     THETARE(2,3) = THETARE(3,2)

      	     ! -- 2. Anti-hermitian part -- !
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
      	     call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

	     temp = THETAIM(3,2)
	     THETAIM(3,2) = -THETAIM(3,1)
	     THETAIM(3,1) = -temp
	     temp = THETAIM(1,1)
	     THETAIM(1,1) = THETAIM(2,2)
	     THETAIM(2,2) = temp
		
      	     THETAIM(1,2) = THETAIM(2,1)
      	     THETAIM(1,3) = THETAIM(3,1)
      	     THETAIM(2,3) = THETAIM(3,2)

             do M=1,3
                do N=1,3
                   THETARE(M,N) = -THETARE(M,N) * DFACTPERI
                   THETAIM(M,N) = -PI * THETAIM(M,N) * abs(DFACTPERI)
                end do
      	     end do
		
          end if  !resonance exists
	
	  do N=1,3
      	     do M=1,3
                WSPEC(M,N) = WSPEC(M,N) + 2. * PI * WPFACT * BETAFACT * &
!                  & cmplx(THETARE(M,N), THETAIM(M,N))
                   & (THETARE(M,N) + zi * THETAIM(M,N))		

      	     end do
          end do

       end if ! skip zero if
       
    end do !harmonic sum

    ! add in extra term in sig33 from the w to u conversion

    WSPEC(3,3)=WSPEC(3,3)+2.*PI*WPFACT*BETAFACT*cmplx(g_33,0.0)

    call WROTATE_AORSA(BETA1,BETA2,WSPEC)

    end subroutine GETNONMAX_SIGMA_AORSA_NEW

!
!*************************************************************************
!


    subroutine getnonmaxsigma_aorsa_newi(w,zspec,aspec,dens,bmag, &
       & k1,xi1,jnxi1,k2,xi2,jnxi2,nbessj,enorm,uparmin,uparmax, &
       & nupar,nuper,uper,upar,dfduper,dfdupar,wspec,ifail, &
	 & l_first, l_interp, kperp_max, nkperp, xkprl0)

!   ---------------------------------------------------------
!   lee's version: interpolates for l_interp = .true
!           does full integrals for l_interp = .false.
!   ---------------------------------------------------------
	
    implicit none


    logical, optional :: l_interp
    logical, optional :: l_first
    integer, optional :: nkperp
    real, optional :: kperp_max
    real :: kperp_max_l
    real :: dk, kperp, p
    integer :: nk, nkperp_l
    integer :: ll, mm, nn
    !nupar, -2 to nkper =2, 0 to nharm=1
    real, save, allocatable, dimension(:,:,:) :: fpint_11, fpint_33
    real, save, allocatable, dimension(:,:) :: fpint_31   ! no nharm
    real, save, allocatable, dimension(:,:) :: pint_11, pint_33
    real, save, allocatable, dimension(:) :: pint_31
    real, intent(in) :: w,zspec,aspec,dens,bmag
    real, dimension(3) :: k1,k2
    integer, intent(in) :: nupar,nuper,nbessj
    real, dimension(nuper) :: xi1,xi2
    real, dimension(nuper, nbessj), intent(in) :: jnxi1,jnxi2
    real, intent(in) :: enorm,uparmin,uparmax
    real, dimension(nuper), intent(in) :: uper
    real, dimension(nupar), intent(in):: upar
    real, dimension(nuper,nupar), intent(in) :: dfduper,dfdupar
    real, save, allocatable, dimension(:,:) ::  dfdthp
    complex, dimension(3,3), intent(inout) :: wspec
    complex :: sig_11, sig_21, sig_22, sig_31, sig_32, sig_33
    complex :: sig_fact
    integer, intent(inout) :: ifail
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: wp2fact=1.745915
    real, parameter:: mpc2=938271998.38
    real, parameter:: c=2.99792458e8
    real, parameter:: pi=3.141592653597932384
    real :: ga_33(nupar)
    real :: w2,wp2,wpfact,wcw,rrp,rrm,rr,irr, du_3, g_33j, xkprl0
    real, save :: g_33, n_save, b_save
    real :: mut0,sqmut0,kpara1, kpara0, npara1, npara0, verp_norm
    real :: isq2,cupar,cuper,nwcw,dfactpar,dfactper, dfactper0, lf0,lnf0uper
    real :: dfactperi, a1, a2, a3, a4, a5,beta1,beta2
    real :: upar0, factor, kperp2, vperp_norm2
    real :: ssgg_11, ssgg_31, ssgg_33, du, temp, kperp_f
    integer :: nharm, iharm, j, k, m, n, dn
    complex :: wref(3,3)
    real :: rerror, ierror, i_fact, r_fact, vperp_norm
    real :: thtrp_11,thtrp_21,thtrp_22,thtrp_31,thtrp_32,thtrp_33
    real :: o_thtrp_11,o_thtrp_21,o_thtrp_22,o_thtrp_31,o_thtrp_32,o_thtrp_33
    real :: e_thtrp_11,e_thtrp_21,e_thtrp_22,e_thtrp_31,e_thtrp_32,e_thtrp_33
    real :: thtip_11,thtip_21,thtip_22,thtip_31,thtip_32,thtip_33
    real :: thtrm_11,thtrm_21,thtrm_22,thtrm_31,thtrm_32,thtrm_33
    real :: o_thtrm_11,o_thtrm_21,o_thtrm_22,o_thtrm_31,o_thtrm_32,o_thtrm_33
    real :: e_thtrm_11,e_thtrm_21,e_thtrm_22,e_thtrm_31,e_thtrm_32,e_thtrm_33
    real :: thtim_11,thtim_21,thtim_22,thtim_31,thtim_32,thtim_33
    real, dimension(nupar) :: temp_11,temp_21,temp_22,temp_31,temp_32,temp_33
    real, dimension(nupar) :: sg_11,sg_21,sg_22,sg_31,sg_32,sg_33




    logical, parameter :: use_ppart6 = .true.
    logical :: is_uniform 
    integer, parameter :: nfxmax = 9
    real*8, dimension(nfxmax) ::  vint
    real*8, dimension(nupar,nfxmax) :: fx
    integer :: nfx
    real*8 :: dx,dh,tol
    integer :: i

!  check to make sure that the v_par mesh is uniform
!  if we need to check par, why not perp as well?

    is_uniform = .true.
    tol = 1.0d-7
    dh = upar(2)-upar(1)
    do i=1,nupar-1
       dx = upar(i+1)-upar(i)
       is_uniform = abs(dx-dh).le. tol*dx
       if (.not.is_uniform) exit
    enddo


!  create a local value of nkperp
    nkperp_l = nkperp !don't change the calling arguement.

    !are the storage arrays allocated?
    !nupar; -2 to nkperp + 2, 0 to nharm + 1 (0 to nbessj -1)
    if(.not. allocated(fpint_11)) then
       if(nkperp_l .lt. 0) then
          nkperp_l = 0  !need to have at least one slot in the nk dimension
       end if !end if set nkperp_1
    
!  allcate storate arrays for interpolation/integration
!  fpint group stores values for interpolation
!  three are needed:  one each for sig_11, sig_f33, and sig_31
!  pint stores the values of the perpendicular integral for the actual k_perp*vperp/omega_c
!  for fpint, the range is:  the parallel mesh; the perpendicular mesh +/- 2 for interpolation;
!  and the harmonic index  (no offset)
!  dfdthp is a common factor for all of sigma that is closely related to dfdtheta

       allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
            &   fpint_31(1:nupar, -2:nkperp_l+2),  &
            &   pint_11(1:nupar, 0:nbessj-1),	&
            &   pint_33(1:nupar, 0:nbessj-1),	&
            &   pint_31(1:nupar),               &
            &   dfdthp(1:nuper,1:nupar))
    else  !check for size
       if((size(fpint_11,2) .ne. nkperp_l + 5)  &
          &   .or. (size(fpint_11,1) .ne. nupar)  &
          &   .or. (size(fpint_11,3) .ne. nbessj)) then
          deallocate(fpint_11,fpint_33,fpint_31,pint_11,pint_33,pint_31)
          print*, 'arrays changed size, reallocated'
          allocate(fpint_11(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_33(1:nupar, -2:nkperp_l+2, 0:nbessj-1),  &
         &   fpint_31(1:nupar, -2:nkperp_l+2),  &
         &   pint_11(1:nupar, 0:nbessj-1),	&
         &   pint_33(1:nupar, 0:nbessj-1),	&
         &   pint_31(1:nupar),			&
         &   dfdthp(1:nuper,1:nupar))
       end if !what size
    end if  !not allocated
!have b or n changed?  if yes, then l_first must be true
    if(((b_save - bmag)**2/bmag**2 + (n_save - dens)**2/dens**2)  &
       &  .gt. 10d-6) then
       l_first = .true.
    end if  !end if b or n changed
    n_save = dens
    b_save = bmag
    

    ifail=0

    w2=w*w    !omega rf squared
    wp2=dens*zspec**2*wp2fact/aspec !plasma frequency squared
    wpfact=wp2/w2  !multiplying factor for sigma
    wcw=bmag*zspec*EOVERMH/aspec/w  !cyclotron frequency/rf frequency
    beta1=atan2(k1(2),k1(1))  !angle of perp k1
    beta2=atan2(k2(2),k2(1))  !angle of perp k2
    
    kpara1 = k1(3)  ! parallel wave vector
    kpara0 = xkprl0 ! no upshift version
    
    npara1 = kpara1 * c/w  ! parallel index of refraction
    npara0 = kpara0 * c/w  ! no upshift version
    
    mut0=0.5*mpc2*aspec/enorm !a factor for converting
    sqmut0=sqrt(mut0)
    
    dfactper  = npara1 / sqmut0   !puts in the velocity nomalization
    dfactper0 = npara0 / sqmut0   !puts in the velocity nomalization
    
    dfactperi = 1.0 / dfactper
    isq2=sqrt(0.5)
    nharm=nbessj-2  !bessel routine starts with array index 1 == j0 
       !and one order higher than nharm is needed
    kperp2 = k1(1)*k1(1)+ k1(2)*k1(2)
    kperp = sqrt(kperp2)
    r_fact = -dfactperi
    i_fact = -pi * abs(dfactperi)

    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  !used in bessel idendities for pependicular velocity integrals.
    vperp_norm = sqrt(vperp_norm2)  ! need to fix kperp = 0 bess argument normalization need to check.

    sig_fact = 2.0 * pi * wpfact
    kperp_max_l = kperp_max  !keep a local copy that can be changed

!      -------------------------------------------------------
!      do extra integral (stix) to convert all elements to "u"
!      do perpendicular integral by inline trapezoidal rule
!      and parallel integral by simpson's rule
!      -------------------------------------------------------
       

!  independent of interpolation or not, we need the extra factor per
!  stix
    if(l_first .or. .not. l_interp) then
       du = (uper(nuper) - uper(1))/(nuper - 1)
          do j = 1, nupar  !par loop for g_33
             g_33j = 0.0
             do k = 2, nuper - 1  !perp loop for g_33
                g_33j = g_33j +            &
             &     upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
             end do !nperp loop for 33 integral

          k = nuper  !do the last point of trapazoidal integratin
          ga_33(j) = g_33j +                       &
   &          0.5 * upar(j)*(uper(k)*dfdupar(k,j)-upar(j)*dfduper(k,j))
       end do  ! end vpar loop all values are filled
       call eqsimpson1d_2(nupar, upar, ga_33, g_33)
       g_33 = g_33 * du
       
    end if !l_first or not l_interp)       

!  the derivative combination in the conductivity can be precomputed
       do j = 1, nupar
          do k = 1 , nuper
             dfdthp(k,j) = upar(j) * dfduper(k,j) - uper(k) * dfdupar(k,j)
	     
	     
!            --------------------------------------------
!            Approximate 2nd term with no upshift, xkprl0
!           ---------------------------------------------	     	     
             dfdthp(k,j) = dfduper(k,j) - dfactper0 * dfdthp(k,j)
	     
	     
!t  his isn't really dfdth now, but a combination that is used
          end do  ! vperp loop
       end do !vpar loop


!      -----------------------------------------------
!      loop over harmonics:  this loop fills integrals
!      -----------------------------------------------

!  logic for dealing with interp/no_interp--
!  for no interp, set kperp to the value we want, and indecies so we execute only once
    if(l_interp) then  !going to interpolate--set up loop
       nk = 0 !starting point for fill loop
       dn = 1 !index increment for fill loop
       dk = kperp_max_l/nkperp_l  !k_perp increment for fill loop
    else
       nk = 1  !start at nk = 1
       dn = 3  !big jump to exit after first loop
       nkperp_l = 1  !one time through loop
       dk = kperp  !k_perp for "exact" calculation
    end if
    if(l_first .or. .not. l_interp) then  !for these cases, we need to do perp integrals
!the loop starts with 0 and increments to nkperp+2 for l_interp, but with 1 for
!one loop for .not. l_interp.
       do  !this loop will either fill the kperp tables, or calculate one value
          if (nk .gt. nkperp_l + 2) exit  !need to fill two past nkperp
             kperp_f = nk*dk
             if(kperp_f .lt. 1.0e-10) then
                kperp_f = 1.0e-8  ! don't let kperp too near zero
             end if
   
   !over all logic comments:
   ! sigma has nine elements, and has a sum over harmonic numbers plus to minus
   ! a max.  when sigma is written for e+ and e-, and has kperp in the x direction,
   ! as is the case here becuase other directions are handled via rotations,
   ! the elements involve the following bessel functions for a harmonic number l
   !              1                        2                        3
   !  1   j(l+1)**2*vperp**2       j(l+1)j(l-1)*vperp**2    j(l+1)j(l)*vperp*vpar
   !  2   j(l+1)j(l-1)*vperp**2    j(l-1)j(l-1)*vperp**2    j(l-1)j(l)*vperp*vpar
   !  1   j(l+1)*j(l)*vperp*vpar   j(l-1)j(l)*vperp*vpar    j(l)j(l)*vpar**2
   !  two types of symmetries and the use of bessel function identities allow us
   !  to reduce the 2*nharm*9 terms to nharm*2
   !  first we note the symmetry in sigma when using +/- coordinates
   !  this reduces the nine elements to six
   !  we can use bessel function identities to get 12, and 13 form 11 and 33
   !  except for nharm = 1
   ! for 12, observe that 2 * l * j(l,z) / z = j(l+1,z) + j(l-1,)
   ! if we square this term (33) and then substract the correct entries for j(l-1), j(l+1)
   ! we can get the cross term.  an extra factor normalization factor must be used to correct
   ! the vperp**2 for the enorm and for kperp that is in the z.
   !  for 13, we do the same thing except use
   ! for harmonics, we note that on changing the sign of l,  l+1=> -(l-1); l=>-l; l-1=>-(l+1)
   ! this allows the following idendifications:
   !  11<=>22 with no sign change--(pairs of bessel functions with the same odd/even
   !  index parity do not change sign on index sign changes
   !  12=21=>12=21 , 33(no changes)
   !  -31=>-32; -32=>32
   !..all symmetries preserved.
   
             iharm = 0  !need 31 for iharm = 0
   !
   ! -- build array with integrand -- !
   ! -- bessel functions are now calculated inside
   
             call wmatprecalc_aorsa(zspec, aspec, enorm, bmag, kperp_f, uper,  &
                 &  nuper, nbessj, 2*nbessj+8,xi1, jnxi1, ifail)
             if(ifail .ne. 0) then
                stop 'bessel function failure in interp intialize'
             end if  !bessel function failure
   
             do j=1,nupar   !start vapr loop for filling iharm = 0
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,1)
                   ssgg_33 = ssgg_33 + lf0 * jnxi1(k,1)
                   lf0 = uper(k) * lf0
                   ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                   ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
                end do !vperp
   
                k = nuper  ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(k,j) * jnxi1(k,1)
                ssgg_33 = ssgg_33 + lf0 * jnxi1(k, 1)
                lf0 = uper(k) * lf0
                ssgg_11 = ssgg_11 + uper(k) * jnxi1(k,1) * lf0
                ssgg_31 = ssgg_31 + upar(j) * jnxi1(k,2) * lf0
   
                ssgg_11 = ssgg_11 * 0.5
                ssgg_31 = ssgg_31 * isq2
   
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
                fpint_31(j,nk) = ssgg_31 * du
             end do !vpar
   !  the 11 and 33 terms are sufficient for all but iharm = 0
   
             do iharm = 1 , nharm + 1
   !do preliminaries for harmonic loop
   
   
   ! -- build array with integrand for the rest of nharms--need one more than nharm!
   
             do j=1,nupar
                ssgg_11 = 0.0
                ssgg_31 = 0.0
                ssgg_33 = 0.0
   
   
                do k = 2, nuper - 1  !do center points for perp integrals
                   lf0 = dfdthp(k,j) * jnxi1(k,iharm+1)**2
                   ssgg_33 = ssgg_33 + lf0
                   ssgg_11 = ssgg_11 + uper(k) * uper(k) * lf0
                end do !vperp
   
   ! do last point for perp integrals--1st point is zero.
   
                lf0 = dfdthp(nuper,j) * jnxi1(nuper,iharm+1)**2
   
                ssgg_33 = ssgg_33 + lf0
   
                ssgg_11 = ssgg_11 + uper(nuper) * uper(nuper) * lf0
                ssgg_11 = ssgg_11 * 0.5
   ! now copy into fipint arrays
                fpint_11(j,nk,iharm) = ssgg_11 * du
                fpint_33(j,nk,iharm) = ssgg_33 * du
             end do !vpar
          end do !iharm for generating integrals
   
   
          nk = nk + dn  !i'm doing the do loop indexing myself to control 1st loop
       end do !fill do loop over kperp (one kperp for exact)
    if(l_interp .and. l_first) then  !fill fpint with extra points for interp
!  use symmetries about origin to add extra terms for interpolation
       fpint_11(1:nupar,-2,0:nharm+1) =   fpint_11(1:nupar,2,0:nharm+1)
       fpint_11(1:nupar,-1,0:nharm+1) =   fpint_11(1:nupar,1,0:nharm+1)
       fpint_33(1:nupar,-2,0:nharm) =   fpint_33(1:nupar,2,0:nharm)
       fpint_33(1:nupar,-1,0:nharm) =   fpint_33(1:nupar,1,0:nharm)
       fpint_31(1:nupar,-2) = - fpint_31(1:nupar,2)
       fpint_31(1:nupar,-1) = - fpint_31(1:nupar,1)
    end if !put in end points at kperp = 0
    end if !if we do perp integrals
    sig_11 = cmplx(0.,0.)
    sig_21 = cmplx(0.,0.)
    sig_22 = cmplx(0.,0.)
    sig_31 = cmplx(0.,0.)
    sig_32 = cmplx(0.,0.)
    sig_33 = cmplx(0.,0.)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!  the remaining code either uses the exact integrals, not l_interp, or interpolates in kperp.	
    du_3 = (upar(nupar)-upar(1))/real(nupar-1)/3.
    if(l_interp) then
       dk = kperp_max_l/nkperp_l
       nk = kperp/dk  + 0.5
       p = kperp/dk - nk
! -- build array with interpolation -- !
       if(.true.) then  !fourth order
          a1 = (p*p-1.)*p*(p-2.)/24.
          a2 = -(p-1.)*p*(p*p-4.)/6.
          a3 = (p*p-1.)*(p*p-4.)/4.
          a4 = -(p+1.)*p*(p*p-4.)/6.
          a5 = (p*p-1.)*p*(p+2.)/24.
          pint_11(1:nupar,0:nharm+1) =   &
            &    fpint_11(1:nupar,nk-2,0:nharm+1) * a1 &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2 &
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3 &
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4 &
            & +  fpint_11(1:nupar,nk+2,0:nharm+1) * a5

          pint_33(1:nupar,0:nharm) =   &
            &    fpint_33(1:nupar,nk-2,0:nharm) * a1	&
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2	&
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3     &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4	&
            & +  fpint_33(1:nupar,nk+2,0:nharm) * a5

          pint_31(1:nupar) =   &
            &    fpint_31(1:nupar, nk-2) * a1   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk  ) * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4   &
            & +  fpint_31(1:nupar, nk+2) * a5

       else	! second order
          a2 = p*(p-1.)/2.
          a3 = (1.-p*p)
          a4 = p*(p+1.)/2.
          pint_11(1:nupar,0:iharm+1) =   &
            & +  fpint_11(1:nupar,nk-1,0:nharm+1) * a2	&
            & +  fpint_11(1:nupar,nk  ,0:nharm+1) * a3	&
            & +  fpint_11(1:nupar,nk+1,0:nharm+1) * a4

          pint_33(1:nupar,0:iharm+1) =   &
            & +  fpint_33(1:nupar,nk-1,0:nharm) * a2   &
            & +  fpint_33(1:nupar,nk  ,0:nharm) * a3   &
            & +  fpint_33(1:nupar,nk+1,0:nharm) * a4

          pint_31(1:nupar) =   &
            & +  fpint_31(1:nupar, nk-1) * a2   &
            & +  fpint_31(1:nupar, nk)   * a3   &
            & +  fpint_31(1:nupar, nk+1) * a4
       end if

    else !integrals are already exactly done for exact--just copy
       pint_11(1:nupar,0:nharm+1) = fpint_11(1:nupar,1,0:nharm+1)
       pint_33(1:nupar,0:nharm) = fpint_33(1:nupar,1,0:nharm)
       pint_31(1:nupar) = fpint_31(1:nupar,1)
    end if

!      the perp integrals have been calculated/interpolated now do
!      parallel, but reuse them for negative harmonics

    do iharm = 0, nharm
!do preliminaries for harmonic loop
       nwcw=real(iharm)*wcw  !l*omega_cyc
       sg_11(1:nupar) = pint_11(1:nupar,iharm + 1)

! do the 31, 32, 21, 22 compression
	
       select case(iharm)
          case(0) !perp blocks are all j1*j1
             sg_22(1:nupar) =   sg_11(1:nupar)  ! no -2 for 0 but they are the same
             sg_21(1:nupar) = - sg_22(1:nupar)     ! 12--22 are anti symmetric
             sg_31(1:nupar) =   pint_31(1:nupar) !this time it has the vpar in it
             sg_32(1:nupar) = - pint_31(1:nupar)
          case default
             sg_22(1:nupar) =   pint_11(1:nupar,iharm - 1)!  look two back for 22
             sg_21(1:nupar) = -0.5*(sg_11(1:nupar)  &
             & +  sg_22(1:nupar))                &
             & +  vperp_norm2*iharm*iharm*pint_33(1:nupar,iharm)  ! use bessel function identity for 21
             sg_32(1:nupar) = sg_31(1:nupar)  !  this is from the previous harmonic
             sg_31(1:nupar) = (sg_11(1:nupar) +  sg_21(1:nupar))*upar(1:nupar)  &
             &  *isq2/iharm/vperp_norm
       end select

       sg_33(1:nupar) = pint_33(1:nupar,iharm) * upar(1:nupar)**2

! now have all the pieces, do the parallel integrals
! -- resonance relation -- !

       rrp=1. - nwcw - dfactper * uparmax
       rrm=1. - nwcw - dfactper * uparmin
	
       if (rrp*rrm.gt.0) then
! -- no resonance here -- !
! -- form the itegrands for the non-singular parallel integral
	
          do j=1,nupar
             rr = 1.- nwcw - dfactper * upar(j)
             irr = 1. / rr
             temp_11(j) = sg_11(j) * irr
             temp_21(j) = sg_21(j) * irr
             temp_22(j) = sg_22(j) * irr
             temp_31(j) = sg_31(j) * irr
             temp_32(j) = sg_32(j) * irr
             temp_33(j) = sg_33(j) * irr
          end do

!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

          e_thtrp_11 = sum(temp_11(2:nupar-1:2))
          e_thtrp_21 = sum(temp_21(2:nupar-1:2))
          e_thtrp_31 = sum(temp_31(2:nupar-1:2))
          e_thtrp_22 = sum(temp_22(2:nupar-1:2))
          e_thtrp_32 = sum(temp_32(2:nupar-1:2))
          e_thtrp_33 = sum(temp_33(2:nupar-1:2))

          o_thtrp_11 = sum(temp_11(3:nupar-1:2))
          o_thtrp_21 = sum(temp_21(3:nupar-1:2))
          o_thtrp_31 = sum(temp_31(3:nupar-1:2))
          o_thtrp_22 = sum(temp_22(3:nupar-1:2))
          o_thtrp_32 = sum(temp_32(3:nupar-1:2))
          o_thtrp_33 = sum(temp_33(3:nupar-1:2))

          thtrp_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrp_11 + 2.*e_thtrp_11))*du_3
          thtrp_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrp_21 + 2.*e_thtrp_21))*du_3
          thtrp_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrp_31 + 2.*e_thtrp_31))*du_3
          thtrp_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrp_22 + 2.*e_thtrp_22))*du_3
          thtrp_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrp_32 + 2.*e_thtrp_32))*du_3
          thtrp_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrp_33 + 2.*e_thtrp_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ --!
	
          sig_11 = sig_11 +  cmplx(thtrp_11,0.0)
          sig_21 = sig_21 +  cmplx(thtrp_21,0.0)
          sig_22 = sig_22 +  cmplx(thtrp_22,0.0)
          sig_31 = sig_31 +  cmplx(thtrp_31,0.0)
          sig_32 = sig_32 +  cmplx(thtrp_32,0.0)
          sig_33 = sig_33 +  cmplx(thtrp_33,0.0)
	
       else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
          upar0 = (1. - nwcw) * dfactperi
	
          if (use_ppart6) then
             fx(1:nupar,1) = sg_11(1:nupar)
             fx(1:nupar,2) = sg_21(1:nupar)
             fx(1:nupar,3) = sg_31(1:nupar)
             fx(1:nupar,4) = sg_22(1:nupar)
             fx(1:nupar,5) = sg_32(1:nupar)
             fx(1:nupar,6) = sg_33(1:nupar)
             nfx = 6
             call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

             thtrp_11 = vint(1)
             thtrp_21 = vint(2)
             thtrp_31 = vint(3)
             thtrp_22 = vint(4)
             thtrp_32 = vint(5)
             thtrp_33 = vint(6)
          else
             call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrp_11)
             call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrp_21)
             call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrp_31)
             call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrp_22)
             call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrp_32)
             call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrp_33)

          endif

! -- 2. anti-hermitian part -- !

          call winterp1d_2(upar,nupar,sg_11,upar0,thtip_11)
          call winterp1d_2(upar,nupar,sg_21,upar0,thtip_21)
          call winterp1d_2(upar,nupar,sg_31,upar0,thtip_31)
          call winterp1d_2(upar,nupar,sg_22,upar0,thtip_22)
          call winterp1d_2(upar,nupar,sg_32,upar0,thtip_32)
          call winterp1d_2(upar,nupar,sg_33,upar0,thtip_33)

          sig_11 = sig_11 +  cmplx(r_fact*thtrp_11,i_fact*thtip_11)
          sig_21 = sig_21 +  cmplx(r_fact*thtrp_21,i_fact*thtip_21)
          sig_22 = sig_22 +  cmplx(r_fact*thtrp_22,i_fact*thtip_22)
          sig_31 = sig_31 +  cmplx(r_fact*thtrp_31,i_fact*thtip_31)
          sig_32 = sig_32 +  cmplx(r_fact*thtrp_32,i_fact*thtip_32)
          sig_33 = sig_33 +  cmplx(r_fact*thtrp_33,i_fact*thtip_33)

       end if

       if(iharm .eq. 0) cycle  !  skip the negative for zero

          nwcw = - nwcw
 ! -- resonance relation -- !

          rrp = 1. - nwcw - dfactper * uparmax
          rrm = 1. - nwcw - dfactper * uparmin
	
          if (rrp*rrm.gt.0) then
! -- no resonance here -- !
	
             do j=1,nupar
                rr = 1.- nwcw - dfactper * upar(j)
                irr = 1. / rr
                temp_11(j) = sg_11(j) * irr
                temp_21(j) = sg_21(j) * irr
                temp_22(j) = sg_22(j) * irr
                temp_31(j) = sg_31(j) * irr
                temp_32(j) = sg_32(j) * irr
                temp_33(j) = sg_33(j) * irr
             end do

	!-- hermitian part of $\theta/(2\pi)$ --  use trap. integration!

             e_thtrm_11 = sum(temp_11(2:nupar-1:2))
             e_thtrm_21 = sum(temp_21(2:nupar-1:2))
             e_thtrm_31 = sum(temp_31(2:nupar-1:2))
             e_thtrm_22 = sum(temp_22(2:nupar-1:2))
             e_thtrm_32 = sum(temp_32(2:nupar-1:2))
             e_thtrm_33 = sum(temp_33(2:nupar-1:2))
	
             o_thtrm_11 = sum(temp_11(3:nupar-1:2))
             o_thtrm_21 = sum(temp_21(3:nupar-1:2))
             o_thtrm_31 = sum(temp_31(3:nupar-1:2))
             o_thtrm_22 = sum(temp_22(3:nupar-1:2))
             o_thtrm_32 = sum(temp_32(3:nupar-1:2))
             o_thtrm_33 = sum(temp_33(3:nupar-1:2))

             thtrm_11 = (temp_11(1)+temp_11(nupar)+2.*(o_thtrm_11 + 2.*e_thtrm_11))*du_3
             thtrm_21 = (temp_21(1)+temp_21(nupar)+2.*(o_thtrm_21 + 2.*e_thtrm_21))*du_3
             thtrm_31 = (temp_31(1)+temp_31(nupar)+2.*(o_thtrm_31 + 2.*e_thtrm_31))*du_3
             thtrm_22 = (temp_22(1)+temp_22(nupar)+2.*(o_thtrm_22 + 2.*e_thtrm_22))*du_3
             thtrm_32 = (temp_32(1)+temp_32(nupar)+2.*(o_thtrm_32 + 2.*e_thtrm_32))*du_3
             thtrm_33 = (temp_33(1)+temp_33(nupar)+2.*(o_thtrm_33 + 2.*e_thtrm_33))*du_3

!-- anti-hermitian part of $\theta/(2\pi)$ is zero!
!  signs and idex variations in the following reflect the negative harmonic
             sig_11 = sig_11 +  cmplx(thtrm_22,0.0)
             sig_21 = sig_21 +  cmplx(thtrm_21,0.0)
             sig_22 = sig_22 +  cmplx(thtrm_11,0.0)
             sig_31 = sig_31 -  cmplx(thtrm_32,0.0)
             sig_32 = sig_32 -  cmplx(thtrm_31,0.0)
             sig_33 = sig_33 +  cmplx(thtrm_33,0.0)


          else

! -- there is a resonance all right -- !
! -- 1. hermitian part -- !
	
             upar0 = dfactperi * (1. - nwcw)


             if (use_ppart6) then

                fx(1:nupar,1) = sg_11(1:nupar)
                fx(1:nupar,2) = sg_21(1:nupar)
                fx(1:nupar,3) = sg_31(1:nupar)
                fx(1:nupar,4) = sg_22(1:nupar)
                fx(1:nupar,5) = sg_32(1:nupar)
                fx(1:nupar,6) = sg_33(1:nupar)

                nfx = 6
                call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint,is_uniform)

                thtrm_11 = vint(1)
                thtrm_21 = vint(2)
                thtrm_31 = vint(3)
                thtrm_22 = vint(4)
                thtrm_32 = vint(5)
                thtrm_33 = vint(6)


             else

                call cauchy_ppart2(upar,nupar,upar0,sg_11,thtrm_11)
                call cauchy_ppart2(upar,nupar,upar0,sg_21,thtrm_21)
                call cauchy_ppart2(upar,nupar,upar0,sg_31,thtrm_31)
                call cauchy_ppart2(upar,nupar,upar0,sg_22,thtrm_22)
                call cauchy_ppart2(upar,nupar,upar0,sg_32,thtrm_32)
                call cauchy_ppart2(upar,nupar,upar0,sg_33,thtrm_33)

             endif

! -- 2. anti-hermitian part -- !
                call winterp1d_2(upar,nupar,sg_11,upar0,thtim_11)
                call winterp1d_2(upar,nupar,sg_21,upar0,thtim_21)
                call winterp1d_2(upar,nupar,sg_31,upar0,thtim_31)
                call winterp1d_2(upar,nupar,sg_22,upar0,thtim_22)
                call winterp1d_2(upar,nupar,sg_32,upar0,thtim_32)
                call winterp1d_2(upar,nupar,sg_33,upar0,thtim_33)


!  signs and idex variations in the following reflect the negative harmonic 
                sig_11 = sig_11 +  cmplx(r_fact*thtrm_22,i_fact*thtim_22)
                sig_21 = sig_21 +  cmplx(r_fact*thtrm_21,i_fact*thtim_21)
                sig_22 = sig_22 +  cmplx(r_fact*thtrm_11,i_fact*thtim_11)
                sig_31 = sig_31 -  cmplx(r_fact*thtrm_32,i_fact*thtim_32)
                sig_32 = sig_32 -  cmplx(r_fact*thtrm_31,i_fact*thtim_31)
                sig_33 = sig_33 +  cmplx(r_fact*thtrm_33,i_fact*thtim_33)

             end if  !resonance exists

         end do !harmonic sum

   ! add in extra term in sig33 from the w to u conversion

         wspec(1,1) = sig_fact *sig_11
         wspec(2,1) = sig_fact *sig_21
         wspec(2,2) = sig_fact *sig_22
         wspec(3,1) = sig_fact *sig_31
         wspec(3,2) = sig_fact *sig_32
         wspec(3,3) = sig_fact * (sig_33 + cmplx(g_33,0.0))
         wspec(1,2) = sig_fact *sig_21
         wspec(1,3) = sig_fact *sig_31
         wspec(2,3) = sig_fact *sig_32

!	 if(l_first) then
!	  print*, 'normalized errors for exact = ', (.not. l_interp), 'calculations'
!	  print*, ' i  j   real    imaginary'
!	  do mm = 1,3
!	     do nn = 1,3
!		  rerror = (real(wspec(mm,nn))- real(wref(mm,nn)))  &
!		   &  /real(wref(mm,nn))
!		   ierror = (imag(wspec(mm,nn))- imag(wref(mm,nn)))  &
!		   &  /imag(wref(mm,nn))
!		   print '(i3, i3, 2e15.3)' , mm, nn, rerror ,ierror
!	     end do
!	  end do
!	 end if

         call wrotate_aorsa(beta1,beta2,wspec)
    return  !finished with orgiginal non iteration version
    end subroutine getnonmaxsigma_aorsa_newi

!
!*************************************************************************
!

  subroutine GETNONMAX_SIGMA_AORSA_NEW1(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL)
    implicit none
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    real  DFDTH(NUPER,NUPAR)
    complex, dimension(3,3), intent(OUT):: WSPEC
    integer, intent(OUT):: IFAIL

    complex:: BETAFACT

    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
	
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPER,3,3):: SSGG
    real, dimension(NUPAR,3,3):: SGG

    real ga_33(nupar)
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_31
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_11

    real, dimension(NUPAR) :: sgg_31_last, sgg_21p
    real, dimension(3,3):: THETARE,THETAIM
    real:: W2,WP2,WPFACT,WCW,RRP,RRM,RR,IRR, du_3, g_33, g_33j
    real:: MUT0,SQMUT0,BETA1,BETA2,KPARA1,NPARA1
    real:: ISQ2,CUPAR,CUPER,NWCW,DFACTPAR,DFACTPER,LF0,LNF0UPER
    real:: DFACTPERI
    real:: UPAR0,SFACT0,SFACTP1,SFACTM1, factor, kperp2, vperp_norm2
    real:: ssgg_11,ssgg_22,ssgg_31,ssgg_33,ssgg_32,ssgg_21,du, temp
    integer:: NHARM,IHARM,NJ,NJP1,NJM1,J,K,M,N

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    DFACTPER = NPARA1 / SQMUT0
    DFACTPERI = 1.0 / DFACTPER
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2
    kperp2 = k1(1)*k1(1)+ k1(2)*k1(2)
    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  ! need to fix kperp = 0
    
    


!      -------------------------------------------------------
!      do extra integral (STIX) to convert all elements to "U"
!      do perpendicular integral by inline trapezoidal rule
!      and parallel integral by simpson's rule
!      -------------------------------------------------------

       do J = 1, NUPAR
	  g_33j = 0.0

          do K = 2, NUPER - 1
             g_33j = g_33j + UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
          end do


	  k = 1
          g_33j = g_33j + 0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))


	  k = nuper
          g_33j = g_33j + 0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))

	
	  du = (uper(nuper) - uper(1))/(nuper - 1)
	  ga_33(j) = g_33j * du

       end do
       
       g_33  = 0.0

!      the uperp integral is done; now do parallel integral by simpson's rule

       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)
       

       do k = 1 , nuper
          do j = 1, nupar
	     DFDTH(k,j) = upar(j) * DFDUPER(k,j) - uper(k) * DFDUPAR(k,j)
	  end do
       end do





!      -------------------
!      Loop over harmonics
!      -------------------

    do IHARM = -1, NHARM
       NWCW=real(IHARM)*WCW  !L*omega_cyc

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)

       ! -- Build array with integrand -- !
       

	  
       do J=1,NUPAR
	  
	  ssgg_11 = 0.0
	  ssgg_31 = 0.0
	  ssgg_33 = 0.0

	  do K = 2, NUPER - 1
	    
             LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
             ssgg_33 = ssgg_33 + JN0XI1(K)**2 *LF0
		 
             LF0 = UPER(K) * JNP1XI1(K) * LF0
             ssgg_11 = ssgg_11 + UPER(K) * JNP1XI1(K)* LF0
             ssgg_31 = ssgg_31 + UPAR(J) * JN0XI1(K) * LF0
          end do
	    
	  ssgg_11 = ssgg_11 * 0.5
	  ssgg_31 = ssgg_31 * ISQ2
	    
	    
	  k = 1

	  LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
	  ssgg_33 = ssgg_33 + 0.5 * JN0XI1(K)**2 *LF0 
	  
          LF0 = UPER(K) * JNP1XI1(K) * LF0
          ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0  
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0


	  k = nuper
	  
	  LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
	  ssgg_33 = ssgg_33 + 0.5 * LF0 * JN0XI1(K)**2
	  
	  LF0 = UPER(K) * JNP1XI1(K) * LF0
	  ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0

       !  end of the kperp loop
       
       
	
	  du = (uper(nuper) - uper(1))/(nuper - 1)
	  sgg_11(j,iharm) =  ssgg_11 * du
	  sgg_21p(j) =  ssgg_33 * du  ! this is a piece of 21
	  sgg_31(j,iharm) =	ssgg_31 * du
	  SGG(j,3,3) =	UPAR(J)**2 * sgg_21p(j)

       end do
       
       
!      the perp integrals have been calculated now do
!      parallel, but reuse them for negative harmonics

       if(iharm .eq. -1) cycle
	 
       du_3 = (UPAR(NUPAR)-UPAR(1))/real(nupar-1)/3.	
	 SGG(1:nupar,1,1) = sgg_11(1:nupar,iharm)
	 
	 ! do the 31--32, 21, 22 compression
	 
	 select case(iharm)  
	 case(0)
            SGG(1:nupar,2,2) = sgg_11(1:nupar,iharm)  ! no -2 for 0 but they are the same
	    SGG(1:nupar,2,1) = - SGG(1:nupar,2,2)     ! 12--21 are anti symmetric
	 case default
	    SGG(1:nupar,2,2) = sgg_11(1:nupar,iharm-2)!  look two back for 22
            SGG(1:nupar,2,1) = -0.5*(SGG(1:nupar,1,1)+ SGG(1:nupar,2,2))&
	      &+ vperp_norm2*iharm*iharm*(sgg_21p(1:nupar))  ! use bessel function identity for 21
	 end select

	 SGG(1:nupar,3,2) = sgg_31(1:nupar,iharm-1)
	 SGG(1:nupar,3,1) = sgg_31(1:nupar,iharm)
	 
	 	
       ! -- Resonance relation -- !
       
       RRP=1. - NWCW - DFACTPER * UPARMAX 
       RRM=1. - NWCW - DFACTPER * UPARMIN 
	 
	 

       if (RRP*RRM.GT.0) then
          ! -- No resonance here -- !
	  
          do J=1,NUPAR
             RR = 1.- NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                end do
             end do
          end do
	  

          !-- Hermitian part of $\Theta/(2\pi)$ --!
	    
	  THETARE(1,1) = (sgg(1,1,1)+sgg(nupar,1,1)+2.*sum(sgg(2:nupar-1,1,1)) +&
	       & 2.*sum(sgg(2:nupar-1:2,1,1)))*du_3
	  THETARE(2,1) = (SGG(1,2,1)+SGG(nupar,2,1)+2.*sum(SGG(2:nupar-1,2,1)) +&
	       & 2.*sum(SGG(2:nupar-1:2,2,1)))*du_3
	  THETARE(3,1) = (sgg(1,3,1)+sgg(nupar,3,1)+2.*sum(sgg(2:nupar-1,3,1)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,1)))*du_3
	  THETARE(2,2) = (sgg(1,2,2)+sgg(nupar,2,2)+2.*sum(sgg(2:nupar-1,2,2)) +&
	       & 2.*sum(sgg(2:nupar-1:2,2,2)))*du_3
          THETARE(3,2) = (sgg(1,3,2)+sgg(nupar,3,2)+2.*sum(sgg(2:nupar-1,3,2)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,2)))*du_3
	  THETARE(3,3) = (sgg(1,3,3)+sgg(nupar,3,3)+2.*sum(sgg(2:nupar-1,3,3)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,3)))*du_3
	
	
	  do J=1,NUPAR ! restore sgg
	     RR = 1.- NWCW - DFACTPER * UPAR(J)
	     
             do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * RR
                end do
             end do
	     
          end do
	    
 	  THETARE(1,2) = THETARE(2,1)
          THETARE(1,3) = THETARE(3,1)
          THETARE(2,3) = THETARE(3,2)

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          ! -- There is a resonance all right -- !
          ! -- 1. Hermitian part -- !
          UPAR0 = (1. - NWCW) * DFACTPERI
	  
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))

          THETARE(1,2) = THETARE(2,1)
          THETARE(1,3) = THETARE(3,1)
          THETARE(2,3) = THETARE(3,2)

          ! -- 2. Anti-hermitian part -- !
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
          call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

          THETAIM(1,2) = THETAIM(2,1)
          THETAIM(1,3) = THETAIM(3,1)
          THETAIM(2,3) = THETAIM(3,2)
	  

          do M=1,3
             do N=1,3
                THETARE(M,N) = -THETARE(M,N) * DFACTPERI
                THETAIM(M,N) = -PI * THETAIM(M,N) * abs(DFACTPERI)
             end do
          end do

       end if

       do N=1,3
          do M=1,3
             WSPEC(M,N)=WSPEC(M,N)+2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
          end do
       end do

	 if(iharm .ne. 0) then  !  skip the negative for zero

	    nwcw = - nwcw

      	  ! -- Resonance relation -- !
	  
	 
          RRP = 1. - NWCW - DFACTPER * UPARMAX 
          RRM = 1. - NWCW - DFACTPER * UPARMIN
	    
          if (RRP*RRM.GT.0) then
      	 ! -- No resonance here -- !
	 
      	  do J=1,NUPAR

             RR = 1. - NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             do M=1,3
                do N = 1,M
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                end do
             end do
      	  end do
		
          !-- Hermitian part of $\Theta/(2\pi)$ --!
	  THETARE(1,1) = (sgg(1,1,1)+sgg(nupar,1,1)+2.*sum(sgg(2:nupar-1,1,1)) +&
	       & 2.*sum(sgg(2:nupar-1:2,1,1)))*du_3
	  THETARE(2,1) = (SGG(1,2,1)+SGG(nupar,2,1)+2.*sum(SGG(2:nupar-1,2,1)) +&
	       & 2.*sum(SGG(2:nupar-1:2,2,1)))*du_3
	  THETARE(3,1) = (sgg(1,3,1)+sgg(nupar,3,1)+2.*sum(sgg(2:nupar-1,3,1)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,1)))*du_3
	  THETARE(2,2) = (sgg(1,2,2)+sgg(nupar,2,2)+2.*sum(sgg(2:nupar-1,2,2)) +&
	       & 2.*sum(sgg(2:nupar-1:2,2,2)))*du_3
          THETARE(3,2) = (sgg(1,3,2)+sgg(nupar,3,2)+2.*sum(sgg(2:nupar-1,3,2)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,2)))*du_3
	  THETARE(3,3) = (sgg(1,3,3)+sgg(nupar,3,3)+2.*sum(sgg(2:nupar-1,3,3)) +&
	       & 2.*sum(sgg(2:nupar-1:2,3,3)))*du_3
	       

		
          do J=1,NUPAR  ! restore sgg
             RR=1.- NWCW - DFACTPER * UPAR(J)
	     do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * RR
            	end do
             end do
      	  end do

	  temp = THETARE(3,2)
	  THETARE(3,2) = -THETARE(3,1)
	  THETARE(3,1) = -temp
		
	  temp = THETARE(1,1)
	  THETARE(1,1) = THETARE(2,2)
	  THETARE(2,2) = temp

 	  THETARE(1,2) = THETARE(2,1)
      	  THETARE(1,3) = THETARE(3,1)
      	  THETARE(2,3) = THETARE(3,2)

      	 !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
      	 THETAIM(1:3,1:3)=0.

          else

      	 ! -- There is a resonance all right -- !
      	 ! -- 1. Hermitian part -- !
	 
	 
      	 UPAR0 = DFACTPERI * (1. - NWCW)

      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))

		 temp = THETARE(3,2)
		 THETARE(3,2) = -THETARE(3,1)
		 THETARE(3,1) = -temp
		 temp = THETARE(1,1)
		 THETARE(1,1) = THETARE(2,2)
		 THETARE(2,2) = temp
		
      	 THETARE(1,2) = THETARE(2,1)
      	 THETARE(1,3) = THETARE(3,1)
      	 THETARE(2,3) = THETARE(3,2)

      	 ! -- 2. Anti-hermitian part -- !
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
      	 call WINTERP1D_2(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

	 temp = THETAIM(3,2)
	 THETAIM(3,2) = -THETAIM(3,1)
	 THETAIM(3,1) = -temp
	 temp = THETAIM(1,1)
	 THETAIM(1,1) = THETAIM(2,2)
	 THETAIM(2,2) = temp
		
      	 THETAIM(1,2) = THETAIM(2,1)
      	 THETAIM(1,3) = THETAIM(3,1)
      	 THETAIM(2,3) = THETAIM(3,2)

         do M=1,3
            do N=1,3
               THETARE(M,N) = -THETARE(M,N) * DFACTPERI
               THETAIM(M,N) = -PI * THETAIM(M,N) * abs(DFACTPERI)
            end do
      	 end do
		

          end if  !resonance exists
	   
	    do N=1,3
      	 do M=1,3
                WSPEC(M,N)=WSPEC(M,N)+2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
      	 end do
          end do

	 end if ! zero skip if
    end do !harmonic sum

    ! add in extra term in sig33 from the w to u conversion

    WSPEC(3,3)=WSPEC(3,3)+2.*PI*WPFACT*BETAFACT*cmplx(g_33,0.0)
    
    call WROTATE_AORSA(BETA1,BETA2,WSPEC)
    
  end subroutine GETNONMAX_SIGMA_AORSA_NEW1

!
!*************************************************************************
!

  subroutine GETNONMAX_SIGMA_AORSA_NEW0(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL)
    implicit none
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    complex, dimension(3,3), intent(OUT):: WSPEC
    integer, intent(OUT):: IFAIL

    complex:: BETAFACT
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPER,3,3):: SSGG
    real, dimension(NUPAR,3,3):: SGG
    real, dimension(3,3):: THETARE,THETAIM
    real:: W2,WP2,WPFACT,WCW,RRP,RRM,RR,IRR
    real:: MUT0,SQMUT0,BETA1,BETA2,KPARA1,NPARA1
    real:: ISQ2,CUPAR,CUPER,NWCW,DFACTPAR,DFACTPER,LF0,LNF0UPER
    real:: UPAR0,SFACT0,SFACTP1,SFACTM1, factor
    integer:: NHARM,IHARM,NJ,NJP1,NJM1,J,K,M,N

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2

    ! -- Loop over harmonics -- !
    do IHARM=-NHARM,NHARM
       NWCW=real(IHARM)*WCW

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)
       JN0XI2(1:NUPER)=SFACT0*JNXI2(1:NUPER,NJ)
       JNP1XI2(1:NUPER)=SFACTP1*JNXI2(1:NUPER,NJP1)
       JNM1XI2(1:NUPER)=SFACTM1*JNXI2(1:NUPER,NJM1)

       ! -- Build array with integrand -- !
       do J=1,NUPAR
          CUPAR=UPAR(J)
          DFACTPAR=NPARA1*CUPAR/SQMUT0
          do K=1,NUPER
             CUPER=UPER(K)
             DFACTPER=NPARA1*CUPER/SQMUT0
             LF0=(1.-DFACTPAR)*DFDUPER(K,J)+DFACTPER*DFDUPAR(K,J)
             LNF0UPER=NWCW*DFDUPER(K,J)*CUPAR+(1.-NWCW)*DFDUPAR(K,J)*CUPER


             SSGG(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
             SSGG(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
             SSGG(K,3,1)=ISQ2*CUPAR*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
!             SSGG(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
             SSGG(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
             SSGG(K,3,2)=ISQ2*CUPAR*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
!             SSGG(K,1,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
!             SSGG(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
             SSGG(K,3,3)=CUPAR*LNF0UPER*JN0XI1(K)*JN0XI2(K)



          end do

          ! -- Integrate of $u_\perp$ -- !
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,1),SGG(J,1,1))
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,1),SGG(J,2,1))
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,1),SGG(J,3,1))
!!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,2),SGG(J,1,2))
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,2),SGG(J,2,2))
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,2),SGG(J,3,2))
!!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,3),SGG(J,1,3))
!!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,3),SGG(J,2,3))
!          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,3),SGG(J,3,3))

          call EQTRAPZ1D(NUPER,UPER,SSGG(1,1,1),SGG(J,1,1))
          call EQTRAPZ1D(NUPER,UPER,SSGG(1,2,1),SGG(J,2,1))
          call EQTRAPZ1D(NUPER,UPER,SSGG(1,3,1),SGG(J,3,1))
          call EQTRAPZ1D(NUPER,UPER,SSGG(1,2,2),SGG(J,2,2))
          call EQTRAPZ1D(NUPER,UPER,SSGG(1,3,2),SGG(J,3,2))
          call EQTRAPZ1D(NUPER,UPER,SSGG(1,3,3),SGG(J,3,3))


!          call sgrate2(NUPER,UPER,SSGG(1,1,1),SGG(J,1,1))
!          call sgrate2(NUPER,UPER,SSGG(1,2,1),SGG(J,2,1))
!          call sgrate2(NUPER,UPER,SSGG(1,3,1),SGG(J,3,1))
!          call sgrate2(NUPER,UPER,SSGG(1,2,2),SGG(J,2,2))
!          call sgrate2(NUPER,UPER,SSGG(1,3,2),SGG(J,3,2))
!          call sgrate2(NUPER,UPER,SSGG(1,3,3),SGG(J,3,3))




       end do

       ! -- Resonance relation -- !
       RRP=1.-NWCW-NPARA1*UPARMAX/SQMUT0
       RRM=1.-NWCW-NPARA1*UPARMIN/SQMUT0

       if (RRP*RRM.GT.0) then
          ! -- No resonance here -- !
          do J=1,NUPAR
             CUPAR=UPAR(J)
             DFACTPAR=NPARA1*CUPAR/SQMUT0
             RR=1.-NWCW-DFACTPAR
             IRR=1./RR
             do N=1,3
                do M=1,3
                   SGG(J,M,N)=SGG(J,M,N)*IRR
                end do
             end do
          end do

          !-- Hermitian part of $\Theta/(2\pi)$ --!
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,1),THETARE(1,1))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,1),THETARE(2,1))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,1),THETARE(3,1))
!          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,2),THETARE(1,2))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,2),THETARE(2,2))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,2),THETARE(3,2))
!          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,3),THETARE(1,3))
!          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,3),THETARE(2,3))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,3),THETARE(3,3))




 	    THETARE(1,2) = THETARE(2,1)
          THETARE(1,3) = THETARE(3,1)
          THETARE(2,3) = THETARE(3,2)

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          ! -- There is a resonance all right -- !
          ! -- 1. Hermitian part -- !
          UPAR0=SQMUT0/NPARA1*(1.-NWCW)
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
!          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,2),THETARE(1,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
!          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,3),THETARE(1,3))
!          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,3),THETARE(2,3))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))

          THETARE(1,2) = THETARE(2,1)
          THETARE(1,3) = THETARE(3,1)
          THETARE(2,3) = THETARE(3,2)

          ! -- 2. Anti-hermitian part -- !
          call WINTERP1D(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
          call WINTERP1D(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
!          call WINTERP1D(UPAR,NUPAR,SGG(1,1,2),UPAR0,THETAIM(1,2))
          call WINTERP1D(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
!          call WINTERP1D(UPAR,NUPAR,SGG(1,1,3),UPAR0,THETAIM(1,3))
!          call WINTERP1D(UPAR,NUPAR,SGG(1,2,3),UPAR0,THETAIM(2,3))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

          THETAIM(1,2) = THETAIM(2,1)
          THETAIM(1,3) = THETAIM(3,1)
          THETAIM(2,3) = THETAIM(3,2)

          do N=1,3
             do M=1,3
                THETARE(M,N)=-THETARE(M,N)*SQMUT0/NPARA1
                THETAIM(M,N)=-PI*THETAIM(M,N)*SQMUT0/abs(NPARA1)
             end do
          end do

       end if

       do N=1,3
          do M=1,3
             WSPEC(M,N)=WSPEC(M,N)+2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
          end do
       end do

    end do

    call WROTATE_AORSA(BETA1,BETA2,WSPEC(1,1))

  end subroutine GETNONMAX_SIGMA_AORSA_NEW0




!
!*************************************************************************
!

  subroutine WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER,UPER, &
       &                   NUPER,NBESSJ,NBESSJ_START,XI,JNXI,IFAIL)

!   ---------------------------------------
!   Nathan's version: optimized for Cray X1
!   ---------------------------------------

    implicit none


    real, intent(IN):: ZSPEC, ASPEC, ENORM,BMAG,KPER
    integer, intent(IN):: NUPER
    real, dimension(NUPER), intent(IN):: UPER
    integer, intent(IN):: NBESSJ,NBESSJ_START
    real, dimension(NUPER), intent(inout):: XI
    real, dimension(NUPER,NBESSJ), intent(inout):: JNXI
    integer, intent(inout):: IFAIL

    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: EVPERAMU=9.314943e8
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real:: WC,MUT0,SQMUT0

    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    WC=BMAG*ZSPEC*EOVERMH/ASPEC
    XI(:)=KPER*UPER(:)*C/SQMUT0/WC

    call BESSJ(XI(1),JNXI(1,1),NUPER,NBESSJ,NBESSJ_START)

  end subroutine WMATPRECALC_AORSA


!
!*************************************************************************
!

  subroutine GETNONMAXSWMAT_AORSA_NEW(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL)
    implicit none
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    complex, dimension(3,3), intent(inout):: WSPEC
    integer, intent(inout):: IFAIL

    complex:: BETAFACT
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
        
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPER,3,3):: SSGG
    real, dimension(NUPAR,3,3):: SGG
    real, dimension(3,3):: THETARE,THETAIM
    real:: W2,WP2,WPFACT,WCW,RRP,RRM,RR,IRR
    real:: MUT0,SQMUT0,BETA1,BETA2,KPARA1,NPARA1
    real:: ISQ2,CUPAR,CUPER,NWCW,DFACTPAR,DFACTPER,LF0,LNF0UPER
    real:: UPAR0,SFACT0,SFACTP1,SFACTM1
    integer:: NHARM,IHARM,NJ,NJP1,NJM1,J,K,M,N

    logical, parameter :: use_ppart6 = .true.
    logical :: is_uniform 
    integer, parameter :: nfxmax = 9
    integer :: nfx 
    real*8, dimension(nfxmax) ::  vint
    real*8, dimension(nupar,nfxmax) :: fx
    real*8 :: dx, dh, tol
    integer :: i

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2

    is_uniform = .true.
    tol = 1.0d-7
    dh = upar(2)-upar(1)
    do i=1,nupar-1
       dx = upar(i+1)-upar(i)
       is_uniform = abs(dx-dh).le. tol*dx
       if (.not.is_uniform) exit
    enddo

    ! -- Loop over harmonics -- !
    do IHARM=-NHARM,NHARM
       NWCW=real(IHARM)*WCW

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)
       JN0XI2(1:NUPER)=SFACT0*JNXI2(1:NUPER,NJ)
       JNP1XI2(1:NUPER)=SFACTP1*JNXI2(1:NUPER,NJP1)
       JNM1XI2(1:NUPER)=SFACTM1*JNXI2(1:NUPER,NJM1)

       ! -- Build array with integrand -- !
       do J=1,NUPAR
          CUPAR=UPAR(J)
          DFACTPAR=NPARA1*CUPAR/SQMUT0
          do K=1,NUPER
             CUPER=UPER(K)
             DFACTPER=NPARA1*CUPER/SQMUT0
             LF0=(1.-DFACTPAR)*DFDUPER(K,J)+DFACTPER*DFDUPAR(K,J)
             LNF0UPER=NWCW*DFDUPER(K,J)*CUPAR+(1.-NWCW)*DFDUPAR(K,J)*CUPER
             SSGG(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
             SSGG(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
             SSGG(K,3,1)=ISQ2*CUPAR*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
             SSGG(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
             SSGG(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
             SSGG(K,3,2)=ISQ2*CUPAR*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
             SSGG(K,1,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
             SSGG(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
             SSGG(K,3,3)=CUPAR*LNF0UPER*JN0XI1(K)*JN0XI2(K)
          end do

          ! -- Integrate of $u_\perp$ -- !
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,1),SGG(J,1,1))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,1),SGG(J,2,1))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,1),SGG(J,3,1))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,2),SGG(J,1,2))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,2),SGG(J,2,2))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,2),SGG(J,3,2))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,1,3),SGG(J,1,3))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,2,3),SGG(J,2,3))
          call EQSIMPSON1D(NUPER,UPER,SSGG(1,3,3),SGG(J,3,3))

       end do

       ! -- Resonance relation -- !
       RRP=1.-NWCW-NPARA1*UPARMAX/SQMUT0
       RRM=1.-NWCW-NPARA1*UPARMIN/SQMUT0

       if (RRP*RRM.GT.0) then
          ! -- No resonance here -- !
          do J=1,NUPAR
             CUPAR=UPAR(J)
             DFACTPAR=NPARA1*CUPAR/SQMUT0
             RR=1.-NWCW-DFACTPAR
             IRR=1./RR
             do N=1,3
                do M=1,3
                   SGG(J,M,N)=SGG(J,M,N)*IRR
                end do
             end do
          end do

          !-- Hermitian part of $\Theta/(2\pi)$ --!
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,1),THETARE(1,1))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,1),THETARE(2,1))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,1),THETARE(3,1))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,2),THETARE(1,2))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,2),THETARE(2,2))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,2),THETARE(3,2))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,1,3),THETARE(1,3))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,2,3),THETARE(2,3))
          call EQSIMPSON1D(NUPAR,UPAR,SGG(1,3,3),THETARE(3,3))

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          ! -- There is a resonance all right -- !
          ! -- 1. Hermitian part -- !
          UPAR0=SQMUT0/NPARA1*(1.-NWCW)

          
          if (use_ppart6) then

           fx(1:nupar,1) = SGG(1:nupar,1,1)
           fx(1:nupar,2) = SGG(1:nupar,2,1)
           fx(1:nupar,3) = SGG(1:nupar,3,1)

           fx(1:nupar,4) = SGG(1:nupar,1,2)
           fx(1:nupar,5) = SGG(1:nupar,2,2)
           fx(1:nupar,6) = SGG(1:nupar,3,2)

           fx(1:nupar,7) = SGG(1:nupar,1,3)
           fx(1:nupar,8) = SGG(1:nupar,2,3)
           fx(1:nupar,9) = SGG(1:nupar,3,3)

           nfx = 9
           call cauchy_ppart6(upar,nupar,upar0,nfx,fx,vint)

           thetare(1,1) = vint(1)
           thetare(2,1) = vint(2)
           thetare(3,1) = vint(3)

           thetare(1,2) = vint(4)
           thetare(2,2) = vint(5)
           thetare(3,2) = vint(6)

           thetare(1,3) = vint(7)
           thetare(2,3) = vint(8)
           thetare(3,3) = vint(9)

          else

          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))

          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,2),THETARE(1,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))

          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,3),THETARE(1,3))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,3),THETARE(2,3))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))

          endif

          ! -- 2. Anti-hermitian part -- !
          call WINTERP1D(UPAR,NUPAR,SGG(1,1,1),UPAR0,THETAIM(1,1))
          call WINTERP1D(UPAR,NUPAR,SGG(1,2,1),UPAR0,THETAIM(2,1))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,1),UPAR0,THETAIM(3,1))
          call WINTERP1D(UPAR,NUPAR,SGG(1,1,2),UPAR0,THETAIM(1,2))
          call WINTERP1D(UPAR,NUPAR,SGG(1,2,2),UPAR0,THETAIM(2,2))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,2),UPAR0,THETAIM(3,2))
          call WINTERP1D(UPAR,NUPAR,SGG(1,1,3),UPAR0,THETAIM(1,3))
          call WINTERP1D(UPAR,NUPAR,SGG(1,2,3),UPAR0,THETAIM(2,3))
          call WINTERP1D(UPAR,NUPAR,SGG(1,3,3),UPAR0,THETAIM(3,3))

          do N=1,3
             do M=1,3
                THETARE(M,N)=-THETARE(M,N)*SQMUT0/NPARA1
                THETAIM(M,N)=-PI*THETAIM(M,N)*SQMUT0/abs(NPARA1)
             end do
          end do

       end if

       do N=1,3
          do M=1,3
             WSPEC(M,N)=WSPEC(M,N)+2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
          end do
       end do

    end do

    call WROTATE_AORSA(BETA1,BETA2,WSPEC(1,1))

  end subroutine GETNONMAXSWMAT_AORSA_NEW

!
!*************************************************************************
!


  subroutine GETNONMAXSWMAT_AORSA(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL)
    implicit none
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    complex, dimension(3,3), intent(inout):: WSPEC
    integer, intent(inout):: IFAIL

    complex:: BETAFACT
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
	
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real:: W2,WP2,WPFACT,WCW,BETA1,BETA2
    real:: MUT0,SQMUT0,ISQ2,NWCW
    real:: SFACT0,SFACTP1,SFACTM1,RRP,RRM,KPARA1
    real:: KPARABS,NPARA1,CUPAR,CUPER,DFACTPAR,DFACTPER
    real:: RR,IRR,LF0,LNF0UPER,UPAR0,DUPARF,UPARMMAX
    real:: CLUPARF,CUPARP,CUPARM,DFACTPARP,DFACTPARM
    real:: IDFDUPARM,IDFDUPERM,IDFDUPARP,IDFDUPERP
    real:: LF0P,LF0M,LNF0PUPER,LNF0MUPER
    real:: DLF0,DLNF0UPER,UDLF0,UDLNF0UPER
    real:: DFDUPER0,DFDUPAR0
    real, dimension(3,3):: THETARE,THETAIM,G
    real, dimension(NUPER,3,3):: SWW
    real, dimension(NUPAR,3,3):: IRRG,GG,UDDGG
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPAR):: LUPARF,LLOGF
    integer:: IHARM,NHARM,NJ,NJP1,NJM1,I,J,J1,K,N,M

    GG = 0.0

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2

    ! -- Loop over harmonics -- !
    do IHARM=-NHARM,NHARM
       NWCW=real(IHARM)*WCW

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)
       JN0XI2(1:NUPER)=SFACT0*JNXI2(1:NUPER,NJ)
       JNP1XI2(1:NUPER)=SFACTP1*JNXI2(1:NUPER,NJP1)
       JNM1XI2(1:NUPER)=SFACTM1*JNXI2(1:NUPER,NJM1)

       ! -- Resonance relation -- !
       RRP=1.-NWCW-NPARA1*UPARMAX/SQMUT0
       RRM=1.-NWCW-NPARA1*UPARMIN/SQMUT0

       if (RRP*RRM.GT.0) then

          !--------------------------------------!
          !-- First case: vector, no resonant  --!
          !-- parallel velocity in the  domain --!
          !--------------------------------------!

          do J=1,NUPAR
             CUPAR=UPAR(J)
             DFACTPAR=NPARA1*CUPAR/SQMUT0
             RR=1.-NWCW-DFACTPAR
             IRR=1./RR

             ! --- Computation of $sww=u_\perp*w$ ---
             do K=1,NUPER
                CUPER=UPER(K)
                DFACTPER=NPARA1*CUPER/SQMUT0
                LF0=(1.-DFACTPAR)*DFDUPER(K,J)+DFACTPER*DFDUPAR(K,J)
                LNF0UPER=NWCW*DFDUPER(K,J)*CUPAR+(1.-NWCW)*DFDUPAR(K,J)*CUPER
                SWW(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,1)=ISQ2*CUPAR*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,2)=ISQ2*CUPAR*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
                SWW(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
                SWW(K,3,3)=CUPAR*LNF0UPER*JN0XI1(K)*JN0XI2(K)
             end do

             !-- Integration over $u_\perp$ --!
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),G(1,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),G(2,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),G(3,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),G(1,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),G(2,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),G(3,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),G(1,3))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),G(2,3))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),G(3,3))

             do N=1,3
                do M=1,3
                   IRRG(J,M,N)=G(M,N)*IRR
                end do
             end do

          end do

          !-- Hermitian part of $\Theta/(2\pi)$ --!
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,1),THETARE(1,1))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,1),THETARE(2,1))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,1),THETARE(3,1))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,2),THETARE(1,2))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,2),THETARE(2,2))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,2),THETARE(3,2))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,3),THETARE(1,3))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,3),THETARE(2,3))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,3),THETARE(3,3))

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          !---------------------------------------------------!
          !-- Second case: there is a resonant              --!
          !-- parallel velocity in the $u_\parallel$ domain --!
          !---------------------------------------------------!

          !-- Resonant parallel velocity --!
          UPAR0=SQMUT0/NPARA1*(1.-NWCW)

          ! --- 1. Hermitian part --- !

          !-- New $u_\parallel$ grid --!
          LUPARF(1)=0.
          LLOGF(1)=0.
          UPARMMAX=max(abs(UPARMIN),abs(UPARMAX))
          DUPARF=(UPARMMAX+ABS(UPAR0))/(NUPAR-1)
          do J1=2,NUPAR
             CLUPARF=DUPARF*(J1-1)
             LUPARF(J1)=CLUPARF
             LLOGF(J1)=CLUPARF*(LOG(ABS(CLUPARF))-1.)
          end do

          do J=1,NUPAR
             CUPARP=UPAR0+LUPARF(J)
             CUPARM=UPAR0-LUPARF(J)
             DFACTPARP=NPARA1*CUPARP/SQMUT0
             DFACTPARM=NPARA1*CUPARM/SQMUT0

             ! --- Computation of $sww=u_\perp*w$ ---
             do K=1,NUPER
                CUPER=UPER(K)
                DFACTPER=NPARA1*CUPER/SQMUT0
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR, &
                     & CUPER,CUPARP,IDFDUPARP)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER, &
                     & CUPER,CUPARP,IDFDUPERP)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR, &
                     & CUPER,CUPARM,IDFDUPARM)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER, &
                     & CUPER,CUPARM,IDFDUPERM)
                LF0P=(1.-DFACTPARP)*IDFDUPERP+DFACTPER*IDFDUPARP
                LF0M=(1.-DFACTPARM)*IDFDUPERM+DFACTPER*IDFDUPARM
                LNF0PUPER=NWCW*IDFDUPERP*CUPARP+(1.-NWCW)*IDFDUPARP*CUPER
                LNF0MUPER=NWCW*IDFDUPERM*CUPARM+(1.-NWCW)*IDFDUPARM*CUPER
                DLF0=LF0P-LF0M
                DLNF0UPER=LNF0PUPER-LNF0MUPER
                UDLF0=CUPARP*LF0P-CUPARM*LF0M
                UDLNF0UPER=CUPARP*LNF0PUPER-CUPARM*LNF0MUPER
                SWW(K,1,1)=0.5*CUPER*DLF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,1)=0.5*CUPER*DLF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,1)=ISQ2*UDLF0*JNP1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,2)=0.5*CUPER*DLF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,2)=0.5*CUPER*DLF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,2)=ISQ2*UDLF0*JNM1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,3)=ISQ2*CUPER*DLNF0UPER*JN0XI1(K)*JNP1XI2(K)
                SWW(K,2,3)=ISQ2*CUPER*DLNF0UPER*JN0XI1(K)*JNM1XI2(K)
                SWW(K,3,3)=UDLNF0UPER*JN0XI1(K)*JN0XI2(K)
             end do

             !-- Integration over $u_\perp$ --!
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),GG(J,1,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),GG(J,2,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),GG(J,3,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),GG(J,1,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),GG(J,2,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),GG(J,3,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),GG(J,1,3))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),GG(J,2,3))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),GG(J,3,3))

          end do

          !-- $gg''*llogf$ --!
          UDDGG(1,1:3,1:3)=0.
          UDDGG(NUPAR,1:3,1:3)=0.
          do J=2,NUPAR-1
             UDDGG(J,1,1)=LLOGF(J)* &
                  & (GG(J+1,1,1)-2.*GG(J,1,1)+GG(J-1,1,1))/DUPARF**2
             UDDGG(J,2,1)=LLOGF(J)* &
                  & (GG(J+1,2,1)-2.*GG(J,2,1)+GG(J-1,2,1))/DUPARF**2
             UDDGG(J,3,1)=LLOGF(J)* &
                  & (GG(J+1,3,1)-2.*GG(J,3,1)+GG(J-1,3,1))/DUPARF**2
             UDDGG(J,1,2)=LLOGF(J)* &
                  & (GG(J+1,1,2)-2.*GG(J,1,2)+GG(J-1,1,2))/DUPARF**2
             UDDGG(J,2,2)=LLOGF(J)* &
                  & (GG(J+1,2,2)-2.*GG(J,2,2)+GG(J-1,2,2))/DUPARF**2
             UDDGG(J,3,2)=LLOGF(J)* &
                  & (GG(J+1,3,2)-2.*GG(J,3,2)+GG(J-1,3,2))/DUPARF**2
             UDDGG(J,1,3)=LLOGF(J)* &
                  &(GG(J+1,1,3)-2.*GG(J,1,3)+GG(J-1,1,3))/DUPARF**2
             UDDGG(J,2,3)=LLOGF(J)* &
                  & (GG(J+1,2,3)-2.*GG(J,2,3)+GG(J-1,2,3))/DUPARF**2
             UDDGG(J,3,3)=LLOGF(J)* &
                  & (GG(J+1,3,3)-2.*GG(J,3,3)+GG(J-1,3,3))/DUPARF**2
          end do

          !-- Real part of $\Theta/(2\pi)$ --!
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,1),THETARE(1,1))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,1),THETARE(2,1))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,1),THETARE(3,1))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,2),THETARE(1,2))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,2),THETARE(2,2))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,2),THETARE(3,2))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,3),THETARE(1,3))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,3),THETARE(2,3))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,3),THETARE(3,3))

          do N=1,3
             do M=1,3
                THETARE(M,N)=-THETARE(M,N)*SQMUT0/NPARA1
             end do
          end do

          ! --- 2. Anti-hermitian part --- !
          do K=1,NUPER
             CUPER=UPER(K)
             DFACTPER=NPARA1*CUPER/SQMUT0
             call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR,CUPER,UPAR0,DFDUPAR0)
             call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER,CUPER,UPAR0,DFDUPER0)

             LF0=NWCW*DFDUPER0+DFACTPER*DFDUPAR0
             LNF0UPER=NWCW*DFDUPER0*UPAR0+(1.-NWCW)*DFDUPAR0*CUPER

             SWW(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
             SWW(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
             SWW(K,3,1)=ISQ2*UPAR0*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
             SWW(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
             SWW(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
             SWW(K,3,2)=ISQ2*UPAR0*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
             SWW(K,1,3)=ISQ2*CUPAR*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
             SWW(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
             SWW(K,3,3)=UPAR0*LNF0UPER*JN0XI1(K)*JN0XI2(K)
          end do

          !--- Integration over $u_perp$ ---!
          call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),G(1,1))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),G(2,1))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),G(3,1))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),G(1,2))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),G(2,2))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),G(3,2))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),G(1,3))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),G(2,3))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),G(3,3))

          !-- Imaginary part of $\Theta/(2\pi)$ --!
          do N=1,3
             do M=1,3
                THETAIM(M,N)=-PI*G(M,N)*SQMUT0/abs(NPARA1)
             end do
          end do

       end if

       do N=1,3
          do M=1,3
             WSPEC(M,N)=WSPEC(M,N)+ &
                  & 2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
          end do
       end do
    end do

    call WROTATE_AORSA(BETA1,BETA2,WSPEC(1,1))

  end subroutine GETNONMAXSWMAT_AORSA

!
!*************************************************************************
!

    subroutine GETNONMAX_SIGMA_AORSA(W, &
       & ZSPEC, ASPEC, DENS, BMAG, &
       & K1, XI1, JNXI1, & 
       & K2, XI2, JNXI2, NBESSJ, & 
       & ENORM, UPARMIN, UPARMAX, &
       & NUPAR, NUPER, UPER, UPAR, & 
       & DFDUPER, DFDUPAR, & 
       & WSPEC, IFAIL)
       
    implicit none
    
    real, intent(IN):: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3):: K1,K2
    integer, intent(IN):: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN):: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN):: JNXI1,JNXI2
    real, intent(IN):: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER,DFDUPAR
    complex, dimension(3,3), intent(inout):: WSPEC
    integer, intent(inout):: IFAIL

    complex:: BETAFACT
    
    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real:: W2,WP2,WPFACT,WCW,BETA1,BETA2
    real:: MUT0,SQMUT0,ISQ2,NWCW
    real:: SFACT0,SFACTP1,SFACTM1,RRP,RRM,KPARA1
    real:: KPARABS,NPARA1,CUPAR,CUPER,DFACTPAR,DFACTPER
    real:: RR,IRR,LF0,LNF0UPER,UPAR0,DUPARF,UPARMMAX
    real:: CLUPARF,CUPARP,CUPARM,DFACTPARP,DFACTPARM
    real:: IDFDUPARM,IDFDUPERM,IDFDUPARP,IDFDUPERP
    real:: LF0P,LF0M,LNF0PUPER,LNF0MUPER
    real:: DLF0,DLNF0UPER,UDLF0,UDLNF0UPER
    real:: DFDUPER0,DFDUPAR0
    real, dimension(3,3):: THETARE,THETAIM,G
    real, dimension(NUPER,3,3):: SWW
    real, dimension(NUPAR,3,3):: IRRG,GG,UDDGG
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPAR):: LUPARF,LLOGF
    integer:: IHARM,NHARM,NJ,NJP1,NJM1,I,J,J1,K,N,M

    GG = 0.0

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERMH/ASPEC/W
    BETA1=ATAN2(K1(2),K1(1))
    BETA2=ATAN2(K2(2),K2(1))
    KPARA1=K1(3)
    NPARA1=KPARA1*C/W
    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    ISQ2=SQRT(0.5)
    NHARM=NBESSJ-2

    ! -- Loop over harmonics -- !
    do IHARM=-NHARM,NHARM
       NWCW=real(IHARM)*WCW

       BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
       NJ=ABS(IHARM)+1
       NJP1=ABS(IHARM+1)+1
       NJM1=ABS(IHARM-1)+1
       SFACT0=real(sign(1,IHARM))**ABS(IHARM)
       SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
       SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)
       JN0XI2(1:NUPER)=SFACT0*JNXI2(1:NUPER,NJ)
       JNP1XI2(1:NUPER)=SFACTP1*JNXI2(1:NUPER,NJP1)
       JNM1XI2(1:NUPER)=SFACTM1*JNXI2(1:NUPER,NJM1)

       ! -- Resonance relation -- !
       RRP=1.-NWCW-NPARA1*UPARMAX/SQMUT0
       RRM=1.-NWCW-NPARA1*UPARMIN/SQMUT0

       if (RRP*RRM.GT.0) then

          !--------------------------------------!
          !-- First case: vector, no resonant  --!
          !-- parallel velocity in the  domain --!
          !--------------------------------------!

          do J=1,NUPAR
             CUPAR=UPAR(J)
             DFACTPAR=NPARA1*CUPAR/SQMUT0
             RR=1.-NWCW-DFACTPAR
             IRR=1./RR

             ! --- Computation of $sww=u_\perp*w$ ---
             do K=1,NUPER
                CUPER=UPER(K)
                DFACTPER=NPARA1*CUPER/SQMUT0
                LF0=(1.-DFACTPAR)*DFDUPER(K,J)+DFACTPER*DFDUPAR(K,J)
                LNF0UPER=NWCW*DFDUPER(K,J)*CUPAR+(1.-NWCW)*DFDUPAR(K,J)*CUPER
                SWW(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,1)=ISQ2*CUPAR*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,2)=ISQ2*CUPAR*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
                SWW(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
                SWW(K,3,3)=CUPAR*LNF0UPER*JN0XI1(K)*JN0XI2(K)
             end do

             !-- Integration over $u_\perp$ --!
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),G(1,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),G(2,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),G(3,1))
           ! call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),G(1,2)) !
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),G(2,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),G(3,2))
           ! call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),G(1,3)) !
           ! call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),G(2,3)) !
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),G(3,3))
	
	     G(1,2) = - G(2,1)
	     G(1,3) =   G(3,1)
	     G(2,3) = - G(3,2)

             do N=1,3
                do M=1,3
                   IRRG(J,M,N)=G(M,N)*IRR
                end do
             end do

          end do

          !-- Hermitian part of $\Theta/(2\pi)$ --!
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,1),THETARE(1,1))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,1),THETARE(2,1))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,1),THETARE(3,1))
       !  call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,2),THETARE(1,2)) !
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,2),THETARE(2,2))
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,2),THETARE(3,2))
       !  call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,1,3),THETARE(1,3)) !
       !  call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,2,3),THETARE(2,3)) !
          call EQTRAPZ1D(NUPAR,UPAR,IRRG(1,3,3),THETARE(3,3))
	
	  THETARE(1,2) = - THETARE(2,1)
	  THETARE(1,3) =   THETARE(3,1)
	  THETARE(2,3) = - THETARE(3,2)

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          !---------------------------------------------------!
          !-- Second case: there is a resonant              --!
          !-- parallel velocity in the $u_\parallel$ domain --!
          !---------------------------------------------------!

          !-- Resonant parallel velocity --!
          UPAR0=SQMUT0/NPARA1*(1.-NWCW)

          ! --- 1. Hermitian part --- !

          !-- New $u_\parallel$ grid --!
          LUPARF(1)=0.
          LLOGF(1)=0.
          UPARMMAX=max(abs(UPARMIN),abs(UPARMAX))
          DUPARF=(UPARMMAX+ABS(UPAR0))/(NUPAR-1)
          do J1=2,NUPAR
             CLUPARF=DUPARF*(J1-1)
             LUPARF(J1)=CLUPARF
             LLOGF(J1)=CLUPARF*(LOG(ABS(CLUPARF))-1.)
          end do

          do J=1,NUPAR
             CUPARP=UPAR0+LUPARF(J)
             CUPARM=UPAR0-LUPARF(J)
             DFACTPARP=NPARA1*CUPARP/SQMUT0
             DFACTPARM=NPARA1*CUPARM/SQMUT0

             ! --- Computation of $sww=u_\perp*w$ ---
             do K=1,NUPER
                CUPER=UPER(K)
                DFACTPER=NPARA1*CUPER/SQMUT0
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR, &
                     & CUPER,CUPARP,IDFDUPARP)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER, &
                     & CUPER,CUPARP,IDFDUPERP)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR, &
                     & CUPER,CUPARM,IDFDUPARM)
                call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER, &
                     & CUPER,CUPARM,IDFDUPERM)
                LF0P=(1.-DFACTPARP)*IDFDUPERP+DFACTPER*IDFDUPARP
                LF0M=(1.-DFACTPARM)*IDFDUPERM+DFACTPER*IDFDUPARM
                LNF0PUPER=NWCW*IDFDUPERP*CUPARP+(1.-NWCW)*IDFDUPARP*CUPER
                LNF0MUPER=NWCW*IDFDUPERM*CUPARM+(1.-NWCW)*IDFDUPARM*CUPER
                DLF0=LF0P-LF0M
                DLNF0UPER=LNF0PUPER-LNF0MUPER
                UDLF0=CUPARP*LF0P-CUPARM*LF0M
                UDLNF0UPER=CUPARP*LNF0PUPER-CUPARM*LNF0MUPER
                SWW(K,1,1)=0.5*CUPER*DLF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,1)=0.5*CUPER*DLF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,1)=ISQ2*UDLF0*JNP1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,2)=0.5*CUPER*DLF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
                SWW(K,2,2)=0.5*CUPER*DLF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
                SWW(K,3,2)=ISQ2*UDLF0*JNM1XI1(K)*JN0XI2(K)*CUPER
                SWW(K,1,3)=ISQ2*CUPER*DLNF0UPER*JN0XI1(K)*JNP1XI2(K)
                SWW(K,2,3)=ISQ2*CUPER*DLNF0UPER*JN0XI1(K)*JNM1XI2(K)
                SWW(K,3,3)=UDLNF0UPER*JN0XI1(K)*JN0XI2(K)
             end do

             !-- Integration over $u_\perp$ --!
             call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),GG(J,1,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),GG(J,2,1))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),GG(J,3,1))
          !  call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),GG(J,1,2)) !
             call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),GG(J,2,2))
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),GG(J,3,2))
          !  call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),GG(J,1,3)) !
          !  call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),GG(J,2,3)) !
             call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),GG(J,3,3))
	
	     GG(J,1,2) = - GG(J,2,1)
	     GG(J,1,3) =   GG(J,3,1)
	     GG(J,2,3) = - GG(J,3,2)

          end do

          !-- $gg''*llogf$ --!
          UDDGG(1,1:3,1:3)=0.
          UDDGG(NUPAR,1:3,1:3)=0.
          do J=2,NUPAR-1
             UDDGG(J,1,1)=LLOGF(J)* &
                  & (GG(J+1,1,1)-2.*GG(J,1,1)+GG(J-1,1,1))/DUPARF**2
             UDDGG(J,2,1)=LLOGF(J)* &
                  & (GG(J+1,2,1)-2.*GG(J,2,1)+GG(J-1,2,1))/DUPARF**2
             UDDGG(J,3,1)=LLOGF(J)* &
                  & (GG(J+1,3,1)-2.*GG(J,3,1)+GG(J-1,3,1))/DUPARF**2
             UDDGG(J,1,2)=LLOGF(J)* &
                  & (GG(J+1,1,2)-2.*GG(J,1,2)+GG(J-1,1,2))/DUPARF**2
             UDDGG(J,2,2)=LLOGF(J)* &
                  & (GG(J+1,2,2)-2.*GG(J,2,2)+GG(J-1,2,2))/DUPARF**2
             UDDGG(J,3,2)=LLOGF(J)* &
                  & (GG(J+1,3,2)-2.*GG(J,3,2)+GG(J-1,3,2))/DUPARF**2
             UDDGG(J,1,3)=LLOGF(J)* &
                  &(GG(J+1,1,3)-2.*GG(J,1,3)+GG(J-1,1,3))/DUPARF**2
             UDDGG(J,2,3)=LLOGF(J)* &
                  & (GG(J+1,2,3)-2.*GG(J,2,3)+GG(J-1,2,3))/DUPARF**2
             UDDGG(J,3,3)=LLOGF(J)* &
                  & (GG(J+1,3,3)-2.*GG(J,3,3)+GG(J-1,3,3))/DUPARF**2
          end do

          !-- Real part of $\Theta/(2\pi)$ --!
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,1),THETARE(1,1))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,1),THETARE(2,1))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,1),THETARE(3,1))
       !  call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,2),THETARE(1,2)) !
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,2),THETARE(2,2))
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,2),THETARE(3,2))
       !  call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,1,3),THETARE(1,3)) !
       !  call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,2,3),THETARE(2,3)) !
          call EQTRAPZ1D(NUPAR,LUPARF,UDDGG(1,3,3),THETARE(3,3))
	
	  THETARE(1,2) = - THETARE(2,1)
	  THETARE(1,3) =   THETARE(3,1)
	  THETARE(2,3) = - THETARE(3,2)

          do N=1,3
             do M=1,3
                THETARE(M,N)=-THETARE(M,N)*SQMUT0/NPARA1
             end do
          end do

          ! --- 2. Anti-hermitian part --- !
          do K=1,NUPER
             CUPER=UPER(K)
             DFACTPER=NPARA1*CUPER/SQMUT0
             call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPAR,CUPER,UPAR0,DFDUPAR0)
             call WINTERP2D(UPER,NUPER,UPAR,NUPAR,DFDUPER,CUPER,UPAR0,DFDUPER0)

             LF0=NWCW*DFDUPER0+DFACTPER*DFDUPAR0
             LNF0UPER=NWCW*DFDUPER0*UPAR0+(1.-NWCW)*DFDUPAR0*CUPER

             SWW(K,1,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNP1XI2(K)*CUPER
             SWW(K,2,1)=0.5*CUPER*LF0*JNP1XI1(K)*JNM1XI2(K)*CUPER
             SWW(K,3,1)=ISQ2*UPAR0*LF0*JNP1XI1(K)*JN0XI2(K)*CUPER
             SWW(K,1,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNP1XI2(K)*CUPER
             SWW(K,2,2)=0.5*CUPER*LF0*JNM1XI1(K)*JNM1XI2(K)*CUPER
             SWW(K,3,2)=ISQ2*UPAR0*LF0*JNM1XI1(K)*JN0XI2(K)*CUPER
             SWW(K,1,3)=ISQ2*CUPAR*LNF0UPER*JN0XI1(K)*JNP1XI2(K)
             SWW(K,2,3)=ISQ2*CUPER*LNF0UPER*JN0XI1(K)*JNM1XI2(K)
             SWW(K,3,3)=UPAR0*LNF0UPER*JN0XI1(K)*JN0XI2(K)
          end do

          !--- Integration over $u_perp$ ---!
          call EQTRAPZ1D(NUPER,UPER,SWW(1,1,1),G(1,1))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,2,1),G(2,1))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,1),G(3,1))
       !  call EQTRAPZ1D(NUPER,UPER,SWW(1,1,2),G(1,2)) !
          call EQTRAPZ1D(NUPER,UPER,SWW(1,2,2),G(2,2))
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,2),G(3,2))
       !  call EQTRAPZ1D(NUPER,UPER,SWW(1,1,3),G(1,3)) !
       !  call EQTRAPZ1D(NUPER,UPER,SWW(1,2,3),G(2,3)) !
          call EQTRAPZ1D(NUPER,UPER,SWW(1,3,3),G(3,3))
	
	  G(1,2) = - G(2,1)
	  G(1,3) =   G(3,1)
	  G(2,3) = - G(3,2)

          !-- Imaginary part of $\Theta/(2\pi)$ --!
          do N=1,3
             do M=1,3
                THETAIM(M,N)=-PI*G(M,N)*SQMUT0/abs(NPARA1)
             end do
          end do

       end if

       do N=1,3
          do M=1,3
             WSPEC(M,N)=WSPEC(M,N)+ &
                  & 2.*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
          end do
       end do
    end do

    call WROTATE_AORSA(BETA1,BETA2,WSPEC(1,1))

  end subroutine GETNONMAX_SIGMA_AORSA

!
!*************************************************************************
!




  subroutine BESSJ(XI,JNXI,LXI,LNBESSJ,LNSBESSJ)
    implicit none
    integer, intent(IN):: LXI,LNBESSJ,LNSBESSJ
    real, dimension(LXI), intent(IN):: XI
    real, dimension(LXI,LNBESSJ), intent(inout):: JNXI
    real, parameter:: BIGNO=1.e10,BIGNI=1.e-10
    real, dimension(LXI):: TOX,AXI,J0XI,J1XI
    real, dimension(LXI):: BJP,BJ,BJM,JSUM
    integer:: MSUM,IX,N,NP
    integer, dimension(LXI):: MCOR

    ! --- Returns the J bessel functions from 0 to LNBESSJ --- !
    ! --- LNSBESSJ is the starting point for the downard recurence --- !
    ! --- XI(LXI) is the argument array and JNXI(LXI,LNBESSJ) is the output ---!

!dir$ preferstream,prefervector
    do IX=1,LXI
    JNXI(IX,1:LNBESSJ)=0.

    TOX(IX)=0.
    AXI(IX)=ABS(XI(IX))
    if (AXI(IX).GT.0.)then
       TOX(IX)=2./AXI(IX)
    endif

    BJP(IX)=0.
    BJ(IX)=1.
    JSUM(IX)=0.
    enddo
!    MSUM=0
    do N=LNSBESSJ,1,-1
!dir$ concurrent,preferstream,prefervector
      do IX=1,LXI
      MSUM=mod((LNSBESSJ-N),2)
       MCOR(IX)=0
       BJM(IX)=REAL(N)*TOX(IX)*BJ(IX)-BJP(IX)
       BJP(IX)=BJ(IX)
       BJ(IX)=BJM(IX)
       if (ABS(BJ(IX)).GT.BIGNO)then
          BJ(IX)=BJ(IX)*BIGNI
          BJP(IX)=BJP(IX)*BIGNI
          JSUM(IX)=JSUM(IX)*BIGNI
          MCOR(IX)=1
       endif
       JSUM(IX)=JSUM(IX)+REAL(MSUM)*BJ(IX)
!       MSUM=1-MSUM
       if (N.LT.LNBESSJ) JNXI(IX,N+1)=BJP(IX)
       do NP=N+2,LNBESSJ
          if (MCOR(IX).EQ.1)then
             JNXI(IX,NP)=JNXI(IX,NP)*BIGNI
          endif
       end do
      enddo
    end do
    do IX=1,LXI

    JSUM(IX)=2.*JSUM(IX)-BJ(IX)
    do N=1,LNBESSJ
       JNXI(IX,N)=JNXI(IX,N)/JSUM(IX)
    end do
    end do

    ! --- Upward Recurrence for x>=n --- !
    call BESSJ0(AXI,J0XI,LXI)
    call BESSJ1(AXI,J1XI,LXI)
    do IX=1,LXI
    BJM(IX)=J0XI(IX)
    BJ(IX)=J1XI(IX)
    do N=1,LNBESSJ
       if (AXI(IX)>=real(N-1))then
          JNXI(IX,N)=BJM(IX)
          BJP(IX)=REAL(N)*TOX(IX)*BJ(IX)-BJM(IX)
          BJM(IX)=BJ(IX)
          BJ(IX)=BJP(IX)
       endif
    end do

    ! --- Zero --- !
       if (XI(IX)==0.) then
          JNXI(IX,1)=1.
          do N=2,LNBESSJ
             JNXI(IX,N)=0.
          end do
       end if
!    end do

    ! --- Revert when needed --- !
    ! --- N.B. J1 does not need to be reverted --- !
!    do IX=1,LXI
       do N=4,LNBESSJ,2
          if (XI(IX).LT.0.) JNXI(IX,N)=-JNXI(IX,N)
       end do
    end do

  end subroutine BESSJ

!
!*************************************************************************
!


  subroutine BESSJ0(XI,J0XI,LXI)
    implicit none
    integer, intent(IN):: LXI
    real, dimension(LXI), intent(IN) :: XI
    real, dimension(LXI), intent(inout) :: J0XI
    real, dimension(LXI) :: YI,AXI,XXI,ZI
    real, dimension(5):: P,Q
    real, dimension(6):: R,S
    P=(/1.0,-0.1098628627e-2,0.2734510407e-4, &
         -0.2073370639e-5,0.2093887211e-6/)
    Q=(/-0.1562499995e-1,0.1430488765e-3, &
         -0.6911147651e-5,0.7621095161e-6,-0.934945152e-7/)
    R=(/57568490574.0,-13362590354.0,651619640.7, &
         -11214424.18,77392.33017,-184.9052456/)
    S=(/57568490411.0,1029532985.0,&
         9494680.718,59272.64853,267.8532712,1.0/)

    where (abs(XI).LT.8.0)
       YI=XI*XI
       J0XI=(R(1)+YI*(R(2)+YI*(R(3)+YI*(R(4)+YI*(R(5)+YI*R(6))))))/ &
            (S(1)+YI*(S(2)+YI*(S(3)+YI*(S(4)+YI*(S(5)+YI*S(6))))))
    elsewhere
       AXI=abs(XI)
       ZI=8./AXI
       YI=ZI*ZI
       XXI=AXI-0.785398164
       J0XI=sqrt(0.636619772/AXI)* &
            (cos(XXI)*(P(1)+YI*(P(2)+YI*(P(3)+YI*(P(4)+YI*P(5))))) &
            - ZI*sin(XXI)*(Q(1)+YI*(Q(2)+YI*(Q(3)+YI*(Q(4)+YI*Q(5))))))
    end where
  end subroutine BESSJ0

!
!*************************************************************************
!



  subroutine BESSJ1(XI,J1XI,LXI)
    implicit none
    integer, intent(IN):: LXI
    real, dimension(LXI), intent(IN):: XI
    real, dimension(LXI), intent(inout):: J1XI
    real, dimension(LXI) :: YI,AXI,XXI,ZI
    real, dimension(5):: P,Q
    real, dimension(6):: R,S
    integer:: IX

    P=(/1.0,0.183105e-2,-0.3516396496e-4, &
         0.2457520174e-5,-0.240337019e-6/)
    Q=(/0.04687499995,-0.2002690873e-3, &
         0.8449199096e-5,-0.88228987e-6,0.105787412e-6/)
    R=(/72362614232.0,-7895059235.0,242396853.1, &
         -2972611.439,15704.48260,-30.16036606/)
    S=(/144725228442.0,2300535178.0,18583304.74, &
         99447.43394,376.9991397,1.0/)

    where (abs(XI).LT.8.0)
       YI=XI*XI
       J1XI=XI*(R(1)+YI*(R(2)+YI*(R(3)+YI*(R(4)+YI*(R(5)+YI*R(6))))))/ &
            (S(1)+YI*(S(2)+YI*(S(3)+YI*(S(4)+YI*(S(5)+YI*S(6))))))
    elsewhere
       AXI=abs(XI)
       ZI=8./AXI
       YI=ZI*ZI
       XXI=AXI-2.356194491
       J1XI=sqrt(0.636619772/AXI)* &
            (cos(XXI)*(P(1)+YI*(P(2)+YI*(P(3)+YI*(P(4)+YI*P(5))))) &
            - ZI*sin(XXI)*(Q(1)+YI*(Q(2)+YI*(Q(3)+YI*(Q(4)+YI*Q(5)))))) &
            *sign(1.0,XI)
    end where

    ! --- Revert if needed --- !
    do IX=1,LXI
       if (XI(IX).LT.0.) J1XI(IX)=-J1XI(IX)
    end do

  end subroutine BESSJ1

!
!*************************************************************************
!


  subroutine EQTRAPZ1D(LNU,LU,ARRAY,INTEGRAL)
    ! -- 1D trapezoidal integration for equally spaced coordinates -- !
    implicit none
    integer, intent(IN):: LNU
    real, dimension(LNU), intent(IN):: LU,ARRAY
    real, intent(inout):: INTEGRAL
    real:: DU

    DU=(LU(LNU)-LU(1))/(LNU-1.)
    INTEGRAL=(sum(ARRAY(2:LNU-1))+0.5*(ARRAY(1)+ARRAY(LNU)))*DU

  end subroutine EQTRAPZ1D

!
!*************************************************************************
!

      subroutine sgrate2(nxmax, x, f, ans)

      implicit none

      integer n, nxmax
      real f(nxmax), x(nxmax), ans

      ans = 0.0
      do n = 1, nxmax-1
         ans = ans + (x(n+1) - x(n)) * (f(n) + f(n+1)) / 2.0
      end do

      return
      end
!
!*************************************************************************
!


  subroutine EQSIMPSON1D(LNU,LU,ARRAY,INTEGRAL)
    ! -- 1D Simpson integration for equally spaced coordinates and odd N -- !
    implicit none
    integer, intent(IN):: LNU
    real, dimension(LNU), intent(IN):: LU,ARRAY
    real, intent(inout):: INTEGRAL

    real:: DU
    integer:: I

    DU=(LU(LNU)-LU(1))/(LNU-1.)
    INTEGRAL=2.*sum(ARRAY(2:LNU-1))+ARRAY(1)+ARRAY(LNU)
    do I=1,(LNU-1)/2
       INTEGRAL=INTEGRAL+2.*ARRAY(2*I)
    end do
    INTEGRAL=INTEGRAL*DU / 3.

  end subroutine EQSIMPSON1D

!
!*************************************************************************
!
  subroutine EQSIMPSON1D_2(LNU,LU,ARRAY,INTEGRAL)
    ! -- 1D Simpson integration for equally spaced coordinates and odd N -- !
    implicit none
    integer, intent(IN):: LNU
    real, dimension(LNU), intent(IN):: LU,ARRAY
    real, intent(inout):: INTEGRAL

    real:: DU

    DU=(LU(LNU)-LU(1))/(LNU-1.)
    INTEGRAL=ARRAY(1)+ARRAY(LNU)+ 2.*sum(array(2:lnu-1))
    integral = integral + 2.*sum(array(2:lnu-1:2))

    INTEGRAL=INTEGRAL*DU / 3.

  end subroutine EQSIMPSON1D_2

!
!*************************************************************************
!


  subroutine WINTERP2D(X,LX,Y,LY,FXY,X0,Y0,FXY0)
    implicit none
    integer, intent(IN):: LX,LY
    real, dimension(LX), intent(IN):: X
    real, dimension(LY), intent(IN):: Y
    real, dimension(LX,LY), intent(IN):: FXY
    real, intent(IN):: X0,Y0
    real, intent(inout):: FXY0
    real:: A0,A1,A2,A3,A4,A5,P,Q
    real:: XMIN,XMAX,YMIN,YMAX,DX,DY
    integer:: I,J

    FXY0=0.
    XMIN=X(1)
    YMIN=Y(1)
    XMAX=X(LX)
    YMAX=Y(LY)
    DX=(XMAX-XMIN)/(LX-1)
    DY=(YMAX-YMIN)/(LY-1)
    if ((X0>=XMIN).and.(X0<=XMAX).and.(Y0>=YMIN).and.(Y0<=YMAX)) then

       I=int(1.+(X0-XMIN)/DX)
       J=int(1.+(Y0-YMIN)/DY)
       P=(X0-X(I))/DX
       Q=(Y0-Y(J))/DY

       if (I.eq.LX) I=LX-1
       if (J.eq.LY) J=LY-1
       if ((I.ne.1).and.(J.ne.1)) then
          A0=0.5*Q*(Q-1.)
          A1=0.5*P*(P-1.)
          A2=1.+P*Q-P*P-Q*Q
          A3=0.5*P*(P-2.*Q+1.)
          A4=0.5*Q*(Q-2.*P+1.)
          A5=P*Q
          FXY0=A0*FXY(I,J-1)+A1*FXY(I-1,J)+A2*FXY(I,J)+A3*FXY(I+1,J) &
               +A4*FXY(I,J+1)+A5*FXY(I+1,J+1)
       else
          A0=(1.-P)*(1.-Q)
          A1=P*(1.-Q)
          A2=Q*(1.-P)
          A3=P*Q
          FXY0=A0*FXY(I,J)+A1*FXY(I+1,J)+A2*FXY(I,J+1)+A3*FXY(I+1,J+1)
       end if
   end if
  end subroutine WINTERP2D

!
!*************************************************************************
!

  SUBROUTINE WROTATE_AORSA(BETA1,BETA2,W)
    !     Apply unitary polarization matrices, eg., Wout = C2*Win*C1
    !     Transforms the perpendicular part of the tensor from right
    !     and left polarized components into fixed real components.
    IMPLICIT NONE
    real:: BETA1,BETA2
    complex :: W(3,3),WW(3,3),C1(3,3),C2(3,3)
    complex :: C1P,C1M,C2P,C2M
    INTEGER I,J,K
    real, parameter :: ZERO = 0.0, ONE = 1.0
    real, parameter :: ROOTHALF = 0.707106781
    C1P = ROOTHALF*EXP(CMPLX(ZERO,BETA1))
    C1M = ROOTHALF*EXP(CMPLX(ZERO,-BETA1))
    C2P = ROOTHALF*EXP(CMPLX(ZERO,BETA2))
    C2M = ROOTHALF*EXP(CMPLX(ZERO,-BETA2))
    C1(1,1) = C1P
    C1(2,1) = C1M
    C1(3,1) = CMPLX(ZERO,ZERO)
    C1(1,2) = C1P*CMPLX(ZERO,-ONE)
    C1(2,2) = C1M*CMPLX(ZERO,ONE)
    C1(3,2) = CMPLX(ZERO,ZERO)
    C1(1,3) = CMPLX(ZERO,ZERO)
    C1(2,3) = CMPLX(ZERO,ZERO)
    C1(3,3) = CMPLX(ONE,ZERO)
    C2(1,1) = C2M
    C2(2,1) = C2M*CMPLX(ZERO,ONE)
    C2(3,1) = CMPLX(ZERO,ZERO)
    C2(1,2) = C2P
    C2(2,2) = C2P*CMPLX(ZERO,-ONE)
    C2(3,2) = CMPLX(ZERO,ZERO)
    C2(1,3) = CMPLX(ZERO,ZERO)
    C2(2,3) = CMPLX(ZERO,ZERO)
    C2(3,3) = CMPLX(ONE,ZERO)
    DO I=1,3
       DO J=1,3
          WW(I,J) = CMPLX(ZERO,ZERO)
          DO K=1,3
             WW(I,J) = WW(I,J) + W(I,K)*C1(K,J)
          enddo
       enddo
    enddo
    DO I=1,3
       DO J=1,3
          W(I,J) = CMPLX(ZERO,ZERO)
          DO K=1,3
             W(I,J) = W(I,J) + C2(I,K)*WW(K,J)
          enddo
       enddo
    enddo
    RETURN
  END subroutine WROTATE_AORSA

!
!*************************************************************************
!


  subroutine CAUCHY_PPART(X,NX,XRES,FX,PINT)
    ! -- Returns principal part of int(x,fx/(x-xres)) -- !
    implicit none
    integer, intent(IN):: NX
    real, dimension(NX), intent(IN):: X,FX
    real, intent(IN):: XRES
    real, intent(inout):: PINT

    real:: CXP,CXM,FXP,FXM,LLOGF,XMAX,XMIN
    real:: CXF,DXF,DX
    real, dimension(NX):: LXF,INTGD,FXPP
    integer:: IX,IXM1,IXP1

    IXM1(IX)=MAX(IX-1,1)
    IXP1(IX)=MIN(IX+1,NX)

    XMAX=X(NX)
    XMIN=X(1)
    DXF=(XMAX+ABS(XRES))/(NX-1)
    DX=(XMAX-XMIN)/(NX-1)

    do IX=1,NX
       CXF=DXF*(IX-1)
       FXPP(IX)=(FX(IXP1(IX))-2.*FX(IX)+FX(IXM1(IX)))/DX**2
       LXF(IX)=CXF
    end do

    INTGD(1)=0.
    do IX=2,NX
       CXF=LXF(IX)
       LLOGF=CXF*(LOG(ABS(CXF))-1.)
       CXP=XRES+CXF
       CXM=XRES-CXF
       call WINTERP1D(X,NX,FXPP,CXP,FXP)
       call WINTERP1D(X,NX,FXPP,CXM,FXM)
       INTGD(IX)=LLOGF*(FXP-FXM)
    end do

    call EQSIMPSON1D(NX,LXF,INTGD,PINT)

  end subroutine CAUCHY_PPART

!
!*************************************************************************
!
  subroutine WINTERP1D(X,LX,FX,X0,FX0)
    implicit none
    integer, intent(IN):: LX
    real, dimension(LX), intent(IN):: X,FX
    real, intent(IN):: X0
    real, intent(inout):: FX0

    real:: FXY0,XMIN,XMAX,DX
    real:: P,A1,A2,A3
    integer:: I

    FX0=0.
    XMIN=X(1)
    XMAX=X(LX)
    DX=(XMAX-XMIN)/(LX-1)
    if ((X0>=XMIN).and.(X0<=XMAX)) then
       I=int(1.+(X0-XMIN)/DX)
       P=(X0-X(I))/DX
       if (I.eq.LX) I=LX-1
       if (I.eq.1) I=2
       A1=0.5*P*(P-1.)
       A2=1.-P*P
       A3=0.5*P*(P+1.)
       FX0=A1*FX(I-1)+A2*FX(I)+A3*FX(I+1)
    end if
  end subroutine WINTERP1D
!
!*************************************************************************
!
  subroutine WINTERP1D_2(X,LX,FX,X0,FX0)
    implicit none
    integer, intent(IN):: LX
    real, dimension(LX), intent(IN):: X,FX
    real, intent(IN):: X0
    real, intent(inout):: FX0

    real:: FXY0,XMIN,XMAX,DX
    real:: P,A1,A2,A3
    integer:: I

    FX0=0.
    XMIN=X(1)
    XMAX=X(LX)
    DX=(XMAX-XMIN)/(LX-1)
    if ((X0>=XMIN).and.(X0<=XMAX)) then
       I=int(1.+(X0-XMIN)/DX)
       P=(X0-X(I))/DX
!       if (I.eq.LX) I=LX-1
!       if (I.eq.1) I=2

       if(i .ne. 1 .and. i .ne. lx)then
          A1=0.5*P*(P-1.)
          A2=1.-P*P
          A3=0.5*P*(P+1.)
          FX0=A1*FX(I-1)+A2*FX(I)+A3*FX(I+1)
       else
          if(i .eq. lx) fx0 = fx(i)
          if(i .eq. 1)  fx0 = fx(i) + (fx(i+1) - fx(i)) * p
       end if

    end if
  end subroutine WINTERP1D_2

!
!*************************************************************************
!



!end module METS2AORSA

