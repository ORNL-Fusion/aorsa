!module METS2AORSA

!  private
!  public:: WMATPRECALC_AORSA,GETNONMAXSWMAT_AORSA,GETNONMAXSWMAT_AORSA_NEW

!contains

!
!*************************************************************************
!

  subroutine WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER,UPER, &
       &                   NUPER,NBESSJ,NBESSJ_START,XI,JNXI,IFAIL)
    implicit none
    real, intent(IN):: ZSPEC, ASPEC, ENORM,BMAG,KPER
    integer, intent(IN):: NUPER
    real, dimension(NUPER), intent(IN):: UPER
    integer, intent(IN):: NBESSJ,NBESSJ_START
    real, dimension(NUPER), intent(OUT):: XI
    real, dimension(NUPER,NBESSJ), intent(OUT):: JNXI
    integer, intent(OUT):: IFAIL

    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EVPERAMU=9.314943e8
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real:: WC,MUT0,SQMUT0

    MUT0=0.5*MPC2*ASPEC/ENORM
    SQMUT0=SQRT(MUT0)
    WC=BMAG*ZSPEC*EOVERAMU/ASPEC
    XI(:)=KPER*UPER(:)*C/SQMUT0/WC

    call BESSJ(XI(1),JNXI(1,1),NUPER,NBESSJ,NBESSJ_START)

  end subroutine WMATPRECALC_AORSA

!
!*************************************************************************
!

  subroutine GETNONMAXSIGMA_AORSA_NEW(W,ZSPEC,ASPEC,DENS,BMAG, &
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
    real, parameter:: EOVERAMU = 9.64853e7
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

    logical, parameter :: use_original  = .false.
    logical, parameter :: use_gemm  = .false.
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
!dir$ ssp_private CAUCHY_PPART2

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
     
    

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERAMU/ASPEC/W
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
    if (use_original) then

        do J = 1, NUPAR
	  g_33j = 0.0

          do K = 2, NUPER - 1
             g_33j = g_33j +                                                 &
     &         UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
          end do


	  k = 1
          g_33j = g_33j +                                                    &
     &      0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))


	  k = nuper
          g_33j = g_33j +                                                    &
     &      0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))

	
	  du = (uper(nuper) - uper(1))/(nuper - 1)
	  ga_33(j) = g_33j * du

       end do

       g_33  = 0.0

!      the uperp integral is done; now do parallel integral by simpson's rule

       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)


       do j = 1, nupar
       do k = 1 , nuper
	     DFDTH(k,j) = upar(j) * DFDUPER(k,j) - uper(k) * DFDUPAR(k,j)
       end do
       end do

       

     else
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

!      the uperp integral is done; now do parallel integral by simpson's rule

       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)



     endif


     if (.not.use_original) then


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

!     -------------------------------
!     precompute LF0kj(k,j) * vectors
!     -------------------------------
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

!	-----------------------------
!	one big happy matrix multiply
!	-----------------------------
        if (use_gemm) then

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

        else

!       LF0JN0N0(1:jdim,1:nbessj) = matmul(                    &
!                         transpose(LF0kj(2:(kdim-1),1:jdim)), &       
! 			JN0N0(2:(kdim-1),1:nbessj) )
!       LF0JN0N0U(1:jdim,1:nbessj) = matmul(               &
!                         transpose(LF0kj(2:(kdim-1),1:jdim)), &
!                         JN0N0U(2:(kdim-1),1:nbessj) )
!       LF0JN0N1U(1:jdim,1:(nbessj-1)) = matmul(               &
!                         transpose(LF0kj(2:(kdim-1),1:jdim)), &
!                         JN0N1U(2:(kdim-1),1:(nbessj-1))  )

         Jmat(1:jdim,1:(3*nbessj-1)) = matmul(                  &
             transpose( LF0kj(2:(kdim-1), 1:jdim)),             &
                        Kmat(2:(kdim-1),1:(3*nbessj-1))  )

        endif
     
       endif
     


      maxabserr = 0.0
      maxrelerr = 0.0


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


       if (use_original) then
       JN0XI1(1:NUPER)=SFACT0*JNXI1(1:NUPER,NJ)
       JNP1XI1(1:NUPER)=SFACTP1*JNXI1(1:NUPER,NJP1)
       JNM1XI1(1:NUPER)=SFACTM1*JNXI1(1:NUPER,NJM1)
       endif

       ! -- Build array with integrand -- !


	
       do J=1,NUPAR
	
          ssgg_11 = 0.0
          ssgg_31 = 0.0
          ssgg_33 = 0.0

          if (use_original) then

          do K = 2, NUPER - 1
	
             LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
             ssgg_33 = ssgg_33 + JN0XI1(K)**2 *LF0
		
             LF0 = UPER(K) * JNP1XI1(K) * LF0
             ssgg_11 = ssgg_11 + UPER(K) * JNP1XI1(K)* LF0
             ssgg_31 = ssgg_31 + UPAR(J) * JN0XI1(K) * LF0
          end do

          else
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

          endif

!

	
          ssgg_11 = ssgg_11 * 0.5
	    ssgg_31 = ssgg_31 * ISQ2
	
	
          if (use_original) then
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


          else
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



          endif

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
	
          if (use_original) then
          do J=1,NUPAR
             RR = 1.- NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                end do
             end do
          end do

          else
          do J=1,NUPAR
             RR = 1.- NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             rr_j(j) = rr
             irr_j(j) = irr
          enddo


!dir$     concurrent
          do mn=1,6
           m = mlist(mn)
           n = nlist(mn)

           do j=1,nupar
            IRR = irr_j(j)
	    SGG(J,M,N) = SGG(J,M,N) * IRR
           enddo

          end do

          endif
	

          !-- Hermitian part of $\Theta/(2\pi)$ --!
	
          if (use_original) then
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
          else

!dir$     concurrent,preferstream,novector
            do mn=1,6
               m = mlist(mn)
               n = nlist(mn)
!dir$ vector             
               sum_even = sum(sgg(2:(nupar-1):2,m,n))
               sum_odd  = sum(sgg(3:(nupar-1):2,m,n))
             
               thetare(m,n) = (sgg(1,m,n) + sgg(nupar,m,n) +    &
                   2.0d0*(sum_odd + sum_even) + 2.0d0*sum_even ) * du_3
            enddo
           endif
	
	
          if (use_original) then
	  do J=1,NUPAR ! restore sgg
	     RR = 1.- NWCW - DFACTPER * UPAR(J)
	
             do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * RR
                end do
             end do
	
          end do
          else

!dir$        concurrent
             do mn=1,6
              m = mlist(mn)
              n = nlist(mn)

              do j=1,nupar
                   SGG(J,M,N) = SGG(J,M,N) * rr_j(j)
              enddo
             end do

          endif
	
 	  THETARE(1,2) = THETARE(2,1)
          THETARE(1,3) = THETARE(3,1)
          THETARE(2,3) = THETARE(3,2)

          !-- Anti-hermitian part of $\Theta/(2\pi)$ --!
          THETAIM(1:3,1:3)=0.

       else

          ! -- There is a resonance all right -- !
          ! -- 1. Hermitian part -- !
          UPAR0 = (1. - NWCW) * DFACTPERI
	
          if (use_original) then
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))
          else
!dir$       concurrent
            do mn=1,6
             m = mlist(mn)
             n = nlist(mn)
             call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,m,n),THETARE(m,n))
            enddo
          endif
            

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
	
          if (use_original) then
      	  do J=1,NUPAR

             RR = 1. - NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             do M=1,3
                do N = 1,M
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                end do
             end do
      	  end do
          else
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

          endif
		
          !-- Hermitian part of $\Theta/(2\pi)$ --!

          if (use_original) then
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
          else
!dir$     concurrent,preferstream,novector

             do mn=1,6
!dir$     vector
               m = mlist(mn)
               n = nlist(mn)
               
               sum_even = sum(sgg(2:(nupar-1):2,m,n))
               sum_odd  = sum(sgg(3:(nupar-1):2,m,n))

               thetare(m,n) = (sgg(1,m,n)+sgg(nupar,m,n) +     &
                  2.0d0*(sum_even + sum_odd) + 2.0d0*sum_even) * du_3
             enddo
           endif
	
               


		
          if (use_original) then
          do J=1,NUPAR  ! restore sgg
             RR=1.- NWCW - DFACTPER * UPAR(J)
	     do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * RR
            	end do
             end do
      	  end do
          else

!dir$        concurrent
             do mn=1,6
               m = mlist(mn)
               n = nlist(mn)

               do j=1,nupar
                   SGG(J,M,N) = SGG(J,M,N)  * rr_j(j)
               enddo

             end do
          endif

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

         if (use_original) then
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
      	 call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))
         else

!dir$       concurrent
            do mn=1,6
              m = mlist(mn)
              n = nlist(mn)
      	      call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,m,n),THETARE(m,n))
            enddo
         endif


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
                WSPEC(M,N)=WSPEC(M,N)+                                &
                2.0d0*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
      	 end do
         end do

	 end if ! zero skip if
    end do !harmonic sum

    ! add in extra term in sig33 from the w to u conversion

    WSPEC(3,3)=WSPEC(3,3)+2.*PI*WPFACT*BETAFACT*cmplx(g_33,0.0)

    call WROTATE_AORSA(BETA1,BETA2,WSPEC)

  end subroutine GETNONMAXSIGMA_AORSA_NEW






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
    complex, dimension(3,3), intent(OUT):: WSPEC
    integer, intent(OUT):: IFAIL

    complex:: BETAFACT
    real, parameter:: EOVERAMU = 9.64853e7
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

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERAMU/ASPEC/W
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
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,1),THETARE(1,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,1),THETARE(2,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,1),THETARE(3,1))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,2),THETARE(1,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,2),THETARE(2,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,2),THETARE(3,2))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,1,3),THETARE(1,3))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,2,3),THETARE(2,3))
          call CAUCHY_PPART2(UPAR,NUPAR,UPAR0,SGG(1,3,3),THETARE(3,3))

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
    complex, dimension(3,3), intent(OUT):: WSPEC
    integer, intent(OUT):: IFAIL

    complex:: BETAFACT
    real, parameter:: EOVERAMU = 9.64853e7
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

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERAMU/ASPEC/W
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

    subroutine GETNONMAXSIGMA_AORSA(W,ZSPEC,ASPEC,DENS,BMAG, &
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
    real, parameter:: EOVERAMU = 9.64853e7
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

    IFAIL=0
    WSPEC(1:3,1:3) = cmplx(0.,0.)
    W2=W*W
    WP2=DENS*ZSPEC**2*WP2FACT/ASPEC
    WPFACT=WP2/W2
    WCW=BMAG*ZSPEC*EOVERAMU/ASPEC/W
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

  end subroutine GETNONMAXSIGMA_AORSA

!
!*************************************************************************
!




  subroutine BESSJ(XI,JNXI,LXI,LNBESSJ,LNSBESSJ)
    implicit none
    integer, intent(IN):: LXI,LNBESSJ,LNSBESSJ
    real, dimension(LXI), intent(IN):: XI
    real, dimension(LXI,LNBESSJ), intent(OUT):: JNXI
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
    real, dimension(LXI), intent(OUT) :: J0XI
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
    real, dimension(LXI), intent(OUT):: J1XI
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
    real, intent(OUT):: INTEGRAL
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
    real, intent(OUT):: INTEGRAL

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
    real, intent(OUT):: INTEGRAL

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
    real, intent(OUT):: FXY0
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
    real, intent(OUT):: PINT

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
    real, intent(OUT):: FX0

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
    real, intent(OUT):: FX0

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


       subroutine intplt1d(x, lx, fx, x0, fx0)

!      -----------------------
!      1D linear interpolation
!      -----------------------

       implicit none

       integer lx, i
       real fx(lx), x(lx), fx0, x0, dx

       dx = (x(lx) - x(1)) / (lx - 1)
       i = int((x0 - x(1)) / dx) + 1

       if(i .eq. lx)then
          fx0 = fx(lx)
       else
          fx0 = fx(i) + (fx(i+1) - fx(i)) * (x0 - x(i)) / dx
       end if

	 return
	 end subroutine intplt1d

!
!*************************************************************************
!



!end module METS2AORSA

