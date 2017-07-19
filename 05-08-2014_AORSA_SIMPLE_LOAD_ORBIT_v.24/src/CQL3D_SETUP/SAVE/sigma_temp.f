!

    subroutine GETNONMAXSIGMA_AORSA_NEW(W,ZSPEC,ASPEC,DENS,BMAG, &
       & K1,XI1,JNXI1,K2,XI2,JNXI2,NBESSJ,ENORM,UPARMIN,UPARMAX, &
       & NUPAR,NUPER,UPER,UPAR,DFDUPER,DFDUPAR,WSPEC,IFAIL, &
	 & l_first, l_interp, kperp_max, nkperp)
	 
    implicit none
    !new declarations for interpolation subroutine
    logical, optional :: l_interp
    logical, optional :: l_first
    integer, optional :: nkperp
    real, optional :: kperp_max
    real :: dk, kperp, p
    integer :: nk
    real, save, dimension(1:200, -1:101, 0:4, 3, 3) :: perp_int   ! fix indices
    ! first is vpar, second is nkerp  (-1 and plus one added for quadratic interp
    !end new declerations
    real, intent(IN) :: W,ZSPEC,ASPEC,DENS,BMAG
    real, dimension(3) :: K1,K2
    integer, intent(IN) :: NUPAR,NUPER,NBESSJ
    real, dimension(NUPER), intent(IN) :: XI1,XI2
    real, dimension(NUPER, NBESSJ), intent(IN) :: JNXI1,JNXI2
    real, intent(IN) :: ENORM,UPARMIN,UPARMAX
    real, dimension(NUPER), intent(IN) :: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN) :: DFDUPER,DFDUPAR
    real, save ::  DFDTH(100,200)  ! put in module
    complex, dimension(3,3), intent(OUT) :: WSPEC
    integer, intent(OUT) :: IFAIL
    complex :: BETAFACT
    real, parameter:: EOVERAMU = 9.64853e7
    real, parameter:: WP2FACT=1.745915
    real, parameter:: MPC2=938271998.38
    real, parameter:: C=2.99792458e8
    real, parameter:: PI=3.141592653597932384
    real, dimension(NUPER):: JN0XI1,JNP1XI1,JNM1XI1
    real, dimension(NUPER):: JN0XI2,JNP1XI2,JNM1XI2
    real, dimension(NUPER,3,3):: SSGG
    real, dimension(NUPAR,3,3):: SGG
    real :: ga_33(nupar)
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_31
    real, dimension(NUPAR,-1:NBESSJ-2) :: sgg_11
    real, dimension(NUPAR) :: sgg_31_last, sgg_21p
    real, dimension(3,3):: THETARE,THETAIM
    real :: W2,WP2,WPFACT,WCW,RRP,RRM,RR,IRR, du_3, g_33j
    real, save :: g_33 
    real :: MUT0,SQMUT0,BETA1,BETA2,KPARA1,NPARA1
    real :: ISQ2,CUPAR,CUPER,NWCW,DFACTPAR,DFACTPER,LF0,LNF0UPER
    real :: DFACTPERI
    real :: UPAR0,SFACT0,SFACTP1,SFACTM1, factor, kperp2, vperp_norm2
    real :: ssgg_11,ssgg_22,ssgg_31,ssgg_33,ssgg_32,ssgg_21,du, temp
    integer :: NHARM,IHARM,NJ,NJP1,NJM1,J,K,M,N
    real :: ssgg_11a,ssgg_33a, ssgg_31a
    

    !fix--need to add bounds checking on stuff
    
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
    kperp = sqrt(kperp2)
    
    vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp2  ! need to fix kperp = 0

    if(.not. l_interp) then


!      -------------------------------------------------------
!      do extra integral (STIX) to convert all elements to "U"
!      do perpendicular integral by inline trapezoidal rule
!      and parallel integral by simpson's rule
!      -------------------------------------------------------

       do J = 1, NUPAR
	 g_33j = 0.0
         do K = 2, NUPER - 1
            g_33j = g_33j +                                                 &
    &         UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
         end do !nperp loop for 33 integral
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

       !  The extra uperp integrals are done; now do parallel integral by simpson's rule
       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)

       !  The derivative combination in the conductivity can be precomputed
       do j = 1, nupar
          do k = 1 , nuper
	       DFDTH(k,j) = upar(j) * DFDUPER(k,j) - uper(k) * DFDUPAR(k,j)
          end do
       end do


       !      -------------------
       !      Loop over harmonics
       !      -------------------

       do IHARM = -1, NHARM
          !do preliminaries for harmonic loop
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
          du = (uper(nuper) - uper(1))/(nuper - 1)
          ! -- Build array with integrand -- !
	
          do J=1,NUPAR
             ssgg_11 = 0.0
             ssgg_31 = 0.0
             ssgg_33 = 0.0
             do K = 2, NUPER - 1  !do center points for perp integrals
                LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
                ssgg_33 = ssgg_33 + JN0XI1(K)**2 *LF0
                LF0 = UPER(K) * JNP1XI1(K) * LF0
                ssgg_11 = ssgg_11 + UPER(K) * JNP1XI1(K)* LF0
                ssgg_31 = ssgg_31 + UPAR(J) * JN0XI1(K) * LF0
      	 end do
             ssgg_11 = ssgg_11 * 0.5
	       ssgg_31 = ssgg_31 * ISQ2
	
 	    k = 1   !do first point for perp integrals
          LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
          ssgg_33 = ssgg_33 + 0.5 * JN0XI1(K)**2 *LF0
          LF0 = UPER(K) * JNP1XI1(K) * LF0
          ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0


	    k = nuper  ! do last point for perp integrals
	    LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
	    ssgg_33 = ssgg_33 + 0.5 * LF0 * JN0XI1(K)**2
	    LF0 = UPER(K) * JNP1XI1(K) * LF0
	    ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0


       !  end of the kperp loop


	
	    du = (uper(nuper) - uper(1))/(nuper - 1)  !could be taken out of loop
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
	    ! -- Form the itegrands for the non-singular parallel integral
	
          do J=1,NUPAR
             RR = 1.- NWCW - DFACTPER * UPAR(J)
             IRR = 1. / RR
             do M=1,3
                do N=1,M
                   SGG(J,M,N) = SGG(J,M,N) * IRR
                end do
             end do
          end do




      	!-- Hermitian part of $\Theta/(2\pi)$ --  use trap. integration!

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

  !       Use symmetry to fill the upper left triangle 	
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
	
	    do N=1,3  !do loop to increment w by harmonic contrib
      	 do M=1,3
                WSPEC(M,N)=WSPEC(M,N)+                                &
                2.0d0*PI*WPFACT*BETAFACT*cmplx(THETARE(M,N),THETAIM(M,N))
      	 end do
          end do

	 end if ! skip if nharm = zero 
    end do !harmonic sum
 
   ! add in extra term in sig33 from the w to u conversion

    WSPEC(3,3)=WSPEC(3,3)+2.*PI*WPFACT*BETAFACT*cmplx(g_33,0.0)
    
	    print*, 'g_33 exact  ', g_33


    call WROTATE_AORSA(BETA1,BETA2,WSPEC)
    
    return  !finished with orgiginal non iteration version
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
    else !going to interpolate
    
       if(l_first) then   !  need to do perp_int's for each vpar and nharm 0->plus
          !   do extra integral (STIX) to convert all elements to "U" form)
          !   do perpendicular integral by inline trapezoidal rule
          du = (uper(nuper) - uper(1))/(nuper - 1)
          do J = 1, NUPAR
	    g_33j = 0.0
             do K = 2, NUPER - 1  !center perp points
                g_33j = g_33j +                                                 &
             &     UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
             end do
	          k = 1
                g_33j = g_33j +                                                    &
             &     0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))
	          k = nuper
                g_33j = g_33j +                                                    &
             &     0.5 * UPAR(J)*(UPER(K)*DFDUPAR(K,J)-UPAR(J)*DFDUPER(K,J))

	
	  
	  ga_33(j) = g_33j * du

       end do

       g_33  = 0.0

!      The extra uperp integral is done; now do parallel integral by simpson's rule
!
       call EQSIMPSON1D_2(NUPAR, UPAR, ga_33, g_33)
	 
	 print*, 'g_33 from init  ',  g_33

!      The derivative combination in the conductivity can be precomputed
       do j = 1, nupar
          do k = 1 , nuper
	       DFDTH(k,j) = upar(j) * DFDUPER(k,j) - uper(k) * DFDUPAR(k,j)
          end do
       end do


!      -------------------
!      loop over kperp space
!      -------------------
       do nk= 0, nkperp+1  !loop to fill interpolation table
           dk = kperp_max/nkperp
	     kperp = dk*nk
	     if(nk == 0) then 
	        kperp = 1.0E-7  ! don't let kperp get too small--I'll put exact in
	     end if
	     vperp_norm2 = mut0*wcw*wcw*w2/c/c/kperp/kperp
	     !need to get Bessel functions at each kperp loop
	     call WMATPRECALC_AORSA(ZSPEC, ASPEC, ENORM, BMAG, KPERP, UPER,  &
           &  NUPER, NBESSJ, 2*NBESSJ+8,XI1, JNXI1, IFAIL)
	     if(ifail .ne. 0) then
	        stop 'bessel function failure in interp intialize'
	     end if

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


          do K = 2, NUPER - 1  !do center points for perp integrals
             LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
             ssgg_33 = ssgg_33 + JN0XI1(K)**2 *LF0
             LF0 = UPER(K) * JNP1XI1(K) * LF0
             ssgg_11 = ssgg_11 + UPER(K) * JNP1XI1(K)* LF0
             ssgg_31 = ssgg_31 + UPAR(J) * JN0XI1(K) * LF0
          end do

          ssgg_11 = ssgg_11 * 0.5
	    ssgg_31 = ssgg_31 * ISQ2
	
	
	    k = 1   !do first point for perp integrals
          LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
          ssgg_33 = ssgg_33 + 0.5 * JN0XI1(K)**2 *LF0
          LF0 = UPER(K) * JNP1XI1(K) * LF0
          ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0


	    k = nuper  ! do last point for perp integrals
	    LF0 = DFDUPER(K,J) - DFACTPER * dfdth(k,j)
	    ssgg_33 = ssgg_33 + 0.5 * LF0 * JN0XI1(K)**2
	    LF0 = UPER(K) * JNP1XI1(K) * LF0
	    ssgg_11 = ssgg_11 + 0.25 *       UPER(K) * JNP1XI1(K)* LF0
          ssgg_31 = ssgg_31 + 0.5 * ISQ2 * UPAR(J) * JN0XI1(K) * LF0


       !  end of the kperp loop


	
	    du = (uper(nuper) - uper(1))/(nuper - 1)  !could be taken out of loop
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
             SGG(1:nupar,2,1) = -0.5*(SGG(1:nupar,1,1)+ SGG(1:nupar,2,2))  &
	         &  + vperp_norm2*iharm*iharm*(sgg_21p(1:nupar))  ! use bessel function identity for 21
	     end select
	     
	     
	     SGG(1:nupar,3,2) = sgg_31(1:nupar,iharm-1)
	     SGG(1:nupar,3,1) = sgg_31(1:nupar,iharm)

	     ! the perpendicular integrals are now saved for each vpar, each iharm
	     perp_int(1:nupar, nk, iharm, 1:3,1:3) = SGG(1:nupar,1:3,1:3)
         end do !harmonic sum
       end do !kperp fill loop
	 ! put in -dk point 2x2, 3,3 are even, xz, yz are odd
	 perp_int(1:nupar, -1, iharm, 1:3,1:3) = perp_int(1:nupar, 1, iharm, 1:3,1:3)
	 perp_int(1:nupar, -1, iharm, 1:3,1:3) = perp_int(1:nupar, 1, iharm, 1:3,1:3)
	 perp_int(1:nupar, -1, iharm, 3,1:2) = perp_int(1:nupar, 1, iharm, 3,1:2)
	 perp_int(1:nupar, -1, iharm, 1:2,3) = perp_int(1:nupar, 1, iharm, 1:2,3)
	 
	 
!---------------------------------------------------------------------------------------------------------------	 
	 
!---------------------------------------------------------------------------------------------------------------	 
	 
!---------------------------------------------------------------------------------------------------------------	 
	 
    else   !else for not here the first time at a point--does real calcualtion
    kperp = sqrt(kperp2)
    dk = kperp_max/nkperp
    nk = kperp/dk
    du_3 = (UPAR(NUPAR)-UPAR(1))/real(nupar-1)/3.
    p = kperp - nk*dk
    

    
!      -------------------
!      Loop over harmonics
!      -------------------
      do IHARM = 0, NHARM
         NWCW=real(IHARM)*WCW  !L*omega_cyc
         BETAFACT=exp(cmplx(0.,IHARM*(BETA1-BETA2)))
         NJ=ABS(IHARM)+1
         NJP1=ABS(IHARM+1)+1
         NJM1=ABS(IHARM-1)+1
         SFACT0=real(sign(1,IHARM))**ABS(IHARM)
         SFACTP1=real(sign(1,IHARM+1))**ABS(IHARM+1)
         SFACTM1=real(sign(1,IHARM-1))**ABS(IHARM-1)

         ! -- Build array with interpolation -- !
         
         SGG(1:nupar,1:3,1:3) =   &
           &   0.5*perp_int(1:nupar,nk-1,iharm,1:3,1:3)*(p*(p-1.))  &
           &   + perp_int(1:nupar,nk,iharm,1:3,1:3)*(1.-p*p)   &
	     &   + 0.5*perp_int(1:nupar,nk+1,iharm,1:3,1:3)*(p*(p+1.))
         ! -- Resonance relation -- !
         RRP=1. - NWCW - DFACTPER * UPARMAX
         RRM=1. - NWCW - DFACTPER * UPARMIN
         if (RRP*RRM.GT.0) then
            ! -- No resonance here -- !
	      ! -- Form the itegrands for the non-singular parallel integral
	
            do J=1,NUPAR
               RR = 1.- NWCW - DFACTPER * UPAR(J)
               IRR = 1. / RR
               do M=1,3
                  do N=1,M
                     SGG(J,M,N) = SGG(J,M,N) * IRR
                  end do
               end do
            end do

      	!-- Hermitian part of $\Theta/(2\pi)$ --  use trap. integration!
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
            !  Use symmetry to fill the upper left triangle 	
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
		  
		  

              !fill symmetry positive iharm real
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
              !fill symetry positive iharm
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

                !fill real, negative iharm, no resonance
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

                  !fill symmetry for real theta and negative iharm
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
			
                  ! fill symmetry positions for negative iharm and imaginary theta
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
	    
	    print*, 'g_33 interp  ', g_33

          call WROTATE_AORSA(BETA1,BETA2,WSPEC)
       end if !interpolation prep, interpolation
    end if !interpolation yes/no
    end subroutine GETNONMAXSIGMA_AORSA_NEW
