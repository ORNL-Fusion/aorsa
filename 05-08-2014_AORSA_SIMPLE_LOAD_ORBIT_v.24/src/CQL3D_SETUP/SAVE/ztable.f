!     Last change:  DNS  20 Jun 2002    3:53 pm
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      MODULE ZTABLE_MOD
!
!        USE KINDS
         IMPLICIT NONE
!
         INTEGER  , PARAMETER :: NTABLES = 2
!
         INTEGER  , PARAMETER :: NZEDS   = 200
!
         INTEGER  , PARAMETER :: NALPS   = 200
!
         REAL :: TWOPI,REALNZEDS,REALNALPS
!
         REAL :: ZED(0:NZEDS),ALP(0:NALPS)
         REAL :: ZFUNC_REAL(0:4,0:NZEDS,0:NALPS,NTABLES), &
                      ZFUNC_IMAG(0:4,0:NZEDS,0:NALPS,NTABLES)
!
         CHARACTER :: ATABLE*10
!
      END MODULE ZTABLE_MOD
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_ZFUNC(ZETA,ALPHA,GAMMA,ZFUNC)
!         USE KINDS
         IMPLICIT NONE
         REAL, INTENT(IN)  :: ZETA,ALPHA,GAMMA
         COMPLEX, INTENT(OUT) :: ZFUNC(0:4)
         REAL    :: ALPHAINV1,ALPHAINV2,ALPHAINV3,ALPHAINV4
         REAL    :: PHASE,COLLISION1,COLLISION2
         COMPLEX :: EYE,PHASE_FACTOR
         COMPLEX :: ZFUNC_NONLOCAL1(0:4),ZFUNC_NONLOCAL2(0:4)
!
         CALL ZTABLE_VALUE(ZETA,ALPHA,ZFUNC,ZFUNC_NONLOCAL1,ZFUNC_NONLOCAL2)
         IF (ALPHA.LE.0.0) RETURN
!
         COLLISION1 = EXP(-0.125*ABS(GAMMA)/ALPHA**3)
         COLLISION2 = COLLISION1**8
!
         EYE = CMPLX(0.0,1.0)
!
         ALPHAINV1 = 1.0/ALPHA
         ALPHAINV2 = ALPHAINV1*ALPHAINV1
         ALPHAINV3 = ALPHAINV2*ALPHAINV1
         ALPHAINV4 = ALPHAINV3*ALPHAINV1
!
         ZFUNC_NONLOCAL1(4) =  ZFUNC_NONLOCAL1(4)                          + &
                              (ZFUNC_NONLOCAL1(3)*EYE)*(4.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL1(2)     *(6.0*ALPHAINV2) - &
                              (ZFUNC_NONLOCAL1(1)*EYE)*(4.0*ALPHAINV3) + &
                               ZFUNC_NONLOCAL1(0)     *         ALPHAINV4
         ZFUNC_NONLOCAL1(3) =  ZFUNC_NONLOCAL1(3)                          + &
                              (ZFUNC_NONLOCAL1(2)*EYE)*(3.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL1(1)     *(3.0*ALPHAINV2) - &
                              (ZFUNC_NONLOCAL1(0)*EYE)*         ALPHAINV3
         ZFUNC_NONLOCAL1(2) =  ZFUNC_NONLOCAL1(2)                          + &
                              (ZFUNC_NONLOCAL1(1)*EYE)*(2.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL1(0)     *         ALPHAINV2
         ZFUNC_NONLOCAL1(1) =  ZFUNC_NONLOCAL1(1)                          + &
                              (ZFUNC_NONLOCAL1(0)*EYE)*         ALPHAINV1
!
         PHASE = ZETA*ALPHAINV1
         PHASE_FACTOR = CMPLX(COS(PHASE),SIN(PHASE))
         ZFUNC_NONLOCAL1 = (ZFUNC_NONLOCAL1*PHASE_FACTOR)*COLLISION1
!
         ALPHAINV1 = 2.0/ALPHA
         ALPHAINV2 = ALPHAINV1*ALPHAINV1
         ALPHAINV3 = ALPHAINV2*ALPHAINV1
         ALPHAINV4 = ALPHAINV3*ALPHAINV1
!
         ZFUNC_NONLOCAL2(4) =  ZFUNC_NONLOCAL2(4)                          + &
                              (ZFUNC_NONLOCAL2(3)*EYE)*(4.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL2(2)     *(6.0*ALPHAINV2) - &
                              (ZFUNC_NONLOCAL2(1)*EYE)*(4.0*ALPHAINV3) + &
                               ZFUNC_NONLOCAL2(0)     *         ALPHAINV4
         ZFUNC_NONLOCAL2(3) =  ZFUNC_NONLOCAL2(3)                          + &
                              (ZFUNC_NONLOCAL2(2)*EYE)*(3.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL2(1)     *(3.0*ALPHAINV2) - &
                              (ZFUNC_NONLOCAL2(0)*EYE)*         ALPHAINV3
         ZFUNC_NONLOCAL2(2) =  ZFUNC_NONLOCAL2(2)                          + &
                              (ZFUNC_NONLOCAL2(1)*EYE)*(2.0*ALPHAINV1) - &
                               ZFUNC_NONLOCAL2(0)     *         ALPHAINV2
         ZFUNC_NONLOCAL2(1) =  ZFUNC_NONLOCAL2(1)                          + &
                              (ZFUNC_NONLOCAL2(0)*EYE)*         ALPHAINV1
!
         PHASE = ZETA*ALPHAINV1
         PHASE_FACTOR = CMPLX(COS(PHASE),SIN(PHASE))
         ZFUNC_NONLOCAL2 = (ZFUNC_NONLOCAL2*PHASE_FACTOR)*COLLISION2
!
         ZFUNC      = ZFUNC + ZFUNC_NONLOCAL1 + ZFUNC_NONLOCAL2
!
         RETURN
      END SUBROUTINE ZTABLE_ZFUNC
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ... same as above, with (ZETA-DELTA) instead of ZETA,
!     and individual values, instead of array.  Used in METS code.
!
      SUBROUTINE GET_ZFUNC_TABLE(ZETA,ALPHA,GAMMA,DELTA, &
                                 ZFUNC0,ZFUNC1,ZFUNC2,ZFUNC3,ZFUNC4)
!         USE KINDS
         IMPLICIT NONE
         REAL, INTENT(IN)  :: ZETA,ALPHA,GAMMA,DELTA
         COMPLEX, INTENT(OUT) :: ZFUNC0,ZFUNC1,ZFUNC2,ZFUNC3,ZFUNC4
         REAL   :: ZETAMDELTA
         COMPLEX :: ZFUNC(0:4)
         ZETAMDELTA = ZETA-DELTA
         CALL ZTABLE_ZFUNC(ZETAMDELTA,ALPHA,GAMMA,ZFUNC)
         ZFUNC0 = ZFUNC(0)
         ZFUNC1 = ZFUNC(1)
         ZFUNC2 = ZFUNC(2)
         ZFUNC3 = ZFUNC(3)
         ZFUNC4 = ZFUNC(4)
         RETURN
      END SUBROUTINE GET_ZFUNC_TABLE
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ... for future use
!
      SUBROUTINE ZTABLE_ZFUNC_PARTS(ZETA,ALPHA,GAMMA, &
                                    ZFUNC_LOCAL,ZFUNC_NONLOCAL1,ZFUNC_NONLOCAL2)
!         USE KINDS
         IMPLICIT NONE
         REAL, INTENT(IN)  :: ZETA,ALPHA,GAMMA
         COMPLEX, INTENT(OUT) :: ZFUNC_LOCAL(0:4), &
                                      ZFUNC_NONLOCAL1(0:4),ZFUNC_NONLOCAL2(0:4)
         REAL    :: COLLISION1,COLLISION2
         CALL ZTABLE_VALUE(ZETA,ALPHA,ZFUNC_LOCAL,ZFUNC_NONLOCAL1,ZFUNC_NONLOCAL2)
         IF (ALPHA.LE.0.0) RETURN
         COLLISION1 = EXP(-0.125*ABS(GAMMA)/ALPHA**3)
         COLLISION2 = COLLISION1**8
         ZFUNC_NONLOCAL1 = COLLISION1*ZFUNC_NONLOCAL1
         ZFUNC_NONLOCAL2 = COLLISION2*ZFUNC_NONLOCAL2
         RETURN
      END SUBROUTINE ZTABLE_ZFUNC_PARTS
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_VALUE(ZETA,ALPHA, &
                              ZFUNC_LOCAL,ZFUNC_NONLOCAL1,ZFUNC_NONLOCAL2)
!
         USE ZTABLE_MOD
!
!         USE KINDS
         IMPLICIT NONE
!
         REAL   , INTENT(IN)  :: ZETA,ALPHA
         COMPLEX, INTENT(OUT) :: ZFUNC_LOCAL(0:4), &
                                      ZFUNC_NONLOCAL1(0:4),ZFUNC_NONLOCAL2(0:4)
!
         INTEGER      :: K
         INTEGER      :: ITABLE
         INTEGER      :: IM,JM,IP,JP
         REAL    :: XM,YM,XP,YP
         REAL    :: ABSALPHA,ABSZETA,ALPVAL,ZEDVAL,ZETASCALE,ZETASIGN,ZSIGN
         REAL    :: ALPHAINV1,ALPHAINV2,ALPHAINV3,ALPHAINV4
         REAL    :: ZREAL,ZIMAG,ZFUNC_SCALE
!
         ABSALPHA  = ABS(ALPHA)
         ALPVAL    = SQRT(ABSALPHA/(1.0+ABSALPHA))
!
         ZETASIGN  = SIGN(1.0,ZETA)
         ZETASCALE = SQRT(1.0+ABSALPHA)
         ABSZETA   = ABS(ZETA)/ZETASCALE
         ZEDVAL    = ABSZETA/(1.0+ABSZETA)
!
         IM = INT(ZEDVAL*REALNZEDS)
         JM = INT(ALPVAL*REALNALPS)
         IP = IM+1
         JP = JM+1
!
         XP = ZEDVAL*REALNZEDS-REAL(IM)
         YP = ALPVAL*REALNALPS-REAL(JM)
         XM = 1.0 - XP
         YM = 1.0 - YP
!
         ITABLE = 1
         DO K=0,4
            ZREAL = XM*YM*ZFUNC_REAL(K,IM,JM,ITABLE) + XM*YP*ZFUNC_REAL(K,IM,JP,ITABLE) + &
                    XP*YM*ZFUNC_REAL(K,IP,JM,ITABLE) + XP*YP*ZFUNC_REAL(K,IP,JP,ITABLE)
            ZIMAG = XM*YM*ZFUNC_IMAG(K,IM,JM,ITABLE) + XM*YP*ZFUNC_IMAG(K,IM,JP,ITABLE) + &
                    XP*YM*ZFUNC_IMAG(K,IP,JM,ITABLE) + XP*YP*ZFUNC_IMAG(K,IP,JP,ITABLE)
            ZSIGN = ZETASIGN**K
            ZREAL = ZREAL*ZSIGN*ZETASIGN
            ZIMAG = ZIMAG*ZSIGN
            ZFUNC_LOCAL(K) = CMPLX(ZREAL,ZIMAG)
            ZFUNC_SCALE    = SQRT(ZETA**2+REAL(K+1)**2*ZETASCALE**2)
            ZFUNC_LOCAL(K) = ZFUNC_LOCAL(K)/ZFUNC_SCALE**(K+1)
            ZFUNC_NONLOCAL1(K) = 0.0
            ZFUNC_NONLOCAL2(K) = 0.0
         ENDDO
!
         IF (ALPHA.LE.0.0) RETURN
!
         ZFUNC_SCALE = 2.0*EXP(-0.0625/ALPHA**2)
         ITABLE = 2
         DO K=0,4
            ZFUNC_NONLOCAL1(K) = ZFUNC_SCALE*CMPLX(0.0,AIMAG(ZFUNC_LOCAL(K)))
            ZREAL = XM*YM*ZFUNC_REAL(K,IM,JM,ITABLE) + XM*YP*ZFUNC_REAL(K,IM,JP,ITABLE) + &
                    XP*YM*ZFUNC_REAL(K,IP,JM,ITABLE) + XP*YP*ZFUNC_REAL(K,IP,JP,ITABLE)
            ZIMAG = XM*YM*ZFUNC_IMAG(K,IM,JM,ITABLE) + XM*YP*ZFUNC_IMAG(K,IM,JP,ITABLE) + &
                    XP*YM*ZFUNC_IMAG(K,IP,JM,ITABLE) + XP*YP*ZFUNC_IMAG(K,IP,JP,ITABLE)
            ZSIGN = ZETASIGN**K
            ZREAL = ZREAL*ZSIGN*ZETASIGN
            ZIMAG = ZIMAG*ZSIGN
            ZFUNC_NONLOCAL2(K) = CMPLX(ZREAL,ZIMAG)
            ZFUNC_NONLOCAL2(K) = ZFUNC_NONLOCAL2(K)
            ZFUNC_LOCAL(K)     = CONJG(ZFUNC_LOCAL(K)) + ZFUNC_NONLOCAL2(K)
            ZFUNC_NONLOCAL2(K) = -CONJG(ZFUNC_NONLOCAL2(K))
         ENDDO
!
         RETURN
!
      END SUBROUTINE ZTABLE_VALUE
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_SETUP( MYID, IFLAG )
!
         USE ZTABLE_MOD
!
!         USE KINDS
         IMPLICIT NONE
!
         INTEGER :: I, MYID, IFLAG
         LOGICAL :: LEXISTING_TABLE
!
         TWOPI = 8.0*ATAN(1.0)
         REALNZEDS = REAL(NZEDS)
         REALNALPS = REAL(NALPS)



!
         ATABLE = 'ZTABLE.TXT'
         INQUIRE (FILE=ATABLE,EXIST=LEXISTING_TABLE)
         IF (LEXISTING_TABLE) THEN
            CALL ZTABLE_READ
         ELSE
            IF ( MYID .EQ. 0 ) THEN

            WRITE (*,*) '*************************************************************'
            WRITE (*,*) '*** File ZTABLE.TXT, the Z-function table, was not found. ***'
            WRITE (*,*) '*** Recomputing now.  This may take a few extra minutes.  ***'
            WRITE (*,*) '*************************************************************'
            CALL ZTABLE_WRITE

            ENDIF

            IFLAG = 1
         ENDIF
!
         RETURN
!
      END SUBROUTINE ZTABLE_SETUP
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_READ
!
         USE ZTABLE_MOD
!
!         USE KINDS
         IMPLICIT NONE
!
         CHARACTER :: AFMT*128
         INTEGER   :: I,J,K,ITABLE
!
         OPEN (UNIT=12,FILE=ATABLE)
!
         DO ITABLE=1,NTABLES
!
            AFMT = "(13X,//)"
            READ (12,AFMT)
!
            DO K=0,4
!
               AFMT = "(13X)"
               READ (12,AFMT)
!
               AFMT = "(/13X)"
               READ (12,AFMT)
               AFMT = "(13X,999(1PE13.5))"
               READ (12,AFMT) (ALP(J),J=0,MIN(NALPS,999))
               AFMT = "(1PE13.5,999(1PE13.5))"
               DO I=0,NZEDS
                  READ (12,AFMT) ZED(I),(ZFUNC_REAL(K,I,J,ITABLE),J=0,MIN(NALPS,999))
               ENDDO
!
               AFMT = "(//13X)"
               READ (12,AFMT)
               AFMT = "(13X,999(1PE13.5))"
               READ (12,AFMT) (ALP(J),J=0,MIN(NALPS,999))
               AFMT = "(1PE13.5,999(1PE13.5))"
               DO I=0,NZEDS
                  READ (12,AFMT) ZED(I),(ZFUNC_IMAG(K,I,J,ITABLE),J=0,MIN(NALPS,999))
               ENDDO
!
               AFMT = "(///)"
               READ (12,AFMT)
!
            ENDDO
!
         ENDDO
!
         CLOSE (UNIT=12)
!
         RETURN
!
      END SUBROUTINE ZTABLE_READ
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_WRITE
!
         USE ZTABLE_MOD
!
!         USE KINDS
         IMPLICIT NONE
!
         INTEGER      :: I,J,K,ITABLE
         REAL    :: ZETA,ALPHA,ZETASCALE,ABSZETA
         REAL    :: ASYMPT(0:4) = (/-1.0,1.0,-2.0,6.0,-24.0/)
         COMPLEX :: ZFUNC
         CHARACTER    :: AFMT*128
!
         DO I=0,NZEDS
            ZED(I) = REAL(I)/REALNZEDS
         ENDDO
!
         DO J=0,NALPS
            ALP(J) = REAL(J)/REALNALPS
         ENDDO
!
         DO ITABLE=1,NTABLES
            DO J=0,NALPS-1
               DO I=0,NZEDS-1
                  DO K=0,4
                     ALPHA = ALP(J)**2/(1.0-ALP(J)**2) * (-1.0)**ITABLE
                     IF (ITABLE.EQ.2.AND.J.EQ.0) ALPHA = 1.0E-10 ! ... must not be exactly zero
                     ZETASCALE = SQRT(1.0+ABS(ALPHA))
                     ABSZETA = ZED(I)/(1.0-ZED(I))
                     ZETA = ZETASCALE*ABSZETA
                     CALL ZTABLE_INTEGRAL(K,ZETA,ALPHA,ZFUNC)
                     IF (ITABLE.EQ.1) THEN
                        ZFUNC = ZFUNC*SQRT(ZETA**2+REAL(K+1)**2*ZETASCALE**2)**(K+1)
                     ENDIF
                     ZFUNC_REAL(K,I,J,ITABLE) = REAL(ZFUNC)
                     ZFUNC_IMAG(K,I,J,ITABLE) = AIMAG(ZFUNC)
                  ENDDO
               ENDDO
               IF (ITABLE.EQ.1) THEN
                  ZFUNC_REAL(:,NZEDS,J,ITABLE) = ASYMPT(:)
                  ZFUNC_IMAG(:,NZEDS,J,ITABLE) = 0.0
               ELSEIF (ITABLE.EQ.2) THEN
                  ZFUNC_REAL(:,NZEDS,J,ITABLE) = 0.0
                  ZFUNC_IMAG(:,NZEDS,J,ITABLE) = 0.0
               ENDIF
            ENDDO
            ZFUNC_REAL(:,:,NALPS,ITABLE) = 2.0*ZFUNC_REAL(:,:,NALPS-1,ITABLE) - &
                                                   ZFUNC_REAL(:,:,NALPS-2,ITABLE)
            ZFUNC_IMAG(:,:,NALPS,ITABLE) = 2.0*ZFUNC_IMAG(:,:,NALPS-1,ITABLE) - &
                                                   ZFUNC_IMAG(:,:,NALPS-2,ITABLE)
         ENDDO
!
         OPEN (UNIT=12,FILE=ATABLE)
!
         DO ITABLE=1,NTABLES
!
           AFMT = "(' == TABLE  ==',//)"
           WRITE (AFMT(12:12),'(I1)') ITABLE
           WRITE (12,AFMT)
!
            DO K=0,4
!
               AFMT = "(' == ZFUNC  ==')"
               WRITE (AFMT(12:12),'(I1)') K
               WRITE (12,AFMT)
!
               AFMT = "(/'REAL_PART')"
               WRITE (12,AFMT)
               AFMT = "('  ZED \ ALP  ',999(1PE13.5))"
               WRITE (12,AFMT) (ALP(J),J=0,MIN(NALPS,999))
               AFMT = "(1PE13.5,999(1PE13.5))"
               DO I=0,NZEDS
                  WRITE (12,AFMT) ZED(I),(ZFUNC_REAL(K,I,J,ITABLE),J=0,MIN(NALPS,999))
               ENDDO
!
               AFMT = "(//'IMAG_PART')"
               WRITE (12,AFMT)
               AFMT = "(' ZED \ ALPHA ',999(1PE13.5))"
               WRITE (12,AFMT) (ALP(J),J=0,MIN(NALPS,999))
               AFMT = "(1PE13.5,999(1PE13.5))"
               DO I=0,NZEDS
                  WRITE (12,AFMT) ZED(I),(ZFUNC_IMAG(K,I,J,ITABLE),J=0,MIN(NALPS,999))
               ENDDO
!
               AFMT = "(///)"
               WRITE (12,AFMT)
!
            ENDDO
!
         ENDDO
!
         CLOSE (UNIT=12)
!
         RETURN
!
      END SUBROUTINE ZTABLE_WRITE
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      SUBROUTINE ZTABLE_INTEGRAL(K,ZETA,ALPHA,ZFUNC)
!         USE KINDS
         IMPLICIT NONE
         INTEGER     , INTENT(IN)  :: K
         REAL        , INTENT(IN)  :: ZETA,ALPHA
         COMPLEX , INTENT(OUT) :: ZFUNC
         INTEGER  , PARAMETER :: NPERIOD = 100  ! ... points per zeta-cycle
         REAL, PARAMETER :: BETA = 4.5
         LOGICAL      :: L_AT_ALPHA_INV,L_ALPHA_POSITIVE
         REAL    :: TOL,TWOPI,ZETASCALE
         REAL    :: DT,HALFDT,T,DTZ
         REAL    :: AMPLITUDE_FACTOR,DERIV_FACTOR,TEST
         COMPLEX :: EYE,PHASE_FACTOR
!
         TWOPI = 8.0*ATAN(1.0)
         EYE = CMPLX(0.0,1.0)
         TOL = EXP(-BETA**2)
         L_ALPHA_POSITIVE = ALPHA.GT.0.0
         ZETASCALE = SQRT(1.0+ABS(ALPHA))
         DT = TWOPI/SQRT(ZETA**2+REAL(K+1)**2*ZETASCALE**2)/NPERIOD
         HALFDT = 0.5*DT
         ZFUNC = CMPLX(0.0,0.0)
         DTZ = HALFDT
         T   = 0.0
         L_AT_ALPHA_INV = .FALSE.
         AMPLITUDE_FACTOR = EXP(-0.25*T**2*(1.0-0.5*ALPHA*T)**2)
         IF (L_ALPHA_POSITIVE) AMPLITUDE_FACTOR = AMPLITUDE_FACTOR - &
                               EXP(-0.0625*(1.0-ALPHA*T)**2*(3.0-ALPHA*T)**2/ALPHA**2 &
                                   -0.0625/ALPHA**2)
         DERIV_FACTOR = T**K
         PHASE_FACTOR = (EYE**(K+1))*EXP(EYE*(ZETA*T))
         ZFUNC = ZFUNC + (DTZ*AMPLITUDE_FACTOR*DERIV_FACTOR)*PHASE_FACTOR
         TEST = AMPLITUDE_FACTOR*MAX(1.0,ABS(DERIV_FACTOR))
         DTZ = DT
         DO WHILE (TEST.GT.TOL.AND.(.NOT.L_AT_ALPHA_INV))
            T = T+DT
            IF (L_ALPHA_POSITIVE) L_AT_ALPHA_INV = (T+HALFDT).GE.(1.0/ALPHA)
            IF (L_AT_ALPHA_INV) DTZ = T+HALFDT-1.0/ALPHA
            IF (L_AT_ALPHA_INV) T   = 1.0/ALPHA-0.5*DTZ
            AMPLITUDE_FACTOR = EXP(-0.25*T**2*(1.0-0.5*ALPHA*T)**2)
            IF (L_ALPHA_POSITIVE) AMPLITUDE_FACTOR = AMPLITUDE_FACTOR - &
                                  EXP(-0.0625*(1.0-ALPHA*T)**2*(3.0-ALPHA*T)**2/ALPHA**2 &
                                      -0.0625/ALPHA**2)
            DERIV_FACTOR = T**K
            PHASE_FACTOR = (EYE**(K+1))*EXP(EYE*(ZETA*T))
            ZFUNC = ZFUNC + (DTZ*AMPLITUDE_FACTOR*DERIV_FACTOR)*PHASE_FACTOR
            TEST = AMPLITUDE_FACTOR*MAX(1.0,ABS(DERIV_FACTOR))
         ENDDO
         IF (.NOT.L_ALPHA_POSITIVE) RETURN
         DTZ = HALFDT
         T   = 0.0
         AMPLITUDE_FACTOR = EXP(-0.25*T**2*(1.0-0.5*ALPHA*T)**2)
         AMPLITUDE_FACTOR = AMPLITUDE_FACTOR - &
                            EXP(-0.0625*(1.0-ALPHA*T)**2*(3.0-ALPHA*T)**2/ALPHA**2 &
                                -0.0625/ALPHA**2)
         DERIV_FACTOR = T**K
         PHASE_FACTOR = (EYE**(K+1))*EXP(EYE*(ZETA*T))
         ZFUNC = ZFUNC + (DTZ*AMPLITUDE_FACTOR*DERIV_FACTOR)*PHASE_FACTOR
         TEST = AMPLITUDE_FACTOR*MAX(1.0,ABS(DERIV_FACTOR))
         DTZ = DT
         DO WHILE (TEST.GT.TOL)
            T = T-DT
            AMPLITUDE_FACTOR = EXP(-0.25*T**2*(1.0-0.5*ALPHA*T)**2)
            AMPLITUDE_FACTOR = AMPLITUDE_FACTOR - &
                               EXP(-0.0625*(1.0-ALPHA*T)**2*(3.0-ALPHA*T)**2/ALPHA**2 &
                                   -0.0625/ALPHA**2)
            DERIV_FACTOR = T**K
            PHASE_FACTOR = (EYE**(K+1))*EXP(EYE*(ZETA*T))
            ZFUNC = ZFUNC + (DTZ*AMPLITUDE_FACTOR*DERIV_FACTOR)*PHASE_FACTOR
            TEST = AMPLITUDE_FACTOR*MAX(1.0,ABS(DERIV_FACTOR))
         ENDDO
!
      END SUBROUTINE ZTABLE_INTEGRAL
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


