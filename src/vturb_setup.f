      SUBROUTINE VTURB_SETUP (VTURB, VTURB_LINE, ND,
     >                        VELO, RADIUS, TAUROSS, T, XMU)
C***********************************************************************
C***  This routine fills the VTURB vector as specified in the CARDS file
C***
C***  This routine is work in progress
C***
C***  called from WRSTART (so far only, planned: ENSURETAUMAX, HYDROSOLVE)
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, INTENT(IN) :: ND
      
      REAL, DIMENSION(ND), INTENT(IN) :: VELO, RADIUS, TAUROSS, T, XMU
      REAL, DIMENSION(ND), INTENT(INOUT) :: VTURB
      
      INTEGER, PARAMETER :: NMAXPAR = 40
      
      CHARACTER*(*) :: VTURB_LINE
      CHARACTER(40), DIMENSION(NMAXPAR) :: ACTPAR
      
      INTEGER, EXTERNAL :: IDX
      
      REAL :: VTURBND, VMICND, VNORM, VMIC_FRAC, VMIC_MAX, DX, X, Q,
     >        DPARVMIC_IN, DPARVMIC_OUT, DUMMY, 
     >        V_INCREASE, TAU_INCREASE, R_INCREASE
      INTEGER :: I, NPAR, L

      LOGICAL, SAVE :: bWARN
      
C***  Physical constants
      REAL, PARAMETER :: PI = 3.141592654 

C***  Local parameters and conversion factors      
      REAL, PARAMETER ::  VMICFRAC_DEFAULT = 0.05
      REAL, PARAMETER ::  FVTURBTOMIC      = SQRT(2.)
      REAL, PARAMETER ::  FVMICTOTURB      = 1./FVTURBTOMIC
            
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      DATA bWARN /.FALSE./  
    

C***  Read number of parameters and put them into ACTPAR array      
      CALL SARGC(VTURB_LINE, NPAR)
      IF (NPAR <= 1) THEN
C***    No VTURB/VMIC line in CARDS file
        DO L=1, ND
          VTURB(L) = 0.
        ENDDO
        IF (.NOT. bWARN) THEN
           WRITE (hCPR,'(A)') 
     >        'VTURB_SETUP> No VTURB or VMIC card has been found!'
           bWARN = .TRUE.
        ENDIF
        RETURN
      ENDIF
      
      
      DO I=1, MIN(NPAR, NMAXPAR)
        CALL SARGV(VTURB_LINE,I,ACTPAR(I))
      ENDDO
      
C***  First parameter: Keyword VTURB or VMIC
C       for VMIC, values are divided by sqrt(2) 
C       to be in line with definition of VMIC in FORMAL      
      IF (ACTPAR(1) == 'VMIC') THEN
        VNORM = 1.
      ELSE
        VNORM = FVTURBTOMIC
      ENDIF
      
C***  Second parameter: standard value for WRSTART initializing     
      READ (ACTPAR(2),'(F20.0)', ERR=99) VTURBND
      VMICND = VNORM * VTURBND
      
      IF (NPAR == 2) THEN
C***    Only two parameters: constant VMIC/VTURB
        DO L=1, ND
          VTURB(L) = FVMICTOTURB * VMICND
        ENDDO      
      ELSE
C***    More than 2 parameters: VMIC/VTURB not constant 
        IF (ACTPAR(3) == 'VELOFRAC') THEN
C***      VELOFRAC option: VMIC is scaled with wind velocity        
          IF (NPAR > 3) THEN
            READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_FRAC
          ELSE 
            VMIC_FRAC = VMICFRAC_DEFAULT
          ENDIF
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * MAX(VMICND, VMIC_FRAC * VELO(L))
          ENDDO
        ELSEIF (ACTPAR(3) == 'HASER') THEN
C***      Similar to MAX option, but with slightly different increasing
C***       motivated by Haser et al. (1998)
          IF (NPAR == 3) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: HASER option needs specified MAX value!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          VMIC_MAX = VNORM * VMIC_MAX
          VMIC_FRAC = (VMIC_MAX - VMICND) / VELO(1)
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * (VMICND + VMIC_FRAC * VELO(L))
          ENDDO
        ELSEIF (ACTPAR(3) == 'HILLIER') THEN
C***      Approach similar to the description used in Hillier's depth-dependent clumping
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: HILLIER option needs two more '//
     >        'parameters, VMIC_MAX and V_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) V_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >                  (VMIC_MAX - VMICND) * EXP(-VELO(L)/V_INCREASE))
          ENDDO
        ELSEIF (ACTPAR(3) == 'EXPTAU') THEN
C***      Exponential VTURB onset (like in Hillier's formula)
C***      but using the Tauross_cont scale
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: EXPTAU option needs two more '//
     >        'parameters, VMIC_MAX and TAU_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) TAU_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            IF (TAUROSS(L) <= 1.E-30) THEN
              VTURB(L) = FVMICTOTURB * VMIC_MAX 
            ELSE
              VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >             (VMIC_MAX - VMICND) * EXP(-TAU_INCREASE/TAUROSS(L)))
            ENDIF
          ENDDO
        ELSEIF (ACTPAR(3) == 'EXPRADIUS') THEN
C***      Exponential VTURB onset (like in Hillier's formula)
C***      but using the Tauross_cont scale
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: EXPRADIUS option needs two more '//
     >        'parameters, VMIC_MAX and R_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) R_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            IF (TAUROSS(L) <= 1.E-30) THEN
              VTURB(L) = FVMICTOTURB * VMIC_MAX 
            ELSE
              VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >         (VMIC_MAX-VMICND) * EXP(-(RADIUS(L)-1.)/(R_INCREASE-1.)))
            ENDIF
          ENDDO
        ELSEIF (ACTPAR(3) == 'MAX') THEN
C***      MAX option: Outer maximum of VMIC/VTURB is given
          IF (NPAR == 3) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: MAX option must have a specified value!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          VMIC_MAX = VNORM * VMIC_MAX
          IF (NPAR == 4) THEN
C***        No further parameters: increase with wind velocity (fraction)
            VMIC_FRAC = VMIC_MAX / VELO(1)
            DO L=1, ND
              VTURB(L) = FVMICTOTURB * MAX(VMICND, VMIC_FRAC * VELO(L))
            ENDDO
          ELSEIF (NPAR == 8) THEN                        
            READ(ACTPAR(6), '(F20.0)', ERR=99) DPARVMIC_IN
            READ(ACTPAR(8), '(F20.0)', ERR=99) DPARVMIC_OUT
            IF (ACTPAR(5) == 'VELO1' .AND. ACTPAR(7) == 'VELO2') THEN
C***          Interpolate between min and max value in specified velocity regime            
              IF (DPARVMIC_IN > DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN >= VELO(1)) .OR. 
     >            (DPARVMIC_OUT <= VELO(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner wind velocity: ', VELO(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer wind velocity: ', VELO(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MAX(DPARVMIC_IN, VELO(ND))
              DPARVMIC_OUT = MIN(DPARVMIC_OUT, VELO(1))
              DX = DPARVMIC_OUT - DPARVMIC_IN
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (VELO(L) <= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (VELO(L) >= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( VELO(L) - DPARVMIC_IN ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
              
            ELSEIF (ACTPAR(5) == 'TAU1' .AND. ACTPAR(7) == 'TAU2') THEN
C***          Interpolate between min and max value in specified TAU regime            
              IF (DPARVMIC_IN < DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN <= TAUROSS(1)) .OR. 
     >            (DPARVMIC_OUT >= TAUROSS(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner Tau: ', TAUROSS(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer Tau: ', TAUROSS(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MIN(DPARVMIC_IN, TAUROSS(ND))
              DPARVMIC_OUT = MAX(DPARVMIC_OUT, TAUROSS(1))
              DX = DPARVMIC_IN - DPARVMIC_OUT
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (TAUROSS(L) >= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (TAUROSS(L) <= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( DPARVMIC_IN - TAUROSS(L) ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
            
            ELSEIF (ACTPAR(5) == 'R1' .AND. ACTPAR(7) == 'R2') THEN
C***          Interpolate between min and max value in specified radius regime            
              IF (DPARVMIC_IN > DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN >= RADIUS(1)) .OR. 
     >            (DPARVMIC_OUT <= RADIUS(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner radius: ', RADIUS(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer radius: ', RADIUS(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MAX(DPARVMIC_IN, RADIUS(ND))
              DPARVMIC_OUT = MIN(DPARVMIC_OUT, RADIUS(1))
              DX = DPARVMIC_OUT - DPARVMIC_IN
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (RADIUS(L) <= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (RADIUS(L) >= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( RADIUS(L) - DPARVMIC_IN ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
            
            ELSE 
              WRITE (hCPR,'(A)') '*** FATAL ERROR: Inconsistent ' //
     >          ' detail options for VTURB/VMIC!'
              GOTO 99
            ENDIF
          ELSE
C***        Format does not seem to be valid          
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: Unknown detail options for VTURB/VMIC!'
            GOTO 99
          ENDIF
          
        ENDIF
      ENDIF
      
      RETURN
      
C***  Error handling

   99 WRITE (hCPR,'(A)') 
     >  '*** FATAL ERROR: Decoding error in VTURB/VMIC line:', 
     >  VTURB_LINE
      STOP 'ERROR IN VTURB_SETUP'
      
      END
