      SUBROUTINE CMFFEOP (XLAMK, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDEXMAX, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT,
     >                    SIGMAACT, OPAFEI, ETAFEI,
     >                    OPAFEION, ETAFEION, T, IVERS_FE_EXPFAC, 
     >                    TEFF, NCHARG, MAXION, bFELASER, bNoIronLaser)

C***********************************************************************
C***  NON-LTE IRON OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  CALLED FROM: COLI and FORMAL -> FORMCMF
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, N, MAXION, LASTFE, INDEXMAX,
     >                       IVERS_FE_EXPFAC, MAXFEACT
      REAL, INTENT(IN) :: TEFF, RSTAR, CLIGHT, XLAMK, XLAM0FE, DXFE, 
     >                    VDOPFE
      LOGICAL, INTENT(IN) :: bNoIronLaser

      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      INTEGER, DIMENSION(LASTFE) :: IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                              INDFEACT, INDRB
      INTEGER, DIMENSION(N) :: IFRBSSTA, IFRBSEND, INDFESACT,
     >                         INDRBS, IFES
      REAL, DIMENSION(ND) :: OPAFE, ETAFE, ENTOT, T
      REAL, DIMENSION(ND, MAXION) :: OPAFEION, ETAFEION
      REAL, DIMENSION(ND,N) :: POPNUM, POPLTE
      REAL, DIMENSION(LASTFE) :: SIGMAACT
      REAL, DIMENSION(ND,LASTFE) :: OPAFEI, ETAFEI
      REAL, DIMENSION(N) :: ELEVEL, WEIGHT
      LOGICAL, DIMENSION(ND) :: bFELASER
      INTEGER, DIMENSION(N) :: NCHARG

      REAL :: SIGMA, SIGMAR, PLASER, ENUP, ENLOW, SUM, G, EMINDU,
     >        DELTA, SIGMA1, SIGMA2, WAVNUM0, WAVNUM03, XINDF,
     >        WEIGHTLU, WAVNUMK, WAVNUMK3
      INTEGER :: L, ION, INDF, INDFSTART, NUP, LOW, IND, INDACT
            
     
C***  Numerical parameters     
      REAL, PARAMETER :: EXPMAX = 500.
     
C***  Physical constants
      REAL, PARAMETER :: C1 = 1.4388        !C1 = H * C / K    ( CM * ANGSTROEM )
      REAL, PARAMETER :: C2 = 3.9724E-16    !C2 = 2 * H * C     ( CGS UNITS )

C***  CALCULATE FE-FREQUECY-INDEX
      XINDF = - LOG10(XLAM0FE/XLAMK) /
     >          LOG10(1. + VDOPFE*DXFE/CLIGHT)

C***  WAVENUMBERS IN KAYSER ( = CM**-1)
      WAVNUMK = 1.E8 / XLAMK
      WAVNUMK3 = WAVNUMK * WAVNUMK * WAVNUMK 

C***  PRESET OPAFE, ETAFE
      DO L=1, ND
         OPAFE(L) = 0.
         ETAFE(L) = 0.
         DO ION=1, MAXION
           OPAFEION(L, ION) = 0.
           ETAFEION(L, ION) = 0.
         ENDDO
      ENDDO
      
C***  LOOP OVER ACTIVE FE-LINES
      DO 10 INDACT=1, MAXFEACT

         IND = INDFEACT(INDACT)
         LOW = IFELOW(IND)
         NUP = IFENUP(IND)
         ION = NCHARG(LOW) + 1
         WEIGHTLU = WEIGHT(LOW) / WEIGHT(NUP)
         WAVNUM0 = ELEVEL(NUP) - ELEVEL(LOW)
         WAVNUM03 = WAVNUM0 * WAVNUM0 * WAVNUM0 

C***    FIND INDICES IN ARRAY 'SIGMAFE'
         INDF = NINT(XINDF-0.5)-IFRBSTA(IND)
         INDFSTART = INDRB(IND)
C***    READ NEIGHBOURING VALUES
         SIGMA1 = SIGMAFE(INDFSTART + INDF)
         SIGMA2 = SIGMAFE(INDFSTART + INDF + 1)
C***    INTERPOLATION
         DELTA = XINDF - FLOAT(INDF+IFRBSTA(IND))
         SIGMA = (1.-DELTA)*SIGMA1 + DELTA*SIGMA2

C***    MULTIPLY WITH RSTAR (DIMENIONS OF OPA AND ETA)
         SIGMAR = SIGMA * RSTAR

ccc     For test use: if "nu-factors" are activated:
ccc     modify SIGMA --> used in FREQUINT for integrartion of XJFEMEAN
ccc     -- the following statement might be de-activated else
ccc         SIGMA = SIGMA * WAVNUM0 / WAVNUMK
ccc    Note!!! SIGMAACT is NOT MANIPULATED in the Version with correct RRATE 
ccc            integration!!

C***    STORE SIGMA FOR LATER USE
         SIGMAACT(INDACT) = SIGMA

         IF (SIGMA <= 0.) CYCLE
         
C***    CALCULATE 'OPA' AND 'ETA' OVER DEPTH INDEX 'L' 
         DO 20 L=1, ND
            ENUP  = POPNUM(L,NUP) * ENTOT(L)
            ENLOW = POPNUM(L,LOW) * ENTOT(L)
C***       Opacity and Emissivity
C***     Different versions according to the IRONLINES-EXPFAC option:
C***     0: OFF  - no exp factor (recommended by Goetz) 
C***     1: TEFF - exp factor with max(Te,Teff) (recommended by wrh, default) 
C***     2: TRADITIONAL - exp factor with Te (standard till 17-Feb-2018)
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              G = WEIGHTLU 
            ELSEIF (IVERS_FE_EXPFAC .EQ. 1) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/MAX(T(L),TEFF))
            ELSEIF (IVERS_FE_EXPFAC .EQ. 2) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/T(L))
            ELSE 
              STOP '*** FATAL INTERNAL ERROR 1 in subr. CMFFEOP'
            ENDIF
            EMINDU = G * ENUP * SIGMAR
            SUM = ENLOW * SIGMAR - EMINDU
            OPAFEI(L, IND) = SUM
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUM03
            ELSE
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUMK3
            ENDIF
ccc         in the preceding stament WAVNUM03 has been used by Goetz

            PLASER = ENUP/WEIGHT(NUP) - ENLOW/WEIGHT(LOW)
            IF (bNoIronLaser .AND. PLASER > 0.) THEN
              !skip level
            ELSEIF (OPAFEI(L,IND) .GT. 0.) THEN
               OPAFE(L) = OPAFE(L) + OPAFEI(L,IND)
               ETAFE(L) = ETAFE(L) + ETAFEI(L,IND)
               OPAFEION(L,ION) = OPAFEION(L,ION) + OPAFEI(L,IND) 
               ETAFEION(L,ION) = ETAFEION(L,ION) + ETAFEI(L,IND) 
                              
              IF (PLASER > 0. .AND. (.NOT. bFELASER(L))) THEN
c                WRITE (hCPR,*) 'FE POP INVERSION AT L = ', L
                bFELASER(L) = .TRUE.
              ENDIF
            
            ENDIF
            
 20      CONTINUE

 10   CONTINUE

      RETURN
      END
