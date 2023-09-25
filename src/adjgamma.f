      SUBROUTINE ADJGAMMA(GAHIST, MAXGAHIST, AG, LASTHYDRO, LASTTAU,
     >                    BTALTER, GAMMAC, GAMMAL, GAMMAR, GAMMAD, 
     >                    BGFIN, GF, BGAMMACFIX, BGAMMALFIX, 
     >                    BGAMMARFIX, BGAMMADFIX, bHDNoAG)
C*************************************************************
C***  Automatic Gamma Adjustment
C***    GAMMA HISTORY has been already shifted!
C*************************************************************

      DIMENSION GAHIST(26,MAXGAHIST), AG(7)
      LOGICAL :: BKONVER1, BKONVER2, BGFIN, BTALTER, bPrintZeroInfo
      LOGICAL :: BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX
      LOGICAL :: bHDNoAG   ! if true, last HD corrections are so large that GAMMAs should be reset to inf (0)
      INTEGER :: LASTHYDRO, LASTTAU

      bPrintZeroInfo = .FALSE.

C*** Explicitely set gamma's win over autogamma
      IF (BGAMMACFIX) GAMMACFIX = GAMMAC 
      IF (BGAMMALFIX) GAMMALFIX = GAMMAL 
      IF (BGAMMARFIX) GAMMARFIX = GAMMAR 
      IF (BGAMMADFIX) GAMMADFIX = GAMMAD 

C***  Number of last STEAL
      IF (GAHIST(2,2) .NE. 0.) THEN
        LST = 3
      ELSE
        LST = 2
      ENDIF

      IF (GAHIST(19,LST) .GT. 0.) THEN
        CORRMAX = ALOG10(GAHIST(19,LST))
      ELSE
        CORRMAX = 999.
      ENDIF

C***  Depth Points not konverged in the last iterations?
C***  BKONVER1 : Last Iteration considered
C***  BKONVER2 : Last 6 Iterations considered
      BKONVER1 = .TRUE.
      BKONVER2 = .TRUE.
      DO I=2, 7
        IF (GAHIST(1,I) .LT. 0. .OR. GAHIST(2,I) .NE. 0.) EXIT
        IF (I .EQ. 2 .AND. GAHIST(25,I) .GT. 0.) BKONVER1 = .FALSE.
        IF (GAHIST(25,I) .GT. 0.) BKONVER2 = .FALSE.
      ENDDO
      IF ((LASTHYDRO > -1) .AND. (LASTHYDRO < 6)) THEN
        BKONVER2 = .FALSE.
        IF (LASTHYDRO < 2) BKONVER1 = .FALSE.
      ENDIF

C***  Last change of Gammas
      LASTCH = 8
      DO I=3, MAXGAHIST
        IF (GAHIST(1,I) .LT. 0. .OR. GAHIST(2,I) .NE. 0.) EXIT
        IF (GAHIST(3,LST) .NE. GAHIST(3,I) .OR. 
     >      GAHIST(4,LST) .NE. GAHIST(4,I) .OR. 
     >      GAHIST(5,LST) .NE. GAHIST(5,I) .OR. 
     >      GAHIST(6,LST) .NE. GAHIST(6,I)) THEN
          LASTCH = I - 2
          EXIT
        ENDIF
      ENDDO

      GAMMAC = GAHIST(3,LST)
      GAMMAL = GAHIST(4,LST)
      GAMMAR = GAHIST(5,LST)
      GAMMAD = GAHIST(6,LST)

      WRITE (0,*) 
      WRITE (0,*) 'Automatic Gamma Adjustment is active'
      WRITE (0,*) '------------------------------------'
      WRITE (0,'(4X,4A11)') 'GAMMAC', 'GAMMAL', 'GAMMAR', 'GAMMAD'
      IF (GAMMAC .GT. 1.E7) THEN
        WRITE (0,'(A4,4(1X,E10.4))')
     >   'Old:', GAMMAC, GAMMAL, GAMMAR, GAMMAD
      ELSE
        WRITE (0,'(A4,4(1X,F10.2))')
     >   'Old:', GAMMAC, GAMMAL, GAMMAR, GAMMAD
      ENDIF

      IF ((GAHIST(1,2) < 0.) .OR. (LASTHYDRO == 2)) THEN
C***  No GAMMA HISTORY AVAILABLE (reset after hydro iteration)
        GAMMAC = AG(1)
        GAMMAL = GAMMAC / AG(5)
        GAMMAR = GAMMAC / AG(6)
        GAMMAD = GAMMAC / AG(7)
C***  Increase Gamma
      ELSEIF (.NOT. BKONVER1 .OR. CORRMAX .GT. AG(4)) THEN
        GAMMAC = GAHIST(3,LST) * 2.
        GAMMAL = GAMMAC / AG(5)
        GAMMAR = GAMMAC / AG(6)
        GAMMAD = GAMMAC / AG(7)
C***  Decrease Gamma
      ELSEIF (BKONVER2 .AND. ((CORRMAX .LE. AG(3) .AND.
     >        (LASTCH .GE. 6 .OR. GAHIST(3,LST) .GE. 100.)) .OR. 
     >        CORRMAX .LE. -2. .OR. 
     >        (CORRMAX .LE. AG(3)+1. .AND. GAHIST(3,LST) .GE. 100.))) 
     >        THEN
        IF (GAHIST(3,LST) .EQ. 0.) THEN
          GAMMAC = 640.
          GAMMAL = GAMMAC / AG(5)
          GAMMAR = GAMMAC / AG(6)
          GAMMAD = GAMMAC / AG(7)
        ELSE
          GAMMAC = GAHIST(3,LST) / 2.
          GAMMAL = GAMMAC / AG(5)
          GAMMAR = GAMMAC / AG(6)
          GAMMAD = GAMMAC / AG(7)
        ENDIF
      ENDIF

C***  Do not decrease below final Gamma
      IF (GAMMAC .LT. AG(2) .AND. GAMMAC .GT. 0.) THEN
        GAMMAC = AG(2)
        GAMMAL = GAMMAC / AG(5)
        IF (GAMMAL .LT. 1.) THEN GAMMAL = 1
        GAMMAR = GAMMAC / AG(6)
        IF (GAMMAR .LT. 1.) THEN GAMMAR = 1
        GAMMAD = GAMMAC / AG(7)
        IF (GAMMAD .LT. 1.) THEN GAMMAD = 1
      ENDIF

C***  Keep within upper-limit Gamma (1st parameter on CARD, if >0)
      IF (GAMMAC .GT. AG(1) .AND. AG(1) .GT. 0.) THEN
        GAMMAC = AG(1)
        GAMMAL = GAMMAC / AG(5)
        GAMMAR = GAMMAC / AG(6)
        GAMMAD = GAMMAC / AG(7)
      ENDIF

C***  Set BGFIN = .TRUE., if final gammas are reached
      BGFIN = GAMMAC .EQ. AG(2)

      IF (GF .GT. 0.) THEN
        write (0,*) 'GAMMAS forced to GF=', gf
        GAMMAC = GF
        GAMMAL = GF
        GAMMAR = GF
        GAMMAD = GF
      ENDIF

C***  Prevent very large GAMMA values as they could crash the code
C***  (and ruin the MODHIST output)
      IF (GAMMAC > 9999.) THEN
        GAMMAC = 0.
        GAMMAL = 0.
        GAMMAR = 0.
        GAMMAD = 0.
      ENDIF      

      IF (BGAMMACFIX) GAMMAC = GAMMACFIX 
      IF (BGAMMALFIX) GAMMAL = GAMMALFIX
      IF (BGAMMARFIX) GAMMAR = GAMMARFIX
      IF (BGAMMADFIX) GAMMAD = GAMMADFIX

C***  If this is the first iteration with temperature equation/correction, 
C***  start with GAMMA=0
      IF ((LASTHYDRO == 1) .OR. (bHDNoAG .AND. LASTHYDRO == 2) .OR.
     >    ((.NOT. BTALTER .AND. LASTHYDRO /= 2) .AND. 
     >            (GAHIST(26,1) > .0 .AND. GAHIST(26,2) == .0)) .OR. 
C***  Go back two STEAL jobs if TCORR ALTERNATE or HD update happened two jobs ago
     >    ((BTALTER .OR. LASTHYDRO == 2) .AND.           
     >            (GAHIST(26,1) > .0 .AND. GAHIST(26,3) == .0)) ) THEN

C***      reset to the upper-limit Gamma (first AUTO GAMMA parameter) 
C***      Note: this replaces the previous version, where all GAMMAs were
C***      reset to zero -- wrh  2-Apr-2020
 
          GAMMAC = AG(1)
          GAMMAL = GAMMAC / AG(5)
          GAMMAR = GAMMAC / AG(6)
          GAMMAD = GAMMAC / AG(7)

          bPrintZeroInfo = .TRUE.
      ENDIF


      IF (GAMMAC .GT. 1.E7) THEN
        WRITE (0,'(A4,4(1X,E10.4))')
     >   'New:', GAMMAC, GAMMAL, GAMMAR, GAMMAD
      ELSE
        WRITE (0,'(A4,4(1X,F10.2))')
     >   'New:', GAMMAC, GAMMAL, GAMMAR, GAMMAD
      ENDIF
      IF (bPrintZeroInfo) THEN
        IF (bHDNoAG) THEN
          WRITE (0,'(A)') 
     >      "STEAL: GAMMA reset due to large HD stratification updates"
        ELSEIF (LASTHYDRO == 1) THEN
          WRITE (0,'(A)') 
     >      "STEAL: First STEAL after HD update -> all GAMMAs = 0"
        ELSE
          WRITE (0,'(A)') 
     >      'STEAL: GAMMA reset due to onset of temperature corrections' 
        ENDIF
      ENDIF

      RETURN
      END
