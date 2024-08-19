      SUBROUTINE FECHECK (XLAMK, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEINDR)
C***********************************************************************
C***  CHECK FOR ACTIVE BOUND-BOUND TRANSITIONS OF GENERIC ION
C***  CALLED FROM COLI
C***********************************************************************

      DIMENSION INDRB(LASTFE), IFRBSTA(LASTFE), IFRBEND(LASTFE)
      DIMENSION INDFEACT(LASTFE)
      LOGICAL BFECHECK, BFEWING

C***  SEARCH FOR LOWEST POSSIBLE TRANSITION

      IND = INDFEACT(1)

C***  CALCULATION OF ACTUAL FREQUENCY INDEX
      XINDF = - ALOG10(XLAM0FE/XLAMK) /
     >          ALOG10(1. + VDOPFE*DXFE/CLIGHT)

 10   IF (IND .GT. LASTFE) THEN
C ***   NO TRANSITIONS ARE ACTIVE ANYMORE
         BFECHECK = .FALSE.
         GOTO 20
      ENDIF

      XINDLOW  = FLOAT(IFRBSTA(IND))

      IF (XINDF .LT. XINDLOW) THEN
C ***   NO ACTIVE TRANSITION
         BFECHECK = .FALSE.
         GOTO 20
      ELSE
C ***   TRASITION ACTIVE?
         XINDUP = FLOAT(IFRBEND(IND))
         IF (XINDF .GT. XINDUP) THEN
C ***      TEST NEXT TRANSITION
            IND = IND + 1
            INDFEACT(1) = IND
            GOTO 10
         ELSE
C ***      TRANSITION ACTIVE!
            BFECHECK = .TRUE.
            GOTO 20
         ENDIF
      ENDIF

 20   CONTINUE

C *** TEST IF RED WING OF ANY LINE IS ACTIVE
      IF ((.NOT. BFECHECK) .AND. (IND .GT. 1)) THEN
C ***   RED WING OF LAST LINE, WHICH HAS BEEN ACTIVE
         XINDUP = FLOAT(IFRBEND(IND-1))
         XINDWING = XINDUP + DFEINDR
         IF (XINDF .LT. XINDWING) THEN
C ***      WING IS ACTIVE!
            BFEWING = .TRUE.
         ELSE
            BFEWING = .FALSE.
         ENDIF
      ELSE
C ***   WING IS ACTIVE, WHEN LINE IS ACTIVE!         
         BFEWING = BFECHECK
      ENDIF

C *** STORE INDICES OF ACTIVE LINES IN ARRAY INDFEACT
      MAXFEACT = 0
      DO IND=INDFEACT(1), LASTFE
         XINDLOW = FLOAT(IFRBSTA(IND))
         IF (XINDF .GE. XINDLOW) THEN
            XINDUP = FLOAT(IFRBEND(IND))
            IF (XINDUP .GE. XINDF) THEN
C ***         TRANSITION ACTIVE !
               MAXFEACT = MAXFEACT + 1
               INDFEACT(MAXFEACT) = IND
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END
