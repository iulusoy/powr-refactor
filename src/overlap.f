      SUBROUTINE OVERLAP (IBLENDS,MAXLAP,LASTIND, LASTINDAUTO,
     >           XLAMZERO,XLAMRED,
     >           XLAMBLUE,XMAX,EINST,NDIM,ELEVEL,INDNUP,INDLOW,NOLAP,
     >           VDOP,NBLENDS, KRUDAUT, EAUTO, NAUTO, EION)

C**********************************************************************
C***  CALLED FROM: STEAL
C***  FILLS THE MATRIX "IBLENDS": ROW 'IND' CONTAINS THE INDICES
C***     OF ALL LOCALLY OVERLAPPING LINES
C**********************************************************************

      DIMENSION IBLENDS(MAXLAP,LASTIND),INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION XLAMZERO(LASTIND), XLAMBLUE(LASTIND), XLAMRED(LASTIND)
      DIMENSION EINST(NDIM,NDIM), ELEVEL(NDIM), NBLENDS(LASTIND)
      DIMENSION KRUDAUT(NAUTO), EAUTO(NAUTO), EION(NDIM)
      LOGICAL NOLAP

C***  VELOCITY OF LIGHT (IN KM/SEC)
      DATA CLIGHT / 2.99792458E5 /

C***  INITIALIZATION AND WAVELENGHTS COMPUTATION
      DO 100 IND=1, LASTINDAUTO
C**      Distinguish normal lines from stabilizing lines (DRTRANSITs)
         IF (IND .LE. LASTIND) THEN
            WLAMZERO = ELEVEL(INDNUP(IND))-ELEVEL(INDLOW(IND))
         ELSE
            WLAMZERO = EION(INDLOW(IND)) - ELEVEL(INDLOW(IND))
     >           + EAUTO(IND-LASTIND)
         ENDIF
         XLAMZERO(IND) = 1.E8 / WLAMZERO
         XLAMBLUE(IND) = XLAMZERO(IND) * (1. - XMAX * VDOP / CLIGHT)
         XLAMRED(IND)  = XLAMZERO(IND) * (1. + XMAX * VDOP / CLIGHT)
         DO 150 LB=1, MAXLAP
  150    IBLENDS(LB,IND) = 0
  100 CONTINUE

C***  LOOP OVER ALL LINES (EXCEPT IRON-LINES)
      DO 200 IND =1, LASTINDAUTO
C***     RUDIMENTAL LINE?
            IF (INDTEST .LE. LASTIND) THEN 
               IF (EINST(LOW,NUP) .EQ. -2.) CYCLE
            ELSE
               IF (KRUDAUT(INDTEST-LASTIND) .EQ. 1) CYCLE
            ENDIF

C***  SPECIAL BRANCH IF OVERLAP OPTION IS NOT ACTIVE
      IF (NOLAP) THEN
         IBLENDS(1,IND) = IND
         NBLENDS(IND)=1
         ELSE

         LB = 0
         DO 300 INDTEST=1, LASTINDAUTO
            LOW = INDLOW(INDTEST)
            NUP = INDNUP(INDTEST)
C***        TEST LINE NOT RUDIMENTAL?
            IF (INDTEST .LE. LASTIND) THEN 
               IF (EINST(LOW,NUP) .EQ. -2.) CYCLE
            ELSE
               IF (KRUDAUT(INDTEST-LASTIND) .EQ. 1) CYCLE
            ENDIF
C***        TEST LINE OVERLAPPING?
            IF ((XLAMBLUE(IND) .LT. XLAMRED (INDTEST)  .AND.
     >            XLAMZERO(IND) .GE. XLAMZERO(INDTEST) )
     >           .OR.
     >           ( XLAMRED (IND) .GT. XLAMBLUE(INDTEST)   .AND.
     >            XLAMZERO(IND) .LE. XLAMZERO(INDTEST) )
     >           ) THEN
               LB = LB + 1
C***           DIMENSION CHECK
               IF (LB .GT. MAXLAP) THEN
                  WRITE (0,'(A)') '*** Dimension insufficient:'
                  WRITE (0,'(A,I5)') '*** MAXLAP =', MAXLAP
                  WRITE (0,'(A)') '*** FATAL ERROR in subr. OVERLAP'
                  STOP 'ERROR'
               ENDIF
               IBLENDS(LB,IND) = INDTEST
            ENDIF
  300    CONTINUE
         NBLENDS(IND)=LB
      ENDIF
  200 CONTINUE

      RETURN
      END
