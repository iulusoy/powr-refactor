      SUBROUTINE DECOMO (LSOPA, LSINT, NOTEMP, LSHTOT, LPLOHTOT, 
     >            MODHIST, BUNLU, BPLOTRTAU1, 
     >            NPLOTOPA, OPTIONPLOTOPA, NPLOTOPADIM)
C*******************************************************************************
C***  DECODING INPUT OPTIONS FOR PROGRAM "COMO"
C*******************************************************************************
 
      CHARACTER(80) :: KARTE
      CHARACTER(80), DIMENSION(NPLOTOPADIM) :: OPTIONPLOTOPA
      PARAMETER (MAXPAR = 3)
      CHARACTER(20), DIMENSION(MAXPAR) :: ACTPAR
      LOGICAL NOTEMP, BUNLU, BPLOTRTAU1
      INTEGER :: NPLOTOPA, NPLOTOPADIM

C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      LSHTOT=-1
      LPLOHTOT=-1
      NCON=0
      NOTEMP=.FALSE.
      BUNLU = .FALSE.
      BPLOTRTAU1 = .FALSE.
      NPLOTOPA = 0

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
 
    8 READ (1,4, END=99) KARTE
    4 FORMAT (A)
      CALL SARGC(KARTE,NPAR)
      IF ( NPAR .LT. 1) GOTO 8
      IF ( NPAR .GT. 20) NPAR = 20
C***  Get actual parameters
      ACTPAR = ''
      DO I=1, NPAR
       CALL SARGV(KARTE,I,ACTPAR(I))
      ENDDO
 
C***  PRINT options
      IF (ACTPAR(1) .EQ. 'PRINT') THEN
C                         =====

         IF(ACTPAR(2) .EQ. 'INT') THEN
C                           ===
            READ (ACTPAR(3), 7) XL
    7       FORMAT (F20.0)
            LSINT=IFIX(XL)
            IF (LSINT.EQ.0) LSINT=1

         ELSEIF (ACTPAR(2) .EQ. 'OPA') THEN
C                                ===
            READ (ACTPAR(3), 7) XL
            LSOPA=IFIX(XL)
            IF (LSOPA.EQ.0) LSOPA=1

         ELSEIF (ACTPAR(2) .EQ. 'HTOTC') THEN
C                               ===========
            READ (ACTPAR(3), 7) XL
            LSHTOT=IFIX(XL)
            IF (LSHTOT.EQ.0) LSHTOT=1
         ENDIF

C***  PLOT options
      ELSEIF (ACTPAR(1) .EQ. 'PLOT') THEN
C                             ====

         IF (ACTPAR(2) .EQ. 'HTOTC') THEN
C                            =====
            LPLOHTOT=1

         ELSEIF (ACTPAR(2) .EQ. 'RTAU1') THEN
C                                =====
            BPLOTRTAU1 = .TRUE.

         ELSEIF (ACTPAR(2) .EQ. 'OPA') THEN
C                                ===
            NPLOTOPA = NPLOTOPA + 1
            IF (NPLOTOPA .LE. NPLOTOPADIM) THEN 
               OPTIONPLOTOPA(NPLOTOPA) = KARTE
            ELSE
               WRITE (*,*) 
     >            'WARNING: more PLOT OPA options than dimensioned' 
               WRITE (0,*) 
     >            'WARNING: more PLOT OPA options than dimensioned' 
            ENDIF
         ENDIF
C     Other options
      ELSE

         IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                            ========
            NOTEMP=.TRUE.
            IF (KARTE(30:40) .NE. ' ') 
     $         CALL DECNOT (MODHIST,MODHIST,KARTE,NOTEMP,'COMO')

         ELSEIF (ACTPAR(1) .EQ. 'UNLUTEC') THEN
C                                =======
            BUNLU=.TRUE.
         ENDIF
      ENDIF

      GOTO 8

C***  END-OF-FILE REACHED:
   99 CONTINUE
      CLOSE (1)
      RETURN

      END
