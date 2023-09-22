      SUBROUTINE DECON (LSOPA,LSINT,IFLUX,JOBMAX,MODHIST, 
     >                  BUNLU, bCLXJC, IVERS, POPMIN)
C***********************************************************************
C***  DECODING INPUT OPTIONS, CALLED FROM WRCONT *******************************
C***********************************************************************

      IMPLICIT NONE

      INTEGER :: I, NPAR, IFLUX, JOBMAX, IVERS,
     >           LSOPA, LSINT
      INTEGER, DIMENSION(1) :: MODHIST

      REAL :: XL, POPMIN
 
      CHARACTER(14), DIMENSION(5) :: ACTPAR
      CHARACTER(80) :: KARTE

      INTEGER, EXTERNAL :: IDX

      LOGICAL :: BUNLU, bCLXJC, bReadPOPMIN

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0 !write to wruniqX.cpr (stderr)

 
C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      IFLUX=-1
      JOBMAX=-1
      BUNLU = .FALSE.
      bCLXJC = .FALSE.
      IVERS = 4
      bReadPOPMIN = (POPMIN < 1.E-99)
      IF (bReadPOPMIN) THEN
        POPMIN = 1.E-25           !default value should be the same as in STEAL -> DECSTE
      ENDIF
 
      OPEN (1,FILE='CARDS', STATUS='UNKNOWN')
      REWIND (1)

      DO !--- Loop over all CARDS lines ---

        READ (1, '(A)', END=1) KARTE

        CALL SARGC(KARTE,NPAR)
        IF ( NPAR .LT. 1) CYCLE
        IF ( NPAR .GT. 5) NPAR = 5

        DO I=1, NPAR
          CALL SARGV(KARTE,I,ACTPAR(I))
        ENDDO


        IF (KARTE(:10) .EQ. 'PRINT FLUX') THEN
C                            ==========
          IFLUX=1
          CYCLE
        ENDIF
        IF (KARTE(:10) .EQ. 'PRINT INT ') THEN
C                            ==========
          !DECODE (80,7,KARTE) XL
          READ (UNIT=KARTE, FMT=7) XL
    7     FORMAT (10X,F10.0)
          LSINT=IFIX(XL)
          IF (LSINT.EQ.0) LSINT=1
          CYCLE
        ENDIF
        IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN
C                            ==========
          !DECODE (80,7,KARTE) XL
          READ (UNIT=KARTE, FMT=7) XL
          LSOPA=IFIX(XL)
          IF (LSOPA.EQ.0) LSOPA=1
          CYCLE
        ENDIF
        IF (KARTE(:7) .EQ. 'JOBMAX=') THEN
C                           =======
          !DECODE (80,3,KARTE) XL
          READ (UNIT=KARTE, FMT=3) XL
    3     FORMAT (7X,F10.9)
          JOBMAX=IFIX(XL)
          CYCLE
        ENDIF
c        IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
cC                           ========
c              NOTEMP=.TRUE.
c              IF (KARTE(30:40) .NE. ' ')
c     $          CALL DECNOT (MODHIST(1),MODHIST,KARTE,NOTEMP,'WRCONT')
c              CYCLE
c              ENDIF
        IF (KARTE(:7) .EQ. 'UNLUTEC') THEN
C                           =======
          BUNLU=.TRUE.
          CYCLE
        ENDIF

        IF (ACTPAR(1) .EQ. 'OB-VERS') THEN
C                           =======
          IF (NPAR .EQ. 2 .OR. ACTPAR(3) .EQ. 'WRCONT') THEN
            READ (ACTPAR(2),'(I10.0)', ERR=90) IVERS
C            write (0,*) 'DECCON: IVERS = ', ivers
          ENDIF
          CYCLE 
        ENDIF
 
        IF (ACTPAR(1)(1:6) == 'POPMIN' .AND. bReadPOPMIN) THEN
C                              ======
          READ (ACTPAR(2),'(F10.0)', ERR=90) POPMIN
          CYCLE
        ENDIF         

        IF (KARTE(:8) == 'COLI-XJC') THEN
C                         ========
          bCLXJC = .TRUE.
          CYCLE
        ENDIF
        
      ENDDO
 
C***  END-OF-FILE REACHED:
   1  CONTINUE
      CLOSE (1)
      RETURN
 
C***  ERROR EXITS **************************************************

   90 WRITE (hCPR,*) '*** DECON: FATAL ERROR WHEN DECODING NUMBER'
      WRITE (hCPR,*) '*** THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (hCPR,*) KARTE(:IDX(KARTE))
      STOP 'ERROR in WRCONT'
      
      END
