      SUBROUTINE COLIHIST (MODHIST, MAXHIST, JOBNUM, BCOLIP,
     >                     ITSTART, ITMAX, BSHORT, IVERS, 
     >                     BEMIX, BEPSGMAXERR, EMIXSTART, bVeloMod)
C***********************************************************************
C***  UPDATING THE MODEL HISTORY FOR PROGRAM COLI
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAXHIST, JOBNUM, IVERS,
     >                       ITSTART, ITMAX
      LOGICAL, INTENT(IN) :: BCOLIP, BSHORT, 
     >                       BEMIX, BEPSGMAXERR, 
     >                       bVeloMod
      REAL, INTENT(IN) :: EMIXSTART
      CHARACTER(MAXHIST*8), INTENT(INOUT) :: MODHIST

      INTEGER :: IDUMMY, IERR
      
      CHARACTER(8) :: BUFFER8
      CHARACTER(16) :: BUFFER16
      CHARACTER(24) :: BUFFER24
  
C***  ENTRY "COLI" OR "COLI+"
C***  - THE LATTER INDICATES NEW COMPUTATION OF EDDIEs
      IF (BCOLIP) THEN
        WRITE(UNIT=BUFFER16,FMT=19) JOBNUM, 'COLI+'
      ELSE
        WRITE(UNIT=BUFFER16,FMT=19) JOBNUM, 'COLI '
      ENDIF
   19 FORMAT ('/',I7,'. ',A5,1X)
      CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,16,BUFFER16)
 
C***  WRITE VERSION NUMBER
      WRITE(UNIT=BUFFER8,FMT=3) IVERS
    3 FORMAT (' (V.',I2,') ')
      CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,BUFFER8)


      IF (BSHORT) THEN
        WRITE(UNIT=BUFFER8,FMT='(A8)') '<BSHORT>'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,BUFFER8)
      ENDIF

      WRITE(UNIT=BUFFER24,FMT=28) ITSTART, ITMAX
   28 FORMAT (' ITSTART, ITMAX : ',I1,1X,I1)
      CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,24,BUFFER24)

      IF (.NOT. BEMIX) THEN
        WRITE(UNIT=BUFFER8,FMT='(A8)') ' NO EMIX'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,BUFFER8)
      ENDIF

      IF (BEPSGMAXERR) THEN
        WRITE(UNIT=BUFFER16,FMT=30) EMIXSTART
   30   FORMAT (' EMIX START',F5.2)
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,16,BUFFER16)
      ENDIF

      IF (bVeloMod) THEN
        BUFFER8 = ' VELOMOD'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,BUFFER8)
      ENDIF

C***  Save Model History
      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)

      RETURN
      END
