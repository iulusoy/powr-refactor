      SUBROUTINE DEC2BETA(KARTE, BETA2, BETA2FRACTION, RONSET)
C***********************************************************************      
C***  Decoding of 2BETALAW card
C***
C***  called by DECSTAR, DECSTE
C***********************************************************************      

      IMPLICIT NONE

      CHARACTER(40) :: TRYPAR
      CHARACTER(40), DIMENSION(20) :: CURPAR
      CHARACTER(100), INTENT(IN) :: KARTE
      
      REAL, INTENT(INOUT) :: BETA2, BETA2FRACTION, RONSET
      
      INTEGER :: NPAR, i, IERR
      
      LOGICAL :: bOldDecode 
      LOGICAL, DIMENSION(2) :: bParamFound

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      
      bOldDecode = .FALSE.
      CALL SARGC (KARTE, NPAR)
      IF (NPAR < 5) THEN
        WRITE (hCPR,'(A)') '*** VELPAR: NOT ENOUGH PARAMETERS'
        STOP '*** FATAL ERROR WHILE DECODING VELPAR CARDS-LINE'
      ENDIF

C***  Default: negative value for R_onset
C***  => means that second beta law starts at R_con
      RONSET = -1.
      
      bParamFound = .FALSE.
      DO i=1, NPAR
        CALL SARGV(KARTE,i,CURPAR(i))
      ENDDO
      IF (NPAR > 2) THEN
        DO i=2, NPAR 
          SELECTCASE (CURPAR(i))
            CASE ('BETA2')
              IF (NPAR >= (i+1)) THEN
                TRYPAR = CURPAR(i+1)                      
                READ (TRYPAR, '(F10.0)', IOSTAT=IERR, ERR=92) BETA2
                IF (IERR == 0) THEN
                  bParamFound(1) = .TRUE.
                ENDIF             
              ENDIF
            CASE ('FRACTION')
              IF (NPAR >= (i+1)) THEN
                TRYPAR = CURPAR(i+1)                      
                READ (TRYPAR, '(F10.0)', IOSTAT=IERR,
     >                                          ERR=93) BETA2FRACTION
                IF (IERR == 0) THEN
                  bParamFound(2) = .TRUE.
                ENDIF             
              ENDIF
            CASE ('ONSET')
C***          optional parameter, allows to start second beta law 
C***            further out than R_con (default)           
              IF (NPAR >= (i+1)) THEN
                TRYPAR = CURPAR(i+1)                      
                READ (TRYPAR, '(F10.0)', IOSTAT=IERR, ERR=94) RONSET
              ENDIF
          ENDSELECT
        ENDDO
      ENDIF
      
      
      DO i=1, 2
C***    One or more parameters have not been found => switch to old decoding
        IF (.NOT. bParamFound(i)) THEN
          WRITE (hCPR,*) '*** DEC2BETA: Old 2BETALAW decoding used'
          bOldDecode = .TRUE.
        ENDIF
      ENDDO
      
      
      IF (bOldDecode) THEN
C***     Old method: Pure position-based decoding, no keyword check      
         CALL SARGV (KARTE,3,TRYPAR)
         READ (TRYPAR, '(F10.0)', ERR=99) BETA2 
         CALL SARGV (KARTE,5,TRYPAR)
         READ (TRYPAR, '(F10.0)', ERR=99) BETA2FRACTION 
      ENDIF      
      
      
      
      RETURN
      

C***  FATAL ERROR CODES      
      
   92 WRITE(hCPR,'(A)') '*** DEC2BETA: CANNOT READ BETA2 IN:'
      WRITE (hCPR,*) KARTE
      STOP 'FATAL ERROR IN DEC2BETA'

   93 WRITE(hCPR,'(A)') '*** DEC2BETA: CANNOT READ FRACTION IN:'
      WRITE (hCPR,*) KARTE
      STOP 'FATAL ERROR IN DEC2BETA'

   94 WRITE(hCPR,'(A)') '*** DEC2BETA: CANNOT READ ONSET IN:'
      WRITE (hCPR,*) KARTE
      STOP 'FATAL ERROR IN DEC2BETA'
      
   99 WRITE (hCPR,*)
     >   '*** DEC2BETA: ERROR WHILE DECODING THE FOLLOWING CARDS-LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'FATAL ERROR IN DEC2BETA'
      
      END
      