      SUBROUTINE DECCOLI (LSOPA,LSINT,MODHIST,
     >             NOLAP,LPOPAB,LPOPABD,LEVDEBUG,LPJNUE,
     >             LPJNUED,LASERV,PARALAS,LPSNUE,LPSNUED,
     $             MAXPLOT, RANGE1, RANGE2, EXLAM1, EXLAM2,MAXEXT,
     $             BLLIST, 
     >             REDISMODE, NEWWRC, BCOLIRAY, BKUDRITZKI,
     >             CLHLP, BITCONT, BPLOT,  
     >             IPLOT, LPLOT, ND, OPC, 
     >             IVERS, BEMIX, EMIXSTART, BEMIXFIX, EMIXFIX, 
     >             IVERS_FE_EXPFAC, BPLOTALPHA, ALPOPT, POPMIN,
     >             XLAM_FINE_START, XLAM_FINE_END, LPLOT_WCHARM, XLP1,
     >             XLP2, GAMMACOLI, GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2,
     >             bKALPHA, bHYDROSOLVE, VELOMODFAK, iTypeAlpha, 
     >             bForceCOLIP, bTDIFFUS, bPLOTRTAU1, bNoIronLaser, 
     >             iHTOTCUT, bMAXEMIX, bALOTri,
     >             bNoNEGEDDIG, bDDVDOP, LPRDH, 
     >             bDEBUG, CUTOPAMEANTHRES, DRLINES_CARD, EMIXMAX)
C*******************************************************************************
C***  DECODES INPUT CARDS FOR MAIN PROGRAM "COLI"
C*******************************************************************************

      !IMPLICIT NONE
 
      CHARACTER(120) :: KARTE, DRLINES_CARD, ALPOPT
      CHARACTER(4) :: REDISMODE
      LOGICAL :: NOLAP, LASERSET, DRNORUD, BLLIST, 
     >           BCOLIRAY, CLHLP, BITCONT, 
     >           BPLOT, BKUDRITZKI,
     >           BEMIX, BEMIXFIX, BPLOTALPHA, bCustomVMOD
      INTEGER :: IVERS_FE_EXPFAC, iHTOTCUT
      INTEGER, DIMENSION(MAXPLOT) :: LPOPAB, LPOPABD, LPJNUE, LPJNUED,
     >                               LPSNUE, LPSNUED
      REAL, INTENT(INOUT) :: VELOMODFAK, POPMIN
      INTEGER, INTENT(INOUT) :: iTypeAlpha
      REAL :: PARALAS, tempREAL, CUTOPAMEANTHRES
      CHARACTER(16), DIMENSION(20) :: ACTPAR
      CHARACTER(8) :: OPC
      LOGICAL, INTENT(INOUT) :: bKALPHA, bForceCOLIP, bTDIFFUS, bDDVDOP
      LOGICAL, INTENT(OUT) :: bHYDROSOLVE, bPLOTRTAU1
      LOGICAL :: bReadPOPMIN, bNoIronLaser, bMAXEMIX, 
     >           bNoNEGEDDIG, bDEBUG, bALOTri
       
      INTEGER :: LSOPA, LSINT, IPLOT, LPLOT, LPRDH
      REAL :: RANGE1, RANGE2


      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


C***  DEFAULT VALUES
      LSOPA   = -1
      LSINT   = -1
      NOLAP   = .FALSE.
      LASERSET= .FALSE.
      LEVDEBUG= 0
      LASERV  = 1
      PARALAS = 0.01
      BLLIST  = .FALSE.
      REDISMODE = 'CONT'
      BCOLIRAY = .FALSE.
      CLHLP = .FALSE.
      BITCONT = .TRUE.
      BPLOT = .FALSE.
      IPLOT = ND / 2
      LPLOT = 10000
      OPC = 'NONE'
      IVERS = 4
      BPLOTALPHA = .FALSE.
      bHYDROSOLVE = .FALSE.
      bTDIFFUS = .TRUE.
      bCustomVMOD = .FALSE.
      iTypeAlpha = 0
      VELOMODFAK = 0.9          !default value for velocity field modification factor
      bReadPOPMIN = (POPMIN < 1.E-99)
      IF (bReadPOPMIN) THEN
        POPMIN = 1.E-25           !default value should be the same as in STEAL -> DECSTE
      ENDIF
      bPLOTRTAU1 = .FALSE.
      bNoIronLaser = .FALSE.
      bALOTri = .FALSE.
      iHTOTCUT = 0
      bNoNEGEDDIG = .FALSE.     !default: Do not reset negative EDDIG
      bDDVDOP = .FALSE.         !default: No depth-dependent VDOP
      LPRDH = 0
      bDEBUG = .FALSE.
      CUTOPAMEANTHRES = 0.      !default: Do not cut out strong lines from mean opacities
    
C***  Version of exponential factor in Fe lines (rates and eta) 
      IVERS_FE_EXPFAC = 1

C***  NUMBER OF REPEAT ITERATIONS BETWEEN TWO COLI+ JOBS
C***   This default must agree with the definition in DECSTE (STEAL)
      NEWWRC=6

C***  EDDIMIX enabled
      BEMIX = .TRUE.
      EMIXSTART = 1.0
      BEMIXFIX = .FALSE.
      EMIXFIX = 1.0
      bMAXEMIX = .FALSE.
      EMIXMAX = -1.0

      DO I=1,MAXPLOT
        LPOPAB (I) = 0
        LPOPABD(I) = 0
        LPJNUE (I) = 0
        LPJNUED(I) = 0
        LPSNUE (I) = 0
        LPSNUED(I) = 0
      ENDDO

      NOPLOT = 0
      NJPLOT = 0
      NSPLOT = 0

      RANGE1 = -1.
      RANGE2 = -1.

      EXLAM1 = -1.
      EXLAM2 = -1.

C***  Default values for extension of the fine grid (for Integration in STEALCL)
      XLAM_FINE_START =  100.
      XLAM_FINE_END   = 2000.

C***  Default-Range which is plotted
      XLP1 = 100.
      XLP2 = 1000.

C***  Special PLOT of WCHARM factors; If this set to zero, no plot is prepared
      LPLOT_WCHARM = 0

c      GAMMACOLI = 0.  !this will ALWAYS lead to a crash
      GAMMACOLI = 1.

C***  Parameters for Unsoeld-Lucy method, needed by FREQUINT
      GAMMAT = 0.
      UNLU_TAUMAX = 1000.
      UNLU_TAUMAX2 = 100.
      
      ALPOPT = ""

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
    6 READ (1,FMT='(A)', END=100) KARTE
 
      CALL SARGC(KARTE,NPAR)
      IF ( NPAR < 1) GOTO 6
      IF ( NPAR > 20) NPAR = 20

      DO I=1, NPAR
        CALL SARGV(KARTE,I,ACTPAR(I))
      ENDDO

      IF (ACTPAR(1) .EQ. 'PRINT') THEN
C                         =====
        IF ( ACTPAR(2) .EQ. 'INTL') THEN
C                            ====
            READ (ACTPAR(3),8, ERR=99) XL
    8       FORMAT (F10.0)
            LSINT=IFIX(XL)
            IF (LSINT.EQ.0) LSINT=1

        ELSE IF ( ACTPAR(2) .EQ. 'OPAL') THEN
C                                 ====
            READ (ACTPAR(3),8, ERR=99) XL
            LSOPA=IFIX(XL)
            IF (LSOPA.EQ.0) LSOPA=1

        ELSE IF ( ACTPAR(2) .EQ. 'LINELIST') THEN
C                                 ========
            BLLIST = .TRUE.
        ELSE IF ( ACTPAR(2) .EQ. 'DELTAH') THEN
C                                 ======
            READ (ACTPAR(3),8, ERR=99) XL
            LPRDH=IFIX(XL)
            IF (LPRDH == 0) LPRDH = 1
        ENDIF
        GOTO 6
      ENDIF

      IF ( ACTPAR(1)(:6) .EQ. 'DRLINE' ) THEN
C                                  ======
         DRLINES_CARD = KARTE
 
      ELSEIF (ACTPAR(1) .EQ. 'CMFDLEV') THEN
C                             =======
        READ(ACTPAR(2),'(I10)', ERR=99) LEVDEBUG
        GOTO 6
      ENDIF

      IF (ACTPAR(1) == 'OPAMEANTHRES') THEN
        READ(ACTPAR(2), '(G10.3)', ERR=99) CUTOPAMEANTHRES
        GOTO 6
      ENDIF

      IF (ACTPAR(1) .EQ. 'LASERV') THEN
C                         ======
        LASERSET = .TRUE.
        READ (ACTPAR(2),'(I10)', ERR=99) LASERV
        IF (NPAR == 3) THEN
          READ (ACTPAR(3),'(G10.1)', ERR=99) PARALAS
        ELSEIF (NPAR > 3) THEN
          DO IPAR=3, NPAR-1
            SELECTCASE (ACTPAR(IPAR))
              CASE ('LINES', 'LINE', 'L')
                IF (NPAR >= (IPAR+1)) THEN
                  READ (ACTPAR(IPAR+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    PARALAS = tempREAL
                  ENDIF                    
                ENDIF
cc              CASE ('KONT', 'CONT', 'K', 'C')
cc                IF (NPAR >= (IPAR+1)) THEN
cc                  READ (ACTPAR(IPAR+1), '(F10.0)', IOSTAT=IERR) tempREAL
cc                  IF (IERR == 0) THEN
cc                    PARALKONT = tempREAL
cc                  ENDIF                    
cc                ENDIF
            ENDSELECT
          ENDDO
        ENDIF
        GOTO 6
        ENDIF

      IF (ACTPAR(1) .EQ. 'EXTEND') THEN
C                         ======
        READ (ACTPAR(2), '(F10.0)', ERR=99) EXLAM1
        READ (ACTPAR(3), '(F10.0)', ERR=99) EXLAM2
        GOTO 6
      ENDIF

      IF (ACTPAR(1) .EQ. 'REDISMODE') THEN
C                         =========
         REDISMODE = ACTPAR(2) 
         GOTO 6
      ENDIF

      IF (ACTPAR(1) .EQ. 'NEWWRC' ) THEN
C                         ======
         READ (ACTPAR(2),'(F10.0)',ERR=99) XL
         NEWWRC=IFIX(XL)
      ENDIF

      IF (ACTPAR(1) .EQ. 'OB-VERS') THEN
C                         =======
         IF (NPAR .EQ. 2 .OR. ACTPAR(3) .EQ. 'COLI') THEN
           READ (ACTPAR(2),'(I10.0)', ERR=99) IVERS
         ENDIF
      ENDIF

      IF (KARTE(:4) == 'VDOP') THEN
C                       ====
         IF (NPAR > 2) THEN
           IF (ACTPAR(3) == 'AUTO') bDDVDOP = .TRUE.
         ENDIF
         GOTO 6
      ENDIF
      
      IF (ACTPAR(1) == 'NOTDIFFUS') THEN
C                       =========
         bTDIFFUS = .FALSE.
      ENDIF
         
C********************  CMF PLOT OPTION  **********************************

      IF (ACTPAR(1) .EQ. 'CMFPLOT') THEN
C                         =======
        IF (NPAR .LT. 2) THEN
         PRINT *,'CMFPLOT: MISSING OPTION(S)'
         GOTO 6
         ENDIF
C***  CMFPLOT RANGE lambda1 lambda2
C***    - WITH THIS OPTION, ONLY ONE PLOT OF EACH QUANTITY (J, KAPPA, S) 
C***     IS POSSIBLE
        IF (ACTPAR(2) .EQ. 'RANGE') THEN
C                           =====
           MAXPLOT = 1
           READ (ACTPAR(3), '(F10.0)', ERR=99) RANGE1
           READ (ACTPAR(4), '(F10.0)', ERR=99) RANGE2
           GOTO 6
           ENDIF

        IF (ACTPAR(2) .EQ. 'OPAB') THEN
C                           ====
         IF (NOPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NOPLOT=NOPLOT+1
         READ(ACTPAR(3),'(I10)', ERR=99) LPOPAB(NOPLOT)
         READ(ACTPAR(4),'(I10)', ERR=99) LPOPABD(NOPLOT)
         GOTO 6
         ENDIF

        IF (ACTPAR(2) .EQ. 'JNUE') THEN
C                           ====
         IF (NJPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NJPLOT = NJPLOT + 1
         READ (ACTPAR(3),'(I10)', ERR=99) LPJNUE(NJPLOT)
         READ (ACTPAR(4),'(I10)', ERR=99) LPJNUED(NJPLOT)
         GOTO 6
         ENDIF

        IF (ACTPAR(2) .EQ. 'SNUE') THEN
C                           ====
         IF (NSPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NSPLOT=NSPLOT+1
         READ (ACTPAR(3), '(I10)', ERR=99) LPSNUE(NSPLOT)
         READ (ACTPAR(4), '(I10)', ERR=99) LPSNUED(NSPLOT)
         GOTO 6
         ENDIF

        ENDIF

        IF (ACTPAR(1) .EQ. 'PLOT') THEN
C                           ====
           IF (ACTPAR(2) == 'ALPHA') THEN
              BPLOTALPHA = .TRUE.
              ALPOPT = KARTE
              GOTO 6
           ENDIF

           IF (ACTPAR(2) == 'RTAU1COLI') THEN
              bPLOTRTAU1 = .TRUE.
              GOTO 6
           ENDIF
           
        ENDIF

        IF (ACTPAR(1) .EQ. 'PURE' .AND. ACTPAR(2) .EQ. 'COLIRAY') THEN
C                           ====                       =======
         BCOLIRAY = .TRUE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'OPC') THEN
C                           ===
         OPC = ACTPAR(2)
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'COLI') THEN
C                           ====
          IF (ACTPAR(2) .EQ. 'OUTPUT' .AND. ACTPAR(3) .EQ. 'ONLY') THEN
            CLHLP = .TRUE.
            GOTO 6
          ENDIF

          IF (ACTPAR(2) .EQ. 'PLOT') THEN
            BPLOT = .TRUE.
            IF (NPAR .GE. 3) THEN
              READ (ACTPAR(3),'(I10)', ERR=99) IPLOT
            ENDIF
            IF (NPAR .GE. 4) THEN
              READ (ACTPAR(4),'(I10)', ERR=99) LPLOT
            ENDIF
            GOTO 6

          ELSE IF (ACTPAR(2) .EQ. 'PLOT_WCHARM') THEN
            READ (ACTPAR(3),'(I10)', ERR=99) LPLOT_WCHARM

          ELSE IF (ACTPAR(2) .EQ. 'UNLUPAR') THEN
            DO IPAR=3, NPAR-1
               IF (ACTPAR(IPAR) .EQ. 'GAMMAT') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) GAMMAT
               ELSE IF (ACTPAR(IPAR) .EQ. 'TAUMAX') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) UNLU_TAUMAX
               ELSE IF (ACTPAR(IPAR) .EQ. 'TAUMAX2') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) UNLU_TAUMAX2
               ENDIF
            ENDDO
          ENDIF
          IF (ACTPAR(2) .EQ. 'COLI+') THEN
C                             =====
            NEWWRC = 1
            bForceCOLIP = .TRUE.
            GOTO 6
          ENDIF
          IF (ACTPAR(2) .EQ. 'PLOTRANGE') THEN
            READ (ACTPAR(3),'(F10.0)', ERR=99) XLP1
            READ (ACTPAR(4),'(F10.0)', ERR=99) XLP2
            GOTO 6
          ENDIF
          IF (ACTPAR(2) .EQ. 'GAMMA') THEN
            READ (ACTPAR(3),'(F10.0)', ERR=99) GAMMACOLI
            GOTO 6
          ENDIF
          IF (ACTPAR(2) == 'KUDRITZKI') THEN
            BKUDRITZKI = .TRUE.
            GOTO 6
          ENDIF
          IF (ACTPAR(2) == 'NONEGACC') THEN
            iHTOTCUT = 1
            GOTO 6
          ENDIF
          IF (ACTPAR(2) == 'NO' .AND. ACTPAR(3)(1:3) == 'NEG'
     >           .AND. ACTPAR(4) == 'EDDIG') THEN
            bNoNEGEDDIG = .TRUE.
            GOTO 6
          ENDIF          
          IF (ACTPAR(2) == 'NOFELASER') THEN
            bNoIronLaser = .TRUE.
            WRITE (0,*) 'COLI: NO IRON LASER ALLOWED'
            GOTO 6
          ENDIF
          IF (ACTPAR(2) == 'DEBUG') THEN
            bDEBUG = .TRUE.
            GOTO 6
          ENDIF
        ENDIF
C***  End of Options beginning with COLI ...

        IF (KARTE(:22) == 'NO CONTINUUM ITERATION') THEN
C                          ======================
         BITCONT = .FALSE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'NO' .AND. 
     >      ACTPAR(2) .EQ. 'EDDIMIX') THEN
C                           =======
         BEMIX = .FALSE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'EDDIMIX' .AND. 
     >      ACTPAR(2) .EQ. 'START') THEN
C                           =======
         READ (ACTPAR(3),'(F10.0)', ERR=99) EMIXSTART
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'EDDIMIX' .AND. 
     >      ACTPAR(2) .EQ. 'FIX') THEN
C                           =======
         BEMIXFIX = .TRUE.
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) EMIXFIX
         ENDIF
         GOTO 6
         ENDIF

        IF (ACTPAR(1) == 'EDDIMIX' .AND.  ACTPAR(2) == 'MAX') THEN
C                         =======                       ===
          bMAXEMIX = .TRUE.
          IF (NPAR .LT. 3) THEN
            WRITE (0,*) 'ERROR: EDDIMIX MAX needs a value'
            GOTO 99
          ENDIF

          READ (ACTPAR(3),'(F10.0)', ERR=99) EMIXMAX 
          GOTO 6
        ENDIF

        IF (ACTPAR(1)(:16) .EQ. 'IRONLINES-EXPFAC') THEN
C                                ================
           IF (ACTPAR(2) .EQ. 'OFF') THEN
              IVERS_FE_EXPFAC = 0
              WRITE (0,*) '*** WARNING: non-standard ' // KARTE
           ELSEIF (ACTPAR(2) == 'TEFF') THEN
              IVERS_FE_EXPFAC = 1
           ELSEIF (ACTPAR(2) == 'TRADITIONAL') THEN
              IVERS_FE_EXPFAC = 2
              WRITE (0,*) '*** WARNING: non-standard ' // KARTE
           ELSE
              WRITE (0,*) '*** WARNING: Invalid option ' // KARTE
           ENDIF
          GOTO 6
        ENDIF


        IF (ACTPAR(1)(1:6) == 'POPMIN' .AND. bReadPOPMIN) THEN
C                              ======
          READ (ACTPAR(2),'(F10.0)', ERR=99) POPMIN
          GOTO 6
        ENDIF         
         
        !Force calculation of force multipliers
        ! Note: Usually this is automatically forced from the HYDRO interval option
        !       By using the FORCEMULTIPLIERS card you can ensure that the calculation
        !       is done at every COLI even if the HYDRO calculations are never performed
        !Full syntax is FORCEMULTIPLIERS but ACTPAR currently only reads up to 14 chars
        IF (ACTPAR(1)(1:9) == 'FORCEMULT') THEN
C                              =========
         bKALPHA = .TRUE.
         IF (.NOT. bHYDROSOLVE .AND. iTypeAlpha == 0) THEN
           iTypeAlpha = 1
         ENDIF
         IF (NPAR > 2) THEN
           DO I=2, NPAR
             IF ((ACTPAR(I) == 'VMOD') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F20.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 VELOMODFAK = tempREAL
                 bCustomVMOD = .TRUE.
               ENDIF               
             ENDIF
             IF ((ACTPAR(I) == 'AT') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F10.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 iTypeAlpha = IFIX(tempREAL)
               ENDIF        
             ENDIF
           ENDDO
         ENDIF
         GOTO 6
        ENDIF

        IF (ACTPAR(1) == 'HYDRO') THEN
C                         =====
         bHYDROSOLVE = .TRUE.
         IF (NPAR > 2) THEN
           DO I=2, NPAR
             IF ((ACTPAR(I) == 'VMOD') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F20.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 VELOMODFAK = tempREAL
                 bCustomVMOD = .TRUE.
               ENDIF               
             ENDIF
             IF ((ACTPAR(I) == 'AT') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F10.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 iTypeAlpha = IFIX(tempREAL)
               ENDIF        
             ENDIF
           ENDDO
         ENDIF
         GOTO 6
        ENDIF

        IF (ACTPAR(1) .EQ. 'XJLAPP' .AND.
     >         ACTPAR(2) .EQ. 'COLI') THEN
          DO I=3, NPAR
            IF (ACTPAR(I)(:3) .EQ. 'TRI') THEN
              bALOTri = .TRUE.
            ENDIF
          ENDDO
        ENDIF

        IF (ACTPAR(1) .EQ. 'XJLAPP' .AND.
     >         ACTPAR(2) .EQ. 'NEW') THEN
C                              ===========
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=99) XLAM_FINE_END
         ENDIF 
        ENDIF

        IF (ACTPAR(1) .EQ. 'XJCAPP' .AND.
     >         ACTPAR(2) .EQ. 'NEW') THEN
C                              ===========
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=99) XLAM_FINE_END
         ENDIF 
        ENDIF

      GOTO 6


C***  END-OF-DATA REACHED: REGULAR EXIT
  100 CONTINUE
      CLOSE(1)
  
C***  Check if Alpha-Type is in valid range:
      IF (iTypeAlpha < 0 .OR. 
     >     (iTypeAlpha > 3  .AND. iTypeAlpha < 10) .OR.
     >     (iTypeAlpha > 12 .AND. iTypeAlpha /= 31)) THEN
        WRITE (hCPR,*) 'INVALID ALPHA TYPE SPECIFIED!'
        STOP 'ERROR in DECCOLI'
      ENDIF

C***  Default is different if response factor is calculated instead of alpha (AT == 3)
      IF (iTypeAlpha == 3 .AND. (.NOT. bCustomVMOD)) THEN
c        VELOMODFAK = 0.99   !more accurate default for response factor calculation
        VELOMODFAK = 0.9   !try this first
      ENDIF
       
C***  IF LASER VERSION NOT EXPLICITELY DECLARED, SET DEFAULT VALUE
      IF (.NOT. LASERSET  .AND. NOLAP) LASERV = 0

      RETURN

C***  ERROR EXIT:  ************************************************
   99 CONTINUE
      WRITE (hCPR,*) 'ERROR WHEN DECODING THE FOLLOWING LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'

      END
