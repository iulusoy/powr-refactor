      SUBROUTINE PLOTACC(PLOTOPT, AGRAV, AMECH, ARAD, APRESS, ACONT, 
     >                   ATHOM, WORKRATIO, VELO, RADIUS, ND, ATMEAN,
     >                   ENTOT, RNE, TAUROSS, RCON, T, TEFF, RSTAR, 
     >                   XMU, VTURB, XMSTAR, Rcritical, bFULLHYDROSTAT, 
     >                   bNoARAD, MODHEAD, JOBNUM, hPLOT, HYSTACC)
C******************************************************************************
C***  DIRECT TRANSFER OF HSUM PLOT
C***  TOTAL (FREQUENCY-INTEGRATED) FLUX versus DEPTH INDEX
C***  for both, continuum, line and sum of these
C******************************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 200
      INTEGER, PARAMETER :: NFINE = 10
      INTEGER, PARAMETER :: NDFINE = NDMAX * NFINE

      INTEGER, INTENT(IN) :: ND, JOBNUM, hPLOT

      REAL, DIMENSION(ND), INTENT(IN) :: AGRAV, AMECH, ARAD, APRESS, T,
     >                                   ACONT, ATHOM, VELO, VTURB, XMU,
     >                                   RADIUS, ENTOT, RNE, TAUROSS
      REAL, INTENT(IN) :: WORKRATIO, RCON, ATMEAN, XMSTAR, RSTAR, TEFF, 
     >                    Rcritical, HYSTACC
     
      LOGICAL, INTENT(IN) :: bFULLHYDROSTAT, bNoARAD

      CHARACTER(110) :: MODHEAD, HEADLINE
      CHARACTER(4) :: NormalizeTo
      CHARACTER(8) :: CENTER, CNORM, TEXT_X2
      CHARACTER(40) :: XTEXT
      CHARACTER(20) :: CUROPT, XAXISMODE, YLABEL
      CHARACTER PLOTOPT*(*)

      REAL, DIMENSION(NDMAX) :: X, Y1, Y2, Y3, Y4, Y5, VSCRATCH, RI
      REAL, DIMENSION(NDFINE) :: APRESSFINE, XFINE, YPFINE
      REAL, DIMENSION(ND) :: ANORM, VMACH, RHO

      INTEGER :: I, L, NPAR, NDIN, Lcand, Ltcand, 
     >           IFINE, NTOTFINE, ITOTFINE
      REAL :: XMIN, XMAX, YMIN, YMAX, YMINVAL, VINT, XINT,
     >        GEDDL, XLSTAR, XLSTARS, RNEINT, Xcrit,
     >        XABST, XTICK, XLENGTH, XOFF, XLIN, XICON, 
     >        Rsonic, Vsonic, Xsonic, Rts, Vts, Xts, XMG,
     >        RFINE, NTFINE, A2, DA2DR, DNTDR, GFINE, DR

      LOGICAL :: bBouretStyle, bDisplayFine

      INTEGER, EXTERNAL :: IDX

      !Physical constants
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)     
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)
      REAL, PARAMETER :: GCONST = 6.6727E-8     !Gravitational Constant (CGS)      

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTACC: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTACC: DIMENSION NDMAX INSUFFICIENT'
         WRITE (hCPR,'(2(A,I4))') '  NDMAX = ', NDMAX, ',   ND = ', ND
         WRITE (hCPR,'(A)') 'PLOTACC: PLOT SKIPPED'
         RETURN
      ENDIF
      
      CENTER = '\CENTER\'

C***  Decode possible options for x-axis mode
      XAXISMODE = 'DEPTHINDEX'
      NormalizeTo = 'AGRAV'
      bDisplayFine = .FALSE.
      bBouretStyle = .FALSE.
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR > 2) THEN
        DO I=3, NPAR
          CALL SARGV (PLOTOPT, I, CUROPT)
          SELECTCASE (CUROPT)
            CASE ('DEPTHINDEX', 'VELOCITY', 'TAUROSS', 'RADIUS',
     >            'L', 'R', 'V', 'ENTOT', 'RHO', 'TAU')
C***          simple keyword option mode            
              XAXISMODE = CUROPT
            CASE ('XAXISMODE', 'X')
              IF (NPAR >= (I+1)) THEN
                CALL SARGV (PLOTOPT, I+1, XAXISMODE)
              ENDIF
            CASE ('GEFF')
              NormalizeTo = 'GEFF'
            CASE ('BOURET')
              XAXISMODE = 'VELOCITY'
              NormalizeTo = 'THOM'
              bBouretStyle = .TRUE.
            CASE ('FINE')
              NormalizeTo = 'GRAV'
              bDisplayFine = .TRUE.
              bBouretStyle = .FALSE.
          ENDSELECT
        ENDDO
      ELSE 
C***    Fallback for unique labeling inside steal.plot / wruniq.plot
        PLOTOPT = PLOTOPT(:IDX(PLOTOPT)) // ' L'
      ENDIF

      IF (bBouretStyle .AND. bNoARAD) THEN
         WRITE (hCPR,'(A)') 'PLOTACC: WARNING - No ARAD calculated yet'
         WRITE (hCPR,'(A)') 'PLOTACC: *** PLOT ACC BOURET SKIPPED ***'
         RETURN
      ENDIF

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

      ANORM(ND) = 0.
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )      
        IF (NormalizeTo == 'GEFF') THEN
          !calculate depth-dependend GEFF 
          IF (bFULLHYDROSTAT) THEN
            !calculate geff with full a_rad
            GEDDL = MIN(ARAD(L)/AGRAV(L), 0.9)
          ELSE
            !calculate geff with 
            RNEINT = 0.5 * (RNE(L) + RNE(L+1))
            GEDDL = 10.**(-4.51) * (RNEINT/ATMEAN) * XLSTARS / XMSTAR
          ENDIF
          ANORM(L) = AGRAV(L) * ( 1. - GEDDL )
          CNORM = 'g&Teff&M'
        ELSEIF (NormalizeTo == 'THOM') THEN
          ANORM(L) = ATHOM(L)
          CNORM = 'g&Telec&M'
        ELSE
          ANORM(L) = AGRAV(L)
          CNORM = 'g'
        ENDIF          
      ENDDO

      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        VMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !v_mach in km/s
        VSCRATCH(L) = VELO(L) - VMACH(L)        
        IF ((VSCRATCH(L) < 0.) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,VMACH,RADIUS,ND)
      ENDIF

      !Find "true" sonic point (incl. VTURB)
      Ltcand = 0
      DO L=1, ND
        VSCRATCH(L) = VELO(L) - SQRT(VMACH(L)**2 + VTURB(L)**2)        
        IF ((VSCRATCH(L) < 0) .AND. (Ltcand == 0)) THEN
          Ltcand = L
        ENDIF
      ENDDO
      IF (Ltcand > 1) THEN
        CALL SPLINPOX(Rts,0.,RADIUS,VSCRATCH,ND,.FALSE.,Ltcand)
        CALL SPLINPOX(Vts,Rts,VELO,RADIUS,ND)
      ENDIF

      DO L=1, ND
        VSCRATCH(L) = VMACH(L)**2 + VTURB(L)**2
      ENDDO
      XMG = GCONST * XMSTAR * XMSUN
      DO L=1, ND
        IF (L > 1) THEN
          DR = RADIUS(L-1) - RADIUS(L)
          DO IFINE=1, NFINE
            ITOTFINE = (L-2)*NFINE + IFINE + 1
            RFINE = RADIUS(L-1) - FLOAT(IFINE)/FLOAT(NFINE) * DR
            CALL SPLINPOX(NTFINE,RFINE,ENTOT,RADIUS,ND,DFDX=DNTDR)
            CALL SPLINPOX(A2,RFINE,VSCRATCH,RADIUS,ND,DFDX=DA2DR)
            APRESSFINE(ITOTFINE) = - A2/NTFINE * DNTDR - DA2DR
            APRESSFINE(ITOTFINE) = APRESSFINE(ITOTFINE)*1.E10/RSTAR
            IF (APRESSFINE(ITOTFINE) > 0.) THEN
              GFINE = XMG/(RFINE*RSTAR)**2
              YPFINE(ITOTFINE) = LOG10(APRESSFINE(ITOTFINE)/GFINE)
            ENDIF
            SELECTCASE (XAXISMODE)
              CASE ('DEPTHINDEX', 'L')
                XFINE(ITOTFINE) = FLOAT(L-1) + FLOAT(IFINE)/FLOAT(NFINE)
              CASE ('RADIUS')
                IF (RFINE > 1.) THEN
                  XFINE(ITOTFINE) = LOG10(RFINE - 1.)
                ELSE
                  XFINE(ITOTFINE) = -100.
                ENDIF
              CASE ('VELOCITY', 'V')
                CALL SPLINPOX(XFINE(ITOTFINE),RFINE,VELO,RADIUS,ND)
                XFINE(ITOTFINE) = XFINE(ITOTFINE)/VELO(1)
              CASE ('ENTOT')
                XFINE(ITOTFINE) = LOG10(NTFINE)
              CASE ('RHO')
                XFINE(ITOTFINE) = LOG10(NTFINE * AMU * ATMEAN)
              CASE ('TAU', 'TAUROSS')
                CALL SPLINPOX(XFINE(ITOTFINE),RFINE,TAUROSS,RADIUS,ND)
            ENDSELECT
          ENDDO
        ELSE
          RFINE = RADIUS(L)
          CALL SPLINPOX(NTFINE,RFINE,ENTOT,RADIUS,ND,DFDX=DNTDR)
          CALL SPLINPOX(A2,RFINE,VSCRATCH,RADIUS,ND,DFDX=DA2DR)
          APRESSFINE(1) = - A2/NTFINE * DNTDR - DA2DR        
          APRESSFINE(1) = APRESSFINE(1)*1.E10/RSTAR
          IF (APRESSFINE(1) > 0.) THEN
            GFINE = XMG/(RFINE*RSTAR)**2
            YPFINE(1) = LOG10(APRESSFINE(1)/GFINE)
          ENDIF
          SELECTCASE (XAXISMODE)
            CASE ('DEPTHINDEX', 'L')
              XFINE(1) = 1.
            CASE ('RADIUS', 'R')
              XFINE(1) = LOG10(RFINE - 1.)
            CASE ('VELOCITY', 'V')
              XFINE(1) = 1.
            CASE ('ENTOT')
              XFINE(1) = LOG10(NTFINE)
            CASE ('RHO')
              XFINE(1) = LOG10(NTFINE * AMU * ATMEAN)
            CASE ('TAU', 'TAUROSS')
              XFINE(1) = -100.
          ENDSELECT
        ENDIF
      ENDDO
      NTOTFINE = (ND-1)*NFINE + 1
      
      IF (XAXISMODE == 'DEPTHINDEX' .OR. XAXISMODE == 'L') THEN 
C***     X-Axis = Depth Index
         XTEXT = CENTER//'DEPTH INDEX L'
         XMIN = 0.
         XMAX = FLOAT(ND)
         TEXT_X2 = ' XMAX '
         XTICK = 5.
         XABST = 10. 
         DO L=1, ND-1
            X(L) = FLOAT(L) + 0.5
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF
         IF (Ltcand > 1 .AND. Ltcand < ND) THEN
           CALL SPLINPOX(Xts,Rts,X,RI,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'VELOCITY' .OR. XAXISMODE == 'V') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v' // CHAR(92) // '8'
         XMIN = 0.
         XMAX = 1.
         TEXT_X2 = ' XMIN '
         XTICK = 0.1
         XABST = 0.5 
         DO L=1, ND-1
           CALL SPLINPOX(VINT, RI(L), VELO, RADIUS, ND)
           X(L) = VINT/VELO(1)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = Vsonic/VELO(1)
         ENDIF
      ELSEIF (XAXISMODE == 'RADIUS' .OR. XAXISMODE == 'R') THEN
         XTEXT = CENTER//'LOG (R/R\*-1)'
         XMIN = LOG10 ( 0.9 * (RI(ND-1) - 1.) )
         XMAX = LOG10 ( RADIUS(1) - 1. )
         XTICK = 0.2
         XABST = 1.0 
         DO L=1, ND-1
           X(L) = LOG10( RI(L) - 1. )
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = LOG10( Rsonic - 1. )
         ENDIF
         IF (Ltcand > 1 .AND. Ltcand < ND) THEN
           Xts = LOG10( Rts - 1. )
         ENDIF
      ELSEIF (XAXISMODE == 'TAUROSS' .OR. XAXISMODE == 'TAU') THEN
         XTEXT = CENTER//'log #t#&TRoss&M'
         XMIN = LOG10(0.4 * TAUROSS(2))
         XMAX = LOG10( TAUROSS(ND) )
         XTICK = 0.1
         XABST = 1.0 
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), TAUROSS, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF
         IF (Ltcand > 1 .AND. Ltcand < ND) THEN
           CALL SPLINPOX(Xts,Rts,X,RI,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'ENTOT' .OR. XAXISMODE == 'N') THEN
         XTEXT = CENTER // 'log(&Rn&N&Ttot&M/cm&H-3&M)'
         XMIN = LOG10(ENTOT(1))
         XMAX = LOG10(ENTOT(ND))
         XTICK = 1.
         XABST = 3.
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), ENTOT, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RADIUS,ND-1)         
         ENDIF      
         IF (Ltcand > 1 .AND. Ltcand < ND) THEN
           CALL SPLINPOX(Xts,Rts,X,RADIUS,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'RHO' .OR. XAXISMODE == 'DENS') THEN
         XTEXT = CENTER // 'log(&R#r#&N/(g cm&H-3&M))'
         DO L=1, ND
           RHO(L) = ENTOT(L) * AMU * ATMEAN
         ENDDO         
         XMIN = LOG10(RHO(1))
         XMAX = LOG10(RHO(ND))
         XTICK = 1.
         XABST = 2.
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), RHO, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF           
         IF (Ltcand > 1 .AND. Ltcand < ND) THEN
           CALL SPLINPOX(Xts,Rts,X,RI,ND-1)         
         ENDIF
      ELSE
         WRITE (hCPR,*) 'PLOTACC: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))
         RETURN
      ENDIF


      CALL JSYMSET ('G2','TRANSFER')

      HEADLINE = 'ACC:'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM

      WRITE (hPLOT, '(A)') '*NEXTPLOT: ' // PLOTOPT
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 5 1.0 0.58 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 7 0.0 0.50 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 9 0.6 0.60 0.6'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 8 0.6 0.8 1.0'
      WRITE (hPLOT, '(A)') '\PEN=1'     !needed to get rid of definecolor pen bug 

      WRITE (hPLOT, '(A)') '\COLOR=3'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0 XMAX 0 0 0'

      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON >= RADIUS(ND) .AND. RCON <= RADIUS(1)) THEN
        CALL SPLINPOX(XICON, RCON, X, RI, ND)
c        IF (XICON >= 0.005) THEN
        IF (XICON >= XMIN + 0.1 * XTICK) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',XICON,' YMIN ',XICON,' YMAX 0. 0. SYMBOL=10'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',XICON,' YMAX 0. -0.2 0.2 -90  R&Tcon'
        ENDIF
      ENDIF
      IF (Lcand > 1 .AND. Lcand < ND) THEN
        IF (ABS(Xsonic) >= 0.005) THEN
          WRITE (hPLOT, '(A,F12.5,A,F12.5,A)') 
     >      '\LINUN ',Xsonic,' YMIN ',Xsonic,
     >                    ' YMAX 0. 0. SYMBOL=20 SIZE=0.05'
          WRITE (hPLOT, '(A,F12.5,A)') 
     >      '\LUNA ',Xsonic,' YMAX 0. -0.2 0.2 -90  R&Tsonic'
        ENDIF
      ENDIF
      IF (Ltcand > 1 .AND. Ltcand < ND) THEN
        IF (ABS(Xts) >= 0.005) THEN
          WRITE (hPLOT, '(A,F12.5,A,F12.5,A)') 
     >      '\LINUN ',Xts,' YMIN ',Xts,
     >                    ' YMAX 0. 0. SYMBOL=20 SIZE=0.025'
          WRITE (hPLOT, '(A,F12.5,A)') 
     >      '\LUNA ',Xts,' YMAX 0. -0.2 0.2 -90  R&Tts'
        ENDIF
      ENDIF
      IF (Rcritical >= RI(ND-1) .AND. Rcritical <= RI(1)) THEN
        CALL SPLINPOX(Xcrit, Rcritical, X, RI, ND-1)
        IF (Xcrit >= 0.005) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',Xcrit,' YMIN ',Xcrit,' YMAX 0. 0. SYMBOL=9'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',Xcrit,' YMAX 0. -0.2 0.2 -90  R&Tcrit'
        ENDIF
      ENDIF
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'

      IF (HYSTACC .GT. 0.) THEN
         WRITE (hPLOT, '(A)') '* HYDROSTATIC INTEGRATION EPSILON'
         WRITE (hPLOT, '(A)') '\COLOR = 8'
         WRITE (hPLOT, '(A,F6.3,1X,F6.3,A,F6.3,A)')
     >    '\RECT ', XICON, ALOG10(1. - HYSTACC), TEXT_X2, ALOG10(1. + HYSTACC),
     >    '  0 0 FILLED'
      ENDIF


      WRITE (hPLOT, '(A)') '\COLOR=1'

      WRITE (hPLOT, '(A,F6.3)') 
     >   'KASDEF LUN XMIN YMIN 1. 1. .5 ' //
     >    '&EWORK Ratio Rad.+Gaspress. / Mech.+Grav. =&N', WORKRATIO  
      WRITE (hPLOT, '(A)') '\INBOX'

      XLENGTH = 2.5
      IF (bBouretStyle) THEN
        XLIN = 2.
        XOFF = XLENGTH + XLIN + 0.5
        
C***    ACC plot in the coloring from Bouret et al. (2012)
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -1. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF, ' M-1.'
     >   // ' 0.3 (a&Tmech&M+g)/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH,' 0 ',XLIN,' -1.5 SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF, ' M-1.5'
     >   // ' 0.3 (a&Tmech&M-a&Tpress&M+g)/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=4'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -2. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-2.0 0.3 a&Trad&M/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=7'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -2.5 SYMBOL=9'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-2.5 0.3 a&Tcont&M/' // CNORM     
      ELSE
        XLIN = 14.
        XOFF = XLENGTH + XLIN + 0.5
        
C***    Standard ACC plot colors           
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -1. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF, ' M-1.'
     >   // ' 0.3 (' // CNORM(:IDX(CNORM)) // '+a&Tmech&M)/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -1.5 SYMBOL=9'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-1.5 0.3 a&Trad&M/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=4'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -2. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-2. 0.3 a&Tmech&M/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=7'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -2.5 SYMBOL=9'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-2.5 0.3 a&Tcont&M/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -3 SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN ' //
     >   'XMIN YMAX ', XOFF, ' M-3.0 0.3 (a&Trad&M+a&Tpress&M)/' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=5'
        WRITE (hPLOT, '(A,2(F4.1,A))') CHAR(92) // 
     >   'LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -3.5 '
     >                                       // ' SYMBOL=9 SIZE=0.05'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN XMIN YMAX ', XOFF, ' M-3.5 0.3 a&Tpress&M/' // CNORM

        WRITE (hPLOT, '(A)') '\COLOR=6'
        WRITE (hPLOT, '(A,2(F4.1,A))') 
     >   '\LINREL XMIN YMAX ',XLENGTH, ' 0 ', XLIN, ' -4.0  SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF,
     >   ' M-4.0 0.3 a&Tthom&M/' // CNORM
      ENDIF
      
      WRITE (hPLOT, '(A)') '\COLOR=1'

C***  HTOT: Conversion into radiation temperatures
C***  Note: HTOT may not be calculated if STEAL is used for OUTPUT ONLY
      YMIN = 0.
      YMAX = 0.
      DO L=1, ND-1
        Y1(L) = (ANORM(L) + AMECH(L))/ANORM(L)
        Y2(L) = ARAD(L)/ANORM(L)
        Y3(L) = AMECH(L)/ANORM(L)
        Y4(L) = ACONT(L)/ANORM(L)
        Y5(L) = APRESS(L)/ANORM(L)
        IF (.NOT. bBouretStyle) THEN
          IF (Y1(L) > .0) Y1(L) = LOG10(Y1(L))
          IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
          IF (Y3(L) > .0) Y3(L) = LOG10(Y3(L))
          IF (Y4(L) > .0) Y4(L) = LOG10(Y4(L))
          IF (Y5(L) > .0) Y5(L) = LOG10(Y5(L))
        ELSE
          Y1(L) = (AGRAV(L) + AMECH(L))/ANORM(L)
        ENDIF
      ENDDO

C***  Restrict blue curve (AMECH) to inbox area
      YMINVAL = MIN(MINVAL(Y2), MINVAL(Y4))
      YMINVAL = MIN(YMINVAL, MINVAL(Y5))
      DO L = 1, ND-1
        NDIN = L
        IF (Y3(L) < YMINVAL) EXIT
      ENDDO

      IF (bBouretStyle) THEN
        YLABEL = 'a/g&Telec&M         '
      ELSE
        YLABEL = 'log&T10&M(a/g)      '
      ENDIF      
      
      CALL PLOTANF (hPLOT,HEADLINE, '&E'//HEADLINE,
     $        XTEXT, CENTER//YLABEL,
     >        0., XMIN, XMAX, XTICK, XABST, 0.,
     >        0., YMIN, YMAX, .1,  0.2, 0.,
     $        X, Y1, ND-1, 5)

      IF (bBouretStyle) THEN
C***    ACC plot in the coloring from Bouret et al. (2012)
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=4') 
        
        DO L=1, ND-1
          Y2(L) = (AMECH(L)-APRESS(L)+AGRAV(L))/ANORM(L)
        ENDDO
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2') 
C***    Also show pure continuum contribution
        CALL PLOTCONS (hPLOT,X,Y4,ND-1,'COLOR=7 SYMBOL=9 SIZE=0.1') 
        
      ELSE
C***    Standard ACC plot colors     
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2 SYMBOL=9 SIZE=0.1') 
        CALL PLOTCONS (hPLOT,X,Y3,NDIN,'COLOR= 4') 
        CALL PLOTCONS (hPLOT,X,Y4,ND-1,'COLOR=7 SYMBOL=9 SIZE=0.1') 

        DO L=1, ND-1
          Y2(L) = (ARAD(L)+APRESS(L))/ANORM(L)
          IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
        ENDDO
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2') 

        CALL PLOTCONS (hPLOT,X,Y5,ND-1,'COLOR=5 SYMBOL=9 SIZE=0.05') 

        DO L=1, ND-1
          Y2(L) = ATHOM(L)/ANORM(L)
          IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
        ENDDO
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=6') 

        IF (bDisplayFine) THEN
          CALL PLOTCONS (hPLOT,XFINE,YPFINE,NTOTFINE,
     >                               'COLOR=6 SYMBOL=9 SIZE=0.05') 
        ENDIF
      ENDIF
      
      RETURN
      END
