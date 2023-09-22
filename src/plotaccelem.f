      SUBROUTINE PLOTACCELEM(PLOTOPT, AGRAV, AMECH, ARAD, ARADELEM,
     >                       ACONTELEM, ARADION, ACONTION, 
     >                       APRESS, ACONT, ATHOM, WORKRATIO, VELO, 
     >                       RADIUS, ND, NATOM, MAXION, ATMEAN, ENTOT, 
     >                       RNE, TAUROSS, RCON, T, TEFF, RSTAR, XMU, 
     >                       XMSTAR, Rcritical, bFULLHYDROSTAT, SYMBOL, 
     >                       ELEMENT, bNoARAD, MODHEAD, JOBNUM, hPLOT)
C***********************************************************************
C***  PLOT ACCELERATION CONTRIBUTIONS FROM DIFFERENT ELEMENTS OR IONS
C***  LOG ACC/G versus any depth-dependent quantity (DEPTH INDEX, R, V, ...)
C***  for both, continuum, line and sum of these
C***
C***  This routine allows a variety of plots
C***  CARDS line syntax:
C***    PLOT ACCELEM [LEADIONS] [CONT|LINE] [RELATIVE|GEFF] [X=key]
C***       LEADIONS: plot leading ions instead of total elements
C***                 (can be restricted to a certain minimum contrib) 
C***       CONT:     restrict to continuum contributions
C***       LINE:     restrict to line contributions
C***       RELATIVE: plot contributions relative to total ARAD
C***                    (respectively to total ALINE, ACONT)
C***       GEFF:     normalize to GEFF instead of GGRAV
C***       X=key:    defines the x-axis scale with key being one 
C***                                             of the following:
C***         DEPTHINDEX|VELOCITY|TAUROSS|RADIUS|ENTOT|RHO|L|R|V|TAU
C***           default is DEPTHINDEX aka L
C***
C***********************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, NATOM, MAXION, JOBNUM, hPLOT

      REAL, DIMENSION(ND), INTENT(IN) :: AGRAV, AMECH, APRESS, ARAD,
     >                                   ACONT, ATHOM, VELO, RADIUS, 
     >                                   T, ENTOT, RNE, XMU, TAUROSS
      REAL, DIMENSION(NATOM, ND-1) :: ARADELEM, ACONTELEM
      REAL, DIMENSION(ND-1, NATOM, MAXION) :: ARADION, ACONTION
      REAL, INTENT(IN) :: WORKRATIO, RCON, ATMEAN, XMSTAR, RSTAR, TEFF, 
     >                    Rcritical
     
      CHARACTER(2), DIMENSION(NATOM) :: SYMBOL
      CHARACTER(10), DIMENSION(NATOM), INTENT(IN) :: ELEMENT
     
      LOGICAL, INTENT(IN) :: bFULLHYDROSTAT, bNoARAD

      CHARACTER(110) :: MODHEAD, HEADLINE
      CHARACTER(5) :: CROMNUM      
      CHARACTER(8) :: CENTER, CNORM, CYOFF, CINT, CCOL
      CHARACTER(10) :: IONELEM
      CHARACTER(40) :: XTEXT, CELEMSYM
      CHARACTER(30) :: CUROPT, NEXTPAR, XAXISMODE, CESSIZE, YLABEL
      CHARACTER PLOTOPT*(*)

      REAL, DIMENSION(NDMAX) :: X, X2, Y1, Y2, Y3, Y4, Y5, 
     >                          RI, VSCRATCH, ANORM, VMACH, RHO
      CHARACTER(20), DIMENSION(NDMAX) :: CLEGEND

      INTEGER :: I, L, NPAR, NDIN, Lcand, ISYM, NA, NAOnly, ION,
     >           LL, LSTART, LEND, ICOL, NCOL, IC, ICONTRIB, IERR
      REAL :: XMIN, XMAX, YMIN, YMAX, YMINVAL, YCUR, YSTEP, YL,
     >        GEDDL, XLSTAR, XLSTARS, RNEINT, Xcrit, XINT, NIONUSED,
     >        XABST, XTICK, XLENGTH, XICON, Rsonic, Vsonic, Xsonic,
     >        RELCONTRIBMIN, tempREAL

      LOGICAL :: bNormalizeToGEFF, bCONT, bLINE, bRELATIVE, bLEADIONS,
     >           bIonUsed

      INTEGER, EXTERNAL :: IDX

C***  Roman numbers for ionization stages
      INTEGER, PARAMETER :: NROMNUM = 27
      CHARACTER(5), DIMENSION(NROMNUM) :: ROMNUM = (/
     >    'I    ', 'II   ', 'III  ', 'IV   ', 'V    ', 
     >    'VI   ', 'VII  ', 'VIII ', 'IX   ', 'X    ',
     >    'XI   ', 'XII  ', 'XIII ', 'XIV  ', 'XV   ',
     >    'XVI  ', 'XVII ', 'XVIII', 'XIX  ', 'XX   ',
     >    'XXI  ', 'XXII ', 'XXIII', 'XXIV ', 'XXV  ',
     >    'XXVI ', 'XXVII'  /)

C***  Plot symbols and sizes for detailed stage curves     
      INTEGER, PARAMETER :: MAXCOL = 6
      INTEGER, DIMENSION(MAXCOL), PARAMETER :: 
     >    ISYMCOLR = (/  8, 7, 5, 3, 2, 1 /)
      INTEGER, PARAMETER :: MAXSYM = 20
      INTEGER, DIMENSION(MAXSYM), PARAMETER ::
     >    NPLOTSYM = (/ 21, 22, 23, 24, 25,
     >                  29, 28,  8, 32, 26,
     >                  27, 14, 18, 30, 33,
     >                  15, 17,  4, 29, 16 /)
      LOGICAL, DIMENSION(MAXSYM), PARAMETER ::
     >    BPERLSYM = (/ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.,
     >                 .FALSE., .TRUE., .FALSE., .TRUE., .FALSE.,
     >                 .FALSE.,.FALSE.,.FALSE.,.FALSE., .FALSE.,
     >                 .FALSE.,.FALSE.,.FALSE.,.FALSE., .FALSE. /)
      REAL, DIMENSION(MAXSYM), PARAMETER ::
     >    SPLOTSYM = (/ 0.5 , 0.5 , 0.5 ,  0.5 , 0.5 ,
     >                  0.3 , 0.5 ,-0.25,  0.5 , 0.3 ,
     >                  0.3 , 0.15, 0.2 ,  0.3 , 0.3 ,
     >                 -0.3 ,-0.3 ,-0.3 , -0.3 , 0.3  /)

     
      
      !Physical constants
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)     
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (bNoARAD) THEN
         WRITE (hCPR,'(A)') 'PLOTACCELEM: Radiative acceleration ' //
     >      'not yet calculated ***'
         WRITE (hCPR,'(3A)') 'PLOTACCELEM: ', PLOTOPT(:IDX(PLOTOPT)),
     >      ' SKIPPED'
         RETURN
      ENDIF      
      
      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTACCELEM: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTACCELEM: DIMENSION INSUFFICIENT'
         WRITE (hCPR,'(A)') 'PLOTACCELEM: PLOT SKIPPED'
         RETURN
      ENDIF

      CENTER = '\CENTER\'

C***  Decode possible options for x-axis mode
      XAXISMODE = 'DEPTHINDEX'
      bNormalizeToGEFF = .FALSE.
      bCONT = .FALSE.
      bLINE = .FALSE.
      bRELATIVE = .FALSE.
      bLEADIONS = .FALSE.
      RELCONTRIBMIN = 0.005   !minimum relative contribution (currently 0.5%)      
      NAOnly = 0
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
C              IF (MINVAL(GEFFL) > 0.) THEN
                bNormalizeToGEFF = .TRUE.
C              ELSE
C                WRITE (hCPR,'(A)') 'PLOTACC: DATA FOR GEFF NOT',
C     >                     ' AVAILABLE - GGRAV IS USED INSTEAD'
C              ENDIF
            CASE ('CONT')
C***          This option switches plot to continuum contributions only            
              bCONT = .TRUE.
              bLINE = .FALSE.
            CASE ('LINE')
C***          This option switches plot to line contributions only            
              bLINE = .TRUE.
              bCONT = .FALSE.
            CASE ('LEADIONS')
C***          Plot all ions from all elements 
C***          contributing more than a defined amount
              bLEADIONS = .TRUE.
              IF (NPAR >= (I+1)) THEN
C***            optional parameter: minimum relative contribution for plotting
                CALL SARGV (PLOTOPT, I+1, NEXTPAR)
                READ (NEXTPAR, '(F15.0)', IOSTAT=IERR) tempREAL
                IF (IERR == 0) THEN
                  RELCONTRIBMIN = tempREAL
                ENDIF                                  
              ENDIF 
            CASE ('ION', 'IONS', 'ELEMENT')
              IF (NPAR >= (I+1)) THEN
                CALL SARGV (PLOTOPT, I+1, IONELEM)
                DO NA=1, NATOM
C***              Allow SYMBOL or full element name as identifier                
                  IF (IONELEM == ELEMENT(NA)) IONELEM = SYMBOL(NA)
                  IF (IONELEM == SYMBOL(NA)) THEN
                    NAOnly = NA
c                    WRITE (0,*) 'ELEMENT FOUND: ', ELEMENT(NAOnly)
                    EXIT
                  ENDIF
                ENDDO
C***            Skip this plot if no element is given or
C***            an unknown or unused element was specified
                IF (NAOnly == 0) THEN
                  WRITE (hCPR,'(A)') 'PLOTACCELEM: Invalid ' //
     >                                 'element specified ***'
                  WRITE (hCPR,'(3A)') 'PLOTACCELEM: ', 
     >                        PLOTOPT(:IDX(PLOTOPT)), ' SKIPPED'
                  RETURN                  
                ENDIF
              ENDIF
            CASE ('RELATIVE')
C***          This option switches plot to relative contributions
              bRELATIVE = .TRUE.
              bNormalizeToGEFF = .FALSE.
          ENDSELECT
        ENDDO
      ELSE
C***    Fallback for unique labeling inside steal.plot / wruniq.plot
        PLOTOPT = PLOTOPT(:IDX(PLOTOPT)) // ' L'
      ENDIF
C      IF (NPAR .GT. 2) CALL SARGV (PLOTOPT, 3, XAXISMODE)

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

      ANORM(ND) = 0.
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )      
        IF (bRELATIVE) THEN
          IF (bCONT) THEN
            ANORM(L) = ACONT(L)
          ELSEIF (bLINE) THEN
            ANORM(L) = ARAD(L)-ACONT(L)
          ELSE
            ANORM(L) = ARAD(L)
          ENDIF
          CNORM = ''
        ELSEIF (bNormalizeToGEFF) THEN
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
          CNORM = '/g&Teff&M'
        ELSE
          ANORM(L) = AGRAV(L)
          CNORM = '/g'
        ENDIF          
      ENDDO

      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        VMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !v_mach in km/s
        VSCRATCH(L) = VELO(L) - VMACH(L)        
        IF ((VSCRATCH(L) < 0) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,VMACH,RADIUS,ND)
      ENDIF
      
      
      IF (XAXISMODE == 'DEPTHINDEX' .OR. XAXISMODE == 'L') THEN 
C***     X-Axis = Depth Index
         XTEXT = CENTER//'DEPTH INDEX L'
         XMIN = 0.
         XMAX = FLOAT(ND)
         XTICK = 5.
         XABST = 10. 
C         DO L=1, ND-2
         DO L=1, ND-1
            X(L) = FLOAT(L) + 0.5
         ENDDO
         IF (Lcand > 1) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND)         
         ENDIF
      ELSEIF (XAXISMODE == 'VELOCITY' .OR. XAXISMODE == 'V') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v' // CHAR(92) // '8'
         XMIN = 0.
         XMAX = 1.
         XTICK = 0.1
         XABST = 0.5 
C         DO L=1, ND-2
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), VELO, RADIUS, ND)
            X(L) = XINT/VELO(1)
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
      ELSEIF (XAXISMODE == 'TAUROSS' .OR. XAXISMODE == 'TAU') THEN
         XTEXT = CENTER//'log #t#&TRoss&M'
         XMIN = LOG10(0.4 * TAUROSS(2))
         XMAX = LOG10( TAUROSS(ND) )
         XTICK = 0.1
         XABST = 1.0 
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), TAUROSS, RADIUS, ND)
           X(L) = LOG10 (XINT)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF
c         IF (Ltcand > 0 .AND. Ltcand < ND) THEN
c           CALL SPLINPOX(Xts,Rts,X,RADIUS,ND-1)         
c         ENDIF
      ELSEIF (XAXISMODE == 'ENTOT') THEN
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
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF      
      ELSEIF (XAXISMODE == 'RHO') THEN
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
      ELSE
         WRITE (hCPR,*) 'PLOTACC: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))
         RETURN
      ENDIF


      CALL JSYMSET ('G2','TRANSFER')

      HEADLINE = 'ACCELEM:'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM
C    9 FORMAT (19X,1HJ,I3)

      WRITE (hPLOT, '(A)') '*NEXTPLOT: ' // PLOTOPT
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\SET_NSETMAX 1000'            
      WRITE (hPLOT, '(A)') '\INSTRUCTION EXPORT'            
C*** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                            
C***  For an easier re-use of the plots in papers and talks,
C***  we store the different element styles in WRPlot array variables

      IF (NAOnly == 0 .AND. .NOT. bLEADIONS) THEN
C***    for the contributions of the different elements      
        WRITE (hPLOT, '(A)') ''
        DO NA=1, NATOM
          WRITE(UNIT=CINT,FMT='(I2)') NA
          WRITE (hPLOT, '(5A)') '\VAR elemname[', TRIM(ADJUSTL(CINT)),
     >                                '] = "', SYMBOL(NA), '"'
        ENDDO 
        WRITE (hPLOT, '(A)') ''
        DO NA=1, NATOM
          WRITE(UNIT=CINT,FMT='(I2)') NA
          WRITE (hPLOT, '(3A)') '\VAR elemcol[', TRIM(ADJUSTL(CINT)),
     >                                '] = 8'
        ENDDO 
        WRITE (hPLOT, '(A)') ''
        DO NA=1, NATOM
          WRITE(UNIT=CINT,FMT='(I2)') NA
          WRITE(UNIT=CELEMSYM,FMT='(A,I2)') 'SYMBOL=', NPLOTSYM(NA)
          WRITE(UNIT=CESSIZE,FMT='(A,F4.2)') 'SIZE=', SPLOTSYM(NA)
          IF (BPERLSYM(NA)) THEN
            WRITE (hPLOT, '(7A)') '\VAR elemsym[', TRIM(ADJUSTL(CINT)),
     >               '] = "', TRIM(CESSIZE), ' ', TRIM(CELEMSYM), '"'
            WRITE (hPLOT, '(3A)') '\VAR elemlin[', TRIM(ADJUSTL(CINT)),
     >               '] = "SYMBOL=0"'
          ELSE
            WRITE (hPLOT, '(7A)') '\VAR elemsym[', TRIM(ADJUSTL(CINT)),
     >               '] = "', TRIM(CESSIZE), ' ', TRIM(CELEMSYM), '"'
            WRITE (hPLOT, '(3A)') '\VAR elemlin[', TRIM(ADJUSTL(CINT)),
     >               '] = "SYMBOL=5"'
          ENDIF
        ENDDO

      ELSE 
        IF (bLEADIONS) THEN
C***    for leading ionization stages of all elements
          ICONTRIB = 0
          NCOL = 1
          ICOL = ISYMCOLR(NCOL)
          WRITE(UNIT=CCOL,FMT='(I2)') ICOL
          DO NA=1, NATOM
            DO ION=1, MAXION

              IF (bCONT) THEN
                bIonUsed = (MAXVAL(ACONTION(1:ND-1,NA,ION)) > 0.)
              ELSEIF (bLINE) THEN
                bIonUsed = (MAXVAL(ARADION(1:ND-1,NA,ION)-
     >                             ACONTION(1:ND-1,NA,ION)) > 0.)
              ELSE
                bIonUsed = (MAXVAL(ARADION(1:ND-1,NA,ION)) > 0.)
              ENDIF
              IF (bIonUsed) THEN
            
                DO L=1, ND-1
                  IF (bCONT) THEN
                    YL = (ACONTION(L,NA,ION))/ANORM(L)
                  ELSEIF (bLINE) THEN
                    YL = (ARADION(L,NA,ION)-ACONTION(L,NA,ION))/ANORM(L)
                  ELSE
                    YL = (ARADION(L,NA,ION))/ANORM(L)
                  ENDIF
                  IF (YL > RELCONTRIBMIN) THEN
              
                    IF (ION <= NROMNUM) THEN
C***                  Get roman number for current ionization stage            
                      CROMNUM = ROMNUM(ION)
                    ELSE 
C***                  Failsafe: use arabic number
                      WRITE(UNIT=CROMNUM,FMT='(I5)') ION
                    ENDIF

C***                symbol choice, legend and line plot              
                    ICONTRIB = ICONTRIB + 1
                    WRITE(UNIT=CINT,FMT='(I2)') ICONTRIB
                    ISYM = ISYM + 1
                    IF (ISYM > MAXSYM) THEN
C***                  change color if we run out of different symbol styles                  
                      ISYM = 1
                      NCOL = NCOL + 1
                      ICOL = ISYMCOLR(NCOL)
                      WRITE(UNIT=CCOL,FMT='(I2)') ICOL
                    ENDIF  
                    CLEGEND(ICONTRIB) = TRIM(SYMBOL(NA)) // '\,'
     >                            // TRIM(ADJUSTL(CROMNUM))
                    WRITE (hPLOT, '(5A)')
     >                '\VAR ionname[', TRIM(ADJUSTL(CINT)), '] = "',
     >                     TRIM(ADJUSTL(CLEGEND(ICONTRIB))), '"'
                    WRITE (hPLOT, '(3A)') 
     >                '\VAR ioncol[', TRIM(ADJUSTL(CINT)), '] = '
     >                                         // TRIM(ADJUSTL(CCOL))          
     
                    EXIT !exit ND loop as soon as we know ion is contributing
                  ENDIF
                ENDDO
                
              ENDIF
            
            ENDDO
          ENDDO
        
          NIONUSED = ICONTRIB
        ELSE            
C***    for ionization stages of a single element (specified in NAOnly)      
          WRITE (hPLOT, '(A)') ''
          ISYM = 0
          DO ION=1, MAXION
C***        Plot only Ions which contribute        
            IF (bCONT) THEN
              bIonUsed = (MAXVAL(ACONTION(1:ND-1,NAOnly,ION)) > 0.)
            ELSEIF (bLINE) THEN
              bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)
     >                           -ACONTION(1:ND-1,NAOnly,ION)) > 0.)
            ELSE
              bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)) > 0.)
            ENDIF
            IF (bIonUsed) THEN
              ISYM = ISYM + 1
              WRITE(UNIT=CINT,FMT='(I2)') ISYM
              IF (ION <= NROMNUM) THEN
C***            Get roman number for current ionization stage            
                CROMNUM = ROMNUM(ION)
              ELSE 
C***            Failsafe: use arabic number
                WRITE(UNIT=CROMNUM,FMT='(I5)') ION
              ENDIF
              WRITE (hPLOT, '(5A)') '\VAR ionname[',TRIM(ADJUSTL(CINT)),
     >                    '] = "', TRIM(ADJUSTL(CROMNUM)), '"'
              WRITE (hPLOT, '(3A)') '\VAR ioncol[', TRIM(ADJUSTL(CINT)),
     >                                '] = 8'
            ENDIF
          ENDDO 
          NIONUSED = ISYM
        ENDIF
        
C***    write line styles into WRplot array variables        
        WRITE (hPLOT, '(A)') ''
        ISYM = 0
        DO IC=1, NIONUSED
          ISYM = ISYM + 1
          IF (ISYM > MAXSYM) ISYM = 1
          WRITE(UNIT=CINT,FMT='(I2)') IC
          WRITE(UNIT=CELEMSYM,FMT='(A,I2)') 'SYMBOL=', NPLOTSYM(ISYM)
          WRITE(UNIT=CESSIZE,FMT='(A,F4.2)') 'SIZE=', SPLOTSYM(ISYM)
          IF (BPERLSYM(ISYM)) THEN
            WRITE (hPLOT, '(7A)') '\VAR ionsym[', TRIM(ADJUSTL(CINT)),
     >               '] = "', TRIM(CESSIZE), ' ', TRIM(CELEMSYM), '"'
            WRITE (hPLOT, '(3A)') '\VAR ionlin[', TRIM(ADJUSTL(CINT)),
     >               '] = "SYMBOL=0"'
          ELSE
            WRITE (hPLOT, '(7A)') '\VAR ionsym[', TRIM(ADJUSTL(CINT)),
     >               '] = "', TRIM(CESSIZE), ' ', TRIM(CELEMSYM), '"'
            WRITE (hPLOT, '(3A)') '\VAR ionlin[', TRIM(ADJUSTL(CINT)),
     >               '] = "SYMBOL=5"'
          ENDIF
        ENDDO
      
      ENDIF

C*** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                            
      WRITE (hPLOT, '(A)') '\INSTRUCTION NO-EXPORT'      
      
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 5 1.0 0.58 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 7 0.0 0.50 0.0'
      IF (bCONT) THEN
        WRITE (hPLOT, '(A)') '\DEFINECOLOR 8 0.6 0.60 0.0'
      ELSEIF (bLINE) THEN
        WRITE (hPLOT, '(A)') '\DEFINECOLOR 8 0.3 0.60 0.5'
      ELSE
        WRITE (hPLOT, '(A)') '\DEFINECOLOR 8 0.5 0.20 0.2'
      ENDIF
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 9 0.6 0.60 0.6'
      WRITE (hPLOT, '(A)') '\PEN=1'     !needed to get rid of definecolor pen bug 
      
      WRITE (hPLOT, '(A)') '\COLOR=3'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0 XMAX 0 0 0'

      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON >= RI(ND) .AND. RCON <= RI(1)) THEN
        CALL SPLINPOX(XICON, RCON, X, RI, ND)
        IF (XICON >= 0.005) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',XICON,' YMIN ',XICON,' YMAX 0. 0. SYMBOL=10'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',XICON,' YMAX 0. -0.2 0.2 -90  R&Tcon'
        ENDIF
      ENDIF
      IF (Lcand > 1) THEN
        IF (Xsonic >= 0.005) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',Xsonic,' YMIN ',Xsonic,
     >                    ' YMAX 0. 0. SYMBOL=20 SIZE=0.05'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',Xsonic,' YMAX 0. -0.2 0.2 -90  R&Tsonic'
        ENDIF
      ENDIF
      IF (Rcritical >= RADIUS(ND) .AND. Rcritical <= RADIUS(1)) THEN
        CALL SPLINPOX(Xcrit, Rcritical, X, RADIUS, ND)
        IF (Xcrit >= 0.005) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',Xcrit,' YMIN ',Xcrit,' YMAX 0. 0. SYMBOL=9'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',Xcrit,' YMAX 0. -0.2 0.2 -90  R&Tcrit'
        ENDIF
      ENDIF
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'
      WRITE (hPLOT, '(A)') '\COLOR=1'

      WRITE (hPLOT, '(A,F6.3)') 
     >   'KASDEF LUN XMIN YMIN 1. 1. .5 ' //
     >    '&EWORK Ratio Rad.+Gaspress. / Mech.+Grav. =&N', WORKRATIO  
      WRITE (hPLOT, '(A)') '\INBOX'

      IF ((NAOnly == 0) .AND. .NOT. bRELATIVE) THEN
        XLENGTH = ( XMAX - XMIN ) / 10.      
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -1. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN ', XLENGTH, ' YMAX 14.5 M-1.'
     >   // ' 0.3 (' // CNORM(2:IDX(CNORM)) // '+a&Tmech&M)' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -1.5 SYMBOL=9'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN ',XLENGTH,' YMAX 14.5 M-1.5 0.3 a&Trad&M' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=4'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -2. SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN ',XLENGTH,' YMAX 14.5 M-2. 0.3 a&Tmech&M' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=7'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -2.5 SYMBOL=9'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN ',XLENGTH,' YMAX 14.5 M-2.5 0.3 a&Tcont&M' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -3 SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN', XLENGTH,
     >   ' YMAX 14.5 M-3.0 0.3 (a&Trad&M+a&Tpress&M)' // CNORM

        WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=5'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LINUN 0. YMAX ',XLENGTH,' YMAX 14. -3.5 SYMBOL=9 SIZE=0.05'
        WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     >   'LUN ',XLENGTH,' YMAX 14.5 M-3.5 0.3 a&Tpress&M' // CNORM

        WRITE (hPLOT, '(A)') '\COLOR=6'
        WRITE (hPLOT, '(A,F4.1,A)') 
     >   '\LINUN 0. YMAX ',XLENGTH,' YMAX 14. -4.0  SYMBOL=5'
        WRITE (hPLOT, '(A,F4.1,A)') '\LUN ', XLENGTH,
     >   ' YMAX 14.5 M-4.0 0.3 a&Tthom&M' // CNORM
      ENDIF
      
      
      
      WRITE (hPLOT, '(A)') '\VAR xleglin = 1.'
      WRITE (hPLOT, '(A)') '\VAR xleglins = 0.5'
      WRITE (hPLOT, '(A)') '\VAR yleglin = 0.'
      WRITE (hPLOT, '(A)') '\VAR xlegol = 0.25'
      WRITE (hPLOT, '(A)') '\VAR xlegols = 0.5'
      WRITE (hPLOT, '(A)') '\VAR xlegotval = 1.5'
      WRITE (hPLOT, '(A)') '\EXPR xlegot = "L+" // $xlegotval'
      WRITE (hPLOT, '(A)') '\VAR xlegcolstep = 2.8'
      WRITE (hPLOT, '(A)') '\VAR tsizeleg = 0.3'
      IF (NAOnly == 0 .AND. .NOT. bLEADIONS) THEN
        WRITE (hPLOT, '(A)') '\VAR ylegstep = 0.5'
        WRITE (hPLOT, '(A)') '\VAR ylegcur = 0.7'
        WRITE(UNIT=CINT,FMT='(I2)') NATOM
        WRITE (hPLOT, '(2A)') '\DO ENDLEG ileg=1,', TRIM(ADJUSTL(CINT))
        WRITE (hPLOT, '(A)') '  \COLOR=$elemcol[$ileg]'
        WRITE (hPLOT, '(A)') '  \VAR curlin = $elemlin[$ileg]'
        WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
        WRITE (hPLOT, '(A)') '  \IF $curlin .EQ. ""'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglin $yleglin '
     >                         // ' $xlegol $ylegcur $elemsym[$ileg]'
        WRITE (hPLOT, '(A)') '  \ELSE'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglin $yleglin '
     >                         // ' $xlegol $ylegcur $elemlin[$ileg]'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglins $yleglin '
     >                         // ' $xlegols $ylegcur $elemsym[$ileg]'
        WRITE (hPLOT, '(A)') '  \ENDIF'
        WRITE (hPLOT, '(A)') '  \LUN SAME SAME $xlegot $ylegcurlab '
     >                             //  ' $tsizeleg $elemname[$ileg]'
        WRITE (hPLOT, '(A)') '  \CALC ylegcur = $ylegcur + $ylegstep'
        WRITE (hPLOT, '(A)') '\LABEL ENDLEG'
        IF (bCONT) THEN
          WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
          WRITE (hPLOT, '(A)') '\LUN XMAX YMIN R1.5 U0 '
     >                             //  ' $tsizeleg &FCONT&f'
        ELSEIF (bLINE) THEN
          WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
          WRITE (hPLOT, '(A)') '\LUN XMAX YMIN R1.5 U0 '
     >                             //  ' $tsizeleg &FLINE&f'
        ENDIF
      ELSEIF (bLEADIONS) THEN

        WRITE (hPLOT, '(A)') '\VAR ylegstep = 0.5'
        WRITE (hPLOT, '(A)') '\VAR ylegcurstart = 0.5'
        WRITE (hPLOT, '(A)') '\VAR ylegcur = $ylegcurstart'
        WRITE (hPLOT, '(A)') '\VAR xlegcol = 1' 
        WRITE(UNIT=CINT,FMT='(I2)') ICONTRIB
        WRITE (hPLOT, '(2A)') '\DO ENDLEG ileg=1,', TRIM(ADJUSTL(CINT))
        WRITE (hPLOT, '(A)') '  \COLOR=$ioncol[$ileg]'
        WRITE (hPLOT, '(A)') '  \VAR curlin = $ionlin[$ileg]'
        WRITE (hPLOT, '(A)') '  \IF $ileg .GT. 30'
        WRITE (hPLOT, '(A)') '    \CALC ilegmod = $ileg - 30'
        WRITE (hPLOT, '(A)') '    \VAR xlegcol = 2'
        WRITE (hPLOT, '(A)') '    \LABEL morecols'
        WRITE (hPLOT, '(A)') '    \IF $ilegmod .GT. 30'
        WRITE (hPLOT, '(A)') '      \CALC ilegmod = $ilegmod - 30'
        WRITE (hPLOT, '(A)') '      \CALC xlegcol = $xlegcol + 1'
        WRITE (hPLOT, '(A)') '      \GOTO morecols'
        WRITE (hPLOT, '(A)') '    \ENDIF'       
        WRITE (hPLOT, '(A)') '    \IF $ilegmod .EQ. 1'
        WRITE (hPLOT, '(A)') '      \VAR ylegcur = $ylegcurstart'       
        WRITE (hPLOT, '(A)') '    \ENDIF'       
        WRITE (hPLOT, '(A)') '  \ENDIF'
        WRITE (hPLOT, '(A)') '  \CALC xlegolcur = $xlegol + ($xlegcol - 1) * $xlegcolstep' 
        WRITE (hPLOT, '(A)') '  \CALC xlegolscur = $xlegols + ($xlegcol - 1) * $xlegcolstep' 
        WRITE (hPLOT, '(A)') '  \CALC xlegotcur = $xlegotval + ($xlegcol - 1) * $xlegcolstep' 
        WRITE (hPLOT, '(A)') '  \EXPR xlegotcur = "L" // $xlegotcur'
        WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
        WRITE (hPLOT, '(A)') '  \IF $curlin .EQ. ""'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglin $yleglin '
     >                         // ' $xlegolcur $ylegcur $ionsym[$ileg]'
        WRITE (hPLOT, '(A)') '  \ELSE'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglin $yleglin '
     >                         // ' $xlegolcur $ylegcur $ionlin[$ileg]'
        WRITE (hPLOT, '(A)') '    \LINREL XMAX YMIN $xleglins $yleglin '
     >                         // ' $xlegolscur $ylegcur $ionsym[$ileg]'
        WRITE (hPLOT, '(A)') '  \ENDIF'
        WRITE (hPLOT, '(A)') '  \LUN SAME SAME $xlegotcur $ylegcurlab '
     >                             //  ' $tsizeleg $ionname[$ileg]'
        WRITE (hPLOT, '(A)') '  \CALC ylegcur = $ylegcur + $ylegstep'
        WRITE (hPLOT, '(A)') '\LABEL ENDLEG'
        IF (bCONT) THEN
          WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
          WRITE (hPLOT, '(A)') '\LUN XMAX YMIN R1.5 U0 '
     >                             //  ' $tsizeleg &FCONT&f'
        ELSEIF (bLINE) THEN
          WRITE (hPLOT, '(A)') '  \EXPR ylegcurlab = "M" // $ylegcur'
          WRITE (hPLOT, '(A)') '\LUN XMAX YMIN R1.5 U0 '
     >                             //  ' $tsizeleg &FLINE&f'
        ENDIF
      
      ELSE  
C***    Special Label for selected element:
        WRITE (hPLOT, '(A)') '\COLOR=8'
        YSTEP = 0.5
        YCUR = -0.7
        WRITE(UNIT=CYOFF,FMT='(F5.2)') YCUR
        WRITE (hPLOT, '(A)') 
     >      '\LUN XMAX YMAX 0.5 ' // TRIM(ADJUSTL(CYOFF))  
     >                 // ' 0.5 ' // ELEMENT(NAOnly)         
c        ISYM = 20
        ISYM = 0
        DO ION=1, MAXION
C***      Plot only Ions which contribute        
          IF (bCONT) THEN
            bIonUsed = (MAXVAL(ACONTION(1:ND-1,NAOnly,ION)) > 0.)
          ELSEIF (bLINE) THEN
            bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)
     >                         -ACONTION(1:ND-1,NAOnly,ION)) > 0.)
          ELSE
            bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)) > 0.)
          ENDIF
          IF (bIonUsed) THEN
            ISYM = ISYM + 1
c            IF (ISYM == 31) ISYM = ISYM + 1
            YCUR = YCUR - YSTEP
            WRITE(UNIT=CYOFF,FMT='(F5.2)') YCUR
            WRITE(UNIT=CELEMSYM,FMT='(A,I2)') 'SYMBOL=', NPLOTSYM(ISYM)
            WRITE(UNIT=CESSIZE,FMT='(A,F4.2)') 'SIZE=', SPLOTSYM(ISYM)
            IF (BPERLSYM(ISYM)) THEN
              WRITE (hPLOT, '(A)') 
     >          '\LINREL XMAX YMAX 1. 0. 0.5 ' // TRIM(ADJUSTL(CYOFF)) 
     >          // ' ' // TRIM(CELEMSYM) // ' ' // TRIM(CESSIZE)
            ELSE
              WRITE (hPLOT, '(A)') 
     >         '\LINREL XMAX YMAX 1. 0. 0.5 ' // TRIM(ADJUSTL(CYOFF)) 
     >          // ' SYMBOL=5  SIZE=0.5'
              WRITE (hPLOT, '(A)') 
     >         '\LINREL XMAX YMAX 0.5 0. 0.75 ' // TRIM(ADJUSTL(CYOFF)) 
     >          // ' ' // TRIM(CELEMSYM) // ' ' // TRIM(CESSIZE)
            ENDIF
            IF (ION <= NROMNUM) THEN
C***          Get roman number for current ionization stage            
              CROMNUM = ROMNUM(ION)
            ELSE 
C***          Failsafe: use arabic number
              WRITE(UNIT=CROMNUM,FMT='(I5)') ION
            ENDIF
            WRITE (hPLOT, '(A)') 
     >        '\LUN SAME SAME L+1.75 M' // TRIM(ADJUSTL(CYOFF)) 
     >            // ' 0.3 ' // TRIM(ADJUSTL(CROMNUM))
          ENDIF
        ENDDO        
      ENDIF
      
      WRITE (hPLOT, '(A)') '\COLOR=1'

C***  HTOT: Conversion into radiation temperatures
C***  Note: HTOT may not be calculated if STEAL is used for OUTPUT ONLY
      YMIN = 0.
      YMAX = 0.
cc      YMIN = -1.
cc      YMAX =  1.
      DO L=1, ND-1
         Y1(L) = (ANORM(L) + AMECH(L))/ANORM(L)
         Y2(L) = ARAD(L)/ANORM(L)
         Y3(L) = AMECH(L)/ANORM(L)
         Y4(L) = ATHOM(L)/ANORM(L)
         Y5(L) = APRESS(L)/ANORM(L)
cc         IF ((Y1(L) .GT. 0.) .AND. (Y2(L) .GT. 0.)) THEN
cc            YMAX = ALOG10(MAX(10**YMAX, 1.05*Y1(L), 1.05*Y2(L)))
cc            YMIN = ALOG10(MIN(10**YMIN, 0.95*Y1(L), 0.95*Y2(L)))
cc         ENDIF
         IF (.NOT. bRELATIVE) THEN
           IF (Y1(L) > .0) Y1(L) = LOG10(Y1(L))
           IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
           IF (Y3(L) > .0) Y3(L) = LOG10(Y3(L))
           IF (Y4(L) > .0) Y4(L) = LOG10(Y4(L))
           IF (Y5(L) > .0) Y5(L) = LOG10(Y5(L))
         ELSE 
           Y1(L) = 1.
         ENDIF
      ENDDO

C***  Restrict blue curve (AMECH) to inbox area
      YMINVAL = MIN(MINVAL(Y2), MINVAL(Y4))
      YMINVAL = MIN(YMINVAL, MINVAL(Y5))
      DO L = 1, ND-1
        NDIN = L
        IF (Y3(L) < YMINVAL) EXIT
      ENDDO
      
      IF (bRELATIVE) THEN
        IF (bCONT) THEN
          YLABEL = 'a&Tc,elem&M/a&Tcont&M'
        ELSEIF (bLINE) THEN
          YLABEL = 'a&Tc,elem&M/a&Tline&M'
        ELSE
          YLABEL = 'a&Trad,elem&M/a&Trad&M'
        ENDIF
      ELSE
        YLABEL = 'log&T10&M(a/g)'
      ENDIF

      CALL PLOTANF (hPLOT,HEADLINE, '&E'//HEADLINE,
     $        XTEXT, CENTER//TRIM(YLABEL),
     >        0., XMIN, XMAX, XTICK, XABST, 0.,
     >        0., YMIN, YMAX, .1,  0.2, 0.,
     $        X, Y1, ND-1, 5)

C***  In the standard version, the absolute terms for 
C***  ARAD, AMECH, ATHOM and APRESS are all displayed. 
C***  For the RELATIVE option, only the normalized ATHOM is shown
C***  together with the elemental contributions calculated further down.
      IF (.NOT. bRELATIVE) THEN
C***    ARAD      
        CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2 SYMBOL=9 SIZE=0.1') 
C***    AMECH (only for values larger than YMIN)
        CALL PLOTCONS (hPLOT,X,Y3,NDIN,'COLOR=4') 
      ENDIF
C***  ATHOM
      IF (NAOnly == 0 .OR. (.NOT. bRELATIVE)) THEN
        CALL PLOTCONS (hPLOT,X,Y4,ND-1,'COLOR=6') 
      ENDIF
      IF (.NOT. bRELATIVE) THEN
C***    APRESS
        CALL PLOTCONS (hPLOT,X,Y5,ND-1,'COLOR=5 SYMBOL=9 SIZE=0.05') 
      ENDIF

      IF (bRELATIVE) THEN
C***    For RELATIVE option, build total checksum of ARAD (or ACONT)
C***    The local Y3 array is re-used for this case
C***    and initialized with the normalized ATHOM, stored in Y4
        DO L=1, ND-1
          IF (NAOnly == 0) THEN
            Y3(L) = Y4(L)
          ELSE 
C***        If we plot the ions of a single element only, 
C***        there is no Thomson contribution in the checksum!
            Y3(L) = 0.
          ENDIF
        ENDDO
      ENDIF
c      ISYM = 20
      ISYM = 0
      IF (NAOnly > 0) THEN
        IF (bRELATIVE) THEN
C***      For relative contributions, normalize with the total
C***      contribution of the selected element
          DO L=1, ND
            IF (bCONT) THEN
              ANORM(L) = ACONTELEM(NAOnly,L)
            ELSEIF (bLINE) THEN
              ANORM(L) = ARADELEM(NAOnly,L)-ACONTELEM(NAOnly,L)
            ELSE 
              ANORM(L) = ARADELEM(NAOnly,L)
            ENDIF
          ENDDO
        ENDIF
        DO ION=1, MAXION
C***      Plot only Ions which contribute        
          IF (bCONT) THEN
            bIonUsed = (MAXVAL(ACONTION(1:ND-1,NAOnly,ION)) > 0.)
          ELSEIF (bLINE) THEN
            bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)
     >                         -ACONTION(1:ND-1,NAOnly,ION)) > 0.)
          ELSE
            bIonUsed = (MAXVAL(ARADION(1:ND-1,NAOnly,ION)) > 0.)
          ENDIF
          IF (bIonUsed) THEN
            ISYM = ISYM + 1
            WRITE(UNIT=CINT,FMT='(I2)') ISYM
c            IF (ISYM == 31) ISYM = ISYM + 1
            DO L=1, ND-1        
              IF (bCONT) THEN
                Y2(L) = ACONTION(L,NAOnly,ION)/ANORM(L)
              ELSEIF (bLINE) THEN
                Y2(L) = (ARADION(L,NAOnly,ION)-ACONTION(L,NAOnly,ION))
     >                                   /ANORM(L)
              ELSE
                Y2(L) = ARADION(L,NAOnly,ION)/ANORM(L)
              ENDIF
              IF (.NOT. bRELATIVE) THEN
C***            Show contributions larger than 10^(-10)
                IF (Y2(L) < .0) THEN
C***              Anormal (inward) acceleration is indicated by values 
C***              below the standard minimum level
                  Y2(L) = -11.
                ELSE
                  Y2(L) = LOG10(MAX(1.E-10,Y2(L)))
                ENDIF
              ELSE 
C***            RELATIVE option: Build total checksum of ARADELEM (or ACONTELEM)
C***            by addding the contribution of the current ion to Y3
                Y3(L) = Y3(L) + Y2(L)
              ENDIF
            ENDDO
            IF (.NOT. BPERLSYM(ISYM)) THEN
              CELEMSYM = 'COLOR=$ioncol[' // TRIM(ADJUSTL(CINT)) // ']'
     >                // ' $ionlin[' // TRIM(ADJUSTL(CINT)) // ']'
              CALL PLOTCONS (hPLOT,X,Y2,ND-1, TRIM(CELEMSYM)) 
            ENDIF
            CELEMSYM = 'COLOR=$ioncol[' // TRIM(ADJUSTL(CINT)) // ']'
     >                // ' $ionsym[' // TRIM(ADJUSTL(CINT)) // ']'
c            WRITE(UNIT=CELEMSYM,FMT='(A,F4.2,A,I2)') 
c     >        'COLOR=8 SIZE=',SPLOTSYM(ISYM),' SYMBOL=', NPLOTSYM(ISYM)
            CALL PLOTCONS (hPLOT,X,Y2,ND-1, TRIM(CELEMSYM)) 
          ENDIF
        ENDDO
      ELSEIF (bLEADIONS) THEN
        ICONTRIB = 0.
        NCOL = 1
        ICOL = ISYMCOLR(NCOL)
        WRITE(UNIT=CCOL,FMT='(I2)') ICOL
        DO NA=1, NATOM
          DO ION=1, MAXION
C***        TODO: How to judge LINE and CONT options??
            
            IF (bCONT) THEN
              bIonUsed = (MAXVAL(ACONTION(1:ND-1,NA,ION)) > 0.)
            ELSEIF (bLINE) THEN
              bIonUsed = (MAXVAL(ARADION(1:ND-1,NA,ION)
     >                           -ACONTION(1:ND-1,NA,ION)) > 0.) 
            ELSE
              bIonUsed = (MAXVAL(ARADION(1:ND-1,NA,ION)) > 0.)
            ENDIF
            IF (.NOT. bIonUsed) CYCLE
            
C***      Plot only the range from the outermost to innerpost point where
C***      the contribution is larger than RELCONTRIBMIN            
            LSTART=0
            LEND=0
            DO L=1, ND-1
              IF (bCONT) THEN
                YL = (ACONTION(L,NA,ION))/ANORM(L)
              ELSEIF (bLINE) THEN
                YL = (ARADION(L,NA,ION)-ACONTION(L,NA,ION))/ANORM(L)
              ELSE
                YL = (ARADION(L,NA,ION))/ANORM(L)
              ENDIF
              IF (YL > RELCONTRIBMIN) THEN
                LSTART = L
                EXIT
              ENDIF
            ENDDO
            DO L=ND-1, 1, -1
              IF (bCONT) THEN
                YL = (ACONTION(L,NA,ION))/ANORM(L)
              ELSEIF (bLINE) THEN
                YL = (ARADION(L,NA,ION)-ACONTION(L,NA,ION))/ANORM(L)
              ELSE
                YL = (ARADION(L,NA,ION))/ANORM(L)
              ENDIF
              IF (YL > RELCONTRIBMIN) THEN
                LEND = L
                EXIT
              ENDIF
            ENDDO
            IF (LSTART > 0) THEN
              LL = 0
              DO L=LSTART, LEND  
                LL = LL+1
                X2(LL) = X(L)
                IF (bCONT) THEN
                  Y2(LL) = (ACONTION(L,NA,ION))/ANORM(L)
                ELSEIF (bLINE) THEN
                  Y2(LL) = (ARADION(L,NA,ION)-ACONTION(L,NA,ION))/ANORM(L)
                ELSE
                  Y2(LL) = (ARADION(L,NA,ION))/ANORM(L)
                ENDIF

                IF (.NOT. bRELATIVE) THEN
                  IF (Y2(LL) < .0) THEN
C***              Anormal (inward) acceleration is indicated by values 
C***              below the standard minimum level
                    Y2(LL) = -11.
                  ELSE
                    Y2(LL) = LOG10(MAX(1.E-10,Y2(LL)))
                  ENDIF
                ENDIF
                
              ENDDO              

C***          symbol choice, legend and line plot              
              ICONTRIB = ICONTRIB + 1
              ISYM = ISYM + 1
              IF (ISYM > MAXSYM) ISYM = 1
              
              WRITE(UNIT=CINT,FMT='(I2)') ICONTRIB
              IF (.NOT. BPERLSYM(ISYM)) THEN
                CELEMSYM = 'COLOR=$ioncol['// TRIM(ADJUSTL(CINT)) //']'
     >                // ' $ionlin[' // TRIM(ADJUSTL(CINT)) // ']'
                CALL PLOTCONS (hPLOT,X2,Y2,LL, TRIM(CELEMSYM)) 
              ENDIF
              CELEMSYM = 'COLOR=$ioncol[' // TRIM(ADJUSTL(CINT)) // ']'
     >                // ' $ionsym[' // TRIM(ADJUSTL(CINT)) // ']'
              
              CALL PLOTCONS (hPLOT,X2,Y2,LL, TRIM(CELEMSYM)) 
            ENDIF
            
          ENDDO
        ENDDO
      ELSE
        DO NA=1, NATOM
          ISYM = ISYM + 1
          WRITE(UNIT=CINT,FMT='(I2)') ISYM
c          IF (ISYM == 31) ISYM = ISYM + 1
          DO L=1, ND-1
            IF (bCONT) THEN
              Y2(L) = (ACONTELEM(NA,L))/ANORM(L)
            ELSEIF (bLINE) THEN
              Y2(L) = (ARADELEM(NA,L)-ACONTELEM(NA,L))/ANORM(L)
            ELSE
              Y2(L) = (ARADELEM(NA,L))/ANORM(L)
            ENDIF
            IF (.NOT. bRELATIVE) THEN
C***        Show contributions larger than 10^(-10)
              IF (Y2(L) < .0) THEN
C***          Anormal (inward) acceleration is indicated by values 
C***          below the standard minimum level
                Y2(L) = -11.
              ELSE
                Y2(L) = LOG10(MAX(1.E-10,Y2(L)))
              ENDIF
            ELSE 
C***          RELATIVE option: Build total checksum of ARAD (or ACONT)
C***          by addding the contribution of the current element to Y3
              Y3(L) = Y3(L) + Y2(L)
            ENDIF
          ENDDO
          IF (.NOT. BPERLSYM(ISYM)) THEN
              CELEMSYM = 'COLOR=$elemcol[' // TRIM(ADJUSTL(CINT)) // ']'
     >                // ' $elemlin[' // TRIM(ADJUSTL(CINT)) // ']'
c            WRITE(UNIT=CELEMSYM,FMT='(A)') 'COLOR=8 SIZE=0.5 SYMBOL=5'
            CALL PLOTCONS (hPLOT,X,Y2,ND-1, TRIM(CELEMSYM)) 
          ENDIF
          CELEMSYM = 'COLOR=$elemcol[' // TRIM(ADJUSTL(CINT)) // ']'
     >                // ' $elemsym[' // TRIM(ADJUSTL(CINT)) // ']'
c          WRITE(UNIT=CELEMSYM,FMT='(A,F4.2,A,I2)') 
c     >        'COLOR=8 SIZE=',SPLOTSYM(ISYM),' SYMBOL=', NPLOTSYM(ISYM)
          CALL PLOTCONS (hPLOT,X,Y2,ND-1, TRIM(CELEMSYM)) 
        ENDDO
      ENDIF
      IF (bRELATIVE .AND. .NOT. bLEADIONS) THEN
C***    Plot total checksum of contributions
C***    (in red for ARAD contribs, green for ACONT contribs)
        IF (bCONT) THEN
          CALL PLOTCONS (hPLOT,X,Y3,ND-1,'COLOR=7 SYMBOL=9 SIZE=0.1')
        ELSE
          CALL PLOTCONS (hPLOT,X,Y3,ND-1,'COLOR=2 SYMBOL=9 SIZE=0.1') 
        ENDIF
      ENDIF
      
      RETURN
      END
