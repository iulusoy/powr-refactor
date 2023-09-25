      SUBROUTINE PLOTFGSTRAT(PLOTOPT, RADIUS, RSTAR, VELO, ND, DOTM4P, 
     >                       AMACH, XMG, CKL, ALPHAL,
     >                       GEDDL, GAMMARAD, F_raw, G_raw,
     >                       ENTOT, TAUROSS,
     >                       MODHEAD, INDOUT, JOBNUM, KANAL, KEY)
C***********************************************************************
C***   plot precalculated (preguessing) values from HDSOLUTION
C***    creates prehydro.plot or fgparts.plot
C***   NOTE: RADIUS and VELO are in CM as in HDSOLUTION and below
C***
C***  called from HDSOLUTION
C***********************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, JOBNUM, INDOUT, KANAL

      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, VELO, AMACH, CKL,
     >                                   ALPHAL, GEDDL, GAMMARAD,
     >                                   F_raw, G_raw, TAUROSS, ENTOT
     
      REAL, DIMENSION(NDMAX) :: RADIUS_LOG, VELO_LOG, AMACH_LOG,
     >                       FMDOT, GMDOT, F_edd, F_a2, F_dadr, FGM,
     >                       G_av, G_k, X, GMDOTr, F_rr, G_rr, RI,
     >                       RHO, RADREL

      REAL, INTENT(IN) :: DOTM4P, XMG, RSTAR
      
      REAL :: DADR, DA2DR, ADUMMY, XMIN, XMAX, XTICK, XABST, 
     >        YMIN, YMAX, YG, YF, XLL, XLR, XSYM, XINT, ATMEAN      
      INTEGER :: I, L, NPAR, FCOLOR, GCOLOR, NDPLOT, NDr, NDr2, hPLOT
      LOGICAL :: bPlotRAW
      
      CHARACTER(100) :: MODHEAD
      CHARACTER(110) :: HEADLINE, PLOTFILENAME
      CHARACTER(40) :: XTEXT
      CHARACTER(20) :: XAXISMODE, CUROPT
      CHARACTER(10) :: CTIME
      CHARACTER(8) :: CENTER, CDATE, CJOB, KEY
      CHARACTER PLOTOPT*(*)


      !Constants
      REAL, PARAMETER :: PI4  = 12.5663706144       !PI4 = 4*PI
      REAL, PARAMETER :: XMSUNPYR = 6.303E25        !1 SOLAR MASS / YR in CGS      
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)
      REAL, PARAMETER :: GCONST = 6.6727E-8     !Gravitational Constant (CGS)

      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOWNPLOT = 2    !current plot

      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTFGSTRAT: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTFGSTRAT: DIMENSION INSUFFICIENT'
         WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDIPMAX = ', NDMAX
         WRITE (hCPR,'(A)') 'PLOTFGSTRAT: PLOTFGSTRAT SKIPPED'
         RETURN
      ENDIF
      
      WRITE (UNIT=CJOB, FMT='(I7)') JOBNUM
      
      IF (KEY == 'TRANSFER') THEN
C***    This plot is part of steal.plot / wruniq.plot      
        CALL JSYMSET ('G2','TRANSFER')      
        hPLOT = KANAL
        WRITE (hPLOT, '(A)') '*NEXTPLOT: ' // PLOTOPT
      ELSE
C***    This plot will be put into its own file      
        hPLOT = hOWNPLOT
        IF (TRIM(KEY) == 'ANALYSE') THEN
          PLOTFILENAME = 'fgparts.'// TRIM(ADJUSTL(CJOB)) //'.plot'
        ELSE
          PLOTFILENAME = 'prehydro.'// TRIM(ADJUSTL(CJOB)) //'.plot'
        ENDIF
        OPEN(hPLOT, FILE=TRIM(PLOTFILENAME), STATUS='UNKNOWN')       
      ENDIF 
            
      XAXISMODE = 'DEPTHINDEX'
      CENTER = '\CENTER\'
      bPlotRAW = (MAXVAL(G_raw) > 0.)

C***  Parse optional parameters      
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
          ENDSELECT
        ENDDO
      ENDIF
      
      HEADLINE = 'HYDRO:'//MODHEAD(13:)

      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A,F6.3)') 
     >   '*\INTERACTIVE'            !auskommentierter interactive-Modus fuer schnellere Analysemglk.
      WRITE (hPLOT, '(A)') '\FONT=HELVET'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR=3 1.0 0.6 0.0'     !orange
      WRITE (hPLOT, '(A)') '\DEFINECOLOR=5 0.0 0.6 0.6'     !gruenblau
      WRITE (hPLOT, '(A)') '\DEFINECOLOR=9 0.6 0.6 0.6'
      WRITE (hPLOT, '(A)') '\PEN=1'

      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0 XMAX 0 0 0 '      !plot zero line in grey
      IF (INDOUT > 0) THEN
        WRITE (hPLOT, '(2(A,I3),A)') 
     >       '\LINUN ', INDOUT, ' YMIN ', INDOUT, ' YMAX 0 0 ' !plot INDOUT line in grey
      ENDIF
      WRITE (hPLOT, '(A)') '\COLOR=1'

      CALL DATE_AND_TIME(CDATE, CTIME)

      WRITE (hPLOT, '(A,I7)') 
     >   '\LUNA XMAX YMAX 0.7 0. 0.4 -90. JOB No.', JOBNUM

      WRITE (hPLOT, '(11A)') 
     >   '\LUNA XMAX YMIN R0.7 0. 0.3 -90. (calculated at ',
     >   CDATE(1:4), '/', CDATE(5:6), '/', CDATE(7:8), ' ', 
     >   CTIME(1:2), ':', CTIME(3:4), ')'


C***  Calculate F (divided by XMG) and G (divided by MDOT)  ***********
C      divided values are dimensionless and therefore easier to plot
            
      DO L=1, ND
        RADREL(L) = RADIUS(L)/RSTAR
        IF (L > 1) THEN
          RI(L-1) = 0.5 * ( RADREL(L-1) + RADREL(L) )      
        ENDIF
        
        CALL SPLINPOX(ADUMMY , RADIUS(L), AMACH, RADIUS, ND, DFDX=DADR)
        DA2DR = 2 * AMACH(L) * DADR
        
        F_edd(L) = (1.-GEDDL(L))
        F_a2(L) = - 2.*AMACH(L)*AMACH(L)*RADIUS(L)/XMG
        F_dadr(L) = RADIUS(L)*RADIUS(L) * DA2DR / XMG
        
        FGM(L) = F_edd(L) + F_a2(L) + F_dadr(L)
        FMDOT(L) = XMG * FGM(L) / DOTM4P
        
        G_av(L) = (1. - (AMACH(L)/VELO(L))**2)
        G_k(L) = - CKL(L)/DOTM4P
        GMDOT(L) = G_av(L) + G_k(L)
        
        F_rr(L) = F_raw(L)/XMG*DOTM4P
      ENDDO      
      
      WRITE(hPLOT,'(2A,F12.6)') '\LUN XMIN YMAX L+1.0 M-1.0 0.2 ',
     >                 'log \M = ', LOG10(DOTM4P * PI4 / XMSUNPYR)
      WRITE(hPLOT,'(A,F12.6)') '\NEXTLUN v\8 = ', VELO(1)/1.E5
     
C***  Boundaries and box settings
      XMIN = 0.
      XMAX = FLOAT(ND) + 1.
      XTICK = 1.                    !Striche
      XABST = 5.                     !Strich mit Wert
      YG = MAXVAL(GMDOT)
      YF = MAXVAL(FGM)
      YMAX = MAX(YG, YF) * 1.5
      YMAX = MAX(1.1, YMAX)
      YMIN = -1. * YMAX             !total minimum usually does not need to be covered
      FCOLOR = 2
      GCOLOR = 4
      NDPLOT = ND

      IF (XAXISMODE == 'DEPTHINDEX' .OR. XAXISMODE == 'L') THEN 
C***     X-Axis = Depth Index
         XTEXT = CENTER//'DEPTH INDEX L'
         XMIN = 0.
         XMAX = FLOAT(ND) + 1.
         XTICK = 1.
         XABST = 10. 
         DO L=1, ND
           X(L) = FLOAT(L)
         ENDDO
      ELSEIF (XAXISMODE == 'VELOCITY' .OR. XAXISMODE == 'V') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v' // CHAR(92) // '8'
         XMIN = 0.
         XMAX = 1.
         XTICK = 0.1
         XABST = 0.5 
         DO L=1, ND
            X(L) = VELO(L)/VELO(1)
         ENDDO
      ELSEIF (XAXISMODE == 'RADIUS' .OR. XAXISMODE == 'R') THEN
         NDPLOT = ND-1
         XTEXT = CENTER//'LOG (R/R\*-1)'
         XMIN = LOG10 ( 0.9 * (RI(ND-1) - 1.) )
         XMAX = LOG10 ( RADREL(1) - 1. )
         XTICK = 0.2
         XABST = 1.0 
         DO L=1, ND-1
           X(L) = LOG10( RADREL(L) - 1. )
         ENDDO
      ELSEIF (XAXISMODE == 'TAUROSS' .OR. XAXISMODE == 'TAU') THEN
c         NDPLOT = ND-1
         XTEXT = CENTER//'log #t#&TRoss&M'
         XMIN = LOG10(0.4 * TAUROSS(2))
         XMAX = LOG10( TAUROSS(ND) )
         XTICK = 0.1
         XABST = 1.0 
         X(1) = XMIN   ! mock value as we cannot set LOG(0)
         DO L=2, ND
           X(L) = LOG10 (TAUROSS(L))
         ENDDO
c         IF (Ltcand > 0 .AND. Ltcand < ND) THEN
c           CALL SPLINPOX(Xts,Rts,X,RADIUS,ND-1)         
c         ENDIF
      ELSEIF (XAXISMODE == 'ENTOT') THEN
         XTEXT = CENTER // 'log(&Rn&N&Ttot&M/cm&H-3&M)'
         XMIN = LOG10(ENTOT(1))
         XMAX = LOG10(ENTOT(ND))
         XTICK = 1.
         XABST = 3.
         DO L=1, ND
           X(L) = LOG10(ENTOT(L))
         ENDDO      
      ELSEIF (XAXISMODE == 'RHO') THEN
         XTEXT = CENTER // 'log(&R#r#&N/(g cm&H-3&M))'
         DO L=1, ND
           RHO(L) = ENTOT(L) * AMU * ATMEAN
         ENDDO         
         XMIN = LOG10(RHO(1))
         XMAX = LOG10(RHO(ND))
         XTICK = 1.
         XABST = 2.
         DO L=1, ND
           X(L) = LOG10(RHO(L))
         ENDDO      
c         IF (Lcand > 1 .AND. Lcand < ND) THEN
c           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
c         ENDIF           
      ELSE
         WRITE (hCPR,*) 'PLOTFGSTRAT: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) TRIM(PLOTOPT)
         RETURN
      ENDIF      
      
      DO L=ND, 1, -1
        IF (GMDOT(L) > YMIN) THEN
          NDr = L+1
          EXIT
         ENDIF
      ENDDO

      DO L=ND, 1, -1
        IF (G_raw(L) > YMIN) THEN
          NDr2 = L+1
          EXIT
         ENDIF
      ENDDO

C***  write reduced arrays for contionous G line plots (full will usually crash WRPlot)
      DO L=1, NDr
        GMDOTr(L) = GMDOT(L)
      ENDDO
      DO L=1, NDr2
        G_rr(L) = G_raw(L)
      ENDDO
      
      XLL = XMAX
      XLR = XLL + XTICK * 1.5
      XSYM = XLL + XTICK * 0.6
                  
      
C***  Explanation of different lines            
      WRITE (hPLOT, '(A)') '\COLOR=2'
      WRITE (hPLOT, '(A,F6.3,A,F6.3,A)')
     >              '\LINUN ',XLL,' YMAX ',XLR,' YMAX 2.0 -1.0'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-1.0 0.4 &2\~'
      WRITE (hPLOT, '(A)') '\NEWSIZELUN 0.34'
      WRITE (hPLOT, '(A)') '\SAMELUN F'
      WRITE (hPLOT, '(A)') '\COLOR=3'
      WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >           '\LINUN ',XLL,' YMAX ',XLR,' YMAX 2.0 -1.8 SYMBOL=5'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-1.8 0.2 &3#G#-Term'
      WRITE (hPLOT, '(A,F6.3,A,F6.3,A)')
     >  '\LINUN ',XLL,' YMAX ',XLR,' YMAX 2.0 -2.6 SYMBOL=10 SIZE=0.05'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-2.6 0.2 &3a&H2&M-Term'
      WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >  '\LINUN ',XLL,' YMAX ',XLR,' YMAX 2.0 -3.4 SYMBOL=20 SIZE=0.02'
      WRITE (hPLOT, '(A)') 
     >           '\LUN XMAX YMAX L3.0 M-3.4 0.2 &3da&H2&M/ds-Term'
      WRITE (hPLOT, '(A)') '\COLOR=4'
      WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >           '\LINUN ',XLL,' YMAX ',XLR,' YMAX 2.0 -5.0'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-5.0 0.4 &4\~'
      WRITE (hPLOT, '(A)') '\NEWSIZELUN 0.34'
      WRITE (hPLOT, '(A)') '\SAMELUN G'
      IF (MAXVAL(ALPHAL) > 0.) THEN
        WRITE (hPLOT, '(A)') '\COLOR=5'
        WRITE (hPLOT, '(A,F6.3,A)') '\SYM ',XSYM,' YMAX 2.0 -5.8 0.25 3'
        WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-5.8 0.2 &5av-Term'
        WRITE (hPLOT, '(A,F6.3,A)') '\SYM ',XSYM,' YMAX 2.0 -6.6 0.25 4'
        WRITE (hPLOT, '(A)') '\LUN XMAX YMAX L3.0 M-6.6 0.2 &5k-Term'
      ENDIF
      WRITE (hPLOT, '(A)') '\COLOR=1'
      WRITE (hPLOT, '(A)') '\INBOX'
      
C***  Plot data curves and symbols      
      CALL    PLOTANFS(hPLOT,'', '&E'//HEADLINE,
     >                 XTEXT,
     >                 '',
     >                 0., XMIN, XMAX, XTICK, XABST, 0.,
     >                 0., YMIN, YMAX, 0.25, 0.5, 0.,
     >                 X,GAMMARAD,NDPLOT, 'COLOR=9 SYMBOL=9 SIZE=0.15')

      CALL PLOTCONS(hPLOT,X,F_edd,NDPLOT,'COLOR=3 SYMBOL=5')
      CALL PLOTCONS(hPLOT,X,F_a2,NDPLOT,'COLOR=3 SYMBOL=10 SIZE=0.1')
      CALL PLOTCONS(hPLOT,X,F_dadr,NDPLOT,'COLOR=3 SYMBOL=20 SIZE=0.02')
      CALL PLOTCONS(hPLOT,X,FGM,NDPLOT,'COLOR=2 SYMBOL=6')

      IF (MAXVAL(ALPHAL) > 0.) THEN
        CALL PLOTCONS(hPLOT,X,G_av,NDPLOT,'COLOR=5 SYMBOL=3')
        CALL PLOTCONS(hPLOT,X,G_k,NDPLOT,'COLOR=5 SYMBOL=4')
      ENDIF
      CALL PLOTCONS(hPLOT,X,GMDOT,NDPLOT,'COLOR=4 SYMBOL=2')
      IF (NDr > 0) THEN
        CALL PLOTCONS(hPLOT,X,GMDOTr,NDr,'COLOR=4 SYMBOL=6 SIZE=0.1')
      ENDIF
      
      IF (bPlotRAW) THEN
        CALL PLOTCONS(hPLOT,X,F_rr,NDPLOT,'COLOR=2 SYMBOL=20 SIZE=0.02')
        CALL PLOTCONS(hPLOT,X,G_rr,NDr2,'COLOR=4 SYMBOL=20 SIZE=0.02')
      ENDIF

      IF (KEY /= 'TRANSFER') CLOSE(hPLOT)
      
      RETURN
      END
      
