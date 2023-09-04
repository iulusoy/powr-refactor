      SUBROUTINE PLOTALPHA (PLOTOPT, ND, RADIUS, ALPHAF, VMACH,
     >                      VELO, TAUROSS, ENTOT, ATMEAN,
     >                      MODHEAD, JOBNUM, bOwn)
C******************************************************************************
C***  PLOT FORCE MULTIPLIER PARAMETER ALPHA VS DEPTH POINT OR LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, DIMENSION(ND) :: RADIUS, ALPHAF, VSCRATCH, RI,
     >                       VELO, VMACH, TAUROSS, ENTOT, RHO
      LOGICAL, INTENT(IN) :: bOwn           !own plot file or part of wruniq.plot?

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2
      CHARACTER(40) :: XTEXT
      CHARACTER(10) :: CTIME
      CHARACTER(8) :: CDATE, CENTER
      CHARACTER(20) :: CUROPT, XAXISMODE      
      CHARACTER PLOTOPT*(*)

      INTEGER :: L, NPAR, I, Lcand
      REAL :: XMIN, XMAX, YMIN, YMAX, Xsonic, Rsonic, Vsonic, ATMEAN,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST

      INTEGER, EXTERNAL :: IDX
     
      !Physical constants
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)     
     
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hPLOT = 2       !write to plot file
 
C***  INITIALIZATION
      IF (bOwn) THEN
        OPEN (hPLOT, FILE='alpha.plot', STATUS='UNKNOWN')
        WRITE (hCPR, '(A)') 'FORCE MULTIPLIERS PLOTTED IN alpha.plot'
      ELSE
        OPEN (hPLOT, FILE='PLOT', STATUS='UNKNOWN')
        CALL JSYMSET ('G2','TRANSFER')
        CALL REMARK ('FORCE MULTIPLIERS TO BE ROUTED')
      ENDIF
        
      XAXISMODE = 'DEPTHINDEX'
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR > 2) THEN
        parloop: DO I=3, NPAR
          CALL SARGV (PLOTOPT, I, CUROPT)
          IF (CUROPT == 'XAXISMODE' .OR. CUROPT == 'X') THEN
            IF (NPAR >= (I+1)) THEN
              CALL SARGV (PLOTOPT, I+1, CUROPT)
            ELSE
              CYCLE parloop
            ENDIF
          ENDIF
          SELECTCASE (CUROPT)
            CASE ('VELOCITY', 'VELO', 'V')
              XAXISMODE = 'VELOCITY'
            CASE ('TAUROSS', 'TAU')
              XAXISMODE = 'TAUROSS'
            CASE ('RADIUS', 'R')
              XAXISMODE = 'RADIUS'
            CASE ('DEPTHINDEX', 'INDEX', 'ND')
              XAXISMODE = 'DEPTHINDEX'
            CASE ('N', 'NTOT', 'ENTOT', 'PARTDENS')
              XAXISMODE = 'ENTOT'
            CASE ('RHO', 'DENS', 'DENSITY')
              XAXISMODE = 'RHO'
          ENDSELECT
        ENDDO parloop
      ENDIF
 
 
      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        IF (L /= ND) THEN
          RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
        ENDIF
        VSCRATCH(L) = VELO(L) - VMACH(L)        
        IF ((VSCRATCH(L) < 0) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,VMACH,RADIUS,ND)
      ENDIF
      
      
C***  HEADER  ------------------------------------------------------
      HEAD1=' FORCE MULTIPLIER PARAMETER ALPHA VERSUS DEPTH INDEX'
      HEAD2      = MODHEAD(13:32)
      HEAD2(22:) = 'FORCE MULTIPLIER ALPHA(L)'

      CENTER = '\CENTER\'
      
C***  X-AXIS: ----------------------------------------------------------
C***  DEPTH POINT INDEX
C***  @TODO: IMPLEMENT RADIUS OPTION

      IF (XAXISMODE == 'DEPTHINDEX') THEN 
C***     X-Axis = Depth Index
         XTEXT = CENTER//'DEPTH INDEX L'
         XMIN = 0.
         XMAX = FLOAT(ND)
         XTICK = 5.
         XABST = 10. 
         DO L=1, ND-1
            X(L) = FLOAT(L) + 0.5
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RADIUS,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'VELOCITY') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v\8'
         XMIN = 0.
         XMAX = 1.
         XTICK = 0.1
         XABST = 0.5 
C         DO L=1, ND-2
         DO L=1, ND-1
            X(L) = 0.5 * (VELO(L) + VELO(L+1)) / VELO(1)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = Vsonic/VELO(1)
         ENDIF
      ELSEIF (XAXISMODE == 'TAUROSS') THEN
         XTEXT = CENTER//'LOG Rosseland Optical Depth #t#'
         XMIN = LOG10(0.9 * TAUROSS(2))
         XMAX = LOG10( TAUROSS(ND) )
         XTICK = 0.1
         XABST = 1.0 
         DO L=1, ND-1
           X(L) = LOG10 (0.5 * (TAUROSS(L) + TAUROSS(L+1)))
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RADIUS,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'RADIUS') THEN
         XTEXT = CENTER//'LOG (R/R\*-1)'
         XMIN = LOG10 ( 0.9 * (RI(ND-1) - 1.) )
         XMAX = LOG10 ( RADIUS(1) - 1. )
         XTICK = 0.1
         XABST = 0.5 
         DO L=1, ND-1
           X(L) = LOG10( RI(L) - 1. )
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = LOG10( Rsonic - 1. )
         ENDIF
      ELSEIF (XAXISMODE == 'ENTOT') THEN
         XTEXT = CENTER // 'log(&Rn&N&Ttot&M/cm&H-3&M)'
         XMIN = LOG10(ENTOT(1))
         XMAX = LOG10(ENTOT(ND))
         XTICK = 1.
         XABST = 3.
         DO L=1, ND-1
           X(L) = LOG10(0.5 * (ENTOT(L) + ENTOT(L+1)))
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RADIUS,ND-1)         
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
           X(L) = LOG10(0.5 * (RHO(L) + RHO(L+1)))
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RADIUS,ND-1)         
         ENDIF      
      ELSE
         WRITE (hCPR,*) 'PLOTACC: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))
         RETURN
      ENDIF
      
      
C***  Y-AXIS:  ---------------------------------------------------------
C***  ALPHA RANGE (AUTO)
      YMAX = 0. 
      YMIN = .0
      YSCALE = 0.
      YTICK=2.5 !unwichtig, da AUTO-Option gesetzt
      YABST=5.  !unwichtig, da AUTO-Option gesetzt

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        Y(L) = ALPHAF(L)
      ENDDO
 
      WRITE (hPLOT, '(A,A)') 'PLOT: ', HEAD1
      WRITE (hPLOT, '(A)') '\INBOX'
      WRITE (hPLOT, '(A)') '\FONT=HELVET'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR=9 0.5 0.5 0.5'
      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0. XMAX 0. 0. 0. SYMBOL=9'
      WRITE (hPLOT, '(A)') '\COLOR=1'

      CALL DATE_AND_TIME(CDATE, CTIME)

      WRITE (hPLOT, '(A,I7)') 
     >   '\LUNA XMAX YMAX 0.7 0. 0.4 -90. JOB No.', JOBNUM

      WRITE (hPLOT, '(11A)') 
     >   '\LUNA XMAX YMIN R0.7 0. 0.3 -90. (calculated at ',
     >   CDATE(1:4), '/', CDATE(5:6), '/', CDATE(7:8), ' ', 
     >   CTIME(1:2), ':', CTIME(3:4), ')'

      CALL PLOTANF (hPLOT,HEAD1,HEAD2
     >        ,XTEXT
     >        ,'\CENTER\force multiplier parameter #a#'
     >        ,0.,XMIN,XMAX,XTICK,XABST,.0
     >        ,YSCALE,0.,0.,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
 
C      CLOSE(hPLOT)
      
      RETURN
      END
