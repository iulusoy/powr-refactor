      SUBROUTINE PLOTRTAU1COLI (hPLOT, K, XLAMK, XLAMKOLD, CMODE, 
     >                          ND, RADIUS, OPAK, OPA, FILLFAC, 
     >                          MAINPRO, MAINLEV, MODHEAD)
C***********************************************************************
C***  PLOT of the radius where tau=1 ("Nick-Knatterton-Diagram")
C***  In contrast to the COMO version, this one uses the full opacity
C***  The plot can be requested by the CARDS line: PLOT RTAU1COLI

C***  called from: COLI
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: hPLOT, ND, K
      REAL, INTENT(IN) :: XLAMK, XLAMKOLD
      
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, OPAK, OPA, FILLFAC
      CHARACTER(1), INTENT(IN) :: CMODE
      CHARACTER(100), INTENT(IN) :: MODHEAD

      CHARACTER(10), DIMENSION(ND) :: MAINPRO, MAINLEV
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER, PARAMETER :: MAXIDENT = 100
      CHARACTER(45), DIMENSION(MAXIDENT) :: IDLINE
 
      REAL :: TAU, TAUOLD, OPAL, OPALP, OPAC, OPACP, Q, DUMMY, RTAU1,
     >        XTICK, XABST, YTICK, YABST, XMIN, XMAX, YBASE,
     >        RTAU1C, TAUC, TAUCOLD, TMPLAM, TMPRTAU1, TMPRTAU1C
      INTEGER :: LTAU1, L, I, IOST

C***  This routine makes use of a few SAVE variables
C***  which must keep its value for the next call
      CHARACTER(10), SAVE :: LASTMAINPRO, LASTMAINLEV
      REAL, SAVE :: RTAU1LAST, RTAU1CLAST
      INTEGER, SAVE :: NIDENT, NDATA

 
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
C***  Temporary file for second dataset            
      INTEGER, PARAMETER :: hRT1C = 150
 
C***  INITIALIZATION
      IF (K == 0) THEN
         RTAU1LAST = 1.0
         RTAU1CLAST = 1.0
         NIDENT=0
         NDATA=0
         OPEN (hRT1C, FILE='rtau1coli.tmp', STATUS='UNKNOWN')
         CALL JSYMSET ('G2','TRANSFER')
         CALL REMARK ('RTAU1-PLOT TO BE ROUTED')

      ENDIF

      TAU=0.
      TAUC=0.
      RTAU1 = 1.
      RTAU1C = 1.
      LTAU1 = ND 

C***  RTAU1 for lines and continua:
      dl: DO L=1, ND-1
         TAUOLD = TAU
         OPAL = OPAK(L) * FILLFAC(L)
         OPALP = OPAK(L+1) * FILLFAC(L+1)
         TAU = TAU + 0.5*(OPAL+OPALP)*(RADIUS(L)-RADIUS(L+1))
         
         IF (TAU >= 1.) THEN
            Q = (1.-TAUOLD) / (TAU-TAUOLD)
            RTAU1 = (1.-Q) * RADIUS(L) + Q * RADIUS(L+1)
c        Do not store LTAU1, but set marks only via continua
c            LTAU1 = L + 1            
            EXIT dl
         ENDIF
      ENDDO dl

C***  RTAU1 for continua only:      
      dlc: DO L=1, ND-1
         TAUCOLD = TAUC
         OPAC = OPA(L) * FILLFAC(L)
         OPACP = OPA(L+1) * FILLFAC(L+1)
         TAUC = TAUC + 0.5*(OPAC+OPACP)*(RADIUS(L)-RADIUS(L+1))
         
         IF (TAUC >= 1.) THEN
            Q = (1.-TAUCOLD) / (TAUC-TAUCOLD)
            RTAU1C = (1.-Q) * RADIUS(L) + Q * RADIUS(L+1)
            WRITE(hRT1C, '(F12.4, 2(1X, F12.4))') 
     >           ALOG10(XLAMK), RTAU1, RTAU1C
            LTAU1 = L + 1            
            EXIT dlc
         ENDIF
      ENDDO dlc
        

      IF (K > 0) THEN
        IF (LASTMAINPRO .NE. MAINPRO(LTAU1) 
     >       .OR. LASTMAINLEV .NE. MAINLEV(LTAU1) ) THEN 
ccc     >       .OR. RTAU1 .LT. RTAU1LAST*0.95) THEN 
          IF (NIDENT .LT. MAXIDENT) THEN
           NIDENT = NIDENT + 1
           WRITE (IDLINE(NIDENT), '(A, 1X, 1PG14.6, 1X, A)')
     >      '\IDENT', XLAMKOLD, '&E'//LASTMAINPRO//' '//LASTMAINLEV 
          ENDIF
        ENDIF
      ENDIF
      LASTMAINPRO = MAINPRO(LTAU1)
      LASTMAINLEV = MAINLEV(LTAU1)
      RTAU1LAST = RTAU1

      IF (CMODE == 'F') THEN 
C***     AUTO-SCALING by WRplot
         XMIN = .0
         XMAX = .0
 
         XTICK = 0.1
         XABST = 1.
         YBASE = 10**(FLOOR(LOG10(RADIUS(1))))
         YTICK = YBASE/50.
         YABST = YBASE/10.
         OPEN (hPLOT, FILE='PLOT', STATUS='UNKNOWN')
         WRITE (hPLOT, '(A)') 'PLOT: RTAU1COLI'          
         WRITE (hPLOT, '(A)') '\SET_NDATMAX=1000000'          
         CALL PLOTANF (hPLOT,'RTAU1COLI', '&E'//MODHEAD
     >        ,'\CENTER\log #l#/\A'
     >        ,'\CENTER\Radius where #t# = 1'
     >        ,0., XMIN, XMAX, XTICK, XABST,.0
     >        ,0., 1.,   RADIUS(1), YTICK, YABST, .0
     >        ,DUMMY, DUMMY, 0, 0)
         BACKSPACE(hPLOT)
         WRITE (hPLOT, '(A)') 'N=? XYTABLE COLOR=4 PEN=5' 
C***     Read continua dataset from temporary file
         REWIND (hRT1C)
         rlc: DO
           READ(hRT1C, '(F12.4, 2(1X, F12.4))', IOSTAT=IOST)
     >       TMPLAM, TMPRTAU1, TMPRTAU1C
           IF (IOST < 0) EXIT rlc
           WRITE(hPLOT, '(F12.4, 1X, F12.4)') TMPLAM, TMPRTAU1C
         ENDDO rlc
         WRITE (hPLOT, '(A)') 'FINISH'
         WRITE (hPLOT, '(A)') 'N=? XYTABLE COLOR=2 PEN=3' 
C***     Read continua + line dataset from temporary file
         REWIND (hRT1C)
         rll: DO
           READ(hRT1C, '(F12.4, 2(1X, F12.4))', IOSTAT=IOST)
     >       TMPLAM, TMPRTAU1, TMPRTAU1C
           IF (IOST < 0) EXIT rll
           WRITE(hPLOT, '(F12.4, 1X, F12.4)') TMPLAM, TMPRTAU1
         ENDDO rll
         WRITE (hPLOT, '(A)') 'FINISH'

         WRITE (hPLOT, '(A)') ' '
         WRITE (hPLOT, '(A)') '* Identification of Edges'
         WRITE (hPLOT, '(A)') '\IDY 10'
         WRITE (hPLOT, '(A)') '\IDSIZE 0.2'
         WRITE (hPLOT, '(A)') '\IDLOG'
         DO I=1, NIDENT
            WRITE (hPLOT, '(A)') IDLINE(I)
         ENDDO
         WRITE (hPLOT, '(A)') 'END'
         
         CLOSE (hRT1C)
      ENDIF



      RETURN
      END
