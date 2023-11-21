      PROGRAM MAINwrstart 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:50 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL wrstart 
      END
      SUBROUTINE ADDHISTENTRY(MODHIST,LAST,MAXHIST,ENTRYLEN,ENTRYSTR)
C***********************************************************************
C***  ADDS A STRING ENTRY TO THE MODEL HISTORY CHARACTER ARRAY
C     ENTRYSTR: new string to add
C     ENTRYLEN: length of ENTRYSTR
C     
C     Note: This routine updates MODHIST and the historical LAST parameter
C           to be compartible with older code that still uses ENDCODE/DECODE 
C           and declares MODHIST as an integer array which is filled with
C           Hollerith constants
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAXHIST, ENTRYLEN      
      CHARACTER(ENTRYLEN), INTENT(IN) :: ENTRYSTR

      INTEGER, INTENT(INOUT) :: LAST
      CHARACTER(8*MAXHIST), INTENT(INOUT) :: MODHIST

      INTEGER :: LASTCHAR, INFROM, INTO, BUFFERINT, CURLAST
      REAL :: ADDBYTES
      CHARACTER(8) :: BUFFER8


      !LAST can be set to -1 => read from MODHIST(1:8)
      IF (LAST < 0) THEN
        BUFFER8 = MODHIST(1:8)
        READ(UNIT=BUFFER8, FMT='(A8)') BUFFERINT
        CURLAST = BUFFERINT
      ELSE
        CURLAST = LAST
      ENDIF

      LASTCHAR = CURLAST * 8
      
      INFROM = LASTCHAR + 1
      INTO = LASTCHAR + ENTRYLEN
      MODHIST(INFROM:INTO) = ENTRYSTR

      ADDBYTES = REAL(ENTRYLEN) / 8.
      CURLAST = CURLAST + INT(ADDBYTES)
      IF (ADDBYTES - REAL(INT(ADDBYTES)) > 0) THEN
        !Entry has a length that is not a multiple of 8, one more byte needed)
        CURLAST = CURLAST + 1
      ENDIF

      !MODHIST(1:8) or MODHIST(1) in integer array definition contains the 
      !currently used length of MODHIST (in Bytes i.e. in CHARs / 8)
      ! The first bytes is for historical written as an integer, 
      !  NOT as a character containing an integer number
      !  THerefore (A8) is used als FORMAT instead of (I8) 
      WRITE(UNIT=BUFFER8, FMT='(A8)') CURLAST          
      MODHIST(1:8)=BUFFER8

      IF (LAST >= 0) THEN
        !Fill LAST only if not called with -1 (otherwise COLI will crash on subroutine call)
        LAST = CURLAST
      ENDIF

      RETURN
      END
      SUBROUTINE APPEND_AUTOLEVELS (N, N_WITH_DRLEVELS, NDIM, MAXIND, 
     $                  MAXAUTO, NAUTO, LOWAUTO, IONAUTO, EAUTO, ELEVEL, 
     $                  LEVEL, EION, WEIGHT, INDLOW, INDNUP, LASTIND,
     >                  LEVUPAUTO, LEVAUTO, WAUTO, 
     >                  NCHARG, IONGRND, NOM)
C******************************************************************************
C***  This subroutine appends autoionization levels 
C***  to the vectors LEVEL, WEIGHT and ELEVEL
C***  Index range: N+1 ... N_WITH_DRLEVELS
C***  THIS SUBROUTINE IS ESSENTIAL TO BE CALLED FROM COLI, 
C**   but also called from STEAL - PRIDAT in order to displaying 
C**   the autoionization levels ("DR-LEVELS") with PRINT DATOM 
C******************************************************************************

      DIMENSION WEIGHT(NDIM), ELEVEL(NDIM), EION(NDIM)
      DIMENSION NCHARG(NDIM), IONGRND(NDIM), NOM(NDIM)
      DIMENSION INDNUP(MAXIND),INDLOW(MAXIND)
      DIMENSION LOWAUTO(MAXAUTO),IONAUTO(MAXAUTO), WAUTO(MAXAUTO)
      DIMENSION EAUTO(MAXAUTO)
      CHARACTER*10 LEVEL(NDIM)
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

C***  Copy vector LEVUPAUTO to LEVAUTO, skipping multiple occurrences 
      NLEVEL_AUTO = 1
      LEVAUTO(1) = LEVUPAUTO(1)
      DO 5 ILEVEL=2, NAUTO
         DO ITEST=1, NLEVEL_AUTO 
            IF (LEVUPAUTO(ILEVEL) .EQ. LEVAUTO(ITEST)) GOTO 5
         ENDDO
         NLEVEL_AUTO = NLEVEL_AUTO + 1
         LEVAUTO(NLEVEL_AUTO) = LEVUPAUTO(ILEVEL)
    5 CONTINUE

C***  Number of levels is increased
      N_WITH_DRLEVELS = N + NLEVEL_AUTO
      IF (N_WITH_DRLEVELS .GT. NDIM) THEN
         WRITE (*,'(A, I5,A,I5)') 
     >    '*** NDIM=', NDIM, ' insufficient: needed', N_WITH_DRLEVELS
      ENDIF

C***  The new AUTO-Levels are appended to LEVEL list
      DO ILEVEL = 1, NLEVEL_AUTO
C***     ... the list of levels
         LEVEL(N+ILEVEL) = LEVAUTO(ILEVEL)
      ENDDO

C***  further level attributes are same as for lower level 
C***    of DRTRANSIT with same upper level name
      DO ILEVEL= N+1, N+NLEVEL_AUTO
         DO INDDR = 1, NAUTO
            IF (LEVEL(ILEVEL) .EQ. LEVUPAUTO(INDDR)) THEN
              NCHARG (ILEVEL) = NCHARG (LOWAUTO(INDDR))
              EION   (ILEVEL) = EION   (LOWAUTO(INDDR))
              IONGRND(ILEVEL) = IONGRND(LOWAUTO(INDDR))
              NOM    (ILEVEL) = NOM    (LOWAUTO(INDDR))

              WEIGHT (ILEVEL) = WAUTO (INDDR)
              ELEVEL (ILEVEL) = EAUTO (INDDR)

              GOTO 10
            ENDIF
         ENDDO
         STOP '*** FATAL INTERNAL ERROR IN subr. append_autolevels'
   10    CONTINUE 
      ENDDO

C***  NEW!! for stabilizing transitions, INDNUP points to the
C***      auto-ionizing level, and NOT to the next-higher ion
      DO 11 INDDR = 1, NAUTO
         DO ILEVEL= N+1, N+NLEVEL_AUTO
            IF (LEVEL(ILEVEL) .EQ. LEVUPAUTO(INDDR)) THEN
               IND = INDDR + LASTIND
               INDNUP (IND) = ILEVEL
               GOTO 11
            ENDIF
         ENDDO
         STOP '*** FATAL INTERNAL ERROR2 IN subr. append_autolevels'
   11 CONTINUE

C***  Energies of the appended auto-ionizing levels now relative 
C***   to the ground level, as for "normal" levels
      DO ILEVEL = N+1, N_WITH_DRLEVELS
         ELEVEL(ILEVEL) = ELEVEL(ILEVEL) + EION(ILEVEL)
      ENDDO

      RETURN
      END 
      SUBROUTINE BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $                    XLAMBDA,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,
     $                    KONTNUP,KONTLOW,LASTKON)
C***********************************************************************
C***  THIS ROUTINE PREPARES AN ARRAY SIGMAKI WITH THE BOUND-FREE CROSS SECTIONS
C***  ( IN CM**2) TO AVOID UNNECCESSARY MULTIPLE CALCULATIONS
C***********************************************************************
 
      DIMENSION EION(NDIM),ELEVEL(NDIM),EINST(NDIM,NDIM)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION XLAMBDA(NF)
      DIMENSION SIGMAKI(NF,LASTKON)
 
C***  LOOP OVER ALL CONTINUUM TRANSITIONS
      DO 8 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
C***  EINST = THRESHOLD CROSS SECTION IN 10**-18 CM**2
      SIGMATH=EINST(LOW,NUP)*1.E-18
C***  EDGE = THRESHOLD ENERGY IN KAYSER *****
      EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)

C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS
      DO 9 K=1,NF
      WAVENUM=1.E8/XLAMBDA(K)
      IF (WAVENUM .LT. EDGE) THEN
         SIGMAKI(K,KON)=.0
      ELSE
         CALL PHOTOCS (SIGMAKI(K,KON),SIGMATH,EDGE,WAVENUM,ALPHA,
     $                 SEXPO,
     $                 ADDCON1, ADDCON2, ADDCON3, 
     $                 IGAUNT,KON)
      ENDIF
    9 CONTINUE

    8 CONTINUE
 
      RETURN
      END
      FUNCTION BNUE(XLAMBDA,T)
C***********************************************************************
C***  PLANCK FUNCTION, LAMBDA IN ANGSTROEM, T IN KELVIN
C***  BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
C***********************************************************************
      DATA C1,C2 / 1.4388E8, 3.9724E+8 /
C***      DATA EXPMAX / 650. /
C***  geaendert wrh  5-Jul-2013 
      DATA EXPMAX / 500. /

C***  PREVENT ERROR FOR NEGATIVE TEMPERATURES
      IF (T .LE. .0) THEN
         BNUE=.0
      ELSE
         W=1./XLAMBDA
         HNUEKT=C1*W/T
C***  PREVENT ERROR (BAD SCALAR ARGUMENT TO ARLIB MATH ROUTINE)
C***  FOR EXTREME WIEN-DOMAIN
         IF (HNUEKT .GT. EXPMAX) THEN
           FAC=C2*W*W*W
           FACLOG=ALOG(FAC)
           ARG= HNUEKT - FACLOG
           IF (ARG .GT. EXPMAX) THEN
            BNUE = EXP(-EXPMAX)
           ELSE  
            BNUE = EXP(-ARG)
           ENDIF  
         ELSE IF (HNUEKT .LT. 1.E-10) THEN
            BNUE=C2/HNUEKT*W*W*W
         ELSE
            BNUE=C2/(EXP(HNUEKT)-1.)*W*W*W
         ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE CLOCK

      STOP 'CLOCK NOT IMPLEMENTED AT DEC/OSF'

      RETURN
      END
      SUBROUTINE CLOSMS (ICHANNEL, IERR)
C************************************************************
C***  ROUTINE VON LARS KOESTERKE      8-Sep-1995 15:50:47
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, IDUMMY, IDUMMY, DUMMY, 
     >              IDUMMY, 'CLOSE', IERR)

      RETURN
      END
      SUBROUTINE CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, VELO, 
     >                         TAUROSS, DENSCON_LINE, RADIUS, T, XMU)

C************************************************************************
C***  This routine allows to specify a density stratification of 
C***  the clumping constrast, DENSCON(L). 
C***  The subroutine interpretes directly the DENSCON-line from the  
C***  CARDS file, which is handed over from DECSTAR via the character
C***  variable DENSCON_LINE. 
C***  Default is depth-independent clumping with DENSCON_FIX, 
C***  which is the second parameter on the DENSCON-line and decoded 
C***  already in Subr. DECSTAR
C************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      REAL, DIMENSION(ND) :: DENSCON, FILLFAC, VELO, DCSMOOTH,
     >                       TAUROSS, RADIUS, T, XMU  !NOTE: Tauross is Tauross_cont
      CHARACTER(10) :: CLUMP_CRIT, CLUMP_CRIT2 
      CHARACTER(20) :: ACTPAR
      CHARACTER*(*) :: DENSCON_LINE

C***  Local arrays      
      INTEGER, PARAMETER :: NDIPMAX = 100
      REAL, DIMENSION(NDIPMAX) :: AMACH, VHELP 
      
      INTEGER :: L, IPAR, ND, NPAR, NDSTART, Lcand, LENPAR
      REAL :: DENSCON_FIX, PAR1, PAR2, PARL, D2, F1, F8, X1, X2, 
     >        DX, X, Q, Rsonic, Vsonic, TAUsonic
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      REAL, PARAMETER :: PI = 3.141592654
      REAL, PARAMETER :: RGAS = 8.3145E7            !GAS CONSTANT in CGS 


      IF (ND > NDIPMAX) THEN
         WRITE (hCPR,'(A)') 'CLUMP_STRUCT: FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'CLUMP_STRUCT: NDIPMAX INSUFFICIENT'
         WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDIPMAX = ', NDIPMAX
         STOP 'FATAL ERROR IN CLUMP_STRUCT'
      ENDIF
      
C***  Decoding the input line with DENSCON specifications
      CALL SARGC (DENSCON_LINE, NPAR)
C***  Clumping not depth-dependent?
      IF (NPAR .LE. 2) THEN
         DO L = 1, ND
         DENSCON(L) = DENSCON_FIX
         FILLFAC(L) = 1. / DENSCON_FIX
         ENDDO
         GOTO 100
      ENDIF

      CALL SARGV (DENSCON_LINE, 3, CLUMP_CRIT) 


      Lcand = 0
      Rsonic = 1.
      Vsonic = VELO(ND)
      DO L=1, ND
C***    Calculate speed of sound for all depth points)
        AMACH(L) = SQRT(RGAS * T(L) / XMU(L)) / 1.E5
        VHELP(L) = VELO(L) - AMACH(L)
        IF ((VHELP(L) < 0.) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO                
C***  Find sonic point parameters
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VHELP,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,AMACH,RADIUS,ND)
        CALL SPLINPOX(TAUsonic,Rsonic,TAUROSS,RADIUS,ND)
      ENDIF

C***  Branch to simulate Hillier's formula, given in Martins et al. 
C***  (2004, A&A 420, 1087; SMC-Paper)       
      IF (CLUMP_CRIT == 'HILLIER') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the sonic speed            
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = PAR1 * Vsonic
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = Vsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
               FILLFAC(L) = F8 + (1.-F8)* EXP(-VELO(L)/PAR1)
               DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF

C***  Branch to simulate Paco's (F. Najarro's) FASTWIND clumping strat.
C***  from Najarro et al. (2009, ApJ 691, 1816)
C***  SYNTAX:   DENSCON [D1] PACO V1 D2 V2
C***      e.g.  DENSCON 10 PACO 100. 4  300.
C***            DENSCON 10 PACO 2.5  1. 2.0  (more typical use)
      IF (CLUMP_CRIT == 'NAJARRO' .OR. CLUMP_CRIT == 'PACO') THEN
         IF (NPAR < 6) GOTO 93
         CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
         LENPAR = LEN_TRIM(ACTPAR)
         IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***        Parameter is interpreted as a fraction of the sonic speed            
            READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
            PAR1 = PAR1 * Vsonic
         ELSEIF (ACTPAR == 'SONIC') THEN
            PAR1 = Vsonic
         ELSE
            READ (ACTPAR, '(F20.0)', ERR=93) PAR1
         ENDIF
         IF (PAR1 <= .0) GOTO 95 

C***     Read outermost clumping factor (reached only in the limit)
         CALL SARGV (DENSCON_LINE, 5, ACTPAR)
         READ (ACTPAR, '(F20.0)', ERR=93) D2

         CALL SARGV (DENSCON_LINE, 6, ACTPAR) 
         LENPAR = LEN_TRIM(ACTPAR)
         IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***        Parameter is interpreted as a fraction of the sonic speed            
            READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR2
            PAR2 = PAR2 * Vsonic
         ELSEIF (ACTPAR == 'SONIC') THEN
            PAR2 = Vsonic
         ELSE
            READ (ACTPAR, '(F20.0)', ERR=93) PAR2
         ENDIF
         IF (PAR2 <= .0) GOTO 95 

         F1 = 1. / DENSCON_FIX
         F8 = 1. / D2
         DO L = 1, ND
            FILLFAC(L) = F1 + (1.-F1)*EXP(-VELO(L)/PAR1) 
     >                      + (F8-F1)*EXP((VELO(L)-VELO(1))/PAR2)
            DENSCON(L) = 1. / FILLFAC(L)
         ENDDO

         GOTO 100
      ENDIF

C***  Exponential clumping onset (like in Hillier's formula)
C***  but using the radius scale
      IF (CLUMP_CRIT == 'EXPRADIUS') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the 
C***          distance of the sonic radius from RSTAR
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = 1. + PAR1 * (Rsonic - 1.)
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = Rsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
               FILLFAC(L) = F8 + (1.-F8)*EXP(-(RADIUS(L)-1.)/(PAR1-1.))
               DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF

C***  Exponential clumping onset (like in Hillier's formula)
C***  but using the Tauross_cont scale
      IF (CLUMP_CRIT == 'EXPTAU') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the 
C***          sonic Tau value
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = PAR1 * TAUsonic
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = TAUsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
              IF (TAUROSS(L) <= 1.E-30) THEN
                FILLFAC(L) = F8
              ELSE
                FILLFAC(L) = F8 + (1.-F8)*EXP(-PAR1/TAUROSS(L))
              ENDIF
              DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF
      
      
C***  Stratification needs enough parameters!
      IF (NPAR .LT. 5) GOTO 98
      CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
      IF (CLUMP_CRIT == 'VELO' .AND. ACTPAR == 'SONIC') THEN
        !Clumping starts at sonic point        
        PAR1 = Vsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .AND. ACTPAR == 'SONIC') THEN
        CALL SPLINPOX(PAR1,Rsonic,TAUROSS,RADIUS,ND)
      ELSEIF (CLUMP_CRIT == 'RADIUS' .AND. ACTPAR == 'SONIC') THEN
        PAR1 = Rsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .OR. CLUMP_CRIT == 'VELO' 
     >                        .OR. CLUMP_CRIT == 'RADIUS') THEN
        READ (ACTPAR, '(F20.0)', ERR=97) PAR1
      ENDIF
      CALL SARGV (DENSCON_LINE, 5, ACTPAR) 
      IF (CLUMP_CRIT == 'VELO' .AND. ACTPAR == 'SONIC') THEN
C***    Clumping ends at sonic point
        PAR2 = Vsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .AND. ACTPAR == 'SONIC') THEN
        CALL SPLINPOX(PAR2,Rsonic,TAUROSS,RADIUS,ND)
      ELSEIF (CLUMP_CRIT == 'RADIUS' .AND. ACTPAR == 'SONIC') THEN
        PAR2 = Rsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .OR. CLUMP_CRIT == 'VELO' 
     >                           .OR. CLUMP_CRIT == 'RADIUS') THEN
        READ (ACTPAR, '(F20.0)', ERR=97) PAR2
      ENDIF

C***  Depth dependent Clumping Factor
         
C***  Clumping increases from VDENS1 to VDENS2


C***  Branch: Clumping scales with Rosseland optical depth
      IF (CLUMP_CRIT .EQ. 'TAU') THEN
         X1 = AMIN1 (PAR1, PAR2)
         X1 = AMAX1 (TAUROSS(1), X1)
         X2 = AMAX1 (PAR1, PAR2)
         X2 = AMIN1 (TAUROSS(ND), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (TAUROSS(L) .LE. X1) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE IF (TAUROSS(L) .GE. X2) THEN
               DENSCON(L) = 1.
            ELSE
               X = PI * (TAUROSS(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = (1.-Q) + Q * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO

C***  Branch: Clumping scales with velocity
      ELSE IF (CLUMP_CRIT .EQ. 'VELO') THEN
         X1 = AMIN1 (PAR1, PAR2)
         X1 = AMAX1 (VELO(ND), X1)
         X2 = AMAX1 (PAR1, PAR2)
         X2 = AMIN1 (VELO(1), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (VELO(L) .LE. X1) THEN
               DENSCON(L) = 1.
            ELSE IF  (VELO(L) .GE. X2) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE
               X = PI * (VELO(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = Q + (1.-Q) * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO

      ELSE IF (CLUMP_CRIT == 'RADIUS') THEN
         X1 = MIN(PAR1, PAR2)
         X1 = MAX(RADIUS(ND), X1)
         X2 = MAX(PAR1, PAR2)
         X2 = MIN(RADIUS(1), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (RADIUS(L) <= X1) THEN
               DENSCON(L) = 1.
            ELSE IF  (RADIUS(L) >= X2) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE
               X = PI * (RADIUS(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = Q + (1.-Q) * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO
         
C***  The LIST criterion allows multiple clumping steps
C***  which can be defined by velocity, optical depth or radius
C***  The first value before the LIST keyword defines the outermost clumping factor
C***  The keyword after LIST defines the criterion for the following list which
C***  can be of arbitrary length, but must follow in (parval dval)-pairs.
C***   
C***  SYNTAX:  DENSCON [Dout] LIST TAU|VELO|RADIUS p1 D1 p2 D2 p3 D3 ...
C***  Example: DENSCON 10 LIST TAU 0.1 1. 0.005 20. 
C***       means that (seen from R_star outwarts)
C***           - clumping starts with D=1. until TAUROSS_cont drops below 0.1
C***           - increases to D=20. until TAUROSS_cont drops below 0.005
C***           - decreases to D=10 in the outermost wind (TAUROSS_cont < 0.005)
      ELSE IF (CLUMP_CRIT == 'LIST') THEN
        DO L=1, ND
          DENSCON(L) = DENSCON_FIX
        ENDDO
        CALL SARGV (DENSCON_LINE, 4, CLUMP_CRIT2) 
        IF (CLUMP_CRIT2 == 'VELO') THEN
          PARL = 0.
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            IF (ACTPAR == 'SONIC') THEN
              PAR1 = Vsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            ENDIF
            IF (PAR1 < PARL) THEN
              WRITE (hCPR,*) 'ERROR: ' 
              WRITE (hCPR,*) 'Velocity steps not in monotonic order'
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop: DO L=NDSTART, 1, -1
              IF (VELO(L) <= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop
              ENDIF
            ENDDO ndloop
          ENDDO
        ELSEIF (CLUMP_CRIT2 == 'TAU') THEN
          PARL = TAUROSS(ND) * 1.1
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            IF (PAR1 > PARL) THEN
              WRITE (hCPR,'(A)') ' ERROR: ' 
              WRITE (hCPR,'(A)') ' Tau steps not in decreasing order'
              WRITE (hCPR,'(2(A,F8.3))') 'Last = ',PARL, ' Next = ',PAR1
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop2: DO L=NDSTART, 1, -1
              IF (TAUROSS(L) >= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop2
              ENDIF
            ENDDO ndloop2
          ENDDO
        ELSEIF (CLUMP_CRIT2 == 'RADIUS') THEN
          PARL = 0.
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            IF (ACTPAR == 'SONIC') THEN
              PAR1 = Rsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            ENDIF
            IF (PAR1 < PARL) THEN
              WRITE (hCPR,*) 'ERROR: ' 
              WRITE (hCPR,*) 'Radius steps not in monotonic order'
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop3: DO L=NDSTART, 1, -1
              IF (RADIUS(L) <= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop3
              ENDIF
            ENDDO ndloop3
          ENDDO
        ELSE
          WRITE (hCPR,*) 'ERROR: ' 
          WRITE (hCPR,*)
     >      'Invalid choice of clumping-structure criterion: ', 
     >       CLUMP_CRIT2
          GOTO 99
        ENDIF
        !Smoothing of clumping grid
        DCSMOOTH(ND) = DENSCON(ND)
        DCSMOOTH(1) = DENSCON(1)
        DO L=ND-1, 2, -1
          DCSMOOTH(L) = DENSCON(L) * 0.6 +
     >          DENSCON(L+1) * 0.2 + DENSCON(L-1) * 0.2
        ENDDO
        DO L=1, ND
          DENSCON(L) = DCSMOOTH(L)
          FILLFAC(L) = 1. / DENSCON(L)
        ENDDO
      ELSE
C***  Error: invalid choice of Clumping Criterion
         GOTO 96
      ENDIF
      
  100 RETURN


C***  ERROR branches *********************************************

  93  WRITE (0,*) 'ERROR: '
      WRITE (0,*) 'Najarro formula needs three additional parameters'
      GOTO 99

  94  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Hillier formula needs one parameter: V_cl ' 
      GOTO 99

  95  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Hillier/Najarro velocity parameters must be > 0'
      GOTO 99

  96  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Invalid choice of Clumping-structure criterion: ', 
     >            CLUMP_CRIT       
      GOTO 99

  97  WRITE (0,*) 'ERROR: ',  ACTPAR, 
     >      ' could not be decodes as floating-point number!'
      GOTO 99    

  98  WRITE (0,*) 'ERROR: Clumping-structure criterion needs two', 
     >            'parameters'
      GOTO 99

  99  WRITE (0,*) 'The error occured in the following input line:'
      WRITE (0,*) DENSCON_LINE
      STOP 'FATAL ERROR detected by SUBROUTINE CLUMP_STRUCT'

      END


      SUBROUTINE CMSSTORE (ICHANNEL, IADRDUMMY, MAXADRDUMMY, 
     >                     NAME, NAME2, X, NDIM, ACTION, IERR)
C*******************************************************************
C***
C***  Cray-MS-STORagE
C***
C***  ROUTINE IS THE ADAPTER BETWEEN THE ROUTINE STORAGE, WHICH
C***    EMULATES MASS-STORAGE, AND THE MS-CALL AT CRAY
C***  IT IS NEEDED BECAUSE AT THE CRAY THE INDEX ARRAYS ARE NOT
C***    GIVEN BY THE CALL OF READMS, WRITMS AND CLOSMS
C***  Version 1.0 :
C***                The actual file is immedeatly closed when a second
C***                file is used
C***  Version 2.0 :
C***                unchanged
C***  Version 3.0 :
C***                The Array IADR is now twodimensional. The actual
C***                file need not to be closed when a new file is used
C***  Version 3.1 :
C***                MAXADR increased from 25600 to 256000 (2000 Records)
C***                wrh 15-Mar-2005 16:33:59
C***  Version 3.2 :
C***                MSMAXCH increased from 3 to 5
C*******************************************************************

      PARAMETER (MSMAXCH = 5)
C***  NOTE: MAXADR MUST BE A MULTIPLE OF 128
C***        (I.E. IADRL IN ROUTINE STORAGE)
      PARAMETER (MAXADR = 2000 * 128)

      CHARACTER*8 STATUS, ACTION
      DIMENSION IFILE(MSMAXCH), IADR(MAXADR,MSMAXCH)

C      SAVE LASTCH, IFILE, NFILE, STATUS, ITRANS
      SAVE

C***  CHECK THE ICHANNEL NUMBER
      IF (ICHANNEL .LE. 0) THEN
        WRITE (0,*) ' NEGATIV ICHANNEL NUMBERS ARE NOT ALLOWED'
        STOP 'ERROR IN CMSSTORE'
      ENDIF

      IF (ICHANNEL .NE. LASTCH) THEN
C***  SEARCH FOR FILE NUMBER
        STATUS = 'UNKNOWN'
        DO I=1, NFILE
          IF (ICHANNEL .EQ. IFILE(I)) THEN
            ITRANS = I
            STATUS = 'KNOWN'
            GOTO 10
          ENDIF
        ENDDO
   10   CONTINUE
C***  CLOSE THE FILE USED LAST. This is now (Vers 3) done by SCLOSE
        IF (LASTCH .GT. 0) THEN
          CALL STORAGE (LASTCH, IADR(1,ITRANSLAST), 
     >                  MAXADR, 'DUMMY', 'DUMMY', X, NDIM, 
     >                  'SCLOSE', 'CRAY', IERR)
        ENDIF
C***  OPEN THE ACTUAL FILE. This is now (Vers 3) done by SOPEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), 
     >                  MAXADR, 'DUMMY', 'DUMMY', 
     >                  X, NDIM, 'SOPEN', 'CRAY', IERR)
        ENDIF
      ENDIF

      IF (ACTION(1:4) .EQ. 'OPEN') THEN
        IF (STATUS .NE. 'UNKNOWN') THEN
          WRITE (0,*) 'DO NOT OPEN AN OPEN FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF
        NFILE = NFILE + 1
        ITRANS = NFILE
        IF (NFILE .GT. MSMAXCH) THEN
          WRITE (0,*) 'DO NOT OPEN MORE FILES THAN MSMAXCH'
          STOP 'ERROR IN CMSSTORE'
        ENDIF
        CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                NAME, NAME2, X, NDIM, 
     >                'OPEN', 'CRAY', IERR)
        IFILE(NFILE) = ICHANNEL
        STATUS = 'KNOWN'

      ELSE IF (ACTION(1:4) .EQ. 'READ') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'READ', 'CRAY', IERR)
        ELSE
          WRITE (0,'(A,I3)') 
     >        'DO NOT READ IN A CLOSED FILE, CHANNEL=', ICHANNEL
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:6) .EQ. 'LENGTH') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'LENGTH', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT READ (LENGTH) IN A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:5) .EQ. 'WRITE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'WRITE', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT WRITE IN A CLOSED FILE'
          WRITE (0,*) 'ICHANNEL=',ICHANNEL
          WRITE (0,'(A5,A8)') 'NAME=',NAME
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:5) .EQ. 'CLOSE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'CLOSE', 'CRAY', IERR)
          STATUS = 'UNKNOWN'
          NFILE = NFILE - 1
          DO I=ITRANS, NFILE
            IFILE(I) = IFILE(I+1)
            DO J=1, MAXADR
              IADR(J,I) = IADR(J,I+1)
            ENDDO
          ENDDO
        ELSE
          WRITE (0,*) 'DO NOT CLOSE A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:6) .EQ. 'CHANGE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'CHANGE', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT CHANGE IN A CLOSED FILE'
          WRITE (0,*) 'ICHANNEL=',ICHANNEL
          WRITE (0,'(A5,A8)') 'NAME=',NAME
          WRITE (0,'(A5,A8)') 'NAME2=',NAME2
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:4) .EQ. 'INFO') THEN
        IF (STATUS .EQ. 'OPEN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  ACTION, 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT DO INFO ON A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE
        WRITE (0,*) ' ACTION ', ACTION( :IDX(ACTION)), ' NOT KNOWN'
        STOP 'ERROR IN CMSSTORE'

      ENDIF

      IF (ACTION(1:5) .NE. 'CLOSE') THEN
        LASTCH = ICHANNEL
        ITRANSLAST = ITRANS
      ELSE
        LASTCH = 0
      ENDIF

      RETURN
      END
      SUBROUTINE COOP (XLAM,ND,T,RNE,POPNUM,POPMIN,ENTOT,RSTAR,
     $                 OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $                 NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,
     $                 EINST,ALPHA,SEXPO,
     $                 ADDCON1, ADDCON2, ADDCON3, 
     $                 IGAUNT,SIGMATHK,SEXPOK,EDGEK,K,NF,SIGMAKI, 
     $                 RADIUS,KONTNUP,KONTLOW,LASTKON,XDATA)
C***********************************************************************
C***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  CALLED FROM VARIOUS PLACES: 
C***                  WRCONT
C***                  COMO
C***                  FORMAL
C***                  FORMAL - FORMCMF
C***                  STEAL - LINPOP - CCORE
C***                  STEAL - LINPOP - OPAROSS
C***                  WRCONT - DIFDTDR - OPAROSS
C***                  COMO - DIFDTDR - OPAROSS
C***                  COLI - DIFDTDR - OPAROSS
C***                  STEAL - TAUSCAL (2x) - OPAROSS
C***                  STEAL - ENSURETAUMAX - TAUSCAL (2x) - OPAROSS
C***  SIMILAR SUBROUTINES: CMFCOOP, COOPFRQ, DCOOP
C***  This version (23-Mar-2007) assumes that KODAT positions
C***  (i.e. KODATIND) give the atomic number (NCORECHARGE)
C***********************************************************************
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR
      PARAMETER ( MAXION = 27 )
C***  MAXIMUM X-RAY DATA
      PARAMETER ( MAXXDAT = 10)
C***  Dimension of the core-charge data locally provided here
      PARAMETER (MAXATOMDIM = 30)

      DIMENSION XDATA(MAXXDAT)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION NOM(N)
      DIMENSION KODAT(MAXATOM)
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION POPNUM(ND,N)
      DIMENSION OPA(ND),ETA(ND),THOMSON(ND)
      CHARACTER(8), DIMENSION(ND) :: IWARN
      DIMENSION RADIUS(ND)
      DIMENSION T(ND),RNE(ND),ENTOT(ND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION SIGMAKI(NF,LASTKON)
      DIMENSION GFF(0:MAXION)
      DIMENSION KODATIND(MAXATOMDIM)
      CHARACTER*10 LEVEL(N),MAINPRO(ND),MAINLEV(ND)
      LOGICAL XRAYS, KSHELL
      
      REAL, PARAMETER :: EXPMIN = -300.

C***  Output of laser warnings for bound-free transitions
      INTEGER, SAVE :: NWARN
      DATA NWARN /0/ ! no warning has been issued yet

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
      DATA C2 / 3.9724E-16 /
C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /
C***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON  ( CM**3 * KELVIN**(3/2) )
      DATA C3 / 2.07E-16 /
C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100, CM**5 )
      DATA CFF / 1.370E-23 /

      W  = 1.E8/XLAM
      W3 = W*W*W
 
C***  K-SHELL ABSORPTION CROSS-SECTIONS PROVIDED BY DATOM FILE ?
      KSHELL = .FALSE.
      DO NA=1, MAXATOM
        DO ISTATE=1, MAXATOM
         IF (SIGMATHK(NA,ISTATE) .NE. .0) THEN
            KSHELL = .TRUE.
            EXIT
         ENDIF
        ENDDO
      ENDDO

C***  X-RAY EMISSION SWITCHED ON BY CARD OPTION ?
      XRAYS = .FALSE.
      IF (XDATA(1) .NE. 0.) XRAYS = .TRUE.

C***  ONLY NEEDED FOR K-SHELL OR XRAY BRANCH
      IF (XRAYS .OR. KSHELL) THEN 
C***     Establish KODAT index for each used element
C***     First find number of used elements
         NATOMMAX = NOM(N)
         IF (MAXATOM .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
            STOP 'ERROR IN COOP'
         ENDIF

C***     Now find for each NA the corresponding KODAT index
         DO NA=1, NATOMMAX
            KODATIND(NA) = 0
            DO J = 1, MAXATOM
               IF (NA .EQ. KODAT(J)) KODATIND(NA) = J 
            ENDDO
            IF (KODATIND(NA) .EQ. 0) THEN
               WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
               STOP 'ERROR IN COOP'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOP'
            ENDIF
         ENDDO
      ENDIF

C***  PARAMETER FOR X-RAY SOURCE, AND SOME PREPARATIONS
      IF (XRAYS) THEN
         XFILL = XDATA(1)
         XRAYT = XDATA(2)
         XMINR = XDATA(3)
         DIFFEMEXP = XDATA(4)
         EXPVAL = -C1*W/XRAYT
         IF (EXPVAL >= EXPMIN) THEN           
           EXPFACXRAY = EXP(EXPVAL)
         ELSE 
           EXPFACXRAY = 0.
         ENDIF
         PRESIGXRAY = CFF / W3 / SQRT(XRAYT)
         IF (XDATA(5) .NE. 0.) THEN
            XFILL2 = XDATA(5)
            XRAYT2 = XDATA(6)
            XMINR2 = XDATA(7)
            EXPVAL = -C1*W/XRAYT2
            IF (EXPVAL >= EXPMIN) THEN           
              EXPFACXRAY2 = EXP(EXPVAL)
            ELSE 
              EXPFACXRAY2 = 0.
            ENDIF
            PRESIGXRAY2 = CFF / W3 / SQRT(XRAYT2)
         ENDIF
C***    Calculate number of free electrons, assuming full ionization
C***    Because of number conservation, this can be done for any depth point
         RNEXRAY = .0
         DO J = 1, N
          RNEXRAY = RNEXRAY + POPNUM(1,J) * KODATIND(NOM(J)) 
         ENDDO
      ENDIF
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      NBFLASER = 0
   55 CONTINUE
      OPAMAX=.0
      OPAL=.0
      ETAL=.0
      IWARN(L)='        '
      TL=T(L)
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL
 
C***  BOUND-FREE  ******************************************************
C***  I = LOW      J = UP
      DO 5 KON=1,LASTKON
      J=KONTNUP(KON)
      I=KONTLOW(KON)
      EDGE=ELEVEL(J)+EION(I)-ELEVEL(I)
      IF (W .LT. EDGE) GOTO 5
 
C***  CALCULATE SIGMA, THE FREQUENCY-DEPENDENT CROSS SECTION
C***  IF ( K .GT. 0 ) IT IS ASSUMED THAT THE BOUND-FREE CROSS SECTIONS SIGMA
C***  HAVE BEEN ALREADY CALCULATED ( ARRAY SIGMAKI )
      IF (K .GT. 0) THEN
            SIGMA=SIGMAKI(K,KON)
            ELSE
            SIGMATH=EINST(I,J)*1.E-18
            CALL PHOTOCS (SIGMA,SIGMATH,EDGE,W,ALPHA,SEXPO,
     >                    ADDCON1, ADDCON2, ADDCON3, 
     >                    IGAUNT,KON)
            ENDIF
            IF (SIGMA < 0.) THEN
              WRITE (0,*) 'BIGWARN: negative BF cross section!!!'
            ENDIF
            
 
C***  RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      WE=C3*RNE(L)*ENTOT(L)/T32
      EXPVAL = C1*(EDGE-W)/TL
      IF (EXPVAL >= EXPMIN) THEN
        G=WEIGHT(I)/WEIGHT(J)*WE*EXP(EXPVAL)
      ELSE 
        G=0.
      ENDIF
      EMINDU=G*POPNUM(L,J)*SIGMA

C***  Set emissivities zero if both levels are equal (=POPMIN)
      IF (POPNUM(L,I) .EQ. POPNUM(L,J)) EMINDU = .0

      SUM=POPNUM(L,I)*SIGMA-EMINDU
C***  LASER WARNING IF STIMULATED EMISSION EXCEEDS ABSORPTION IN THIS TRANSITION
C***  IF TOTAL CONT. OPA WAS < 0 AT FIRST TRIAL. THIS BF TRANSITION IS SKIPPED
cc      IF (SUM.LT. .0 .AND. NBFLASER .EQ. 1) THEN
      IF (SUM.LT. .0) THEN
         IWARN(L)='*       '
         IF (NWARN .EQ. 0)  THEN
            WRITE (0, 90) L, LEVEL(I), LEVEL(J)
   90       FORMAT ('*** WARNING FROM Subr. COOP: ',
     >       'LASERING BOUND-FREE CONTINUA SUPPRESSED',
     >       /,'*** THIS OCCURED FOR THE FIRST TIME AT DEPTH INDEX', I7, 
     >       /,'*** BETWEEN LEVELS ', A10, ' AND ', A10 )
             NWARN = NWARN + 1
         ENDIF
      ELSE
         OPAL=OPAL+SUM
         ETAL=ETAL+EMINDU
      ENDIF
      IF(SUM .LT. OPAMAX) GOTO 5
        OPAMAX=SUM
        MAINPRO(L)='BOUND-FREE'
        MAINLEV(L)=LEVEL(I)
    5 CONTINUE

C***  K-SHELL IONISATION  **********************************************
      IF (KSHELL) THEN      
C***     LOOP OVER ALL LEVELS 
         LASTNOMJ = -1
         LASTISTATE = -1
         DO 6 J=1,N
            NOMJ = NOM(J)
            ISTATE = NCHARG(J) + 1
C***        ARE THERE K-SHELL-DATA FOR CURRENT ELEMENT?
            IF (SIGMATHK(NOMJ,ISTATE) .EQ. 0.) GOTO 6

C***        IS RADIATION HARDER THAN K-SHELL EDGE ? 
            IF (W .LT. EDGEK(NOMJ,ISTATE)) GOTO 6

C***        K-SHELL IONIZATION NEEDS IONS WITH AT LEAST 3 ELECTRONS LEFT
            IF (KODATIND(NOMJ) - NCHARG(J) .LT. 3) THEN  
               WRITE (0,*) 'UNEXPECTED INCONSISTENCY WITH K-SHELL DATA'
               STOP 'ERROR in Subr. COOP'
            ENDIF

            IF (LASTNOMJ .NE. NOMJ .OR. LASTISTATE .NE. ISTATE) THEN
               CALL KSIGMA (SIGMAK, SIGMATHK(NOMJ,ISTATE), 
     >                   EDGEK(NOMJ,ISTATE), W, SEXPOK(NOMJ,ISTATE))
               LASTNOMJ = NOMJ
               LASTISTATE = ISTATE
            ENDIF

c            IF (POPNUM(L,J) < 1.1 * POPMIN) CYCLE
            SUM = POPNUM(L,J) * SIGMAK
            OPAL = OPAL + SUM
            IF (SUM .GT. OPAMAX) THEN
               OPAMAX=SUM
               MAINPRO(L) = 'K-SHELL'
               MAINLEV(L) = LEVEL(J)(:2)
               WRITE (MAINLEV(L)(4:5), '(I2)') ISTATE 
            ENDIF
    6    CONTINUE
      ENDIF

C***  X-RAY EMISSION -------------------------------------------------
      IF (XRAYS) THEN
C***     X-RAY SOURCE: FREE-FREE BREMSSTRAHLUNG 
         IF ((XFILL .GT. 0.) .AND. (RADIUS(L) .GE. XMINR)) THEN
            DO I=1, N
               NCHARI = KODATIND(NOM(I)) 
               SIGMAFF= PRESIGXRAY * FLOAT(NCHARI * NCHARI)
               SUM = RNEXRAY * ENTOT(L) * POPNUM(L,I) * SIGMAFF * XFILL
C***           For the differential emission measure option, 
C***           the exponential factor is replaced by a function 
C***           that depends on the DEM exponent - see documentation
               IF (DIFFEMEXP .EQ. .0) THEN
                  FDEM = EXPFACXRAY
               ELSE IF (DIFFEMEXP .EQ. 1.5) THEN
                  XXMAX = C1*W/XRAYT 
                  XXMIN = C1*W/TL
                  FDEM  = (EXPFACXRAY - EXP(-XXMIN)) / XXMAX
C*                normalization
                  FDEM = FDEM * 0.5 / ((XRAYT/TL)**0.5 - 1.)
               ELSE IF (DIFFEMEXP .EQ. 2.5) THEN
                  XXMAX = C1*W/XRAYT 
                  XXMIN = C1*W/TL
                  FDEM  = EXPFACXRAY - EXP(-XXMIN)
     >                  + XXMAX * EXPFACXRAY - XXMIN * EXP(-XXMIN)
                  FDEM = FDEM / (XXMAX*XXMAX)
C*                normalization
                  FDEM = FDEM * 1.5 / ((XRAYT/TL)**1.5 - 1.)
               ELSE
                  WRITE (0,*) '*** XRAY DIFFERENTIAL EMISSION MEASURE:' 
                  WRITE (0,*) '*** INVALID EXPONENT:', DIFFEMEXP
                  STOP '*** FATAL ERROR IN SUBR. COOP'
               ENDIF
               EMINDU = SUM * FDEM
               SUM = SUM - EMINDU
               OPAL = OPAL + SUM
               ETAL = ETAL + EMINDU 
               IF (SUM .GE. OPAMAX) THEN
                  OPAMAX=SUM
                  MAINPRO(L)='XRAYSOURCE'
                  MAINLEV(L)=LEVEL(I)
               ENDIF
            ENDDO
         ENDIF
         IF ((XFILL2 .GT. 0.) .AND. (RADIUS(L) .GE. XMINR2)) THEN
            DO I=1, N
               NCHARI = KODATIND(NOM(I)) 
               SIGMAFF= PRESIGXRAY2 * FLOAT(NCHARI * NCHARI)
               SUM = RNEXRAY * ENTOT(L) * POPNUM(L,I) * SIGMAFF *XFILL2
               EMINDU = SUM * EXPFACXRAY2
               SUM = SUM - EMINDU
               OPAL = OPAL + SUM
               ETAL = ETAL + EMINDU 
               IF (SUM .GE. OPAMAX) THEN
                  OPAMAX=SUM
                  MAINPRO(L)='XRAYSOURCE'
                  MAINLEV(L)=LEVEL(I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
 
C***  FREE-FREE  *******************************************************
C***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      GFF(0)=.0
      DO 10 ION=1, MAXION
      CALL GAUNTFF (GIII,ION,XLAM,TL)
      GFF(ION)=GIII*FLOAT(ION*ION)
   10 CONTINUE

C***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
      PRESIG=CFF/W3/ROOTTL
      EXPVAL = -C1*W/TL
      IF (EXPVAL >= EXPMIN) THEN
        EXPFAC=EXP(EXPVAL)
      ELSE 
        EXPFAC = 0.
      ENDIF
      DO 3 I=1,N
      NCHARI=NCHARG(I)
      IF (NCHARI .GT. MAXION) THEN
         WRITE (0,*) '*** ERROR COOP: MAXION TOO SMALL'
         STOP '*** ERROR IN COOP'
      ENDIF
      SIGMAFF=PRESIG*GFF(NCHARI)
c      IF (POPNUM(L,I) < 1.1*POPMIN) GOTO 3
      SUM=RNE(L)*ENTOT(L)*POPNUM(L,I)*SIGMAFF
      EMINDU=SUM*EXPFAC
      SUM=SUM-EMINDU
      OPAL=OPAL+SUM
      ETAL=ETAL+EMINDU
      IF (SUM.LT.OPAMAX) GOTO 3
      OPAMAX=SUM
      MAINPRO(L)='FREE-FREE'
      MAINLEV(L)=LEVEL(I)
    3 CONTINUE

C***  If total true continuum opacity is negative: 
C***  re-do the whole calculation, but skip lasering bound-free transitions 
C***  Note: in the version with strict suppression (wrh  4-Apr-2019)
C***        this condition should never be met, and is therefore
C***        commented
cc      IF (OPAL .LE. .0 .AND. NBFLASER .EQ. 0) THEN
cc         NBFLASER = 1 
cc         GOTO 55
cc      ENDIF

C***  ATENTION NOTICE: LINE IS IMPORTANT FOR ALL ETAL MODIFIED BEFORE 
      ETAL=ETAL*C2*W3
 
C***  THOMSON SCATTERING ***********************************************
      SUM=RNE(L)*SIGMAE
      IF (SUM.LT.OPAMAX) GOTO 4
      MAINPRO(L)='THOMSON'
      MAINLEV(L)='ELECTRON'
    4 OPAL=OPAL+SUM
C***  THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
      THOMSON(L)=SUM/OPAL
 
      OPA(L)=OPAL*ENTOT(L)*RSTAR
      ETA(L)=ETAL*ENTOT(L)*RSTAR
    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
      SUBROUTINE COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $                   WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $                   RNEL,TL,SIGMAFF,MAXION,RL,XDATA,
     $                   SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $                   KONTNUP,KONTLOW,LASTKON,OPATHOM)
C***********************************************************************
C***  NON-LTE CONTINUOUS OPACITY AT CURRENT DEPTH FOR ALL FREQUENCIES
C***  NOTE: ONLY TRUE OPACITY, WITHOUT THOMSON SCATTERING TERM.
C***  THOMSON OPACITY IS PREPARED AS AN EXTRA VARIABLE: OPATHOM
C***  This version (23-Mar-2007) assumes that KODAT positions
C***  (i.e. KODATIND) give the atomic number (NCORECHARGE)
C***  Called from: STEAL --> LINPOP --> COMA
C***          and: WRSTART --> GREY --> OPAGREY
C***********************************************************************
 
      DIMENSION NOM(N)
      REAL, DIMENSION(MAXATOM, MAXION) :: SIGMATHK, SEXPOK, EDGEK
      DIMENSION XDATA(10)
      DIMENSION KODAT(MAXATOM)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N),EN(N)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),NFEDGE(LASTKON)
      DIMENSION OPAC(NF),ETAC(NF),XLAMBDA(NF),EXPFAC(NF)
      DIMENSION SIGMAKI(NF,LASTKON),SIGMAFF(NF,0:MAXION)
      LOGICAL XRAYS, KSHELL
C***  Dimension of the core-charge data locally provided here
      PARAMETER ( MAXATOMDIM = 30)
      DIMENSION KODATIND(MAXATOMDIM)

C***  Output of laser warnings for bound-free transitions
      DATA NWARN /0/ ! no warning has been issued yet
      SAVE NWARN
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
      DATA C2 / 3.9724E-16 /
C***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      DATA C3 / 2.07E-16 /
C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION (ALLEN PAGE 100)
      DATA CFF / 1.370E-23 /
C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /


C***  K-SHELL ABSORPTION CROSS-SECTIONS PROVIDED BY DATOM FILE ?
      KSHELL = .FALSE.
      DO NA=1, MAXATOM
        DO ISTATE=1, MAXATOM
         IF (SIGMATHK(NA,ISTATE) .NE. .0) THEN
            KSHELL = .TRUE.
            EXIT
         ENDIF
        ENDDO
      ENDDO

C***  PARAMETER FOR X-RAY SOURCE
      XRAYS = .FALSE.
      IF (XDATA(1) .NE. 0.) THEN
        XRAYS = .TRUE.
        XFILL = XDATA(1)
        XRAYT = XDATA(2)
        XMINR = XDATA(3)
        DIFFEMEXP = XDATA(4)
        IF (XDATA(5) .NE. 0.) THEN
           XFILL2 = XDATA(5)
           XRAYT2 = XDATA(6)
           XMINR2 = XDATA(7)
        ENDIF
      ENDIF

      T32=TL*SQRT(TL)
 
ccC***  Safety against negative opacities from lasering bound-free trans.
ccC***  Note: In contrast to the oher opacity routines, in COOPFRQ
ccC***   all lasering b-f continue are suppressed *at all frequencies* 
ccC**    if a negative "true" opacity is encountered at any frequency 
cc      NBFLASER = 0
cc   55 CONTINUE

      DO 10 K=1,NF
      OPAC(K)=.0
      ETAC(K)=.0
   10 CONTINUE
 
C***  BOUND-FREE  ******************************************************

C***  LOOP OVER ALL CONTINUUM TRANSITIONS
      DO 5 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON) 
      EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
      IF (C1*EDGE/TL .GT. 700.) THEN
        WE = 0.
      ELSE
        EXPEDGE=EXP(C1*EDGE/TL)
        WE=C3*RNEL*ENTOTL/T32 *WEIGHT(LOW)/WEIGHT(NUP)*EXPEDGE
      ENDIF
      NFLOW=NFEDGE(KON)

C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS WITH XLAMBDA(K) < EDGE
      DO 11 K=1,NFLOW
      SIGMA=SIGMAKI(K,KON)
C***  RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      G=WE*EXPFAC(K)
      EMINDU=G * SIGMA * EN(NUP)

C***  Set emissivities zero if both levels are equal (=POPMIN)
      IF (EN(LOW) .EQ. EN(NUP)) EMINDU = .0

      SUM=    EN(LOW)*SIGMA-EMINDU

C***  LASER bf continua are skipped!  3-Feb-2016
cc      IF (SUM .LT. .0 .AND. NBFLASER .EQ. 1) THEN
C***    Note: the above version, which makes use of NBLASER. has been
C***          replaced by the more radical version that any lasering b-f
C***          transition is immediately disregarded - wrh  4-Apr-2019

cc      IF (SUM .LT. .0) THEN
cc         IF (NWARN .EQ. 0)   WRITE (0, 90) 
cc   90    FORMAT ('*** WARNING FROM Subr. COOPFRQ: ',
cc     >   'LASERING BOUND-FREE CONTINUA SUPPRESSED')
cc         NWARN = NWARN + 1
cc      ELSE
      IF (SUM .GT. .0) THEN
         OPAC(K)=OPAC(K)+SUM
         ETAC(K)=ETAC(K)+EMINDU
      ENDIF
   11 CONTINUE
 
    5 CONTINUE

C***  ONLY NEEDED FOR K-SHELL OR XRAY BRANCH
      IF (XRAYS .OR. KSHELL) THEN 
C***     Establish KODAT index for each used element
C***     First find number of used elements
         NATOMMAX = NOM(N)
         IF (MAXATOM .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
            STOP 'ERROR IN COOPFRQ'
         ENDIF
C***     Now find for each NA the corresponding KODAT index
C***     in order to know the core charge
         DO NA=1, NATOMMAX
            KODATIND(NA) = 0
            DO II = 1, MAXATOM
               IF (NA .EQ. KODAT(II)) KODATIND(NA) = II
            ENDDO
            IF (KODATIND(NA) .EQ. 0) THEN
               WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
         ENDDO
      ENDIF

C***  K-SHELL IONISATION  **********************************************
      IF (KSHELL) THEN
C***  LOOP OVER ALL LEVELS 
      DO 6 J=1,N
         NOMJ = NOM(J)
         ISTATE = NCHARG(J) + 1
C***     ARE THERE K-SHELL-DATA FOR CURRENT ELEMENT?
         IF (SIGMATHK(NOMJ,ISTATE) .EQ. 0.) GOTO 6

C***     K-SHELL IONIZATION NEEDS IONS WITH AT LEAST 3 ELECTRONS LEFT
         IF (KODATIND(NOMJ) - NCHARG(J) .LT. 3) THEN  
            WRITE (0,*) 'UNEXPECTED INCONSISTENCY WITH K-SHELL DATA'
            STOP 'ERROR in Subr. COOPFRQ'
         ENDIF

         WK = 1.E8 / EDGEK(NOMJ,ISTATE)
         DO K=1,NF
C***        IS RADIATION TOO SOFT FOR K-SHELL-IONISATION? 
            IF (XLAMBDA(K) .GT. WK) EXIT
            W = 1.E8 / XLAMBDA(K)
            CALL KSIGMA (SIGMAK, SIGMATHK(NOMJ,ISTATE), 
     >                   EDGEK(NOMJ,ISTATE), W, SEXPOK(NOMJ,ISTATE))
  
            SUM = EN(J) * SIGMAK
            OPAC(K) = OPAC(K) + SUM
         ENDDO
    6 CONTINUE
      ENDIF

C***  X-RAY SOURCE: FREE-FREE BREMSSTRAHLUNG HeIII ASSUMED
      IF (XRAYS) THEN
      IF ((XFILL .GT. 0.) .AND. (RL .GE. XMINR)) THEN
        DO 15,K=1,NF
          W = 1.E8/XLAMBDA(K)
          W3 = W * W * W
          EXPFACXRAY = EXP(-C1*W/XRAYT)
          SUM    = 2. * ENTOTL * 4. * CFF / W3 / SQRT(XRAYT)
C***           For the differential emission measure option,   
C***           the exponential factor is replaced by a function  
C***           that depends on the DEM exponent - see documentation  
               IF (DIFFEMEXP .EQ. .0) THEN
                  FDEM = EXPFACXRAY
               ELSE IF (DIFFEMEXP .EQ. 1.5) THEN
                  XXMAX = C1*W/XRAYT
                  XXMIN = C1*W/TL
                  FDEM  = (EXPFACXRAY - EXP(-XXMIN)) / XXMAX
C*                normalization             
                  FDEM = FDEM * 0.5 / ((XRAYT/TL)**0.5 - 1.)
               ELSE IF (DIFFEMEXP .EQ. 2.5) THEN
                  XXMAX = C1*W/XRAYT
                  XXMIN = C1*W/TL
                  FDEM  = EXPFACXRAY - EXP(-XXMIN)
     >                  + XXMAX * EXPFACXRAY - XXMIN * EXP(-XXMIN)
                  FDEM = FDEM / (XXMAX*XXMAX)
C*                normalization    
                  FDEM = FDEM * 1.5 / ((XRAYT/TL)**1.5 - 1.)
               ELSE
                  WRITE (0,*) '*** XRAY DIFFERENTIAL EMISSION MEASURE:'
                  WRITE (0,*) '*** INVALID EXPONENT:', DIFFEMEXP
                  STOP '*** FATAL ERROR IN SUBR. COOPFRQ'
               ENDIF
               EMINDU = SUM * FDEM
C*             In the DEM branch, the opacity is set to zero 
               IF (DIFFEMEXP .EQ. .0) THEN
                 OPAX   = (SUM-EMINDU) * XFILL
                 OPAC(K)=OPAC(K)+OPAX
               ENDIF
          ETAX   = EMINDU * XFILL
          ETAC(K)=ETAC(K)+ETAX
   15   CONTINUE
      ENDIF
      IF ((XFILL2 .GT. 0.) .AND. (RL .GE. XMINR2)) THEN
        DO 16,K=1,NF
          W = 1.E8/XLAMBDA(K)
          W3 = W * W * W
          EFACTX = EXP(-C1*W/XRAYT2)
          SUM    = 2. * ENTOTL * 4. * CFF / W3 / SQRT(XRAYT2)
          EMINDU = SUM * EFACTX
          OPAX   = (SUM-EMINDU) * XFILL2
          ETAX   = EMINDU * XFILL2
          OPAC(K)=OPAC(K)+OPAX
          ETAC(K)=ETAC(K)+ETAX
 16     CONTINUE
      ENDIF
      ENDIF


C***  FREE-FREE  *******************************************************
C***  Loop 3 is most time consuming in WRSTART (53 percent
C***  Tested on 30-Jan-1997 19:52:38, Lars
      DO 30 K=1,NF
      W=1.E8/XLAMBDA(K)
      W3=W*W*W
      DO 3 I=1,N
      SUM=RNEL*ENTOTL*EN(I)*SIGMAFF(K,NCHARG(I))
      EMINDU=SUM*EXPFAC(K)
      SUM=SUM-EMINDU
      OPAC(K)=OPAC(K)+SUM
      ETAC(K)=ETAC(K)+EMINDU
    3 CONTINUE
      ETAC(K)=ETAC(K)*C2*W3*ENTOTL*RSTAR
      OPAC(K)=OPAC(K)*ENTOTL*RSTAR
   30 CONTINUE

ccC***  If total true continuum opacity is negative at at least one freq.:
ccC***  re-do the whole calculation, but skip lasering bound-free trans.
cc      OPACMIN = OPAC(1) 
cc      DO K=2, NF
cc         IF (OPAC(K) .LT. OPACMIN) OPACMIN = OPAC(K) 
cc      ENDDO
C***  Note: in the version with strict suppression (wrh  4-Apr-2019)
C***        this condition should never be met, and is therefore
C***        commented
cc      IF (OPACMIN .LE. .0 .AND. NBFLASER .EQ. 0) THEN
cc         NBFLASER = 1
cc         GOTO 55
cc      ENDIF

C***  THOMSON SCATTERING OPACITY (FREQUENCY INDEPENDENT! ) ************
      OPATHOM = RNEL * SIGMAE * ENTOTL * RSTAR
 
      RETURN
      END
      SUBROUTINE COUNT(J, ISTR)
C **********************************************************************
C *** Writes a 3-digit integer number with leading zeros (001, 002 ...)
C ***  CALLED BY FEDAT
C **********************************************************************

      CHARACTER*3 ISTR

ccc   test: restore the error of Leindecker's version!
ccc      if (j .eq. 100) then
ccc         istr = '099'
ccc         return
ccc      endif 

C *** COUNTER FOR LEVEL NAMES
      IF (J .LT. 0 .OR. J .GT. 999) THEN
         STOP 'FATAL ERROR IN SUBR. COUNT'
      ELSE IF (J .LT. 10) THEN
         WRITE (ISTR,'(A2,I1)') '00', J
      ELSE IF (J .LT. 100) THEN
         WRITE (ISTR,'(A1,I2)') '0', J
      ELSE
         WRITE (ISTR,'(I3)') J
      ENDIF
    
      RETURN
    
      END
      SUBROUTINE DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,SIGMATHK,
     $                  SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND, KEYCBF,
     >                  ROUTINE, INDEXMAX, NFEREADMAX, MAXFEIND,
     >                  LASTFE, SIGMAFE, INDRB, INDRF,
     >                  IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY, 
     >                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >                  LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

c!!!!!! Folgende parameter wurden entfernt:
C!!!    CBFC, BOUND, EINSTINT, COCOFE, NCOMAX, NCO
c!!!    folgende Parameter sind neu: MAXFEIND, FEDUMMY
C!!!    umbenannt wurden: NMAX -> NFEREADMAX, DUMMY -> FEDUMMY


C*******************************************************************************
C***  READS ATOMIC DATA FROM TAPE4=DATOM  **************************************
C***  The decoded elements are indexed in sequence of their occurrence
C***   (J = 1 ... NATOM)
C***  The meaning of vector KODAT is weird:
C***  Each chemical element is assigned a position in this vector
C***    (originally by an arbitrary definition; 
C***    in the version from March 2007, the index is equal to 
C***    the corecharge (NZ, atomic number, Kernladungszahl)
C***  KODAT(NZ) contains the index J under which this element was found
C***    in the DATOM file; unused elements have KODAT(NZ)=0 
C*******************************************************************************
 
      INTEGER, INTENT(IN) :: NDIM, MAXIND, MAXION, MAXKONT, MAXAUTO

      INTEGER, DIMENSION(NDIM) :: NCHARG, IONGRND, MAINQN, NOM
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,ADDCON1,ADDCON2,ADDCON3
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      DIMENSION LOWAUTO(MAXAUTO),WAUTO(MAXAUTO),EAUTO(MAXAUTO)
     $         ,AAUTO(MAXAUTO),IONAUTO(MAXAUTO),KRUDAUT(MAXAUTO)
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(NDIM),LEVUP,LEVLOW, LEVION
      CHARACTER*10 ELEMENT(MAXATOM),NEWELE
      CHARACTER*8 IGLOW, ICBF, IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      CHARACTER*4 CEY,KEYCBB(MAXIND)
      CHARACTER*3 KRUDI,DRRUDI
      CHARACTER*2 SYMBOL(MAXATOM), KSHELLSYM
      CHARACTER(LEN=*), INTENT(IN) :: ROUTINE

      LOGICAL :: BFEMODEL
 
      DO 15 NA=1,MAXATOM
         DO ISTAGE=1, MAXION
           SIGMATHK(NA,ISTAGE)=.0
           SEXPOK  (NA,ISTAGE)=.0
           EDGEK   (NA,ISTAGE)=.0
         ENDDO
   15 KODAT(NA)=0
      DO 6 I=1,NDIM
      EION(I)=.0
      IONGRND(I)=-1
      ALTESUM(1,I)=-1.
C***  INITIALIZE TRANSTION MATRIX TO DETECT MISSING LINE TRANSITIONS
      DO 6 J=1,NDIM
    6 EINST(I,J)=-99.
    
      DO IND=1,MAXIND
       INDNUP(IND)=0
       INDLOW(IND)=0
       COCO(1,IND)=.0
       COCO(2,IND)=.0
       COCO(3,IND)=.0
       COCO(4,IND)=.0
       KEYCBB(IND)='    '
      ENDDO

      DO KON=1,MAXKONT
       KONTNUP(KON)=-1
       KONTLOW(KON)=0
       IGAUNT(KON)= ' '
      ENDDO

      DO 96 J=1,MAXAUTO
      WAUTO(J)=0.0
      EAUTO(J)=0.0
      AAUTO(J)=0.0
      IONAUTO(J)=-1
      KRUDAUT(J)=0
   96 LOWAUTO(J)=0
      NATOM=0
      N=0
      IND=0
      KONT=0
      NAUTO=0
      LEVSEQ=0

      BFEMODEL = .FALSE.
 
      OPEN (4, FILE='DATOM', STATUS='OLD')
    1 READ(4,2,END=3) KARTE
    2 FORMAT(A)

      IF (KARTE(:1) .EQ. '*' .OR. KARTE(:1) .EQ. ' ') GOTO 1
      IF (KARTE(:10) .EQ. 'ELEMENT   ' ) GOTO 5
      IF (KARTE(:10) .EQ. 'LEVEL     ' ) GOTO 10
      IF (KARTE(:10) .EQ. 'LINE      ' ) GOTO 20
      IF (KARTE(:10) .EQ. 'CONTINUUM ' ) GOTO 30
      IF (KARTE(:10) .EQ. 'K-SHELL   ' ) GOTO 40
      IF (KARTE(:10) .EQ. 'LTESUM    ' ) GOTO 50
      IF (KARTE(:10) .EQ. 'DRTRANSIT ' ) GOTO 60
      CALL REMARK ('UNRECOGNIZED DATA INPUT')
      GOTO 990
      
C***  ELEMENTS ---------------------------------------------------------
    5 CONTINUE
C***  DECODED ELEMENT IS ALREADY KNOWN
      NEWELE=KARTE(13:22)
      DO 19 NA=1,NATOM
      IF (NEWELE .EQ. ELEMENT(NA)) GOTO 1
   19 CONTINUE

C***  If DATOM was called with parameter ROUTINE = 'NOIRON, 
C****   the 'ELEMENT GENERIC' card is ignored. 
C***    This feature was introduced to facilitate 
C***    Sonja's NEWFORMAL_CARDS program
      IF (ROUTINE(1:6) == 'NOIRON' .AND. NEWELE == 'GENERIC') GOTO 1

C***  NEW ELEMENT DECODED:
      LEVSEQ=0
      NATOM=NATOM+1
      IF (NATOM .GT. MAXATOM) THEN
         CALL REMARK ('DATOM: MORE ELEMENTS THAN DIMENSIONED')
         GOTO 990
      ENDIF
      READ (KARTE,9) ELEMENT(NATOM),SYMBOL(NATOM),ATMASS(NATOM),
     $                    STAGE(NATOM)
    9 FORMAT (12X,A10,2X,A2,4X,F6.2,3X,F5.0)

      CALL FINDCHARGE (ELEMENT(NATOM), NZ)
      KODAT(NZ) = NATOM      

      IF (NZ .EQ. 0) THEN
         WRITE (0,*) 'UNKNOWN ELEMENT DECODED: ', ELEMENT(NATOM)
         GOTO 990
      ENDIF

C***  "GENERIC" MODEL ATOM OF IRON GROUP ELEMENTS DECODED
      IF  (NZ .EQ. 26) THEN
         BFEMODEL = .TRUE.
C***    DECODE INPUT CARD AGAIN
 109     FORMAT (12X,A10,2X,A2,4X,I6,2X,I6)
         READ (KARTE,109) ELEMENT(NATOM),SYMBOL(NATOM),IONLOW,IONTOP
C***    COMPLETE IRON-DATA IS READ IN FROM MASS STORAGE FILE 'FEDAT'
         NFIRSTFE = N + 1

C***     Iron line indices are arranged *behind* the DRTRANSITs
         LASTINDAUTO = LASTIND + NAUTO

         CALL FEDAT (ROUTINE, INDEXMAX, NFEREADMAX, IONLOW, IONTOP,
     &               MAXATOM, NDIM, MAXIND, MAXKONT, NATOM,
     &               N, LASTFE, LASTKON, LASTINDAUTO, MAXFEIND,
     &               EINST, SIGMAFE, INDRB, INDRF, IFENUP, 
     &               IFELOW, INDNUP, INDLOW, KONTNUP, KONTLOW,
     &               LEVEL, ELEMENT, SYMBOL, ATMASS, STAGE,
     &               ELEVEL, WEIGHT, EION, NCHARG, NOM, KODAT,
     &               NFIRST, NLAST, IFRBSTA, IFRBEND, FEDUMMY,
     &               VDOPFE, DXFE, XLAM0FE, SIGMAINT, KEYCBB)
C***  Fill vector with ionization energies
         DO I=NFIRSTFE+1, N-1
            IF (EION(I) .EQ. .0) EION(I) = EION(I-1)
         ENDDO
      ENDIF
      GOTO 1
 
C***  LEVELS -----------------------------------------------------------
   10 N=N+1
      IF (LEVSEQ .NE. 0) THEN
         CALL REMARK ('DATOM: LEVEL CARD OUT OF SEQUENCE')
         GOTO 990
      ENDIF
      IF(N .GT. NDIM) THEN
         CALL REMARK ('DATOM : MORE LEVELS THEN DIMENSIONED (NDIM)')
         GOTO 990
      ENDIF
      IF (NATOM .NE. 0) NOM(N)=NATOM

      READ (KARTE,11,ERR=985) LEVEL(N),NCHARG(N),NW,ELEVEL(N),E,MAINQN(N)
   11 FORMAT(12X,A10,1X,I2,1X,I4,2F10.0,1X,I2)

      WEIGHT(N)=FLOAT(NW)

C***  If EION is empty, it will be repeated from last level
C***  after checking that the latter belongs to the same element and ion
      IF (E .EQ. 0.0) THEN
         IF (N .GT. 1) THEN    
            IF (NOM   (N-1) .NE. NOM   (N) .OR. 
     >          NCHARG(N-1) .NE. NCHARG(N) ) THEN
C***            Setting missing ionization energy to zero;
C***              this may be OK, if it is the highest level of an element
C***              and will be checked later (in the continuum block)
                  EION(N) = .0
            ELSE
                EION(N) = EION(N-1)
            ENDIF
         ELSE
            CALL REMARK ('ERROR: FIRST LEVEL WITHOUT IONIZATION ENERGY')
            GOTO 990
         ENDIF 
      ELSE
         EION(N) = E
C***  If EION is given redundantly, it must be identical  
         IF (N .GT. 1) THEN    
            IF (NOM   (N-1) .EQ. NOM   (N) .AND. 
     >          NCHARG(N-1) .EQ. NCHARG(N) .AND.
     >          EION  (N-1) .NE. E ) THEN
                  CALL REMARK ('ERROR: DIFFERENT IONIZATION ENERGIES')
                  GOTO 990
            ENDIF
         ENDIF
      ENDIF

C***  If level energy entry is empty, it is calculated from MainQN 
C***  by RYDBERG's formula
C***  If both are not given, this may be a mistake. Unfortunately,
C***  old DATOM files use these empty files for He III, H II as meaning
C***  ELEVEL = 0.0 -- therefore no error check can be made here
      IF (KARTE(31:40) .EQ. '          ') THEN 
         IF (KARTE(52:53) .NE. '  ') THEN
            F=FLOAT(MAINQN(N))
            ELEVEL(N) = (1.-1./F/F)*EION(N)
         ELSE
            ELEVEL(N) = .0
         ENDIF
      ENDIF
      GOTO 1
 
C***  LINE TRANSITIONS  ------------------------------------------------
   20 READ (KARTE,21) LEVUP,LEVLOW,AUPLOW,KRUDI,CEY,CO1,CO2,CO3,CO4
   21 FORMAT(10X,A10,2X,A10,G10.0,2X,A3,A4,1X,4G7.0)
      LEVSEQ=1
C***  FIND UPPER INDEX
      DO 22 J=1,N
      NUP=J
      IF (LEVEL(J).EQ.LEVUP ) GOTO 23
   22 CONTINUE
      CALL REMARK ('UPPER LINE LEVEL NOT FOUND')
      GOTO 990
C***  FIND LOWER INDEX
   23 DO 24 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 25
   24 CONTINUE
      CALL REMARK ('LOWER LINE LEVEL NOT FOUND')
      WRITE (0,*) 'LOWER LEVEL = ',LEVLOW
      GOTO 990
   25 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK ('ERROR: LINE BETWEEN DIFFERENT ELEMENTS')
            GOTO 990
         ENDIF
         IF (EION(LOW) .EQ. .0) THEN
            CALL REMARK ('ERROR: MISSING IONIZATION ENERGY')
            GOTO 990
         ENDIF

      ENDIF
      IF (NCHARG(NUP) .NE. NCHARG(LOW)) THEN
         CALL REMARK ('LINE BETWEEN DIFFERENT IONIZATION STAGES')
         GOTO 990
      ENDIF
      IF (NUP.LE.LOW) THEN
         CALL REMARK ('LINE TRANSITION INDICES WRONG')
         GOTO 990
      ENDIF
      IF (ELEVEL(NUP) < ELEVEL(LOW)) THEN
         CALL REMARK ('LINE ERROR: UPPER LEVEL NOT IN FRONT')
         GOTO 990
      ENDIF
C***  CORRECT LINE TRANSITION DETECTED:
      IND=IND+1

C***  ERROR STOP
      IF (IND .GT. MAXIND) THEN
         CALL REMARK ('ERROR: IND .GT. MAXIND')
         GOTO 990
      ENDIF
      INDNUP(IND)=NUP
      INDLOW(IND)=LOW
      EINST(NUP,LOW) = AUPLOW
      KEYCBB(IND)=CEY
      COCO(1,IND)=CO1
      COCO(2,IND)=CO2
      COCO(3,IND)=CO3
      COCO(4,IND)=CO4
C***  RUDIMENTAL TRANSITIONS ARE MARKED BY -2. IN THE TRANSPOSED
C***  MATRIX ELEMENT  EINST(LOW,NUP)
      IF (KRUDI.NE.'   ') EINST(LOW,NUP)=-2.
      LASTIND=IND
      GOTO 1
 
C***  CONTINUUM TRANSITIONS    -----------------------------------------
   30 READ (KARTE,31,ERR=990) 
     >    LEVLOW, SIGMA, ALPLOW, SLOW, IGLOW, ICBF, LEVION
   31 FORMAT (10X,A10,3G10.0,1X,A8,1X,A8,1X,A10)
      LEVSEQ=1
C***  FIND LOWER INDEX
      DO 34 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 35
   34 CONTINUE
      CALL REMARK ('LOWER CONTINUUM LEVEL NOT FOUND')
      GOTO 990
   35 CONTINUE
C***  FIND UPPER INDEX
      IF (LEVION .EQ. '          ') THEN
          EMIN=999999.
          NUP=0
          DO 39 J=1, N
          IF ((NATOM .GT. 1) .AND. (NOM(J) .NE. NOM(LOW))) GOTO 39
          IF ((NCHARG(J) .EQ. NCHARG(LOW)+1).AND.(ELEVEL(J) .LT. EMIN)) 
     $       THEN
             NUP=J
             EMIN=ELEVEL(J)
          ENDIF 
   39     CONTINUE
          IF (NUP .NE. 0) GOTO 33
      ELSE
          DO 32   J=1, N
          NUP=J
          IF ((LEVEL(J) .EQ. LEVION).AND.(NCHARG(J) .EQ. NCHARG(LOW)+1))
     $        GOTO 33
   32     CONTINUE
      ENDIF

      WRITE (0,*) 'ERROR: UPPER CONTINUUM LEVEL ', LEVION, 
     $         ' NOT FOUND'
      GOTO 990

   33 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK ('CONTINUUM BETWEEN DIFFERENT ELEMENTS')
            GOTO 990
         ENDIF
      ENDIF
C***  CORRECT CONTINUUM TRANSITION DETECTED:
      KONT=KONT+1
C***  ERROR STOP
      IF (KONT .GT. MAXKONT) THEN
         CALL REMARK ('ERROR: MORE CONTINUA THAN DIMENSIONED (MAXKONT)')
         GOTO 990
      ENDIF
      KONTNUP(KONT)=NUP
      KONTLOW(KONT)=LOW
      EINST(LOW,NUP)=SIGMA
      ALPHA(KONT)=ALPLOW
      SEXPO(KONT)=SLOW
      IGAUNT(KONT)=IGLOW
C***  KEYCBF(KONT)=ICBF
      IF (IGLOW .EQ. 'PIKB12' .OR. IGLOW .EQ. 'OPAPROIX') THEN
C*** READ FURTHER INFORMATION (ADDCON1-3) IN FOLLOWING LINE
  121   READ(4,122,END=3) KARTE
  122   FORMAT(A)
  130   READ (KARTE,131) ADD1, ADD2, ADD3
  131   FORMAT (21X, 3G10.0)
      ELSE
        ADD1 = 0.
        ADD2 = 0.
        ADD3 = 0.
      ENDIF
      ADDCON1(KONT) = ADD1
      ADDCON2(KONT) = ADD2
      ADDCON3(KONT) = ADD3
      LASTKON=KONT
      GOTO 1
 
C***  K-SHELL IONISATION
   40 CONTINUE
C***  DOES AKTUELL ELEMENT FIT TO K-SHELL-DATA ?
      IF (SYMBOL(NATOM) .NE. KARTE(16:17)) GOTO 982

C***  Are the K-shell data split for ionization stages?
      IF (KARTE(18:20) .NE. '   ') THEN
         READ (KARTE(18:20),'(I3)',ERR=980) ISTAGE 
         IF (ISTAGE .LT. 1) GOTO 981
         IF (ISTAGE .GT. NZ-2) GOTO 981
         READ (KARTE,41,ERR=983) SIGMATHK(NATOM,ISTAGE), 
     >                           SEXPOK  (NATOM,ISTAGE)
   41    FORMAT (20X,2G10.0)
C***     Special option: energy might be given in electron Volts
         IF (KARTE(52:53) .EQ. 'EV') THEN
            READ (KARTE(44:51), '(F8.0)', ERR=983) EVOLT
            EDGEK(NATOM,ISTAGE) = EVOLT * 1.E8 / 12397.7
         ELSE 
            READ (KARTE(44:53), '(F10.0)', ERR=983) EDGEK(NATOM,ISTAGE)
         ENDIF
      ELSE
         READ (KARTE,42,ERR=983) SIGMATHK(NATOM,1), 
     >                           SEXPOK  (NATOM,1),
     >                           EDGEK   (NATOM,1)
   42    FORMAT (20X,2G10.0,3X,F10.0)
         DO ISTAGE=2, NZ-2
            SIGMATHK(NATOM,ISTAGE) = SIGMATHK(NATOM,1)
            SEXPOK  (NATOM,ISTAGE) = SEXPOK  (NATOM,1)
            EDGEK   (NATOM,ISTAGE) = EDGEK   (NATOM,1)
         ENDDO
      ENDIF

      LEVSEQ=1
      GOTO 1

C***  SUM OF TRANSITIONS TO UPPER LEVELS WHICH ARE ASSUMED TO BE IN LTE
   50 CONTINUE
      WRITE (0,*) 'ERROR: the LTESUM branch is no longer supported' 
      WRITE (0,*) 'ERROR: the LTESUM branch is no longer supported' 
      STOP 'FATAL ERROR reported from subr. DATOM'

c      READ (KARTE,51) LEVLOW,IRANGE,ASUM,COEFF1,COEFF2
c   51 FORMAT (10X,A10,1X,A8,1X,G9.0,1X,F7.0,1X,F7.0)
c      LEVSEQ=1
cC***  FIND LOWER INDEX
c      DO 52 J=1,N
c      LOW=J
c      IF (LEVEL(J) .EQ. LEVLOW) GOTO 53
c   52 CONTINUE
c      CALL REMARK ('LOWER LTESUM LEVEL NOT FOUND')
c      GOTO 990
c   53 CONTINUE
c      ALTESUM(1,LOW)=ASUM
c      ALTESUM(2,LOW)=COEFF1
c      ALTESUM(3,LOW)=COEFF2
c      ENCODE (8,54,ALTESUM(4,LOW)) IRANGE
c      WRITE  (8,54,ALTESUM(4,LOW)) IRANGE
c   54 FORMAT (A8)
c      GOTO 1

C***  AUTOIONIZATION AND DIELECTRONIC RECOMBINATION  -------------------
C***  DRTRANSIT line encountered
   60 NAUTO=NAUTO+1
      IF (NAUTO .GT. MAXAUTO) THEN
         WRITE (0,*) '*** NAUTO .GT. MAXAUTO'
         GOTO 990
      ENDIF
C***  Iron line indices are arranged at the end of the range.
C***  Therefore, no DRTRANSIT lines may occur after the GENERIC element
      IF (LASTFE .GT. 0) THEN
         WRITE (0,*) '*** DRTRANSIT line found behind GENERIC element'
         GOTO 990
      ENDIF
      READ (KARTE,61,ERR=986) 
     >     LEVLOW,LEVUP,NW,EAUTO(NAUTO),AAUTO(NAUTO),LEVION,DRRUDI
   61 FORMAT(10X,A10,2X,A10,1X,I4,1X,F10.0,1X,G10.0,2X,A10,1X,A3)
      LEVUPAUTO(NAUTO) = LEVUP
C***  wrstart may perform some consistency checks
      IF (ROUTINE .EQ. 'WRSTART') THEN
         IF (NW .LE. 0) GOTO 987
         IF (KARTE(33:33) .NE. ' ') GOTO 988
         IF (KARTE(38:38) .NE. ' ') GOTO 988
      ENDIF
      WAUTO(NAUTO)=FLOAT(NW)
      LEVSEQ=1
C***  FIND LOWER INDEX
      DO 64 J=1,N
         LOW=J
         IF (LEVEL(J) .EQ. LEVLOW ) GOTO 65
   64 CONTINUE
      CALL REMARK ('LOWER LEVEL FOR DIEL. RECOMBINATION NOT FOUND')
      STOP 'DRLOWER'
   65 CONTINUE
C***  FIND INDEX OF PARENT ION
      IF (LEVION .EQ. '          ') GOTO 69
      DO 68 J=1, N
      NUP=J
      IF ((LEVEL(J) .EQ. LEVION) .AND. (NCHARG(LOW)+1 .EQ. NCHARG(J)))
     $    GOTO 63
   68 CONTINUE
      PRINT *, 'ERROR: PARENT ION FOR DR ', LEVION, 
     $         'NOT FOUND'
      CALL REMARK ('PARENT ION FOR DR NOT FOUND')
      GOTO 990
   63 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK 
     $        ('DIELECTRONIC RECOMBINATION BETWEEN DIFFERENT ELEMENTS')
            STOP 'DRATOMS'
         ENDIF
      ENDIF
   69 CONTINUE
      IF (LEVION .EQ. '          ') THEN
         IONAUTO(NAUTO)=0
      ELSE
         IONAUTO(NAUTO)=NUP
      ENDIF
      LOWAUTO(NAUTO)=LOW

C***  RUDIMENTAL DR-TRANSITIONS 
C***  in the current DRTANSIT data there are no entries in that column 
C***  in the future, this could be used to override the global 
C***  setting via the DRLEVELS option in CARDS
      IF (DRRUDI .NE. '   ') KRUDAUT(NAUTO) = 1

      GOTO 1

C***  END OF INPUT DATA REACHED  ---------------------------------------
    3 CLOSE (4)
 
C***  OLD DATOM FILE (CONTAINING ONLY A HELIUM MODEL ATOM) RECOGNIZED
      IF (NATOM .EQ. 0) THEN
         NATOM=1
         ELEMENT(1)='HELIUM    '
         SYMBOL(1)='HE'
         ATMASS(1)=4.
         STAGE(1)=3.
         KODAT(1)=1
         DO I=1,N
            NOM(I)=1
         ENDDO
      ENDIF

C***  SOME CHECKS OF THE INPUT DATA (only in WRSTART, else skipped)) 
      If (ROUTINE .NE. 'WRSTART') GOTO 78

      IF (N .EQ. 0) THEN
            CALL REMARK ('NO ENERGY LEVELS RECOGNIZED')
            STOP '*** ERROR detected in subr. DATOM'
      ENDIF

C***  ALL ELEMENTS ARE CHECKED ONE BY ONE FOR EXISTENCE OF ANY LEVEL
      DO 84 NA=1,NATOM
         DO 86 I=1,N
            IF (NOM(I) .EQ. NA) GOTO 84
   86    CONTINUE
         CALL REMARK ('ERROR: ELEMENT WITHOUT ANY LEVEL DECODED')
         STOP '*** ERROR detected in subr. DATOM'
   84 CONTINUE

C***  LEVELS ARE CHECKED FOR CORRECT ELEMENT MEMBERSHIP
      DO 85 J=1,N
         IF (LEVEL(J)(:2) .NE. SYMBOL(NOM(J))) THEN
            IF (LEVEL(J)(:1) .NE. SYMBOL(NOM(J))(:1) 
     >         .OR. SYMBOL(NOM(J))(2:2) .NE. ' ') THEN 
               CALL REMARK ('WRONG ELEMENT MEMBERSHIP OF LEVELS')
               GOTO 990
            ENDIF
         ENDIF
   85 CONTINUE

C***  TRANSITIONS ARE CHECKED FOR COMPLETENESS
      DO 7 I=1,N
         DO 7 J=1,N
            IF (NOM(I) .NE. NOM(J)) GOTO 7
            IF (ELEMENT(NOM(I)) .EQ. 'GENERIC') GOTO 7
            IF (NCHARG(I) .NE. NCHARG(J)) GOTO 8
            IF (I.LE.J) GOTO 7
            IF (EINST(I,J) .LE. -99.  ) THEN
            CALL REMARK ('LINE TRANSITION MISSING OR f-VALUE IS GREATER
     >-OR-EQUAL TO 99.0')
            WRITE (0,*) 'LEVEL-NO.: ',I,J
            WRITE (0,*) 'LEVEL-NAMES: "',LEVEL(I),'"   "',LEVEL(J),'"'
     >                  ,' f-VALUE: ', EINST(I,J)
            STOP '*** ERROR detected in subr. DATOM'
            ENDIF
            GOTO 7

C***     Continuum transition
    8    IF (I .GE. J) GOTO 7
C***     CHARGES MUST DIFFER BY 1
         IF (NCHARG(I)+1 .NE. NCHARG(J)) GOTO 7
C***     THERE MUST BE AT LEAST ONE CONTINUUM TRANSITION FOR EACH LOWER LEVEL
         DO 77 KON=1, LASTKON
            IF (KONTLOW(KON) .EQ. I) GOTO 7
   77    CONTINUE
         IF (EINST(I,J) .LT. .0 ) THEN
            CALL REMARK ('CONTINUUM TRANSITION MISSING')
            WRITE (0,'(4A)') 'LOWER LEVEL: ', LEVEL(I), 
     >                     '  UPPER LEVEL: ', LEVEL(J)
            STOP '*** ERROR detected in subr. DATOM'
            ENDIF
    7 CONTINUE
C***  Checks for completeness finished

C***  Consistency check for DRTRANSIT lines:
      IF (NAUTO .GT. 0) THEN
        DO I=1, NAUTO
         DO J=1, NAUTO
            IF (LEVUPAUTO(I) .EQ. LEVUPAUTO(J)) THEN
               IF (WAUTO(I) .NE. WAUTO(J)) THEN
                  WRITE (0,'(A,2I4)') '*** ERROR in DRTRANSIT data: ' //
     >                        'different statistical weights for ' //
     >                        'level ' // LEVUPAUTO(I),  
     >                        INT(WAUTO(I)), INT(WAUTO(J))  
                  GOTO 989
               ENDIF
               IF (EAUTO(I) .NE. EAUTO(J)) THEN
                  WRITE (0,'(A,2X,2F10.2)') 
     >                '*** ERROR in DRTRANSIT data: different energies ' //
     >                'for level ' // LEVUPAUTO(I),  
     >                EAUTO(I), EAUTO(J)  
                  GOTO 989
               ENDIF
            ENDIF
         ENDDO
        ENDDO
      ENDIF

   78 CONTINUE
C***  End of data checks only performed in WRSTART *********************

C***  GENERATE VECTORS NFIRST, NLAST: FIRST AND LAST LEVEL OF EACH ELEMENT
      DO 90 NA=1,NATOM
      IF (NA .EQ. 1) THEN
          NFIRST(NA)=1
      ELSE
          NFIRST(NA)=NLAST(NA-1)+1
      ENDIF
      IF (NA .LT. NATOM) THEN
          NLAST(NA)= ISRCHEQ(N,NOM(1),1,NA+1) - 1
      ELSE
          NLAST(NA)=N
      ENDIF
   90 CONTINUE

C***  GENERATION OF VECTOR IONGRND: DEFAULT LEVEL FOR IONIZATION (LOWEST
C***  LEVEL OF PARENT ION)
      DO 92 J=1, N
      IONGRND(J)=0
      EMIN=999999.
      DO 92 I=1, N
      IF ((NOM(I) .EQ. NOM(J)) .AND. (NCHARG(I) .EQ. NCHARG(J)+1) .AND.
     $    (ELEVEL(I) .LT. EMIN)) THEN
         EMIN=ELEVEL(I)
         IONGRND(J)=I
      ENDIF
   92 CONTINUE

C***  Convert f-values into EINSTEIN coeficients A_up,low
C***  NEGATIVE LINE-CARD ENTRIES INDICATE OSCILLATOR STRENGTHS
      DO 66 IND=1,LASTIND
         NUP=INDNUP(IND)
         LOW=INDLOW(IND)
         AUPLOW=EINST(NUP,LOW)
         IF (AUPLOW .GE. 0.0) GOTO 66
         WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
         EINST(NUP,LOW)=-0.6669*WAVENUM*WAVENUM*AUPLOW*WEIGHT(LOW)/
     /               WEIGHT(NUP)
   66 CONTINUE

C****************
C***  DRTRANSITs 
C****************

C***  CHECK for max. number of line transitions incl. DRTRANSITS
      IF (LASTIND+NAUTO+LASTFE > 99999) THEN
        WRITE (0,*)
     >     '*** MORE THAN 99999 LINE TRANSITIONS ENCOUNTERED ***'
        WRITE (0,*) 'This is not compatible with the encoding of the'
        WRITE (0,'(A)') ' line index in the MODEL file variables'
     >     // ' XJLnnnnn.'
        STOP '*** FATAL ERROR IN DATOM'
      ENDIF

C***  ASSIGNMENT OF DEFAULT IONIZATION LEVEL (GROUND STATE OF PARENT ION)
C***  FOR DIELECTRONIC RECOMBINATION TRANSITIONS
C***  NOTE: ASSUMPTION IS THAT ALL DOUBLY EXCITED STATES AUTOIONIZE
C***        INTO THE GROUND STATE OF THE PARENT ION
      DO 97 I=1, NAUTO
         IF (IONAUTO(I) .EQ. 0) IONAUTO(I)=IONGRND(LOWAUTO(I))
         IF (IONAUTO(I) .NE. IONGRND(LOWAUTO(I))) STOP 'IONAUTO'
C***     Check that the stabilizing transitions have positive wavelength
         LOW=LOWAUTO(I)
C***     WAVENUMBER OF STABILIZING TRANSITION
         WSTABIL = EION(LOW) - ELEVEL(LOW) + EAUTO(I)
         IF (WSTABIL .LE. .0) THEN
           WRITE (0,*)  '*** INCONSISTENCY IN DRTRANSIT DATA (DATOM):'
           WRITE (0,*)  '*** STABILIZING LINE HAS NEGATIVE WAVELENGTH'
           WRITE (0,*)  '*** Transition ', LEVEL(LOWAUTO(I)), ' - ', 
     >               LEVEL(IONAUTO(I))
           STOP '*** FATAL ERROR DETECTED BY SUBR. DATOM'
         ENDIF
   97 CONTINUE

C***  Conversion of oscillator strength f (indicated by neg-sign)
C***   into Aup-low
C***   Bug removed (wrong: EAUTO(LOW) ) wrh 10-Apr-2003 18:47:27
      DO 67 J=1,NAUTO
      AAUTOJ=AAUTO(J)
      IF (AAUTOJ .LT. 0.0) THEN
         LOW=LOWAUTO(J)
         WAVENUM=EION(LOW)-ELEVEL(LOW)+EAUTO(J)
         AAUTO(J)=-0.6669*WAVENUM*WAVENUM*AAUTOJ*WEIGHT(LOW)/WAUTO(J)
      ENDIF
   67 CONTINUE

C***  The DRTRANSITs are arranged behind the regular line transitions in 
C***      the vectors INDLOW, INDNUP
C***      Here, INDNUP is the index of the next-higher ground level
C***      The upper level names are stored in LEVUPAUTO(NAUTO) 
      DO IND=1, NAUTO
        INDLOW(LASTIND+IND) = LOWAUTO(IND)
        INDNUP(LASTIND+IND) = IONAUTO(IND)
      ENDDO

C***  Append the auto-ionizing levels to those vectors that specify
C***    the energy levels (index range N+1 ... N_WITH_DRLEVELS)
      N_WITH_DRLEVELS = N 
      IF (NAUTO .GT. 0)
     >    CALL APPEND_AUTOLEVELS (N, N_WITH_DRLEVELS, NDIM, MAXIND,
     >            MAXAUTO, NAUTO, LOWAUTO, IONAUTO, EAUTO, ELEVEL,
     $            LEVEL, EION, WEIGHT, INDLOW, INDNUP, LASTIND,
     >            LEVUPAUTO, LEVAUTO, WAUTO, NCHARG, IONGRND, NOM)

      RETURN

C***  ERROR exits **************************

  980 WRITE (0,*) 'Ionization Stage must be an integer number'
      GOTO 990

  981 WRITE (0,*) 'Ionization Stage outside valid range'
      GOTO 990

  982 WRITE (0,*) 'K-SHELL-DATA DO NOT FIT TO CURRENT ELEMENT'
      GOTO 990

  983 WRITE (0,*) 'K-SHELL DATA COULD NOT BE DECODED AS NUMBERS'
      GOTO 990

  985 WRITE (0,*) 'ERROR WHEN DECODING LEVEL CARD'
      GOTO 990

  986 WRITE (0,*) 'ERROR WHEN DECODING DRTRANSIT CARD'
      GOTO 990

  987 WRITE (0,*) 'stat. weight < 0 read from DRTRANSIT card'
      GOTO 990

  988 WRITE (0,*) 'Non-blank entry falls into column gap'
      GOTO 990

C***  ERROR BRANCH ********************************************
  989 WRITE (0,'(A,2I4)') '*** You must provide DRTRANSIT data ' //
     >                    'ion the corrected version (after June 2023)'
      GOTO 990


  990 WRITE (0,*) 'The Error was detected when decoding DATOM line:'
      WRITE (0,*) KARTE
      STOP 'ERROR detected by Subr. DATOM'

      END
      FUNCTION DBNUEDT (XLAMBDA,T)
C**********************************************************************
C***  DERIVATIVE OF PLANCK FUNCTION WITH RESPECT TO TEMPERATURE
C***  XLAMBDA IN ANGSTROEM, T IN KELVIN  
C***  BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ      
C**********************************************************************

C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
      DATA C1,C2 / 1.4388E8, 3.9724E+8 /

C***  PREVENT ERROR FOR NEGATIVE TEMPERATURES 
      IF (T .LE. .0) THEN
         DBNUEDT=.0
         ELSE
         HNUEKT=C1/XLAMBDA/T
C*** PREVENT ERROR AR004 (BAD SCALAR ARGUMENT TO ARLIB MATH ROUTINE)
C*** FOR EXTREME WIEN-DOMAIN
C***                                AND FOR RALEIGH-JEANS DOMAIN

C***   CRAY VERSION
C!!!         IF (HNUEKT .GT. 5000. .OR. HNUEKT .LT. 1.E-6) THEN
C***   DEC-VERSION
         IF (HNUEKT .GT. 700. .OR. HNUEKT .LT. 1.E-6) THEN
            DBNUEDT=0.
         ELSE
            X2=XLAMBDA*XLAMBDA
            X4=X2*X2
            T2=T*T
            EXFAC=EXP(HNUEKT)
            DBNUEDT=C1*C2/X4/T2/(EXFAC + 1./EXFAC - 2.)
            ENDIF
         ENDIF

      RETURN
      END
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
            SUBROUTINE DECFREQ (XLAMBDA, NF, NFDIM, TREF, OLDFGRID, KEY, 
     >                    MODOLD, XLAMBLUE)
C*******************************************************************************
C***  PURPOSE: READS THE FREQUENCY GRID XLAMBDA(K) (WAVELENGTHS IN A) 
C***       EITHER: INPUT FROM TAPE6 = FGRID
C***       ALTERNATIVELY: IF (OLDFGRID) : TAKEN FROM OLD MODEL FILE
C***       - NOTE: ONLY THE HAND-MADE POINTS ARE DEFINED HERE, 
C***               FURTHER POINTS (EDGES!) ARE INSERTED BY SUBR. FGRID
C***  CALLING TREE: WRSTART - FGRID - DECFREQ  
C*******************************************************************************
 
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
      CHARACTER KARTE*80, ACTPAR*20
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      DIMENSION XLAMBDA(NFDIM)

      REAL, INTENT(IN) :: XLAMBLUE
      LOGICAL OLDFGRID
 
C***  IF NO REFERENCE TEMPERATURE IS SPECIFIED, TEFF IS DEFAULT
      TREF=TEFF

C***  BRANCH FOR OLDFGRID OPTION:READ FROM CHANNEL 9 = OLD MODEL FILE 
      IF (OLDFGRID) THEN

      CALL OPENMS (9, IDUMMY, IDUMMY, 1, IERR)
      IERR=1
      CALL READMS (9, MODOLD, 13, 'MODHEAD ', IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK (' OLD MODEL FILE NOT AVAILABLE FOR READING FGRID')
         PRINT *,     ' OLD MODEL FILE NOT AVAILABLE FOR READING FGRID'
         STOP 'ERROR'
         ENDIF
      CALL READMS (9, NF     , 1 , 'NF      ', IERR)
C***  ARRAY BOUND CHECK
      IF (NF .GT. NFDIM) THEN
         CALL REMARK (' OLD MODEL HAS TOO MANY FREQUENCY POINTS')
         PRINT *,      'OLD MODEL HAS TOO MANY FREQUENCY POINTS'
         STOP 'ERROR'
         ENDIF
      CALL READMS (9, XLAMBDA, NF, 'XLAMBDA ', IERR)
      CALL READMS (9, KEY    , NF, 'KEY     ', IERR)
      CALL CLOSMS (9, IERR)

C***  REMOVE ALL ENTRIES WHICH ARE NOT HAND-MADE FREQUENCY POINTS 
C***    (RECOGNIZED BY BLANK ENTRY IN ARRAY "KEY")
      KK = 0
      DO 10 K=1, NF
        IF (KEY(K) .NE. '        ') GOTO 10
        KK = KK + 1
        XLAMBDA(KK) = XLAMBDA(K)
   10 CONTINUE
      NF = KK 

      GOTO 7
      ENDIF
C***  END OF THE BRANCH with OLD FGRID ****************************

C*** If BLUEMOST is given in CARDS, use this number instead of reading FGRID 
      IF (XLAMBLUE .GT. .0) THEN
         NF = 1
         XLAMBDA(1) = XLAMBLUE
         RETURN
      ENDIF


C***  ELSE: DECODING INPUT CARDS FROM TAPE 7 = FGRID
      NF=0
      OPEN (7, FILE='FGRID', STATUS='UNKNOWN')
    1 READ (7,6, END=7) KARTE
    6 FORMAT(A)

      IF (KARTE(:1) .EQ.'*' ) THEN
            PRINT 2,KARTE
    2       FORMAT (1X,A)
            GOTO 1
            ENDIF
      IF (KARTE(:4) .EQ. 'TREF' ) THEN
            CALL SARGV (KARTE, 2, ACTPAR)
            READ (ACTPAR,'(F20.0)', ERR=99) TREF
            GOTO 1
            ENDIF
      NF=NF+1
      IF(NF.GT.NFDIM) THEN
            CALL REMARK ('FREQUENCY DIMENSION NFDIM INSUFFICIENT')
            STOP 'ERROR in SUBR. DECFREQ'
            ENDIF
            READ (KARTE,'(F20.0)', ERR=99) XLAMBDA(NF)
    5 FORMAT(F10.0)
      GOTO 1
 
C*********** ERROR STOPS ***********************************************

    7 IF(NF .LE. 0) THEN
            CALL REMARK ('*** NO FREQUENCY SCALE ENCOUNTERED')
            CALL REMARK ('You must either use OLD FGRID')
            CALL REMARK ('or provide a file FGRID')
            CALL REMARK ('or specify BLUEMOST=x.x in the CARDS file')
            STOP '*** FATAL ERROR in SUBR. DECFREQ'
            ENDIF
 
      DO 21 K=2,NF
      IF((XLAMBDA(K-1)-XLAMBDA(K))*(XLAMBDA(1)-XLAMBDA(NF)).LE..0) THEN
            CALL REMARK ('WAVELENGTH SCALE OUT OF SEQUENCE:')
            WRITE (0,'(F20.2)') XLAMBDA(K)
            STOP 'ERROR in SUBR. DECFREQ'
            ENDIF
   21 CONTINUE

      RETURN

C*************
   99 WRITE (0,*) '*** File FGRID: CANNOT DECODE FLOATING POINT NUMBER!'
      WRITE (0,*) 'The error occured in the following line:'
      WRITE (0,*) KARTE
      STOP 'ERROR in SUBR. DECFREQ'

      END
      SUBROUTINE DECSTAR (MODHEAD, FM, RSTAR, VDOP, RMAX, TTABLE, 
     >                    NDDIM, OLDTEMP, MAXATOM, NATOM, ABXYZ, KODAT,
     >                    VPLOT, ATMASS, XDATA, MAXXDAT, OLDFGRID, 
     >                    THIN, ThinCard, GLOG, GEFFLOG, GEDD, 
     >                    bSaveGEFF, XMSTAR, WRTYPE, TAUMAX, TAUACC,
     >                    bTauFix, BTWOT, TFAC, DENSCON_LINE, BLACKEDGE, 
     >                    bOLDSTART, RadiusGridParameters, XMDOT, XLOGL,
     >                    RTRANS, BTAUR, DENSCON_FIX, MASSORIGIN, 
     >                    LRTinput, ELEMENT, iOldStratification, 
     >                    bHYDROSOLVE, GEddFix, bOldMdot, VoldMod, 
     >                    bOVTauMax, MLRELATION, NC, VTURBND, 
     >                    VTURB_LINE, iOLDRAD, RADGAMMASTART, 
     >                    fHYDROSTART, bFULLHYDROSTAT, bGAMMARADMEAN, 
     >                    cVEXTEND, bGREYSTART, GEFFKEY, POPMIN, 
     >                    XLAMBLUE, LTESTART, bOLDMODEL, XMDOTold, 
     >                    XMSTARold, RSTARold, VNDold, VFINALold, 
     >                    bForceDCUpdate, bOLDJ, bOVTauCut, bDCSCALE,
     >                    TEFFold, XLOGLold, MDOTINPUT, DRLINES_CARD)
C***********************************************************************
C***  DECODES INPUT CARDS, CALLED FROM WRSTART
C***********************************************************************
 
      IMPLICIT NONE

C***  COMMON/VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS TO FUNCTION WRVEL
C***  The second line has additional parameters for the 2-beta-law 
C***     -- wrh  6-Apr-2006 17:21:57
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, RONSET,
     >     bSMOCO, VCON, VOFF

C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE TO: 
C***                    WRSTART, GREY, PRIMOD, DECFREQ, DTDR, JSTART 
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
     
C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      REAL :: FM, RSTAR, VDOP, RMAX, RTRANS, TEFF, 
     >          TMIN, TMODIFY, TFAC, RONSET,
     >          VMIN, VFINAL, BETA, VPAR1, VPAR2, VCON,
     >          RCON, HSCALE, BETA2, BETA2FRACTION, VOFF,
     >          VPAR1_2, VPAR2_2, fHYDROSTART, RSTARold,
     >          XMDOT, XMSTAR, XMSTARG, XLOGL, XMSTARold,
     >          ABUND, ABREST, DENSUM, SUM, XMDOTold, VNDold,
     >          XHY, YHE, XC, XO, tempREAL, POPMIN, 
     >          VFINALold, XLAMBLUE, XMDOTTRANS,
     >          XLOGLold, TEFFOLD

      INTEGER :: NDDIM, MAXATOM, NATOM, MAXXDAT,
     >           MASSORIGIN, NPAR, IDXA, IDX,
     >           KHY, KHE, KC, KO, IERR, MFORM, NC, 
     >           LcalcCond, RcalcCond
      INTEGER, INTENT(OUT) :: GEddFix, LRTinput, iOldStratification, 
     >                        iOLDRAD, MDOTINPUT

      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS
      REAL, DIMENSION(MAXXDAT) :: XDATA
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      INTEGER, DIMENSION(8) :: DTVALUES
      LOGICAL :: LTESTART, TTABLE, SPHERIC, 
     >           OLDTEMP, VPLOT, OLDFGRID, bOLDMODEL, bOLDJ,
     >           ABMASS,ABNUMB,THIN, BTWOT, RMAX_IN_RSUN, BTAUR
      LOGICAL, INTENT(OUT) :: bSaveGEFF, bOLDSTART, bGREYSTART,
     >                        bTauFix, bForceDCUpdate, bOVTauCut,
     >                        bHYDROSOLVE, bOldMdot, bOVTauMax,
     >                        bFULLHYDROSTAT, bGAMMARADMEAN, bDCSCALE
      LOGICAL :: bCalcGLOGfromM, bLRTcomplete, bSMOCO
      CHARACTER(100) :: MODHEAD
      CHARACTER(8)  :: DAT, TIM, GEFFKEY, cVEXTEND
      CHARACTER(80) :: KARTE, ACTPAR, NEXTPAR, ThinCard
      CHARACTER(120) :: DRLINES_CARD
      CHARACTER(2)  :: WRTYPE
      CHARACTER(20) :: ACTPAR2
      CHARACTER(10) :: SYS, NODE, REL, VER, MACH
      CHARACTER(33) :: HOST
      CHARACTER*(*) :: DENSCON_LINE, VTURB_LINE
      CHARACTER(10), DIMENSION(NATOM) :: ELEMENT
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters     !contains all RADIUS-GRID CARDS (for subr. RGRID)
      CHARACTER(9) :: MLRELATION

      REAL :: BLACKEDGE, DENSCON_FIX
      REAL :: GLOG, GEFFLOG, GEDD, q, QLOG, RADGAMMASTART
      REAL :: TAUMAX, TAUACC, VoldMod, VTURBND

      !Laufvariablen und co.
      INTEGER :: I, K, NA, NZ, IPAR,
     >           IFOUNDELEMENT

      !Konstanten
      REAL, PARAMETER :: RSUN = 6.96E10     !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: TEFFSUN = 5780.    !EFFECTIVE TEMPERATURE OF THE SUN
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


C***  DEFAULT VALUES **********
      XDATA(1)=.0
      LTESTART=.FALSE.
      DO i=1, 3, 1
        RadiusGridParameters(i) = ' '
      ENDDO
      TEFF=-1.
      RSTAR=-1.
      RMAX=-1.
      XLOGL = -99.
      XMSTAR = -99.
      GLOG = -99.
      RTRANS = -99.
      QLOG = -99.
      VDOP=-1.
      XMDOT= -99.
      XMDOTTRANS= -99.
      RMAX_IN_RSUN = .FALSE.
      DO NA=1,MAXATOM
        ABXYZ(NA)=0.
      ENDDO
      TTABLE=.FALSE.
      SPHERIC=.TRUE.
      OLDTEMP=.FALSE.
      OLDFGRID = .FALSE.
      TMODIFY=0.0
      TMIN=.0
      VPLOT=.FALSE.
      ABMASS=.FALSE.
      ABNUMB=.FALSE.
      THIN=.FALSE.
      TAUMAX = 0.
      TAUACC = 1.E-4
      bTauFix = .FALSE.
      BTWOT = .FALSE.
      DENSCON_FIX = 1.
      DENSCON_LINE = ' '
      VTURB_LINE = ' '
      DRLINES_CARD = ''
      BLACKEDGE = .0
      BTAUR = .FALSE.
      BETA2FRACTION = .0
      VCON = 1.0      
      bSMOCO = .FALSE.
      bForceDCUpdate = .FALSE.
      bDCSCALE = .FALSE.
      
      q = -1.0
      GEDD = -1.0
      XLAMBLUE = -1. 
C***  Number of core-rays
      NC = 4
C***  determines if GEDD is directly specified (=2) or implied (=1):
      GEddFix = 0          
C***  nedative default value indicates that GEFF has not been set in CARDS
      GEFFLOG = -1.0
      bSaveGEFF = .FALSE.    !Determines if GEFFLOG is fixed in the MODEL file
      bCalcGLOGfromM = .FALSE.
      bOLDSTART = .FALSE.   !true if OLDSTART CARDS option has been set
      bGREYSTART = .FALSE.
      bOLDJ = .FALSE.       !true if radiation field XJC from old model is used 
      iOldStratification = 0  !> 0 if stratification (ND, VELO, RADIUS, DENSCON) from old model should be kept (@todo: also MDOT?)
      bHYDROSOLVE = .FALSE. !true if hydro iteration will be done in STEAL
      bOldMdot = .FALSE.
      bOVTauMax = .FALSE.   !true if option TAUMAX has been set on OLD STRATIFICATION card
      bOVTauCut = .FALSE.
      cVEXTEND = ' '
      fHYDROSTART = -99.    !starting from OLD MODEL using v infered from hydro (fraction > 0.)
      VoldMod = 1.
      bFULLHYDROSTAT = .FALSE.
      bGAMMARADMEAN = .FALSE.
      iOLDRAD = 0
      RADGAMMASTART = -99.
      POPMIN = 1.E-25           !default value should be the same as in STEAL -> DECSTE
      GEFFKEY = '        '
C***  Mass-Luminosity relation default: 
      MFORM = 2     ! 1 = Langer (1989); 2 = Graefener et al. (2011) 
      WRTYPE = ''   ! default: automatic choice of type for M-L relation
C***  LcalcCont = Codenumber for the conditions to calculate the Luminosity
C***    0 =
C***    1 = Mass input given
C***    2 = Eddington-Gamma input given
C***    3 = 1 + 2 = Mass and Eddington-Gamma input given 
C***                -> L can be determined if R is given
      LcalcCond = 0
     
C***  RcalcCond = Codenumber for the calculation of the radius
C***    0 = not yet defined
C***    1 = direct input
C***    2 = LRT relation
C***    4 = inferred from M and log g
      RcalcCond = 0

C***  MASSORIGIN = Codenumber for source of the stellar mass (for PRIPARAM)
C***    0 = Mass-luminosity relation   
C***    1 = input   
C***    2 = log g   
C***    3 = log geff plus Eddington Gamma
      MASSORIGIN = 0
      VTURBND = .0

C***  MDOTINPUT: origin of the mass-loss rate
C      MDOTINPUT = -1 : not specified --> ERROR
C      MDOTINPUT = 1 :  MDOT specified 
C      MDOTINPUT = 2 :  RTRANS specified 
C      MDOTINPUT = 3 :  MDOTTRANS specified (def. Graefener & Vink 2013) 
C      MDOTINPUT = 4 :  LOG Q specified
      MDOTINPUT = -1 


C***  CONSTRUCT MODEL HEADER  *********
      IF (OPSYS .EQ. 'DEC/UNIX') THEN
        CALL MY_DATE (DAT)
C        CALL DATE_AND_TIME(VALUES=DTVALUES)
C        WRITE(UNIT=DAT,FMT='(I4,"/",I2,"/",I2)') 
C     >          DTVALUES(1), DTVALUES(2), DTVALUES(3)
        CALL MY_CLOCK(TIM)
      ELSE IF (OPSYS .EQ. 'CRAY') THEN
        !CALL DATE (DAT)
        CALL DATE_AND_TIME(VALUES=DTVALUES)
        WRITE(UNIT=DAT,FMT='(I4,"/",I2,"/",I2)') 
     >          DTVALUES(1), DTVALUES(2), DTVALUES(3)
        !CALL DATE_AND_TIME(DATE=DAT)
        CALL CLOCK(TIM)
      ELSE
        WRITE (0,*) 'OPSYS NOT RECOGNIZED'
        WRITE (0,*) 'OPSYS=', OPSYS
        STOP 'ERROR IN SUBR. DECSTAR'
      ENDIF
      MODHEAD='MODEL START'
      MODHEAD(13:)=DAT
      MODHEAD(25:)=TIM

C***  READ, ECHO AND DECODE EACH INPUT CARDS ********
C***  AND ADD INFORMATION ABOUT HOST (SUBR. UNAME: SYSTEM CALL) ********
C!!!      CALL UNAME(SYS, NODE, REL, VER, MACH)
          sys = ' '
      HOST = OPSYS
      IF (SYS .EQ. 'craSH') THEN
         HOST(10:) = ' M92 2/64 (craSH Kiel)'
      ELSE IF (SYS .EQ. 'sn5208') THEN
         HOST(10:) = ' EL 4/1024 (craSHi Kiel)'
      ELSE IF (SYS .EQ. 'sn422') THEN 
         HOST(10:) = '/216 (crax Berlin)'
      ELSE IF (SYS .EQ. 'sn1619') THEN
         HOST(10:) = '2E/264 (cray Berlin)'
      ELSE IF (SYS .EQ. 'DEC/UNIX') THEN
         HOST      = 'DEC/UNIX at Potsdam'
      ELSE
         HOST(11:) = '(unknown)'
      ENDIF

C***  CHECK IF SYSTEM IS INSTALLED CORRECTLY
      IF (SYS .EQ. 'DEC/UNIX' .AND. SYS .NE. OPSYS) THEN
        WRITE (0,*) 'SYSTEM IS NOT INSTALLED CORRECTLY'
        WRITE (0,'(3A)') 'SYS, OPSYS=', SYS, OPSYS
      ENDIF

      PRINT 1, HOST
    1 FORMAT ('1>>> ECHO OF INPUT CARDS <<<', 53X, 'HOST: ', A33, /, 
     >        1X, 27('-'), 52X, 7('='), //)
      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND (1)

      !-- Zeilenweises Abarbeiten der CARDS-Datei ----
   10 READ (1,11, END=100) KARTE
   11 FORMAT (A)

      PRINT 2,KARTE
    2 FORMAT (1X,A)
 

      IF (KARTE .EQ. '') GOTO 10
      CALL SARGV (KARTE, 1, ACTPAR)
      CALL SARGC (KARTE, NPAR)

      IF (KARTE(:5) .EQ. 'TABLE') THEN
C                         =====
            TTABLE=.TRUE.
      
      ELSE IF (ACTPAR(:6) .EQ. 'DRLINE') THEN
C                               ======  
            DRLINES_CARD = KARTE

      ELSE IF ((KARTE(:4) == 'THIN') .OR.
     >         (KARTE(:9) == 'HYDROSTAT')) THEN
C                              ====
            THIN=.TRUE.
            ThinCard = KARTE
            IF (NPAR > 2) THEN
              DO IPAR=3, NPAR
                CALL SARGV(ThinCard,IPAR,ACTPAR) 
                IF (ACTPAR(1:4) == 'FULL') THEN
                  bFULLHYDROSTAT = .TRUE.
                ENDIF
                IF (ACTPAR == 'MEAN') THEN
                  bGAMMARADMEAN = .TRUE.
                ENDIF
              ENDDO
            ENDIF
            
      ELSE IF (KARTE(:5) .EQ. 'PLANE') THEN
C                              =====
            SPHERIC=.FALSE.

C***  Option GREY_TAUSCALE_START for compatibility with versions before
C***  23-Jan-2016
      ELSE IF (KARTE(:9) == 'GREYSTART' .OR. 
     >         KARTE(:9) == 'GREY_TAUS'     ) THEN
C                            =========
            bGREYSTART = .TRUE.
            
      ELSE IF (KARTE(:8) == 'OLDSTART') THEN
C                            ========
            bOLDSTART = .TRUE.
            IF (NPAR > 1) THEN  
              DO IPAR=2, NPAR
                CALL SARGV (KARTE, IPAR, ACTPAR)
C                IF (ACTPAR == 'DEPART') THEN
C                  write info here!
C                ENDIF
C                IF (ACTPAR == 'TAU') THEN
C                  write info here!
C                ENDIF
ccc             possible further values: DEPART, TAU (only read by ADAPTER currently)               
              ENDDO
           ENDIF
           IF (.NOT. bOLDMODEL) GOTO 95
              
      ELSE IF (KARTE(:10) == 'HYDROSTART') THEN
C                             ==========
         fHYDROSTART = 1.
        IF (NPAR > 1) THEN
          CALL SARGV (KARTE,2,ACTPAR)
          READ (ACTPAR, '(F10.0)', ERR=99) fHYDROSTART
          IF (fHYDROSTART > 1. .OR. fHYDROSTART < 0.) THEN
            WRITE (hCPR,'(A)') 'INVALID DAMPING VALUE FOR HYDROSTART'
            GOTO 92
          ENDIF
        ENDIF

      ELSE IF (KARTE(:5) .EQ. 'OLD T') THEN
C                              =====
           OLDTEMP=.TRUE.
           IF (KARTE(:9) .EQ. 'OLD T TAU') BTAUR = .TRUE.
           IF (.NOT. bOLDMODEL) GOTO 95
            

      ELSE IF (KARTE(:9) .EQ. 'OLD FGRID') THEN
C                              =========
           OLDFGRID = .TRUE.
           IF (.NOT. bOLDMODEL) GOTO 95
           
      ELSE IF ((KARTE(:9) .EQ. 'OLD STRAT')
C                               =========
     >     .OR. (KARTE(:5) .EQ. 'OLD V')) THEN
C                                ======   
        IF (.NOT. bOLDMODEL) GOTO 95
        iOldStratification = 1
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (NPAR > 2) THEN
          DO i=3, NPAR 
            CALL SARGV (KARTE, i, ACTPAR)
            SELECTCASE (ACTPAR)
              CASE ('VFINAL', 'VINF')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    VoldMod = tempREAL
                    !store abolute value as negative number, modificator as positive
                    VoldMod = -1. * VoldMod
                  ENDIF
                ENDIF
              CASE ('MOD', 'VMOD', 'MODFAK')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    !store abolute value as negative number, modificator as positive
                    VoldMod = tempREAL
                  ENDIF
                ENDIF
              CASE ('GRID')
                iOldStratification = 2
              CASE ('TAU') 
                cVEXTEND = 'TAU'
              CASE ('EXTRAP') 
                cVEXTEND = 'EXTRAP'
              CASE ('SCALE', 'STRETCH')
                cVEXTEND = 'STRETCH'
              CASE ('DCNEW', 'DENSCON-UPDATE')
                bForceDCUpdate = .TRUE.
              CASE ('DCSCALE', 'DENSCON-SCALE')
                bDCSCALE = .TRUE.
              CASE ('MDOT', 'KEEPMDOT')
                bOldMdot = .TRUE.
              CASE ('TAUCUT')
                bOVTauCut = .TRUE.
              CASE ('TAUMAX')
C***            revised Mar 2019: still iterate for TAUMAX by shifting V up/down to "hit" TAUMAX
                bOVTauMax = .TRUE.
            ENDSELECT
          ENDDO
        ENDIF
            
      ELSE IF (KARTE(:7) .EQ. 'TMODIFY') THEN
C                              =======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F20.0)',ERR=92) TMODIFY

      ELSE IF (KARTE(:4) .EQ. 'TMIN') THEN
C                              ====
         IF (KARTE(:10) /= 'TMIN-START') THEN
           WRITE (hCPR, '(A)') '*** WARNING: Deprecated TMIN card'
     >       // ' used: Please use TMIN-START card instead!'
         ENDIF
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F20.0)',ERR=92) TMIN

      ELSE IF (KARTE(:8) .EQ. 'LTESTART') THEN
C                              ========
         LTESTART=.TRUE.
         IF (NPAR >= 3) THEN 
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'BLACKEDGE') THEN
               CALL SARGV (KARTE, 3, ACTPAR)
               READ (ACTPAR, '(F10.0)', ERR=92) BLACKEDGE
            ENDIF
         ENDIF

      ELSE IF (KARTE(:6) .EQ. 'PLOT V') THEN
C                              ======
            VPLOT=.TRUE.

      ELSE IF (KARTE(:5) .EQ. 'RGRID') THEN
C                              =====
            RadiusGridParameters(1) = KARTE

      ELSE IF (ACTPAR .EQ. 'RTRANS') THEN
C                           ======
         CALL SARGV (KARTE, 2, ACTPAR)
         IDXA = IDX(ACTPAR)
         IF (ACTPAR(IDXA-2:IDXA) .EQ. 'DEX') THEN
             READ (ACTPAR(:IDXA-3), '(F20.0)',ERR=92) RTRANS
             RTRANS = 10.**RTRANS
         ELSE
            READ (ACTPAR, '(F20.0)',ERR=92) RTRANS
         ENDIF
         MDOTINPUT = 2

      ELSE IF (KARTE(:5) .EQ. 'LOG Q') THEN
C                              =====
        CALL SARGV (KARTE, 3, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) QLOG
        MDOTINPUT = 4
      
      ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
C                           =====
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'OLD') THEN
          IF (.NOT. bOLDMODEL) GOTO 95
          RSTAR = RSTARold/RSUN
          WRITE (hCPR,'(A,F8.4)')
     >      'RSTAR TAKEN FROM OLD MODEL: ', RSTAR
        ELSE 
          READ (ACTPAR, '(F20.0)',ERR=92) RSTAR
        ENDIF
        RcalcCond = 1
 
      ELSE IF (KARTE(:5) .EQ. 'LOG L') THEN
C                              =====
        CALL SARGV (KARTE, 3, ACTPAR)
        IF (ACTPAR == 'OLD') THEN
          IF (.NOT. bOLDMODEL) GOTO 95
          XLOGL = XLOGLold
          WRITE (hCPR,'(A,F8.4)')
     >      'LUMINOSITY TAKEN FROM OLD MODEL: log L/Lsun = ', XLOGL
        ELSE 
          READ (ACTPAR, '(F20.0)',ERR=92) XLOGL
        ENDIF
  
      ELSE IF (KARTE(:6) .EQ. 'VELPAR') THEN
C                              ======
         VMIN = VNDold
         VFINAL = VFINALold
         CALL DECVELPAR(KARTE, VFINAL, VMIN, BETA, RMAX, VOFF)         
         
      ELSE IF (ACTPAR .EQ. '2BETALAW') THEN
         CALL DEC2BETA(KARTE, BETA2, BETA2FRACTION, RONSET)
c         CALL SARGV (KARTE,3,ACTPAR)
c         READ (ACTPAR, '(F10.0)', ERR=99) BETA2 
c         CALL SARGV (KARTE,5,ACTPAR)
c         READ (ACTPAR, '(F10.0)', ERR=99) BETA2FRACTION 
 
      ELSE IF (KARTE(:4) .EQ. 'VDOP') THEN
C                              ====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) VDOP
 
      ELSE IF (KARTE(:8) .EQ. 'HEADLINE') THEN
C                              ========
        MODHEAD(35:) = KARTE(10:)
 
      ELSE IF (ACTPAR .EQ. 'MDOT') THEN
C                           ====
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'OLD') THEN
          IF (.NOT. bOLDMODEL) GOTO 95
          IF (XMDOTold < -100.) GOTO 96
          XMDOT = XMDOTold
          WRITE (hCPR,'(A,F8.4)')
     >      'MASS LOSS RATE TAKEN FROM OLD MODEL: ', XMDOT
        ELSE 
          READ (ACTPAR, '(F20.0)',ERR=92) XMDOT
        ENDIF
        MDOTINPUT = 1

      ELSE IF (ACTPAR(1:5) == 'MDOTT' .OR. ACTPAR == 'MDTRANS') THEN
C                              =====                  =======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) XMDOTTRANS
        MDOTINPUT = 3

      ELSE IF (ACTPAR .EQ. 'TEFF') THEN
C                           ====
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'OLD') THEN
          IF (.NOT. bOLDMODEL) GOTO 95
          TEFF = TEFFold
          WRITE (hCPR,'(A,F8.4)')
     >      'TSTAR (in CARDS: TEFF) TAKEN FROM OLD MODEL: ', TEFF
        ELSE
          IDXA = IDX(ACTPAR)
          IF (ACTPAR(IDXA-2:IDXA) == 'DEX') THEN
            READ (ACTPAR(:IDXA-3), '(F20.0)',ERR=92) TEFF
            TEFF = 10.**TEFF
          ELSE
            READ (ACTPAR, '(F20.0)',ERR=92) TEFF
          ENDIF
        ENDIF

      ELSE IF (ACTPAR .EQ. 'NCORE') THEN
C                           =====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(I3)', ERR=92) NC
        
      ELSE IF (KARTE(:5) .EQ. 'MSTAR') THEN
C                              =====
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR .NE. '?') THEN
          IF (MASSORIGIN .NE. 0) THEN
            IF (RcalcCond > 0) THEN
              WRITE (0, '(A)') '*** DOUBLE DEFINITION OF STELLAR MASS'
              GOTO 92
            ELSE
              RcalcCond = 4
            ENDIF
          ELSE
            MASSORIGIN = 1
            LcalcCond = LcalcCond + 1
          ENDIF

          IF (ACTPAR == 'OLD') THEN
            IF (.NOT. bOLDMODEL) GOTO 95
            XMSTAR = XMSTARold
            WRITE (hCPR,'(A,F8.4)')
     >        'STELLAR MASS TAKEN FROM OLD MODEL: ', XMSTAR
          ELSE 
            READ (ACTPAR, '(F20.0)', ERR=92) XMSTAR
          ENDIF

        ENDIF
 
      ELSE IF ((KARTE(:10) .EQ. 'LOG G_GRAV') .OR.
     >         (KARTE(:9) .EQ. 'LOG GGRAV') .OR.
     >         (KARTE(:8) .EQ. 'LOG GRAV')) THEN
C                              =====
        CALL SARGV (KARTE, 3, ACTPAR)
        IF (ACTPAR .NE. '?') THEN
          IF (MASSORIGIN .NE. 0) THEN
            IF (RcalcCond > 0) THEN
              WRITE (0, '(A)') '*** DOUBLE DEFINITION OF STELLAR MASS'
              GOTO 92
            ELSE
              RcalcCond = 4
            ENDIF
          ELSE
            MASSORIGIN = 2
          ENDIF
          READ (ACTPAR, '(F20.0)', ERR=92) GLOG
        ENDIF
 
      ELSE IF ((KARTE(:9) == 'LOG G_EFF') .OR.
C                               =========
     >         (KARTE(:8) == 'LOG GEFF') .OR.
C                               ========
     >         ((KARTE(:5) == 'LOG G') .AND. 
C                              =====
        !note: for backward compartibility log g is always interpreted as log geff
        !      To avoid confusion LOG G is considered as deprecated, instead you
        !      should always write LOG GEFF or LOG GGRAV
     >          ((KARTE(6:6) == ' ') .OR. (KARTE(6:6) == '=') .OR.
     >           (KARTE(6:6) == ',') .OR. (KARTE(6:6) == ':') ) ) 
     >        ) THEN

        CALL SARGV (KARTE, 3, ACTPAR)
        IF (ACTPAR /= '?') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) GEFFLOG
          bSaveGEFF = .TRUE.
        ENDIF
        !Warning if deprecated syntax is used
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'G') THEN
          WRITE (hCPR, '(A)') '*** WARNING: Deprecated syntax LOG G '
     >      // ' used: Please use LOG GEFF (or LOG GGRAV) instead!'
        ENDIF
        IF (bSaveGEFF) THEN
          !Check if interpretation keyword exists
          GEFFKEY = 'AUTO    '
          CALL SARGV (KARTE, 4, ACTPAR)
          IF (NPAR > 4 .AND. ACTPAR == 'RADFORCE') THEN
            CALL SARGV (KARTE, 5, ACTPAR)
            SELECTCASE(ACTPAR)
              CASE ('RAD', 'FULL')
                GEFFKEY = 'RAD     '
              CASE ('THOMSON', 'ELECTRON', 'THOM', 'E', 'e')
                GEFFKEY = 'THOM    '
            ENDSELECT
          ELSEIF (NPAR >= 4) THEN
            WRITE (hCPR, '(A)') '*** ERROR: Invalid or missing '
     >        // ' parameter on LOG GEFF card!'
            WRITE (hCPR, '(A)') '    Use RADFORCE=FULL or '
     >        // 'RADFORCE=ELECTRON to specify the intention'
     >        // ' of the LOG GEFF card!'
            GOTO 92          
          ENDIF
        ENDIF

      ELSE IF (ACTPAR == 'EDDINGTON-GAMMA') THEN
C                         ===============
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR /= 'AUTO') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) GEDD
          bSaveGEFF = .TRUE.
          GEddFix = 2
          LcalcCond = LcalcCond + 2
        ENDIF

      ELSE IF (ACTPAR == 'QION-START') THEN

        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR /= 'AUTO') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) q
        ENDIF
 
      ELSE IF (ACTPAR == 'RADGAMMA-START') THEN
C                         ==============
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'CONT') THEN
C***      Approximate Gamma_rad via continuum opacities        
          iOLDRAD = 1
        ELSEIF (ACTPAR == 'OLD') THEN
C***      Take Gamma_rad from old MODEL
          IF (.NOT. bOLDMODEL) GOTO 95
          iOLDRAD = 2
          IF (NPAR > 2) THEN
            CALL SARGV (KARTE, 3, ACTPAR)
            !Use MEAN value only? (for WRSTART)
            IF (ACTPAR == 'MEAN') bGAMMARADMEAN = .TRUE.
          ENDIF
        ELSE 
C***      Take explicit given value        
          READ (ACTPAR, '(F20.0)', ERR=92) RADGAMMASTART        
        ENDIF  

      ELSE IF (KARTE(:7) == 'MLANGER') THEN
C                            =======
         MFORM = 1
      ELSE IF (KARTE(:6) == 'MGOETZ') THEN
C                            ======
         MFORM = 2
      ELSE IF (ACTPAR .EQ. 'TAUMAX') THEN
C                           ======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) TAUMAX
        IF (NPAR > 2) THEN
          DO i=3, NPAR 
            CALL SARGV (KARTE, i, ACTPAR)
            SELECTCASE (ACTPAR)
              CASE ('FIX')
                bTauFix = .TRUE.
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF
                ENDIF
C              CASE ('MIN')             !not used in wrstart
C                bTauStrict = .FALSE.
              CASE ('EPS', 'ACC')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF                  
                ENDIF           
              CASE ('REPS', 'RELEPS', 'RELACC')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL * TAUMAX
                  ENDIF                  
                ENDIF           
C              CASE ('REDUCE')          !not used in wrstart
C                IF (NPAR >= (i+1)) THEN
C                  CALL SARGV (KARTE, i+1, NEXTPAR)
C                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
C                  IF (IERR == 0) THEN
C                    ReduceTauCorrections = tempREAL
C                  ELSE
C                    ReduceTauCorrections = 0.5
C                  ENDIF                  
C                ENDIF
            ENDSELECT
          ENDDO
        ENDIF


      ELSE IF (ACTPAR == 'TAUFIX') THEN
C                         ======
        bTauFix = .TRUE.
        IF (NPAR > 1) THEN
          CALL SARGV (KARTE, 2, ACTPAR)
          IF (ACTPAR /= 'STRICT') THEN
            READ (ACTPAR, '(F10.0)', ERR=92) TAUACC
          ELSEIF (NPAR > 2) THEN
            CALL SARGV (KARTE, 3, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=92) TAUACC                                     
          ENDIF
        ENDIF

      ELSEIF (ACTPAR == 'HYDRO') THEN
C                        =====
        bHYDROSOLVE = .TRUE.
        
      ELSE IF (ACTPAR .EQ. 'TWOTEMP') THEN
C                           ======
        BTWOT = .TRUE.
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) TFAC

      ELSE IF (ACTPAR .EQ. 'DENSCON') THEN
C                           ======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) DENSCON_FIX
        DENSCON_LINE = KARTE
        
      ELSE IF (ACTPAR .EQ. 'WRTYPE' .OR. ACTPAR .EQ. 'STARTYPE') THEN
        CALL SARGV (KARTE, 2, WRTYPE)
        IF (WRTYPE .NE. 'OB' .AND. WRTYPE .NE. 'WN' .AND.
     >      WRTYPE .NE. 'WC') THEN
            WRITE (0, *) '*** ERROR: Invalid choice of WRTYPE' 
            WRITE (0, *) '*** Valid types are: OB, WN, WC'
            GOTO 92
        ENDIF
        
      ELSE IF (ACTPAR == 'VTURB') THEN
        VTURB_LINE = KARTE
        IF (NPAR .GT. 1) THEN
           CALL SARGV (KARTE, 2, ACTPAR2)
           READ (ACTPAR2,'(F20.0)',ERR=98) VTURBND
        ELSE
           GOTO 97
        ENDIF
        
      ELSE IF (ACTPAR == 'VMIC') THEN
        VTURB_LINE = KARTE
        IF (NPAR .GT. 1) THEN
           CALL SARGV (KARTE, 2, ACTPAR2)
           READ (ACTPAR2,'(F20.0)',ERR=98) VTURBND
           VTURBND = VTURBND / SQRT(2.)
        ELSE
           GOTO 97
        ENDIF
        
      ELSE IF (ACTPAR == 'POPMIN') THEN
C                         ======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) POPMIN
        
      ELSE IF (ACTPAR(:8) == 'BLUEMOST') THEN
C                             ========
        CALL SARGV (KARTE, 2, ACTPAR2)
        READ (ACTPAR2, '(F20.0)',ERR=92) XLAMBLUE

C*********** SPECIAL BLOCK FOR X-RAYS  **************************
      ELSE IF (ACTPAR .EQ. 'XRAY') THEN
C                           ====
        DO I=2, NPAR-1
          CALL SARGV (KARTE, I, ACTPAR)
          IF (ACTPAR .EQ. 'XFILL') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(1)
          ELSEIF (ACTPAR .EQ. 'XRAYT') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(2)
          ELSEIF (ACTPAR .EQ. 'XRMIN') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(3)
C         Differential Emission Measure for X-rays
          ELSEIF (ACTPAR .EQ. 'DIFF-EM-EXP') THEN
           CALL SARGV (KARTE, I+1, ACTPAR)
           READ (ACTPAR, '(F20.0)', ERR=92) XDATA(4)
           IF (XDATA(4) .NE. 1.5 .AND. XDATA(4) .NE. 2.5) THEN
              WRITE (0,*) '*** SORRY, only values 1.5 or 2.5 permitted'
              WRITE (0,*) '*** as exponents for diff. emmission measure'
              GOTO 92
           ENDIF
C         Or second component
          ELSEIF (ACTPAR .EQ. 'XFILL2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(5)
          ELSEIF (ACTPAR .EQ. 'XRAYT2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(6)
          ELSEIF (ACTPAR .EQ. 'XRMIN2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(7)
          ENDIF
        ENDDO
C       Check for consistency
        IF (XDATA(4) .GT. 0. .AND. XDATA(5) .GT. 0.) THEN
           WRITE (0,*) '*** SORRY, DIFF-EM-EXP does NOT allow 2nd XFILL'
           GOTO 92         
        ENDIF

C***  BLACK EDGE IN THE RADIATION FIELD FOR NEW-START
      ELSE IF (ACTPAR .EQ. 'JSTART') THEN
         LTESTART = .FALSE.
         IF (NPAR .GE. 3) THEN 
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'BLACKEDGE') THEN
               CALL SARGV (KARTE, 3, ACTPAR)
               READ (ACTPAR, '(F10.0)', ERR=92) BLACKEDGE
            ELSEIF (ACTPAR == 'OLD') THEN
               bOLDJ = .TRUE.
               IF (.NOT. bOLDMODEL) GOTO 95
            ENDIF
         ENDIF

      ELSE IF (KARTE(:5) == 'OLD J') THEN
         LTESTART = .FALSE.
         bOLDJ = .TRUE.
         IF (.NOT. bOLDMODEL) GOTO 95

      ELSE IF (ACTPAR .EQ. 'SPECIAL_OUTER_POINTS') THEN
         RadiusGridParameters(2) = KARTE

      ELSE IF (ACTPAR .EQ. 'SPECIAL_INNER_POINTS') THEN
         RadiusGridParameters(3) = KARTE

      ELSE IF (ACTPAR .EQ. 'RMAX_IN_RSUN') THEN
         RMAX_IN_RSUN = .TRUE.

C*********** END OF SPECIAL BLOCK FOR X-RAYS  **************************

      ELSE IF (ACTPAR .EQ. 'HELIUM') THEN
         WRITE (0,*) 'NOT ALLOWED TO SPECIFY THE HELIUM ABUNDANCE!'
         GOTO 92

      ELSE
C***  Check if the card refers to an element abundance
         CALL FINDCHARGE (ACTPAR, NZ)
C***     NZ > 0 means that an element of this name is known
         IF (NZ .GT. 0) THEN
            IFOUNDELEMENT = 0
            DO K=1, NATOM
               IF (ACTPAR .EQ. ELEMENT(K)) THEN
                  IFOUNDELEMENT = 1
                  CALL SARGV (KARTE, 2, ACTPAR2)
                  READ (ACTPAR2, '(G20.0)', ERR=99) ABUND
                  IF (ABUND .LT. 0.) THEN
                     WRITE (0,*) 'NEGATIVE ABUNDANCE ENCOUNTERED'
                     GOTO 92
                  ENDIF
                  ABXYZ(K)=ABUND
                  IF (KARTE(21:30) .NE. ' ') THEN
                     ABMASS=.TRUE.
                  ELSE
                     ABNUMB=.TRUE.
                  ENDIF
                  EXIT
               ENDIF
            ENDDO

            IF (IFOUNDELEMENT .EQ. 0) THEN
               WRITE (0,*) 
     >             'ABUNDANCE GIVEN, BUT ELEMENT NOT FOUND IN DATOM'  
               GOTO 92
            ENDIF
         ENDIF


      ENDIF

      GOTO 10

 
  100 CONTINUE
      CLOSE (1)
C******************************************************************

C***  ERROR STOP IN CASE OF MIXED TYPES (BY NUMBER / MASS FRACTION) OF
C***  ABUNDANCE VALUES
      IF (ABMASS .AND. ABNUMB) THEN
        WRITE (0,*) 
     >  'ALL abundances must be given EITHER by mass OR by number!'
         GOTO 92
      ENDIF

C***  COMPUTATION OF THE RELATIVE ABUNDANCE (BY NUMBER OR MASS FRACTION)
C***  OF HELIUM
      ABREST=0.
      DO NA=1,MAXATOM
        IF (KODAT(NA) .GT. 0) THEN
          ABREST=ABREST+ABXYZ(KODAT(NA))
        ENDIF
      ENDDO

      IF (ABREST .GT. 1.) THEN
         WRITE(hCPR,'(A)') 'REL. ABUNDANCES add up to more than 100%'
         WRITE(hCPR,*) ' SUM of non-HE elements: ', ABREST
         STOP 'ERROR detected in DECSTAR'
      ENDIF

C***  Helium abundance is ALWAYS the complement to 100% 
      ABXYZ(KODAT(2)) = 1.-ABREST

C***  CHECK OF UNNECESSARY ATOMIC DATA 
      DO NA=1,MAXATOM
        IF (KODAT(NA) .GT. 0) THEN
          IF (ABXYZ(KODAT(NA)) .EQ. 0.) THEN
            WRITE (0,*) 'ELEMENT ', ELEMENT(KODAT(NA)), 
     >                  ' has ZERO abundance'
            STOP 'ERROR detected in DECSTAR'
          ENDIF
        ENDIF
      ENDDO
                 
      
C***  CONVERTING (POSSSIBLE) MASS FRACTIONS INTO FRACTIONAL ABUNDANCES
C***  BY NUMBER
      IF (ABMASS) THEN
         DENSUM=0.0
         DO NA=1,NATOM
           DENSUM=DENSUM+ABXYZ(NA)/ATMASS(NA)
         ENDDO
         DO NA=1,NATOM
           ABXYZ(NA)=ABXYZ(NA)/ATMASS(NA)/DENSUM
         ENDDO
      ENDIF

C***  Calculate Helium and Carbon Mass fraction from the number fractions
      SUM = .0
      DO NA=1, NATOM
        SUM = SUM + ABXYZ(NA)*ATMASS(NA)
      ENDDO
C***  Hydrogen mass fraction
      KHY = KODAT(1)
      IF (KHY == 0) THEN
        XHY = .0
      ELSE
        XHY = ABXYZ(KHY)*ATMASS(KHY) / SUM
      ENDIF
C***  Helium mass fraction: 
      KHE = KODAT(2)
      IF (KHE .EQ. 0) THEN 
        YHE = .0
      ELSE
        YHE = ABXYZ(KHE)*ATMASS(KHE) / SUM
      ENDIF
C***  Carbon mass fraction; 
      KC = KODAT(6)
      IF (KC .EQ. 0) THEN
        XC = .0
      ELSE
        XC = ABXYZ(KC)*ATMASS(KC) / SUM
      ENDIF
C***  Oxygen mass fraction;
      KO = KODAT(8)
      IF (KO .EQ. 0) THEN
        XO = .0
      ELSE
        XO = ABXYZ(KO)*ATMASS(KO) / SUM
      ENDIF
      
C***  Automatic definition of type (OB, WN or WC)
C***   this type may be used for the M-L relation
      IF (WRTYPE .EQ. '') THEN
         IF (XHY .GE. 0.45) THEN 
           WRTYPE = 'OB' 
         ELSE IF (XC .GE. 0.1) THEN 
           WRTYPE = 'WC' 
         ELSE
           WRTYPE = 'WN' 
         ENDIF
      ENDIF      

      IF (q < 0.) THEN
C***  Start approximation for q (number of free electrons per mass unit)
        !@todo: Set Q depending on type of star or s.th. like that
        !q = n_e / (AMU n_i) aber diese Groessen sind erst spaeter bekannt
        !using the following approximations:
        IF (XHY >= 0.5) THEN
          q = 0.88     !H-rich Of/WN Star (fully ionized solar comp.)
        ELSEIF (XHY >= 0.1) THEN
          q = 0.39    !WNL
        ELSEIF (XC >= 0.2) THEN
          IF (XO > 0.1) THEN
            q = 0.37  !WO
          ELSE
            q = 0.21  !WC 
          ENDIF
        ELSE
          q = 0.25    !WNE
        ENDIF
      ENDIF

C***  CHECK OF MISSING SPECIFICATIONS

C***  Rstar has not been specified in the CARDS
C***    but can be inferred from LOG GGRAV and MSTAR            
      IF (RSTAR < .0 .AND. XMSTAR > 0. .AND. GLOG > 0.) THEN
        RSTAR = SQRT(GCONST * XMSTAR * XMSUN / 10.**GLOG)
        RSTAR = RSTAR / RSUN
        RcalcCond = 4
      ENDIF

      bLRTcomplete = .FALSE.
      IF (TEFF < .0) THEN
        IF ((RSTAR > 0.) .AND. (XLOGL /= -99)) THEN
          !TEFF has not set in the CARDS file but can be inferred from L and R
          !Note: RSTAR should be still in RSUN at this point
          TEFF = ( 10**(XLOGL) / RSTAR**2 )**(0.25) * TEFFSUN
          WRITE (hCPR, *) 'NOTE: Temperature was infered from',
     >               ' radius and given luminosity: ', TEFF
          bLRTcomplete = .TRUE.
          LRTinput = 3
        ELSEIF ((RSTAR > 0.) .AND. (LcalcCond == 3)) THEN
          !In the rare case that MSTAR and EDDINGTON-GAMMA are given, L can be calculated:
          XLOGL = LOG10( GEDD * XMSTAR / ( 10**(-4.51) * q ) )
          TEFF = ( 10**(XLOGL) / RSTAR**2 )**(0.25) * TEFFSUN
          WRITE (hCPR, *) 'NOTE: Temperature was infered from',
     >              ' radius and calculated luminosity: ', TEFF
          bLRTcomplete = .TRUE.
          LRTinput = 5
        ELSE
          !Not enough information to calculate TEFF
          WRITE(hCPR,'(A)') 'EFFECTIVE TEMPERATURE NOT SPECIFIED'
          STOP 'ERROR'
        ENDIF
      ENDIF
      
      IF (RadiusGridParameters(1) .EQ. ' ') THEN
            WRITE(hCPR,'(A)') 'RADIUS GRID NOT SPECIFIED'
            STOP 'ERROR'
      ENDIF


      IF (.NOT. bLRTcomplete) THEN
        !This is only called if TEFF has been set in the CARDS file
        IF (RSTAR < .0) THEN
          !Rstar has not been specified in the CARDS and could not be inferred from log g and M
          IF (XLOGL .NE. -99.) THEN
            RSTAR = 10.**(0.5*XLOGL) / (TEFF/TEFFSUN)**2
            RcalcCond = 2
            LRTinput = 2
          ELSEIF (LcalcCond == 3) THEN
            !In the rare case that MSTAR and EDDINGTON-GAMMA are given, L can be calculated:
            XLOGL = LOG10( GEDD * XMSTAR / ( 10**(-4.51) * q ) )
            RSTAR = 10.**(0.5*XLOGL) / (TEFF/TEFFSUN)**2
            RcalcCond = 2
            LRTinput = 4
          ELSE
              WRITE(hCPR,'(A)') 
     >          'NEITHER RSTAR NOR LOG L OR EDDINGTON-GAMMA SPECIFIED'
              STOP 'ERROR'
          ENDIF
        ELSE 
          !Rstar has been preset in the CARDS
          IF (XLOGL .NE. -99.) THEN
              !TEFF, Rstar und L gesetzt => System ueberbestimmt
              WRITE(hCPR,'(A)') 'RSTAR AND LOG L SPECIFIED BOTH'
              WRITE(hCPR,'(A)') 'BUT ONLY ONE OF THEM IS INDEPENDENT!'
              STOP 'ERROR'
          ELSE
              XLOGL = 2 * ALOG10(RSTAR) + 4 * ALOG10(TEFF/TEFFSUN)
              LRTinput = 1
          ENDIF
        ENDIF
      ENDIF

      IF (MDOTINPUT < 0) THEN
            WRITE(hCPR,'(A)') 
     >        'NEITHER MDOT NOR RTRANS NOR LOG Q SPECIFIED'
            STOP 'ERROR'
      ENDIF

C***  Calculate XMDOT from transformed mass-loss rate
      IF (MDOTINPUT == 3 .AND. XMDOT == -99.) THEN
         XMDOT = XMDOTTRANS + ALOG10(VFINAL / 1000.)
     >                      - 0.5 * ALOG10(DENSCON_FIX)
     >                      - 0.75 * ( XLOGL - 6. )
      ENDIF
      IF (RTRANS /= -99.  .AND. QLOG /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, LOG Q AND RTRANS, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF (XMDOT /= -99.  .AND. QLOG /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, LOG Q AND MDOT, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF (RTRANS /= -99.  .AND. XMDOT /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, MDOT AND RTRANS, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF ( QLOG /= -99. .AND. RTRANS == -99.) THEN
         !Rt aus log Q berechnen
         RTRANS = (10**(-4.) / 2500.)**(2./3.) *
     >      (VFINAL * (10**(QLOG))**(2.) )**(-1./3.)
      ENDIF
      IF ( XMDOT == -99.) THEN
         !Mdot aus Rt-Def. berechnen
         XMDOT = ALOG10(VFINAL / 2500.) - 4.0 + (3./2.) * 
     >      ALOG10(RSTAR/RTRANS) - 0.5 * ALOG10(DENSCON_FIX)
      ENDIF

      IF (RMAX < .0) THEN
         WRITE(hCPR,'(A)') 'RMAX NOT SPECIFIED'
         STOP 'ERROR'
      ENDIF
      
      IF (RMAX_IN_RSUN) THEN
          RMAX = RMAX / RSTAR
          IF (RMAX < 1.) THEN
            WRITE(hCPR,'(A)') '*** FATAL ERROR: RMAX < RSTAR ***'
            WRITE(hCPR,'(A,G12.5)') ' RSTAR/SUN = ', RSTAR/RSUN
            WRITE(hCPR,'(A,G12.5)') ' RMAX/SUN  = ', RMAX/RSUN
            STOP 'ERROR IN DECSTAR'
          ENDIF
      ENDIF

      IF (RMAX < 1.) THEN
         WRITE(hCPR,'(A)') 'RMAX LESS THAN RSTAR'
         STOP 'ERROR'
      ENDIF
      
      IF (VDOP < .0) THEN
         WRITE(hCPR,'(A)') 'VDOP NOT SPECIFIED'
         STOP 'ERROR'
      ENDIF

      IF (TAUMAX <= .0) THEN
         WRITE(hCPR,'(A)') '**** WARNING: TAUMAX NOT SPECIFIED ****'
         WRITE(hCPR,'(A)')
     >     'Beware that PoWR might crash during JSTART->SPLINPO ',
     >     'if vmin is not high enough.'
      ENDIF
 
C***  OPTION 'OLDTEMP' OVERWRITES OTHER OPTIONS:
      IF (OLDTEMP) THEN
ccc        De-activated by ansander, 17-03-2021 
ccc         TTABLE=.FALSE.
         SPHERIC=.FALSE.
         TMIN=.0
      ENDIF

C***  OPTION 'TTABLE' OVERWRITES OTHER OPTIONS:
      IF (TTABLE) THEN
         SPHERIC=.FALSE.
         TMIN=.0
      ENDIF
      

 
C***  COMPUTATION OF THE MASS FLUX
      IF (FM .NE. .0) FM= 10.**FM
      IF (XMDOT .NE. .0) FM= 10.**(XMDOT+3.02) /RSTAR/RSTAR

     
      !Warning if full a_rad integration should be performed, but no good start is given
      IF (bFULLHYDROSTAT .AND. (iOldStratification == 0)
     >      .AND. (iOLDRAD == 0) .AND. RADGAMMASTART < 0.) THEN
        WRITE (hCPR,'(A)') 'WARNING: HYDROSTATIC INTEGRATION'
     >                          // ' will have a poor start. ' 
     >      // 'Set RADGAMMA-START card to improve the situation!'
      ENDIF
      
      !System g/M, geff, GEDD ueberbestimmt        
      IF ((GEDD >= 0.) .AND. (GEFFLOG > 0.) .AND. 
     >    ((MASSORIGIN == 1) .OR. (MASSORIGIN == 2))) THEN
        WRITE(hCPR,'(A)') 'LOG G, GEDD AND LOG GEFF SPECIFIED'
        STOP 'ERROR detected in DECSTAR'
      ENDIF

      

      IF ((GEFFLOG > 0.) .AND. (MASSORIGIN == 0)) THEN
         !Spezialfall: Nur LOG GEFF (und ggf. EDDINGTON GAMMA) gegeben, aber LOG G oder MSTAR nicht
         MASSORIGIN = 3
         IF (GEDD >= 0.) THEN
            GLOG = GEFFLOG - ALOG10 ( 1. - GEDD )          
         ELSEIF (RADGAMMASTART >= 0.) THEN
            GLOG = GEFFLOG - ALOG10 ( 1. - RADGAMMASTART )          
         ELSE
C***        (CAUTION: This should not be used if FULL HYDROSTATIC INTEGRATION is performed)
            XMSTAR = 10.**( GEFFLOG - 4.4371 ) * RSTAR**2. 
     >                 + 10.**(-4.51) * q * (10.**XLOGL)
            bCalcGLOGfromM = .TRUE.
         ENDIF
      ELSEIF ( (GEFFLOG > 0.) .AND.
     >         ((MASSORIGIN == 1) .OR. (MASSORIGIN == 2)) ) THEN
        !System stark bestimmt, Gamma implizit fest (=> keine Berechnung ueber q im weiteren Code)
        GEddFix = 1
        WRITE(hCPR,*)
     >    "WARNING: Eddington Gamma has been implicitly fixed"
        WRITE(hCPR,*)
     >    ' (Remove either g_eff, g_grav or Mstar from CARDS',
     >    ' file if you want to avoid this.)'
      ENDIF

C***  STELLAR RADIUS IN CM
      RSTAR=RSTAR*RSUN


C***  LOG G or MSTAR might not be specified; in that case, MSTAR is 
C***     calculated from the mass-luminosity-relation
      IF (MASSORIGIN == 0) THEN
C***     M-L relations: 
C***        OB type: Goetz (mandatory) - H-burner
C***        WN type: Goetz (default) or Langer '89 (optional) - He burner
C***        WC type: Langer (mandatory)
         IF (WRTYPE == 'WC' .OR. 
     >      (WRTYPE == 'WN' .AND. MFORM == 1) .OR.
     >       WRTYPE == 'WC') THEN
            CALL MLANGER (XMSTARG, TEFF, RSTAR, YHE, WRTYPE)         
            XMSTAR = XMSTARG / XMSUN
            MLRELATION='Langer'
         ELSE
            CALL MGOETZ (XMSTAR, TEFF, RSTAR, XHY, WRTYPE)
            XMSTARG = XMSTAR * XMSUN
            MLRELATION='Graefener'
         ENDIF
      ENDIF

      IF ((MASSORIGIN >= 2) .AND. (.NOT. bCalcGLOGfromM)) THEN
         !LOG G vorgegeben oder aus LOG GEFF und GEDD berechnet
         XMSTARG = 10.**GLOG * RSTAR * RSTAR / GCONST 
         XMSTAR = XMSTARG / XMSUN
      ELSEIF (RcalcCond .NE. 4) THEN
         !Masse direkt vorgegeben
         bCalcGLOGfromM = .TRUE.
      ENDIF

      IF (bCalcGLOGfromM) THEN         
         XMSTARG = XMSTAR * XMSUN          
         GLOG = ALOG10(GCONST * XMSTARG / RSTAR / RSTAR)
      ENDIF

C***  Calculate GEFFLOG or GEDD or both (if not specified from input)
      IF (GEDD >= 0.) THEN    
C***    EDDINGTON-GAMMA was specified:
        IF (GEDD > 1.) THEN
          WRITE(hCPR,*) "DECSTAR: EDDINGTON-GAMMA may not exceed 1"
          STOP 'FATAL ERROR in DECSTAR'
        ENDIF
        GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
        GEDDFIX = 2         !Gamma explizit fest
      ELSEIF (GEFFLOG >= 0.) THEN
C***    EDDINGTON-GAMMA implicitely specified (from LOG GEFF and LOG GRAV)     
        GEDD = 1. - 10**( GEFFLOG - GLOG )
      ELSEIF (RADGAMMASTART >= 0.) THEN
        GEFFLOG = ALOG10( (10.**GLOG) * (1. - RADGAMMASTART) ) 
      ELSE
C***    EDDINGTON-GAMMA must be estimated (using Thomson only)
C***    CAREFUL: This branch is also called if RADGAMMASTART: OLD is used
C***    Real GEFFLOG is calculated in PREP_GAMMARAD instead
        GEDD = 10**(-4.51) * q * (10.**XLOGL) / XMSTAR
        IF (GEDD >= 1.) THEN
          WRITE(hCPR,*) 'GEDD:', GEDD
          WRITE(hCPR,*) "DECSTAR: Initial GEDD guessing failed"
          WRITE(hCPR,*) "-- Star is beyond the Eddington limit!"
          WRITE(hCPR,'(A,F7.2)') '  Suggestion: Set MSTAR > ', 
     >                   10**(-4.51) * q * (10**XLOGL)
          STOP 'ERROR in DECSTAR'
        ENDIF
        GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
      ENDIF      

      IF (GEddFix > 0) THEN
        WRITE(hCPR,*) "Fixed Eddington Gamma value: ", GEDD
        IF (GEDD < 0.) THEN
          WRITE(hCPR,*) "DECSTAR: unphysical negative Eddington Gamma"
          STOP 'FATAL ERROR in DECSTAR'
        ENDIF
      ENDIF
  
C***  Consistency check between GEFF keyword and HYDROSTATIC INTEGRATION
      IF (THIN .AND. TRIM(GEFFKEY) /= '') THEN
        IF (TRIM(GEFFKEY) == 'AUTO') THEN
          WRITE(hCPR,*) "DECSTAR: Deprecated Syntax for LOG GEFF"
          WRITE(hCPR,*) 'When using the HYDROSTATIC INTEGRATION card, ' 
     >            // 'the intention of LOG GEFF should be specified.'
          WRITE(hCPR,*) 'Add RADFORCE=ELECTRON or RADFORCE=FULL as a '
     >            // 'keyword after LOG GEFF!'                  
          STOP 'ERROR in DECSTAR'
        ELSEIF (bFULLHYDROSTAT .AND. TRIM(GEFFKEY) /= 'RAD') THEN
          WRITE(hCPR,*) "DECSTAR: Inconsistent meaning for LOG GEFF"
          WRITE(hCPR,*) 'When using the FULL option on the HYDROSTATIC' 
     >      // ' INTEGRATION card, '
          WRITE(hCPR,*) 'LOG GEFF must also be flagged with '
     >      // ' RADFORCE=FULL.'
          STOP 'ERROR in DECSTAR'
        ELSEIF (.NOT. bFULLHYDROSTAT .AND. TRIM(GEFFKEY) /= 'THOM') THEN
          WRITE(hCPR,*) "DECSTAR: Inconsistent meaning for LOG GEFF"
          WRITE(hCPR,*) 'When using the HYDROSTATIC INTEGRATION card' 
     >      // ' without FULL option, '
          WRITE(hCPR,*) 'LOG GEFF must be flagged with '
     >      // ' RADFORCE=ELECTRON.'
          STOP 'ERROR in DECSTAR'
        ENDIF
      ENDIF
  
  
      IF (iOldStratification > 0) THEN
        THIN = .FALSE.      !Ignore HYDROSTATIC INTEGRATION in WRSTART if OLD STRAT is used
      ENDIF
  
      RETURN


C***  ERROR EXITS ****************************************************

   92 CONTINUE
      WRITE (0,*)'DECSTAR: ERROR WHILE DECODING THE FOLLOWING LINE:'
      WRITE (0,*) KARTE
      STOP 'ERROR'

   95 WRITE (hCPR,'(A)') 'DECSTAR: OLD MODEL REQUIRED BUT NOT FOUND!'
      WRITE (hCPR,'(A)') 'AN OLD MODEL IS REQUIRED DUE TO THE LINE: '
      WRITE (hCPR,'(A)') KARTE(:IDX(KARTE))
      STOP 'FATAL ERROR IN DECSTAR'

   96 WRITE (hCPR,'(A)') 'DECSTAR: OLD MODEL DOES NOT CONTAIN MDOT!'
      WRITE (hCPR,'(A)') 'AN OLD MODEL IS REQUIRED DUE TO THE LINE: '
      WRITE (hCPR,'(A)') KARTE(:IDX(KARTE))
      STOP 'FATAL ERROR IN DECSTAR'
      
   97 WRITE (0,'(A)') '*** ERROR: PARAMETER MISSING'
      WRITE (0,'(A)') 'THE ERROR OCCURED IN THE FOLLOWING LINE: '
      WRITE (0,'(A)') KARTE(:IDX(KARTE))
      STOP 'ERROR IN DECSTAR'

   98 WRITE (0,'(A)') 
     >  '*** ERROR: THE FOLLOWING STRING COULD NOT BE DECODED AS A '
     > // 'FLOATING POINT NUMBER:', ACTPAR2(:IDX(ACTPAR2)), 
     >    'THE ERROR OCCURED IN THE FOLLOWING LINE:', KARTE
      STOP 'ERROR IN DECSTAR'

   99 WRITE (0,'(A)') 
     >  '*** ERROR: THE FOLLOWING STRING COULD NOT BE DECODED AS A '
     > // 'FLOATING POINT NUMBER:', ACTPAR(:IDX(ACTPAR)), 
     >    'THE ERROR OCCURED IN THE FOLLOWING LINE:', KARTE
      STOP 'ERROR IN DECSTAR'

      END
      SUBROUTINE DECVELPAR(KARTE, VFINAL, VMIN, BETA, RMAX)
C***  Decodes the VELPAR line form the CARDS file      

      IMPLICIT NONE

      CHARACTER(40) :: TRYPAR
      CHARACTER(40), DIMENSION(20) :: CURPAR
      CHARACTER(100) :: KARTE
      REAL :: VFINAL, VMIN, BETA, RMAX

      INTEGER :: NPAR, i, IERR
      
      LOGICAL :: bOldDecode 
      LOGICAL, DIMENSION(4) :: bParamFound

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
     
      bOldDecode = .FALSE.
      CALL SARGC (KARTE, NPAR)
      IF (NPAR < 5) THEN
        WRITE (hCPR,'(A)') '*** VELPAR: NOT ENOUGH PARAMETERS'
        STOP '*** FATAL ERROR WHILE DECODING VELPAR CARDS-LINE'
      ENDIF

      !New decoding => allows flexible format and modern syntax
      bParamFound = .FALSE.
      DO i=1, NPAR
        CALL SARGV(KARTE,i,CURPAR(i))
      ENDDO
      IF (NPAR > 2) THEN
        DO i=2, NPAR 
         SELECTCASE (CURPAR(i))
          CASE ('VFINAL')          
            IF (NPAR >= (i+1)) THEN
              TRYPAR = CURPAR(i+1)
              IF (TRYPAR == '(KM/S)') THEN
                IF (NPAR >= (i+2)) THEN
                  TRYPAR = CURPAR(i+2)
                ELSE
                  GOTO 92
                ENDIF
              ENDIF
              READ (TRYPAR, '(F10.0)', IOSTAT=IERR, ERR=92) VFINAL
              IF (IERR == 0) THEN
                bParamFound(1) = .TRUE.
              ENDIF      
            ENDIF
          CASE ('VMIN')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) VMIN
              IF (IERR == 0) THEN
                bParamFound(2) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('BETA')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) BETA
              IF (IERR == 0) THEN
                bParamFound(3) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('RMAX')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) RMAX
              IF (IERR /= 0) THEN
                WRITE(hCPR,'(A)') '*** DECVELPAR: CANNOT READ RMAX'
                STOP 'ERROR'
              ELSE
                bParamFound(4) = .TRUE.
              ENDIF                  
            ENDIF
         ENDSELECT
        ENDDO
      ENDIF

      DO i=1, 4 
        !One or more parameters have not been found => switch to old decoding
        IF (.NOT. bParamFound(i)) THEN
          WRITE (hCPR,*) '*** DECVELPAR: Old VELPAR decoding used'
          bOldDecode = .TRUE.
        ENDIF
      ENDDO

      IF (bOldDecode) THEN
         READ (KARTE,19,ERR=99) VFINAL,VMIN,BETA,RMAX
   19    FORMAT(22X,F7.0,5X,F6.0,5X,F4.0,5X,F6.0)
      ENDIF
      
      RETURN
      
C***  FATAL ERROR CODES      
      
   92 WRITE(hCPR,'(A)') '*** DECVELPAR: CANNOT READ VFINAL IN:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'
      
   99 WRITE (hCPR,*)
     >   'DECVELPAR: ERROR WHILE DECODING THE FOLLOWING CARDS-LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'
      
      END
            FUNCTION DELTAGR(R)
C***  DIFFERENCE OF VELOCITY GRADIENTS FROM BOTH VELOCITY LAWS: 
C***  THE BETA LAW (OUTER REGION) AND THE EXPONENTIAL LAW (INNER REGION)
C***  THIS ROUTINE SERVES FOR ESTIMATING THE CONNECTION POINT BETWEEN BOTH
C***  REGIONS AND IS CALLED FROM SUBROUTINE INITVEL, MAIN PROGRAM WRSTART
C***  BRANCH FOR BETA .LE. 0 (SQRT-LOG-LAW) REMOVED (wrh 14-Aug-2009) 

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE, 
     >       BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  New parametrization of beta law, wrh  5-Mar-2001 
C***  Two-beta-law: no change, since no contribution from second term
      RP = VPAR2 + R

      GRADOUT = VPAR1*BETA*(1.-1./RP)**(BETA-1.) /(RP*RP)      

C***  Prevent overflow of exp
      ARG = AMIN1 (100., (R-1.)/HSCALE)       
      GRADIN = VMIN * EXP(ARG) / HSCALE

      DELTAGR= GRADOUT - GRADIN 

      RETURN
      END
      FUNCTION DELTAGRTHIN(R,ND,RADIUS,VELO)
C***  DIFFERENCE OF VELOCITY GRADIENTS FROM BOTH VELOCITY LAWS: 
C***  THE BETA LAW (OUTER REGION) AND THE EXPONENTIAL LAW (INNER REGION)
C***  THIS ROUTINE SERVES FOR ESTIMATING THE CONNECTION POINT BETWEEN BOTH
C***  REGIONS AND IS CALLED FROM SUBROUTINE INITVEL, MAIN PROGRAM WRSTART
C***  BRANCH FOR BETA .LE. 0 (SQRT-LOG-LAW) REMOVED (wrh 14-Aug-2009) 

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      REAL :: DELTAGRTHIN

      INTEGER :: ND

      REAL, DIMENSION(ND) :: RADIUS, VELO

      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2,
     >        R, RP, DVDR, GRADIN, GRADOUT, Vdummy

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE, 
     >       BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C      COMMON/THINGRAD/ GRADIN

      CALL SPLINPOX(Vdummy, R, VELO, RADIUS, ND, DFDX=DVDR)
C      GRADIN = DVDR

C***  New parametrization of beta law, wrh  5-Mar-2001 
C***  Two-beta-law: no change, since no contribution from second term
      RP = VPAR2 + R

      GRADOUT = VPAR1*BETA*(1.-1./RP)**(BETA-1.) /(RP*RP)      

      DELTAGRTHIN = GRADOUT - DVDR 

      RETURN
      END
      SUBROUTINE DTDR(N,R,T,TPRIME)
C*******************************************************************************
C***  THIS SUBROUTINE HELPS TO CALCULATE THE TEMPERATURE STRUCTURE OF A
C***  SPHERICAL, GREY ATMOSPHERE.
C***  IT IS A FORMAL PARAMETER OF THE IMSL-ROUTINE"DVERK" AND CALLED ONLY
C***  FROM THERE.
C***  THE DIFFERENTIAL EQUATION IS WRITTEN: Y' = F(X,Y)
C***  HERE:  DT/DR = FUNCTION OF (R,T)
C***  THE EDDINGTON FACTOR F IS ASSUMED TO BE 1/3 INSIDE THE RADIUS R1, I.E.
C***  WHERE TAU=1.
C***  OUTSIDE THE RADIUS R13, I.E. WHERE TAU=1/3, WE ASSUME GEOMETRICAL DILUTION
C***  DILUTION OF THE RADIATION FIELD AS FOR A UNIFORM RADIATING SPHERE
C***  OF RADIUS R23
C***  F IS INTERPOLATED IN THE INTERJACENT REGION
C*******************************************************************************
 
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
C***  OPAROSS FROM SUBROUTINE GREY
      COMMON /COMDTDR/ OPAMEAN,R23COM,R1COM,R13COM
 
      IF (R .LE. R1 COM) THEN
C***     OPTICALLY THICK REGION
          FEDDI=.33333333333333
          BFAC=0.
      ELSE IF (R .LT. R13COM) THEN
C***  INTERPOLATED REGION AROUND R23, I.E. BETWEEN TAU=1/3 AND TAU=1
C***  F1 = FEDDI AT R1
      F1=0.333333333333
C***  F2 = FEDDI AT R13
      R2=R13COM*R13COM
         SMUE=SQRT(1.-R23COM*R23COM/R2)
      F2=((SMUE+1.)*SMUE+1.)/3.
C***  F3 = FPRIME AT R1
      F3=0.
C***  F4 = FPRIME AT R13
      F4=(2.+1./SMUE)*R23COM*R23COM/R2/R13COM/3.
C***  POLYNOMIAL COEFFICIENTS, SEE SUBROUTINE CUBIC
      D1=1./(R13COM-R1COM)
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
C***  EVALUATION OF THE POLYNOMIAL
      FEDDI=P1*(R-R1COM)**3 + P2*(R-R1COM)
     +     +P3*(R13COM-R)**3+P4*(R13COM-R)
C***  ... AND ITS DERIVATIVE
      FPRIME=3.*P1*(R-R1COM)**2 + P2
     -      -3.*P3*(R13COM-R)**2-P4
          BFAC=(3.*FEDDI-1.)/R+FPRIME
      ELSE
          R2=R*R
         SMUE=SQRT(1.-R23COM*R23COM/R2)
          FEDDI=((SMUE+1.)*SMUE+1.)/3.
          FPRIME=(2.+1./SMUE)*R23COM*R23COM/R2/R/3.
          BFAC=(3.*FEDDI-1.)/R+FPRIME
      ENDIF
      TT=TEFF/T
      TT3=TT*TT*TT
      TPRIME=-BFAC*T/4./FEDDI-OPAMEAN/R/R*TEFF/16./FEDDI*TT3
      RETURN
      END
      SUBROUTINE FEDAT (ROUTINE, INDEXMAX, NFEREADMAX, IONLOW, IONTOP,
     &                  MAXATOM, NDIM, MAXIND, MAXKONT, NATOM,      
     &                  N, LASTFE, LASTKON, LASTINDAUTO, MAXFEIND, 
     &                  EINST, SIGMAFE, INDRB, INDRF, IFENUP, 
     &                  IFELOW, INDNUP, INDLOW, KONTNUP, KONTLOW,
     &                  LEVEL, ELEMENT, SYMBOL, ATMASS, STAGE,
     &                  ELEVEL, WEIGHT, EION, NCHARG, NOM, KODAT,
     &                  NFIRST, NLAST, IFRBSTA, IFRBEND, FEDUMMY,
     &                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, KEYCBB)

c!!!!!! Folgende parameter wurden entfernt: 
C!!!    CBFC, BOUND, EINSTINT, COCOFE, NCOMAX, NCO
c!!!    folgende Parameter sind neu: MAXFEIND, FEDUMMY
C!!!    umbenannt wurden: NMAX -> NFEREADMAX, DUMMY -> FEDUMMY

C **********************************************************************
C ***
C *** CALLED BY: SUBROUTINE DATOM
C ***    READS ALL RELEVANT ATOMIC DATA FOR IRON GROUP LINE BLANKETING
C ***    FROM A MASS-STORAGE FILE CREATED BY THE IRON-PACKAGE (TAPE 21)
C ***
C **********************************************************************

      IMPLICIT NONE

C***  Local dimensions:
C***  Maximum number of bound-bound transitions within one ion
c      PARAMETER ( NBBMAX  = 400 )
      INTEGER, PARAMETER :: NBBMAX = 999  !taken from CREATMS (from blanket program)
C***  NIONMAX must cover the number of Iron Ionization stages in the FEDAT file
      INTEGER, PARAMETER :: NIONMAX = 27
C***  NROMMAX: Roman Numbers encoded: cf. DATA-Statement for ROMNUM 
      INTEGER, PARAMETER :: NROMMAX = 27
C***  Local arrays
      INTEGER, DIMENSION(NIONMAX) :: NB, NTRA, NTRB
      INTEGER, DIMENSION(NBBMAX) :: NX_A, NX_B
      CHARACTER(2) :: CIONBUFFER
      CHARACTER(5), DIMENSION(NROMMAX) :: ROMNUM 
      CHARACTER(16), DIMENSION(NBBMAX) :: NAMARRAY
      CHARACTER(16) :: NAMBUFFER
      CHARACTER(LEN=8) :: NAME
      CHARACTER(LEN=24) :: GENER
      CHARACTER(LEN=40) :: IONNAME(NIONMAX)
      CHARACTER(LEN=3) :: ISTR
      
C***  Maximum number of superlevels within the same ion
      INTEGER, PARAMETER :: MAXLEVEL = 100
      CHARACTER(LEN=8), DIMENSION(MAXLEVEL) :: LEVNAMES

      LOGICAL :: BEXTEND

      INTEGER, INTENT(IN) :: IONLOW, IONTOP, NATOM, NDIM,
     >                       MAXIND, MAXATOM, MAXKONT, MAXFEIND, 
     >                       INDEXMAX, NFEREADMAX, LASTINDAUTO
      INTEGER, INTENT(INOUT) :: N, LASTFE, LASTKON
      
C***  Formal Parameters
      CHARACTER(LEN=*)  ROUTINE
      CHARACTER(LEN=2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(LEN=10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(LEN=10), DIMENSION(NDIM) :: LEVEL
      
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      INTEGER, DIMENSION(NDIM) :: NCHARG, NOM
      REAL, DIMENSION(NDIM) :: EION, ELEVEL, WEIGHT
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW

C***  The following arrays are ONLY used in STEAL, WRSTART 
C***  -> not necessary to define these arrays and MAXIND in other calls 
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      CHARACTER(LEN=4), DIMENSION(MAXIND) :: KEYCBB

C***  IRON-SPECIFIC ARRAYS; MAXFEIND = MAX. NUMBER OF IRON SUPERLINES
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFENUP, IFELOW,
     >                                IFRBSTA, IFRBEND
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      
      REAL :: SIGMAUL, SIGMALU, SIGMAINTCUR,
     >        WLOW, XF, XFOLD, DSIGMA, XLAM, XLOGSTEP, DXFE, 
     >        XLAM0FE, VDOPFE, XBAK, WNUP, XLAMCMIND
      INTEGER :: INDEX, INDEXS, LOWION, NUPION, I, J, K, NION,
     >           NDATA, LEVCOUNT, ILAM, NUP, LOW, INDSTA,
     >           INDEND, INDSELF, IND, IFREQRBSTA, IFREQRBEND,
     >           NOLD, KONT, IADR, MAXADR, IERR, NTEST, IBAK,
     >           NREADIN, IERRLEVNAM, KZEROS, KZEROE
      
      LOGICAL :: bFEINFO, bFEULSEP, bINTRA

C***  Constants:
      REAL, PARAMETER :: CLIGHT = 2.99792458E10   !C in cm/s
      REAL, PARAMETER :: PI8 = 25.1327412288      !PI8 = 8*PI
      REAL, PARAMETER :: FSIG = 2.6540E-2         !PI*e^2/m/c in CGS units
      
C***  Roman Numbers for Level names
      DATA ROMNUM / 'I....', 'II...', 'III..', 'IV...', 'V....', 
     >              'VI...', 'VII..', 'VIII.', 'IX...', 'X....',
     >              'XI...', 'XII..', 'XIII.', 'XIV..', 'XV...',
     >              'XVI..', 'XVII.', 'XVIII', 'XIX..', 'XX...',
     >              'XXI..', 'XXII.', 'XXIII', 'XXIV.', 'XXV..',
     >              'XXVI.', 'XXVII' /

      CALL OPENMS(21,IADR,MAXADR,1,IERR)

C *** READ GENERAL PARAMETERS OF GENERIC ION
      CALL READMS (21, GENER, 3, 'GENERIC ', IERR)
      READ (GENER(15:16), '(F2.0)') STAGE(NATOM)
      READ (GENER(17:24), '(F8.5)') ATMASS(NATOM)
C***  NION = Number of Iron Ionization stages in the FEDAT Data file 
      CALL READMS (21, NION,    1,      'NION    ', IERR)
      IF (NION .GT. NIONMAX .OR. IONTOP .GT. NIONMAX) THEN
         WRITE (0,'(A)') '*** Local dimension NIONMAX insufficient !' 
         WRITE (0,*) NION, IONTOP, NIONMAX
         STOP            '*** ERROR STOP IN FEDAT *********'
      ENDIF

C *** READ PROPERTIES AND NO. OF TRANSITIONS FOR ALL IONIZATION STAGES
      CALL READMS (21, IONNAME, NION*5, 'IONNAME ', IERR)

      CALL READMS (21, NTRA,    NION,   'NTRA_A  ', IERR)
      CALL READMS (21, NTRB,    NION,   'NTRA_B  ', IERR)

C *** READ PARAMETERS OF FREQUENCY GRID
      CALL READMS (21, VDOPFE,  1,      'VDOPP   ', IERR)
      CALL READMS (21, DXFE,    1,      'FSTEP   ', IERR)
      CALL READMS (21, XLAM0FE, 1,      'XLAMNULL', IERR)

C *** INITIALIZE COUNTERS
      INDEX    = 0
      LEVCOUNT = N
      IND      = 0
      KONT     = LASTKON


C *** READ NUMBER OF SUPERLEVELS PER IONIZATION STAGE
        CALL READMS (21, NB, NION, 'NLEV    ', IERR)

C **********************************************************************
C *** LOOP OVER IONIZATION STAGES
C **********************************************************************
      NOLD = N
      DO 10 I=IONLOW, IONTOP

C***    IONLOW == 0 is a placeholder for using the full levels for neutral Fe
        IF (I == 0) CYCLE 

        IF (IONTOP .GT. NION+1) THEN
           WRITE (0, '(A,I3)') 
     >         '*** IONTOP stage requested from DATOM file:', IONTOP 
           WRITE (0, '(A,I3)') 
     >         '*** Highest stage +1 available in FEDAT:', NION+1  
           STOP '*** ERROR DETECTED BY FEDAT' 
        ENDIF

C***    Note: A further ionization stage NION+1, which is not in the data, 
C***        is added if requested. This EXTEND stage has only one level
        BEXTEND = I .EQ. NION+1

C ***   REDUCTION TO 1 LEVEL FOR HIGHEST AND LOWEST IONISATIION STAGE
        IF ((I .EQ. IONLOW).OR.(I .EQ. IONTOP) .OR. BEXTEND) NB(I) = 1

        IF (N+NB(I) .GT. NDIM) THEN
           WRITE (0,'(A)') '*** Dimension NDIM insufficient !' 
           WRITE (0,'(A,I4)') '*** Present value NDIM = ', NDIM
           WRITE (0,'(A,I4)') '*** Required value = ', N+NB(I)
           WRITE (0,*)  I, N, NB(I)
           STOP            '*** ERROR STOP IN FEDAT *********'
        ENDIF

        IF (NB(I) .GT. MAXLEVEL) THEN
           WRITE (0,'(A)') '*** Dimension MAXLEVEL insufficient !' 
           WRITE (0,'(A)') '*** max. number of superlevels in one ion' 
           WRITE (0,'(A,I4)') '*** Present value = ', MAXLEVEL
           STOP            '*** ERROR STOP IN FEDAT *********'
        ENDIF

        IF (.NOT.BEXTEND) THEN

C***      Check dimension: Max number of bound-bound transitions within present ion
           IF (NTRA(I) .GT. NBBMAX .OR. NTRB(I) .GT. NBBMAX) THEN
              WRITE (0,'(A)') '*** Dimension NBBMAX insufficient !' 
              STOP '*** ERROR STOP IN FEDAT'
           ENDIF

C***       STORE MEAN ENERGIES AND STATISTICAL WEIGHTS        
           NAME = 'ELEV' // IONNAME(I)(2:4) 
           CALL READMS (21, ELEVEL(LEVCOUNT+1), NB(I), NAME, IERR)

           NAME = 'WEIG' // IONNAME(I)(2:4) 
           CALL READMS (21, WEIGHT(LEVCOUNT+1), NB(I), NAME, IERR)

C***       Levelnames with parity (if used) - new since 18-Jan-2016
           NAME = 'LEVN' // IONNAME(I)(2:4)
           CALL READMS (21, LEVNAMES, NB(I), NAME, IERRLEVNAM)

C***       READ NUMBER OF DATA-POINTS PER CROSS-SECTION (PRESENT ION I)
           NAME = 'N' // IONNAME(I)(2:4) // '_A  '
           CALL READMS (21, NX_A, NTRA(I), NAME, IERR)

           NAME = 'N' // IONNAME(I)(2:4) // '_B  '
           CALL READMS (21, NX_B, NTRB(I), NAME, IERR)

        ELSE
           ELEVEL(LEVCOUNT+1) = 0.
           WEIGHT(LEVCOUNT+1) = 1.
        ENDIF

        
C***  CREATE SUPERLEVEL NAMES, READ CHARGES AND IONIZATION ENERGIES      
        DO J=1,NB(I)
          N = LEVCOUNT+J
          IF (.NOT. BEXTEND) THEN
            READ(IONNAME(I)(2:3),'(I2)') NCHARG(N)        
            IF (J .EQ. 1) READ(IONNAME(I)(17:24),'(F8.0)') EION(N)
          ELSE
C***        For the extended Level no IONNAME exists
            READ(IONNAME(I-1)(2:3),'(I2)') NCHARG(N)
            NCHARG(N)=NCHARG(N)+1
            IF (J .EQ. 1) EION(N)=0.
          ENDIF
          NOM(N) = NATOM


          IF (IERRLEVNAM /= -10 .AND. .NOT. BEXTEND) THEN
C***        level names are already prepared in the FEDAT file
            LEVEL(N) = SYMBOL(NATOM) // LEVNAMES(J)
          ELSE
C***        use default level names if no predefined names are available
            IF (NCHARG(N)+1 .GT. NROMMAX) THEN
              WRITE (0,'(A)') '*** Roman Number for Ion. stage not known' 
              STOP            '*** ERROR STOP IN FEDAT *********'
            ENDIF
            LEVEL(N) = SYMBOL(NATOM) // ' ' // ROMNUM(NCHARG(N)+1) // '.'
            WRITE (LEVEL(N)(9:10),'(I2)') J
            IF (J .LE. 9) LEVEL(N)(9:9) = '.'
          ENDIF 
        ENDDO

        IF (N .GT. NOLD) THEN
          NFIRST(NATOM) = NOLD+1
          NLAST(NATOM) = N
        ENDIF
                
        NAME = 'A' // IONNAME(I)(2:4) // 'NAM ' 
        CALL READMS (21, NAMARRAY, 2*NTRA(I), NAME, IERR)
                
C**********************************************************************
C***    STORE RBB TRANSITION-DATA IN ONE-DIMENSIONAL ARRAY >>SIGMAFE<<
C***    LOOP OVER ALL BOUND-BOUND TRANSITIONS
        DO 20 J=1,NTRA(I)

C***      NO BB-TRANSITIONS FOR HIGHEST AND LOWEST IONISATION STAGES
          IF ((I .GE. IONTOP).OR.(I .LE. IONLOW)) GOTO 20

C***      READ LEVEL NUMBERS ASSOCIATED WITH TRANSITION          
          CIONBUFFER = NAMARRAY(J)(3:4)
          READ (UNIT=CIONBUFFER,FMT='(I2)') LOWION 
          CIONBUFFER = NAMARRAY(J)(6:7)
          READ (UNIT=CIONBUFFER,FMT='(I2)') NUPION
          LOW=LOWION+LEVCOUNT
          NUP=NUPION+LEVCOUNT

C***      OMIT TRANSITION IF LOW=NUP
          IF (LOW .EQ. NUP) GOTO 20

C***      BB-TRANSITION INDEX FOR IRON LINES (STARTING FROM 1)
          IND = IND + 1        
          IF (IND. GT. MAXFEIND) THEN
           WRITE (0,'(A)') '*** Dimension MAXFEIND insufficient !' 
           STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF     

C***      CREATE POINTER TO STARTING INDEX OF RBB TRANSITION-DATA         
          INDSTA = INDEX+1
          NDATA = NX_A(J)
          IF (NDATA+2 .GT. NFEREADMAX) THEN
              WRITE (0,'(A)') '*** Dim. NFEREADMAX insufficient for b-b!' 
              WRITE (0,'(A, I10)') 
     >          '*** dimensioned: NFEREADMAX = ', NFEREADMAX
              WRITE (0,'(A, I10)') 
     >          '*** required   : NFEREADMAX = ', NDATA+2
              STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          INDEND = INDSTA+NDATA
          IF (INDEND .GE. INDEXMAX) THEN 
              WRITE (0,'(A)') '*** Dimension INDEXMAX insufficient !' 
              WRITE (0,'(A, I10)') 
     >          '*** dimensioned: INDEXMAX =', INDEXMAX
              WRITE (0,'(A, I10)') 
     >          '*** required   : INDEXMAX =', INDEND
              STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          
C ***   STORE LEVEL-NUMBERS IN INDEX-ARRYS          
          IFENUP(IND)=NUP
          IFELOW(IND)=LOW
          
C ***   READ TRANSITION DATA          
          CALL COUNT(J, ISTR)
          NAME = 'A' // IONNAME(I)(2:4) // ISTR 
          CALL READMS (21, FEDUMMY, NDATA+2, NAME, IERR)

C ***   STORE FREQUENCY INDICES
          IFREQRBSTA = - INT(FEDUMMY(2))
          IFREQRBEND = - INT(FEDUMMY(1))

C ***   STORE CROSS-SECTIONS IN ARRAY >>SIGMAFE<<          
          XLOGSTEP = ALOG10(1. + VDOPFE*1.E5*DXFE/CLIGHT)
          SIGMAINTCUR = 0.
          KZEROS = 0
          KZEROE = 0
          DO K=1,NDATA
C***        Calculation of Lambda and Nu in cgs
             ILAM = IFREQRBSTA + K - 1
             XLAM = XLAM0FE*1.E-8 * 10.**(ILAM*XLOGSTEP)
             XF   = CLIGHT / XLAM
C ***       CROSS-SECTION
               SIGMAFE(INDEX+K) = FEDUMMY(NDATA-K+3) 
C***           Determine zero cross section regions at start and end
               IF (SIGMAFE(INDEX+K) <= 0.) THEN
                 IF (KZEROS >= 0) KZEROS = KZEROS + 1
                 KZEROE = KZEROE + 1
               ELSE 
C***             Non-zero cross section found:
                 IF (K == 1) THEN
C***               The first entry in the array is already non-zero.
C***               Therefore stop all further increasements of KZEROS by
C***               setting the counter to a negative value:
                   KZEROS = -1
                 ELSEIF (KZEROS > 0) THEN
C***               To stop further increasing of KZEROS after the first time
C***               this has occured, multiply the result with -1.
                   KZEROS = -1. * KZEROS
                 ENDIF
C***             Reset KZEROE since the formet part with zero cross-section
C***             was definately not at the end of the cross section array.
                 KZEROE = 0
               ENDIF
 
C ***          INTEGRATION OF EINSTEIN-COEFFICIENT AND NORM FOR SIGMA
               IF (K .GT. 1) THEN
                 DSIGMA = (SIGMAFE(INDEX+K)+SIGMAFE(INDEX+K-1))/2.
                 SIGMAINTCUR = SIGMAINTCUR + DSIGMA * (XFOLD - XF)
ccc the following statement is deactivated in libcr_cl version 16-Feb-1999
C***           NU^2 - Term is now accounted for
ccc             XNU = XLAMCMIND/XLAM
ccc                XNUMID = (XNU+XNUOLD)/2.
ccc                XNUMID2= XNUMID*XNUMID
ccc    IMPORTANT: --- consitent change required in CMFFEOP
ccc               --- new linking of both, COLI *and* STEALCL 
ccc                SIGMAINT(IND) = SIGMAINT(IND)
ccc     >                          + XNUMID2 * DSIGMA * (XFOLD - XF)
             ENDIF

             XFOLD = XF
ccc             XNUOLD= XNU
          ENDDO


C***      Band-Band transition 
          INDRB(IND) = INDSTA
            
          IFRBSTA(IND) = IFREQRBSTA
          IFRBEND(IND) = IFREQRBEND
C***      Remove empty cross sections regions from pointer range
C***      (leaves only at maximum one zero entry at beginning and end)
!             IF (ABS(KZEROS) > 0.) THEN
!               DO K=1, ABS(KZEROS)-1
!                 IF (SIGMAFE(INDEX+K) > 0.) STOP 'FATAL: KREDUCINGS FAILED!'
!               ENDDO
! c              WRITE (0,*) 'STA: ', IFREQRBSTA, KZEROS
!               IFRBSTA(IND) = IFREQRBSTA + ABS(KZEROS) - 1
!               INDRB(IND) = INDSTA + ABS(KZEROS) - 1
!             ENDIF
!             IF (ABS(KZEROE) > 0.) THEN
!               DO K=NDATA, NDATA-ABS(KZEROE)+1, -1
!                 IF (SIGMAFE(INDEX+K) > 0.) STOP 'FATAL: KREDUCINGE FAILED!'
!               ENDDO
! c              WRITE (0,*) 'END: ', IFREQRBEND, KZEROE
!               IFRBEND(IND) = IFREQRBEND - ABS(KZEROE) + 1
!             ENDIF
            
            
          SIGMAINT(IND) = SIGMAINTCUR                  
                    
          XLAMCMIND = 1./(ELEVEL(NUP) - ELEVEL(LOW))
          WLOW = WEIGHT(LOW)
          WNUP = WEIGHT(NUP)

          EINST(NUP,LOW) = SIGMAINTCUR * 
     >                       PI8*WLOW/WNUP/(XLAMCMIND*XLAMCMIND)
            
                    
C***      enhance index for next cross section reading     
          INDEX = INDEND

 20    CONTINUE
 
 21    CONTINUE

C **********************************************************************
C *** STORE RBF TRANSITION-DATA IN ARRAY >>EINST<< 
C *** LOOP OVER ALL BOUND-FREE TRANSITIONS
        NAME = 'B' // IONNAME(I)(2:4) // 'NAM'
        CALL READMS (21, NAMARRAY, 2*NTRB(I), NAME, IERR)

        DO 30 J=1, NTRB(I)

C ***   NO BF-TRANSITION FOR HIGHEST IONISATION STAGE
          IF (I .EQ. IONTOP) GOTO 30

C ***   READ LEVEL NUMBERS ASSOCIATED WITH TRANSITION          
          READ (NAMARRAY(J)(6:7),'(I2)') LOWION
          LOW = LOWION + LEVCOUNT
          NUP = LEVCOUNT + NB(I) + 1

C ***   MODEL ION REDUCED TO ONE LEVEL?
          IF (LOWION .GT. NB(I)) GOTO 30
          
C ***   INCREASE CONTINUUM-INDEX (ADDING UP TO "NORMAL" CONTINUA)
          KONT = KONT + 1
          IF (KONT. GT. MAXKONT) THEN
           WRITE (0,'(A)') '*** Dimension MAXKONT insufficient !' 
           STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF     

C ***   STORE LEVEL-NUMBERS IN INDEX-ARRAYS          
          KONTNUP(KONT) = NUP
          KONTLOW(KONT) = LOW
          
C ***   CREATE POINTER AT STARTING INDICES OF RBF TRANSITION-DATA         
          IF (J.EQ.1) THEN
            INDRF(KONT) = INDEX + 1
          ELSE
            INDRF(KONT) = NX_B(J-1) + 2 + INDRF(KONT-1)
          ENDIF
          
C ***   READ TRANSITION DATA          
          CALL COUNT(J, ISTR)
          NAME = 'B' // IONNAME(I)(2:4) // ISTR 
C ***   THE LAST INDEX IS THE COMPOSED CROSS SECTION AT THE EDGE 
          NREADIN = NX_B(J)+3
          IF (NREADIN .GT. NFEREADMAX) THEN
             WRITE (0,'(A)') '*** Dim. NFEREADMAX insufficient for b-f!' 
             STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          CALL READMS (21, FEDUMMY, NREADIN, NAME, IERR)

C ***     STORE COMPOUND THRESOLD CROSS-SECTIONS IN "EINST" [10**-18 CM**2]
          EINST(LOW, NUP) = FEDUMMY(NREADIN) * 1.E18

C ***     STORE CROSS-SECTIONS IN ARRAY >>SIGMAFE<<          
          DO K=1, NX_B(J)+2
            SIGMAFE(INDEX+K) = FEDUMMY(K)
          ENDDO
          INDEX = INDEX+NX_B(J) + 2

 30     CONTINUE
C***    END-OF LOOP OVER BOUND-FREE TRANSITIONS  **********

C***    ACTUALIZE LEVEL-COUNTER        
        LEVCOUNT = LEVCOUNT + NB(I)

 10   CONTINUE
C *** END OF LOOP OVER IONIZATION STAGES
C **********************************************************************

C ***  REORDER BB-TRANSITIONS TO INCREASING FREQUECY INDEX 'IFRBSTA'
 66   NTEST = 0
      DO K=2, IND
         IF (IFRBSTA(K) .LT. IFRBSTA(K-1)) THEN
            NTEST = 1
             
             IBAK = IFRBSTA(K)
             IFRBSTA(K) = IFRBSTA(K-1)
             IFRBSTA(K-1) = IBAK
             
             IBAK = IFRBEND(K)
             IFRBEND(K) = IFRBEND(K-1)
             IFRBEND(K-1) = IBAK
             
             IBAK = INDRB(K)
             INDRB(K) = INDRB(K-1)
             INDRB(K-1) = IBAK
             
             IBAK = IFENUP(K)
             IFENUP(K) = IFENUP(K-1)
             IFENUP(K-1) = IBAK

             IBAK = IFELOW(K)
             IFELOW(K) = IFELOW(K-1)
             IFELOW(K-1) = IBAK
             
             XBAK = SIGMAINT(K)
             SIGMAINT(K) = SIGMAINT(K-1)
             SIGMAINT(K-1) = XBAK
             
          ENDIF
       ENDDO 

       IF (NTEST .EQ. 1) GOTO 66
C ***  END OF BUBBLESORT

C *** SAVE ARRAY-LENGHTS
      LASTKON  = KONT
      LASTFE   = IND

C***  Dimension check
      IF (LASTINDAUTO+LASTFE .GT. MAXIND) THEN
         WRITE (0,'(A)') '*** Dimension MAXIND insufficient !' 
         WRITE (0,'(A,I5)') '*** Available: MAXIND = ', MAXIND
         WRITE (0,'(A,I5)') 
     >            '*** Required : MAXIND = ', LASTINDAUTO+LASTFE
         STOP            '*** ERROR STOP IN FEDAT *********'
      ENDIF
      
C***  Append the Superline indices to the line-transition vectors 
      DO K=1, LASTFE
         INDNUP(LASTINDAUTO+K) = IFENUP(K)
         INDLOW(LASTINDAUTO+K) = IFELOW(K)
C***     set valid keyword for bound-bound collision rates
         KEYCBB(LASTINDAUTO+K) = 'NULL'
      ENDDO

      CALL CLOSMS(21,IERR)
            
      RETURN
      END
      SUBROUTINE FGRID (NFDIM, NF, XLAMBDA, FWEIGHT, KEY, NOM, SYMBOL, 
     $                  N, NCHARG, ELEVEL, EION, EINST, NDIM,
     $                  EDGEK, KODAT, MAXATOM, MAXION,
     $                  INDNUP, INDLOW, LASTIND, KONTNUP, KONTLOW, 
     >                  LASTKON, OLDFGRID, NF2, XLAMBDA2, 
     >                  VDOP, XLAMBLUE)
 
C***********************************************************************
C***  GENERATION OF THE FREQUENCY GRID AND INTEGRATION WEIGHTS
C***  INCLUDING PREDEFINED FREQUENCY POINTS (FROM TAPE6)
C***  AND THE CONTINUUM FREQUENCY POINTS (NO LINE FREQUENCY POINTS)
C***  XLAMBDA = CORRESPONDING WAVELENGTH POINTS IN ANGSTROEMS
C***********************************************************************
 
C***  ATTENTION: The following common block is not complete 
C                since only VFINAL is used in this routine
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE

      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, FWEIGHT
      REAL, DIMENSION(N) :: ELEVEL, EION
      REAL, DIMENSION(MAXATOM,MAXION) :: EDGEK
      INTEGER, DIMENSION(N) :: NCHARG, NOM
      INTEGER, DIMENSION(LASTIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(LASTKON) :: KONTNUP, KONTLOW
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      LOGICAL :: OLDFGRID, BADD
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8) :: NAME, NEDGE, CKEY
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      CHARACTER(100) :: CFORMAT, MODOLD
      REAL, INTENT(IN) :: XLAMBLUE
      
      REAL :: WINGpos, WINGneg
      INTEGER :: KOMIT, KLASTK
 
C***  Constants
      REAL, PARAMETER :: STEBOL = 1.8046E-5     !STEFAN-BOLTZMANN CONSTANT / PI  (ERG/CM**2/S/STERAD/KELVIN**4)
      REAL, PARAMETER :: CLIGHT = 2.99792458E18 !SPEED OF LIGHT IN ANGSTROEM / SECOND

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
C***  BEGIN OF THE LISTING
      WRITE (hOUT,FMT='(A,//,20X,A,/)') '1',
     >        '2. FREQUENCY POINTS AND INTEGRATION WEIGHTS'
 
C***  TAKE THE FREQUENCY GRID (HAND-MADE WAVELENGTH POINTS) 
C***      EITHER FROM TAPE6 = FGRID, OR FROM OLD MODEL FILE
      CALL       DECFREQ (XLAMBDA, NF, NFDIM, TREF, OLDFGRID, KEY, 
     >                    MODOLD, XLAMBLUE)

C***  The FREQUENCY GRID IS TESTED FOR TEFF AND FOR 2.0 * TEFF;
C***  THE SECOND VALUE IS THE EXPECTED TEMPERATURE AT TAUROSS = 20
      TREF2 = 2.0 * TREF
 
C***  SET KEYWORD ARRAY TO BLANKS
      DO K=1,NF
        KEY(K) = '        '
      ENDDO
 
C***  ADDITION OF THE REDMOST LINE FREQUENCY POINT (MINUS EPSILON!)
 
C***  FIND THE REDMOST LINE CENTER FREQUENCY POINT
      IF (XLAMBDA(NF) .GT. XLAMBDA(1)) THEN
          REDMOST=XLAMBDA(NF)
      ELSE
          REDMOST=XLAMBDA(1)
      ENDIF
      DO IND=1,LASTIND
        J=INDNUP(IND)
        I=INDLOW(IND)
        WLENG=1.E8/(ELEVEL(J)-ELEVEL(I))
        IF (WLENG .GT. REDMOST) THEN
            REDMOST=WLENG
            INDRED=IND
            IRED=I
            JRED=J
        ENDIF
      ENDDO
      IND=INDRED
      WLENG = REDMOST * (1. + VFINAL / 100000.)
C***  Minimum redmost wavelength is 100000. A
      IF (WLENG .LT. 100000.) THEN
        WLENG = 100000.
        WRITE (UNIT=NAME, FMT='(A8)') 'REDMIN  '
      ELSE
        IF (EINST(IRED,JRED) .EQ. -2.) THEN
            WRITE (UNIT=NAME, FMT='(A4,I4)') 'RUD ', IND
        ELSE
            WRITE (UNIT=NAME, FMT='(A4,I4)') 'LINE', IND
        ENDIF
      ENDIF

      IF (NF .GT. 1) THEN 
         CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WLENG)
         IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NAME)
         IF (K.LT.0) KEY(-K)=NAME
      ELSE
C***  If FGRID gives only one lambda point, 
C***    the redmost lambda becomes the second point  
         NF = 2
         XLAMBDA(2) = WLENG
         KEY(2) = NAME
      ENDIF 
 
C***  INSERTION OF CONTINUUM EDGES
C***  One point each is inserted before and after the edge
C***  at minus and plus 0.3*VDOP  
      EDGEWIDTH = 0.3 * VDOP * 1.E13 / CLIGHT   

      DO KON=1,LASTKON
        J=KONTNUP(KON)
        I=KONTLOW(KON)
        WLENG=1.E8/(ELEVEL(J)+EION(I)-ELEVEL(I))
        WPLUS = (1. + EDGEWIDTH) * WLENG
        CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WPLUS)
        IF (SYMBOL(NOM(I)) /= 'G ') THEN
          WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >       'EDGE+', SYMBOL(NOM(I)),NCHARG(I)+1
        ELSE 
          WRITE (UNIT=NEDGE, FMT='(A5,A1,I2)') 
     >       'EDGE+', SYMBOL(NOM(I))(1:1),NCHARG(I)+1
        ENDIF
        IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
        WMINUS = (1. - EDGEWIDTH) * WLENG
        CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WMINUS)
        IF (SYMBOL(NOM(I)) /= 'G ') THEN
          WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >      'EDGE-', SYMBOL(NOM(I)),NCHARG(I)+1
        ELSE 
          WRITE (UNIT=NEDGE, FMT='(A5,A1,I2)') 
     >      'EDGE-', SYMBOL(NOM(I)),NCHARG(I)+1
        ENDIF
        IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
      ENDDO
 
C***  INSERTION OF K-SHELL-IONISATION EDGES

      DO NZ=1,MAXATOM
         KOJ=KODAT(NZ)
         IF  (KOJ == 0) CYCLE
         DO ISTATE=1, NZ-2
            IF (EDGEK(KOJ,ISTATE) .EQ. .0) CYCLE
            WLENG = 1.E8 / EDGEK(KOJ,ISTATE)
            WPLUS = (1. + EDGEWIDTH) * WLENG
            CALL SEQUIN (NFDIM, NF, XLAMBDA, K, WPLUS)
            WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >          'K-ED+', SYMBOL(KOJ), ISTATE
            IF (K .GT. 0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
            WMINUS = (1. - EDGEWIDTH) * WLENG
            CALL SEQUIN (NFDIM,NF,XLAMBDA,K,WMINUS)
            WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >          'K-ED-', SYMBOL(KOJ), ISTATE
            IF (K .GT. 0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
         ENDDO
      ENDDO

      BTOT = STEBOL * TREF**4.
      BTOT2 = STEBOL * TREF2**4.

      NF2 = NF - 2
      DO K=1, NF2
        XLAMBDA2(K) = XLAMBDA(K+1)
      ENDDO

C***  Entry-Point for additional frequencies
      IADD = 0
  100 CONTINUE

C***  Calculate RELMAX, and the Index of its Maximum (KMAX) for both
C***  Temperatures
      KMAX1 = 0
      KMAX2 = 0
      RELMAX1 = 0.
      RELMAX2 = 0.
      DO K=1,NF-1
        XLAM1 = XLAMBDA(K)
        XLAM2 = XLAMBDA(K+1)
        W = (1./XLAMBDA(K) - 1./XLAMBDA(K+1)) * CLIGHT
        RELCONT1 = 0.5 * W * 
     >             (BNUE(XLAM1,TREF) + BNUE(XLAM2,TREF)) / BTOT
        IF (RELCONT1 .GT. RELMAX1) THEN
          KMAX1 = K
          RELMAX1 = RELCONT1
        ENDIF
        RELCONT2 = 0.5 * W * 
     >             (BNUE(XLAM1,TREF2) + BNUE(XLAM2,TREF2)) / BTOT2
        IF (RELCONT2 .GT. RELMAX2) THEN
          KMAX2 = K
          RELMAX2 = RELCONT2
        ENDIF
      ENDDO

C***  Insert additional frequency points
C***  First: to ensure small relative contributions
      BADD = .FALSE.
      XNF = FLOAT(NF)
      IF (RELMAX1 .GT. RELMAX2) THEN
        IMAX = KMAX1
        RELMAX = RELMAX1
      ELSE
        IMAX = KMAX2
        RELMAX = RELMAX2
      ENDIF
      IF (NF .GT. NFDIM-1) THEN
        WRITE (0,*) 
     >        'WARNING : NFDIM insufficient for automatic Continuum', 
     >        '  Frequency-Point adjustment'
        WRITE (0,'(A, I5)') 'NFDIM=', NFDIM
        WRITE (*,*) 
     >        'WARNING : NFDIM insufficient for automatic Continuum', 
     >        '  Frequency-Point adjustment'
        WRITE (*,'(A, I5)') 'NFDIM=', NFDIM
      ELSE
C***  Criterium set to 2.5 percent
        IF (RELMAX .GT. 0.025) THEN
C***    Red: Put New Point in the Intervall IMAX and IMAX+1
          IF (IMAX .LT. NF) THEN
            WNEWR = (XLAMBDA(IMAX) + XLAMBDA(IMAX+1)) / 2.
            WRITE (UNIT=CKEY, FMT='(A8)') KEY(IMAX+1)
            CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WNEWR)
            IADD = IADD + 1
            NAME = 'ADD     '
            IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)
            BADD = .TRUE.
          ENDIF
        ENDIF
      ENDIF
      IF (BADD) GOTO 100

C***  Second: To ensure small Delta-Lambda / Lambda steps
      IADD_DL = 0
C***  Define Treshhold to distiguish between intervalls
C***    Treshhold is the first Continuum frequency point redwards starting value
ccc      TRESH_START  = 2000.
      TRESH_START2 = 227.83774
      TRESH_START3 = 504.259
      THRESH_IRMID = 30000.  ! 3 micron
      CRIT1 = 0.1
      CRIT2 = 1.0
      IF (TRESH_START2 .LT. XLAMBDA(1) .OR. 
     >    TRESH_START2 .GT. XLAMBDA(NF-1)) THEN
        WRITE (0,'(2A,3(G12.5,1X))') 
     >              'Troubles in determining TRESHHOLD: ', 
     >              'TRESH_START2, XLAMBDA(1), XLAMBDA(NF-1)=', 
     >               TRESH_START2, XLAMBDA(1), XLAMBDA(NF-1)
      ENDIF
      IF (TRESH_START3 .LT. XLAMBDA(1) .OR. 
     >    TRESH_START3 .GT. XLAMBDA(NF-1)) THEN
        WRITE (0,'(2A,3(G12.5,1X))') 
     >              'Troubles in determining TRESHHOLD: ', 
     >              'TRESH_START3, XLAMBDA(1), XLAMBDA(NF-1)=', 
     >               TRESH_START3, XLAMBDA(1), XLAMBDA(NF-1)
      ENDIF

  110 CONTINUE
      BADD = .FALSE.
C***  Find first Step larger than criterion
      DO K=1, NF-1
        IF (XLAMBDA(K) >= THRESH_IRMID) THEN
          CRIT = CRIT2
        ELSEIF (XLAMBDA(K) .GE. TRESH_START3) THEN
          CRIT = 0.1 * CRIT1
        ELSE
          IF (XLAMBDA(K) .LE. TRESH_START2 .AND. 
     >        (XLAMBDA(K)-20.) .GT. 0.) THEN
            F2 = ( (XLAMBDA(K)-20.)**8 / 6.4E16 ) + 1
C!!!        write (0,*) '!!!', K, XLAMBDA(K), F2
          ELSE
            F2 = 1.
          ENDIF
          IF (XLAMBDA(K) .LE. TRESH_START3 .AND. 
     >        (XLAMBDA(K)-300.) .GT. 0.) THEN
            F3 = ( (XLAMBDA(K)-300.)**8 / 6.4E16 ) + 10
C!!!        write (0,*) '!!!', K, XLAMBDA(K), F3
          ELSE
            F3 = 1.
          ENDIF
          CRIT = CRIT1 / (F2 * F3)
        ENDIF
        XL_MID = 0.5 * (XLAMBDA(K) + XLAMBDA(K+1))
        DLL = (XLAMBDA(K+1) - XLAMBDA(K)) / XL_MID
        IF (DLL .GT. CRIT) THEN
          KINS = K
          BADD = .TRUE.
          EXIT
        ENDIF
      ENDDO

C***  Insert new Point in Interval KINS, KINS+1
      IF (BADD) THEN
        IF (NF .GT. NFDIM-1) THEN
          WRITE (0,*) 
     >          'WARNING : NFDIM insufficient for automatic Continuum', 
     >          '  Frequency-Point adjustment'
          WRITE (0,'(A, I5)') 'NFDIM=', NFDIM
          WRITE (*,*) 
     >          'WARNING : NFDIM insufficient for automatic Continuum', 
     >          '  Frequency-Point adjustment'
          WRITE (*,'(A, I5)') 'NFDIM=', NFDIM
        ELSE
          IF (IMAX .LT. NF) THEN
            WRITE (UNIT=CKEY, FMT='(A8)') KEY(IMAX+1)
            CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XL_MID)
            IADD_DL = IADD_DL + 1
            NAME = 'ADD_DL  '
            IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)
          ENDIF
        ENDIF
        GOTO 110
      ENDIF

C***  New 17-Aug-2015: 
C***  insert exact wavelengths of Smith ubv monochromatic magnitudes
C***  at 3650., 4270., 5160. Ang
      XLAMSMITH = 3650.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH u '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMSMITH = 4270.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH b '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMSMITH = 5160.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH v '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMMASSEY = 6000.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMMASSEY)
      NAME = 'MASSEY r'
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

C***  All additional Frequency Points are now inserted

C***  Now OMIT frequency points which are too close!
C***  Criterion is 1.0 * VDOP
      KLAST=2
  200 DO K=KLAST, NF
        IF (XLAMBDA(K)/XLAMBDA(K-1) .LT. (1.+EDGEWIDTH) ) THEN
          WRITE (0,'(A,F10.2,2X,A)') 
     >    'FGRID: Cont. frequency point omitted: ', XLAMBDA(K), KEY(K) 
          NF = NF - 1
          KLAST=K
          DO KK=K, NF
            XLAMBDA(KK) = XLAMBDA(KK+1)
            KEY    (KK) = KEY    (KK+1)
          ENDDO
          GOTO 200
        ENDIF
      ENDDO

C***  FREQUENCY INTEGRATION WEIGHTS ACCORDING TO THE TRAPEZOIDAL RULE **********
      NFM=NF-1
      DO K=2,NFM
        FWEIGHT(K)=.5*(1./XLAMBDA(K-1) - 1./XLAMBDA(K+1))*CLIGHT
      ENDDO
      FWEIGHT(1)=.5*(1./XLAMBDA(1) - 1./XLAMBDA(2))*CLIGHT
      FWEIGHT(NF)=.5*(1./XLAMBDA(NFM) - 1./XLAMBDA(NF))*CLIGHT
 
C***  RENORMALIZATION TO RETAIN THE EXACT INTEGRAL SUM OF PLANCKS FUNCTION
C***   AT REFERENCE TEMPERATURE TREF
      SUM=.0
      DO K=1,NF
        SUM=SUM+FWEIGHT(K)*BNUE(XLAMBDA(K),TREF)
      ENDDO
      RENORM=BTOT/SUM
      DO K=1,NF
        FWEIGHT(K)=FWEIGHT(K)*RENORM
      ENDDO
 
C***  OUTPUT
C***  CONTINUATION OF THE LISTING
      IF (OLDFGRID) THEN
        WRITE (hOUT, FMT='(A,A)') ' FGRID TAKEN FROM OLD MODEL: ',MODOLD
      ENDIF

      WRITE (hOUT,*) 
      WRITE (hOUT,'(A)') 
     >      'NOTE : Continuum frequencies inserted'
      WRITE (hOUT,'(A,I4)')
     >      '       Relative Contribution Criterion: ', IADD
      WRITE (hOUT,'(A,I4)')
     >      '       Delta-Lambda / Lambda Criterion: ', IADD_DL

      PRINT 18
   18 FORMAT(//,1X,
     > '  NR   LAMBDA/ANGSTROEM    KEY     WAVENUMBER(KAYSER)  ',
     > ' FREQUENCY(HERTZ)        WEIGHT    REL.INTEGRAL CONTRIBUTION',/)

      CFORMAT = '(1X,I5,F15.2,2X,A8,1X,F15.2,2E22.6,F8.3,3X,F8.3)'
      DO K=1,NF
        XLAM=XLAMBDA(K)
        RELCO=NF*FWEIGHT(K)*BNUE(XLAM,TREF)/BTOT
        RELCO2=NF*FWEIGHT(K)*BNUE(XLAM,TREF2)/BTOT2
        WRITE (hOUT, FMT=CFORMAT)  K, XLAM, KEY(K), 1.E8/XLAM,
     >             CLIGHT/XLAM, FWEIGHT(K), RELCO, RELCO2       
      ENDDO

      WRITE (hOUT, FMT='(//,A,F10.6,A,F8.0,A)')
     >    ' RENORMALIZATION FACTOR :', RENORM,
     >    '     REFERENCE TEMPERATURE :', TREF, ' K'

      IF (NF > 9999) THEN
        WRITE (0,*)
     >     '*** MORE THAN 9999 COARSE FREQUENCY POINTS ENCOUNTERED ***'
        WRITE (0,*) 'This is not compatible with the encoding of the'
        WRITE (0,'(A)') ' frequency index in the MODEL file variables'
     >     // ' XJCnnnn, WJCnnnn and EDDInnnn.'
        STOP 'FATAL ERROR IN FGRID'
      ENDIF
      
      RETURN
      END
      SUBROUTINE FINDCHARGE (THISNAME, NZ)
C*****************************************************************
C***  This subroutine searches for THISNAME in the list of 
C***  chemical elements and returns NZ = core charge 
C***  (Ordnungszahl im Periodischen System)
C***  If THISNAME is not found, NZ=0 is returned.
C***  Note that GENERIC is in the position of IRON 
C***  Called from: DATOM, DECSTAR
C*****************************************************************

      CHARACTER THISNAME*(*)
      PARAMETER (MAXELEM = 26)
      CHARACTER*10 ELEMNAME(MAXELEM)
      DATA ELEMNAME /'HYDROGEN  ', 'HELIUM    ', 'LITHIUM   ', 
     2               'BERYLLIUM ', 'BORON     ', 'CARBON    ', 
     3               'NITROGEN  ', 'OXYGEN    ', 'FLUORINE  ',
     4               'NEON      ', 'SODIUM    ', 'MAGNESIUM ', 
     5               'ALUMINIUM ', 'SILICON   ', 'PHOSPHORUS',
     6               'SULFUR    ', 'CHLORINE  ', 'ARGON     ', 
     7               'POTASSIUM ', 'CALCIUM   ', 'SCANDIUM  ', 
     8               'TITANIUM  ', 'VANADIUM  ', 'CHROMIUM  ', 
     9               'MANGANESE ', 'GENERIC   '/

      NZ = 0
      DO K = 1, MAXELEM
         IF (THISNAME .EQ. ELEMNAME(K)) THEN
           NZ = K
           EXIT
         ENDIF
      ENDDO

      END
      SUBROUTINE GAUNTFF (GIII,NCHARGE,XLAM,TEMP)
C***********************************************************************
C***  THE MODULE COMPUTES THE G-III FACTOR FOR THERMAL BREMSSTRAHLUNG
C***  GIII = GAUNT FACTOR
C***  NCHARGE = REST CHARGE OF ION ( =1 FOR H+ AND HE+, =2 FOR HE++ )
C***  XLAM    = WAVELEGTH IN ANGSTROEM
C***  TEMP    = TEMPERATURE IN KELVIN
C***  THE GAUNT FACTORS ARE COMPILED FROM THE PUPLICATIONS OF BERGER (1956)
C***  APJ 124,P550  AND KARZAS AND LATTER (1961) APJ SUPPL 6,P167 (FIG 3,4,5)
C***  TO GET THE HERE TABULATED VALUES BICUBIC SPLINE INTERPOLATION WAS USED
C***  --------  TABULATED RANGE FOR NCHARGE = 1  -----------------------
C***  TEMPERATURE (KELVIN)             :  1577  ...  157 700
C***  WAVELENGTH (ANGSTROEM)           :   120  ...  1.2 E6
C***********************************************************************
 
      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP

      PARAMETER ( MAXION = 27 )
 
      REAL A(21,21)
      DIMENSION ZLOG(MAXION)
 
      DATA      ZLOG / 0.     , 0.30103, 0.47712, 0.60206, 0.69897,
     >                 0.77815, 0.84510, 0.90309, 0.95424, 1.0    ,
     >                 1.04139, 1.07918, 1.11394, 1.14613, 1.17609, 
     >                 1.20412, 1.23045, 1.25527, 1.27875, 1.30103,
     >                 1.32222, 1.34232, 1.36173, 1.38021, 1.39794,
     >                 1.41497, 1.43136 /
      DATA ((A(I,J),I=1,21),J=1,6) /
     $ 1.331, 1.274, 1.232, 1.200, 1.177, 1.158, 1.143, 1.130, 1.120,
     * 1.115, 1.112, 1.108, 1.103, 1.100, 1.100, 1.101, 1.098, 1.090,
     $ 1.080, 1.069, 1.057,
     $ 1.412, 1.347, 1.299, 1.258, 1.214, 1.173, 1.145, 1.132, 1.126,
     * 1.122, 1.118,1.114, 1.110, 1.106, 1.102, 1.099, 1.098, 1.096,
     * 1.090, 1.078, 1.060,
     $ 1.503, 1.418, 1.350, 1.289, 1.230, 1.178, 1.146, 1.133, 1.129,
     * 1.125, 1.121, 1.117, 1.113, 1.109, 1.104, 1.100, 1.100, 1.100,
     $1.098, 1.084, 1.062,
     $ 1.601, 1.489, 1.394, 1.311, 1.239, 1.182, 1.148, 1.135, 1.131,
     $ 1.128, 1.123, 1.119, 1.115, 1.111, 1.106, 1.102, 1.102, 1.105,
     $ 1.102, 1.087, 1.064,
     $ 1.704, 1.566, 1.443, 1.338, 1.254, 1.193, 1.157, 1.142, 1.136,
     $ 1.131, 1.126, 1.121, 1.118, 1.114, 1.109, 1.105, 1.105, 1.104,
     $ 1.103, 1.087, 1.064,
     $ 1.807, 1.652, 1.509, 1.386, 1.288, 1.218, 1.176, 1.156, 1.147,
     * 1.139, 1.131, 1.126, 1.122, 1.119, 1.114, 1.110, 1.109, 1.107,
     $ 1.101, 1.085, 1.063 /
      DATA ((A(I,J),I=1,21),J=7,12) /
     $ 1.911, 1.750, 1.600, 1.466, 1.352, 1.264, 1.207, 1.179, 1.165,
     $ 1.153, 1.142, 1.129, 1.126, 1.126, 1.121, 1.116, 1.111, 1.106,
     $ 1.096, 1.081, 1.060,
     $ 2.016, 1.849, 1.702, 1.562, 1.431, 1.322, 1.248, 1.209, 1.188,
     $ 1.172, 1.158, 1.142, 1.136, 1.135, 1.130, 1.122, 1.114, 1.104,
     $ 1.091, 1.076, 1.059,
     $ 2.121, 1.947, 1.790, 1.644, 1.500, 1.375, 1.288, 1.241, 1.216,
     $ 1.196, 1.179, 1.158, 1.149, 1.145, 1.137, 1.128, 1.116, 1.103,
     $ 1.090, 1.076, 1.061,
     $ 2.226, 2.046, 1.876, 1.713, 1.561, 1.429, 1.335, 1.281, 1.250,
     $ 1.225, 1.204, 1.177, 1.164, 1.157, 1.144, 1.132, 1.119, 1.105,
     $ 1.091, 1.078, 1.064,
     $ 2.331, 2.146, 1.970,1.802, 1.644,  1.507, 1.407, 1.343, 1.301,
     $ 1.265, 1.234, 1.201, 1.182, 1.172, 1.145, 1.138, 1.122, 1.105,
     $ 1.090, 1.076, 1.061,
     $ 2.500, 2.307, 2.123, 1.947, 1.766, 1.623, 1.513, 1.435, 1.374,
     $ 1.318, 1.268, 1.229, 1.204, 1.189, 1.168, 1.147, 1.126, 1.105,
     $ 1.085, 1.067, 1.046 /
      DATA ((A(I,J),I=1,21),J=13,18) /
     $ 2.659, 2.460, 2.268,2.084, 1.907, 1.755, 1.633, 1.539, 1.458,
     $ 1.379, 1.309, 1.262, 1.229, 1.208, 1.185, 1.158, 1.132, 1.106,
     $ 1.081, 1.055, 1.026,
     $ 2.809, 2.604, 2.407, 2.217, 2.035, 1.873, 1.740, 1.634, 1.538,
     $ 1.442, 1.356, 1.300, 1.259, 1.230, 1.202, 1.172, 1.141, 1.112,
     $ 1.082, 1.049, 1.012,
     $ 2.953, 2.743, 2.541, 2.345, 2.156, 1.973, 1.830, 1.714, 1.610,
     $ 1.506, 1.410, 1.344, 1.294, 1.255, 1.221, 1.187, 1.155, 1.124,
     $ 1.090, 1.050, 1.005,
     $ 3.091, 2.878, 2.670, 2.469, 2.276, 2.090, 1.923, 1.797, 1.683,
     $ 1.573, 1.471, 1.394, 1.333, 1.285, 1.242, 1.204, 1.170, 1.136,
     $ 1.099, 1.054, 1.001,
     $ 3.261, 3.042, 2.829, 2.622, 2.421, 2.227, 2.039, 1.895, 1.766,
     $ 1.646, 1.537, 1.451, 1.378, 1.319, 1.268, 1.224, 1.184, 1.146,
     $ 1.104, 1.055, 0.997,
     $ 3.432, 3.208, 2.989, 2.776, 2.569, 2.368, 2.175, 2.008, 1.858,
     $ 1.725, 1.610, 1.514, 1.429, 1.358, 1.299, 1.246, 1.199, 1.155,
     $ 1.108, 1.054, 0.993 /
      DATA ((A(I,J),I=1,21),J=19,21) /
     $ 3.593, 3.365, 3.142, 2.924, 2.711, 2.505, 2.305, 2.127, 1.956,
     $ 1.811, 1.689, 1.583, 1.485, 1.403, 1.335, 1.271, 1.216, 1.166,
     $ 1.114, 1.057, 0.994,
     $ 3.747, 3.515, 3.289, 3.067, 2.850, 2.638, 2.432, 2.233, 2.055,
     $ 1.903, 1.776, 1.660, 1.549, 1.455, 1.374, 1.301, 1.239, 1.183,
     $ 1.127, 1.066, 1.001,
     $ 3.895, 3.661, 3.431, 3.205, 2.984, 2.768, 2.558, 2.354, 2.152,
     $ 2.000, 1.873, 1.743, 1.620, 1.513, 1.415, 1.336, 1.269, 1.209,
     $ 1.151, 1.083, 1.018 /
 
      DATA NHELP / 0 /

C***  LOGARITHM OF TEMP AND XLAM
      TEMPLOG=ALOG10(TEMP)
      XLAMLOG=ALOG10(XLAM)
      GOTO 1
 
C***  ENTRY FOR CALLING WITH LOGARITHMIC ARGUMENTS
      ENTRY GFFLOG (GIII,NCHARGE,XLOG,TLOG)
      XLAMLOG=XLOG
      TEMPLOG=TLOG
    1 CONTINUE
 
      IF (NCHARGE .LE. 0 .OR. NCHARGE .GT. MAXION) THEN
         CALL REMARK ('NCHARGE OUTSIDE VALID RANGE')
         STOP 'ERROR in Subr. GAUNTFF'
         ENDIF
 
C***  AT THE FIRST CALL WITHIN A MAIN PROGRAMM, THE COUNTERS ARE SET ZERO.
C***  THESE COUNTERS INDICATE HOW OFTEN THIS ROUTINE HAS BEEN CALLED
C***  WITH PARAMETERS OUTSIDE OF THE TABULATED RANGE.
      IF (NHELP .EQ. 0) THEN
         NHELP=1
         NTUP=0
         NTLOW=0
         NFUP=0
         NFLOW=0
      ENDIF
 
C***  SCALE TEMP AND XLAM WITH THE REST IONIC CHARGE
      XZLOG=XLAMLOG+2.*ZLOG(NCHARGE)
      TZLOG=TEMPLOG-2.*ZLOG(NCHARGE)
 
C***  CALCULATE TABULAR INDICES (BROKEN NUMBERS)
      AT=10.*TZLOG    -30.98364
      AF=31.3944035-5.*XZLOG
 
      MT=AT
C***  LOWER TEMPERATURE BOUNDARY
      IF (MT.LT.1) THEN
         MT=1
         AT=1.
         NTLOW=NTLOW+1
         ENDIF
C***  UPPER TEMPERATURE BOUNDARY
      IF (MT.GT.20) THEN
         MT=20
         AT=20.
         NTUP=NTUP+1
         ENDIF
 
      MF=AF
C***  LOWER FREQUENCY BOUNDARY
      IF (MF.LT.1) THEN
         MF=1
         AF=1.
         NFLOW=NFLOW+1
         ENDIF
C***  UPPER FREQUENCY BOUNDARY
      IF (MF.GT.20) THEN
         GIII=1.
C***     THIS APPROXIMATION BECOMES WORSE IN THE X-RAY REGION.
C***     BEYOND NU=10**17 BETTER VALUES SHOULD BE CALCULATED
         IF (MF .GT. 24) NFUP=NFUP+1
         RETURN
         ENDIF
 
C***  INTERPOLATION
      XF=AF-MF
      GF1=A(MF,MT)+XF*(A(MF+1,MT)-A(MF,MT))
      GF2=A(MF,MT+1)+XF*(A(MF+1,MT+1)-A(MF,MT+1))
      GIII=GF1+(AT-MT)*(GF2-GF1)

      IF (GIII < 0.) THEN
        WRITE (0,*) "XLAM=", XLAM, " TEMP=", TEMP
        WRITE (0,*) "TZLOG=", TZLOG, " ZLOG(NCHARGE)=", ZLOG(NCHARGE)
        WRITE (0,*) "MF=", MT, " MT=", MT
        WRITE (0,*) "GF1=", GF1, " GF2=", GF2
        WRITE (0,*) "AT=", AT, " AT-MT=", AT-MT
        WRITE (0,*) "negativer Gauntfaktor: ", GIII
        STOP "------- FATAL ERROR: SOMETHING IS TOTALLY WRONG -------"
      ENDIF
 
c      giii = 1.

      RETURN
      END
      SUBROUTINE GEOMESH (RADIUS, INCRIT,P,Z,ND,NDDIM,NP,NPDIM,RMAX, 
     >                    RadiusGridParameters, bNoRGrid, NC)
C***********************************************************************
C***  THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN R, P AND Z
C****   called only by WRSTART (wrh, goetz) and STEAL->HYDROSOLVE (wrh)
C***********************************************************************
  
      IMPLICIT NONE

      LOGICAL bNoRGrid
      INTEGER :: ND, NDDIM, NP, NPDIM, NC
      REAL, DIMENSION(NDDIM) :: R
      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT      
      REAL, DIMENSION(ND) :: RADIUS, P
      REAL, DIMENSION(2) :: Z
      REAL :: RMAX, RR, PJ, PJPJ
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters
      
      INTEGER I, J, L, JMAX
      
      IF (.NOT. bNoRGrid) THEN
        CALL RGRID (NDDIM,ND,RADIUS,INCRIT,RMAX, 
     >              RadiusGridParameters)
      ENDIF

      CALL PGRID (NPDIM,NP,ND,RADIUS,P, NC)
      DO L=1,ND
        RR=RADIUS(L)*RADIUS(L)
        JMAX=NP+1-L
        DO J=1,JMAX
          PJ=P(J)
          PJPJ=PJ*PJ
          I=(J-1)*ND+L
          Z(I)=SQRT(RR-PJPJ)
        ENDDO
      ENDDO

      RETURN
      END
      SUBROUTINE GRADIFF  (ND,VELO,GRADI,RADIUS)
C***********************************************************************
C***  FOR THE VELOCITY FIELD GIVEN BY VECTOR VELO(L), THE GRADIENTS GRADI(L)
C***  ARE COMPUTED BY LINEAR INTERPOLATION BETWEEN THE NEIGHBORING POINTS
C***********************************************************************
      DIMENSION VELO(ND),GRADI(ND),RADIUS(ND)
      NDM=ND-1
      DO 1 L=2,NDM
      GRADI(L)=(VELO(L+1)-VELO(L-1))/(RADIUS(L+1)-RADIUS(L-1))
      IF (L.EQ.2) GRADI(1)=GRADI(2)
      IF (L.EQ.NDM) GRADI(ND)=GRADI(NDM)
    1 CONTINUE

      RETURN
      END
      SUBROUTINE GREY (ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR,
     $            ALPHA,SEXPO,
     $            ADDCON1, ADDCON2, ADDCON3, 
     >            IGAUNT,POPNUM,TAUROSS,R23,TEXIST,NDIM,N,
     $            LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,KODAT,
     $            ABXYZ,NOM,NFIRST,NLAST,NATOM,EXPFAC,SIGMAKI,NFEDGE,
     $            OPAC,ETAC,SIGMAFF,MAXION,MAXATOM,SIGMATHK,EDGEK,
     $            SEXPOK,KONTNUP,KONTLOW,LASTKON,XDATA, DENSCON, 
     >            FILLFAC)

C***********************************************************************
C***  COMPUTATION OF THE TEMPERATURE STRUCTURE
C***********************************************************************
C***  IF (TEXIST=.TRUE.), THE TEMPERATURE STRUCTURE IS NOT CALCULATED BUT
C***  IS ASSUMED TO BE EXISTING.
C***  ELSE IF (SPHERIC=.TRUE.) THE TEMPERATURE STRUCTURE IS CALCULATED
C***  AS FOR A SPHERICAL, GREY ATMOSPHERE. THE EDDINGTON FACTOR F IS
C***  APPROXIMATED BY ASSUMING A GEOMETRICALLY DILUTED RADIATION FIELD
C***  EMERGING FROM A SPHERE OF RADIUS R23 (I.E. THE RADIUS WHERE
C***  TAU-ROSSELAND= 2/3).
C***  THE SOLUTION OF THE 1. MOMENT EQUATION IS OBTAINED BY USING THE
C***  IMSL-ROUTINE DVERK. DTDR IS A USER-PROVIDED COEFFICIENT FUNCTION.
C***  ELSE
C***  THE TEMPERATURE STRUCTURE IS CALCULATED AS FOR A PLANE-PARALLEL,
C***  GREY ATMOSPHERE.
C***  ENDIF
C***
C***  IF (TMODIFY .NE. .0), THE CALCULATED TEMPERATURE STRUCTURE IS FINALLY
C***  MODIFIED BY A FACTOR R**TMODIFY .
C***  IF (TMIN .GT. 0), THE TEMPERATURE STRUCTURE IS MODIFIED NOT TO FALL
C***  BELOW TMIN.
C***  
C***  Clumping: LTE popnumbers is calculated with clump (electron) density, 
C***  Rosseland opacity is calculated with clump density, and then
C***  scaled down with the filling factor
C***********************************************************************
 
      COMMON /COMDTDR/ OPAMEAN,R23COM,R1COM,R13COM
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      INTEGER, PARAMETER :: ITER23MAX = 50

      DIMENSION XLAMBDA(NF),FWEIGHT(NF),EXPFAC(NF)
      DIMENSION SIGMAKI(NF,LASTKON)
      DIMENSION T(ND),RADIUS(ND),RNE(ND),ENTOT(ND),TAUROSS(ND)
      DIMENSION NCHARG(NDIM),ENLTE(NDIM),EION(NDIM),ELEVEL(NDIM)
      DIMENSION POPNUM(ND,N)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),NFEDGE(LASTKON)
      DIMENSION NOM(N)
      DIMENSION KODAT(NATOM),ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION COMVEC(24),WORKSP(1,9)
      LOGICAL TEXIST, SPHERIC
      CHARACTER*10 LEVEL(N)
      EXTERNAL DTDR

c*** tiefenabh. clumping nach goetz...
      DIMENSION DENSCON(ND), FILLFAC(ND)

C***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
C***  SIGMAKI(K,KON) IN CM**2
      CALL BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $              XLAMBDA,ALPHA,SEXPO,
     $              ADDCON1, ADDCON2, ADDCON3, 
     $              IGAUNT,
     $              KONTNUP,KONTLOW,LASTKON)

C***  PRE-CALCULATE FREQUENCY INDICES OF IONIZATION EDGES
      DO 30 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
      EDGELAM=1.E8/(ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW))
      NFEDGE(KON)=ISRCHFGT(NF,XLAMBDA,1,EDGELAM) - 1
   30 CONTINUE

C***  ABUNDANCES OF HELIUM AND HYDROGEN FOR THE MODIFIED START APPROXIMATION
            ABHE=0.
            ABH=0.
            IF (KODAT(1) .GT. 0) ABH =ABXYZ(KODAT(1))
            IF (KODAT(2) .GT. 0) ABHE=ABXYZ(KODAT(2))
C***  changed according to WR-Memo.dir/070327.txt --  wrh 26-Sep-2008 

C***  FIRST RGRID-POINT:  R=RMAX (L=1)
      TAUROSS(1)=.0
      IF (.NOT. TEXIST) THEN
         IF (SPHERIC) THEN
            T(1)=TEFF*0.7071068/SQRT(RADIUS(1))
            ELSE
            T(1)=TEFF*0.8112
            ENDIF
         ENDIF
 
C***  ITERATION LOOP  **************************************************
C***  MAIN PURPOSE IS THE ITERATION OF R23 IN THE SPHERICAL CASE
C***  START VALUE OF THE RADIUS R23  (ONLY USED IN THE SPHERICAL CASE)
      R23COM=1.
      R1COM=1.
      R13COM=1.
      ITER=0
  100 ITER=ITER+1
      DTMAX=.0
C***  INITIALIZATION OF VARIABLES USED IN THE IMSL-SUBROUTINE DVERK
      TOL=0.0001
      IND=1
      NW=1
 
C***  LOOP OVER ALL DEPTH POINTS  *******************************************
c      write (*,*) 'iter=',iter,'!!!!!!!!!!!!!!'
      DO 10 L=1,ND
      TL=T(L)
C***  Clump density!
      ENTOTL = ENTOT(L) * DENSCON(L) 
      IF (ITER .GT. 1) TOLD=T(L+1)
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY
      TP=TL
      IF (TP.LT.TMIN) TP=TMIN
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT POINT L
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY
      RNEL=RNE(L)
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
      IF (ITER .EQ. 1 .AND. TP .LT. 10000.) THEN
            T32=TP*SQRT(TP)
            RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP)*T32/2.07E-16
     $            *ABHE/ENTOTL)
            ENDIF

C***  Iteration of electron density
    3 ENE=RNEL*ENTOTL
      CALL       LTEPOP (N,ENLTE,TP,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
      RNEOLD=RNEL
      RNEL=.0
      DO 2 J=1,N
    2 RNEL=RNEL+NCHARG(J)*ENLTE(J)
      RNEDIF=RNEL-RNEOLD
      IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.001) GOTO 3

C***  STORE LTE POPNUMBERS TO BE WRITTEN AT THE MODEL FILE (START APPROXIMAT9ON@
      DO 5 J=1,N
    5 POPNUM(L,J)=ENLTE(J)
C*** ELECTRON DENSITY ALSO STORED
      RNE(L)=RNEL
 
      IF (L .EQ. ND) GOTO 10
 
      CALL OPAGREY (OPARL,ENLTE,TL,RNEL,ENTOTL,RSTAR,N,
     $              NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,NOM,
     $              EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,MAXION,
     $              SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $              KONTNUP,KONTLOW,LASTKON,RADIUS(L),XDATA)
      OPARL = OPARL * FILLFAC(L)
 
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY  AT POINT L+1
C***  IN THE FIRST ITERATION USING T(L)
C***  IN THE FOLLOWING ITERATIONS USING TOLD
      IF (ITER .EQ. 1) THEN
         TL1=TL
         ELSE
         TL1=TOLD
         ENDIF
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY
         TP1=TL1
         IF (TP1.LT.TMIN) TP1=TMIN
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY
      RNEL=RNE(L+1)
C***  Clump density!
      ENTOTL1 = ENTOT(L+1) * DENSCON(L+1)
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
      IF (ITER .EQ. 1 .AND. TP1 .LT. 1.E4) THEN
            T32=TP1*SQRT(TP1)
           RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP1)*T32/2.07E-16
     $            *ABHE/ENTOTL1)
            ENDIF
   13 ENE=RNEL*ENTOTL1
      CALL       LTEPOP (N,ENLTE,TP1,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
      RNEOLD=RNEL
      RNEL=.0
      DO 12 J=1,N
   12 RNEL=RNEL+NCHARG(J)*ENLTE(J)
      RNEDIF=RNEL-RNEOLD
      IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.001) GOTO 13
 
      CALL OPAGREY (OPARL1,ENLTE,TL1,RNEL,ENTOTL1,RSTAR,N,
     $              NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,NOM,
     $              EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,MAXION,
     $              SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $              KONTNUP,KONTLOW,LASTKON,RADIUS(L),XDATA)
      OPARL1 = OPARL1 * FILLFAC(L+1)
 
C***  ARITHMETIC MEAN OF OPARL AND OPARL1
      OPAMEAN=0.5*(OPARL+OPARL1)
      TAUROSS(L+1)=OPAMEAN*(RADIUS(L)-RADIUS(L+1))+TAUROSS(L)
 
      IF (.NOT. TEXIST) THEN
         IF (SPHERIC) THEN
            RL=RADIUS(L)
            RL1=RADIUS(L+1)
c            IF (OPSYS .EQ. 'CRAY') THEN
c              CALL DVERK(1,DTDR,RL,TL,RL1,TOL,IND,COMVEC,NW,WORKSP,IER)
c            ELSE
C!!!              CALL RUKU (RL, RL1, TL, 10000)
              CALL RUKU (RL, RL1, TL, 100)
c            ENDIF
c            stop 'stop in grey nach first step!!!!!!!!!'
c            write (*,*) 'grey: rl1, tl1=',rl1,tl1
            T(L+1)=TL
            ELSE
C***        HOPF FUNCTION, C.F. UNSOELD P. 138
            Q=0.6940-0.1167*EXP(-1.9720*TAUROSS(L+1))
            T(L+1)=TEFF*(0.75*(TAUROSS(L+1)+Q))**0.25
            ENDIF
         ENDIF
 
C***  MAXIMUM TEMPERATURE CORRECTION
      IF (ITER .GT. 1) DTMAX=AMAX1(DTMAX,ABS(TOLD-T(L+1)))
 
   10 CONTINUE
c      stop 'stop in grey'
C*****************************************************************************
 
C***  CALCULATE RADIUS R23 WHERE TAUROSS=2/3
      TAU23=0.666666666666
      IF (TAUROSS(ND) .LT. TAU23) THEN
         R23=1.
         ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND)
         ENDIF
      R23COM=R23
      TAU1=1.
      IF (TAUROSS(ND) .LT. TAU1  ) THEN
      R1COM=1.
         ELSE
         CALL LIPO (R1COM,TAU1 ,RADIUS,TAUROSS,ND)
         ENDIF
      TAU13=0.333333333333
      IF (TAUROSS(ND) .LT. TAU13 ) THEN
      R13COM=1.
         ELSE
         CALL LIPO (R13COM,TAU13,RADIUS,TAUROSS,ND)
         ENDIF
 
      IF (ITER .LE. 1 .OR. 
     >    (DTMAX .GT. 10. .AND. ITER .LE. ITER23MAX)) GOTO 100
C      IF (DTMAX .GT. 10.) GOTO 100
      IF (ITER .GT. ITER23MAX) THEN
        WRITE (0,*) 'Max Number of Iterations exceeded in Subr. GREY'
      ENDIF
 
      IF (.NOT. TEXIST) THEN
         IF (T(ND) .LE. TEFF ) THEN
            T(ND)=TEFF
            T(ND-1)=TEFF
            ENDIF
         ENDIF
 
C***  MODIFY THE TEMPERATURE STRATIFICATION BY A FACTOR RADIUS**TMODIFY
C***  AND/OR BY REQUIRING A MINIMUM TEMPERATURE TMIN
      DO 4 L=1,ND
      IF (TMODIFY .EQ. .0 .AND. T(L) .GT. TMIN ) GOTO 4
      IF (TMODIFY .NE. .0) T(L)=T(L)*RADIUS(L)**TMODIFY
      IF (T(L) .LT. TMIN) T(L)=TMIN
 
C***     CALCULATE LTE POP.NUMBERS WITH MODIFIED TEMPERATURE]
            TL=T(L)
            T32=TL*SQRT(TL)
            ENTOTL = ENTOT(L) * DENSCON(L)
            RNEL=RNE(L)
C***     MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
            IF (TL.LT.1.E4) THEN
            RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TL)*T32/2.07E-16
     $            *ABHE/ENTOTL)
            ENDIF
   23       ENE=RNEL*ENTOTL
      CALL       LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
            RNEOLD=RNEL
            RNEL=.0
            DO 24 J=1,N
   24       RNEL=RNEL+NCHARG(J)*ENLTE(J)
            RNEDIF=RNEL-RNEOLD
            IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.00001)
     $         GOTO 23
            DO 25 J=1,N
   25       POPNUM(L,J)=ENLTE(J)
            RNE(L)=RNEL
    4 CONTINUE
 
      RETURN
      END
      SUBROUTINE HYSTHDRUKU(RSTART, REND, VSTART, VEND, 
     >                      RADIUS, GEFFL, A2SUM, ND, RSTAR)
C***********************************************************************
C***  4th-order Runge-Kutta iteration for the equation of motion
C***  which has to be solved to obtain the velocity field in the
C***  quasi-hydrostatic regime
C***  
C***  called by VELTHIN  (for each depth point!)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
     
C***  Steps for Runge-Kutta integration     
      INTEGER, PARAMETER :: IRUKUSTEP = 50  
     
      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: RSTAR, RSTART, REND, VSTART
      REAL, INTENT(OUT) :: VEND
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, GEFFL, A2SUM
     
      REAL :: DR, DRHALF, DELTAR, RIN, RL, X1, X2, X3, X4,
     >        H, GEFFI, A2SUMI, DVDR, DA2DR, VIN, VL, R
     
      INTEGER :: I
     
      DELTAR = REND - RSTART 
      
      DR = DELTAR / IRUKUSTEP
      DRHALF = DR / 2.

      RL = RSTART
      VL = VSTART * 1.E5
                 
      DO I = 1, IRUKUSTEP

        RIN = RL      
        VIN = VL
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR = (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X1 = DR * DVDR

        RIN = RL + DRHALF
        VIN = VL + X1/2.
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR = (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X2 = DR * DVDR

        RIN = RL + DRHALF
        VIN = VL + X2/2.
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR =  (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X3 = DR * DVDR

        RIN = RL + DR
        VIN = VL + X3
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR =  (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X4 = DR * DVDR

        RL = RL + DR
        VL = VL + ((X1/2.) + X2 + X3 + (X4/2.)) / 3.
      ENDDO
      
      VEND = VL / 1.E5

      RETURN
      END
      SUBROUTINE HYSTRUKU(RSTART, REND, DELTAB, 
     >                    RADIUS, GEFFL, A2SUM, ND, H0, RSTAR)
C***********************************************************************
C***  simple Runge-Kutta iteration for the differential equation
C***  which has to be solved to obtain the solution of the transformed
C***  hydrostatic equation
C***  
C***  called by VELTHIN  (for each depth point!)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
     
C***  Steps for Runge-Kutta integration     
      INTEGER, PARAMETER :: IRUKUSTEP = 10  
     
      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: H0, RSTAR, RSTART, REND
      REAL, INTENT(OUT) :: DELTAB
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, GEFFL, A2SUM
     
      REAL :: DR, DRHALF, DELTAR, RIN, RL, X1, X2, X3, X4,
     >        H, GEFFI, A2SUMI, DBDR
     
      INTEGER :: I
     
      DELTAR = REND - RSTART 
      
      DR = DELTAR / IRUKUSTEP
      DRHALF = DR / 2.

      RL = RSTART
      DELTAB = 0.
         
      DO I = 1, IRUKUSTEP
        RIN = RL        
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X1 = DR * DBDR

        RIN = RL + DRHALF
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X2 = DR * DBDR

        RIN = RL + DR
        IF (RIN > RADIUS(1)) THEN
C***      In the rare case that this routine runs up to
C***      the outer boundary, we have to avoid that due to
C***      numerical reasons RIN is slightly larger than RMAX
C***      and thus SPLINPOX calls would fail.
          A2SUMI = A2SUM(1)
          GEFFI = GEFFL(1)
        ELSE
          CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
          CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        ENDIF
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X3 = DR * DBDR

        RL = RL + DR
        DELTAB = DELTAB + ((X1/2.) + 2. * X2 + (X3/2.)) / 3.
      ENDDO

      RETURN
      END

       FUNCTION IDX(TEXT)
C***  Returns the number of non-blank characters of string TEXT
       CHARACTER*(*) TEXT
       DO IDX=LEN(TEXT),1,-1
          IF(TEXT(IDX:IDX).NE.' ') RETURN
       END DO
       IDX=0
       RETURN
       END
      SUBROUTINE INITVEL (RMAX, TEFF, GEFFLOG, RSTAR, XMASS,
     >                    VTURB, bHScaleOnly, bHydroStat)
C*******************************************************************************
C***  INITIALIZATION OF THE VELOCITY-FIELD PARAMETERS
C***  CALLED FROM: WRSTART; STEAL - ENSURETAUMAX
C***  New parametrization of beta law, wrh  5-Mar-2001 
C***  2BETA-LAW, 6-Apr-2006
C*******************************************************************************

      IMPLICIT NONE

      REAL, INTENT(IN) :: RMAX, TEFF, GEFFLOG, RSTAR, XMASS
      LOGICAL, INTENT(INOUT) :: bHydroStat,    !returns .FALSE. if no hydrostatic domain is encountered
     >                          bHScaleOnly    !if true only HSCALE is calculated in this routine

      REAL :: VFINAL, VMINCAND, BETA, BETA2, BETA2FRACTION, VMIN,
     >        VPAR1, VPAR2, HSCALE, VPAR1_2, VPAR2_2, RCON, VCON, 
     >        VTURB

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >            BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      REAL, EXTERNAL :: DELTAGR
 
      REAL :: VMAX, Q, S, T, RCONOLD, RCONMIN, RCONMAX, RGRADMAX, DXX
      INTEGER :: ITER

      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: AMU   = 1.66E-24   !ATOMIC MASS UNIT (GRAMM)

C***  DETERMINATION OF THE SCALE HEIGHT (HYDROSTATIC EQ.)
C***  XMASS= MEAN MASS IN AMU (e.g.  2.0 FOR HE II, 1.33 FOR HE III )
C***  Turbulence pressure added 13-Mar-2014
      HSCALE = (BOLTZ*TEFF/(XMASS*AMU) + (VTURB*1.E5)**2)
     >                   /10.**GEFFLOG /RSTAR

      IF (bHScaleOnly) RETURN     
      
C***  In case of 2BETA-lwa, the second beta must be larger than 1
C***  in order to make sure that the contribution at the connection
C***  point vanishes 
      IF (BETA2FRACTION > 0) THEN
         IF (BETA2 .LT. 1.) THEN       
           WRITE (0, '(A)') '*** SECOND BETA MUST BE > 1' 
           STOP             '*** FATAL ERROR IN INITVEL'
         ENDIF
      ENDIF

      VMAX = VFINAL * (1.-BETA2FRACTION)
      Q = (VMIN/VMAX)**(1./BETA)
      
ccc   HSCALE=(BOLTZ*TEFF) / (XMASS*AMU*RSTAR*10.**GEFFLOG)
ccc   WRITE (0,*) " HSCALE: ", HSCALE
 

C***  ITERATIVE DETERMINATION OF THE CONNECTION POINT RCON
C***  BY REQUIRING A CONTINUOUS AND SMOOTH CONNECTION
 
C**   Initialization
      RCON=1.
      ITER = 0

C---  RCON Iteration loop ------------------------------------------
      DO
        RCONOLD=RCON
        ITER = ITER + 1
      
C        WRITE (0, '(A, I2, A, F12.6, A, F12.6)') 
C     >      'ITER=', ITER, '  VPAR2=', VPAR2, '  RCON=', RCON 
        IF (ITER >= 99) THEN
          WRITE (0, '(A)') '*** ITERATION FOR VPARs NOT CONVERGED!' 
          STOP             '*** FATAL ERROR IN INITVEL'
        ENDIF

C***    CHOSE THE BETA-LAW PARAMETERS SUCH THAT BOTH VELOCITIES AGREE
C***    AT THE CONNECTION POINT
        VCON = VMIN * EXP((RCON-1.)/HSCALE)
        CALL INITVELBETAPAR(RMAX, RCON, VCON, .TRUE.)
C        Q = (VCON/VMAX)**(1./BETA)
C        S = (RMAX - 1. + RCON) / 2.
C        T = (RMAX - Q * RCON) / (1.-Q) - RMAX * RCON
C        VPAR2 = SQRT(S*S + T) - S
C        VPAR1 = VMAX / (1.-1./(VPAR2+RMAX))**BETA
C        IF (BETA2FRACTION > 0.) THEN
C          VPAR2_2 = 1. - RCON   
C          VPAR1_2 = BETA2FRACTION * 
C     >             VFINAL/(1.-1./(VPAR2_2+RMAX))**BETA2
C        ENDIF

C***    Now it is additionally demanded that the velocity gradient
C***    is also smooth. The FUNCTION DELTAGR gives the difference 
C***    of the gradients from the outer minus the inner law 

C***    CHECK IF AT THE INNER BOUNDARY THE EXPONENTIAL LAW GRADIENT
C***    EXCEEDS ALREADY THE BETA-LAW GRADIENT
C***    Note: this check may only be performed in the first iteration
        IF (ITER <= 1) THEN
          IF (DELTAGR(1.) .LT. .0) THEN
            RCON=1.
            CALL REMARK('INITVEL: NO HYDROSTATIC DOMAIN ENCOUNTERED')
            bHydroStat = .FALSE.
            RETURN
          ENDIF
        ENDIF
 
C***    LOWER GUESS FOR RCON
        RCONMIN=AMAX1(1.,1.-VPAR2+1.E-10)
C***    Maximum of gradient in beta laws:
        RGRADMAX = 1. + (BETA-1.)/2. - VPAR2


C***    Loop: INCREASE RCONMIN STEPWISE TO ENSURE DELTAGR(RCONMIN) .GE. 0.
        !@Check: Kann dies eine Endlosschleife werden?
        DO
          IF (DELTAGR(RCONMIN) >= .0) EXIT
          DXX = 0.1 * HSCALE
          RCONMIN = RCONMIN + DXX
          IF (RCONMIN > RGRADMAX) THEN
            WRITE (0, '(A)') '*** INITVEL: NO CONNECTION POINT FOUND'
            STOP             '*** ERROR IN INITVEL'
          ENDIF
cc           WRITE (0, '(A, F12.6, A, F12.6)') 
cc     >           'RCONMIN INCREASED TO', RCONMIN, 
cc     >           '  DELTAGR=', DELTAGR(RCONMIN)
        ENDDO
 
C***    UPPER GUESS FOR RCON
        RCONMAX=RCONMIN
C**     Loop: 
        DO
           RCONMAX=RCONMAX+HSCALE
           IF (DELTAGR(RCONMAX) < .0) EXIT
           IF (RCONMAX >= RMAX) THEN
             RCON=RMAX
             CALL REMARK('INITVEL: NO BETA-LAW DOMAIN ENCOUNTERED')
             RETURN
           ENDIF
        ENDDO
 
C***    Now find new RCON where DELTAGR is zero:
        CALL REGULA (DELTAGR,RCON,.0,RCONMIN,RCONMAX,1.E-8)

        IF (ABS(RCON-RCONOLD) <= 1.E-7) EXIT
        
      ENDDO
C---  End of RCON Iteration loop --------------------------------------
 
      RETURN
      END
      SUBROUTINE INITVELBETAPAR (RMAX, RCONin, VCON, 
     >                           bUseMaxForPar1)
C*******************************************************************************
C***  RAW INITIALIZATION OF THE VELOCITY-FIELD PARAMETERS
C***  unifies code parts that were previously in 
C***   INITVEL, VELTHIN and ENSURETAUMAX
C*******************************************************************************

      IMPLICIT NONE

      REAL RMAX, RCONin, VCON
      LOGICAL bUseMaxForPar1, bPAR1AtCon

      REAL :: VFINAL, VMAX, BETA, BETA2, BETA2FRACTION, VMIN,
     >        VPAR1, VPAR2, HSCALE, VPAR1_2, VPAR2_2, RCON,
     >        Q, S, T, RM

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >            BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

       bPAR1AtCon = .NOT. bUseMaxForPar1

      RCON = RCONin

      VMAX = VFINAL * (1.-BETA2FRACTION)

      IF (BETA == .0) THEN
        Q = (Vcon/VMAX)**2
        VPAR2 = (RMAX**Q / RCON)**(1./(Q-1.))
        VPAR1 = VMAX / SQRT(ALOG(RMAX/VPAR2))
      ELSEIF (BETA < 0.) THEN
        Q = (VCON/VMAX)**(1./ABS(BETA))
        VPAR2 = (Q - 1.) / (Q/RMAX - 1./RCON)
        IF (bPAR1AtCon) THEN
          VPAR1 = Vcon / (1.-VPAR2/RCON)**(ABS(BETA))
        ELSE
          VPAR1 = VMAX / (1.-VPAR2/RMAX)**(ABS(BETA))         !INITVEL approach
        ENDIF
        IF (BETA2FRACTION > 0.) THEN
          VPAR2_2 = RCON
          VPAR1_2 = BETA2FRACTION * 
     >      VFINAL / (1.-VPAR2_2/RMAX)**(ABS(BETA2))     !note: 2BETA has always RMAX here
        ENDIF
      ELSE
        Q = (VCON/VMAX)**(1./BETA)
        S = (RMAX - 1. + RCON) / 2.
        T = (RMAX - Q * RCON) / (1.-Q) - RMAX * RCON
        VPAR2 = SQRT(S*S + T) - S
        IF (bPAR1AtCon) THEN
          RM = RCON
          IF (VPAR2 + RCON <= 1.) RM = 1. + ABS(VPAR2) + 1.E-10
          VPAR1 = Vcon/(1.-1./(VPAR2+RM))**BETA             !VELTHIN approach
        ELSE
          VPAR1 = VMAX / (1.-1./(VPAR2+RMAX))**BETA         !INITVEL approach
        ENDIF
        IF (BETA2FRACTION > 0.) THEN
          VPAR2_2 = 1. - RCON   
          VPAR1_2 = BETA2FRACTION * 
     >             VFINAL/(1.-1./(VPAR2_2+RMAX))**BETA2     !note: 2BETA has always RMAX here
        ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE INSTALL
C***********************************************************************
C***  Called from Main Programs
C***********************************************************************

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C*********************************************
C***  Hier Hauptschalter:
C***  INST = 1 : Cray
C***  INST = 2 : Potsdam DEC/UNIX
C***  INST = 3 : Potsdam SGI Origin 2000
C*********************************************

C     vvvvvvvv
      INST = 2
C     ^^^^^^^^

C******  Cray  **************************************
      IF (INST .EQ. 1) THEN

      OPSYS = 'CRAY'

C******  DEC  **************************************
      ELSE IF (INST .EQ. 2) THEN

      OPSYS = 'DEC/UNIX'

C******  SGI  **************************************
      ELSE IF (INST .EQ. 3) THEN

      OPSYS = 'SGI'

C****** ERROR  **************************************
      ELSE
      STOP 'ERROR IN INSTALL'
      ENDIF

      RETURN
      END
      FUNCTION ISRCHEQ(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHEQ, ISRCHNE - Searches a vector for the first element equal or
C***       not equal to a target

C***  SYNOPSIS
C***       index = ISRCHEQ (n, x, incx, target)

C***       index = ISRCHNE (n, x, incx, target)

C***  IMPLEMENTATION
C***       Cray PVP systems

C***  DESCRIPTION
C***       ISRCHEQ searches a real or integer vector for the first element that
C***       is equal to a real or integer target.

C***       ISRCHNE searches a real or integer vector for the first element that
C***       is not equal to a real or integer target.

C***       These functions have the following arguments:

C***       index  Integer.  (output)
C***              Index of the first element equal or not equal to target.
C***              If target is not found, n+1 is returned.
C***              If n <= 0, 0 is returned.

C***       n      Integer.  (input)
C***              Number of elements to be searched.

C***       x      Real or integer array of dimension  (n-1)*|incx|+1.  (input)
C***              Array x contains the vector to be searched.

C***       incx   Integer.  (input)
C***              Increment between elements of the searched array.

C***       target Real or integer.  (input)
C***              Value for which to search in the array.

C***  The Fortran equivalent code for ISRCHEQ is as follows:

      INTEGER X(*), TARGET

      J=1
      ISRCHEQ=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).EQ.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHEQ=I

      RETURN
      END

      FUNCTION ISRCHFGT(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHFGT - Searches a vector for the first element greater 
C***       to a target

C***  SYNOPSIS
C***  SEE ISRCHEQ

      REAL X(*), TARGET

      J=1
      ISRCHFGT=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).GT.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFGT=I

      RETURN
      END

      FUNCTION ISRCHNE(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHEQ, ISRCHNE - Searches a vector for the first element equal or
C***       not equal to a target

C***  SYNOPSIS
C***       index = ISRCHEQ (n, x, incx, target)

C***       index = ISRCHNE (n, x, incx, target)

C***  IMPLEMENTATION
C***       Cray PVP systems

C***  DESCRIPTION
C***       ISRCHEQ searches a real or integer vector for the first element that
C***       is equal to a real or integer target.

C***       ISRCHNE searches a real or integer vector for the first element that
C***       is not equal to a real or integer target.

C***       These functions have the following arguments:

C***       index  Integer.  (output)
C***              Index of the first element equal or not equal to target.
C***              If target is not found, n+1 is returned.
C***              If n <= 0, 0 is returned.

C***       n      Integer.  (input)
C***              Number of elements to be searched.

C***       x      Real or integer array of dimension  (n-1)*|incx|+1.  (input)
C***              Array x contains the vector to be searched.

C***       incx   Integer.  (input)
C***              Increment between elements of the searched array.

C***       target Real or integer.  (input)
C***              Value for which to search in the array.

C***  The Fortran equivalent code for ISRCHEQ is as follows:

      INTEGER X(*), TARGET

      J=1
      ISRCHNE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).NE.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHNE=I

      RETURN
      END

      SUBROUTINE JSTART (NF,XLAMBDA,KEY,ND,R,T,XJC,XJL,ELEVEL, N,EINST,
     $                   NDIM, INDNUP, INDLOW, LASTINDALL, LASTIND,
     >                   NAUTO, WSTABIL, R23, TAUROSS, 
     >                   LTESTART, BLACKEDGE, bOLDJ, XJCold, XLAMBDAold, Rold,
     >                   NFold, NDold, KRUDAUT)
C******************************************************************************
C***  START APPROXIMATION OF THE RADIATION FIELD
C***  XJC: CONTINUUM RADIATION FIELD
C***  XJL: LINE RADIATION FIELD
C******************************************************************************
 
      INTEGER, INTENT(IN) :: NFold, NDold, LASTIND, LASTINDALL, NAUTO, ND, NF
      REAL, DIMENSION(NF) :: XLAMBDA
      CHARACTER(8), DIMENSION(NF) :: KEY
      REAL, DIMENSION(ND) :: R, T, XJC, XJL, TAUROSS
      REAL, DIMENSION(NDIM) :: ELEVEL
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      INTEGER, DIMENSION(LASTINDALL) :: INDNUP, INDLOW
      INTEGER, DIMENSION(NAUTO) :: KRUDAUTO
      REAL, DIMENSION(NAUTO) :: WSTABIL
      REAL, DIMENSION(NDold) :: Rold
      REAL, DIMENSION(NFold) :: XLAMBDAold
      REAL, DIMENSION(NDold,NFold) :: XJCold
      CHARACTER(8) :: NAME
      REAL, DIMENSION(4) :: XJSP, TAUSP

      LOGICAL :: LTESTART, bOLDJ

C***  Find temperature at tau = 2/3
      CALL LIPO (T23, R23, T, R, ND)

C***  Find depth index inside tau = 1
      L1  = ISRCHFGT(ND,TAUROSS, 1, 1.5)
      IF (L1 .EQ. 0) L1 = ND-1

C***  Find depth index outside tau = 1/3
      L13 = ISRCHFGT(ND,TAUROSS, 1, 0.333333)
      IF (L13 .GT. 1) L13 = L13 - 1
      IF (L13 .EQ. 0) L13 = 1

C***  Construct vector with 4 elements for spline interpolation
C***   in order to smooth over the tau=1 discontinuity 
      IF (L13 .EQ. 1 .OR. L1 .GE. (ND-1)) THEN 
         NSP = 0
      ELSE
         NSP = 4
         ISP1 = L13 - 1       
         ISP2 = L13       
         ISP3 = L1       
         ISP4 = L1+1       
         TAUSP(1) = TAUROSS(ISP1)
         TAUSP(2) = TAUROSS(ISP2)     
         TAUSP(3) = TAUROSS(ISP3)     
         TAUSP(4) = TAUROSS(ISP4)     
      ENDIF

C***  CONTINUUM RADIATION FIELD XJC  *****************************************
C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS
C***  ( DEPTH VEKTOR AT EACH FREQUENCY POINT )
      DO 6 K=1,NF
      XLAM=XLAMBDA(K)
 
      IF (bOLDJ) THEN
C***     BRANCH FOR XJC FROM OLD MODEL
         DO L=1, ND
C***       Determine the corresponding depth point in the old model
C***       (This could be optimized with a more precise interpolation)
           lofind: DO LL=1, NDold
             IF (R(L) >= Rold(LL)) THEN
               Lold = LL
               EXIT lofind
             ENDIF
           ENDDO lofind
           WN = 1.E8/XLAM
           IF (XLAM < XLAMBDAold(1)) THEN
             XJC(L) = MIN(XJCold(L,1), BNUE(XLAM, T(L)))
           ELSEIF (XLAM > XLAMBDAold(NFold)) THEN 
             XJC(L) = MIN(XJCold(L,NFold), BNUE(XLAM, T(L))) 
           ELSE
             CALL XRUDI (XJC(L),WN,XJCold,XLAMBDAold,NDold,NFold,Lold)
           ENDIF
         ENDDO
      ELSEIF (LTESTART) THEN
C***     BRANCH FOR LTESTART
         DO L=1,ND
           XJC(L)=BNUE(XLAM,T(L))
           IF (R(L) > R23 .AND. XLAM < BLACKEDGE) THEN
             XJC(L) = XJC(L) * EXP(-20.*(2./3. - TAUROSS(L)))
           ENDIF
         ENDDO   
      ELSE
C***     BRANCH FOR GEOMETRICAL DILUTION OF BLACKBODY FIELD 
C***       AT PHOTOSPHERIC RADIUS
         BPHOT = BNUE(XLAM,T23)
cc         IF (XLAM .LT. BLACKEDGE) BPHOT = .0

         DO 20 L=1,ND
            IF (R23 .LT. R(L)) THEN
               W=0.5
               ARG=1. - R23*R23 / (R(L)*R(L))
               IF (ARG .GT. .0) W=0.5*(1.-SQRT(ARG))
               XJC(L) = BPHOT * W
cc naechste Zeile testweise
               IF (XLAM .LT. BLACKEDGE) XJC(L) = XJC(L) * 
     >             EXP(-20.*(2./3. - TAUROSS(L)))
            ELSE
               XJC(L)=BNUE(XLAM,T(L))
            ENDIF
   20    CONTINUE

C***  replaced by a new version with Spline interpolation, wrh  3-Mar-2006
         IF (NSP .EQ. 4) THEN
            XJSP(1) = XJC(ISP1)
            XJSP(2) = XJC(ISP2)
            XJSP(3) = XJC(ISP3)
            XJSP(4) = XJC(ISP4)
            DO L = L13+1, L1-1
               CALL SPLINPO (XJC(L), TAUROSS(L), XJSP, TAUSP, NSP)
            ENDDO
         ENDIF
      ENDIF
 
      WRITE (UNIT=NAME, FMT='(A3,I4,A1)') 'XJC', K, ' '
      CALL WRITMS(3,XJC,ND,NAME,-1,0, IERR)
    6 CONTINUE
C******************************************************************************
 
C***  LINE RADIATION FIELD XJL  ***********************************************
C***  ( DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND )
      DO 99 IND=1, LASTINDALL
C***    calc XLAM: Branch for normal lines and for superlevels
        IF (IND .LE. LASTIND .OR. IND .GT. LASTIND+NAUTO) THEN
           J=INDNUP(IND)
           I=INDLOW(IND)
           IF (EINST(I,J) .EQ. -2.) CYCLE
           XLAM=1.E8/(ELEVEL(J)-ELEVEL(I))
        ELSE
           INDDR = IND - LASTIND
           IF (KRUDAUT(INDDR) .EQ. 1) CYCLE
           XLAM = 1.E8 / WSTABIL(INDDR)
        ENDIF

C***    THIS VERSION: SAME APPROXIMATION AS FOR CONTINUUM
        IF (bOLDJ) THEN
C***      Branch for old radiation field: Just interpolate XJC from above
          DO L=1,ND
C***        Determine the corresponding depth point in the old model
C***        (This could be optimized with a more precise interpolation)
            lofind2: DO LL=1, NDold
              IF (R(L) >= Rold(LL)) THEN
                Lold = LL
                EXIT lofind2
              ENDIF
            ENDDO lofind2
            WN = 1.E8/XLAM
            CALL XRUDI (XJL(L),WN,XJCold,XLAMBDAold,NDold,NFold,Lold)
          ENDDO        
        ELSEIF (LTESTART) THEN
C***      BRANCH FOR LTESTART
          DO L=1,ND
            XJL(L)=BNUE(XLAM,T(L))
            IF (R(L) > R23 .AND. XLAM < BLACKEDGE) THEN
              XJL(L) = XJL(L) * EXP(-20.*(2./3. - TAUROSS(L)))
            ENDIF
          ENDDO
        ELSE
C***    BRANCH FOR GEOMETRICAL DILUTION OF BLACKBODY FIELD
          BPHOT = BNUE(XLAM,T23)

cc         IF (XLAM .LT. BLACKEDGE) BPHOT = .0

          DO 9 L=1,ND
            IF (R23 .LT. R(L)) THEN
               W=0.5
               ARG=1. - R23*R23 / (R(L)*R(L))
               IF (ARG .GT. .0) W=0.5*(1.-SQRT(ARG))
               XJL(L)=BPHOT*W
cc naechste Zeile testweise
               IF (XLAM .LT. BLACKEDGE) XJL(L) = XJL(L) * 
     >             EXP(-20.*(2./3. - TAUROSS(L)))
            ELSE
               XJL(L)=BNUE(XLAM,T(L))
            ENDIF
    9     CONTINUE

C***     Test output for the smoothing 
C***     For activation, set INDPLO=137 for HeII Lyman-alpha, for instance
         INDPLO = 0
         IF (IND .EQ. INDPLO) THEN
          OPEN (2, FILE='PLOT')
          NSPPLOT = L1 - L13 + 3
          CALL PLOTANFS (2,' ', ' ',
     $        'TAUROSS',
     $        'Radiation intensity',
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     >        TAUROSS(L13-1), XJL(L13-1), NSPPLOT, 
     >        'SYMBOL=2 SIZE=0.2 COLOR=4')
          CALL JSYMSET ('G2','TRANSFER')
         ENDIF

C***     replaced by a new version with Spline interpolation, wrh  3-Mar-2006
         IF (NSP .EQ. 4) THEN
            XJSP(1) = XJL(ISP1)
            XJSP(2) = XJL(ISP2)
            XJSP(3) = XJL(ISP3)
            XJSP(4) = XJL(ISP4)
            DO L = L13+1, L1-1
               CALL SPLINPO (XJL(L), TAUROSS(L), XJSP, TAUSP, NSP)
            ENDDO
         ENDIF

C***     Continue testplot (now after smoothing)
         IF (IND .EQ. INDPLO)
     >    CALL PLOTCONS (2, 
     >        TAUROSS(L13-1), XJL(L13-1), NSPPLOT, 
     >        'SYMBOL=8 SIZE=0.2 COLOR=2')

      ENDIF

C***  WRITE RADIATION FIELD TO THE MODEL FILE 
      IF (IND <= 9999) THEN
        WRITE (UNIT=NAME, FMT='(A3,I4,A1)') 'XJL',IND,' '
      ELSE
        WRITE (UNIT=NAME, FMT='(A3,I5)') 'XJL', IND
      ENDIF
      CALL WRITMS (3,XJL,ND,NAME,-1,0, IERR)
   99 CONTINUE
C***********************************************************************

      RETURN
      END

      SUBROUTINE JSYMSET (SYMBOL,WORD)
C***********************************************************************
C***  USER -DEFINED UNICOS VERSION OF THE COS LIBRARY ROUTINE JSYMSET
C***  OPENS A FILE NAMED AS THE SPECIFIED COS SYMBOL AND
C***  WRITES THE CHARACTER STRING "WORD" INTO THIS FILE
C***********************************************************************
      CHARACTER SYMBOL*(*), WORD*(*)

      OPEN (66,FILE=SYMBOL, STATUS='UNKNOWN')
      WRITE (66,'(A)') WORD
      CLOSE (66)

      RETURN
      END
      SUBROUTINE KSIGMA(SIGMAK,SIGMATHK,EDGEK,WAVENUM,SEXPOK)
C******************************************************************
C***  CLCULATE KSIGMA, THE FREQUENZ-DEPENDENT K-SHELL CROSS SECTION
C***  CALLED FROM SOBROUTINE COOP
C******************************************************************

     
      X=EDGEK/WAVENUM
      
      SIGMAK = SIGMATHK * 1.E-18 * X ** SEXPOK
      
      RETURN
      END
      SUBROUTINE LIPO (F,X,FI,XI,N)
C***********************************************************************
C***  LINEAR INTERPOLATION: FOR GIVEN X, FIND F(X) FROM A TABLE FI(XI)
C ------ INDEXSUCHE DURCH BISECTION ------
C***********************************************************************

      DIMENSION FI(N),XI(N)

      NA=1
      A=XI(1)
      NB=N
      B=XI(N)
      IF((X-A)*(X-B).GT..0) STOP 'ERROR IN SUBR. LIPO'
    1 IF((NB-NA).EQ.1) GOTO 2
      NH=(NA+NB)/2
      H=XI(NH)
      IF((X-A)*(X-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
 
    2 P=(X-A)/(B-A)
      F=P*FI(NB)+(1.-P)*FI(NA)
      RETURN
      END
      SUBROUTINE LOADOLDJ(XJCold, XLAMBDAold, NDold, NFold, 
     >                    NDDIM, NFDIM, hOldMODEL)
C***********************************************************************
C***  read out the XJC radiation field from an old model
C***
C***  called from WRSTART
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDold, NDDIM, NFDIM, hOldMODEL
      INTEGER, INTENT(INOUT) :: NFold
      REAL, DIMENSION(NFDIM) :: XLAMBDAold
      REAL, DIMENSION(NDDIM * NFDIM) :: XJCold

      CHARACTER(8) :: NAME
      INTEGER :: K, IERR


      CALL READMS (hOldMODEL, NFold     ,     1, 'NF      ', IERR)
      CALL READMS (hOldMODEL, XLAMBDAold, NFold, 'XLAMBDA ', IERR)
      DO K=1, NFold
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL READMS(hOldMODEL, XJCold(1+NDold*(K-1)), NDold, NAME, IERR)
      ENDDO


      RETURN
 
      END


        SUBROUTINE LOWERCASE (TEXT)
C ***
C ******************************************************************************
C ***   converts string TEXT to lower case.
C ******************************************************************************

      CHARACTER*(*) TEXT, TAB(2)*26
      DATA TAB /'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 
     >	        'abcdefghijklmnopqrstuvwxyz'/

      DO L = 1, LEN(TEXT)
         J = INDEX( TAB(1), TEXT(L:L) )
         IF( J .NE. 0 ) TEXT(L:L) = TAB(2)(J:J)
      ENDDO

      RETURN
      END
      SUBROUTINE LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
C***********************************************************************
C***  POPULATION NUMBERS IN THERMODYNAMIC EQUILIBRIUM FOR ALL ELEMENTS
C***********************************************************************
 
      DIMENSION ENLTE(N), WEIGHT(N),NCHARG(N),EION(N),ELEVEL(N),NOM(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1/1.4388/
C***  C2 = FACTOR IN SAHA EQ.   ( C.F. MIHALAS P.113 )
      DATA C2/2.07E-16/
 
      T32=TL*SQRT(TL)
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
      DO 9 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      ENLTE(NFIRNA)=1.
      DO 1 J=NFIRNA+1,NLANA
      IF(NOM(J) .NE. NOM(NFIRNA)) STOP'LTEPOP:WRONG ELEMENT MEMBERSHIP'
      IF (NCHARG(J) .EQ. NCHARG(J-1) ) THEN
C***     BOLTZMANN FACTOR
         ENLTE(J) = EXP(C1*(ELEVEL(J-1)-ELEVEL(J))/TL) 
     >             * WEIGHT(J)/WEIGHT(J-1) * ENLTE(J-1)

      ELSE IF (NCHARG(J) .EQ. NCHARG(J-1)+1 ) THEN
C***     SAHA FACTOR
         ENLTE(J) = EXP(C1*(ELEVEL(J-1)-ELEVEL(J)-EION(J-1))/TL) 
     >       * T32/ENE/C2 * WEIGHT(J)/WEIGHT(J-1) * ENLTE(J-1)
      ELSE
         STOP 'LTEPOP: INVALID CHARGE DIFFERENCE'
      ENDIF
    1 CONTINUE
 
C***  NORMALIZATION
      SUM=.0
      DO 2 J=NFIRNA,NLANA
        ENLTE(J) = MAX(1.E-200,ENLTE(J))
    2   SUM = SUM+ENLTE(J)
      SUM = SUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
    3   ENLTE(J) = ENLTE(J)/SUM
 
    9 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
      SUBROUTINE MGOETZ (XMSTAR, TEFF, RSTAR, XHY, WRTYPE)     
C***********************************************************************
C***  Calculates Stellar Mass from the Luminosity (mass-luminosity relation)
C***  according to Graefener et al. (2011: A&A 535, 56)
C***  Notes: 
C***  - masses are for chemically homogeneous stars, and thus upper limits
C***  - two different M-L relations: 
C***       WRTYPE='WN': He-burning stars (H-free)
C***       WRTYPE='OB': Hydrogen-burning stars
C***  - resulting mass XMSTAR is in units of MSUN  
C***********************************************************************

      IMPLICIT NONE

      REAL, INTENT(OUT) :: XMSTAR
      REAL, INTENT(IN) :: TEFF, RSTAR, XHY
      
      REAL :: F1, F2, F3, F4, F5, F6, F7, F8, F9, F1r, F2r, F3r, F4r,
     >        XLSTAR, XLSTARS, XLOGL, XLOGM, Radikant
      
      !Constants
      REAL, PARAMETER :: STEBOL = 5.6705E-5   !STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: PI4 = 12.5663706144  !PI4 = 4*PI
      REAL, PARAMETER :: XLSUN = 3.85E33      !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33     !Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      CHARACTER WRTYPE*2

C***  Fit parameters for hydrogen-burning stars, Eqs. 11 and 12, Table A.1  
      F1=4.026 
      F2=4.277
      F3=-1.0
      F4=25.48
      F5=36.93
      F6=-2.792
      F7=-3.226
      F8=-5.317
      F9=1.648

C***  Fit parameters for helium-burning stars, Eq. 13, Table A.1
      F1r=3.997
      F2r=-1.0
      F3r=25.83
      F4r=-3.268

            
      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN
      XLOGL = ALOG10(XLSTARS)

C***  Hydrogen-burning star 
      IF (WRTYPE == 'OB') THEN
        IF (XHY .LT. 0.1) GOTO 90
        Radikant = F4 + F5 * XHY + F6 * XHY**2 + (F7+ F8*XHY) * XLOGL
        XLOGM = ( F1+F2*XHY + F3*SQRT( Radikant ) ) / ( 1.+F9*XHY )

C***  Helium-burning star
      ELSEIF (WRTYPE == 'WN') THEN
        Radikant = F3r + F4r * XLOGL
        XLOGM = F1r + F2r * SQRT( Radikant )
      ELSE
        GOTO 91
      ENDIF
      
      XMSTAR = 10.**XLOGM

      RETURN

C***  ERROR BRANCHES
   90 WRITE (0,*) '*** ERROR: M-L relation only valid for X_H > 0.1' 
      GOTO 99

   91 WRITE (0,*) '*** ERROR: M-L relation (Graefener version) not '
     >             // 'valid for type ' // WRTYPE 

   99 STOP '*** INTERNAL ERROR in Subr. MGOETZ'

      END
      SUBROUTINE MLANGER (XMSTAR, TEFF, RSTAR, YHE, WRTYPE) 
C***  Calculates Stellar Mass for a H-free Star
C***  from the Luminosity (mass-luminosity relations)
C***  Note: the resulting mass XMSTAR is given in gramms !!  

      CHARACTER*2 WRTYPE
      
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      DATA STEBOL / 5.6705E-5 /
C***  PI4 = 4*PI
      DATA PI4 / 12.5663706144 /
C***  XLSUN = Solar Luminosity (CGS-Units)
      DATA XLSUN / 3.85E33 /
C***  XMSUN = Solar Mass (g)
      DATA XMSUN / 1.989E33 /

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

C***  Mass-Luminosity relation from Langer N., A&A 210, 93 (1989)


      XLOGL = ALOG10(XLSTARS)

      IF (WRTYPE .EQ. 'WC') THEN
         A0 = -0.487870 * YHE
         A0 = A0 + 2.971463
         A1 = +0.434909 * YHE
         A1 = A1 + 2.771634
         A2 = -0.093793 * YHE
         A2 = A2 - 0.487209
         P = A1/A2
         Q = (A0 - XLOGL)/A2
         PQS = P*P/4. - Q
         IF (PQS .GE. 0.) THEN
cc            XLOGMP = -P/2. + SQRT(PQS)
            XLOGMM = -P/2. - SQRT(PQS)
            XMSTARS = 10.**XLOGMM
            XMSTAR = XMSTARS * XMSUN
         ELSE
            WRITE (0,*) '*** ERROR: Langer-Formula failed!'
            STOP '*** ERROR IN SUBROUTINE MLANGER'
         ENDIF
      ELSE IF (WRTYPE .EQ. 'WN') THEN
         XLOGM = - 0.158206 - 0.053868*XLOGL + 0.055467*XLOGL*XLOGL
         XMSTARS = 10.**XLOGM
         XMSTAR = XMSTARS * XMSUN
      ELSE
         WRITE (0,*) '*** ERROR: ILLEGAL WRTYPE: ', WRTYPE
         STOP '*** ERROR IN SUBROUTINE MLANGER'
      ENDIF

c      WRITE (0, *) '*** MASS FROM LUMINOSITY ***'
C      WRITE (0, *) XLOGMP, XLOGMM
C      WRITE (0, *) XLSTARS, 10**(A0 + A1*XLOGMM + A2*XLOGMM*XLOGMM)
c      WRITE (0, *) 'L=', XLSTARS, 'L_SUN'
c      WRITE (0, *) 'M=', XMSTARS, 'M_SUN'
      

C***  Zeta(M) correcture from Heger A., A&A 315, 421

C      ZETA = 1.8548 - 0.039711*MSTAR +6.1946E-4*MSTAR*MSTAR


C***  Mass-Radius relation from Langer N., A&A 210, 93 (1989)

C      B0 = -0.114128 * YHE
C      B0 = B0 - 0.497381
C      B1 = +0.208654 * YHE
C      B1 = B1 + 0.371859
C      B2 = -0.156504 * YHE
C      B2 = B2 + 0.156504

C      RTH = B0 + B1*XLOGM + B2*XLOGM2

      RETURN
      END
      SUBROUTINE MY_CLOCK (STRING)

      CHARACTER*8 STRING

      CALL TIME(STRING)

      RETURN
      END
      SUBROUTINE MY_DATE (STRING)

      CHARACTER STRING*8, STRING2*10

      CALL DATE(STRING2)

C***  Create STRING as yy/mm/dd
      STRING = '  /  /  '

C***  Year
      STRING(1:2) = STRING2(8:9)

C***  Month
      IF      (STRING2(4:6) .EQ. 'Jan') THEN
        STRING(4:5) = '01'
      ELSE IF (STRING2(4:6) .EQ. 'Feb') THEN
        STRING(4:5) = '02'
      ELSE IF (STRING2(4:6) .EQ. 'Mar') THEN
        STRING(4:5) = '03'
      ELSE IF (STRING2(4:6) .EQ. 'Apr') THEN
        STRING(4:5) = '04'
      ELSE IF (STRING2(4:6) .EQ. 'May') THEN
        STRING(4:5) = '05'
      ELSE IF (STRING2(4:6) .EQ. 'Jun') THEN
        STRING(4:5) = '06'
      ELSE IF (STRING2(4:6) .EQ. 'Jul') THEN
        STRING(4:5) = '07'
      ELSE IF (STRING2(4:6) .EQ. 'Aug') THEN
        STRING(4:5) = '08'
      ELSE IF (STRING2(4:6) .EQ. 'Sep') THEN
        STRING(4:5) = '09'
      ELSE IF (STRING2(4:6) .EQ. 'Oct') THEN
        STRING(4:5) = '10'
      ELSE IF (STRING2(4:6) .EQ. 'Nov') THEN
        STRING(4:5) = '11'
      ELSE IF (STRING2(4:6) .EQ. 'Dec') THEN
        STRING(4:5) = '12'
      ELSE
        WRITE (0,*) 'MONTH NOT RECOGNIZED'
        WRITE (*,*) 'STRING2(4:6)=',STRING2(4:6)
        STOP 'ERROR IN MY_DATE'
      ENDIF

C***  day
      STRING(7:8) = STRING2(1:2) 


      RETURN
      END
      SUBROUTINE OPAGREY (OPARL,EN,TL,RNEL,ENTOTL,RSTAR,N,
     $                   NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,
     $                   NOM,EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,
     $                   MAXION,SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $                   KONTNUP,KONTLOW,LASTKON,RL,XDATA)
C***********************************************************************
C***  --- ONLY CALLED FROM SUBR. GREY ---
C***  FAST VERSION (VECTORIZED "COOPFRQ"!!!) OF ORIGINAL SUBR. OPAROSS:
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT DEPTH POINT L
C***  FOR GIVEN POPNUMBERS
C***********************************************************************
 
      DIMENSION EN(N)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF),EXPFAC(NF),OPAC(NF),ETAC(NF)
      DIMENSION SIGMAFF(NF,0:MAXION)
 
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /

C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN P.100 )
      DATA CFF / 1.370E-23 /

C***  SIGMAE = ELECTRON SCATTERING CROSS SECTION ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /

C***  PRE-CALCULATE EXPONENTIAL FACTORS AND FREE-FREE CROSS SECTIONS 
C***  FOR THE TEMPERATURE OF THE CURRENT DEPTH POINT
      TLOG=ALOG10(TL)
      ROOTTL=SQRT(TL)
      DO 1 K=1,NF
      XLAM=XLAMBDA(K)
      W=1.E8/XLAM
      EXPFAC(K)=EXP(-C1*W/TL)
      W3=W*W*W
      PRESIG=CFF/W3/ROOTTL
      XLAMLOG=ALOG10(XLAM)
      SIGMAFF(K,0)=0.0
      DO 1 ION=1,MAXION
      CALL GFFLOG (GIII,ION,XLAMLOG,TLOG)
      SIGMAFF(K,ION)=PRESIG*FLOAT(ION*ION)*GIII
    1 CONTINUE

      CALL COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $             WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $             RNEL,TL,SIGMAFF,MAXION,RL,XDATA,
     $             SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $             KONTNUP,KONTLOW,LASTKON,OPATHOM)

C***  FREQUENCY INTEGRATION
C***  NOTE: FOR NUMERICAL REASONS, THE NORMALIZATION CONSTANT Q IS
C**         ALSO INTEGRATED NUMERICALLY!
      SUM=.0
      Q=.0
C***  THOMSON SCATTERING OPACITY
      OPATH=RNEL*SIGMAE*ENTOTL*RSTAR
      DO 2 K=1,NF
C***  DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
      DBDT=DBNUEDT(XLAMBDA(K),TL)
C***  CONSIDERING THE THOMSON OPACITY
      OPAC(K)=OPAC(K)+OPATH
      SUM=SUM+DBDT*FWEIGHT(K)/OPAC(K)
      Q=Q+DBDT*FWEIGHT(K)
    2 CONTINUE
 
      OPARL=Q/SUM      

      RETURN
      END
      SUBROUTINE OPAROSS (OPARL,EN,TL,RNEL,ENTOTL,RSTAR,NDIM,N,
     $                   LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                   ALPHA,SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     $                   MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,RL,
     $                   KONTNUP,KONTLOW,LASTKON,POPMIN)
C***********************************************************************
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT DEPTH POINT L
C***  FOR GIVEN POPNUMBERS
C***********************************************************************
 
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION EN(NDIM)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
      CHARACTER*10 LEVEL(N),MAINPRO,MAINLEV
      REAL :: POPMIN
 
C***  FREQUENCY INTEGRATION ***
C***  NOTE: FOR NUMERICAL REASONS, THE NORMALIZATION CONSTANT Q IS
C**     ALSO INTEGRATED NUMERICALLY!
      SUM=.0
      Q=.0
C***  USE INTEGER VARIABLE 'LYMP' INSTEAD OF INTEGER CONSTANT "1" 
C***  (OLD VERSION) TO PREVENT RUNTIME ABORT ON Y-MP SYSTEMS
C***  (cft77 5.0.5.0 OR segldr 8.0c ERROR???)
C***  [ U. WESSOLOWSKI, 29-SEP-1993 20:39:49 ]
      LYMP = 1
      DO 1 K=1,NF
      XLAM=XLAMBDA(K)
      CALL COOP (XLAM,LYMP,TL,RNEL,EN,POPMIN,ENTOTL,RSTAR,OPA,ETA,
     >           THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     >           NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,
     $          ELEVEL,EION,EINST,ALPHA,SEXPO,
     $          ADDCON1, ADDCON2, ADDCON3, 
     $          IGAUNT,SIGMATHK,SEXPOK,EDGEK,0,DUMMY,DUMMY,RL,
     $          KONTNUP,KONTLOW,LASTKON,0.)
C***  DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
      DBDT=DBNUEDT(XLAM,TL)
      SUM=SUM+DBDT*FWEIGHT(K)/OPA
      Q=Q+DBDT*FWEIGHT(K)
    1 CONTINUE
 
      OPARL=Q/SUM      

      RETURN
      END
      SUBROUTINE OPENMS(ICHANNEL, IADR, MAXADR, IFNAME, IERR)
C************************************************************
C***  ROUTINE BY LARS KOESTERKE      8-Sep-1995 15:49:02
C************************************************************

      INTEGER :: IFNAME
      CHARACTER(8) :: FNAME, CSTAT     !IMPORTANT: only the first 7 characters can be used!

      IF ((IFNAME /= 0) .AND. (IFNAME /= 1)) THEN
        WRITE(UNIT=FNAME, FMT='(A8)') IFNAME
      ELSE
        FNAME = '        '
      ENDIF
      CSTAT = 'AUTO'

      CALL CMSSTORE (ICHANNEL, IADR, MAXADR, CSTAT, FNAME, 
     >              DUMMY, IDUMMY, 'OPEN', IERR)

      RETURN
      END
      SUBROUTINE PGRID (NPDIM,NP,ND,R,P, NC)
C***********************************************************************
C***  GRID OF IMPACT-PARAMETER POINTS
C***********************************************************************

      DIMENSION  P(NPDIM),R(ND)

C***  NC = NUMBER OF CORE-INTERSECTING RAYS
      NP=ND+NC
      IF (NP.GT.NPDIM) STOP 'PGRID: TOO MANY IMPACT-PARAMETER POINTS'

C***  CORE RAYS EQUALLY SPACED
      D=1./FLOAT(NC)
      DO 1 J=1,NC
    1 P(J)=(J-1)*D
      DO 2 L=1,ND
      J=NP+1-L
    2 P(J)=R(L)

      RETURN
      END
      SUBROUTINE PHOTOCS (SIGMA,SIGMATH,EDGE,WAVENUM,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,KON)
C***********************************************************************
C***  CALCULATES SIGMA(NUE), THE FREQUENCY DEPENDENT PHOTO CROSS SECTION
C***  THIS ROUTINE IS ONLY CALLED FROM BFCROSS, COOP AND CMFCOOP
C***********************************************************************

      CHARACTER(8), DIMENSION(1) :: IGAUNT
      DIMENSION ALPHA(1),SEXPO(1)
      DIMENSION ADDCON1(1), ADDCON2(1), ADDCON3(1)
C***  THE FOLLOWING DATA ARE FOR MIHALAS' GAUNT FACTOR FIT ( HY AND HE II, N=1) 
      DATA A0,A1,A2,A3,AM1,AM2 /
     > 1.2302628,  -3.1927214E-2, 8.9105122E-4, -1.1544111E-5,
     > -0.50812150, 0.10631895 /
      DATA IWARN_PIKB12 / 0 / 
 
      X=EDGE/WAVENUM
      XINV=1./X

C***  The calling programms should not call PHOTOCS beyond the ionisation edge,
C***  althogh this is not a catastrophy. 
C***  The followong warning might be switched off if annoying.
      IF (WAVENUM .LT. EDGE) THEN
         WRITE (0,'(a)') '*** WARNING: PHOTOCS CALLED OUTSIDE EDGE FREQUENCY'
         SIGMA = .0
         RETURN
      ENDIF


C***  VARIABLE IGAUNT IS MISUSED TO CARRY THE KEYWORD FOR NEW
C***  PHOTOIONIZATION CROSS SECTIONS: 
C***  'KOESTER': KOESTER ET AL. 1985, A+A 149, 423
C***   =======
C***             FIT COEFFICIENTS FOR THE MODIFIED (]) FORMULA: 
C***                     SIGMATH = F(A0)
C***                  ALPHA(KON) = A1
C***                  SEXPO(KON) = A2
      IF (IGAUNT(KON) .EQ. 'KOESTER') THEN
          XLN=ALOG(1.E8/WAVENUM)
          XLN2=XLN*XLN
          X0LN=ALOG(1.E8/EDGE)
          X0LN2=X0LN*X0LN
          SIGMA=SIGMATH*X**ALPHA(KON)*EXP(SEXPO(KON)*(XLN2-X0LN2))
C***  'BUTLER..': K. BUTLER (MUNICH), PRIVATE COMMUNICATION
C***   ========
C***  '......12': FORMULA 12 FROM PROGRAM DETAIL
C***              MODIFIED FORM OF THE SEATON FORMULA
      ELSE IF (IGAUNT(KON) .EQ. 'BUTLER12') THEN
          SIGMA=SIGMATH*X**(ALPHA(KON)+SEXPO(KON)*ALOG(X))
C***  'DETAILN3': EXTENDED VERSION (6 COEFFICIENTS!!!) OF FORMULA 12
C***   ========   (SEE "BUTLER12")
      ELSE IF (IGAUNT(KON) .EQ. 'DETAILN3') THEN
          NIII4NO=INT(SEXPO(KON))
C***      PREVENT LARGE PHOTOIONIZATION CROSS SECTIONS FOR SPECIFIED 
C***      N III LEVELS AT X-RAY FREQUENCIES (XLAMBDA <= 30A)
          IF ((NIII4NO .EQ. 4 .OR. NIII4NO .EQ. 5 .OR. NIII4NO .EQ. 6)
     $        .AND. (X .LT. 0.05)) THEN
              SIGMA=0.0
          ELSE
              SUMI36=PHOTON3(X,NIII4NO)
              SIGMA=SIGMATH*X**(ALPHA(KON)+SUMI36)
          ENDIF

      ELSE IF (IGAUNT(KON) .EQ. 'PIKB12  ') THEN
C***          XLN = ALOG(X)
C***          SIGMA = SIGMATH*X**(ALPHA(KON)+XLN*(SEXPO(KON)+XLN*(
C***     >            ADDCON1(KON)+XLN*(ADDCON2(KON)+XLN*ADDCON3(KON)))))
C***      changed by wrh 14-Mar-2005 11:43:35
C***      Corresponding data are in our files for C II and C III
C***      The origin of this fit-formula is unclear. It is rubbish:
C***      The exponent of X:=nue_o/nue is a polynomial of ln(x). 
C***      Taking ln of the whole formula, gives:
C***      ln(SIGMA) = ln (SIGMATH) + ALPHA*ln(X) + SEXPO*(ln(X))**2
C***                               + ADDCON1*(ln(X))**3 ...
C***      As X approaches zero for high frequencies, 
C***      ln(x) becomes -infinity, and the exponent goes towards plus 
C***      infinity when the highest non-zero coefficient (mostly ADDCON1) 
C***      is NEGATIVE - which is the case in all (but one) data!
C***      Therefore, I replace this formula by a hydrogenic slope. 
C***      Corresponding data should be avoided. Therefore, a warning is 
C***      issued when PIKB12 is used for the first time. 
              SIGMA = SIGMATH * X**3
              IF (IWARN_PIKB12 .EQ. 0) THEN
                 IWARN_PIKB12= IWARN_PIKB12 + 1
                 WRITE (0, *) '*** WARNING issued from PHOTOCS:'
                 WRITE (0, *) '*** Obsolete Formula PIKB12 ' 
     >            //'replaced by hydrogenic slope ***'
              ENDIF

C***  'OPAPROIX': ATOMIC DATA FOR OPACITY CALCULATIONS (OPACITY PROJECT):
C***  ==========  IX. The lithium isoelectronic sequence
C***              (Peach, Saraph & Seaton 1988, 
C***               J. Phys. B: At. Mol. Opt. Phys. 21, 3669)
      ELSE IF (IGAUNT(KON) .EQ. 'OPAPROIX') THEN
          XLOG = ALOG10(XINV)
          SIGMA = SIGMATH*XINV**(ALPHA(KON)+XLOG*(SEXPO(KON)+
     +                                      XLOG*ADDCON1(KON)))
      ELSE
C***  DEFAULT - OLD VERSION: SEATON / HYDROGENIC
C***  =======
          IF (ALPHA(KON) .NE. .0) THEN
C***      ALPHA(KON) DEFINED: SEATON FORMULA
            SIGMA=SIGMATH * X**SEXPO(KON)*(ALPHA(KON)+(1.-ALPHA(KON))*X)
          ELSE
C***      ALPHA(KON) NOT DEFINED: HYDROGENIC EXPONENT NUE**(-3)
            SIGMA=SIGMATH*X*X*X
          ENDIF
      ENDIF
 
C***  BOUND-FREE GAUNT FACTORS ARE CALCULATED DEPENDING ON KEYWORD IGAUNT
 
C***  POLYNOMIAL FIT FOR GII(N=1) FROM MIHALAS 1967 APJ 149,P187
C***  THIS MIHALAS GAUNT FACTOR IS ONLY VALID FOR GROUND STATE N=1 ]
      IF (IGAUNT(KON) .EQ. 'MIHALAS' ) THEN
            IF (X .GT. 0.055) THEN
C***           GOOD ACCURACY ONLY FOR WAVELENGTHS (OF HYDROGEN!) FROM
C***           THRESHOLD TO 50 A:
               GAUNT=A0+X*(AM1+X*AM2)+(A1+(A2+A3*XINV)*XINV)*XINV
            ELSE
               GAUNT=1.
            ENDIF
            SIGMA=SIGMA*GAUNT
      ENDIF
 
C***  GII FIT FROM SEATON (1960), REP. PROG. PHYS. 23, P. 313
C***  using FORMULA (2.4) from page 316
C***     note that DEN := ( n * (u + 1) )^(-2/3)
C***  THIS FORMULA IS VALID FOR ALL HYDROGENIC LEVELS,
C***  I.E.  H I, HE II, ..., C IV, N V, O VI
C***  VARIABLE SEXPO IS MISUSED TO CARRY THE MAIN QUANTUM NUMBER
      IF (IGAUNT(KON) .EQ. 'SEATON' ) THEN
            IF (SEXPO(KON) .LE. .0) THEN
               CALL REMARK ('MAIN QUANTUM NUMBER UNDEFINED')
               STOP 'ERROR IN SUBR. PHOTOCS : SEXPO0'
               ENDIF
            U=XINV - 1.
            DEN=(X/SEXPO(KON))**0.666666666666
            GAUNT=1. + 0.1728 * (U-1.) * DEN -
     -            0.0496 * (U*(U+1.333333333333)+1.) * DEN * DEN
            SIGMA=SIGMA*GAUNT
            ENDIF
 
C***  PREVENT NEGATIVE PHOTOIONIZATION CROSS SECTIONS:
      IF (SIGMA .LT. 0.0)  SIGMA = 0.0

C***  PREVENT LARGE PHOTOIONIZATION CROSS SECTIONS:
      IF (SIGMA .GT. 10.*SIGMATH)  THEN
         WRITE (0,*) 'WARNING : VERY HIGH PHOTOIONISATION CROSS ',
     >               'SECTION DETECTED; SET TO ZERO'

cc        write (0,'(a,2(e15.8,2x))') 'sigma, sigmath=', sigma, sigmath
cc        write (0,'(a,a8)') 'igaunt=',igaunt(kon)
cc        write (0,*) 'kon, wavenum=',kon, wavenum
cc        write (0,*) 'alpha, sexpo...=', alpha(KON), sexpo(KON), 
cc     >              addcon1(KON), addcon2(KON), addcon3(KON)
cc        write (0,*) 'ERROR IN SUBR. PHOTOCS : SIGMA'
cc        STOP 'ERROR IN SUBR. PHOTOCS : SIGMA'

         SIGMA = 0.

      ENDIF

      RETURN
      END
      FUNCTION PHOTON3 (X,NO)
C***********************************************************************
C***  CALCULATION OF PART OF THE SUM (I= 3 TO 6) IN THE EXTENDED VERSION
C***  OF FORMULA 12 (SEE DESCRIPTION OF PROGRAM "DETAIL") FOR RBF CROSS
C***  SECTIONS
C***  ---  CALLED FROM PHOTOCS: IF (IGAUNT(LOW) .EQ. 'DETAILN3')  ---
C***  NECESSARY FOR PHOTOIONISATION CROSSSECTIONS OF 10 LEVELS OF THE
C***  QUARTET SYSTEM IN NITROGEN III
C***  THE LAST 4 COEFFICIENTS ARE STORED IN DATA STATEMENTS!!!
C***********************************************************************

C***  NUMBER OF CONSIDERED LEVELS:
      PARAMETER ( N3QUART = 10 )

      DIMENSION AI(3:6,N3QUART)
C***  DETAIL LEVEL:     A4P1  ==>  N 32P2P4.2               (1)
      DATA (AI(I,1),I=3,6) /
     $ 6.37367E-2, 0.23060, 0.24562, 6.99181E-2 /
C***  DETAIL LEVEL:     A4S1  ==>  N 32P3S4.6               (2)
      DATA (AI(I,2),I=3,6) /
     $ -0.22891, -0.13430, -3.79655E-3, 5.93594E-3 /
C***  DETAIL LEVEL:     A4P2  ==>  N 33S'P412               (3)
      DATA (AI(I,3),I=3,6) /
     $ -0.56510, -0.30266, -0.10805, 3.69164E-3 /
C***  DETAIL LEVEL:     A4D1  ==>  N 33P'D415               (4)
      DATA (AI(I,4),I=3,6) /
     $ -0.50098, -0.32843, -0.15584, -3.08269E-2 /
C***  DETAIL LEVEL:     A4S2  ==>  N 33P'S417               (5)
      DATA (AI(I,5),I=3,6) /
     $ -0.55232, -0.40950, -0.20078, -3.97182E-2 /
C***  DETAIL LEVEL:     A4P3  ==>  N 33P'P418               (6)
      DATA (AI(I,6),I=3,6) /
     $ -0.52906, -0.41945, -0.21380, -4.24112E-2 /
C***  DETAIL LEVEL:     A4F1  ==>  N 33D'F422               (7)
      DATA (AI(I,7),I=3,6) /
     $ 4.06244E-2, -5.55455E-2, 4.59197E-2, 1.38044E-2 /
C***  DETAIL LEVEL:     A4D2  ==>  N 33D'D423               (8)
      DATA (AI(I,8),I=3,6) /
     $ 7.02743E-2, 2.16001E-2, 7.56214E-2, 1.59807E-2 /
C***  DETAIL LEVEL:     A4P4  ==>  N 33D'P424               (9)
      DATA (AI(I,9),I=3,6) /
     $ 4.44634E-2, 8.28837E-2, 9.21828E-2, 1.73687E-2 /
C***  DETAIL LEVEL:     A4P5  ==>  N 34S'P431               (10)
      DATA (AI(I,10),I=3,6) /
     $ -0.13676, 0.17396, 0.12893, 2.99775E-2 /

C***  ERROR STOP:
      IF ((NO .LT. 1) .OR. (NO .GT. N3QUART)) STOP 'ERROR'

C***  CALCULATION OF PART OF THE SUM (I= 3 TO 6):
      XLN=ALOG(X)
      PHOTON3=XLN*(AI(3,NO)+XLN*(AI(4,NO)+XLN*(AI(5,NO)+AI(6,NO)*XLN)))

      RETURN
      END
      SUBROUTINE PLOTANF (KANAL,NPLOT,NHEAD,NX,NY,
     $ B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6,
     $ X,Y,N,ISYMBOL )
C***********************************************************************
C***  DIESE ROUTINE BEREITET EIN NEUES PLOTFILE ZUR UEBERTRAGUNG VOR
C***  DIE ANGEGEBENEN WERTE BESCHREIBEN DEN PLOTKASTEN
C***  Use NPLOT = '' in CALL PLOTANF to skip the "PLOT:" line
C***********************************************************************
      CHARACTER*(*) NPLOT,NHEAD,NX,NY
      INTEGER, EXTERNAL :: IDX

      IF (IDX(NPLOT) > 0) THEN
        !Some plots insert their own "PLOT:"-line to allows KASDEF commands
        ! before the coordinate box. Therefore only write "PLOT:"-line if not empty!
        WRITE (KANAL,1) NPLOT
      ENDIF
      WRITE (KANAL, '(A)') 'KASDEF FONT=HELVET'
      WRITE (KANAL,2) NHEAD
      WRITE (KANAL,3) NX
      WRITE (KANAL,4) NY
    1 FORMAT (' PLOT   :',A)
    2 FORMAT (' HEADER :',A)
    3 FORMAT (' X-ACHSE:',A)
    4 FORMAT (' Y-ACHSE:',A)

C***  Special branch activating Lars' AUTO axes: set XMIN = XMAX and/or YMIN = YMAX!
      IF (B2 == B3 .AND. C2 == C3) THEN
         WRITE (KANAL, 5)
    5   FORMAT (5X,
     > 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     >  / ,' X: AUTO',/,' Y: ')
      ELSEIF (B2 == B3) THEN
         WRITE (KANAL, 6) C1,C2,C3,C4,C5,C6
    6   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: AUTOX',/,' Y: ',6(G12.6,1X))
      ELSE IF (C2 .EQ. C3) THEN
        WRITE (KANAL,7) B1,B2,B3,B4,B5,B6
    7   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(G12.6,1X),/,' Y: AUTO')
      ELSE 
        WRITE (KANAL,8) B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6
    8   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(G12.6,1X),/,' Y: ',6(G12.6,1X))
      ENDIF

C***  No Dataset opened for Zero Length
      IF (N .GT. 0) THEN
        CALL PLOTTAB (KANAL,X,Y,N,ISYMBOL)
      ELSE
        WRITE (KANAL, '(A)') 'END'
      ENDIF

      RETURN
      END
      SUBROUTINE PLOTANFS (KANAL,NPLOT,NHEAD,NX,NY,
     $ B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6,
     $ X,Y,N,STYLE)
C***********************************************************************
C***  DIESE ROUTINE BEREITET EIN NEUES PLOTFILE ZUR UEBERTRAGUNG VOR
C***  DIE ANGEGEBENEN WERTE BESCHREIBEN DEN PLOTKASTEN
C***********************************************************************
      CHARACTER*(*) NPLOT,NHEAD,NX,NY,STYLE

      WRITE (KANAL,1) NPLOT
      WRITE (KANAL, '(A)') 'KASDEF FONT=HELVET'
      IF (N .GT. 100000) WRITE (KANAL,'(A,I7)') 
     >    'KASDEF SET_NDATMAX ', N
      WRITE (KANAL,2) NHEAD
      WRITE(KANAL,3) NX
      WRITE (KANAL,4) NY
    1 FORMAT (' PLOT   :',A)
    2 FORMAT (' HEADER :',A)
    3 FORMAT (' X-ACHSE:',A)
    4 FORMAT (' Y-ACHSE:',A)

C***  Special branch activating Lars' AUTO axes: set XMIN = XMAX!
      IF (B2 .EQ. B3) THEN
         WRITE (KANAL, 6)
    6   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: AUTO',/,' Y: ')
      ELSE IF (C2 .EQ. C3) THEN
        WRITE (KANAL,7) B1,B2,B3,B4,B5,B6
    7   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(G13.6,1X),/,' Y: AUTO')
      ELSE
        WRITE (KANAL,5) B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6
    5   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(1PG13.6,1X),/,' Y: ',6(G13.6,1X))
      ENDIF

      CALL PLOTTABS (KANAL,X,Y,N,STYLE)

      RETURN
      END
      SUBROUTINE PLOTCON (KANAL,X,Y,N,ISYMBOL)
C***********************************************************************
C***  DIESE ROUTINE DIENST ZUM EINTRAGEN EINER 2. ODER 3. FUNKTION IN DEN
C***  BEGONNENEN PLOT ( PLOTCON ) , BZW. AB ENTRY PLOTTAB ZUM EINTRAGEN
C***  DER ERSTEN FUNKTION
C***********************************************************************
      DIMENSION X(N),Y(N)

C*** Backspace last Line (current ENDE)
      BACKSPACE KANAL

      ENTRY PLOTTAB (KANAL,X,Y,N,ISYMBOL)

C***  CHECK FOR SMALL X- AND Y-VALUES
      DO I=1, N
        IF (ABS(X(I)) .LT. 1.E-30) X(I) = 0.
        IF (ABS(Y(I)) .LT. 1.E-30) Y(I) = 0.
      ENDDO

      WRITE (KANAL,1) N,ISYMBOL
    1 FORMAT (' N=',I5,'   PLOTSYMBOL=',I3)

      DO 3 I=1,N,5
      J=MIN0( 5,N-I+1)
      WRITE (KANAL,2) (X(I+L-1),L=1,J)
      WRITE (KANAL,2) (Y(I+L-1),L=1,J)
    2 FORMAT (1X, 5(G15.8,1X))
    3 CONTINUE

      WRITE (KANAL,4)
    4 FORMAT (' ENDE')

      RETURN
      END
      SUBROUTINE PLOTCONS (KANAL,X,Y,N,STYLE)
C***********************************************************************
C***  DIESE ROUTINE DIENST ZUM EINTRAGEN EINER 2. ODER 3. FUNKTION IN DEN
C***  BEGONNENEN PLOT ( PLOTCON ) , BZW. AB ENTRY PLOTTAB ZUM EINTRAGEN
C***  DER ERSTEN FUNKTION
C***********************************************************************
      DIMENSION X(N),Y(N)
      CHARACTER STYLE*(*)

C*** Backspace last Line (current ENDE)
      BACKSPACE KANAL

      ENTRY PLOTTABS (KANAL,X,Y,N,STYLE)

C***  CHECK FOR SMALL X- AND Y-VALUES
      DO I=1, N
        IF (ABS(X(I)) .LT. 1.E-30) X(I) = 0.
        IF (ABS(Y(I)) .LT. 1.E-30) Y(I) = 0.
      ENDDO

      WRITE (KANAL,1) N, STYLE
    1 FORMAT (' N=', I7, 1X, A)

      DO 3 I=1,N,5
      J=MIN0( 5,N-I+1)
      WRITE (KANAL,2) (X(I+L-1),L=1,J)
      WRITE (KANAL,2) (Y(I+L-1),L=1,J)
    2 FORMAT (1X, 5(G15.8,1X))
    3 CONTINUE

      WRITE (KANAL,4)
    4 FORMAT (' ENDE')

      RETURN
      END
      SUBROUTINE PLOTV (ND ,RADIUS, VELO, MODHEAD, JOBNUM)
C******************************************************************************
C***  DIRECT PLOT TRANSFER OF THE VELOCITY STRATIFICATION
C***              V(R) VERSUS LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, DIMENSION(ND) :: RADIUS, VELO

      !Velocity law common block
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE, 
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER :: KANAL, L
      REAL :: XMIN, XMAX, YMIN, YMAX, VCON,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
 
C***  INITIALIZATION
      KANAL=2
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('V-STRATIFICATION TO BE ROUTED')
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' WR VELOCITY STRATIFICATION V(R) VERSUS LOG(R/R*-1)'
      HEAD2      = MODHEAD(13:32)
      HEAD2(22:) = 'VELOCITY STRATIFICATION v(r)'

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX=2.9
      XMIN=-3.0
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
      YMAX = ANINT(VFINAL/100.+1.)
      YMIN = .0
      YSCALE = 15./(YMAX-YMIN)
      YTICK=2.5
      YABST=5.

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.)
        Y(L) = VELO(L) / 100.
      ENDDO
 
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'
C***  BORDER BETWEEN INNER/OUTER VELOCITY LAW
      IF ((RCON > RADIUS(ND)) .AND. (RCON <= RADIUS(1))) THEN
        CALL SPLINPOX(VCON, RCON, VELO, RADIUS, ND)
        WRITE (KANAL, '(A,2F10.3,A)')
     >    '\SYM ', LOG10(RCON-1.), VCON / 100., ' 0. 0. 0.3 4'
      ENDIF
      WRITE (KANAL, '(A,F10.3)') '*RCON=', RCON
      DO L=1, ND-1
        IF (MOD(L, 10) == 0) THEN
          WRITE (KANAL,41) X(L), X(L), X(L), L
   41     FORMAT ('\LINREL ', F7.3, ' YMAX 0. -0.5', /,
     >            '\LINREL ', F7.3, ' YMIN 0.  0.5', /,
     >            '\LUN    ', F7.3, ' YMIN -0.2 0.7 0.3 ', I3)
        ENDIF
      ENDDO

      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,'\CENTER\log (r/R&T*&M - 1)'
     >        ,'\CENTER\v(r) / [100 km/s]'
     >        ,XSCALE,0.,0.,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
      CALL PLOTCONS (KANAL,X,Y,ND-1,'SYMBOL=8 SIZE=0.1') 
 
      RETURN
      END
      SUBROUTINE PLOTVGRAD (ND ,RADIUS, VELO, MODHEAD, JOBNUM, 
     >                      RMIN, RMAX)
C******************************************************************************
C***  PLOT OF THE VELOCITY GRADIENT STRATIFICATION
C***              DV/DR(R) VERSUS LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, INTENT(IN) :: RMIN, RMAX
      REAL, DIMENSION(ND) :: RADIUS, VELO, VELO2, GRAD1, GRAD2

      !Velocity law common block
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON / VELPAR / VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      REAL, EXTERNAL :: WRVEL

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER :: hPLOT, L, NDBG
      REAL :: XMIN, XMAX, YMIN, YMAX, Vdummy, RP, RM,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
 
      !calculate Gradients
      DO L=1, ND
        IF (VELO(L) > 0.) THEN
          CALL SPLINPOX(Vdummy,RADIUS(L),VELO,RADIUS,ND,DFDX=GRAD1(L))
        ELSE 
          GRAD1(L) = 0.
        ENDIF
      ENDDO

      !GRAD2 = Beta law gradient
      DO L=1, ND
        IF (RADIUS(L) + VPAR2 <= 1.) THEN
          EXIT 
        ENDIF
        RP = VPAR2 + RADIUS(L)
        GRAD2(L) = VPAR1*BETA*(1.-1./RP)**(BETA-1.) /(RP*RP) 
        NDBG = L
      ENDDO

C***  INITIALIZATION
      hPLOT=2
      OPEN (hPLOT, FILE='vgrad.plot', STATUS='UNKNOWN')
C      CALL JSYMSET ('G2','TRANSFER')
C      CALL REMARK ('V-STRATIFICATION TO BE ROUTED')
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' WR VELOCITY GRADIENT STRATIFICATION DV/DR(R)' //
     >          ' VERSUS LOG(R/R*-1)'
      IF (MODHEAD(:5) /= 'DEBUG') THEN
        HEAD2      = MODHEAD(13:32)
        HEAD2(22:) = 'VELOCITY GRADIENT STRATIFICATION dv/dr(r)'
      ELSE
        HEAD2 = 'VELOCITY GRADIENT STRATIFICATION dv/dr(r)'
      ENDIF

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX= LOG10(RADIUS(1)-1.)
      XMIN= LOG10( (RADIUS(ND-1)-1.) * 0.9 ) 
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY GRADIENT RANGE
C*** ( 1. + (BETA-1.)/2. - VPAR2 ) * 1.5
      YMAX = 1.5 * MAXVAL(GRAD2(1:NDBG))   !150% of maximum beta law gradient
      YMIN = -9.0
      YSCALE = 15./(YMAX-YMIN)
      YTICK= INT( (YMAX-YMIN) / 20. )
      YABST=5. * YTICK

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.)
C        Y(L) = VELO(L) / 100.
        Y(L) = GRAD1(L)
      ENDDO
 
      WRITE (hPLOT, '(A,A)') 'PLOT: ', HEAD1
      WRITE (hPLOT, '(A)') '\INBOX'

C***  BORDER BETWEEN INNER/OUTER VELOCITY LAW
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=9'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RCON-1.),
     >    ' YMIN ',LOG10(RCON-1.),' YMAX 0. 0. SYMBOL=10'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RCON-1.),' YMAX 0. -0.2 0.2 -90  R&Tcon&M'
      ENDIF
      IF (RMIN > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=5'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RMIN-1.),
     >    ' YMIN ',LOG10(RMIN-1.),' YMAX 0. 0. SYMBOL=9'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RMIN-1.),' YMAX 0. -0.2 0.2 -90  min.'
      ENDIF
      IF (RMAX > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=6'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RMAX-1.),
     >    ' YMIN ',LOG10(RMAX-1.),' YMAX 0. 0. SYMBOL=9'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RMAX-1.),' YMAX 0. -0.2 0.2 -90  max.'
      ENDIF
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'
      WRITE (hPLOT, '(A)') '\COLOR=1'


      CALL PLOTANF (hPLOT,HEAD1,HEAD2
     >        ,'\CENTER\log (r/R&T*&M - 1)'
     >        ,'\CENTER\dv/dr(r) / [km/s/Rstar]'
     >        ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
 
      NDBG = MIN(NDBG, ND-1)
      CALL PLOTCONS (hPLOT,X,GRAD2,NDBG,'COLOR=2 SYMBOL=9 SIZE=0.1') 

  
      CLOSE(hPLOT)

      RETURN
      END
      SUBROUTINE PREP_DRLINES (DRLINES_CARD, NAUTO, KRUDAUT, EAUTO)
C******************************************************************
C***  Prepares which of the DRTRANSITs (i.e. stabilizing transitions
C***  from auto-ionisation levels) should be treated as: 
C***    - "normal" lines:     KRUDAUT = 0 
C***    - "rudimental" lines: KRUDAUT = 1 
C***  This depends on the CARDS line DRLINES
C***  called from: WRSTART, COLI, STEAL
C******************************************************************
      CHARACTER*120 DRLINES_CARD
      CHARACTER*4 ACTPAR
      DIMENSION KRUDAUT(NAUTO), EAUTO(NAUTO)

C***  Default: no DRLINES_CARD was encountered
      IF (DRLINES_CARD .EQ. '') THEN
C***     Only lines that are truely auto-ionizing are set rudimantel         
         DO IND=1, NAUTO
            IF (EAUTO(IND) < .0) THEN
               KRUDAUT(IND) = 0
            ELSE
               KRUDAUT(IND) = 1
            ENDIF
         ENDDO
      ELSE
         CALL SARGC (DRLINES_CARD, NPAR)
         IF (NPAR .LT. 2) GOTO 97
         CALL SARGV (DRLINES_CARD, 2, ACTPAR)

         IF (ACTPAR .EQ. 'ALL') THEN
C***                      ===
C***        All lines set non-rudimental:
            DO IND=1, NAUTO
               KRUDAUT(IND) = 0
            ENDDO
          
         ELSEIF (ACTPAR .EQ. 'NONE') THEN
C***                          ====
C***        All lines set rudimental:
         write (0,*) 'in PREP_DRLINES: parameter NONE detected'
            DO IND=1, NAUTO
               KRUDAUT(IND) = 1
            ENDDO

         ELSE
            GOTO 98
         ENDIF 
      ENDIF   

      RETURN 

C***  ERROR exits
   97 WRITE (0,*) '*** ERROR: DRLEVELS needs one parameter:' 
      WRITE (0,*) '*** allowed are: ALL or NONE'  
      GOTO 99

   98 WRITE (0,*) '*** ERROR: DRLEVELS has invalid parameter:'  
      WRITE (0,*) '*** allowed are: ALL or NONE'  
      GOTO 99

   99 WRITE (0,*) '*** the error occured on the following CARDS line:'
      WRITE (0,*) DRLINES_CARD(:IDX(DRLINES_CARD)) 
      STOP 'Fatal ERROR detected by subroutine PREP_DRLINES'

      END
      SUBROUTINE PREP_GAMMARAD (bOLDRAD, bFULLHYDROSTAT, GEDD, 
     >      IGEddFix, GAMMARAD, NDold, ARAD, GLOG, RSTAR,
     >      RADIUSold, RSTARold, XMGold, TAURold, RCONold, GEFFLOG,
     >      STAPEL, ATMEAN, XLOGL, XMSTAR, RADGAMMASTART, 
     >      GEDDRAD, bOldStratification, bGAMMARADMEAN, bSaveGEFF)
C***********************************************************************
C***  This subroutine prepares GAMMARAD and related quentities for
C***  the hydrostatic equation
C***  Called from: WRSTART
C***********************************************************************

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NDold, IGEddFix
      REAL, INTENT(IN) :: XMGold, XLOGL, STAPEL, ATMEAN,
     >                    RSTARold, RCONold, RADGAMMASTART
      
      REAL, DIMENSION(NDold), INTENT(IN) :: RADIUSold, TAURold, ARAD
      REAL, DIMENSION(NDold), INTENT(OUT) :: GAMMARAD
      REAL, DIMENSION(NDold) :: ARADL 

      REAL, INTENT(INOUT) :: GEDD, GEFFLOG, GLOG, XMSTAR
      REAL, INTENT(OUT) :: GEDDRAD

      LOGICAL, INTENT(IN) :: bOLDRAD, bFULLHYDROSTAT, bSaveGEFF,
     >                       bOldStratification, bGAMMARADMEAN

      INTEGER :: L
      REAL :: GEDDRADMEAN, TAUNORM, WTAU, RSTAR

      !Physical constants
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      GEDDRAD = 0.      
      
      IF (bOLDRAD .AND. bFULLHYDROSTAT) THEN
        !calculate effective g via ARAD from OLD MODEL

        IF (GEDD > .0 .AND. IGEddFix == 2) THEN
          !Gamma value has been fixed via CARDS 
          DO L=1, NDold
            GAMMARAD(L) = GEDD
          ENDDO

        ELSE

C***      Interpolate ARAD from interstices to ARADL on full depth points  
          DO L=1, NDold
            IF (L==1) THEN
              ARADL(L) = ARAD(1)
            ELSEIF (L==NDold) THEN
              ARADL(L) = ARAD(NDold-1)
            ELSE
              ARADL(L) = 0.5 * (ARAD(L-1) + ARAD(L))
            ENDIF
            GAMMARAD(L) = ARADL(L)*(RADIUSold(L)*RSTARold)**2/XMGold
C***        restrict GAMMARAD 
            GAMMARAD(L) = MIN(GAMMARAD(L), 0.9)
          ENDDO

          !Calculate mean value, needed for g-g_eff-Relation
          ! if RADGAMMASTART: OLD is used (cannot be done in DECSTAR)
          GEDDRADMEAN = GAMMARAD(NDold) * EXP(-TAURold(NDold))
          TAUNORM = 0.
          DO L=NDold-1, 1, -1
            IF (RADIUSold(L) > RCONold) EXIT  
            WTAU = (TAURold(L+1) - TAURold(L)) * EXP(-TAURold(L))
            GEDDRADMEAN = GEDDRADMEAN + GAMMARAD(L) * WTAU
            TAUNORM = TAUNORM + WTAU
          ENDDO
          IF (TAUNORM > 0.) THEN
            GEDDRADMEAN = GEDDRADMEAN / TAUNORM
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDDRADMEAN))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = LOG10((10.**GLOG)*(1.-GEDDRADMEAN))
          ENDIF
          
          IF (bGAMMARADMEAN) THEN
            !Use mean value of the inner part instead of individual values (CARDS option)
            DO L=1, NDold
              GAMMARAD(L) = GEDDRADMEAN
            ENDDO
          ENDIF
        ENDIF

      ELSEIF (GEDD > 0. .OR. RADGAMMASTART >= 0.) THEN
C***    for new models only g_THOMSON is available, 
C***    unless RADGAMMA-START is specified
        IF (IGEddFix == 0) THEN        
          IF (RADGAMMASTART > 0) THEN
            GEDD = RADGAMMASTART
          ELSE
            GEDD = 10.**(-4.51) * (STAPEL/ATMEAN) *(10.**XLOGL) /XMSTAR
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDD))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
          ENDIF
        ENDIF
        IF (RADGAMMASTART >= 0.) THEN
          DO L=1, NDold
            GAMMARAD(L) = RADGAMMASTART
          ENDDO
        ENDIF
      ENDIF

      IF (bFULLHYDROSTAT) THEN
        IF (bOLDRAD) THEN
          GEDDRAD = GEDDRADMEAN
        ELSEIF (RADGAMMASTART >= 0.) THEN
          GEDDRAD = RADGAMMASTART
        ELSE
          GEDDRAD = GEDD      !poor start with Thomson only
        ENDIF        
      ENDIF
      
C***  Only if not copying the old stratification
      IF (.NOT. bOldStratification) THEN
        IF (bFULLHYDROSTAT) THEN
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Full Gamma (ND) = ', GEDDRAD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG

        ELSE
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Eddington Gamma = ', GEDD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG
        ENDIF
      ENDIF
      
      RETURN 
      END
      SUBROUTINE PRICOMP (NDIM,EINST,N,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     $                   STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     $                   INDLOW,INDNUP,NAUTO,LOWAUTO,
     $                   EAUTO,KONTNUP,KONTLOW,LASTKON,XMASS,KRUDAUT)
C***********************************************************************
C***  PRINTOUT OF THE CHEMICAL COMPOSITION OF THE WR MODEL ATMOSPHERE
C***********************************************************************
 
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION NCHARG(N),NOM(N)
      DIMENSION ABXYZ(NATOM)
      DIMENSION ATMASS(NATOM),STAGE(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION LOWAUTO(NAUTO),EAUTO(NAUTO),KRUDAUT(NAUTO)
      CHARACTER*80 KARTE
      CHARACTER*10 ELEMENT(NATOM)

      PARAMETER ( MAXION = 24 )
      CHARACTER*4 ROMION(0:MAXION)
      CHARACTER*2 SYMBOL(NATOM)


      DATA (ROMION(I),I=0,MAXION) /'   I','  II',' III','  IV','   V',
     >                             '  VI',' VII','VIII','  IX','   X',
     >                             '  XI',' XII','XIII',' XIV','  XV',
     >                             ' XVI','XVII',' 18 ',' XIX','  XX',
     >                             ' XXI','XXII',' 23 ','XXIV',' XXV'/

 
      PRINT 1
    1 FORMAT (1H1,//,
     $        10X,'C H E M I C A L  C O M P O S I T I O N',
     $         4X,'(FILE: "DATOM")',/,
     $        10X,38('='),/)

C***  DECODING FIRST INPUT CARDS (COMMENTS) FROM TAPE 4 = DATOM
      OPEN (4, FILE='DATOM', STATUS='UNKNOWN')
    9 READ (4,6, END=99) KARTE
    6 FORMAT(A)
      IF (KARTE(:1) .EQ. '*') THEN
C***     EXPLICIT END ("*END") FOR INPUT CARDS IS GIVEN:
         IF (KARTE(2:4) .EQ. 'END') GOTO 99
         PRINT 2, KARTE
    2    FORMAT (1X,A)
         GOTO 9
      ENDIF
   99 CLOSE (4)
      PRINT 22
   22 FORMAT (//)

C***  CALCULATION OF THE MEAN ATOMIC WEIGHT "ATMEAN";
C***  SEARCH FOR THE INDEX OF MAIN ELEMENT "HELIUM"
      ATMEAN=0.0
      NAHE=0
      DO 10 NA=1, NATOM
      ATMEAN=ATMEAN+ABXYZ(NA)*ATMASS(NA)
      IF (SYMBOL(NA) .EQ. 'HE') NAHE=NA
   10 CONTINUE
      IF (NAHE .GT. 0) THEN
         PRINT 11
   11    FORMAT
     >   (1X,' REL. ABUNDANCES:',12X,' BY NUMBER ',5X,'  BY MASS  ',
     >    11X,' REL. TO HELIUM:',12X,' BY NUMBER ',5X,'  BY MASS  ')
      ELSE
         PRINT 12
   12    FORMAT
     >   (1X,' REL. ABUNDANCES:',12X,' BY NUMBER ',5X,'  BY MASS  ')
      ENDIF

      DO 19 NA=1,NATOM
      ABNA=ABXYZ(NA)
      ATMASNA=ATMASS(NA)
      FRACM=ABNA*ATMASNA/ATMEAN
      IF (NAHE .GT. 0) THEN
         ABHE=ABXYZ(NAHE)
         PRINT 13, SYMBOL(NA),ABNA,FRACM,
     >             SYMBOL(NA),ABNA/ABHE,ABNA*ATMASNA/(ABHE*ATMASS(NAHE))
      ELSE
         PRINT 13, SYMBOL(NA),ABNA,FRACM
      ENDIF
   13 FORMAT (23X,A2,5X,1PE11.4,5X,E11.4,32X,A2,5X,E11.4,5X,E11.4)
   19 CONTINUE
 
      NION=1
      DO 20 NLEV=2,N
      IF ((NOM(NLEV) .NE. NOM(NLEV-1)) .OR.
     $    (NCHARG(NLEV) .NE. NCHARG(NLEV-1))) NION=NION+1
   20 CONTINUE
      NDR=0
      NLTE=0
      DO 26 IND=1,NAUTO
      IF (EAUTO(IND) .GT. 0.) NDR=NDR+1
      IF (EAUTO(IND) .LT. 0.) NLTE=NLTE+1
   26 CONTINUE
      IF ((NDR+NLTE) .NE. NAUTO) STOP 'NAUTO'
      PRINT 21, NATOM,NION,N,LASTKON,LASTIND,NDR,NLTE
   21 FORMAT (///,1X,' STATISTICS:',5X,I3,'  ELEMENTS',/,18X,I3,'  IONS'
     $        ,/,17X,I4,'  LEVELS'
     $        ,/,17X,I4,'  CONTINUUM TRANSITIONS'
     $        ,/,16X,I5,'  LINE TRANSITIONS'
     $        ,/,17X,I4,'  STABIL. TRANSITIONS FROM AUTOIONIZING '
     $                 ,'STATES'
     $        ,/,17X,I4,'  TRANSITIONS FROM LTE LEVELS',///)
 
      PRINT 31
   31 FORMAT (1X,' INDEX',3X,'ELEMENT',6X,'ATOMIC MASS',3X,'IONS',3X,
     $        'MAIN ION',/,1X,52('-'))
      DO 39 NA=1,NATOM
      NAION=1
      DO 38 NLEV=NFIRST(NA)+1,NLAST(NA)
      IF (NCHARG(NLEV) .NE. NCHARG(NLEV-1)) NAION=NAION+1
   38 CONTINUE
      ISTAGE = INT(STAGE(NA))-1
      IF (ISTAGE .GT. MAXION) THEN
         WRITE (0,*) '*** MAXION too small ***************'
         STOP        '*** FATAL ERROR IN SUBR. PRICOMP ***'
      ENDIF
      PRINT 32,NA,ELEMENT(NA),ATMASS(NA),NAION,SYMBOL(NA),
     $         ROMION(ISTAGE)
   32 FORMAT (3X,I2,5X,A10,5X,F6.2,7X,I2,5X,A2,1X,A4)
   39 CONTINUE

      PRINT 50, ATMEAN, XMASS
   50 FORMAT (///, '  MEAN ATOMIC MASS:', F6.2, 5X, 
     >        'MEAN PARTICLE MASS:', F6.2)

 
      PRINT 41
   41 FORMAT (///,1X,' ELEMENT',3X,'ION',3X,'CHARGE',3X,'LEVELS',3X,
     $        'CONTINUA',3X,
     $        'LINES',1X,'(RUD.)',3X,
     $        'LINES(LTE)',1X,'(RUD.)',3X,
     $        'STAB.LINES',1X,'(RUD.)',
     $        /,1X,99('-'))
      DO 49 NA=1,NATOM
      N1=NFIRST(NA)
   44 NCH1=NCHARG(N1)
      IF (NCH1 .EQ. NCHARG(NLAST(NA))) THEN
          NION=NLAST(NA)-N1+1
      ELSE
          NION= ISRCHNE(NLAST(NA)-N1+1,NCHARG(N1+1),1,NCH1)
      ENDIF
      NKONT=0
      NLINE=0
      NRUD1=0
      NLTE=0
      NRUD2=0
      NDR=0
      NRUD3=0
      DO 43 KON=1,LASTKON
      LOW=KONTLOW(KON)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) 
     $    NKONT=NKONT+1
   43 CONTINUE
      DO 45 IND=1,LASTIND
      LOW=INDLOW(IND)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) THEN
          NLINE=NLINE+1
          IF (EINST(LOW,INDNUP(IND)) .EQ. -2.) NRUD1=NRUD1+1
      ENDIF
   45 CONTINUE
      DO 46 IND=1,NAUTO
      LOW=LOWAUTO(IND)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) THEN
         IF (EAUTO(IND) .LT. 0.) THEN
            NLTE=NLTE+1
            IF (KRUDAUT(IND) .EQ. 1) NRUD2=NRUD2+1
         ENDIF
         IF (EAUTO(IND) .GT. 0.) THEN
            NDR=NDR+1
            IF (KRUDAUT(IND) .EQ. 1) NRUD3=NRUD3+1
         ENDIF
      ENDIF
   46 CONTINUE
      PRINT 42, SYMBOL(NA),ROMION(NCH1),NCH1,NION,NKONT,NLINE,NRUD1,
     >          NLTE,NRUD2,NDR,NRUD3
   42 FORMAT (4X,A2,6X,A4,4X,I2,7X,I2,8X,I2,6X,I4,3X,I4,8X,I3,6X,I3,
     >        7X,I4,6X,I3)
      N1=N1+NION
      IF (N1 .LE. NLAST(NA)) GOTO 44
   49 CONTINUE

      RETURN
      END
      SUBROUTINE PRIMOD (ND,RADIUS,INCRIT,ENTOT,T,VELO,GRADI,NP,OLDTEMP,
     >          MODHEAD,JOBNUM,MODOLD,JOBNOLD,TTABLE,TAUROSS,R23,
     >          TEFFOLD, THIN, ITTAU, MAXITTAU, RCON, BTWOT, MODOLD2, 
     >          JOBNOLD2, TFAC, BETA, VPAR1, VPAR2, DENSCON, BETA2,
     >          BETA2FRACTION, HSCALE, bNoDetails,bStartCall, VTURB)
C*******************************************************************************
C***  PRINTOUT OF THE DEPTH-DEPENDENT MODEL SPECIFICATIONS
C***  CALLED FROM WRSTART, STEAL
C*******************************************************************************
 
      IMPLICIT NONE
       
      REAL, DIMENSION(ND), INTENT(IN) :: 
     >      RADIUS, ENTOT, T, VELO, GRADI, TAUROSS, DENSCON
      REAL, INTENT(IN) :: 
     >      RCON, BETA, VPAR1, VPAR2, BETA2, BETA2FRACTION, HSCALE,
     >      R23, TEFFOLD, TFAC, VTURB
      INTEGER, INTENT(IN) ::
     >      ND, NP, JOBNUM, JOBNOLD, JOBNOLD2, ITTAU, MAXITTAU
      LOGICAL, INTENT(IN) ::
     >      OLDTEMP,
     >      TTABLE,
     >      BTWOT,
     >      THIN,         !THIN-Wind option set
     >      bNoDetails,   !Do not print table with values for all depth points
     >      bStartCall    !if false, informations only known by WRSTART are skipped
      CHARACTER(8), DIMENSION(ND), INTENT(IN) :: 
     >      INCRIT
      CHARACTER(100), INTENT(IN) :: 
     >      MODHEAD, MODOLD, MODOLD2

      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

      REAL :: RL1, RL2, DENS, TEFF, TMIN, TMODIFY
      INTEGER :: L
      LOGICAL :: SPHERIC
      CHARACTER(40) :: MESSAGE
     
      IF (MODHEAD /= ' ') THEN
        PRINT 1,MODHEAD,JOBNUM
    1   FORMAT (1X,  A  ,20X,'JOB NO.',I7,
     $       ///,20X,'M O D E L   S P E C I F I C A T I O N',/,20X,
     $  38('='))
      WRITE (*,*)
      ENDIF
      
      IF (.NOT. bNoDetails) THEN    !details option
        WRITE (*,*)
        WRITE (*,'(20X,A)') '   DEPTH-DEPENDENT QUANTITIES'
        WRITE (*,*)
        WRITE (*,'(3A)')
     >   ' DEPTH     R-1     RMAX-R    CRITERION     VELOCITY        ',
     >   'GRADIENT      LOG NUMB. DENS.   EL. TEMPERATURE ',
     >   'TAU-ROSS.  CLUMPING'
        WRITE (*,'(3A)')
     >   ' INDEX                                      (KM/S)     (KM/',
     >   'S PER RSTAR)  (ATOMS PER CM+3)        (KELVIN)   ',
     >   '(IN LTE)   DENSCON'
        WRITE (*,*)

C***    LOOP OVER ALL DEPTH POINTS **************************************
        DO L=1,ND
          RL1=RADIUS(L)-1.
          RL2=RADIUS(1)-RADIUS(L)
          DENS=ALOG10(ENTOT(L))
          PRINT 3, L, RL1, RL2, INCRIT(L), VELO(L), GRADI(L), DENS,
     >         T(L), TAUROSS(L), DENSCON(L)
        ENDDO
    3   FORMAT(I3, 5X, G10.3, G10.3, 3X, A8, 2X, G11.4, F14.2, F19.3,
     >       F19.0, F12.3, F10.3)
C**********************************************************************
      ENDIF !end of details option

      IF (bStartCall) PRINT 4,NP
    4 FORMAT (/,' NUMBER OF IMPACT-PARAMETER POINTS  NP =',I3)

      PRINT 9, R23
    9 FORMAT (' RADIUS WHERE TAU-ROSSELAND = 2/3: R23 =',F7.3)


C***  MESSAGES CONCERNING THE TEMPERATURE STRUCTURE
      IF (TTABLE) PRINT *,'TEMPERATURE STRUCTURE FROM TABLE INPUT'

      IF (bStartCall .AND. SPHERIC) PRINT 7,TEFF
    7 FORMAT (' TEMPERATURE STRUCTURE AS FOR A SPHERICAL, ',
     $   'GREY LTE ATMOSPHERE (APPROXIMATELY) WITH TEFF=',F7.0)

      IF (OLDTEMP) THEN
        PRINT 11, MODOLD, JOBNOLD
   11   FORMAT (' TEMPERATURE STRUCTURE TAKEN FROM THE FOLLOWING ',
     >     'OLD MODEL:',/,1X,A,5X,'AFTER JOB NO.',I3)
        IF (BTWOT) THEN
          IF (TFAC .GE. 0. .AND. TFAC .LE. 1.) THEN
            IF (TFAC .EQ. 0.5) THEN 
              MESSAGE = ' INTERPOLATION'
            ELSE IF (TFAC .LT. 0.5) THEN 
              MESSAGE = ' INTERPOLATION, MODEL 1 STRESSED'
            ELSE
              MESSAGE = ' INTERPOLATION, MODEL 2 STRESSED'
            ENDIF
          ELSE
            IF (TFAC .LT. 0.) THEN
              MESSAGE = ' EXTRAPOLATION, MODEL 1 STRESSED'
            ELSE
              MESSAGE = ' EXTRAPOLATION, MODEL 2 STRESSED'
            ENDIF
          ENDIF
          WRITE (*,'(A,/,A,A,A,I3)') 
     >           ' SECOND MODEL : ',' ', MODOLD2, 
     >           'AFTER JOB NO.', JOBNOLD2
          WRITE (*,'(A,F4.1,A)')
     >           ' COMBINE FACTOR = ', TFAC, MESSAGE
        ENDIF
      ENDIF

      IF (OLDTEMP .AND.  TEFF .NE. TEFFOLD) PRINT 12, TEFF, TEFFOLD
   12 FORMAT (' ... SCALED WITH TEFF (NEW/OLD) =', F7.0, ' / ', F7.0)

      IF (bStartCall .AND.
     >      .NOT. TTABLE .AND. .NOT. SPHERIC .AND. .NOT. OLDTEMP)
     $   PRINT 8,TEFF
    8    FORMAT (' TEMPERATURE STRUCTURE AS FOR A ',
     $           'PLANE-PARALLEL, GREY LTE ATMOSPHERE WITH TEFF=',F7.0)
          
      IF (TMODIFY .NE. .0) PRINT 5,TMODIFY
    5 FORMAT (' TEMPERATURE STRATIFICATION MODIFIED: TMODIFY=',F6.3)
 
      IF (TMIN .GT. .0) PRINT 10, TMIN
   10 FORMAT (' MINIMUM TEMPERATURE SPECIFIED: TMIN=',F7.0)



C***  MESSAGES CONCERNING THE VELOCITY FIELD 
      IF (RCON == 1.) THEN
         PRINT *, 'VELOCITY FIELD: NO HYDROSTATIC DOMAIN ENCOUNTERED'
      ELSEIF (RCON > RADIUS(1) .AND. bStartCall) THEN
         PRINT *, 'VELOCITY FIELD: OLD STRATIFICATION USED'
      ELSE 
         PRINT '(2A,G10.3,A,1PG12.3,A)', 
     >         ' VELOCITY FIELD: HYDROSTATIC DOMAIN ',
     >         'BELOW RADIUS RCON = 1 + ', RCON-1. , '  (HSCALE = ',
     >         HSCALE,')'
         IF (THIN) PRINT *,'HYDROSTATIC EQ. WITH TEMPERATURE',
     >                     '-DEPENDENT SCALE HEIGHT IS USED'
         IF (VTURB .GT. .0) PRINT '(A, F6.1)',
     >      ' TURBULENCE PRESSURE TAKEN INTO ACCOUNT, VTURB=', VTURB
         ENDIF
      IF (ABS(BETA) .GT. .0) PRINT '(A,F4.1,A,3PG12.6,1(A,1PG12.3))',
     >  ' OUTER PART: BETA=', ABS(BETA), ' LAW, VPAR1=', VPAR1
     >  ,' , VPAR2=', VPAR2
      IF (BETA2FRACTION .GT. .0) PRINT '(13X,A,F4.2,A,F5.2)',
     >  '2-BETA-LAW:  BETA2FRACTION=', BETA2FRACTION,
     >    '  BETA2=', BETA2

      IF (ITTAU .GT. 1) THEN
       IF (ITTAU .LT. MAXITTAU) THEN
        PRINT '(A,I3,A)', 
     >    ' VMIN AUTOMATICALLY ADJUSTED TO SPECIFIED TAUMAX (', 
     >    ITTAU, ' ITERATIONS)'
        ELSE                  
        PRINT '(A,I3,A)', 
     >   '*** WARNING: AUTOMATICAL ADJUSTMENT OF VMIN NOT CONVERGED IN',
     >   ITTAU, ' ITERATIONS ! ***'
        ENDIF
       ENDIF
 
      PRINT *, ' '

      RETURN
      END
      SUBROUTINE PRIPARAM (MODHEAD, TEFF, RSTAR, XMDOT, XLOGL, RTRANS,
     >           VFINAL, VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG, GEDD, 
     >           GEddFix, RMAX, XMSTAR, WRTYPE, MASSORIGIN, LRTinput, 
     >           ND, MLRELATION, VTURB, MDOTINPUT)
C*******************************************************************************
C***  PRINTOUT OF THE STELLAR PARAMETERS
C***  CALLED FROM WRSTART
C***   Note: Some of these values can chance significantly during the iteration
C***    process, use results of STEAL->PRINTMODELSUMMARY to obtain final values
C*******************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, MASSORIGIN, GEddFix, LRTinput
      INTEGER, INTENT(IN) :: MDOTINPUT
      REAL, DIMENSION(ND), INTENT(IN) :: DENSCON, FILLFAC
      REAL, INTENT(IN) :: TEFF, RSTAR, XMDOT, GEDD,
     >                    VFINAL, VDOP, GLOG, GEFFLOG, RMAX, XMSTAR
      REAL, INTENT(INOUT) :: XLOGL, RTRANS
      CHARACTER(100), INTENT(IN) :: MODHEAD
      CHARACTER(2), INTENT(IN) :: WRTYPE
      CHARACTER(9) :: MLRELATION
      CHARACTER(30) :: FORMSTR

      REAL :: RSTARSU, VTURB

      REAL, PARAMETER :: RSUN = 6.96E10     !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: TEFFSUN = 5780.    !SOLAR EFFECTIVE TEMPERATURE
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


      RSTARSU = RSTAR/RSUN                  !STELLAR RADIUS IN SOLAR UNITS
 
      WRITE (hOUT,1) MODHEAD(:32), MODHEAD(33:)
    1 FORMAT (///,80('*'),/,'*',/,
     >        '* FUNDAMENTAL PARAMETERS',/,
     >        '* ======================',/,'*',/,
     >        '* ',A,/,'* ',A,/'*')

      SELECTCASE (LRTinput) 
        CASE (1)  !TEFF & RSTAR given => L calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM RSTAR AND TEFF)'
        CASE (2)  !TEFF & L given => RSTAR calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (CALCULATED FROM LUMINOSITY AND TEFF)'
        CASE (3)  !L & RSTAR given => TEFF calculated
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (CALCULATED FROM RSTAR AND LUMINOSITY)'
        CASE (4)  !MSTAR & GEDD & TEFF given => L & RSTAR calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM MSTAR AND EDDINGTON_GAMMA)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (CALCULATED FROM LUMINOSITY AND TEFF)'
        CASE (5)  !MSTAR & GEDD & RSTAR given => TEFF & L calculated
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM MSTAR AND EDDINGTON_GAMMA)'
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (CALCULATED FROM RSTAR AND LUMINOSITY)'
        CASE DEFAULT
          WRITE(hCPR,'(A,I5)') 'INVALID VALUE FOR LRTinput: ', LRTinput
          STOP 'FATAR ERROR detected in DECSTAR'
      ENDSELECT
            

C***  calculation of Rtrans if not done in DECSTAR
      IF (RTRANS .LE. .0)  
     >    RTRANS = RSTARSU * ( VFINAL / 2500. 
     >          / 10.**(XMDOT+0.5*ALOG10(DENSCON(1))+4.0))**(2./3.)
      IF (RTRANS > 999.) THEN
         FORMSTR = '(A,F7.0,A,F6.3,A)'
      ELSE
         FORMSTR = '(A,F7.3,A,F6.3,A)'

      ENDIF

      SELECTCASE (MDOTINPUT) 
        CASE (1)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (INPUT)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM MDOT)'

        CASE (2)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM RTRANS)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (INPUT)'

        CASE (3)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM MDTRANS)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM MDTRANS)'

        CASE (4)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM LOG Q)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM LOG Q)'

      ENDSELECT 
     

      WRITE (hOUT,'(A,I7,A)') '* VFINAL  =', NINT(VFINAL), ' KM/S'
      WRITE (hOUT,'(A,I7,A)') '* VDOP    =', NINT(VDOP)  , ' KM/S'


c***  tiefenabh. clumping...
C***  aufpassen, evtl. ueber den L-Bereich??
      IF (DENSCON(1) .NE. 1.) THEN 
         WRITE (hOUT,6) DENSCON(1), FILLFAC(1)
    6    FORMAT ('* DENSCON(1) =', F7.2, '    FILLFAC(1) =', F6.4)  
      ELSE 
         WRITE (hOUT,'(A)') '* MODEL WITHOUT CLUMPING'
      ENDIF

      IF (MASSORIGIN == 0) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]  CALCULATED FROM:'  
         WRITE (hOUT,'(A)') '*   MASS-LUMINOSITY ' //
     >                   'RELATION FOR TYPE = ' // WRTYPE
     >                   // ' from ' // MLRELATION
         WRITE (hOUT,'(A,F6.2,A,A,F5.2,A,A,F5.2)')
     >          '*   MASS =', XMSTAR, ' M_SUN',  
     >                ' -- CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
      ELSEIF (MASSORIGIN == 1) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* MSTAR   =', XMSTAR, ' M_SUN   (INPUT)'  
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F5.2, A,/,A,F5.2,A,A,F5.2)')
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >         '*   CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD 
         ELSE
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A,/,A,F5.2,A,A,F5.2,A)')
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)',
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' AND EDDINGTON_GAMMA = ', GEDD , ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A,/,A,F5.2,A,A,F5.2,A)')
     >         '* EDDINGTON_GAMMA = ', GEDD, '   (INPUT)',
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' AND LOG G_EFF = ', GEFFLOG, ' [CGS]'
           ENDIF
         ENDIF
      ELSEIF (MASSORIGIN == 2) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]   (INPUT)'  
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F6.2, A,/,A,F5.2,A,A,F5.2)') 
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >         '*   CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
         ELSE
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A,/,A,F6.2,A,A,F5.2,A)') 
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)',
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >                ' AND EDDINGTON_GAMMA = ', GEDD, ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A,/,A,F6.2,A,A,F5.2,A)') 
     >         '* EDDINGTON_GAMMA = ', GEDD, '   (INPUT)',
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >                ' AND LOG G_EFF = ', GEFFLOG,' [CGS]'
           ENDIF
         ENDIF
      ELSEIF (MASSORIGIN == 3) THEN
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_EFF =', GEFFLOG, ' [CGS]   (INPUT)'  
           WRITE (hOUT,'(A,F5.2, A, A, F5.2)') 
     >         '*   CALCULATED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
           WRITE (hOUT,'(A,F6.2, A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN'
         ELSE
           WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]   (INPUT)'  
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)'           
             WRITE (hOUT,'(A,F6.2, A, A, F5.2,A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN',
     >                ' AND EDDINGTON_GAMMA = ', GEDD, ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A)') 
     >         '*  EDDINGTON_GAMMA = ', GEDD, '   (INPUT)'           
             WRITE (hOUT,'(A,F6.2, A, A, F5.2,A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN',
     >                ' AND LOG G_EFF = ', GEFFLOG, ' [CGS]'
           ENDIF
        ENDIF
      ELSE
         STOP 'PRIPARAM: ERROR - UNKNOWN CODE FOR MASSORIGIN'
      ENDIF

      IF (VTURB > .0) THEN
        WRITE (hOUT,'(A,F7.2,A)') '* VTURB   =', VTURB, ' KM/S'
      ENDIF

      WRITE (hOUT,'(A,F7.2,A,F8.2,A)') '* RMAX    =', RMAX, ' RSTAR =', 
     >                          RMAX*RSTARSU, ' R_SUN'


      WRITE(hOUT,80) 
   80 FORMAT ('*',/,80('*'),//)


      RETURN
      END
      SUBROUTINE PRIXDAT(XDATA,MAXXDAT)
C***  CALLED BY WRSTART
C***  PRINT OUT THE X-RAY-DATA

      DIMENSION XDATA(MAXXDAT)

C***  K-BOLTZMANN IN ERG/K
      DATA BOLZK/1.38E-16/

C***  ATOMIC-UNIT IN GR
      DATA AU/2.23E-24/

C***  MEAN PARTICLE MASS OF FULLY IONISIZED HELIUM
      COMION = 4 / 3 * AU
      
      
C***  VECTOR XDATA DEVIDED IN ITS PARTITIONS
      XFILL = XDATA(1)
      XRAYT = XDATA(2)
      XRMIN = XDATA(3)

      XFILL2 = XDATA(5)
      XRAYT2 = XDATA(6)
      XRMIN2 = XDATA(7)

      
C***  PRINTOUT THE X-RAY-INFORMATION
      PRINT 5
    5 FORMAT(/,10X,'INFORMATIONS ABOUT SHOCK-WAVE PARAMETERS : ')

C***  EXIT IF XFILL=0 NO X-RAY-EMISSION
      IF ((XFILL .EQ. 0.)) THEN
        PRINT 10
 10     FORMAT(10X,'NO X-RAY-DATA')
        RETURN

      ELSE

       IF (XFILL .GT. 0.) THEN
        XVELOC = SQRT(3 * BOLZK / COMION * XRAYT)/100000
        PRINT 15,XRAYT,XVELOC,XRMIN,XFILL
 15     FORMAT(/,10X,'X-RAY TEMPERATURE IN K:',E10.3,
     $         /,10X,'THAT MEANS A PARTICLE VELOCITY IN KM/SEC:',F10.2,
     $         /,10X,'MINIMUM RADIUS FOR SHOCKWAVE IN RSTARS:',F6.2,
     $  /,10X,'RELATIV NUMBER OF PARTICLE MAKING BREMSSTRAHLUNG:',E11.4,
     $         /)
       ENDIF

       IF (XDATA(4) .NE. .0) WRITE (*,16) XDATA(4)
 16     FORMAT (10X,'DIFFERENTIAL EMISSION MEASURE EXPONENT:', F5.2)
        WRITE (*,*)

       IF (XFILL2 .GT. 0.) THEN
        XVELOC2 = SQRT(3 * BOLZK / COMION * XRAYT2)/100000
        PRINT 17,XRAYT2,XVELOC2,XRMIN2,XFILL2
 17     FORMAT(/,10X,'X-RAY TEMPERATURE2 IN K:',E10.3,
     $         /,10X,'THAT MEANS A PARTICLE VELOCITY IN KM/SEC:',F10.2,
     $         /,10X,'MINIMUM RADIUS2 FOR SHOCKWAVE IN RSTARS:',F6.2,
     $  /,10X,'RELATIV NUMBER OF PARTICLE MAKING BREMSSTRAHLUNG:',E11.4,
     $         /)
       ENDIF

      ENDIF

      RETURN

      END
      SUBROUTINE READMS(ICHANNEL, X, NDIM, NAME, IERR)
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, X, NDIM, 
     >              'READ', IERR)

      RETURN
      END
      SUBROUTINE READOLDT (ICH, ND, NDDIM, T, RADIUS,
     $                     TOLD, ROLD, MODOLD, JOBNOLD, TEFF, TEFFOLD,
     >                     TAURCONT, TAURCONTOLD, BTAUR)
C**********************************************************************
C***  CALLED FROM: WRSTART
C***  READS TEMPEREATURE STRUCTURE FROM OLD MODEL FILE, IF REQUESTED
C***  AS THE RADIUS GRIDS MAY DEVIATE, THE NEW TEMPERATURE STRUCTURE
C***    IS OBTAINED BY INTERPOLATION
C**********************************************************************

      DIMENSION T(ND), RADIUS(ND), TOLD(ND), ROLD(ND)
      DIMENSION TAURCONT(ND), TAURCONTOLD(NDDIM)
      LOGICAL BTAUR, BTAUR_INTERPO

      BTAUR_INTERPO = BTAUR

C***  READ FROM CHANNEL ICH = OLD MODEL FILE ****************************
      CALL OPENMS (ICH, IDUMMY, IDUMMY, 1, IERR)
      IERR=1
      CALL READMS (ICH,MODOLD,13,'MODHEAD ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK (' OLD MODEL FILE NOT AVAILABLE')
         PRINT *,     ' OLD MODEL FILE NOT AVAILABLE'
         STOP 'ERROR'
         ENDIF
      CALL READMS (ICH,JOBNOLD,     1, 'JOBNUM  ' , IERR)
      CALL READMS (ICH,NDOLD  ,     1, 'ND      ' , IERR)
      CALL READMS (ICH,TEFFOLD,     1, 'TEFF    ' , IERR)
      CALL READMS (ICH,TOLD   , NDOLD, 'T       ' , IERR)
      CALL READMS (ICH,ROLD   , NDOLD, 'R       ' , IERR)

C***  ARRAY BOUND CHECK
      IF (NDOLD .GT. NDDIM) THEN
         CALL REMARK (' OLD MODEL HAS TOO MANY DEPTH POINTS')
         PRINT *,      'OLD MODEL HAS TOO MANY DEPTH POINTS'
         STOP 'ERROR'
      ENDIF

C***  OLD T TAU: CHECK IF TAURCONT IS ON MODEL OR FALLBACK NEEDED
      IF (BTAUR) THEN
         CALL READMS (ICH,TAURCONTOLD, NDOLD, 'TAURCONT' , IERR)
         IF (IERR .EQ. -10) THEN
            WRITE (0,*) 
     >       '*** WARNING: TAURCONT not on MODEL file, take TAUROSS'
            CALL READMS (ICH,TAURCONTOLD, NDOLD, 'TAUROSS ' , IERR)
            IF (IERR .EQ. -10) THEN
               WRITE (0,*) 
     >          '*** WARNING: TAUROSS not on MODEL file,' 
     >          // ' TAU option in OLD T disabled'
               BTAUR_INTERPO = .FALSE.
            ENDIF        
         ENDIF
      ENDIF
      
      CALL CLOSMS (ICH, IERR)


C***  INTERPOLATION ***************************************************
cc      BTAUR_INTERPO = BTAUR_INTERPO .AND. (TAURCONT(ND) .NE. 0.)
      BTAUR_INTERPO = BTAUR_INTERPO .AND. (TAURCONT(ND) > 0.)

      IF (BTAUR_INTERPO) THEN
         WRITE (0,*) 'Interpolation of Temperature on Tau-Grid'

cC***     scale TAUROLD
c         DO L=1, NDOLD
c            TAURCONTOLD(L) = 
c     >         TAURCONTOLD(L) * (TAURCONT(ND)+1.E-8)/TAURCONTOLD(NDOLD)
c         ENDDO
c Auskommentiert, und testweise durch das nachstehende IF ersetzt:
c     wrh  4-Mar-2019

         DO L=1, ND
            TAURL = TAURCONT(L)
            IF (TAURL .GT. TAURCONTOLD(NDOLD)) THEN
               T(L) = TOLD(NDOLD)
            ELSE             
               CALL LIPO (T(L),TAURL,TOLD,TAURCONTOLD,NDOLD)
            ENDIF
         ENDDO

      ELSE

         WRITE (0,*) 'Interpolation of Temperature on Radius-Grid'
         DO 1 L=1, ND
            IF (RADIUS(L) .GT. ROLD(1)) THEN
               T(L)=TOLD(1)
            ELSE
               CALL LIPO (T(L),RADIUS(L),TOLD,ROLD,NDOLD)
            ENDIF
    1    CONTINUE 

      ENDIF

C***  SCALE TEMPERATURE WITH THE RATIO TEFF ( NEW / OLD )
      IF (BTAUR_INTERPO) THEN
      IF (TEFF .NE. TEFFOLD) THEN
         Q = TEFF / TEFFOLD
         DO L=1, ND
C***  NOTE: THIS DOES NO LONGER GUARANTEE TMIN
C***        and can lead to non-monotonic situations
c           FTAU = EXP(-TAUROSS(L)*2.)
C           FTAU = 1.-MIN(TAUROSS(L), 1.)
c           T(L) = T(L) * (FTAU + (1.-FTAU) * Q)
           T(L) = T(L) * Q
         ENDDO
      ENDIF

      ELSE
      IF (TEFF .NE. TEFFOLD) THEN
         Q = TEFF / TEFFOLD
         DO 2 L=1, ND
    2    T(L) = T(L) * Q
         ENDIF

      ENDIF

      RETURN
      END
      SUBROUTINE REGULA(F,X,Y,X1,X2,EPS)
C***********************************************************************
C***  THIS ROUTINE CALCULATES THE SOLUTION X OF F(X)=Y IN THE INTERVAL
C***  (X1,X2) ! METHOD: 
C***  -------- REGULA FALSI --------
C***  USING BISECTION STEPS TO GUARANTEE CONVERGENCE, PRECISION IN X: EPS
C***********************************************************************
      LOGICAL BI
 
      BI=.TRUE.
      A=X1
      B=X2
      FA=F(A)-Y
      IF(FA.EQ..0) GOTO 3
      FB=F(B)-Y
      IF(FB.EQ..0) GOTO 4
      IF (FA*FB.GT..0) THEN
         WRITE (0,*) '*** INVALID ARGUMENTS WHEN CALLING REGULA ***'
         WRITE (0,10) A, B
   10    FORMAT ('INTERVAL:     A  =', E15.5, 5X, '  B  =', E15.5)
         WRITE (0,11) FA, FB
   11    FORMAT ('FUNCTION:   F(A) =', E15.5, 5X, 'F(B) =', E15.5)
         CALL TRBK
         STOP 'ERROR'
         ENDIF

    1 D=A-B
      IF(ABS(D).LT.EPS) GOTO 5
      BI=.NOT.BI
      IF(BI) X=A-FA*D/(FA-FB)
      IF(.NOT.BI) X=.5*(A+B)
      FX=F(X)-Y
      IF(FX.EQ..0) RETURN
      IF(FX*FA.GT..0) GOTO 2

      B=X
      FB=FX
      GOTO 1

    2 A=X
      FA=FX
      GOTO 1

    3 X=A
      RETURN

    4 X=B
      RETURN

    5 X=(A+B)*.5
      RETURN

      END
      SUBROUTINE REMARK (STRING)

      CHARACTER*(*) STRING

      WRITE (0,'(A)') STRING( :IDX(STRING))

      RETURN
      END
      SUBROUTINE RGRID (NDDIM,ND,R     ,INCRIT,RMAX, 
     >                  RadiusGridParameters)

C***********************************************************************
C***  SPECIFICATION OF THE GRID OF DEPTH-POINTS   ******************************
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: NDDIM
      INTEGER, INTENT(INOUT) :: ND 
      REAL, INTENT(IN) :: RMAX
      REAL, DIMENSION(NDDIM), INTENT(INOUT) :: R
      CHARACTER(8), DIMENSION(NDDIM), INTENT(INOUT) :: INCRIT
      CHARACTER(80), DIMENSION(3), INTENT(IN) :: RadiusGridParameters


      INTEGER :: ND2, ND3, NSPEC
      INTEGER :: N, NH, NH2, NDACT
      INTEGER :: N1, N2, ND1, ND1M, NH1
      INTEGER :: I, K, L, M, IERR
      REAL :: QSPEC, QSPEC_OUT
      REAL :: SUM, RI, DELR, Q, RL
      REAL :: RHO, RHOOLD, RHONEW, RHOLAST
      REAL :: V, VL, VLAST, TAUL
      REAL :: XND, XND2, XND3, DLOGTAU, DLOGR, DLOGR2

      !COMMON /COMRPAR/ RPAR
      REAL, DIMENSION(1000) :: RHELP, TAU
      !CHARACTER RPAR*80
      REAL :: SPEC_IN_SPACING, SPEC_OUT_SPACING
      INTEGER :: NPAR, NSPEC_IN, NSPEC_OUT
      CHARACTER(80) :: RPAR, ACTPAR
      CHARACTER(40), DIMENSION(20) :: CURPAR

      LOGICAL :: bOldDecode
      LOGICAL, DIMENSION(4) :: bParamFound

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      REAL, EXTERNAL :: WRVEL                 !own function for velocity field

 
C***  DECODE THE INPUT OPTION CARD WHICH SPECIFIES THE RADIUS GRID
      !RadiusGridParameters contains up to 3 lines from CARDS:
      ! RGRID, SPECIAL_OUTER_POINTS and SPECIAL_INNER_POINTS
      !RPAR contains the current CARDS line
      RPAR = RadiusGridParameters(1)
      !DECODE (80,7,RPAR ) XND,XND2,XND3,DLOGTAU
      
      bOldDecode = .FALSE.
      CALL SARGC (RPAR, NPAR)
      IF (NPAR < 5) THEN
        WRITE (hCPR,'(A)') '*** RGRID: NOT ENOUGH PARAMETERS'
        STOP '*** FATAL ERROR WHILE DECODING RGRID CARDS-LINE'
      ENDIF

      !New decoding => allows flexible format and modern syntax
      bParamFound = .FALSE.
      DO i=1, NPAR
        CALL SARGV(RPAR,i,CURPAR(I))
      ENDDO
      IF (NPAR > 2) THEN
        DO i=2, NPAR 
         SELECTCASE (CURPAR(i))
          CASE ('ND','NDTOTAL','ND_TOTAL','NDTOT','ND_TOT')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND
              IF (IERR == 0) THEN
                ND=IFIX(XND)
                bParamFound(1) = .TRUE.
              ENDIF      
            ENDIF
          CASE ('ND2','NDDENS','ND_DENS','DENS','NDRHO','ND_RHO','RHO')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND2
              IF (IERR == 0) THEN
                ND2=IFIX(XND2)
                bParamFound(2) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('ND3','NDVELO','ND_VELO','NDV','VELO','V')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND3
              IF (IERR == 0) THEN
                ND3=IFIX(XND3)
                bParamFound(3) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('DLOGTAU')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) DLOGTAU
              IF (IERR /= 0) THEN
                WRITE(hCPR,'(A)') '*** RGRID: CANNOT READ DLOGTAU'
                STOP 'ERROR'
              ELSE
                bParamFound(4) = .TRUE.
              ENDIF                  
            ENDIF
         ENDSELECT
        ENDDO
      ENDIF

      DO i=1, 4 
        !One or more parameters have not been found switch to old decoding
        IF (.NOT. bParamFound(i)) THEN
          WRITE (hCPR,*) '*** RGRID: Old RGRID decoding used'
          bOldDecode = .TRUE.
        ENDIF
      ENDDO

      IF (bOldDecode) THEN
        !Failsafe: old decoding with fixed format
        READ (UNIT=RPAR, FMT=7, ERR=92) XND,XND2,XND3,DLOGTAU
    7   FORMAT (10X,F5.0,4X,F5.0,4X,F5.0,8X,F10.0 )
        ND=IFIX(XND)
        ND2=IFIX(XND2)
        ND3=IFIX(XND3)
      ENDIF
      
C***  Decode input options for special points
      
      !default values for Special Outer and Inner Points
      NSPEC_IN         = 4
      SPEC_IN_SPACING  = 2.
      NSPEC_OUT        = 0
      SPEC_OUT_SPACING = 2.
      
      RPAR = RadiusGridParameters(2)
      CALL SARGC (RPAR, NPAR)
      IF (NPAR >= 2) THEN         
         CALL SARGV (RPAR, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=92) NSPEC_OUT
         IF (NPAR > 2) THEN
           CALL SARGV (RPAR, 3, ACTPAR)
           IF (ACTPAR == 'SPACING') THEN
             CALL SARGV (RPAR, 4, ACTPAR)
             READ (ACTPAR, '(F10.0)', ERR=92) SPEC_OUT_SPACING
           ENDIF
         ENDIF
      ENDIF
      RPAR = RadiusGridParameters(3)
      CALL SARGC (RPAR, NPAR)
      IF (NPAR >= 2) THEN
         CALL SARGV (RPAR, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=92) NSPEC_IN
         IF (NPAR > 2) THEN
           CALL SARGV (RPAR, 3, ACTPAR)
           IF (ACTPAR == 'SPACING') THEN
             CALL SARGV (RPAR, 4, ACTPAR)
             READ (ACTPAR, '(F10.0)', ERR=92) SPEC_IN_SPACING
           ENDIF
         ENDIF
      ENDIF      
      
      
      NSPEC=NSPEC_IN
      QSPEC=1. / SPEC_IN_SPACING
      QSPEC_OUT=1. / SPEC_OUT_SPACING
      
      IF (ND > NDDIM) THEN
          CALL REMARK ('ND .GT. NDIM')
          STOP 'ERROR'
      ENDIF
      ND1=ND-ND2-ND3-NSPEC - NSPEC_OUT
      IF (ND1 < 2) THEN
          CALL REMARK ('ND1 .LT. 2')
          STOP 'ERROR'
      ENDIF

 
C***  ND1 POINTS EQUALLY SPACED IN LOG TAU
 
C*** FIRST, A TABLE OF 1000 ENTRIES: RHELP(I), TAU(I) , IS ESTABLISHED
      N=1000
      NH1=800
      NH2=200
C***  THE RADIUS POINTS RHELP ARE SPACED LOGARITHMICALLY
      DLOGR=ALOG(RMAX)/FLOAT(NH1-1)
      SUM=.0
      RHOOLD=1./RMAX/RMAX/WRVEL(RMAX)
      TAU(1)=.0
      RHELP(1)=RMAX
      !Hilfsfeingitter aus 1000 (800+200) Punkten, feineres Spacing innen
      DO I=2,NH1-1
        RI=EXP((NH1-I)*DLOGR)
        RHELP(I)=RI
        RHONEW=1./RI/RI/WRVEL(RI)
        DELR=RI-RHELP(I-1)
        SUM=SUM+.5*DELR*(RHOOLD+RHONEW)
        RHOOLD=RHONEW
        TAU(I)=SUM
      ENDDO
      DLOGR2=LOG(RI)/FLOAT(NH2+1)
      DO I=NH1,N
        RI=EXP((N-I)*DLOGR2)
        IF (I == N) RI=1.
        RHELP(I)=RI
        RHONEW=1./RI/RI/WRVEL(RI)
        DELR=RI-RHELP(I-1)
        SUM=SUM+.5*DELR*(RHOOLD+RHONEW)
        RHOOLD=RHONEW
        TAU(I)=SUM
      ENDDO
 
      !Gitterpunkte nach TAU-Kriterium vergeben (gleichm. Spacing in DLOGTAU)
      ND1M=ND1-1
      DO L=2,ND1M
        TAUL=TAU(N)*10**(DLOGTAU*(FLOAT(L-2)/FLOAT(ND1-2)-1.))
        CALL LIPO (R(L),TAUL,RHELP,TAU,N)
        INCRIT(L)='TAU     '
      ENDDO
      R(1)=RMAX
      INCRIT(1)='RMAX    '
      R(ND1)=1.
 
C***  INSERT ND2 POINTS BY DENSITY CRITERION
      N2=ND1+ND2-1
      DO K=ND1,N2
        Q=1.
        RHOLAST=1./R(1)/R(1)/WRVEL(R(1))
        DO L=2,K
          RL=R(L)
          RHO=1./RL/RL/WRVEL(RL)
          IF (RHO/RHOLAST >= Q) THEN 
            Q=RHO/RHOLAST
            M=L
          ENDIF
          RHOLAST=RHO
        ENDDO
        CALL SHIFT (R,M,K)
        R(M)=.5*(R(M-1)+R(M+1))
        CALL SHIFTSTRING (INCRIT,M,K)
        INCRIT(M)='DENSITY '
      ENDDO
 
C***  INSERT ND3 POINTS BY VELOCITY CRITERION
      N1=ND1+ND2
      N2=ND1+ND2+ND3-1
      DO K=N1,N2
        Q=.0
        VLAST=WRVEL(R(1))
        DO L=2,K
          VL=WRVEL(R(L))
          IF (VLAST-VL >= Q) THEN
            Q=VLAST-VL
            M=L
          ENDIF
          VLAST=VL
        ENDDO
        CALL SHIFT (R,M,K)
        R(M)=.5*(R(M-1)+R(M+1))
        CALL SHIFTSTRING (INCRIT,M,K)
        INCRIT(M)='VELOCITY'
      ENDDO
 
C***  Insert NSPEC-OUT Points near the Outer Boundary
      DO I=1, NSPEC_OUT
        NDACT = ND - NSPEC_OUT + I - 2
        CALL SHIFT (R, 2, NDACT)
        CALL SHIFTSTRING (INCRIT, 2, NDACT)
C!!!        R(2) = 0.5 * (R(1) + R(3))
        R(2) = QSPEC_OUT * (R(3) - R(1)) + R(1)
        INCRIT(2)='SPECIAL '
      ENDDO
 
C***  INSERT NSPEC POINTS NEAR THE INNER BOUNDARY
      DO L=ND-NSPEC,ND-1
        R(L)=1.+(R(L-1)-1.)* QSPEC
        INCRIT(L)='SPECIAL '
      ENDDO
 
C***  INNER BOUNDARY
      R(ND)=1.
      INCRIT(ND)='INNER-B.'
 
C***  A special smoothing in radius is applied in order to get
C***    equally spaced depth points. This is done on the logarithmic 
C***    radius grid.         wrh+lars 17-Dec-1997 14:39:34
      DO L=2+NSPEC_OUT, ND-NSPEC_IN-1
        RHELP(L) = 0.25 * 
     >    (ALOG10(R(L-1)) + ALOG10(R(L+1)) + 2.*ALOG10(R(L)))
      ENDDO
      DO L=2+NSPEC_OUT, ND-NSPEC_IN-1
        R(L) = 10.**RHELP(L)
      ENDDO

      RETURN

      !Error message in case of CARDS line decoding fail
   92 WRITE (hCPR,*)
     >   'RGRID: ERROR WHILE DECODING THE FOLLOWING CARDS-LINE:'
      WRITE (hCPR,*) RPAR
      STOP 'ERROR'
      
      END
      SUBROUTINE RUKU(RSTART, REND, T, ISTEP)

c      write (*,*) 'rstart, rend, t, istep:',rstart, rend, t, istep

      H = (REND - RSTART) / ISTEP
      H2 = H / 2.

c      write (*,*) 'h,h2=',h,h2

      TL = T
      RL = RSTART

c      stop 'stop in ruku'

      DO I = 1, ISTEP

        RIN = RL
        TIN = TL
        CALL DTDR(1, RIN, TIN, TPRIME)
        X1 = H * TPRIME

        RIN = RL + H2
        TIN = TL + X1/2.
        CALL DTDR(1, RIN, TIN, TPRIME)
        X2 = H * TPRIME

        RIN = RL + H2
        TIN = TL + X2/2.
        CALL DTDR(1, RIN, TIN, TPRIME)
        X3 = H * TPRIME

        RIN = RL + H
        TIN = TL + X3
        CALL DTDR(1, RIN, TIN, TPRIME)
        X4 = H * TPRIME

        RL = RL + H
        TL = TL + ((X1/2.) + X2 + X3 + (X4/2.)) / 3.

c        write (*,*) 'rl, tl=',rl,tl

      ENDDO

      T = TL

      RETURN
      END
      SUBROUTINE SARGC(TEXT,N)

C**   Diese Funktion ist aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt die Anzahl der Argumente in einem String.
C**   Das Parsen des Strings uebernimmt die Routine sargp.
C**   Dort sind auch die syntaktischen Regeln beschrieben.

      CHARACTER*(*) TEXT
      INTEGER AS,AE,N,I

      I=0 				! Do not look for any Argument
      CALL SARGP(TEXT,N,I,AS,AE)

      RETURN

      END
      SUBROUTINE SARGP(TEXT,N,I,AS,AE)

C**	Die Subroutine sargp zerlegt einen String.
C**	Ermittelt werden:
C**	n : die Anzahl der Argumente und, 
C**     	falls i in [1..n],
C**	as und ae : Start- und Endindex des i-ten Arguments.
C**	Geparsed wird nach folgenden Regel:
C**	Leerzeichen werden grunsaetzlich nicht beachtet
C**	(Ausnahmen siehe unten).
C**	Argumente werden durch Leerzeichen oder Komma oder '=' oder ':'
C**     getrennt.
C**	Wird ein Argument durch Leerzeichen und Komma getrennt,
C**     so gilt dies als eine Trennung.
C**	Folgen zwei Kommata ohne Argument, also hoechstens durch
C**	Leerzeichen getrennt, gilt dies als Leerargument.
C**	Zeichen zwischen zwei doppelten Anfuehrungszeichen gelten
C**	als ein Argument, auch wenn Leerzeichen oder Kommata
C**	enthalten sind. In diesem Fall werden die Anfuehrungszeichen
C**	als nicht zum Argument gehoerig betrachtet.
C**	Interpretaion der Rueckgabewerte:
C**	n: immer die Anzahl der Argumente
C**	as=-1 -> Es wurde kein i-tes Argument gefunden (i nicht in [1..n])
C**	as=0  -> Das i-te Argument war leer (z.B. ",,")
C**	sonst : text(as:ae) = i-tes Argument

        CHARACTER*(*) TEXT
        INTEGER N                  ! Out: Anzahl der Argumente
        INTEGER I                  ! In : gesuchtes Argument
        INTEGER AS,AE              ! Out: erste u. letzte Position des 
                                   !         Arguments im String

	INTEGER TL,TI
        INTEGER STATE

	TL=LEN(TEXT)

	state=0
C	state=0 -> kein Argument aktiv
C	state=1 -> Argumentende gefunden, naechstes Komma
C			ist   k e i n   Leerargument
C	state=2 -> normales Argument aktiv
C	state=3 -> Argument in '"' aktiv


        N=0
        AS=-1
	AE=TL

	DO 100 TI=1,TL

C** Go here looking for next start of an argument
        IF (STATE .EQ. 0) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                            N=N+1
                            IF (N .EQ. I) THEN
                                AS=0
                            ENDIF
                            GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
			   IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 1) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			   STATE=0
			   GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
                           IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 2) THEN
           IF (TEXT(TI:TI) .EQ. ' ') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI
			  GOTO 100
	   ENDIF
	   IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			  STATE=0
                          IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
           GOTO 100
        ELSE 
C***  ! IF (STATE .EQ. 3)
           IF (TEXT(TI:TI) .EQ. '"') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
        ENDIF

100	CONTINUE

	RETURN
        END
      SUBROUTINE SARGV(TEXT,I,ARGTEXT)

C**   Diese Funktion ist fast aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt das i-te Argument in einem String.
C**   Die Regeln zur Argumenttrennung sind in sargp beschrieben.
C**   Ist das i-te Argument nicht vorhanden, wird argtext nicht veraendert.

      CHARACTER*(*) TEXT,ARGTEXT
      INTEGER I
      INTEGER N,AS,AE

      CALL SARGP(TEXT,N,I,AS,AE)

      IF (AS .EQ. -1) GOTO 10
      IF (AS .EQ. 0) THEN
             ARGTEXT=' '
      ELSE
             ARGTEXT=TEXT(AS:AE)
      ENDIF

10    RETURN

      END
      SUBROUTINE SECOND

      STOP 'SECOND NOT IMPLEMENTED AT DEC/OSF'

      RETURN
      END
      SUBROUTINE SEQUIN (NDIM,N,X,K,XNEW)
C***********************************************************************
C***  X(K) MUST BE A STRICKTLY MONOTONIC SEQUENCE
C***   THIS SUBROUTINE INSERTS THE POINT XNEW AT SUITABLE POSITION K
C***  NO INSERTION IN CASE OF EXACT COINCIDENCE (K:= NEGATIVE INDEX )
C***********************************************************************

      DIMENSION X(NDIM)
 
C***  INDEX ESTIMATE BY BISECTION
      NA=1
      A=X(1)
      NB=N
      B=X(N)
      IF (XNEW.EQ.A) GOTO 9
      IF (XNEW.EQ.B) GOTO 10
      IF((XNEW-A)*(XNEW-B).GE..0) GOTO 4
    1 IF(NB-NA .EQ. 1) GOTO 2
      NH=(NA+NB)/2
      H=X(NH)
      IF (XNEW.EQ.H) GOTO 11
      IF((XNEW-A)*(XNEW-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
 
C***  CASES OF EXACT COINCIDENCE
    9 K=-1
      RETURN
   10 K=-N
      RETURN
   11 K=-NH
      RETURN
 
      ENTRY SEQUINE  (NDIM,N,X,K,XNEW)
C***  ENTRY POINT IF THE INSERTION INDEX K IS GIVEN *********************
      IF (K.EQ.N+1) GOTO 8
      IF (K.LE.0 .OR. K.GT.N) THEN
            CALL REMARK ('ENTRY SEQUINE: WRONG INDEX GIVEN')
            STOP 'ERROR'
            ENDIF
      NB=K
 
C***  INSERTION OF XNEW AT INDEX NB
    2 IF(N+1 .GT. NDIM) THEN
            WRITE (0,*) 'ARRAY DIMENSION TOO SMALL'
            WRITE (0,*) 'NDIM=', NDIM
            STOP 'ERROR'
            ENDIF
      K=NB
    7 DO 5 J=K,N
      I=N+K-J
    5 X(I+1)=X(I)
    8 X(K)=XNEW
      N=N+1
      RETURN
 
C***  XNEW LIES OUTSIDE X(1) ... X(N)
    4 IF(N+1 .GT. NDIM) THEN
            WRITE (0,*) 'ARRAY DIMENSION TOO SMALL'
            WRITE (0,*) 'NDIM=', NDIM
            STOP 'FATAL ERROR IN SUBR. SEQUIN'
            ENDIF
      IF((XNEW-A)*(B-A).GT..0) GOTO 6
       K = 1
       GOTO 7
    6 N=N+1
      K=N
      X(N)=XNEW
      RETURN
      END
      SUBROUTINE SHIFTREAL (RealArray,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***    uses INTERFACE to overload shift routine to handle
C***     arrays of types real, integer and string(character)
C***********************************************************************
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L, ND
      REAL, DIMENSION(ND), INTENT(INOUT) :: RealArray
      INTEGER :: II, I

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFTREAL, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO II=L,ND
        I=ND+L-II
        RealArray(I+1)=RealArray(I)
      ENDDO

      RETURN
      END SUBROUTINE SHIFTREAL
      SUBROUTINE SHIFTSTRING (StringArray,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***    uses INTERFACE to overload shift routine to handle
C***     arrays of types real, integer and string(character)
C***********************************************************************

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L, ND
      CHARACTER(*), DIMENSION(ND), INTENT(INOUT) :: StringArray
      INTEGER :: II, I

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFTSTRING, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO II=L,ND
        I=ND+L-II
        StringArray(I+1)=StringArray(I)
      ENDDO

      RETURN
      END SUBROUTINE SHIFTSTRING
      SUBROUTINE SPLINPO (F, X, FI, XI, N)
C***********************************************************************
C***  CUBIC SPLINE INTERPOLATION, READY-FOR-USE
C***  XI(I), FI(I)  TABLE WHICH DEFINES THE FUNCTION TO BE INTERPOLATED
C***  X             ARGUMENT FOR WHICH THE FUNCTION VALUE IS REQUIRED
C***  FX            RESULTING FUNCTION VALUE
C***  THE RESULTING FUNCTION IS A PIECEWISE CUBIC INTERPOLATION 
C***     POLYNOMIAL WITH CONTINUOUS DERIVATIVE. EXTREMA CAN ONLY OCCUR 
C***     AT GIVEN MESHPOINTS (MONOTONIC VERSION AFTER M. STEFFEN)
C***********************************************************************

      DIMENSION XI(N), FI(N)

C***  CHECK FOR STRICTLY MONOTONIC ORDER
      DN = XI(N) - XI(1)
      DO 1 L=2, N
      DX = XI(L) - XI(L-1)
      IF (DX*DN .LE. .0) THEN
         WRITE (0, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3    FORMAT (' *** BAD USE OF SUBROUTINE SPLINPO:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
         CALL TRBK
         STOP 'ERROR'
         ENDIF
    1 CONTINUE

C***  FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
      DO 4 I=2, N
      L = I
      IF ( (X-XI(L-1)) * (X-XI(L)) .LE. .0) GOTO 2
    4 CONTINUE

      CALL REMARK('BAD USE OF SUBR. SPLINPO - X OUTSIDE TABLE')
      WRITE (0,'(A,G12.5,5X,A,G12.5,5X,A,I3,A,G12.5)') 
     >         ' X=', X,' X(1)=', XI(1),' X(', N, ')=', XI(N)
      CALL TRBK
      STOP 'ERROR'

    2 CONTINUE

C***  DETERMINATION OF THE COEFFICIENTS P1, P2, P3, P4 (CF. SUBR. CUBIC)
 
C***  SET UP THE COEFFICIENT MATRIX
      D1=1./(XI(L)-XI(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,N)
C***  FUNCTION TO BE INTERPOLATED: FI
      F1 = FI(L-1)
      F2 = FI(L)

      IF (.TRUE.) THEN
C***     Standard version: zentrierte Ableitungen an den Stuetzstellen. Das 
C***     ist bei nicht aequidistanten Stuetzstellen fragwuerdig, verringert 
C***     aber andererseits das Ueberschwingen insbesondere wenn man nicht 
C***     MONO verwendet 
         F3 = (FI(L) - FI(LA)) / (XI(L) - XI(LA))
         F4 = (FI(LB) - FI(L-1)) / (XI(LB) - XI(L-1))

      ELSE
C***  Alternative Version wrh 16-May-2002 13:39:29
C***  Statt der zentrierten Ableitung werden gewichtete Mittel der
C***  Ableitungen der angrenzenden Intervalle genommen. Das entspricht 
C***  der Ableitung eines Parabel-Fits durch drei Punkte.  
         FS0 = (FI(L) - FI(L-1)) / (XI(L) - XI(L-1))
         IF (L .EQ. 2) THEN 
            F3 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L) - XI(L-2))
            FSM = (FI(L-1) - FI(L-2)) / (XI(L-1) - XI(L-2))
            F3 = P * FSM + (1.-P) * FS0
         ENDIF
         IF (L .EQ. N) THEN 
            F4 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L+1) - XI(L-1))
            FSP = (FI(L+1) - FI(L)) / (XI(L+1) - XI(L))
            F4 = P * FSP + (1.-P) * FS0
         ENDIF
      ENDIF

C***  SET TRUE FOR MONO OPTION
      IF (.TRUE.) THEN

ccc   Diese bis heute (3-Sep-2002) verwendete Version erscheint mir 
ccc   merckwuerdig und an mehreren Stellen fehlerhaft! Ich lasse sie 
ccc   aus Dokumentationsgruenden hier stehen. Nachfolgend dann eine 
ccc   Version nach heutiger Erkenntnis. wrh  
c       S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
c       IF (LA .NE. L-2 .AND. LB .NE. L) THEN
c         S3 = S4
c         S5 = S4
c       ELSE
c         IF (LA .EQ. L-2) THEN
c           S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
c         ELSE
c           S3 = 1.4 * S4 - 0.5 * F3
c         ENDIF
c         IF (LB .EQ. L+1) THEN
c           S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
c         ELSE
c           S5 = 1.5 * S4 - 0.5 * F4
c         ENDIF
c       ENDIF

          S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
C***   We are not in the first interval:
          IF (LA .NE. L-2) THEN
            S3 = S4
          ELSE
            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
          ENDIF
C***   We are not in the last interval:
          IF (LB .NE. L+1) THEN
             S5 = S4
          ELSE
             S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
          ENDIF

       F3 = (SIGN(1.0,S3)+SIGN(1.0,S4))*MIN(ABS(S3),ABS(S4),0.5*ABS(F3))
       F4 = (SIGN(1.0,S4)+SIGN(1.0,S5))*MIN(ABS(S4),ABS(S5),0.5*ABS(F4))

      ENDIF
 
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
 

C***  EVALUATION OF THE INTERPOLATION POLYNOMIAL
      DXM = X - XI(L-1)
      DX  = XI(L) - X
      F = (P1 * DXM * DXM + P2 ) * DXM
     >  + (P3 * DX  * DX  + P4 ) * DX

      RETURN
      END
      SUBROUTINE SPLINPOX(F, X, FI, XI, N, SAFE, LPARAM, DFDX, D2FD2X)
C***********************************************************************
C***  CUBIC SPLINE INTERPOLATION, READY-FOR-USE
C***  XI(I), FI(I)  TABLE WHICH DEFINES THE FUNCTION TO BE INTERPOLATED
C***  X             ARGUMENT FOR WHICH THE FUNCTION VALUE IS REQUIRED
C***  FX            RESULTING FUNCTION VALUE
C***  THE RESULTING FUNCTION IS A PIECEWISE CUBIC INTERPOLATION 
C***     POLYNOMIAL WITH CONTINUOUS DERIVATIVE. EXTREMA CAN ONLY OCCUR 
C***     AT GIVEN MESHPOINTS (MONOTONIC VERSION AFTER M. STEFFEN)
C
C     Unified version implementing all features from classic routines
C       SPLINPO, SPLINPO_FAST and SPLINP from Goetz branch (A. Sander, Jan 2012)
C
C     Due to the optional arguments (e.g. for the derivatives), all main routines
C     need to have the following interface block:
C
C      INTERFACE SPLINPO
C        SUBROUTINE SPLINPO(F, X, FI, XI, N, DFDX, D2FD2X)
C          INTEGER, INTENT(IN) :: N          
C          REAL, DIMENSION(N), INTENT(IN) :: XI, FI
C          REAL, INTENT(OUT) :: F
C          REAL, INTENT(IN) :: X
C          LOGICAL, INTENT(IN), OPTIONAL :: SAFE
C          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
C          REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X
C        END SUBROUTINE
C      END INTERFACE SPLINPO
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      
      REAL, DIMENSION(N), INTENT(IN) :: XI, FI
      REAL, INTENT(OUT) :: F
      REAL, INTENT(IN) :: X
      LOGICAL, INTENT(IN), OPTIONAL :: SAFE
      INTEGER, INTENT(IN), OPTIONAL :: LPARAM
      REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X

      REAL :: DN, DX, DXM, FS0,
     >        D1, D2, D3, D23, H11, H12, H13, H14,
     >        H21, H22, H23, H24, H31, H32, H33, H34,
     >        H41, H42, H43, H44,
     >        F1, F2, F3, F4, FSM, FSP, S3, S4, S5,
     >        P, P1, P2, P3, P4

      INTEGER :: L, I, LA, LB

      LOGICAL :: bFoundInterval, bSafeMode

      !File and channel handles
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqXX.cpr

      SAVE      !this is required due to the _SAME_X entry point

      IF (present(SAFE)) THEN
        bSafeMode = SAFE
      ELSE
        bSafeMode = .TRUE.              !safe mode is default
      ENDIF

C***  CHECK FOR STRICTLY MONOTONIC ORDER (safe mode only)
      IF (bSafeMode) THEN
        DN = XI(N) - XI(1)
        DO L=2, N
          DX = XI(L) - XI(L-1)
          IF (DX*DN <= .0) THEN
            WRITE (hCPR, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3       FORMAT (' *** BAD USE OF SUBROUTINE SPLINPOX:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
            CALL TRBK
            STOP 'ERROR'
          ENDIF
        ENDDO

C***    FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
        bFoundInterval = .FALSE.
        DO I=2, N
          L = I
          IF ( (X-XI(L-1)) * (X-XI(L)) <= .0) THEN
            bFoundInterval = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF (.NOT. bFoundInterval) THEN
          CALL REMARK('BAD USE OF SUBR. SPLINPOX - X OUTSIDE TABLE')
          WRITE (hCPR,'(A,G12.5,5X,A,G12.5,5X,A,I3,A,G12.5)') 
     >         ' X=', X,' X(1)=', XI(1),' X(', N, ')=', XI(N)
          CALL TRBK
          STOP 'ERROR'
        ENDIF
        IF (present(LPARAM)) THEN
          IF ((LPARAM > 0) .AND. (LPARAM /= L)) THEN
            !If LPARAM has been specified (and is > 0), it must have the correct value
            WRITE (hCPR,'(A, I5)') '**** LPARAM =', LPARAM
            WRITE (hCPR,'(A, I5)') '**** CORRECT L =', L
            STOP ' *** ERROR IN SUBR. SPLINPO - WRONG INDEX LPARAM'
          ENDIF
        ENDIF

      ENDIF

      !Note: Unlike in other languages, the following two IF statements
      !      cannot be combined, because the second check would be made
      !      even if present(LPARAM) already returns FALSE. Therere it
      !      would cause a crash if LPARAM is not set 
      IF (present(LPARAM)) THEN
        IF (LPARAM >= 0) THEN
          !Preset L if LPARAM has been specified and is inside 2..N
          L = LPARAM
        ENDIF
      ENDIF

C***  DETERMINATION OF THE COEFFICIENTS P1, P2, P3, P4 (CF. SUBR. CUBIC)
 
C***  SET UP THE COEFFICIENT MATRIX
      D1=1./(XI(L)-XI(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,N)

C      WRITE (hCPR,*) 'Debug: L-2=', L-2, '  LA=', LA

C***  Entry point for a subsequent call with same x point, 
C***     but different function 
C      ENTRY SPLINPOX_SAME_X (F, X, FI, XI, N, SAFE)

C***  FUNCTION TO BE INTERPOLATED: FI
      F1 = FI(L-1)
      F2 = FI(L)

      IF (.TRUE.) THEN
C***     Standard version: zentrierte Ableitungen an den Stuetzstellen. Das 
C***     ist bei nicht aequidistanten Stuetzstellen fragwuerdig, verringert 
C***     aber andererseits das Ueberschwingen insbesondere wenn man nicht 
C***     MONO verwendet 
         F3 = (FI(L) - FI(LA)) / (XI(L) - XI(LA))
         F4 = (FI(LB) - FI(L-1)) / (XI(LB) - XI(L-1))

      ELSE
C***  Alternative Version wrh 16-May-2002 13:39:29
C***  Statt der zentrierten Ableitung werden gewichtete Mittel der
C***  Ableitungen der angrenzenden Intervalle genommen. Das entspricht 
C***  der Ableitung eines Parabel-Fits durch drei Punkte.  
         FS0 = (FI(L) - FI(L-1)) / (XI(L) - XI(L-1))
         IF (L == 2) THEN 
            F3 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L) - XI(L-2))
            FSM = (FI(L-1) - FI(L-2)) / (XI(L-1) - XI(L-2))
            F3 = P * FSM + (1.-P) * FS0
         ENDIF
         IF (L == N) THEN 
            F4 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L+1) - XI(L-1))
            FSP = (FI(L+1) - FI(L)) / (XI(L+1) - XI(L))
            F4 = P * FSP + (1.-P) * FS0
         ENDIF
      ENDIF

C***  SET TRUE FOR MONO OPTION
      IF (.TRUE.) THEN

ccc   Diese bis heute (3-Sep-2002) verwendete Version erscheint mir 
ccc   merckwuerdig und an mehreren Stellen fehlerhaft! Ich lasse sie 
ccc   aus Dokumentationsgruenden hier stehen. Nachfolgend dann eine 
ccc   Version nach heutiger Erkenntnis. wrh  
c       S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
c       IF (LA .NE. L-2 .AND. LB .NE. L) THEN
c         S3 = S4
c         S5 = S4
c       ELSE
c         IF (LA .EQ. L-2) THEN
c           S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
c         ELSE
c           S3 = 1.4 * S4 - 0.5 * F3
c         ENDIF
c         IF (LB .EQ. L+1) THEN
c           S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
c         ELSE
c           S5 = 1.5 * S4 - 0.5 * F4
c         ENDIF
c       ENDIF

          S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
C***   We are not in the first interval:
          IF (LA /= L-2) THEN
            S3 = S4
          ELSE
            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
          ENDIF
C***   We are not in the last interval:
          IF (LB /= L+1) THEN
             S5 = S4
          ELSE
             S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
          ENDIF

       F3 = (SIGN(1.0,S3)+SIGN(1.0,S4))*MIN(ABS(S3),ABS(S4),0.5*ABS(F3))
       F4 = (SIGN(1.0,S4)+SIGN(1.0,S5))*MIN(ABS(S4),ABS(S5),0.5*ABS(F4))

      ENDIF
 
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
 

C***  EVALUATION OF THE INTERPOLATION POLYNOMIAL
      DXM = X - XI(L-1)
      DX  = XI(L) - X
      F = (P1 * DXM * DXM + P2 ) * DXM
     >  + (P3 * DX  * DX  + P4 ) * DX

C***  Calculation of derivatives (optional)
      !added on 06.10.2011 to provide the same functionality as goetz
      IF (present(DFDX)) THEN
        DFDX = 3. * P1 *  DXM * DXM + P2 - 3. * P3 * DX * DX - P4
      ENDIF
      IF (present(D2FD2X)) THEN
        D2FD2X = 6. * (P1 *  DXM + P3 * DX)
      ENDIF

      RETURN
      END
      SUBROUTINE STAMP (OPSYS, PROGNAME, TIM1)
C**************************************************************
C***  Writes a time stamp in channel 0 (cpr-file)
C***    and the elapsed CPU time into standard-out
C**************************************************************

      CHARACTER OPSYS*(*), PROGNAME*(*), TIM1*(*), TIM2*10 
c      CHARACTER OPSYS*(*), PROGNAME*(*), TIM1*10, TIM2*10 
      REAL*4 DTIME, ETIME, TARRAY(2)



      IF (OPSYS .EQ. 'CRAY') THEN
        CALL CLOCK(TIM2)
        WRITE (0,'(A,A8,2X,3A,F8.1)')
     >      'My Wallclock: ', PROGNAME, TIM1, TIM2, 
     >      'CPU-sec.:', SECOND()
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', SECOND(), ' sec'
      ELSEIF (OPSYS .EQ. 'DEC/UNIX') THEN
           CALL DATE_AND_TIME (DATE=TIM1, TIME=TIM2)
ccc        CALL TIME(TIM2)
ccc        ET = ETIME(TARRAY)
ccc        DTIME gives decimals:
        ET = DTIME(TARRAY)
        WRITE (0,'(A,A8,2X,A,1X,A,1X,A,1X,F8.1)')
     >      'My Wallclock: ', PROGNAME, 
     >       TIM1(1:4) // '/' // TIM1(5:6) // '/' // TIM1(7:8),  
     >       TIM2(1:2) // ':' // TIM2(3:4) // ':' // TIM2(5:6),  
     >      'CPU-sec.:', ET
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', ET, ' sec'
      ELSEIF (OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM2)
        ET = ETIME(TARRAY)
        WRITE (0,'(A,A8,2X,2(A8,1X),A,3X,F8.1)')
     >      'My Wallclock: ', PROGNAME, TIM1, TIM2, 
     >      'CPU-sec.:', ET
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', ET, ' sec'
      ENDIF
 
      RETURN
      END
      SUBROUTINE STORAGE (ICHANNEL, IADR, MAXADR, 
     >                    CNAME,    !variable name (or file status if opening a file)
     >                    CNAME2,   !used for: 2nd var (CHANGE), format (MSINFO), kind (WRITE), file name (OPEN)
     >                    X, NDIM, ACTION, MODE, IERR)

C******************************************************************
C***  MASS-STORAGE EMULATOR FOR DEC/OSF1 BY LARS KOESTERKE
C***  VERSION 2.0     21-Jul-1997 13:41:16
C***    New Feautures : INFO-D
C***  VERSION 3.0     15-May-1998 14:16:39
C***    IADR is now stored for several MS-Files 
C***    SOPEN and SCLOSE are new Actions
C******************************************************************

C***  IRECL : RECORD LENGTH OF THE FILE OPENENED BY FORTRAN. THE VARIABLE 
C***                     LENGTH IS 4 BYTE
C***   => changed with compiler parameter "assume byterecl" to be compartible with gfortran
C***  IADRL = IRECL / 2, BECAUSE ALL VARIABLES, EVEN THE CHARACTERS, SHOULD 
C***                     HAVE 8 BYTE

C***  Variable-names have 8 Characters. Note that the Character ^ is not valid
C***    because it is used to transport blanks from the MSFILE plotutility
C***    mcplot.com to msinfo.com (INFO-D)

C***  For INTEL Compiler, 
C***    unless the compiler option "-assume byterecl" is set:
ccc      PARAMETER (IRECL = 256)
C***  For gfortran compiler:
      PARAMETER (IRECL = 1024)

      PARAMETER (IADRL = 128)

      DIMENSION IADR(MAXADR), ISCRATCH(IADRL), SCRATCH(IADRL)
      DIMENSION X(NDIM)
      INTEGER(KIND=1), DIMENSION(IADRL*8) :: SCRATCHBYTE   !for msinfo only - works just with 8 byte defaults

C***  MODE IS NOT USED SO FAR

      CHARACTER(1) :: CKIND
      CHARACTER(7) :: CSTATUS, FN
      CHARACTER(8) :: ACTION, MODE, KSTR, KINDSTR, CNAME, CNAME2
      CHARACTER(10) :: FMTSTR, CACTION

      REAL :: RDUMMY, RSCRATCH
      INTEGER :: IKIND, IDEFKIND, NDIMR, LFMT, IFMTLEN, INTSCRATCH

      LOGICAL BEXIST, BVARKN, BDIMEQ, BNEWKINDINFO, bWRINT

      INTEGER, EXTERNAL :: IDX  !function to obtain the (non-blank) length of a string

C      SAVE IRECL2, NIND2, NIND2U, NINDEX, NREC, NVAR, LASTCH
      SAVE

C***  IWARN = 0  :  No Warning
C***          1  :  Warnings when Reading unknown Variables
C***          2  :  Verbose output
      IWARN = 0

C***  Substitute ^ in blanks in CNAME and CNAME2
      IF (ACTION .EQ. 'INFO-D') THEN
        DO I=1,8
          IF (CNAME(I:I)  .EQ. '^') CNAME(I:I)  = ' '
          IF (CNAME2(I:I) .EQ. '^') CNAME2(I:I) = ' '
        ENDDO
C***  Check for the substring '^' in CNAME and CNAME2. This is not longer
C***  allowed
      ELSEIF (ACTION .EQ. 'READ' .OR. 
     >        ACTION .EQ. 'WRITE' .OR. 
     >        ACTION .EQ. 'LENGTH' .OR. 
     >        ACTION .EQ. 'CHANGE') THEN
        DO I=1,8
          IF (CNAME(I:I)  .EQ. '^' .OR. 
     >        CNAME2(I:I) .EQ. '^') THEN
            WRITE (0,*) 'The Substring ^ is not allowed in the Names', 
     >        CNAME, CNAME2
            STOP 'ERROR in Subr. STORAGE'
          ENDIF
        ENDDO
      ENDIF

      IF (MODE(1:4) .EQ. 'CRAY') THEN
        READ(UNIT=CNAME, FMT='(A8)') NAME
        READ(UNIT=CNAME2, FMT='(A8)') NAME2
c    1   FORMAT (A8)
      ELSE
        WRITE (0,*) 'MODE NICHT CRAY'
        WRITE (0,*) 'MODE=',MODE(1:4)
        STOP 'ERROR IN STORAGE'
      ENDIF

      IF (IWARN .GE. 2) THEN
        WRITE(0,'(A,A8,A,I2,2X,A,A)') 
     >        'storage: action=',action,' ICHANNEL=',ICHANNEL, 
     >        'NAME=',NAME
      ENDIF

C***  CHECK THE ICHANNEL NUMBER
      IF (ICHANNEL .LE. 0) THEN
        WRITE (0,*) ' NEGATIVE ICHANNEL NUMBERS ARE NOT ALLOWED'
        STOP 'ERROR IN STORAGE'
      ENDIF

C***  CHECK FOR MINIMUM LENGTH OF THE INDEX ARRAY IADR
      IF (MAXADR .LT. IRECL) THEN
        WRITE (0,*) ' DIMENSION (MAXADR) OF INDEX ARRAY IADR SMALLER',
     >              ' THAN RECORD LENGTH (IRECL) OF MASS-STORAGE FILE'
        STOP 'ERROR IN ROUTINE STORAGE'
      ENDIF

C***  NUMBER OF INDEX RECORDS (NIND) CLAIMED BY MAXADR
      NIND = (MAXADR - 1) / IADRL + 1
c      write (*,*) 'STORAGE : action, nind=', action,nind

C***                        ====
      IF (ACTION( :4) .EQ. 'OPEN') THEN
C***                        ====

C***  FILE-NAME ON DISK: fort.<ICHANNEL>
        IF (IDX(CNAME2) > 0) THEN
          FN = CNAME2(1:7)
        ELSEIF (ICHANNEL < 10) THEN
          WRITE(FN, '("fort.",I1,1X)') ICHANNEL
        ELSE
          WRITE(FN, '("fort.",I2   )') ICHANNEL
        ENDIF

C***  CHECK IF FILE DOES EXIST
        INQUIRE (FILE=FN, EXIST=BEXIST)
        CACTION = 'READWRITE'
        IF (CNAME == 'AUTO') THEN
          IF (BEXIST) THEN
            CSTATUS = 'OLD'
          ELSE
            CSTATUS = 'NEW'
          ENDIF
        ELSEIF (CNAME == 'READ') THEN
          CSTATUS = 'OLD'
          CACTION = 'READ'
        ELSE
          CSTATUS = CNAME(1:7)
          IF ((CSTATUS == 'REPLACE') .OR. (CSTATUS(1:3) == 'NEW')) THEN
            BEXIST = .FALSE.
          ENDIF
        ENDIF

        OPEN (UNIT=ICHANNEL,
     >        FILE=FN,
     >        ACCESS='DIRECT',
     >        FORM='UNFORMATTED',
     >        RECL=IRECL,
     >        STATUS=CSTATUS,
     >        ACTION=CACTION,
     >        ERR=90,
     >        IOSTAT=IOS)

C***  READ IADR IF FILE EXISTS
C***    READ FIRST RECORD
        IF (BEXIST) THEN
          READ (UNIT=ICHANNEL, REC=1, ERR=91, IOSTAT=IOS) ISCRATCH
          DO I=1, IADRL
            IADR(I) = ISCRATCH(I)
          ENDDO

C***  INTERPRET THE FIRST ELEMENTS
C***          IRECL2   : RECORD LENGTH OF THE EXISTING FILE
C***          NIND2    : NUMBER OF INDEX RECORDS OF THE EXISTING FILE
C***          NIND2U   : NUMBER OF INDEX RECORDS USED IN THE EXISTING FILE
C***          NINDEX   : NUMBER OF INDICES IN INDEX RECORD
C***          NREC     : TOTAL NUMBER OF RECORDS
C***          NVAR     : NUMBER OF VARIABLES
C***          THE INDEX ARRAY (IADR) IS USED FOR THE STORAGE OF THE 
C***            INFORMATION ABOUT THE ARRAYS STORED IN THE FILE
C***            FIVE ENTRYS ARE USED FOR EACH ARRAY
C***          IADR(11) : NAME OF THE ARRAY
C***          IADR(12) : NUMBER OF FIRST RECORD
C***          IADR(13) : NUMBER OF RECORDS
C***          IADR(14) : NUMBER OF VARIABLES
C***          IADR(15) : DATA TYPE and KIND (was UNUSED until May 2012)
C***          IADR(16) : etc. (next entry!)
          IRECL2 = IADR(1)
          NIND2  = IADR(2)
          NIND2U = IADR(3)
          NINDEX = IADR(4)
          NREC   = IADR(5)
          NVAR   = IADR(6)

C***  CONSISTENCY CHECKS

C***    COMPARE RECORD LENGTH
C***      the following is no longer fatal and therefore not reported
c          IF (IRECL .NE. IRECL2) THEN
c            WRITE (0,*) ' RECORD LENGTH OF FILE (IRECL2) AND',
c     >                  ' ROUTINE (IRECL) DO NOT MATCH'
c            WRITE (0,'(A7,I4,A8,I4)') ' IRECL=',IRECL, ' IRECL2=',IRECL2
c            WRITE (0,'(A)') '- If compiled with "ifort", make sure that'
c     >                   // ' the option "-assume byterecl" was used!'
cC!            STOP 'ERROR IN STORAGE'
c          ENDIF

C***    COMPARE NUMBER OF INDEX RECORDS
          IF (NIND .EQ. NIND2) THEN
          ELSE IF (NIND .GT. NIND2) THEN
cc  message suppressed! wrh 12-Aug-2008 13:58:44
cc            WRITE (0,*) 'INFO from STORAGE: ' // FN 
cc     >       // ' has less Index-Records than dimensioned > expanded'
          ELSE
            WRITE (0,*) 'Number of Index-Records in File is greater'
            WRITE (0,*) 'than in Index-Array'
C            WRITE (0,*) ' MIN IS TAKEN?'
            WRITE (0,'(A6,I4,A6,I4)') 'NIND2=', NIND2, ' NIND=', NIND
            WRITE (0,'(A,I3)') 'CHANNEL=',ICHANNEL
            STOP 'ERROR IN STORAGE'
          ENDIF

C***    READ THE REST OF THE INDEX ARRAY
          DO I=2, NIND2U
            READ (UNIT=ICHANNEL, REC=I) ISCRATCH
            DO J=1, IADRL
              IADR((I-1)*IADRL + J) = ISCRATCH(J)
            ENDDO
          ENDDO

        ELSE
C***  NEW FILE OPENED
          IRECL2 = IRECL
          NIND2  = (MAXADR - 1) / IADRL + 1
          NIND2U = 1
          NINDEX = 0
          NREC   = NIND2
          NVAR   = 0
        ENDIF

C***                        ====
      ELSE IF (ACTION( :5) .EQ. 'SOPEN') THEN
C***                        ====
          IRECL2 = IADR(1)
          NIND2  = IADR(2)
          NIND2U = IADR(3)
          NINDEX = IADR(4)
          NREC   = IADR(5)
          NVAR   = IADR(6)

C***                             =====
      ELSE IF (ACTION( :5) .EQ. 'WRITE') THEN
C***                             =====
C***  TRANSFORM SIZE IF KIND IS DIFFERENT 
        NDIMR = NDIM
        IF (IDX(CNAME2) > 0) THEN
          IF ( (CNAME2(1:1) == 'I') .OR.
     >         (CNAME2(1:1) == 'i') .OR.
     >         (CNAME2(1:1) == 'R') .OR.
     >         (CNAME2(1:1) == 'r') ) THEN
            !kind calculations only required for data types which can have different kinds (integer, real)
            KINDSTR = ''
            DO J=2, 8
              SELECTCASE(CNAME2(J:J))
                CASE ('1':'9', '0')
                  KINDSTR = TRIM(ADJUSTL(KINDSTR)) // CNAME2(J:J)
                CASE DEFAULT
                  EXIT
              ENDSELECT
            ENDDO
            READ(UNIT=KINDSTR, FMT='(I8)') IKIND       !transform number part from kind string into integer
            !Get default kind for used type
            IF ((CNAME2(1:1) == 'I') .OR. (CNAME2(1:1) == 'i')) THEN
              IDEFKIND = KIND(IDEFKIND)
            ELSE
              IDEFKIND = KIND(RDUMMY)
            ENDIF
            !rescale dimension with used kind
            NDIMR = NDIMR * IKIND / IDEFKIND
          ENDIF
        ENDIF

C***  CHECK IF VARIABLE IS KNOWN
        BVARKN = .FALSE.
        INUM = 0
        BNEWKINDINFO = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
c          WRITE (0,'(A9,A8,1x,a8,I6)') 'VARIABLE=',NAME,IADR(1+INDEX),
c     >    IADR(4+INDEX)
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            INUM   = IADR(4+INDEX)
            IF ((IDX(CNAME2) > 0) .AND. (NAME2 /= IADR(5+INDEX))) THEN
              !new identification string => can only be used if total byte size matches
              BNEWKINDINFO = .TRUE.
              IF (IADR(INDEX+5) == -1) THEN
                !no format set so far
                BDIMEQ = (NDIMR == INUM)                   
              ELSE
                WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)  !transform from integer into character
                KINDSTR = ''
                DO J=2, 8
                  SELECTCASE(KSTR(J:J))
                    CASE ('1':'9', '0')
                      KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                    CASE DEFAULT
                      EXIT
                  ENDSELECT
                ENDDO
                READ(UNIT=KINDSTR, FMT='(I8)') IKIND      !transform number part from kind string into integer
                IDEFKIND = KIND(RDUMMY)
                BDIMEQ = (NDIMR == NDIM * IKIND / IDEFKIND)
              ENDIF
            ELSEIF (INUM .EQ. NDIM) THEN
              BDIMEQ = .TRUE.
            ELSE
              BDIMEQ = .FALSE.
            ENDIF
            IF (BDIMEQ) THEN
              BVARKN = .TRUE.
              IFIRST = IADR(2+INDEX)
C***          EXIT THIS LOOP
              GOTO 40
            ELSE
              WRITE (0,*) ' INCONSISTENCE: VARIABLE FOUND BUT WITH',
     >                    ' DIFFERENT ARRAY LENGTH'
              WRITE (0,'(A9,A8)') 'VARIABLE=',NAME
              WRITE (0,'(A21,I6)') 'DIMENSION OF ARRAY : ', NDIM
              WRITE (0,'(A21,I6)') 'DIMENSION IN FILE  : ', INUM
              STOP 'ERROR IN STORAGE (ACTION : WRITE)'
            ENDIF
          ENDIF
        ENDDO
   40   CONTINUE

C***    NUMBER OF RECORDS USED FOR THE ARRAY WHICH WILL BE STORED
        NIND3 = (NDIMR - 1) / IADRL + 1
        IF (.NOT. BVARKN) THEN
C***    UPDATE INDEX AND APPEND NAME OF THE ARRAY
          NVAR = NVAR + 1
          IFIRST = NREC + 1
          NREC = NREC + NIND3
          NIND2U = (10 + NVAR*5 - 1) / IADRL + 1
C          IF (NIND2U .GT. NIND2) THEN
C            WRITE (0,*) ' INDEX ARRAY IS NOT LARGE ENOUGH TO',
C     >                  ' RECEPT A NEW ARRAY'
C            WRITE (0,*) 'NIND2=', NIND2
C            STOP 'ERROR IN STORAGE'
C          ENDIF
          NINDEX = NINDEX + 1
          INDEX = 10 + (NINDEX-1)*5
          IADR(INDEX+1) = NAME
          IADR(INDEX+2) = IFIRST
          IADR(INDEX+3) = NIND3
          IADR(INDEX+4) = NDIM
          IF (IDX(CNAME2) > 0) THEN
            IADR(INDEX+5) = NAME2
          ELSE
            IADR(INDEX+5) = -1
          ENDIF
        ELSEIF (BNEWKINDINFO) THEN
          IADR(INDEX+4) = NDIM
          IADR(INDEX+5) = NAME2
        ENDIF
        
C!!!  OLD VERSION (MIT UMKOPIEREN DES GESAMTEN ARRAYS)
C***    AS LONG AS THE ARRAY FILLS THE NEXT RECORD COPY IT TO SCRATCH
C***    THE LAST RECORD IS, IF NECESSARY, FILLD UP WITH 0.
          NREST = NDIMR
          DO I=1, NIND3
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)
              SCRATCH(J) = X(J+INDEX)
            ENDDO
            DO J=NREST+1, IADRL
              SCRATCH(J) = 0.
            ENDDO
            WRITE(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=93, IOSTAT=IOS) 
     >            SCRATCH
            NREST = NREST - IADRL
          ENDDO
C!!!          DO I=1, NIND3
C!!!            INDEX = (I-1) * IADRL
C!!!            WRITE(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=93, IOSTAT=IOS) 
C!!!     >            X(INDEX+1)
C!!!          ENDDO

C***                             ====
      ELSE IF (ACTION( :4) .EQ. 'READ') THEN
C***                             ====
C***  CHECK IF VARIABLE IS KNOWN
        NDIMR = NDIM
        BVARKN = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
c      write (0,'(A,2A8)') 'testname=',IADR(1+INDEX), name
c      write (0,'(A,2I)') 'testname=',IADR(1+INDEX), name
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            IF (IADR(5+INDEX) /= -1) THEN
              !Rescale NDIM (or copy NDIMR) if different KIND specified
              WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)  !transform from integer into character
              KINDSTR = ''
              IF ( (KSTR(1:1) == 'I') .OR.
     >             (KSTR(1:1) == 'i') .OR.
     >             (KSTR(1:1) == 'R') .OR.
     >             (KSTR(1:1) == 'r') ) THEN
                !kind calculations only required for data types which can have different kinds (integer, real)
                DO J=2, 8
                  SELECTCASE(KSTR(J:J))
                    CASE ('1':'9', '0')
                      KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                    CASE DEFAULT
                      EXIT
                  ENDSELECT
                ENDDO
                READ(UNIT=KINDSTR, FMT='(I8)') IKIND      !transform number part from kind string into integer
                !Get default kind for used type
                IF ((KSTR(1:1) == 'I') .OR. (KSTR(1:1) == 'i')) THEN
                  IDEFKIND = KIND(IDEFKIND)
                ELSE
                  IDEFKIND = KIND(RDUMMY)
                ENDIF
                !rescale dimension with used kind
                NDIMR = NDIMR * IKIND / IDEFKIND
              ENDIF
            ENDIF
            IF (IADR(4+INDEX) .EQ. NDIMR) THEN
              BVARKN = .TRUE.
            ELSEIF (IADR(4+INDEX) .GT. NDIMR) THEN
              BVARKN = .TRUE.
              IF (IWARN .GE. 2) THEN
                WRITE (0,*) ' WARNING from Subroutine STORAGE'
                WRITE (0,*) ' Inconsistence: Variable found but with',
     >                      ' different Array Length (Action: READ)'
                WRITE (0,'(A9,A8)') 'Variable=',NAME
                WRITE (0,'(A21,I6)') 'Dimension of Array : ', NDIMR
                INUM   = IADR(4+INDEX)
                WRITE (0,'(A21,I6)') 'Dimension in File  : ', INUM
              ENDIF
            ELSE
              BVARKN = .TRUE.
              INUM   = IADR(4+INDEX)
              IF (IWARN .GE. 1) THEN
                WRITE (0,*) ' WARNING from Subroutine STORAGE'
                WRITE (0,*) ' Inconsistence: Variable found but with',
     >                      ' different Array Length (Action: READ)'
                WRITE (0,'(A9,A8)') 'Variable=',NAME
                WRITE (0,'(A21,I6)') 'Dimension of Array : ', NDIMR
                WRITE (0,'(A21,I6)') 'Dimension in File  : ', INUM
C!!!                STOP 'ERROR IN STORAGE (ACTION: READ)'
              ENDIF
              NDIMR = MIN(NDIMR,INUM)
            ENDIF
C***        EXIT THIS LOOP
            GOTO 45
          ENDIF
        ENDDO
   45   CONTINUE
        IF (.NOT. BVARKN) THEN
          IF (IWARN .GE. 1) THEN
            WRITE (0,'(A,A8,A)') 
     >        ' WARNING from STORAGE: Variable ', NAME, 'not found'
          ENDIF
          IERR = -10
          RETURN
C!!!          STOP 'ERROR in STORAGE (ACTION: READ)'
        ELSE
          IERR = 0
          IFIRST = IADR(2+INDEX)
          INUM   = IADR(4+INDEX)
        ENDIF

C!!!  ALTE VERSION (MIT UMKOPIEREN DES GESAMTEN ARRAYS)
C***    AS LONG AS THE ARRAY FILLS THE NEXT RECORD COPY IT TO SCRATCH
C***    THE LAST RECORD IS, IF NECESSARY, FILLED UP WITH 0.
          NIND3 = (NDIMR-1) / IADRL + 1
          NREST = NDIMR
          DO I=1, NIND3
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >            SCRATCH
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)
              X(J+INDEX) = SCRATCH(J)
            ENDDO
            NREST = NREST - IADRL
          ENDDO

C!!!          NIND3 = (NDIM-1) / IADRL + 1
C!!!          DO I=1, NIND3-1
C!!!            INDEX = (I-1) * IADRL
C!!!            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
C!!!     >            X(INDEX+1)
C!!!          ENDDO
C!!!          INDEX = (NIND3-1) * IADRL
C!!!          READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
C!!!     >          SCRATCH
C!!!          NREST = NDIM - ((NIND3-1) * IADRL)
C!!!          DO J=1, NREST
C!!!            X(J+INDEX) = SCRATCH(J)
C!!!          ENDDO


C***  New option, returns in NDIM the array lenth of the variable
C***  wrh 10-Aug-2007 
C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'LENGTH') THEN
C***                             =====
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
          IF (IADR(1+INDEX) .EQ. NAME) THEN 
             NDIM = IADR(4+INDEX) 
             GOTO 14
          ENDIF
        ENDDO
        WRITE (0,*) 'ERROR: COULD NOT FIND LENGTH OF MS-VARIABLE ', 
     >              CNAME
        STOP 'FATAL ERROR in subr. STORAGE (action: LENGTH)'
   14   CONTINUE
C***                             =====
      ELSE IF (ACTION( :5) .EQ. 'CLOSE') THEN
C***                             =====

C***  STORE IADR(1..6)
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR
C***  WRITE INDEX RECORD TO FILE
        DO I=1, NIND2U
          DO J=1, IADRL
            ISCRATCH(J) = IADR((I-1)*IADRL + J)
          ENDDO
          WRITE (UNIT=ICHANNEL, REC=I, ERR=92, IOSTAT=IOS) ISCRATCH
        ENDDO
        CLOSE(ICHANNEL)

C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'SCLOSE') THEN
C***                             =====
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR

C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'CHANGE') THEN
C***                             =====
        BVARKN = .FALSE.
        DO I=1, IADR(6)
          INDEX = 10 + (I-1)*5
c          write (0,'(a,i8,a,a8,a,a8)') 'index=',index, 
c     >                            'nameold=',IADR(1+INDEX),
c     >                            'testname=',cname
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            BVARKN = .TRUE.
            GOTO 10
          ENDIF
        ENDDO
   10   CONTINUE
        IF (BVARKN) THEN
          IADR(INDEX+1) = NAME2
        ELSE
          WRITE (0,*) 'WARNING: CHANGE: Variable not found'
c          WRITE (0,'(A9,A8,1x,a8)') 'VAR1,2',NAME,NAME2
        ENDIF

C***                             =====
      ELSE IF (ACTION( :4) .EQ. 'INFO' .AND. 
     >         ACTION( :6) .NE. 'INFO-D') THEN
C***                             ====
C        write (*,*) 'test---------------'
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR
        WRITE(*,'(I9,4X,A)') IRECL2, ': RECORD LENGTH OF THE FILE'
        WRITE(*,'(I9,4X,A)') NIND2,  ': NUMBER OF INDEX RECORDS'
        WRITE(*,'(I9,4X,A)') NIND2U, ': NUMBER OF INDEX RECORDS USED'
        WRITE(*,'(I9,1X,A)') NINDEX, ': NUMBER OF INDICES IN INDEX REC'
        WRITE(*,'(I9,1X,A)') NREC,   ': TOTAL NUMBER OF RECORDS'
        WRITE(*,'(I9,4X,A)') NVAR,   ': NUMBER OF VARIABLES'

        IF (ACTION( :6) .EQ. 'INFO-L') THEN
C***                          ------
          DO I=1, NVAR
            INDEX = 10 + (I-1)*5
            IF (IADR(5+INDEX) == -1) THEN
              KSTR = 'default'
            ELSE
              WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)
            ENDIF
            WRITE (*,50) I, IADR(1+INDEX), 
     >                   IADR(2+INDEX), 
     >                   IADR(3+INDEX), IADR(4+INDEX), KSTR
   50       FORMAT (' * Variable Nr.', I6, 2X, 'Name=', A, 2X, 
     >              'First=', I6, 2X, 'Num-Rec=', I6, 2X, 
     >              'Num-Var=', I8, 2X, 'Vartype=', A)

          ENDDO
        ENDIF
      ELSE IF (ACTION( :6) .EQ. 'INFO-D') THEN
C***                             ------
C***  CHECK IF VARIABLE IS KNOWN
        BVARKN = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
          ID = IDX(CNAME)
c        write (0,'(4(a,1x))') 'INFO-D :', IADR(1+INDEX), NAME, ':'
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            BVARKN = .TRUE.
C***        EXIT THIS LOOP
            GOTO 44
          ENDIF
        ENDDO
   44   CONTINUE
        IF (.NOT. BVARKN) THEN
          WRITE (0,'(3A)') 'Variable ',CNAME, ' not known'
          STOP 'ERROR WHEN ACTION = INFO-D'
        ENDIF
        IERR = 0
        IFIRST = IADR(2+INDEX)
        INUM   = IADR(4+INDEX)

C***    Kurzuebersicht ueber die Variable mit der Nummer I
        INDEX = 10 + (I-1)*5
        IF (IADR(5+INDEX) == -1) THEN
          KSTR = 'default'
          CKIND = ' '
        ELSE
          WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)
          SELECTCASE(KSTR(1:1))
            CASE ('i','I') 
              CKIND = 'I'
            CASE ('r','R') 
              CKIND = 'R'
            CASE ('a','A')
              CKIND = 'A'
            CASE ('c','C') 
              CKIND = 'C'
            CASE DEFAULT
              CKIND = ' '
          ENDSELECT
          IF (CKIND /= ' ') THEN
            KINDSTR = ''
            DO J=2, 8
              SELECTCASE(KSTR(J:J))
                CASE ('1':'9', '0')
                  KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                CASE DEFAULT
                  EXIT
              ENDSELECT
            ENDDO
            READ(UNIT=KINDSTR, FMT='(I8)') IKIND       !transform number part from kind string into integer
          ENDIF
        ENDIF
        WRITE (*,50) I, IADR(1+INDEX), 
     >               IADR(2+INDEX), 
     >               IADR(3+INDEX), IADR(4+INDEX), KSTR

C***    Ausgabe der Variable mit der Nummer I
        NIND3 = (INUM-1) / IADRL + 1
        NREST = INUM
        LFMT = IDX(CNAME2)
        bWRINT = .FALSE.
        IF ((CNAME2(1:4) == 'AUTO') .OR. (CNAME2(1:4) == 'auto')) THEN
          IF (CKIND == 'I') THEN
            IFMTLEN = INT(  LOG10(2. ** FLOAT(IKIND * 8)) + 2. )
            WRITE(UNIT=FMTSTR, FMT='(I7,A1)') IFMTLEN, ')'
            FMTSTR = '(I' // TRIM(ADJUSTL(FMTSTR))
          ELSEIF (CKIND == 'R') THEN
            FMTSTR = '(G16.6)'
          ELSEIF ((CKIND == 'A') .OR. (CKIND == 'C')) THEN
            FMTSTR = '(A)'
          ELSE
            WRITE (*,*) '* WARNING: No vartype specified,',
     >                    ' auto-format is not recommended'
            SELECTCASE (CNAME(1:1)) !guessing based on first letter (implicit fortran)
              CASE ('I':'N','i':'n')
                FMTSTR = '(I16)'
              CASE DEFAULT
                FMTSTR = '(G16.6)'
            ENDSELECT
          ENDIF  
        ELSEIF (       (CNAME2(1:1) /= '(') 
     >           .AND. (CNAME2(LFTM:LFTM) /= ')')       ) THEN
          FMTSTR = '(' // TRIM(ADJUSTL(CNAME2)) // ')'
        ELSE
          FMTSTR = CNAME2
        ENDIF
        FMTSTR = ADJUSTL(FMTSTR)
        IF (FMTSTR(1:2) == '(I') bWRINT = .TRUE.
        WRITE (*,*) 'N=?'
        DO I=1, NIND3
          IF (CKIND /= ' ') THEN
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >              SCRATCHBYTE
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST), IKIND
              WRITE (*,FMT=TRIM(FMTSTR)) (SCRATCHBYTE(JJ), JJ=J, J+IKIND-1)
            ENDDO
          ELSE
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >              SCRATCH
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)    
              IF (bWRINT) THEN
                INTSCRATCH = TRANSFER(SCRATCH(J), INTSCRATCH)
                WRITE (*,FMT=TRIM(FMTSTR)) INTSCRATCH
              ELSE 
                WRITE (*,FMT=TRIM(FMTSTR)) SCRATCH(J)
              ENDIF               
            ENDDO
          ENDIF
          NREST = NREST - IADRL
        ENDDO
        WRITE (*,*) 'FINISH'

      ELSE
        WRITE (0,*) ' ACTION ', ACTION( :IDX(ACTION)), ' NOT KNOWN'
        STOP 'ERROR IN STORAGE'

      ENDIF

      RETURN

C********** Error Stops ****************************************

   90 WRITE (0,*) ' ERROR WHEN OPENING MASS-STORAGE FILE'
      GOTO 99

   91 WRITE (0,*) ' ERROR WHEN READING MASS-STORAGE FILE (LABEL=91)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   92 WRITE (0,*) ' ERROR WHEN WRITING MASS-STORAGE FILE (LABEL=92)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   93 WRITE (0,*) ' ERROR WHEN WRITING MASS-STORAGE FILE (LABEL=93)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   94 WRITE (0,*) ' ERROR WHEN READING MASS-STORAGE FILE (LABEL=94)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99


99    WRITE (0,'(A,I4)') 'Fortran Channel:', ICHANNEL
      WRITE (0,'(A,I4)') 'IOS=',IOS
      STOP 'ERROR IN SUBROUTINE STORAGE'

      END
      SUBROUTINE TABREAD (ND,RADIUS,VELO,GRADI,T,ITAB)
C***********************************************************************
C***  READS TEMPERATURE AND/OR VELOCITY STRUCTURE FROM FILE TAPE8=TABLE ********
C***  ITAB = 0: NO TABLE INPUT (DEFAULT)
C***  ITAB = 1: INPUT OF TABULATED TEMPERATURE STRUCTURE T(R)
C***  ITAB = 2: INPUT OF TABULATED VELOCITY FIELD V(R)
C***  ITAB = 3: INPUT OF T(R) AND V(R)
C***********************************************************************
 
      REAL RADIUS(ND),VELO(ND),GRADI(ND),T(ND)
      REAL R(70),V(70),G(70),TH(70),IP(70)
      CHARACTER KARTE*80
 
      DO 10 L=1,70
   10 IP(L)=0
 
      OPEN (8,FILE='TABLE', STATUS='UNKNOWN')
    1 READ(8,2, END=99) KARTE
    2 FORMAT(A)

      IF (KARTE(:1) .EQ. '*' ) GOTO 1
      IF (KARTE(:5).EQ.'TABLE') THEN
         DECODE (80,11,KARTE) NRP,VTAB,TTAB
   11    FORMAT (15X,I3,7X,F7.1,7X,F8.0)
         IF (NRP.LE.0) THEN
            CALL REMARK ('NO TABLE INPUT')
            STOP 'ERROR'
            ENDIF
         IF (NRP.GT.70) THEN
            CALL REMARK ('NR TAB PTS .GT. 70')
            STOP 'ERROR'
            ENDIF
         GOTO 1
      ENDIF
      IF (KARTE(:7).EQ.'ENDGRID') GOTO 30
 
C***  DECODE INPUT
      DECODE (80,12,KARTE) L,R(L),V(L),G(L),TH(L)
   12 FORMAT (I5,4F10.5)
      IP(L)=1
      GOTO 1
   30 CONTINUE
      CLOSE (8)
C***  END INPUT
 
      PRINT 31, VTAB,TTAB
   31 FORMAT (1H1,//,20X,' TABULAR INPUT OF THE VELOCITY AND/OR TEMPERATURE
     $ URE STRUCTURE',/,20X,63(1H=),//,10X,
     $ 'VELOCITY UNIT =',F8.1,10X,'TEMPERATURE UNIT=',F9.0,//,10X,
     $                                    '  NO   RADIUS    VELOCITY   GRAD
     $ RADIENT   TEMPERATURE',//)
      DO 33 L=1,NRP
      PRINT 32, L,R(L),V(L),G(L),TH(L)
   32 FORMAT (10X,I4,4F11.5)
      IF (IP(L).NE.1) THEN
         CALL REMARK ('GRID POINT MISSING')
         STOP 'ERROR'
         ENDIF
   33 CONTINUE
 
C***  INTERPOLATION OF THE TABLE GRID
C***  TEMPERATURE STRUCTURE T(R)
      IF (TTAB.GT.0.0) THEN
      ITAB=1
      DO 51 L=1,ND
      REF=RADIUS(L)
      CALL LIPO (F,REF,TH,R,NRP)
   51 T(L)=F*TTAB
      END IF
C***  VELOCITY FIELD V(R)
      IF (VTAB.GT.0.0) THEN
      IF (ITAB .GT. 0) THEN
         ITAB=3
      ELSE
         ITAB=2
      ENDIF
      DO 41 L=1,ND
      REF=RADIUS(L)
      CALL LIPO (F,REF,V,R,NRP)
      VELO(L)=F*VTAB
      CALL LIPO (F,REF,G,R,NRP)
      GRADI(L)=F*VTAB
   41 CONTINUE
      END IF
 
      RETURN

C***  ERROR EXIT:
   99 CONTINUE
      CALL REMARK ('NO ENDGRID FOUND')
      STOP 'ERROR'

      END
      SUBROUTINE TAUSCAL (RSTAR,ND,RADIUS,RNE,
     $                   ENTOT,T,POPNUM,NDIM,N,EN,LEVEL,NCHARG,WEIGHT,
     $                   ELEVEL,EION,EINST,ALPHA,SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT,NOM,NF,
     >                   XLAMBDA,FWEIGHT,TAUTHOM,TAUROSS,MAXATOM,MAXION,
     >                    SIGMATHK,SEXPOK,EDGEK,KODAT,KONTNUP,KONTLOW,
     >                    LASTKON, DENSCON, FILLFAC, POPMIN)
C***********************************************************************
C***  CALCULATION OF THE NLTE OPTICAL DEPTH SCALES (ROSSELAND, THOMSON)
C***   for continuum opacities
C***
C***   called from 
C***     STEAL -> HYDROSOLVE -> TAUSCAL (3x)
C***     STEAL -> HYDROSOLVE -> HDSOLUTION -> TAUSCAL 
C***     STEAL -> HYDROSOLVE -> HYDROVELO -> CALCTAUCONTMAX -> TAUSCAL 
C***     STEAL -> ENSURETAUMAX -> TAUSCAL (2x)
C***     STEAL -> TAUSCAL (2x)
C***     WRSTART -> TAUSCAL
C***  
C***  note: this routine updates the EN vector
C***********************************************************************
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NDIM, N, MAXATOM, MAXION, LASTKON, NF
      REAL, INTENT(IN) :: RSTAR

      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, RNE, ENTOT, T
      REAL, DIMENSION(ND), INTENT(INOUT) :: TAUTHOM, TAUROSS
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPNUM
      REAL, DIMENSION(N), INTENT(INOUT) :: EN
      
      REAL, DIMENSION(N), INTENT(IN) ::  WEIGHT, ELEVEL, EION
      REAL, DIMENSION(N, N), INTENT(IN) :: EINST
      REAL, DIMENSION(LASTKON), INTENT(IN) :: 
     >                ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM, MAXION), INTENT(IN) :: 
     >                SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(NF), INTENT(IN) :: FWEIGHT, XLAMBDA
      CHARACTER(10), DIMENSION(N), INTENT(IN) :: LEVEL(NDIM)

      INTEGER, DIMENSION(N), INTENT(IN) :: NCHARG, NOM
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: 
     >                IGAUNT, KONTNUP, KONTLOW
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: KODAT

C*** tiefenabh. clumping nach goetz
      REAL, DIMENSION(ND), INTENT(IN) :: DENSCON, FILLFAC
 
      INTEGER :: L, I
      REAL :: TL, RL, ENTOTDENSL, RNEL, OPAMEAN, DR, RM1, 
     >        ENELM1, ENEL, OPARM1, OPARL, ENEMEAN, POPMIN
      LOGICAL :: bUseOPAROSSpre

C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      REAL, PARAMETER :: SIGMAE = 0.6652E-24
 
C***  INITIALIZATION
      TAUTHOM(1) = 0.0
      TAUROSS(1) = 0.0

C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO L=1,ND
        RL=RADIUS(L)
        TL=T(L)
C***    Calculate opacities for clump density
        ENTOTDENSL=ENTOT(L) * DENSCON(L)
        RNEL=RNE(L)
        ENEL=RNEL*ENTOTDENSL
        DO I=1,N
          EN(I)=POPNUM(L,I)
        ENDDO
        !note that OPAROSS only calculates continuum opacities
        !(this is used in steal for taumax fixing)
        CALL OPAROSS (OPARL,EN,TL,RNEL,ENTOTDENSL,RSTAR,NDIM,N,
     >                  LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     >                  ALPHA,SEXPO,
     >                  ADDCON1, ADDCON2, ADDCON3, 
     >                  IGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     >                  MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >                  RL,KONTNUP,KONTLOW,LASTKON,POPMIN)
ccc        WRITE(0,*) 'L, OPARL ', L, OPARL
        IF (OPARL <= 0.) THEN
          WRITE(0,*) " WARNING: negative Opacity ", OPARL
C          WRITE(0,*) '   taking absolute value... '
C          OPARL = ABS(OPARL)
C***       Absolute of negative value might still be desasterous
C***       -> use Thomson only (recalculated here)
           OPARL = SIGMAE * ENEL
        ENDIF
        OPARL = OPARL * FILLFAC(L) !downscale opacity with filling factor
        IF (L > 1) THEN
          DR = RM1 - RL
          ENEMEAN    = 0.5 * (ENEL+ENELM1)
          TAUTHOM(L) = TAUTHOM(L-1) + 
     >                  FILLFAC(L) * SIGMAE * ENEMEAN * DR * RSTAR
          OPAMEAN    = 0.5*(OPARL+OPARM1) 
          TAUROSS(L) = TAUROSS(L-1)+OPAMEAN*DR
        ENDIF
        RM1=RL
        ENELM1=ENEL
        OPARM1=OPARL
      ENDDO
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
      FUNCTION TRADFUN (XLAMBDA,XJ)
C***********************************************************************
C***  RADIATION TEMPERAURE IN KELVIN, FROM XJ = J-NUE (CGS UNITS)
C***  AND XLAMBDA IN ANGSTROEM
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
C***  Version improved to prevent overflow for almost-zero XJ
C***        wrh 14-Mar-2005 17:15:26
C***********************************************************************

      DATA C1,C2 / 1.4388E8, 3.9724E+8 /
C***  Threshold for Taylor expansion of ALOG(1+X)
      DATA EPS / 1.E-14 /

C***  Zero Trad for negative J
cc      IF (XJ .LE. .0) THEN 
cc wegen Absturz:
      IF (XJ .LE. 1.E-100) THEN 
         TRADFUN=.0
         RETURN
      ENDIF

      W=1./XLAMBDA
      W3 = W * W * W

      IF (C2 * W3 / EPS .LT. XJ) THEN 
C***  Equivalent to:      IF (X .LT. 1.E-14) THEN
C***     SERIES EXPANSION: LN(1+X) = X
         TRADFUN = C1 * XLAMBDA * XLAMBDA * XJ / C2
      ELSE
         X = C2 * W3 / XJ
         TRADFUN=C1*W/ALOG(1.+X)
      ENDIF

      RETURN
      END
      SUBROUTINE TRBK

      WRITE (*,*) 'TRACE BACK FACILITY (TRBK) IS NOT AVAILABLE AT',
     >            ' DEC/OSF1 AT THE MOMENT'

      CALL ABORT()
      RETURN
      END
      FUNCTION VELOBETA(R)
C***********************************************************************
C***  VELOCITY FIELD from analytic (beta) law
C***  BETA < 0: different fine parametrization of beta law with abs(beta)
C***  BETA = 0: SWITCH TO SQRT-LOG-LAW
C***
C***  called by WRVEL, DELTAGRTHIN, VELTHIN
C***********************************************************************

      IMPLICIT NONE

      REAL :: VELOBETA
      REAL, INTENT(IN) :: R
      
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, RM 
      LOGICAL :: bSMOCO

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >               BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

        IF (BETA == 0.) THEN
          VELOBETA = VPAR1 * SQRT(LOG(R/VPAR2))
        ELSEIF (BETA < 0.) THEN
          VELOBETA = VPAR1 * (1.-VPAR2/R)**(ABS(BETA))
          IF (BETA2FRACTION > .0) THEN 
             VELOBETA = VELOBETA + 
     >                VPAR1_2 * (1.-VPAR2_2/R)**(ABS(BETA2))
          ENDIF
        ELSE
C***     New parametrization of beta law, wrh  5-Mar-2001 
          RM = (VPAR2+R)
          IF (RM <= 1.) RM = 1. + 1.E-10
          VELOBETA = VPAR1 * (1.-1./RM)**BETA
          IF (BETA2FRACTION > .0) THEN 
             VELOBETA = VELOBETA + 
     >                VPAR1_2 * (1.-1./(VPAR2_2+R))**BETA2   
          ENDIF
        ENDIF

      RETURN
      
      END
      SUBROUTINE VELTHIN (T,       !Electron Temperature per depth point
     >                    RADIUS,  !Radius (in units of Rstar) per depth point
     >                    VELO,    !Velocity (in km/s) per depth point (modified here)
     >                    ND,      !Total number of depth points
     >                    RSTAR,   !Rstar in Rsun
     >                    RMAX,    !Rmax in Rstar 
     >                    GEFFL,   !g_eff per depth point
     >                    XMU,     !mean particle mass per depth point
     >                    VTURB,   !turbulence velocity (hydrostatic part only) in km/s
     >                    ThinCard,    !CARDS line with HYDROSTATIC INTEGRATION (for parameters)
     >                    CRITERION    !return velocity criterion
     >                   )
C*******************************************************************************
C***  INITIALIZATION OF THE VELOCITY-FIELD VECTOR VELO(L) 
C***  IN THE INNER PART, THE CORRECT HYDROSTATIC EQUATION IS INTEGRATED,
C***  ACCOUNTING FOR THE DEPTH-DEPENDENT SCALE HEIGHT ACCORDING TO THE
C***  GIVEN TEMPERATURE STRATIFICATION AND THE R**2 - SCALING OF GRAVITY
c***  MAJOR REVISION by Andreas Sander (Feb. 2012)
C***  detailed documentation of the first part is given in Sander et al. (2015)
C*******************************************************************************
  
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: RSTAR, RMAX, VTURB
      !Note: XMU is now an array, earlier versions used only XMASS = XMU(ND)
      REAL, DIMENSION(ND), INTENT(IN) :: T, RADIUS, XMU, GEFFL
      !Do not use INTENT(OUT) for VELO because outer part is read later on
      REAL, DIMENSION(ND), INTENT(INOUT) :: VELO     
      CHARACTER(8), DIMENSION(ND), INTENT(OUT) :: CRITERION
      CHARACTER(80), INTENT(IN) :: ThinCard

      REAL, DIMENSION(4) ::  RIP, VIP
      CHARACTER(40), DIMENSION(20) :: CURPAR

      INTEGER, PARAMETER :: NDTHINMAX = 100.            
      REAL, DIMENSION(NDTHINMAX) :: VBETA, VHELP, RHELP, 
     >                              A2SUM, VELOMINUSSONIC

      INTEGER :: I, J, IMIN, ICON, ICONOLD, K, L, MAXIT, IMAXBETA,
     >           NONMONO, NEXTRA, NSameRcon, LOszil, IE, NPAR, IERR, 
     >           LL, Lguess, LEnforcerNeeded, INMS, INME, NDHELP
      REAL :: DGP, GRAD1, GRAD2, BPL, DR, H, RM, TMID, Q, HCONST, Vcon,
     >        VA, VMAX, S, ST, GRADIN, Vdummy, VINT, VCUR, VFINAL, DXX,
     >        VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE, FPLAST, FPMAX, H0,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, RCONOLD, ROUTMAX,
     >        DVDRcon, VDIFF, RCONMIN, RCONMAX, RGRADMAX, BAD, VADD,
     >        RFA, RFB, RFX, RFFA, RFFB, RFFX, RFD, EXPO, RSTEP, RLAST,
     >        RINT, RIN, XMUCUR, Hlast, RBETAMIN, RCONSTART, P, VAM,
     >        GEFF, GEFFR, A2SUML, A2, DA2DR, BETACON, VOFF, RONSET,
     >        RCONMAXguess, RCONMAXnow, tempREAL, fsonic, DVDR, VSOUND,
     >        VMINUSA, VMINUSAM, DELTAB, VL,
     >        GRADRCONMIN, GRADRCONMAX
      LOGICAL :: bConFound, bConPointConv, bRFbi, bSmoothVelo, bDEBUG,
     >           bUseSonicConnectionCriterion, bForceMonotonic, bRUKU,
     >           bFailsafe, bPrintVelo, bFullHD, bSMOCO, bTestMono, 
     >           bNewSonic

      REAL, EXTERNAL :: DELTAGRTHIN, VELOBETA

      !Common block VELPAR needed to transfer velocity parameters
      ! This routine updates VPAR1, VPAR2, RCON
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                 BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      !Numerical constants
      INTEGER, PARAMETER :: MAXCONIT = 99    !Maximum iterations for connection point search (was 20)
      REAL, PARAMETER :: EXPMAX = 500        !Maximum exponent for velocity scaling (prevent overflow)
      REAL, PARAMETER :: RCONACC = 1.E-6     !Accuracy for RCON
      REAL, PARAMETER :: RFACC = 1.E-12      !Accuracy for Regula Falsi Call
      REAL, PARAMETER :: EPSRCONMIN = 1.E-12 !Offset from beta law singularity

      REAL, DIMENSION(MAXCONIT) :: RCONHIST !History with calculated RCON values


      !Physical Constants
      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: AMU = 1.66E-24     !ATOMIC MASS UNIT (GRAMM)
      REAL, PARAMETER :: RGAS = BOLTZ / AMU !Gas Constant (CGS) = BOLTZ / AMU

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      !Defaults for Switches
      bUseSonicConnectionCriterion = .FALSE.
      bForceMonotonic = .TRUE.
      bSmoothVelo = .FALSE.
      bFailsafe = .FALSE.
      LEnforcerNeeded = 0
      bPrintVelo = .FALSE.        !Debug option: Print resulting VELO vector
      bFullHD = .FALSE.
      bRUKU = .TRUE.
      fsonic = 1.

      IF (ND > NDTHINMAX) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'VELTHIN: NDTHINMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'ND = ', ND,
     >                         ', NDTHINMAX = ', NDTHINMAX
        STOP 'FATAL ERROR IN VELTHIN'      
      ENDIF
          
c      IF (BETA < .0) THEN
c        CALL REMARK ('VELTHIN: BETA OPTION INVALID')
c        STOP 'ERROR'
c      ENDIF
      
      !Read additional CARDS line parameters (if existing)
      CALL SARGC (ThinCard, NPAR)
      DO i=1, NPAR
        CALL SARGV(ThinCard,i,CURPAR(I))
      ENDDO
      IF (NPAR >= 2) THEN
        DO i=2, NPAR 
          SELECTCASE(CURPAR(i))
            CASE ('NONMONO', 'NONMONOTONIC') 
              bForceMonotonic = .FALSE.
            CASE ('SMOOTH')
              bSmoothVelo = .TRUE.
            CASE ('SONIC')
              bUseSonicConnectionCriterion = .TRUE.
              IF (NPAR >= (i+1)) THEN
                READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                IF (IERR == 0) THEN
                  fsonic = tempREAL
                ENDIF                          
              ENDIF
            CASE ('RUKU')
              bRUKU = .TRUE.
            CASE ('NORUKU')
              bRUKU = .FALSE.
            CASE ('FULLHD')
              bFullHD = .TRUE.
            CASE ('PRINTV')
              bPrintVelo = .TRUE.
          ENDSELECT
        ENDDO
      ENDIF
       
  20  IF (bFailsafe) THEN
        !failsafe branch (only reached via jump)
        bUseSonicConnectionCriterion = .TRUE. 
      ENDIF
      
     
C***  LOWER PART OF VELOCITY LAW: 
C***  INTEGRATION OF THE HYDROSTATIC EQUATION

C***  ANSATZ    p = p0 f(r) exp(-(r-1)/HSCALE)
C***            rho = AMU 
C***   
C***            f(r) is obtained numerically by integrating (f = FP in the code)
C***            df/dr = f (1/HSCALE - 1/H)

C***  FIRST STEP: INTEGRATE   dp/dr = -1/H p
C***  WITH DEPTH-DEPENDEND SCALE HEIGHT H = HCONST * T(I)
C***  1. STEP: SOLVE DIFF. EQ. FOR THE CORRECTION FACTOR TO THE
C***  EXPONENTIAL LAW WITH CONSTANT SCALE HEIGHT HSCALE
C***  - THIS FACTOR IS STORED IN VECTOR VELO

      DO L=1, ND
        !initialize VELO vector with VFINAL (just in case)
        VELO(L) = VFINAL
      ENDDO
      
      H0 = HSCALE 

      IMIN = 1
      BPL = 0.
      VELO(ND) = VMIN

C***  Calculate isothermal sound speed a (squared) 
C        and add turbulence contribution
      DO L=1, ND
        A2  = RGAS * T(L) / XMU(L)
        A2SUM(L) = A2 + VTURB*VTURB*1.E10        
      ENDDO
      
C***  Initialize the loop at l=ND
      GEFFR = GEFFL(ND) / (RADIUS(ND) * RADIUS(ND))
      H     = A2SUM(ND) / GEFFR / RSTAR

      hydrostat: DO L=ND-1,1,-1        
        DR = RADIUS(L) - RADIUS(L+1)         


        GEFFR = GEFFL(L) / (RADIUS(L) * RADIUS(L))
        
        IF (bFullHD) THEN
C***      use full hydrodynamic formula  
          IF (bRUKU) THEN
c            WRITE (hCPR,'(A,3(2X,G15.8))') 'A/V= ', 
c     >        SQRT(A2SUM(L+1))/VELO(L+1)/1.E5, 
c     >        SQRT(A2SUM(L+1))/1.E5, VELO(L+1)
            IF (SQRT(A2SUM(L+1))/VELO(L+1)/1.E5 - 1. < 1.) THEN
              RCONMAX = RADIUS(L+1)
              IMIN = L+1
              EXIT hydrostat
            ELSE
              CALL HYSTHDRUKU(RADIUS(L+1), RADIUS(L), VELO(L+1), VL, 
     >                        RADIUS, GEFFL, A2SUM, ND, RSTAR)
              VELO(L) = VL
            ENDIF
          ELSE 
            CALL SPLINPOX(A2SUML, RADIUS(L), A2SUM, RADIUS, ND, 
     >                                               DFDX=DA2DR)
            DVDR = (GEFFR - 2. * A2SUML/RADIUS(L) + DA2DR) /
     >             ( A2SUML/VELO(L+1) - VELO(L+1) )
            IF (DVDR < 0.) THEN
              RCONMAX = RADIUS(L+1)
              IMIN = L+1
              EXIT hydrostat
            ENDIF
            VELO(L) = DVDR * DR
          ENDIF
        ELSE
C***      uses hydrostatic formula with exp split        
          Hlast = H        ! from last loop, i.e. at L+1
          
          H     = A2SUM(L) / GEFFR / RSTAR

          H = 0.5 * (H + Hlast)
          IF (bRUKU) THEN
            CALL HYSTRUKU (RADIUS(L+1), RADIUS(L), DELTAB, 
     >                   RADIUS, GEFFL, A2SUM, ND, H0, RSTAR)
            BPL = BPL + DELTAB
c          WRITE (0,*) 'BLA: ', DR * (1./H0 - 1./H), DELTAB
          ELSE
            BPL = BPL + DR * (1./H0 - 1./H)
          ENDIF
          
          EXPO = (RADIUS(L)-1.)/H0 - BPL

C***    prevent overflow because of exponential growth of VELO
          IF (EXPO .LT. EXPMAX) THEN
            VELO(L) = VMIN * A2SUM(L)/A2SUM(ND) * EXP(EXPO)/RADIUS(L)**2
            IMIN = L 
            RCONMAX = RADIUS(IMIN)
C*         end the static part one point after VFINAL
            IF (VELO(L+1) >= VFINAL) THEN
              VELO(L) = MIN(10.*VFINAL, VELO(L))
              EXIT hydrostat 
            ENDIF
          ELSE
            VELO(L) = VMIN*A2SUM(L)/A2SUM(ND) * EXP(EXPMAX)/RADIUS(L)**2
            VELO(L) = MIN (VELO(L), VFINAL)
            EXIT hydrostat
          ENDIF
        ENDIF
      ENDDO hydrostat

C***  Define the index of the connection point 
C***  (only if defined by the sonic speed, other wise this will happen later)
C***  - find index interval (ICON-1, ICON) where VELO first time exceeds f*VSOUND 
      IF (bUseSonicConnectionCriterion) THEN
        IF (bFullHD) THEN
C***      For FULLHD option take last good point of integration        
C***      and check if we can still go lower from there
          ICON = IMIN
          IF (ICON < ND) THEN
            VSOUND = SQRT(RGAS * T(ICON) / XMU(ICON)) * 1.E-5   !speed of sound in km/s
            VMINUSAM = VELO(ICON) - fsonic*VSOUND
            DO L=ICON+1, ND
              VSOUND = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5       !speed of sound in km/s
              VMINUSA = VELO(L) - fsonic*VSOUND
              IF (VMINUSA < 0. .OR. L == ND) THEN
                ICON = L
                EXIT
              ELSE 
                VMINUSAM = VMINUSA
              ENDIF
            ENDDO
          ENDIF
          IF (bDEBUG) THEN
            WRITE (0,*) 'DEBUG: ICON, IMIN = ', ICON, IMIN
          ENDIF          
        ELSE
          ICON = ND+1  ! No hydrostatic domain!
          DO L=ND, IMIN, -1
              IF (L .LT. ND) VMINUSA = VMINUSAM
              VSOUND = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5       !speed of sound in km/s
              VMINUSAM = VELO(L) - fsonic*VSOUND
              IF (VMINUSAM > .0) EXIT
              ICON = L
          ENDDO
        ENDIF
      ENDIF

C***  Ensure monotonic increase of the velocity
C***  note: this applies only to the hydrostatic part, i.e. from index IMIN
      IF (.NOT. bForceMonotonic) GOTO 16
   15 Continue
C***  Copy all points to help vectors, omitting those which are non-monotonic      
      LL = 0
      DO L=IMIN, ND
        IF (L .EQ. ND) THEN
           bTestMono = .TRUE.
        ELSE
           bTestMono =  (VELO(L) > VELO(L+1)) 
        ENDIF
        IF (bTestMono) THEN
          LL = LL + 1
          VHELP(LL) = VELO(L)
          RHELP(LL) = RADIUS(L)
        ELSEIF (LEnforcerNeeded == 0) THEN
           LEnforcerNeeded = L
        ENDIF
      ENDDO
      NDHELP = LL

      IF (NDHELP .EQ. ND-IMIN+1) GOTO 16 ! all points are (now) monotonic
C***  If not, then replace all points by their interpolated value      
C***    (trivially, points which had been copied will not change)
      DO L=ND-1, IMIN, -1
        IF (RADIUS(L) < RHELP(1)) THEN
c          CALL SPLINPO (VELO(L), RADIUS(L), VHELP, RHELP, NDHELP)
          CALL SPLINPOX(VELO(L), RADIUS(L), VHELP, RHELP, NDHELP)
        ELSE 
C***       non-monotonic point(s) might include IMIN, thus shorteneing 
C***       the range in which interpolation is possible. 
C***       Linear extrapolation is applied for these point(s) 
           IF (NDHELP >= 2) THEN
             VELO(L) = VHELP(1) + (RADIUS(L)-RHELP(1)) * 
     >               (VHELP(2)-VHELP(1))/(RHELP(2)-RHELP(1))  
           ELSE ! extreme case that only one good point exists
             VELO(L) = VELO(L+1) * 1.001
           ENDIF
        ENDIF
      ENDDO
      GOTO 15 ! Check once more for monotony, just in case

   16 CONTINUE

      RCONMAX = RADIUS(IMIN+1)+0.9*(RADIUS(IMIN)-RADIUS(IMIN+1))

      DO I=ND, 1, -1
        RHELP(ND-I+1) = RADIUS(I)
        VHELP(ND-I+1) = VELO(I)
      ENDDO
      
      VMAX = VFINAL * (1.-BETA2FRACTION)

C***  ITERATIVE SEARCH FOR THE CONNECTION POINT INDEX ICON
C***   WHERE THE DERIVATIVES ARE CONTINUOUS
C***   - GIVEN THAT POINT, THE PARAMETERS OF THE ANALYTIC LAW ARE
C***     RE-ADJUSTED SUCH THAT THE VELOCITIES ARE CONTINUOUS.
C***  THE SEARCH PROCEEDS INWARDS.
      
      IF (bUseSonicConnectionCriterion) THEN
C***    Sonic point is between indices (ICON,ICON-1) 
        bConPointConv = .TRUE.
        MAXIT = 1
        IF (bFullHD) THEN
          VCON = VELO(ICON)
          RCON = RADIUS(ICON)
        ELSEIF (ICON .LE. ND) THEN
C*      linear interpolation for RCON
          P = VMINUSA / (VMINUSA - VMINUSAM) 
          RCON = P * RADIUS(ICON-1) + (1.-P) * RADIUS(ICON)
          VA  = VELO(ICON)   - VMINUSA
          VAM = VELO(ICON-1) - VMINUSAM
          VCON = VA + (RADIUS(ICON-1) - RCON) *
     >      (RADIUS(ICON-1)-RADIUS(ICON))/(VAM - VA) 

c          write (*,*) 'VCON, VA, ICON=', VCON, VA, ICON
cc          CALL SPLINPO (VCON, RCON, 
cc     >     VELO(ICON-1:ND), RADIUS(ICON-1:ND), ND-ICON+2) 
        ELSE
           RCON = 1.
           Vcon = VMIN
        ENDIF

      ELSE
        bConPointConv = .FALSE.
        RCON = 1.
        MAXIT = MAXCONIT
        Vcon = VMIN
      ENDIF
            
      RCONSTART = RCON
      RCONOLD = RCON
      NSameRcon = 0
c      WRITE (0,*) 'VELTHING: Searchin RCON...'

      conpointloop: DO K=1, MAXIT !------------------------------------------

        RCONOLD = RCON
        IF (bDEBUG) WRITE (hCPR,*) ' conpointloop IT, =', K

        IF (.NOT. bUseSonicConnectionCriterion) THEN
          !Obtain VCON via spline interpolation on integration result
          IF (RCON > RADIUS(IMIN)) THEN
            WRITE (hCPR,*) "WARNING RCON=", RCON, " > ", RADIUS(IMIN)
          ENDIF
          IF (RCON < RHELP(1)) THEN
            WRITE (hCPR,*) ' RCON limited to minimum value: ', 
     >                             RHELP(1)            
            WRITE (hCPR,*) ' RHELP-: Rcon=', Rcon,' Rhelp(1)=', RHELP(1)
            RCON = RHELP(1)
            Vcon = VMIN
          ELSEIF (RCON >= RHELP(ND-IMIN+1)) THEN
            WRITE (hCPR,*) ' RCON limited to maximum value: ', 
     >                             RHELP(ND-IMIN+1)
C            WRITE (hCPR,*) ' RHELP+: Rcon=',Rcon,' Rhelp(+)=', 
C     >                             RHELP(ND-IMIN)
            RCON = RHELP(ND-IMIN+1)
            Vcon = VHELP(ND-IMIN+1)
          ELSE
            CALL SPLINPOX(Vcon,RCON,VHELP,RHELP,ND-IMIN+1,DFDX=DVDRcon)
          ENDIF
          !Failsafe-Check for Vcon (in case of bad SPLINPOX result)
          IF (Vcon < VMIN) THEN
            Vcon = VMIN
          ENDIF
        ENDIF

C***    CALCULATION OF VELOCITY-FIELD PARAMETERS (ANALYTIC LAW)
        IF (bUseSonicConnectionCriterion .AND. bNewSonic 
     >                                      .AND. BETA /= 0.) THEN
C***      In this branch v(Rmax) = v_inf is not guaranteed
C***      Possibly add iteration loop to ensure v_inf again
C***      TODO: Extend this for 2-beta-laws
          CALL SPLINPOX(Vdummy, RCON, VELO, RADIUS, ND, DFDX=DVDRcon)
          WRITE (0,*) 'DEBUG:  DVDRcon = ', DVDRcon
          VPAR1 = VFINAL - Vcon
          DO
          IF (BETA > 0) THEN
            VOFF = Vcon - (DVDRcon/BETA/VPAR1)**(BETA/(BETA-1.))
            VPAR2 = 1./(1.- ((Vcon - VOFF)/(VPAR1 - VOFF))**(1./BETA) ) 
     >                                   - RCON
          ELSE
            VOFF = Vcon - (RCON*DVDRcon/BETA/VPAR1)**(BETA/(BETA-1.))
            VPAR2 = RCON*(1.- ((Vcon - VOFF)/(VPAR1 - VOFF))**(1./BETA))
          ENDIF          
C***        Iteration to ensure VMAX 
c            WRITE (0,*) 'DEBUG: VINF, WRVEL(RMAX)', VFINAL, VELOBETA(RMAX)
            IF (ABS(VFINAL - VELOBETA(RMAX)) < 1.) EXIT
            IF (VFINAL - VELOBETA(RMAX) > 0) THEN
              VPAR1 = 1.1 * VPAR1
            ELSE 
              VPAR1 = 0.9 * VPAR1
            ENDIF
          ENDDO
        ELSE
          CALL INITVELBETAPAR(RMAX, RCON, Vcon, .TRUE.)
        ENDIF
c        CALL INITVELBETAPAR(RMAX, RCON, Vcon, .FALSE.)
c        WRITE (0,*) 'VCON = ', VCON, VELOBETA(RCON)

        IF (bUseSonicConnectionCriterion) THEN
          !no connection point loop needed if sonic criterion is used
          IF (BDEBUG) 
     >       CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 1., 1.)
          EXIT conpointloop
        ENDIF

C***    Define minimum radius for the BETA law       
        IF (BETA < 0.) THEN
          RBETAMIN = VPAR2 + EPSRCONMIN
        ELSE
          RBETAMIN = 1.-VPAR2+EPSRCONMIN  !Beta law cannot be calculated inside of this point
        ENDIF
        
C***    Initialize RCON and ICON in the first run of this loop        
        IF (K <= 1) THEN
          RIN = MAX(RBETAMIN,1.)
          IF (DELTAGRTHIN(RIN,ND,RADIUS,VELO) < .0) THEN
C***        Inner gradient is already larger than outer gradient at innermost point:
C***        => no hydrostatic domain          
            RCON=RIN
            ICON = ND
            IF (BDEBUG) 
     >         CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 1., 1.)
            EXIT conpointloop
          ENDIF

          IF ((IMIN >= ND) .OR. RCONMAX > RMAX) THEN          
            RCON=RMAX
            WRITE (hCPR,'(A)') 
     >          'VELTHIN: NO BETA-LAW DOMAIN ENCOUNTERED'
            RETURN
          ENDIF
        ENDIF

C***    LOWER BOUNDARY FOR RCON         
        RCONMIN=MAX(RCONSTART, RBETAMIN) !Lower guess for RCON
        !Radius of maximum of gradient in beta laws
        RGRADMAX = 1. + (BETA-1.)/2. - VPAR2    

        CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, RCONMIN, RGRADMAX)

C***    Loop: INCREASE RCONMIN STEPWISE TO ENSURE DELTAGRTHIN(RCONMIN) .GE. 0.
C***    DELTAGRTHIN returns GRAD_BETA - GRAD_HYST
C***    In the outer wind this value should be below zero (i.e. GRAD_HYST > GRAD_BETA)
C***    However, for BETA > 1, there can also be an inner range where 
C***      GRAD_HYST > GRAD_BETA, which must be avoided to find the correct solution.
C***    We therefore increase RCONMIN until we get out of this region.
        GRADRCONMIN = DELTAGRTHIN(RCONMIN,ND,RADIUS,VELO)
        DO WHILE (GRADRCONMIN < .0) 
          DXX = 0.1 * HSCALE        
          RCONMIN = RCONMIN + DXX
          GRADRCONMIN = DELTAGRTHIN(RCONMIN,ND,RADIUS,VELO)
          IF (RCONMIN > RGRADMAX) THEN
            WRITE (hCPR,'(A)') '*** VELTHIN: NO CONNECTION POINT FOUND'
            WRITE (hCPR,'(A)') '*** using sonic connection criterion...'
            bFailsafe = .TRUE.
            GOTO 20            
C            STOP               '*** ERROR IN VELTHIN'
          ENDIF
        ENDDO

C***    RCONMIN has now been increased such that we avoid the second solution for 
C***    cases with BETA > 1
C***    RCONMAX denotes the maximum radius where we have a hydrostatic solution
C***    IF RCONMAX > RCONMIN we are fine, but for RCONMIN > RCONMAX we cannot find 
C***    a solution
        
        RCONMAXguess = MAX(1., RBETAMIN)           !RCONMAXguess must be initialized
        Lguess = ND
        conmaxmin: DO L=ND-1, 1, -1
          IF (RADIUS(L) > RCONMIN) THEN
            Lguess = L
            RCONMAXguess = RADIUS(L)
            EXIT conmaxmin
          ENDIF
        ENDDO conmaxmin


C***    Try to find a good range for RCONMAX    
C       DELTAGRTHIN = BETALAWGRADIENT - HYSTGRADIENT
C***    We loop outwards through the grid until we reach the end of the 
C***     tabulated (hydrostatic) velocity field (index IMIN)
        conmaxmax: DO L=Lguess, IMIN, -1
          RCONMAXguess = RADIUS(L)
          GRADRCONMAX = DELTAGRTHIN(RCONMAXguess,ND,RADIUS,VELO)
          IF (GRADRCONMAX < .0) THEN
C***        Hydrostatic gradient is (again) larger than BETA law gradient
C***        => we have reached the outer limit for RCON
            EXIT conmaxmax
          ENDIF
          IF (L == IMIN) THEN
C***        We have never reached a part where we can find a connection point          
            RCONMAXguess = RCONMAX
          ENDIF
        ENDDO conmaxmax
        GRADRCONMAX = DELTAGRTHIN(RCONMAXguess,ND,RADIUS,VELO)

        
C***    Update vgrad (debug) plot with final boundaries
        CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 
     >                 RCONMIN, RCONMAXguess)
        
C***    Check if RCON search boundaries exclude each other
C***      or if there is no sign change between both boundaries
        IF (RCONMIN > RCONMAXguess .OR. 
     >           (GRADRCONMAX > 0. .AND. GRADRCONMIN > 0.)) THEN
C***         => no hydrostatic domain        
          RCON = 1.
          ICON = ND
          EXIT conpointloop
        ELSEIF (GRADRCONMAX < 0. .AND. GRADRCONMIN < 0.) THEN
C***         => no wind domain        
          RCON = 1.1 * RMAX
          ICON = 1
CCC       TODO: Check if ICON = 0 would be possible and better
          EXIT conpointloop
        ENDIF                        


C***    Now find new RCON where DELTAGRTHIN is zero:        
        !------------- Regula Falsi -------------------------
 
        bRFbi=.TRUE.
        RFA=RCONMIN
        RFB=RCONMAXguess
        RFFA = DELTAGRTHIN(RFA, ND, RADIUS, VELO)
        RFFB = DELTAGRTHIN(RFB, ND, RADIUS, VELO)
        IF (RFFA == 0.) THEN
          RCON = RFA
        ELSEIF (RFFB == 0.) THEN
          RCON = RFB
        ELSEIF (RFFA*RFFB > 0) THEN
          WRITE (hCPR,*) '*** INVALID ARGUMENTS FOR REGULA FALSI ***'
          WRITE (hCPR,FMT=10) RFA, RFB
   10     FORMAT ('INTERVAL:     A  =', E15.5, 5X, '  B  =', E15.5)
          WRITE (hCPR,FMT=11) RFFA, RFFB
   11     FORMAT ('FUNCTION:   F(A) =', E15.5, 5X, 'F(B) =', E15.5)
          WRITE (hCPR,'(A)') '*** using sonic connection criterion...'
          bFailsafe = .TRUE.
          GOTO 20            
C          CALL TRBK
C          STOP 'VELTHIN: ERROR IN REGULA FALSI'
        ELSE
          rfloop: DO !-  -  -  -  -  -  -  -  -  -  -
            RFD = RFA - RFB
            IF (ABS(RFD) < RFACC) THEN
              RFX = 0.5 * (RFA + RFB)
              EXIT rfloop
            ENDIF
            bRFbi = .NOT. bRFbi
            IF (bRFBi) THEN
              RFX = RFA - RFFA*RFD/(RFFA-RFFB)
            ELSE
              RFX = 0.5 * (RFA + RFB)
            ENDIF
            RFFX = DELTAGRTHIN(RFX, ND, RADIUS, VELO)
            IF (RFFX == 0.) EXIT rfloop
            IF (RFFX * RFFA > 0.) THEN
              RFA = RFX
              RFFA = RFFX
            ELSE
              RFB = RFX
              RFFB = RFFX
            ENDIF
          ENDDO rfloop !-  -  -  -  -  -  -  -  -  -  -
          RCON = RFX
        ENDIF
          
        !-------------- end of Regula Falsi -----------------

CC        WRITE (0,*) ' RCON found = ', RCON
        IF (RCON > RCONMAX) THEN
          WRITE (0,*) 'Error: RCON larger than RCONMAX '
          STOP "*** FATAL ERROR IN VELTHIN"        
        ENDIF
        RCONHIST(K) = RCON
C        WRITE (hCPR,*) ' THIN-IT:', K, ' RCON=', RCON

        IF (ABS(RCON-RCONOLD) <= RCONACC) THEN
          bConPointConv = .TRUE.
          EXIT conpointloop
        ELSEIF (K >= 2) THEN
          !Check if solution is oscillating
          NSameRcon = 0
          DO I=K-1, 1, -1
            IF (ABS(RCON-RCONHIST(I)) <= RCONACC) THEN
              NSameRcon = NSameRcon + 1
              !Cancel if more than 10 oscillations
              IF (NSameRcon > 10) THEN
                DO J=K-1, 1, -1
                  IF (ABS(RCON-RCONHIST(J)) <= RCONACC) THEN
                    LOszil = K - J
                    EXIT
                  ENDIF
                ENDDO
                !RCON = (RCONOLD + RCON) / 2.
                RCON = 0.
                DO J=K, K-LOszil+1, -1
                  RCON = RCON + RCONHIST(J)
                ENDDO
                RCON = RCON / FLOAT(LOszil)
                bConPointConv = .TRUE.
                BAD = ABS( MAXVAL(RCONHIST(K-LOszil+1:K))
     >                   - MINVAL(RCONHIST(K-LOszil+1:K)) )
     >            / ABS( 0.5 * ( MAXVAL(RCONHIST(K-LOszil+1:K))
     >                   + MINVAL(RCONHIST(K-LOszil+1:K)) ) )
C                CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1337)
                WRITE (hCPR,'(A,A,I2,A,F12.8,A)')
     >              'VELTHIN: WARNING - OSCILLATING SOLUTION',
     >              ' (LENGTH: ', LOszil, '  BADNESS: ', BAD, ')'
C                WRITE (hCPR,'(A,F12.8)')
C     >              '  averaged RCON: ', RCON
                EXIT conpointloop
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ENDDO conpointloop !------------------------------------------

      IMAXBETA = 0
      
      IF (ICON /= ND) THEN
        IF (.NOT. bConPointConv) THEN
          WRITE (hCPR,'(A)') 'VELTHIN: CONVERGENCE PROBLEMS'
          WRITE (hCPR,*) ABS(RCON-RCONOLD), ' ACC:',RCONACC
          STOP 'ERROR'
        ENDIF

        IF (.NOT. bUseSonicConnectionCriterion) THEN
          iconloop: DO I=ND, 1, -1  
            IF (RADIUS(I) >= RCON) THEN
              ICON = I
              EXIT iconloop
            ENDIF
          ENDDO iconloop
        ENDIF

C        WRITE (hCPR,FMT='(A,I2)') " ICON=", ICON
        IMAXBETA = ICON-1
      ELSE
        !no hydrostatic domain: use beta-law only
        IMAXBETA = ND
      ENDIF
      
c      IF (BETA > 0.) THEN
        !Calculate beta law velocities for comparison
        DO I=1, ND
          RM = RADIUS(I)
          IF (RM <= RCON) THEN
            IMAXBETA = I-1
            EXIT
          ELSE
            VBETA(I) = VELOBETA(RM)
            IF (bDEBUG) THEN
               WRITE (hCPR,FMT='(A,I2,A,F12.6)') 
     >           "DEBUG: VELOBETA(",I,")=",VBETA(I)
            ENDIF
          ENDIF
        ENDDO
  
c      ENDIF


C***  DEFINITION OF VELO(I) - OVERWRITING THE OUTER PART WITH THE
C***      BETA OR SQRT(ALOG(R)) LAW
      CRITERION = "HYSTINT "
C      WRITE (hCPR,*) " IMAXBETA=", IMAXBETA
      DO I=1, IMAXBETA
        IF (ABS(BETA) <= .0) THEN
          VELO(I) = VPAR1*SQRT(ALOG(RADIUS(I)/VPAR2))
          CRITERION(I) = "SQRT    "
        ELSE
          VELO(I) = VBETA(I)
          IF (BETA2FRACTION > .0) THEN 
            CRITERION(I) = "2BETA   "
          ELSE
            CRITERION(I) = "BETA    "
          ENDIF
        ENDIF
      ENDDO

C***  Ensure monotonic field at connection point     
C***  (This can be an issue in case of oscillating solutions)
      IF (IMAXBETA > 1 .AND. IMAXBETA < ND 
     >      .AND. VELO(IMAXBETA) < VELO(IMAXBETA+1)) THEN
        INME = IMAXBETA+1        
        nmcloop: DO J=IMAXBETA-1, 1, -1
          IF (VELO(J) > VELO(IMAXBETA+1)) THEN
            INMS = J
            EXIT nmcloop
          ENDIF
        ENDDO nmcloop
        IF (INMS > 2) THEN
          RIP(1) = RADIUS(INME)
          RIP(2) = RADIUS(INMS)
          RIP(3) = RADIUS(INMS-1)
          RIP(4) = RADIUS(INMS-2)          
          VIP(1) = VELO(INME)
          VIP(2) = VELO(INMS)
          VIP(3) = VELO(INMS-1)
          VIP(4) = VELO(INMS-2)          
        ELSE
          WRITE (hCPR,*) 'NONMONOTONIC BETA-PART IN WIND VELOCITY'
          STOP 'FATAL ERROR IN VELTHIN'
        ENDIF
        DO J=INMS+1, INME, 1
          CALL SPLINPOX(VELO(J), RADIUS(J), VIP, RIP, 4)
        ENDDO        
      ENDIF      
      
C***  Smooth velocity field if enforced (re-uses VHELP)
      IF (bSmoothVelo) THEN
        DO L=2, ND-1
          VHELP(L) = 0.25 * 
     >      (LOG10(VELO(L-1)) + 2. * LOG10(VELO(L)) + LOG10(VELO(L+1)))        
        ENDDO
        DO L=2, ND-1
          VELO(L) = 10**VHELP(L)
        ENDDO
      ENDIF

      !Debug Output
      IF (bPrintVelo) THEN
        DO I=1, ND
          WRITE (hCPR,FMT='(A,G10.3,A,I2,A,F12.6,A,A)') 
     >      "R-1=", RADIUS(I)-1.,
     >      '  VELTHIN-VELO(',I,')=',VELO(I),'  ',CRITERION(I)
        ENDDO
        WRITE (hCPR,FMT='(A,F10.5)') ' *** VELTHIN: RCON = ', RCON
        IF (BETA2FRACTION > .0) THEN 
          WRITE (hCPR,FMT='(A,4(2X,G15.8))') 'BETA VPARs ', 
     >          VPAR1, VPAR2, VPAR1_2, VPAR2_2
        ELSE 
          WRITE (hCPR,FMT='(A,2(2X,G15.8))') 'BETA VPARs ', VPAR1, VPAR2
        ENDIF
      ENDIF
C      CALL SPLINPOX(Vcon,RCON,VELO,RADIUS,ND)
      IF (IMAXBETA < 1) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: NO BETA-LAW DOMAIN ENCOUNTERED'
      ENDIF
      IF (IMAXBETA >= ND-1) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: NO HYDROSTATIC DOMAIN ENCOUNTERED'
      ENDIF
      IF (bForceMonotonic .AND. LEnforcerNeeded > ICON) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: MONOTONIC STRATIFICATION ENFORCED'
      ENDIF

      RETURN
      END


      SUBROUTINE VTURB_SETUP (VTURB, VTURB_LINE, ND,
     >                        VELO, RADIUS, TAUROSS, T, XMU)
C***********************************************************************
C***  This routine fills the VTURB vector as specified in the CARDS file
C***
C***  This routine is work in progress
C***
C***  called from WRSTART (so far only, planned: ENSURETAUMAX, HYDROSOLVE)
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, INTENT(IN) :: ND
      
      REAL, DIMENSION(ND), INTENT(IN) :: VELO, RADIUS, TAUROSS, T, XMU
      REAL, DIMENSION(ND), INTENT(INOUT) :: VTURB
      
      INTEGER, PARAMETER :: NMAXPAR = 40
      
      CHARACTER*(*) :: VTURB_LINE
      CHARACTER(40), DIMENSION(NMAXPAR) :: ACTPAR
      
      INTEGER, EXTERNAL :: IDX
      
      REAL :: VTURBND, VMICND, VNORM, VMIC_FRAC, VMIC_MAX, DX, X, Q,
     >        DPARVMIC_IN, DPARVMIC_OUT, DUMMY, 
     >        V_INCREASE, TAU_INCREASE, R_INCREASE
      INTEGER :: I, NPAR, L

      LOGICAL, SAVE :: bWARN
      
C***  Physical constants
      REAL, PARAMETER :: PI = 3.141592654 

C***  Local parameters and conversion factors      
      REAL, PARAMETER ::  VMICFRAC_DEFAULT = 0.05
      REAL, PARAMETER ::  FVTURBTOMIC      = SQRT(2.)
      REAL, PARAMETER ::  FVMICTOTURB      = 1./FVTURBTOMIC
            
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      DATA bWARN /.FALSE./  
    

C***  Read number of parameters and put them into ACTPAR array      
      CALL SARGC(VTURB_LINE, NPAR)
      IF (NPAR <= 1) THEN
C***    No VTURB/VMIC line in CARDS file
        DO L=1, ND
          VTURB(L) = 0.
        ENDDO
        IF (.NOT. bWARN) THEN
           WRITE (hCPR,'(A)') 
     >        'VTURB_SETUP> No VTURB or VMIC card has been found!'
           bWARN = .TRUE.
        ENDIF
        RETURN
      ENDIF
      
      
      DO I=1, MIN(NPAR, NMAXPAR)
        CALL SARGV(VTURB_LINE,I,ACTPAR(I))
      ENDDO
      
C***  First parameter: Keyword VTURB or VMIC
C       for VMIC, values are divided by sqrt(2) 
C       to be in line with definition of VMIC in FORMAL      
      IF (ACTPAR(1) == 'VMIC') THEN
        VNORM = 1.
      ELSE
        VNORM = FVTURBTOMIC
      ENDIF
      
C***  Second parameter: standard value for WRSTART initializing     
      READ (ACTPAR(2),'(F20.0)', ERR=99) VTURBND
      VMICND = VNORM * VTURBND
      
      IF (NPAR == 2) THEN
C***    Only two parameters: constant VMIC/VTURB
        DO L=1, ND
          VTURB(L) = FVMICTOTURB * VMICND
        ENDDO      
      ELSE
C***    More than 2 parameters: VMIC/VTURB not constant 
        IF (ACTPAR(3) == 'VELOFRAC') THEN
C***      VELOFRAC option: VMIC is scaled with wind velocity        
          IF (NPAR > 3) THEN
            READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_FRAC
          ELSE 
            VMIC_FRAC = VMICFRAC_DEFAULT
          ENDIF
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * MAX(VMICND, VMIC_FRAC * VELO(L))
          ENDDO
        ELSEIF (ACTPAR(3) == 'HASER') THEN
C***      Similar to MAX option, but with slightly different increasing
C***       motivated by Haser et al. (1998)
          IF (NPAR == 3) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: HASER option needs specified MAX value!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          VMIC_MAX = VNORM * VMIC_MAX
          VMIC_FRAC = (VMIC_MAX - VMICND) / VELO(1)
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * (VMICND + VMIC_FRAC * VELO(L))
          ENDDO
        ELSEIF (ACTPAR(3) == 'HILLIER') THEN
C***      Approach similar to the description used in Hillier's depth-dependent clumping
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: HILLIER option needs two more '//
     >        'parameters, VMIC_MAX and V_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) V_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >                  (VMIC_MAX - VMICND) * EXP(-VELO(L)/V_INCREASE))
          ENDDO
        ELSEIF (ACTPAR(3) == 'EXPTAU') THEN
C***      Exponential VTURB onset (like in Hillier's formula)
C***      but using the Tauross_cont scale
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: EXPTAU option needs two more '//
     >        'parameters, VMIC_MAX and TAU_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) TAU_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            IF (TAUROSS(L) <= 1.E-30) THEN
              VTURB(L) = FVMICTOTURB * VMIC_MAX 
            ELSE
              VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >             (VMIC_MAX - VMICND) * EXP(-TAU_INCREASE/TAUROSS(L)))
            ENDIF
          ENDDO
        ELSEIF (ACTPAR(3) == 'EXPRADIUS') THEN
C***      Exponential VTURB onset (like in Hillier's formula)
C***      but using the Tauross_cont scale
          IF (NPAR < 5) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: EXPRADIUS option needs two more '//
     >        'parameters, VMIC_MAX and R_INCREASE!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          READ(ACTPAR(5), '(F20.0)', ERR=99) R_INCREASE
          VMIC_MAX = VNORM * VMIC_MAX
          DO L=1, ND
            IF (TAUROSS(L) <= 1.E-30) THEN
              VTURB(L) = FVMICTOTURB * VMIC_MAX 
            ELSE
              VTURB(L) = FVMICTOTURB * (VMIC_MAX -
     >         (VMIC_MAX-VMICND) * EXP(-(RADIUS(L)-1.)/(R_INCREASE-1.)))
            ENDIF
          ENDDO
        ELSEIF (ACTPAR(3) == 'MAX') THEN
C***      MAX option: Outer maximum of VMIC/VTURB is given
          IF (NPAR == 3) THEN
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: MAX option must have a specified value!'
            GOTO 99
          ENDIF
          READ(ACTPAR(4), '(F20.0)', ERR=99) VMIC_MAX
          VMIC_MAX = VNORM * VMIC_MAX
          IF (NPAR == 4) THEN
C***        No further parameters: increase with wind velocity (fraction)
            VMIC_FRAC = VMIC_MAX / VELO(1)
            DO L=1, ND
              VTURB(L) = FVMICTOTURB * MAX(VMICND, VMIC_FRAC * VELO(L))
            ENDDO
          ELSEIF (NPAR == 8) THEN                        
            READ(ACTPAR(6), '(F20.0)', ERR=99) DPARVMIC_IN
            READ(ACTPAR(8), '(F20.0)', ERR=99) DPARVMIC_OUT
            IF (ACTPAR(5) == 'VELO1' .AND. ACTPAR(7) == 'VELO2') THEN
C***          Interpolate between min and max value in specified velocity regime            
              IF (DPARVMIC_IN > DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN >= VELO(1)) .OR. 
     >            (DPARVMIC_OUT <= VELO(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner wind velocity: ', VELO(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer wind velocity: ', VELO(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MAX(DPARVMIC_IN, VELO(ND))
              DPARVMIC_OUT = MIN(DPARVMIC_OUT, VELO(1))
              DX = DPARVMIC_OUT - DPARVMIC_IN
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (VELO(L) <= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (VELO(L) >= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( VELO(L) - DPARVMIC_IN ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
              
            ELSEIF (ACTPAR(5) == 'TAU1' .AND. ACTPAR(7) == 'TAU2') THEN
C***          Interpolate between min and max value in specified TAU regime            
              IF (DPARVMIC_IN < DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN <= TAUROSS(1)) .OR. 
     >            (DPARVMIC_OUT >= TAUROSS(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner Tau: ', TAUROSS(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer Tau: ', TAUROSS(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MIN(DPARVMIC_IN, TAUROSS(ND))
              DPARVMIC_OUT = MAX(DPARVMIC_OUT, TAUROSS(1))
              DX = DPARVMIC_IN - DPARVMIC_OUT
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (TAUROSS(L) >= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (TAUROSS(L) <= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( DPARVMIC_IN - TAUROSS(L) ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
            
            ELSEIF (ACTPAR(5) == 'R1' .AND. ACTPAR(7) == 'R2') THEN
C***          Interpolate between min and max value in specified radius regime            
              IF (DPARVMIC_IN > DPARVMIC_OUT) THEN
                DUMMY = DPARVMIC_IN
                DPARVMIC_IN = DPARVMIC_OUT
                DPARVMIC_OUT = DUMMY
              ENDIF
              IF ((DPARVMIC_IN >= RADIUS(1)) .OR. 
     >            (DPARVMIC_OUT <= RADIUS(ND))) THEN
                WRITE(hCPR,'(A)') '*** FATAL ERORR: VMIC/VTURB MAX '
     >            // 'boundaries not within model boundaries!'
                WRITE(hCPR,'(A,F10.5)') 'Inner radius: ', RADIUS(ND)
                WRITE(hCPR,'(A,F10.5)') 'Outer radius: ', RADIUS(1)
                STOP '*** ERROR in subroutine VTURB_SETUP'
              ENDIF
              DPARVMIC_IN = MAX(DPARVMIC_IN, RADIUS(ND))
              DPARVMIC_OUT = MIN(DPARVMIC_OUT, RADIUS(1))
              DX = DPARVMIC_OUT - DPARVMIC_IN
C***          Interpolation between MIN and MAX value:
              DO L=1, ND
                IF (RADIUS(L) <= DPARVMIC_IN) THEN
                  VTURB(L) = FVMICTOTURB * VMICND
                ELSEIF (RADIUS(L) >= DPARVMIC_OUT) THEN
                  VTURB(L) = FVMICTOTURB * VMIC_MAX
                ELSE
                  X = PI * ( RADIUS(L) - DPARVMIC_IN ) / DX
                  Q = 0.5 + 0.5 * COS(X)
                  VTURB(L) = FVMICTOTURB *
     >                           ( Q * VMICND  + (1.-Q) * VMIC_MAX )
                ENDIF
              ENDDO
            
            ELSE 
              WRITE (hCPR,'(A)') '*** FATAL ERROR: Inconsistent ' //
     >          ' detail options for VTURB/VMIC!'
              GOTO 99
            ENDIF
          ELSE
C***        Format does not seem to be valid          
            WRITE (hCPR,'(A)') 
     >        '*** FATAL ERROR: Unknown detail options for VTURB/VMIC!'
            GOTO 99
          ENDIF
          
        ENDIF
      ENDIF
      
      RETURN
      
C***  Error handling

   99 WRITE (hCPR,'(A)') 
     >  '*** FATAL ERROR: Decoding error in VTURB/VMIC line:', 
     >  VTURB_LINE
      STOP 'ERROR IN VTURB_SETUP'
      
      END
      SUBROUTINE WRITMS(ICHANNEL, !file handle 
     >                         X, !(fortran) variable that should be saved
     >                      NDIM, !size of variable (length of array)
     >                      NAME, !index name in mass-storage file (must be <= 8 characters)
     >                  IKINDSTR, !kind string (use -1 for default kind)
     >                    IDUMMY, !unused integer dummy
     >                      IERR)
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************

      IMPLICIT NONE

      CHARACTER(8) :: BUFFER8, NAME
      INTEGER :: IKINDSTR, ICHANNEL, IDUMMY, IERR, NDIM
      REAL :: X

      !IKINDSTR contains a string containing type and format
      ! first character defines the type, numbers following specifiy the kind
      ! NOTE: for compartibility issues IKINDSTR must be read as integer by the subroutine
      IF (IKINDSTR /= -1) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') IKINDSTR    !transform from integer into character
      ELSE
        BUFFER8 = '        '
      ENDIF

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, BUFFER8, X, NDIM, 
     >              'WRITE', IERR)

      RETURN
      END
C***  MAIN PROGRAM WRSTART  ****************************************************
      SUBROUTINE WRSTART
C***********************************************************************
C***  THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
C***  CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
C***  IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
C***    AND (IF NOT TAKEN FROM OLD MODEL) THE FREQUENCY GRID (FILE FGRID)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
 
C***  DEFINE ARRAY DIMENSIONS
      INTEGER, PARAMETER :: MAXATOM =          26 
      INTEGER, PARAMETER :: NDIM    =         2560 
      INTEGER, PARAMETER :: NFDIM   =  2*NDIM + 400 + 1000  !wrstart needs higher NFDIM buffer due to FGRID setup
      INTEGER, PARAMETER :: MAXIND  =        40000 
      INTEGER, PARAMETER :: MAXFEIND  =       2500 
      INTEGER, PARAMETER :: MAXKONT =      NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR =         NDIM 
      INTEGER, PARAMETER :: NDDIM   =           89 
      INTEGER, PARAMETER :: NPDIM   =           94 
      INTEGER, PARAMETER :: MAXHIST =         4000 
      INTEGER, PARAMETER :: MAXXDAT =           10 
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION  =   27 

C***  COMMON /VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS (HERE USED: VMIN)
      COMMON /VELPAR/ VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
C***  OLD RADIUS AND VELOCITY ARE TRANSFERRED TO FUNCTION WRVEL 
C***     BY SPECIAL COMMON BLOCKS:
      COMMON /COMRADI/ OLDRADI(NDDIM)
      COMMON /COMVELO/ NEWVELO, NDv, OLDVELO(NDDIM)

      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT, VCRIT, VCRITold      
      !INCRIT in common block has been replaced by IDUMP as it is never used there

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 3200
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO, WSTABIL

      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND,
     >                            INDRBS, IFRBSSTA, IFRBSEND, IFES
      INTEGER, DIMENSION(MAXKONT) :: IGAUNT, KONTNUP, KONTLOW, 
     >                               KEYCBF, NFEDGE
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(NDDIM) :: IDUMP
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ENLTE
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      REAL, DIMENSION(NDDIM) :: RADIUS, ENTOT, T, TOLD, VELO, GRADI,
     >                          ROLD, RNE, XJC, XJL, TAUROSSOLD, 
     >                          TAURCONT, TAURCONTOLD,
     >                          RHO, XMU, VELOold, RADIUSold, GEFFL,
     >                          ARAD, APRESS, VTEMP, DR, RI, OLDGRADI,
     >                          GAMMARAD, EN, TAUTHOM, TAUGREY, 
     >                          VTURB, VMIC, GRSTATIC,
C***  PROVIDE VARIABLES FOR THE READ OF THE SECOND MODEL FOR TEMPERATURE 
C***    INTERPOLATION
     >                          T2, RADIUS2, TOLD2, ROLD2, TEFFOLD2
      REAL, DIMENSION(NPDIM) :: P
      REAL, DIMENSION(NDDIM,NPDIM) :: Z
      REAL, DIMENSION(NDDIM,NFDIM) :: XJCold
      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, FWEIGHT,
     >                          EXPFAC, OPAC, ETAC, XLAMBDAold
      REAL, DIMENSION(NFDIM,MAXKONT) :: SIGMAKI
      REAL, DIMENSION(NFDIM,0:MAXION) :: SIGMAFF
      REAL, DIMENSION(4,NDIM) :: ALTESUM
      REAL, DIMENSION(4,MAXIND) :: COCO
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NDDIM,NDIM) :: POPNUM

      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(MAXION) :: TFEEXC

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C      INTEGER, PARAMETER :: INDEXMAX = 4000000, NFEREADMAX = 100000     !standard vd100
C      INTEGER, PARAMETER :: INDEXMAX = 7000000, NFEREADMAX = 300000     !standard vd100
C      INTEGER, PARAMETER :: INDEXMAX = 5000000, NFEREADMAX = 150000     !vd50
c      INTEGER, PARAMETER :: INDEXMAX = 15000000, NFEREADMAX = 400000    !vd20
      INTEGER, PARAMETER :: INDEXMAX = 100000000, NFEREADMAX = 600000    !custom hydro
c      INTEGER, PARAMETER :: INDEXMAX = 100000000, NFEREADMAX = 600000    !xxl branch
C      INTEGER, PARAMETER :: INDEXMAX = 12000000, NFEREADMAX = 300000    !Goetz       
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      LOGICAL :: bFEULSEP

      INTEGER :: N, ND, NDold, JOBNUM, NATOM, LASTIND, IDX, IERR,
     >           NAUTO, IDUMMY, MAXITTAU, ITTAU, IERRVELO, LASTINDALL,
     >           LASTKON, LASTKDR, LASTFE, LAST, IFLAG, L, J,
     >           NA, NF, NF2, MASSORIGIN, NP, JOBNOLD, JOBNOLD2,
     >           NEXTK, LRTinput, NC, NDv, NFold, LASTSELF, 
     >           LOW, INDDR, MDOTINPUT, N_WITH_DRLEVELS

      REAL :: GLOG, GEFFLOG, GEDD, TROLD, VMINOLD, VMINOLD2, FM,
     >        qpLINK_USER, BLACKEDGE, VA, RCSAVE, Vfac, GFLSAV,
     >        VDOPFE, DXFE, XLAM0FE, XMSTARold, DTAUCDTAUTH,
     >        TAUMAX, TAUACC, VFINAL, XMASS, ATMEAN, STAPEL, AMIN,
     >        RSTAR, RMAX, VDOP, XMSTAR, TFAC, XMDOT, XLOGL, RTRANS,
     >        DENSCON_FIX, VMIN, TEFF, TEFFOLD, TMIN, TMIN2, TS, RL,
     >        R23, TAU23, OLDRADI, OLDVELO, TROLD2, RCON, TOTOUT, BETA,
     >        VPAR1, VPAR2, BETA2, BETA2FRACTION, HSCALE, VMINhydro,
     >        VPAR1_2, VPAR2_2, TMODIFY, SPHERIC, XMDOTold, VoldMod,
     >        RCONold, fHYDROSTART, XMG, XMGold, VTURBND, RADEXP,
     >        dummy, GRADLAST, GRADIL, STEPDAMP, RINT, PL, PLP,
     >        RHOINT, VINMAX, VMINcard, GAMMAL, RSTARold, TL, ENE,
     >        RADGAMMASTART, GEDDRAD, GEDDPRINT, RCRIT, POPMIN,
     >        XLAMBLUE, VNDold, VFINALold, RMAXnew,
     >        VL, DVDR, NORMINERTIA, DTDRIN_OLD, DV, RRMAXDIFF,
     >        DCINF_OLD, XLOGLold
  
      CHARACTER(MAXHIST*8) :: MODHIST
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters

      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(100) :: MODHEAD, MODOLD, MODOLD2
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(10), DIMENSION(MAXAUTO) :: LEVUPAUTO, LEVAUTO
      CHARACTER(8), DIMENSION(NDDIM) :: VELOCRITERION
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2) :: WRTYPE
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(120) :: DENSCON_LINE, ThinCard,VTURB_LINE,DRLINES_CARD
      CHARACTER(6) :: BUFFER6
      CHARACTER(8) :: BUFFER8, GEFFKEY, CVEXTEND
      CHARACTER(9) :: MLRELATION
      CHARACTER(144) :: BUFFER144

      LOGICAL :: TTABLE, OLDTEMP, TEXIST, VPLOT, THIN, NEWVELO, BTWOT,
     >           OLDFGRID
      LOGICAL :: BFEMODEL
      LOGICAL :: BTAUR
      LOGICAL :: bTauFix,             !true if fix option has been set in the TAUMAX CARDS line
     >           bHydroStat,          !true if hydrostatic domain encountered in velocity law
     >           bNoRGrid,            !Parameter for GEOMESH, defines if radius grid should be (re)made or not
     >           bNoDetails,          !decides whether PRIMOD should print all values per depth point or not
     >           bSaveGEFF,            !if true geff is used for hydrostatic scale height instead of g
     >           bOLDSTART,           !true if OLDSTART CARDS option has been set
     >           bThinImprove,        !true if THIN artistic has already been run in this iteration
     >           bOldMdot,            !true if MDOT option has been set in OLD V line
     >           bHYDROSOLVE,         !true if hydrodynamic consistent model should be obtained
     >           bOVTauMax,           !true if TAUMAX option has been set on OLD STRATIFICATION card
     >           bOVTauCut,           !true if TAUCUT option has been set on OLD V line
     >           bHScaleOnly,         !determines if INITVEL calulates only HSCALE or full static law
     >           bFULLHYDROSTAT,      !if true, VELTHIN uses full GAMMA instead of EDDINGTON GAMMA
     >           bGAMMARADMEAN,       !if true, a mean value for GAMMARAD is used instead of individuals
     >           bGREYSTART,
     >           bNDfirst,
     >           bSMOCO,
     >           bOLDMODEL,           !true is old MODEL file exists
     >           LTESTART,
     >           bForceDCUpdate,
     >           bOLDJ,               !true if XJC from old MODEL file is used  
     >           bDCSCALE             !true if old DENSCON stratification should be scaled to new D_inf
      INTEGER :: GEddFix,             ! > 0 if GEDD should kept fixed in all calculations
     >           iOldStratification,  ! > 0 if OLD STRATIFICATION option has been set in the CARDS file
     >           iOLDRAD              !use continuum approximation (1) or old radiation field (2)
                                      ! for HYDROSTATIC INTEGRATION if possible

      REAL, EXTERNAL :: WRVEL

C***  Tiefenabhaengiges Clumping nach Goetz Graefener
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC, 
     >                          DENSCON_OLD, FILLFAC_OLD

C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE FROM SUBR. DECSTAR
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  CONTROL PARAMETER FOR TABULATED INPUT OF T(R) AND V(R)
      INTEGER :: ITAB = 0
       
      !Konstanten
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: AMU = 1.66E-24         !atomic mass unit (constant)
      REAL, PARAMETER :: CLIGHT = 2.99792458E10 !SPEED OF LIGHT IN CM / SECOND
      REAL, PARAMETER :: RSUN = 6.96E10         !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: BOLTZ = 1.38E-16       !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)
      REAL, PARAMETER :: GCONST = 6.670E-8      !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)     
      REAL, PARAMETER :: TEFFSUN = 5780.        !SOLAR EFFECTIVE TEMPERATURE
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !(new) MODEL file
      INTEGER, PARAMETER :: hOldMODEL = 9   !old MODEL file (used in OLDSTART)
      INTEGER, PARAMETER :: hOldMODEL2 = 8  !second old "model" file (used if TWOTEMP cards option is set)

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      bHydroStat = .TRUE.      !default value: hydrostatic domain is achieved


C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> WRSTART started: Program Version from ', 
     >                LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL

      IF (OPSYS .EQ. 'CRAY') THEN
        CALL CLOCK(TIM1)
c      ELSE
c        CALL TIME(TIM1)
      ENDIF

      !Put dummy stuff into COMMON block data to ensure storage reservation
      NDv = NDDIM
      NEWVELO = .FALSE.
      DO L=1, NDDIM
        OLDRADI(L) = 1337.
        OLDVELO(L) = 1338.
      ENDDO
      
C***  JOB NUMBER OF THIS JOB
      JOBNUM=0

C***  INITIALIZE SOME VARIABLES (TODT 04.05.2010)
      TROLD = 0.
      RCRIT = -1.
      DTDRIN_OLD = -1.

C***  READ ATOMIC DATA FROM FILE "DATOM"
      CALL       DATOM (NDIM,N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND,KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'WRSTART', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL,
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

C***  Check for possible old model and read values that
C***  might be needed before DECSTAR due to CARDS options
      CALL OPENMS (hOldMODEL, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (hOldMODEL, RSTARold,    1, 'RSTAR   ', IERR)
      bOLDMODEL = (IERR /= -10)
      IF (bOLDMODEL) THEN
        CALL READMS (hOldMODEL, TEFFOLD,     1, 'TEFF    ', IERR)
        XLOGLold = ALOG10( (RSTARold/RSUN)**2 * (TEFFOLD/TEFFSUN)**4 )
        CALL READMS (hOldMODEL, XMDOTold,    1, 'XMDOT   ', IERR)
        IF (IERR == -10) THEN
          XMDOTold = -999.
        ENDIF
        CALL READMS (hOldMODEL, XMSTARold,   1, 'XMSTAR  ', IERR)
        CALL READMS (hOldMODEL, NDold  ,     1, 'ND      ', IERR)
        CALL READMS (hOldMODEL, DTDRIN_OLD,  1, 'DTDRIN  ', IERR)
        CALL READMS (hOldMODEL, VELOold, NDold, 'VELO    ', IERRVELO)
        VFINALold = VELOold(1)
        VNDold = VELOold(NDold)
      ELSE 
        IERRVELO = -10
        VNDold = -1.
        VFINALold = -1.
      ENDIF
      CALL CLOSMS (hOldMODEL, IERR)
 
C***  DECODING INPUT DATA
C***  VDENS1,VDENS2,DENSCON_FIX und CLUMP_CRIT nach goetz...
      CALL DECSTAR (MODHEAD, FM, RSTAR, VDOP, RMAX, TTABLE, NDDIM,
     >              OLDTEMP, MAXATOM, NATOM, ABXYZ, KODAT, VPLOT, 
     >              ATMASS, XDATA, MAXXDAT, OLDFGRID, THIN, ThinCard,
     >              GLOG, GEFFLOG, GEDD, bSaveGEFF, XMSTAR, WRTYPE, 
     >              TAUMAX, TAUACC, bTauFix, BTWOT, TFAC, DENSCON_LINE, 
     >              BLACKEDGE, bOLDSTART, RadiusGridParameters, XMDOT,
     >              XLOGL, RTRANS, BTAUR, DENSCON_FIX, MASSORIGIN,
     >              LRTinput, ELEMENT, iOldStratification, bHYDROSOLVE,
     >              GEddFix, bOldMdot, VoldMod, bOVTauMax, MLRELATION,
     >              NC, VTURBND, VTURB_LINE,
     >              iOLDRAD, RADGAMMASTART, fHYDROSTART, 
     >              bFULLHYDROSTAT, bGAMMARADMEAN, CVEXTEND, bGREYSTART,
     >              GEFFKEY, POPMIN, XLAMBLUE, LTESTART, 
     >              bOLDMODEL, XMDOTold, XMSTARold, RSTARold, VNDold,
     >              VFINALold, bForceDCUpdate, bOLDJ, bOVTauCut, 
     >              bDCSCALE, TEFFOLD, XLOGLold, MDOTINPUT, 
     >              DRLINES_CARD)

C***  GENERATION OF THE CONTINUOUS FREQUENCY GRID
      CALL       FGRID (NFDIM, NF, XLAMBDA, FWEIGHT, KEY, NOM, SYMBOL, 
     $                  N, NCHARG, ELEVEL, EION, EINST, NDIM,
     $                  EDGEK,KODAT,MAXATOM, MAXION,
     $                  INDNUP, INDLOW, LASTIND, KONTNUP, KONTLOW, 
     >                  LASTKON, OLDFGRID, NF2, XLAMBDA2, 
     >                  VDOP, XLAMBLUE)
 
C***  ATMEAN = MEAN ATOMIC WEIGHT ( IN AMU )
      !What does ABXYZ contain?
      ATMEAN=0.
      DO NA=1,NATOM
          ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO

C***  STAPEL: NUMBER OF FREE ELECTRONS PER ATOM
C***   S T A R T   A P P R O X I M A T I O N  
      !STAPEL = RNE start approx
      STAPEL=0.
      DO NA=1,NATOM
         STAPEL = STAPEL + ABXYZ(NA) * (STAGE(NA)-1.)
      ENDDO

C***  MEAN MASS PER PARTICLE (ATOMS AND ELECTRONS) IN AMU
C***  -  HAND INPUT, AS IN OLDER PROGRAM VERSIONS, IS IGNORED
      XMASS = ATMEAN / (1. + STAPEL)
 
C***  LOOP FOR AUTOMATICAL ADJUSTMENT OF MAX. ROSSELAND DEPTH (OPTINAL)
      ITTAU = 0

      VMINOLD  = 0.
      VMINOLD2 = 0.

C***  REQUIRED ACCURACY (IN TAU): 
      !TAUACC = 1.E-4  !(This is now a CARDS option with this default value)
C***  MAXIMUM NUMBER OF ITERATIONS
      MAXITTAU = 30
      
      IF (fHYDROSTART > 0.) THEN
        !Hydrostart needs old stratification
        iOldStratification = 2
        THIN = .FALSE.      
      ENDIF
      
      VMINcard = VMIN
      IF (bOLDJ) THEN
        CALL OPENMS (hOldMODEL, IDUMMY, IDUMMY, 1, IERR)
        CALL LOADOLDJ(XJCold, XLAMBDAold, NDold, NFold, 
     >                NDDIM, NFDIM, hOldMODEL)
        CALL CLOSMS (hOldMODEL, IERR)
      ENDIF
      IF (bOLDSTART .OR. iOLDRAD > 0 .OR. iOldStratification > 0) THEN
        !use old vmin from old MODEL instead of CARDS value if OLDSTART has been set
        CALL OPENMS (hOldMODEL, IDUMMY, IDUMMY, 1, IERR)
        CALL READMS (hOldMODEL, ARAD, NDold-1,   'ARAD    ', IERR)
        CALL READMS (hOldMODEL,RADIUSold,NDold,'R       ',IERR)
        IF (IERR == -10 .AND. iOldStratification > 0) THEN
          WRITE (hCPR, '(A)') 'Error: Cannot find old radius grid!'
          STOP '*** FATAL ERROR IN WRSTART'
        ENDIF
        CALL READMS (hOldMODEL,RCONold,NDold,'RCON    ',IERR)
        CALL READMS (hOldMODEL,TAUROSSOLD,NDold,'TAUROSS ',IERR)
        CALL READMS (hOldMODEL,TAURCONTOLD, NDold, 'TAURCONT', IERR)
        IF (IERR == -10) THEN
          WRITE (hCPR, '(A)') 'Warning: Could not find TAURCONT in '
     >      // 'old MODEL file. Using TAUROSS instead.'
          DO L=1, NDold
            TAURCONTOLD(L) = TAUROSSOLD(L)
          ENDDO
        ENDIF
        CALL READMS (hOldMODEL,DENSCON_OLD,NDold,'DENSCON ',IERR)
        DCINF_OLD = DENSCON_OLD(1)
        DO L=1,NDold
          IF (DENSCON_OLD(L) <= 0. ) THEN
              IF (DENSCON_OLD(1) <= 0.) THEN
                CALL REMARK (
     >            'Zero or negative depth-dep. clumping in OLD model!')
                STOP 'Error in reading depth-dep. clumping '
     >            // 'from OLD model during WRSTART'
              ENDIF 
              DENSCON_OLD(L) = DENSCON_OLD(1)
          ENDIF  
          IF (bDCSCALE) THEN
C***        OLD V DCSCALE option:
C***        scale old clumping stratification to match a new D_inf
            DENSCON_OLD(L) = (DENSCON_OLD(L) - 1.) 
     >        * (DENSCON_FIX - 1.)/(DCINF_OLD - 1.) + 1.
          ENDIF
          FILLFAC_OLD(L) = 1. / DENSCON_OLD(L)
        ENDDO
        
        XMGold = XMSTARold * GCONST * XMSUN
        IF (IERR /= -10) THEN
          IF (IERRVELO /= -10) THEN
            IF (iOldStratification > 0) THEN
              VMIN = VELOold(NDold)
              CALL READMS (hOldMODEL,VCRITold,NDold,'VCRIT   ' ,IERR)
              IF (IERR == -10) THEN
                VCRITold = '        '
              ENDIF
              WRITE (hCPR,'(A)') 'Using old velocity field...'
              DO L=1, NDold
                OLDRADI(L) = RADIUSold(L)
                OLDVELO(L) = VELOold(L)
              ENDDO              
              NDv = NDold
              VFINAL = VELOold(1)
              IF (VoldMod /= 1.) THEN
                IF (VoldMod < 0.) THEN
                  VFINAL = ABS(VoldMod)
                  VoldMod = VFINAL / VELOold(1)
                ENDIF
                !now scale with given facor
                DO L=1, NDold
                  OLDVELO(L) = VoldMod * VELOold(L)
                ENDDO
                VFINAL = OLDVELO(1)
                VMIN = OLDVELO(NDold)
                WRITE (hCPR,'(A,F8.2)') 
     >            'Adjusted old velocity field to VFINAL = ', VFINAL
              ENDIF
              IF (LEN_TRIM(CVEXTEND) > 0) THEN
C***            CVEXTEND options allow you to change
C***            (and especially enlarge) RMAX when using OLD V 
C***            We do this on a logarithmic scale of r/RSTAR in order 
C***            to get a suitable stretch on the whole range.
                RADEXP = LOG10(RMAX)/LOG10(OLDRADI(1))
                DO L=1, NDold
                  RADIUSold(L) = RADIUSold(L)**RADEXP
                ENDDO
C***            If RCONold is outside of radius grid, we need to ensure this now#
                IF (RCONold > OLDRADI(1)) THEN
                  RCONold = 1.1 * RADIUSold(1) + 1.
                ENDIF      
C***            @TODO: HOW TO HANDLE OLD V TAU?         
                IF (CVEXTEND == 'STRETCH') THEN
C***              In the STRETCH option, we scale the old velocity
C***              with regards to the radius scale, i.e. the radius
C***              vector is changed but the velocity field vector is left 
C***              untouched
                  IF (RCONold <= OLDRADI(1)) THEN
                    RCONold = RCONold**RADEXP
                  ENDIF
                  DO L=1, NDold
                    OLDRADI(L) = OLDRADI(L)**RADEXP
                  ENDDO
                  WRITE (hCPR,'(A,F8.2)') 'Old velocity field '
     >              // ' is stretched to a new RMAX by exponentiation '
     >              // ' of the radius grid to the power ', RADEXP
                ELSEIF (CVEXTEND == 'EXTRAP') THEN
C***              In the EXTRAP option the velocity field is simply
C***              continued outwards using the outermost gradient
                  CALL SPLINPOX(VL, OLDRADI(1), OLDVELO, OLDRADI, NDold,
     >                                     DFDX=DVDR)
ccc               TEST: Outermost point should have intended v_inf 
                  OLDRADI(1) = RMAX
                  DO L=1, NDold
                    IF (RADIUSold(L) > OLDRADI(1)) THEN
C***                  linear extrapolation of v(r)
                      VELOold(L) = OLDVELO(1) 
     >                              + DVDR * (RADIUSold(L)-OLDRADI(1)) 
                    ELSEIF (RADIUSold(L) < OLDRADI(NDold)) THEN
C***                  This should only happen due to numerical inaccuracy
C***                    and thus never more than once!
                      VELOold(L) = OLDVELO(NDold)
                    ELSE 
                      CALL SPLINPOX(VELOold(L), RADIUSold(L), 
     >                                          OLDVELO, OLDRADI, NDold)
                    ENDIF
                  ENDDO
C***              Copy into OLD*-vectors which are later used for
C***              estabilishing the new grid in GEOMESH->RGRID->WRVEL                  
                  DO L=1, NDold
                    OLDRADI(L) = RADIUSold(L)
                    OLDVELO(L) = VELOold(L)
                  ENDDO              
C***              RCONold is not changed since we only extrapolate                   
                  WRITE (hCPR,'(A,F8.2)') 'Old velocity field '
     >              // ' is extended to a new RMAX by extrapolation. '
                ENDIF
C***            Ensure outer boundary is exactly at desired RMAX
                OLDRADI(1) = RMAX
                RADIUSold(1) = RMAX
              ELSEIF (RADIUSold(1) < RMAX) THEN
C***            if this triggers, it is sometimes just due to numerical accuracy
C***            TODO: Check by how many digits RMAX is given in the CARDS
C***                  and test if the old value differs by less than this accuracy
                RRMAXDIFF = ABS(1. - RADIUSold(1)/RMAX)
                IF (RRMAXDIFF < 1.E-9) THEN
C***              We attribute this to numerical accuracy: Set Rold(1) = RMAX
                  OLDRADI(1) = RMAX
                  RADIUSold(1) = RMAX
                  WRITE (hCPR,'(A,F12.5)') 
     >              '*** WARNING: RMAX was adjusted by ', RRMAXDIFF                      
                ELSE 
                  RMAXnew = RMAX
                  RMAX = RADIUSold(1)
                  WRITE (hCPR,'(A,F12.5)') 
     >              '*** ERROR: OLD V forces lower RMAX = ', RMAX     
                  WRITE (hCPR,'(2(A,F12.5))') 
     >              '*** Current requested RMAX is = ', RMAXnew, 
     >                  ' which differs by ', RRMAXDIFF     
                  WRITE (hCPR,'(A)') '*** Please restart the model with'
     >              // ' one of the following options:'
                  WRITE (hCPR,'(A,F12.5)') '  1) Adjust the RMAX option'
     >              // ' on the VELPAR card to the lower value RMAX = ', 
     >              RMAX
                  WRITE (hCPR,'(A)') '     (Make sure that the'
     >              // ' RMAX_IN_RSUN option is switched off)'
                  WRITE (hCPR,'(A)') '  2) Add the STRETCH option to'
     >              // ' the OLD STRATIFICATION (or OLD V) card'
                  WRITE (hCPR,'(A)') '     (This will stretch the old'
     >              // ' velocity field to the new RMAX)'
                  WRITE (hCPR,'(A)') '  3) Add the EXTRAP option to the'
     >              // ' OLD STRATIFICATION (or OLD V) card'
                  WRITE (hCPR,'(A)') '     (This will extrapolate the'
     >              // ' old velocity field to the new RMAX)'
                  STOP '*** FATAL ERROR IN WRSTART' 
                ENDIF 
              ENDIF
C***          We Use old RCON if inside the old radius grid
C***          However, this can only be done after the GEOMESH call,
C***          otherwise WRVEL would fail as we do not know the old 
C***          BETA law parameters. Therefore ensure interpolation 
C***          in WRVEL by setting RCON > RMAX
              RCON = 1.5 * RADIUSold(1)
C***          TODO: check if the following lines still make sense!!!
              IF (bOldMdot) THEN
                CALL READMS (hOldMODEL,XMDOTold,1,'XMDOT   ',IERR)
                IF (VoldMod /= 1.) THEN
                  XMDOTold = XMDOTold + LOG10(VoldMod)
                ENDIF
                IF (IERR /= -10) THEN
                  XMDOT = XMDOTold 
                  FM= 10.**(XMDOT+3.02) * (RSUN/RSTAR)**2
                  RTRANS = 0.
                ENDIF
              ENDIF
              IF (bOVTauMax .OR. fHYDROSTART > 0.) THEN
                CALL READMS (hOldMODEL, T, NDold,     'T       ', IERR)
                CALL READMS (hOldMODEL, RNE, NDold,   'RNE     ', IERR)
                DO L=1, NDold
                  XMU(L) = ATMEAN / (1. + RNE(L))
                ENDDO
              ENDIF
              IF (bOVTauMax) THEN
                CALL READMS (hOldMODEL, RCSAVE, 1, 'RCSAVE  ', IERR)
                IF (IERR == -10) THEN
                  !use sonic point instead of critical point
                  sploop: DO L=NDold, 1, -1
                    VA = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5
                    IF (OLDVELO(L) > VA) THEN
                      RCSAVE = OLDRADI(L) 
                      EXIT sploop
                    ENDIF
                  ENDDO sploop                         
                ENDIF
              ENDIF
              IF (fHYDROSTART > 0.) THEN
                CALL READMS (hOldMODEL, ENTOT, NDold, 'ENTOT   ', IERR)
                CALL READMS (hOldMODEL, OLDGRADI, NDold, 'GRADI   ', IERR)
                DO L=1, NDold
                  RHO(L) = ENTOT(L) * AMU * ATMEAN
                ENDDO
                CALL READMS (hOldMODEL, VMIC, NDold, 'VMIC    ', IERR)
                IF (IERR == -10) THEN                
                  CALL READMS (hOldMODEL,VTURB(NDold),1,'VTURB   ',IERR)
                  IF (IERR == -10) THEN
                    VTURB(1:NDold) = 0.
                  ELSE
                    DO L=1, NDold-1
                      VTURB(L) = VTURB(NDold)
                    ENDDO
                  ENDIF
                ELSE
                  DO L=1, NDold
                    VTURB(L) = VMIC(L) / SQRT(2.)
                  ENDDO
                ENDIF
                DO L=1, NDold-1
                  RI(L) = 0.5 * ( OLDRADI(L) + OLDRADI(L+1) )
                  DR(L) = OLDRADI(L) - OLDRADI(L+1) 
                  PL  = RHO(L)*(RGAS*T(L)/XMU(L) + (VTURB(L)*1.E5)**2.)
                  PLP = RHO(L+1)*(RGAS*T(L+1)/XMU(L+1)
     >                                          +(VTURB(L+1)*1.E5)**2.)
                  RHOINT = 0.5 * ( RHO(L) + RHO(L+1) )
                  APRESS(L) = - (PL - PLP) / (DR(L) * RSTAR) / RHOINT
                ENDDO       
                
                VINMAX = MAX(VMINcard, 10 * OLDVELO(NDold))
                VINMAX = VINMAX * 1.E5
                RINT = RI(NDold-1) 
                XMG = GCONST * XMSTAR * XMSUN
                CALL SPLINPOX (GRADIL, RINT, OLDGRADI, OLDRADI, NDold)
                RINT = RINT * RSTAR
                GRADIL = GRADIL * 1.E5 / RSTAR
                VTEMP(ND) = ( ARAD(ND-1)+APRESS(ND-1)-XMG/(RINT**2) )
     >                              /GRADIL 
                IF ( VTEMP(ND) < 0.) THEN
                  VTEMP(ND) = MAX(VMINcard, 10 * OLDVELO(ND))
                  VTEMP(ND) = VTEMP(ND) * 1.E5
                ENDIF
                WRITE (hCPR,*) VTEMP(ND) /1.E5, GRADIL/1.E5*RSTAR
                DO WHILE (VTEMP(ND) > VINMAX .OR. VTEMP(ND) < 0.)
                  IF (VTEMP(ND) < 1.E-10) THEN
                    GRADIL = GRADIL / 2.
                  ELSE
                    GRADIL = GRADIL * 1.2
                  ENDIF
                  VTEMP(ND) = ( ARAD(ND-1)+APRESS(ND-1)-XMG/(RINT**2) )
     >                             / GRADIL
                  WRITE (hCPR,*) VTEMP(ND) /1.E5, GRADIL/1.E5*RSTAR
                ENDDO
                DO L=NDold-1, 2, -1
                  IF (GRADIL > 0) THEN
                    GRADLAST = GRADIL
                  ENDIF
                  VTEMP(L) = VTEMP(L+1) + GRADIL * DR(L) * RSTAR
                  RINT = RI(L) * RSTAR
                  GRADIL = ( ARAD(L)+APRESS(L)-XMG/(RINT**2) ) /VTEMP(L)
                  IF (GRADIL <= 0.) THEN
                    GRADIL = GRADLAST * 0.75
                  ENDIF
                  GRADIL = MAX(0., GRADIL)
C                  WRITE (hCPR,*) 'L, v (km/s) :', L, VTEMP(L) /1.E5
                ENDDO
                VTEMP(1) = VTEMP(2) + GRADIL * DR(1) * RSTAR
                
                !Daempfung mit fHYDROSTART
                DO L=1, NDold
                  VELO(L) = fHYDROSTART * VTEMP(L)/1.E5
     >                       + (1. - fHYDROSTART) * OLDVELO(L)
                  WRITE (hCPR,*) OLDRADI(L), VELO(L)
                ENDDO
                
                !Glaettung
                DO L=2, NDold-1
                  VTEMP(L) = 0.25 * ( LOG10(VELO(L-1)) 
     >              + 2. * LOG10(VELO(L)) + LOG10(VELO(L+1)) )
                ENDDO
                DO L=2, NDold-1
                  VELO(L) = 10.**VTEMP(L)
                ENDDO
                
                !Kopieren nach OLDVELO
                DO L=1, NDold
                  OLDVELO(L) = VELO(L)
                ENDDO
                RCON = RADIUS(1) * 1.5
                WRITE (hCPR,*) '**STARTING with a hydrodynamic approach'
                WRITE (hCPR,*) '**DAMPING FACTOR ', fHYDROSTART
                
              ENDIF
            ENDIF
          ELSEIF (iOldStratification > 0) THEN
            !No stratification data found => create new one
            WRITE (hCPR,*) '**WARNING: OLD STRATIFICATION NOT FOUND **'
            iOldStratification = 0
          ENDIF
        ENDIF
        CALL CLOSMS (hOldMODEL, IERR)
        WRITE (hCPR,*) 'Oldstart: Old model has been read!'
      ELSE 
        WRITE (hCPR,*) 'Fresh start: No old model is used!'
        NDold = 1  !must be set for DIMENSION-statements in PREP_GAMMARAD
      ENDIF
      
C***  Preparation of GAMMARAD for the hydrostatic equation
      CALL PREP_GAMMARAD (iOLDRAD, bFULLHYDROSTAT, GEDD,
     >      GEddFix, GAMMARAD, NDold, ARAD, GLOG, RSTAR,
     >      RADIUSold, RSTARold, XMGold, TAUROSSOLD, RCONold, GEFFLOG,
     >      STAPEL, ATMEAN, XLOGL, XMSTAR, RADGAMMASTART, 
     >      GEDDRAD, iOldStratification, bGAMMARADMEAN, bSaveGEFF )      
      
      bNDfirst = .TRUE.
      XMG = GCONST * XMSTAR * XMSUN
      VMINhydro = -99.      !init negative => not yet calculated
      
C    ------- Main TAU iteration loop starts here ----------------------
      tauit: DO
          ITTAU = ITTAU + 1
C          WRITE (hCPR,*) " ITTAU=", ITTAU

C***      CHECK WETHER NEW VMIN IS INBETWEEN 0. AND VFINAL
C          IF (VMIN < 0.) VMIN = 1.E-4
          IF (VMIN > VFINAL) VMIN = 0.9 * VFINAL
C         IF (VMIN == VMINOLD) VMIN = 0.9 * VMIN        
          IF (VMINhydro > 0. .AND. VMIN > VMINhydro) THEN
            VMIN = VMINhydro
            WRITE (hCPR,*) '*** WARNING: No hydrodynamic solution ***'
          ENDIF

C***      INITIALISATION OF THE VELOCITY-FIELD PARAMETERS
C***      Turbulence pressure added, 13-Mar-2014
          IF (iOldStratification == 0) THEN
            bHScaleOnly = (THIN .AND. (ITTAU > 1))
            CALL INITVEL (RMAX,TEFF,GEFFLOG,RSTAR,XMASS,
     >                    VTURBND, bHScaleOnly, bHydroStat)
            IF (ITTAU == 1 .OR. (.NOT. THIN)) NEWVELO=.TRUE.
          ELSE            
            IF (iOldStratification == 2) THEN
              ND = NDold
              DO L=1, ND
                RADIUS(L) = OLDRADI(L)
                VELO(L) = OLDVELO(L)                
              ENDDO
C***           Unused old version of OVTauMax option:                
c              IF (bOVTauMax .AND. XMDOTold /= XMDOT) THEN
cc                Vfac = VELOold(ND) / VMIN - 1.
c                Vfac = 10**(XMDOT - XMDOTold) - 1.
c                VELO = VELOold
c                vfacloop: DO L=ND, 1, -1
cc                  Vfac = Vfac * EXP(FLOAT(L-ND) / 5.) + 1.
c                  
c                  CALL SPLINPOX(VL, RADIUS(L),VELO,RADIUS,ND,DFDX=DVDR)
c                  
c                  NORMINERTIA = VL * DVDR / XMG * RADIUS(L) * RADIUS(L)
c     >                               * RSTAR * 1.E10
cc                  WRITE (0,*) 'NORMINERTIA = ', L, NORMINERTIA
c                  
c                  Vfac = (VELOold(ND) / VMIN - 1.)
c     >                       * MAX(1.-NORMINERTIA,0.) + 1.
c                  
cc                  IF (RADIUS(L) >= RCSAVE) THEN
cc                    EXIT vfacloop
cc                  ENDIF
c                  VELO(L) = VELO(L) / Vfac
c                  OLDVELO(L) = VELO(L)
c                ENDDO vfacloop
c              ENDIF
            ENDIF
            NEWVELO=.FALSE.
          ENDIF

C***  OLD V WITH TAUMAX option: Initialize OLD..-variables with shifted
C***  velocity field for interpolation on new grid.
          IF (bOVTauMax .AND. VMIN /= VELOold(NDold)) THEN
              DV = VMIN - VELOold(NDold)
              DO L=1, NDold
                OLDVELO(L) = VELOold(L) + DV 
                OLDRADI(L) = RADIUSold(L)
              ENDDO
              WRITE (hCPR,*) 'OVTauMAX: DV = ', DV
          ENDIF
            


C***      LOOP FOR IMPROVED HYDROSTATIC EQUATION ("THIN WIND" OPTION)
C***       (MUST BE DONE TWICE, FIRST TO ESTABLISH A RADIUS MESH AND THE 
C***        TEMPERATURE STRUCTURE, SECOND FOR THE EXACT HYDROSTATIC EQ.)

          bThinImprove = .FALSE.
c          WRITE (0,*) 'ITTAU loop: ', ITTAU

C         ------- THIN WIND loop ----------------------
          thinwind: DO

C***  GENERATION OF THE RADIUS GRID, P-GRID, Z-GRID, AND ND
            IF (iOldStratification > 1) THEN
              bNoRGrid = .TRUE.
              IF (fHYDROSTART > 0.) THEN
                INCRIT = 'HDSTART '
              ELSE
                INCRIT = 'OLDMODEL'
              ENDIF
            ELSE
              bNoRGrid = .FALSE.
            ENDIF
c            WRITE (hCPR,* ) 'VMIN now = ', VMIN
c            WRITE (hCPR,* ) 'HSCALE now = ', HSCALE
c            DO L=1, NDv
c              WRITE (0,*) 'OLD: ', OLDRADI(L), OLDVELO(L)
c            ENDDO
            CALL GEOMESH (RADIUS,INCRIT,P,Z,ND,NDDIM,NP,NPDIM,RMAX, 
     >                    RadiusGridParameters, bNoRGrid,
     >                    NC, XMG, RSTAR, RCRIT)
            IF (bNDfirst) THEN
              WRITE (hCPR,*) ' ND, NDold: ', ND, NDold
              IF (bGREYSTART) THEN
                WRITE (hCPR,*) 'TAUROSS scale used: GREY'
              ELSE
                WRITE (hCPR,*) 'TAUROSS scale used: CONT'
              ENDIF
              bNDfirst = .FALSE.
            ENDIF

C***  START APPROXIMATION FOR EL. DENSITY PUT INTO ARRAY (NEEDS ND)
            DO L=1, ND
              RNE(L) = STAPEL
              XMU(L) = ATMEAN / (1. + RNE(L))
            ENDDO

C***        TAUROSS must be initialized in the first run            
            IF (iOldStratification > 0 .AND. ITTAU == 1) THEN
              DO L=1, ND
                IF (bOVTauMax) THEN
C***              If we iterate OLD V on Tau, we cannot use the Tau-Scale from the old model
C***              but instead will start with OLD DENSCON and then iterate with the new DENSCON
C***              from ITTAU > 1 onwards
                  TAURCONT(L) = -99.
                ELSE
C***              If we do not iterate OLD V on TauMax, there is just one run, so we need a 
C***              Tau-Scale for the DENSCON setup.
C***              (This option should be used if there is no serious change of TAUMAX)
                  IF (L == 1) THEN
C***                Numerical interpolation on the edge might give negative value, ensure zero
                    TAURCONT(L) = 0.
                  ELSE
                    CALL SPLINPOX(TAURCONT(L), RADIUS(L),
     >                             TAURCONTOLD, OLDRADI, NDold)
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
            
   
C***  READ TEMPERATURE STRUCTURE FROM OLD MODEL, IF REQUESTED
C     (note: REQUIRES GEOMESH first for interpolation)      
            IF (OLDTEMP) THEN
              CALL READOLDT (hOldMODEL, ND, NDDIM,T,RADIUS,
     >                      TOLD,ROLD,MODOLD,JOBNOLD,TEFF,TEFFOLD,
     >                      TAURCONT, TAURCONTOLD, BTAUR)
              IF (BTWOT) THEN
                CALL READOLDT (hOldMODEL2, ND,NDDIM,T2,RADIUS,
     >                        TOLD2,ROLD2,MODOLD2,JOBNOLD2,TEFF,
     >                        TEFFOLD2, TAURCONT, TAURCONTOLD, BTAUR)
                IF (TMIN > 6000.) THEN
                  TMIN2 = TMIN
                ELSE
                  TMIN2 = 6000.
                ENDIF
                DO L=1, ND
                  TS = T(L)
                  T(L) = T(L) - TFAC*(T(L) - T2(L))
                  IF (T(L) > 1.2 * TS) THEN
                    T(L) = 1.2 * TS
                  ELSE IF (T(L) < 0.8 * TS) THEN
                    T(L) = 0.8 * TS
                  ENDIF
                ENDDO
              ENDIF
            ENDIF

C***  READ TEMPERATURE OR VELOCITY (OR BOTH) FROM FILE 'TABLE'
C***  READING GRADI PREVENTS ODD FEATURES DUE TO NUMERICAL DIFF IN GRADIFF
            IF (TTABLE) THEN
              CALL TABREAD (ND,RADIUS,VELO,GRADI,T,ITAB)
            ENDIF
 
C***  ENTOT = TOTAL NUMBER DENSITY OF ALL ATOMS
            DO L=1,ND
              RL = RADIUS(L)
              IF (.NOT.TTABLE .OR. ITAB < 2) THEN
                IF (NEWVELO) THEN
                  VELO(L) = WRVEL(RL)
                ELSEIF (RL > OLDRADI(1)) THEN
                  VELO(L) = OLDVELO(1)
                ELSE
                  CALL SPLINPOX(VELO(L), RL, OLDVELO, OLDRADI, NDv)
                ENDIF
              ENDIF
cc              WRITE (0,*) 'L, R, V ', L, RL, VELO(L)
              RHO(L) = FM / RL / RL / VELO(L) / 1.E5
              ENTOT(L) = RHO(L) / AMU / ATMEAN
cc              WRITE (0,*) 'L, R, ENTOT ', L, RL, ENTOT(L)*RNE(L)
            ENDDO            

            IF (.NOT. TTABLE .OR. ITAB < 2) THEN
              CALL       GRADIFF  (ND,VELO,GRADI,RADIUS)
            ENDIF
C***  Definition of the Depth-dependent Clumping Factor
            IF (iOldStratification > 0 .AND. 
     >              (.NOT. bForceDCUpdate .OR. TAURCONT(ND) < 0.)) THEN
              DO L=1, ND
                IF (RADIUS(L) > RADIUSold(1)) THEN
                  DENSCON(L) = DENSCON_OLD(1)
                  FILLFAC(L) = FILLFAC_OLD(1)
                ELSEIF (RADIUS(L) < RADIUSold(ND)) THEN
                  DENSCON(L) = DENSCON_OLD(ND)
                  FILLFAC(L) = FILLFAC_OLD(ND)
                ELSE                
                  CALL SPLINPOX(DENSCON(L), RADIUS(L), 
     >                          DENSCON_OLD, RADIUSold, NDold)
                  CALL SPLINPOX(FILLFAC(L), RADIUS(L), 
     >                          FILLFAC_OLD, RADIUSold, NDold)
                ENDIF
              ENDDO
            ELSE
c              WRITE (0,*) 'DEBUG before CLUMP_STRUCT ---------------'
c              WRITE (0,*) 'L, RADIUS(L), VELO(L), TAU(L)'
c              DO L=1, ND
c                WRITE (0,'(I3,3(1X,G15.5))') L, RADIUS(L), VELO(L), TAUROSS(L)
c              ENDDO
c              WRITE (0,*) 'DEBUG before CLUMP_STRUCT ---------------'
              CALL CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, 
     >                           VELO, TAURCONT, DENSCON_LINE, 
     >                           RADIUS, T, XMU)
c              WRITE (0,*) 'DEBUG after CLUMP_STRUCT ~~~~~~~~~~~~~~~'
c              WRITE (0,*) 'L, RADIUS(L), VELO(L), TAU(L)'
c              DO L=1, ND
c                WRITE (0,*) L, DENSCON(L)
c              ENDDO
c              WRITE (0,*) 'DEBUG after CLUMP_STRUCT ~~~~~~~~~~~~~~~'
            ENDIF
            CALL VTURB_SETUP (VTURB, VTURB_LINE, ND, 
     >                        VELO, RADIUS, TAURCONT, T, XMU)

C***  ITAB = 2: ONLY TABULATED INPUT OF V(R), I.E. T(R) MUST BE CALCULATED
            IF (ITAB .EQ. 2) TTABLE=.FALSE.
            TEXIST=TTABLE .OR. OLDTEMP

C***  TEMPERATURE STRATIFICATION AND INITIAL POPNUMBERS (LTE)         
            CALL GREY(ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR,
     >                ALPHA,SEXPO,
     >                ADDCON1, ADDCON2, ADDCON3, 
     >                IGAUNT,POPNUM,TAUGREY,R23,TEXIST,NDIM,N,
     >                LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,
     >                KODAT,ABXYZ,NOM,NFIRST,NLAST,NATOM,
     >                EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,
     >                MAXION,MAXATOM,SIGMATHK,EDGEK,
     >                SEXPOK,KONTNUP,KONTLOW,LASTKON,XDATA, DENSCON,
     >                FILLFAC)
            IF (.NOT. bGREYSTART) THEN
C***          improve TAUROSS scale by using OPAROSS instead of using OPAGREY
C               This has an influence on the velocity and the radiation field 
C               (JSTART uses TAUROSS scale, note that R23 is always from GREY)
              CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >                     T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, 
     >                     WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO,
     >                     ADDCON1, ADDCON2, ADDCON3, 
     >                     IGAUNT, NOM, NF,
     >                     XLAMBDA, FWEIGHT,
     >                     TAUTHOM, TAURCONT,
     >                     MAXATOM, MAXION, SIGMATHK, SEXPOK, EDGEK, 
     >                     KODAT, KONTNUP, KONTLOW, LASTKON, 
     >                     DENSCON, FILLFAC, POPMIN
     >        ) 
            ELSE
              DO L=1, ND
                TAURCONT(L) = TAUGREY(L)
              ENDDO
            ENDIF

            DO L=1,ND
              !Added XMU in WRSTART for use in VELTHIN and to avoid problems 
              !  with first steal run if hydro or TAUFIX is enabled
              !  XMU is recalculated here because RNE has been changed by GREY              
              XMU(L) = ATMEAN / (1. + RNE(L))              
              IF (bFULLHYDROSTAT) THEN
                IF (iOLDRAD == 2) THEN
                  IF (RADIUS(L) > RADIUSold(1)) THEN
                    GAMMAL = GAMMARAD(1)
                  ELSEIF (RADIUS(L) < RADIUSold(NDold)) THEN
                    GAMMAL = GAMMARAD(NDold)
                  ELSE
                    CALL SPLINPOX(GAMMAL, RADIUS(L), 
     >                           GAMMARAD, RADIUSold, NDold)
                  ENDIF
                ELSEIF (iOLDRAD == 1 .AND. .NOT. bGREYSTART) THEN
                  IF (L == ND) THEN
                    DTAUCDTAUTH = (TAURCONT(ND) - TAURCONT(ND-1)) /
     >                               (TAUTHOM(ND) - TAUTHOM(ND-1))
                  ELSEIF (L == 1) THEN
                    DTAUCDTAUTH = (TAURCONT(2) - TAURCONT(1)) /
     >                               (TAUTHOM(2) - TAUTHOM(1))
                  ELSE 
                    DTAUCDTAUTH = (TAURCONT(L+1) - TAURCONT(L-1)) /
     >                               (TAUTHOM(L+1) - TAUTHOM(L-1))
                  ENDIF
                  GAMMAL = DTAUCDTAUTH * GEDD
                ELSEIF (RADGAMMASTART >= 0.) THEN
                  GAMMAL = RADGAMMASTART
                ELSE 
                  GAMMAL = GEDD
                ENDIF
                GEFFL(L) = (10.**GLOG) * (1. - GAMMAL)
              ELSE
                !either Thompson only or fixed => no depth-dependent value
                GEFFL(L) = (10.**GLOG) * (1. - GEDD)     
              ENDIF
            ENDDO

            IF (bHYDROSOLVE) THEN
              AMIN = SQRT(RGAS * T(ND) / XMU(ND))  !sound speed at inner boundary
              IF (VMIN*1.E5 >= AMIN) THEN
                !inner boundary values for v must be less than sound speed
                ! otherwise G~ would never cross zero
                VMINhydro = AMIN / SQRT(2.) / 1.E5
                CYCLE tauit
              ENDIF
            ENDIF
            
C***  NEW VELOCITY FIELD FROM INTEGRATION OF HYDROSTATIC EQUATION
            IF ((ITAB .LT. 2) .AND. THIN) THEN
              CALL VELTHIN(T, RADIUS, VELO, ND, RSTAR, RMAX,
     >                  GEFFL, XMU, VTURB, ThinCard, VELOCRITERION)

C***          STORE RADIUS AND VELO FOR WRVEL, BEFORE A NEW RADIUS 
C***            GRID WILL BE DEFINED
              DO L=1, ND
                OLDRADI(L) = RADIUS(L)
                OLDVELO(L) = VELO(L)
              ENDDO
              NDv = ND
              NEWVELO = .FALSE.

              DO L=1, ND
                RL = RADIUS(L)
                RHO(L) = FM / RL / RL / VELO(L) / 1.E5
                ENTOT(L) = RHO(L) / AMU / ATMEAN
              ENDDO

              IF (.NOT. bThinImprove) THEN
                bThinImprove = .TRUE.
                CYCLE thinwind
              ELSE
                EXIT thinwind
              ENDIF
              
              EXIT
            ELSE
              EXIT
            ENDIF
            
          ENDDO thinwind
C         ------- THIN WIND loop ends here ----------------------


C***  AUTOMATICAL ADJUSTMENT OF VMIN (OPTIONAL: IF TAUMAX SPECIFIED)
          IF (iOldStratification > 0 .AND. (fHYDROSTART > 0.)) THEN
            XMDOT = LOG10(FM) - 2. * LOG10(RSUN/RSTAR) - 3.02
            IF (TAUMAX > 0.) THEN
              
              WRITE (hCPR,FMT='(A,I3,A,F8.4,A,F9.4,A)') 
     >          ' **  TAUMAX ITERATION: ITTAU=', ITTAU, 
     >          '   TAUROSS=', TAURCONT(ND), 
     >          '  LOG MDOT = ', XMDOT, ' [Msun/yr]'
              
              STEPDAMP = 1./ ( 1 +  LOG10( FLOAT(ITTAU) ) )
              
              IF (TAURCONT(ND) > TAUMAX+TAUACC) THEN
                FM = FM * 10**(-0.1 * STEPDAMP)
                IF (ITTAU < MAXITTAU) CYCLE tauit
              ENDIF

              IF (TAURCONT(ND) < TAUMAX-TAUACC) THEN
                FM = FM * 10**(0.1 * STEPDAMP)
                IF (ITTAU < MAXITTAU) CYCLE tauit
              ENDIF            
           ENDIF
            
           EXIT tauit 
            
          ELSEIF (iOldStratification > 0 .AND. 
     >         (.NOT. bOVTauMax) .AND. (.NOT. bOVTauCut)) THEN
c          ELSEIF (iOldStratification > 0) THEN
C***        Use old RCON in case of OLD STRATIFICATION 
C***        if inside of both, old and new grid
            IF (RCONold < RADIUSold(1) .AND. RCONold < RADIUS(1)) THEN
              RCON = RCONold
            ENDIF          
            WRITE (hCPR,FMT='(A,F10.4,A,F11.6)') 
     >          ' *** OLD STRATIFICATION USED: TAURCONT=', TAURCONT(ND), 
     >          '   VMIN=', VMIN
            EXIT tauit
          ELSEIF (TAUMAX > .0) THEN
            IF (bOVTauMax .AND. ITTAU == 1) THEN
              WRITE (hCPR,FMT='(A,F10.4)') 
     >          ' *** OLD STRATIFICATION USED: will be'
     >          // ' shifted to match the desired TAUMAX of ', TAUMAX
              IF (bForceDCUpdate) THEN
C***            We need at least 2 iterations if we update the DENSCON-Scale
C***            and iterate on TAUMAX
                CYCLE
              ENDIF
            ENDIF

            IF (bOVTauCut) THEN
C***          OLD V TAUCUT option: Only use the part of OLD V 
C***                               which is in the Tau-scale of the new model
              IF ( TAURCONT(ND) < (TAUMAX - TAUACC) .AND. 
     >              (ABS(TAURCONT(ND))-TROLD) > TAUACC ) THEN
C***            Common block for velo field has not to be filled with the new
C***            radius grid and the v(r)-part old the old velocity field up to
C***            the new TAURCONT(ND)
C***            @todo: sensible iteration system (derived TAU should not change) 
                NDv = ND
                DO L=1, NDv
                  CALL SPLINPOX(VL, TAURCONT(L), VELOold, TAURCONTOLD, NDold)
                  OLDVELO(L) = VL
                  OLDRADI(L) = RADIUS(L)
                ENDDO          
                VMIN = OLDVELO(NDv)
                TROLD = TAURCONT(ND)
                WRITE (hCPR,FMT='(A,I3,A,F10.4)') 
     >            ' *** OLD V TAUCUT Iteration: IT=', ITTAU, 
     >            '  TAURCONT=', TAURCONT(ND) 
                CYCLE      
              ELSE
                WRITE (hCPR,FMT='(A,F10.4,A,F11.6)') 
     >            ' *** OLD V WITH CUT AT: TAURCONT=', TAURCONT(ND), 
     >            '   VMIN=', VMIN
                EXIT tauit
              ENDIF
            ENDIF

            IF (VMIN == VMINhydro) THEN
              TAUMAX = REAL(CEILING(TAURCONT(ND)))   !taumax = next larger integer 
              WRITE (hCPR,*) '*** TAUMAX has been adjusted for hydro'
              WRITE (hCPR,FMT='(A,F8.2)') '*** New TAUMAX = ', TAUMAX
            ENDIF
          
            !Enforce monotonic TAUROSS scale
            DO L=2, ND
              IF (TAURCONT(L) < TAURCONT(L-1)) THEN
                TAURCONT(L) = TAURCONT(L-1) + 1.E-10
              ENDIF
            ENDDO

            WRITE (hCPR,FMT=77) ITTAU, TAURCONT(ND), VMIN
   77       FORMAT (' *** TAUMAX ITERATION: ITTAU=', I3,
     >          '   TAURCONT=', F10.4, '  VMIN=', G12.4)

            !MAKE BISECTION IN EVERY THIRD STEP
            VMINOLD2 = VMINOLD
            VMINOLD = VMIN
            TROLD2 = TROLD
            TROLD = TAURCONT(ND)
C***        Check how effective this really is            
            IF ((ABS(TAURCONT(ND))-TAUMAX) > TAUACC  .AND.  
     >        ITTAU/10*10 == ITTAU .AND. ITTAU > 2) THEN
C     >         ITTAU .GT. 2) THEN
C***            VMIN = 0.5 * (VMINOLD + VMINOLD2)
                VMIN = (VMINOLD2*(TROLD-TAUMAX) - VMINOLD * 
     >          (TROLD2-TAUMAX)) /  (TROLD - TROLD2)

                IF (ITTAU < MAXITTAU) CYCLE
            ENDIF

C***        TAUROSS TOO LARGE:
            IF (TAURCONT(ND) > TAUMAX+TAUACC) THEN
              CALL SPLINPOX(VMIN, TAUMAX, VELO, TAURCONT, ND)
              IF (ITTAU < MAXITTAU) CYCLE
            ENDIF

C***        TAUROSS TOO SMALL:
            IF (TAURCONT(ND) < TAUMAX-TAUACC) THEN
              VMIN = VMIN * (TAURCONT(ND) / TAUMAX)**.5
              IF (ITTAU < MAXITTAU) CYCLE
            ENDIF

          ENDIF

          EXIT !Exit the TAU iteration loop if no cycling criteria was met
        
      ENDDO tauit
C     ------- Main TAU iteration loop ends here ----------------------


C***  VELOCITY GRADIENT BY DIFFERENTIATION (final calculation)
      IF (ITAB < 2) THEN
        CALL       GRADIFF  (ND,VELO,GRADI,RADIUS)
      ENDIF

      !Store VELOCRITERION as two-character vector in model file
      VCRIT = '        '
      IF (iOldStratification > 1) THEN
        !if old v is used, take also old VCRIT as starting approach
        ! (note: this might not be exact if taumax changes signifcantly!)
        DO L=1, ND
          VCRIT(L) = VCRITold(L)
        ENDDO
      ELSEIF (iOldStratification == 1) THEN
        DO L=1, ND
          VCRIT(L) = 'IP      '
        ENDDO
      ELSEIF (TTABLE .AND. ITAB > 1) THEN
        DO L=1, ND
          VCRIT(L) = 'TAB     '
        ENDDO
      ELSE     
        DO L=1, ND
          !Note: static identifier "ST" is already set above if used
          SELECTCASE (VELOCRITERION(L))
            CASE ('HYSTINT ')
              VCRIT(L) = 'HS      '
            CASE ('BETA    ')
              VCRIT(L) = 'B       '
            CASE ('2BETA   ')
              VCRIT(L) = '2B      '
            CASE ('SQRT    ')
              VCRIT(L) = 'R       '
            CASE DEFAULT
              IF (RADIUS(L) > RCON) THEN
                IF (BETA2FRACTION > 0.) THEN
                  VCRIT(L) = '2B      '
                ELSE
                  VCRIT(L) = 'B       '
                ENDIF
              ENDIF
          ENDSELECT
        ENDDO
      ENDIF
      
C***  MODEL-FILE: MASS STORAGE FILE IN NAME-INDEX MODE ****************
      CALL OPENMS (3,IDUMMY,IDUMMY,1, IERR)
      IFLAG = -1  !'default' flagging by default ;)
      CALL WRITMS (3,ND,1,         'ND      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RADIUS,ND,    'R       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NP,1,         'NP      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,P,NP,         'P       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,Z,ND*NP,      'Z       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,ENTOT,ND,     'ENTOT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,DENSCON,ND,   'DENSCON ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,TEFF,1,       'TEFF    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,GLOG,1,       'GLOG    ', IFLAG, IDUMMY, IERR)
      IF (bSaveGEFF) THEN
        !Save geff only if this value should be fixed in the following
        !calculations (meaning geff or EddGamma was given in CARDS)
        CALL WRITMS (3,GEFFLOG,1,    'GEFFLOG ', IFLAG, IDUMMY, IERR)
      ELSE
        GFLSAV = -1. * GEFFLOG  !negative value indicating flexible value
        CALL WRITMS (3, GFLSAV,1,    'GEFFLOG ', IFLAG, IDUMMY, IERR)
      ENDIF
      CALL WRITMS (3,GEFFKEY,1,    'GEFFKEY ', IFLAG, IDUMMY, IERR)
      IF (GEddFix > 0) THEN
        CALL WRITMS (3,GEDD,1,       'GEDD    ', IFLAG, IDUMMY, IERR)
      ELSE
        CALL WRITMS (3,-1.0   ,1,    'GEDD    ', IFLAG, IDUMMY, IERR)
      ENDIF
      IF (GEDDRAD > 0. .AND. iOldStratification == 0) THEN
        CALL WRITMS (3,GEDDRAD,1,    'GEDDRAD ', IFLAG, IDUMMY, IERR)
      ENDIF
      CALL WRITMS (3,T,ND,         'T       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VELO,ND,      'VELO    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VMIN,1,       'VMIN    ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,VFINAL,1,     'VFINAL  ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,HSCALE,1,     'HSCALE  ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,BETA,1,       'BETA    ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,BETA2,1,      'BETA2   ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,BETA2FRACTION,1,'B2FRAC  ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,VPAR1,1,      'VPAR1   ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,VPAR2,1,      'VPAR1   ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,VPAR1_2,1,    'VPAR12  ', IFLAG, IDUMMY, IERR)
!       CALL WRITMS (3,VPAR2_2,1,    'VPAR22  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,GRADI,ND,     'GRADI   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,TAURCONT,ND,  'TAURCONT', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RSTAR,1,      'RSTAR   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RCON,1,       'RCON    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,INCRIT,ND,    'INCRIT  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VCRIT, ND,    'VCRIT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMSTAR,1,     'XMSTAR  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMDOT,1,      'XMDOT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RHO,ND,       'RHO     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VDOP,1,       'VDOP    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NF,1,         'NF      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XLAMBDA,NF,   'XLAMBDA ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NF2,1,        'NF2     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XLAMBDA2,NF2, 'XLAMBD2 ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,FWEIGHT,NF,   'FWEIGHT ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,KEY,NF,       'KEY     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,POPNUM,ND*N,  'POPNUM  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,POPNUM,ND*N,  'POPLTE  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,LEVEL, N,     'LEVEL   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RNE,ND,       'RNE     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMU,ND,       'XMU     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,ABXYZ,NATOM,  'ABXYZ   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,MODHEAD,13,   'MODHEAD ', IFLAG, IDUMMY, IERR)
      NEXTK=1
      CALL WRITMS (3,NEXTK,1,      'NEXTK   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,JOBNUM,1,     'JOBNUM  ', IFLAG, IDUMMY, IERR)
      BUFFER6 = 'UNDEF.'
      READ(UNIT=BUFFER6, FMT='(A6)') TOTOUT
      CALL WRITMS (3,TOTOUT,1,     'TOTOUT  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XDATA,MAXXDAT,'XDATA   ', IFLAG, IDUMMY, IERR)
      DO L=1, ND
        VMIC(L) = VTURB(L) * SQRT(2.)
      ENDDO
      CALL WRITMS (3,VMIC, ND,     'VMIC    ', IFLAG, IDUMMY, IERR)
      IF (THIN) THEN
        DO L=1, ND
          GRSTATIC(L) = 1. - GEFFL(L) / (10.**GLOG)
        ENDDO
        CALL WRITMS (3,GRSTATIC, ND, 'GRSTATIC', IFLAG, IDUMMY, IERR)
      ENDIF
      IF (OLDTEMP .AND. iOldStratification > 0 
     >                          .AND. DTDRIN_OLD > 0.) THEN
        CALL WRITMS (3, DTDRIN_OLD, 1, 'DTDRIN  ', IFLAG, IDUMMY, IERR)
      ENDIF

C***  WRITE THE BEGINNING OF THE MODEL HISTORY      
      MODHIST(9:32) = '/      0. WRSTART       '
      LAST=4    !3 * 8 +  offset (8)
      WRITE(UNIT=BUFFER8, FMT='(A8)') LAST 
      MODHIST(1:8)=BUFFER8 !First entry contains the number of used integers (=chars / 8)


      IF (OLDTEMP) THEN
         WRITE(UNIT=BUFFER144, FMT=2) MODOLD, JOBNOLD
    2    FORMAT ('T(R) FROM OLD MODEL: ',A,'  AFTER JOB NO.',I4,4X)
         CALL ADDHISTENTRY(MODHIST, LAST, MAXHIST, 144, BUFFER144)
         
         IF (BTWOT) THEN
           WRITE (UNIT=BUFFER144, FMT=22) TFAC, MODOLD2, JOBNOLD2
   22      FORMAT ('SECOND (TFAC=',F4.1',) : ',A,
     >             '  AFTER JOB NO.',I4,4X)
           CALL ADDHISTENTRY(MODHIST, LAST, MAXHIST, 144, BUFFER144)
         ENDIF
      ENDIF

      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ', IFLAG, IDUMMY, IERR)
 
C***  START APPROXIMATION FOR THE RADIATION FIELD
      LASTINDALL = LASTIND + NAUTO + LASTFE
C***  prepare the rudimental settings for DRTRANSIT lines
      IF (NAUTO > 0) THEN
         CALL PREP_DRLINES (DRLINES_CARD, NAUTO, KRUDAUT, EAUTO)
C***     WAVENUMBER OF STABILIZING TRANSITIONS
         DO INDDR=1, NAUTO
            LOW=LOWAUTO(INDDR)
            WSTABIL(INDDR) = EION(LOW) - ELEVEL(LOW) + EAUTO(INDDR)
         ENDDO
      ENDIF
      CALL JSTART (NF,XLAMBDA,KEY,ND,RADIUS,T,XJC,XJL,ELEVEL,
     >             N,EINST,NDIM,INDNUP,INDLOW, LASTINDALL, LASTIND,
     >             NAUTO, WSTABIL, R23, TAUGREY, 
     >             LTESTART, BLACKEDGE, bOLDJ, XJCold, XLAMBDAold,
     >             RADIUSold, NFold, NDold, KRUDAUT)
 
      CALL CLOSMS (3, IERR)
C**********************************************************************

C***  PRINTOUT OF THE FUNDAMENTAL STELLAR PARAMETERS
      IF (bFULLHYDROSTAT) THEN
        GEDDPRINT = GEDDRAD
      ELSE
        GEDDPRINT = GEDD
      ENDIF
      CALL PRIPARAM (MODHEAD, TEFF, RSTAR, XMDOT, XLOGL, RTRANS, 
     >        VFINAL, VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG, GEDDPRINT, 
     >        GEddFix, RMAX, XMSTAR, WRTYPE, MASSORIGIN, LRTinput, ND, 
     >        MLRELATION, VTURB(ND), MDOTINPUT)



C***  PRINTOUT OF THE CHEMICAL COMPOSITION
      CALL PRICOMP (NDIM, EINST, N, NCHARG, NOM, NATOM, ABXYZ, ATMASS,
     $              STAGE, NFIRST, NLAST, ELEMENT, SYMBOL, LASTIND,
     $              INDLOW, INDNUP, NAUTO, LOWAUTO,
     $              EAUTO, KONTNUP, KONTLOW, LASTKON, XMASS, KRUDAUT)

C***  PRINTOUT OF THE X-RAY DATA
      CALL PRIXDAT(XDATA,MAXXDAT)

C***  PRINTOUT OF VARIOUS MODEL SPECIFICATIONS
      bNoDetails = .FALSE.  !print all details at the start
      CALL PRIMOD (ND,RADIUS,INCRIT,ENTOT,T,VELO,GRADI,NP,OLDTEMP,
     $             MODHEAD,JOBNUM,MODOLD,JOBNOLD,TTABLE,TAURCONT,R23,
     $             TEFFOLD,THIN, ITTAU, MAXITTAU, RCON, BTWOT, MODOLD2, 
     >             JOBNOLD2, TFAC, BETA, VPAR1, VPAR2, DENSCON,
     >             BETA2, BETA2FRACTION, HSCALE,
     >             bNoDetails,.TRUE., VTURB(ND))
 
C***  PLOT OF THE VELOCITY LAW
      IF (VPLOT) CALL PLOTV (ND,RADIUS,VELO,MODHEAD,JOBNUM)
      

      CALL JSYMSET ('G0','0')

      CALL STAMP (OPSYS, 'WRSTART', TIM1)

      STOP 'O.K.'
      END
      FUNCTION WRVEL(R)
C***********************************************************************
C***  VELOCITY FIELD
C***  NEWVELO = TRUE: ANALYTIC LAW;
C***  NEWVELO = FALSE: BELOW RCON: INTERPOLATION IN TABLE OLDVELO
C***            OVER OLDRADI
C***            (NEWVELO is called after VELTHIN routine has been passed)
C***  BETA .LE. 0: SWITCH TO SQRT-LOG-LAW
C***********************************************************************
      
      REAL, INTENT(IN) :: R
      INTEGER :: NDVAL

C***  COMMON/VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS TO FUNCTION WRVEL
C***  The second line has additional parameters for the 2-beta-law
C***     -- wrh  6-Apr-2006 17:21:57
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  (OLD) RADIUS AND VELOCITY ARE TRANSFERRED TO FUNCTION WRVEL 
C***     BY SPECIAL COMMON BLOCKS:
      COMMON /COMRADI/ OLDRADI(2)
      COMMON /COMVELO/ NEWVELO, NDVAL, OLDVELO(2)
      LOGICAL NEWVELO
      
      REAL, EXTERNAL :: VELOBETA

      IF (R .LE. 1.) THEN
         WRVEL=VMIN
      ELSE IF (R. LT. RCON) THEN
         IF (NEWVELO) THEN
            WRVEL = (VMIN) * EXP(MIN(100.,(R-1.)/HSCALE))
         ELSE
           IF (R < OLDRADI(NDVAL)) THEN
             WRVEL = OLDVELO(NDVAL)
           ELSEIF (R > OLDRADI(1)) THEN
             WRVEL = OLDVELO(1)
           ELSE
             CALL SPLINPO(WRVEL, R, OLDVELO, OLDRADI, NDVAL)
           ENDIF
         ENDIF
      ELSE
C***    Calculate velocity via (double-)beta law      
        WRVEL = VELOBETA(R)
      ENDIF


      IF (WRVEL.GT.VFINAL) WRVEL=VFINAL

      RETURN

      END
      SUBROUTINE XRUDI (XJ,WAVENUM,XJC,XLAMBDA,ND,NF,L)
C***********************************************************************
C***  INTERPOLATION OF THE CONTINUUM RADIATION FIELD AT WAVENUM
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE
C***********************************************************************
      DIMENSION XJC(ND,NF),XLAMBDA(NF)

      WLENG=1.E8/WAVENUM
      NA=1
      A=XLAMBDA(1)
      NB=NF
      B=XLAMBDA(NF)
      IF ((WLENG-A)*(WLENG-B) .GT. .0) THEN
         WRITE (0,*) 'RUDIMENTAL TRANSITION OUTSIDE WAVELENGTH GRID',
     >               ' AT ', WLENG, ' ANGSTROEM'
         WRITE (0,*) '*** FATAL ERROR in SUBROUTINE XRUDI'
         STOP '*** FATAL ERROR in SUBROUTINE XRUDI'
         ENDIF
   10 IF ( NB-NA .EQ. 1) GOTO 12
      NH=(NA+NB)/2
      H=XLAMBDA(NH)
      IF ((WLENG-A)*(WLENG-H) .GT. .0) GOTO 13
      NB=NH
      B=H
      GOTO 10
   13 NA=NH
      A=H
      GOTO 10
   12 P=(WLENG-A)/(B-A)
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE
      TRADA=TRADFUN(XLAMBDA(NA),XJC(L,NA))
      TRADB=TRADFUN(XLAMBDA(NB),XJC(L,NB))
      TRAD=P*TRADB+(1.-P)*TRADA
      XJ=BNUE(WLENG,TRAD)
      RETURN
      END
