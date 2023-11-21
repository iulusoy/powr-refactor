      PROGRAM MAINformal 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:25 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL formal 
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
      SUBROUTINE BACKJC (XJC, ND, NF, XJCIND, XLAMBDA, 
     >                   XLAMREF, RADIUS)

C***  CALCULATE THE CONTINUUM INTENSITIES AT THE LINE CENTER FREQUENCIES
C***    BY INTERPOLATION IN THE CURRENT MODEL-FILE DATA, I.E. XJC AS OB
C***    BY THE LAST EXECUTION OF "COMO"
C***    NOTE: EVEN IF THE SUBROUTINE ELIMIN HAS BEEN CALLED, THE
C***          VECTOR XJCIND IS NOW OVERWRITTEN

      DIMENSION XJC(ND,NF), XLAMBDA(NF), XJCIND(ND)
      DIMENSION RADIUS(ND)

C***  FIND THE CONTINUUM FREQUENCY INTERVAL OF PRESENT LINE
      KB = ISRCHFGT(NF,XLAMBDA,1,XLAMREF)
      KA = KB - 1
      IF (KB .EQ. 1 .OR. KB .EQ. NF) THEN
        KA = KB
      ELSE IF (KB .LT. 1 .OR. KB .GT. NF) THEN
        WRITE (0,*) 'SUBR. BACKJC: UNEXPECTED LINE FREQUENCY'
        WRITE (0,*) 'XLAMREF=',XLAMREF
        WRITE (0,*) 'BOUNDARIES: ', XLAMBDA(1), XLAMBDA(NF)
        STOP 'ERROR'
      ENDIF
C***  CALCULATE INTERPOLATION WEIGHTS
      XLAMA = XLAMBDA(KA)
      XLAMB = XLAMBDA(KB)
      IF (XLAMB .NE. XLAMA) THEN
        Q = (XLAMREF - XLAMA) / (XLAMB-XLAMA)
      ELSE
        Q = 0.
      ENDIF

C***  INTERPOLATION AT ALL DEPTH POINTS
C***   NOTE: RADIUS**2 * XJC IS EXPRESSED AS TRAD
      DO 1 L=1, ND

        R2 = RADIUS(L) * RADIUS(L)

        XJCA = XJC(L,KA) * R2
        XJCB = XJC(L,KB) * R2

        TRADA = TRADFUN(XLAMA, XJCA)
        TRADB = TRADFUN(XLAMB, XJCB)
        TRAD = (1.-Q) * TRADA + Q * TRADB
        XJCIND(L) = BNUE(XLAMREF,TRAD) / R2

c        if (xjcind(l) .eq. 0) then
c          write(*,*) 'backjc:',xjca,xjcb,trada,
c     >                tradb,xjcind(l)
c        endif

    1 CONTINUE

      RETURN
      END
      SUBROUTINE BANDWIDTH (ND, RADIUS, OPAC, OPAL, LINPRO, AVOIGT, 
     >                      NDDIM,NBL, MAXLAP, XMAX, XMAXBROAD, XMAXLIN,
     >                      PHITAB, NFDIMPHITAB, 
     >                      DXMAX, LEVELNUP, LEVELLOW,
     >                      BDD_VDOP, GRIEMPAR, VDOP, ALN, bHYDROGEN,
     >                      TAUMAX, TAUMINBROAD)
C***********************************************************************
C***  Called from: STARKBROAD
C***  This subroutine estimates the necessary line bandwidth
C***  such that the *static* optical depth drops under a specified limit 
C***********************************************************************
      INTEGER, INTENT(IN) :: NFDIMPHITAB, ND, NBL, MAXLAP
      REAL, INTENT(IN) :: TAUMINBROAD

      REAL, DIMENSION(-NFDIMPHITAB:NFDIMPHITAB, ND) :: PHITAB
      REAL, DIMENSION(MAXLAP,NDDIM) :: AVOIGT, GRIEMPAR
      REAL, DIMENSION(ND) :: RADIUS, OPAC, OPAL
      CHARACTER(8) :: LINPRO
      CHARACTER(10) :: LEVELNUP, LEVELLOW
      LOGICAL :: BDD_VDOP, bHYDROGEN

      REAL, EXTERNAL :: PHIHOLTSMARK
      
      REAL, PARAMETER :: CLIGHT  = 2.99792458E5    ! c in km / s
      
C***  Optical Depth limit: the smaller this value, the larger the bandwidth 
C***  (Jan 2017: Default now in formal.f, can be adjusted via FORMAL_CARDS)
      TAUCRIT = TAUMINBROAD

C***  Branch for pressure-broadening  (depth-dependent VOIGT profiles) 
      IF (     LINPRO .EQ. 'BRD-HeI ' .OR.
     >         LINPRO .EQ. 'VOIGT   ' .OR.
     >         LINPRO .EQ. 'Q-STARK ') THEN     

C***     TAU = integral OPAL * PHI(X) * dr
C***     the Voigt profile is radius-dependent (i.e. inside the sum)

         X = XMAX

    3    CONTINUE
         PHI = VOIGTH(AVOIGT(NBL,1), X)
         TAU  = 0.5 * (RADIUS(1) - RADIUS(2)) * OPAL(1) * PHI
         TAUC = 0.5 * (RADIUS(1) - RADIUS(2)) * OPAC(1)

         DO L=2, ND-1
            PHI = VOIGTH(AVOIGT(NBL,L), X)
            TAU = TAU + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
            TAUC = TAUC + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAC(L)
C***        Integrate only up to a point 
C***        where the continuum opacity exceeds the specified TAUMAX     
            IF (TAUC > TAUMAX) EXIT
         ENDDO

         IF (TAUC < TAUMAX) THEN
           PHI = VOIGTH(AVOIGT(NBL,ND), X) 
           TND = 0.5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHI
           TAU = TAU + TND
         ENDIF
         
         IF (TAU .GT. TAUCRIT) THEN
            X = X + DXMAX
            GOTO 3
         ENDIF

         XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)

C***  Branch for pressure-broadening  (depth-dependent Holtsmark profiles) 
      ELSE IF (LINPRO == 'L-STARK ') THEN     
         
         X = XMAX
         DLAMDOPREL = (VDOP / CLIGHT)
         
    4    CONTINUE
         
C***     Approximative treatment for linear stark broadening of hydrogenic ions
C***     using a joined profile function of a Doppler and a Holtsmark profile
C***     (Note: This is inferior to precise tables!)
C***     Details:
C***     - The parameter GRIEMPAR has been precalculated in LINSTARK
C***     - Normalized profile function is calculated in PHIHOLTSMARK
C***     - BETA and BETADOP need to be provided in Angstroem/cm, which is
C***         guaranteed by the definition of GRIEMPAR in LINSTARK
         BETA = GRIEMPAR(NBL,1) * (EXP(ALN*X) - 1.)
         BETADOP = GRIEMPAR(NBL,1) * DLAMDOPREL
         PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
                  
         TAU  = 0.5 * (RADIUS(1)   - RADIUS( 2)) * OPAL( 1) * PHI
         TAUC = 0.5 * (RADIUS(1)   - RADIUS( 2)) * OPAC( 1)

         DO L=2, ND-1
            BETA = GRIEMPAR(NBL,L) * (EXP(ALN*X) - 1.)
            BETADOP = GRIEMPAR(NBL,L) * DLAMDOPREL
            PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
            TAU = TAU + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
C***        Integrate only up to a point 
C***        where the continuum opacity exceeds the specified TAUMAX     
            TAUC = TAUC + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
            IF (TAUC > TAUMAX) EXIT
         ENDDO

         IF (TAUC < TAUMAX) THEN
            BETA = GRIEMPAR(NBL,ND) * (EXP(ALN*X) - 1.)
            BETADOP = GRIEMPAR(NBL,ND) * DLAMDOPREL
            PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
            TND = 0.5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHI
            TAU = TAU + TND
         ENDIF
         
         IF (TAU > TAUCRIT) THEN
            X = X + DXMAX
            IF (X < 1./DLAMDOPREL) THEN
C***          Increase BANDWIDTH if X < C/VDOP
              GOTO 4
            ELSE 
C***          Cut profile, but write warning into CPR and OUT            
              WRITE (0,'(4A)') 'WARNING: BANDWIDTH CANNOT'
     >          //  ' FULLY COVER L-STARK PROFILE FOR ',              
     >                    LEVELNUP, ' - ', LEVELLOW
              WRITE (*,'(4A)') 'WARNING: BANDWIDTH CANNOT'
     >          //  ' FULLY COVER L-STARK PROFILE FOR ',              
     >                    LEVELNUP, ' - ', LEVELLOW
            ENDIF
         ENDIF
         
         XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)
      
C***  Branch for tabulated STARK profiles (H I and He II)
      ELSE IF (LINPRO(:3) .EQ. 'BRD') THEN

         KMIN = INT(XMAX/DXMAX)

C***     Note: only the bue wing is tested
         DO K=KMIN, NFDIMPHITAB
C***        TAU = integral OPAL * PHI(X,r) * dr
            X = K * DXMAX

            TAU = .5 * (RADIUS(1) - RADIUS(2)) * OPAL(1) * PHITAB(K,1)
            TAUC= .5 * (RADIUS(1) - RADIUS(2)) * OPAC(1)
            DO L=2, ND-1
               TAU = TAU + 0.5 * (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L)
     >             * PHITAB(K,L)
               TAUC= TAUC + 0.5 * (RADIUS(L-1) - RADIUS(L+1)) * OPAC(L)
               IF (TAUC > TAUMAX) EXIT
            ENDDO
            IF (TAUC < TAUMAX)
     >        TAU = TAU +
     >        .5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHITAB(K,ND)
            IF (TAU .LT. TAUCRIT) GOTO 2

         ENDDO

C***     Print a warning if TAU is still larger than TAUCRIT at the end 
C***       of the tabulated profile
         WRITE (0,'(5A)') 'WARNING: PROFILE NOT FULLY COVERED BY ',
     >                    'TABLE DIM. NFDIMPHITAB: ', 
     >                    LEVELNUP, ' - ', LEVELLOW

C***     XMAX is raised to the frequecy which was required 
C***       to reach TAU < TAUCRIT (at max NFDIMPHITAB*DXMAX) 
    2    XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)

C***  Error branch
      ELSE
         WRITE (0,*) '*** UNKNOWN PROFILE TYPE: ', LINPRO
         STOP ' *** FATAL ERROR IN SUBROUTINE BANDWIDTH'
      ENDIF


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
      SUBROUTINE CMFFEOP (XLAMK, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDEXMAX, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT,
     >                    SIGMAACT, OPAFEI, ETAFEI,
     >                    OPAFEION, ETAFEION, T, IVERS_FE_EXPFAC, 
     >                    TEFF, NCHARG, MAXION, bFELASER, bNoIronLaser)

C***********************************************************************
C***  NON-LTE IRON OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  CALLED FROM: COLI and FORMAL -> FORMCMF
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, N, MAXION, LASTFE, INDEXMAX,
     >                       IVERS_FE_EXPFAC, MAXFEACT
      REAL, INTENT(IN) :: TEFF, RSTAR, CLIGHT, XLAMK, XLAM0FE, DXFE, 
     >                    VDOPFE
      LOGICAL, INTENT(IN) :: bNoIronLaser

      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      INTEGER, DIMENSION(LASTFE) :: IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                              INDFEACT, INDRB
      INTEGER, DIMENSION(N) :: IFRBSSTA, IFRBSEND, INDFESACT,
     >                         INDRBS, IFES
      REAL, DIMENSION(ND) :: OPAFE, ETAFE, ENTOT, T
      REAL, DIMENSION(ND, MAXION) :: OPAFEION, ETAFEION
      REAL, DIMENSION(ND,N) :: POPNUM, POPLTE
      REAL, DIMENSION(LASTFE) :: SIGMAACT
      REAL, DIMENSION(ND,LASTFE) :: OPAFEI, ETAFEI
      REAL, DIMENSION(N) :: ELEVEL, WEIGHT
      LOGICAL, DIMENSION(ND) :: bFELASER
      INTEGER, DIMENSION(N) :: NCHARG

      REAL :: SIGMA, SIGMAR, PLASER, ENUP, ENLOW, SUM, G, EMINDU,
     >        DELTA, SIGMA1, SIGMA2, WAVNUM0, WAVNUM03, XINDF,
     >        WEIGHTLU, WAVNUMK, WAVNUMK3
      INTEGER :: L, ION, INDF, INDFSTART, NUP, LOW, IND, INDACT
            
     
C***  Numerical parameters     
      REAL, PARAMETER :: EXPMAX = 500.
     
C***  Physical constants
      REAL, PARAMETER :: C1 = 1.4388        !C1 = H * C / K    ( CM * ANGSTROEM )
      REAL, PARAMETER :: C2 = 3.9724E-16    !C2 = 2 * H * C     ( CGS UNITS )

C***  CALCULATE FE-FREQUECY-INDEX
      XINDF = - LOG10(XLAM0FE/XLAMK) /
     >          LOG10(1. + VDOPFE*DXFE/CLIGHT)

C***  WAVENUMBERS IN KAYSER ( = CM**-1)
      WAVNUMK = 1.E8 / XLAMK
      WAVNUMK3 = WAVNUMK * WAVNUMK * WAVNUMK 

C***  PRESET OPAFE, ETAFE
      DO L=1, ND
         OPAFE(L) = 0.
         ETAFE(L) = 0.
         DO ION=1, MAXION
           OPAFEION(L, ION) = 0.
           ETAFEION(L, ION) = 0.
         ENDDO
      ENDDO
      
C***  LOOP OVER ACTIVE FE-LINES
      DO 10 INDACT=1, MAXFEACT

         IND = INDFEACT(INDACT)
         LOW = IFELOW(IND)
         NUP = IFENUP(IND)
         ION = NCHARG(LOW) + 1
         WEIGHTLU = WEIGHT(LOW) / WEIGHT(NUP)
         WAVNUM0 = ELEVEL(NUP) - ELEVEL(LOW)
         WAVNUM03 = WAVNUM0 * WAVNUM0 * WAVNUM0 

C***    FIND INDICES IN ARRAY 'SIGMAFE'
         INDF = NINT(XINDF-0.5)-IFRBSTA(IND)
         INDFSTART = INDRB(IND)
C***    READ NEIGHBOURING VALUES
         SIGMA1 = SIGMAFE(INDFSTART + INDF)
         SIGMA2 = SIGMAFE(INDFSTART + INDF + 1)
C***    INTERPOLATION
         DELTA = XINDF - FLOAT(INDF+IFRBSTA(IND))
         SIGMA = (1.-DELTA)*SIGMA1 + DELTA*SIGMA2

C***    MULTIPLY WITH RSTAR (DIMENIONS OF OPA AND ETA)
         SIGMAR = SIGMA * RSTAR

ccc     For test use: if "nu-factors" are activated:
ccc     modify SIGMA --> used in FREQUINT for integrartion of XJFEMEAN
ccc     -- the following statement might be de-activated else
ccc         SIGMA = SIGMA * WAVNUM0 / WAVNUMK
ccc    Note!!! SIGMAACT is NOT MANIPULATED in the Version with correct RRATE 
ccc            integration!!

C***    STORE SIGMA FOR LATER USE
         SIGMAACT(INDACT) = SIGMA

         IF (SIGMA <= 0.) CYCLE
         
C***    CALCULATE 'OPA' AND 'ETA' OVER DEPTH INDEX 'L' 
         DO 20 L=1, ND
            ENUP  = POPNUM(L,NUP) * ENTOT(L)
            ENLOW = POPNUM(L,LOW) * ENTOT(L)
C***       Opacity and Emissivity
C***     Different versions according to the IRONLINES-EXPFAC option:
C***     0: OFF  - no exp factor (recommended by Goetz) 
C***     1: TEFF - exp factor with max(Te,Teff) (recommended by wrh, default) 
C***     2: TRADITIONAL - exp factor with Te (standard till 17-Feb-2018)
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              G = WEIGHTLU 
            ELSEIF (IVERS_FE_EXPFAC .EQ. 1) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/MAX(T(L),TEFF))
            ELSEIF (IVERS_FE_EXPFAC .EQ. 2) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/T(L))
            ELSE 
              STOP '*** FATAL INTERNAL ERROR 1 in subr. CMFFEOP'
            ENDIF
            EMINDU = G * ENUP * SIGMAR
            SUM = ENLOW * SIGMAR - EMINDU
            OPAFEI(L, IND) = SUM
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUM03
            ELSE
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUMK3
            ENDIF
ccc         in the preceding stament WAVNUM03 has been used by Goetz

            PLASER = ENUP/WEIGHT(NUP) - ENLOW/WEIGHT(LOW)
            IF (bNoIronLaser .AND. PLASER > 0.) THEN
              !skip level
            ELSEIF (OPAFEI(L,IND) .GT. 0.) THEN
               OPAFE(L) = OPAFE(L) + OPAFEI(L,IND)
               ETAFE(L) = ETAFE(L) + ETAFEI(L,IND)
               OPAFEION(L,ION) = OPAFEION(L,ION) + OPAFEI(L,IND) 
               ETAFEION(L,ION) = ETAFEION(L,ION) + ETAFEI(L,IND) 
                              
              IF (PLASER > 0. .AND. (.NOT. bFELASER(L))) THEN
c                WRITE (hCPR,*) 'FE POP INVERSION AT L = ', L
                bFELASER(L) = .TRUE.
              ENDIF
            
            ENDIF
            
 20      CONTINUE

 10   CONTINUE

      RETURN
      END
      SUBROUTINE CMFSET_FORMAL(Z,
     $          ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,ETA,
     $          PP,BCORE,DBDR,RADIUS,XIMINUS,DXI,OPAK,ETAK,
     >          OPAFE, ETAFE, BWITHLINES)
C***********************************************************************
C***  THIS SUBROUTINE IS TO SET UP THE ARRAY ELEMENTS FOR THE CMF FORMALISM
C***********************************************************************

      DIMENSION S(ND),OPA(ND),ETA(ND),VA(ND),VB(ND)
      DIMENSION TA(ND),TB(ND),TC(ND),UB(ND),GA(ND),H(ND)
      DIMENSION PP(ND),Z(ND),OPAK(ND),ETAK(ND)
      DIMENSION OPAFE(ND), ETAFE(ND)
      LOGICAL BWITHLINES
 
      LZ=LMAX-1
  
C***  OUTER BOUNDARY CONDITION  -  FIRST ORDER
      IF (BWITHLINES) THEN
         AK=OPA(1)+OPAK(1)+OPAFE(1) ! PHIK*OPAL(1)
         AZ=0.5*(AK+OPA(2)+OPAK(2)+OPAFE(2)) ! PHIK*OPAL(2))
      ELSE
         AK=OPA(1)
         AZ=0.5*(AK+OPA(2))
      ENDIF  
      TAUZ=Z(1)-Z(2)
      DX=PP(1)/AK
      S(1)=XIMINUS+DX*DXI
      TC(1)=1./(AK*TAUZ)
      TB(1)=TC(1)+DX+1.
      UB(1)=DX
      VB(1)=0.0
C***  FOR G AND H, THE MATRIX ELEMENT S ARE NOT DIFFERENT FROM INNER POINTS
      DTZM=1./(AZ*TAUZ)
      DXZM=(PP(1)+PP(2))/AZ/2.
      DAZM=DTZM/(1.+DXZM)
      DBZM=DXZM/(1.+DXZM)
      GA(1)=-DAZM
      H(1)=DBZM
      IF(LZ.LT.2) GOTO 2
 
C***  NON-BOUNDARY POINTS
      DO 1 L=2,LZ
         IF (BWITHLINES) THEN
            AK=OPA(L)+OPAK(L)+OPAFE(L) ! PHIK*OPAL(L)
            EK=ETA(L)+ETAK(L)+ETAFE(L) ! PHIK*ETAL(L)
            AZ=0.5*(AK+OPA(L+1)+OPAK(L+1)+OPAFE(L+1)) ! PHIK*OPAL(L+1))
            AZM=0.5*(AK+OPA(L-1)+OPAK(L-1)+OPAFE(L-1)) ! PHIK*OPAL(L-1))
         ELSE
            AK=OPA(L)
            EK=ETA(L) 
            AZ=0.5*(AK+OPA(L+1)) 
            AZM=0.5*(AK+OPA(L-1))
         ENDIF
         S(L)=EK/AK
         TAU=0.5*(Z(L-1)-Z(L+1))
         TAUZ=Z(L)-Z(L+1)
         TAUZM=Z(L-1)-Z(L)
         DT=1./(AK*TAU)
         DTZ=1./(AZ*TAUZ)
         DTZM=1./(AZM*TAUZM)
         DX=PP(L)/AK
         DXZ=(PP(L)+PP(L+1))*0.5/AZ
         DXZM=(PP(L)+PP(L-1))*0.5/AZM
         DAZ=DTZ/(1.+DXZ)
         DAZM=DTZM/(1.+DXZM)
         DBZ=DXZ/(1.+DXZ)
         DBZM=DXZM/(1.+DXZM)
         TA(L)=DT*DAZM
         TC(L)=DT*DAZ
         TB(L)=TA(L)+TC(L)+DX+1.
         UB(L)=DX
         VA(L)=-DT*DBZM
         VB(L)=DT*DBZ
         GA(L)=-DAZ
         H(L)=DBZ
C     DTZM=DTZ
C     DAZM=DAZ
C     DBZM=DBZ
    1 CONTINUE
 
    2 L=LMAX
      IF(LMAX.LT.ND) GOTO 4
      
C***  INNER BOUNDARY CONDITION (CORE RAYS)  -  ONLY TO FIRST ORDER
      IF (BWITHLINES) THEN
         AK=OPA(ND)+OPAK(ND)+OPAFE(ND) ! PHIK*OPAL(ND)
      ELSE
         AK=OPA(ND)
      ENDIF  
      S(ND)=BCORE+DBDR*Z(ND)/AK
      TAUZ=Z(ND-1)-Z(ND)
      DT=1./(TAUZ*AK)
      DX=PP(L)/AK
      TA(L)=DT
      TB(L)=DT+DX+1.
      UB(L)=DX
      VA(L)=0.0
      RETURN
 
C***  INNER BOUNDARY CONDITION (NON-CORE RAYS)  -  SECOND ORDER
    4 CONTINUE
      IF (BWITHLINES) THEN      
         AK=OPA(L)+OPAK(L)+OPAFE(L) ! PHIK*OPAL(LMAX)
         EK=ETA(L)+ETAK(L)+ETAFE(L) ! PHIK*ETAL(L)
      ELSE
         AK=OPA(L)
         EK=ETA(L)
      ENDIF        
      S(L)=EK/AK
      TAUZ=Z(LZ)
      DT=1./(AK*TAUZ)
      DX=PP(L)/AK
      DA=DT/(1.+DX)
      DB=DX/(1.+DX)
      TA(L)=2*DT*DA
      TB(L)=TA(L)+DX+1.
      UB(L)=DX
      VA(L)=-2*DT*DB
      RETURN
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
      SUBROUTINE CONVOLGAUSS_FLEX (YORIG, YCONV, NDAT, DX, VDOPDU, 
     >                             SIGMA_SQRD)
C**********************************************************************
C***  CONVOLUTION OF PROFILE WITH GAUSS-Function exp(-x*x/(sigma^2))
C**********************************************************************
      DIMENSION YORIG(NDAT), YCONV(NDAT)

C**** LOOP OVER ALL DATA POINTS *******************************
      DO L=1, NDAT

C***    CENTRAL POINT
        WEIGHTSUM = 1.
        YCONV(L) = YORIG(L) 

C***    LOOP OVER STEPS TO THE LEFT
        K = L
510     K = K - 1
        X = (L-K) * DX
        IF (K .LT. 1 .OR. X .GT. 4.5*VDOPDU) GOTO 610
        WEIGHT = EXP(-X*X/ SIGMA_SQRD)
        WEIGHTSUM = WEIGHTSUM + WEIGHT
        YCONV(L) = YCONV(L) + YORIG(K) * WEIGHT
        GOTO 510

610     K = L

C***    LOOP OVER STEPS TO THE RIGHT
620     K = K + 1
        X = (K-L) * DX
        IF (K .GT. NDAT .OR. X .GT. 4.5*VDOPDU) GOTO 720
        WEIGHT = EXP(-X*X/ SIGMA_SQRD)
        WEIGHTSUM = WEIGHTSUM + WEIGHT
        YCONV(L) = YCONV(L) + YORIG(K) * WEIGHT
        GOTO 620

C*** NACHNORMIERUNG
720     YCONV(L) = YCONV(L) / WEIGHTSUM

      ENDDO
C*****END OF LOOP ******************************************************

C***  The original date are finally overwritten by the convolved data
      DO L=1, NDAT
         YORIG(L) = YCONV(L)
      ENDDO



      RETURN
      END
      SUBROUTINE CONVOLOPAFE(YORIG, YSCRATCH, NDDIM, ND, NFL, DX, 
     >                       VDOP, VDOPFE, DD_VDOPDUFE)
C**********************************************************************
C***  CONVOLUTION OF PROFILE WITH GAUSS-Function exp(-x*x/(sigma^2))
C***
C***  called from FORMCMF
C**********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: NDDIM, NFL, ND

      REAL, DIMENSION(ND), INTENT(IN) :: DD_VDOPDUFE
      REAL, DIMENSION(NDDIM, NFL), INTENT(INOUT) :: YORIG
      REAL, DIMENSION(NFL), INTENT(INOUT) :: YSCRATCH

      INTEGER, PARAMETER :: KRELMAX = 250000
      INTEGER, PARAMETER :: NDCONVMAX =  200

      REAL, DIMENSION(KRELMAX) :: DOPWEIGHT
      INTEGER, DIMENSION(NDCONVMAX) :: KRELINDEX
      
      REAL, INTENT(IN) :: VDOP, VDOPFE, DX
      
      REAL :: SIGMA_SQRD, WEIGHT, WEIGHTSUM, X
      INTEGER :: L, KL, K, KREL, KR
      
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      IF (ND > NDCONVMAX) THEN
        WRITE (hCPR,'(A)') 'CONVOLOPAFE: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'CONVOLOPAFE: NDHDMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDCONVMAX = ', NDCONVMAX
        STOP 'FATAL ERROR IN FORMAL->CONVOLOPAFE'
      ENDIF
            
C***  LOOP OVER ALL DEPTH POINTS      
      DO L=1, ND
C***    Do not perform any convolution if current intended FE VDOP
C***      is below FEDAT VDOPFE
        IF (DD_VDOPDUFE(L) <= VDOPFE/VDOP) CYCLE
        
        SIGMA_SQRD = (DD_VDOPDUFE(L))**2 - (VDOPFE/VDOP)**2
        
C***    Prepare doppler weight vector        
        X = DX
        KREL = 0
        DO WHILE (X < 4.5*DD_VDOPDUFE(L))
          KREL = KREL + 1
          IF (KREL > KRELMAX) THEN
            WRITE (hCPR,'(A)') 'CONVOLOPAFE: FATAL ERROR ******'
            WRITE (hCPR,'(A)') 'CONVOLOPAFE: KRELMAX INSUFFICIENT'
            STOP 'FATAL ERROR IN FORMAL->CONVOLOPAFE'          
          ENDIF
          
          DOPWEIGHT(KREL) = EXP(-X*X/ SIGMA_SQRD)
          X = X + DX
        ENDDO
C***    Store maximum relative index
        KRELINDEX(L) = KREL
                  
C****   LOOP OVER ALL DATA POINTS **************************************
        DO KL=1, NFL

C***    CENTRAL POINT
          WEIGHTSUM = 1.
          YSCRATCH(KL) = YORIG(L,KL) 

C***      Loop over relative steps to the left          
          DO KR = 1, KRELINDEX(L)
            K = KL - KR
            IF (K < 1) EXIT
            YSCRATCH(KL) = YSCRATCH(KL) + YORIG(L,K) * DOPWEIGHT(KR)
            WEIGHTSUM = WEIGHTSUM + DOPWEIGHT(KR)
          ENDDO

C***      Loop over relative steps to the right          
          DO KR = 1, KRELINDEX(L)          
            K = KL + KR
            IF (K > NFL) EXIT
            YSCRATCH(KL) = YSCRATCH(KL) + YORIG(L,K) * DOPWEIGHT(KR)
            WEIGHTSUM = WEIGHTSUM + DOPWEIGHT(KR)
          ENDDO
          
C***      Normalization of Doppler profile
C***      (This cannot be done with DOPWEIGHT vector as they would
C***       not cover the boundary situations)
          YSCRATCH(KL) = YSCRATCH(KL) / WEIGHTSUM

        ENDDO
C*****  END OF DATAPOINT LOOP ******************************************

C***  The original data is overwritten with the convolved data
        DO KL=1, NFL
          YORIG(L,KL) = YSCRATCH(KL)
        ENDDO

      ENDDO

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
C**********************************************************************
C***  This subroutine copies the specified SECOND MODEL to fort.4
C***     and checks that the DATOM files are identical 
C***  Called from: FORMAL
C**********************************************************************
      SUBROUTINE COPY_SECONDMODEL 
     >           (SECONDMODEL_PATH, IGNORE_DIFF, BIRONLINES) 

      CHARACTER SECONDMODEL_PATH*(*), TESTLINE*80
      LOGICAL IGNORE_DIFF, BIRONLINES, SECONDMODEL_EXISTS

C***  Test if SECONDMODEL_PATH exists, i.e. if output of ls not empty
      CALL SYSTEM  ('echo `ls -d ' // 
     >         SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) //
     >          ' `   > secmod_path ')

      OPEN (44, FILE='secmod_path', STATUS='OLD', ERR=990)
      READ (44, '(A)') TESTLINE
      IF (TESTLINE .NE. '' ) THEN
         WRITE (0,*)
     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' was found'
      ELSE
         WRITE (0,*) '*** ERROR: ' //
     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' not found!'
         GOTO 999
      ENDIF
      CLOSE (44)

ccc   Helge's suggestion with INQUIRE (does not work!??)
ccc      INQUIRE (FILE='SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))', 
ccc     >         EXIST=SECONDMODEL_EXISTS)
ccc      IF (SECONDMODEL_EXISTS) THEN
ccc         WRITE (0,*)
ccc     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' was found'
ccc      ELSE
ccc         WRITE (0,*) '*** ERROR: ' //
ccc     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' not found!'
ccc         GOTO 999
ccc      ENDIF

C***  Note: mass-storage files cannot be opened by name; hence copy first
      CALL SYSTEM ('cp ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/MODEL fort.4')
      CALL SYSTEM ('chmod u+w fort.4')

      WRITE (0,*) 'Copying secondmodel: ', 
     >             SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))

C******************************************************************
C***  Check if DATOM files are identical:
      CALL SYSTEM ('diff ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/DATOM DATOM | wc -c > DATOM.diff')
      OPEN (44, FILE='DATOM.diff', STATUS='OLD', ERR=991)
      READ (44, '(A)') TESTLINE
      CLOSE (44)
      IF (TESTLINE .NE. '0') THEN
         WRITE (0,*) 
     >    '*** ERROR: DATOM files of first and second model differ'
         IF (IGNORE_DIFF) THEN
            WRITE (0,*) 
     >      '*** NO ABORT since option IGNORE_DIFF has been requested'
         ELSE
            WRITE (0,'(A,/,A)') 
     >         '*** If you are sure that this difference is ' 
     >         // 'not significant:', '*** use the option: IGNORE_DIFF'
            GOTO 999
         ENDIF
      ELSE
         WRITE (0,*) 'DATOM files of first and second model agree'
      ENDIF
C*******************************************************************

C******************************************************************
C***  Check if FEDAT_FORMAL files differ:
      IF (BIRONLINES) THEN
       CALL SYSTEM ('diff ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/FEDAT_FORMAL fort.21 | wc -c > FEDAT.diff')
       OPEN (45, FILE='FEDAT.diff', STATUS='OLD', ERR=992)
       READ (45, '(A)') TESTLINE
       CLOSE (45)
       IF (TESTLINE .NE. '0') THEN
         WRITE (0,*) '*** IMPORTANT WARNING:  ********************'
         WRITE (0,*) 
     >     '*** ERROR: FEDAT files of first and second model differ'
         WRITE (0,*) '*** Make sure that both FEDAT files have same ' //
     >                    'superlevel structure!'
         WRITE (0,*) '*** In particular, check for parity splitting!'
         WRITE (0,*) '********************************************'
         IF (IGNORE_DIFF) THEN
            WRITE (0,*) 
     >      '*** NO ABORT since option IGNORE_DIFF has been requested'
         ELSE
            WRITE (0,*) 
            WRITE (0,'(A,/,A)') 
     >         '*** If you are sure that this difference is ' 
     >         // 'not significant:', '*** use the option: IGNORE_DIFF'
            GOTO 999
         ENDIF
       ELSE
         WRITE (0,*) 'FEDAT_FORMAL files of first and second model agree'
       ENDIF
      ENDIF
C*******************************************************************

  100 CONTINUE
      RETURN

C**********************************************************************
C***  ERROR BRANCHES 
C**********************************************************************

  990 WRITE (0,*) '*** Internal ERROR when opening secmod_exist'
      GOTO 999

  991 WRITE (0,*) '*** Internal ERROR when opening DATOM.diff'
      GOTO 999

  992 WRITE (0,*) '*** Internal ERROR when opening FEDAT.diff'
      GOTO 999

  999 STOP '*** FATAL ERROR detected by subr. COPY_SECONDMODEL' 

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
      SUBROUTINE DECFORM (KARTE, LSOPA, RANGERED, RANGEBLUE,
     $                   VSINI,DISP, REDIS, BWESEX, 
     $                   LSPRO, LSDWL, FIN, VDOP, IDTERM, 
     $                   JPFIRST_ORIG,JPLAST_ORIG,IVERSION,
     >                   TAUMAX, XMAX, DXMAX, PATH_VCSSB,
     $                   IDENT,OSMIN,BROAD, MODHEAD,JOBNUM,
     >                   STRING1,MAXSTRI,NSTRING,ABSWAV,BCONT,FUNIT, 
     >                   LINELIST,
     >                 BIRONLINES, BNOCONT, SPZ1, SPZ2, MANIPOP_OPTIONS, 
     >                 XUNIT, BCALIBRATED, FREQIN, MACROCLUMPLINE, 
     >                 PATH_LEMKE_DAT,BAIRWAVELENGTHSET,
     >                 BAIRWAVELENGTH, MAXSPZ, TRANSDWLLINE, RCOROTLINE,
     >                 DD_VDOP_LINE, BPLOTVDOP,
     >                 BMICROTURB, LIMB_LINE,
     >                 VDOPPLOT_LINE, bBIGBANDLIMIT, 
     >                 NOWIND_LINE, TAUMINBROAD, 
     >                 bDDOPAFORMCMF, bDDFECONVOL, 
     >                 LPHISTA_ORIG, LPHIEND_ORIG,
     >                 BVSINI_AT_RCOROT, DX_WINDROT, SECONDMODEL_LINE)
C***********************************************************************
C***  DECODING INPUT CARDS FOR PROGRAM "FORMAL" FROM FTO5 = $IN (DEFAULT)
C***  THIS ROUTINE IS LEFT WHEN A "LINE" OR "BLEND" OPTION IS FOUND,
C***  AND ENTERED AGAIN WHEN THE LINE HAS BEEN COMPUTED.
C***********************************************************************
 
      LOGICAL FIN, REDIS, IDENT, BMICROTURB
      LOGICAL BVSINI_AT_RCOROT
      LOGICAL BROAD, BAIRWAVELENGTHSET,BAIRWAVELENGTH
      LOGICAL ABSWAV, BCONT, LINELIST, BIRONLINES, BNOCONT, BCALIBRATED
      LOGICAL :: BTRANSVELO, BTRANSOD,  BPLOTVDOP, 
     >           bBIGBANDLIMIT, bDDOPAFORMCMF, bDDFECONVOL
      CHARACTER*132 STRING1(0:MAXSTRI), FUNIT*(*)
      CHARACTER KARTE*(*), ACTPAR*80, MODHEAD*100, MACROCLUMPLINE*(*)
      CHARACTER ACTPAR2*80
      CHARACTER TRANSDWLLINE*(*), RCOROTLINE*(*), DD_VDOP_LINE*(*)
      CHARACTER LIMB_LINE*(*), VDOPPLOT_LINE*(*), NOWIND_LINE*(*)
      CHARACTER*20 FREQIN 
      CHARACTER*256 PATH_VCSSB, PATH_LEMKE_DAT
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ), MANIPOP_OPTIONS, XUNIT
      CHARACTER*(*) SECONDMODEL_LINE

      REAL :: TAUMINBROAD
      
C***  Set Defaults
      DO I=1, MAXSPZ
         SPZ1(I) = ''
         SPZ2(I) = ''
      ENDDO
      NSPZ = 0
      RANGERED  = .0
      RANGEBLUE = .0

      JPFIRST_ORIG = 0
      JPLAST_ORIG  = 999

      LPHISTA_ORIG = 0
      LPHIEND_ORIG = 999

      SECONDMODEL_LINE = ''

C***  BEFORE ANY NEW DECODED 'LINE'-CARD OR 'BLEND'-BLOCK:
C***  FIRST, CLEAR THE STRING ARRAY FOR COMMENTS WHICH ARE TO BE 
C***  WRITTEN TO FILE 'PLOT' (KANAL = 1)
      NSTRING = 0
      DO 11 NSTRI=0,MAXSTRI
   11 STRING1(NSTRI) = ' '

    1 READ (2,2, IOSTAT=IOSTAT) KARTE
    2 FORMAT (A)

C*** END OF INPUT DATA REACHED
      IF (IOSTAT .LT. 0 .OR. KARTE(:3) .EQ. 'END') THEN
            FIN=.TRUE.
            RETURN
            ENDIF

      CALL SARGC(KARTE,NPAR)
      IF (NPAR. LE. 0) GOTO 1
      IF (KARTE(1:1) .EQ. '*') GOTO 1
      CALL SARGV (KARTE, 1, ACTPAR)

C*** NOTE: 'LINE'-CARDS OR 'BLEND'-BLOCKS ARE FURTHER DECODED IN SUBR. PREFORM
      IF ( KARTE(:4) .EQ. 'LINE' .OR. KARTE(:5) .EQ. 'BLEND') RETURN
C                          ====                       =====

      IF (ACTPAR .EQ. 'RANGE') THEN
C                          =====
         IF (NPAR .LT. 3) THEN
            WRITE (0,*) '*** ERROR: RANGE option needs two parameters'
            GOTO 99
         ENDIF 
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) RANGEBLUE
         CALL SARGV (KARTE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) RANGERED
         IF (RANGEBLUE .EQ. RANGERED) THEN
            WRITE (0,*) '*** ERROR: RANGE has zero length'
            GOTO 99
         ELSEIF (RANGEBLUE .GT. RANGERED) THEN
            SCRATCH   = RANGEBLUE
            RANGEBLUE = RANGERED
            RANGERED  = SCRATCH
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'PLOT' ) THEN
C                               ====
            CALL SARGV(KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'VDOP') THEN
               BPLOTVDOP = .TRUE.
               VDOPPLOT_LINE = KARTE
            ENDIF

      ELSE IF (ACTPAR(:4) .EQ. 'LIMB') THEN
C                               ====
            CALL SARGV(KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'OFF') THEN
               LIMB_LINE = "OFF"
            ELSE
               LIMB_LINE = KARTE
            ENDIF

      ELSE IF (KARTE(:9) .EQ. 'TRANS DWL') THEN
C                              =========
            IF (LSDWL.LT.0) LSDWL=2
            IF (LSDWL.EQ.1) LSDWL=4            
C** TRANSDWLLINE holds TRANS DWL x-axis options. See routine "tradwl" for description.           
            TRANSDWLLINE = KARTE
            
      ELSE IF ( KARTE(:9) .EQ. 'PRINT DWL') THEN
C                               ========
            IF (LSDWL.LT.0) LSDWL=1
            IF (LSDWL.EQ.2) LSDWL=4

      ELSE IF ( ACTPAR .EQ. 'VSINI' ) THEN
C                            =====
            CALL SARGV (KARTE,2,ACTPAR2)
            READ (ACTPAR2, '(F20.0)', ERR=98) VSINI
            DX_WINDROT = 1. ! Default
            DO IPAR=3, NPAR
              CALL SARGV (KARTE,IPAR,ACTPAR2)
              IF (ACTPAR2     .EQ. 'RSTAR') BVSINI_AT_RCOROT = .FALSE. 
              IF (ACTPAR2(:4) .EQ. 'RCOR' ) BVSINI_AT_RCOROT = .TRUE. 
              IF (ACTPAR2(:2) .EQ. 'DX') THEN
                 CALL SARGV (KARTE,IPAR+1,ACTPAR2)
                 READ (ACTPAR2, '(F20.0)', ERR=98) DX_WINDROT
                 IF (DX_WINDROT .LE. .0) THEN
                    WRITE (0,*) '*** ERROR: DX must be positive!'
                    GOTO 99
                 ENDIF
              ENDIF
            ENDDO
              
C** RCOROTLINE contains RCOROT and unit description. See routine "rotation_prep" for description.            
      ELSE IF ( ACTPAR .EQ. 'RCOROT' ) THEN
C                            ======
         RCOROTLINE = KARTE
      ELSE IF (KARTE(:4) .EQ. 'DISP') THEN
C                              ====
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) DISP

      ELSE IF (ACTPAR .EQ. 'SECONDMODEL') THEN
C                           ===========
         SECONDMODEL_LINE(1+IDX(SECONDMODEL_LINE):)
     >         = ' ' // KARTE(12:IDX(KARTE))

      ELSE IF ( ACTPAR .EQ. 'LPHISTA' ) THEN
C                            =======      
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=98) LPHISTA_ORIG

      ELSE IF ( ACTPAR .EQ. 'LPHIEND' ) THEN
C                            =======      
      
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=98) LPHIEND_ORIG
         
      ELSE IF ( ACTPAR .EQ. 'JPFIRST' ) THEN
C                               ======
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I3)', ERR=98) JPFIRST_ORIG

      ELSE IF ( ACTPAR .EQ. 'JPLAST' ) THEN
C                               ======
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I3)', ERR=98) JPLAST_ORIG

      ELSE IF ( ACTPAR == "NOWIND" ) THEN
C                          ======      
         NOWIND_LINE = KARTE 

      ELSE IF ( KARTE(:9) .EQ. 'PRINT PRO') THEN
C                               ========
            LSPRO=1

      ELSE IF ( KARTE(:10) .EQ. 'PRINT OPAL' ) THEN
C                                ==========
            CALL SARGV(KARTE, 3, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=98) XL
            LSOPA=1
            IF (XL .GT. 1.) LSOPA=IFIX(XL)

      ELSE IF ( ACTPAR .EQ. 'VDOP' ) THEN
C                               ====
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) VDOP 
         DD_VDOP_LINE = KARTE
      ELSE IF ( ACTPAR .EQ. 'VMIC' ) THEN
C                               ====
         BMICROTURB = .TRUE.
         DD_VDOP_LINE = KARTE

      ELSE IF ( KARTE(:8) .EQ. 'IVERSION' ) THEN
C                               ========
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'Z'  ) THEN 
               IVERSION=1
            ELSE IF (ACTPAR .EQ. 'TAU') THEN 
               IVERSION=0
            ELSE
               IVERSION=0 
               WRITE (0, '( A)') '************************************'
               WRITE (0, '(2A)') 'WARNING: INVALID IVERSION: ', ACTPAR
               WRITE (0, '( A)') 'TAU VERSION WILL BE APPLIED INSTEAD!'
               WRITE (0, '( A)') '************************************'
            ENDIF

      ELSE IF ( KARTE(:5) .EQ. 'REDIS' ) THEN
C                               =====
            REDIS = .TRUE.

      ELSE IF ( KARTE(:7) .EQ. 'NOREDIS' ) THEN
C                               =======
            REDIS = .FALSE.

      ELSE IF ( KARTE(:5) .EQ. 'CALIB' ) THEN
C                               =====
            BCALIBRATED = .TRUE.
            
            CALL SARGC (KARTE, NPAR)
            IF (NPAR .GT. 1) CALL SARGV (KARTE, 2, FUNIT)
                
      ELSE IF ( KARTE(:7) .EQ. 'NOCALIB' ) THEN
C                               =======
            BCALIBRATED = .FALSE. 

      ELSE IF ( KARTE(:4) .EQ. 'CONT' ) THEN
C                               =====
            BCONT = .TRUE.
            
            CALL SARGC (KARTE, NPAR)
            IF (NPAR .GT. 1) CALL SARGV (KARTE, 2, FUNIT)

      ELSE IF ( KARTE(:6) .EQ. 'NOCONT' ) THEN
C                               =====
            BCONT = .FALSE.

      ELSE IF ( KARTE(:5) .EQ. 'IDENT' ) THEN
C                               =====
            IDENT = .TRUE.
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'OSMIN') THEN
                CALL SARGV (KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F10.0)', ERR=98) OSMIN
            ENDIF

      ELSE IF ( KARTE(:7) .EQ. 'NOIDENT' ) THEN
C                               =======
            IDENT = .FALSE.

      ELSE IF ( KARTE(:6) .EQ. 'BWESEX' ) THEN
C                               ======
            CALL SARGV(KARTE, 2, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=98) BWESEX

      ELSE IF ( KARTE(:5) .EQ. 'BROAD' .OR. 
C                               =====
     >          KARTE(:8) .EQ. 'ALLBROAD') THEN
C                               ========
            BROAD = .TRUE.

      ELSE IF ( KARTE(:7) .EQ. 'NOBROAD' .OR.
C                               =======
     >         KARTE(:10) .EQ. 'NOALLBROAD' ) THEN
C                               ==========
            BROAD = .FALSE.

      ELSE IF ( KARTE(:8) .EQ. 'ABS WAVE' ) THEN
C                               ========
            ABSWAV = .TRUE.

      ELSE IF ( KARTE(:8) .EQ. 'REL WAVE' ) THEN
C                               ========
            ABSWAV = .FALSE.

      ELSE IF ( KARTE(:8) .EQ. 'LISTONLY' ) THEN
C                               ========
            LINELIST = .TRUE.

      ELSE IF ( KARTE(:6) .EQ. 'STRING' ) THEN
C                               ======
         CALL SARGC (KARTE, NPAR)
         IF (NPAR .LE. 2) THEN 
            WRITE (*,*) 'STRING option has not enough parameters'
            WRITE (0,*) 'STRING option has not enough parameters'
            STOP '*** ERROR detected by Subr. DECFORM'
         ENDIF
         CALL SARGV (KARTE, 2, ACTPAR)
         IF (ACTPAR .NE. 'COMMENT') THEN
            WRITE (*,*) 'STRING COMMENT is the only supported option'
            WRITE (0,*) 'STRING COMMENT is the only supported option'
            STOP '*** ERROR detected by Subr. DECFORM'
         ENDIF
         CALL SARGREST (KARTE,NPAR, 3, ISTART, IEND)         
         IEND = MIN0(IEND, 126)
         IF (NSTRING .EQ. 0) THEN
            WRITE (STRING1(0),'(A)') 'PLOT: ' // KARTE(ISTART:IEND)
            FREQIN=KARTE(ISTART:IEND)
         ENDIF
         NSTRING = NSTRING+1
         IF (NSTRING .LE. MAXSTRI) THEN
               WRITE (STRING1(NSTRING),'(A)') KARTE(ISTART:IEND)
         ELSE
            WRITE (*,*) 'WARNING: too many string options - ignored:'
            WRITE (*,*) KARTE
            WRITE (0,*) 'WARNING: too many string options - ignored:'
            WRITE (0,*) KARTE
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'TAUMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) TAUMAX 

      ELSE IF ( ACTPAR .EQ. 'XMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) XMAX 

      ELSE IF ( ACTPAR .EQ. 'DXMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) DXMAX 

      ELSE IF ( ACTPAR == 'TAUBROAD' ) THEN
C                          ========
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) TAUMINBROAD

      ELSE IF ( ACTPAR .EQ. 'NO-IRONLINES' ) THEN
C                            ===========
              BIRONLINES =.FALSE.

      ELSE IF ( ACTPAR .EQ. 'IRONLINES' ) THEN
C                            =========
              BIRONLINES = .TRUE.

      ELSE IF ( ACTPAR .EQ. 'WAVELENGTH' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, ACTPAR)
         IF ( ACTPAR .EQ. 'AIR') THEN
              BAIRWAVELENGTHSET=.TRUE.
              BAIRWAVELENGTH=.TRUE.
         ELSE IF ( ACTPAR .EQ. 'VACUUM' ) THEN
              BAIRWAVELENGTHSET=.TRUE.
              BAIRWAVELENGTH=.FALSE.
         ELSE
              WRITE (0,*) '*** INVALID OPTION'
              WRITE (0,*) '*** Allowed are only AIR or VACUUM'
              GOTO 99
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'NO-MODELCONT' ) THEN
C                            ===========
              BNOCONT=.TRUE.

      ELSE IF ( ACTPAR .EQ. 'SET_POP_ZERO' ) THEN
C                            ============
         CALL SARGV (KARTE, 2, ACTPAR)
         NSPZ = NSPZ + 1
         IF (NSPZ .GT. MAXSPZ) THEN
            WRITE (0,'(2A, I4)') 
     >       '*** WARNING: more SET_POP_ZERO commands than ',
     >       'dimensioned - MAXSPZ =', MAXSPZ
         ELSE
            SPZ1(NSPZ) = ACTPAR
            IF (NPAR .GE. 4) THEN
              CALL SARGV (KARTE, 3, ACTPAR)
              IF (ACTPAR .EQ. 'EXCEPT') THEN
                CALL SARGV (KARTE, 4, ACTPAR)
                SPZ2(NSPZ) = ACTPAR
              ELSE
                WRITE (0,'(A,A10)') 'Keyword wrong: ACTPAR = ', ACTPAR
                WRITE (0,'(2A)') 'KARTE = ', KARTE( :IDX(KARTE))
                STOP 'ERROR in Subr. DECFORM'
              ENDIF
            ENDIF
         ENDIF
      ELSE IF ( ACTPAR .EQ. 'MANIPOP_OPTIONS' ) THEN
C                            ===============
         CALL SARGP (KARTE, NPAR, 2, ISTART, IEND)
         IF (ISTART .GT. 0) MANIPOP_OPTIONS = KARTE(ISTART:)

      ELSE IF ( ACTPAR .EQ. 'XUNIT' ) THEN
C                            =====
         CALL SARGV (KARTE, 2, XUNIT)

      ELSE IF ( ACTPAR .EQ. 'MACROCLUMP') THEN
C                            ==========
         MACROCLUMPLINE = KARTE

      ELSE IF ( ACTPAR .EQ. 'PATH_VCSSB' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, PATH_VCSSB)

      ELSE IF ( ACTPAR .EQ. 'PATH_LEMKE_DAT' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, PATH_LEMKE_DAT)

      ELSE IF ( ACTPAR == 'NO-BIGBANDCUT' ) THEN
C                          =============
         bBIGBANDLIMIT = .FALSE.
         
      ELSE IF ( ACTPAR == 'BIGBANDCUT' ) THEN
C                          ==========
         bBIGBANDLIMIT = .TRUE.
         
      ELSE IF ( ACTPAR == 'NO-DDVDOPREDIS' ) THEN
C                          ==============
C***     Switches off depth-dependent opacity profile functions in FORMCMF
C***     (use for comparison calculations with models before Feb 2017)
         bDDOPAFORMCMF = .FALSE.
         
      ELSE IF ( ACTPAR == 'DDVDOPREDIS' ) THEN
C                          ===========
         bDDOPAFORMCMF = .TRUE.
         
      ELSE IF ( ACTPAR == 'NO-FECONVOL' ) THEN
C                          ===========
C***     Switches off depth-dependent iron opacity convolution
C***     (use for comparison calculations with models before Feb 2017)
         bDDFECONVOL = .FALSE.
         
      ELSE IF ( ACTPAR == 'FECONVOL' ) THEN
C                          ========
         bDDFECONVOL = .TRUE.
         
      ENDIF
 
      GOTO 1

C***  Error branches
   98 WRITE (0,*) '*** ERROR when decoding parameter as number'

   99 WRITE (0,*) '*** The error occured in the following line:'
      WRITE (0,*) KARTE(:IDX(KARTE))
      STOP ' *** FATAL ERROR DETECTED BY SUBR. DECFORM'

 
      END
 

      SUBROUTINE DIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR,DTDR,TEFF,NOTEMP)
C***********************************************************************
C***  CALLED FROM: WRCONT, ETL
C***  GIVES THE PLANCK FUNCTION, BCORE, AND ITS RADIUS-DERIVATIVE, DBDR, AT
C***  THE INNER BOUNDARY FROM ...
C***  IF (NOTEMP) ... THE GIVEN TEMPERATURE STRATIFICATION
C***  IF (.NOT. NOTEMP) ... "DTDR", THE DERIVATIVE OF THE TEMPERATURE WITH
C***                    RESPECT TO THE RADIUS, AS CALCULATED IN SUBR. DIFDTDR
C***                    TO YIELD THE CORRECT TOTAL FLUX  H = L/(4*PI*RSTAR)**2
C***  IN DIFFUSION APPROXIMATION, THEN THE INCIDENT INTENSITY WILL BE
C***      IPLUS = BCORE + DBDR * Z / X
C***  Z = MUE, X = OPACITY
C***********************************************************************
 
      DIMENSION T(ND),RADIUS(ND)
      LOGICAL NOTEMP
 
C***  Prevent Zeros in the radiation field XJC in case of 
c     small (X-ray) wavelength together with low temperatures
      BCORE=MAX(BNUE(XLAM,T(ND)),1.E-286)
 
      IF (NOTEMP) THEN
         BNDM=BNUE(XLAM,T(ND-1))
         DBDR=(BCORE-BNDM)/(RADIUS(ND-1)-1.)
      ELSE
C***     DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
         DBDT=DBNUEDT(XLAM,T(ND))
         DBDR=DBDT*DTDR
      ENDIF
 
      RETURN
      END
      SUBROUTINE DRTRANS (XLAM,LINE,LOW,NUP,INDLAP,XLAMLAP,
     >                    DELXLAP,NBLINE,MAXLAP,INDLOW,INDNUP,LASTIND,
     >                    MAXIND,LEVEL,WEIGHT,ELEVEL,N,EINST,NDIM,
     >                    POPNUM,T,ND,ALN,VDOP,EION,ENTOT,RNE,
     >                    MAXSUBL,NSUBLOW,BROAD,LINPRO,AVOIGT, 
     >                    DENSCON,NMOD, MAXMOD, NDDIM, 
     >                    MAINQN, NCHARG, NOM, IND_ORIGLEV)
C***********************************************************************
C***  SUBROUTINE FOR DRTRANSIT HANDLING OF MAIN PROGRAM "FORMAL"
C***  CALLED FROM SUBR. PREFORM IN CASE OF DECODED OPTION "DRTRANSIT"
C***  INPUT OPTIONS: "/AUTONIVEAU" - DATA FOR UPPER (AUTOIONIZATION) LEVEL
C***                 "/LOWERLEVEL" - DATA FOR LOWER SUBLEVEL
C***                 "/ADDLINE"    - DATA FOR HANDLING OF ADDITIONAL LINE
C***                 "-DRTRANS"    - END OF DRTRANSIT INPUT 
C***  ACTION: 1. READ INPUT DATA (LEVELS, LINES) FOR DRTRANSIT HANDLING
C***          2. CALCULATE POPNUMBER FROM SAHA-BOLTZMANN FORMULA
C***             (I.E. LTE POPNUMBER)
C***          3. SPLIT LOWER LEVEL INTO LTE-LEVELS (OPTIONAL)
C***********************************************************************

      DIMENSION ND(MAXMOD)
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD)
      DIMENSION WEIGHT(N),ELEVEL(N),EION(N), MAINQN(N), NOM(N)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM), NCHARG(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION POPNUM(NDDIM,NDIM,NMOD)
      REAL, DIMENSION (NDDIM,NMOD) :: T, ENTOT, RNE, DENSCON
      DIMENSION NSUBLOW(MAXSUBL)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(N), LEV, LEVUP,LEVLOW
      CHARACTER*8 LINPRO(MAXLAP), LINPROBL
      LOGICAL BROAD

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

C***  DEFAULTS FOR DRTRANSIT HANDLING:
      NBLSAVE=NBLINE
C***  NLOW: NUMBER OF /LOWERLEVELS from splitting LEVELS
      NLOW=0
C***  NSUBLOW: POINTER TO ADDITIONAL LEVELS IN THE ORIGINAL ARRAYS
      DO 10 I=1,MAXSUBL
   10 NSUBLOW(I)=0

C***  1. READ INPUT FOR DRTRANSIT HANDLING OF SPECIFIED LINE  ----------
    1 READ (2,'(A)') KARTE
      IF (KARTE(1:1) .EQ. '*') GOTO 1
      IF (KARTE      .EQ. '' ) GOTO 1

      IF (KARTE(:8) .EQ. '-DRTRANS') GOTO 20
C                         ========

C************************************************************************
      IF (KARTE(:11) .EQ. '/AUTONIVEAU') THEN
C                          ===========
         READ (KARTE,3) LEV, NW, ELEV
    3    FORMAT (12X,A10,1X,I4,1X,F20.0)

C***     FIND INDEX of AUTOLEVEL (IF ALREADY EXISTING):
         JNUP=0
C***     a) ALREADY EXISTING ?
         DO J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV) .OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
  202             FORMAT (' >>>>> SUBR. DRTRANS: WARNING: ', //, 
     >                    'INCONSISTENT ASSIGNMENT OF LEVEL ', A)
  203             FORMAT ('ELEVEL=',2(F9.1,1X),F7.4)
                  WRITE (0,202) LEV
                  WRITE (0,203) ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,202) LEV
                  WRITE (*,203) ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ENDIF
               JNUP=J
               EXIT
            ENDIF
         ENDDO

C***     b) LEVEL NOT YET DEFINED -> append to the list 
         IF (JNUP .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *, ' >>>>> SUBR. DRTRANS: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('DRTRANS: N GREATER THAN NDIM')
               STOP 'NDIM exceeded in DRTRANS'
            ENDIF

            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
C*          The AUTO level belongs to the lower ion - copy from there:
            MAINQN(N)= MAINQN(LOW)
            NCHARG(N)= NCHARG(LOW)
            NOM(N)   = NOM(LOW)
            EION(N)  = EION(LOW)
            IND_ORIGLEV(N) = LOW

C***        POPNUMBER OF THE NEW AUTOIONIZATION LEVELS: RELATIVE TO THE 
C***        GIVEN UPPER STATE by SAHA-BOLTZMANN (LTE ratio)
C***        Note: AUTONIVEAU population is *not* subtracted from 
C***              population of the parent level! --> inconsistent! ллллллллллллллл
            WAUTO= WEIGHT(N)
            WNUP = WEIGHT(NUP)
            EAUTO= ELEVEL(N)
            EILOW= EION(LOW)
            ENUP = ELEVEL(NUP)

            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  TL = T(L,IMOD)
                  SQRTL = SQRT(TL)
                  ENTOTL = ENTOT(L,IMOD) * DENSCON(L,IMOD)
                  RNEL = RNE(L,IMOD)
                  IPOP = L+ (NUP-1)*ND(IMOD)
                  POPNUP = POPNUM(IPOP,1,IMOD)
                  BOLTZ = EXP(-C1*(EAUTO-EILOW-ENUP)/TL)
                  SAHAPHI = WAUTO/WNUP*CI/TL/SQRTL * BOLTZ

ccc wrong versions of BOLTZ:  wrh, 26-Oct-2021
cc???     >           BOLTZ = EXP(-C1*(EAUTO-EILOW)/TL)
cc???     >           BOLTZ = EXP(-C1*(ELOW-EILOW-ENUP)/TL)

                  POPAUTO = POPNUP * ENTOTL * RNEL * SAHAPHI
                  IPOP = L+ (N-1)*ND(IMOD)
                  POPNUM(IPOP,1,IMOD) = POPAUTO

               ENDDO
            ENDDO 
         ENDIF

C************************************************************************
      ELSE IF (KARTE(:11) .EQ. '/LOWERLEVEL') THEN
C                               ===========
         N=N+1
         IF (N .GT. NDIM) THEN
            PRINT *, ' >>>>> SUBR. DRTRANS: ERROR STOP (N .GT. NDIM)'
            CALL REMARK ('DRTRANS: N GREATER THAN NDIM')
            STOP 'NDIM1'
         ENDIF
         NLOW=NLOW+1
         IF (NLOW .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. DRTRANS: ERROR STOP (NLOW .GT. MAXSUBL)'
            CALL REMARK ('DRTRANS: NLOW GREATER THAN MAXSUBL')
            STOP 'NLOW'
         ENDIF
         NSUBLOW(NLOW)=N
         READ (KARTE,3) LEVEL(N),NW,ELEVEL(N)
         WEIGHT(N)=FLOAT(NW)

C*       The sublevel stems from the lower level - copy from there:
         MAINQN(N)= MAINQN(LOW)
         NCHARG(N)= NCHARG(LOW)
         NOM(N)   = NOM(LOW)
         EION(N)  = EION(LOW)
         IND_ORIGLEV(N) = LOW

C************************************************************************
      ELSE IF (KARTE(:8) .EQ. '/ADDLINE') THEN
C                              ========

C***     First: if the LOWERLEVEL is to be split:
C***     BOLTZMANN (LTE)
         IF (NLOW .GT. 0) THEN
            DO IMOD=1, NMOD
               DO 30 L=1, ND(IMOD)
                  IPOP = L+ (NSUBLOW(1)-1)*ND(IMOD)
                  POPNUM(IPOP,1,IMOD) = 1.
                  DO J=2,NLOW
                    IPOP   = L+ (NSUBLOW(J  )-1)*ND(IMOD)
                    IPOPM1 = L+ (NSUBLOW(J-1)-1)*ND(IMOD)
                    POPNUM(IPOP,1,IMOD) = 
     >               EXP(C1*(ELEVEL(NSUBLOW(J-1))
     -               -ELEVEL(NSUBLOW(J)))/T(L,IMOD))*WEIGHT(NSUBLOW(J))
     /               /WEIGHT(NSUBLOW(J-1))*POPNUM(IPOPM1,1,IMOD)
                  ENDDO

C***              NORMALIZATION
                  SUM=0.
                  DO J=1,NLOW     
                     IPOP = L+ (NSUBLOW(J)-1)*ND(IMOD)
                     SUM = SUM + POPNUM(IPOP,1,IMOD)
                  ENDDO

                  IPOP   = L+ (LOW-1)*ND(IMOD)
                  SUM = SUM / POPNUM(IPOP,1,IMOD)

                  DO J=1,NLOW
                     IPOP = L+ (NSUBLOW(J)-1)*ND(IMOD)
                     POPNUM(IPOP,1,IMOD) 
     >                 = POPNUM(IPOP,1,IMOD) / SUM
                  ENDDO

   30          CONTINUE
            ENDDO
         ENDIF

C***     Now executing the ADDLINE actions
         LASTIND=LASTIND+1

         IF (LASTIND .GT. MAXIND) THEN
            PRINT *,
     >        ' >>>>> SUBR. DRTRANS: ERROR STOP (LASTIND .GT. MAXIND)'
            CALL REMARK ('DRTRANS: LASTIND GREATER THAN MAXIND')
            STOP 'MAXIND'
         ENDIF
         IF (NBLINE+1 .GT. MAXLAP) THEN
            PRINT *,
     >          ' >>>>> SUBR. DRTRANS: ERROR STOP (NBLINE .GT. MAXLAP)'
            CALL REMARK ('DRTRANS: NBLINE GREATER THAN MAXLAP')
            STOP 'MAXLAP'
         ENDIF

         READ (KARTE,23) LEVUP,LEVLOW,AUPLOW
   23    FORMAT (9X,A10,2X,A10,G10.0)

C***     FIND UPPER INDEX:
         DO 24 J=1,N
         JNUP=J
         IF (LEVEL(J) .EQ. LEVUP) GOTO 25
   24    CONTINUE

   90    FORMAT ('DRTRANS: UPPER LINE LEVEL NOT FOUND: ', A10)
         WRITE (0, 90)  LEVUP
         STOP 'UPPER'

C***     FIND LOWER INDEX:
   25    DO 26 J=1,N
         JLOW=J
         IF (LEVEL(J) .EQ. LEVLOW) GOTO 27
   26    CONTINUE

   91    FORMAT ('DRTRANS: LOWER LINE LEVEL NOT FOUND: ', A10)
         WRITE (0, 91) LEVLOW
         STOP 'LOWER'

   27    INDNUP(LASTIND)=JNUP
         INDLOW(LASTIND)=JLOW

C***     Check for correct sequence UP - LOW
         IF (ELEVEL(JNUP) .LE. ELEVEL(JLOW)) THEN
            WRITE (0,'(A)') '*** FATAL ERROR in FORMAL_CARDS:'
            WRITE (0,'(A)') '*** E(UP) must be higher than E(LOW) !'
            WRITE (0,'(A)') '*** error occured in the following line:'
            WRITE (0,'(A)') KARTE(:IDX(KARTE))
            STOP '*** ERROR in subroutine DRTRANS'
         ENDIF

C***     default wavelength 
         XLAMSUB=1.E8/(ELEVEL(JNUP)-ELEVEL(JLOW))

C***     read wavelength if given, covert from AIR to VAC
C***     read optional VOIGT parameter only if BROAD=.TRUE.
          CALL READ_LINECARD_PARAMETERS (KARTE(43:), XLAMSUB,
     >            BROAD, LINPROBL, AVOIGTBL)
          IF (BROAD .AND. (LINPROBL .EQ. '')) LINPROBL = 'DRTRANS '

C***     If first line in blend block: it defines reference lambda
         IF (NBLINE .EQ. 0) XLAM = XLAMSUB

C***     In case AUPLOW < 0 it means f-Value --> convert
         IF (AUPLOW .LT. 0.) THEN
            EINST(JNUP,JLOW)=-6.669E15/XLAMSUB/XLAMSUB*WEIGHT(JLOW)
     /                                      /WEIGHT(JNUP)*AUPLOW
         ELSE
            EINST(JNUP,JLOW)=AUPLOW
         ENDIF

C***     Insert current ADDLINE in the list sorted by increasing DELTAX
         CALL INSERT_LINE (LINE, LASTIND, NBLINE, INDLAP,
     >                   XLAMLAP,DELXLAP,
     $                   XLAMSUB, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C************************************************************************
      ELSE
         WRITE (0,*) 'UNRECOGNIZED INPUT CARD IN SUBR. DRTRANS!'
         WRITE (0,*) 'KARTE=', KARTE( :IDX(KARTE))

      ENDIF
C************************************************************************

      GOTO 1


   20 CONTINUE
C***  END OF INPUT FOR DRTRANSIT HANDLING  -----------------------------


C***  NO SUBLINES DECODED: RETURN WITHOUT ANY DR-TRANSITION
      IF (NBLINE .EQ. NBLSAVE) THEN
        WRITE (0,*) '*** WARNING: DRTRANSIT BLOCK WITHOUT LINES!'
      ENDIF

   99 CONTINUE


      RETURN
      END
      SUBROUTINE ELIMIN (XLAM,EMFLUX,FLUXIN,U,Z,
     $          A,B,C,W,BX,WX,XJC,RADIUS,P,BCORE,DBDR,
     $          OPA,ETA,THOMSON,EDDI,ND,NP,NPDIM,ENTOT,K,
     $          IWARN, IWARN2, ST, BELIFI, IVERS)
C***  FEAUTRIER SCHEME FOR CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY
C***  TAPE7 = MASS STORAGE FILE FOR FEAUTRIER MATRICES
C***  LAST PARAMETER -1 IN WRITMS IS VERY IMPORTANT ]
C***  OTHERWISE "MASS STORAGE LIMIT" EXCEEDED BECAUSE OLD RECORDS ARE
C***  NOT OVERWRITTEN .
C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY ]
      DIMENSION XJC(ND),RADIUS(ND),OPA(ND),ETA(ND),THOMSON(ND)
      DIMENSION EDDI(3,ND)
      DIMENSION P(NPDIM),BX(NPDIM,NPDIM),WX(NPDIM),B(2)
      DIMENSION U(ND,NP),Z(ND,NP)
      
C***  To avoid underflows define minimum allowed positive real/intensity/flux
      REAL, PARAMETER :: RTINY = 1.E-300

C***  To Store Feautrier Matrices ST(94*95,89)
      DIMENSION ST((NPDIM+1)*NPDIM,ND)
      LOGICAL BELIFI
      CHARACTER*4 CKEY

C***  Operating system:                    
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C***  GAUSS-ELIMINATION
      DO 1 L=1,ND
      CALL       SETUP (L,A,B,C,W,JMAX,ND,NP,NPDIM,OPA,ETA,THOMSON,
     $          XLAM,Z,RADIUS,BCORE,DBDR,XIMINUS,ENTOT,K, IVERS)
      IF (L.EQ.1) GOTO 6
      CALL MDMV (A,BX,JMAX,NPDIM)
      CALL MSUB (B,BX,JMAX,NPDIM)
      CALL MDV (A,WX,JMAX)
      CALL VADD (W,WX,JMAX)
    6 CONTINUE
      IF (XLAM .LE. 4. .AND. OPSYS(1:4) .EQ. 'CRAY') THEN
        CKEY = 'OWN'
      ELSE
        CKEY = 'FREE'
      ENDIF
      CALL INV (JMAX,NPDIM,B,CKEY)
      CALL MVV (WX,B,W,JMAX,JMAX,NPDIM)
      IF (L.EQ.ND) GOTO 2
      CALL MVMD (BX,B,C,JMAX,JMAX-1,NPDIM)
C***  COMPRESSING THE MATRIX BX  AND VECTOR WX  INTO THE RANGE OF B AND C
      DO 7 J=1,JMAX
      JC=1+(J-1)*JMAX
    7 CALL EQUAL (B (JC)  ,BX(1,J),JMAX)
      CALL EQUAL (B (JMAX*JMAX+1)  ,WX,JMAX)
      LL=JMAX*(JMAX+1)
      IF (BELIFI) THEN
        CALL WRITMS (7,B ,LL,L,-1, IDUMMY, IERR)
      ELSE
        DO I=1, LL
          ST(I,L) = B(I)
        ENDDO
      ENDIF
    1 CONTINUE
 
C***  BACK SUBSTITUTION
C***  RECENT WX IS THE FEAUTRIER-INTENSITY U AT THE INNER BOUNDARY
    2 CALL MOMENT0 (ND,RADIUS,ND,JMAX,Z,WX,XJC(ND),.FALSE.)
      CALL MOMENT1 (RADIUS(ND),JMAX,P,WX,H)
      HPLUS=BCORE/2. + DBDR/3./OPA(ND)
      CALL MOMENT2 (RADIUS(ND),JMAX,P,WX,XK)
C***  EDDI(1,L) IS THE EDDINGTON FACTOR  F = K / J
C***  EDDI(2,L) IS THE SPHERICITY FACTOR Q
C***  EDDI(3,L) IS THE EDDINGTON FACTOR H / J  (ONLY AT THE BOUNDARIES)
C***  EDDI(3,ND-1) IS THE OUTWARD FLUX HPLUS AT THE INNER BOUNDARY
      IF (XJC(ND) .GT. .0) THEN
         EDDI(1,ND)=XK/XJC(ND)
      ELSE
         EDDI(1,ND)=1./3.
         IWARN = IWARN + 1
      ENDIF

      EDDI(2,ND)=1.
      RRQ=1.
      IF (XJC(ND) .GT. .0) THEN
         EDDI(3,ND)=H/XJC(ND)
      ELSE
         EDDI(3,ND)=0.5
         IWARN = IWARN + 1
      ENDIF
 
      EDDI(3,ND-1)=HPLUS
      FLUXIN=4*(HPLUS-H)
      IF (EDDI(1,ND) .LT. 0.01) THEN
        IWARN = IWARN + 1
        EDDI(1,ND) = 0.01
      ENDIF
      IF (EDDI(1,ND) .GT. 1.) THEN
        IWARN2 = IWARN2 + 1
        EDDI(1,ND) = 1.0
      ENDIF
      FL=3.-1./EDDI(1,ND)
      DO 5 J=1,JMAX
        IF (ABS(WX(J)) < RTINY) WX(J) = RTINY
    5 U(ND,J)=WX(J)
      CALL EQUAL (A,WX,JMAX)
      NC2=NP-ND+2
 
C***  L = ND-1 ... 1
      DO 4 JMAX=NC2,NP
      L=NP+1-JMAX
      RL=RADIUS(L)
      LL=JMAX*(JMAX+1)
      IF (BELIFI) THEN
        CALL READMS (7,B ,LL,L, IERR)
      ELSE
        DO I=1, LL
          B(I) = ST(I,L)
        ENDDO
      ENDIF
C***  DECOMPRESSING THE MATRIX BX AND VECTOR WX  OUT OF B (AND C)
      DO 8  J=1,JMAX
      JC=1+(J-1)*JMAX
    8 CALL EQUAL (BX(1,J),B(JC)  ,JMAX)
      CALL EQUAL (WX,B (JMAX*JMAX+1)  ,JMAX)
      CALL MVV (W,BX,A,JMAX,JMAX-1,NPDIM)
      CALL VADD (WX,W,JMAX)
C***  WX(J) IS THE FEAUTRIER-INTENSITY U AT RADIUS R(L)
      DO 3 J=1,JMAX
        IF (ABS(WX(J)) < RTINY) WX(J) = RTINY
    3 U(L,J)=WX(J)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,WX,XJC(L),.FALSE.)
c      WRITE (0,*) ' L=', L
c      WRITE (0,*) ' RL=', RL
c      WRITE (0,*) ' JMAX=', JMAX
c      WRITE (0,*) ' P=', P
c      WRITE (0,*) ' WX=', WX
c      WRITE (0,*) ' XK=', XK
      IF (ABS(XK) < RTINY) XK = RTINY
      CALL MOMENT2 (RL,JMAX,P,WX,XK)
      IF (XJC(L) .GT. .0) THEN
         EDDI(1,L)=XK/XJC(L)
      ELSE
         EDDI(1,L)=1./3.
         IWARN = IWARN + 1
      ENDIF

      IF (EDDI(1,L) .LT. 0.01) THEN
        IWARN = IWARN + 1
        EDDI(1,L) = 0.01
      ENDIF
      IF (EDDI(1,L) .GT. 1.) THEN
        IWARN2 = IWARN2 + 1
        EDDI(1,L) = 1.
      ENDIF
C***  THIS IS AN INGENIOUS (;) RECURSION FORMULA FOR THE SPHERICITY FACTOR ]
      RLP=RADIUS(L+1)
      FLP=FL
      FL=3.-1./EDDI(1,L)
      RRQ=RRQ *EXP(FL-FLP)*(RL/RLP)**((FLP*RL-FL*RLP)/(RL-RLP))
      EDDI(2,L)=RRQ/RL/RL
    4 CALL EQUAL (A,WX,JMAX)

      CALL MOMENT1 (RADIUS(1), NP, P, WX, H)

C***  CORRECT THE EDDINGTON FLUX AT THE OUTER BOUDARY
C***    FOR A (POSSIBLE) INCIDENT RADIATION:
      H = H - 0.5 * XIMINUS

C***  EMFLUX IS THE EMERGENT FLUX, BUT RELATED TO THE INNER RADIUS
      EMFLUX=4.*H*RADIUS(1)*RADIUS(1)

C***  BOUNDARY EDDINGTON FACTOR
      IF (XJC(1) .GT. .0) THEN
         EDDI(3,1)=H/XJC(1)
      ELSE
         EDDI(3,1)=0.5
         IWARN = IWARN + 1
      ENDIF
 
      RETURN
      END
      SUBROUTINE EQUAL (A,B,N)
C***********************************************************************
C***  VEKTOREN   A := B
C***********************************************************************

      DIMENSION A(N),B(N)
      DO 1 I=1,N
    1 A(I)=B(I)

      RETURN
      END
      SUBROUTINE FECHECK (XLAMK, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEINDR)
C***********************************************************************
C***  CHECK FOR ACTIVE BOUND-BOUND TRANSITIONS OF GENERIC ION
C***  CALLED FROM COLI
C***********************************************************************

      INTEGER, DIMENSION(LASTFE) :: INDRB, IFRBSTA, IFRBEND, INDFEACT
      INTEGER :: LASTFE, MAXFEACT
      LOGICAL :: BFECHECK, BFEWING

C***  SEARCH FOR LOWEST POSSIBLE TRANSITION

      IND = INDFEACT(1)

C***  CALCULATION OF ACTUAL FREQUENCY INDEX
      XINDF = - ALOG10(XLAM0FE/XLAMK) /
     >          ALOG10(1. + VDOPFE*DXFE/CLIGHT)

 10   IF (IND .GT. LASTFE) THEN
C ***   NO TRANSITIONS ARE ACTIVE ANYMORE
         BFECHECK = .FALSE.
         GOTO 20
      ENDIF

      XINDLOW  = FLOAT(IFRBSTA(IND))

C**** Cross-check: Is the first transition INDFEACT(1) an active transition?
C****              If not, set INDFEACT(1) to first active transition!
C****              If no transition is active, set BFECHECK = .FALSE.
      IF (XINDF .LT. XINDLOW) THEN
C ***   NO ACTIVE TRANSITION
         BFECHECK = .FALSE.
         GOTO 20
      ELSE
C ***   TRASITION ACTIVE?
         XINDUP = FLOAT(IFRBEND(IND))
         IF (XINDF .GT. XINDUP) THEN
C ***      TEST NEXT TRANSITION
            IND = IND + 1
            INDFEACT(1) = IND
            GOTO 10
         ELSE
C ***      TRANSITION ACTIVE!
            BFECHECK = .TRUE.
            GOTO 20
         ENDIF
      ENDIF

 20   CONTINUE

C  *** Improved Wing test: all former transitions are checked ***
C *** TEST IF RED WING OF ANY LINE IS ACTIVE
      IF ((.NOT. BFECHECK) .AND. (IND .GT. 1)) THEN
C ***   RED WING OF LAST LINE, WHICH HAS BEEN ACTIVE
         BFEWING = .FALSE.
         DO INDW=1, IND-1
           XINDUP = FLOAT(IFRBEND(INDW-1))
           XINDWING = XINDUP + DFEINDR
           IF (XINDF < XINDWING) THEN
             BFEWING = .TRUE.
             EXIT
           ENDIF
         ENDDO
      ELSE
C ***   WING IS ACTIVE, WHEN LINE IS ACTIVE!         
         BFEWING = BFECHECK
      ENDIF
      
C *** STORE INDICES OF ACTIVE LINES IN ARRAY INDFEACT
      MAXFEACT = 0
      DO IND=INDFEACT(1), LASTFE
         XINDLOW = FLOAT(IFRBSTA(IND))
         IF (XINDF .GE. XINDLOW) THEN
            XINDUP = FLOAT(IFRBEND(IND))
            IF (XINDUP .GE. XINDF) THEN
C ***         TRANSITION ACTIVE !
               MAXFEACT = MAXFEACT + 1
               INDFEACT(MAXFEACT) = IND
            ENDIF
C         ELSE           !these two lines
C            EXIT        !are taken from Goetz (can be added for speed, since INDEFEACT is sorted by wavelength)
         ENDIF
      ENDDO

      
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

      FUNCTION FIERFC(X)
C**********************************************************************
C***  FIRST REPEATED INTEGRAL OF THE ERROR FUNCTION (COMPLEMENT)
C***  - RATIONAL APPROXIMATION FROM ABRAMOWITZ AND STEGUN, P. 299
C***  - NOT VALID FOR 
C**********************************************************************
      DIMENSION A(3)
      DATA A / 0.3480242, -0.0958798, 0.7478556 /
      DATA WPIINV / 0.564189583549 /

      IF ( X .LT. 0. ) THEN
        WRITE (0,*) 'FUNCTION FIERFC CALLED WITH NEGATIVE ARGUMENT'
        STOP 'ERROR'
      ENDIF
      T = 1. / (1. + 0.47047 * X)
      POLY = ((A(3) * T + A(2)) * T + A(1)) * T
      FIERFC = (WPIINV - X*POLY) * EXP(-X*X)

      RETURN
      END
      SUBROUTINE FILTERFUNCTIONS (FILTERNAME, NFILT, FLAM, FILT)
C***********************************************************************
C***  Photometric filter functions (JOHNSON, STROEMGREN, etc.)
C***  REFERENCES: Heber et al. (1984: A&A 130, 119)
C***  https://ui.adsabs.harvard.edu/abs/1984A%26A...130..119H/
C***  In this form, newly written by wrh 25-May-2021
C***    but based on previous code that was part of subr. PRIJOST
C***  New filters added: 2MASS J, H, Ks
C***  18.05.21 by Helge: Gaia Gbp, Grp, G ; 
C***                     corrected Ks2MASS #points
C***
C*** Input : FILTERNAME
C*** Output: Vectors of length NFILT (MAXNFILT = 200)
C***              FLAM (wavelength points in Ang)
C***              FILT (filter function at wavelength points) 
C***          The filter functions are *not* yet normalized
C***
C*** Caveat: The wavelength grid is defined here as 
C***             FLAM(I) = FBLUE + I*FINCR
C***             i.e. not with (I-1)*FINCR as one would expect.
C***             This seems to be intentional and correct (wrh 28-May-2021).
C***             Therefore, FBLUE is not the first lambda point FLAM(1),  
C***             but FBLUE = FLAM(1) - FINCR
C***
C*** Called from: WRCONT - PRICOLR - PRIJOST
C***              STEAL  - PRICOLR - PRIJOST
C***              FORMAL - LIMBDARK_OUTPUT 
C***********************************************************************

      IMPLICIT NONE

      REAL UJ(14),BJ(22),VJ(25),US(28),VS(29),BS(29),YS(29)
      REAL J2MASS(194), H2MASS(90), Ks2MASS(90)
      REAL Gaia3Gbp(177), Gaia3Grp(152), Gaia3G(145)

      INTEGER NFILT, I

C***  Filter function: blue start, wavelength increment
      REAL FBLUE, FINCR, FILTMAX

      INTEGER, PARAMETER :: MAXNFILT = 200
      REAL, DIMENSION(MAXNFILT) :: FLAM, FILT
      CHARACTER FILTERNAME*(*)

C***  NAMES OF THE CONSIDERED COLOR BANDS
C***  1 - 3: JOHNSON U,B,V; 
C***  4 - 7: STROEMGREN u,v,b,y; 
C***  8-10: 2MASS J, H, K
C***  11-13: Gaia Gbp, Grp, G


C*** JOHNSON FILTER FUNCTIONS

      DATA UJ/0.000,0.025,0.250,0.680,1.137,1.650,2.006,2.250,2.337,
     -        1.925,0.650,0.197,0.070,0.000/


      DATA BJ/0.000,0.006,0.080,0.337,1.425,2.253,2.806,2.950,3.000,
     -        2.937,2.780,2.520,2.230,1.881,1.550,1.275,0.975,0.695,
     -        0.430,0.210,0.055,0.000/

      DATA VJ/0.020,0.175,0.900,1.880,2.512,2.850,2.820,2.625,2.370,
     -        2.050,1.720,1.413,1.068,0.795,0.567,0.387,0.250,0.160,
     -        0.110,0.081,0.061,0.045,0.028,0.017,0.007/

C***  STROEMGREN FILTER FUNCTIONS

      DATA US/00.00,00.60,05.92,11.80,17.80,23.80,29.12,33.32,36.08,
     -        38.12,39.28,39.44,39.08,38.52,37.08,35.32,33.08,30.08,
     -        26.16,21.68,16.92,11.76,07.68,04.68,02.16,00.88,00.40,
     -        00.00/
      DATA VS/00.00,00.16,00.36,00.84,01.48,02.28,03.04,04.80,07.84,
     -        13.00,20.00,29.88,40.00,47.28,49.32,48.00,43.52,37.20,
     -        28.20,18.12,11.12,06.72,04.00,02.68,02.00,01.40,00.72,
     -        00.36,00.00/
      DATA BS/00.00,00.40,00.96,01.60,02.32,03.60,05.00,08.00,12.40,
     -        20.00,30.20,40.40,45.88,46.88,45.04,38.20,27.84,17.40,
     -        11.00,06.88,03.92,02.52,01.52,01.24,00.92,00.72,00.40,
     -        00.20,00.00/
      DATA YS/00.00,00.68,01.68,02.72,04.00,06.96,10.20,15.20,23.40,
     -        33.00,41.00,45.40,47.92,50.12,52.20,51.32,45.44,34.80,
     -        24.80,15.48,10.00,06.60,04.16,02.52,01.36,00.96,00.80,
     -        00.40,00.00/

C***  2MASS filter functions
      DATA J2MASS 
     >  / 0.01, 0.03, 0.07, 0.12, 0.17, 0.22, 0.29, 0.46, 0.72, 1.06,
     >    1.63, 2.37, 3.10, 4.00, 4.60, 5.12, 5.48, 5.90, 6.41, 6.92,
     >   13.91,26.52,33.48,35.80,37.29,36.78,34.31,26.32,24.66,26.12,
     >   27.93,27.33,26.32,27.71,44.26,38.07,30.16,45.85,32.36,24.99,
     >   47.76,39.56,28.07,27.86,29.52,31.19,29.10,23.27,20.32,34.74,
     >   52.05,61.06,62.19,65.68,69.35,72.59,74.94,78.85,82.18,80.37,
     >   79.51,80.48,82.33,83.67,83.62,80.73,75.45,72.90,70.99,70.19,
     >   69.97,70.14,70.32,70.45,70.36,70.23,70.10,70.69,71.98,73.01,
     >   72.38,71.47,71.54,77.41,83.78,87.42,90.60,92.96,94.49,95.63,
     >   96.38,96.49,95.86,94.72,93.50,92.05,89.34,86.63,82.76,79.19,
     >   77.10,74.91,72.13,69.25,67.30,66.27,65.28,64.56,63.98,63.98,
     >   64.09,64.19,64.34,64.55,64.76,65.54,66.89,68.42,71.06,73.88,
     >   76.06,77.85,80.44,78.86,75.83,74.50,87.53,92.00,86.24,84.83,
     >   85.55,94.71,97.37,92.10,88.14,77.90,40.89,48.88,84.07,79.79,
     >   69.15,59.86,44.99,25.47,12.88, 2.69, 0.40, 1.62, 3.21, 2.82,
     >   0.73, 0.33, 0.21, 0.09, 1.26, 3.18, 1.89, 0.28, 0.03, 0.02,
     >   0.00, 0.00, 0.01, 0.09, 0.25, 0.23, 0.26, 4.44, 5.93, 1.31,
     >   3.13, 2.91, 1.24, 0.17, 0.39, 0.61, 0.78, 0.70, 0.39, 0.11,
     >   0.05, 0.04, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04,
     >   0.04, 0.03, 0.02, 0.01 /

      DATA H2MASS 
     >  / 0.00, 0.02, 0.03, 0.04, 0.06, 0.10, 0.15, 0.21, 0.26,0.38,
     >    0.55, 0.71, 1.33, 2.74, 6.21,10.79,15.98,22.54,33.06,43.67,
     >   53.37,62.91,71.02,77.02,82.66,87.22,90.52,92.33,92.80,92.41,
     >   88.90,86.32,86.96,89.05,91.02,92.31,91.54,91.26,92.10,92.41,
     >   92.12,92.30,92.39,92.76,93.77,94.55,95.62,97.38,99.29,99.69,
     >   99.86,98.65,96.96,95.46,94.48,93.54,92.73,94.07,96.14,98.00,
     >   98.61,99.00,99.09,98.97,98.71,98.04,96.74,92.44,88.07,82.49,
     >   74.67,65.22,46.83,31.70,39.98,22.79,13.57, 9.72, 6.03, 1.17,
     >    1.60, 0.63, 0.03, 0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.00 /

      DATA Ks2MASS 
     > /  0.02, 0.05, 0.24, 0.52, 0.87, 1.30, 2.20, 3.77, 6.48,10.76,
     >   17.06,22.33,20.39,23.14,26.12,29.60,36.93,35.95,44.71,59.87,
     >   69.29,75.65,75.22,72.74,68.95,65.16,64.05,68.44,71.10,72.46,
     >   78.62,81.12,81.97,82.26,84.31,86.60,87.46,86.31,87.93,90.12,
     >   90.96,91.78,92.40,92.68,92.67,92.21,91.29,90.42,91.19,90.83,
     >   86.50,86.43,90.72,94.53,96.32,97.82,98.63,98.48,97.83,97.12,
     >   96.83,98.04,98.25,98.21,97.10,96.87,98.64,97.22,95.79,94.14,
     >   89.55,90.72,80.14,62.72,50.22,41.48,31.24,20.55,14.93,10.36,
     >    7.30, 5.06, 3.68, 2.60, 1.26, 0.87, 0.33, 0.22, 0.26, 0.16 /

C***  Gaia EDR3 filter functions

      DATA Gaia3Gbp
     > /   0.13, 0.80, 2.77, 6.30,10.74,15.39,19.92,23.63,25.19,24.14
     >   ,21.67,18.96,17.88,18.46,20.27,22.31,23.40,23.39,22.85,22.24
     >   ,21.52,20.68,19.69,18.58,17.61,16.77,16.26,16.43,17.30,18.91
     >   ,21.35,24.57,28.15,31.93,35.69,39.22,42.46,45.22,47.67,49.74
     >   ,51.51,53.01,54.33,55.38,56.27,57.05,57.69,58.23,58.73,59.09
     >   ,59.50,59.82,60.12,60.47,60.84,61.11,61.44,61.79,62.12,62.37
     >   ,62.52,62.69,62.74,62.79,62.78,62.71,62.71,62.70,62.69,62.78
     >   ,62.80,62.90,62.99,63.14,63.26,63.42,63.46,63.55,63.50,63.55
     >   ,63.57,63.54,63.46,63.38,63.39,63.31,63.18,63.15,63.11,63.07
     >   ,63.07,63.00,63.03,62.98,63.00,62.96,62.98,62.96,62.83,62.73
     >   ,62.58,62.39,62.28,62.22,62.10,62.01,62.01,62.01,61.97,62.04
     >   ,62.18,62.25,62.26,62.34,62.27,62.15,61.96,61.77,61.58,61.37
     >   ,61.18,60.90,60.65,60.53,60.56,60.81,61.21,61.79,62.54,63.22
     >   ,63.83,64.25,64.30,64.27,64.08,63.60,62.96,62.30,61.68,61.31
     >   ,61.17,61.28,61.57,62.05,62.75,63.58,64.66,65.69,66.32,66.73
     >   ,66.90,66.77,66.42,65.77,64.87,63.80,62.43,60.93,59.28,57.40
     >   ,55.50,53.36,51.13,48.55,45.34,41.17,35.66,28.85,21.34,14.19,
     >     8.35, 4.30, 1.94, 0.78, 0.30, 0.13, 0.09 /

      DATA Gaia3Grp
     > /   0.12, 0.57, 2.45, 9.31,25.12,46.68,64.20,72.62,72.58,68.54
     >   ,68.06,69.02,69.67,69.94,70.31,70.86,71.45,72.14,72.61,73.00
     >   ,73.15,73.20,73.05,72.96,73.03,73.23,73.37,73.56,73.76,73.71
     >   ,73.75,73.56,73.40,73.21,73.11,73.18,73.33,73.71,73.99,74.29
     >   ,74.35,74.37,74.19,73.99,73.89,73.88,73.92,74.02,73.95,73.90
     >   ,73.63,73.19,72.74,72.35,71.84,71.38,70.99,70.64,70.46,70.40
     >   ,70.31,70.38,70.40,70.37,70.24,70.06,69.96,69.52,69.09,68.69
     >   ,68.15,67.64,67.17,66.70,65.92,65.27,64.66,63.76,62.90,62.03
     >   ,61.18,60.15,59.22,58.24,57.09,56.23,55.42,54.40,53.27,52.33
     >   ,51.32,50.41,49.30,48.22,47.14,45.90,44.73,43.55,42.37,41.11
     >   ,39.94,38.70,37.44,36.19,34.92,33.70,32.45,31.15,29.90,28.61
     >   ,27.34,26.10,24.85,23.67,22.43,21.25,20.12,19.02,17.96,16.90
     >   ,15.91,14.95,14.00,13.10,12.15,11.16,10.16, 9.14, 8.16, 7.22,
     >     6.32, 5.48, 4.70, 4.01, 3.40, 2.87, 2.41, 2.01, 1.67, 1.38,
     >     1.14, 0.93, 0.76, 0.61, 0.49, 0.39, 0.31, 0.24, 0.19, 0.15,
     >     0.11, 0.10 /

      DATA Gaia3G
     > /   0.39, 3.13, 7.22,10.42,12.21,12.72,12.40,11.71,11.11,10.85
     >   ,11.20,13.21,17.93,24.77,31.61,37.27,41.78,45.39,48.29,50.62
     >   ,52.56,54.17,55.53,56.73,57.80,58.73,59.57,60.32,61.02,61.65
     >   ,62.26,62.76,63.27,63.72,64.18,64.58,64.99,65.38,65.72,66.12
     >   ,66.42,66.80,67.11,67.38,67.69,67.95,68.23,68.49,68.73,68.92
     >   ,69.18,69.45,69.68,69.86,70.09,70.32,70.50,70.74,70.91,71.08
     >   ,71.19,71.36,71.49,71.57,71.66,71.73,71.75,71.77,71.81,71.77
     >   ,71.67,71.53,71.38,71.23,71.02,70.72,70.39,70.03,69.61,69.08
     >   ,68.58,68.01,67.33,66.58,65.77,64.94,64.05,63.10,62.03,60.94
     >   ,59.76,58.55,57.25,55.97,54.56,53.11,51.60,50.17,48.57,46.98
     >   ,45.35,43.81,42.14,40.49,38.74,37.05,35.40,33.73,32.02,30.42
     >   ,28.76,27.13,25.57,23.98,22.45,20.95,19.52,18.09,16.75,15.43
     >   ,14.17,12.98,11.83,10.75, 9.71, 8.76, 7.84, 6.98, 6.19, 5.45,
     >   4.78, 4.16, 3.60, 3.09, 2.62, 2.21, 1.85, 1.52, 1.25, 1.00,
     >   0.80, 0.63, 0.49, 0.37, 0.27 /

      NFILT = 0

      SELECTCASE (FILTERNAME)

      CASE ('U', 'JOHN U') 
         NFILT = 14
         FBLUE = 2800.
         FINCR = 100.
         FILT(1:NFILT) = UJ
         
      CASE ('B', 'JOHN B') 
         NFILT = 22 
         FBLUE = 3400.
         FINCR = 100.
         FILT(1:NFILT)  = BJ

      CASE ('V', 'JOHN V') 
         NFILT = 25
         FBLUE = 4700.
         FINCR = 100.
         FILT(1:NFILT)  = VJ

      CASE ('u', 'STROEM U') 
         NFILT = 28
         FBLUE = 3125.
         FINCR = 25.
         FILT(1:NFILT)  = US

      CASE ('v', 'STROEM V') 
         FBLUE = 3725.
         FINCR = 25.
         NFILT = 29
         FILT(1:NFILT) = VS

      CASE ('b', 'STROEM B') 
         FBLUE = 4325.
         FINCR = 25.
         NFILT = 29
         FILT(1:NFILT)  = BS

      CASE ('y', 'STROEM Y') 
         NFILT = 29
         FBLUE = 5000.
         FINCR = 25.
         FILT(1:NFILT)  = YS

      CASE ('J_2MASS') 
         NFILT = 194
         FBLUE = 10610.
         FINCR = 20.
         FILT(1:NFILT)  = J2MASS

      CASE ('H_2MASS') 
         NFILT = 90
         FBLUE = 14150.
         FINCR = 50.
         FILT(1:NFILT)  = H2MASS

      CASE ('Ks_2MASS') 
         NFILT = 90
         FBLUE = 19275.
         FINCR = 50.
         FILT(1:NFILT)  = Ks2MASS

      CASE ('Gaia_Gbp') 
         NFILT = 177
C wrong!         FBLUE = 3275.
         FBLUE = 3255.
         FINCR = 20.
         FILT(1:NFILT)  = Gaia3Gbp

      CASE ('Gaia_Grp') 
         NFILT = 152
C wrong!         FBLUE = 6160.
         FBLUE = 6130.
         FINCR = 30.
         FILT(1:NFILT)  = Gaia3Grp

      CASE ('Gaia_G') 
         NFILT = 145
C wrong!         FBLUE = 3280.
         FBLUE = 3230.
         FINCR = 50.
         FILT(1:NFILT)  = Gaia3G

      ENDSELECT

      IF (NFILT .EQ. 0) THEN
         WRITE (0,'(A)') '*** FATAL INTERNAL ERROR:'
         WRITE (0,'(A)') '*** Invalid Filter name: ' // FILTERNAME
         STOP '*** ERROR in subr. FILTERFUNCTIONS'

      ENDIF

C***  WAVELENGTH GRID of selected filter
      DO I=1, NFILT
         FLAM(I) = FBLUE + I * FINCR
      ENDDO

C******** Test Plot Facility
      IF (.FALSE.) THEN 
C*       Normalization: Maximum = 1
         FILTMAX = MAXVAL (FILT(1:NFILT))
         DO I=1, NFILT
           FILT(I) = FILT(I) / FILTMAX
         ENDDO

         CALL PLOTANFS (17, FILTERNAME, FILTERNAME
     >        ,'\CENTER\#l# [\A]'
     >        ,'\CENTER\Transmission'
     >        , .0 , .0 , .0 , .0 , .0 , .0
     >        , .0 , .0 , .0 , .0 , .0 , .0
     >        , FLAM, FILT, NFILT, 'SYMBOL=5 COLOR=2')
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
      SUBROUTINE FINDIND (IND, INDSTR, LEVEL, N, INDNUP,INDLOW,LASTIND, 
     >                    BFINDERR, LOW, NUP)
C***********************************************************************
C***  FIND LINE-INDEX "IND" FROM GIVEN LEVEL-NAMES
C***  - ACTION: IF INDEX-FIELD "INDSTR" CONTAINS A "?":
C***    READ NEXT INPUT LINE, WHICH MUST HAVE THE FOLLOWING FORMAT:
C***    UPPERLEVEL=.......... LOWERLEVEL=..........
C***    ELSE: READ IND FROM "INDSTR"
C***  - CALLED FROM: PREFORM < DECFORM < FORMAL
C***********************************************************************

      CHARACTER INDSTR*4, KARTE*80, ERRMES*80, LEVEL(N)*10
      DIMENSION INDNUP(LASTIND), INDLOW(LASTIND)
      LOGICAL BFINDERR

      BFINDERR = .FALSE.

      IF (INDSTR(1:1) .NE. '?' .AND. INDSTR(2:2) .NE. '?' .AND.
     $    INDSTR(3:3) .NE. '?' .AND. INDSTR(4:4) .NE. '?') THEN

C***     INDEX IS ASSUMED TO BE GIVEN AS NUMBER IN THE "INDSTR" FIELD
         READ (INDSTR, '(I4)', ERR=110) IND
         NUP = INDNUP(IND)
         LOW = INDLOW(IND)
         IF (IND .LE. 0 .OR. IND .GT. LASTIND) THEN
            ERRMES = 'LINE INDEX OUT-OF-RANGE: ' // INDSTR
            GOTO 100
         ENDIF

      ELSE

C***     READ NEXT INPUT LINE WHICH MUST GIVE THE LEVEL NAMES
         READ (2, '(A)') KARTE
         IF (KARTE( 1:10) .NE. 'UPPERLEVEL' .OR.
     $       KARTE(23:32) .NE. 'LOWERLEVEL') THEN
            ERRMES = 'INCORRECT INPUT CARD: ' // KARTE(:43)
            GOTO 100
         ENDIF

C***     FIND UPPER INDEX
         DO 2 I=2, N
            IF (KARTE(12:21) .EQ. LEVEL(I)) THEN
               NUP = I
               GOTO 3
            ENDIF
    2    CONTINUE
         ERRMES = 'UPPER LEVEL NOT FOUND: ' // KARTE(12:21)
         GOTO 100
    3    CONTINUE

C***     FIND LOWER INDEX
         DO 4 I=1, N
         IF (KARTE(34:43) .EQ. LEVEL(I)) THEN
            LOW = I
            GOTO 5
            ENDIF
    4    CONTINUE
         ERRMES = 'LOWER LEVEL NOT FOUND: ' // KARTE(34:43)
         GOTO 100
    5    CONTINUE

C***     FIND LINE INDEX
         IF (LOW .EQ. NUP) THEN
C***        lower and upper level have the same index
            IND = -1
            GOTO 7
         ELSE
            DO 6 I=1, LASTIND
            IF (INDLOW(I) .EQ. LOW .AND. INDNUP(I) .EQ. NUP) THEN
               IND = I
               GOTO 7
            ENDIF
    6       CONTINUE
         ENDIF

         WRITE (0,'(A)') 'LINE INDEX NOT FOUND FROM THE FOLLOWING CARD:'
         ERRMES = KARTE(:IDX(KARTE))
         GOTO 100

    7    CONTINUE
      ENDIF

      GOTO 200

C***  ERROR EXIT: 
  110 WRITE (0,'(A)') '*** ERROR: LINE Index entry invalid:'
      ERRMES = KARTE(:IDX(KARTE))
    
  100 CONTINUE
      BFINDERR = .TRUE.
      PRINT 101, ERRMES
  101 FORMAT ('THIS LINE or MULTIPLET SKIPPED! ', /, 
     >        ' NON-FATAL ERROR ',
     $          'DETECTED BY SUBROUTINE FINDIND:', /, A,/)

  200 RETURN

      END        
      SUBROUTINE FINDLDR (LOW, NUP, LEVEL, NOM, NCHARG, N, BFINDERR)
C***********************************************************************
C***  FIND LEVEL-INDICES FROM GIVEN LEVEL-NAMES
C***  - READ NEXT INPUT LINE, WHICH MUST HAVE THE FOLLOWING FORMAT:
C***    UPPERLEVEL=.......... LOWERLEVEL=..........
C***  - CALLED FROM: FORMAL > PREFORM
C***********************************************************************

      DIMENSION NOM(N), NCHARG(N)
      CHARACTER KARTE*80, ERRMES*80, LEVEL(N)*10
      LOGICAL BFINDERR

      BFINDERR = .FALSE.

C***  READ NEXT INPUT LINE WHICH MUST GIVE THE LEVEL NAMES
      READ (2, '(A)') KARTE
      IF (KARTE( 1:10) .NE. 'UPPERLEVEL' .OR.
     $    KARTE(23:32) .NE. 'LOWERLEVEL') THEN
            ERRMES = 'INCORRECT INPUT CARD: ' // KARTE(:43)
            GOTO 100
      ENDIF

C***  FIND UPPER INDEX
      DO 2 I=2, N
         IF (KARTE(12:21) .EQ. LEVEL(I)) THEN
            NUP = I
            GOTO 3
         ENDIF
    2 CONTINUE
      ERRMES = 'UPPER LEVEL NOT FOUND: ' // KARTE(12:21)
      GOTO 100

    3 CONTINUE

C***  FIND LOWER INDEX
      DO 4 I=1, N
         IF (KARTE(34:43) .EQ. LEVEL(I)) THEN
            LOW = I
            GOTO 5
         ENDIF
    4 CONTINUE
      ERRMES = 'LOWER LEVEL NOT FOUND: ' // KARTE(34:43)
      GOTO 100

    5 CONTINUE

C***  CHECK FOR CORRECT IONIZATION STAGES: CHARGE DIFFERENCE = 1
      IF ((NOM(LOW) .NE. NOM(NUP)) .OR. 
     >    (NCHARG(LOW) .NE. NCHARG(NUP)-1)) THEN
         ERRMES = '*** ERROR: INVALID COMBINATION OF LEVELS'
         GOTO 100
      ENDIF

      RETURN

C***  ERROR EXIT
  100 BFINDERR = .TRUE.

      PRINT 101, ERRMES
      WRITE (0,101) ERRMES
  101 FORMAT (/, 1X, 80('*'), /, ' NON-FATAL ERROR ON INPUT OPTIONS, ',
     $        'DETECTED BY SUBROUTINE FINDLDR:', /,
     $        1X, A, /, 1X, 80('*'), /)
      PRINT *, KARTE(:IDX(KARTE))

      RETURN

      END        
      SUBROUTINE FOLR (N, X, IX, A, IA)

      DIMENSION X(2), A(2)

      IF (IX .NE. 1 .OR. IA .NE. 1) THEN
        STOP 'ERROR IN FOLR: INCREMENTS NOT ALLOWED'
      ENDIF

      DO I=2, N
        A(I) = A(I) - X(I)*A(I-1)
      ENDDO

      RETURN
      END
C***  MAIN PROGRAM FORMAL  ****************************************************
      SUBROUTINE FORMAL
 
C*******************************************************************************
C***  FORMAL INTEGRAL IN THE OBSERVERS FRAME, YIELDING EMERGENT FLUX PROFILES
C*******************************************************************************
 
      IMPLICIT NONE
      CHARACTER(1) :: onechar
 
C***  DEFINE ARRAY DIMENSIONS
 
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM =          26 

      INTEGER, PARAMETER :: MAXLAP  =           15000
      INTEGER, PARAMETER :: MAXSUBL =           15000
      INTEGER, PARAMETER :: MAXIND  = 45000 + MAXSUBL
      INTEGER, PARAMETER :: MAXFEIND  =       1500 
C***  NDIM = Number-of-original-levels + MAXSUBLevels (without 2*)
      INTEGER, PARAMETER :: NDIM    =  1560 + MAXSUBL
      INTEGER, PARAMETER :: NFLDIM  =          300000   
      INTEGER, PARAMETER :: NFODIM  =          250000   
      INTEGER, PARAMETER :: MAXMOD  =               2
      INTEGER, PARAMETER :: NFDIM   =    2*NDIM + 400 
      INTEGER, PARAMETER :: MAXKONT =         NFDIM/2 
      INTEGER, PARAMETER :: NDDIM   =             100 
      INTEGER, PARAMETER :: NDADDIM =    2*NDDIM + 10 
      INTEGER, PARAMETER :: NPDIM   =             130 
      INTEGER, PARAMETER :: MAXXN   =            4000
      INTEGER, PARAMETER :: NREDMAX =            6000
      INTEGER, PARAMETER :: MAXSTRI =              20 
      INTEGER, PARAMETER :: MAXXDAT =              10 
      INTEGER, PARAMETER :: NPHIMAX =             500 

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 
      
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5    !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5    !xxl
      
C***  ARRAYS FOR TREATMENT OF LINE OVERLAPS:
      REAL, DIMENSION(MAXLAP) :: XLAMLAP, DELXLAP 
      INTEGER, DIMENSION(MAXLAP) :: INDLAP, IPOINTERPHITAB
      REAL, DIMENSION(NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL 
      REAL, DIMENSION(NDDIM,MAXLAP) :: OPALRAY, ETALRAY
      REAL, DIMENSION(MAXXN,MAXLAP) :: OPALFIN, ETALFIN
      REAL, DIMENSION(MAXLAP,NDDIM,MAXMOD) :: AVOIGT, GRIEMPAR

C***  ARRAYS FOR MULTIPLET HANDLING:
      INTEGER, DIMENSION(MAXSUBL) :: NSUBLOW, NSUBNUP

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 2850
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

C***  VECTORS FOR USE IN SUBR. ZONEINT
      REAL, DIMENSION(MAXXN) :: TAU, DTAU, WTAU, TAUC, DTAUC, WTAUC,
     >                          XCMFFINE, POROLENGTHFINE

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ENLTE
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST, KODATIND
      REAL, DIMENSION(NFDIM, MAXMOD) :: XLAMBDA
      REAL, DIMENSION(NDDIM, MAXMOD) :: OPA, ETA
      REAL, DIMENSION(NDDIM) :: ETANOTH, THOMSON, DELW, ADELW,
     >                          XJCIND, OPAFEFT, ETAFEFT
      INTEGER, DIMENSION(NDDIM) :: IWARN
      REAL, DIMENSION(NDDIM):: VEC_SECMOD
      REAL, DIMENSION(NDDIM,MAXMOD) :: RADIUS, ENTOT, T, RNE,
     >                                  VELO, VDU, GRADI, TAURCONT
      REAL, DIMENSION(NDDIM) :: RADIUS_MERGED 
      REAL, DIMENSION(NDDIM,MAXMOD) :: VDU_ORIG, 
     >               T_ORIG, RNE_ORIG, ENTOT_ORIG
      REAL, DIMENSION(MAXMOD) :: RCON, XMDOT, POPMIN
      REAL, DIMENSION(MAXATOM,MAXMOD) :: ABXYZ
      REAL, DIMENSION(NDDIM,NFDIM,MAXMOD) :: XJC
      REAL, DIMENSION(NDDIM,NPDIM) :: U, UK
      REAL, DIMENSION(3,NDDIM) :: EDDI ! For dimensioning, used in FORMCMF
      REAL, DIMENSION(NDADDIM) :: OPARAY, ETARAY, ETACRAY, 
     >                                    ZRAY, XCMF, RRAY
      REAL, DIMENSION(MAXXN) :: OPAFINE, ETAFINE, OPAFC, ETAFC,
     >                          SFINE, CSFINE, ZFINE
      REAL, DIMENSION(NFODIM) :: PROFILE, CPROFILE, DLAM, PROFILEC
      REAL, DIMENSION(MAXION) :: TFEEXC
      REAL, DIMENSION(NPDIM) :: A,WE, WX
      REAL, DIMENSION(NPDIM,NPDIM) :: BX
      LOGICAL, DIMENSION(NDIM,NDDIM,MAXMOD) :: ZERO_RATES

C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY]
C***  @TODO: GET RID OF THIS CONDITION ASAP!!!
      REAL, DIMENSION(NPDIM,NPDIM) :: B
      REAL, DIMENSION(NPDIM) :: C

      
      REAL, DIMENSION(NPDIM,MAXMOD) :: PGRID
      REAL, DIMENSION(NPDIM) :: PGRID_ORIG, PGRID_MERGED, EMINT_P
      INTEGER NP_ORIG, NDNP_ORIG
      REAL, DIMENSION(NDDIM,NPDIM,MAXMOD) :: ZGRID, ZGRID_ORIG
      REAL, DIMENSION(NDDIM,NPDIM) :: ZGRID_MERGED

      REAL, DIMENSION(NDDIM,NDIM,MAXMOD) :: POPNUM, POPNUM_ORIG
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO, 
     >                            ADDCON1, ADDCON2, ADDCON3
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW 
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW, MULTIIND


C***  ARRAYS FOR CMF (ELECTRON SCATTERING REDISTRIBUTION OPTION)
      REAL, DIMENSION(NDDIM) :: TA, TB, TC, UB, GA, H, QQ, SCMF,
     >                          VA, VB, PP, W0
      REAL, DIMENSION(NDDIM,NPDIM) :: V
      REAL, DIMENSION(NDDIM,NFLDIM) :: XJNUE, ETANCK, THOMCK,
     >                                 OPAK, ETAK
      REAL, DIMENSION(NDDIM,NFLDIM,MAXMOD) :: ETACK, OPACK, ETACCK
      REAL, DIMENSION(NFLDIM,MAXMOD) :: BCORE, DBDR
      REAL, DIMENSION(NFLDIM) :: XPLOT, YPLOT, YSCRATCH

C***  ARRAYS TO HANDLE DEPTH-DEPENDENT REDISTRIBUTION INTEGRAL
      REAL, DIMENSION(NDDIM) :: VDUEL, WREDI0
      REAL, DIMENSION(NREDMAX,NDDIM) :: WREDI

C***  X-RAY DATA
      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK

C***  Density Contrast
      REAL, DIMENSION(NDDIM,MAXMOD) :: DENSCON, FILLFAC 
      REAL, DIMENSION(NDDIM) :: ENTOTDENS, POROLENGTH
      REAL, DIMENSION(NDADDIM) :: POROLENGTHRAY

C***  Line profiles for pressure broadening
C***  NLDIMPHITAB = max. number of H + He lines in one spectral range 
C***  NFDIMPHITAB = max. number of frequency points in one half-profile 
      INTEGER, PARAMETER :: NLDIMPHITAB = 100
      INTEGER, PARAMETER :: NFDIMPHITAB = 2000
      REAL, DIMENSION 
     > (-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB,MAXMOD) :: PHITAB
      REAL, DIMENSION(2*NFDIMPHITAB+1) :: PHISCRATCH
      REAL, DIMENSION(MAXLAP) :: XMAXLIN
      CHARACTER(132), DIMENSION(0:MAXSTRI) :: STRING1
      CHARACTER(100), DIMENSION(MAXMOD) :: MODHEAD
      CHARACTER(100) :: MODHEADT
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2) :: EL_MIN
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8), DIMENSION(MAXLAP) :: LINPRO
      CHARACTER(10) :: XUNIT
      CHARACTER(20) :: FREQIN 
      CHARACTER(80) :: MANIPOP_OPTIONS, DD_VDOP_LINE, NOWIND_LINE, 
     >                 VDOPPLOT_LINE, TRANSDWLLINE, RCOROTLINE, 
     >                 MACROCLUMPLINE
      CHARACTER(256) :: KARTE, PATH_VCSSB, PATH_LEMKE_DAT

      LOGICAL :: PLOT, FIN, REDIS, IDENT
      LOGICAL :: BROAD, bBIGBANDLIMIT
      LOGICAL :: ABSWAV, BCONT, BNOCONT
      LOGICAL :: CORE, LINELIST, BCALIBRATED  
      LOGICAL :: BDD_VDOP, BPLOTVDOP, BPLOTVMIC, BMICROTURB
      INTEGER :: INDCUT, ND_MERGED, NP_MERGED

      CHARACTER(8) :: FUNIT
C***  IND_ORIGLEV(LEVEL) = index of "mother level" of multiplets 
C***                       (identical for normal lines)      
      INTEGER, DIMENSION(NDIM) :: IND_ORIGLEV

C***  For rotation    
      INTEGER, DIMENSION(NPDIM) :: NPHI
      REAL, DIMENSION(NPHIMAX,NPDIM) :: PHIWEIGHT, PHIARR    
 
C***  Variables needed for SECONDMODEL_DEFINE
      CHARACTER SECONDMODEL_LINE*1000, SECONDMODEL_PATH*400
      LOGICAL SECONDMODEL_CHANGED, IGNORE_DIFF
      REAL, DIMENSION(2) :: SECMOD_RRANGE

      REAL, DIMENSION(NPHIMAX) :: PHI_VEC

C***  VARIABLES READ BY FORMOSA
      INTEGER, DIMENSION(MAXMOD) :: ND, NP, NF, JOBNUM
      REAL, DIMENSION(MAXMOD) :: RSTAR, VDOP_MODEL, TEFF
      REAL, DIMENSION(NDDIM, MAXMOD) :: VMIC_MODEL

C***  FOR depth-dependent VDOP:
      REAL :: VMICFRAC_DEFAULT, VDOP, VDOP_FIRSTMOD

C*** vectors in km/s, Doppler units, and squared Doppler units are prepared
C     to avoid recalculation in long loops.
      REAL, DIMENSION(NDDIM,MAXMOD) :: 
     >                       DD_VMIC, DD_VMICDU, DD_VMICDU_ORIG
      REAL, DIMENSION (NDDIM, MAXATOM, MAXMOD) :: 
     >     DD_VDOP, DD_VDOPDU, DD_VDOP_ORIG, DD_VDOPDU_ORIG
      REAL, DIMENSION (MAXXN, MAXATOM) :: DD_VDOPDU_FINE_NORMFAC,
     >                  DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE
      REAL, DIMENSION(NDADDIM, MAXATOM) :: DD_VDOPDU_RAY

C*** Variables for limb-darkening plot
      REAL, DIMENSION(NFODIM) :: WEIGHT_LIMBDARK
      CHARACTER(132) :: LIMB_LINE
      LOGICAL :: BLIMB

      REAL :: DXLAM, DX_WINDROT, WEIGHTSUM

      REAL, DIMENSION(NPDIM) :: SUMFILT, SUMINT
      
      INTEGER, DIMENSION(MAXLAP) :: IND_ELLINE
      
C***  FROM PREPRAY -> OBSFRAM
      REAL BCOREL, DBDRL

C***  Intersection points with secondmodel domain
      REAL, DIMENSION(2, NPDIM, NPHIMAX) :: ZINTER

C***  To Store Feautrier Matrices ST(94*95,89)
      REAL, DIMENSION((NPDIM+1)*NPDIM,NDDIM) :: ST
      LOGICAL :: BELIFI

C***  SPZ is for the Input Option SET_POP_ZERO
C***  - can be sed to suppress specific lines 
      INTEGER, PARAMETER :: MAXSPZ = 10
      CHARACTER(10), DIMENSION(MAXSPZ) :: SPZ1, SPZ2

C***  Iron data variables      
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFELOW, IFENUP, INDFEACT
      REAL, DIMENSION(MAXFEIND) :: SIGMAACT, SIGMAACTUL, SIGNU3ACT,
     >                             SIGMAINT, SIGMAINTUL
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
ccc     SIGMAFEUL, SIGNU3FE unused - from Andreas
ccc      REAL, DIMENSION(INDEXMAX) :: SIGMAFE, SIGMAFEUL, SIGNU3FE
      INTEGER, DIMENSION(NDIM) :: INDRBS, IFRBSSTA, IFRBSEND, 
     >                            INDFESACT, IFES
      REAL, DIMENSION(NDDIM,NFLDIM,MAXMOD) :: OPAFE, ETAFE
      REAL, DIMENSION(NDDIM,MAXFEIND) :: OPAFEI, ETAFEI
      LOGICAL, DIMENSION(NDDIM) :: bFELASER
      LOGICAL :: BFECHECK, BFEMODEL, BIRONLINES, BFEWARNING,
     >           BAIRWAVELENGTHSET, BAIRWAVELENGTH, 
     >           bDDOPAFORMCMF, bDDFECONVOL, BVSINI_AT_RCOROT

      REAL :: VSINI, TAUMAX, XMAXMIN, DXMAX, VDOPFE,
     >        DXFE, XLAM0FE, VMAX, OSMIN, RMAX, XMAX, YMAX,
     >        XMAXBROAD, VSIDU, RANGEBLUE, VMAXDU,
     >        RANGERED, RCOROT, XLAM, XLAMLN,
     >        DISP, BWESEX, ALN, XRANGERED, XRANGEBLUE, XLAP,
     >        FREMAX, FREMIN, DUMMY, XCMFRED, XCMFBLUE, DXCMF, FNUEC,
     >        XLAMREF, DXOBS, FINCRI, BEGLAM, XOBS0, XO, PJPJ, 
     >        PWEIGHT, XLAM2, DISMO, FACLOG, FAC, CEMINT, EMINT, 
     >        ATMEAN, XMUL, AMACHL, TAUMINBROAD
      INTEGER LSOPA, LSDWL, LSPRO, IFIRST, IERR, IDUMMY,
     >        LPHISTA, LPHIEND, IMOD, JPFIRST, JPLAST, KANAL1, LINE,
     >        LASTIND, LASTFE, KWORDS, K,
     >        L, N, I, J, NMOD, NATOM, NAUTO, LASTKON, NA, NDN, NFL,
     >        NSTRING, IDTERM, NLPHITAB, IVERSION, NPHIROT_DEFAULT,
     >        NBLINE, NMULTI, NBL, IND, LBREF, IRANGEBLUE, IRANGERED,
     >        LOW, NUP, INDREF, IWARNJ0, JP, NFOBS, IFOBR,
     >        LPHI, IRAY, LTOT, NSTRI, LASTSELF, IVDOPSTATUS,
     >        JPFIRST_ORIG, JPLAST_ORIG, LPHISTA_ORIG, LPHIEND_ORIG, 
     >        NBFIRST, NBLAST, N_WITH_DRLEVELS 

      INTEGER, EXTERNAL :: IDX, ISRCHFLE, ISRCHFLT, ISRCHFGT, ISRCHFGE
     
C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  Physical constants      
      REAL, PARAMETER :: CLIGHT = 2.99792458E5      ! IN KM/SECOND
      REAL, PARAMETER :: CLIGHT2 = 2.99792458E18    ! in Angstroem/sec
      REAL, PARAMETER :: PARSEC = 3.08561E18    !PARSEC IN CM
      REAL, PARAMETER :: CONMV = 48.64 !Calibration CONSTANT FOR F-NUE IN MV
      REAL, PARAMETER :: PI = 3.14159 
      REAL, PARAMETER :: WPIINV = 0.564189583549       

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT   = 6    !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR   = 0    !write to wruniqX.cpr (stderr)

C***  Initialize values (only DATA statements here!)
      DATA LSOPA,LSDWL,LSPRO,IFIRST /-1, -1, -1, 1/
      DATA VSINI / 0. /
      DATA BVSINI_AT_RCOROT / .TRUE. /
C***  X-UNITS OF PLOT IN ANGSTROEM (ALTERNATIVELY: MICROMETER)
      DATA XUNIT / 'ANGSTROEM ' /

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  For progress-bar
      INTEGER NRAYDONE, NPHISUM, IMPATIENCE
C********* END OF DECLARATIONS *****************************************

C***  Write Link Data (Program Version) to CPR file
      WRITE (hCPR,'(2A)') '>>> FORMAL started: Program Version from '
     >                 ,LINK_DATE
      WRITE (hCPR,'(4A)') '>>> created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL
ccc not in gfortran      CALL TIME(TIM1)

C***  XN = NUMBER OF INTEGRATION INTERVALS ACROSS THE SCATTERING ZONE (REAL)
C***  TAUMAX = MAXIMUM OPTICAL DEPTH WHERE INTEGRATION IS TRUNCATED
C***  XMAX   = COMOVING-FRAME BANDWIDTH OF THE SCATTERING ZONE
C***  DXMAX  = MAXIMUM CMF FREQUENCY STEP IN THE SCATTERING ZONE
C***  TAUMINBROAD = MININUM TAU VALUE REQUIRED FOR EXTENDING THE 
C***                BROADENING CALCULATION RANGE IN SUBR. BANDWIDTH
C***  DEFAULT VALUES - MAY BE CHANGED BY INPUT OPTIONS!
      TAUMAX = 10.0
      XMAXMIN = 3.5
      DXMAX  = 0.3
      TAUMINBROAD = 0.01       !was 0.1 before Jan 2017 
      MANIPOP_OPTIONS = ' '      
      NMOD = 1

C***  Initialize "STRING COMMENT"
      FREQIN = ' '

C***  Initialize NOWIND option
      NOWIND_LINE = 'NONE'

C***  Initialize LIMBDARKENIG option
      LIMB_LINE = 'OFF'

C***  Initialize MACROCLUMP option
      MACROCLUMPLINE = 'NONE'

C***  Initialize porosity length
      DO L=1, NDDIM
         POROLENGTH(L) = .0
      ENDDO

C***  Initialize BELIFI; To Store Feautrier Matrices in Memory is now default
      BELIFI = .FALSE.

C***  OPEN CARDS-FILE
      OPEN (UNIT=2, FILE='CARDS', STATUS='UNKNOWN')

      CALL DATOM (NDIM, N, LEVEL, NCHARG, WEIGHT, 
     >            ELEVEL, EION, MAINQN,
     >            EINST, ALPHA, SEXPO,
     >            ADDCON1, ADDCON2, ADDCON3, 
     >            IGAUNT, COCO, KEYCBB, ALTESUM,
     >            INDNUP, INDLOW, LASTIND, MAXIND, MAXATOM, NATOM,
     >            ELEMENT, SYMBOL, NOM, KODAT, ATMASS, STAGE,
     >            SIGMATHK, SEXPOK, EDGEK, NFIRST,
     >            NLAST, NAUTO, MAXAUTO, LOWAUTO, WAUTO, EAUTO, AAUTO,
     >            IONAUTO, KRUDAUT, KONTNUP, KONTLOW, LASTKON, MAXKONT,
     >            IONGRND, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'FORMAL', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS)

      BIRONLINES = BFEMODEL    ! Default as found in DATOM
    
C***  KODATIND gives the Core Charge Number of an element
C***  This is needed to identify H and He for pressure broadening
      DO NA=1, NATOM
         KODATIND(NA) = 0
         DO J = 1, MAXATOM
            IF (NA .EQ. KODAT(J)) KODATIND(NA) = J
         ENDDO
         IF (KODATIND(NA) .EQ. 0) THEN
            WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
            STOP 'ERROR IN FORMAL'
         ENDIF
         IF (KODATIND(NA) .GT. MAXATOM) THEN
            WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
            STOP 'ERROR IN FORMAL'
         ENDIF
      ENDDO

C***  READING OF THE MODEL FILE
      IMOD=1
      CALL FORMOSA (ND(IMOD), RADIUS(1,IMOD), NP(IMOD), 
     >          PGRID(1,IMOD), ZGRID(1,1,IMOD),
     >          ENTOT(1,IMOD), RNE(1,IMOD), ABXYZ(1,IMOD), NATOM, 
     >          T(1,IMOD), VELO(1,IMOD), NF(IMOD),
     >          XLAMBDA(1,IMOD), GRADI(1,IMOD),
     >          POPNUM_ORIG(1,1,IMOD), 
     >          RSTAR(IMOD), VDOP_MODEL(IMOD), VMIC_MODEL(1,IMOD),
     >          JOBNUM(IMOD), N, NDDIM, NPDIM, NFDIM, 
     >          MODHEAD(IMOD), TEFF(IMOD),
     >          MAXXDAT, XDATA, XJC(1,1,IMOD), IMOD, 
     >          DENSCON(1,IMOD), FILLFAC(1,IMOD), TAURCONT(1,IMOD),
     >          POPMIN(IMOD), ZERO_RATES(1,1,IMOD), RCON(IMOD), NDIM, 
     >          XMDOT(IMOD) )

      VMAX = VELO(1,1)
      RMAX = RADIUS(1,1)
C***  Default (if not specified otherwise in FORMAL_CARDS):
      VDOP = VDOP_MODEL(1)

      NDN = ND(IMOD) * N
      DO I=1, NDN
         POPNUM(I,1,IMOD) = POPNUM_ORIG(I,1,IMOD)
      ENDDO

      CALL POPMIN_NULLING (ZERO_RATES(1,1,IMOD), POPNUM(1,1,IMOD), 
     >                     POPMIN(IMOD), ND(IMOD), N)

C***  The no. of core-intersecting impact parameters might be modified 
C***  in subr. ROTATION_PREP; 

C***  The radius-grid might be changed by subr. SECONDMODEL_PREP
C***  as well as bu subr. NOWIND
C***  therefore, P and Z are cloned and restored after one RANGE 
C***  has been completed 
      NP_ORIG = NP(1) 
      DO JP=1, NP_ORIG
         PGRID_ORIG(JP) = PGRID(JP,1)
      ENDDO

      ZGRID_ORIG = ZGRID

C***  Printout of Model Parameters

      WRITE (*,'(2A)') 'MODHEAD=', MODHEAD(IMOD)
      CALL PRI_PAR (TEFF(IMOD), RSTAR(IMOD), 
     >         VELO(1,IMOD), DENSCON(1,IMOD), XMDOT(IMOD))

C***  Open file for plots of the wind-rotation geometric grid
      OPEN (65, FILE='windgrid.plot', STATUS='UNKNOWN')

C***  OPEN FILE 1 = 'PLOT' FOR DIRECT PLOT TRANSFER
      KANAL1 = 1
      CALL JSYMSET ('G1','TRANSFER')
      CALL REMARK ('PLOT DATA TO BE ROUTED')
      OPEN (KANAL1, FILE='PLOT', STATUS='UNKNOWN')

C***  MASS STORAGE ON FILE 7 FOR THE FEAUTRIER MATRICES (CONT.)
      CALL OPENMS (7, IDUMMY, IDUMMY, 0, IERR)
 
C***  DEFAULT OPTIONS
      BROAD = .FALSE.
      IDENT = .TRUE.
      OSMIN = 0.05
      REDIS = .TRUE.
      BCONT = .FALSE.
      BWESEX = 1.
      FIN=.FALSE.
      DISP= .0
      BDD_VDOP = .FALSE.
      RCOROTLINE = 'DEFAULT'
      VMICFRAC_DEFAULT = 0.05
      IVERSION=0
      ABSWAV=.TRUE.
      LINELIST = .FALSE.
      BNOCONT = .FALSE.
      BCALIBRATED = .FALSE.
      FUNIT = 'FLAM10PC'
      BAIRWAVELENGTHSET=.FALSE.
      BAIRWAVELENGTH=.FALSE.
      BPLOTVDOP = .FALSE.
      BMICROTURB = .FALSE.
      BLIMB = .FALSE.
      DD_VDOP_LINE = ''
      VDOPPLOT_LINE = ''
      PATH_VCSSB     = 'default'
      PATH_LEMKE_DAT = 'default'
      bBIGBANDLIMIT = .TRUE.
      bDDOPAFORMCMF = .TRUE.
      bDDFECONVOL = .TRUE.
      IVDOPSTATUS = 1
      IGNORE_DIFF = .FALSE. 

C***  LOOP FOR EVERY SPECTRAL RANGE TO BE SYNTHESIZED    ---------------------
    1 CONTINUE

C***  DECFORM reads from FORMAL_CARDS all options till a new
C***  spectral range is opened with a LINE, BLEND or RANGE option
      CALL DECFORM (KARTE, LSOPA, RANGERED, RANGEBLUE,
     >              VSINI, DISP, REDIS, BWESEX,
     >              LSPRO, LSDWL, FIN, VDOP, IDTERM, 
     >              JPFIRST_ORIG, JPLAST_ORIG,
     >              IVERSION, TAUMAX, XMAXMIN, DXMAX, PATH_VCSSB,
     >              IDENT, OSMIN, BROAD, MODHEAD(1), JOBNUM(1),
     >              STRING1, MAXSTRI, NSTRING, ABSWAV, BCONT, FUNIT,
     >              LINELIST, 
     >              BIRONLINES, BNOCONT, SPZ1, SPZ2, MANIPOP_OPTIONS, 
     >              XUNIT, BCALIBRATED, FREQIN, MACROCLUMPLINE, 
     >              PATH_LEMKE_DAT, BAIRWAVELENGTHSET, BAIRWAVELENGTH, 
     >              MAXSPZ, TRANSDWLLINE, RCOROTLINE, 
     >              DD_VDOP_LINE, BPLOTVDOP, BMICROTURB,  
     >              LIMB_LINE, VDOPPLOT_LINE, bBIGBANDLIMIT, 
     >              NOWIND_LINE, TAUMINBROAD, 
     >              bDDOPAFORMCMF, bDDFECONVOL, 
     >              LPHISTA_ORIG, LPHIEND_ORIG,
     >              BVSINI_AT_RCOROT, DX_WINDROT, SECONDMODEL_LINE)

C***  If no more ranges requested, terminate the program
      IF (FIN) GOTO 20

      XMAX = XMAXMIN   
      
      IF (.NOT. BFEMODEL) BIRONLINES = .FALSE.

      WRITE (0,'(/,A)') '----------- Starting with range '//FREQIN

C***  If JPFIRST and/or JPLAST were specified: 
C***     make sure that theu fall into the allowed range
      JPFIRST = MAX(1,       JPFIRST_ORIG)
      JPFIRST = MIN(JPFIRST, NP(1)-1)
      JPLAST  = MIN(NP(1)-1, JPLAST_ORIG)
      JPLAST  = MAX(JPLAST, JPFIRST)

      LPHISTA = 1
      LPHIEND = 1


C***  If iron lines are disabled by the option NO-IRONLINES :
C***  WARNING issued here to the cpr-File, and by Subr. PRIPRO to formal.out 
      IF (.NOT. BIRONLINES  .AND. BFEMODEL) then
        BFEWARNING = .TRUE.
        WRITE(0,'(A)') 'WARNING: MODEL CONTAINS IRON LINES THAT ARE '
     >              // 'SUPPRESSED BY THE OPTION "NO-IRON LINES"'
      ELSE 
        BFEWARNING = .FALSE.
      ENDIF    

C***  Preparation of depth-dependent ("DD") VDOP. 
C***  DD is activatd by specifying VMIC 
C***  If DD is not activated: 
C***     VDOP is constant, taken from MODEL file, or
C***     overwritten by FORMAL_CARDS option VDOP 
C***  If DD is activated: VDOP refers to the smallest Doppler-broadening
C***     of any line at any depth for ensuring sufficient resolution

C***  VDOP_STRUCT must be REPEATED FOR THE SECOND MODE (different T)
      IMOD = 1
      CALL VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, 
     >   DD_VDOP(1,1,IMOD), VDOP, VMIC_MODEL(1,IMOD), VELO(1,IMOD), 
     >   T(1,IMOD), ND(IMOD), NDDIM, NATOM, MAXATOM,  
     >   DD_VDOPDU(1,1,IMOD), 
     >   VMICFRAC_DEFAULT, ATMASS, XMAX, XMAXMIN,
     >   SYMBOL, VDOPFE, BMICROTURB, BIRONLINES,
     >   DD_VMIC(1,IMOD), TAURCONT, RADIUS(1,IMOD), EL_MIN,  
     >   DD_VMICDU(1,IMOD), bDDFECONVOL, IVDOPSTATUS)

C***  Note for Multi-Model mode: all dimensionless frequencies
C***  refer to VDOP of model 1
      VSIDU  = VSINI / VDOP
      VMAXDU = VMAX / VDOP

C*******************************************************************
C***  Preparation in case of SECONDMODEL
C*******************************************************************
      IF (SECONDMODEL_LINE .NE. '')
     >  CALL SECONDMODEL_DEFINE (SECONDMODEL_LINE, NMOD, 
     >        SECONDMODEL_PATH, SECONDMODEL_CHANGED, 
     >        RMAX, SECMOD_RRANGE, IGNORE_DIFF)

C***  Error stop if DWL-Plot requested but multi-model active 
      IF (LSDWL .GT. 0 .AND. NMOD .GT. 1) THEN
        WRITE (0,*) 'DWL-PLOT NOT POSSIBLE FOR MORE THAN ONE MODEL'
        STOP 'ERROR detected in main program FORMAL'
      ENDIF

C***  Macroclumping 
C***  Note: in case of SECOND MODEL, the same POROLENGTH vector will be
C***        used there
      IF (MACROCLUMPLINE(:10) .EQ. 'MACROCLUMP') THEN
         CALL PREPMACROCLUMP (MACROCLUMPLINE, DENSCON, VELO, RADIUS,
     >            TAURCONT, ND, POROLENGTH)
      ENDIF

      IF (SECONDMODEL_CHANGED .AND. NMOD .EQ. 2) THEN
        
        CALL COPY_SECONDMODEL 
     >        (SECONDMODEL_PATH, IGNORE_DIFF, BIRONLINES)
        WRITE (0,*) 'SECOND MODEL COPIED'
        IMOD=2
        CALL FORMOSA (ND(IMOD), RADIUS(1,IMOD), NP(IMOD), 
     >          PGRID(1,IMOD), ZGRID(1,1,IMOD),
     >          ENTOT(1,IMOD), RNE(1,IMOD), ABXYZ(1,IMOD), NATOM, 
     >          T(1,IMOD), VELO(1,IMOD), NF(IMOD),
     >          XLAMBDA(1,IMOD), GRADI(1,IMOD),
     >          POPNUM_ORIG(1,1,IMOD),
     >          RSTAR(IMOD), VDOP_MODEL(IMOD), VMIC_MODEL(1, IMOD),
     >          JOBNUM(IMOD), N, NDDIM, NPDIM, NFDIM, 
     >          MODHEAD(IMOD), TEFF(IMOD),
     >          MAXXDAT, XDATA, XJC(1,1,IMOD), IMOD,
     >          DENSCON(1,IMOD), FILLFAC(1,IMOD), TAURCONT(1,IMOD),
     >          POPMIN(IMOD), ZERO_RATES(1,1,IMOD), RCON(IMOD), NDIM,
     >          XMDOT(IMOD) )

        VMAX = MAX (VELO(1,1),VELO(1,2))

        NDN = ND(IMOD) * N
        DO I=1, NDN
           POPNUM(I,1,IMOD) = POPNUM_ORIG(I,1,IMOD)
        ENDDO

        CALL POPMIN_NULLING (ZERO_RATES(1,1,IMOD), POPNUM(1,1,IMOD), 
     >                       POPMIN(IMOD), ND(IMOD), N)

C***    Printout of Model Parameters (second model)
        WRITE (*,'(2A)') 'SECOND MODEL READ FROM ', 
     >      SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))
        WRITE (*,'(2A)') 'MODHEAD=', MODHEAD(IMOD)
        CALL PRI_PAR (TEFF(IMOD), RSTAR(IMOD), 
     >         VELO(1,IMOD), DENSCON(1,IMOD), XMDOT(IMOD))

C***    Preparation of depth-dependent ("DD") VDOP. 
C***    This is repeated here for the SECOND MODEL (different T(r)!)
        VDOP_FIRSTMOD = VDOP
        CALL VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, 
     >   DD_VDOP(1,1,IMOD), VDOP, VMIC_MODEL(1,IMOD), VELO(1,IMOD), 
     >   T(1,IMOD), ND(IMOD), NDDIM, NATOM, MAXATOM,  
     >   DD_VDOPDU(1,1,IMOD), 
     >   VMICFRAC_DEFAULT, ATMASS, XMAX, XMAXMIN,
     >   SYMBOL, VDOPFE, BMICROTURB, BIRONLINES,
     >   DD_VMIC(1,IMOD), TAURCONT, RADIUS(1,IMOD), EL_MIN, 
     >   DD_VMICDU(1,IMOD), bDDFECONVOL, IVDOPSTATUS)

C***     VDOP should not be overwritten by SECOND MODEL
         VDOP = VDOP_FIRSTMOD 
      ENDIF
C***  End of preparations in case of SECONDMODEL *******************

C***  Manipulation of Popnumbers, in order to simulate the emission of a hot 
C***    component (wrh,  5-Feb-2001)
      IF (MANIPOP_OPTIONS .NE. ' ') THEN
         IF (NMOD .EQ. 1) THEN 
            CALL MANIPOP (ENLTE, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >              ABXYZ, NFIRST, NLAST, NATOM, POPNUM, RNE, ENTOT, 
     >              N, ND, MANIPOP_OPTIONS, DENSCON, T, level)
         ELSE
            WRITE (0,*) '*** MANIPOP and SECONDMODEL options ',
     >               'cannot be combined!'
            STOP ' *** ABORT in FORMAL ***'
         ENDIF
      ENDIF

C***  Manipulation of Popnumbers: Zero setting for specified levels
      DO IMOD=1, NMOD
         CALL SET_POP_ZERO (LEVEL, NDIM, N, SPZ1, SPZ2, POPNUM(1,1,IMOD),
     >                 ND(IMOD), NMOD, MAXSPZ)
      ENDDO

C***  INTRODUCING DIMENSIONLESS VELOCITY UNITS
      DO IMOD=1, NMOD
        DO L=1, ND(IMOD)
          VDU (L,IMOD)=VELO(L,IMOD) / VDOP
        ENDDO
      ENDDO

C*******************************************************************
C***  Reading of FORMAL-CARDS is now continued for the spectral range
C***  (current single LINE, or till -BLEND closes the range)
C***
C***  For the SECOND-MODEL mode, the atomic date in the FORMAL_CARDS
C***  are for both models; the DATOM data must be identical anyhow.
C***  Subr. MULTIPLE - MULTISPLI adds additional SUBLEVELS with relative  
C***  LTE population numbers to the POPNUM array; this requires T(L) 
C***  and thus an internal loop over IMOD. 
C***  Moreover PREFORM might encounter a DRTRANSIT option which adds 
C***  auto-ionizing levels, requiring T, ENTOT and RNE and 
C***  also an internal loop IMOD=1, NMOD 
C*******************************************************************

      CALL PREFORM (KARTE, N, ELEVEL, LINE, INDLOW, INDNUP, LASTIND,
     >                CLIGHT, VDOP, INDLAP, XLAMLAP, DELXLAP, ALN, 
     >                XLAM, NBLINE, MAXLAP, MAXIND, MAXATOM, 
     >                LEVEL, WEIGHT, EINST, NDIM, POPNUM,
     >                T, ND, NOM, NCHARG, EION, ENTOT, RNE,
     >                MAXSUBL, NSUBLOW, NSUBNUP, BROAD, 
     >                LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD, DENSCON, 
     >                MAINQN, MULTIIND, 
     >                NMULTI, DD_VDOP, NATOM,
     >                IND_ORIGLEV)

C***  The reference wavelenght XLAM is the wavelength of the first line
C***  specified in the current BLEND range. If the BLEND is void of
C***  lines, XLAM is taken as the middle of the range. 
      IF (NBLINE == 0) THEN
        XLAM = (RANGERED + RANGEBLUE) / 2.
      ENDIF

C***  X-Scale: Frequency in harmonic Doppler units referring to XLAM
      XLAMLN = ALOG (XLAM)

C***  Restrict wavelength range if RANGE option given,
C***  remove lines outside range
      IF (RANGERED .NE. RANGEBLUE) THEN
         XRANGERED  = ALOG (RANGERED  / XLAM) / ALN
         XRANGEBLUE = ALOG (RANGEBLUE / XLAM) / ALN

         IRANGERED = ISRCHFLE(NBLINE,XLAMLAP,1,RANGERED)
         NBLINE = NBLINE+1-IRANGERED
         DO I=1, NBLINE
            INDLAP(I) =  INDLAP(I+IRANGERED-1)
            XLAMLAP(I) = XLAMLAP(I+IRANGERED-1)
            DELXLAP(I) = DELXLAP(I+IRANGERED-1)
            LINPRO(I) =  LINPRO(I+IRANGERED-1)
            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  AVOIGT(I,L,IMOD) =  AVOIGT(I+IRANGERED-1,L,IMOD)
               ENDDO
            ENDDO 
         ENDDO
         IF (RANGEBLUE .GT. XLAMLAP(NBLINE)) THEN 
            IRANGEBLUE = ISRCHFLE(NBLINE,XLAMLAP,1,RANGEBLUE)
            NBLINE = IRANGEBLUE - 1
         ENDIF
         IF (NBLINE .LT. 1) THEN
            WRITE (0,*) '*** WARNING: RANGE CONTAINS NO LINES'
         ENDIF
      ELSE
         XRANGERED  = .0
         XRANGEBLUE = .0
      ENDIF

C***  ONLY LINELIST TO BE GENERATED -- NO CALCULATION!
      IF (LINELIST) THEN
        WRITE (*,*) ' ******** LISTONLY OPTION, NO CALCULATION ! *****'
C*      If RANGE option is used together with LISTONLY,
C*      the definition of FREMIN, FREMAX would be missing
        IF (RANGERED .NE. RANGEBLUE) THEN
           FREMIN = XRANGERED
           FREMAX = XRANGEBLUE
        ENDIF
        GOTO 10
      ENDIF

C***  To each line will be assigned the CMF bandwidth that must be
C***  covered; this will be done in Subr. STARKBROAD -> BANDWIDTH
C***  and stored in the vector XMAXLIN
C***  XMAXBOAD is initialzed here and might be incremented in BANDWIDTH
      XMAXBROAD = XMAX

C***  LOOP OVER DETECTED AND ALL BLENDING LINES  -----------------------
      DO 25 NBL=1,NBLINE
        IND=INDLAP(NBL)
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)

C***    IND_ELLINE(LEVEL) gives corresponding element index = NA.         
        IND_ELLINE(NBL) = NOM(LOW)
        NA = NOM(LOW)
        XLAP=XLAMLAP(NBL)
C***    Initialize bandwidth of all lines to XMAX for the case that
C***    STARKBROAD is not called (no broadening requested)
        XMAXLIN(NBL) = XMAX
        DO IMOD=1, NMOD
          CALL LIOP (EINST(NUP,LOW), WEIGHT(LOW), WEIGHT(NUP), LOW,NUP,
     >               ND(IMOD), XLAP, ENTOT(1,IMOD), POPNUM(1,1,IMOD),
     >               RSTAR(IMOD), OPAL(1,NBL,IMOD), ETAL(1,NBL,IMOD),
     >               VDOP)
        ENDDO


C***    OPTION: PRINT LINE OPACITIES ***************************************
        IF (LSOPA.GT.0) THEN
C***      Not allowed for Multi-Model Model !
C***      Not accounting for clumping!
          IF (NMOD .GT. 1) THEN
            WRITE (0,*) 'ERROR : PRINTING OF LINE OPACITIES NOT ',
     >                  'POSSIBLE,'
            WRITE (0,*) '        IF MORE THAN ONE MODEL IS SPECIFIED'
            STOP 'ERROR IN SUBR. FORMAL'
          ENDIF

C***    CALL COOP AND BACKJC ONLY FOR THE USE OF OPA, ETA, THOMSON, XJCIND
C***    IN THE ROUTINE PRIOPAL
          CALL COOP (XLAP, ND, T(1,1), RNE(1,1),
     >               POPNUM(1,1,1), ENTOT(1,1), RSTAR(1),
     >               OPA(1,1), ETANOTH, THOMSON, IWARN, MAINPRO, 
     >               MAINLEV, NOM, KODAT,
     >               NDIM, N, MAXATOM, LEVEL, NCHARG, WEIGHT, 
     >               ELEVEL, EION, EINST,
     >               ALPHA, SEXPO,
     >               ADDCON1, ADDCON2, ADDCON3,
     >               IGAUNT, SIGMATHK, SEXPOK, EDGEK, 
     >               0, NF(1), DUMMY,
     >               RADIUS, KONTNUP, KONTLOW, LASTKON,XDATA)
          CALL BACKJC (XJC(1,1,1), ND, NF(1), 
     >                 XJCIND, XLAMBDA(1,1), 
     >                 XLAP, RADIUS)
          DO 117 L=1, ND(1)
            ETA(L,1) = ETANOTH(L) + 
     >                    OPA(L,1) * THOMSON(L) * XJCIND(L)
  117     CONTINUE

          CALL PRIOPAL(KARTE, XLAP, ND, 
     >                   OPA(1,1), OPAL(1,NBL,1), 
     >                   ETA(1,1), ETAL(1,NBL,1),
     >                   RADIUS, JOBNUM(1), 
     >                   LSOPA, MODHEAD(1))
        ENDIF 
C***    END OF PRINT-OPACITY-BLOCK *************************************

   25 CONTINUE
C***  END LOOP OVER ALL LINES IN BLEND -------------------------

C***********************************************************************
C***  FORMCMF: Co-moving frame calculation of the frequeny
C***           redistribution by electron scattering
C***********************************************************************

      DO IMOD=1, NMOD
         IF (IMOD .EQ. 2) WRITE (0,'(A)') 
     >      'NEXT: FORMCMF for SECOND MODEL'

         VMAX = AMAX1(VMAX, VELO(1,IMOD))
         CALL FORMCMF (TEFF(IMOD), U, UK, NDDIM, NFLDIM,
     >        ZGRID(1,1,IMOD), OPA(1,IMOD), ETA(1,IMOD), ETANOTH, 
     >        ETACK(1,1,IMOD), ETACCK(1,1,IMOD), 
     >        OPACK(1,1,IMOD), ETANCK,
     >        XCMFRED, XCMFBLUE, DXCMF, XMAX, RMAX, VMAXDU,
     >        THOMSON, THOMCK, OPAL(1,1,IMOD),
     >        ETAL(1,1,IMOD),  
     >        RADIUS(1,IMOD), ND(IMOD), NP(IMOD), PGRID(1,IMOD),
     >        TA, TB, TC, UB, GA, H, QQ, SCMF, V, VA, VB, PP,
     >        BCORE(1,IMOD), DBDR(1,IMOD), VDU(1,IMOD), 
     >        GRADI(1,IMOD), VDOP,
     >        W0, NBLINE, DELXLAP, XRANGERED, XRANGEBLUE,
     >        XJNUE, OPAK, ETAK,
     >        XPLOT, YPLOT, XJC(1,1,IMOD), NFDIM, XJCIND,
     >        LINE, MODHEAD(IMOD), JOBNUM(IMOD), NREDMAX, WREDI,
     >        WREDI0, VDUEL, BWESEX,
     >        XLAMBDA(1,IMOD), XLAM, T(1,IMOD), 
     >        RNE(1,IMOD), POPNUM(1,1,IMOD), 
     >        ENTOT(1,IMOD), RSTAR(IMOD), MAINPRO,
     >        MAINLEV, NOM, KODAT, NDIM, N, 
     >        MAXATOM, LEVEL, NCHARG,
     >        WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO, ADDCON1,
     >        ADDCON2, ADDCON3, IGAUNT, 
     >        SIGMATHK, SEXPOK, EDGEK, NF(IMOD),
     >        KONTNUP, KONTLOW, LASTKON, XDATA, REDIS, LBREF, INDREF,
     >        NPDIM, A, B, C, WE, BX, WX, EDDI, IWARNJ0,
     >        IWARN, FNUEC, LSDWL, ST, BELIFI, 
     >        DENSCON(1,IMOD), FILLFAC(1,IMOD), ENTOTDENS, 
     >        bBIGBANDLIMIT,
C *** IRON: Additional Parameter used in SUBROUTINE FECHECK called from FORMCMF
     >        INDRB, IFRBSTA, IFRBEND, LASTFE,
     >        CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >        INDFEACT, BFECHECK,
C *** IRON: Additional Parameter used in SUBROUTINE CMFFEOP called from FORMCMF
     >        SIGMAFE, OPAFE(1,1,IMOD), ETAFE(1,1,IMOD), 
     >        IFENUP, IFELOW, BIRONLINES,
     >        SIGMAACT, OPAFEI, ETAFEI, NFL, BAIRWAVELENGTHSET, 
     >        BAIRWAVELENGTH, VSIDU,
     >        DD_VDOPDU(1,1,IMOD), NATOM, YSCRATCH,
     >        MAXLAP, IND_ELLINE, bDDOPAFORMCMF, bDDFECONVOL)   

      ENDDO
      WRITE (0,*) 'FORMCMF finished for '//FREQIN
C***********************************************************************

C***  A second model, if involved, might have a different radius grid. 
C***  All relevant quantities are now interpolated to the grid of the 
C***  main model. 
C***  VDU and DD_VDOPDU are not only specific for the current RANGE
C***  and therefore saved before being scaled, and will be restored 
C***  when the current range is done. 
      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         VDU_ORIG       = VDU
         DD_VDOP_ORIG   = DD_VDOP
         DD_VDOPDU_ORIG = DD_VDOPDU
         DD_VMICDU_ORIG = DD_VMICDU
         T_ORIG         = T
         ENTOT_ORIG     = ENTOT
         RNE_ORIG       = RNE
      ENDIF

      IF (NMOD .EQ. 2) THEN
         CALL MERGE_RGRID (RADIUS, NDDIM, MAXMOD, ND, RADIUS_MERGED,
     >                     ND_MERGED, PGRID, NP, NPDIM, PGRID_MERGED, 
     >                     NP_MERGED, ZGRID_MERGED, SECMOD_RRANGE)
C***     Note: Number of impact parameters NP might have been enhanced 
C***           by subr. MERGE_RGRID
         JPFIRST = MAX(1,           JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST,     NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)

      ELSE
C***  "Normal" case without SECOND MODEL 
         ND_MERGED     = ND(1)
         NP_MERGED     = NP(1)
         RADIUS_MERGED = RADIUS(1:ND(1), 1)
         PGRID_MERGED  = PGRID(1:NP(1), 1)
         ZGRID_MERGED  = ZGRID(1:NDDIM, 1:NPDIM,1)
      ENDIF

C***  Default (no WINDROT, no SECONDMODEL)
      DO JP = 1, NP_MERGED
            NPHI(JP) = 1
            PHIWEIGHT(1,JP) = 1.
      ENDDO


C***  NOWIND option: omit (part of) the wind in the formal integral;
C***     this subroutine modifies RADIUS_MERGED etc. 
      IF (NOWIND_LINE .NE. 'NONE') THEN 
         CALL NOWIND (NOWIND_LINE, RCON, NATOM, ATMASS, ABXYZ,
     >                ND, RNE, T, RADIUS, VELO, VDOP, TAURCONT, 
     >                ND_MERGED, RADIUS_MERGED, PGRID_MERGED, 
     >                NP_MERGED, ZGRID_MERGED, IERR)
         IF (IERR .GT. 0) GOTO 2
C***     Note: Number of impact parameters NP might have been enhanced 
C***           by subr. MERGE_RGRID
         JPFIRST = MAX(1,           JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST,     NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)
      ENDIF

C*******************************************************************
C***  Preparation in case of wind rotation
C*******************************************************************
      IF (VSINI .GT. 0.) THEN
         WRITE(0,'(A,F7.1)') 'Wind Rotation: VSINI [km/s] =', VSINI
C***     ROTATION_PREP evaluates RCOROT, inserts more core-intersecting 
C***     points in PGRID_MERGED, re-calculates ZGRID_MERGED, 
C***     and establishes / modifies the array of azimuthal angles PHIARR 
        
         CALL ROTATION_PREP (RCOROTLINE, RADIUS, VELO, TAURCONT, 
     >                    ND, RCOROT, NPHI, VSIDU, XMAX, 
     >                    ND_MERGED, RADIUS_MERGED,
     >                    NP_MERGED, PGRID_MERGED, ZGRID_MERGED,
     >                    NDDIM, NPDIM, NPHIMAX, PHIWEIGHT, PHIARR, 
     >                    BVSINI_AT_RCOROT, DX_WINDROT)        

C***  If JPFIRST and/or JPLAST were specified: 
C***     make sure that theu fall into the allowed range
C***     Note: Number of impact parameters NP_MERGED might have been enhanced 
C***           by subr. ROTATION_PREP
         JPFIRST = MAX(1,       JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST, NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)

         CALL PLOT_WINDROT_GRID (PGRID_MERGED, NPDIM, JPFIRST, JPLAST, 
     >      LPHISTA_ORIG, LPHIEND_ORIG, NPHI, 
     >      PHIARR, NPHIMAX, XPLOT, YPLOT)

          IF (LPHISTA_ORIG .GT. 0 .OR. LPHIEND_ORIG .LT. 999) 
     >       WRITE (0,'(A,/,A, 2I4)') 
     >       '**** TEST RUN WITH RESTRICTED RANGE OF ANGLE INTEGRAL:', 
     >       '**** LPHISTA, LPHIEND =',  LPHISTA_ORIG, LPHIEND_ORIG

      ENDIF

      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         CALL RESCALE_SECMOD (NDDIM, NFLDIM, MAXLAP, MAXMOD, MAXATOM,
     >       ND, NATOM, NFL, NBLINE, RADIUS, RADIUS_MERGED, ND_MERGED,
     >       VDU, OPAL, ETAL, ETACK, ETACCK, OPACK, OPAFE, ETAFE,
     >       DD_VDOPDU, DD_VDOP, DD_VMICDU, T, ENTOT, RNE, VEC_SECMOD, 
     >       RSTAR, NMOD)
      ENDIF

      IF (NMOD .EQ. 2) THEN
         CALL SECONDMODEL_PREP (ZINTER, NPHI, PGRID_MERGED, NP_MERGED, 
     >        NPDIM, NPHIMAX, PHIARR, PHIWEIGHT, PHI_VEC, 
     >        SECONDMODEL_LINE, 
     >        JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG)

C***  Test output in case of a single ray
         IF (JPFIRST == JPLAST .AND. LPHISTA_ORIG == LPHIEND_ORIG) THEN
            WRITE (0,'(A,I3,A,F8.3,A,I3,A,F8.3)') 
     >      '*** Single Ray at P(', JPFIRST, ')=', 
     >      PGRID_MERGED(JPFIRST), 
     >      ' and PHI(', LPHISTA_ORIG, ')=',
     >      PHIARR(LPHISTA_ORIG,JPFIRST)*180./PI
            WRITE (0,'(A,2F8.3)') 'Intersection points with second model:'
     >       // ' Z1, Z2 =', ZINTER(1,JPFIRST,LPHISTA_ORIG), 
     >                       ZINTER(2,JPFIRST,LPHISTA_ORIG)
         ENDIF
      ENDIF

C*************************************************************
C***  Line broadening preparation (if requested)
C************************************************************************
      IF (BROAD) THEN 
         DO IMOD=1, NMOD
C***        initialize counter for line-broadening profile table
            NLPHITAB = 0
            DO NBL=1, NBLINE
               IND=INDLAP(NBL)
               LOW=INDLOW(IND)
               NUP=INDNUP(IND)
               NA = NOM(LOW)
               CALL STARKBROAD (KODATIND, NOM, 
     >          PHITAB(-NFDIMPHITAB,1,1,IMOD), 
     >          NFDIMPHITAB, NLDIMPHITAB, NLPHITAB,
     >          ND_MERGED, RADIUS_MERGED, OPA, OPAL(1,NBL,1), LINPRO(NBL), 
     >          AVOIGT, NBL, MAXLAP, IPOINTERPHITAB(NBL),
     >          XMAX, XMAXBROAD, XMAXLIN(NBL), NDDIM, DXMAX, PATH_VCSSB, 
     >          LOW, NUP, PHISCRATCH, XLAMLAP(NBL), ALN, T(1,IMOD), 
     >          ENTOT(1,IMOD), RNE(1,IMOD), LEVEL, MAINQN, 
     >          NCHARG, POPNUM(1,1,IMOD), N, VDOP, ELEVEL, EION, 
     >          PATH_LEMKE_DAT, DD_VDOP(1,NA,IMOD), 
     >          DD_VDOPDU(1,NA,IMOD), 
     >          IND_ORIGLEV, BDD_VDOP, DD_VMICDU(1,IMOD), 
     >          GRIEMPAR, TAUMAX, TAUMINBROAD, IMOD, MAXMOD, EINST, NDIM)

C***          Ensure XMAX covers maximum broadening  
              XMAX = MAX(XMAX, XMAXBROAD)

            ENDDO
         ENDDO
       ENDIF
C************************************************************************

ccc   **** test output *****************************
      IF (.false.) THEN
      do imod=1, nmod
         write (0,*) 'GRIEMPAR of MODEL', imod
         do nbl=1, nbline
            do l=1, nd(imod)
               if (griempar(nbl,l,imod) .gt. .0) then
                  onechar = '+'
               elseif (griempar(nbl,l,imod) .lt. .0) then
                  onechar = '-'
               else
                  onechar = '0'
               endif
               write (0,'(a1,$)') onechar
            enddo
            write (0,'(1x,I2,2x,A)') nbl, linpro(nbl)
         enddo
      enddo
      ENDIF
ccc   **** end of test output *****************************

C***  DEFINING THE OBSERVER'S FRAME FREQUENCY RANGE (FREMIN, FREMAX) 
C***  Note: this range must be smaller than covered by FORMCMF
      IF (RANGERED .NE. RANGEBLUE) THEN
         FREMIN = XRANGERED
         FREMAX = XRANGEBLUE
         RANGERED  = .0
         RANGEBLUE = .0
      ELSE
         FREMIN = XCMFRED  + 1.1*VMAXDU
         FREMAX = XCMFBLUE - 1.1*VMAXDU
      ENDIF

      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'XCMFBLUE=',XCMFBLUE," (",XLAM * EXP( XCMFBLUE * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'FREMAX  =',FREMAX," (" ,XLAM * EXP( FREMAX * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'FREMIN  =',FREMIN," (",XLAM * EXP( FREMIN * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'XCMFRED =',XCMFRED," (",XLAM * EXP( XCMFRED * ALN )," Ang)"

      IF (LSDWL.GT.0) THEN
          CALL PRIDWL (VELO(1,1), GRADI(1,1), VDOP,
     >                 POPNUM(1,1,1), ENTOT(1,1),
     >                 RADIUS, RSTAR(1), OPAL(1,1,1),
     >                 EINST, WEIGHT, FNUEC, XLAP, 
     >                 LOW, NUP, ND, NDIM, DELW,
     >                 LSDWL, OPA(1,1), THOMSON, ADELW,
     >                 JOBNUM(1), MODHEAD(1))
        GOTO 501
      ENDIF
           
C***  Prepare LOOP FOR EVERY OBSERVER'S FRAME-FREQUENCY  -----------------------
C***
C***  Define DXOBS = frequency spacing in Doppler units
      DXOBS = DXMAX
      NFOBS = IFIX((FREMAX-FREMIN)/DXOBS) + 2
      IF (NFOBS .GT. NFODIM) THEN
            WRITE (0,*) '********************************************'
            WRITE (0,*) 'ERROR: NFODIM TOO SMALL, REQUIRED:', NFOBS
            WRITE (0,*) '********************************************'
            WRITE (*,*) '********************************************'
            WRITE (*,*) 'ERROR: NFODIM TOO SMALL, REQUIRED:', NFOBS
            WRITE (*,*) '********************************************'
          STOP 'FATAL ERROR detected in FORMAL'
          ENDIF

C**  DEFINING ZERO-POINT AND INCREMENT OF THE OBSERVER'S FRAME FREQUENCY
      XOBS0 = FREMAX + DXOBS

      WRITE (0,'(1X,A,I6)')
     >      'Number of Observers-frame frequencies : ', NFOBS

      DO K = 1, NFOBS
          PROFILE(K)=0.
         CPROFILE(K)=0.
         PROFILEC(K)=0.
      ENDDO

C***  Decode LIMBDARKENING options and 
C***  prepare limb-darkening weights if requested
      IF (LIMB_LINE .NE. 'OFF') 
     >   CALL LIMBDARK_PREP (LIMB_LINE, XOBS0, DXOBS, NFOBS, XLAM, ALN, 
     >                      WEIGHT_LIMBDARK, BLIMB)

C***  Prepare total number of rays (for progress-bar)
      NPHISUM = 0
      DO JP= JPFIRST, JPLAST
         NPHISUM = NPHISUM + NPHI(JP)
      ENDDO
      IMPATIENCE = MAX (1, NPHISUM * NFOBS / 50)

      WRITE (0,'(A,I4,I4,I4)') 'JPFIRST, JPLAST, NP:', 
     >                          JPFIRST, JPLAST, NP_MERGED
      WRITE (0,'(A,I6)') 'Total number of azimuth-angle points: ',
     >               NPHISUM
      WRITE (0,'(A,/,A)') 'Progress Bar for this range:', 
     > '0% 10% 20%  30%  40%  50%  60%  70%  80%  90% 100%'

      NRAYDONE = 0


C***  LOOP OVER IMPACT PARAMETERS *********************************
      DO 14 JP=JPFIRST, JPLAST
         PJPJ = PGRID_MERGED(JP) * PGRID_MERGED(JP)
C***     IRAY = starting position of depth vector with impact parameter JP
C***            within the z-array if used as 1-dimensional field 
         IRAY=ND_MERGED * (JP-1) + 1

C***     for limb-darkening function
         EMINT_P(JP) = .0

C***     Prepare Loop over azimuth angles 
C***      - if no wind rotation and no second-model: only 1 point) 
C***      - if LPHISTA or LPHIEND are specified, make sure that 
C***        that they fall into the allowed range at current JP
          IF  (LPHISTA_ORIG .GT. 0 .OR. LPHIEND .LT. 999) THEN
              LPHISTA = MAX(1,       LPHISTA_ORIG)
              LPHISTA = MIN(LPHISTA, NPHI(JP))
              LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
              LPHIEND = MAX(LPHIEND, LPHISTA)
          ELSE
              LPHISTA = 1
              LPHIEND = NPHI(JP)
          ENDIF

C***     Normalization of PHIWEIGHT to unity
C***       i.e.phi integral always returns the average over the P-ring
C***      (or part of that ring, if LPHISTA or LPHIEND are restricted)

         WEIGHTSUM = .0
         DO LPHI=1, NPHI(JP)
             WEIGHTSUM = WEIGHTSUM + PHIWEIGHT(LPHI,JP)
         ENDDO

         DO LPHI=1, NPHI(JP)
             PHIWEIGHT(LPHI,JP) = PHIWEIGHT(LPHI,JP) / WEIGHTSUM
         ENDDO

C***     Loop over azimuth angles ***********************************
         DO 21 LPHI = LPHISTA, LPHIEND

C***         Loop over obsframe-frequencies  ************************
             DO 13 K=1, NFOBS 

             XO=XOBS0 - K*DXOBS
             XLAMREF = EXP (XLAMLN + XO * ALN)
             DLAM(K)= XLAMREF - XLAM

             CALL PREPRAY (ZGRID_MERGED(IRAY,1), COS(PHIARR(LPHI,JP)),
     >                PGRID_MERGED, ND_MERGED, NDDIM, NP_MERGED, JP,XO, 
     >                LTOT, PWEIGHT, 
     >                CORE, VDU, 
     >                RADIUS_MERGED, OPAL, ETAL, RRAY, OPARAY, ETARAY, 
     >                ETACRAY, OPALRAY, ETALRAY, NFLDIM, ZRAY, XCMF,
     >                NDADDIM, NBLINE, MAXLAP, REDIS, ETACK, ETACCK, OPACK, 
     >                     XCMFRED, XCMFBLUE, DXCMF, BCORE, BCOREL, 
     >                     DBDR, DBDRL, K, POROLENGTH, POROLENGTHRAY,
     >                     RCOROT, VSIDU, DD_VDOPDU_RAY, DD_VDOPDU, NATOM, 
     >                     NBFIRST, NBLAST, 
     >                     MAXATOM, BVSINI_AT_RCOROT, NMOD, MAXMOD, 
     >                     ZINTER(1,JP,LPHI), XMAX, DELXLAP)
C***         If no p-integral is calculated (test case), omit
C***         integration weight p*dp
             IF (JPFIRST .EQ. JPLAST) PWEIGHT = 1.

             CALL OBSFRAM (LTOT, CORE, XMAX, XMAXLIN, 
     >               EMINT, CEMINT, BCOREL, DBDRL,
     >               TAUMAX, PJPJ, ZFINE,  
     >               OPAFINE, OPAFC, OPALFIN, ETAFINE,
     >               ETAFC, ETALFIN, SFINE, CSFINE, RRAY, 
     >               OPARAY, OPALRAY, ETARAY, ETACRAY, 
     >               ETALRAY, ZRAY, XCMF, MAXXN,
     >               NDADDIM, NDDIM, DELXLAP, NBLINE,
     >               IVERSION, LINPRO, AVOIGT,
     >               TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >               XCMFFINE, POROLENGTHFINE, MAXLAP, DXMAX, 
     >               BIRONLINES, OPAFE, ETAFE, NFLDIM, 
     >               XCMFBLUE, XCMFRED, DXCMF, POROLENGTHRAY,
     >               PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >               RADIUS_MERGED, ND_MERGED, 
     >               DD_VDOPDU_RAY, NATOM, NBFIRST, NBLAST,
     >               IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >               DD_VDOPDU_FINE_SQRD, 
     >               DD_VDOPDU_FINE,
     >               GRIEMPAR, KODAT, VDOP, 
     >               INDCUT, ZINTER(1,JP,LPHI), NMOD, MAXMOD, K)

           PROFILEC(K)= PROFILEC(K) +  EMINT*PHIWEIGHT(LPHI,JP)*PWEIGHT
           CPROFILE(K)= CPROFILE(K) + CEMINT*PHIWEIGHT(LPHI,JP)*PWEIGHT
           IF (BLIMB) EMINT_P(JP) = EMINT_P(JP) + 
     >                EMINT*PHIWEIGHT(LPHI,JP)*WEIGHT_LIMBDARK(K)

C***       Progress bar
           NRAYDONE = NRAYDONE + 1
           IF ( (NRAYDONE / IMPATIENCE) * IMPATIENCE .EQ. NRAYDONE) THEN
              WRITE (0,'(A,$)') 'X'
           ENDIF

   13      CONTINUE
C***       End of frequency loop ********************************************

   21     CONTINUE
C***      End of phi-loop ***********************************************

   14 CONTINUE
C***  End of p-loop *************************************************

      WRITE (0,'(/,A)') 'Range done: ' // FREQIN 
C*************************************************

C***  Normalization: 
C***   PROFILEC = calibrated flux
C***   PROFILE  = normalized flux
C***  CPROFILE  = contimuum  flux
      DO K=1, NFOBS
         PROFILE(K) = PROFILEC(K) / CPROFILE(K)
      ENDDO
 
C***  Transformation of the wavelength scale to air, if appropriate
      IF (BAIRWAVELENGTH) THEN
        DO K=1, NFOBS
          XO=XOBS0 - K*DXOBS
          XLAMREF = EXP (XLAMLN + XO * ALN)
          XLAM2   = XLAMREF * XLAMREF
          DLAM(K) = DLAM(K) - XLAMREF*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
        ENDDO
      ENDIF

C***  OUTPUT FOR THE DETECTED LINE: 
C***  WRITE BUFFERED STRING ARRAY WITH INDIVIDUAL COMMENTS 
C***  (OPTION 'STRING') TO FILE 'PLOT' (KANAL = 1)
      DO 11 NSTRI=0, MIN0(NSTRING,MAXSTRI)
        WRITE (KANAL1,'(A)') STRING1(NSTRI)(:IDX(STRING1(NSTRI)))
        IF (NSTRI .GT. 0) 
     >   WRITE (*     ,'(A)') STRING1(NSTRI)(:IDX(STRING1(NSTRI)))
   11 CONTINUE
      LOW=INDLOW(LINE)
      NUP=INDNUP(LINE)
      IF (LSDWL .LE. 1) THEN
        CALL TRAPLO (1, PROFILE, DLAM, NFOBS,
     >               LINE, FREQIN, MODHEAD(1), JOBNUM(1),
     >               DISP, N, NBLINE, IDENT, OSMIN,
     >               XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >               ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >               .FALSE., FUNIT, NMOD, 
     >               LPHISTA, LPHIEND, NPHI,  
     >               XUNIT, .FALSE., BAIRWAVELENGTH, 
     >               RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
      
      IF (BCONT) THEN 
           WRITE(0,'(A,A)') 'Plotting Continuum flux'
           WRITE (KANAL1,'(A)') 'PLOT: CONTINUUM '//FREQIN
           WRITE (KANAL1,'(A)') '* CONTINUUM '//FREQIN
          CALL TRAPLO (1, CPROFILE, DLAM, NFOBS, LINE,
     >          FREQIN, MODHEAD(1), JOBNUM(1), DISP, N, NBLINE, 
     >          IDENT, OSMIN,
     >          XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >          ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >          BCONT, FUNIT, NMOD, 
     >          LPHISTA, LPHIEND, NPHI, 
     >          XUNIT, BCALIBRATED, BAIRWAVELENGTH,
     >          RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
C***    Plotten des kompletten, nichtrektifizierten spektrums
      IF (BCALIBRATED) THEN
          WRITE(0,'(A,A)') 'Plotting flux-calibrated spectrum'
          WRITE (KANAL1,'(A)') 'PLOT: CALIBRATED '//FREQIN
          WRITE (KANAL1,'(A)') '* CALIBRATED '//FREQIN
          CALL TRAPLO (1, PROFILEC, DLAM, NFOBS, LINE,
     >         FREQIN, MODHEAD(1), JOBNUM(1),
     >         DISP, N, NBLINE, 
     >         IDENT, OSMIN,
     >         XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >         ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >         .TRUE., FUNIT, NMOD, 
     >         LPHISTA, LPHIEND, NPHI, 
     >         XUNIT, BCALIBRATED, BAIRWAVELENGTH,
     >         RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
      

  501 CONTINUE
      IF (LSDWL .GT. 1)
     > CALL TRADWL (1, DELW, ADELW, VELO(1,1), ENTOT(1,1), ND,
     >              LINE, MODHEAD(1), JOBNUM(1),
     >              LEVEL(NUP), XLAM, XPLOT, TAURCONT, RADIUS,
     >              TRANSDWLLINE)
   10 CONTINUE
 
C***  PRINTOUT OF PROFILE TABLE ONLY WHEN NOT "DWL" MODE ACTIVE
      IF (LSDWL.LE.1) THEN
         CALL PRIPRO (XLAM, VDOP, NFOBS, PROFILE,
     >         XOBS0, DXOBS, JOBNUM(1),
     >         VSINI, MODHEAD(1), DLAM, LSPRO, IFIRST,
     >         NPHI, LPHISTA, LPHIEND,
     >         DXMAX, TAUMAX, XMAX, TAUMINBROAD, JPFIRST, JPLAST, 
     >         PGRID_MERGED,  
     >         PWEIGHT, ELEVEL, NDIM, INDNUP, INDLOW,
     >         LASTIND, INDLAP, MAXLAP, FREMAX, FREMIN,
     >         NBLINE, RMAX, XLAMLAP,
     >         WEIGHT, LEVEL, EINST, VMAX,
     >         IVERSION, LINE, REDIS,
     >         BWESEX, BROAD, LINPRO, AVOIGT, BNOCONT, NMOD, 
     >         DENSCON, FILLFAC, BFEWARNING, SPZ1,
     >         SPZ2,MANIPOP_OPTIONS,ND, 
     >         MULTIIND, NMULTI, BAIRWAVELENGTH, MAXSPZ, RCOROT,
     >         BDD_VDOP, DD_VDOP_LINE, BMICROTURB, EL_MIN, 
     >         IVDOPSTATUS, BVSINI_AT_RCOROT)
      ENDIF

      IF (BPLOTVDOP) THEN
        CALL PLOTVDOP (RADIUS, TAURCONT, VELO, DD_VDOP, 
     >                 ND, NATOM, BPLOTVDOP, SYMBOL, BMICROTURB,
     >                 DD_VMIC, DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD(1),
     >                 BIRONLINES)
      ENDIF

      IF (BLIMB) 
     >   CALL LIMBDARK_OUTPUT (KANAL1, LIMB_LINE, NFOBS, EMINT_P,
     >                            NP_MERGED, PGRID_MERGED)

C***  Entry in case of error which enforced skipping the  current range
    2 CONTINUE

C***  The no. of core-intersecting impact parameters might have been modified 
C***  in subr. ROTATION_PREP; therefore, P and Z were cloned and 
C***  are now restored for the potential next RANGE

      NP(1) = NP_ORIG 
      DO JP=1, NP_ORIG
         PGRID(JP,1) = PGRID_ORIG(JP)
      ENDDO

      ZGRID = ZGRID_ORIG

C***  In case of a second model, several vectors of both models have been
C***  saved before being rescaled to the merged radius grid 
C***  by subr. RESCALE_SECMOD, and must now be restored for 
C***  a possible further range
      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         VDU       = VDU_ORIG
         DD_VDOP   = DD_VDOP_ORIG
         DD_VDOPDU = DD_VDOPDU_ORIG
         DD_VMICDU = DD_VMICDU_ORIG
         T         = T_ORIG
         ENTOT     = ENTOT_ORIG
         RNE       = RNE_ORIG
      ENDIF


      GOTO 1
C***  ENDLOOP   (RANGE COMPLETED - continue with next RANGE) -----------

C***********************************************************************
      
 20   CLOSE (1)
 
      CALL CLOSMS (7, IERR)

C***  CLOSE THE CARDS-FILE  
      CLOSE (2)

C***  CLOSE THE FILE for plots of the wind-rotation geometry grid 
      CLOSE (65)

      CALL JSYMSET ('G0', '0')

      CALL STAMP (OPSYS, 'FORMAL', TIM1)

      STOP 'O.K.'
      END
 

       SUBROUTINE FORMCMF (TEFF, U, UK, NDDIM, NFLDIM,
     $        Z, OPA, ETA, ETANOTH, ETACK, ETACCK, OPACK, ETANCK, 
     $        XCMFRED, XCMFBLUE, DXCMF, XMAX, RMAX, VMAXDU, 
     $        THOMSON, THOMCK, OPAL, ETAL, 
     $        RADIUS, ND, NP, P,
     $        TA, TB, TC, UB, GA, H, QQ, S, V, VA, VB, PP,
     $        BCORE, DBDR, VDU, GRADI, VDOP,
     $        W0, NBLINE, DELXLAP, XRANGERED, XRANGEBLUE,
     $        XJNUE, OPAK, ETAK,
     $        XPLOT, YPLOT, XJC, NFDIM, XJCIND,
     $        LINE, MODHEAD, JOBNUM, NREDMAX, WREDI, 
     $        WREDI0, VDUEL, BWESEX,
     $        XLAMBDA, XLAM, T, RNE, POPNUM, ENTOT, RSTAR, MAINPRO, 
     $        MAINLEV, NOM, KODAT, NDIM, N, MAXATOM, LEVEL, NCHARG, 
     $        WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO, ADDCON1, 
     $        ADDCON2, ADDCON3, IGAUNT, SIGMATHK, SEXPOK, EDGEK, NF, 
     $        KONTNUP, KONTLOW, LASTKON, XDATA, REDIS, LBREF, INDREF, 
     $        NPDIM, A, B, C, W, BX, WX, EDDI, IWARNJ0, 
     $        IWARN, FNUEC, LSDWL, ST, BELIFI, 
     >        DENSCON, FILLFAC, ENTOTDENS, bBIGBANDLIMIT,
C *** IRON: Additional Parameter used in SUBROUTINE FECHECK called from FORMCMF
     >        INDRB, IFRBSTA, IFRBEND, LASTFE,
     >        CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >        INDFEACT, BFECHECK,
C *** IRON: Additional Parameter used in SUBROUTINE CMFFEOP called from FORMCMF
     >        SIGMAFE, OPAFE, ETAFE, 
     >        IFENUP, IFELOW, BIRONLINES,
     >        SIGMAACT, OPAFEI, ETAFEI,NFL,BAIRWAVELENGTHSET, 
     >        BAIRWAVELENGTH, VSIDU,
     >        DD_VDOPDU, NATOM, YSCRATCH, MAXLAP, IND_ELLINE, 
     >        bDDOPAFORMCMF, bDDFECONVOL)

C***********************************************************************
C***  CALLED FROM: FORMAL
C***  LINE RADIATION TRANSFER IN THE COMOVING FRAME, ACCOUNTING ITERATIVELY
C***     FOR PHOTON REDISTRIBUTION BY ELECTRON SCATTERING.
C***     THE RESULTING FREQUENCY-DEPENDING CONTINUUM EMISSIVITY IS STORED
C***     IN THE ARRAY ETACK(L,K). 
C***     THE ORIGINAL CONTINUUM EMISSIVITY IS STORED IN THE ARRAY
C***     ETACCK(L,K). THIS ARRAY IS USED IN PREPRAY TO CALCULATE THE RAY
C***     EMISSIVITY OF THE CONTINUUM
C***     THE CMF FREQUENCY OF INDEX K IS: XK = XCMFBLUE - (K-1) * DXCMF
C***     K RANGES FROM 1 TO NFL. L IS THE DEPTH INDEX.
C***
C***  The main part runs twice: 
C***     1) without lines (BWITHLINES = .FALSE.)
C***     2) with    lines (BWITHLINES = .TRUE.)
C***
C***  Output arrays:
C***     ETACCK(L,K) = continuum emissivity incl. redistribution, no lines
C***     ETACK (L,K) = continuum emissivity incl. redistributed line photons 
C***     OPACK (L,K) = continuum opacity    incl. Thomson term 
C***
C***  Only onternally used:
C***     ETAK  (L,K) = true emissivity: continuum, in 2) + sum of lines 
C***     ETANCK(L,K) = internal storage of emissivity without Thomson
C***********************************************************************

      INTEGER, INTENT(IN) :: ND, NATOM, NDDIM, NFLDIM, MAXLAP

      INTEGER, DIMENSION(MAXLAP) :: IND_ELLINE
      REAL, DIMENSION(ND) :: RADIUS, OPA, ETA, ETANOTH, THOMSON
      REAL, DIMENSION (NDDIM, NATOM), INTENT(IN) :: DD_VDOPDU
      DIMENSION ETACK(NDDIM,NFLDIM), OPACK(NDDIM,NFLDIM)
      DIMENSION ETACCK(NDDIM,NFLDIM)
      DIMENSION ETANCK(NDDIM,NFLDIM),THOMCK(NDDIM,NFLDIM)
      DIMENSION BCORE(NFLDIM), DBDR(NFLDIM)
      DIMENSION XLAMBDA(NFDIM), XJC(ND,NF)
      DIMENSION T(NDDIM), ENTOT(NDDIM), ENTOTDENS(NDDIM)
      DIMENSION W0(ND)
      DIMENSION U(ND,NP), UK(ND,NP), V(ND,NP), Z(ND,NP)
      DIMENSION PP(ND), VDU(ND), GRADI(ND)
      DIMENSION XJNUE(NDDIM,NFLDIM)
      DIMENSION OPAL(NDDIM,NBLINE), ETAL(NDDIM,NBLINE)
      DIMENSION DELXLAP(NBLINE)
      DIMENSION OPAK(NDDIM,NFLDIM), ETAK(NDDIM,NFLDIM)
      DIMENSION XPLOT(NFLDIM), YPLOT(NFLDIM), XJCIND(ND)
      DIMENSION INDFEACT(LASTFE) 
      REAL, DIMENSION(NFLDIM) :: YSCRATCH
C***  ARRAYS TO HANDLE DEPTH-DEPENDENT REDISTRIBUTION INTEGRAL
      DIMENSION VDUEL(NDDIM), WREDI0(NDDIM), WREDI(NREDMAX,NDDIM)
      REAL, DIMENSION(NDDIM, NFLDIM) :: OPAFE, ETAFE
      REAL, DIMENSION(NDDIM, NFLDIM) :: OPACFEK, ETACFEK
      INTEGER, DIMENSION(MAXATOM) :: KODAT

      CHARACTER*80  HEADER, XTEXT, YTEXT
      CHARACTER MODHEAD*100
      LOGICAL REDIS, bBIGBANDLIMIT, 
     >        bDDOPAFORMCMF, bDDFECONVOL
      LOGICAL BFECHECK, BFEWING, BIRONLINES, BWITHLINES
      LOGICAL BAIRWAVELENGTHSET, BAIRWAVELENGTH
      
      INTEGER :: NAFE

      DIMENSION QQ(ND), S(ND)

      DIMENSION DENSCON(ND),FILLFAC(ND)

C***  1 / SQRT(PI)
      DATA WPIINV / 0.564189583549 /

C***  TEST PLOT FACILITY: J-NUE AT DEPTH POINT LJPLOT (0 = DISABLED)
      LJPLOT = 0

C-----------------------------------------------------------------
C***  MAXIMUM NUMBER OF ITERATION
      MAXITER = 100
C-----------------------------------------------------------------
C***  EPSILON CONVERGENCE-CRITERION
      EPS = 0.001
C-----------------------------------------------------------------
C***  MAXIMUM HALF-BANDWIDTH FOR CMF LINE TRANSFER IN DOPPLER UNITS
      CMFBAND = 4.5
C-----------------------------------------------------------------
C***  SPACING OF CMF FREQUENCY POINTS
      DXCMF = 0.3
C-----------------------------------------------------------------
C***  HALF BANDWIDTH OF ELECTRON REDISTRIBUTION INTEGRAL
C***     IN UNITS OF THE ELECTRON DOPPLER VELOCITY
      BWES = 2. * BWESEX
C-----------------------------------------------------------------
C***  THE FREQUENCY RANGES (BOTH, CMF AND OBS. FRAME) ARE EXTENDED
C***      BY ONE HALF BANDWIDTH, MULTIPLIED WITH THE FOLLOWING FACTORS:
      BWEXRED = 2.
      BWEXBLU = 1.5
C-----------------------------------------------------------------
 
C***  Switch for adding lines (incl. iron bands) - else only continuum 
      BWITHLINES = .FALSE.
      WRITE (0,'(A)') "FORMCMF: Continuum Thomson redistribution"

C***  CALCULATE THE RANDOM VELOCITY OF THE ELECTRONS IN DOPPLER UNITS
C***       ATOMIC MASS OF HELIUM
           ATMASS = 4.
      DO 27 L=1, ND
        VTHHE2 = 0.01651 * T(L) / ATMASS
        VMICRO2 = VDOP * VDOP - VTHHE2
        VTHEL2 = 30.3165 * T(L)
        VDUEL(L) = SQRT (VMICRO2 + VTHEL2) / VDOP
   27 CONTINUE
C----------------------------------------------------------------------
C***  Prepare Clump density
      DO L=1, ND
         ENTOTDENS(L) = ENTOT(L) * DENSCON(L)
      ENDDO

      ALN = ALOG (1. - VDOP / CLIGHT)
      XLAMLN = ALOG (XLAM)

C***  DEFINE CMF FREQUENCY RANGE
C***  Since for the continuum, electron redistribution is now calculated
C***  in any case,, the range is always extended independent from the
C***  REDIS option - wrh  3-Jul-2020  
C***  The extension for electron reditribution is also returned to
C***  the calling FORMAL
cc      IF (REDIS) THEN
        FREMINEX = BWEXRED * BWES * VDUEL(ND)
        FREMAXEX = BWEXBLU * BWES * VDUEL(ND)
cc      ELSE
cc        FREMINEX = .0
cc        FREMAXEX = .0
cc      ENDIF

C***  The Band (XCMFRED,XBLUE) must cover the whole range needed in 
C***      OBSFRAM --> ZONEINT
C***  It is not clear to me why *twice* VDU(1) is needed (wrh 27-May-2008)
      BIGGERBAND =  AMAX1(XMAX,CMFBAND)
C***  Range should not grow by more than 0.2c due to XMAX (Tomer, 8.1.2015)
C***  This range cutting can be suppressed by the CARDS option NO-BIGBANDCUT     
      IF (bBIGBANDLIMIT) THEN
        BIGGERBAND = AMIN1(BIGGERBAND, 0.2 * CLIGHT / VDOP)
      ENDIF  
      XCMFBLUE = DELXLAP(NBLINE) + BIGGERBAND + FREMAXEX + 2.1*VMAXDU
     >                                                   + VSIDU
      XCMFRED  = DELXLAP(1)      - BIGGERBAND - FREMINEX - 2.1*VMAXDU
     >                                                   - VSIDU

C***  The range may be extended by the RANGE option (but not shrinked)
      IF (XRANGERED .NE. XRANGEBLUE) THEN
         XCMFBLUE = AMAX1 (XCMFBLUE, XRANGEBLUE+ BIGGERBAND + 
     >                              FREMAXEX + 1.1*VMAXDU) 
         XCMFRED  = AMIN1 (XCMFRED , XRANGERED - BIGGERBAND - 
     >                              FREMAXEX - 1.1*VMAXDU) 
         XLAMBDAC = XLAM * ( EXP( XRANGEBLUE * ALN ) 
     >               +        EXP( XRANGERED  * ALN ) ) / 2. 
      ELSE
         XLAMBDAC = XLAM * ( EXP( XCMFBLUE * ALN ) 
     >               +        EXP( XCMFRED  * ALN ) ) / 2. 
      ENDIF

C***  Calculate Iron OPAs/ETAs for Air, lower limit = 3200 Ang
      BAIRWAVELENGTH = ( (.NOT. BAIRWAVELENGTHSET) .AND. 
     >      ( XLAMBDAC .GE. 3200. .AND. XLAMBDAC .LE. 10000. ) )
     >  .OR. ( BAIRWAVELENGTHSET .AND. BAIRWAVELENGTH ) 

      NFL = 1 + (XCMFBLUE - XCMFRED) / DXCMF

      IF (NFL .GT. NFLDIM) THEN
         CALL REMARK ('NFLDIM INSUFFICIENT')
         WRITE (0,*) 'CMFBLUE, RED: ', XCMFBLUE, XCMFRED
         WRITE (0,*) 'XMAX, TOMERTHRESH: ', XMAX, 0.2 * CLIGHT / VDOP
   99    FORMAT (' NFLDIM GIVEN:', I8, ', BUT REQIERED:', I8)
         WRITE (0,99) NFLDIM, NFL
         PRINT 99, NFLDIM, NFL
         STOP 'ERROR'
         ENDIF

C*********************************************************
C***  DEFINE WEIGHTS FOR REDISTRIBUTION INTEGRAL
C***     THE REDISTRIBUTION FUNCTION IS THE ANGLE-AVERAGED
C***      REDISTRIBUTION FUNCTION FOR ELECTRON SCATTERING, R(e,A)
C***      (cf. MIHALAS: STELLAR ATMOSPHERES, P. 432)
C***      R(e,a) IS CODED AS FUNCTION FIERFC
C*********************************************************
C***      NUMBER OF INTEGRATION STEPS TO ONE SIDE : NREDI
      NREDI = NINT (BWES * 2. * VDUEL(ND) / DXCMF)
      IF (NREDI .GT. NREDMAX) THEN
        CALL REMARK ('DIMENSION NREDMAX INSUFFICIENT')
        WRITE (0,*) 'NREDI=',NREDI, '  NREDMAX=',NREDMAX
        STOP 'ERROR'
      ENDIF
      DO 28 L=1, ND
        WREDI0(L) = FIERFC(0.)
        SUM = WREDI0(L)
        DO 40 I=1,NREDI
          X = I * DXCMF/ (2. * VDUEL(L))
          WREDI(I,L) = FIERFC(X)
          SUM = SUM + 2. * WREDI(I,L)
   40   CONTINUE
C***  NORMALISATION
        WREDI0(L) = WREDI0(L) / SUM
        DO 41 I=1,NREDI
   41     WREDI(I,L) = WREDI(I,L) / SUM
   28 CONTINUE
 
C***  VERSION OF OUTER BOUNDARY CONDITION : NO INCIDENT RADIATION
      XIMINUS=0.
C***  FREQUENCY DERIVATIVE OF IMINUS
      DXI=0.

C***  IRON: WIDTH OF RED LINE WING IN IRON-DELTAX-UNITS
      IF (BIRONLINES) THEN
         DFEDUMMY= 1.
      ENDIF
  
C***  IRON: SET LOWEST LINE TO '1', NO ACTIVE FE-LINES 
      INDFEACT(1) = 1
      MAXFEACT    = 0
      BFECHECK    = .FALSE.
      BFEWING     = .FALSE.
      VDOPMAX = MAXVAL(DD_VDOPDU)

C***  CALCULATION OF OPAK, ETAK = ALL BLENDING LINES
      WRITE (0,*) 'Calculating total opacities...'
      DO 12 K=1,NFL
        XK = XCMFBLUE - (K-1) * DXCMF
        XLAMREF = EXP (XLAMLN+ XK * ALN)
        DO 15 L=1,ND
          OPAK(L,K) = 0.
          ETAK(L,K) = 0.
   15   CONTINUE
C***  LOOP EINSCHRAENKEN?
        DO 14 IBLEND=1,NBLINE
          XKBL = XK - DELXLAP(IBLEND)
C***      Worst Doppler width appoximation          
          XKLBMAXDU = XKBL*XKBL / VDOPMAX/VDOPMAX
C***      Do not evaluate PHI for more than 10 Dopplerwidths          
          IF (XKLBMAXDU > 4.5) CYCLE
          IF (bDDOPAFORMCMF) THEN
            NA = IND_ELLINE(IBLEND)          
          ELSE
            PHIX = EXP(-XKBL*XKBL) * WPIINV
          ENDIF
          DO 13 L=1,ND
            IF (bDDOPAFORMCMF) THEN
              XKBLL = XKBL / DD_VDOPDU(L,NA)
              PHIX = EXP(-XKBLL*XKBLL) * WPIINV / DD_VDOPDU(L,NA)
            ENDIF
            OPAK(L,K) = OPAK(L,K) + OPAL(L,IBLEND) * PHIX
            ETAK(L,K) = ETAK(L,K) + ETAL(L,IBLEND) * PHIX
              
   13     CONTINUE
   14   CONTINUE

       
        IF (BIRONLINES) THEN
C***       IRON: CHECK FOR ACTIVE BOUND-BOUND TRANSITIONS OF GENERIC ION
           CALL FECHECK (XLAMREF, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEDUMMY)

C***       IRON: NON-LTE IRON OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
           CALL CMFFEOP (XLAMREF, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE(1,K), ETAFE(1,K), INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT,
     >                    SIGMAACT, OPAFEI, ETAFEI, T, 1, TEFF) 

        ENDIF

C***    Inner Boundary values at each CMF frequency (diffusion approx.):
C***    BCORE = B_nue    DBDR = dB/dr  
C***    BCOREL and DBDRL re stored in the vectors BCORE and DBDR. 
C***    The local variables are used later in CALL BACKJCU
        CALL DIFFUS (XLAMREF, T, RADIUS, ND,
     >               BCOREL, DBDRL, DUMMY, DUMMY, .TRUE.)
        BCORE(K) = BCOREL
        DBDR(K) = DBDRL
 
        CALL COOP (XLAMREF,ND,T,RNE,POPNUM,ENTOTDENS,RSTAR,
     $             OPA,ETANOTH,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $             NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $             ALPHA,SEXPO,
     $             ADDCON1, ADDCON2, ADDCON3,
     $             IGAUNT,SIGMATHK,SEXPOK,EDGEK,0,NF,DUMMY,
     $             RADIUS,KONTNUP,KONTLOW,LASTKON,XDATA)
C***    Opacity is scaled down from clump to average value
        DO L=1, ND
          OPA(L) = OPA(L) * FILLFAC(L)
          ETANOTH(L) = ETANOTH(L) * FILLFAC(L)
        ENDDO

C***    U IS NEEDED AS BLUE WING BOUNDARY FOR THE RADIATION TRANSFER;
C***    IN CASE OF DWL PLOT (LSDWL > 0) FNUEC IS ALSO NEEDED
        IF (K .EQ. 1 .OR. LSDWL .GT. 0) THEN
C+++       abarnisk 06/2004: "ELIMIN" ersetzt hier "BACKJCU"-Aufruf 
C***       Note: EDDI is output parameter, only needs dimensioned space
           CALL ELIMIN (XLAMREF,FNUEC,DUMMY2,U,Z,A,B,C,W,BX,WX,XJCIND,
     $          RADIUS,P,BCORE,DBDR,OPA,ETANOTH,THOMSON,EDDI,ND,NP,NPDIM,
     $          ENTOT,INDREF, IDUMMY, ST, BELIFI, 0)
           IF (LSDWL .GT. 0) GOTO 31
        ELSE
            CALL BACKJC (XJC, ND, NF, XJCIND, XLAMBDA, XLAMREF, RADIUS)
        ENDIF

C***    FILL FREQUENCY-INDEXED ARRAYS 
        DO 17 L=1,ND
          ETANCK(L,K) = ETANOTH(L)
          OPACK (L,K) = OPA(L)
          THOMCK(L,K) = THOMSON(L)
C***      Initialize ETACK with coherent scattering
          ETACK(L,K) = ETANOTH(L) + OPA(L) * THOMSON(L) * XJCIND(L)
   17   CONTINUE

   12 CONTINUE

      IF (BIRONLINES .AND. bDDFECONVOL) THEN
C***      Convolve iron line opacities and emissivities with
C***       the correct DD_VDOPDU (if larger than VDOPFE)
        WRITE (0,*) 'Convolving iron line opacities with specified'
     >              // ' VDOP stratification...'
        NAFE = KODAT(26)
        CALL CONVOLOPAFE(OPAFE, YSCRATCH, NDDIM, ND, NFL, DXCMF,
     >                   VDOP, VDOPFE, DD_VDOPDU(1,NAFE))
        CALL CONVOLOPAFE(ETAFE, YSCRATCH, NDDIM, ND, NFL, DXCMF,
     >                   VDOP, VDOPFE, DD_VDOPDU(1,NAFE))
      ENDIF   
      WRITE(0,'(A)') 'STARTING FORMCMF WITH CONTINUUM ONLY'
   
C*** Originally, NOREDIS would cause to exit the routine here, adding 
C*** only coherent Thomson term to ETA
C*** Now ETA and J are calculated iteratively, with the continuum Thomson 
C*** term redistributed
C*** Tomer, 27/5/14
     
 301  CONTINUE 
C*** ITERATION LOOP =================================================
      DO 30 ITER=1, MAXITER

C***  SET XJNUE = 0.
      DO 16 L=1,ND
      DO 16 K=1,NFL
   16   XJNUE(L,K) = 0.


C***  LOOP OVER IMPACT PARAMETERS
      DO 11 JP=1,NP-1

      LMAX=MIN0(NP+1-JP,ND)
      LZ=LMAX-1

C***  *****************************************************************
C***  The following formalism is a difference scheme for the ray-by-ray 
C***  radiative transfer, as descibed in my Thesis (Hamann 1979). 
C***  U and V vectors are the intensity-like and flux-like Feautrier
C***  variables
C***  - In program COLI we have meanwhile replaced this method by a
C***    short-characteristic integration, which guarentees positive
C***    intensities. One might consider to apply shortchar here, too
C***    (wrh, 16-Apr-2013)
C***  *****************************************************************

C***  INITIALIZE UK = FEAUTRIER U AT BLUE WING
      DO 6 L=1,LMAX
    6   UK(L,JP) = U(L,JP)

C***  INITIALIZE FEAUTRIER V AS DERIVATIVE OF U
       DO 5 L=1,LZ
        X=0.5*(OPACK(L,1)+OPACK(L+1,1))
        V(L,JP)=(UK(L+1,JP)-UK(L,JP)) / (X*(Z(L,JP)-Z(L+1,JP)))
    5  CONTINUE

C***  GENERATE WEIGHTS FOR INTEGRATION OF THE 0.-MOMENT Jnue
      CALL GENW0 (JP, LMAX, ND, NP, Z, RADIUS, P, W0)

C***  PP(L) = VELOCITY GRADIENT, PROJECTED ON THE PRESENT RAY J, /DELTAX
C***  GDU = GRADI, CONVERTED IN DOPPLER UNITS
C***  VDU = VELO IN DOPPLER UNITS, WAS CONVERTED ALREADY IN FORMAL 
      DO 4 L=1,LMAX
      RL = RADIUS(L)
      Y  = Z(L,JP) / RL
      YY = Y*Y
      GDU = GRADI(L) / VDOP
    4 PP(L) = (YY * GDU + (1.-YY) * VDU(L)/RL) / DXCMF


C***  LOOP FOR ALL (ORIGINAL) FREQUENCY POINTS
C***  THE CURRENT INDEX K REFERS TO THE RESULTING U,V , WHERE ALL QUANTITIES
C***  ARE TAKEN. THE FREQUENCY DERIVATIVE STRETCHES BACK TO K-1

      DO 1 K=1,NFL
         
C+++ abarnisk 06/2004: Verbesserung des Startwertes fuer U (bzw. UK) durch
C+++ mehrfaches "auf der Stelle iterieren" von "CMFSET_FORMAL" mit OPACK(1,1)
C+++ ETACK(1,1). Laenge der Schleife willkuerlich gewaehlt!
         
         IF (K.EQ.1) THEN 
            DO 311 I=1,100      ! oder 50,20,10???
               CALL CMFSET_FORMAL(Z(1,JP),
     $           ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,
     $           OPACK(1,K),ETACK(1,K),PP,BCORE(K),DBDR(K),
     $           RADIUS,XIMINUS,DXI,OPAK(1,K),ETAK(1,K),
     >           OPAFE(1,K),ETAFE(1,K),BWITHLINES)
               CALL VMALV (VA,VB,V(1,JP),QQ,LMAX)
               DO 501 L=1, LMAX
 501              QQ(L) = QQ(L) + S(L)
                  CALL MDV (UB,UK(1,JP),LMAX)
                  CALL VADD (UK(1,JP),QQ,LMAX)
                  CALL INVTRI (TA,TB,TC,UK(1,JP),LMAX)
                  CALL MDV (H,V(1,JP),LZ)
                  CALL GMALU (GA,UK(1,JP),S,LMAX)
                  CALL VADD (V(1,JP),S,LZ)
 311           CONTINUE
               GOTO 2
         ENDIF
            
            CALL CMFSET_FORMAL(Z(1,JP),
     $           ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,
     $           OPACK(1,K),ETACK(1,K),PP,BCORE(K),DBDR(K),
     $           RADIUS,XIMINUS,DXI,OPAK(1,K),ETAK(1,K),
     >           OPAFE(1,K),ETAFE(1,K),BWITHLINES)
            
      CALL VMALV (VA,VB,V(1,JP),QQ,LMAX)

C***  THE FOLLOWING INLINING SAVES CPU TIME:
C***   (NOTE: IT SEEMS THAT INLINING OF THE OTHER CALLS WITH MATRIX COLUMNS
C***          AS FORMAL PARAMETERS DOES NOT SAVE TIME.)
C!!      CALL VADD (QQ,S,LMAX)
      DO 50 L=1, LMAX
   50 QQ(L) = QQ(L) + S(L)

      CALL MDV (UB,UK(1,JP),LMAX)

      CALL VADD (UK(1,JP),QQ,LMAX)

      CALL INVTRI (TA,TB,TC,UK(1,JP),LMAX)
C***  NOW U IS THE FIELD AT THE NEW INDEX K
      CALL MDV (H,V(1,JP),LZ)
      CALL GMALU (GA,UK(1,JP),S,LMAX)

      CALL VADD (V(1,JP),S,LZ)
 
2     CONTINUE
C***  end of the Feautrier difference scheme ************************

C***  ADDING THE NEW UK(L) TO THE ANGEL-AVERAGED INTENSITY XJNUE
      DO 20 L=1,LMAX
20    XJNUE(L,K) = XJNUE(L,K) + UK(L,JP) * W0(L)

    1 CONTINUE
   11 CONTINUE

C***  UPDATING THE FREQUENZ DEPENDENT CONTINUUM EMISSIVITY
C***  ANGLE AVERAGED FREQUENCY-REDISTRIBUTION
C***  CORMAX = MAX. RELATIVE CORRECTION OF ETA
      CORMAX = 0.
      DO 18 L=1,ND
        DO 18 K=1,NFL
C***  INITIALIZING OF REDISTRIBUTION INTEGRAL
          OPATH = OPACK(L,K) * THOMCK(L,K)
          ETACKNEW = ETANCK(L,K) + OPATH * XJNUE(L,K) * WREDI0(L)
          DO 19 I=1,NREDI
            KPI = K+I
            IF (KPI .GT. NFL) KPI = NFL
            KMI = K-I
            IF (KMI .LT.1 ) KMI = 1
            ETACKNEW = ETACKNEW + OPATH * (XJNUE(L,KPI) +
     >                 XJNUE(L,KMI)) * WREDI(I,L)
   19     CONTINUE

            RELCOR = ABS (ETACKNEW / ETACK(L,K) - 1.) 
            IF (RELCOR .GT. CORMAX) CORMAX = RELCOR
          ETACK(L,K) = ETACKNEW
   18 CONTINUE

C***  TEST PLOT FACILITY ********************************************
      IF (LJPLOT .GT. 0) THEN

C***  DEFINE THE Y-VECTOR
      DO K=1, NFL
         YPLOT(K) = XJNUE(LJPLOT,K) 
      ENDDO

      IF (ITER .EQ. 1) THEN
C***  DEFINE X-VECTOR: Wavelength
      DO  K=1, NFL
        XK = XCMFBLUE - (K-1) * DXCMF
          XPLOT(K) = EXP (XLAMLN+ XK * ALN)
      ENDDO
      XTEXT   = '\CENTER\#l# [\A]'

C***  Label Y-AXIS
      YTEXT = 'INTENSITY J&T#n#&M'
      WRITE (HEADER, 10) BWITHLINES, LJPLOT, MODHEAD(13:33)
   10 FORMAT ('J&T#n#&M:  WITHLINES=', L1, '  L=', I2,
     >        '   MODEL ', A )

        CALL PLOTANFS (1, HEADER, HEADER, XTEXT, YTEXT,
     $             0., 0., 0., 0., 0., 0.,
     $             0., 0., 0., 0., 0., 0.,
     $             XPLOT, YPLOT, NFL, 'COLOR=2')
      ELSE
        ISYMBOL = 5
        IF (ITER .EQ. 2) ISYMBOL = 9
        IF (ITER .EQ. 3) ISYMBOL = 10
        CALL PLOTCON (1,XPLOT,YPLOT,NFL,ISYMBOL)
      ENDIF
      ENDIF
C********************************************************************

      WRITE (0,'(A,I3,2X,A,F11.5)') 
     >      'FORMCMF : Iteration-Nr.:', ITER, 'Cormax=', CORMAX

C***  Check for convergence 
      IF (CORMAX .LT. EPS .AND. ITER .GT. 1) THEN
         WRITE(0,'(A,F6.4)')'Converged: Cormax <',  EPS
         IF (.NOT.BWITHLINES) THEN 
C***        This was the first run without lines
C***        -> save the continuum emissivity (whole array)
         DO L=1, ND
            DO K=1, NFL
               ETACCK(L,K) = ETACK(L,K)
            ENDDO
         ENDDO
C***     If NOREDIS=TRUE, then skip the Thomson redistribution of lines 
C***     and exit routine
            IF (.NOT. REDIS) THEN 
                WRITE(0,'(A)') 'Line redistribution due to electron '
     >                         // ' scattering skipped! (NOREDIS)'
                GOTO 31
            ENDIF
            BWITHLINES = .TRUE.
            WRITE(0,'(A)')'NOW FORMCMF WITH LINES'
            GOTO 301
         ELSE
C***     This was the second run **with** lines
            GOTO 31
         ENDIF
      ENDIF

C***  Check for divvergence 
      IF (CORMAX .GT. 1.E6 .AND. ITER .GT. 1 .AND. BWITHLINES) THEN
         WRITE (0,'(A)')  
     >   '*** WARNING: Iteration of electron-scattering redistribution',
     >   ' terminated early because of apparent divergence'
         WRITE (*,'(2A)')  
     >   '*** WARNING: Iteration of electron-scattering redistribution',
     >   ' terminated early because of apparent divergence'
          GOTO 31

         STOP '*** FATAL ERROR in FORMCMF'
      ENDIF

C***  END OF ITERRATION LOOP ===================================
30    CONTINUE
      WRITE (*,'(A,I3,2A)')  'MAX. NUMBER (',MAXITER,') OF ITERATIONS ',
     >                       'IN FORMCMF INSUFFICIENT'
      WRITE (0,'(A,I3,2A)')  'MAX. NUMBER (',MAXITER,') OF ITERATIONS ',
     >                       'IN FORMCMF INSUFFICIENT'

31    CONTINUE

      RETURN
      END
 
 





      SUBROUTINE FORMOSA (ND, R, NP, P, Z, ENTOT, RNE, ABXYZ, NATOM, 
     >                    T, VELO, NF, XLAMBDA, GRADI,
     >                    POPNUM, RSTAR, VDOP, VMIC,
     >                    JOBNUM, N, NDDIM, NPDIM, NFDIM,
     >                    MODHEAD, TEFF, MAXXDAT, XDATA, XJC, IMOD, 
     >                    DENSCON, FILLFAC, TAURCONT, POPMIN, 
     >                    ZERO_RATES, RCON, NDIM, XMDOT)
C***********************************************************************
C***  CALLED FROM FORMAL, READS THE MODEL FILE
C***********************************************************************
 
      DIMENSION XDATA(MAXXDAT), ABXYZ(NATOM)
      DIMENSION XJC(2),EDDI(2),xlambda(2)
      CHARACTER*8 NAME, MODHEAD*(*)
      DIMENSION DENSCON(ND), FILLFAC(ND), TAURCONT(ND), VMIC(ND)
      LOGICAL, DIMENSION(NDIM*NDDIM) :: ZERO_RATES
      REAL :: RCON

C***  Local arrays:      
      INTEGER, PARAMETER :: NDMAX = 200
      REAL, DIMENSION(NDMAX) :: VTURB
           
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
           
      M = IMOD - 1

      CALL OPENMS(3+M, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS(3+M, ND,1,          'ND      ', IERR)
      IF (ND .GT. NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      IF (ND > NDMAX) THEN
          WRITE (hCPR,'(A)') 'FORMOSA: FATAL ERROR ******'
          WRITE (hCPR,'(A)') 'FORMOSA: LOCAL DIMENSION FOR VTURB INSUFFICIENT'
          STOP 'ERROR IN FORMOSA'
      ENDIF
      CALL READMS (3+M, R,ND,         'R       ', IERR)
      CALL READMS (3+M, NP,1,         'NP      ', IERR)
      IF (NP.GT.NPDIM) THEN
            CALL REMARK ('TOO MANY IMPACT-PARAMETER POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3+M, P,NP,         'P       ', IERR)
      CALL READMS (3+M, Z,ND*NP,      'Z       ', IERR)
      CALL READMS(3+M, ENTOT,ND,      'ENTOT   ', IERR)
      CALL READMS(3+M, TAURCONT,ND,   'TAURCONT', IERR)
      CALL READMS(3+M, POPMIN,1,      'POPMIN  ', IERR)
      DENSCON(1:ND) = 0.0
      CALL READMS(3+M,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               WRITE(0,*)'Error in depth-dep. clumping during FORMAL!'
               STOP 'Error in depth-dep. clumping during FORMALLARGECL!'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
      CALL READMS (3+M, RNE,ND,       'RNE     ', IERR)
      CALL READMS (3+M, T,ND,         'T       ', IERR)
      CALL READMS (3+M, VELO,ND,      'VELO    ', IERR)
      CALL READMS (3+M, GRADI,ND,     'GRADI   ', IERR)
      CALL READMS (3+M, NF,1,         'NF      ', IERR)
      IF (NF.GT.NFDIM) THEN
         CALL REMARK ('TOO MANY FREQUENCY POINTS')
         STOP 'ERROR'
      ENDIF
      CALL READMS (3+M, XLAMBDA,   NF,'XLAMBDA ', IERR)
      CALL READMS (3+M, POPNUM,ND*N,  'POPNUM  ', IERR)
      CALL READMS (3+M, RSTAR,1,      'RSTAR   ', IERR)
      CALL READMS (3+M, VDOP,1,       'VDOP    ', IERR)
      CALL READMS (3+M, JOBNUM,1,     'JOBNUM  ', IERR)
      CALL READMS (3+M, MODHEAD,13,   'MODHEAD ', IERR)
      CALL READMS (3+M, TEFF,1,       'TEFF    ', IERR)
      CALL READMS (3+M, XMDOT,1,      'XMDOT   ', IERR)

C***  Read microturbulence vector from main iteration      
      CALL READMS (3,VMIC, ND,  'VMIC    ', IERR)
      IF (IERR == -10) THEN
C*      old MODEL file => only one VTURB value 
        CALL READMS(3,VTURB(ND),1, 'VTURB   ', IERR)
C*      very old MODEL file => neither VMIC nor VTURB        
        IF (IERR == -10) VTURB(ND) = -99.
        DO L=1, ND
C*        convert VTURB values to VMIC 
          IF (L /= ND) VTURB(L) = VTURB(ND)
          VMIC(L) = VTURB(L) * SQRT(2.)
        ENDDO
      ELSE
C*      Indicate non-existing VMIC      
        DO L=1, ND
          VMIC(L) = -99.
        ENDDO
      ENDIF
      
C***  Flags for the POPMIN levels
      CALL READMS (3+M,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
C*    Default if variable does not exist yet
      IF (IERR .EQ. -10) THEN
        DO I=1, N*ND
          ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF

      CALL READMS(3+M, RCON,1,        'RCON    ', IERR)
      
C***  'XDATA' IN FORMAL NOT NECESSARY
C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3+M, XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK ('ERROR WHEN READING XDATA FROM MODEL FILE')
C***     XFILL EQ 0.
         XDATA(1) = 0.  
      ENDIF

C***  READ ALL CONTINUUM INTENSITIES
      ND3=3*ND
      DO 15 K=1,NF
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL READMS(3+M, XJC(1+ND*(K-1)),ND,NAME, IERR)
   15 CONTINUE
 
C***  ABXYZ only needed for calling LTEPOP in MANIPOP, and for PRI_PAR
C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3+M, ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF

C***  NOT EXISTING RECORD 'ABXYZ': DEFAULT IS AN ATOMIC DATA FILE "DATOM"
C***  CONTAINING "HELIUM" AS THE ONLY ELEMENT
      IF (IERR .EQ. -10) THEN
         IF (NATOM .EQ. 1) THEN
            ABXYZ(1)=1.
         ELSE
            CALL REMARK ('NOT EXISTING RECORD ABXYZ')
            STOP 'ERROR'
         ENDIF
      ENDIF
 
      CALL CLOSMS(3+M,  IERR)

      RETURN
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
      SUBROUTINE GENW0 (JP,LMAX,ND,NP,Z,RADIUS,P,W0)
C***********************************************************************
C***  AT GIVEN IMPACT-PARAMETER INDEX JP, THIS SUBROUTINE CALCULATES THE
C***  WEIGHTS FOR THE ANGLE INTEGRATION OF THE 0. MOMENT, W0(L),
C***  AT ALL DEPTH POINTS L.
C***  THE ANGLE INTEGRATION IS PERFORMED IN THE VARIABLE "Z".
C***   - CALLED FROM: SUBROUTINE CMFRAY
C***********************************************************************
 
      DIMENSION W0(ND), Z(ND,NP),P(NP),RADIUS(ND)
 
C***  LOOP OVER ALL DEPTH POINTS
      DO 10 L=1,LMAX

      RL=RADIUS(L)
      RL2=RL+RL
      IF (JP.GT.1) GOTO 8

C***  FIRST STEP IF JP=1
      B=Z(L,1)
      A=Z(L,2)
      W0(L)=(B-A)/RL2
      GOTO 10

    8 IF (L.EQ.LMAX .AND. JP.GT.(NP-ND) ) GOTO 9

C***  INTERMEDIATE STEP
      A=Z(L,JP+1)
      B=Z(L,JP)
      C=Z(L,JP-1)
      W0(L)=(C-A)/RL2
      GOTO 10

C***  LAST STEP, IMPLYING Z(L,JMAX)=0
    9 B=Z(L,JP-1)
      W0(L)=B/RL2

   10 CONTINUE

      RETURN
      END
      SUBROUTINE GMALU (GA,U,V,LMAX)
C***********************************************************************
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
C***********************************************************************
      DIMENSION GA (LMAX),U(LMAX),V(LMAX)

      LZ=LMAX-1
      DO 1 L=1,LZ
    1 V(L)=GA(L)*(U(L)-U(L+1))

      RETURN
      END
      SUBROUTINE HORNER (X, Y, KPLUS1, Z)
C***  Polynomial Y(X) 
C***  The first coefficient belongs to the highest exponent!
      DIMENSION Z(KPLUS1)

      Y = 0.

      DO I=1, KPLUS1-1
        Y = (Y + Z(I)) * X
      ENDDO

      Y = Y + Z(KPLUS1)

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
      SUBROUTINE INSERT_LINE 
     >                 (LINE, INDBL, NBLINE, INDLAP, XLAMLAP, DELXLAP,
     $                   XLAMBL, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C*******************************************************************************
C***  This subroutine inserts an additional line in the list
C***  such that it stays sorted by ascending DELTAX
C***  Various line parameters are inserted accordingly in
C***  their respective arrays:
C***     INDLAP, XLAMLAP, LINPRO, AVOIGT 
C***  CALLED FROM FORMAL - PREFORM  
C*******************************************************************************

      IMPLICIT NONE

      INTEGER I, LOCDXGT, NDDIM, MAXLAP, MAXMOD, L, INDBL, NBL, LINE
      INTEGER ISRCHFGT, IMOD, NMOD, NBLINE, NPAR
      REAL ALN, AVOIGTBL, XLAMBL, XLAM, DELTAX 

      INTEGER, DIMENSION (MAXMOD) :: ND
      INTEGER, DIMENSION (MAXLAP) :: INDLAP
      REAL,    DIMENSION (MAXLAP) :: XLAMLAP, DELXLAP
      REAL, DIMENSION (MAXLAP,NDDIM,MAXMOD) :: AVOIGT
      CHARACTER KARTE*80
      CHARACTER*8 LINPRO(MAXLAP), PROF, PROFDEFAULT, LINPROBL
      INTEGER, PARAMETER :: NPARMAX = 4 
      CHARACTER*10 ACTPAR(NPARMAX)

C***  If this is the very first line: initialize the vectors
      IF (NBLINE .EQ. 0) THEN
         DELTAX = .0
         LOCDXGT = 1
         LINE = INDBL
      ELSE
         DELTAX = ALOG(XLAMBL/XLAM) / ALN
         LOCDXGT=ISRCHFGT(NBLINE,DELXLAP,1,DELTAX)
         IF (LOCDXGT .LE. NBLINE) THEN
            CALL SHIFT (INDLAP,LOCDXGT,NBLINE)
            CALL SHIFT (XLAMLAP,LOCDXGT,NBLINE)
            CALL SHIFT (DELXLAP,LOCDXGT,NBLINE)
            CALL SHIFT (LINPRO,LOCDXGT,NBLINE)
            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  CALL SHIFT (AVOIGT(1,L,IMOD),LOCDXGT,NBLINE)
               ENDDO
            ENDDO
          ENDIF
      ENDIF

      NBLINE=NBLINE+1

      IF (NBLINE .GE. MAXLAP) THEN
               PRINT *,'*** ERROR: MAX. NUMBER OF BLEND'
     $                ,' COMPONENTS EXCEEDED: MAXLAP=', MAXLAP
               WRITE (0,*) '*** ERROR: MAX. NUMBER OF BLEND'
     >                     ,' COMPONENTS EXCEEDED: MAXLAP=', MAXLAP
               STOP '*** FATAL ERROR detected by subd. INSERT_LINE' 
      ENDIF

      INDLAP(LOCDXGT)=INDBL
      XLAMLAP(LOCDXGT)=XLAMBL
      DELXLAP(LOCDXGT)=DELTAX
      LINPRO(LOCDXGT)=LINPROBL
      DO IMOD=1, NMOD
         AVOIGT(LOCDXGT,1,IMOD)=AVOIGTBL
      ENDDO

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
      SUBROUTINE INV (N,NDIM,A,CKEY)
C**********************************************************************
C***  MATRIX INVERSION BY CALLING THE CFT-LIBRARY SUBR. MINV
C***  A      = MATRIX
C***  N      = RANK OF SUBMATRIX (UPPER-LEFT CORNER) TO BE INVERTED
C***  NDIM   = ROW DIMENSION OF TOTAL MATRIX
C***  NMAX   = MAXIMUM VALUE OF NDIM
C**********************************************************************

c      PARAMETER (NMAX=2000)
      PARAMETER (NMAX=4000)
      DIMENSION A(NDIM,NDIM)
      CHARACTER*4 CKEY

C***  Dimensum must be greater-equal 2*NMAX
C***  If DEC-DXML Routines are used, it should be 64*NPDIM
C***  Here NPDIM=94 is assumed ==> 94 * 64 = 6016
C***  Or, for solving the statistical equations NDIM=300 ==> 300 * 64 = 19200
C      PARAMETER (NDIMSC = 19200)
      PARAMETER (NDIMSC = NMAX * 64)
      DIMENSION SCRATCH(NDIMSC)

C***  Array for special for the Use in DEC-DXML-Routines 
C***  DGETRF and DGETRI
      DIMENSION IPIV(NMAX)

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C***  OUTPUT of matrix dimension parameters (testing only)
ccc      write (0,*) '**** test INV> N=', N, 'NDIM=', NDIM

      IF (N .GT. NDIM) THEN
         PRINT *,'ERROR IN SUBROUTINE INV: N .GT. NDIM'
         CALL REMARK ('ERROR IN SUBROUTINE INV: N .GT. NDIM')
         STOP 'ERROR'
         ENDIF
      IF (NDIM .GT. NMAX) THEN
         WRITE (0,'(A)') '*** FATAL ERROR: NDIM .GT. NMAX'
         WRITE (0,'(A, I4, A, I4)') 'NDIM=', NDIM, ' .GT. NMAX=', NMAX
         STOP 'ERROR IN SUBROUTINE INV'
         ENDIF

C***  disable owninv for testing
      IF (CKEY(1:3) .EQ. 'OWN' .OR. CKEY(1:4) .EQ. 'OWNL') THEN
        CALL OWNINV(N, NDIM, A, CKEY) 
      ELSE
          if (ndim .gt. NMAX) then
            write (0,*) '*** obsolete program branch'
            write (0,*) '*** check DIMENSION SCRATCH'
            stop '*** FATAL ERROR in SUBR. INV' 
          endif 
          CALL DGETRF(N, N, A, NDIM, IPIV, INFO)
          CALL DGETRI(N, A, NDIM, IPIV, SCRATCH, NDIMSC, INFO)
      ENDIF

C***  Nur zum Loesen ginge auch folgendes
C***     Inverse Matrix stuende dann aber nicht zur Verfuegung
C***     Wert von Parameter 8 (LDB) ist unklar
C!!!          CALL DGETRS('N', N, N, A, NDIM, IPIV, B, NDIM, INFO)
C!!!          WRITE (0,*) 'INFO, SCRATCH(1)=', INFO, SCRATCH(1)

      RETURN
      END
      SUBROUTINE INVTRI (A,B,C,Q,N)
C***********************************************************************
C***  TRIDIAGONALE MATRIX   -A(L)   B(L)  -C(L)
C***  RECHTE SEITE Q(L)
C***  LOESUNG AUF VEKTOR Q(L)
C***  ACHTUNG -- AUCH C(L) WIRD VERAENDERT --
C***********************************************************************

      DIMENSION A(N),B(N),C(N),Q(N)
      DIMENSION AHELP(200),BHELP(200)

      IF (N .GT. 200) STOP '*** ERROR: dimension overflow in INVTRI'
      CI=C(1)/B(1)
      C(1)=CI
      QI=Q(1)/B(1)
      Q(1)=QI
      NM=N-1
      IF (N.EQ.1) RETURN
      IF (N.EQ.2) GOTO 3
      DO 1 I=2,NM
      AI=A(I)
      H=B(I)-AI*CI
      CI=C(I)/H
      C(I)=CI
      QI=(Q(I)+QI*AI)/H
    1 Q(I)=QI
    3 QI=(Q(N)+QI*A(N))/(B(N)-CI*A(N))
      Q(N)=QI
C**  THE BACKWARD ELIMINATION MAY BE SPEEDED UP BY CRAY VECTOR ROUTINE FOLR
      DO 4 L=2,N
      AHELP(L)=-C(N+1-L)
    4 BHELP(L)= Q(N+1-L)
      BHELP(1)=Q(N)
      CALL FOLR(N,AHELP,1,BHELP,1)
      DO 5 L=1,N
    5 Q(N+1-L)=BHELP(L)
      RETURN
C***  DEAD BRANCH: RECURSION BY HAND (CAN NOT BE AUTO-VECTORIZED)
      DO 2 I=1,NM
      L=N-I
      QI=Q(L)+C(L)*QI
    2 Q(L)=QI
      RETURN
      END
      FUNCTION ISAMAX (N, X, INC)

      ISAMAX = IDAMAX (N, X, INC)

c      DIMENSION X(N)

c      XMAX = ABS(X(1))
c      IMAX = 1

c      DO I=2*INC, N, INC
c        IF (ABS(X(I)) .GT. XMAX) THEN
c          XMAX = ABS(X(I))
c          IMAX = I
c        ENDIF
c      ENDDO

c      ISAMAX = IMAX

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

      FUNCTION ISRCHFGE(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHFGE - Searches a vector for the first element greater 
C***       than or equal to a target

C***  SYNOPSIS
C***  SEE ISRCHEQ

      REAL X(*), TARGET

      J=1
      ISRCHFGE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).GE.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFGE=I

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

      FUNCTION ISRCHFLE(N,X,INCX,TARGET)

C***  NAME
C***     ISRCHFLE searches a real vector for the first element that is less
C***     than or equal to a real target.

C***  SYNOPSIS
C***  SEE ISRCHLE

      REAL X(*), TARGET

      J=1
      ISRCHFLE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).LE.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFLE=I

      RETURN
      END

      FUNCTION ISRCHFLT(N,X,INCX,TARGET)

C***  NAME
C***     ISRCHFLT searches a real vector for the first element that is less
C***     than a real target.

C***  SYNOPSIS
C***  SEE ISRCHLT

      REAL X(*), TARGET

      J=1
      ISRCHFLT=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).LT.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFLT=I

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
      FUNCTION KHOLTSMARK(IQNLOW, IQNUP, IZ)
C***********************************************************************      
C***  returns the Holtsmark wing broadening coefficient K_low,up
C***  (required for the treatment of linear Stark broadening)
C***
C***  IZ = NUCLEAR CHARGE 
C***       (note: this is sometimes called "ion charge" in the literature)
C***
C***  table and calculation taken from H. Griem, 1960, ApJ 132, p.883
C***   (The Griem calculations are an approximative method 
C***    for the results described in Unsoeld, Page 320ff)
C***
C***  Please note that K_low,up is NOT dimensonless, but instead has
C***   the unit of [ c*h^7 / (m^3 * e^9) ]  = cm^(3/2) * s * g^(-1/2)
C***   [All constant units have to evaluated in cgs]
C***
C***  called by LINSTARK
C***********************************************************************      

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: IQNLOW, IQNUP, IZ
      REAL :: KHOLTSMARK
      INTEGER :: L1, U1

C**** Table 1 from Griem for KHOLTSMARK 
C        ( -1  = no valid entry )
      REAL, DIMENSION(2:18,1:4) :: KTAB = 
     >   RESHAPE(  (/
     >      0.293E-3,       -1.,     -1.,    -1.,
     >      0.560E-3,  0.143E-1,     -1.,    -1.,
     >      0.940E-3,  0.188E-1,   0.163,    -1.,
     >      1.430E-3,  0.262E-1,   0.174,    0.98,
     >      2.040E-3,  0.356E-1,   0.213,    0.91,
     >      2.750E-3,  0.470E-1,   0.267,    1.02,
     >      3.580E-3,  0.600E-1,   0.332,    1.20,
     >      4.510E-3,  0.720E-1,   0.406,    1.42,
     >      5.600E-3,  0.880E-1,   0.49 ,    1.68,
     >      6.700E-3,  1.100E-1,   0.58 ,    1.96,
     >      8.000E-3,  1.290E-1,   0.68 ,    2.28,
     >      9.400E-3,  1.520E-1,   0.80 ,    2.63,
     >     10.800E-3,  1.760E-1,   0.92 ,    3.01,
     >     12.400E-3,  2.010E-1,   1.04 ,    3.41,
     >     14.100E-3,  2.280E-1,   1.18 ,    3.84,
     >     16.000E-3,  2.570E-1,   1.33 ,    4.31,
     >     17.900E-3,  2.880E-1,   1.48 ,    4.8
     >  /) , (/ 17, 4 /), ORDER=(/ 2, 1 /) )
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      
      
C***  STOP if called with invalid quantum numbers      
      IF (IQNLOW <= 0) THEN
        WRITE (hCPR,*) 'ERROR: IQNLOW <= 0'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF
      IF (IQNUP <= 0) THEN
        WRITE (hCPR,*) 'ERROR: IQNUP <= 0'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF
C***  STOP if called with reserve QN roles
      IF (IQNLOW >= IQNUP) THEN
        WRITE (hCPR,*) 'ERROR: IQNLOW >= IQNUP NOT ALLOWED'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF

C***  Determine if table values can be used
      IF (IQNUP <= 18 .AND. IQNLOW <= 4) THEN
c        KHOLTSMARK = KTAB(IQNLOW, IQNUP)
        KHOLTSMARK = KTAB(IQNUP, IQNLOW)
      ELSE
C***    asymptotic formula for everything else
C        (note that Hubeny incorrectly writes 5.5 * 1.E-4 in his 1994 paper, 
C         but recalculating Griem's equations leads to 5.5 * 1.E-5 as he writes in his 1960 paper)
        KHOLTSMARK = 5.5 * 1.E-5 
     >        * (IQNUP*IQNLOW)**4 / (IQNUP*IQNUP - IQNLOW*IQNLOW)
      ENDIF
      
C***  Correction for hydrogenic ions beyond hydrogen
C***   after Hubeny & Mihalas (Book, 2015), page 260, Eq. (8.151)
C***   see also Eq. (33) in Griem (1960)
C***  Note: This factor can either be put on K_low,up or Holtsmark normal field strength F0
C***        SYNSPEC (formal integral for TLUSTY models) puts it in F0, Griem puts it on K_low,up
      KHOLTSMARK = KHOLTSMARK / (FLOAT(IZ)**5)
      
      
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
      SUBROUTINE LIMBDARK_OUTPUT (KANAL, LIMB_LINE, NFOBS, EMINT_P, 
     >                            NP_MERGED, PGRID_MERGED) 
C******************************************************************************
C***  PLOT OF THE LIMB DARKENING 
C******************************************************************************
      IMPLICIT NONE
      CHARACTER :: LIMB_LINE*(*), LIMBFILENAME*100
      INTEGER, INTENT(IN) :: NP_MERGED, NFOBS
      REAL, DIMENSION(NP_MERGED) :: EMINT_P, PGRID_MERGED, TAURCONT

C***  Local variables
      REAL, DIMENSION(NP_MERGED) :: XPLOT, YPLOT
      CHARACTER(100) :: HEAD1, HEAD2, STYLE, ACTPAR,  
     >                 XAXISSTR, YAXISSTR
      INTEGER :: KANAL, JP, NPAR, IPAR, IDX, NP_CUT, ISRCHFGT, ISRCHFLT
      INTEGER :: LIMB_UNIT
      REAL PCUT, EMINT_MIN

C***  Defaults: 
      LIMBFILENAME = ''
      XAXISSTR = '\CENTER\p [R&T*&M]'
      YAXISSTR = 
     >   '\CENTER\I&T#n#&M [erg cm&H-2&M s&H-1&M Hz&H-1&M sterad&H-1&M]'

      DO JP = 1, NP_MERGED
         XPLOT(JP) = PGRID_MERGED(JP)
         YPLOT(JP) = EMINT_P(JP)
      ENDDO

      NP_CUT = NP_MERGED


C***  Search for plot options 
      CALL SARGC(LIMB_LINE, NPAR)

      DO IPAR=2, NPAR
         CALL SARGV(LIMB_LINE, IPAR, ACTPAR)

         IF (ACTPAR(:4) .EQ. 'NORM') THEN
C                             ====
             DO JP=1, NP_MERGED
                YPLOT(JP) = YPLOT(JP) / EMINT_P(1)
             ENDDO
                YAXISSTR = '\CENTER\I&T#n#&M / I&T#n#&M(p=0)'

         ELSEIF (ACTPAR .EQ. 'MU') THEN
C                             ==
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFGT(NP_MERGED,PGRID_MERGED,1, 1.) - 1 )
            DO JP=1, NP_CUT
               XPLOT(JP) = SQRT(1. - PGRID_MERGED(JP)**2)
            ENDDO
            XAXISSTR = '\CENTER\#m#'

C***     Plot ends before impact parameter PCUT
         ELSEIF (ACTPAR .EQ. 'PCUT') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, ACTPAR)
            READ (ACTPAR, '(F12.0)',ERR=100) PCUT
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFGT(NP_MERGED,PGRID_MERGED,1,PCUT) - 1)

         ELSEIF (ACTPAR .EQ. 'MIN') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, ACTPAR)
            READ (ACTPAR, '(F12.0)',ERR=100) EMINT_MIN
            EMINT_MIN = EMINT_MIN * EMINT_P(1)
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFLT(NP_MERGED, EMINT_P, 1, EMINT_MIN) -1)

         ELSEIF (ACTPAR .EQ. 'FILE') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, LIMBFILENAME)

         ENDIF
      ENDDO
     

  110 CONTINUE
C***  HEADER  ------------------------------------------------------
      HEAD1 = 'LIMB DARKENING'
      HEAD2 = LIMB_LINE(:IDX(LIMB_LINE))

      WRITE (KANAL, '(A,A)') 'PLOT: ', LIMB_LINE
      WRITE (KANAL, '(A)') '\INBOX'

      CALL PLOTANFS (KANAL,HEAD1,HEAD2
     >        ,XAXISSTR
     >        ,YAXISSTR
     >        , 0., 0., 0., 0., 0.,.0
     >        , 0., 0., 0., 0., 0.,.0
     >        , XPLOT ,YPLOT, NP_CUT, 'COLOR=2')
      
      WRITE (0,'(A)') 'LIMBDARKENING plot written'

C***  ASCII file output (optional)
      IF (LIMBFILENAME .NE. '') THEN
         LIMB_UNIT = 52
         OPEN (UNIT=LIMB_UNIT, FILE=LIMBFILENAME, ACTION="write", 
     >            status="unknown")
         WRITE (LIMB_UNIT, '(A)') '# (1) ' // XAXISSTR(9:IDX(XAXISSTR)) 
         WRITE (LIMB_UNIT, '(A)') '# (2) ' // YAXISSTR(9:IDX(YAXISSTR)) 
         DO JP=1, NP_CUT
            WRITE (LIMB_UNIT, '(1P,G14.4,1X,G14.4)') 
     >             XPLOT(JP), YPLOT(JP)
         ENDDO
         CLOSE (UNIT=LIMB_UNIT)
      ENDIF

      RETURN

C***  ERROR branch ************************************************
  100 WRITE (0, 101) LIMB_LINE(:IDX(LIMB_LINE)) 
  101 FORMAT  ("*** LIMBDARK: SYNTAX ERROR with plot options" /
     >         "The error occured in the following line:" /, A, /,
     >         "Possible plot options:" /
     >         "LIMBDARK ... [PCUT=x.x] [MIN=x.x] [MU] [NORM]" /
     >         'Non-fatal error detected in subr. LIMBDARK_OUTPUT' /
     >         ' --> LIMBDARK plot option(s) ignored' )
      GOTO 110

      END
      SUBROUTINE LIMBDARK_PREP (LIMB_LINE, XOBS0, DXOBS, NFOBS,
     >                          XLAM, ALN, WEIGHT_LIMBDARK, BLIMB)

C***************************************************************************
C*** This routine decodes the LIMBDARKENING options
C*** and prepares the vector WEIGHT_LIMBDARK for the
C*** integration/extraction of the intensities I_nu(p) 
C*** The options allow to extract for either
C*** (1) a specific wavelength, 
C*** (2) a wavelenght range, or 
C*** (3) a broad-band color
C*** Syntax:
C*** (1) LIMB LAM x.x [plotoptions]
C***     provides I_nu(p) at given wavelength x.x 
C*** (2) LIMB RANGE x.x y.y [plotoptions]
C***     provides I_nu(p) averaged over the given wavelength x.x to y.y 
C*** (3) LIMB BAND U|B|V|u|b|v|y [plotoptions]
C***     provides I_nu(p) averaged over the specified filter 
C************************************************************************** 

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NFOBS 
      CHARACTER LIMB_LINE*(*), SELECTEDBAND*8 
      REAL, INTENT(IN) :: DXOBS, XOBS0, XLAM, ALN
      REAL WEIGHT_LIMBDARK(NFOBS)

C***  For the filter funcions:
      INTEGER, PARAMETER :: MAXNFILT = 200
      INTEGER NFILT
      REAL, DIMENSION(MAXNFILT) :: FLAM, FILT

      REAL XLAMTEST, FILTVAL

C***  LOCAL VARIABLES:
      CHARACTER ACTPAR*80
      REAL DUMMY, XLAM_LIMBDARK, XLAM1, XLAM2, XK, XK1, XK2, P, Q 
      REAL XOBS1, XLAMMIN, XLAMMAX, XBLUE, XRED, WEIGHTSUM, DNUEDX
      INTEGER :: I, J, L, NPAR, PLIPO, K, K1, K2, IDX
      LOGICAL BLIMB, BDEFINED
      
      BLIMB = .TRUE.

C***  frequency range is from XOBS0 to XOBS1
C***  X... are dimensionless frequencies
      XBLUE = -(XOBS0 - DXOBS)
      XRED  = -(XOBS0 - DXOBS*NFOBS)
C***  current LAMBDA  range
      XLAMMIN = XLAM * EXP(-XBLUE*ALN)
      XLAMMAX = XLAM * EXP(-XRED *ALN)

ccc      write (0,*) 'XBLUE=', XBLUE
ccc      write (0,*) 'XRED =', XRED
ccc      write (0,*) 'XLAMMIN=', XLAMMIN
ccc      write (0,*) 'XLAMMAX=', XLAMMAX


      DO K=1, NFOBS
         WEIGHT_LIMBDARK(K) = .0
      ENDDO

C***  Check: only one BAND or LAMBDA or RANGE can be specified !
      BDEFINED = .FALSE.

C***  Decode the line LIMB_LINE specified in FORMAL_CARDS
      CALL SARGC(LIMB_LINE, NPAR)
      IF ((NPAR .LT. 2) .OR. (NPAR .GT. 12)) GOTO 100

      DO I=2, NPAR
         CALL SARGV(LIMB_LINE, I, ACTPAR)

C*******************************************************************
         IF (ACTPAR(:3) .EQ. 'LAM') THEN
C                             ===
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, ACTPAR)
            READ (ACTPAR, '(F12.0)', ERR=100) XLAM_LIMBDARK

C***        Skip if XLAM_LIMBDARK-XLAMMIN not covered:
            IF ((XLAM_LIMBDARK-XLAMMIN)*(XLAM_LIMBDARK-XLAMMAX) .GE. .0) 
     >         GOTO 200

C*          Find non-integer K-index of requested XLAM_LIMBDARK
            XK   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM_LIMBDARK)) / ALN )
     >             / DXOBS

            K = INT(XK)
            Q = XK - FLOAT(K)
            P = 1. - Q
            WEIGHT_LIMBDARK(K)   = P
            WEIGHT_LIMBDARK(K+1) = Q
               

ccc               write (0,*) 'XLAM_LIMBDARK lies inside:', XLAM_LIMBDARK
ccc               xlamtest = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
ccc               write (0,*) 'K, lam =', K,   XLAMtest, 1-q
ccc               xlamtest = XLAM * EXP ((XOBS0 - (K+1)*DXOBS) * ALN)  
ccc               write (0,*) 'K, lam =', K+1, XLAMtest, q

C*******************************************************************
         ELSEIF (ACTPAR .EQ. 'RANGE') THEN
C                             =====
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+2 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, ACTPAR)
            READ(ACTPAR, '(F12.0)', ERR=100) XLAM1
            CALL SARGV(LIMB_LINE, I+2, ACTPAR)
            READ(ACTPAR, '(F12.0)', ERR=100) XLAM2

C***        Sort lam1 < lam2 
            IF (XLAM1 .GT. XLAM2) THEN
                DUMMY = XLAM1
                XLAM1 = XLAM2
                XLAM2 = DUMMY
            ENDIF   

C***        Skip RANGE if not covered:
            IF ((XLAM1-XLAMMIN)*(XLAM1-XLAMMAX) .GE. .0) GOTO 201 
            IF ((XLAM2-XLAMMIN)*(XLAM2-XLAMMAX) .GE. .0) GOTO 201 

C*          Find non-integer K-index of requested XLAM1
            XK1   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM1)) / ALN )
     >             / DXOBS

            K1 = INT(XK1)
            Q = XK1 - FLOAT(K1)
            P = 1. - Q
            WEIGHT_LIMBDARK(K1)   = P*P / 2.
            WEIGHT_LIMBDARK(K1+1) = P * (Q + 1.) / 2.
               
C*          Find non-integer K-index of requested XLAM2
            XK2   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM2)) / ALN )
     >             / DXOBS

            K2 = INT(XK2)
            Q = XK2 - FLOAT(K2)
            P = 1. - Q
            WEIGHT_LIMBDARK(K2)   = Q * (P + 1.) / 2.
            WEIGHT_LIMBDARK(K2+1) = Q*Q / 2.
            DO K=K1, K2
               IF (K .EQ. K1 .OR. K .EQ. K1) THEN 
                  WEIGHT_LIMBDARK(K) = 0.5
               ELSE
                  WEIGHT_LIMBDARK(K) = 1.
               ENDIF
            ENDDO

C***        Integration is done over \nu. Unlike \Delta x, \Delta\nu is 
C***        not constant but propto \nu (cf. manpowr) 
            DO K=K1, K2+1
               DNUEDX = EXP(ALN)**(XOBS0-K*DXOBS)
               WEIGHT_LIMBDARK(K)   = WEIGHT_LIMBDARK(K) * DNUEDX
            ENDDO

C*******************************************************************
         ELSE IF (ACTPAR .EQ. 'BAND') THEN
C                              ====
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, SELECTEDBAND)
            
            CALL FILTERFUNCTIONS (SELECTEDBAND, NFILT, FLAM, FILT)

C***        Skip RANGE if filter not covered:
            XLAM1 = FLAM(1)
            XLAM2 = FLAM(NFILT)
            IF ((XLAM1-XLAMMIN)*(XLAM1-XLAMMAX) .GE. .0) GOTO 202 
            IF ((XLAM2-XLAMMIN)*(XLAM2-XLAMMAX) .GE. .0) GOTO 202 

C*          Find non-integer K-index of first filter-wavelength:
            XK1   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM1)) / ALN )
     >             / DXOBS

            K1 = INT(XK1)
            Q = XK1 - FLOAT(K1)
            P = 1. - Q
            WEIGHT_LIMBDARK(K1)   = P*P / 2.
            WEIGHT_LIMBDARK(K1+1) = P * (Q + 1.) / 2.
               
C*          Find non-integer K-index of last filter-wavelength:
            XK2   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM2)) / ALN )
     >             / DXOBS

            K2 = INT(XK2)
            Q = XK2 - FLOAT(K2)
            P = 1. - Q
            WEIGHT_LIMBDARK(K2)   = Q * (P + 1.) / 2.
            WEIGHT_LIMBDARK(K2+1) = Q*Q / 2.
            DO K=K1, K2
               IF (K .EQ. K1 .OR. K .EQ. K1) THEN 
                  WEIGHT_LIMBDARK(K) = 0.5
               ELSE
                  WEIGHT_LIMBDARK(K) = 1.
               ENDIF
            ENDDO

C***        Integration is done over \nu. Unlike \Delta x, \Delta\nu is 
C***        not constant but propto \nu (cf. manpowr) 
ccc            DO K=K1, K2+1
            DO K=K1+1, K2
               XLAMTEST = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
               CALL SPLINPO (FILTVAL, XLAMTEST, FILT, FLAM, NFILT)
               DNUEDX = EXP(ALN)**(XOBS0-K*DXOBS)
               WEIGHT_LIMBDARK(K) = WEIGHT_LIMBDARK(K) * DNUEDX *FILTVAL
            ENDDO

         ENDIF
      ENDDO

      IF (.NOT. BDEFINED) GOTO 95

C***  Normalization of frequency integral
      WEIGHTSUM = .0
      DO K=1, NFOBS
         WEIGHTSUM = WEIGHTSUM + WEIGHT_LIMBDARK(K)
      ENDDO

      DO K=1, NFOBS
         WEIGHT_LIMBDARK(K) = WEIGHT_LIMBDARK(K) / WEIGHTSUM
      ENDDO

C***  Test output
      if (.false.) then 
         do k=1, nfobs
            if (WEIGHT_LIMBDARK(K) .EQ. .0) CYCLE
            xlamtest = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
            write (0,*) 'K, lam =', K,   XLAMtest, WEIGHT_LIMBDARK(K)
         enddo
      endif


      RETURN

C******************************************************************
C***  Skipping LIMBDARK if LAMBDA, RANGE or FILTER is not covered 
C***  by the present wavelength range

  200 WRITE (0,'(A, F8.0,A)') '*** WARNING: ' // 
     >      'LIMBDARK: LAM=', XLAM_LIMBDARK,  
     >      ' lies outside the current RANGE'
      GOTO 210

  201 WRITE (0,'(A, 2F8.0,A)') '*** WARNING: ' // 
     >      'LIMBDARK: interval LAM=', XLAM1, XLAM2,  
     >      ' not covered by current RANGE'
      GOTO 210

  202 WRITE (0,'(A, F8.0,A, F8.0, A)') '*** WARNING: ' // 
     >      'LIMBDARK: FILTER ' // SELECTEDBAND(:IDX(SELECTEDBAND)) //
     >      ' (', XLAM1, ' - ', XLAM2, 
     >             ') not covered by current RANGE'
      GOTO 210

  210 CONTINUE
      BLIMB = .FALSE.
      WRITE (0,'(A, F10.0,A,F10.0)') '*** --> LIMBDARKENING skipped ' //
     >            '- current range is ', XLAMMIN, ' to ', XLAMMAX
      RETURN 

C***  ERROR BRANCHES ************************************************

   90 WRITE (0,'(A)') 
     >  '*** ERROR: only one LAMBDA, RANGE or BAND can be selected'
      GOTO 110

   95 WRITE (0,'(A)') 
     >  '*** ERROR: either LAMBDA, RANGE or BAND must be selected'
      GOTO 110

  100 WRITE(0,'(A)') "*** ERROR: Wrong syntax"
      GOTO 110

  110 CONTINUE
      WRITE(0,'(A)') "The error occured in the following line"
      WRITE(0,'(A)') LIMB_LINE(:IDX(LIMB_LINE))
      WRITE(0,'(A)') "Possible syntax:"
      WRITE(0,'(A)') 'LIMB LAM x.x [plotoptions]'
      WRITE(0,'(A)') 'LIMB RANGE x.x y.y [plotoptions]'
      WRITE(0,'(A)') 'LIMB BAND U|B|V|u|b|v|y [plotoptions]'

      WRITE(0,'(A)')  'Fatal error in subroutine LIMBDARK_PREP'
      WRITE(0,'(A)')  '--> LIMB DARKENING not calculaated (skipped)'
      
      BLIMB = .FALSE.
      RETURN

      END
      SUBROUTINE LINSTARK (GRIEMPAR, XNE, NUCCHARG, 
     >                     MAINQNLOW, MAINQNNUP, XLAM0,
     >                     LINPRO, LEVELLOW, LEVELNUP)
C**********************************************************************
C***  Linear Stark effect, prepartion routine for GRIEMPAR vector
C***     motivated by the approach from Hubeny et al. (1994)
C***  
C***  GRIEMPAR             - conversion factor for Griem's beta parameter
C***                               = lamda_0/(F_0 * K_low,up)
C***  XNE                  - electron density
C***  NUCCHARG             - nuclear charge 
C***                           (often misleadingly labeled ion charge)
C***  MAINQNLOW, MAINQNNUP - main quantum numbers of lower and upper level 
C***  XLAM0                - transition wavelength in Angstroem
C***
C***  returns GRIEMPAR in unit of:  Angstroem/cm
C***    (i.e. needs to be mulitplied with _relative_ frequency difference)
C***
C***  written by ansander (Dec 2015)
C***
C***  called from STARKBROAD
C**********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAINQNLOW, MAINQNNUP, NUCCHARG
      
      CHARACTER(8) :: LINPRO
      CHARACTER(10) :: LEVELLOW, LEVELNUP
      CHARACTER(40) :: CQNERRDETAIL

      REAL, EXTERNAL :: KHOLTSMARK
      
      REAL :: GRIEMPAR, F0, T, XNE, XLAM0, HOLTSK
            
C***  Holtsmark field strength normalization factor      
C***  1.25 = 2 * Pi * (4/15)^(2/3) * e 
      REAL, PARAMETER :: FNORM = 1.25E-9            ! in cm^(3/2) * g^(1/2) * s^(-1)
      

C***  The main quantum numbers of both levels need to be known for this method!
      IF (MAINQNLOW <= 0 .AND. MAINQNNUP <= 0) THEN
        CQNERRDETAIL = '  (both main quantum numbers not known) '
      ELSEIF (MAINQNLOW <= 0) THEN
        CQNERRDETAIL = '  (lower quantum number not known)      '
      ELSEIF (MAINQNNUP <= 0) THEN
        CQNERRDETAIL = '  (upper quantum number not known)      '
      ENDIF

C***  Fallback to quadratic Stark effect in case of unknown quantum numbers      
      IF (MAINQNLOW <= 0 .OR. MAINQNNUP <= 0) THEN
           LINPRO = 'Q-STARK '
           GRIEMPAR = .0
           WRITE (0,'(5A)') '*** WARNING: LINSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP,  CQNERRDETAIL
           WRITE (*,'(5A)') '*** WARNING: LINSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP,  CQNERRDETAIL
           RETURN
      ENDIF 

C***  F0 = Holtsmark normal field strength
      F0 = FNORM * XNE**(2./3.)
C***  HOLTSK is Griems K_low,up correction factor, generalized with
C***    a charge factor to approximate for all hydrogenic ions
      HOLTSK = KHOLTSMARK(MAINQNLOW, MAINQNNUP, NUCCHARG)

C***  GRIEMPAR contains the transition-dependent quantities which can
C***  be precalculated before the wavelength integration
C***  The factor XLAM0 in the nominator is neeed to convert the relative wavelengths used
C***  in the integration in ZONEINT back into absolute wavelengths in Angstroem which
C***  are required by the PHISTARK routine.
      GRIEMPAR = XLAM0 / (F0 * HOLTSK)
      
              
      RETURN
      END
      SUBROUTINE LIOP (EINST,WEIGHTI,WEIGHTJ,I,J,
     $           ND,XLAM,ENTOT,POPNUM,RSTAR,OPAL,ETAL,VDOP)
C***********************************************************************
C***  CALCULATES THE LINE OPACITY OPAL AND THE EMISSIVITY ETAL FOR ONE LINE
C***  AT ALL DEPTH POINTS
C***  PHYSICAL DIMENSION OF OPAL : PER RSTAR AND PER DELTA-NUE-DOPPLER
C***    ( THE LATTER WILL BE CANCELLED OUT BY THE NORMALISED PROFILE FUNCTION )
C***  GIVEN QUANTITIES : EINST = EINSTEIN COEFFICIENT A(UP-LOW)  (PER SECOND)
C***                     WEIGHTI, WEIGHTJ = STATISTICAL WEIGHTS (UP, LOW)
C***                     POPNUM(L,J) = RELATIVE POPULATION NUMBERS
C***                     ENTOT(L) = TOTAL NUMBER DENSITY
C***                     I = LOW,  J = UP
C***  called from:
C***       COLI -> CHECK_LINES -> LIOP
C***      STEAL -> LINPOP -> LCORE -> LIOP
C***      STEAL -> LINPOP -> COMA -> SETXJL -> LIOP
C***      STEAL -> LINPOP -> COMA -> SETXJLCF -> LIOP
C***      STEAL -> LINPOP -> COMA -> SETXJFINE -> LIOP
C***     FORMAL -> LIOP
C***********************************************************************
 
      DIMENSION ENTOT(ND),OPAL(ND),ETAL(ND),POPNUM(ND,1)
      REAL :: ENI, ENJ, VDOP

C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  PI8 = 8*PI
      DATA PI8 /25.1327412288 /

C***  FOR EXACT COMPATIBILITY, CALCULATE C3 = H*C/(4*PI) :
      C3=C2/PI8

C***  WAVELENGTH IN CENTIMETER
      XLAMCM=XLAM/1.E8

C***  DND= DELTA-NUE-DOPPLER  (HERTZ)
      DND=VDOP/XLAM*1.E13

      EMINDU=EINST*XLAMCM*XLAMCM/PI8*RSTAR
      ABSORP = EMINDU*WEIGHTJ/WEIGHTI
      EMSPON=C3*RSTAR*EINST/XLAMCM

      DO L=1,ND
         ENI=POPNUM(L,I)*ENTOT(L)
         ENJ=POPNUM(L,J)*ENTOT(L)
C***     Set emissivities zero if both levels are equal (=POPMIN)
         IF (ENI .EQ. ENJ) ENJ=.0
         OPAL(L)=(ENI*ABSORP-ENJ*EMINDU)/DND
         ETAL(L)=ENJ*EMSPON/DND
      ENDDO

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
      SUBROUTINE MACROCLUMP (OPA, POROLENGTH, CORRFAC_MACROCLUMP)
C*****************************************************************
C***  For a given opacity OPA and porosity length, this routine 
C***  returns the correction factor: OPA_EFF = OPA * CORRFAC_MACROCLUMP 
C***  due to clumps of the optical thickness TAUCLUMP = OPA * POROLENGTH
C*****************************************************************

      IF (POROLENGTH .EQ. .0) THEN
         CORRFAC_MACROCLUMP = 1.
         RETURN
      ENDIF
 
      TAUCLUMP = OPA * POROLENGTH
 
      IF (TAUCLUMP .GT. 1.E-10) THEN
         CORRFAC_MACROCLUMP = (1. - EXP(-TAUCLUMP)) / TAUCLUMP
      ELSE IF (TAUCLUMP .GT. .0) THEN
         CORRFAC_MACROCLUMP = 1. - TAUCLUMP * 0.5
      ELSE
         CORRFAC_MACROCLUMP = 1.
      ENDIF
 
      RETURN
      END
      SUBROUTINE MANIPOP (ENLTE, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >              ABXYZ, NFIRST, NLAST, NATOM, POPNUM, RNE, ENTOT, 
     >              N, ND, MANIPOP_OPTIONS, DENSCON, T, level)
C************************************************************************
C***  Manipulation of popnumbers, in order to simulate the emission of 
C***  a hot component. 
C************************************************************************

      DIMENSION RNE(ND), ENTOT(ND), POPNUM(ND,N), ENLTE(N), T(ND)
      PARAMETER ( MAXPAR = 6 )
      CHARACTER MANIPOP_OPTIONS*(*), PAR(MAXPAR)*20, level(N)*(*)
      LOGICAL BTCOOL
      DIMENSION DENSCON(ND)

      BTCOOL = .FALSE.

C***  THOT    = Temperature of the hot component (optional 'TCOOL') 
C***  DENSHOT = Density of the hot component, relative to the 
C***            density of the cool component including the clumping factor
C***  HOTMASS = Mass fraction of the hot material

      CALL SARGC (MANIPOP_OPTIONS, NPAR)      
      IF (NPAR .EQ. 0) RETURN

      DO I=1, NPAR
         CALL SARGV (MANIPOP_OPTIONS, I, PAR(I))
      ENDDO       

      THOT    = -1. 
      DENSHOT = -1.
      HOTMASS = -1.
      DO I=1, NPAR-1
        IF (PAR(I) .EQ. 'THOT'   )  THEN
           IF (PAR(I+1) .EQ. 'TCOOL') THEN
              BTCOOL = .TRUE.
           ELSE
              READ (PAR(I+1),'(G20.0)') THOT
           ENDIF
        ENDIF
        IF (PAR(I) .EQ. 'DENSHOT') READ (PAR(I+1),'(G20.0)') DENSHOT
        IF (PAR(I) .EQ. 'HOTMASS') READ (PAR(I+1),'(G20.0)') HOTMASS
      ENDDO
      IF (THOT .LE. 0. .AND. .NOT. BTCOOL) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER THOT'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 
      IF (DENSHOT .LE. 0.) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER DENSHOT'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 
      IF (HOTMASS .LE. 0.) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER HOTMASS'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 

      DO L=1, ND
        ENE = RNE(L) * ENTOT(L) * DENSCON(L) * DENSHOT 
        IF (BTCOOL) THOT = T(L)
        
C***    No manipulation if 'HOT' component not hotter than the cool one 
        IF (THOT .LT. T(L)) CYCLE

        CALL LTEPOP (N,ENLTE,THOT,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)

        DO I=1, N
           POPNUM(L,I) = POPNUM(L,I) + HOTMASS * ENLTE(I)
        ENDDO

cccc    test output
c        if (l .eq. 30 ) then
c           write (0, '(A,I3)') '*** L =', L
c           do i=1,n
c            write (0, '(A,X,G14.4)') level(i), enlte(i)
c           enddo
c        endif

      ENDDO

      RETURN
      END
      SUBROUTINE MDMV (A,B,JMAX,NP)
C***********************************************************************
C***  MATRIX (DIAGONAL)  A  *  MATRIX (VOLL)  B
C***  ERGEBNIS-MATRIX UEBERSCHREIBT  B
C***********************************************************************

      DIMENSION A(NP),B(NP,NP)

      DO 1 I=1,JMAX
      AI=A(I)
      DO 1 K=1,JMAX
    1 B(I,K)=B(I,K)*AI
      RETURN
      END
      SUBROUTINE MDV (A,W,N)
C***********************************************************************
C*** MATRIX A (DIAGONAL)  *  VEKTOR W
C***  ERGEBNIS-VEKTOR UEBERSCHREIBT  W
C***********************************************************************

      DIMENSION  A(N),W(N)

      DO 1 I=1,N
    1 W(I)=A(I)*W(I)

      RETURN
      END
         SUBROUTINE MERGE_RGRID (RADIUS, NDDIM, MAXMOD, ND, 
     >      RADIUS_MERGED, ND_MERGED, PGRID, NP, NPDIM, PGRID_MERGED,
     >                     NP_MERGED, ZGRID_MERGED, SECMOD_RRANGE)
C******************************************************
C***  The RADIUS-Grids of each model are optimized; 
C***  In order to secure the numerical accuracy of the integration
C***  in case of the combination of two models, these grids must
C***  be combined. 
C******************************************************

      DIMENSION RADIUS(NDDIM,MAXMOD), RADIUS_MERGED(NDDIM)
      DIMENSION ND(MAXMOD)
      DIMENSION PGRID(NPDIM), PGRID_MERGED(NPDIM) 
      DIMENSION ZGRID_MERGED(NDDIM*NPDIM)
      REAL, DIMENSION(2) :: SECMOD_RRANGE

      RMAX = MIN (RADIUS(1,1), RADIUS(1,2))
      RADIUS_MERGED(1) = MIN (RADIUS(1,1), RADIUS(1,2))

      IF (RADIUS(1,1) .NE. RADIUS(1,2)) THEN
         WRITE (0,'(A)')      'Note:  RMAX of the two models differ:'
         WRITE (0,'(A,F6.1)') '         Main model: RMAX=', RADIUS(1,1) 
         WRITE (0,'(A,F6.1)') '       Second model: RMAX=', RADIUS(1,2) 
         WRITE (0,'(A,F6.1)') '--> Minimum adopted: RMAX=', RMAX 
      ENDIF

C***  The "density" of RADIUS-Points is compared between the two models;
C***  the merged-grid points are taken from that model where the points
C***  are denser. 
C***  At radii outside the SECMOD_RRANGE, only model 1 (the main model) 
C***  is taken; this can happen when the second-model region has the
C***  SHAPE od a SHERE

      L = 2
C***  RADIUS grids are restricted to values below the new RMAX 
      L1= ISRCHFLT(ND(1), RADIUS(1,1), 1, RADIUS_MERGED(1))
      L2= ISRCHFLT(ND(2), RADIUS(1,2), 1, RADIUS_MERGED(1))

    1 CONTINUE
      IF (RADIUS(L1-1,1) - RADIUS(L1+1,1) .LT.
     >    RADIUS(L2-1,2) - RADIUS(L2+1,2) .OR. 
     >    (RADIUS(L1,1)-SECMOD_RRANGE(1)) * 
     >    (RADIUS(L1,1)-SECMOD_RRANGE(2)) .GE. .0) THEN
          RADIUS_MERGED(L) = RADIUS(L1,1)
ccc          write (0,*) 'from 1: L, RADIUS_MERGED =', L, RADIUS_MERGED(L)
      ELSE
          RADIUS_MERGED(L) = RADIUS(L2,2)
ccc          write (0,*) 'from 2: L, RADIUS_MERGED =', L, RADIUS_MERGED(L)
      ENDIF


      L1 = ISRCHFLT(ND(1), RADIUS(1,1), 1, RADIUS_MERGED(L))
      L2 = ISRCHFLT(ND(2), RADIUS(1,2), 1, RADIUS_MERGED(L))
      L  = L  + 1
      IF (L .GT. NDDIM) THEN
         WRITE (0,*) 'NDDIM insufficient for merging RADIUS-Grids'
         STOP 'FATAL ERROR in Subr. MERGE_RGRID'
      ENDIF
      IF (L1 .LT. ND(1) .AND. L2 .LT. ND(2)) GOTO 1

      RADIUS_MERGED(L) = 1.
      ND_MERGED = L

      WRITE (0,'(A)') 'MERGED RADIUS GRID CONSTRUCTED'
      WRITE (0,'(A,I4)') 'MODEL 1: ND =', ND(1)
      WRITE (0,'(A,I4)') 'MODEL 2: ND =', ND(2)
      WRITE (0,'(A,I4)') 'Merged : ND =', ND_MERGED

C***  Construct PGRID_MERGED on basis of the new RADIUS_MERGED
C***  Note: PGRID has actually two vectors for each model, 
C***        but only the first one is used here to maintain the core-rays 
C***        which might have been changed by ROTATION_PREP

      NCORE = NP - ND(1)
      DO JP = 1, NCORE
         PGRID_MERGED(JP) = PGRID(JP)
      ENDDO

      NP_MERGED = NCORE + ND_MERGED
      IF (NP_MERGED .GT. NPDIM) THEN
         WRITE (0,'(A)') '*** DIMENSION NPDIM insufficient!'
         WRITE (0,'(A,I4)') '*** presently: NPDIM =', NPDIM
         WRITE (0,'(A,I4)') '*** needed   : NPDIM =', NP_MERGED
         STOP '*** FATAL ERROR IN SUBR. MERGE_RGRID'
      ENDIF

      DO L=1,ND_MERGED
         JP = NP_MERGED + 1 - L
         PGRID_MERGED(JP) = RADIUS_MERGED(L)
      ENDDO


C***  Z-array must be updated with new RADIUS_MERGED and the new
C***        PGRID_MERGED
C***  Note: the Z array is filled in the upper-left corner 
C***  The impact-parameters P are taken from the first (main) model,
C***  where the core-rys might have been increased in ROTATION_PREP
      DO L = 1, ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED + 1 - L
        DO JP = 1, JMAX
          PJ = PGRID_MERGED(JP)
          PJPJ = PJ * PJ
          I=(JP-1)*ND_MERGED+L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I) = SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I) = .0
          ENDIF
        ENDDO
      ENDDO

C***  test output
c      jp = 1
c      LMAX=MIN0(NP_MERGED+1-JP,ND_MERGED)
c      do l=1, lmax
c         I=(JP-1)*ND_MERGED+L
c         write (0,*) 'L, ZGRID_MERGED =', L, ZGRID_MERGED(i)
c      enddo 

      RETURN
      END
      SUBROUTINE MOMENT0 (ND,R,L,JMAX,Z,XJ,XJMEAN,MODE)
C***********************************************************************
C***  INTEGRATION OF THE ZERO-MOMENT OF THE RADIATION FIELD (MEAN INTENSITY)
C***   IF ( MODE = .TRUE.    ) THE INTEGRATION WEIGHTS ARE GENERATED (VECTOR XJ)
C***   ELSE : XJ IS CONSIDERED AS ANGLE-DEPENDENT INTENSITY AND THE
C***          INTEGRATION IS PERFORMED ( RESULT XJMEAN)
C***  RADIUS-MESH R, ACTUAL INDEX L, AND Z-MESH Z(L,J) ARE GIVEN
C***  WEIGHTS ARE ACCORDING TO TRAPEZOIDAL RULE IN Z=SQRT(R*R-P*P)
C***********************************************************************
 
      DIMENSION R(ND),Z(ND,JMAX),XJ(JMAX)
      LOGICAL MODE
      RL2=2.*R(L)
C***  FIRST STEP
      ZJ=Z(L,1)
      ZNEXT=Z(L,2)
      IF (MODE) XJ(1)=(ZJ-ZNEXT)/RL2
      IF (.NOT.MODE) XJMEAN=(ZJ-ZNEXT)*XJ(1)
C***  MIDDLE STEPS
      DO 1 J=3,JMAX
      ZLAST=ZJ
      ZJ=ZNEXT
      ZNEXT=Z(L,J)
      IF (MODE) XJ(J-1)=(ZLAST-ZNEXT)/RL2
      IF (.NOT.MODE) XJMEAN=XJMEAN+XJ(J-1)*(ZLAST-ZNEXT)
    1 CONTINUE
C***  LAST STEP, IMPLYING Z(L,JMAX)=.0
      IF (MODE) GOTO 2
      XJMEAN=XJMEAN+XJ(JMAX)*ZJ
      XJMEAN=XJMEAN/RL2
      RETURN
    2 XJ(JMAX)=ZJ/RL2
      RETURN
      END
      SUBROUTINE MOMENT1 (R,NP,P,U,H)
C***********************************************************************
C***  CALCULATES AN ANGLE INTEGRAL H OF THE RADIATION FIELD U
C***  BESIDES OF THE OUTER BOUNDARY, THIS IS NOT THE 1. MOMENT H,
C***  BUT RATHER AN INTENSITY-LIKE QUANTITY
C***  INTEGRATION WITH TRAPEZOIDAL RULE, WEIGHTS P * DP
C***    -  VECTORIZING VERSION, REPLACED 26-MARCH-1991
C***********************************************************************
 
      DIMENSION U(NP),P(NP)

C***  FIRST POINT
      A = P(1)
      B = P(2)
      W = (B - A) * (B + 2. * A)
      H = W * U(1)

C***  INNER POINT
      DO 2 J=2,NP-1
      A = P(J-1)
      B = P(J)
      C = P(J+1)
      W = (A + B + C) * (C - A)
      H = H + W * U(J)
    2 CONTINUE

C***  LAST POINT
      A = P(NP-1)
      B = P(NP)
      W = (B - A) * (2. * B + A)
      H = H + W * U(NP)
      H = H / (R * R * 6.)

      RETURN
      END
      SUBROUTINE MOMENT2 (R,JMAX,P,U,XK)
C***********************************************************************
C***  INTEGRATION OF THE 2. MOMENT XK OF THE RADIATION FIELD U
C***  FEAUTRIER-INTENSITY U(J), IMPACT PARAMETER MESH P(J)
C***  AND RADIUS POINT R ARE GIVEN.
C***  WEIGHTS ARE ACCORDING TO TRAPEZOIDAL RULE IN Z*Z*DZ, Z=SQRT(R*R-P*P)
C***********************************************************************

      DIMENSION P(JMAX),U(JMAX)

c     DON'T change this back, fp consistency needed!
c      RR=R*R
C***  FIRST STEP, IMPLYING P(1)=0
      Z=R
      ZQ=R*R
      PJ=P(2)
      ZNQ=R*R-PJ*PJ
      ZNQ=MAX(ZNQ,0.)
      ZNEXT=SQRT(ZNQ)
      W=Z*(3*ZQ-ZNQ)-ZNEXT*(ZQ+ZNQ)
      XK=W*U(1)
C***  MIDDLE STEPS
      DO 1 J=3,JMAX
      ZLAST=Z
      ZLQ=ZQ
      Z=ZNEXT
      ZQ=ZNQ
      PJ=P(J)
      ZNQ=R*R-PJ*PJ
      ZNQ=MAX(ZNQ,0.)
      ZNEXT=SQRT(ZNQ)
      W=Z*(ZLQ-ZNQ)+ZLAST*(ZLQ+ZQ)-ZNEXT*(ZQ+ZNQ)
    1 XK=XK+W*U(J-1)
C***  LAST STEP, IMPLYING P(JMAX)=R
      W=Z*ZQ
      XK=XK+W*U(JMAX)
c     old:
c      XK=XK/R/RR/12.
      XK=XK/R/R/R/12.

      RETURN
      END
      SUBROUTINE MSUB (A,B,JMAX,NP)
C***********************************************************************
C***  A := A - B
C***********************************************************************

      DIMENSION A(NP,NP),B(NP,NP)

      DO 1 K=1,JMAX
      DO 1 I=1,JMAX
    1 A(I,K)=A(I,K)-B(I,K)

      RETURN
      END
      SUBROUTINE MULTIPLE (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP,
     >                    DELXLAP,NBLINE,MAXLAP,INDLOW,INDNUP,LASTIND,
     >                    MAXIND,LEVEL,WEIGHT,ELEVEL,N,EINST,NDIM,
     >                    POPNUM,T,ND,ALN,VDOP, 
     >                    MAXSUBL,NSUBLOW,NSUBNUP,BROAD,LINPRO,AVOIGT,
     >                    NMOD, MAXMOD, NDDIM, 
     >                    MAINQN, NCHARG, EION, NOM, IND_ORIGLEV)
C***********************************************************************
C***  SUBROUTINE FOR MULTIPLET HANDLING OF MAIN PROGRAM "FORMAL"
C***  CALLED FROM SUBR. PREFORM IN CASE OF DECODED OPTION "MULTIPLET"
C***  INPUT OPTIONS: "/LOWER"   - DATA FOR LOWER SUBLEVEL
C***                 "/UPPER"   - DATA FOR UPPER SUBLEVEL
C***                 "/SUBLINE" - DATA FOR SUBLINE HANDLING
C***                 "-MULTI"   - END OF MULTIPLET INPUT 
C***  ACTION: 1. READ INPUT DATA (LEVELS, LINES) FOR MULTIPLET HANDLING
C***          2. ARRANGE SUBLINES IN A SEQUENCE OF INCREASING FREQUENCIES
C***          3. SPLIT ORIGINAL ENERGY LEVELS INTO SPECIFIED SUBLEVELS
C***             WITH POPNUMBERS CALCULATED FROM THE BOLTZMANN FORMULA
C***             (I.E. LTE POPNUMBERS)
C***********************************************************************

C***  DEFINE FORTRAN CHANNEL FOR LOGFILE MESSAGES
      PARAMETER ( LOGCHAN = 0 )
    
      DIMENSION ND(NMOD)
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,NMOD)
      DIMENSION WEIGHT(N), ELEVEL(N), MAINQN(N), NOM(N)
      DIMENSION NCHARG(N), EION(N)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION POPNUM(NDDIM,NDIM,NMOD),T(NDDIM,NMOD)
      DIMENSION NSUBLOW(MAXSUBL),NSUBNUP(MAXSUBL)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(N), LEV, LEVUP, LEVLOW
      CHARACTER*8 LINPRO(MAXLAP)
      LOGICAL BROAD, BAIR, BVAC

C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /
      
C***  DEFAULTS FOR MULTIPLET HANDLING:
      NBLSAVE=NBLINE
C***  NLOW, NNUP: NUMBER OF LOWER, UPPER SUBLEVELS
      NLOW=0
      NNUP=0
C***  NSUBLOW, NSUBNUP: POINTER TO SUBLEVELS IN THE ORIGINAL ARRAYS
      DO 10 I=1,MAXSUBL
      NSUBLOW(I)=0
   10 NSUBNUP(I)=0

C***  1. READ INPUT FOR MULTIPLET HANDLING OF SPECIFIED LINE  ----------
    1 READ (2,2) KARTE
      IF (KARTE(1:1) .EQ. '*' .OR. KARTE .EQ. ' ') GOTO 1
    2 FORMAT (A)
      IF (KARTE(:6) .EQ. '-MULTI') GOTO 20
C                         ======
      IF (KARTE(:6) .EQ. '/LOWER') THEN
C                         ======
         READ (KARTE,3) LEV,NW,ELEV
    3    FORMAT (12X,A10,1X,I4,1X,F20.0)

C***     FIND LOWER INDEX:
         JLOW=0
         DO 200 J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV).OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
                  WRITE (0,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (0,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (*,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ELSE
                 JLOW=J
               ENDIF
            ENDIF
 200     CONTINUE

C***     LEVEL NOT YET DEFINED
         IF (JLOW .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *,
     >            ' >>>>> SUBR. MULTIPLE: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('MULTIPLE: N GREATER THAN NDIM')
               STOP 'NDIM1'
            ENDIF
            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
            MAINQN(N)=MAINQN(LOW)
            NCHARG(N)=NCHARG(LOW)
            NOM(N) = NOM(LOW)
            EION(N)  =EION(LOW)
            IND_ORIGLEV(N) = LOW
            JLOW=N
         ENDIF

         NLOW=NLOW+1
         IF (NLOW .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. MULTIPLE: ERROR STOP (NLOW .GT. MAXSUBL)'
            CALL REMARK ('MULTIPLE: NLOW GREATER THAN MAXSUBL')
            STOP 'NLOW'
         ENDIF
         NSUBLOW(NLOW)=JLOW

      ELSE IF (KARTE(:6) .EQ. '/UPPER') THEN
C                              ======
         READ (KARTE,3) LEV, NW, ELEV 

C***     FIND UPPER INDEX:
         JNUP=0
         DO 201 J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV).OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
                  WRITE (0,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (0,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (*,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ELSE
                 JNUP=J
               ENDIF
            ENDIF
 201     CONTINUE

C***     LEVEL NOT YET DEFINED
         IF (JNUP .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *,
     >            ' >>>>> SUBR. MULTIPLE: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('MULTIPLE: N GREATER THAN NDIM')
               STOP 'NDIM2'
            ENDIF
            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
            MAINQN(N)=MAINQN(NUP)
            NCHARG(N)=NCHARG(NUP)
            NOM(N) = NOM(NUP)
            EION(N)  =EION(NUP)
            IND_ORIGLEV(N) = NUP
            JNUP=N
         ENDIF

         NNUP=NNUP+1
         IF (NNUP .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. MULTIPLE: ERROR STOP (NNUP .GT. MAXSUBL)'
            CALL REMARK ('MULTIPLE: NNUP GREATER THAN MAXSUBL')
            STOP 'NNUP'
         ENDIF
         NSUBNUP(NNUP)=JNUP

      ELSE IF (KARTE(:8) .EQ. '/SUBLINE') THEN
C                              ========
         LASTIND=LASTIND+1
         IF (LASTIND .GT. MAXIND) THEN
            PRINT *,
     >        ' >>>>> SUBR. MULTIPLE: ERROR STOP (LASTIND .GT. MAXIND)'
            CALL REMARK ('MULTIPLE: LASTIND GREATER THAN MAXIND')
            STOP 'MAXIND'
         ENDIF
         IF (NBLINE+1 .GT. MAXLAP) THEN
            PRINT *,
     >          ' >>>>> SUBR. MULTIPLE: ERROR STOP (NBLINE .GT. MAXLAP)'
            CALL REMARK ('MULTIPLE: NBLINE GREATER THAN MAXLAP')
            STOP 'MAXLAP'
         ENDIF

         READ (KARTE,23) LEVUP,LEVLOW,AUPLOW
   23    FORMAT (9X,A10,2X,A10,G10.0)

C***     FIND UPPER INDEX:
         DO 24 J=1,N
         JNUP=J
         IF (LEVEL(J) .EQ. LEVUP) GOTO 25
   24    CONTINUE

   90    FORMAT ('*** ERROR: UPPER LINE LEVEL NOT FOUND: ', A10)
         WRITE (LOGCHAN, 90)  LEVUP
         STOP 'ERROR STOP IN SUBR. MULTIPLE'

C***     FIND LOWER INDEX:
   25    DO 26 J=1,N
         JLOW=J
         IF (LEVEL(J) .EQ. LEVLOW) GOTO 27
   26    CONTINUE

   91    FORMAT ('*** ERROR: LOWER LINE LEVEL NOT FOUND: ', A10)
         WRITE (LOGCHAN, 91) LEVLOW
         STOP 'ERROR STOP IN SUBR. MULTIPLE'

   27    INDNUP(LASTIND)=JNUP
         INDLOW(LASTIND)=JLOW
C***     default wavelength 
         XLAMSUB=1.E8/(ELEVEL(JNUP)-ELEVEL(JLOW))

C***     read wavelength if given, covert from AIR to VAC
C***     VOIGT parameter?
         CALL SARGC (KARTE(43:),NPAR)
         IF (NPAR .GT. 0) 
     >     CALL READ_LINECARD_PARAMETERS (KARTE(43:), XLAMSUB,
     >            BROAD, LINPROBL, AVOIGTBL)

         IF (NBLINE .EQ. 0) XLAM = XLAMSUB 

         IF (AUPLOW .LT. 0.) THEN
            EINST(JNUP,JLOW)=-6.669E15/XLAMSUB/XLAMSUB*WEIGHT(JLOW)
     /                                      /WEIGHT(JNUP)*AUPLOW
         ELSE
            EINST(JNUP,JLOW)=AUPLOW
         ENDIF


C***     Insert current SUBLINE in the list sorted by increasing DELTAX
         CALL INSERT_LINE (LINE, LASTIND, NBLINE, INDLAP, XLAMLAP,DELXLAP,
     $                   XLAMSUB, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

      ELSE
         WRITE (0,*) 'UNRECOGNIZED INPUT CARD IN SUBR. MULTIPLE!'
         WRITE (0,*) 'KARTE=', KARTE( :IDX(KARTE))
      ENDIF

      GOTO 1

   20 CONTINUE
C***  END OF INPUT FOR MULTIPLET HANDLING  -----------------------------

C***  NO SUBLINES DECODED
      IF (NBLINE .EQ. NBLSAVE) THEN
         WRITE (0,'(A)') 
     >    ' >>>>> SUBR. MULTIPLE: WARNING - NO SUBLINE IN MULTIPLET:'
         WRITE (0,'(1X,A,1X,A)') LEVEL(NUP) , LEVEL(LOW)
         WRITE (*,'(A)') 
     >    ' >>>>> SUBR. MULTIPLE: WARNING - NO SUBLINE IN MULTIPLET:'
         WRITE (*,'(1X,A,1X,A)') LEVEL(NUP) , LEVEL(LOW)
         GOTO 99
      ENDIF

C***  3. SPLITTING OF THE ENERGY LEVELS: BOLTZMANN (LTE POPNUMBERS)
      DO IMOD=1, NMOD
        CALL MULTSPLI(ND(IMOD), N, NSUBLOW, MAXSUBL, NLOW, POPNUM(1,1,IMOD), 
     >                ELEVEL, 
     >                T(1,IMOD), WEIGHT, NNUP, NSUBNUP, LOW, NUP, 
     >                LEVEL)
      ENDDO


C***  Copy the constant AVOIGT values to second_model 
C***  (all lines, entry at L=1)
      DO IMOD=2, NMOD
         DO NBL = 1, NBLINE
          AVOIGT(NBL,1,IMOD) = AVOIGT(NBL,1,1)
         ENDDO
      ENDDO

   99 CONTINUE
      RETURN
      END
      SUBROUTINE MULTSPLI(ND, N, NSUBLOW, MAXSUBL, NLOW, POPNUM, 
     >                 ELEVEL, T, WEIGHT, NNUP, NSUBNUP, LOW, NUP, 
     >                 LEVEL)

      DIMENSION NSUBLOW(MAXSUBL), NSUBNUP(MAXSUBL)
      DIMENSION POPNUM(ND, N), T(ND), ELEVEL(N), WEIGHT(N)
      CHARACTER*10 LEVEL(N) 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
      NWEIGHTLOWREST = 0
      NWEIGHTNUPREST = 0      
C***  A. SPLITTING THE LOWER ENERGY LEVEL:
      IF (NLOW .GT. 0) THEN
C***     Check if sum of sublevel weights conforms with weight of original level
         WEIGHTSUM = WEIGHT(NSUBLOW(1))
          ELOWREST = WEIGHT(NSUBLOW(1)) * ELEVEL(NSUBLOW(1))
         DO J=2, NLOW
            WEIGHTSUM = WEIGHTSUM + WEIGHT(NSUBLOW(J))
            ELOWREST = ELOWREST + WEIGHT(NSUBLOW(J))*ELEVEL(NSUBLOW(J))
         ENDDO
         IF (NINT(WEIGHTSUM) .NE. NINT(WEIGHT(LOW))) THEN
            WRITE (*,'(A)') '*** WARNING: Inconsistency of stat. weights'
            WRITE (*,'(A, I3, A, I3)') 'LOWERLEVEL: ' // LEVEL(LOW)
     >           // '    WEIGHT: ', NINT(WEIGHT(LOW)),  
     >              '  =!= Sum of sublevel weights: ', NINT(WEIGHTSUM)
            WRITE (0,'(A)') 'WARNING: Inconsistency of WEIGHTS '
     >           // 'for level ' // LEVEL(LOW) // '  -- see output file'
            ! construct rest level:
            NWEIGHTLOWREST = NINT(WEIGHT(LOW) - WEIGHTSUM)
            ELOWREST=(ELEVEL(LOW) *WEIGHT(LOW)-ELOWREST)/NWEIGHTLOWREST
         ENDIF  

         DO 30 L=1,ND
         POPNUM(L,NSUBLOW(1))=1.
         DO 35 J=2,NLOW
   35    POPNUM(L,NSUBLOW(J))=EXP(C1*(ELEVEL(NSUBLOW(J-1))
     -         -ELEVEL(NSUBLOW(J)))/T(L))*WEIGHT(NSUBLOW(J))
     /         /WEIGHT(NSUBLOW(J-1))*POPNUM(L,NSUBLOW(J-1))
C***  NORMALIZATION
         SUM=0.
         IF (NWEIGHTLOWREST .GT. 0) THEN
            SUM=EXP(C1*(ELEVEL(NSUBLOW(1))-ELOWREST)/T(L))
     >           *NWEIGHTLOWREST/WEIGHT(NSUBLOW(1))*1.
            IF (L .EQ. 1) THEN
               WRITE (*,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
               WRITE (0,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
            ENDIF
         ENDIF
         DO 36 J=1,NLOW     
   36    SUM=SUM+POPNUM(L,NSUBLOW(J))
         SUMINV = POPNUM(L,LOW) / SUM
         DO 37 J=1,NLOW
   37    POPNUM(L,NSUBLOW(J))=POPNUM(L,NSUBLOW(J)) * SUMINV
   30    CONTINUE
      ENDIF

C***  B. SPLITTING THE UPPER ENERGY LEVEL: BOLTZMANN (LTE)
      IF (NNUP .GT. 0) THEN
         WEIGHTSUM = WEIGHT(NSUBNUP(1))
          ENUPREST = WEIGHT(NSUBNUP(1)) * ELEVEL(NSUBNUP(1))         
         DO J=2, NNUP
            WEIGHTSUM = WEIGHTSUM + WEIGHT(NSUBNUP(J))
            ENUPREST = ENUPREST + WEIGHT(NSUBNUP(J))*ELEVEL(NSUBNUP(J))
         ENDDO
C***     Check if sum of sublevel weights conforms with weight of original level
         IF (NINT(WEIGHTSUM) .NE. NINT(WEIGHT(NUP))) THEN
            WRITE (*,'(A)') '*** WARNING: Inconsistency of stat. weights'
            WRITE (*,'(A, I3, A, I3)') 'UPPERLEVEL: ' // LEVEL(NUP)
     >           // '    WEIGHT: ', NINT(WEIGHT(NUP)),  
     >              '  =!= Sum of sublevel weights: ', NINT(WEIGHTSUM)
            WRITE (0,'(A)') 'WARNING: Inconsistency of WEIGHTS '
     >           // 'for level ' // LEVEL(NUP) // '  -- see output file'
            NWEIGHTNUPREST = NINT(WEIGHT(NUP) - WEIGHTSUM)
            ENUPREST = (ELEVEL(NUP)*WEIGHT(NUP)-ENUPREST)/NWEIGHTNUPREST
         ENDIF  
         DO 40 L=1,ND
         POPNUM(L,NSUBNUP(1))=1.
         DO 45 J=2,NNUP
   45    POPNUM(L,NSUBNUP(J))=EXP(C1*(ELEVEL(NSUBNUP(J-1))
     -         -ELEVEL(NSUBNUP(J)))/T(L))*WEIGHT(NSUBNUP(J))
     /         /WEIGHT(NSUBNUP(J-1))*POPNUM(L,NSUBNUP(J-1))
C***     NORMALIZATION
         SUM=0.
         IF (NWEIGHTNUPREST .GT. 0) THEN
            SUM=EXP(C1*(ELEVEL(NSUBNUP(1))-ENUPREST)/T(L))
     >           *NWEIGHTNUPREST/WEIGHT(NSUBNUP(1))*1.
            IF (L .EQ. 1) THEN 
               WRITE (*,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
               WRITE (0,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
            ENDIF
         ENDIF         
         DO 46 J=1,NNUP
   46    SUM=SUM+POPNUM(L,NSUBNUP(J))
         SUMINV = POPNUM(L,NUP) / SUM
         DO 47 J=1,NNUP
   47    POPNUM(L,NSUBNUP(J))=POPNUM(L,NSUBNUP(J)) * SUMINV
   40    CONTINUE
      ENDIF

      RETURN
      
C********** Error branches ***************************************************
      
      END
      SUBROUTINE MVMD (BX,B,C,JMAX,JMM,NP)
C***********************************************************************
C***  MATRIX (VOLL)  B  *  MATRIX (DIAGONAL)  C
C***  ERGEBNIS-MATRIX  BX(VOLL)
C*** AKTUELLES FORMAT BX(JMAX,JMM)=B(JMAX,JMAX)*C(JMAX,JMM)
C***  WOBEI DIE UEBERZAEHLIGEN ZEILEN DER DIAGONALMATRIX C VERSCHWINDEN
C***********************************************************************

      DIMENSION BX(NP,NP),B(NP,NP),C(NP)

      DO 1 K=1,JMM
      CK=C(K)
      DO 1 I=1,JMAX
    1 BX(I,K)=B(I,K)*CK

      RETURN
      END
      SUBROUTINE MVV (WX,B,W,JMAX,JMM,NP)
C***********************************************************************
C***  MATRIX (VOLL)  B  *  VEKTOR W
C***  ERGEBNIS-VEKTOR  WX
C***  AKTUELLES FORMAT  WX(JMAX) = B(JMAX,JMM) * W(JMM)
C***********************************************************************

      DIMENSION WX(NP),B(NP,NP),W(NP)

      DO 1 I=1,JMAX
      WXI = .0
      DO 2 K=1,JMM
    2 WXI=WXI+B(I,K)*W(K)
    1 WX(I)=WXI

      RETURN
      END
      SUBROUTINE NEWPOL2(Y,X,YY,XX,ND)

C **** INTERPOLATION DURCH GEMITTELTE NEWTONPOLYNOME 2. GRADES.
C ***
C ***   
C ***
      IMPLICIT NONE
      INTEGER ND
      REAL Y(ND), X(ND), YY, XX

      INTEGER I, IS
      REAL Z, C21, C22, C23, C31, C32

C *** find XX
      DO I = 1, ND
          IF( XX .lt. X(I) ) GO TO 1
      end do
1     CONTINUE

      IS = I
      IF( IS .EQ. 1 ) THEN		! no extrapolation
	  YY = Y(1)
      ELSE IF( IS .GT. ND ) THEN	! still no extrapolation
	  YY = Y(ND)
      ELSE IF( IS .EQ. 2 ) THEN		! second point has problems
          C22 = (Y(2)-Y(1))/ (X(2)-X(1))
          C23 = (Y(3)-Y(2))/ (X(3)-X(2))
          C32 = (C23-C22)/ (X(3)-X(1))
          Z = XX - X(1)
          YY = Y(1) + C22*Z + C32*Z* (XX-X(2))
      ELSE IF( IS .LT. ND ) THEN	! normal interpolation
          C21 = (Y(IS-1)-Y(IS-2))/ (X(IS-1)-X(IS-2))
          C22 = (Y(IS)-Y(IS-1))/ (X(IS)-X(IS-1))
          C23 = (Y(IS+1)-Y(IS))/ (X(IS+1)-X(IS))
          C31 = (C22-C21)/ (X(IS)-X(IS-2))
          C32 = (C23-C22)/ (X(IS+1)-X(IS-1))
          Z = XX - X(IS-1)
          YY = Y(IS-1) + C22*Z + 0.5* (C31+C32)*Z* (XX-X(IS))
      ELSE ! IS .EQ. ND 		! some trouble again
          C21 = (Y(ND-1)-Y(ND-2))/ (X(ND-1)-X(ND-2))
          C22 = (Y(ND)-Y(ND-1))/ (X(ND)-X(ND-1))
          C31 = (C22-C21)/ (X(ND)-X(ND-2))
          Z = XX - X(ND-1)
          YY = Y(ND-1) + C22*Z + C31*Z* (XX-X(ND))
      END IF

      END
      SUBROUTINE NOWIND (NOWIND_LINE, RCON, NATOM, ATMASS, ABXYZ, 
     >                   ND, RNE, T, RADIUS, VELO, VDOP, TAURCONT, 
     >                   ND_MERGED, RADIUS_MERGED, PGRID_MERGED,
     >                   NP_MERGED, ZGRID_MERGED, IERR)

C***********************************************************************
C*** "NOWIND" CARD allows to calculate the mergent spectrum 
C*** as if the (outer part of) the wind would not exist
C***     The wind is either cut due to VELO, RADIUS or TAU 
C***  In case of SECONDMODEL, the criteria refer to the the first model 
C***  Note that the _MERGED input parameters are over-written
C***********************************************************************

      IMPLICIT NONE
      CHARACTER NOWIND_LINE*(*), ACTPAR*20
      INTEGER NPAR, NATOM, ND, L, NA, IDX, IERR
      INTEGER ND_MERGED, NP_MERGED, NCORE, JP, JMAX, I
      INTEGER INDCUT, ISRCHFLT

      REAL RCUT, RCON, XMUL, VSONIC, VDOP, VCUT, ATMEAN, TAUCUT
      REAL RR, PJ, PJPJ
      REAL ATMASS(NATOM), ABXYZ(NATOM)
      REAL RNE(ND), T(ND), RADIUS(ND), VELO(ND), TAURCONT(ND)
      REAL RADIUS_MERGED(ND_MERGED), PGRID_MERGED(NP_MERGED)
      REAL ZGRID_MERGED(ND_MERGED*NP_MERGED)

      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)

C****************************************************************
C***  First part: defining RCUT based on structure of first MODEL 
C****************************************************************

      CALL SARGC (NOWIND_LINE, NPAR)

C***  No parameters (default): cutting at connection radius 
      IF (NPAR == 1) THEN                  
         RCUT = RCON
         GOTO 20
      ENDIF

      CALL SARGV(NOWIND_LINE, 2, ACTPAR)
      IF (ACTPAR == 'OFF') THEN
         NOWIND_LINE = 'NONE'
         RETURN

      ELSEIF (ACTPAR == 'RADIUS') THEN
      IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=91) RCUT
         GOTO 20
      ELSEIF (ACTPAR == 'VDOP') THEN
         VCUT = VDOP
         GOTO 10

      ELSEIF (ACTPAR == 'TAU') THEN
      IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=91) TAUCUT
         GOTO 25

C***  For compatibility with older syntax, allow keywords
C***      VDOP or SONIC without preceeding VEL0=
      ELSEIF (ACTPAR == 'VDOP') THEN
         VCUT = VDOP
         GOTO 10

      ELSEIF (ACTPAR == 'SONIC') THEN
C***     Need to calculate sonic velocity  
         ATMEAN = 0.
         DO NA=1, NATOM
            ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
         ENDDO
         DO L=1, ND
            XMUL = ATMEAN / (1. + RNE(L))
            VSONIC = SQRT(RGAS * T(L)/ XMUL) * 1.E-5
            IF (VELO(L) < VSONIC) EXIT
         ENDDO
         VCUT = VSONIC
         GOTO 10

      ELSEIF (ACTPAR == 'VELO') THEN
         IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         IF (ACTPAR == 'VDOP') THEN
            VCUT = VDOP
         ELSEIF (ACTPAR == 'SONIC') THEN
C***        Need to calculate sonic velocity  
            ATMEAN = 0.
            DO NA=1, NATOM
               ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
            ENDDO
            DO L=1, ND
               XMUL = ATMEAN / (1. + RNE(L))
               VSONIC = SQRT(RGAS * T(L)/ XMUL) * 1.E-5
               IF (VELO(L) < VSONIC) EXIT
            ENDDO
            VCUT = VSONIC

         ELSE
            READ (ACTPAR, '(F10.0)', ERR=91) VCUT
         ENDIF
         GOTO 10

C***  Assume that the second parameter is the cut velocity value
      ELSE
         READ (ACTPAR, '(F10.0)', ERR=91) VCUT
         GOTO 10
      ENDIF

C***  velocity VCUT was specified
   10 CONTINUE
C***  Test if VCUT is within valid range
      IF (VCUT .GE. VELO(1)) THEN
         WRITE (0,*) '**** ERROR: VCUT .GE. max. velocity'
         GOTO 99
      ENDIF 
      IF (VCUT .LE. VELO(ND)) THEN
         WRITE (0,*) '**** ERROR: VCUT .LE. min. velocity'
         GOTO 99
      ENDIF 
      CALL SPLINPO (RCUT,   VCUT, RADIUS,   VELO, ND) 
      CALL SPLINPO (TAUCUT, VCUT, TAURCONT, VELO, ND) 
      GOTO 30

C***  If RADIUS was specified:
   20 CONTINUE
C***  Test if RCUT is within valid range
      IF (RCUT .GE. RADIUS(1)) THEN
         WRITE (0,*) '**** ERROR: RCUT .GE. RMAX'
         GOTO 99
      ENDIF 
      IF (RCUT .LE. 1.) THEN
         WRITE (0,*) '**** ERROR: RCUT .LE. 1.'
         GOTO 99
      ENDIF 
      CALL SPLINPO (VCUT,   RCUT, VELO,     RADIUS, ND) 
      CALL SPLINPO (TAUCUT, RCUT, TAURCONT, RADIUS, ND) 
      GOTO 30

C***  If TAU was specified:
   25 CONTINUE
C***  Test if TAUCUT is within valid range
      IF (TAUCUT .LE. TAURCONT(1)) THEN
         WRITE (0,*) '**** ERROR: TAUCUT .LE. TAURCONT(1)'
         GOTO 99
      ENDIF 
      IF (TAUCUT .GE. TAURCONT(ND)) THEN
         WRITE (0,*) '**** ERROR: TAUCUT .GE. TAUMAX'
         GOTO 99
      ENDIF 
      CALL SPLINPO (VCUT, TAUCUT, VELO,     TAURCONT, ND) 
      CALL SPLINPO (RCUT, TAUCUT, RADIUS,   TAURCONT, ND) 
      GOTO 30

    
   30 Continue
      WRITE(0,80) NOWIND_LINE(:IDX(NOWIND_LINE)), VCUT, RCUT, TAUCUT
      WRITE(*,80) NOWIND_LINE(:IDX(NOWIND_LINE)), VCUT, RCUT, TAUCUT
   80 FORMAT (A,/,'Thus removing layers above:',/,
     >        '     VCUT   =', F9.3, ' km/s',/,
     >        '     RCUT   =', F9.3, ' Rstar',/,
     >        '     TAUCUT =', F9.3)

C****************************************************************
C***  Second part: usining RCUT to shorten RADIUS_MERGED 
C****************************************************************

      IF (RCUT .GE. RADIUS_MERGED(1)) THEN
         WRITE (0,*) '**** ERROR: RCUT .GE. RMAX of second model'
         GOTO 99
      ENDIF

C***  Now define the new RADIUS_MERGED vector
      NCORE = NP_MERGED - ND_MERGED
      RADIUS_MERGED(1) = RCUT
      INDCUT =  ISRCHFLT(ND_MERGED, RADIUS_MERGED, 1, RCUT) - 1
      DO L= INDCUT+1, ND_MERGED
         RADIUS_MERGED(L-INDCUT+1) = RADIUS_MERGED(L)
      ENDDO  
      WRITE (0,'(A,I4)') 
     >     'Number of radial points (total wind) ND=', ND_MERGED
      ND_MERGED = ND_MERGED - INDCUT + 1
      NP_MERGED = NP_MERGED - INDCUT + 1
      WRITE (0,'(A,I4)') 
     >     'Reduced by NOWIND to                 ND=', ND_MERGED

C***  Re-establish PGRID and ZGRID according to the new RADIUS vector
      DO L=1, ND_MERGED
         JP = NP_MERGED + 1  - L
         PGRID_MERGED(JP) = RADIUS_MERGED(L)
      ENDDO

      DO L = 1, ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED + 1 - L
        DO JP = 1, JMAX
          PJ = PGRID_MERGED(JP)
          PJPJ = PJ * PJ
          I = (JP-1)*ND_MERGED + L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I) = SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I) = .0
          ENDIF
        ENDDO
c        write (0,*) L, RADIUS_MERGED(L)
      ENDDO

c      do jp=1, np_merged 
c        write (0,*) jp, pgrid_MERGED(jp)
c      enddo

      RETURN

C***  Error branches
   90 WRITE (0,*) '*** ERROR: This option requieres a value'
      GOTO 99

   91 WRITE (0,*) '*** ERROR when decoding floating-point number'
      GOTO 99

   99 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) NOWIND_LINE(:IDX(NOWIND_LINE))
      WRITE (0,*) 'Therefore, this spectral range must be skipped!'
      IERR = 1
      RETURN

      END
      SUBROUTINE OBSFRAM (LTOT,CORE,XMAX,XMAXLIN,EMINT,CEMINT,
     >                   BCOREL,DBDRL,TAUMAX,PJPJ, 
     >                   ZFINE, OPAFINE,OPAFC,OPALFIN,
     >                   ETAFINE,ETAFC,ETALFIN,SFINE,CSFINE,RRAY,
     >                   OPARAY,OPALRAY,
     >                   ETARAY,ETACRAY,ETALRAY, ZRAY,XCMF,
     >                   MAXXN,NDADDIM,NDDIM,DELXLAP,NBLINE,IVERSION, 
     >                   LINPRO,AVOIGT,TAU,TAUC,DTAU,DTAUC,WTAU,WTAUC,
     >                   XCMFFINE, POROLENGTHFINE, MAXLAP, DXMAX,
     >                   BIRONLINES, OPAFE, ETAFE, NFLDIM, 
     >                   XCMFBLUE, XCMFRED, DXCMF, POROLENGTHRAY, 
     >                PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >                RADIUS, ND, DD_VDOPDU_RAY, NATOM, NBFIRST, NBLAST,
     >                IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >                DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE,
     >                GRIEMPAR, KODAT, VDOP, INDCUT, 
     >                ZINTER, NMOD, MAXMOD, KINDEX)

C***********************************************************************
C***  INTEGRATION OF THE EMERGENT INTENSITY IN THE OBSERVER'S FRAME
C***********************************************************************
  
      INTEGER, INTENT(IN) :: NATOM
      REAL, DIMENSION(NDADDIM, NATOM), INTENT(IN) :: DD_VDOPDU_RAY
      REAL, DIMENSION(MAXXN, NATOM) :: DD_VDOPDU_FINE, DD_VDOPDU_FINE_SQRD, 
     >                                 DD_VDOPDU_FINE_NORMFAC
      INTEGER, DIMENSION (MAXLAP) :: IND_ELLINE
      DIMENSION XCMF(NDADDIM)
      DIMENSION RRAY(NDADDIM), OPARAY(NDADDIM)
      DIMENSION ETARAY(NDADDIM), ETACRAY(NDADDIM)
      DIMENSION ZRAY(NDADDIM)
      REAL, DIMENSION(NDADDIM,NBLINE) :: OPALRAY, ETALRAY
      DIMENSION DELXLAP(NBLINE)
      
      DIMENSION ZFINE(MAXXN), OPAFINE(MAXXN)
      DIMENSION OPAFC(MAXXN), OPALFIN(MAXXN)
      DIMENSION ETAFINE(MAXXN), ETAFC(MAXXN)
      DIMENSION ETALFIN(MAXXN), SFINE(MAXXN)
      DIMENSION CSFINE(MAXXN)

      LOGICAL CORE, BIRONLINES, BTAUMAX_REACHED
      CHARACTER LINPRO(MAXLAP)*8
      REAL, DIMENSION(MAXLAP) :: XMAXLIN

C***  WPI = SQRT(PI)
      DATA WPI /1.772454/

C***  INITIALIZATION
      EMINT = .0
      CEMINT = .0
      TAUSUM = .0
      TAUSUMC = .0
      BTAUMAX_REACHED = .FALSE.

      LRED = 1
      LBLUE = LTOT

      CALL ZONEINT (LRED, LBLUE, EMINT, CEMINT, TAUSUM, TAUSUMC,
     >                PJPJ, ZFINE, KINDEX,
     >                OPAFINE, OPAFC, OPALFIN, ETAFINE, ETAFC,
     >                ETALFIN, SFINE, CSFINE,
     >                RRAY, ZRAY, 
     >                XCMF, OPARAY,
     >                OPALRAY, ETARAY,
     >                ETACRAY, ETALRAY,
     >                MAXXN, LTOT, DELXLAP, NBLINE, NBFIRST, NBLAST, 
     >                NDADDIM, IVERSION, LINPRO, AVOIGT,
     >                TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >                XCMFFINE, POROLENGTHFINE, TAUMAX, XMAX, XMAXLIN,
     >                DXMAX, BIRONLINES, OPAFE, ETAFE, NDDIM, NFLDIM, 
     >                XCMFBLUE, XCMFRED, DXCMF, CORE, POROLENGTHRAY, 
     >                PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >                RADIUS, ND, MAXLAP, DD_VDOPDU_RAY, NATOM, 
     >                IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >                DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE, 
     >                GRIEMPAR, KODAT, VDOP, ZINTER, NMOD, MAXMOD, 
     >                BTAUMAX_REACHED)
 
C***  FOR CORE RAYS, ADD INCIDENT RADIATION
      IF (CORE .AND. .NOT. BTAUMAX_REACHED) THEN
         X = OPARAY(LTOT)  ! X means kappa

C***     Without lines: PLUSIC
         PLUSIC = BCOREL + DBDRL * ZRAY(LTOT) / X

C***     With lines: add line opacities
         DO 45 NLOC = NBFIRST, NBLAST
            XLBLTOT=XCMF(LTOT)-DELXLAP(NLOC)
            IF (XLBLTOT .GT. 300) THEN
               PHI = .0
            ELSE
               PHI=EXP(-XLBLTOT*XLBLTOT)/WPI
            ENDIF            
            X=X+PHI*OPALRAY(LTOT,NLOC)
   45    CONTINUE
         PLUSI  = BCOREL + DBDRL * ZRAY(LTOT) / X
         EMINT  = EMINT  + PLUSI  * EXP(-TAUSUM)
         CEMINT = CEMINT + PLUSIC * EXP(-TAUSUMC)
      ENDIF
 
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
      SUBROUTINE OWNINV (N, NDIM, A, CKEY)
C*******************************************************************************
C***  MATRIX INVERSION (CRAY FORTRAN)
C***  MATRIX A
C***  N = RANK OF SUBMATRIX (LEFT UPPER BLOCK) TO BE INVERTED
C***  NDIM = ROW DIMENSION OF TOTAL MATRIX
C***  N .LE. NDIM .LE. NMAX
C***  Most time consuming routine in STEAL 
C***  in WRSTART-Job : (69 percent), (Loop 80: 66 percent)
C***  Most time consuming routine in WRCONT
C***                   (80 percent), (Loop 80: 76 percent)
C***  Tested on 30-Jan-1997 20:54:25, Lars
C***
C***  Now a Normalization of the rows and columns is possible (BNORM = .TRUE.)
C***  Row-Norm is Quadratic Norm
C***  Column-Norm is Maximum Norm
C***
C***
C*******************************************************************************

      LOGICAL BNORM

      PARAMETER (BNORM = .TRUE.)
C      PARAMETER (BNORM = .FALSE.)
C      PARAMETER (NMAX = 560)  old version before split
      PARAMETER (NMAX = 2002)

      DIMENSION A(NDIM,NDIM),IK(NMAX),JK(NMAX)
      DIMENSION CNORM(NMAX), RNORM(NMAX)
      CHARACTER*4 CKEY

C***  OUTPUT of matrix dimension parameters (testing only)
C      WRITE (0,*) 'OWNINV> N=', N, 'NDIM=', NDIM

      IF (N .GT. NDIM) THEN
            CALL REMARK ('N .GT. NDIM')
            STOP 'ERROR IN OWNINV'
            ENDIF
      IF (NDIM .GT. NMAX) THEN
            CALL REMARK ('NDIM .GT. NMAX')
            WRITE (0,'(A,i4,1x,a,i4)') 'NDIM=',NDIM, 'NMAX=',NMAX
            STOP 'ERROR IN OWNINV'
            ENDIF

C*** New Branch to normalize Columns and Rows
C*** Calculation of the Row-Norm (Quadratic Norm)
      IF (BNORM) THEN
        DO K=1, N
          R = 0.
          DO L=1, N
            XM = A(L,K)
            IF (ABS(XM) .GT. 1.E100) THEN
              R = 1.
              EXIT
            ENDIF
            R = R + XM*XM
          ENDDO
          RNORM(K) = 1. / SQRT(R)
        ENDDO

C*** Calculation of the Column-Norm (Maximum Norm)
        DO L=1, N
          C = A(L,1)
          DO K=2, N
            XM = ABS(A(L,K))
            IF (XM .GT. C) C = XM
          ENDDO
          CNORM(L) = 1. / C
        ENDDO

C*** Applikation of CNORM and RNORM
        DO K=1, N
          DO L=1, N
            A(L,K) = A(L,K) * CNORM(L) * RNORM(K)
          ENDDO
        ENDDO
      ENDIF

      DO 100  K=1,N
C
C     SUCHE MAXIMALES MATRIXELEMENT
      AMAX=0.
      ABSAMAX=0.
      DO 30 J=K,N
        L=N+1-K
        IMAX=ISAMAX(L,A(K,J),1)+K-1
        IF(ABSAMAX.GT.ABS(A(IMAX,J))) GOTO 30
        AMAX=A(IMAX,J)
        ABSAMAX=ABS(AMAX)
        IK(K)=IMAX
        JK(K)=J
   30 CONTINUE
C
C***  WARNING IN CASE OF SINGULARITY (PIVOT ELEMENT = 0 )
C***  Exit the Routine with CKEY = 'SING'
      IF (AMAX .EQ. .0) THEN
        CALL REMARK ('Subr. OWNINV: Singularity discovered')
        IF (CKEY .EQ. 'OWNL') THEN
          CKEY = 'SING'
          RETURN
        ELSE
          STOP 'ERROR in Subr. OWNINV'
        ENDIF
      ENDIF
C
C     VERTAUSCHEN DER ZEILEN I UND K
      I=IK(K)
      IF(I.EQ.K) GOTO 51
      DO 50 J=1,N
      SAVE =A(K,J)
      A(K,J)=A(I,J)
   50 A(I,J)=-SAVE
C
C     VERTAUSCHEN DER SPALTEN J UND K
   51 J=JK(K)
      IF(J.EQ.K) GOTO 61
      DO 60 I=1,N
      SAVE = A(I,K)
      A(I,K)=A(I,J)
   60 A(I,J)=-SAVE
C
C     DIVISION DER SPALTE K DURCH AMAX
   61 DO 70  I=1,N
      A(I,K)=-A(I,K)/AMAX
   70 CONTINUE
      A(K,K)=0.
C
C     UMFORMEN DER ZEILEN UND SPALTEN
C***  ELEMENTARE UMFORMUNG:  ZEILE I = ZEILE I + A(I,K) * ZEILE K
      DO 81 I=1,N
        IF (I .EQ. K) GOTO 81
        AIK=A(I,K)
        DO 80 J=1,N
          A(I,J)=A(I,J)+AIK*A(K,J)
   80 CONTINUE
C!!!  Tried to replace Loop 80
C!!!  It works, but it is slower!
c      CALL DGEMA ('N', 'N', 1, N, 1., A(I,1), NDIM, AIK, A(K,1), NDIM, 
c     >            A(I,1), NDIM)
   81 CONTINUE
C!!!  Tried to replace Loops 80 and 81 by
C!!!  It works, but it is also slower!
c      FORALL (I=1:N, J=1:N, I/=K) A(I,J) = A(I,J) + A(I,K)*A(K,J)


C
C***  SPALTE K: DIVISION DUCH AMAX
      DO 90 J=1,N
      A(K,J)=A(K,J)/AMAX
  90  CONTINUE
C
C***  DIAGONALELEMENT
      A(K,K)=1./AMAX
  100 CONTINUE
C
C
C     ZEILEN UND SPALTEN RUECKTAUSCHOPERATIONEN
      DO 130 L=1,N
      K=N-L+1
      J=IK(K)
      IF(J.LE.K) GOTO 111
C***  VERTAUSCHEN DER SPALTEN J UND K
      DO 110 I=1,N
       SAVE=A(I,K)
      A(I,K)=-A(I,J)
  110 A(I,J)=SAVE
  111 I=JK(K)
      IF(I .LE. K) GOTO 130
C***  VERTAUSCHEN DER ZEILEN I UND K
      DO 120 J=1,N
      SAVE=A(K,J)
      A(K,J)=-A(I,J)
  120 A(I,J)=SAVE
  130 CONTINUE
C

C*** Applikation of RNORM and CNORM
      IF (BNORM) THEN
        DO K=1, N
          DO L=1, N
            A(L,K) = A(L,K) * RNORM(L) * CNORM(K)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
      FUNCTION PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
C***********************************************************************
C***  called from STARKHOLTZMARK
C***
C***  The profile is obtained by connecting a doppler core to an
C***   asymptotic Holtsmark function
C***
C***  BETA = DELTALAM / F_0 / K_low,up           
C***  BETADOP = c * DELTANUED / F_0 / K_low,up
C***
C***  Attention: BETA and BETADOP are not truly dimensionless, but
C***             instead need to be provided in Angstroem/cm
C***
C***
C***  concept taken from Hubeny et al. 1994, A&A 282, 151 (Appendix B)
C***********************************************************************

      IMPLICIT NONE
      
      LOGICAL, INTENT(IN) :: bHYDROGEN

      REAL :: PHIHOLTSMARK
      REAL :: BETADCRIT, FHOLTSMARK, X0, X1, DX, C, BETASTAR
      REAL, INTENT(INOUT) :: BETA
      REAL, INTENT(IN) :: BETADOP
      
      REAL, PARAMETER :: XEPS = 1.E-3
      
      REAL, PARAMETER :: BETADCRITH  = 5.82         !for Hydrogen only      
c      REAL, PARAMETER :: BETADCRITHE = 3.66         !for all other hydrogenic ions
      REAL, PARAMETER :: BETADCRITHE = 3.67         !for all other hydrogenic ions
      
C***  Edmonds fitting parameters      
      REAL, PARAMETER :: A0 = 0.07209481
      REAL, PARAMETER :: A1 = 0.4796232
      REAL, PARAMETER :: A2 = -0.5758228        !Attention: The minus sign is erronously missing in Hubeny et al (1994)

C***  Stark broadening: empirical branch limits
C      Original set from Hubeny et al. (1994)
c      REAL, PARAMETER :: BL1 = 1.14
c      REAL, PARAMETER :: BL2 = 11.4
C      Revised set from Hubeny's Synspec code
C       better connection between branches
      REAL, PARAMETER :: BL1 = 1.52
      REAL, PARAMETER :: BL2 = 8.325        
      
      REAL, PARAMETER :: WPI         = 1.77245385   !sqrt(Pi)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      
C***  Ensure that beta is always a positive wavelength difference      
      BETA = ABS(BETA)      
      IF (BETADOP < 0.) STOP 'FATAL ERROR: BETADOP < 0'
      
      IF (bHYDROGEN) THEN
        BETADCRIT = BETADCRITH
        FHOLTSMARK = 1.
      ELSE
        BETADCRIT = BETADCRITHE
        FHOLTSMARK = 1./2.
      ENDIF
           
           
C**********************************************************************            
      IF (BETADOP >= BETADCRIT) THEN     
C***    Doppler broadening is dominant in the line core
C***       => join Doppler profile with Holtsmark profile in the wing
C***          switch is done at BETA > BETASTAR
C***            with BETASTAR denoting the position where both profile
C***            functions have the same value. With x = BETA/BETADOP
C***            this leads to the Equation x^2 - 2.5 * ln(x) - C = 0

C***    Calculation of the constant C
        C = 1.5 * LOG(BETADOP) - LOG(3.*FHOLTSMARK*WPI)

C***    BETASTAR is determined via a short iteration using a
C***    Newton-Raphson approach to find x where
C***       x^2 - 2.5 * ln(x) - C = 0
        
C***    N-R starting value is chosen depending on C (taken from Hubeny)        
        IF (C > 1.26) THEN
          X0 = SQRT(C) * ( 1. + 1.25 * LOG(C) / (4.*C - 5) )
        ELSE
          X0 = SQRT(C + 0.28)
        ENDIF

        DO
          X1 = X0 - (X0*X0 - 2.5*LOG(X0) - C) / (2.*X0 - 2.5/X0)
          DX = 1. - X0/X1
          IF (ABS(DX) < XEPS) THEN
C***        iteration converged          
            BETASTAR = X1 * BETADOP
            EXIT
          ELSE
            X0 = X1
            IF (X0 <= 0) THEN
              WRITE (hCPR,*) 'PROBLEM WITH LINEAR STARK BROADENING'
              WRITE (hCPR,*) 'CANNOT DETERMINE CORE-WIND CONNECTION'
              WRITE (hCPR,*) 'FATAL (X0 < 0),  X0 = ', X0
              WRITE (hCPR,*) 'BETADOP > BETACRIT', BETADOP, BETADCRIT
              WRITE (hCPR,*) 'C, bHYDROGEN = ', C, bHYDROGEN
              STOP 'FATAL ERROR IN PHIHOLTSMARK'
            ENDIF 
          ENDIF
        ENDDO
C***    -- end of N-R iteration --        
        
        
        IF (BETA <= BETASTAR) THEN
C***      Doppler profile in the core
          PHIHOLTSMARK = 1./WPI/BETADOP * EXP( -1.* (BETA/BETADOP)**2 )
        ELSE
C***      Holtsmark profile in the wing
          PHIHOLTSMARK = 3. * FHOLTSMARK * BETA**(-5./2.)
        ENDIF

        
C**********************************************************************      
      ELSE
C***    BETADOP < BETADCRIT => Stark broadening is dominant everywhere
C                               (Doppler broadening is neglected)
c        WRITE (hCPR,*) 'Dominant Stark broadening'

C***    Fitting according to calculations of Edmonds et al. (1967)
        IF (BETA < BL1) THEN
          PHIHOLTSMARK = 0.08
        ELSEIF (BETA < BL2) THEN
          PHIHOLTSMARK = A0 * EXP(A1 * LOG(BETA) + A2 * (LOG(BETA))**2 )
        ELSE
          PHIHOLTSMARK = 3. * BETA**(-5./2.)
        ENDIF
C***    PHI needs to be scaled with FHOLTSMARK-factor (different for Hydrogen)        
        PHIHOLTSMARK = PHIHOLTSMARK * FHOLTSMARK

      ENDIF
C**********************************************************************      

C***  Profile normalization
C***  (except for depth-dependent VDOP, similar to VOIGTH)
      PHIHOLTSMARK = BETADOP * PHIHOLTSMARK 

      
C***  Fatal error for fatal results
C***  (indicating that something must be very wrong in the calling parameters)      
      IF (PHIHOLTSMARK < 0) STOP 'FATAL ERROR: PHIHOLTSMARK < 0'
      
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
      SUBROUTINE PLOT_SECONDMODEL_GRID (P, NP, NPDIM, NPHI, PHIARR,
     >           NPHIMAX, JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG, 
     >           ZINTER)
C******************************************************************
C**   Plot of the ray positions, seen as facing the stellar disk 
C******************************************************************
      DIMENSION P(NPDIM), NPHI(NPDIM) 
      DIMENSION PHIARR(NPHIMAX,NPDIM)
      PARAMETER ( MAXPLOT = 1000 )
      DIMENSION XPLOT(MAXPLOT), YPLOT(MAXPLOT)
      DIMENSION ZINTER (2,NPDIM, NPHIMAX)

      WRITE (66,*) 'PLOT: Grid for second-model geometry'
      WRITE (66,*) '\NOBOX'
      WRITE (66,*) '\OFS 4 2'
      IF (2*NP .GT. 99) WRITE (66,'(A,I4)') '\SET_NSETMAX ',  2*NP+1  
    
      SCALE = 9./ P(JPLAST)
      DO JP=1, JPLAST
         WRITE (66,'(A,F10.4,A)') '\ARC 0 0 0 0 ', SCALE*P(JP), ' 0 360'
         WRITE (66,'(A,F10.4,A, I3)') 
     >             '\LUN 0 0 ', SCALE*P(JP), ' 1 .2 ', JP
      ENDDO

      XPLOT(1)=.0
      YPLOT(1)=.0

      CALL PLOTANFS (66, '', '', '', '', 
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > XPLOT, YPLOT, 1, 'SYMBOL=8 SIZE=-0.05 COLOR=1')

      DO JP=JPFIRST, JPLAST
         IF (NPHI(JP) .GT. MAXPLOT) THEN
            WRITE (*,*) 'ERROR: INSUFFICIENT DIMENSION' 
            STOP '*** ERROR in subr. PLOT_SECONDMODEL_GRID'
         ENDIF
         LPHISTA = MAX(1,       LPHISTA_ORIG)
         LPHISTA = MIN(LPHISTA, NPHI(JP))
         LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
         LPHIEND = MAX(LPHIEND, LPHISTA)

C***     Red dots: intersecting with SECONDMODEL domain
         NCOUNT = 0
         DO LPHI= LPHISTA, LPHIEND
            IF (ZINTER(1,JP,LPHI) .NE. ZINTER(2,JP,LPHI)) THEN
               NCOUNT = NCOUNT + 1 
               XPLOT(NCOUNT) = P(JP) * COS(PHIARR(LPHI,JP))
               YPLOT(NCOUNT) = P(JP) * SIN(PHIARR(LPHI,JP))
            ENDIF
         ENDDO
         IF (NCOUNT .GT. 0) CALL PLOTCONS (66, XPLOT, YPLOT, 
     >            NCOUNT, 'SYMBOL=8 SIZE=-0.1 COLOR=2')

C***     Blue dots: not intersecting with SECONDMODEL domain
         NCOUNT = 0
         DO LPHI= LPHISTA, LPHIEND
            IF (ZINTER(1,JP,LPHI) .EQ. ZINTER(2,JP,LPHI)) THEN
               NCOUNT = NCOUNT + 1 
               XPLOT(NCOUNT) = P(JP) * COS(PHIARR(LPHI,JP))
               YPLOT(NCOUNT) = P(JP) * SIN(PHIARR(LPHI,JP))
            ENDIF
         ENDDO

         IF (NCOUNT .GT. 0) CALL PLOTCONS (66, XPLOT, YPLOT, 
     >            NCOUNT, 'SYMBOL=8 SIZE=-0.1 COLOR=4')

      ENDDO

      RETURN
      END
      SUBROUTINE PLOTVDOP (RADIUS, TAUROSS, VELO, DD_VDOP, ND, NATOM, 
     >                     BPLOTVDOP, SYMBOL, BMICROTURB, DD_VMIC,
     >                     DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD, 
     >                     BIRONLINES)
C******************************************************************************
C***  Plot of the Dopper-broadening velocity versus v(r) or tauross
C***  Called from: FORMAL
C******************************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NATOM
      REAL, DIMENSION(ND,NATOM) :: DD_VDOP 
      REAL, DIMENSION(ND) :: RADIUS, DD_VMIC, TAUROSS, VELO
      REAL, DIMENSION(ND) :: X, Y
      CHARACTER*60 HEAD1, HEAD2, STYLE
      CHARACTER*(*) DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD
      CHARACTER*50 ACTPAR1, ACTPAR2, XAXSTR, YAXSTR
      CHARACTER*2 SYMBOL(NATOM)
      INTEGER :: KANAL, L, NA, NPAR1, NPAR2
      REAL :: XMIN, XMAX, YMIN, YMAX, VCON,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
     
      LOGICAL :: BPLOTVDOP, BMICROTURB, BIRONLINES
 
C***  INITIALIZATION
      KANAL=1
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
 
C***  HEADER  ------------------------------------------------------
      HEAD1='VDOP'
      HEAD2= 'M' //  MODHEAD(12:32) // ' Depth-dependent VDOP'

! C***  X-AXIS: ----------------------------------------------------------
! C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
!       XMAX=2.9
!       XMIN=-3.0
!       XSCALE = 22./(XMAX-XMIN)
!       XTICK=0.5
!       XABST=1.
!       
!       
! C***  Y-AXIS:  ---------------------------------------------------------
! C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
!       YMAX = MAXVAL(DD_VDOP)
!       YMAX = YMAX + 0.05*YMAX
!       YMIN = .0
!       YSCALE = 15./(YMAX-YMIN)
!       YTICK=2.5
!       YABST=5.
!         
C*** First data set
      CALL SARGC(VDOPPLOT_LINE, NPAR1) 
      IF (NPAR1 .EQ. 3) THEN
        CALL SARGV (VDOPPLOT_LINE, 3, ACTPAR1)
      ELSE
        ACTPAR1 = ''
      ENDIF
      CALL SARGC(DD_VDOP_LINE, NPAR2)
      IF (NPAR2 .GE. 5) THEN
        CALL SARGV (DD_VDOP_LINE, 5, ACTPAR2)
      ELSE
        ACTPAR2 = ''
      ENDIF
      IF (ACTPAR1 .EQ. '') THEN
        IF ((ACTPAR2 .EQ. '') .OR. (ACTPAR2(:4) .EQ. 'VELO')) THEN
            DO L=1, ND-1
                X(L) = VELO(L)   
            ENDDO
            XAXSTR = '\CENTER\Wind velocity [km/s]'
        ELSE IF (ACTPAR2(:3) .EQ. 'TAU') THEN
            DO L=1, ND-1
                X(L) = TAUROSS(L)
            ENDDO
            XAXSTR = '\CENTER\#t#'
        ELSE IF (ACTPAR2(:1) .EQ. 'R') THEN
            DO L=1, ND-1
            X(L) = ALOG10(RADIUS(L)-1.) 
            ENDDO
            XAXSTR = '\CENTER\log (r/R&T*&M - 1)'
        ELSE 
            WRITE(0,*) "VDOP Plot version not known"
            STOP "Fatal error in subroutine PLOTVDOP"
        ENDIF
      ELSE IF (ACTPAR1 .EQ. 'R') THEN
        DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.) 
        ENDDO
        XAXSTR = '\CENTER\log (r/R&T*&M - 1)'
      ELSE IF (ACTPAR1 .EQ. 'VELO') THEN
        DO L=1, ND-1
            X(L) = VELO(L)   
        ENDDO
        XAXSTR = '\CENTER\Wind velocity [km/s]'
      ELSE IF (ACTPAR1 .EQ. 'TAU') THEN
        DO L=1, ND-1
            X(L) = TAUROSS(L)
        ENDDO
        XAXSTR = '\CENTER\#t#'
      ELSE
        WRITE(0,*) "VDOP Plot version not known"
        STOP "Fatal error in subroutine PLOTVDOP"
      ENDIF
      DO L=1, ND-1
        Y(L) = DD_VDOP(L,1)
      ENDDO
      YAXSTR = '\CENTER\Doppler broadening velocity VDOP [km/s]'

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX=MAXVAL(X(:ND-1))
      XMIN=MINVAL(X(:ND-1))
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
      YMAX = MAXVAL(Y(:ND-1))
      YMIN = 0.
      YSCALE = 15./(YMAX-YMIN)
      YTICK=2.5
      YABST=5.
      
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'
      WRITE (KANAL, '(A)') '\FONT=HELVET'
      WRITE (KANAL, '(A)') '\PENDEF=3'
      WRITE (KANAL, '(A)') '\DEFINECOLOR 3 = .0 .5 .0'

      
C*** Legend
      IF (BMICROTURB) THEN
       WRITE (KANAL, '(A)') '\COLOR=1'
       WRITE (KANAL, '(A)') "\LINREL XMAX YMAX  1 0 0.5 -1"
       WRITE (KANAL, '(A, A)') "\LUN XMAX YMAX L1.7 M-1 0.3 ", SYMBOL(1)
       DO NA=2, NATOM
         WRITE (KANAL, '(A, I2)') '\COLOR=',NA 
         WRITE (KANAL, '(A, F10.5)') 
     >         "\LINREL XMAX YMAX  1 0 0.5", -1.-0.6*(NA-1)
         WRITE (KANAL, '(A, A)') "\NEXTLUN ", SYMBOL(NA)
       ENDDO
      
       WRITE (KANAL, '(A)') '\COLOR=1'
       WRITE (KANAL, '(A, F10.5, A)') "\LINREL XMAX YMAX  1 0 0.5", 
     >                      -1.-0.6*NATOM, " SYMBOL=9 SIZE=0.07"
        WRITE (KANAL, '(A)') "\NEXTLUN v&Tmic&M"
      ELSE
        IF (BIRONLINES) THEN
         WRITE (KANAL, '(A,I3)') '\COLOR=', NATOM
         WRITE (KANAL, '(A)') "\LINREL XMAX YMAX  1 0 0.5 -1"
         WRITE (KANAL, '(A, A)') "\LUN XMAX YMAX L1.7 M-1 0.3 VDOP-FEDAT"
        ENDIF
      ENDIF

      WRITE (KANAL, '(A)') '\COLOR=1'
      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,XAXSTR
     >        ,YAXSTR
     >        ,XSCALE,0.,0.,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
      DO NA=2, NATOM
        IF (.NOT. BMICROTURB .AND. .NOT. SYMBOL(NA) .EQ. 'G ') CYCLE
        DO L=1, ND-1
            Y(L) = DD_VDOP(L,NA)
        ENDDO
        STYLE = 'SYMBOL=5 COLOR='
        WRITE(STYLE(15:16), '(I2)') NA
        CALL PLOTCONS (KANAL,X,Y,ND-1, STYLE) 
      ENDDO
      IF (BMICROTURB) THEN
        DO L=1, ND-1
            Y(L) = DD_VMIC(L)
        ENDDO
        STYLE = 'SYMBOL=9 COLOR=1 SIZE=0.07'
        CALL PLOTCONS (KANAL,X,Y,ND-1, STYLE) 
      ENDIF

C***  Set back to initial values (for next blend-range):
      BPLOTVDOP = .FALSE.
      
      RETURN
      END
      SUBROUTINE PLOT_WINDROT_GRID (P, NPDIM, JPFIRST, JPLAST, 
     >      LPHISTA_ORIG, LPHIEND_ORIG, NPHI, 
     >      PHIARR, NPHIMAX, XPLOT, YPLOT)
C******************************************************************
C**   Plot of the ray positions, seen as facing the stellar disk 
C******************************************************************
      DIMENSION P(NPDIM), NPHI(NPDIM) 
      DIMENSION PHIARR(NPHIMAX,NPDIM)
      DIMENSION XPLOT(NPHIMAX), YPLOT(NPHIMAX)

      WRITE (65,*) 'PLOT: Grid for wind rotation'
      WRITE (65,*) '\NOBOX'
      NSET = 2 + JPLAST - JPFIRST
      IF (NSET .GT. 99) WRITE (65,'(A,I4)') '\SET_NSETMAX ', NSET  
    
      RMAX = 8.
      SCALE = 12./P(JPLAST)
      DO JP=JPFIRST, JPLAST
         WRITE (65,'(A,F10.4,A)') '\ARC 0 0 0 0 ', SCALE*P(JP), ' 0 180'
      ENDDO

      XPLOT(1)=.0
      YPLOT(1)=.0
      CALL PLOTANFS (65, '', '', '', '', 
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > SCALE,     .0, P(JPLAST), 1., 10., 0.,
     > XPLOT, YPLOT, 1, 'SYMBOL=8 SIZE=-0.05 COLOR=1')

      DO JP=JPFIRST, JPLAST
C***     If specified, the plot is restricted to the points 
C***     between LPHISTA and LPHIEND
         LPHISTA = MAX(1,       LPHISTA_ORIG)
         LPHISTA = MIN(LPHISTA, NPHI(JP))
         LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
         LPHIEND = MAX(LPHIEND, LPHISTA)

         DO LPHI=LPHISTA, LPHIEND
            XPLOT(LPHI) = P(JP) * COS(PHIARR(LPHI,JP))
            YPLOT(LPHI) = P(JP) * SIN(PHIARR(LPHI,JP))
         ENDDO

         CALL PLOTCONS (65, XPLOT(LPHISTA), YPLOT(LPHISTA), 
     >            1+LPHIEND-LPHISTA, 'SYMBOL=8 SIZE=-0.1 COLOR=2')
      ENDDO

      RETURN
      END
      SUBROUTINE POLYFIT (XFIT, YFIT, WFIT, M, KPLUS1, 
     >                    X, A, B, SCRATCH, ATEST, BTEST, DTEST)
C******************************************************************************
C***  THIS SUBROUTINE CALCULATES THE POLYNOM KOEFFICIENTS FOR AN
C***  POLYNOMINAL INTERPOLATION OF RANK (KPLUS1-1) FOR M POINTS
C***  GIVEN BY XFIT AND YFIT AND WEIGHTED WITH WFIT.
C***  SCRATCH, ATEST, BTEST AND DTEST ARE USED FOR THE ROUTINE LINSOL
C******************************************************************************

      DIMENSION XFIT(M),YFIT(M),WFIT(M)
      DIMENSION X(KPLUS1), A(KPLUS1,KPLUS1), B(KPLUS1)
      DIMENSION SCRATCH(2*KPLUS1), ATEST(KPLUS1,KPLUS1)
      DIMENSION BTEST(KPLUS1), DTEST(KPLUS1)
      CHARACTER*4 CKEY

      CKEY = 'OWN'

      IF (M .LT. KPLUS1) THEN
         CALL REMARK ('KPLUS1 EXCEEDS M')
         STOP 'ERROR IN POLYFIT'
         ENDIF

C***  NEW BRANCH TO AVOID NAG-CALL E02ADF
C***  Probably (!) this is only made for cubic degree
      IF (KPLUS1 .NE. 4) THEN
         WRITE (0,*) 'INVALID DEGREE OF POLYNOMIAL FIT: .NE. 4'  
         STOP '****** ERROR IN POLYFIT **********'
      ENDIF

      DO I=0, 6
        SCRATCH(I+1) = 0.
      ENDDO
      DO I=0, 3
        DTEST(I+1) = 0.
      ENDDO
      DO J=1, M
        XF = 1.
        DO I=0, 6
          SCRATCH(I+1) = SCRATCH(I+1) + WFIT(J) * XF
          IF (I .LT. KPLUS1) DTEST(I+1) = DTEST(I+1) + WFIT(J) * 
     >                       XF * YFIT(J)
          XF = XF * XFIT(J)
        ENDDO
      ENDDO

      DO I=1, KPLUS1
        B(I) = DTEST(KPLUS1+1-I)
      ENDDO

      A(1,1) = SCRATCH(7)
      A(2,1) = SCRATCH(6)
      A(1,2) = SCRATCH(6)
      A(3,1) = SCRATCH(5)
      A(2,2) = SCRATCH(5)
      A(1,3) = SCRATCH(5)
      A(4,1) = SCRATCH(4)
      A(3,2) = SCRATCH(4)
      A(2,3) = SCRATCH(4)
      A(1,4) = SCRATCH(4)
      A(4,2) = SCRATCH(3)
      A(3,3) = SCRATCH(3)
      A(2,4) = SCRATCH(3)
      A(4,3) = SCRATCH(2)
      A(3,4) = SCRATCH(2)
      A(4,4) = SCRATCH(1)

      CALL INV (KPLUS1, KPLUS1, A, CKEY)
      CALL VMF (X, B, A, KPLUS1, KPLUS1)

C!!!      CALL LINSOL (X, A, B, KPLUS1, KPLUS1, 
C!!!     >             SCRATCH, ATEST, BTEST, DTEST)

      RETURN
      END
      SUBROUTINE POPMIN_NULLING (ZERO_RATES, POPNUM, POPMIN, ND, N)
C***********************************************************************
C***  Sets all POPNUMS that were flagged by ZERO_RATES to zero 
C***  Reason: these level populations have been set to POPMIN in steal 
C***      but should be set to 0.0 in all radiative-transfer programs
C***      in oder to avoid any artificial contributions to the 
C***      emissivities and opacities 
C***  Called from: WRCONT, COMO, COLI, FORMAL  
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND, N
      REAL, INTENT(IN) :: POPMIN

      REAL,    DIMENSION(ND, N) :: POPNUM
      LOGICAL, DIMENSION(N, ND) :: ZERO_RATES
    
      INTEGER :: L, J

      DO L=1, ND
        DO J=1, N
          IF (ZERO_RATES(J,L) .OR. POPNUM(L,J) < 1.1 * POPMIN) THEN
            POPNUM(L,J) = .0
          ENDIF
        ENDDO        
      ENDDO

      RETURN
      END
      SUBROUTINE PREFORM (KARTE,N,ELEVEL,LINE,INDLOW,INDNUP,LASTIND,
     $                   CLIGHT, VDOP, INDLAP, XLAMLAP, DELXLAP, ALN,
     $                   XLAM, NBLINE, MAXLAP, MAXIND, MAXATOM, 
     $                   LEVEL, WEIGHT, EINST, NDIM, POPNUM,
     >                   T, ND, NOM, NCHARG, EION, ENTOT, RNE, MAXSUBL,
     $                   NSUBLOW, NSUBNUP, BROAD, LINPRO,
     >                   AVOIGT, NMOD, NDDIM, MAXMOD, DENSCON, MAINQN,
     >                   MULTIIND, NMULTI, DD_VDOP, NATOM, IND_ORIGLEV)
C*******************************************************************************
C***  This subroutine reads atomic data from FORMAL_CARDS   
C***  CALLED FROM MAIN PROGRAM "FORMAL" 
C***  after DECFORM has encountered the input-line 
C***      BLEND (i.e. beginning of a BLEND-Block)
C***   or LINE ... (calculation of one isolated line)
C***   This subroutine is left when encountering the input-line
C***      -BLEND (i.e. end of a BLEND-Block)
C*******************************************************************************
 
C***  POPNUM enters, because MULTIPLET is inserting additional (LTE) levels
C***  The same holds for the the DRTRANSIT option
C***  Hence, all entering parmeters must be indexed with IMOD
      DIMENSION POPNUM(NDDIM, NDIM, MAXMOD)
      DIMENSION ND(MAXMOD)
      REAL, DIMENSION (NDDIM,MAXMOD) :: DENSCON, ENTOT, RNE, T
      DIMENSION DD_VDOP(NDDIM, MAXATOM, MAXMOD) 
        
      INTEGER, DIMENSION(NDIM) :: NOM
      DIMENSION ELEVEL(N)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND), MULTIIND(LASTIND)
      INTEGER, DIMENSION(MAXSUBL) :: NSUBLOW, NSUBNUP
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM)
      CHARACTER KARTE*80, INDSTR*4, KARTE2*80
      CHARACTER*8 LINPRO(MAXLAP), LINPROBL
      LOGICAL BLEND, BROAD, BFINDERR, BAIR, BVAC

C***  4 * PI
      DATA PI4 / 12.56637062 /

      ALN = ALOG (1. - VDOP / CLIGHT)

      BFINDERR = .FALSE.
      NMULTI = 0
      XLAM = -1.
      NBLINE = 0

C***  SAVE NUMBER OF ORIGINAL LEVELS AND LINES:
C***  MULTIPLET HANDLING (IN SUBR. MULTIPL) MAY INCREASE THESE NUMBERS 
      NNEW=N
      INDNEW=LASTIND

C***  levels which aren't the product 
C***   of "splitting" are their own "original levels" 
      DO NLEV=1, N
        IND_ORIGLEV(NLEV) = NLEV
      ENDDO

C***  DEFAULTS:
      BLEND=.FALSE.

      DO NBL = 1, MAXLAP
         LINPRO(NBL) = ''
         AVOIGT(NBL,1,1) = -1.
      ENDDO

C***  Check if the current CARD-line was "BLEND"  -------------------
C***    otherwise, it must have been an isolated LINE card -> re-read 
      IF (KARTE(:5) .EQ. 'BLEND') THEN 
C                         =====
         BLEND = .TRUE.
      ELSE
         BACKSPACE (2)
      ENDIF

C******************************************************************
C***  Begin of loop over all LINE, MULTIPLET or DRTRANSIT blocks
C******************************************************************

C***  Read next line, omit comment lines
    4 READ (2,'(A)',END=100) KARTE
      IF (KARTE(:1) .EQ. '*'  .OR.  KARTE(:1) .EQ. ' '
     >   .OR. KARTE(:1) .EQ. '.' ) GOTO 4
       
C***  Remove leading "+" signs (historical relic)
      IF (KARTE(:5) .EQ. '+LINE') KARTE=KARTE(2:)
      IF (KARTE(:10) .EQ. '+MULTIPLET') KARTE=KARTE(2:)
      IF (KARTE(:10) .EQ. '+DRTRANSIT') KARTE=KARTE(2:)

      IF (KARTE(:6) .EQ. '-BLEND') THEN
         IF (.NOT. BLEND) 
     >       WRITE (0,*) '*** WARNING: "-BLEND" encountered ',
     >                     'albeit no BLEND block was opened before'
         IF (NBLINE .EQ. 0)
     >       WRITE (0,*) '*** WARNING: BLEND-Block contained no lines'
         GOTO 20
       ENDIF

C******************************************************************
      IF (KARTE(:4) .EQ. 'LINE') THEN
C                         ====

         CALL SARGV (KARTE, 2, INDSTR)
         CALL FINDIND (INDBL,INDSTR,LEVEL,N,INDNUP,INDLOW,LASTIND
     >                ,BFINDERR, LOWBL, NUPBL)
         IF (BFINDERR) GOTO 4

C***     default wavelength 
         XLAMBL=1.E8/(ELEVEL(NUPBL)-ELEVEL(LOWBL))

C***     read wavelength if given, covert from AIR to VAC
         CALL SARGREST (KARTE, NPAR, 3, IRESTSTART, IRESTEND)
         IF (NPAR .GT. 0)
     >     CALL READ_LINECARD_PARAMETERS (KARTE(IRESTSTART:), XLAMBL,
     >            BROAD, LINPROBL, AVOIGTBL)

C***     In case that this line is the first one encounetred, 
C***        its wavelength serves as reference value:
         IF (NBLINE .EQ. 0) XLAM = XLAMBL

         CALL INSERT_LINE (LINE, INDBL, NBLINE, INDLAP, XLAMLAP, DELXLAP,
     $                   XLAMBL, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN, 
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C***     for isolated line: finish
         IF (.NOT. BLEND) GOTO 20

C*****************************************************************
      ELSEIF (KARTE(:9) .EQ. 'MULTIPLET') THEN
C                             ---------

         CALL SARGV (KARTE, 2, INDSTR)
         CALL FINDIND (INDBL,INDSTR,LEVEL,N,INDNUP,INDLOW,LASTIND
     >                ,BFINDERR, LOW, NUP)

C***     If LEVELS were not found: Skip the whole multiplet
         IF (BFINDERR) THEN
   44       READ (2,'(A)') KARTE
            IF (KARTE(:10) .EQ. '-MULTIPLET') GOTO 4
            GOTO 44
         ENDIF

C***    Store LINE index for checking if this transition was split
            NMULTI = NMULTI + 1
            MULTIIND(NMULTI) = INDBL

            CALL MULTIPLE (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP,
     >                 DELXLAP, NBLINE,MAXLAP,INDLOW,INDNUP,INDNEW,MAXIND,
     >                 LEVEL,WEIGHT,ELEVEL,NNEW,EINST,NDIM,
     >                 POPNUM, T, ND, ALN, VDOP,
     >                 MAXSUBL, NSUBLOW, NSUBNUP, BROAD, LINPRO, 
     >                 AVOIGT, NMOD, MAXMOD, NDDIM, 
     >                 MAINQN, NCHARG, EION, NOM, IND_ORIGLEV)



      ELSE IF (KARTE(:9) .EQ. 'DRTRANSIT') THEN
C                              ---------

         CALL FINDLDR (LOW, NUP, LEVEL, NOM, NCHARG, N, BFINDERR)

C***     If LEVELS were not found: Skip the whole DRTRANSIT block
         IF (BFINDERR) THEN
  55       READ (2,'(A)') KARTE
            IF (KARTE(:10) .EQ. '-DRTRANSIT') GOTO 4
            GOTO 55
         ENDIF

         CALL DRTRANS (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP, DELXLAP,
     >                 NBLINE, MAXLAP, INDLOW, INDNUP, INDNEW, MAXIND,
     >                 LEVEL, WEIGHT, ELEVEL, NNEW, EINST, NDIM,
     >                 POPNUM, T, ND,
     >                 ALN, VDOP, EION, ENTOT, RNE,
     >                 MAXSUBL, NSUBLOW, BROAD, LINPRO, 
     >                 AVOIGT, DENSCON, NMOD, MAXMOD, NDDIM, 
     >                 MAINQN, NCHARG, NOM, IND_ORIGLEV) 

C***  ERROR branch *******************************
      ELSE
         WRITE (0,*) 'ERROR: INVALID OPTION FOUND'
         WRITE (0,*) KARTE(:IDX(KARTE))
      ENDIF

      GOTO 4
C***  END OF BLEND BLOCK  ****************************************


   20 CONTINUE

      RETURN

  100 WRITE (0,*) '*** UNEXPECTED E-O-F in FORMAL_CARDS ***'
      WRITE (0,*) '-BLEND missing?'
      RETURN

      END
      SUBROUTINE PREPMACROCLUMP (MACROCLUMPLINE, DENSCON, VELO, RADIUS, 
     >            TAUROSS, ND, POROLENGTH)
C******************************************************************** 
C***  Decodes the MACROCLUMP option from input, and creates a vector
C***  over depth index  POROLENGTH(L)  
C******************************************************************** 

      CHARACTER MACROCLUMPLINE*(*), ACTPAR*20, ACTPAR2*20
      DIMENSION DENSCON(ND), VELO(ND), RADIUS(ND), POROLENGTH(ND)
      DIMENSION TAUROSS(ND)

      DATA PI / 3.141592654 /

C***  Defaults
      CLUMP_SEP0 = 1.
      TAU1 = -1.
      TAU2 = -1.
      VELO1 = -1.
      VELO2 = -1.

      CALL SARGC (MACROCLUMPLINE, NPAR)
      DO IPAR=2, NPAR, 2
         CALL SARGV (MACROCLUMPLINE, IPAR, ACTPAR)
         IF (NPAR .LT. IPAR+1) GOTO 90
         CALL SARGV (MACROCLUMPLINE, IPAR+1, ACTPAR2)
         IF (ACTPAR .EQ. 'CLUMP_SEP') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) CLUMP_SEP0
         ELSE IF (ACTPAR .EQ. 'TAU1') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) TAU1
         ELSE IF (ACTPAR .EQ. 'TAU2') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) TAU2
         ELSE IF (ACTPAR .EQ. 'VELO1') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) VELO1
         ELSE IF (ACTPAR .EQ. 'VELO2') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) VELO2
         ELSE
            GOTO 92
         ENDIF      
      ENDDO

      IF ((TAU1.NE.-1. .OR.  TAU2.NE.-1.) .AND. 
     >   (VELO1.NE.-1. .OR. VELO2.NE.-1.)) GOTO 93


      DO L=1, ND
         POROLENGTH(L) = DENSCON(L)**(2./3.) * CLUMP_SEP0 
     >       * (VELO(L) / VELO(1) * RADIUS(L)**2)**(1./3.)
      ENDDO


C***  Optional suppression of macroclumping at large Rosseland depth
C***    or at low velocities. Soft switch between two boundaries  

C***  Soft switch with TAU-Roseland parameters
      IF (TAU1 .NE. -1. .OR.  TAU2 .NE. -1.) THEN
         IF (TAU2 .LT. .0) TAU2 = TAU1
         IF (TAU1 .LT. .0) TAU1 = TAU2

C***   Put TAU1, TAU2 in increasing order
         IF (TAU2 .LT. TAU1) THEN
            TEMP = TAU1
            TAU1 = TAU2
            TAU2 = TEMP
         ENDIF 

         DO L=1, ND
            IF (TAUROSS(L) .LT. TAU1) CYCLE 
            IF (TAUROSS(L) .GT. TAU2) THEN
               POROLENGTH(L) = .0 
            ELSE
               W = (TAU2 - TAUROSS(L)) / (TAU2 - TAU1) 
               Q = 0.5 - 0.5 * COS(PI*W)
               POROLENGTH(L) = POROLENGTH(L) * Q 
            ENDIF
         ENDDO

C***  Soft-switch with VELOCITY parameters
      ELSEIF (VELO1 .NE. -1. .OR.  VELO2 .NE. -1.) THEN
         IF (VELO2 .LT. .0) VELO2 = VELO1
         IF (VELO1 .LT. .0) VELO1 = VELO2

C***   Put VELO1, VELO2 in increasing order
         IF (VELO2 .LT. VELO1) THEN
            TEMP = VELO1
            VELO1 = VELO2
            VELO2 = TEMP
         ENDIF 

         DO L=1, ND
            IF (VELO(L) .GT. VELO2) CYCLE 
            IF (VELO(L) .LT. VELO1) THEN
               POROLENGTH(L) = .0 
            ELSE
               W = (VELO2 - VELO(L)) / (VELO2 - VELO1) 
               Q = 0.5 + 0.5 * COS(PI*W)
               POROLENGTH(L) = POROLENGTH(L) * Q 
            ENDIF
         ENDDO

      ENDIF

C*** Output
      WRITE (*,'(A)') MACROCLUMPLINE(:IDX(MACROCLUMPLINE))
      WRITE (*,'(A)') 'PARAMETERS:'
      WRITE (*,'(A, G10.2)') 'CLUMP_SEP = ', CLUMP_SEP0

      IF (TAU1 .GE. .0) THEN 
         CALL LIPO (V1, TAU1, VELO, TAUROSS, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' TAU1 = ', TAU1, '  ---> VELO1= ', V1
      ENDIF

      IF (TAU2 .GE. .0) THEN
         CALL LIPO (V2, TAU2, VELO, TAUROSS, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' TAU2 = ', TAU2, '  ---> VELO2= ', V2
      ENDIF

      IF (VELO1 .GE. .0) THEN
         CALL LIPO (T1, VELO1, TAUROSS, VELO, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' VELO1= ', VELO1, '  ---> TAU1 = ', T1
      ENDIF

      IF (VELO2 .GE. .0) THEN
         CALL LIPO (T2, VELO2, TAUROSS, VELO, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' VELO2= ', VELO2, '  ---> TAU2 = ', T2
      ENDIF

cc      WRITE (*,*) 
cc      DO L=1, ND
cc         write (*,'(I3, G10.2)') l, porolength(l)
cc      ENDDO

      MACROCLUMPLINE = 'DONE'

      RETURN

C**** ERROR BRANCHES *******************
   90 WRITE (0,*) 
     >    'FATAL ERROR: MACROCLUMP OPTION HAS TOO FEW PARAMETERS'
      GOTO 99

   91 WRITE (0,*) 'FATAL ERROR: CANNOT DECODE ',ACTPAR2(:IDX(ACTPAR2)),
     >       ' AS FLOATING POINT NUMBER'
      GOTO 99

   92 WRITE (0,*) 'FATAL ERROR: UNRECOGNIZED PARAMETER: ',
     >       ACTPAR(:IDX(ACTPAR))
      GOTO 99

   93 WRITE (0,*) 'FATAL ERROR: TAU and VELO range cannot ',
     >      'be used simultaneously!'
      GOTO 99

   99 WRITE (0,*) 'THE ERROR OCCURRED IN THE FOLLOWING LINE:'
      WRITE (0,*) MACROCLUMPLINE
      STOP 'ERROR DETECTED IN SUBR. PREPMACROCLUMP'

       END
      SUBROUTINE PREPRAY (Z, COSPHI, P, ND, NDDIM, NP, JP, XO, LTOT, 
     >                   PWEIGHT, CORE, VDU, R,
     $                   OPAL, ETAL, RRAY, OPARAY, ETARAY,
     $                   ETACRAY, OPALRAY, ETALRAY, NFLDIM, ZRAY,
     $                   XCMF, NDADDIM, NBLINE, MAXLAP, 
     $                   REDIS, ETACK, ETACCK, OPACK, 
     $                   XCMFRED, XCMFBLUE, DXCMF,
     $                   BCORE,BCOREL,DBDR,DBDRL, K,
     >                   POROLENGTH, POROLENGTHRAY,
     >                   RCOROT, VSIDU, 
     >                   DD_VDOPDU_RAY, DD_VDOPDU, NATOM, NBFIRST,
     >                   NBLAST, MAXATOM, 
     >                   BVSINI_AT_RCOROT, NMOD, MAXMOD, 
     >                   ZINTER, XMAX, DELXLAP)
C***********************************************************************
C***  DEFINING WHOLE RAYS INCLUDING THE BACKWARD HALF-SPHERE
C***  All names of the rat vectors end on ...RAY
C***   Note: XCMF should also be called XCMFRAY! 
C***  THE CONTINUUM EMISSIVITY IS CALCULATED TWICE: 
C***  ETARAY:  EMISSIVITY WITH REDISTRIBUTION FROM ETACK
C***  ETACRAY: EMISSIVITY OF THE CONTINUUM FROM ETACCK
C***  If a second model is active (i.e. NMOD=2), the quantities that are
C***  copied to the RAY vectors are taken from IMOD=2 if Z falls in the
C***  interval specified by ZINTER   
C***********************************************************************

      INTEGER, INTENT(IN) :: NATOM      
      REAL, DIMENSION (NDDIM, MAXATOM, MAXMOD) :: DD_VDOPDU
      REAL, DIMENSION (NDADDIM, NATOM) :: DD_VDOPDU_RAY

      REAL, DIMENSION (NDADDIM) :: ZRAY, XCMF, RRAY, 
     >                             OPARAY, ETARAY, ETACRAY
      REAL, DIMENSION (NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL
      REAL, DIMENSION (NDADDIM,MAXLAP) :: OPALRAY, ETALRAY
      DIMENSION DELXLAP(NBLINE)

      DIMENSION VDU(NDDIM,MAXMOD), R(ND)
      DIMENSION Z(ND),P(NP), ZINTER(2)
      REAL, DIMENSION (NDDIM,NFLDIM,MAXMOD) :: ETACK, ETACCK, OPACK
      REAL, DIMENSION(NFLDIM,MAXMOD) :: BCORE, DBDR
      DIMENSION POROLENGTH(ND), POROLENGTHRAY(NDADDIM)
      LOGICAL CORE, REDIS
C***  BCOROT: local boolean, .true. if in co-rotating domain
      LOGICAL :: BCOROT, BVSINI_AT_RCOROT

C***  in case VSINI refers to RSTAR: compensate for the scaling with 1/RCOROT
      IF (BVSINI_AT_RCOROT) THEN
        VROTPROJ1 = - P(JP) * COSPHI * VSIDU   
      ELSE 
        VROTPROJ1 = - P(JP) * COSPHI * VSIDU * RCOROT
      ENDIF

C***  Check if ray intersects core 
      LMAX = MIN0(NP+1-JP,ND)
      CORE=LMAX .EQ. ND    

C***  LTOT = number of grid points along ray including back hemisphere
      IF (CORE) THEN
        LTOT = LMAX
      ELSE
        LTOT = 2*LMAX - 1
      ENDIF

C***  LL runs over both hemispheres, 
C***     while L is only within the range for the front hemisphere
      DO 1 LL = 1, LTOT
        IF (LL .LE. LMAX) THEN
           L = LL
           ZRAY(LL) = Z(L)
        ELSE
           L = 2 * LMAX - LL
           ZRAY(LL) = -Z(L)
        ENDIF

C***    Check if current point lies in second-model domain
        IMOD=1
        IF (NMOD .EQ. 2) THEN
           IF ( (ZRAY(LL)-ZINTER(1))*(ZRAY(LL)-ZINTER(2)) .LT. .0) THEN
              IMOD=2
           ENDIF
        ENDIF

        RRECIP = 1. / R(L)
        VRAD = VDU(L,IMOD)
        VRADPROJ = VRAD*ZRAY(LL)*RRECIP
        BCOROT = (R(L) .LE. RCOROT)
        IF (BCOROT) THEN
          VROTPROJ = VROTPROJ1 / RCOROT
        ELSE
          VROTPROJ = VROTPROJ1 * RCOROT * RRECIP * RRECIP
        ENDIF

        XCMF(LL) = XO - (VRADPROJ + VROTPROJ)

        RRAY(LL) = R(L)
        POROLENGTHRAY(LL) = POROLENGTH(L) 
        IF (XCMF(LL) .GT. XCMFBLUE) THEN
          XC = XCMFBLUE
        ELSE IF (XCMF(L) .LE. XCMFRED) THEN
          XC = XCMFRED
        ELSE
          XC = XCMF(LL)
        ENDIF
        KCMF = INT ( (XCMFBLUE - XC) / DXCMF) + 1
        Q = (XCMFBLUE - XCMF(LL))/DXCMF - (KCMF-1)
        PW= 1.-Q
        ETARAY(LL)  = PW*ETACK (L,KCMF,IMOD) + Q*ETACK (L,KCMF+1,IMOD)
        ETACRAY(LL) = PW*ETACCK(L,KCMF,IMOD) + Q*ETACCK(L,KCMF+1,IMOD)
        OPARAY(LL)  = PW*OPACK (L,KCMF,IMOD) + Q*OPACK (L,KCMF+1,IMOD)

        DO NA=1, NATOM
           DD_VDOPDU_RAY(LL,NA) = DD_VDOPDU(L, NA, IMOD)
        ENDDO
    1 CONTINUE    

C***  FIND THE FIRST AND LAST INVOLVED LINE FOR THAT RAY
C***  due to rotation and/or second-model, XCMF might not be monotonic
      XCMFMIN = MINVAL (XCMF(1:LTOT))
      XCMFMAX = MAXVAL (XCMF(1:LTOT))
      NBFIRST = ISRCHFGT (NBLINE, DELXLAP, 1, XCMFMIN-XMAX)
      NBLAST  = ISRCHFGE (NBLINE, DELXLAP, 1, XCMFMAX+XMAX) - 1

C***  Note: this loop over LL cannot be unified with the above one, 
C***        bacause it needs NBFIRST, NBLAST to be defined before      
      DO LL = 1, LTOT
        IF (LL .LE. LMAX) THEN
           L = LL
        ELSE
           L = 2 * LMAX - LL
        ENDIF
        DO NBL=NBFIRST, NBLAST
            OPALRAY(LL,NBL)=OPAL(L,NBL,IMOD)
            ETALRAY(LL,NBL)=ETAL(L,NBL,IMOD)
        ENDDO
      ENDDO

C***  STORE THE INTERPOLATED BCORE AND DBDR IN LOCAL VARIABLES 
C***  WHEN RAY IS A CORE RAY
C***  Since last depth index was LMAX, IMOD has still the correct value
      IF (CORE) THEN
        BCOREL = PW * BCORE(KCMF,IMOD) + Q * BCORE(KCMF+1,IMOD)
        DBDRL  = PW * DBDR (KCMF,IMOD) + Q * DBDR (KCMF+1,IMOD)
      ENDIF
      
C***  PWEIGHT = WEIGHT p*dp at current P(JP) FOR THE FLUX INTEGRAL
      IF (JP .EQ. 1) THEN
         PWEIGHT = P(2)*P(2)/3.
      ELSE 
         PWEIGHT = (P(JP-1)+P(JP)+P(JP+1))*(P(JP+1)-P(JP-1))/3.
      ENDIF

      RETURN
      END
      SUBROUTINE PRIDWL (VELO,GRADI,VDOP,POPNUM,ENTOT,R,RSTAR,OPAL,
     $ EINST,WEIGHT,FNUEC,XLAM,LOW,NUP,ND,NDIM,DELW,LSDWL,OPA,THOMSON,
     $ ADELW,JOBNUM,MODHEAD)
C**********************************************************************
C***   PRINT OUT OF THE LINE EMISSION CONTRIBUTION PER RADIUS
C**********************************************************************

      PARAMETER ( MUNUM = 64 )

      REAL VELO(ND),GRADI(ND),DELW(ND),POPNUM(ND,NDIM),ENTOT(ND),R(ND),
     $ EINST(NDIM,NDIM),WEIGHT(NDIM),OPAL(ND),OPA(ND),THOMSON(ND)
     $ ,ADELW(ND)
      LOGICAL CORE
      CHARACTER MODHEAD*100

      DATA PI / 3.141592654 /
C***  H IN  ERG/SEC
      DATA H / 6.6261965E-27 /
C***  C IN ANGSTROEM
      DATA C / 2.99792456211E+18 /
C***  CONV=1/2HC   CONVERTS A TO B FOR I IN ERG/CM2SECHZ, L IN ANG
      DATA CONV / 2.5170103E-9 /
 
C***  OPTICAL DEPTH WHERE THE CONTINUUM IS DEFINED
      T1=0.66666666666666
 
C***  RT1 = CONTINUUM RADIUS
      TAUC= 0.
      OPAT=OPA(1)
      L=1
C***  DO WHILE (TAUC.LT.T1 .AND. L.LT.ND)
   10 CONTINUE
         L=L+1
         OPAT1=OPAT
         TAUC1=TAUC
         OPAT=OPA(L)
         TAUC=(OPAT+OPAT1)/2.*(R(L-1)-R(L))+TAUC1
      IF ((TAUC .LT. T1) .AND. (L .LT .ND)) GOTO 10
C***  END DO
      RT1=R(L-1)+(R(L)-R(L-1))/(TAUC-TAUC1)*(T1 -TAUC1)
      ENTOT1=ENTOT(L-1)+(ENTOT(L)-ENTOT(L-1))/(TAUC-TAUC1)*(T1 -TAUC1)
      LT1=L

C***  HEADER OUTPUT
      IF (LSDWL.EQ.1.OR.LSDWL.EQ.4)
     $   PRINT 2, MODHEAD,JOBNUM,XLAM,T1,RT1,ALOG10(ENTOT1)
    2 FORMAT (1X,A,20X,'JOB NR.',I3,//,10X,
     $ 'CONTRIBUTION TO THE EQUIVALENT WIDTH AS FUNCTION OF RADIUS',/
     $ 10X,'                    LINE ',F10.2,'A',
     $ 24X,'STAR RADIUS IS DEFINED AT OPTICAL DEPTH TAU-CONT =',F6.2,/,
     $ 10X,58('-'),14X,'I.E. AT R=',1PE10.3,', LOG(ENTOT)=',0PF6.2,//,
     $ 4X,'NR',11X,'R',7X,'TAU-C',7X,'TAU-S',7X,'TAU-P',6X,'TAU-MU',7X,
     $ 'BETAC',2X,'BETA-BETAC',7X,'DW-EM',7X,'DW(R)',7X,'IW(R)',/)
 
C***  INITIALIZATION
      DO 11 L=1,ND
      ADELW(L)=0.0
   11 DELW(L)=0.0
      RL1=R(1)
      DELW1=0.
      RINTW1=0.
      ADELW1=0.
      ENORM=0.0
      TAUC1=0.
      OPAT1=0.
C***  FXLAM IS THE ASTROPHYSICAL FLUX
      HNU=H/XLAM*C
      FXLAM=PI*FNUEC*C/XLAM/XLAM/HNU
C***  XIC = FLUX/PI
      XIC=FNUEC
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      RL=R(L)
      VELOL=VELO(L)
      GRADIL=GRADI(L)
      SIGMAL=RL/VELOL*GRADIL-1.
      OPALL=OPAL(L)
      DRL=(RL1-RL)
      OPAT=OPA(L)
      TAUC=TAUC1+(OPAT+OPAT1)/2.*DRL
      IF (L .LT. LT1) THEN
         XMUT1=SQRT(1.-RT1*RT1/RL/RL)
      ELSE
         XMUT1=0.
      END IF
      TSOBS=OPALL*RL/VELOL*VDOP
      TSOBP=OPALL/GRADIL*VDOP
      DBS=(1.-EXP(-TSOBS))/TSOBS
      DBP=(1.-EXP(-TSOBP))/TSOBP

C***  INTEGRATION OVER MU  **********
      CORE=.FALSE.
      DXMU1=(1.+XMUT1)/(FLOAT(MUNUM)-1.)
      DXMU2=(1.-XMUT1)/(FLOAT(MUNUM)-1.)
      XMUMAX=SQRT(1.-1./RL/RL)
C***  BETA-BETAC ( INTEGRAL -1 TO  MU(T1): BETA(MU) EXP(-TAUCMU) DMU )
      BETA=0.0
      BETAXI=0.0
      XMU=-1.
      TSOBMU1=TSOBS/(1.+SIGMAL)
      DBMUTC1=(1.-EXP(-TSOBMU1))/TSOBMU1
      ETCMU1=EXP(-TAUC)
      DO 21 K=2,MUNUM
      XMU=XMU+DXMU1
      XMU2=XMU*XMU
      TSOBMU=TSOBS/(1.+XMU2*SIGMAL)
      DBMUTC=(1.-EXP(-TSOBMU))/TSOBMU
      BETA=BETA+(DBMUTC1+DBMUTC)/2.*DXMU1
      ETCMU=EXP(-TAUCMU(XMU,L,R,OPA,ND))
      BETAXI=BETAXI+(DBMUTC1*ETCMU1+DBMUTC*ETCMU)/2.*DXMU1
      ETCMU1=ETCMU
   21 DBMUTC1=DBMUTC
C***  BETAC ( INTEGRAL MU(T1) TO +1: BETA(MU) EXP(-TAUCMU) DMU )
      BETAC=0.0
      XMU=XMUT1
      TSOMUT=TSOBS/(1.+XMUT1*XMUT1*SIGMAL)
      DBMUTC1=(1.-EXP(-TSOMUT))/TSOMUT
      ETCMU1=EXP(-TAUCMU(XMU,L,R,OPA,ND))
      DO 22 K=2,MUNUM
      XMU=XMU+DXMU2
      XMU2=XMU*XMU
      TSOBMU=TSOBS/(1.+XMU2*SIGMAL)
      DBMUTC=(1.-EXP(-TSOBMU))/TSOBMU
      BETAC=BETAC+(DBMUTC1+DBMUTC)/2.*DXMU2
      IF (.NOT. CORE) THEN
         IF (XMU .GT. XMUMAX) THEN
C***        LAST STEP FOR INTEGRATION OF BETAXI
            CORE=.TRUE.
            XMU2=XMUMAX*XMUMAX
            TSOBMU=TSOBS/(1.+XMU2*SIGMAL)
            ETCMU=EXP(-TAUCMU(XMUMAX,L,R,OPA,ND))
            DBMUTCX=(1.-EXP(-TSOBMU))/TSOBMU*ETCMU
            BETAXI=BETAXI+(DBMUTC1*ETCMU1+DBMUTCX)/2.*
     *                    (XMUMAX-(XMU-DXMU2))
            GOTO 22
         ENDIF
         ETCMU=EXP(-TAUCMU(XMU,L,R,OPA,ND))
         BETAXI=BETAXI+(DBMUTC1*ETCMU1+DBMUTC*ETCMU)/2.*DXMU2
         ETCMU1=ETCMU
      ENDIF
   22 DBMUTC1=DBMUTC

C***  CORRECT DEFINITIONS OF BETA, BETAC: CASTOR (1970), MNRAS 149, 111
      BETA1C=BETA/2.
      BETAC=BETAC/2.
      BETA=BETA1C+BETAC

C***  CONTRIBUTION IN ANGSTROEM/RSTAR
      XJBETC=BETAC*XIC
      FACT=RL*RL*ENTOT(L)*EINST(NUP,LOW)/FXLAM*RSTAR*EXP(-TAUC)
C***  EMISSION
      DELEM=POPNUM(L,NUP)*BETA1C*FACT
C***  ABSORBTION
      DELAB=CONV*XLAM*XLAM*XLAM*WEIGHT(NUP)*FACT*
     $ (POPNUM(L,LOW)/WEIGHT(LOW)-POPNUM(L,NUP)/WEIGHT(NUP))*XJBETC
      DELWL=DELEM-DELAB
      DELW(L)=DELWL
      RINTW=RINTW1+(DELW1+DELW(L))/2.*DRL
C***  CALCULATION OF XI-HILLIER (1987, APJ SUPPL. SER. 63, 965):
      ADELWL=POPNUM(L,NUP)*ENTOT(L)*RL*RL*RL*BETAXI
      ADELW(L)=ADELWL
C***  NORMALIZATION CONSTANT FOR XI-HILLIER:
C***  INTEGRAL OF XI-HILLIER OVER D(LOG R) IS PROPORTIONAL TO LINE ENERGY
      ENORM=ENORM+(ADELW1/RL1+ADELW(L)/RL)/2.*DRL
      DELW1=DELWL
      RINTW1=RINTW
      ADELW1=ADELWL
      RL1=RL
      TAUC1=TAUC
      OPAT1=OPAT
C***  OUTPUT:
      IF (LSDWL.EQ.1.OR.LSDWL.EQ.4) PRINT 3,
     $   L,RL,TAUC,TSOBS,TSOBP,TSOMUT,BETAC,BETA1C,DELEM,DELW(L),RINTW
    3 FORMAT (I6,1P10E12.3)
    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
C**** NORMALIZATION OF "ADELW" (TOTAL LINE EMISSION)
      DO 90 L=1,ND
   90 ADELW(L)=ADELW(L)/ENORM
 
      RETURN
      END
      SUBROUTINE PRIOPAL(KARTE,XLAM,ND,OPA,OPAL,ETA,ETAL,R,JOBNUM,LSOPA,
     $          MODHEAD)
C***********************************************************************
C***  PRINTOUT OF THE LINE OPACITIES ETC.
C***********************************************************************

      DIMENSION OPA(ND),OPAL(ND),ETA(ND),ETAL(ND),R(ND)
      CHARACTER MODHEAD*100, KARTE*80
C***  WPI = SQRT(PI)
      DATA WPI /1.772454/
      DATA IHELP / 0 /
 
C***  Print header only once at the beginning
      
      IF (IHELP .EQ. 0) GOTO 1
      IHELP=1
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1H1,1X,  A  , 4X,'AFTER JOB NO.',I3,//,10X,
     $ 'LINE OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,44('-'),//,
     $ ' LINE       DEPTH LINE OPACITY LINE/CONT.   TOTAL OPT.DEPTH  ',
     $ 'R (TAU=1)  LINE EMISS.   LINE SOURCE F.    TOTAL SOURCE F.',/,
     $ '            INDEX  (PER RSTAR)               (LINE CENTER)   ',
     $ '                          TRAD/KELVIN       TRAD/KELVIN   ',/)

    1 PRINT 3
    3 FORMAT (1H )
      TAU=.0
      RTAU1=.0
      DO 6 L=1,ND
      OPALC=OPAL(L)/WPI
      ETALC=ETAL(L)/WPI
      OPATOT=OPA(L)+OPALC
      ETATOT=ETA(L)+ETALC
      IF (L.EQ.ND) GOTO 7
      TAUOLD=TAU
      OPATOTP=OPA(L+1)+OPAL(L+1)/WPI
      TAU=TAU+0.5*(OPATOT+OPATOTP)*(R(L)-R(L+1))
      IF( TAUOLD.GE.1. .OR. TAU.LT.1. )  GOTO 7
      Q=(1.-TAUOLD)/(TAU-TAUOLD)
      RTAU1=(1.-Q)*R(L)+Q*R(L+1)
    7 IF(((L-1)/LSOPA)*LSOPA.NE.(L-1) .AND. L.NE.ND) GOTO 6
      S=ETALC/OPALC
      TRADLIN=TRADFUN (XLAM,S)
      S=ETATOT/OPATOT
      TRADTOT=TRADFUN (XLAM,S)
      IF (L.EQ.ND) GOTO 8
      PRINT 5,      L,OPALC,OPALC/OPA(L),TAUOLD,   ETALC,TRADLIN,TRADTOT
    5 FORMAT ( 10X, I6,1P,3E14.2,    11X,    E14.3,0P,2F16.0)
    6 CONTINUE
    8 PRINT 9,KARTE,L,OPALC,OPALC/OPA(L),TAU,RTAU1,ETALC,TRADLIN,TRADTOT
    9 FORMAT (1X,A9,I6,1P,3E14.2,0P,F11.3,1P,E14.3,0P,2F16.0)
      RETURN
      END
      SUBROUTINE PRI_PAR (TEFF, RSTAR, VFIN, DENSCON, XMDOT)
C     ************************************************************
C     * Short model parameter printout                           *
C     *  This routine is only called from FORMAL                 *
C     ************************************************************

      IMPLICIT NONE

      REAL, PARAMETER :: RSUN = 6.96E10 ! solar radius in cm
      REAL DENSCON, RT, RSTAR, VFIN, TEFF, RSTAR_S, XMDOT

      RSTAR_S = RSTAR / RSUN
      RT = (VFIN/2500. * 1.E-4 / 10.**XMDOT /SQRT(DENSCON))**(2./3.) 
     >      * RSTAR_S
      
      WRITE (*,*)
      WRITE (*,*) '--------------------------'
      WRITE (*,*) 'Model Parameters'
      WRITE (*,*) 
      WRITE (*,'(A, F7.0)') ' T* = ', TEFF
      WRITE (*,'(A, 1PG12.3,0P,A,F7.3,A)') 
     >      ' Rt = ', RT, ' = ', ALOG10(RT), 'dex'
      WRITE (*,'(A, F7.4)') ' R* = ', RSTAR_S
      WRITE (*,'(A, F7.3,A)') ' Md = ', XMDOT, 'dex'
      WRITE (*,'(A, F7.1)') ' V8 = ', VFIN
      WRITE (*,'(A, F7.2)') ' D(L=1)  = ', DENSCON
      WRITE (*,*) '--------------------------'
      WRITE (*,*) 

      RETURN
      END
      SUBROUTINE PRIPRO (XLAM,VDOP,NFOBS,PROFILE,XOBS0,DXOBS,JOBNUM,
     $        VSINI,MODHEAD,DLAM,LSPRO,IFIRST,NPHI,LPSTA,LPEND,DXMAX,
     $        TAUMAX,XMAX,TAUMINBROAD,JFIRST,JLAST,
     >        P, PWEIGHT, ELEVEL,NDIM,INDNUP,INDLOW,
     $        LASTIND,INDLAP,MAXLAP,FREMAX,FREMIN,NBLINE,RMAX,XLAMLAP,
     $        WEIGHT,LEVEL,EINST,VMAX,IVERSION,LINE, REDIS,
     $        BWESEX,BROAD,LINPRO,AVOIGT, BNOCONT, NMOD, 
     >        DENSCON, FILLFAC, 
     >        BFEWARNING, SPZ1, SPZ2, MANIPOP_OPTIONS, ND,
     >        MULTIIND, NMULTI, BAIRWAVELENGTH, MAXSPZ, RCOROT,
     >        BDD_VDOP, DD_VDOP_LINE, BMICROTURB, EL_MIN, 
     >        IVDOPSTATUS, BVSINI_AT_RCOROT)
C***********************************************************************
C***  PRINTOUT OF THE EMERGENT LINE PROFILES
C***********************************************************************

      DIMENSION PROFILE(NFOBS),DLAM(NFOBS),P(JFIRST)
      DIMENSION ELEVEL(NDIM), WEIGHT(NDIM), EINST(NDIM,NDIM)
      DIMENSION INDNUP(LASTIND), INDLOW(LASTIND), MULTIIND(LASTIND)
      DIMENSION INDLAP(NBLINE), XLAMLAP(NBLINE)
      DIMENSION AVOIGT(NBLINE)
      DIMENSION JOBNUM(NMOD)

      LOGICAL SPIKE, REDIS, BROAD, BNOCONT, BFEWARNING, BAIRWAVELENGTH
      LOGICAL BDD_VDOP, BMICROTURB, BVSINI_AT_RCOROT
      CHARACTER LEVEL(NDIM)*10, OVERLAP*30, SUBLINE*1, EL_MIN*2
      CHARACTER*100 MODHEAD(NMOD), ACTPAR
      CHARACTER*(*) DD_VDOP_LINE
      CHARACTER LINPRO(NBLINE)*8, AVOIGTSTR*8
      CHARACTER REDIMES*30, BROAMES*6
      CHARACTER TXHALF*8, TXKM*8
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ), MANIPOP_OPTIONS

      DIMENSION DENSCON(ND,NMOD),FILLFAC(ND,NMOD)

      INTEGER :: KSPIKE
      REAL :: SPIKELAM

C***  CLIGHT = VELOCITY OF LIGHT IN KM/SEC
      DATA CLIGHT /2.99792458E+5/
C***  C4 = 1 / (SIGMA-CLASSIC * 8 PI)     (IN ANGSTROEM**2)
      DATA C4 /1.499E-16/
C***  CFACS: ALL CONSTANT FACTORS USED IN THE CALCULATION OF THE ABSOLUTE
C***  LINE EMISSION (ABLIEM) IN UNITS OF THE SOLAR LUMINOSITY LSUN=3.82E33 5I7/S
C***  CFACS = 1.E8*PI*4.*PI*C/LSUN
      DATA CFACS /3.098E-14/
      DATA PI / 3.14159 /

C***  logarithmic doppler unit (here negative!
      ALN = ALOG (1. - VDOP / CLIGHT)

 
C***  OUTPUT OF HEADER: ONE TIME (VARIABLE "IFIRST" IS SET TO ZERO AFTERWARDS)
      IF (IFIRST .EQ. 1 .OR. LSPRO .GT. 0 .OR. NMOD .GT. 1) THEN
        WRITE (*,*)
        DO IMOD=1, NMOD
          PRINT 4, MODHEAD(IMOD)( :IDX(MODHEAD(IMOD))),JOBNUM(IMOD)
        ENDDO
    4   FORMAT (1X, A, 4X, 'AFTER JOB NO.', I7)
      ENDIF

      XLAMMAX = XLAM * EXP(ALN*FREMIN)
      XLAMMIN = XLAM * EXP(ALN*FREMAX)
      WRITE (*,'(/,A,2F10.2)') 
     >     ' WAVELENGTH RANGE [Angstroem]:', XLAMMIN, XLAMMAX

      IF ( BAIRWAVELENGTH ) THEN
       WRITE (*,'(A,/)') ' WAVELENGTHS REFER TO AIR'
      ELSE
       WRITE (*,'(A,/)') ' WAVELENGTHS REFER TO VACUUM'
      ENDIF


      WRITE (*,'(A,I12)')
     >      ' Total number of lines considered = ', NBLINE
      WRITE (*,'(A,F5.2,A,F5.2)') 
     >      ' Clumping parameters: DENSCON(1,1) = ', DENSCON(1,1), 
     >      '   FILLFAC(1,1) = ', FILLFAC(1,1) 
      IF (MANIPOP_OPTIONS .NE. ' ') WRITE (*,'(2A)') 
     > ' Hot LTE component added: MANIPOP_OPTIONS = ', MANIPOP_OPTIONS


cc      IF (IFIRST .EQ. 1 .OR. LSPRO .GT. 0 .OR.
cc     >  (LPSTA .NE. 1 .OR. JFIRST .NE. 1 .OR. JLAST .EQ. 1) .OR. 
cc     >  NMOD .GT. 1) THEN
        PRINT 1, DXMAX, TAUMAX, XMAX, TAUMINBROAD, NFOBS
    1   FORMAT(' NUMERICAL PARAMETERS: DXMAX=',F4.2,', TAUMAX=',F4.1,
     >        ' , XMAX=',F8.1,', TAUMINBROAD=',F6.3,' , NFOBS=',I6)

        PRINT 111, NPHI, LPSTA, LPEND, JFIRST, JLAST
  111   FORMAT (' AZIMUTH ANGLES NPHI=',I4, ', LPHI FROM',I4,' TO',I4, 
     >  /, ' IMPACT PARAMETER INDEX JP FROM',I4,' TO',I4 )

      IF (IDX(DD_VDOP_LINE) .NE. 0)
     > WRITE(*,'(A, A)') ' VDOP was specified by following line: ', 
     >                  DD_VDOP_LINE

      IF (BDD_VDOP) THEN
        WRITE(*, '(A)') ' VDOP is depth-dependent'
        IF (BMICROTURB) THEN
           WRITE(*, '(A)') ' VDOP depends on element, since VMIC was '
     >                  // 'specified: VDOP^2 = vmic^2 + vtherm^2'
        ENDIF
      ENDIF

      WRITE (*,*) 
     >   'Criterion that defined VDOP for wavelength resolution:'  
      IF (IVDOPSTATUS == 1) THEN
         WRITE (*,*) 'VDOP from MODEL file' 
      ELSEIF (IVDOPSTATUS == 2) THEN
         WRITE (*,*) 'Minimum from all Doppler line widths, element = '
     >               // EL_MIN 
      ELSEIF (IVDOPSTATUS == 3) THEN
         WRITE (*,*) 'specified as VDOP in FORMAL_CARDS'
      ELSEIF (IVDOPSTATUS == 4) THEN
         WRITE (*,*) 'VDOPFE from FEDAT_FORMAL, NO-FECONVOL requested'
      ENDIF

      IF (VSINI .GT. .0 .AND. RCOROT .GT. 1.) THEN
        IF (BVSINI_AT_RCOROT) THEN
           WRITE (*,'(A,G14.4)') 'VSINI refers to RCOROT=', RCOROT 
        ELSE
           WRITE (*,*) 'VSINI refers to RSTAR'
        ENDIF
      ENDIF

      IFIRST=0

      WRITE (*,*)
      DO ISPZ=1, MAXSPZ 
        IF (SPZ1(ISPZ) .EQ. '') EXIT
        WRITE (*,'(1X,A,A10,A,A10,A)') 
     >   'Levels starting with >', SPZ1(ISPZ), 
     >   '<, have been set to Zero, except those starting with >', 
     >   SPZ2(ISPZ), '<'
      ENDDO
      WRITE (*,*)

C***  PRINT HEADER FOR ONE BLEND BLOCK (INDEX OF REFERENCE LINE: >LINE<)
      IF (REDIS) THEN
        WRITE (REDIMES, '(F3.1)') BWESEX
        REDIMES = 'ELECTRON REDISTR. BWESEX=' // REDIMES
      ELSE
        REDIMES = ' '
      ENDIF
      IF (BROAD) THEN
        BROAMES = 'BROAD.'
      ELSE
        BROAMES = ' '
      ENDIF
      PRINT 13, LINE, XLAM, VSINI, RCOROT, VDOP, REDIMES, BROAMES
   13 FORMAT (/, 132('-'),/,'LINE', I5, '   AT', F10.2, ' A',
     $ 3X, 'VSINI/(KM/S)=', F5.0, 
     > 3X, 'RCOROT/RSTAR=', F5.1,
     $ 3X, 'VDOP/(KM/S)=', F5.1, 3X, A, 3X, A,
     $ /, 132('-'), / )

C***  LOOP OVER ALL BLEND COMPONENTS
      DO 14 NBL = NBLINE ,1 ,-1
      IND = INDLAP(NBL)
      LOW = INDLOW(IND)
      NUP = INDNUP(IND)
      XLAMI=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
      F = C4 * XLAMI * XLAMI * 
     >    EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***  MARK SUBLINES BY AN ASTERISK
      IF (IND .LE. LASTIND) THEN
         SUBLINE = ' '
      ELSE
         SUBLINE = '*'
      ENDIF
      IF (LINPRO(NBL) .EQ. 'VOIGT   ') THEN
         WRITE (AVOIGTSTR,'(1PG8.2)') AVOIGT(NBL)
      ELSE
         AVOIGTSTR = ' '
      ENDIF

      XLAMBL = XLAMLAP(NBL)
      IF ( BAIRWAVELENGTH ) THEN
         XLAM2   = XLAMBL * XLAMBL
         XLAMBL = XLAMBL - XLAMBL*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
      ENDIF

         PRINT 6, IND, SUBLINE, XLAMBL, XLAMI,  
     $    LEVEL(LOW), LEVEL(NUP), F, LINPRO(NBL), AVOIGTSTR
    6    FORMAT ('LINE', I5, A1, '  AT', F10.2, ' A',
     $      '  (From dE:', F10.2, ' A)', 3X,
     $      'LEVELS: ', A, ' (LOW) - ', A, ' (UP)', 3X, 'F= ', 
     >      1P,G10.3, 1X,A8,1X,A)
   14 CONTINUE


C***  PRINT WARNINGS ABOUT POSSIBLE FURTHER BLENDS **********************

C***  ESTABLISH BANDWIDTH OF ONE LINE:
      FMIN = XMAX + VMAX / VDOP * SQRT(1.-1./RMAX/RMAX) 
      FMAX = XMAX + VMAX / VDOP

      NADDBL=0
      DO 30 IND=1, LASTIND

C***  CHECK WHETHER THIS LINE IS ALREADY TAKEN INTO ACCOUNT
      DO 25 NBL=1,NBLINE
      IF ((IND .EQ. INDLAP(NBL)) .OR. (IND .EQ. LINE)) GOTO 30
   25 CONTINUE

C***  CHECK WHETHER THIS LINE IS SPLIT INTO A MULTIPLET
      DO 27 IM=1, NMULTI
      IF (IND .EQ. MULTIIND(IM)) GOTO 30
   27 CONTINUE

      LOW=INDLOW(IND)
      NUP=INDNUP(IND)
      XLAMI = 1.E8 / (ELEVEL(NUP)-ELEVEL(LOW))
C***  logarithmic scale
      DELTAX = ALOG (XLAMI / XLAM) / ALN

C***  BLEND, IF DELTAX IN THE INTERVAL [ FREMIN-FMAX) , FREMAX+FMIN ]:
      IF ((DELTAX-FREMIN+FMAX)*(DELTAX-FREMAX-FMIN) .LT. .0) THEN
          F = C4 * XLAMI * XLAMI * 
     *        EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***  ATTENTION: NO WARNING ISSUED IF F-VALUE IS VERY SMALL!
          IF (F .LE. 1.E-8) GOTO 30

          NADDBL=NADDBL+1
          IF (NADDBL .EQ. 1) PRINT 28
   28     FORMAT (/,'WARNING - FURTHER LINE BLENDS EXISTING:')

          IF (XLAMI .LT. XLAMMIN) THEN
             DIFF = XLAMI - XLAMMIN
             WRITE (OVERLAP,'(A,F11.1,A)') 
     >             '(BELOW RANGE BY', DIFF, ' A)'
          ELSE IF (XLAMI .GT. XLAMMAX) THEN 
             DIFF = XLAMI - XLAMMAX
             WRITE (OVERLAP,'(A,F11.1,A)')
     >             '(ABOVE RANGE BY', DIFF, ' A)'
          ELSE
             OVERLAP = '(INSIDE RANGE                                )'
          ENDIF
          PRINT 29, IND, XLAMI, OVERLAP, LEVEL(LOW), LEVEL(NUP), F
   29     FORMAT (2X, 'LINE', I5, ' AT', F10.2,' A',  
     $        4X, A, 
     $        4X, 'LEVELS: ', A, ' (LOW) - ', A, ' (UP)',
     $        4X, 'F= ', 1PG10.3) 
      ENDIF
   30 CONTINUE


C***  PRINT LIST OF TRANSITION WHICH ARE SPLIT IN MULTIPLETS ***************

      IF (NMULTI .GT. 0) THEN
      WRITE (*,'(/,A)') 
     >       'THE FOLLOWING TRANSITIONS WERE SPLIT INTO SUBLINES:'  
      DO IND=1, LASTIND
        DO 36 IM=1, NMULTI
          IF (IND .EQ. MULTIIND(IM)) GOTO 37
   36   CONTINUE
C       Index was not found in MULTIIND:
        CYCLE

   37   CONTINUE
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)
        XLAMI = 1.E8 / (ELEVEL(NUP)-ELEVEL(LOW))

        PRINT 38, IND, XLAMI, LEVEL(LOW), LEVEL(NUP)
   38   FORMAT (2X, 'LINE', I5, ' AT', F10.2,' A',  
     $        4X, 'LEVELS: ', A, ' (LOW) - ', A, ' (UP)')
      ENDDO
      ENDIF


C**************************************
C***  Print Line profile (if requested)
C**************************************

      IF (LSPRO.GT.0) PRINT 8
    8 FORMAT (/,11X,'INDEX    ',  
     $ 'DELTA-NUE           DELTA LAMBDA          LAMBDA        ',
     $ '      RELATIVE         EQUIV.WIDTH',/,
     $          11X,' (K)     ',
     $ '(DOPPLER UNITS)     (ANGSTROEM)         (ANGSTROEM)     ',
     $ '        FLUX           (ANGSTROEM)',/)
 
C***  START VALUES:
      KHALF=0
      KMAX = 1
      PMAX=1.
      PHALF=1.
      EQWI=.0
      ASYM=.0
      ABEQWI=0.
      EMEQWI=0.
      DLAM2=(DLAM(2)-DLAM(1))/2.
      REFDIFF=0.05
      SPIKE=.FALSE.
 
C***  LOOP OVER ALL FREQUENCY POINTS -----------------------------------
      X=XOBS0+DXOBS
      IF (LSPRO.GT.0) PRINT 2, 1,X,DLAM(1),DLAM(1)+XLAM,PROFILE(1),EQWI
    2 FORMAT ( 11X,I4,F13.2,F20.2,F20.2,F20.3,F20.3)
      DO 3 K=2, NFOBS
      X=XOBS0+K*DXOBS
      EQWI=EQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ASYM=ASYM+(PROFILE(K-1)*DLAM(K-1)+PROFILE(K)*DLAM(K))*DLAM2
      IF (PROFILE(K).LT.1.) THEN
         ABEQWI=ABEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ELSE
         EMEQWI=EMEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ENDIF
      IF (PROFILE(K).GT.PMAX) THEN
         KMAX=K
         PMAX=PROFILE(K)
      ENDIF
      IF (LSPRO.GT.0) PRINT 2,K,X,DLAM(K),DLAM(K)+XLAM,PROFILE(K),EQWI
      IF (K .EQ. NFOBS) GOTO 3
C***  SPIKE DETECTOR:
      PROFILK=PROFILE(K)
      DIFF1=PROFILK-PROFILE(K-1)
      DIFF2=PROFILK-PROFILE(K+1)
      IF ((DIFF1 .GT. REFDIFF) .AND. (DIFF2 .GT. REFDIFF)
     $     .AND. (DIFF1*DIFF2 .GT. 0.0)) THEN
             SPIKE =.TRUE.
             KSPIKE=K
             SPIKELAM = DLAM(K)+XLAM
      ENDIF
    3 CONTINUE
C***  ENDLOOP ----------------------------------------------------------

C***  OUTPUT: USED METHOD OF INTEGRATION IN SUBR. ZONEINT
      IF (IVERSION .EQ. 1) PRINT 21
   21 FORMAT (/,31X, '=====  CAUTION: INTEGRATION WAS PERFORMED IN Z ',
     $          ' (IVERSION = 1)  =====')
      IF (IVERSION .EQ. 2) PRINT 22
   22 FORMAT (/,31X, '=====  CAUTION: INTEGRATION VERSION: OLDDTAU ',
     $          ' (IVERSION = 2)  =====')

C***  LOOP OUTPUT:
      IF (SPIKE) PRINT 20, SPIKELAM, KSPIKE
   20 FORMAT (/,24X, '===== SPIKE at ', F20.2,
     >          ' (K=',I3,') IN LINE PROFILE DETECTED',
     >          ' - WARNING FOR QUICK LOOK ONLY =====')

      IF (BNOCONT) PRINT 26, ' '
   26 FORMAT (/,24X,'+++++ WARNING: FORMAL INTEGRAL OF AN OLD MODEL',A1,
     >        'CALCULATED WITHOUT INTERPOLATED XJC +++++')


      IF (LSPRO.LE.0)  PRINT 10, EQWI, EMEQWI, ABEQWI
   10 FORMAT (/,5X, 'EQUIVALENT WIDTHS: TOTAL:',F10.3,' A',
     $        10X, 'EMISSION:', F10.3, ' A',
     $        10X, 'ABSORPTION:', F10.3, ' A', /, 30X, 12('=') )
 
C***  OUTPUT OF INFORMATION ABOUT THE LINE PROFILE
      PHALF=PMAX/2.+0.5
      XPEAK=DLAM(KMAX)+XLAM
      DO 23 K=KMAX+1, NFOBS
        IF (PROFILE(K).GT.PHALF) KHALF=K
   23 CONTINUE
      IF (((KHALF .EQ. 0).OR.(KHALF .EQ. NFOBS)).AND.(.NOT. SPIKE)) THEN
         PRINT 9
    9    FORMAT(/,28X,'===== ATTENTION: ABSORPTION PROFILE! =====',/)
         ENDIF

      IF (KHALF .LT. 1 .OR. KHALF .GE. NFOBS) THEN
         TXHALF = ' UNDEF. '
         TXKM   = ' UNDEF. '
         ELSE
         DENOM = PROFILE(KHALF) - PROFILE(KHALF+1)
         IF (DENOM .NE. 0.) THEN 
            XHALF = DLAM(KHALF) - (DLAM(KHALF)-DLAM(KHALF+1)) 
     >              * (PROFILE(KHALF)-PHALF) / DENOM
            XKM = XHALF / XLAM * CLIGHT
            WRITE (TXHALF, '(F8.2)') XHALF
            WRITE (TXKM  , '(F8.2)') XKM
            ENDIF
         ENDIF

         PRINT 7, PMAX, XPEAK, TXHALF, TXKM
    7    FORMAT(/,5X, 'PEAK INTENSITY:', F7.3, ' AT ', F8.2,' A',
     $          10X, 'HALF WIDTH AT HALF MAXIMUM (RED WING):', A8,
     $          ' A  =', A8, ' KM/SEC', /, 20X, '=======')


C***  CALCULATION OF DIMENSIONLESS MOMENTS (LINE PROFILE)
C***  SOURCE: CASTOR ET AL. 1981, MNRAS 194, 547
      FMOM = CLIGHT / XLAM / VMAX
      W0   = -FMOM * EQWI
      W1   = FMOM * FMOM * ASYM
C>>>>      PRINT 17, W0, W1
   17    FORMAT(/,10X,'DIMENSIONLESS MOMENTS (CASTOR ',
     $          'ET AL. 1981):   W0 = ',F9.5,/,56X,'W1 = ',F9.5)
 
      IF (JFIRST .EQ. JLAST)
     > WRITE (*,'(A,I3,A,F10.4,A,F10.4,A)') 
     $   '          INTENSITY PROFILE AT P(',JFIRST,') = ', P(JFIRST), 
     $   '    FLUX integration weight p*dp=',  PWEIGHT, ' not applied'
 
      IF (BFEWARNING) THEN
      WRITE(*,*)' '
        WRITE(*,901)'WARNING:  ACTIVE IRON-LINES DETECTED BUT NOT 
     >           ACCOUNTED IN SPECTRA'
 901  FORMAT(A)
      ENDIF

      WRITE(*,'(///)')


      RETURN
      END
 
      SUBROUTINE QUADSTARK (GAMMAQUAD, T, XNE, NCHARG, ELEVEL, EION, 
     >                      LINPRO, LEVELLOW, LEVELNUP)
C**********************************************************************
C***  Quadratic Stark effect
C***     following Cowley (1971: Observatory 91, 139)
C***  GAMMAQUAD     - damping parameter - Lorentz profil
C***  T             - temperature
C***  XNE           - electron density
C***  ELEVEL        - Excitation energy of the upper level  
C***  NCHARG        - Ion charge
C***   Note: Z in Cowley's formula is "Charge seen by the active electron"  
C***           i.e. Z = NCHARG + 1
C**********************************************************************

      CHARACTER*10 LEVELLOW, LEVELNUP
      CHARACTER*8  LINPRO

C***  Rydberg Energy (infinite mass) in Kaiser
      DATA RYD / 109 737.3 / 
C** The constant is: 8/sqrt(3) * pi  * hbar^2 * m_e^(-3/2) * k_B^(-1/2)  in cgs (See also Eqs. (11) & (23) in Freudenstein+77,  224)
      DATA COWLEY_FAC / 5.0E-5 /
C***  The effective quantum number is calculated with reference to the
C***  ionization energy. However, EION may be not known if there is no higher 
C***  ion in the data than the considered one. 
      IF (EION .LE. .0) THEN
           LINPRO = 'VOIGT   '
           GAMMAQUAD = .0
           WRITE (0,*) '*** WARNING: QUADSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP, 
     >           '  (ionization energy not known)'
           RETURN
      ENDIF 

C***  "distance from continuum" in kaiser
      DELTAE = EION - ELEVEL


C***    WARNING !!!!!!!!
C***    eigentlich muesste man die Energiedifferenz bei doppelt 
C***    angeregten Zustaenden zu dem entsprechenden angeregten Level des 
C***    oberen Ions bilden, siehe Cowley. Diese Information zu beschaffen
C***    erfordert aber einiges Nachdenken.

C***  If energy is close or above 100 Kayser: skip QUADSTARK
C***  use VOIGT as fallback 
      IF (DELTAE .LT. 100) THEN
           LINPRO = 'VOIGT   '
           GAMMAQUAD = .0
           WRITE (0,*) '*** WARNING: QUADSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP, 
     >           '  (upper level too high)'
           RETURN
      ENDIF 
        
C***  primary effective quantum number n* of the upper level
      EFFQN = (NCHARG+1) * SQRT (RYD / DELTAE)

C**   This Gamma is the FWHM of the Lorentz profile 
C***  (see Freudenstein+77,  224, 2w = GAMMAQUAD = FWHM)
C**   Note: It might be more accurate to use T = const = 10kK, see 
C***  discussion by Cowley + Freudenstein. -- Tomer, 18.11.2014
      GAMMAQUAD = COWLEY_FAC * XNE / SQRT(T) * (EFFQN**2/(NCHARG+2))**2

      RETURN
      END

      SUBROUTINE READ_H_STARKDATA (PATH_LEMKE_DAT,LOW,NUP,LEMKEFOUND,       
     >      STKTA_DWS, NWS, STKTA_TS, NTS,
     >      STKTA_ES,  NES, STKTA_PS, NPS,
     >      NWS_DIM, NTS_DIM, NES_DIM, NPS_DIM)
C*****************************************************************
C     Preparatory subroutine for stark-broadening of hydrogen lines
C     This routines reads for the specified line the broadening table 
C      LEMKE_HI.DAT (default path: /home/corona/wrh/work/wrdata/) 
C      from Simbad, created by Michael Lemke
C      see -->  Lemke M.: 1997, A&A Suppl. 122, 285
C     The table comprises the first4 spectral series
C      (Lyman, Balmer, Bracket, Paschen) fur upper levels till n=22
C
C     Content of this subroutine:
C     The file is opened, the block for the line specified by the 
C      principle quantum numbers (nup,low) is searched, and the 
C      data arrays are read and stored:
C       - One vector each for the three axes 
C          ("scaled wavelength", temperature, el. density)
C       - One matrix (here as vector) with the three indices   
C         for the broadening profile
C     * Note: The profile is NOT NORMALIZED 
C     The subroutine returns LEMKEFOUND=1 when the line was found, 
C         and LEMKEFOUND=0 else. 
C******************************************************************

      CHARACTER*256 PATH_LEMKE_DAT, FILENAME 
      CHARACTER*3 LOCAL_SPECIES

      DIMENSION STKTA_DWS(NWS), STKTA_TS (NTS), STKTA_ES (NES)
      DIMENSION STKTA_PS (NPS)

      LOGICAL STKTA_QHALF      !If TRUE, on half of profile tabulated
      LOGICAL STKTA_WSCA       ! Lambda (in A) o
C***  Note: Both are always true in the used table

      IF (PATH_LEMKE_DAT .EQ. 'default') THEN
         FILENAME = '/home/corona/wrh/work/wrdata/LEMKE_HI.DAT'
      ELSE
         FILENAME = 
     >    PATH_LEMKE_DAT(:IDX(PATH_LEMKE_DAT)) // '/' // 'LEMKE_HI.DAT'
      ENDIF
      KANAL=13


C***  Fortran channel
      LUSTK = 80

      LOCAL_SPECIES=' '
      OPEN(UNIT=LUSTK,FILE=FILENAME, FORM='FORMATTED',
     &               STATUS='OLD',ACTION='READ',ERR=99)
      DO 
           READ(LUSTK,*,END=90) LOCAL_SPECIES(1:3),
     &               LOW_STRK,NUP_STRK,STKTA_QHALF,STKTA_WSCA

            READ(LUSTK,*,END=90) NWS, NTS, NES, NPS
C*** check if arrays are dimension sufficiently large
            IF (NWS .GT. NWS_DIM) GOTO 91
            IF (NTS .GT. NTS_DIM) GOTO 92
            IF (NES .GT. NES_DIM) GOTO 93
            IF (NPS .GT. NPS_DIM) GOTO 94

            READ(LUSTK,*, END=90) (STKTA_DWS(I),I=1,NWS)
            READ(LUSTK,*, END=90) (STKTA_TS(I), I=1,NTS)
            READ(LUSTK,*, END=90) (STKTA_ES(I), I=1,NES)
            READ(LUSTK,*, END=90) (STKTA_PS(I), I=1,NPS)

C***  Check if wanted profile was found, then exit reading
          IF (LOCAL_SPECIES .EQ. 'HI ' .AND.
     >        LOW_STRK .EQ. LOW .AND. NUP_STRK .EQ. NUP) EXIT

      END DO

      CLOSE(UNIT=LUSTK)
      LEMKEFOUND = 1
      RETURN


C***  Set flag if this line was not found in the data table
   90 CONTINUE
      LEMKEFOUND = 0
      RETURN

C*******************************************************
   99 WRITE (0,*) '*** ERROR: cannot open file ', 
     >        FILENAME(:IDX(FILENAME))
      WRITE (0,*) '*** ERROR: Broadening data for H I not available'
      GOTO 100

   91 WRITE (0,*) 'DIMENSION NWS_DIM TOO SMALL:',
     >  ' Needed: ', NWS 
      GOTO 100

   92 WRITE (0,*) 'DIMENSION NTS_DIM TOO SMALL:',
     >  ' Needed: ', NTS
      GOTO 100

   93 WRITE (0,*) 'DIMENSION NES_DIM TOO SMALL:',
     >  ' Needed: ', NES 
      GOTO 100

   94 WRITE (0,*) 'DIMENSION NPS_DIM TOO SMALL:',
     >  ' Needed: ', NPS 
      GOTO 100


  100 STOP '*** FATAL ERROR in READ_H_STARKDATA'

      END
      SUBROUTINE READ_LINECARD_PARAMETERS (KARTE_REST, XLAMBL, 
     >             BROAD, LINPROBL, AVOIGTBL)
C*************************************************************************
C***  This subroutine reads additional (optional) parameters :
C***  - Wavelength (first parameter, no keyword)
C***  - 'AIR' or 'VAC'  (second parameter if wavelenth given; defaults)
C***  - VOIGT (radiation-damping on for this line only)
C***  - VOIGT-parameter a (optional, otherwise calculated from lifetimes)
C***  CALLED FROM FORMAL - PREFORM 
C*************************************************************************

      IMPLICIT NONE

      INTEGER NPAR, I, IDX
      REAL XLAMBL, XLAM2, AVOIGTBL

      CHARACTER*80 KARTE_REST
      CHARACTER*20 ACTPAR
      CHARACTER*8  LINPROBL

      LOGICAL BAIR, BVAC, BROAD

      BAIR = .FALSE.
      BVAC = .FALSE.
      LINPROBL = ''
      AVOIGTBL = -1.
 
C***  Count supplementary parameters
      CALL SARGC (KARTE_REST,NPAR)

C***  Input empty? (Should not happen!)
      IF (NPAR .LE. 0) RETURN

C***  reference wavelength explicitely given as 1st parameter?
      CALL SARGV(KARTE_REST,1,ACTPAR)
      IF (ACTPAR(1:1) .LT. 'A') THEN 
         READ (ACTPAR,'(F10.0)', ERR=100) XLAMBL

C***     Wavelength might be followed by keyword "AIR" or "VAC" 
         IF (NPAR .GE. 2) THEN
            CALL SARGV (KARTE_REST, 2, ACTPAR)
            IF (ACTPAR .EQ. 'AIR') BAIR = .TRUE.
            IF (ACTPAR .EQ. 'VAC') BVAC = .TRUE.
         ENDIF

C***     All calculations will be done for vacuum wavelengths
C***     Therefore, explicitely given wavelenths are converted if appropriate 
         IF (BAIR .OR. (.NOT. BVAC .AND.
     >       XLAMBL .GT. 2000. .AND. XLAMBL .LT. 20000.)) THEN
             XLAM2  = XLAMBL * XLAMBL
             XLAMBL = XLAMBL * (1.0 + 2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
          ENDIF
      ENDIF

C***  Voigt parameter definition by input (in fact outdated!) 
      IF (BROAD) THEN
         DO 12 I=1, NPAR
            CALL SARGV (KARTE_REST, I, ACTPAR)
            IF (ACTPAR .EQ. 'VOIGT') THEN
               IF (I .LT. NPAR) THEN
                  CALL SARGV (KARTE_REST, I+1, ACTPAR)
cc                  IF ((ACTPAR(1:1) .LT. 'A'))
                  READ (ACTPAR,'(F20.0)', ERR=110) AVOIGTBL
                  LINPROBL = 'VOIGT   '
               ELSE
                  WRITE (0,*) '*** WARNING: VOIGT keyword without ' //
     >            'value ignored' 
               ENDIF
            ENDIF
   12    CONTINUE
      ENDIF

      RETURN

C***  ERROR Branches ******************************************

  100 WRITE (0,*) '*** ERROR: Wavelength cannot be decoded as number'
      GOTO 130

  110 WRITE (0,*) '*** ERROR: AVOIGT value cannot be decoded as number'
      GOTO 130

  130 WRITE (0,*) '*** The problematic string is: ' // 
     >      KARTE_REST(:IDX(KARTE_REST))
      STOP '*** FATAL ERROR detected bu subr. READ_LINECARD_PARAMETERS'

      END
      SUBROUTINE READMS(ICHANNEL, X, NDIM, NAME, IERR)
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, X, NDIM, 
     >              'READ', IERR)

      RETURN
      END
      SUBROUTINE REMARK (STRING)

      CHARACTER*(*) STRING

      WRITE (0,'(A)') STRING( :IDX(STRING))

      RETURN
      END
C**********************************************************************
C***  If a second model is employed for a part of the atmosphere: 
C***    the radius points of both models were combined to a new 
C***    RADIUS_MERGED grid.
C***  This subroutine now interpolates all relevent arrays to the  
C***     new RADIUS_MERGED grid 
C***  Note: this routine is also abused after NOWIND for restricting 
C***        the radius grid, irrespective whether NMOD=1 or 2  
C***
C***  Called from: FORMAL
C**********************************************************************

      SUBROUTINE RESCALE_SECMOD (NDDIM, NFLDIM, MAXLAP, MAXMOD, MAXATOM, 
     >   ND, NATOM, NFL, NBLINE, RADIUS, RADIUS_MERGED, ND_MERGED,
     >   VDU, OPAL, ETAL, ETACK, ETACCK, OPACK, OPAFE, ETAFE, 
     >   DD_VDOPDU, DD_VDOP, DD_VMICDU, T, ENTOT, RNE, VEC_SECMOD, 
     >   RSTAR, NMOD)

      DIMENSION ND(MAXMOD), RSTAR(MAXMOD), VEC_SECMOD(NDDIM)  
      REAL, DIMENSION (NDDIM,MAXATOM,MAXMOD) :: DD_VDOPDU, DD_VDOP 
      REAL, DIMENSION (NDDIM, MAXMOD) :: 
     >          RADIUS, VDU, T, ENTOT, RNE, DD_VMICDU      
      REAL, DIMENSION (NDDIM,NFLDIM,MAXMOD) :: ETACK, ETACCK, OPACK,
     >                                         OPAFE, ETAFE
      REAL, DIMENSION (NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL
      REAL, DIMENSION (NDDIM) :: RADIUS_MERGED
      DIMENSION QRSTAR(2)
      REAL, PARAMETER :: RSUN = 6.96E10 ! solar radius in cm

ccc   temporary for test plots:
      DIMENSION XPLOT(100), YPLOT(100)

C***  In each of the models, opacities and emissivities are per  
C***  RSTAR of the respctive model, which might differ.
C***  The flux output (Subr. TRAPLO) multiplies the flux with 
C***  the square of RSTAR(1), i.e. adopting the radius of the first
C***  (main) model. Hence, the opacities and emisivities of the 
C***  second model must be re-scaled by the factor QRSTAR(2)

      QRSTAR(1) = 1.

      IF (NMOD .GT. 1) THEN
       QRSTAR(2) = RSTAR(1) / RSTAR(2)

       IF (ABS(QRSTAR(2)-1.0) .GT. 0.001 ) THEN
        WRITE (0,'(A)')      'WARNING: second model has different RSTAR'  
        WRITE (0,'(A,F8.3)') '         RSTAR(1) =', RSTAR(1)/RSUN   
        WRITE (0,'(A,F8.3)') '         RSTAR(2) =', RSTAR(2)/RSUN   
        WRITE (0,'(A,F8.3)') 'second model scaled by factor ', QRSTAR(2)  
        WRITE (0,'(A)')      'Note that this also affects log g and Mdot'  
       ENDIF
      ENDIF
  
      DO IMOD= 1, NMOD

      VEC_SECMOD = VDU(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (VDU(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

C***  Note: There is only one VMIC option in FORMAL_CARDS which
C***        holds for both models, but since there are parameters 
C***        which refer tp RADIUS or VELO, vmic(r) might differ
C***        between the first and the second model
      VEC_SECMOD = DD_VMICDU(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (DD_VMICDU(1,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = T(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (T(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = ENTOT(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (ENTOT(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = RNE(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (RNE(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      DO NA=1, NATOM

         VEC_SECMOD = DD_VDOPDU(1:ND(IMOD),NA,IMOD) 
         CALL TRANSFORM_RGRID (DD_VDOPDU(1,NA,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = DD_VDOP(1:ND(IMOD),NA,IMOD) 
         CALL TRANSFORM_RGRID (DD_VDOP(1,NA,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))
      ENDDO

      DO NBL=1, NBLINE
     
         VEC_SECMOD = OPAL(1:ND(IMOD),NBL,IMOD) * QRSTAR(IMOD) 
         CALL TRANSFORM_RGRID (OPAL(1,NBL,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

C******* Test plot facility
        IF (nbl .eq. 0) then
C          note: after rescale_secmod, all vectors are over radius of imod=1
c        do l=1, ND(2)
        do l=1, ND(1)
           xplot(L) = alog10(RADIUS(L,1)-.999)
c           yplot(L) = alog10(opal(L,nbl,1))
           yplot(L) = opal(L,nbl,1)
c           yplot(L) = VEC_SECMOD(L)
        enddo

        CALL PLOTANFS (77, '', '', '', '',
     >        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
c     >        XPLOT, YPLOT, ND(2), 'SYMBOL=8 SIZE=-0.05 COLOR=2')
     >        XPLOT, YPLOT, ND(1), 'SYMBOL=8 SIZE=-0.05 COLOR=2')

        do l=1, ND(1)
           xplot(L) = alog10(RADIUS(L,1)-.999)
c           yplot(L) = alog10(opal(L,nbl,2))
           yplot(L) = opal(L,nbl,2)
c           yplot(L) = vdu(L,2)
        enddo

        CALL PLOTCONS (77, XPLOT, YPLOT, 
     >                ND(1), 'SYMBOL=5 COLOR=4')

      endif
***************************************************************


         VEC_SECMOD = ETAL(1:ND(IMOD),NBL,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETAL(1,NBL,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      ENDDO

C***  Loop over frequencies
      DO K=1, NFL

         VEC_SECMOD = ETACK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETACK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = ETACCK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETACCK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = OPACK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (OPACK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = OPAFE(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (OPAFE(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = ETAFE(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETAFE(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      ENDDO
C***  End of loop over frequencies

      ENDDO
C***  End of loop over IMOD (model 1, 2)

C******* Test plot facility
      k = 1
      IF (k .eq. 1) then
        do l=1, ND_MERGED
           xplot(L) = alog10(RADIUS_MERGED(L)-.999)
           yplot(L) = alog10(etacck(L,k,1))
        enddo

        CALL PLOTANFS (77, '', '&2mod1 rescaled  &4mod2 rescaled' , 
     >        'log ETACCK', 'log (radius - 0.999)',
     >        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     >        XPLOT, YPLOT, ND_MERGED, 'SYMBOL=8 SIZE=-0.05 COLOR=2')

        do l=1, ND_MERGED
            xplot(L) = alog10(RADIUS_MERGED(L)-.999)
            if (etacck(L,k,2) .gt. .0) then
               yplot(L) = alog10(etacck(L,k,2))
            else
               yplot(L) = .0
            endif
        enddo

        CALL PLOTCONS (77, XPLOT, YPLOT, 
     >                ND_MERGED, 'SYMBOL=5 COLOR=4')

      endif
***************************************************************


      RETURN

C**********************************************************************
C***  ERROR BRANCHES 
C**********************************************************************


      END
      SUBROUTINE ROTATION_PREP (KARTE, RADIUS, VELO, TAURCONT, 
     >                    ND, RCOROT, NPHI, VSIDU, XMAX,
     >                    ND_MERGED, RADIUS_MERGED, 
     >                    NP_MERGED, PGRID_MERGED, ZGRID_MERGED, 
     >                    NDDIM, NPDIM, NPHIMAX, PHIWEIGHT, PHIARR,
     >                    BVSINI_AT_RCOROT, DX_WINDROT)

C*************************************************************************     
C***  Prepares wind-rotation formalism 
C***  This subroutine can be skipped if wind-rotation is not requested, 
C***    i.e. if no VSINI > .0 has been specified
C***  In the first part, the input line is analyzed on basis of model1
C*************************************************************************   

      IMPLICIT NONE     
      CHARACTER KARTE*(*), ACTPAR*20, LINE*200, ACTPAR2*40
      REAL, DIMENSION(NDDIM) ::  RADIUS, VELO, TAURCONT
      REAL RADIUS_MERGED(NDDIM)
      INTEGER, DIMENSION(NPDIM) :: NPHI, LPHIEND
      REAL, DIMENSION(NPDIM) :: PGRID_MERGED
      REAL, DIMENSION(NDDIM*NPDIM) :: ZGRID_MERGED
      REAL, DIMENSION(NPHIMAX,NPDIM) :: PHIWEIGHT, PHIARR
      LOGICAL BVSINI_AT_RCOROT      

C***  because of IMPLICIT NONE:
      REAL PI, RCOROT, RCOROT_V, RCOROT_TAU, XMAX, VSIDU, DPMAX
      REAL DX_WINDROT, RR, COSPHI, PJ, PJPJ, VSIDU_AT_P, DELTACOSPHI
      INTEGER IDX, NCORE, NCORE_ADD, NCORE_NEW, NP_MERGED, NP_NEW
      INTEGER NDDIM, NPDIM, NPHIMAX, ND, ND_MERGED, JP, L, I, J, JMAX 
      INTEGER LPHI, MAXNPHI, NPHISUM, MAXJP, NPHI_MEAN

      DATA PI / 3.14159265358979 /
      
C***********************************************************************
C***  Part 1:interpret the RCOROT line; 
C***    this is done on basis of MODEL1 
C***  alternative parameters are:
C***    RSTAR = x.x
C***    TAUROSS = x.x
C***    VELOKMS = x.x 
C***********************************************************************
      CALL SARGV(KARTE, 1, ACTPAR)
C***  In every case, RCOROT in all units (velocity, optical depth, 
C***    star radius) is interpolated from the given one.      
      IF (ACTPAR .NE. 'RCOROT') THEN
        WRITE (0,'(2A)') 'RCOROT not specified in FORMAL_CARDS, ', 
     >            'adopting the default: RCOROT RSTAR=1.0' 
        RCOROT = 1.        
        CALL SPLINPO (RCOROT_V, RCOROT, VELO, RADIUS, ND)
        CALL SPLINPO (RCOROT_TAU, RCOROT, TAURCONT, RADIUS, ND)   
      ELSE
        CALL SARGV(KARTE, 2, ACTPAR)
        IF (ACTPAR .EQ. 'VELOKMS') THEN            
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT_V
                CALL SPLINPO (RCOROT, RCOROT_V, RADIUS, VELO, ND)
                CALL SPLINPO (RCOROT_TAU, RCOROT_V, TAURCONT, VELO, ND)
        ELSE IF (ACTPAR .EQ. 'TAUROSS') THEN
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT_TAU
                CALL SPLINPO (RCOROT, RCOROT_TAU, RADIUS, TAURCONT, ND)
                CALL SPLINPO (RCOROT_V, RCOROT_TAU, VELO, TAURCONT, ND)
        ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT
                CALL SPLINPO (RCOROT_V, RCOROT, VELO, RADIUS, ND)
                CALL SPLINPO (RCOROT_TAU, RCOROT, TAURCONT, RADIUS, ND)
        ELSE
                WRITE(0,*) '*** ERROR: RCOROT specification invalid'
                WRITE(0,*) '*** Allowed keywords:  VELOKMS | RSTAR |', 
     >                     'TAUROSS'
                GOTO 99
        ENDIF
  
       WRITE (0, '(A)') 'Corotation radius specified by input line:'
       WRITE (0, '(A)') KARTE(:IDX(KARTE))
       WRITE (0, '(A,G15.5)') 
     >      'Corresponding radius [in Rstar]: ', RCOROT
       WRITE (0, '(A,G15.5)') 
     >      'Corresponding velocity: [in km/s]: ', RCOROT_V
       WRITE (0, '(A,G15.5,/)') 
     >      'Corresponding optical depth (continuum Rosseland mean):  ',  
     >      RCOROT_TAU  

      ENDIF



C***********************************************************************
C***  Part 2:
C***  Add more impact parameters if necessary to resolve rotation 
C***    profile with DX_WINDROT Doppler units
C***  This id done on basis of the _MERGED arrays, which may differ
C***    from model 1 in case of SECOND_MODEL
C***********************************************************************

C***  XMAX should be big enough to contain potential rotational broadening 
      XMAX = AMAX1(XMAX, VSIDU)

C***  First, express Delta-P in Doppler units
      DPMAX = VSIDU / DX_WINDROT 
      IF (BVSINI_AT_RCOROT) DPMAX = DPMAX / RCOROT 

C***  Requested number of core rays
      NCORE_NEW = INT(DPMAX) + 1  

C***  Was the original number of core-rays sufficient?
      NCORE = NP_MERGED - ND_MERGED
      IF (NCORE_NEW .LE. NCORE) THEN
        WRITE (0,'(A,I4)') 
     >   'Number of core impact parameters is sufficient: NCORE=', NCORE
        GOTO 1
      ENDIF

C***  Not sufficient: determine the necessary number NCORE_NEW
      NCORE_ADD =  NCORE_NEW - NCORE
      NP_NEW = NP_MERGED + NCORE_ADD
      IF (NP_NEW .GT. NPDIM) THEN
          WRITE(0,'(A,I4)') '*** ERROR: insufficient NPDIM=', NPDIM
          WRITE(0,'(A,I4)') 'Needed for wind rotation: NPDIM=', NP_NEW
          STOP '*** FATAL ERROR in subroutine ROTATION_PREP'
      ENDIF

      WRITE (0,'(A,I4,A,I4)') 
     >       'Number of core impact parameters enhanced from ', NCORE,
     >       ' to', NCORE_NEW

C***  Shift the non-core impact parameters 
      DO JP=NP_MERGED, NCORE, -1
         PGRID_MERGED(JP+NCORE_ADD) = PGRID_MERGED(JP)
      ENDDO

C***  Calculate equally spaced core impact parameters 
      DO JP=2, NCORE_NEW
         PGRID_MERGED(JP) = (JP-1)/FLOAT(NCORE_NEW)
      ENDDO
         
      NP_MERGED = NP_NEW

C***  Z-array must be updated with new P-vector
      DO L=1,ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED +1 -L
        DO J=1,JMAX
          PJ=PGRID_MERGED(J)
          PJPJ=PJ*PJ
          I=(J-1)*ND_MERGED+L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I)=SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I)=.0
          ENDIF
        ENDDO
      ENDDO
C***  End of branch for insertion of more core impact parameters 

    1 CONTINUE

C***********************************************************************
C***  Determine the number of phi angles NPHI
C***    needed to resolve line profiles with DX_WINDROT Doppler units
C***  Note: NPHI is individual for each impact parameter, depending on
C***    the maximum rotation of the shell with radius P(JP)  
C***********************************************************************
      MAXNPHI = 0
      NPHISUM = 0
      DO JP=1, NP_MERGED 
        IF (PGRID_MERGED(JP) .LE. RCOROT) THEN
          VSIDU_AT_P = VSIDU * PGRID_MERGED(JP)
        ELSE 
          VSIDU_AT_P = VSIDU * RCOROT * RCOROT / PGRID_MERGED(JP)
        ENDIF
        IF (BVSINI_AT_RCOROT) VSIDU_AT_P = VSIDU_AT_P / RCOROT 

        NPHI(JP) = 2*INT(VSIDU_AT_P/DX_WINDROT) + 1
        NPHI(JP) = MAX(2, NPHI(JP))
        IF (NPHI(JP) .GT. MAXNPHI) THEN
           MAXNPHI = NPHI(JP)
           MAXJP   = JP
        ENDIF
        NPHISUM = NPHISUM + NPHI(JP)
ccc        WRITE(0,*) JP, P(JP), NPHI(JP)
      ENDDO  
      
      WRITE(0,'(A,I4,A,I4,A,F10.3)') 
     > 'Maximum number of phi angles:', MAXNPHI, 
     >     ' at P(', MAXJP, ') =', PGRID_MERGED(MAXJP)

      NPHI_MEAN = NINT(FLOAT(NPHISUM)/FLOAT(NP_MERGED))
      WRITE(0,'(A,I4,/)') 
     > 'Average number of phi angles per impact parameter:', NPHI_MEAN 

      DO JP=1, NP_MERGED  

C***     Define angles: equally spaced in cosine
         DELTACOSPHI = 2. / (NPHI(JP)-1)
         DO LPHI = 1, NPHI(JP)
           COSPHI = 1. - DELTACOSPHI*(LPHI-1) 
           COSPHI = MAX(COSPHI, -1.)
           PHIARR(LPHI,JP) = ACOS(COSPHI)
         ENDDO

C***    Calculate the angle integration weights by trapezoidal rule  
C***      but omitting the factors 0.5 at all weights
         PHIWEIGHT(1,      JP) = PHIARR(2,JP) - PHIARR(1,JP)
         DO LPHI=2, NPHI(JP)-1
           PHIWEIGHT(LPHI,JP) = PHIARR(LPHI+1,JP) - PHIARR(LPHI-1,JP)
         ENDDO
         PHIWEIGHT(NPHI(JP),JP) = PHIARR(NPHI(JP),JP)
     >                           - PHIARR(NPHI(JP)-1,JP)

C***     Normalization: since phi = 0 ... pi, the sum of 
C***       PHIWEIGHT is normalized to 2*pi
         DO LPHI=1, NPHI(JP)
            PHIWEIGHT(LPHI,JP) =  PHIWEIGHT(LPHI,JP) / (2.*PI) 
         ENDDO

      ENDDO

      RETURN 
      
   99 WRITE (0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE (0,'(A)') '*** The error occured in the following line:'
      WRITE (0,'(A)') KARTE(:IDX(KARTE))
      STOP '*** Fatal Error in rotation_prep'
   
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
      SUBROUTINE SARGREST (TEXT, N, I, IFIRST, LAST)
C**********************************************************************
C***  Ermittelt den Beginn des i-ten Arguments und das nichtleere Ende 
C***  gesamten restlichen Strings. 
C***  Balancierte "..." am Anfang und Ende werden entfernt, sonst 
C***  bleiben " erhalten. Faengt des Argument mit einem Trennzeichen
C***  an (=,:/), so muss (!!!, sonst Fehlerabbruch!) der Reststring in 
C***  "..." eingeschlossen werden, ebenso natuerlich wenn der 
C***  Reststring mit einem Blank beginnen soll. An spaeterer Position sind
C***  Trennzeichen ohne Wirkung. 
C**********************************************************************

      CHARACTER TEXT*(*)

      CALL SARGP (TEXT, N, I, IFIRST, ILAST)
      LAST = IDX(TEXT)

C***  Leerer Reststring
      IF (IFIRST .EQ. -1) THEN
         IFIRST = 1
         LAST   = 1

C***  Der Reststring  beginnt mit einem Trennzeichen ohne " davor
      ELSE IF (IFIRST .EQ. 0) THEN
         WRITE (6, '(A)') 'Error when parsing the following string:'
         WRITE (6, '(A)') TEXT
         WRITE (6, '(A, I3, A)') 'Argument ', I, 
     >                        ' begins with delimiter -> use "..."'
         STOP '>>> ERROR IN SUBROUTINE SARGREST <<<'
      ENDIF

      IF (IFIRST .GT. 1) THEN
C***  Reststring begins with "
         IF (TEXT(IFIRST-1:IFIRST-1) .EQ. '"') THEN 
C***     Remove balanced closing quote if present, else restore leading quote 
            IF (TEXT(LAST:LAST) .EQ. '"' .AND. LAST .GT. IFIRST) THEN
               LAST = LAST - 1
               ELSE
               IFIRST = IFIRST - 1
            ENDIF
         ENDIF
      ENDIF

c      print *, 'Test: ', text
c      print *, 'Test: ', ifirst 
c      print *, 'Test: ', text(ifirst:ifirst)
c      print *, 'Test: ', last
c      print *, 'Test: ', text(last:last)

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
      FUNCTION SDOT (N, X, IX, Y, IY)

      DIMENSION X(N), Y(N)

      SUM = 0.

      DO I=0, N-1
        SUM = SUM + (X(IX+I) * Y(IY+I))
      ENDDO

      SDOT = SUM

      RETURN
      END
      SUBROUTINE SECOND

      STOP 'SECOND NOT IMPLEMENTED AT DEC/OSF'

      RETURN
      END
      SUBROUTINE SECONDMODEL_DEFINE (SECONDMODEL_LINE, NMOD, 
     >             SECONDMODEL_PATH, SECONDMODEL_CHANGED,
     >             RMAX, SECMOD_RRANGE, IGNORE_DIFF)
C******************************************************************
C***  Reads the pathname from the SECONDMODEL_LINE, and checks whether
C***  this name has changed or the option is swiched OFF
C***  Moreover, the radius-range in which the radius-grid needs to be 
C***  refined is estimated here:
C***  For CONE: 
C***      SECMOD_RRANGE = [RMAX, 1]
C***  For SPHERE:
C***      SECMOD_RRANGE = [DSPHERE+RSPHERE, DSPHERE-RSPHERE]
C******************************************************************

      CHARACTER SECONDMODEL_LINE*(*), ACTPAR*100, ACTPAR2*100
      CHARACTER SECONDMODEL_PATH*(*), NEW_SECONDMODEL_PATH*400
      CHARACTER SHAPE*10
      LOGICAL SECONDMODEL_CHANGED, IGNORE_DIFF
      DIMENSION SECMOD_RRANGE(2)

      DATA NEW_SECONDMODEL_PATH / 'UNDEFINED' /

C***  Secondmodel active
      NMOD = 2

C***  Decode parameters from input line
      CALL SARGC (SECONDMODEL_LINE, NPAR)
      DO IPAR=1, NPAR
         CALL SARGV (SECONDMODEL_LINE, IPAR, ACTPAR)

         IF (ACTPAR .EQ. 'SHAPE') THEN
C*                        =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,SHAPE)
            ENDIF

         ELSEIF (ACTPAR .EQ. 'RSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) RSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'DSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'PATH') THEN
C*                            ====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1, 
     >                     NEW_SECONDMODEL_PATH)
            ENDIF

         ELSEIF (ACTPAR(:11) .EQ. 'IGNORE_DIFF') THEN
C*                                 ===========
            IGNORE_DIFF = .TRUE.

         ELSEIF (ACTPAR .EQ. 'OFF') THEN
C*                            ===
            NMOD = 1
         ENDIF
      ENDDO         

C***  Check for mandatory parameters
      IF (NEW_SECONDMODEL_PATH .EQ. 'UNDEFINED') THEN
         GOTO 901
      ELSEIF (SHAPE .EQ. 'UNDEFINED') THEN
         GOTO 92
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         IF (RSPHERE .EQ. -999.) GOTO 97
         IF (DSPHERE .EQ. -999.) GOTO 971
         IF (RSPHERE .GE. DSPHERE+RMAX) GOTO 972
      ENDIF

      IF (SHAPE .EQ. 'CONE') THEN
         SECMOD_RRANGE (1) = RMAX 
         SECMOD_RRANGE (2) = 1. 
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         SECMOD_RRANGE (1) = MIN(RMAX, DSPHERE+RSPHERE) 
         SECMOD_RRANGE (2) = MAX(1.  , DSPHERE-RSPHERE) 
      ENDIF

C***  Check if PATH has changed
      IF (NEW_SECONDMODEL_PATH .NE. SECONDMODEL_PATH) THEN
         SECONDMODEL_CHANGED = .TRUE.
         SECONDMODEL_PATH = NEW_SECONDMODEL_PATH
      ELSE
         SECONDMODEL_CHANGED = .FALSE.
      ENDIF

      RETURN

C*********************************************************************
C***  ERROR branches  ************************************************
C*********************************************************************
   90 WRITE (0,*) '*** ERROR: Option ', ACTPAR(:IDX(ACTPAR)), 
     >            ' needs a value (keyword)'
      GOTO 99

  901 WRITE (0,*) '*** ERROR: PATH to the SECONDMODEL is not defined'
      GOTO 99

   92 WRITE (0,*) '*** ERROR: mandatory parameter SHAPE is missing'
      GOTO 99

   97 WRITE (0,*) '*** ERROR: mandatory parameter RSPHERE is missing'
      GOTO 99

  971 WRITE (0,*) '*** ERROR: mandatory parameter DSPHERE is missing'
      GOTO 99

  972 WRITE (0,*) '*** ERROR: SPHERE covers the whole atmosphere'
      GOTO 99

   98 WRITE (0,*) '*** ERROR: Parameter cannot be decoded as number:'
      WRITE (0,*) '*** ERROR: ', ACTPAR(:IDX(ACTPAR)), '=', 
     >                           ACTPAR2(:IDX(ACTPAR2))
      GOTO 99


   99 STOP '*** FATAL ERROR in subroutine SECONDMODEL_DEFINE'

      END
      SUBROUTINE SECONDMODEL_PREP (ZINTER, NPHI, P, NP, NPDIM, NPHIMAX, 
     >     PHIARR, PHIWEIGHT, PHI_VEC, SECONDMODEL_LINE, 
     >     JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG)
C******************************************************************
C***  CALCULATES THE POINTS OF INTERSECTION FOR EITHER:
C***  - THE CONE geometry:
C***     THETA IS THE OPENING ANGLE OF THE CONE
C***     CONEI is its INCLINATION angle
C***  - THE SPHERE geometry:
C***    RSPHERE= radius of the sphere
C***    DSPHERE= radial distance of the sphere's center from the origen
C***    DELTASPHERE= the sphere-center's elevation angle 
C***    ALPHASPHERE= meridian angle to the sphere center (see Manual)
C***    
C***  The pre-existing set of angle points is enlarged to cover the
C***  cone, but not wasting points outside
C******************************************************************

      DIMENSION ZINTER(2, NPDIM, NPHIMAX), P(NP)
      DIMENSION PHIARR(NPHIMAX,NPDIM), PHIWEIGHT(NPHIMAX,NPDIM)
      DIMENSION NPHI(NPDIM)
      DIMENSION PHI_VEC(NPHIMAX)
      CHARACTER SECONDMODEL_LINE*(*), ACTPAR*100, ACTPAR2*100
      CHARACTER SHAPE*10

      PARAMETER (NPHIMAX_TEMP = 5000)
      DIMENSION PHI_VEC_TEMP(NPHIMAX_TEMP)
      LOGICAL BPHI_ORIG(NPHIMAX_TEMP), BPHI_ORIG_LAST

      DATA PI / 3.14159265358979 /

C***  Defaults
      DATA PHI_REFINE / 1. /
      DATA SHAPE / 'UNDEFINED' /
      DATA THETA_DEG / -999. /
      DATA CONEI_DEG / -999. /
      DATA DELTA_DEG / -999. /
      DATA ALPHA_DEG / -999. /
      DATA RSPHERE   / -999. /
      DATA DSPHERE   / -999. /

      RMAX = P(NP)

C***  Decode parameters from input line
      CALL SARGC (SECONDMODEL_LINE, NPAR)
      DO IPAR=1, NPAR
         CALL SARGV (SECONDMODEL_LINE, IPAR, ACTPAR)

         IF (ACTPAR .EQ. 'SHAPE') THEN
C*                        =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,SHAPE)
            ENDIF

         ELSEIF (ACTPAR .EQ. 'THETA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) THETA_DEG
            ENDIF
         
         ELSEIF (ACTPAR .EQ. 'CONEI') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) CONEI_DEG
            ENDIF
         
         ELSEIF (ACTPAR .EQ. 'DELTA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DELTA_DEG
            ENDIF

         ELSEIF (ACTPAR .EQ. 'ALPHA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) ALPHA_DEG
            ENDIF

         ELSEIF (ACTPAR .EQ. 'RSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) RSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'DSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'PHI_REFINE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) PHI_REFINE
            ENDIF
         ENDIF
      ENDDO         

C***  Check for mandatory parameters
      IF (SHAPE .EQ. 'UNDEFINED') THEN
         GOTO 92
      ELSEIF  (SHAPE .EQ. 'CONE') THEN
         IF (THETA_DEG .EQ. -999.) GOTO 93
         THETA = THETA_DEG * PI / 180.
         IF (CONEI_DEG .EQ. -999.) GOTO 94
         CONEI = CONEI_DEG * PI / 180.

         IF (THETA_DEG .GE. 89.) THEN
            WRITE (0,'(A,F6.1,A)') '***** ERROR: ' //
     >        'Cone opening angle given: THETA=', THETA_DEG, ' degrees'
            WRITE (0,'(A)') '***** ERROR: ' //
     >        'Cone must have opening angle THETA .le. 89 degrees'
            GOTO 99
         ENDIF

         IF (CONEI_DEG - THETA_DEG .LE. .0) THEN
            WRITE (0,'(A, /, A, F6.1, A, F6.1, A)') '***** ERROR: ' //
     >       ' Invalid choice of geometry:', 'THETA=', THETA_DEG,
     >       'deg,  CONE-Inclination=', CONEI_DEG, 'deg'
            WRITE (0,'(A)') '***** ERROR: ' //
     >       'Cone MUST be seen from the side'
            GOTO 99
         ENDIF

      ELSEIF  (SHAPE .EQ. 'SPHERE') THEN
         IF (DELTA_DEG .EQ. -999.) GOTO 95
         DELTA = DELTA_DEG * PI / 180.
         IF (ALPHA_DEG .EQ. -999.) GOTO 96
         ALPHA = ALPHA_DEG * PI / 180.
         IF (RSPHERE .EQ. -999.) GOTO 97
         IF (DSPHERE .EQ. -999.) GOTO 971
         IF (RSPHERE .GE. DSPHERE+RMAX) GOTO 972

      ELSE
         GOTO 91
      ENDIF

C***  Preparation for the CONE case
      IF (SHAPE .EQ. 'CONE') THEN
         SINI = SIN(CONEI)
         COSI = COS(CONEI)
         COTMINUS = 1. / TAN(CONEI-THETA) 
         COTPLUS  = 1. / TAN(CONEI+THETA) 

C***  Preparation for the SPHERE case
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         SIND = SIN (DELTA)
         COSD = COS (DELTA)
         SINA = SIN (ALPHA)
         COSA = COS (ALPHA)
         RSPHERE2 = RSPHERE**2
      ENDIF

      RMAX2 = RMAX * RMAX
      MAXNPHI = 0
      NPHISUM = 0

C***  Loop over impact parameters
      DO JP=1, NP-1

C***  The angle-points will be established  

C***  Only 1 angle point for JP=1 (center)
         IF (JP .EQ. 1) THEN
            PHI_VEC_TEMP(1) = .0
            NPHI_TEMP = 1
            GOTO 10
         ENDIF         

C***  Preparation of angle points
         IF (NPHI(JP) .EQ. 1) THEN
            PHI_VEC(1) = .0
            PHI_VEC(2) = 2.*PI
            NPHI_JP = 2
         ELSE

C**     Copy pre-existing phi points (from wind rotation)
            NPHI_JP = NPHI(JP)
            DO LPHI=1, NPHI_JP
               PHI_VEC(LPHI) = PHIARR (LPHI,JP)
            ENDDO

C**         mirror the list to the southern hemisphere 
            IF ((2 * NPHI_JP - 1) .GT. NPHIMAX) THEN
               WRITE (0,200) NPHIMAX
               WRITE (*,200) NPHIMAX
  200    FORMAT ('*** ERROR: More PHI points needed than dimensioned',
     >          /, 'NPHIMAX= ', I6)
               STOP '*** FATAL ERROR in subr. SECONDMODEL_PREP'
            ENDIF
            DO LPHI=1, NPHI_JP-1         
               PHI_VEC(2*NPHI_JP-LPHI) = 2.*PI -  PHI_VEC(LPHI)
            ENDDO
            NPHI_JP = 2 * NPHI_JP - 1
         ENDIF

C***     Now we add a fine grid of angle points, spaced by the requested
C***     angular resolution DELTA_PHI

C***     Resolution in cone
         IF (SHAPE .EQ. 'CONE') THEN
            NPHI_PER_CONE = NINT(10 * PHI_REFINE)
            DELTA_PHI = 2. * THETA / NPHI_PER_CONE

C***     Resolution in sphere
         ELSEIF (SHAPE .EQ. 'SPHERE') THEN
            NPHI_ACROSS_SPHERE = NINT(10 * PHI_REFINE)
            DELTA_PHI = MIN(2.*PI, RSPHERE / P(JP)) / NPHI_ACROSS_SPHERE

         ELSE
            STOP '*** ERROR: Invalid SHAPE in SECONDMODEL_PREP'
         ENDIF

         PHI_VEC_TEMP(1) = PHI_VEC(1)
         LPHI = 1
         LPHI_TEMP = 1
         DO  
            PHI_NEXT = PHI_VEC_TEMP(LPHI_TEMP) + DELTA_PHI
            LPHI_TEMP = LPHI_TEMP + 1

C***        Error stop if vector PHI_VEC_TEMP too short
            IF (LPHI_TEMP .GT. NPHIMAX_TEMP) THEN
              WRITE (0,201) NPHIMAX_TEMP
              WRITE (*,201) NPHIMAX_TEMP
  201         FORMAT ('*** ERROR: More PHI points needed than' // 
     >                ' dimensioned', /, 'NPHIMAX_TEMP= ', I6)
              STOP '*** FATAL INTERNAL ERROR in subr. SECONDMODEL_PREP'
            ENDIF

            IF (PHI_VEC(LPHI+1) .LT. PHI_NEXT .OR.
     >          PHI_NEXT .GE. 2*PI) THEN
               LPHI = LPHI + 1
               PHI_VEC_TEMP(LPHI_TEMP) = PHI_VEC(LPHI) 
               BPHI_ORIG(LPHI_TEMP) = .TRUE.
            ELSE
               PHI_VEC_TEMP(LPHI_TEMP) = PHI_NEXT
               BPHI_ORIG(LPHI_TEMP) = .FALSE.
            ENDIF
            IF (LPHI .GE. NPHI_JP) EXIT
         ENDDO

         NPHI_TEMP = LPHI_TEMP

   10    CONTINUE
C***     Angle points are now in PHI_VEC_TEMP

C********Loop over all angle points at current JP to find intersections **
         DO LPHI_TEMP=1, NPHI_TEMP
            PHI = PHI_VEC_TEMP (LPHI_TEMP)
            X0  = P(JP) * COS(PHI)
            X02 = X0 * X0
            Y0  = P(JP) * SIN(PHI)
            Y02 = Y0 * Y0 

C***        Intersection points with cone
            IF (SHAPE .EQ. 'CONE') THEN
              ZM = 0.5 * Y0 * (COTPLUS + COTMINUS)
              AAXIS = 0.5 * Y0 * (COTPLUS - COTMINUS)
              AAXIS2 = AAXIS * AAXIS
              BAXIS2 = (Y02 + ZM**2) * (TAN(THETA))**2
              IF (X02 .GE. BAXIS2) THEN
                 Z1 = .0
                 Z2 = .0
              ELSE
                 TERM2 = AAXIS2 * (1- X02/BAXIS2)
                 TERM2 = SQRT(TERM2)
                 Z1 = ZM + TERM2
                 Z2 = ZM - TERM2
              ENDIF

C***        Intersection points with sphere
            ELSEIF (SHAPE .EQ. 'SPHERE') THEN
              XM = DSPHERE * COSD * SINA
              YM = DSPHERE * SIND
              ZM = DSPHERE * COSD * COSA
              DX2 = (X0-XM)**2
              DY2 = (Y0-YM)**2
              TERM2 = RSPHERE2 - DX2 - DY2
              IF (TERM2 .LE. .0) THEN
                 Z1 = .0
                 Z2 = .0
              ELSE
                 TERM2 = SQRT(TERM2)
                 Z1 = ZM + TERM2
                 Z2 = ZM - TERM2
              ENDIF
      
            ELSE
              STOP '*** Internal ERROR in SECONDMODEL_PREP: unknown SHAPE'
            ENDIF

C**         Clipping the intersection line at the RMAX sphere
            ZMAX2 = RMAX2 - X02 - Y02
            ZMAX2 = MAX (.0, ZMAX2)
            ZMAX = SQRT (ZMAX2)
            IF (Z2 .GE. ZMAX .OR. Z1 .LE. -ZMAX) THEN
C**            Intersection interval is entirely outside RMAX sphere
               Z1 = .0
               Z2 = .0
            ELSE
               Z1 = MIN (Z1,  ZMAX)
               Z2 = MAX (Z2, -ZMAX)
            ENDIF

C**         Core rays:
            IF (P(JP) .LT. 1.) THEN
C**            Intersection lines behind the stellar disc are obscured
               ZMIN = SQRT (1. - P(JP)*P(JP))
               IF (Z1 .LE. ZMIN) THEN
                  Z1 = .0
                  Z2 = .0
               ELSE
                  Z2 = ZMIN
               ENDIF
            ENDIF

C***        First angle point is always kept
            IF (LPHI_TEMP .LE. 1) THEN
               ZINTER(1, JP, 1) = Z1
               ZINTER(2, JP, 1) = Z2
               PHIARR(1,JP) = PHI_VEC_TEMP(1)
               LPHI = 1
               BPHI_ORIG_LAST = .TRUE.

            ELSE
C***        Only if last OR current phi intersects -> increase counter; 
C***        else overwrite last point
C***        "original" angle points are also kept and not overwritten
C***        additionally, first non-intersection points are not overwritten
              DZ = Z2 - Z1
              DZLAST = ZINTER(2, JP, LPHI) - ZINTER(1, JP, LPHI)
              IF (DZ .NE. .0 .OR. DZLAST .NE. .0 .OR. BPHI_ORIG_LAST) 
     >           LPHI = LPHI + 1
              ZINTER(1, JP, LPHI) = Z1
              ZINTER(2, JP, LPHI) = Z2
              PHIARR(LPHI,JP) = PHI_VEC_TEMP(LPHI_TEMP)
C**           If last point was inside cone, do not overwrite it
              BPHI_ORIG_LAST= BPHI_ORIG(LPHI_TEMP) .OR. (DZLAST .NE. .0)
            ENDIF

         ENDDO ! phi-loop, index LPHI_TEMP -----------------------------

         NPHI(JP) = LPHI
         IF (NPHI(JP) .GT. MAXNPHI) THEN
           MAXNPHI = NPHI(JP)
           MAXJP   = JP
         ENDIF
         NPHISUM = NPHISUM + NPHI(JP)

      ENDDO ! p-loop --------------------------------------------------

C***  Output of data cube for visualization with gnuplot
      OPEN (20, FILE='secondmodel.dat',
     >      STATUS='UNKNOWN')
      DO J=1, NP
         DO LPHI = 1, NPHI(J)
            X0 = P(J)/P(NP) * COS(PHIARR(LPHI,J))
            Y0 = P(J)/P(NP) * SIN(PHIARR(LPHI,J))
            IF (ZINTER(1,J,LPHI) .NE. .0 .OR.
     >          ZINTER(1,J,LPHI) .NE. .0) THEN
                  WRITE (20, *) X0, Y0, ZINTER(1,J,LPHI)/P(NP)
                  WRITE (20, *) X0, Y0, ZINTER(2,J,LPHI)/P(NP)
                  WRITE (20, *) '   '
                  WRITE (20, *) '   '
            ENDIF
         ENDDO
      ENDDO
      CLOSE (20)

C***  plot of the secondmodel's geometrical (phi,p) grid
      OPEN (66, FILE='secondmodel.plot', STATUS='UNKNOWN')

      CALL PLOT_SECONDMODEL_GRID (P, NP, NPDIM, NPHI, PHIARR,
     >      NPHIMAX, JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG, 
     >      ZINTER)

      CLOSE (66)

C***  Statistical output

      WRITE (*,'(A)') 'SECONDMODEL invoked; specifications:'
      WRITE (*,'(A,$)') 'SHAPE=', SHAPE
      IF (SHAPE .EQ. 'CONE') WRITE (*, '(A,F6.1,A,F6.1)') 
     >         '  THETA=', THETA_DEG, '  CONEI=', CONEI_DEG
      IF (SHAPE .EQ. 'SPHERE') WRITE (*, '(4(A,F6.1))') 
     >         '  RSPHERE=', RSPHERE, '  DSPHERE=', DSPHERE,
     >         '  DELTASPHERE=', DELTASPHERE_DEG, 
     >         '  ALPHASPHERE=', ALPHASPHERE_DEG

      WRITE (0,'(/,A,A)') 
     >       'SECONDMODEL invoked; parameters: SHAPE=', SHAPE
      WRITE (0,'(A,I4,A,I4,A,F10.3)')
     > 'Maximum number of phi angles:', MAXNPHI,
     >     ' at P(', MAXJP, ') =', P(MAXJP)

      NPHI_MEAN = NINT(FLOAT(NPHISUM)/FLOAT(NP))
      WRITE(0,'(A,I4,/)')
     > 'Average number of phi angles per impact parameter:', NPHI_MEAN

C***  Calculate angle integration weights
C***    Note: Normalization of PHIWEIGHT will be done in program FORMAL
cccc    note: this part is identical in ROTATION_PREP 
      DO JP=1, NP-1
        IF (NPHI(JP) .EQ. 1) THEN
            PHIWEIGHT(1,JP) = 1.
        ELSE 
C***       Calculate the angle integration weights by trapezoidal rule
C***       but omitting the factors 0.5 at all weights
           PHIWEIGHT(1,      JP) = PHIARR(2,JP) - PHIARR(1,JP)
           DO LPHI=2, NPHI(JP)-1
             PHIWEIGHT(LPHI,JP) = PHIARR(LPHI+1,JP) - PHIARR(LPHI-1,JP)
           ENDDO
           PHIWEIGHT(NPHI(JP),JP) = PHIARR(NPHI(JP),JP)
     >                           - PHIARR(NPHI(JP)-1,JP)
        ENDIF

      ENDDO


      RETURN

C*********************************************************************
C***  ERROR branches  ************************************************
C*********************************************************************

   90 WRITE (0,*) '*** ERROR: Option ', ACTPAR(:IDX(ACTPAR)), 
     >            ' needs a value (keyword)'
      GOTO 99

   91 WRITE (0,*) '*** ERROR: Invalid keyword after option SHAPE: ', 
     >            SHAPE(:IDX(SHAPE))
      GOTO 99

   92 WRITE (0,*) '*** ERROR: mandatory parameter SHAPE is missing'
      GOTO 99

   93 WRITE (0,*) '*** ERROR: mandatory parameter THETA is missing'
      GOTO 99

   94 WRITE (0,*) '*** ERROR: mandatory parameter CONEI is missing'
      GOTO 99

   95 WRITE (0,*) '*** ERROR: mandatory parameter DELTA is missing'
      GOTO 99

   96 WRITE (0,*) '*** ERROR: mandatory parameter ALPHA is missing'
      GOTO 99

   97 WRITE (0,*) '*** ERROR: mandatory parameter RSPHERE is missing'
      GOTO 99

  971 WRITE (0,*) '*** ERROR: mandatory parameter DSPHERE is missing'
      GOTO 99

  972 WRITE (0,*) '*** ERROR: SPHERE covers the whole atmosphere'
      GOTO 99

   98 WRITE (0,*) '*** ERROR: Parameter cannot be decoded as number:'
      WRITE (0,*) '*** ERROR: ', ACTPAR(:IDX(ACTPAR)), '=', 
     >                           ACTPAR2(:IDX(ACTPAR2))
      GOTO 99


   99 STOP '*** FATAL ERROR in subroutine SECONDMODEL_PREP'

      END
      SUBROUTINE SET_POP_ZERO(LEVEL, NDIM, N, SPZ1, SPZ2, POPNUM, 
     >                        ND, NMOD, MAXSPZ)
C************************************************************************
C***  Some selected Popnums are set to Zero 
C***  Called by FORMALCL
C************************************************************************


      DIMENSION POPNUM(ND,N,NMOD)
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ)
      CHARACTER*10 LEVEL(NDIM)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT   = 6    !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR   = 0    !write to wruniqX.cpr (stderr)

      INTEGER, EXTERNAL :: IDX

      DO ISPZ = 1, MAXSPZ

      ID_SPZ1 = IDX(SPZ1(ISPZ))
      ID_SPZ2 = IDX(SPZ2(ISPZ))
      IF (ID_SPZ1 <= 0) RETURN

c      write (0,*) 'SPZ1=', SPZ1(ISPZ)
c      write (0,*) 'SPZ2=', SPZ2(ISPZ)


      DO I=1, N
        IF ((LEVEL(I)(1:ID_SPZ1) /= SPZ1(ISPZ)(1:ID_SPZ1)) .OR.
     >      (ID_SPZ2 > 0 .AND.  
     >       LEVEL(I)(1:ID_SPZ2) == SPZ2(ISPZ)(1:ID_SPZ2))) THEN
          CYCLE
        ENDIF
        WRITE (0,'(A,A,A)') 'SET_POP_ZERO: Popnums for Level ',
     >                              LEVEL(I), ' have been nulled.'
        DO L=1, ND
          DO M=1, NMOD
            POPNUM(L,I,M) = 0.
          ENDDO
        ENDDO
      ENDDO

      ENDDO ! loop over all input cards

      RETURN
      END
      SUBROUTINE SETUP (L,A,B,C,W,JMAX,ND,NP,NPDIM,OPA,ETA,THOMSON,
     $          XLAM,Z,RADIUS,BCORE,DBDR,XIMINUS,ENTOT,K, IVERS)
C***********************************************************************
C***  ANGLE-DEPENDENT CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY: 
C***  SET UP THE COEFFICINT MATRICES FOR A FEAUTRIER SCHEME:
C***     A (DIAGONAL), B (FULL), C (DIAGONAL) AND W (VECTOR)
C***  CALLED FROM: SUBROUTINE ELIMIN
C***********************************************************************

      DIMENSION A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM)
      DIMENSION RADIUS(ND),OPA(ND),ETA(ND),THOMSON(ND)
      DIMENSION Z(ND,NP)
      LOGICAL PLOT

      JMAX=NP+1-L
      JMM=JMAX-1
 
C***  EVERY L = 1 ... ND
      X=OPA(L)
      G=-X*THOMSON(L)
      ETAL=ETA(L)
 
C***  MEAN INTENSITY INTEGRATION WEIGHTS FROM SUBROUTINE MOMENT0 (VEKTOR W)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,W,DUMMY,.TRUE.)
      DO 1 J=1,JMAX
      WJG=W(J)*G
      DO 1 JS=1,JMAX
    1 B(JS,J)=WJG
      DO 3 J=1,JMAX
    3 W(J)=ETAL
 
      IF(L.EQ.1) GOTO 9
      IF(L.EQ.ND) GOTO 10
 
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1
      XP=(X+OPA(L+1))/2.
      XM=(X+OPA(L-1))/2.
      DO 2 J=1,JMM
      ZLPLUS=Z(L+1,J)
      ZLJ=Z(L,J)
      ZLMIN=Z(L-1,J)
      DT=2./(ZLMIN-ZLPLUS)
      DTM=XM*(ZLMIN-ZLJ)
      DTP=XP*(ZLJ-ZLPLUS)
      A(J)=DT/DTM
      C(J)=DT/DTP
    2 B(J,J)=B(J,J)+A(J)+C(J)+X
 
C     LAST ROW OF BLOCK, J=JMAX
      ZLMIN=Z(L-1,JMAX)
      DT=ZLMIN*XM
      A(JMAX)=2.*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      RETURN
 
C***  OUTER BOUNDARY CONDITION (AT L=1) 
    9 CONTINUE
      
C***  DIFFERENT VERSIONS FOR IMINUS
      IF (IVERS .EQ. 0) THEN
         XIMINUS=0.

      ELSE IF (IVERS .EQ. 1) THEN
         TAU = RADIUS(1) * X
         SBOUND = ETAL / X
         XIMINUS = SBOUND * (1.-EXP(-TAU))
   
      ELSE IF (IVERS .EQ. 4) THEN
C***     TESTPLOT FACILITY: DELETE "C" + SET "K" !!!
C***     PROGRAM             WRCONT: K = CONTINUUM INDEX
C***               ETL, CMF, FORMAL: K = LINE INDEX
         PLOT=.FALSE.
C         PLOT = K .EQ. 13
         CALL SFIT (ND, OPA, ETA, ENTOT, XLAM, SBOUND, PLOT)
         TAU = RADIUS(1) * X
         IF (TAU .LT. 1.E-3) THEN
            XIMINUS = SBOUND * TAU / 3.
         ELSE
            EXPTB = EXP(-TAU)
            FAK1 = 2./TAU
            FAK2 = FAK1/TAU
            XIMINUS = SBOUND * (1.-FAK1+FAK2-FAK2*EXPTB)
         ENDIF

      ELSE
         PRINT *,'INVALID VERSION OF OUTER BOUNDARY CONDITION'
         STOP 'ERROR'
         ENDIF

      XP=(X+OPA(2))/2.

      DO 8 J=1,JMM
      ZLPLUS=Z(2,J)
      ZLJ=Z(1,J)
      DT=XP*(ZLJ-ZLPLUS)
C***  MODIFICATION FOR NONZERO INCIDENT RADIATION  FROM TRUNCATED LAYERS
      W(J)=ETAL + XIMINUS * X *2./DT
      C(J)=2.*X/DT/DT
      B(JMAX,J)=.0
    8 B(J,J)=B(J,J)+C(J)+X+2.*X/DT
      B(JMAX,JMAX)=X
      W(JMAX)=XIMINUS * X
      RETURN
 
C***  INNER BOUNDARY CONDITION    L = ND
   10 XM=(X+OPA(ND-1))/2.
      DO 14 J=1,JMM
      ZLMIN=Z(ND-1,J)
      ZLJ=Z(ND,J)
      DT=XM*(ZLMIN-ZLJ)
      A(J)=2.*X/DT/DT
      B(J,J)=B(J,J)+A(J)+X+2.*X/DT
      PLUSI=BCORE+DBDR*ZLJ/X
   14 W(J)=ETAL+PLUSI*2.*X/DT
      A(JMAX)=2.*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      W(JMAX)=ETAL
      RETURN
      END
C******************************************************************************
      SUBROUTINE SFIT (ND, OPA, ETA, ENTOT, XLAM, SBOUND, PLOT)
C******************************************************************************
C***  THIS SUBROUTINE PREPARES SBOUND = SOURCE FUNCTION AT OUTER BOUNDARY,
C***  OBTAINED BY A SMOOTH EXTRAPOLATION. AT THE OUTERMOST MBOUND
C***  DEPTH POINTS, THE LINE SOURCE FUNCTION IS CONVERTED INTO RADIATION
C***  TEMPERATURES, WHICH ARE THEN FITTED BY A WEIGHTED LEAST-SQUARE FIT
C***  POLYNOMIAL. 
C***  INDEPENDETN VATIABLE: LOG ENTOT(L) - LOG ENTOT(L=1)
C***  WEIGHTS: AS APPROPRIATE FOR AN INTEGRATION (TRAPEZOIDAL RULE). THUS 
C***  THE CROWDED POINTS CLOSE TO THE BOUNDARY ARE NOT HEAVILY WEIGHTED.
C***
C***  THE APPLIED SUBROUTINES E02ADF AND E02AEF ARE FROM THE NAG FORTRAN LIBRARY
C***  THE FOLLOWING PARAMETERS MUST BE SUITABLE CHOOSEN:
C***     MBOUND = NUMBER OF DEPTH POINTS ACCOUNTED FOR IN THE EXTRAPOLATION
C***     KPLUS1 = DEGREE OF FIT POLYNOMIAL + 1 (E.G. QUADRATIC ... KPLUS1=3)
C******************************************************************************

      PARAMETER ( MBOUND = 30 )
      PARAMETER ( KPLUS1 =  4 )
C***  Note: Other degrees are not supported by SUBR. POLYFIT !

      DIMENSION OPA(ND), ETA(ND), ENTOT(ND)
      DIMENSION XFIT(MBOUND),YFIT(MBOUND),WFIT(MBOUND),SDEV(KPLUS1)
      DIMENSION AFIT(KPLUS1,KPLUS1), WORK1(3*MBOUND),WORK2(2*KPLUS1)
      DIMENSION A(KPLUS1,KPLUS1), B(KPLUS1)
      LOGICAL PLOT

      DIMENSION ATEST(KPLUS1,KPLUS1), BTEST(KPLUS1), DTEST(KPLUS1)
      DIMENSION SCRATCH(2*KPLUS1), XKOEFF(KPLUS1)

C***  Operating system:                    
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      IF (MBOUND .GT. ND) THEN
         WRITE (0,*) 'MBOUND EXCEEDS ND'
         WRITE (0,*) 'MBOUND, ND=', MBOUND, ND
         PRINT *, 'MBOUND EXCEEDS ND'
         STOP 'ERROR'
         ENDIF

C***  DEFINE THE TABLE OF (X,Y)-VALUES
      A1=ALOG10(ENTOT(1))
      DO 10 L=1,MBOUND
      XFIT(L)=ALOG10(ENTOT(L))-A1
      IF (OPA(L) .LE. .0) THEN
         YFIT(L)=.0
         ELSE
         S=ETA(L)/OPA(L)
         TRAD=TRADFUN(XLAM,S)
         YFIT(L)=TRAD
         ENDIF
   10    CONTINUE

C***  DEFINE WEIGHTS
      WFIT(1)=XFIT(2)-XFIT(1)
      DO 11 L=2,MBOUND-1
      WFIT(L)=XFIT(L+1)-XFIT(L-1)
   11 CONTINUE
      WFIT(MBOUND)=XFIT(MBOUND)-XFIT(MBOUND-1)

C***  Former CRAY branch deleted wrh  1-Jun-2021
        CALL POLYFIT (XFIT, YFIT, WFIT, MBOUND, KPLUS1,
     >                XKOEFF, A, B, SCRATCH, ATEST, BTEST, DTEST)

C***   As the X scale is defined such that X=0 for RMAX, 
C***      the value of the polynomial is simply the last coefficient
ccc        CALL HORNER (XFIT(1), TRADBOUND, KPLUS1, XKOEFF)
        TRADBOUND = XKOEFF(KPLUS1)

C***  If the fit polynomial has a minimum inside the interval, 
C***     the minimum value is taken instead of the value at x=0
C***  Note: This branch is also restricted to KPLUS1=4 !!
cccC***    First test: negative derivative at x=0
ccc        IF (XKOEFF(3) .GT. .0) GOTO 16
        IF (XKOEFF(KPLUS1) .LE. 0.) GOTO 16
ccc     Fix added: Skip for XKOEFF(1) == 0    (ansander, 06-Oct-2016)
        IF (ABS(XKOEFF(1)) <= 1.E-20) GOTO 16
        PHALB =    XKOEFF(2) / (3.*XKOEFF(1))
        Q     =    XKOEFF(3) / (3.*XKOEFF(1))
C***    Second test: Existence of (real) solution for f'=0
        WURZ = PHALB*PHALB - Q
        IF (WURZ .LT. .0) GOTO 16 
        WURZ = SQRT(WURZ)

C***    First solution of quadratic equation
        X1 = -PHALB + WURZ
C***    Test if X1 inside interval 0, X(MBOUND)
        IF (X1 .LT. .0) GOTO 15
        IF (X1 .GT. XFIT(MBOUND)) GOTO 15
C***    Test, if minimum, not maximum (f" > 0)
        F2STRICH = 2. * XKOEFF(2) + 6. * XKOEFF(1) * X1  
        IF (F2STRICH .LT. .0) GOTO 15
C***    X1 is minimum!        
        CALL HORNER (X1, TRADMIN, KPLUS1, XKOEFF)
        IF (TRADMIN .LT. TRADBOUND) TRADBOUND = TRADMIN

   15   CONTINUE
C***    Now the same for second solution of quadratic equation
        X1 = -phalb - WURZ
C***    Test if X1 inside interval 0, X(MBOUND)
        IF (X1 .LT. .0) GOTO 16
        IF (X1 .GT. XFIT(MBOUND)) GOTO 16
C***    Test, if minimum, not maximum (f" > 0)
        F2STRICH = 2. * XKOEFF(2) + 6. * XKOEFF(1) * X1  
        IF (F2STRICH .LT. .0) GOTO 16
C***    X1 is minimum!        
        CALL HORNER (X1, TRADMIN, KPLUS1, XKOEFF)
        IF (TRADMIN .LT. TRADBOUND) TRADBOUND = TRADMIN

   16   CONTINUE

      SBOUND=BNUE(XLAM,TRADBOUND)

C***  TEST PLOT, IF REQUESTED:
      IF (PLOT) CALL SFITPLO (ENTOT,ETA,OPA,ND,MBOUND,AFIT,
     >                        KPLUS1,XLAM,XKOEFF, TRADBOUND)

      RETURN
      END
      SUBROUTINE SFITPLO (ENTOT,ETAL,OPAL,ND,MBOUND,AFIT,KPLUS1,XLAM, 
     >                    XKOEFF, TRADBOUND)

C************************************************************************
C***  TEST PLOT OF SELECTED LINE SOURCE FUNCTION (RADIATION TEMPERATURE)
C***  VERSUS NUMBER DENSITY (AS DEPTH PARAMETER)
C***  AND COMPARISON WITH THE INTERPOLATION POLYNOM CALCULATED FOR
C***  VERSION 4 OF THE OUTER BOUNDARY CONDITIONS
C***  CALLED FROM: SUBROUTINE SFIT
C************************************************************************

      DIMENSION ENTOT(ND),ETAL(ND),OPAL(ND)
      DIMENSION AFIT(KPLUS1,KPLUS1), XKOEFF(KPLUS1)
      DIMENSION X(100), Y(100), A(20)
C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      IF (ND .GT. 100) THEN
         CALL REMARK ('INSUFFICIENT DIMENSION')
         PRINT *, 'INSUFFICIENT DIMENSION'
         STOP 'ERROR'
         ENDIF

      DO 1 L=1, ND
      X(L)=ALOG10(ENTOT(L))
      IF (OPAL(L) .LE. 0) THEN
         TRAD=.0
         ELSE
         S=ETAL(L)/OPAL(L)
         TRAD=TRADFUN(XLAM,S)
         ENDIF
      Y(L)=TRAD/1000.
    1 CONTINUE

      XMIN=7.
      XMAX=16.
      XSCALE=2.5
      XTICK=1.
      XABST=3.

      YMIN=0.
      YMAX=80.
      YSCALE=0.
      YTICK=5.
      YABST=10.

      OPEN (2, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
      CALL PLOTANF (2
     $ ,'TESTPLOT: SOURCE FUNCTION AND FIT CURVE'
     $ ,'TESTPLOT: SOURCE FUNCTION AND FIT CURVE'
     $ ,'LOG OF NUMBER DENSITY / (CM**-3)'
     $ ,'TRAD / kK'
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,0.
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,0.
     $ ,X,Y,ND,2)

C***  INTERPOLATION POLYNOM
      IFAIL=0
C***  COEFFICIENT VECTOR A
      DO 3 I=1,KPLUS1
      A(I)=AFIT(KPLUS1,I)
    3 CONTINUE
      X1=ALOG10(ENTOT(1))
      XM=ALOG10(ENTOT(MBOUND))
      DO 2 L=1,MBOUND
      XL=ALOG10(ENTOT(L))
      XX= XL - X1

      CALL HORNER (XX, TRAD, KPLUS1, XKOEFF)

      Y(L)=TRAD/1000.
    2 CONTINUE

      CALL PLOTCON (2,X,Y,MBOUND,5)

      Y(1) = TRADBOUND / 1000.
      CALL PLOTCONS (2,X1,Y(1),1, 'SYMBOL=8 COLOR=2')

      RETURN
      END
      SUBROUTINE SHIFT (R,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***********************************************************************
      DIMENSION R(ND)

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFT, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO 1 II=L,ND
      I=ND+L-II
    1 R(I+1)=R(I)

      RETURN
      END
	FUNCTION SOFBET(B,P,N,M)
C
C  Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for
C  quasistatic charged particles.  The alpha and beta lines of the first
C  three series are explicitly included. All other cases use the H18 
C  profile. Profiles are normalised to full oscillator strength. Method 
C  is based on Griem (1960, ApJ 132, 883).
C
C  By Deane Peterson and Bob Kurucz.
C
C  STORAGE FOR CORRECTIONS (P,BETA,IND),(P,IND),(P,IND)
C

      DIMENSION PROPBM(5,15,7),C(5,7),D(5,7)
      DIMENSION PP(5),BETA(15)
      DIMENSION PROB1(75),PROB2(75),PROB3(75),PROB4(75),PROB5(75)
      DIMENSION PROB6(75),PROB7(75)
      DIMENSION C1(5),C2(5),C3(5),C4(5),C5(5),C6(5),C7(5)
      DIMENSION D1(5),D2(5),D3(5),D4(5),D5(5),D6(5),D7(5)
      EQUIVALENCE (PROPBM(1,1,1),PROB1(1)),(PROPBM(1,1,2),PROB2(1))
      EQUIVALENCE (PROPBM(1,1,3),PROB3(1)),(PROPBM(1,1,4),PROB4(1))
      EQUIVALENCE (PROPBM(1,1,5),PROB5(1)),(PROPBM(1,1,6),PROB6(1))
      EQUIVALENCE (PROPBM(1,1,7),PROB7(1))
      EQUIVALENCE (C(1,1),C1(1)),(C(1,2),C2(1)),(C(1,3),C3(1))
      EQUIVALENCE (C(1,4),C4(1)),(C(1,5),C5(1)),(C(1,6),C6(1))
      EQUIVALENCE (C(1,7),C7(1))
      EQUIVALENCE (D(1,1),D1(1)),(D(1,2),D2(1)),(D(1,3),D3(1))
      EQUIVALENCE (D(1,4),D4(1)),(D(1,5),D5(1)),(D(1,6),D6(1))
      EQUIVALENCE (D(1,7),D7(1))
      SAVE PROPBM,C,D,PP,BETA
C
C  Lyman alpha
C

      DATA PROB1/
     1-.980,-.967,-.948,-.918,-.873,-.968,-.949,-.921,-.879,-.821,
     2-.950,-.922,-.883,-.830,-.764,-.922,-.881,-.830,-.770,-.706,
     3-.877,-.823,-.763,-.706,-.660,-.806,-.741,-.682,-.640,-.625,
     4-.691,-.628,-.588,-.577,-.599,-.511,-.482,-.484,-.514,-.568,
     5-.265,-.318,-.382,-.455,-.531,-.013,-.167,-.292,-.394,-.478,
     6 .166,-.056,-.216,-.332,-.415, .251, .035,-.122,-.237,-.320,
     7 .221, .059,-.068,-.168,-.247, .160, .055,-.037,-.118,-.189,
     8 .110, .043,-.022,-.085,-.147/
      DATA C1 /-18.396, 84.674,-96.273,  3.927, 55.191/
      DATA D1 / 11.801,  9.079, -0.651,-11.071,-26.545/
C
C  Lyman beta
C

      DATA PROB2/
     1-.242, .060, .379, .671, .894, .022, .314, .569, .746, .818,
     2 .273, .473, .605, .651, .607, .432, .484, .489, .442, .343,
     3 .434, .366, .294, .204, .091, .304, .184, .079,-.025,-.135,
     4 .167, .035,-.082,-.189,-.290, .085,-.061,-.183,-.287,-.374,
     5 .032,-.127,-.249,-.344,-.418,-.024,-.167,-.275,-.357,-.420,
     6-.061,-.170,-.257,-.327,-.384,-.047,-.124,-.192,-.252,-.306,
     7-.043,-.092,-.142,-.190,-.238,-.038,-.070,-.107,-.146,-.187,
     8-.030,-.049,-.075,-.106,-.140/
      DATA C2 / 95.740, 18.489, 14.902, 24.466, 42.456/
      DATA D2 / -6.665, -7.136,-10.605,-15.882,-23.632/
C
C  Balmer alpha
C

      DATA PROB3/
     1-.484,-.336,-.206,-.111,-.058,-.364,-.264,-.192,-.154,-.144,
     2-.299,-.268,-.250,-.244,-.246,-.319,-.333,-.337,-.336,-.337,
     3-.397,-.414,-.415,-.413,-.420,-.456,-.455,-.451,-.456,-.478,
     4-.446,-.441,-.446,-.469,-.512,-.358,-.381,-.415,-.463,-.522,
     5-.214,-.288,-.360,-.432,-.503,-.063,-.196,-.304,-.394,-.468,
     6 .063,-.108,-.237,-.334,-.409, .151,-.019,-.148,-.245,-.319,
     7 .149, .016,-.091,-.177,-.246, .115, .023,-.056,-.126,-.189,
     8 .078, .021,-.036,-.091,-.145/
      DATA C3 /-25.088,145.882,-50.165,  7.902, 51.003/
      DATA D3 /  7.872,  5.592, -2.716,-12.180,-25.661/
C
C  Balmer beta
C

      DATA PROB4/
     1-.082, .163, .417, .649, .829, .096, .316, .515, .660, .729,
     2 .242, .393, .505, .556, .534, .320, .373, .394, .369, .290,
     3 .308, .274, .226, .152, .048, .232, .141, .052,-.046,-.154,
     4 .148, .020,-.094,-.200,-.299, .083,-.070,-.195,-.299,-.385,
     5 .031,-.130,-.253,-.348,-.422,-.023,-.167,-.276,-.359,-.423,
     6-.053,-.165,-.254,-.326,-.384,-.038,-.119,-.190,-.251,-.306,
     7-.034,-.088,-.140,-.190,-.239,-.032,-.066,-.103,-.144,-.186,
     8-.027,-.048,-.075,-.106,-.142/
      DATA C4 / 93.783, 10.066,  9.224, 20.685, 40.136/
      DATA D4 / -5.918, -6.501,-10.130,-15.588,-23.570/
C
C  Paschen alpha
C

      DATA PROB5/
     1-.819,-.759,-.689,-.612,-.529,-.770,-.707,-.638,-.567,-.498,
     2-.721,-.659,-.595,-.537,-.488,-.671,-.617,-.566,-.524,-.497,
     3-.622,-.582,-.547,-.523,-.516,-.570,-.545,-.526,-.521,-.537,
     4-.503,-.495,-.496,-.514,-.551,-.397,-.418,-.448,-.492,-.547,
     5-.246,-.315,-.384,-.453,-.522,-.080,-.210,-.316,-.406,-.481,
     6 .068,-.107,-.239,-.340,-.418, .177,-.006,-.143,-.246,-.324,
     7 .184, .035,-.082,-.174,-.249, .146, .042,-.046,-.123,-.190,
     8 .103, .036,-.027,-.088,-.146/
      DATA C5 /-19.819, 94.981,-79.606,  3.159, 52.106/
      DATA D5 / 10.938,  8.028, -1.267,-11.375,-26.047/
C
C  Paschen beta
C

      DATA PROB6/
     1-.073, .169, .415, .636, .809, .102, .311, .499, .639, .710,
     2 .232, .372, .479, .531, .514, .294, .349, .374, .354, .279,
     3 .278, .253, .212, .142, .040, .215, .130, .044,-.051,-.158,
     4 .141, .015,-.097,-.202,-.300, .080,-.072,-.196,-.299,-.385,
     5 .029,-.130,-.252,-.347,-.421,-.022,-.166,-.275,-.359,-.423,
     6-.050,-.164,-.253,-.325,-.384,-.035,-.118,-.189,-.252,-.306,
     7-.032,-.087,-.139,-.190,-.240,-.029,-.064,-.102,-.143,-.185,
     8-.025,-.046,-.074,-.106,-.142/
      DATA C6 /111.107, 11.910,  9.857, 21.371, 41.006/
      DATA D6 / -5.899, -6.381,-10.044,-15.574,-23.644/
C
C  Balmer 18
C

      DATA PROB7/
     1 .005, .128, .260, .389, .504, .004, .109, .220, .318, .389,
     2-.007, .079, .162, .222, .244,-.018, .041, .089, .106, .080,
     3-.026,-.003, .003,-.023,-.086,-.025,-.048,-.087,-.148,-.234,
     4-.008,-.085,-.165,-.251,-.343, .018,-.111,-.223,-.321,-.407,
     5 .032,-.130,-.255,-.354,-.431, .014,-.148,-.269,-.359,-.427,
     6-.005,-.140,-.243,-.323,-.386, .005,-.095,-.178,-.248,-.307,
     7-.002,-.068,-.129,-.187,-.241,-.007,-.049,-.094,-.139,-.186,
     8-.010,-.036,-.067,-.103,-.143/
      DATA C7 /511.318,  1.532,  4.044, 19.266, 41.812/
      DATA D7 / -6.070, -4.528, -8.759,-14.984,-23.956/
      DATA PP/0.,.2,.4,.6,.8/
      DATA BETA/1.,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.310,
     1          7.943,10.,12.59,15.85,19.95,25.12/
C

      IF(B.GT.500.) THEN
C
C  Very large B
C

        B2=B*B
        SOFBET=(1.5/SQRT(B)+27./B2)/B2
        RETURN
      END IF
C
C  Other cases
C

      CORR=1.
      B2=B*B
      SB=SQRT(B)
      INDX=7
      MMN=M-N
      IF(N.LE.3.AND.MMN.LE.2) INDX=2*(N-1)+MMN
C
C  Determine relevant Debye range
C

      IM=MIN(INT(5.*P)+1,4)
      IP=IM+1
      WTPP=5.*(P-PP(IM))
      WTPM=1.-WTPP
      IF(B.LE.25.12) THEN
C

        JP=2
   1    IF(B.GT.BETA(JP).AND.JP.LT.15) THEN
          JP=JP+1
          GO TO 1
        END IF
        JM=JP-1
C

        WTBP=(B-BETA(JM))/(BETA(JP)-BETA(JM))
        WTBM=1.-WTBP
        CBP=PROPBM(IP,JP,INDX)*WTPP+PROPBM(IM,JP,INDX)*WTPM
        CBM=PROPBM(IP,JM,INDX)*WTPP+PROPBM(IM,JM,INDX)*WTPM
        CORR=1.+CBP*WTBP+CBM*WTBM
C
C  Get approximate profile for the inner part
C

        PR1=0.
        PR2=0.
        WT=AMAX1(MIN(0.5*(10.-B),1.),0.)
        IF(B.LE.10.) PR1=8./(83.+(2.+0.95*B2)*B)
        IF(B.GE.8.)  PR2=(1.5/SB+27./B2)/B2
        SOFBET=(PR1*WT+PR2*(1.-WT))*CORR
      ELSE
C
C  Asymptotic part for medium B's
C

        CC=C(IP,INDX)*WTPP+C(IM,INDX)*WTPM
        DD=D(IP,INDX)*WTPP+D(IM,INDX)*WTPM
        CORR=1.+DD/(CC+B*SB)
        SOFBET=(1.5/SB+27./B2)/B2*CORR
      END IF
C

      RETURN
      END

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
      SUBROUTINE SPLINPO_FAST (F, X, FI, XI, N, LPARAM, SAFE)
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
      LOGICAL SAFE
C***  For using the entry point, all variables must be static!
      SAVE

C***  CHECK FOR STRICTLY MONOTONIC ORDER
C***  (SKIPPED IF UNSAVE-MODE IS REQUESTED) 
      IF (SAFE) THEN
        DN = XI(N) - XI(1)
        DO L=2, N
           DX = XI(L) - XI(L-1)
           IF (DX*DN .LE. .0) THEN
              WRITE (0, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3         FORMAT (' *** ERROR IN SUBROUTINE SPLINPO_FAST:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
              STOP 'ERROR'
           ENDIF
        ENDDO

C***  FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
        DO I=2, N
           L = I
           IF ( (X-XI(I-1)) * (X-XI(I)) .LE. .0) GOTO 2
        ENDDO
        STOP ' *** ERROR IN SUBR. SPLINPO - X OUTSIDE TABLE'
    2   CONTINUE
        IF (LPARAM .NE. L) THEN
           WRITE (0,'(A, I5)') '**** LPARAM =', LPARAM
           WRITE (0,'(A, I5)') '**** CORRECT L =', L
           STOP ' *** ERROR IN SUBR. SPLINPO - WRONG INDEX LPARAM'
        ENDIF
      ENDIF

C***  In the unsafe mode, current index is taken from the call
      L = LPARAM

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

C***  Entry point for a subsequent call with same x point, 
C***     but different function 
      ENTRY SPLINPO_FAST_SAME_X (F, X, FI, XI, N, SAFE)

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
c            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
c         ELSE
c            S3 = 1.4 * S4 - 0.5 * F3
c         ENDIF
c         IF (LB .EQ. L+1) THEN
c            S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
c         ELSE
c            S5 = 1.5 * S4 - 0.5 * F4
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
      SUBROUTINE STARKBROAD (KODATIND, NOM, PHITAB,
     >    NFDIMPHITAB, NLDIMPHITAB, NLPHITAB, ND, RADIUS, 
     >    OPAC, OPAL, LINPRO, AVOIGT, NBL, MAXLAP, IPOINTERPHITAB,
     >    XMAX, XMAXBROAD, XMAXLIN, NDDIM, DXMAX, PATH_VCSSB, LOW, NUP,
     >    PHISCRATCH, XLAMNBL, ALN, T, ENTOT, RNE, LEVEL, MAINQN, 
     >    NCHARG, POPNUM, N, VDOP, ELEVEL, EION, PATH_LEMKE_DAT,
     >    DD_VDOP, DD_VDOPDU, IND_ORIGLEV, BDD_VDOP,
     >    DD_VMICDU, GRIEMPAR, TAUMAX, TAUMINBROAD, IMOD, MAXMOD, 
     >    EINST, NDIM)


C***********************************************************************      
C***  This subroutine prepares line broadening. 
C***  For H I and He II lines, the profiles 
C***  for each depth point are stored in array PHITAB.
C***  In ZONEINT, the profile fuction is interpolated from that table.
C***  
C***  For all other lines, the VOIGT function is used, 
C***  and the depth-dependent parameter AVOIGT is prepared here.
C***
C***  called from FORMAL, separately for each line (NBL) in the blend
C***********************************************************************      

      DIMENSION PHITAB(-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB)
      CHARACTER(8) :: LINPRO
      CHARACTER*256 PATH_VCSSB, PATH_LEMKE_DAT
      
      INTEGER, INTENT(IN) :: NDDIM
C* take only vector!
      REAL, DIMENSION(ND), INTENT(IN) :: DD_VDOP, DD_VDOPDU, DD_VMICDU
      DIMENSION KODATIND(2), NOM(2), MAINQN(2), NCHARG(2) 
      DIMENSION ELEVEL(2), EION(2)
      DIMENSION T(ND), ENTOT(ND), RNE(ND)
      DIMENSION IND_ORIGLEV(MAXLAP)
      DIMENSION POPNUM(NDDIM,N)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD), GRIEMPAR(MAXLAP,NDDIM,MAXMOD)
      DIMENSION EINST(NDIM,NDIM)

      CHARACTER*10 LEVEL(NDIM)
C***  Arrays for the Lemke tables 
      PARAMETER (NWS_DIM = 54)
      PARAMETER (NTS_DIM = 7)
      PARAMETER (NES_DIM = 17)
      PARAMETER (NPS_DIM = NWS_DIM * NTS_DIM * NES_DIM)
      DIMENSION STKTA_DWS(NWS_DIM), STKTA_TS (NTS_DIM)
      DIMENSION STKTA_ES (NES_DIM), STKTA_PS (NPS_DIM)
      
      LOGICAL BRESONANCE, BDD_VDOP, bHYDROGEN
C***  1 / SQRT(PI)
      DATA WPIINV / 0.564189583549 /
C***  4 * PI
      DATA PI4 / 12.56637062 /
C***  Clight in km/s
      DATA CLIGHT / 299792.458 /
C***  Prepare the profile functions for pressure broadening 
        
C***  Line is resonance line if: 
C**     1) mother lower level is the first to appear in DATOM, 
C***  or   
C***    2) level previous to mother lower level is of a different ion,
C***  or
C***    3) level previous to mother lower level is of different element

      IF (IND_ORIGLEV(LOW) .EQ. 1) THEN
        BRESONANCE = .TRUE.
      ELSE
        BRESONANCE = 
     >          (NCHARG(IND_ORIGLEV(LOW)) .NE. NCHARG(IND_ORIGLEV(LOW)-1)) 
     >          .OR.
     >          (NOM(IND_ORIGLEV(LOW)) .NE. NOM(IND_ORIGLEV(LOW)-1))
      ENDIF 
      
      bHYDROGEN = .FALSE.

C***  Pre-set AVOIGT with radiation-damping 
C***    (will not be used in all cases)

C***  IF NO EXPLICITE VALUE OF GAMMA: SUM OVER EINSTEIN COEFFICIENTS
      IF (AVOIGT(NBL,1,1) .LT. 0.) THEN
         AV=.0
C***     UPPER LEVEL:
         DO I=1, NUP-1
            IF (EINST(NUP,I) .GT. 0.) AV=AV+EINST(NUP,I)
         ENDDO
C***     LOWER LEVEL:
         DO I=1, LOW-1
            IF (EINST(LOW,I) .GT. 0.) AV=AV+EINST(LOW,I)
         ENDDO

         AVOIGT(NBL,1,1) = AV
      ENDIF     

C***  TRANSFORMATION:  GAMMA ---> A (VOIGT FUNCTION)
C***  Note: AVOIGT becomes depth dependent if VDOP is depth dependent.
      AV = AVOIGT(NBL,1,1)
      DO L=1, ND
         DNUED = 1.E13 * DD_VDOP(L) / XLAMNBL
            AVOIGT(NBL,L,IMOD) = AV / PI4 / DNUED
      ENDDO

C***  Choose the appropriate broadening mode LINPRO

C***  If LINPRO has been set already to VOIGT -> keep that setting
      IF (LINPRO .EQ. 'VOIGT   ') THEN

C***  If LINPRO='DRTRANS' was set as default by subr. DRTRANS
      ELSE IF (LINPRO .EQ. 'DRTRANS ') THEN
         LINPRO = ''
      ELSE IF (KODATIND(NOM(LOW)) .EQ. 1) THEN
C***     HYDROGEN BROADENING        
         LINPRO = 'BRD-H   '
         bHYDROGEN = .TRUE.
      ELSE IF (KODATIND(NOM(LOW)) .EQ. 2) THEN
C***     HELIUM BROADENING
         IF (NCHARG(LOW) .EQ. 0) THEN
            LINPRO = 'BRD-HeI '
         ELSE IF (NCHARG(LOW) .EQ. 1) THEN
            LINPRO = 'BRD-HeII'
         ELSE
            STOP '*** INTERNAL ERROR IN STARKPROF: NCHARG'
         ENDIF
      ELSE
C***     ALL OTHER ELEMENTS
         NELECTRON = KODATIND(NOM(LOW)) - NCHARG(LOW)
C***     Linear stark broadening for hydrogenic ions
c           (only possible if MAINQN is different!)
         IF (((NELECTRON == 1) .OR. (NELECTRON == 3)
     >        .OR. (NELECTRON == 11) .OR. (NELECTRON == 19))
     >        .AND. MAINQN(LOW) /= MAINQN(NUP)) THEN      
             LINPRO = 'L-STARK '

C***  Quadratic stark broadening, except for resonance lines  
         ELSEIF (.NOT. BRESONANCE) THEN   
            LINPRO = 'Q-STARK '
         ELSEIF (BRESONANCE) THEN   
            LINPRO = 'VOIGT   '
         ENDIF
      ENDIF
      
C******************************************************************
C********* H I ***************************************************
C******************************************************************
      IF (LINPRO .EQ. 'BRD-H   ') THEN
         NLPHITAB = NLPHITAB + 1
         IPOINTERPHITAB = NLPHITAB
         IF (NLPHITAB .GT. NLDIMPHITAB) THEN
            WRITE (0,*) '*** DIMENSION FOR LINE-BROADENING ',
     >                  ' TABLES (H and He) INSUFFICIENT'
            WRITE (0,'(A,I4)') 
     >                  'Present value: NLDIMPHITAB=', NLDIMPHITAB
            STOP ' *** FATAL ERROR'
         ENDIF

C***     Read Lemke table for the current line
         CALL READ_H_STARKDATA 
     >       (PATH_LEMKE_DAT, MAINQN(LOW), MAINQN(NUP), LEMKEFOUND,
     >       STKTA_DWS, NWS, STKTA_TS, NTS,
     >       STKTA_ES,  NES, STKTA_PS, NPS,
     >       NWS_DIM, NTS_DIM, NES_DIM, NPS_DIM)

         DO L=1, ND
            XNE = ENTOT(L) * RNE(L)

C***        Setting the profile function
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               X = K * DXMAX
               XLAMK = XLAMNBL * EXP(ALN*X)
               IF (LEMKEFOUND .EQ. 1) THEN
                   PHITAB(K, L, NLPHITAB) = STARK_HI_LEMKE
     >              (XLAMK, XLAMNBL, T(L), XNE,
     >               STKTA_DWS, NWS, STKTA_TS, NTS,
     >               STKTA_ES,  NES, STKTA_PS)
               ELSE
                  PHITAB(K, L, NLPHITAB) = STARKHI
     >             (MAINQN(LOW), MAINQN(NUP), XLAMK, XLAMNBL, T(L), XNE,
     >              LEVEL(LOW), LEVEL(NUP)) 
               ENDIF 
            ENDDO

C***        Normalization
            SUM = .0
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               SUM = SUM + PHITAB(K,L,NLPHITAB)
            ENDDO
            FAKNORM = 1. / (DXMAX * SUM)
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               PHITAB(K,L,NLPHITAB) = PHITAB(K,L,NLPHITAB)*FAKNORM
            ENDDO

C***   We realized that the Lemke tables already include thermal broadening. 
C***   Therefore, the routine VDOP_STRUCT now provides the array 
C***   DD_VMICDU (vmic) in addition to DD_VDOPDU. ~~Tomer, 7.5.2015
C***   CONVOLGAUUS_FLEX convolves with gaussian corresponding to an arbitrary VDOP
            NDAT = 2 * NFDIMPHITAB + 1
            DD_VMICDU_SQRD = DD_VMICDU(L)**2
            CALL CONVOLGAUSS_FLEX (PHITAB(-NFDIMPHITAB,L,NLPHITAB), 
     >          PHISCRATCH, NDAT, DXMAX, DD_VDOPDU(L), DD_VMICDU_SQRD)        
         ENDDO

C******************************************************************
C********* He I ***************************************************
C******************************************************************
       ELSE IF (LINPRO .EQ. 'BRD-HeI ') THEN
          
C***   Add Stark broadening by electron impact
          DO L=1, ND
             DLAMD = XLAMNBL * DD_VDOP(L) / CLIGHT
             XNE = ENTOT(L) * RNE(L)
             CALL STARKDAMP_HEI (GAMMAHE1, SHIFT, T(L), XNE, 
     >                 XLAMNBL, NUP, LOW, LINPRO, LEVEL)
             IF (LINPRO .EQ. 'Q-STARK ') EXIT
C***            Add this GAMMA to the "a" (=Voigt function parameter)
                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + GAMMAHE1/DLAMD
             ENDDO

ccc          !!! the follwing part gives nonsense results; 
ccc          !!! Is the returned broadening parameter HEWID in diffent units???  
ccc          !!! So far, I must omit this effect         wrh 11-Jul-2008 
ccc          Attention: in case of a SECOND MODEL, POPNUM would be needed to
ccc            be scaled to te second RADIUS grid
c
cccC***      Add Stark broadening by neutral-helium impact
ccc             DO L=1, ND
cccC***            Number density of neutral helium
ccc                SUM = .0
ccc                DO J = 1, N
ccc                   IF (KODATIND(NOM(J)) .EQ. 2) THEN
ccc                      IF (NCHARG(J) .EQ. 0) SUM = SUM + POPNUM(L,J)
ccc                   ENDIF
ccc                ENDDO
ccc                HE1FRC = ENTOT(L) * SUM
ccc                CALL STARKDAMP_HEI_NEUTRAL (NUP, LOW, T(L),
ccc     >                              HE1FRC, LINPRO, HEWID,
ccc     >                              MAINQN, LEVEL)
ccc                IF (LINPRO .EQ. 'VOIGT   ') EXIT
ccc                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + HEWID/DLAMD
ccc             ENDDO

             
C******************************************************************
C********* He II **************************************************
C******************************************************************
          ELSE IF (LINPRO .EQ. 'BRD-HeII') THEN
             CALL STARKHEIIPREP 
     >          (MAINQN, NUP, LOW, LINPRO, LEVEL, PATH_VCSSB)
C***         if line not in broadening table, LINPRO was changed to L-STARK
             IF (LINPRO .EQ. 'BRD-HeII') THEN
                NLPHITAB = NLPHITAB + 1
                IPOINTERPHITAB = NLPHITAB
                IF (NLPHITAB .GT. NLDIMPHITAB) THEN
                   WRITE (0,*) '*** DIMENSION FOR LINE-BROADENING ',
     >                      ' TABLES (H and He) INSUFFICIENT'
                   WRITE (0,'(A,I4)') 
     >                      'Present value: NLDIMPHITAB=', NLDIMPHITAB
                   STOP ' *** FATAL ERROR'
                ENDIF

                DO L=1, ND
                   XNE = ENTOT(L) * RNE(L)

C***            Setting the profile function
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                      X = K * DXMAX
                      XLAMK = XLAMNBL * EXP(ALN*X)
                      PHITAB(K, L, NLPHITAB) = STARKHEII
     >                        (XLAMK, XLAMNBL, T(L), XNE) 
                   ENDDO


C***            Normalization
                   SUM = .0
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                      SUM = SUM + PHITAB(K,L,NLPHITAB)
                   ENDDO
                   FAKNORM = 1. / (DXMAX * SUM)
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                       PHITAB(K,L,NLPHITAB) = 
     >                      PHITAB(K,L,NLPHITAB)*FAKNORM
                   ENDDO

C***               CONVOLGAUUS_FLEX convolves with gaussian corresponding 
C***               to an arbitrary VDOP
                   NDAT = 2 * NFDIMPHITAB + 1
                   DD_VDOPDU_SQRD = DD_VDOPDU(L)**2
                   CALL CONVOLGAUSS_FLEX
     >                  (PHITAB(-NFDIMPHITAB,L,NLPHITAB), PHISCRATCH, 
     >                   NDAT, DXMAX, DD_VDOPDU(L), DD_VDOPDU_SQRD)  
                ENDDO
             ENDIF         
          ENDIF

C******************************************************************
C********* all other elements beside H and He: quadratic Stark  ***
C***       also fallback branch for He I if no table found
C******************************************************************
          IF (LINPRO == 'L-STARK ') THEN
C***         Linear Stark broadening for hydrogenic ions          
C            GRIEMPAR is prepared here
             DO L=1, ND
                XNE = ENTOT(L) * RNE(L)
                CALL LINSTARK (GRIEMPAR(NBL,L,IMOD), XNE, 
     >             KODATIND(NOM(LOW)),
     >             MAINQN(LOW), MAINQN(NUP), XLAMNBL, 
     >             LINPRO, LEVEL(LOW), LEVEL(NUP))
             ENDDO
          ENDIF
          
C***      Q-STARK is also fallback if no QNs are available for L-STARK
          IF (LINPRO .EQ. 'Q-STARK ') THEN

C***     Add Stark broadening by electron impact (Quadratic Stark Effect)
             DO L=1, ND
                DNUED = 1.E13 * DD_VDOP(L) / XLAMNBL
                XNE = ENTOT(L) * RNE(L)
                CALL QUADSTARK (GAMMAQUAD, T(L), XNE, 
     >             NCHARG(NUP), ELEVEL(NUP), EION(NUP), LINPRO, 
     >             LEVEL(LOW), LEVEL(NUP))
                IF (LINPRO .EQ. 'VOIGT   ') EXIT
C***            Add this GAMMA to the "a" (=Voigt function parameter)
                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + GAMMAQUAD/DNUED
             ENDDO
C***         No preset of VOIGT for DRTRANSIT lines
             IF (NCHARG(LOW) .NE. NCHARG(NUP)) LINPRO = ''
          ENDIF

C***  Bandwidth XMAXBROAD is only calculated for the first MODEL:
      IF (LINPRO .NE. '' .AND. IMOD .EQ. 1) 
     >    CALL BANDWIDTH (ND, RADIUS, OPAC, OPAL, LINPRO, AVOIGT, NDDIM, 
     >                NBL, MAXLAP, XMAX, XMAXBROAD, XMAXLIN,
     >                PHITAB(-NFDIMPHITAB, 1, NLPHITAB), 
     >                NFDIMPHITAB, DXMAX, LEVEL(NUP), LEVEL(LOW),
     >                BDD_VDOP, GRIEMPAR, VDOP, ALN, 
     >                bHYDROGEN, TAUMAX, TAUMINBROAD)

      RETURN
      END

	SUBROUTINE STARKDAMP_HEI (GAMMAHE1, SHIFT, T, XNE, XLAM, 
     >                NUP, LOW, LINPRO, LEVEL)
C ***
C ******************************************************************************
C ***      
C ***	ABSORPTIONSKOEFFIZIENT ISOLIERTER HELIUMLINIEN.	
C ***   NACH GRIEM, PHYS.REV. 125,177 (1962)
C ***	STARK-PARAMETER NACH BENNETT/GRIEM, UNIVERS. MARYLAND REPORT (71)
C ***
C ***   T       > R	Temperature [K]
C ***   XNE     > R	Electron density [cm**-3]
C ***   K       > I	Line index, see below
C ***   XLAM    > R	Central wavelength of line [е]
C ***
C ***   Modifications:
C ***
C ***   28-OCT-1997 22:58:34.11		ML / Stegaurach.  Copied 4026 & 4471
C ***                                   from BCSS.FOR.  Verify correctness...
C ***   23-May-2008 --> only gam4 berechnung / 
C***             short fast version without bcss data tables
C ***   		older version linhe1.f
C ***
C ******************************************************************************

      CHARACTER LINPRO*8, LEVEL(2)*10
      REAL TEM(4), ALF(4,16), W(4,16), SHF(4,16), EP(16)
      REAL BK, AN0, HM, PI, HC, C

      parameter( bk = 1.38046E-16, an0 = 8.85255E-13,
     >             hm = 1.6732E-24, pi = 3.14159, hc = 1.986182E-16,
     >             c = 2.997929E+10 )

      real c1, c2, c3, c4, c5, c6, c7, c8  ! DATA doesn't allow expressions
      parameter( c1 = 0.172*5.62, c2 = 0.193*5.62,
     >             c3 = 0.218*5.62, c4 = 0.249*5.62,
     >             c5 = 0.107*5.62, c6 = 0.119*5.62,
     >             c7 = 0.134*5.62, c8 = 0.154*5.62 ) 

      DATA TEM / 5000., 10000., 20000., 40000. /

******************************************************************
C**** LINES CALCULATED
C**** K = 1    LINE 4922    NOT USED, SEE HELIUM5
C**** K = 2    LINE 4713
C**** K = 3    LINE 4437
C**** K = 4    LINE 4388    NOT USED SEE HELIUM2
C**** K = 5    LINE 4121
C**** K = 6    LINE 3964    BLENDED WITH H. CALCULATED ONLY IF H.LT. .01
C**** K = 7    LINE 3888
C**** K = 8    LINE 5875    CALCULATED IF KHEM .GT. 7
C**** K = 9    LINE 6678    AS 5875
C**** K = 10   LINE 7065    AS 5875
C**** K = 11   LINE 5015   
C**** K = 12   LINE 5048
C**** K = 13   LINE 2945    CALCULTED IF KHEM .GT. 12
C**** K = 14   LINE 2829    AS 2945
C *** K = 15   LINE 4026    Griem 1974, copied from BCSS.FOR -- shift ignored
C *** K = 16   LINE 4471    BCS 1974, JQSRT 14, 1025; BCSS.FOR -- no shift
C ***                       forbidden component ignored.
C***********************************************************************

C***  In order to find the proper "line index" K, we will compare 
C***  XLAM with the approximate wavelengths:
      DIMENSION XLAMK(16)
      DATA XLAMK / 4922., 4713., 4437., 4388., 
     2             4121., 3964., 3888., 5875., 
     3             6678., 7065., 5015., 5048., 
     4             2945., 2829., 4026., 4471.  / 

      DATA ALF / 
     1           0.683, 0.773, 0.885, 1.023,
     2		 0.115, 0.103, 0.095, 0.092,
     3		 0.199, 0.184, 0.177, 0.179,
     4		 1.159, 1.321, 1.527, 1.782,
     5		 0.171, 0.155, 0.145, 0.141,
     6		 0.275, 0.290, 0.311, 0.341,
     7	         0.075, 0.070, 0.067, 0.067,
     8		 0.064, 0.061, 0.059, 0.059,
     9		 0.146, 0.157, 0.169, 0.181,
     Z		 0.067, 0.060, 0.055, 0.052,
     1           0.154, 0.160, 0.169, 0.180,
     2           0.135, 0.123, 0.117, 0.116,
     3		 0.204, 0.195, 0.195, 0.203,
     4		 0.285, 0.276, 0.279, 0.294,
     >           c1, c2, c3, c4,		! see PARAMETER above
     >           c5, c6, c7, c8/

      DATA W / 
     1         2.30,  1.96,  1.63,  1.35,
     2	       0.343, 0.394, 0.438, 0.459,
     3	       1.41,  1.57,  1.65,  1.62,
     4	       6.13,  5.15,  4.24,  3.45,
     5	       0.787, 0.900, 0.987, 1.02,
     6	       1.030, 0.966, 0.877, 0.776,
     7	       0.102, 0.112, 0.117, 0.117,
     8	       0.159, 0.170, 0.176, 0.177,
     9	       0.423, 0.386, 0.349, 0.318,
     Z	       0.180, 0.208, 0.235, 0.254,
     1         0.378, 0.359, 0.334, 0.306,
     2         0.625, 0.704, 0.755, 0.760,
     3	       0.808, 0.857, 0.857, 0.812,
     4	       1.790, 1.870, 1.840, 1.720,
     5         4.040, 3.490, 2.960, 2.470,
     6         1.460, 1.269, 1.079, 0.898/

      DATA SHF / 
     1            1.020,  0.772,  0.584,  0.439,
     2		  0.403,  0.416,  0.391,  0.335,
     3		  1.51,   1.43,    1.24,  0.995,
     4		  2.52,   1.87,    1.38,   1.01,
     5		  0.900,  0.989,  0.811,  0.673,
     6		 -0.647,  -.504,  -.374,  -.266,
     7		  .0743,  .0603,  .0463,  .0347,
     8		 -.0880, -.0552, -.0255, -.0050,
     9		  0.275,  0.233,  0.196,  0.161,
     Z		  0.215,  0.231,  0.227,  0.203,
     1           -0.250, -0.200, -0.152, -0.111,
     2            0.700,  0.685,  0.611,  0.504,
     3		  0.521,  0.411,  0.311,  0.231,
     4		  1.070,  0.830,  0.620,  0.455,
     5            0.000,  0.000,  0.000,  0.000,
     6            0.000,  0.000,  0.000,  0.000/

      DATA EP / 171135., 169087., 171135., 171135., 169087., 166278.,
     A		159856., 169087., 171135., 169087., 166278., 171135.,
     A          159856., 159856., 169087., 169087./

C***  Find K index
      DO KK = 1, 16
         IF (XLAMK(KK) .LT. 0.999*XLAM) CYCLE 
         IF (XLAMK(KK) .GT. 1.001*XLAM) CYCLE 
         K = KK
         GOTO 1
      ENDDO

      WRITE (*,11)  LEVEL(LOW), LEVEL(NUP) 
      WRITE (0,11)  LEVEL(LOW), LEVEL(NUP) 
   11 FORMAT ('*** WARNING from subr. STARKDAMP_HEI: ', 
     >        'no tabulated data for ', A, ' - ', A, 
     >        '--> using Q-STARK instead')

      GAMMAHE1 = 0.0
      LINPRO   = 'Q-STARK '
      RETURN

    1 CONTINUE

ccc test output
ccc      WRITE (0,'(A,2F8.2)') '*** STARKDAMP_HEI: ', XLAMK(K),  XLAM 

C**** Interpolation der Daten
      CALL NEWPOL2(W(1,K), TEM, WE, T, 4 )
      WE = WE * XNE * 1.E-16
      CALL NEWPOL2(ALF(1,K), TEM, ALFI, T, 4 )
      ALFI = ALFI * 1.E-4 * XNE**0.25
      CALL NEWPOL2(SHF(1,K), TEM, D, T, 4 )
      D = D * XNE * 1.E-16
      VV = SQRT(8./PI*BK*T*(1./(4.*HM) + 1./HM/16.))	
      SIG = WE * (0.75/PI/XNE)**(1./3.) / VV / 1.E-8 * C / (XLAM**2)

C**** GESAMTDAEMPFUNG UND VERSCHIEBUNG
      B = ALFI**(8./9.) * SIG**(-1./3.)
      GAMMAHE1 = WE*(1. + 1.36*B)

      IF (D .LT. 0.)  B = - B
      SHIFT = D + WE*2.36*B

C***  Beziehung zwischen gam4 und gamrad
C***  gam = gammahe1 + gamrad/4./pi  
C***  a = gam/dop      

      RETURN
      END
      FUNCTION STARKHEII (WAVE, WAVEH, T, XNE)
C**********************************************************************
C***  Called from: SUBROUTINE STARKBROAD
C***  Stark-Broadening for helium (He II) lines
C***  The data table for the current line must have been loaded before 
C***     to the common block VCSSBDAT by Subr. STARKHEIIPREP
C***  Normalization is arbitrary
C**********************************************************************

      PARAMETER (MPDIM     = 50)     ! Max. # OF PROFILE POINTS IN vcssb TABLE
      PARAMETER (MTEMPDIM  =  6)     ! Max. # OF TEMP POINTS IN vcssb TABLE
      PARAMETER (MNEDIM    = 11)     ! Max. # OF DENSITY POINTS IN vcssb TABLE

      COMMON /VCSSBDAT/ MP, MNE, MTEMP, ALOG_NEL(MNEDIM),
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC,
     >                   SVCS(MPDIM,MTEMPDIM,MNEDIM)

C***  Definition of DEX function
      EXP10( X ) = EXP( 2.30258509299405 * X )

      PHI = .0

      ALOG_T = ALOG10( T )
      XNE0    = EXP10( ALOG_NEL(1) )

C***  Interpolation interval, weights in temperature grid
      BTEMP = ( ALOG_T - ALOG_T0 ) / ALOG_T_INC + 1.
      ITEMP = MAX0 ( MIN0(INT(BTEMP), MTEMP-1), 1 )
      WTTEMP = BTEMP - FLOAT(ITEMP)

C***  Interpolation interval, weights in electron-density grid
      ALOG_NE = ALOG10( XNE )
      DO INE = 1, MNE
         IF (ALOG_NEL(INE) .GT. ALOG_NE) EXIT
         INELAST = INE
      END DO

      INE = INELAST - 1
      INE = MAX0 (MIN0(INE, MNE-1), 1)
      WTNE = (ALOG_NE - ALOG_NEL(INE))/(ALOG_NEL(INE+1) - ALOG_NEL(INE))

      IF( WTNE .LT. 0. ) WTNE = 0.    ! NO EXTRAPOLATION
      IF( WTNE .GT. 1. ) WTNE = 1.

C***  Find wavelength interval index: (ialp, ialp+1)
      DL = ABS(WAVE - WAVEH)
C***  Center point
      IF (DL .LT. 1.E-10) THEN
         IALP = 1
         WT   = .0
      ELSE
         F0 = 1.25E-9 * AMAX1(XNE,XNE0)**(2./3.)   ! FIX AT TABLE LIMIT
         ALOG_ALPHA = ALOG10( DL/F0 )
         BALP = (ALOG_ALPHA - ALOG_ALPHA0) / ALOG_ALPHA_INC + 1.
         IALP = INT(BALP)
         IALP = MAX0 (IALP,1)
         WT = BALP - FLOAT(IALP)
      ENDIF

C***  Set profile to zero outside the tabulated wavelength range
      IF (IALP .GE. MP) GOTO 10

      PROFI = (1. - WTNE) * (1. - WTTEMP) * SVCS(IALP,ITEMP  ,INE  ) +
     +        (1. - WTNE) * WTTEMP        * SVCS(IALP,ITEMP+1,INE  ) +
     +           WTNE     * (1. - WTTEMP) * SVCS(IALP,ITEMP  ,INE+1) +
     +           WTNE     * WTTEMP        * SVCS(IALP,ITEMP+1,INE+1)

      PROFIP= (1. - WTNE) * (1. - WTTEMP) * SVCS(IALP+1,ITEMP  ,INE  ) +
     +        (1. - WTNE) * WTTEMP        * SVCS(IALP+1,ITEMP+1,INE  ) +       
     +           WTNE     * (1. - WTTEMP) * SVCS(IALP+1,ITEMP  ,INE+1) +
     +           WTNE     * WTTEMP        * SVCS(IALP+1,ITEMP+1,INE+1)

      PROF = (1. - WT) * PROFI + WT * PROFIP


      PHI = EXP10(PROF)


   10 STARKHEII = PHI
      RETURN
      END

      SUBROUTINE STARKHEIIPREP 
     >   (MAINQN, NUP, LOW, LINPRO, LEVEL, PATH_VCSSB)
C*******************************************************************
C***  Called from: SUBROUTINE STARKBROAD 
C***  This Subr. reads the line-broadening tables for He II
C***  The tables are searched for the transition low-nup
C*******************************************************************
      CHARACTER LINPRO*8
      CHARACTER*10 LEVEL(2)
      CHARACTER*256 FILENAME
      CHARACTER*256 PATH_VCSSB
      DIMENSION MAINQN(2)

      PARAMETER (MPDIM     = 50)     ! Max. # OF PROFILE POINTS IN vcssb TABLE
      PARAMETER (MTEMPDIM  =  6)     ! Max. # OF TEMP POINTS IN vcssb TABLE
      PARAMETER (MNEDIM    = 11)     ! Max. # OF DENSITY POINTS IN vcssb TABLE

      COMMON /VCSSBDAT/ MP, MNE, MTEMP, ALOG_NEL(MNEDIM), 
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC, 
     >                   SVCS(MPDIM,MTEMPDIM,MNEDIM)

      IF (PATH_VCSSB .EQ. 'default') THEN  
         FILENAME = '/home/corona/wrh/work/wrdata/VCSSB.DAT'
      ELSE
         FILENAME = PATH_VCSSB(:IDX(PATH_VCSSB)) // '/' // 'VCSSB.DAT'
      ENDIF
      KANAL=13

      MAINQNLOW = MAINQN(LOW)
      MAINQNNUP = MAINQN(NUP)
C***  Check if Principle Quantum numbers are known
      IF (MAINQNLOW .LE. 0) GOTO 93
      IF (MAINQNNUP .LE. 0) GOTO 94

      OPEN(KANAL, FILE=FILENAME(:IDX(FILENAME)), 
     >  STATUS = 'OLD', ERR=99 )

C***  Loop over the lines ***************************************
      DO 
        READ(KANAL, *, ERR=95, END=98) NL, NU
        READ(KANAL, *, ERR=95, END=96) MNE, MTEMP, MP
C***    CHECK DIMENSIONS
	IF( MP     .GT. MPDIM    .OR.
     >      MTEMP  .GT. MTEMPDIM .OR.
     >      MNE    .GT. MNEDIM   ) GOTO 97

        READ (KANAL,*,ERR=95) (ALOG_NEL(J), J = 1, MNE),
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC

        READ (KANAL,*,ERR=95) (( (SVCS(K,M,J), K = 1, MP),
     >                           M = 1, MTEMP), J = 1, MNE)
C***    Leave loop if the line was found, otherwise read next line

cc        write (0,'(A,4I3)') 'NL, LOW, NU, NUP =', NL, LOW, NU, NUP

        IF (NL .EQ. MAINQNLOW .AND. NU .EQ. MAINQNNUP) GOTO 1
      END DO

   1  CONTINUE

C***  Regular exit: line data found
      GOTO 101

C***  Error branches *******************************************

   93 WRITE (0,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(LOW)
      WRITE (*,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(LOW)
      LINPRO = 'Q-STARK '
      GOTO 100

   94 WRITE (0,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(NUP)
      WRITE (*,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(NUP)
      LINPRO = 'Q-STARK '
      GOTO 100

   95 WRITE (0,*) '*** READ ERROR on file VCSSB.DAT'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'

   96 WRITE (0,*) '*** ERROR: E-O-F reached prematurely on VCSSB table'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'


   97 WRITE (0,*) '*** ERROR: INSUFFICIENT DIMENSION FOR VCSSB TABLE'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'

   98 WRITE (0,'(6A)') '*** WARNING: no Stark broadening table found ',
     >   'for ', LEVEL(LOW), ' - ', LEVEL(NUP),
     >   '  -> using L-STARK approximation instead'
      WRITE (*,*) '*** WARNING: no Stark broadening table found ',
     >   'for ', LEVEL(LOW), ' - ', LEVEL(NUP), 
     >   '  -> using L-STARK approximation instead'
      LINPRO = 'L-STARK '
      GOTO 100

   99 WRITE (0,*) '*** ERROR: No Stark data found for He II '
      WRITE (0,*) '*** Cannot open file: ', FILENAME(:IDX(FILENAME))
      WRITE (0,*) '*** You may want to define a different path by'
      WRITE (0,*) '*** using the FORMAL_CARDS option PATH_VCSSB = ...'
      WRITE (*,*) '*** ERROR: No Stark data found for He II '
      WRITE (*,*) '*** Cannot open file: ', FILENAME(:IDX(FILENAME))
      WRITE (*,*) '*** You may want to define a different path by'
      WRITE (*,*) '*** using the FORMAL_CARDS option PATH_VCSSB = ...'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'


  100 CONTINUE
  101 CLOSE (KANAL)

      RETURN
      END
  
      FUNCTION STARKHI(N,M,WAVE,WAVEH,T,XNE, LEVELLOW, LEVELNUP)
C********************************************************************
C     Stark-Broadening for hydrogen (H I) lines
C     XNE - electron density
C     M - upper , N - lower principle quantum number
C
C     The unmodified Code returns a stark-profile for lyman-series
C
C     the modified code returns stark-profiles for lyman, paschen and 
C     bracket-series
C     mario parade - !!! using a vcs-approximation-function vcse1f !
C
C     Returns the Stark broadened line profile.  The code follows Griem's 
C     theories (mostly 1960 ApJ, 132, 883) with corrections to approximate 
C     the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.
C
C     Area normalised to unity with frequency.
C
C     by Deane Peterson & Bob Kurucz.
C     (adapted, corrected and comments added by PB)
C************************************************************************

      REAL K
      DIMENSION Y1WTM(2,2),XKNMTB(4,3)
      LOGICAL LYMANALF
      CHARACTER*10 LEVELLOW, LEVELNUP
      SAVE

C  Knm constants as defined by Griem (1960, ApJ 132, 883) for the long 
C  range Holtsmark profile (due to ions only). Lyman and Balmer series 
C  are from VCS, higher series from elsewhere.

      DATA XKNMTB/0.0001716, 0.0090190, 0.1001000, 0.5820000,
     1            0.0005235, 0.0177200, 0.1710000, 0.8660000,
     2            0.0008912, 0.0250700, 0.2230000, 1.0200000/
  
C    new data - mario parade
C
C    DATA KNMTAB/.000356,	.000523,	.00109,		.00149,
C		 .00225,	.0125,		.0177,		.028,
C    		1.0348,		.0493,		.124,		.171,
C		 .223,		.261,		.342,		.683,
C		 .866,		1.02,		1.19,		1.46/
C     DATA FSTARK/.1387,.07910,.02126,.01394,.006462,.004814,.002779,
C    1 .002216,.001443,.001201,.3921,.1193,.03766,.02209,.01139,
C    2 .008036,.005007,.003850,.002658,.002151,.6103,.1506,.04931,
C    3 .02768,.01485,.01023,.006588,.004996,.003542,.002838,.8163,.1788,
C    4 .05985,.03189,.01762,.01196,.007825,.005882,.004233,.003375/

      DATA Y1WTM/1.E18, 1.E17, 1.E16, 1.E14/
      DATA N1/0/, M1/0/

      PARAMETER (CLIGHT = 2.9979258E18)
      PARAMETER (PI = 3.14159265359, SQRTPI = 1.77245385)
      PARAMETER (H = 6.62618E-27)  !Planck in cgs
      PARAMETER (K = 1.38066E-16)  !Boltzmann in cgs

C  Variables depending on conditions

C***  Check if Principle Quantum numbers are known
      IF (N .LE. 0) GOTO 93
      IF (M .LE. 0) GOTO 94
	
cc unused! wrh 19-Feb-2010      HE1FRC=1.0
      T4 = T/10000.
      T43 = T4**0.3
cc unused! wrh 19-Feb-2010     T3NHE = T43*HE1FRC ! fraction H/He
      XNE16 = XNE**0.1666667
cc      PRINT *,(XNE16)
      PP = XNE16*0.08989/SQRT(T) ! the shielding parameter  austausch F=-->F1 mario
      F1 = (XNE16*XNE16*XNE16*XNE16)*1.25E-9      ! Holtsmark normal field strength
cc      PRINT *,(PP)
      Y1B = 2./(1.+0.012/T*SQRT(XNE/T))
      Y1S = T43/XNE16
      C1D = F1*78940./ T
      C2D = F1**2/5.96E-23/XNE
      GCON1 = 0.2+0.09*SQRT(T4)/(1.+XNE/1.E13)
      GCON2 = 0.2/(1.+XNE/1.E15)
cc	PRINT *,(N),(M),(WAVE),(XNE16),(C1D),(C2D)
C

      DELW = WAVE-WAVEH
      FREQNM = CLIGHT/WAVEH
      FREQ = CLIGHT/WAVE
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C

      IF((N.NE.N1).OR.(M.NE.M1)) THEN  
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
C
C  Knm constants not tabulated from approximate asymptotic expression 
C  (Griem 1960 eqn 33) where 1./(1.+.13/FLOAT(MMN)) appears to be a 
C  correction factor to match to the tables.
C

         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5E-5/GNM*XMN2/(1.+.13/FLOAT(MMN))
         END IF
	
	
C
C  Some weighting factors which relate to y1, which is the velocity at 
C  which the minimum impact parameter (where second order perturbation 
C  theory breaks down) and the Lewis cutoff (the limit of validity of 
C  the impact approximation) are the same.
C

         IF(M.EQ.2) THEN
            Y1NUM = 550.
         ELSE IF (M.EQ.3) THEN
            Y1NUM = 380.
         ELSE
            Y1NUM = 320.
         END IF
         IF (MMN.LE.2 .AND. N.LE.2) THEN
            Y1WHT = Y1WTM(N,MMN)
         ELSE IF (MMN.LE.3) THEN
            Y1WHT = 1.E14
         ELSE
            Y1WHT = 1.E13
         END IF
C

         C1CON = XKNM/WAVEH*GNM*XM2MN2
         C2CON = (XKNM/WAVEH)**2
cc	 PRINT *,(C1CON),(C2CON)


      ENDIF

C  Compute line profile
C
C  PRQS is the quasistatic ion contribution
C  FNS  is the quasistatic electron contribution rel to PRQS
C  F    is the impact electron contribution
C
C  First compute the width of the impact electron profile roughly Griem
C  (1967, ApJ 147, 1092) eqn for w.
C

      WTY1 = 1./(1.+XNE/Y1WHT)
      Y1SCAL = Y1NUM*Y1S*WTY1+Y1B*(1.-WTY1)
      C1 = C1D*C1CON*Y1SCAL
      C2 = C2D*C2CON
      G1 = 6.77*SQRT(C1)
      BETA = ABS(DELW)/F1/XKNM
      Y1 = C1*BETA
      Y2 = C2*BETA**2

      IF ((Y2.LE.1.E-4).AND.(Y1.LE.1.E-5)) THEN
         GAM = G1*AMAX1(0.,0.2114+LOG(SQRT(C2)/C1))*(1.-GCON1-GCON2)
      ELSE
         GAM = G1*(0.5*EXP(-MIN(80.,Y1))+VCSE1F(Y1)-0.5*VCSE1F(Y2))*
     *            (1.-GCON1/(1.+(90.*Y1)**3)-GCON2/(1.+2000.*Y1))
         IF (GAM.LE.1.E-20) GAM = 0.
      END IF

C
C  Compute individual quasistatic and impact profiles.
C  PP - shielding parameter
C
	
      PRQS = SOFBET(BETA,PP,N,M) 
      IF (GAM.GT.0.) THEN
         F = GAM/PI/(GAM*GAM+BETA*BETA)
      ELSE
	 F = 0.
      ENDIF
C
C  Fraction of electrons which count as quasistatic. A fit to eqn 8 
C  (2nd term) of Griem (1967, ApJ 147, 1092).

      P1 = (0.9*Y1)**2
      FNS = (P1+0.03*SQRT(Y1))/(P1+1.)
C
C  DBETA (=dBeta/dfreq) changes the area normalisation. 
C  DSQRT(WAVE/WAVEH) corrects the long range part to dfreq**-5/2
C  asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).

cc      DBETA = CLIGHT/FREQ/FREQ/XKNM/F1
cc    richtiger mit Zentralwellenlaenge:
      DBETA = WAVEH/FREQNM/XKNM/F1
      STARKHI = (PRQS*(1.+FNS)+F)*DBETA * SQRT(WAVE/WAVEH)

C
C  The red wing is multiplied by the Boltzmann factor to roughly account
C  for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
C  absorption case.  If emission do for DEL.GT.0.
C

      IF (DEL.LT.0.) STARKHI = STARKHI * EXP(-ABS(H*DEL)/K/T)
C
      RETURN

C***  Error branches *******************************************

   93 WRITE (0,*) '*** ERROR: Missing MAINQN for lower level, '
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW      
      WRITE (*,*) '*** ERROR: Missing MAINQN for lower level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW
      GOTO 100

   94 WRITE (0,*) '*** WARNING: Missing MAINQN for upper level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW
      WRITE (*,*) '*** WARNING: Missing MAINQN for upper level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW

      GOTO 100

  100 STOP '*** FATAL ERROR detected by Subr. STARKHI'

      END

      FUNCTION STARK_HI_LEMKE (XLAM, XLAMCENTER, T, ED,
     >      STKTA_DWS, NWS, STKTA_TS, NTS, 
     >      STKTA_ES,  NES, STKTA_PS)       
C*****************************************************************
C     Function to obtain stark profile for a tabulated transition
C         as a function of electron density, temperature and wavelength
C input:
C     ED      : electron density per cubic centimeter
C     T       : temperature in K 
C     XLAM    : wavelength for which the profile should be calculated
C     XLAMCENTER : line center wavelength
C
C To prepare the work of this function, the broadening table for the 
C    considered line must have been stored into STKTA_PS by calling 
C    SUBROUTINE READ_H_STARKDATA (see there for more comments)
C The data array STKTA_PS has indices::   
C     NWS            Number of "scaled frequency"  points tabulated 
C     NTS            Number of temperature points tabulated 
C     NES            Number of electron densitiy points tabulated 
C
C output: stark profile at xlam in arbitrary units (!!)
C         - must be re-normalized after!
C******************************************************************

      DIMENSION STKTA_DWS(NWS), STKTA_TS (NTS), STKTA_ES (NES)
      DIMENSION STKTA_PS (NWS,NTS,NES)

      ELOG = ALOG10(ED)                                    
      TLOG = ALOG10(T) 
      DLAM = XLAM - XLAMCENTER

C***  Attention, only for symmetric profiles !!!!!!!
      DLAM_LOCAL = ABS(DLAM)

C***  The tables are over a re-scaled frequency 
C***   called "delta alpha" in the paper Lemke M.: 1997, A&A Suppl. 122, 285 
      DLAM_LOCAL = DLAM_LOCAL / (1.25E-9 * (10**(ELOG*2.0/3.0)))

C***  Find T interval (INDX_T, INDX_T+1)
      TLOG = MAX(STKTA_TS(1),MIN(STKTA_TS(NTS),TLOG))
      DO I = 2, NTS
         INDX_T=I-1
         IF (TLOG .LT. STKTA_TS(I)) EXIT
      END DO
      T_WGT=STKTA_TS(INDX_T+1)-STKTA_TS(INDX_T)
      T_WGT =(TLOG-STKTA_TS(INDX_T))/T_WGT

C***  Find electron density interval (INDX_ED, INDX_ED+1)
      ELOG = MAX(STKTA_ES(1),MIN(STKTA_ES(NES),ELOG))
      DO I = 2, NES                                     
         INDX_ED = I-1                                           
         IF (ELOG .LT. STKTA_ES(I)) EXIT
      END DO
      ED_WGT = (ELOG-STKTA_ES(INDX_ED))/
     >                 (STKTA_ES(INDX_ED+1)-STKTA_ES(INDX_ED))

C***  Find wavelength index interval (ML, ML+1)
      DLAM_LOCAL = MAX(STKTA_DWS(1),MIN(STKTA_DWS(NWS),DLAM_LOCAL))
      DO I = 2, NWS
         ML = I-1
         IF (DLAM_LOCAL .LT. STKTA_DWS(I)) EXIT
      END DO
      MLP = ML + 1
      WL_WGT = (DLAM_LOCAL - STKTA_DWS(ML)) 
     >       / (STKTA_DWS(MLP) - STKTA_DWS(ML))

C---------------------------------------------
C***  First at wavelength index ML
C***  Interpolation in temperature
      ST00 = STKTA_PS(ML,INDX_T  ,INDX_ED)
      ST01 = STKTA_PS(ML,INDX_T+1,INDX_ED)
      ST0  = T_WGT * (ST01 - ST00) + ST00

      ST00 = STKTA_PS(ML,INDX_T  ,INDX_ED+1)
      ST01 = STKTA_PS(ML,INDX_T+1,INDX_ED+1)
      ST1  = T_WGT * (ST01 - ST00) + ST00

C***  Interpolation in density
      STARK0 = ED_WGT*(ST1-ST0)+ST0
C---------------------------------------------
C***  second at wavelength index ML+1 alias MLP
C***  Interpolation in temperature
      ST00 = STKTA_PS(MLP,INDX_T  ,INDX_ED)
      ST01 = STKTA_PS(MLP,INDX_T+1,INDX_ED)
      ST0  = T_WGT * (ST01 - ST00) + ST00

      ST00 = STKTA_PS(MLP,INDX_T  ,INDX_ED+1)
      ST01 = STKTA_PS(MLP,INDX_T+1,INDX_ED+1)
      ST1  = T_WGT * (ST01 - ST00) + ST00

C***  Interpolation in density
      STARK1 = ED_WGT*(ST1-ST0)+ST0
C---------------------------------------------

C***  Interpolation in wavelength
      STARK = WL_WGT * (STARK1 - STARK0) + STARK0

      STARK_HI_LEMKE = 10.**STARK

      RETURN
      END
      FUNCTION STARKHOLTSMARK(XI, VDOP, VDOPDU_FINE, GRIEMPAR,
     >                        NLOC, MAXLAP, ZFINE, LRIP, RFINERIP, 
     >                        PR, QR, RADIUS, ND, PJPJ, bHYDROGEN)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***   that are broadened by the linear Stark effect
C****  
C***  GRIEMPAR = convertion factor for Griem's beta parameter 
C***              (has been multiplied with VDOP in STARKBROAD)
C***
C***  
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      INTEGER, INTENT(IN) :: ND

      REAL, DIMENSION(MAXLAP,ND) :: GRIEMPAR
      REAL, DIMENSION(ND) :: RADIUS
      
      REAL, INTENT(IN) :: VDOP, VDOPDU_FINE
      REAL, INTENT(INOUT) :: RFINERIP
      
      REAL :: BETA, BETADOP, FGRIEM, DLAMDOPREL
      
      LOGICAL :: bHYDROGEN

      REAL, PARAMETER :: CLIGHT  = 2.99792458E5    ! c in km / s
      
      
C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

         IF (RFINERIP .LE. RADIUS(LRIP)) THEN
            IF (RFINERIP .GE. RADIUS(LRIP+1)) THEN
               QR = (RADIUS(LRIP) - RFINERIP)/
     >              (RADIUS(LRIP) - RADIUS(LRIP+1))
               PR = 1. - QR
            ELSE
               LRIP = LRIP + 1
               GOTO 200
            ENDIF
         ELSE
            LRIP = LRIP - 1
            GOTO 200
         ENDIF
      ENDIF
C***  End of radius interpolation weights: PR, QR established

      FGRIEM = PR * GRIEMPAR(NLOC,LRIP) + QR * GRIEMPAR(NLOC,LRIP+1)
      
      DLAMDOPREL = VDOP / CLIGHT
      ALN = LOG(1. - DLAMDOPREL)
      BETA = FGRIEM * ABS( EXP(ALN*XI) - 1. )
      BETADOP = FGRIEM * DLAMDOPREL
      
      STARKHOLTSMARK = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)

      RETURN 
      END
      FUNCTION STARKPROF (XI, IPOINTERPHITAB, ZFINE, LRIP, RFINERIP, 
     >          PR, QR, PHITAB, NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, 
     >          NDDIM, PJPJ, DXMAX)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***  (used for H I and He II)
C***  The profile is obtained by interpolation in the table PHITAB
C***  Table PHITAB must have been prepared by SUBROUTINE STARKBROAD
C***  The index IPOINTERPHITAB is the pointer to the line index in that table  
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      DIMENSION PHITAB(-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB)
      DIMENSION RADIUS(ND)

C***  Prepare frequency interpolation
      XIK = XI / DXMAX
      K = INT(XIK)

C***  If frequency outside the bandwidth: PHI = 0.0
      IF (K .LE. -NFDIMPHITAB .OR. K .GE. NFDIMPHITAB) THEN
         STARKPROF = .0
         RETURN
      ENDIF

      Q = ABS( XIK - AINT(XIK) )
      IF (XI .LT. .0) THEN
         KK = K-1
      ELSE
         KK = K+1
      ENDIF

C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

ccc for tests only
ccc         IF (LRIP .LT. 1 .OR. LRIP .GT. ND) THEN
ccc            WRITE (0,*) 'RFINERIP=', RFINERIP
ccc            STOP '*** INTERNAL ERROR IN STARKPROF'
ccc         ENDIF


         IF (RFINERIP .LE. RADIUS(LRIP)) THEN
            IF (RFINERIP .GE. RADIUS(LRIP+1)) THEN
               QR = (RADIUS(LRIP) - RFINERIP)/
     >              (RADIUS(LRIP) - RADIUS(LRIP+1))
               PR = 1. - QR
            ELSE
               LRIP = LRIP + 1
               GOTO 200
            ENDIF
         ELSE
            LRIP = LRIP - 1
            GOTO 200
         ENDIF
      ENDIF
C***  End of radius interpolation weights: PR, QR established

      IP = IPOINTERPHITAB
      PHITABK   = PR * PHITAB(K ,LRIP,IP) + QR * PHITAB(K ,LRIP+1,IP)
      PHITABKK  = PR * PHITAB(KK,LRIP,IP) + QR * PHITAB(KK,LRIP+1,IP)

      STARKPROF = (1.-Q) * PHITABK + Q * PHITABKK

      RETURN
      END
      FUNCTION STARKVOIGT (XI, AVOIGT, NLOC, MAXLAP, 
     >     ZFINE, LRIP, RFINERIP, PR, QR, RADIUS, ND, PJPJ)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***   that are represented by a Voigt function with depth-dependent
C***   Voigt parameter AVOIGT
C***  The profile is obtained by interpolation in AVOIGT and call
C***  of the VOIGTH function. 
C***  Table AVOIGT must have been prepared by SUBROUTINE STARKBROAD
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      DIMENSION AVOIGT(MAXLAP,ND)
      DIMENSION RADIUS(ND)

C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

ccc for tests only
cc         IF (LRIP .LT. 1 .OR. LRIP .GT. ND) THEN
cc            WRITE (0,*) 'RFINERIP=', RFINERIP
cc            STOP '*** INTERNAL ERROR IN STARKPROF'
cc         ENDIF


         IF (RFINERIP .LE. RADIUS(LRIP)) THEN
            IF (RFINERIP .GE. RADIUS(LRIP+1)) THEN
               QR = (RADIUS(LRIP) - RFINERIP)/
     >              (RADIUS(LRIP) - RADIUS(LRIP+1))
               PR = 1. - QR
            ELSE
               LRIP = LRIP + 1
               GOTO 200
            ENDIF
         ELSE
            LRIP = LRIP - 1
            GOTO 200
         ENDIF
      ENDIF
C***  End of radius interpolation weights: PR, QR established

      AV   = PR * AVOIGT(NLOC,LRIP) + QR * AVOIGT(NLOC,LRIP+1)

      STARKVOIGT = VOIGTH(AV,XI)

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
      FUNCTION TAUCMU (XMU,L,R,OPA,ND)
C***********************************************************************
C***  CALCULATION OF CONTINUUM OPTICAL DEPTH AT RADIUS R(L)
C***  AS FUNCTION OF MU (ANGLE COSINE)
C***  CORE RAYS: CONTINUUM OPTICAL DEPTH TO INNER BOUNDARY
C***********************************************************************
      DIMENSION R(ND),OPA(ND)

      RL=R(L)
      XSIN2=1.-XMU*XMU
      XSIN=SQRT(XSIN2)
      RMIN=RL*XSIN
      RMIN2=RMIN*RMIN
      TCMU=0.0

C***  INTEGRATION OF TAUC(MU) INSIDE RL (I.E. FOR RADII .LT. RL)
      IF (XMU .GT. 0.0) THEN
         RLMU=RL
         OPALM=OPA(L)
         SXOLD=0.0
         LMAX=ISRCHFLT(ND,R(1),1,RMIN)-1
C***     INTEGRATION FROM RL TO R(LMAX)
         DO 11 LMU=L+1,LMAX
         RLMU1=R(LMU)
         OPALM1=OPA(LMU)
         SX=RL*XMU-SQRT(RLMU1*RLMU1-RMIN2)
         DS=SX-SXOLD
         TCMU=TCMU+(OPALM+OPALM1)/2.*DS
         RLMU=RLMU1
         OPALM=OPALM1
   11    SXOLD=SX
C***     STEP FROM R(LMAX) TO RMIN
         IF (LMAX .GE. ND) THEN

C***        CORE RAY
            TAUCMU=TCMU
            RETURN

         ELSE
            OPARMIN=OPALM+(OPA(LMAX+1)-OPALM)/(RLMU-R(LMAX+1))*
     *                    (RLMU-RMIN)
            SX=RL*XMU
            DS=SX-SXOLD
            TCMU=TCMU+(OPALM+OPARMIN)/2.*DS
         ENDIF
C***     FACTOR 2. FOR SYMMETRICAL RUN OF MU-RAY INSIDE RL
         TCMU=2.*TCMU
      ENDIF

C***  INTEGRATION OF TAUC(MU) OUTSIDE RL (I.E. FOR RADII .GT. RL)
      RLMU=RL
      OPALM=OPA(L)
      SXOLD=0.0
      DO 12 LMU=L-1,1,-1
      RLMU1=R(LMU)
      OPALM1=OPA(LMU)
      SX=RL*XMU+SQRT(RLMU1*RLMU1-RMIN2)
      DS=SX-SXOLD
      TCMU=TCMU+(OPALM+OPALM1)/2.*DS
      RLMU=RLMU1
      OPALM=OPALM1
   12 SXOLD=SX

      TAUCMU=TCMU

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
      SUBROUTINE TRADWL (KANAL,DELW,ADELW,VELO,ENTOT,ND,LINE,MODHEAD,
     $                   JOBNUM,LEVEL,XLAM, XPLOT, TAUROSS, RADIUS,
     $                   TRANSDWLLINE)
C***********************************************************************
C***  DIRECT TRANSFER OF PLOT DATA: DEPTH OF LINE FORMATION ("DWL-PLOT")
C***  Formal Cards Syntax: TRANS DWL <VELO>,<VELOLOG>,<TAUROSS>,<RSTAR>,
C***  <RSTARLOG>. If empty, default = <DENSITY> 
C***********************************************************************
 
      DIMENSION DELW(ND),ADELW(ND),VELO(ND),ENTOT(ND), XPLOT(ND), TAUROSS(ND), RADIUS(ND)
      CHARACTER*60 HEADER, MODHEAD, XTEXT, YTEXT, XAXIS, TRANSDWLLINE*(*), ACTPAR
      INTEGER NPAR
      CHARACTER*10 LEVEL

C**********************************************************************
C***  DEFINE HEADER LINE                                            ***
C**********************************************************************

      CALL SARGC(TRANSDWLLINE, NPAR)
      IF (NPAR .EQ. 2) THEN
        ACTPAR = 'DENSITY'
      ELSE 
        CALL SARGV (TRANSDWLLINE, 3, ACTPAR)
      ENDIF      
      HEADER         = ' '
      WRITE (HEADER (:8), '(4HLINE,I4)') LINE
      HEADER (11:12) = LEVEL (1:2)
      LAM=NINT(XLAM)
      WRITE (HEADER(13:17),'(I5)') LAM
      HEADER (20:25) = 'MODEL'
C***  DATE AND TIME
      HEADER (27:44) = MODHEAD (15:32)
      HEADER (47:49) = 'J >'
      WRITE (HEADER(51:54),'(I4)') JOBNUM
      
C* Y-AXIS SCALING AND PARAMETERS 
C***  SCALING OF Y-AXIS
      YMIN=.0
      YMAX=.0
      DO 4 L=1,ND
      IF (ADELW(L) .LT. YMIN) YMIN=ADELW(L)
    4 IF (ADELW(L) .GT. YMAX) YMAX=ADELW(L)
      YMIN=AINT(YMIN-0.9999999999)
      YMAX=AINT(YMAX+0.9999999999)
      YSCALE=15./(YMAX-YMIN)
      YTICK=YMAX/4.
      YABST=YTICK*2.0
      YTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'NORMALIZED  XI - HILLIER'
C* Plot over radial velocity     
      IF (ACTPAR .EQ. 'VELO') THEN
        XMIN = 0.
        XMAX = VELO(1)
        XSCALE = 0.
        XTICK = 100.
        XABST = 300.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Radial Velocity (km/s)'        
        DO L=1, ND
          XPLOT(L) = VELO(L)
        ENDDO       
C* Plot over log(velocity)        
      ELSE IF (ACTPAR .EQ. 'VELOLOG') THEN
        XMIN = 0.
        XMAX = ALOG10(VELO(1))
        XSCALE = 0.
        XTICK = 0.25
        XABST = 0.5
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Logarithmic Radial Velocity (Log(v/km/s))'          
        DO L=1, ND
          XPLOT(L) = ALOG10(VELO(L))
        ENDDO
C** OD = plot over optical depth         
      ELSE IF (ACTPAR .EQ. 'TAUROSS') THEN
        XMIN = 0.
        XMAX = TAUROSS(1)
        XSCALE = 0.
        XTICK = 1.
        XABST = 3.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'optical depth'          
        DO L=1, ND
          XPLOT(L) = TAUROSS(L)
        ENDDO
C* Plot over radial distance        
      ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
        XMIN = 0.
        XMAX = RADIUS(1)
        XSCALE = 0.
        XTICK = 100.
        XABST = 300.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Radial Distance (R*)'        
        DO L=1, ND
          XPLOT(L) = RADIUS(L)
        ENDDO   
C* Plot over log(radius)
      ELSE IF (ACTPAR .EQ. 'RSTARLOG') THEN
        XMIN = 0.
        XMAX = ALOG10(RADIUS(1))
        XSCALE = 0.
        XTICK = 0.25
        XABST = 0.5
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'logarithmic Radial Distance (Log(R/R*))'        
        DO L=1, ND
          XPLOT(L) = ALOG10(RADIUS(L))
        ENDDO        
C***  Standard: Plot over log of number density        
      ELSE IF (ACTPAR .EQ. 'DENSITY') THEN
        XMIN = 7.
        XMAX = 16.
        XSCALE = 2.5
        XTICK = 1.0
        XABST = 3.0
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'LOG OF NUMBER DENSITY / (CM**-3)'    
        DO L=1, ND
          XPLOT(L) = ALOG10(ENTOT(L))
        ENDDO
      ELSE
        WRITE (0,'(A)') 'DWL x-axis not known'
        WRITE (0,'(A)') 'Syntax: TRANS DWL <VELO>/<VELOLOG>/<TAUROSS>/<RSTAR>/<RSTARLOG>. <> = <DENSITY>' 
        STOP 'ERROR IN SUBROUTINE TRADWL'
      ENDIF
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,YTEXT
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,XPLOT,ADELW,ND,5)


C*** weiterer Plot
ccc  zur Zeit disabled wegen ungueltiger Ausgabe (suffix U in KASDEF ARC)
      return

      XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Depth Point Number'
C***  Plot over depth index
      DO L=1, ND
         XPLOT(L) = FLOAT(L)
      ENDDO 
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,YTEXT
     $ ,0., 0., 80., 5., 10., 0.
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,XPLOT,ADELW,ND,5)
 
      WRITE (KANAL,*) 'PLOT   : DWL-Plot'
      WRITE (KANAL,*) 'KASDEF INBOX'
      WRITE (KANAL,*) 'KASDEF NOBOX'
      WRITE (KANAL,*) 'KASDEF LUN XMIN YMAX 0. D0.1 0.3 ',
     >      'MARKED: 0.1, 3/4*MAX, MAX'
      WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.7 0.7 0.7'
      WRITE (KANAL,*) 'KASDEF COLOR=2'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0.  2.0U 0. 180. FILLED'
      WRITE (KANAL,*) 'KASDEF COLOR=1'
C***  FIND MAXIMUM
      LMAX = 1
      DO L=2, ND
        IF (ADELW(L) .GT. ADELW(LMAX)) LMAX = L
      ENDDO
      R = 2. + FLOAT(LMAX) / 5.
      ADMAX = ADELW(LMAX)
      ADMAX2 = ADMAX * 3. / 4.
      WRITE (KANAL,*) 'KASDEF PEN = 8'
      WRITE (KANAL,'(A,F4.1,A)') 
     >      'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
      WRITE (KANAL,*) 'KASDEF PEN = 1'
C***  FIND 0.1
      DO L=2, ND
        IF ((ADELW(L-1)-0.1)*(ADELW(L)-0.1) .LE. 0.) THEN
          R = 2. + FLOAT(L) / 5.
          WRITE (KANAL,*) 'KASDEF PEN = 8'
          WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.7 0.7 0.7'
          WRITE (KANAL,*) 'KASDEF COLOR=2'
          WRITE (KANAL,'(A,F4.1,A)') 
     >          'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
          WRITE (KANAL,*) 'KASDEF COLOR=1'
          WRITE (KANAL,*) 'KASDEF PEN = 1'
        ENDIF
        IF ((ADELW(L-1)-ADMAX2)*(ADELW(L)-ADMAX2) .LE. 0.) THEN
          R = 2. + FLOAT(L) / 5.
          WRITE (KANAL,*) 'KASDEF PEN = 8'
          WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.5 0.5 0.5'
          WRITE (KANAL,*) 'KASDEF COLOR=2'
          WRITE (KANAL,'(A,F4.1,A)') 
     >          'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
          WRITE (KANAL,*) 'KASDEF PEN = 1'
        ENDIF
      ENDDO      

      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 3.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 23. 0. M0. U-0.1 0.2 65'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 4.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 24. 0. M0. U-0.1 0.2 60'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 5.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 25. 0. M0. U-0.1 0.2 55'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 6.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 26. 0. M0. U-0.1 0.2 50'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 7.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 27. 0. M0. U-0.1 0.2 45'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 8.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 28. 0. M0. U-0.1 0.2 40'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 9.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 29. 0. M0. U-0.1 0.2 35'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 10.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 30. 0. M0. U-0.1 0.2 30'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 11.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 31. 0. M0. U-0.1 0.2 25'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 12.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 32. 0. M0. U-0.1 0.2 20'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 13.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 33. 0. M0. U-0.1 0.2 15'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 14.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 34. 0. M0. U-0.1 0.2 10'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 15.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 35. 0. M0. U-0.1 0.2 5'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 16.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 36. 0. M0. U-0.1 0.2 1'

      XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Depth Points'
c      YTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
c     >        'NORMALIZED  XI - HILLIER'
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,' '
     $ ,0., 0., 40., 25., 25., 0. 
     $ ,0., 0., 30., 20., 20., 0.
     $ ,XDUM, XDUM, 1, 8)

      RETURN
      END
      SUBROUTINE TRANSFORM_RGRID (VEC_NEW, N_NEW, VEC_OLD, N_OLD, 
     >                            R_NEW, R_OLD)
C********************************************************************
C***  Transformation of vector VEC_OLD (N_OLD) over radius-grid R_OLD
C***               into vector VEC_NEW (N_NEW) over radius-grid R_NEW
C***  If the new grid exceeds the old one: extrapolation as constant
C*********************************************************************

      DIMENSION VEC_NEW(N_NEW), R_NEW(N_NEW)
      DIMENSION VEC_OLD(N_OLD), R_OLD(N_OLD)

      DATA IWARN / 0 /

      IF (IWARN .EQ. 0) THEN 
         IF ( R_NEW(1) .GT. R_OLD(1) ) THEN 
         WRITE (0,'(A,/, A, F8.2,A,F8.2)') 
     >        '*** WARNING issued by subr. TRANSFORM_RGRID: ',  
     >        ' Extrapolation needed from ', R_NEW(1), ' TO ', R_OLD(1)
            IWARN = 1
         ENDIF
      ENDIF

      DO L=1, N_NEW
         IF (R_NEW(L) .GT. R_OLD(1)) THEN
            VEC_NEW(L) = VEC_OLD(1)
         ELSE
            CALL LIPO (VEC_NEW(L), R_NEW(L), VEC_OLD, R_OLD, N_OLD)
         ENDIF
      ENDDO

      RETURN
      END
      SUBROUTINE TRAPLO (KANAL,PROFILE,DLAM,NFOBS,LINE,
     >                   FREQIN,MODHEAD,JOBNUM,
     $                   DISP,N,NBLINE,IDENT,OSMIN,
     >                   XLAM,XLAMLAP,INDLAP,INDLOW,INDNUP,LEVEL,
     >                   ELEVEL,WEIGHT,EINST,NDIM,ABSWAV,BCONT,FUNIT, 
     >                   NMOD, LPSTA, LPEND, NPHI, XUNIT, 
     >                   BCALIBRATED, BAIRWAVELENGTH, RSTAR, 
     >                   KODATIND, MAXATOM, NOM)

C**********************************************************************
C***               DIRECT TRANSFER OF PLOT DATA                     ***
C**********************************************************************
 
      DIMENSION DLAM(NFOBS),PROFILE(NFOBS),DLAMPLO(NFOBS)
      DIMENSION XLAMLAP(NBLINE),INDLAP(NBLINE)
      DIMENSION INDLOW(N),INDNUP(N),ELEVEL(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION JOBNUM(NMOD)
      DIMENSION KODATIND(MAXATOM), NOM(NDIM)

      CHARACTER(*) :: FREQIN    !Name des RANGE-Bereiches
      CHARACTER(60) :: XCHAR 
      CHARACTER(65) :: YCHAR
      CHARACTER(100) :: HEADER
      CHARACTER(*) :: MODHEAD
      CHARACTER*10 LEVEL(N), FUNIT*(*)
      CHARACTER*1 BACKSLASH
      CHARACTER*8 IDSTR
      CHARACTER LAMBDASUBSCRIPT*7
      CHARACTER*(*) XUNIT
      LOGICAL AUTOMAX, IDENT, ABSWAV, BCONT, BCALIBRATED
      LOGICAL BAIRWAVELENGTH, BIDENT_ON
 
      REAL, PARAMETER :: PI = 3.14159   !Pi
      REAL, PARAMETER :: C4 = 1.499E-16 !C4 = 1 / (SIGMA-CLASSIC * 8 PI)     (IN ANGSTROEM**2)

C***  CLIGHT = SPEED OF LIGHT in Angstroem per second
      DATA CLIGHT2 / 2.99792458E18 /

C***  CALIBRATION CONSTANT TO USE FOR THE CONVERSION OF F-NUE TO MV
      DATA CONMV / 48.64 /

C***  PARSEC IN CM
      DATA PARSEC / 3.08561E18 /


      BACKSLASH = CHAR(92)

C**********************************************************************
C***  DEFINE HEADER LINE                                            ***
C**********************************************************************
 
      HEADER         = ' '
      HEADER(1:20) = FREQIN(1:20)

      IF (BCONT .AND..NOT. BCALIBRATED) HEADER(23:31) = 'Continuum'
      IF (BCONT .AND. BCALIBRATED     ) HEADER(23:31) = 'Absolut'

      HEADER (36:36) = 'M'
C***  DATE AND TIME
      WRITE (HEADER(37:),'(7A2)') MODHEAD(13:14),
     >       MODHEAD(16:17),MODHEAD(19:20),'  ',MODHEAD(25:26),
     >       MODHEAD(28:29),MODHEAD(31:32)
      HEADER (53:55) = 'J >'
      WRITE (HEADER (56:61),'(I6)') JOBNUM(1)

C***  Subscript to Lambda: Air or Vacuum
      IF (BAIRWAVELENGTH) THEN
         LAMBDASUBSCRIPT = '&Tair&M'
      ELSE
         LAMBDASUBSCRIPT = '&Tvac&M'
      ENDIF
 
C**********************************************************************
C***                SCALING OF Y-AXIS                               ***
C**********************************************************************

      IF (BCONT .OR. BCALIBRATED) THEN
                IF    (FUNIT .NE. 'MV'        
     >           .AND. FUNIT .NE. 'FLAM10PC'
     >           .AND. FUNIT .NE. 'LOGFLAM'
     >           .AND. FUNIT .NE. 'LOGFNUE'       ) THEN        
                    FUNIT = 'LOGFNUE'
                    WRITE(0,'(A)') '*** WARNING: ' //
     >                 'Flux units invalid: taking LOGFNUE instead'
                ENDIF

C***     Convert flux units if appropriate
         IF (FUNIT .EQ. 'LOGFLAM') THEN
            DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) = ALOG10 (PROFILE(K) * CLIGHT2 / 
     >                     (DLAM(K)+XLAM)**2)
               ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'LOGFNUE') THEN
            DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) =
     >               ALOG10 (PROFILE(K))
               ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'MV') THEN
             DISMO = -2.5 * ALOG10 (RSTAR*RSTAR /
     >                         100. / (PARSEC * PARSEC))
             DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) =
     >               -2.5 * ALOG10 (PI*PROFILE(K)) + DISMO - CONMV
 
              ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'FLAM10PC') THEN
               FAC = PI * RSTAR*RSTAR / (PARSEC*PARSEC) / 100.
               DO K=1, NFOBS
                  PROFILE(K) = PROFILE(K) * CLIGHT2 / 
     >                         (DLAM(K)+XLAM)**2 * FAC 
               ENDDO 
          ENDIF
      ENDIF

C***  Invoke automatic scaling of Y-axis by WRplot
C***  Note: the formar YMAX option is now disabled (wrh 9-Jun-2014)
      YMIN=0.0
      YMAX=0.0

      IF (.NOT.((BCONT .OR. BCALIBRATED) .AND. FUNIT .EQ. 'MV')) 
     >   GOTO 1
C***  Manual scaling of Y axis for MV (negative direction!)
         YMIN = PROFILE(1)    
         YMAX = PROFILE(1)    
         DO K=2,NFOBS
            IF (PROFILE(K) .GT. YMIN) YMIN=PROFILE(K) 
            IF (PROFILE(K) .LT. YMAX) YMAX=PROFILE(K) 
         ENDDO
         YDIFF = YMAX - YMIN
         YMIN = YMIN - 0.05 * YDIFF 
         YMAX = YMAX + 0.05 * YDIFF 

         IF (YMAX .LT. -100.) THEN
           YMAX=-100
         ELSEIF (YMIN .GT. 100.) THEN
           YMIN=100.
         ENDIF

         YTICK = FLOAT (NINT(YDIFF*20.))/400.
         YABST = 10.*YTICK

    1 YSCALE=0.

C**********************************************************************
C***                SCALING OF X-AXIS                               ***
C**********************************************************************

C***  X-Units in micrometer (default: Angstroem)
      IF (XUNIT(:1) .EQ. 'M') THEN
         XUNITFAC = 1.E-4
      ELSE
         XUNITFAC = 1.
      ENDIF

      XLAMPLO = .0
      IF (ABSWAV) XLAMPLO = XLAM 

      DO 110 I=1,NFOBS
        DLAMPLO(I) = ( DLAM(I) + XLAMPLO ) * XUNITFAC
  110 CONTINUE

C**********************************************************************
C***  BREITE DES PLOTS IN CM                                        ***
C**********************************************************************

      BREITE= 20.
      IF (DISP .EQ. .0) THEN

C**********************************************************************
C***  DISPERSION NICHT SPEZIFIZIERT: PLOT UEBER DAS GERECHNETE      ***
C***                                 INTERVALL                      ***
C**********************************************************************

        XMIN = DLAMPLO(1)
        XMAX = DLAMPLO(NFOBS)

      ELSE

C**********************************************************************
C***  DISPERSION GEMAESS DER INPUT OPTION                           ***
C**********************************************************************

        XMAX= 0.5*BREITE*DISP + XLAMPLO
        XMIN=-0.5*BREITE*DISP + XLAMPLO

      ENDIF

      XSCALE=0.
 
C**********************************************************************
C***  BESCHRIFTUNGEN UND TICK-MARKS                                 ***
C**********************************************************************
      IF (XMAX-XMIN .LE. 1.) THEN
         XABST=.1
         XTICK=.02
      ELSE IF (XMAX-XMIN .LE. 3.) THEN
         XABST=.5
         XTICK=.1
      ELSE IF (XMAX-XMIN .LE. 10.) THEN
         XABST=1.
         XTICK=0.2
      ELSE IF (XMAX-XMIN .LE. 30.) THEN
         XABST=5.
         XTICK=1.
      ELSE IF (XMAX-XMIN .LE. 100) THEN
         XABST=10.
         XTICK=2.
      ELSE IF (XMAX-XMIN .LE. 300.) THEN
         XABST=50.
         XTICK=10.
      ELSE IF (XMAX-XMIN .LE. 1000) THEN
         XABST=100.
         XTICK=20.
      ELSE IF (XMAX-XMIN .LE. 3000.) THEN
         XABST=500.
         XTICK=100.
      ELSE IF (XMAX-XMIN .LE. 10000) THEN
         XABST=1000.
         XTICK=200.
      ELSE
         XABST= 5000.
         XTICK= 1000.
      ENDIF
 
C**********************************************************************
C     XMIN=XTICK*(IFIX(DLAM(  1  )/XTICK)-1.)                       ***
C     XMAX=XTICK*(IFIX(DLAM(NFOBS)/XTICK)+1.)                       ***
C**********************************************************************
 
C**********************************************************************
C***  WRITE MODEL HEADERS VERTICALLY
C**********************************************************************
      WRITE (KANAL,'(A)') 'PLOT: ' // HEADER
      XPOS = 1.5

C      DO IMOD=1, NMOD
        WRITE (KANAL, '(A,F3.1,2A)') 
     >    '\LUNA XMAX YMAX ', XPOS, ' 0. 0.3 -90. '
C        WRITE (KANAL, '(A)') 
C     >    '\> &E' // MODHEAD(IMOD)( :IDX(MODHEAD(IMOD)))
        WRITE (KANAL, '(A)') 
     >    '\> &E' // MODHEAD(:IDX(MODHEAD))
        XPOS = XPOS - 0.50
C      ENDDO

C**********************************************************************
C     IF (IDENT): LINE IDENTIFICATION 
C**********************************************************************

      IF (IDENT) THEN
          WRITE (KANAL, '(A,G8.2)') 
     >          '\LAB 17.5 15.3 0.2 IDENT: &Rf&N >= ',OSMIN
          WRITE (KANAL,*) '\ID_INBOX'
          WRITE (KANAL,*) '\IDSIZE=0.2'
          DO 105 NBL=NBLINE,1,-1
          IF (BAIRWAVELENGTH) THEN
             XLAMREF = XLAMLAP(NBL)
             XLAM2   = XLAMREF * XLAMREF
             XLAMREF = XLAMREF - XLAMREF*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
          ELSE
             XLAMREF = XLAMLAP(NBL)
          ENDIF
          DXLAM = XLAMREF - XLAM + XLAMPLO
          DXLAM = DXLAM * XUNITFAC
          IND=INDLAP(NBL)
          NUP=INDNUP(IND)
          LOW=INDLOW(IND)
          XLAMI=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
          F = C4 * XLAMI * XLAMI * 
     >        EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***      IDENT shall not be commented for H, He or if f-value > OSMIN
          BIDENT_ON = (F .GT. OSMIN) .OR.( KODATIND(NOM(LOW)) .LE. 2)
          IF (BIDENT_ON) THEN
              IDSTR = ' \IDENT '
          ELSE
              IDSTR = '*\IDENT '
          ENDIF
          WRITE (KANAL,101) IDSTR, DXLAM, LEVEL(NUP), LEVEL(LOW)
  101     FORMAT (A8, G14.7,' &E',A10, ' - ', A10,'&N')

  105     CONTINUE
      ENDIF

      IF (XUNIT(:1) .EQ. 'M') THEN
       IF (ABSWAV) THEN
         XCHAR  = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#l#'
     >     // LAMBDASUBSCRIPT // '  / #m#m'
       ELSE
         XCHAR = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#D l#'
     >     // LAMBDASUBSCRIPT // '  / #m#m'
       ENDIF
      ELSE
       IF (ABSWAV) THEN
         XCHAR  = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#l#' 
     >     // LAMBDASUBSCRIPT // '  / \A'
       ELSE
         XCHAR = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#D l#'
     >     // LAMBDASUBSCRIPT // '  / \A'
       ENDIF
      ENDIF

C***  Y-AXIS descriptor
      IF (BCONT .OR. BCALIBRATED) THEN
        IF (FUNIT .EQ. 'MV  ') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 'M&T#l#&M [mag]'
        ELSE IF (FUNIT .EQ. 'FLAM10PC') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 
     >     'f&T#l#&M / (erg cm&H-2&M s&H-1&M ' //
     >     BACKSLASH // 'A&H-1&M) at 10 pc'
        ELSE IF (FUNIT .EQ. 'LOGFLAM') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 
     >     'log F&T#l#&M / (erg cm&H-2&M s&H-1&M ' //
     >     BACKSLASH // 'A&H-1&M)'
        ELSE IF (FUNIT .EQ. 'LOGFNUE') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH //
     >     'log F&T#n#&M / (erg cm&H-2&M s&H-1&M Hz&H-1&M)'
        ELSE
          YCHAR = 'units undefined -- internal error'
        ENDIF
      ELSE
C***  else: normalized profile
         YCHAR= BACKSLASH // 'CENTER' // BACKSLASH // 'rel. Flux'
         WRITE (KANAL, '(A)') BACKSLASH // 'LINUN XMIN 1 XMAX 1' 
      ENDIF

      CALL PLOTANFS (KANAL,HEADER,HEADER
     $ ,XCHAR,YCHAR
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,DLAMPLO,PROFILE,NFOBS,'COLOR=2')
      RETURN
      END
      SUBROUTINE TRBK

      WRITE (*,*) 'TRACE BACK FACILITY (TRBK) IS NOT AVAILABLE AT',
     >            ' DEC/OSF1 AT THE MOMENT'

      CALL ABORT()
      RETURN
      END
      SUBROUTINE VADD (A,B,N)
C***********************************************************************
C***  VECTOR ADDITION  A = A + B
C***********************************************************************
      DIMENSION A(N),B(N)

      DO 1 I=1,N
    1 A(I)=A(I)+B(I)

      RETURN
      END
	FUNCTION VCSE1F(X)
C
C  E1 function calculator for VCS approximation. It's rough, but 
C  arranged to be fast. X must be >=0.
C
C  From Atlas9 & Kurucz code (mario)
C  Approximation der vidal-cooper-smith daten
C

      VCSE1F=0.0
      IF(X.LE.0.0) RETURN
      IF(X.LE.0.01) THEN
        VCSE1F=-LOG(X)-0.577215+X
      ELSE IF(X.LE.1.0) THEN
        VCSE1F=-LOG(X)-0.57721566+X*(0.99999193+X*(-0.24991055+
     +                            X*(0.05519968+X*(-0.00976004+
     +                            X*0.00107857))))
      ELSE IF(X.LE.30.) THEN
        VCSE1F=(X*(X+2.334733)+0.25062)/(X*(X+3.330657)+
     +         1.681534)/X*EXP(-X)
      END IF
C

      RETURN
      END
      SUBROUTINE VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, DD_VDOP, VDOP, 
     >                        VMIC_MODEL, VELO, T, 
     >                        ND, NDDIM, NATOM, MAXATOM, DD_VDOPDU, 
     >                        VMICFRAC_DEFAULT, ATMASS, 
     >                        XMAX, XMAXMIN, 
     >                        SYMBOL, VDOPFE, 
     >                        BMICROTURB, BIRONLINES,
     >                        DD_VMIC, TAUROSS, RADIUS, EL_MIN, 
     >                        DD_VMICDU, bDDFECONVOL, IVDOPSTATUS)

C********************************************************************************
C***  Called from formal, 
C***  This routine prepares the depth-dependent DD_VDOP and further,
C***  related, depth-dependent arrays ("DD_..."). 
C***
C***  DD_VDOP(L,NA) returns the Dopplerbroadening -velocity 
C***  at depth point L of atom NA. 
C***  If no depth-dependence is specified, neither for VDOP not VMIC,
C***  then DD_VDOP(L,NA) = VDOP = constant.
C***
C***  Otherwise, 
C***  DD_VDOP(L,NA) = MAX(VDOP, SQRT (Vmic(L)^2 + Vtherm(L,NA)^2))
C***
C***  The squared thermal velocity is given by 2 * k_B * T[L] / m[NA]
C***  
C***  The stratification of Vmic(L) can be specified by the VMIC 
C***  line in FORMAL_CARDS.
C***  Current options: VDOP X.X [VELOFRAC|OUTERVMIC|VMICCONST [X.X]].
C***  In the first two versions, Vmic[L] = f*VELO[L], where the user either
C***  gives f explicity (VELOFRAC) or implicitly 
C***  (OUTERVMIC - f = OUTERVMIC/VELO(1))
C***  In the last version, VMIC is assumed to be constant.
C***  If no value is  stated: 

C***  VDOP is the MINIMUM value of DD_VDOP (relevant for the frequency
C***  spacing to resolve the narrowest lines).  
C***
C***  In case of a SECONDMODEL, this subroutine is caaled a second time, 
C***   with the corresponding parts of the arrays addressed by the CALL.
C***
C********************************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND, NDDIM, NATOM, MAXATOM
      INTEGER NA_VDOPMIN
      INTEGER  :: L, NA, NPAR, I, NA_MIN, ISRCHEQ, IVDOPSTATUS
      REAL, DIMENSION(NATOM), INTENT(IN) :: ATMASS

      REAL, DIMENSION(ND), INTENT(IN) :: VELO, T, TAUROSS, RADIUS,
     >                                   VMIC_MODEL
      REAL, DIMENSION(ND) :: DD_VMIC, DD_VMICDU

      REAL, INTENT(IN) :: VMICFRAC_DEFAULT, XMAXMIN 
      REAL :: VDOP, VMIC_FRAC, VDOP_FRAC, DD_VMIC_SQRD, VDOPMIN, MASS 
      REAL :: THERMVELO_SQRD, XMAX, Q, DUMMY
      REAL :: VMIC_MIN, VMIC_MAX, VDOP_MAX, DIST_IN, DIST_OUT, DX, X
      CHARACTER*2 SYMBOL(NATOM)
      REAL, INTENT(IN) :: VDOPFE
      REAL, DIMENSION(NDDIM, MAXATOM):: 
     >       DD_VDOP, DD_VDOPDU
      LOGICAL :: BDD_VDOP, BMICROTURB, BIRONLINES, bDDFECONVOL
      CHARACTER ACTPAR*20, DD_VDOP_LINE*(*), EL_MIN*2
      INTEGER, EXTERNAL :: IDX
C***  WPIINV = 1. / SQRT(PI)
      REAL, PARAMETER :: WPIINV = 0.564189583549 
      REAL, PARAMETER :: PI = 3.141592654 
      REAL VDOPNOIRONMAX
C*** atomic mass unit in g
      REAL, PARAMETER :: AMU = 1.6605387E-24 
C*** IMPORTANT! K_B is given here in units km^2 g s^-2 K^1 to obtain 
C*** a velocity in [km/s]!
      REAL, PARAMETER :: BOLTZK = 1.3807E-26
      

C***  VMIC VERSION => VDOP^2 = VMIC^2 + VTHERM^2 ************
C***  Note: DD_VDOP(L,NA) depends on element (NA)!
      IF (BMICROTURB) THEN
        WRITE(0,'(A)') "VDOP is depth-dependent"
        BDD_VDOP = .TRUE.
        WRITE(0,'(A)') "VDOP^2 = vmic^2 + vtherm^2"
        CALL SARGC (DD_VDOP_LINE, NPAR)
        IF ((NPAR .LT. 2) .OR. (NPAR .GT. 8)) GOTO 102
        CALL SARGV (DD_VDOP_LINE, 2, ACTPAR)
        IF (ACTPAR == 'MODEL') THEN
C***      use VMIC as specified in the MODEL file
          IF (MINVAL(VMIC_MODEL) < 0.) THEN
C***        FATAL ERROR if no (useful) VMIC is stored in the MODEL
            WRITE (0, '(A)') '*** ERROR: No VMIC stored in MODEL file'
            WRITE (0, '(A)') '*** VMIC MODEL is invalid in this case!'
            STOP '*** Fatal error in vdop_struct'            
          ENDIF
          DO L=1, ND
            DD_VMIC(L) = VMIC_MODEL(L)
          ENDDO
        ELSE
C***      VMIC was specified via FORMAL_CARDS
          READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_MIN
C***      VMIC(L) is never smaller than specified by user
          WRITE(0,'(A,F10.5)') "Inner microturbulence:", VMIC_MIN
          IF (NPAR .EQ. 2) THEN
            WRITE(0,'(A,F10.2)') 'VMIC is constant:', VMIC_MIN
            DO L=1, ND
              DD_VMIC(L) = VMIC_MIN
            ENDDO
          ELSE 
            CALL SARGV (DD_VDOP_LINE, 3, ACTPAR)
C***        In this branch: vmic(L) = VMIC_FRAC * VELO(L)
            IF (ACTPAR .EQ. 'VELOFRAC') THEN
                WRITE(0,'(A)') "VMIC is fraction of wind velocity"
                IF (NPAR .EQ. 3) THEN
                    WRITE(0,'(A,F10.5)') 
     >               "fraction set to default value:", VMICFRAC_DEFAULT
                    VMIC_FRAC = VMICFRAC_DEFAULT
                ELSE IF (NPAR .EQ. 4) THEN
                    CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                    READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_FRAC
                    WRITE(0,'(A,F10.5)') 
     >               "fraction set by user to", VMIC_FRAC
                ELSE
                    GOTO 102
                ENDIF
                DO L=1, ND
                    DD_VMIC(L) = AMAX1(VMIC_MIN, VMIC_FRAC*VELO(L))
                ENDDO
C***        In this branch, the outer VMIC is specified by user
            ELSE IF (ACTPAR .EQ. 'MAX') THEN
                CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_MAX
                WRITE(0,'(A,F6.1)') "Outer microturbulence:", VMIC_MAX
C***            If no further arguments are given, 
C***              VMIC is calculated as VMIC(L) = VMIC_FRAC * VELO(L)
                IF (NPAR .EQ. 4) THEN
                    WRITE(0,'(A)') 
     >              "microturbulence = fraction of wind velocity"
                    VMIC_FRAC = VMIC_MAX/VELO(1)
                    WRITE(0,'(A,F10.5)') 
     >               "fraction implied from maximum microturbulence:", 
     >               VMIC_FRAC
                    DO L=1, ND
                        DD_VMIC(L) = AMAX1(VMIC_MIN, VMIC_FRAC*VELO(L))
                    ENDDO
C***            Otherwise, inner and outer interpolation points 
C***            and method specified by user:
                ELSE IF (NPAR .EQ. 8) THEN
                    CALL SARGV (DD_VDOP_LINE, 5, ACTPAR)
C***                1) Interpolation on velocity
                    IF (ACTPAR .EQ. 'VELO1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'VELO2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF 
                        IF ((DIST_IN .GE. VELO(1)) .OR. (DIST_OUT .LE. VELO(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner wind velocity: ", VELO(ND)
                            WRITE(0,'(A,F10.5)') "Outer wind velocity: ", VELO(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
C***                    interpolation boundaries may not exceed model boundaries
                        DIST_IN = AMAX1(DIST_IN, VELO(ND))
                        DIST_OUT = AMIN1(DIST_OUT, VELO(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between wind velocity", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
C***                    Interpolation takes place here
                        DO L=1, ND
                            IF (VELO(L) .LE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (VELO(L) .GE. DIST_OUT) THEN
                                DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (VELO(L) - DIST_IN ) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
C***                2) Interpolation on Rosseland tau (analog to last block)
                    ELSE IF (ACTPAR .EQ. 'TAU1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'TAU2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN < DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF 
                        IF ((DIST_IN .LE. TAUROSS(1)) .OR. (DIST_OUT .GE. TAUROSS(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner tau: ", TAUROSS(ND)
                            WRITE(0,'(A,F10.5)') "Outer tau: ", TAUROSS(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
                        DIST_IN = AMIN1(DIST_IN, TAUROSS(ND))
                        DIST_OUT = AMAX1(DIST_OUT, TAUROSS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between Rosseland tau", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_IN - DIST_OUT
                        DO L=1, ND
                            IF (TAUROSS(L) .GE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (TAUROSS(L) .LE. DIST_OUT) THEN
                                DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (DIST_IN - TAUROSS(L)) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = 
     >                            Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
C***                3) Interpolation on Radius (analog to last block)
                    ELSE IF (ACTPAR .EQ. 'R1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'R2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF                    
                        IF ((DIST_IN .GE. RADIUS(1)) .OR. (DIST_OUT .LE. RADIUS(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner radius: ", RADIUS(ND)
                            WRITE(0,'(A,F10.5)') "Outer radis: ", RADIUS(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
                        DIST_IN = AMAX1(DIST_IN, RADIUS(ND))
                        DIST_OUT = AMIN1(DIST_OUT, RADIUS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between radii", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            IF (RADIUS(L) .LE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (RADIUS(L) .GE. DIST_OUT) THEN
                                 DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (RADIUS(L) - DIST_IN ) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = 
     >                            Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
                    ELSE
                        GOTO 102
                    ENDIF
C***                End of interpolation branches 
                ELSE
                    GOTO 102
                ENDIF
C***            End of VMICMAX branch
            ELSE
                GOTO 102
            ENDIF
          ENDIF
C***      End of depth-dependent VMIC
        ENDIF
C***    VMIC(L) is now defined for every case
C***    VDOP(L,NA) may now be filled via 
C***      VDOP^2(L,NA) = VMIC^2(L) + VTHERM^2(L,NA)
        DO NA=1, NATOM  
           MASS = ATMASS(NA)* AMU
           DO L=1, ND               
C***        v-thermal squared, in km^2/s^2
            THERMVELO_SQRD = 2 * BOLTZK * T(L) / MASS  
            DD_VMIC_SQRD = DD_VMIC(L) * DD_VMIC(L)
            DD_VDOP(L,NA) = SQRT(THERMVELO_SQRD + DD_VMIC_SQRD)
           ENDDO
        ENDDO
*********** END OF VMIC VERSION ************
      ELSE 
C***    VDOP VERSION => No microturbulence or thermal motion 
C***    Note: DD_VDOP(L,NA) does not depend on the element here 
C***          (except for iron), i.e. VDOP(L,1) = VDOP(L,2) = ... 
C***    Note2: This branch is almost completely analog to the 
C***           VMIC branch - (see detailed comments there)

        CALL SARGC (DD_VDOP_LINE, NPAR)
        IF (NPAR .GT. 8) GOTO 101   
        IF (NPAR .LE. 2) THEN
            WRITE(0,'(A)') "VDOP not depth-dependent"
            DO L = 1,ND
                DO NA=1, NATOM
                    DD_VDOP(L, NA) = VDOP
                ENDDO
            ENDDO
        ELSE 
            BDD_VDOP = .TRUE.
            WRITE(0,'(A)') "VDOP is depth-dependent"
            CALL SARGV (DD_VDOP_LINE, 3, ACTPAR)
            IF (ACTPAR .EQ. 'VELOFRAC') THEN
                WRITE(0,'(A)') "VDOP = fraction of wind velocity"
                IF (NPAR .EQ. 3) THEN
                    WRITE(0,'(A,F10.5)') 
     >                 "fraction set to default value:", VMICFRAC_DEFAULT
                    VDOP_FRAC = VMICFRAC_DEFAULT
                ELSE IF (NPAR .EQ. 4) THEN 
                    CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                    READ (ACTPAR, '(F20.0)', ERR = 101) VDOP_FRAC
                    WRITE(0,'(A,F10.5)') 
     >                  "fraction set by user to", VDOP_FRAC
                ELSE
                    GOTO 101
                ENDIF
                DO L=1, ND
                    DO NA = 1, NATOM
                        DD_VDOP(L, NA) = AMAX1(VDOP, VDOP_FRAC*VELO(L))
                    ENDDO
                ENDDO
            ELSE IF (ACTPAR .EQ. 'MAX') THEN
                CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 101) VDOP_MAX
                WRITE(0,'(A,F6.1)') "Maximum VDOP:", VDOP_MAX
                IF (NPAR .EQ. 4) THEN
                    WRITE(0,'(A)') "VDOP = fraction of wind velocity"
                    VDOP_FRAC = VDOP_MAX/VELO(1)
                    WRITE(0,'(A,F10.5)') 
     >                 "fraction implied from maximum VDOP:", VDOP_FRAC
                    DO L=1, ND
                        DO NA = 1, NATOM
                          DD_VDOP(L,NA) = MAX(VDOP, VDOP_FRAC*VELO(L))
                          IF (SYMBOL(NA) .EQ. 'G ') THEN
                            IF (bDDFECONVOL) THEN
                              DD_VDOP(L,NA) = MAX(VDOPFE, DD_VDOP(L,NA))
                            ELSE 
                              DD_VDOP(L,NA) = VDOPFE
                            ENDIF
                          ENDIF
                        ENDDO
                    ENDDO
                ELSE IF (NPAR .EQ. 8) THEN 
                    CALL SARGV (DD_VDOP_LINE, 5, ACTPAR)
                    IF (ACTPAR .EQ. 'VELO1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'VELO2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        IF ((DIST_IN .GE. VELO(1)) .OR. (DIST_OUT .LE. VELO(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner wind velocity: ", VELO(ND)
                            WRITE(0,'(A,F10.5)') "Outer wind velocity: ", VELO(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
                        DIST_IN = AMAX1(DIST_IN, VELO(ND))
                        DIST_OUT = AMIN1(DIST_OUT, VELO(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "VDOP interpolated between wind velocity", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            DO NA = 1, NATOM
                                IF (VELO(L) .LE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (VELO(L) .GE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (VELO(L) - DIST_IN ) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                    DD_VDOP(L,NA) = 
     >                                Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                       ENDDO
                    ELSE IF (ACTPAR .EQ. 'TAU1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'TAU2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN < DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        IF ((DIST_IN .LE. TAUROSS(1)) .OR. (DIST_OUT .GE. TAUROSS(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner tau: ", TAUROSS(ND)
                            WRITE(0,'(A,F10.5)') "Outer tau: ", TAUROSS(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
                        DIST_IN = AMIN1(DIST_IN, TAUROSS(ND))
                        DIST_OUT = AMAX1(DIST_OUT, TAUROSS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)')
     >                   "VDOP interpolated between Rosseland tau", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_IN - DIST_OUT
                        DO L=1, ND
                            DO NA=1, NATOM
                                IF (TAUROSS(L) .GE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (TAUROSS(L) .LE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (DIST_IN - TAUROSS(L)) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                DD_VDOP(L,NA) = 
     >                              Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                        ENDDO
                    ELSE IF (ACTPAR .EQ. 'R1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'R2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        IF ((DIST_IN .GE. RADIUS(1)) .OR. (DIST_OUT .LE. RADIUS(ND))) THEN
                            WRITE(0,'(A,F10.5)') "Interpolation boundaries not in model boundaries!"
                            WRITE(0,'(A,F10.5)') "Inner radius: ", RADIUS(ND)
                            WRITE(0,'(A,F10.5)') "Outer radius: ", RADIUS(1)
                            STOP "Fatal error in subroutine VDOP_STRUCT"
                        ENDIF
                        DIST_IN = AMAX1(DIST_IN, RADIUS(ND))
                        DIST_OUT = AMIN1(DIST_OUT, RADIUS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "VDOP interpolated between radii", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            DO NA=1, NATOM
                                IF (RADIUS(L) .LE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (RADIUS(L) .GE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (RADIUS(L) - DIST_IN ) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                    DD_VDOP(L,NA) = 
     >                                 Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                        ENDDO
                    ELSE
                        GOTO 101
                    ENDIF
                ELSE
                    GOTO 101
                ENDIF
            ELSE
                GOTO 101
            ENDIF
        ENDIF
C*** VMIC array is filled even when not specified so that convolution works in subr. STARKBROAD
      DO L=1, ND
        DD_VMIC(L) = DD_VDOP(L,1)
      ENDDO
*********** END OF VDOP VERSION ************
      ENDIF
C*** Generic (Symbol 'G') has its own VDOP velocity VDOPFE from the iron file
C*** if NO-IRONLINES is given in FORMAL_CARDS, the Doppler velocity of FE is set 
C*** to the largest value in the array to avoid redundent increase of resolution
      VDOPNOIRONMAX = MAXVAL(DD_VDOP)
      DO NA=1, NATOM 
        IF (SYMBOL(NA) .EQ. 'G') THEN
          DO L=1, ND               
            IF (BIRONLINES .AND. bDDFECONVOL) THEN
               DD_VDOP(L,NA) = MAX(VDOPFE, DD_VDOP(L,NA))
            ELSEIF (BIRONLINES) THEN
               DD_VDOP(L,NA) = VDOPFE
            ELSE
               DD_VDOP(L,NA) = VDOPNOIRONMAX
            ENDIF
          ENDDO
        ENDIF
      ENDDO

C*** Final VDOP: minimum of matrix DD_VDOP(L,NA) (**Including iron!!)
C*** NOTE1: It could be that the final VDOP is therefore different 
C***        than stated by user!
C*** NOTE2: It could be that the final VDOP is *SMALLER* than stated 
C***        by user if VDOPFE < VDOP_USER!
      VDOP = MINVAL(DD_VDOP(:ND,:NATOM))
      WRITE(0,'(A,F6.1)') "Minimum VDOP:", VDOP
      IF (VDOP .LE. 0) THEN
        WRITE (0,'(A)') 'Invalid VDOP=', VDOP
        STOP '*** FATAL ERROR in subr. VDOP_STRUCT'
      ENDIF

C***  EL_MIN = element with narrowest lines (smallest DD_VDOP)
C***  -> Handed to PRIPRO for output 
      VDOPMIN = MINVAL(DD_VDOP(:ND,1))  
      NA_VDOPMIN = 1
      DO NA=2, NATOM
        IF (MINVAL(DD_VDOP(:ND,NA)) .LT. VDOPMIN) THEN 
           NA_VDOPMIN = NA  
           VDOPMIN = MINVAL(DD_VDOP(:ND,NA))  
        ENDIF
      ENDDO
      EL_MIN = SYMBOL(NA_VDOPMIN)

      WRITE(0,'(A,F6.1,A)') "Resolution corresponds to VDOP= ", VDOP,
     >                      ' (see output file for reason)'

C***  Give a warning if the minimum VDOP is smaller than VDOPFE      
C***  (In this case the iron lines might be broader than intended)
      IF (VDOP < VDOPFE) THEN
        WRITE (0,1) VDOP, VDOPFE
        WRITE (6,1) VDOP, VDOPFE
    1   FORMAT ('*** WARNING: FEDAT resolution is not fine enough.', /,
     >   '*** WARNING: min(VDOP) = 'F7.2, ', VDOPFE = ', F7.2, /, 
     >   '*** WARNING: Iron lines are therefore more Doppler-broadened',
     >   ' than others' )
      ENDIF

C***  Analog arrays in Doppler units
      DO L=1,ND
          DO NA=1,NATOM
                DD_VDOPDU(L,NA) = DD_VDOP(L,NA) / VDOP
                IF (.NOT. BMICROTURB) THEN
                  DD_VMIC(L) = DD_VDOP(L,1)
                ENDIF
                DD_VMICDU = DD_VMIC(L) / VDOP
          ENDDO
      ENDDO

C***  XMAX = minimum integration interval of lines - determined by maximum Doppler velocity!
      XMAX = AMAX1(XMAX, XMAXMIN*MAXVAL(DD_VDOPDU))   
         
C***  IVDOPSTATUS contains a code number, which criterion has caused 
C***  the VDOP that finally defines the wavelength resolution.
C***  This code number will lead to a corresponding message by PRIPRO 
C***  1: VDOP from MODEL file (default)
C***  2: minimum of all Doppler linewidths (accounting for Vmic) 
C***  3: specified as VDOP in FORMAL_CARDS
C***  4: VDOPFE from FEDAT file, when NO-FECONVOL requested
      IF (BMICROTURB) THEN
         IVDOPSTATUS = 2
      ELSEIF (IDX(DD_VDOP_LINE) .GT. 0) THEN
         IVDOPSTATUS = 3
      ENDIF
C***  minimum VDOP is not from ion, if NOIRINLINES or VDOPFE not minimum
      IF (BIRONLINES .AND. VDOPFE .LE. VDOP 
     >          .AND. .NOT. BDDFECONVOL) THEN      
         IVDOPSTATUS = 4
      ENDIF 

      RETURN
      
C***  ERROR branches *********************************************


 101  WRITE(0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE(0,'(A)') '*** The error occured in the following line:'
      WRITE(0,'(A)') DD_VDOP_LINE
      WRITE(0,'(A)') 'Version not known! possible versions are:' 
      WRITE(0,'(A)') 'VDOP X.X'
      WRITE(0,'(A)') 'VMIC X.X VELOFRAC [X.X]' 
      WRITE(0,'(A)') 'VDOP X.X MAX X.X ' // 
     >          '[TAU1 | R1 | VELO1 X.X TAU2 | R2 | VELO2 X.X]'
      STOP '*** Fatal error in vdop_struct'
      
 102  WRITE (0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE (0,'(A)') '*** The error occured in the following line:'
      WRITE (0,'(A)') DD_VDOP_LINE
      WRITE(0,'(A)') 'Version not known! possible versions are:' 
      WRITE(0,'(A)') 'VMIC X.X'    
      WRITE(0,'(A)') 'VMIC X.X VELOFRAC [X.X]' 
      WRITE(0,'(A)') 'VMIC X.X MAX X.X ' // 
     >       '[TAU1 | R1 |VELO1 X.X TAU2 | R2 | VELO2 X.X]'
      STOP '*** Fatal Error in vdop_struct'

      END


      SUBROUTINE VMALV (VA,VB,V,Q,LMAX)
C***********************************************************************
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
C***********************************************************************
      DIMENSION VA(LMAX),VB(LMAX),V(LMAX),Q(LMAX)

      LZ=LMAX-1
      Q(1)=VB(1)*V(1)
      DO 1 L=2,LZ
      Q(L)=VA(L)*V(L-1)+VB(L)*V(L)
    1 CONTINUE
      Q(LMAX)=VA(LMAX)*V(LZ)

      RETURN
      END
      SUBROUTINE VMF (V2,V1,A,N,NDIM)
C***********************************************************************
C***  MULTIPLICATION VECTOR = VECTOR * MATRIX (FULL)  --  V2 = V1 * A
C***********************************************************************
      DIMENSION V1(NDIM),V2(NDIM),A(NDIM,NDIM)

      DO 1 J=1,N
      SUM=.0
      DO 2 I=1,N
    2 SUM=SUM+V1(I)*A(I,J)
    1 V2(J)=SUM

      RETURN
      END
      FUNCTION VOIGTH(A,V)
C***********************************************************************
C***  normalized (!!!) VOIGT function H(a,v)
C***  ----- nach DETLEF KOESTER -----
C***********************************************************************
      DATA WPI /1.77245385/

      V=ABS(V)
      VV=V*V
      IF ((A.LT.0.2) .AND. (V.GT.5.0)) GO TO 70
      IF ((A.GE.0.2) .AND. (A.GT.1.4 .OR. (A+V).GT.3.2)) GO TO 80
      H0=EXP(-VV)
      H2=(1.-2.*VV)*H0
      IF (V.GT.2.4) GO TO 30
      IF (V.LE.1.3) GO TO 20
      H1=(-.220416*VV+1.989196*V-6.61487)*VV+9.39456*V-4.4848
      GO TO 40
   20 H1=(.42139*VV-2.34358*V+3.28868)*VV-.15517*V-1.1247
      GO TO 40
   30 H1=((-.0032783*VV+.0429913*V-.188326)*VV+.278712*V+.55415)/
     *   (VV-1.5)
   40 CONTINUE
      IF (A.GE..2) GO TO 60
   50 H=((H2*A+H1)*A+H0)
      VOIGTH=H/WPI
      RETURN

   60 HH1=H1+H0*1.12838
      HH2=H2+HH1*1.12838-H0
      HH3=(1.-H2)*0.37613-HH1*0.66667*VV+HH2*1.12838
      HH4=(3.*HH3-HH1)*0.37613+H0*0.66667*VV*VV
      H=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0)*
     *  (((-.122727278*A+.532770573)*A-.96284325)*A+.979895032)
      VOIGTH=H/WPI
      RETURN

   70 H=((2.12/VV+.8463)/VV+.5642)*A/VV
      VOIGTH=H/WPI
      RETURN

   80 AA=A*A
      U=(AA+VV)*1.4142
      UU=U*U
      H=((((AA-10.*VV)*AA*3.+15.*VV*VV)/UU+3.*VV-AA)/UU+1.)*A*.79788/U
      VOIGTH=H/WPI
      RETURN

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
      SUBROUTINE ZONEINT (LRED, LBLUE, EMINT, CEMINT, TAUSUM, TAUSUMC,
     >                   PJPJ, ZFINE, KINDEX,
     >                   OPAFINE, OPAFC, OPALFIN, ETAFINE, ETAFC,
     >                   ETALFIN, SFINE, CSFINE, 
     >                   RRAY, ZRAY, XCMF, 
     >                   OPARAY, OPALRAY, ETARAY, ETACRAY, ETALRAY,
     >                   MAXXN, LTOT, DELXLAP, NBLINE, NLOBLFI, NLOBLLA, 
     >                   NDADDIM, IVERSION, LINPRO, AVOIGT,
     >                   TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >                 XCMFFINE, POROLENGTHFINE, TAUMAX, XMAX, XMAXLIN,
     >                 DXMAX, BIRONLINES, OPAFE, ETAFE, NDDIM, NFLDIM, 
     >                 XCMFBLUE, XCMFRED, DXCMF, CORE, POROLENGTHRAY,
     >                 PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB, 
     >                 RADIUS, ND, MAXLAP, DD_VDOPDU_RAY, NATOM, 
     >                 IND_ELLINE,  
     >                 DD_VDOPDU_FINE_NORMFAC,
     >                 DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE,
     >                 GRIEMPAR, KODAT, VDOP, ZINTER, NMOD, MAXMOD, 
     >                 BTAUMAX_REACHED)

C***********************************************************************
C***  CALLED FROM: OBSFRAM
C***  FORMAL INTEGRATION ACROSS THE SCATTERING ZONE
C***  IVERSION = 0: INTEGRATION IN TAU - improved version 2-SEP-1998
C***  IVERSION = 1: INTEGRATION IN Z
C***********************************************************************

      IMPLICIT NONE

      REAL, DIMENSION(NDADDIM, NATOM), INTENT(IN) :: DD_VDOPDU_RAY
      REAL, DIMENSION(MAXXN, NATOM) :: DD_VDOPDU_FINE_SQRD, 
     >                                 DD_VDOPDU_FINE_NORMFAC, DD_VDOPDU_FINE
      INTEGER, DIMENSION (MAXLAP) :: IND_ELLINE
      REAL OPALRAY(NDADDIM,NBLINE), ETALRAY(NDADDIM,NBLINE)
      REAL OPALFIN(MAXXN,NBLINE), ETALFIN(MAXXN,NBLINE)
      REAL OPARAY(LTOT), ETARAY(LTOT), ETACRAY(LTOT)
      REAL RRAY(LTOT), ZRAY(LTOT), XCMF(LTOT), DELXLAP(NBLINE)  
      INTEGER IPOINTERPHITAB(NBLINE)
      REAL, DIMENSION(MAXLAP,NDDIM,MAXMOD) :: AVOIGT, GRIEMPAR
      REAL OPAFINE(MAXXN),ETAFINE(MAXXN),SFINE(MAXXN),ZFINE(MAXXN)
      REAL OPAFC(MAXXN), ETAFC(MAXXN), CSFINE(MAXXN)
      REAL TAU(MAXXN),DTAU(MAXXN),WTAU(MAXXN)
      REAL TAUC(MAXXN),DTAUC(MAXXN),WTAUC(MAXXN), XCMFFINE(MAXXN)
      REAL OPAFE(NDDIM,NFLDIM,MAXMOD), ETAFE(NDDIM,NFLDIM,MAXMOD)
      CHARACTER(8), DIMENSION(NBLINE) :: LINPRO 
      REAL, DIMENSION(NBLINE) :: XMAXLIN
      REAL POROLENGTHRAY(LTOT), POROLENGTHFINE(MAXXN)
      INTEGER, DIMENSION(NATOM) :: KODAT
      REAL ZINTER(2)      
      REAL, DIMENSION
     > (-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB,MAXMOD) :: PHITAB

      LOGICAL LASER, BIRONLINES, CORE, bHYROGEN, DEBUG
      LOGICAL BTAUMAX_REACHED

C***  next declarations needed only with IMPLICIT NONE
      REAL P, Q, WPIINV, XI, EPS, EPS2, DXFINE, DXMAX, DZFINE, RFINE
      REAL XCMFBLUE, XKCMF, XCMFRED, PJPJ, TEMP, DXCMF, WIF_P, WIF_M
      REAL ETAFEL, OPAFEL, OPAFEI, OPAFELM, ETAFELM,  ETAFEI
      REAL EMINT, CEMINT, RFINERIP, PHI, XI_VDOP, CORRFAC_MACROCLUMP
      REAL TAUSUM, TAUSUMC, WB, EXPTAUP, PR, QR, EXPTAU, EXPTAUM 
      REAL EXPTAUC, EXPTAUCP, EXPTAUCM, DZ, W, WAC, WBC
      REAL TAUMAX, VDOP, XMAX
      REAL RADIUS(NDDIM)
      REAL STARKHOLTSMARK, STARKPROF, SDOT, STARKVOIGT, WA

      INTEGER L, NFINE, LRIP, LMAX, LTOT, LRED, LBLUE, NPOINT, IPOINT 
      INTEGER I, NA, NATOM, KCMF_M, KCMF_P, NLOC, NLOBLFI, NLOBLLA
      INTEGER LL, IMOD, NMOD, IVERSION, LAST, KINDEX, ND, MAXXN, NBLINE
      INTEGER MAXMOD, NDADDIM, NDDIM, NFLDIM, NFDIMPHITAB, NLDIMPHITAB
      INTEGER MAXLAP

C***  WPIINV = 1. / SQRT(PI)
      DATA WPIINV / 0.564189583549 /
C***  The following two parameters are to make the new TAUVERSION 
C***     laser-resistant. 
      DATA EPS  / 1.E-10 /
      DATA EPS2 /  0.5   /
      SAVE EPS, EPS2

C***  DEBUG option; slows down the code!
      DEBUG = .FALSE. 

      LRIP = 1   ! initialization for depth interpolation (STARKPROF)
    
C***  LMAX is required because OPAFE, ETAFE was not mirrored at symmetry plane
      IF (CORE) THEN
         LMAX = LTOT
      ELSE
         LMAX = (LTOT+1) / 2
      ENDIF


C***  Speed up the formal: ZFINE only established at KINDEX=1
C***  All other interpolations are done (for each KINDEX) subsequently
      IF (KINDEX .EQ. 1) THEN

         DXFINE = DXMAX
      
C**      Loop on sub-intervals [L-1,L] of XCMF rough array.
C**      Refines every sub interval whose frequency-range is bigger than DXMAX 

C**      NFINE  - total number of inserted points in fine array
         NFINE = 0
      
      DO L=LRED+1,LBLUE
C***    NPOINT = number of fine points in current sub-interval [L-1,L] 
C***              including ZRAY(L-1), but excluding ZRAY(L)
        
C***    DXFINE defines the maximum Doppler shift that is allowed over an
C***    integration step along the ray. DXFINE must be small enough to resolve
C***    a Gaussian profile (0.3 Doppler units by default). When the *local*
C***    Doppler width is larger (e.g. by depth-dependent microturbulence or 
C***    Doppler width), this resolution can be relaxed according to the 
C***    finest resolution needed, as derived from the array DD_VDOPDU_RAY.
        DXFINE = DXMAX * MINVAL(DD_VDOPDU_RAY(L,1:NATOM))
        NPOINT = MAX(1, NINT(ABS(XCMF(L-1)-XCMF(L))/DXFINE))
    
        IF (NFINE+NPOINT .GE. MAXXN) THEN
            WRITE (0,'(A,I5,A,I5)') 
     >       'DIMENSION MAXXN=', MAXXN, ' INSUFFICIENT', MAXXN 
            STOP 'ERROR IN ZONEINT'
        ENDIF

        DZFINE = (ZRAY(L-1) - ZRAY(L)) / FLOAT(NPOINT)
        DO IPOINT=1, NPOINT
           ZFINE(NFINE+IPOINT) = ZRAY(L-1) - (IPOINT-1) * DZFINE
        ENDDO        
        NFINE = NFINE + NPOINT
      
      ENDDO
C***  End of Loop over the coarse intervals [L,L+1]
C***  Add endpont of last interval
      NFINE = NFINE+1
      ZFINE(NFINE) = ZRAY(LBLUE)

      ENDIF 
C**************************************************************

c      if (nfine .le. 1) write (0,*) '!!!! nfine .le 1 !!!!' 

c      DO i=1, nfine
c         write (0,*) 'i, zfine =', i, zfine(i)
c      enddo
c
c      stop 'test'

ccccccccccccccccccc begin steinbruch ccccccccccccccccccccccc      DXFINE = DXMAX
      
C** Interpolates the other arrays (ETA, OPA...) 

      L=2 ! meaning that ZRAY(L-1) > ZFINE(I) > ZRAY(L)
      DO I=1, NFINE
         IF (ZFINE(I) .LT. ZRAY(L)) L = L + 1
         CALL SPLINPO_FAST
     >        (XI, ZFINE(I), XCMF, ZRAY, LTOT, L, DEBUG)
            XCMFFINE(I) = XI
c            write (0,'(A,I3,F10.5,I3,F10.5,I5)') 'i, zfine, L, xcmf =', 
c     >           i, zfine(i), L, xi, KINDEX

C***     INTERPOLATION WEIGHT, LINEARLY IN RADIUS R
         RFINE=SQRT(PJPJ+ZFINE(I)*ZFINE(I))
         IF (RRAY(L) .NE. RRAY(L-1)) THEN
           P=(RFINE-RRAY(L-1))/(RRAY(L)-RRAY(L-1))
         ELSE
           P = 0.5
         ENDIF
         Q=1.-P

         POROLENGTHFINE(I) = 
     >      Q * POROLENGTHRAY(L-1) + P * POROLENGTHRAY(L)
          
C***     Interpolation of depth dependent VDOP arrays (linearly in radius)
         DO NA=1, NATOM
           DD_VDOPDU_FINE(I, NA) = 
     >         Q * DD_VDOPDU_RAY(L-1, NA) + P * DD_VDOPDU_RAY(L, NA)
           DD_VDOPDU_FINE_SQRD(I,NA) = 
     >         DD_VDOPDU_FINE(I,NA) * DD_VDOPDU_FINE(I,NA)
C***  This array contains the normalization factors for the gaussians
           DD_VDOPDU_FINE_NORMFAC(I,NA) = WPIINV / DD_VDOPDU_FINE(I,NA)
          ENDDO
         
C***      INTERPOLATION WEIGHTS MODIFIED: LINEAR INTERPOLATION IN R**2
          Q=Q*RRAY(L-1)*RRAY(L-1)/(RFINE*RFINE)
          P=P*RRAY(L  )*RRAY(L  )/(RFINE*RFINE)

          TEMP = Q * OPARAY(L-1) + P * OPARAY(L)
          OPAFINE(I) = TEMP
          OPAFC  (I) = TEMP

          ETAFINE(I) = Q * ETARAY (L-1) + P * ETARAY (L)
          ETAFC  (I) = Q * ETACRAY(L-1) + P * ETACRAY(L)

          
C***  Add Iron Opacities
          IF (BIRONLINES) THEN
C***        Find CMF-Frequency Index of the current fine frequency XI
C***        and interpolate OPAFE, ETAFE for that frequency
C***        at depth points (L-1) and (L)

C***        Frequency Indices and Interpolation Weights
            XKCMF  = (XCMFBLUE -  XI) / DXCMF + 1. 
            IF (XCMFBLUE .LT. XI) THEN
                WRITE (0,*) '*** INTERNAL ERROR DETECTED IN ZONEINT:'
                WRITE (0,*) '*** XCMFBLUE .LT. XI'
                WRITE (0,*) '*** XCMFBLUE=', XCMFBLUE
                WRITE (0,*) '*** XI      =', XI
                WRITE (0,*) '*** Frequency range in FORMCMF too small??'
                STOP        '*** FATAL ERROR ***'
            ENDIF  
            IF (XCMFRED .GT. XI) THEN            
                WRITE (0,*) '*** INTERNAL ERROR DETECTED IN ZONEINT:'
                WRITE (0,*) '*** XCMFRED .GT. XI'
                WRITE (0,*) '*** XCMFRED=', XCMFRED
                WRITE (0,*) '*** XI      =', XI
                WRITE (0,*) '*** Frequency range in FORMCMF too small??'
                STOP        '*** FATAL ERROR ***'
            ENDIF   
            KCMF_M = INT(XKCMF)
            KCMF_P = KCMF_M + 1
            WIF_P  = XKCMF - FLOAT(KCMF_M)
            WIF_M  = 1. - WIF_P

C***        Opacities at XKCMF at (L)
            IF (L .LE. LMAX) THEN
                LL = L
            ELSE
                LL = 2 * LMAX - L
            ENDIF

C***     Check if current point lies in second-model domain
            IMOD=1
            IF (NMOD .EQ. 2) THEN
               IF ((ZFINE(I)-ZINTER(1))*(ZFINE(I)-ZINTER(2)) .LT. .0)
     >            IMOD=2
            ENDIF

            OPAFEL  = WIF_M * OPAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * OPAFE(LL,KCMF_P,IMOD)
            ETAFEL  = WIF_M * ETAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * ETAFE(LL,KCMF_P,IMOD)
C***     Opacities at XKCMF at (L-1)
            IF (L-1 .LE. LMAX) THEN
                LL = L-1
            ELSE
                LL = 2 * LMAX - L-1
            ENDIF
            OPAFELM = WIF_M * OPAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * OPAFE(LL,KCMF_P,IMOD)
            ETAFELM = WIF_M * ETAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * ETAFE(LL,KCMF_P,IMOD)
            
C***        Final Radius-Interpolation
            OPAFEI = Q * OPAFELM + P * OPAFEL
            ETAFEI = Q * ETAFELM + P * ETAFEL
C***        Adding up Iron to the total Opacity/Emissivity
            OPAFINE(I) = OPAFINE(I) + OPAFEI
            ETAFINE(I) = ETAFINE(I) + ETAFEI
         ENDIF
C***     -- end of iron branch ---

         DO NLOC=NLOBLFI,NLOBLLA
            OPALFIN(I,NLOC) = Q * OPALRAY(L-1,NLOC) +
     >                          P * OPALRAY(L,NLOC)
            ETALFIN(I,NLOC) = Q * ETALRAY(L-1,NLOC) +
     >                          P * ETALRAY(L,NLOC)
         ENDDO
      ENDDO
C***  End of Loop over all fine depth points

C***  Macroclumping correction for the continuum opacity
      DO I = 1, NFINE
       CALL MACROCLUMP (OPAFC(I), POROLENGTHFINE(I), CORRFAC_MACROCLUMP)
       OPAFC(I) = OPAFC(I) * CORRFAC_MACROCLUMP
       ETAFC(I) = ETAFC(I) * CORRFAC_MACROCLUMP
      ENDDO
      
C******************************************************************
C***  Z version
C***  NOTE: This version has no problems with LASER effects.
C***        The weight function exp(-tau) is NOT incorporated in the 
C***        quadrature weights. This implies a loss of accuracy.
C******************************************************************
      IF (IVERSION .EQ. 1) THEN
C***  Opacities and optical depths ...
      DO I=1, NFINE
C***     This setting of RFINE assures that L-interpolation weigts
C***     have not been calculated yet (for Subr. STARKPROF)
         RFINERIP = .0

C***     Check if current point lies in second-model domain
         IMOD=1
         IF (NMOD .EQ. 2) THEN
            IF ((ZFINE(I)-ZINTER(1))*(ZFINE(I)-ZINTER(2)) .LT. .0) 
     >         IMOD=2
         ENDIF

C***     ADD LINE PLUS CONTINUUM OPACITIES
         DO NLOC=NLOBLFI,NLOBLLA
            NA = IND_ELLINE(NLOC) !NA ist the element index corresponding to the line
            XI = XCMFFINE(I) - DELXLAP(NLOC)
C***        Max. bandwidth holds for all profile types!
            IF (ABS(XI) .GT. XMAXLIN(NLOC)) CYCLE
            IF (LINPRO(NLOC) .EQ. '        ') THEN
               PHI = EXP(-XI*XI / DD_VDOPDU_FINE_SQRD(I,NA))
               PHI = PHI * DD_VDOPDU_FINE_NORMFAC(I,NA)
               OPAFINE(I)=OPAFINE(I)+PHI*OPALFIN(I,NLOC)
               ETAFINE(I)=ETAFINE(I)+PHI*ETALFIN(I,NLOC)

            ELSE IF (LINPRO(NLOC) .EQ. 'VOIGT   '
     >          .OR. LINPRO(NLOC) .EQ. 'BRD-HeI '
     >          .OR. LINPRO(NLOC) .EQ. 'Q-STARK ') THEN
C***           XI Argument of VOIGTH needs to account for VDOP
               XI_VDOP = XI/DD_VDOPDU_FINE(I,NA)
C***           STARKVOIGT interpolates appropriate value for AV 
C***             as function of radius and for primary and second model, 
C***             then calls the Voigt function VOIGTH 
               PHI = STARKVOIGT (XI_VDOP, AVOIGT(1,1,IMOD), NLOC, MAXLAP, 
     >               ZFINE(I), LRIP, RFINERIP, PR, QR, RADIUS, ND, 
     >               PJPJ)
C***           Profile is normalized appropriately
               PHI = PHI / DD_VDOPDU_FINE(I,NA)
               OPAFINE(I) = OPAFINE(I) + PHI * OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI * ETALFIN(I,NLOC)
            ELSE IF (LINPRO(NLOC) == 'L-STARK ') THEN
C***           Stark broadening: depth-dependent HOLTSMARK PROFILE
               bHYROGEN = (NA == KODAT(1))
C***           STARKHOLTSMARK interpolates appropriate value for GRIEMPAR
               PHI = STARKHOLTSMARK (XI, VDOP, DD_VDOPDU_FINE(I,NA),
     >                               GRIEMPAR(1,1,IMOD), NLOC, MAXLAP, 
     >                               ZFINE(I), LRIP, RFINERIP, PR, QR,
     >                               RADIUS, ND, PJPJ, bHYROGEN)
               OPAFINE(I) = OPAFINE(I) + PHI * OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI * ETALFIN(I,NLOC)
            ELSE IF (LINPRO(NLOC)(:3) .EQ. 'BRD') THEN
C***              PRESSURE BROADENING: Hydrogen and Helium / tabulated
C*** If depth-dependent VDOP active, data is tabulated appropriately in STARKBROAD 
               PHI = STARKPROF (XI, IPOINTERPHITAB(NLOC), ZFINE(I),
     >                  LRIP, RFINERIP, PR, QR,
     >                  PHITAB(-NFDIMPHITAB,1,1,IMOD),
     >                  NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, NDDIM,
     >                  PJPJ, DXMAX)
               OPAFINE(I) = OPAFINE(I) + PHI*OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI*ETALFIN(I,NLOC)
            ELSE
               WRITE (0, '(2A)') 
     >          '*** UNDEFINED LINE PROFILE TYPE: ', LINPRO(NLOC)
               STOP '*** FATAL ERROR IN ZONEINT'
            ENDIF
         ENDDO
 
C***  Macroclumping correction for the total opacity
      CALL MACROCLUMP (OPAFINE(I), POROLENGTHFINE(I), 
     >                  CORRFAC_MACROCLUMP)
      OPAFINE(I) = OPAFINE(I) * CORRFAC_MACROCLUMP
      ETAFINE(I) = ETAFINE(I) * CORRFAC_MACROCLUMP

C***  ESTABLISH DTAU(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAU(I) = OPTICAL DEPTH SCALE
         IF (I .EQ. 1) THEN
            DTAU (1) = 0.
            DTAUC(1) = 0.
            TAU  (1) = TAUSUM
            TAUC (1) = TAUSUMC
         ELSE
            DTAU (I) = 0.5 * 
     >              (OPAFINE(I-1)+OPAFINE(I)) * (ZFINE(I-1)-ZFINE(I))  
            DTAUC(I) = 0.5 * 
     >              (OPAFC  (I-1)+OPAFC  (I)) * (ZFINE(I-1)-ZFINE(I))

            TAU (I) = TAU (I-1) + DTAU (I)
            TAUC(I) = TAUC(I-1) + DTAUC(I)
         ENDIF
      ENDDO   !end of fine frequency loop

      TAUSUM  = TAU (NFINE)
      TAUSUMC = TAUC(NFINE)


C***      S means: EMISSIVITY * EXP(-TAU)
          DO 14 I=1, NFINE
            IF (TAU(I) .LT. 700.) THEN
               SFINE(I) = ETAFINE(I) * EXP(-TAU (I))
            ELSE
              SFINE(I) = 0.
            ENDIF
            IF (TAUC(I) .LT. 700.) THEN
              CSFINE(I) = ETAFC  (I) * EXP(-TAUC(I))
            ELSE
              CSFINE(I) = 0.
            ENDIF
   14     CONTINUE

          WTAU (1) = 0.5 * (ZFINE(1) - ZFINE(2))
          WTAUC(1) = WTAU(1)
          DO 7 I = 2, NFINE-1
            WTAU (I) = 0.5 * (ZFINE(I-1) - ZFINE(I+1))
            WTAUC(I) = WTAU(I)
    7     CONTINUE
          WTAU (NFINE) = 0.5 * (ZFINE(NFINE-1) - ZFINE(NFINE))
          WTAUC(NFINE) = WTAU(NFINE)

C***  INTEGRATION SUM, USING CRAY VECTOR FUNCTION SDOT (SCALAR PRODUCT)
      CEMINT = CEMINT + SDOT(NFINE,CSFINE,1,WTAUC,1)
       EMINT =  EMINT + SDOT(NFINE, SFINE,1,WTAU ,1)



C******************************************************************
C***  TAU version (NEW! 2-SEP-1998, wrh)
C***        This version combines the advantages of both the (old) TAU 
C***        and the Z version, as it is accurate AND laser-resistant. 
C***        The weight function exp(-tau) is incorporated in the 
C***        quadrature weights for better accuracy (equivalent to the 
C***        TAU version) 
C***        In LASER intervals, the normal Z integration is restored. 
C***        This is still possible, because the divison by opa (S=ETA/OPA)
C***        is not performed in advance, but ETA is kept and 1/OPA is 
C***        implied in the weights. For this switching of the integretion 
C***        methods, it is necessary to apply the endpoint-terms exp(-tau)
C***        in the integration weights each time, although they cancel out 
C***        in case of subsequent tau-version intervals
C******************************************************************
      ELSE IF (IVERSION .EQ. 0) THEN
         EPS = 1.E-10

C***     LINE + CONTINUUM ***************************************
 
C***  ESTABLISH DTAU(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAU(I) = OPTICAL DEPTH SCALE
C***     AND CALCULATE INTEGRATION WEIGHTS EXP(-TAU) DTAU
C***     THE LOOP RUNS FOR I=0 FIRST IN ORDER TO PREPARE LINE OPACITIES AT I=1

         TAU(1) = TAUSUM
         EXPTAU = EXP(-TAU(1))
         DO 5 I = 0, NFINE-1
C***        This setting of RFINE assures that L-interpolation weigts 
C***        have not been calculated yet (for Subr. STARKPROF)
            RFINERIP = .0

C***        Check if current point lies in second-model domain
            IMOD=1
            IF (NMOD .EQ. 2) THEN
               IF ((ZFINE(I+1)-ZINTER(1))*(ZFINE(I+1)-ZINTER(2)) .LT. .0) 
     >            IMOD=2
            ENDIF

C***        ADD LINE PLUS CONTINUUM OPACITIES
            DO NLOC=NLOBLFI,NLOBLLA
               NA = IND_ELLINE(NLOC) !NA ist the element index corresponding to the line
               XI = XCMFFINE(I+1) - DELXLAP(NLOC)
C***           Max. bandwidth holds for all profile types!
cc               IF (ABS(XI) .GT. XMAX) CYCLE
               IF (ABS(XI) .GT. XMAXLIN(NLOC)) CYCLE
               IF (LINPRO(NLOC) .EQ. '        ') THEN
C***              DOPPLER BROADENING ONLY: GAUSS PROFILE
                  PHI = EXP(-XI*XI / DD_VDOPDU_FINE_SQRD(I+1,NA))
                  PHI = PHI * DD_VDOPDU_FINE_NORMFAC(I+1,NA)
                  OPAFINE(I+1) = OPAFINE(I+1) 
     >                         + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) 
     >                         + PHI * ETALFIN(I+1,NLOC)

               ELSE IF (LINPRO(NLOC) .EQ. 'VOIGT   '
     >             .OR. LINPRO(NLOC) .EQ. 'BRD-HeI '
     >             .OR. LINPRO(NLOC) .EQ. 'Q-STARK ') THEN
                  XI_VDOP = XI/DD_VDOPDU_FINE(I+1,NA)
                  PHI = STARKVOIGT (XI_VDOP, AVOIGT(1,1,IMOD), NLOC, MAXLAP, 
     >                  ZFINE(I+1), LRIP, RFINERIP, PR, QR, RADIUS, ND, 
     >                  PJPJ)
C***              Profile is normalized appropriately
                  PHI = PHI / DD_VDOPDU_FINE(I+1,NA)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI * ETALFIN(I+1,NLOC)

               ELSE IF (LINPRO(NLOC) == 'L-STARK ') THEN
C***              Stark broadening: depth-dependent HOLTSMARK PROFILE
                  bHYROGEN = (NA == KODAT(1))
C***              STARKHOLTSMARK interpolates appropriate value for GRIEMPAR
                  PHI = STARKHOLTSMARK (XI, VDOP, DD_VDOPDU_FINE(I+1,NA),
     >                              GRIEMPAR(1,1,IMOD), NLOC, MAXLAP, 
     >                              ZFINE(I+1), LRIP, RFINERIP, PR, QR,
     >                              RADIUS, ND, PJPJ, bHYROGEN)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI * ETALFIN(I+1,NLOC)
               ELSE IF (LINPRO(NLOC)(:3) .EQ. 'BRD') THEN
C***              PRESSURE BROADENING: H I and He II / tabulated
C*** If depth-dependent VDOP active, data is tabulated appropriately in STARKBROAD 
                  PHI = STARKPROF (XI, IPOINTERPHITAB(NLOC), ZFINE(I+1), 
     >                  LRIP, RFINERIP, PR, QR, 
     >                  PHITAB(-NFDIMPHITAB,1,1,IMOD),
     >                  NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, NDDIM,
     >                  PJPJ, DXMAX)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI*OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI*ETALFIN(I+1,NLOC)
               ELSE
                  WRITE (0, '(2A)') 
     >             '*** UNDEFINED LINE PROFILE TYPE: ', LINPRO(NLOC)
                  STOP '*** FATAL ERROR IN ZONEINT'
               ENDIF
            ENDDO

C***  Macroclumping correction for the total opacity
           CALL MACROCLUMP (OPAFINE(I+1), POROLENGTHFINE(I+1), 
     >                  CORRFAC_MACROCLUMP)
           OPAFINE(I+1) = OPAFINE(I+1) * CORRFAC_MACROCLUMP
           ETAFINE(I+1) = ETAFINE(I+1) * CORRFAC_MACROCLUMP
         
            IF (I .EQ. 0) GOTO 5

            DTAU (I+1) = 0.5 * 
     >              (OPAFINE(I)+OPAFINE(I+1)) * (ZFINE(I)-ZFINE(I+1))
            TAU (I+1) = TAU (I) + DTAU (I+1)


C***        WB COMES FROM THE INTERVAL I-1, I (ZERO for first interval!)
            IF (I .EQ. 1) THEN
               WB = .0
               EXPTAUP  = EXP(-TAU (2))
            ELSE
               EXPTAUM  = EXPTAU
               EXPTAU   = EXPTAUP
               EXPTAUP  = EXP(-TAU(I+1))
C***           LASER condition for that interval (LINE)
               LASER = ABS(DTAU (I)) .LT. EPS  .OR. 
     >              OPAFINE(I-1) .LT. EPS2 * OPAFC(I-1)    .OR.
     >              OPAFINE(I  ) .LT. EPS2 * OPAFC(I)
               IF (LASER) THEN 
                  DZ = ZFINE(I-1) - ZFINE(I)
                  WB =  EXPTAU  * DZ * 0.5
               ELSE            
                  WB = (-EXPTAU - (EXPTAU  - EXPTAUM ) / DTAU (I)) 
     >                 / OPAFINE(I)
               ENDIF
            ENDIF

C***        WA COMES FROM THE INTERVAL I, I+1
C***        LASER condition for that interval (LINE)
            LASER = ABS(DTAU (I+1)) .LT. EPS  .OR. 
     >           OPAFINE(I+1) .LE. EPS2 * OPAFC(I+1)    .OR.
     >           OPAFINE(I  ) .LE. EPS2 * OPAFC(I  )
            IF (LASER) THEN
               DZ = ZFINE(I) - ZFINE(I+1)
               WA  = EXPTAU  * DZ * 0.5
            ELSE
               WA  = (EXPTAU + (EXPTAUP - EXPTAU ) / DTAU (I+1)) 
     >               /  OPAFINE(I)
            ENDIF
          
            W = WA  + WB
            EMINT = EMINT + ETAFINE(I) * W 

            IF (TAU(I+1) .GT. TAUMAX) THEN
               LAST = I + 1
               BTAUMAX_REACHED = .TRUE.
               GOTO 51
            ENDIF 

    5    CONTINUE

C***     Final Interval N-1, N (or LAST-1, LAST, if TAUMAX was reached)

         LAST = NFINE
   51    CONTINUE


C***     LASER condition for that interval (LINE)
         LASER = ABS(DTAU (LAST)) .LT. EPS  .OR. 
     >           OPAFINE(LAST-1) .LT. EPS2 * OPAFC(LAST-1)  .OR.
     >           OPAFINE(LAST  ) .LT. EPS2 * OPAFC(LAST  )
         IF (LASER) THEN 
            DZ = ZFINE(LAST-1) - ZFINE(LAST)
            W = EXPTAUP * DZ * 0.5
         ELSE
            W = (-EXPTAUP - (EXPTAUP-EXPTAU)/DTAU(LAST) ) 
     >          / OPAFINE(LAST) 
         ENDIF

         EMINT = EMINT + ETAFINE(LAST) * W 
         TAUSUM  = TAU(LAST)


C***     ONLY CONTINUUM  *******************************************

C***  ESTABLISH DTAUC(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAUC(I) = OPTICAL DEPTH SCALE
C***         Note: This calculation must hurry one step ahead
C***     AND CALCULATE INTEGRATION WEIGHTS EXP(-TAU) DTAU


         TAUC (1) = TAUSUMC
         EXPTAUC  = EXP(-TAUC(1))

         DO 55 I = 1, NFINE-1

            DTAUC(I+1) = 0.5 * 
     >              (OPAFC(I)+OPAFC(I+1)) * (ZFINE(I)-ZFINE(I+1))
            TAUC(I+1) = TAUC(I) + DTAUC(I+1)

C***        WB COMES FROM THE INTERVAL I-1, I (ZERO for the first interval!)
            IF (I .EQ. 1) THEN
               WBC = .0
               EXPTAUCP = EXP(-TAUC(2))
            ELSE
               EXPTAUCM = EXPTAUC
               EXPTAUC  = EXPTAUCP
               EXPTAUCP = EXP(-TAUC(I+1))
C***           LASER condition for that interval (Continuum)
               LASER = ABS(DTAUC (I)) .LT. EPS  .OR. 
     >              OPAFC(I-1) .LT. .0  .OR.   
     >              OPAFC(I  ) .LT. .0    
               IF (LASER) THEN 
                  DZ = ZFINE(I-1) - ZFINE(I)
                  WBC = EXPTAUC * DZ * 0.5
               ELSE            
                  WBC = (-EXPTAUC - (EXPTAUC - EXPTAUCM) / DTAUC(I)) 
     >                  / OPAFC  (I)
               ENDIF
            ENDIF

C***        WA COMES FROM THE INTERVAL I, I+1
C***        LASER condition for that interval (Continuum)
            LASER = ABS(DTAUC (I+1)) .LT. EPS  .OR. 
     >           OPAFC(I+1) .LT. 0.  .OR.   
     >           OPAFC(I  ) .LT. 0.    
            IF (LASER) THEN
               DZ = ZFINE(I) - ZFINE(I+1)
               WAC = EXPTAUC * DZ * 0.5
            ELSE
               WAC = (EXPTAUC + (EXPTAUCP - EXPTAUC) / DTAUC(I+1)) 
     >               / OPAFC  (I)
            ENDIF
          
            W = WAC + WBC
            CEMINT = CEMINT + ETAFC(I) * W 

         IF (TAUC(I+1) .GT. TAUMAX) THEN
            LAST = I + 1
            BTAUMAX_REACHED = .TRUE.
            GOTO 52
         ENDIF 

   55    CONTINUE

C***     Final Interval N-1, N (or LAST, if TAUMAX is reached)

         LAST = NFINE

   52    CONTINUE

C***     LASER condition for that interval (Continuum)
         LASER = ABS(DTAUC (LAST)) .LT. EPS  .OR. 
     >           OPAFC(LAST-1) .LE. 0.  .OR.   
     >           OPAFC(LAST  ) .LE. 0.    

         IF (LASER) THEN 
            DZ = ZFINE(LAST-1) - ZFINE(LAST)
            W = EXPTAUCP * DZ * 0.5
         ELSE
            W = (-EXPTAUCP - (EXPTAUCP - EXPTAUC) / DTAUC(LAST)) 
     >          / OPAFC(LAST) 
         ENDIF
         CEMINT = CEMINT + ETAFC(LAST) * W 
         TAUSUMC = TAUC(LAST)

C***********************************************************************
      ELSE
          STOP 'IVERSION INVALID IN ZONEINT'
      ENDIF

      RETURN
      END
