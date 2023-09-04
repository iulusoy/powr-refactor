      SUBROUTINE FREQUINT (K, FWEIGHTL, XLAMK, XLAMKOLD,    
     >                   XLAMBDA, DXMAX, VELO, GRADI,
     >                   XJL, XJLOLD, XHL, XHLOLD, 
     >                   OPAKOLD, NF, ND, NDDIM, NATOM,
     >                   KONTACT, KONTAUP, DFKONT, LIND, PWEIGHT, 
     >                   MAXLIN, INDFEACT, SIGMAACT, MAXFEACT, LASTFE,
     >                   OPAK, ETAK, OPAFEI, ETAFEI, IFENUP, IFELOW,
     >                   DJDS, DJDS_OLD, DJDSraw, 
     >                   bALOTri, DJDSU, DJDSL,
     >                   WEIGHT, N, POPNUM, XKLOLD,XNLOLD,
     >                   OPAKNOTH, ETAKNOTH, OPAKNOTHO, ETAKNOTHO,
     >                   EDDIFO, EDDIGO, EPSGO, RADIUS, BCOLIRAY,
     >                   ETANOTH, OPA, THOMSON, T, ETAKOLD,
     >                   QLFOLD, QLHOLD, OPAKHOLD,
     >                   EDDIHOUTOLD, EDDIHINOLD, XHID,
     >                   FWEIGHT, OPAO, THOMSONO,
     >                   OPAKOLDELEM, OPACOLDELEM, ETAKOLDELEM, 
     >                   ETACOLDELEM, OPAKOLDION, OPACOLDION, 
     >                   ETAKOLDION, ETACOLDION, OPAFE,
C***  Temporary vectors
     >                   SUMJ, SUMJW, SUMDJDSC, SUMDJDSCW, 
     >                   SUMOPA, SUMOPAW, SUMETA, SUMETAW, SC, SCO, 
C***             OUTPUT:
     >                   XJCINT, FWTEST, XJLMEAN, XJFEMEAN, XLAMAPPMEAN, HTOTL,
     >                   HTOTMINUSND, HTOTND, HTOTNDCOR,
     >                   HTOTNDS, ARAD, ACONT, ATHOM, ARADELEM, 
     >                   ACONTELEM, ARADION, ACONTION,
     >                   iIgnoreK, XJTOTL, XKTOTL, XNTOTL, 
     >                   WFELOW, WFENUP, FTCOLI, FTCOLIB,
     >                   DSDFELOW, DSDFENUP, WJC, WJC_MIN, 
     >                   DSDSC, DJDSC, DJDSCOLD, OPAKINT, ETAKINT,
     >                   DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
C***       Unsoeld-Lucy:
     >                   OPASMEAN, OPASMEANTC, SMEAN, QFJMEAN, OPAJMEAN,
     >                   OPAJMEANTC, OPAPMEAN, QOPAHMEAN, HMEAN, 
     >                   EDDIHOUTJMEAN, XHOM, HTOTOUTMINUS,
     >                   RADIUS2, OPC, 
     >                   FERATUL, FERATLU, ELEVEL, 
     >                   EMCOLI, FTFE, IVERS_FE_EXPFAC, LPLOT_WCHARM, 
     >                   GAMMACOLI, OPAROSS, OPAROSSCONT, OPALAMBDAMEAN, 
     >                   GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2, TEFF,
C*** die folgenden SKALAREN Parameter werden ausgereicht, um den
C*** Compiler-Bug seit Compaq Tru64 UNIX V5.0A (Ende 2000) zu umgehen
     >                   XNUEK, XNUEKOLD, XNUEKOLDOLD, VDOPUNIT,
     >                   IFRBSTA, IFRBEND, VDOPFE, DXFE, XLAM0FE, 
     >                   NCHARG, MAXION, POPMIN, 
     >                   TAULAST, XHLOLDTEST, HTOTLTEST, iHTOTCUT, 
     >                   HTOTCUT, bOSKIPLAST, XLAMLASTOSKIP, XHBFSKIP, 
     >                   LPRDH, bDEBUG, CUTOPAMEANTHRES)
C***  CLIGHT = SPEED OF LIGHT IN ANGSTROEM/SECOND
      DATA CLIGHT / 2.99792458E18 /
      DATA PI8 /25.1327412288 /
C***  C1 = H * C / K        (CM * KELVIN)
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /

      DATA WJMAX / 0.9999999999 /

      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: WPI =  1.772454        !WPI = SQRT(PI)
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5  !Stephan-Boltzmann constant (CGS) / Pi
      REAL, PARAMETER :: CLIGHTKM = 2.9979E5    !CLIGHT = SPEED OF LIGHT IN KILOMETER/SECOND
      
      REAL, PARAMETER :: EXPMAX = 500.

      INTEGER, INTENT(IN) :: ND, NDDIM, N, NF, MAXLIN, MAXFEACT, LASTFE,
     >                       NATOM, LPRDH
      REAL, INTENT(IN) :: TEFF
      
      REAL, DIMENSION(ND) :: RADIUS, RADIUS2, XJL, XJLOLD, XHL, XHLOLD,
     >                       XKLOLD, XNLOLD, HTOTL, XJTOTL, XKTOTL, 
     >                       XNTOTL, ARAD, ACONT, ATHOM, VELO, GRADI,
     >                       OPAKOLD, ETAKOLD, OPAKFEOLD, OPAFE,
     >                       OPAROSSCONT, T, OPARFLUX, OPARCFLUX,
     >                       OPAK, ETAK, OPAKNOTH, ETAKNOTH, ETANOTH,
     >                       OPA, THOMSON, OPAO, THOMSONO, OPAKNOTHO,
     >                       ETAKNOTHO, DSDETA, DSDOPA, DJDS, DJDO,
     >                       DJDS_OLD, DJDO_OLD, DSDSC, DJDSC, DJDSCOLD,
     >                       FTCOLI, SC, SCO, TAULAST,
     >                       OPAKNOFENOTHO, ETAKNOFENOTHO, FTCOLIB,
     >                       OPAKFEFTOLD, ETAKFEFTOLD, OPACTOT, OPALTOT,
     >                       XHLOLDTEST, HTOTLTEST, HTOTCUT, XHBFSKIP,
     >                       XLAMLASTOSKIP, DJDSraw
      INTEGER, DIMENSION(MAXLIN) :: LIND
      REAL, DIMENSION(NDDIM, MAXLIN) :: PWEIGHT, OPAL
      INTEGER, DIMENSION(MAXFEACT) :: INDFEACT
      REAL, DIMENSION(MAXFEACT) :: SIGMAACT

      REAL, DIMENSION(ND, NF) :: XJCINT, OPAKINT, ETAKINT, WJC, WJC_MIN
      REAL, DIMENSION(NF) :: XLAMBDA, FWTEST, FWEIGHT, EMCOLI
      INTEGER, DIMENSION(ND) :: iIgnoreK
      REAL, DIMENSION(ND, LASTFE) :: XJFEMEAN
      REAL, DIMENSION(NDDIM,MAXLIN) :: XJLMEAN, XLAMAPPMEAN,
     >                                 XLAMAPPUMEAN, XLAMAPPLMEAN

      INTEGER, DIMENSION(LASTFE) :: IFENUP, IFELOW
      REAL, DIMENSION(ND,LASTFE) :: OPAFEI, ETAFEI, WFELOW, WFENUP,
     >                              DSDFELOW, DSDFENUP, FTFE
      REAL, DIMENSION(N) :: WEIGHT
      REAL, DIMENSION(ND,N) :: POPNUM

      REAL, DIMENSION(NATOM, ND) :: OPAKELEM, OPAKOLDELEM,
     >                              OPACELEM, OPACOLDELEM,
     >                              ETAKELEM, ETAKOLDELEM,
     >                              ETACELEM, ETACOLDELEM,
     >                              OPAROSSELEM, OPATOTELEM
      REAL, DIMENSION(ND, NATOM, MAXION) :: OPAKION, OPAKOLDION,
     >                                      OPACION, OPACOLDION,
     >                                      ETAKION, ETAKOLDION,
     >                                      ETACION, ETACOLDION
      REAL, DIMENSION(NATOM, ND-1) :: ARADELEM, ACONTELEM  
      REAL, DIMENSION(ND-1, NATOM, MAXION) :: ARADION, ACONTION

C***  Unsoeld-Lucy
      DIMENSION OPASMEAN(ND), SMEAN(ND), QFJMEAN(ND), ALAMBDAMEAN(ND)
      DIMENSION OPAJMEAN(ND), OPAPMEAN(ND)
      REAL, DIMENSION(ND) :: QOPAHMEAN, QOPAHMEANL, HMEAN
      DIMENSION OPASMEANTC(ND), OPAJMEANTC(ND) 
      LOGICAL, DIMENSION(ND) :: bOSKIPLAST

      DIMENSION QLFOLD(ND), QLHOLD(ND), OPAKHOLD(ND)
      DIMENSION EDDIFO(ND), EDDIGO(ND), EPSGO(ND)
      LOGICAL BCOLIRAY, BNEWINTERVAL, BNEWINTERVALOLD, BEMPTY
      REAL, DIMENSION(ND) :: SUMJ, SUMJW, SUMDJDSC, SUMDJDSCW, 
     >                       SUMOPA, SUMOPAW, SUMETA, SUMETAW      
      CHARACTER(LEN=8) OPC
      REAL, DIMENSION(LASTFE,ND) :: FERATLU, FERATUL
      DIMENSION ELEVEL(LASTFE)
C***  ALI Acceleration terms
      REAL, DIMENSION(ND) :: OPAROSS, OPALAMBDAMEAN

      INTEGER, INTENT(IN) :: iHTOTCUT, IVERS_FE_EXPFAC
      LOGICAL :: BWCELAB, BCUTOFF, BCUTOFF2, bALOTri,
     >           bIgnoreK, bOSKIP, bDEBUG, bKSKIPMEAN
      DIMENSION IFRBSTA(LASTFE), IFRBEND(LASTFE)
      REAL :: VDOPFE, DXFE, XLAM0FE, bInBand, TL, XHLUSE, 
     >        HTOTMINUSND, HTOTND, HTOTNDCOR, XHI, XHID, 
     >        EDDIHINTOLD, HTOTNDS, HCURL, AF, CF,
     >        DELTAHNUE, VDR, DHMECH, TCUT, XHLCUT, DXMAX,
     >        CUTOPAMEANTHRES

      INTEGER, DIMENSION(N) :: NCHARG


      
      SAVE      

cccccc!!!!!!!!!! LTE test option: set radiation field to Blackbody!
ccc   Consistent change AT TWO OTHER PLACES in this subroutine !
ccc      dimension xjlsave(100), xjloldsave(100)
ccc      DO L=1, ND
ccc      xjlsave(l) = xjl(l)      
ccc      xjloldsave(l) = xjlold(l)      
ccc      XJL(L)    = amax1 ( 1.e-100, BNUE(XLAMK   , T(L)) * RADIUS2(L))
ccc      if (k .gt. 0) XJLOLD(L) = 
ccc     >      amax1 ( 1.e-100, BNUE(XLAMKOLD, T(L)) * RADIUS2(L))
ccc      ENDDO
cccccccccccccccccccccccccccccccccccc

C***  FIXED Switch for the 'elaborated' Calculation of WJC
      BWCELAB = .TRUE.

C***  Application of GAMMACOLI to DJDS
C***  --------------------------------
      DJDS = MIN(DJDS, WJMAX)
      DJDS = MAX(DJDS, 0.)

C***  store undamped DJDS       
      DJDSraw = DJDS
      
C***  DJDS -> TAU
      DJDS = -LOG(1.-DJDS)
C***  DJDS -> TAU/GAMMA
      DJDS = DJDS / GAMMACOLI
C***  Back Transformation
      DJDS = 1.-EXP(-DJDS)

C********************************************************************
C***  Lines: Integration of the scattering integrals XJLMEAN
C********************************************************************

C***  ADDING THE NEW XJL TO THE MEAN INTENSITY XJLMEAN (NORMAL LINES)
      DO NL=1, MAXLIN
        IF (LIND(NL) .EQ. 0) CYCLE
        DO L=1, ND
C***       Note: The factor FWEIGHTL was missing here since we divide
C***             by the integrated PWEIGHTs afterwards and can assume
C***             that FWEIGHTL is roughly constant over a line
C***             This has now been add, but its effect is minor
           XJLMEAN(L,NL) = XJLMEAN(L,NL) + XJL(L)*PWEIGHT(L,NL)*FWEIGHTL
C***       Freq.-integrated LAMBDASTAR operator for current line NL
C***       The ALO at the current frequency is weighted with the
C***       contribution of the considered line opacity (excluding iron) 
C***       (at the current freq.) over the sum of all line and true cont. contrib.
           OPALFRAC = (OPAK(L) - OPAFE(L) - OPA(L)) /OPAKNOTH(L)
           XLAMAPPMEAN(L,NL) = XLAMAPPMEAN(L,NL) 
     >       + DJDS(L) * OPALFRAC * PWEIGHT(L,NL)*FWEIGHTL
           IF (bALOTri) THEN
             IF (L < ND) THEN
               OPALFRAC = (OPAK(L+1) - OPAFE(L+1) - OPA(L+1)) 
     >           /OPAKNOTH(L+1)
               XLAMAPPUMEAN(L,NL) = XLAMAPPUMEAN(L,NL) 
     >           + DJDSU(L) * OPALFRAC * PWEIGHT(L,NL)*FWEIGHTL
             ENDIF
             IF (L > 1) THEN
               OPALFRAC = (OPAK(L-1) - OPAFE(L-1) - OPA(L-1)) 
     >           /OPAKNOTH(L-1)
               XLAMAPPLMEAN(L,NL) = XLAMAPPLMEAN(L,NL) 
     >           + DJDSL(L) * OPALFRAC * PWEIGHT(L,NL)*FWEIGHTL
             ENDIF
           ENDIF
        ENDDO
      ENDDO

C***  ADDING THE NEW XJL TO THE MEAN INTENSITY XJFEMEAN (IRON LINES)
C***     AND CALCULATION OF DIAGONAL WEIGHTS

C***  DERIV. OF S WITHOUT THOMSON TERMS
      DSDETA = 0.
      DSDOPA = 0.
      WHERE (OPAKNOTH .GT. 0.)
         DSDETA = 1. / OPAKNOTH
         DSDOPA = -ETAKNOTH / (OPAKNOTH * OPAKNOTH)
      END WHERE

      XLAMKCM = XLAMK * 1.E-8
      WAVK = 1. / XLAMKCM
      WAVK3 = WAVK * WAVK * WAVK  
      DO INDACT=1, MAXFEACT
         IND = INDFEACT(INDACT)
         SIGMA = SIGMAACT(INDACT)

         NUP = IFENUP(IND)
         LOW = IFELOW(IND)
         IF (LOW .EQ. NUP) CYCLE

         WLU = WEIGHT(LOW)/WEIGHT(NUP)
         WAV0 = ELEVEL(NUP) - ELEVEL(LOW)
         XLAM0CM = 1. / WAV0

         IF (SIGMA <= 0.) CYCLE
         
         DO L=1, ND
           POPNUP = POPNUM(L,NUP)
           POPLOW = POPNUM(L,LOW)

C***       Integration of XJFEMEAN (now with superlvl freq. weighting)
           XJFEMEAN(L,IND) = XJFEMEAN(L,IND)
     >                         + WAV0/WAVK * XJL(L) * FWEIGHTL * SIGMA
            
C***       CORRECTION FOR ENERGY EQUATION
            FTFE(L,IND) = FTFE(L,IND)
     >            + (ETAFEI(L,IND) - OPAFEI(L,IND)*XJL(L)/RADIUS2(L))
     >              * FWEIGHTL

C***        Radiative rates of iron superlines
            XJLPLAIN = XJL(L) / RADIUS2(L)

C***        optional de-activation of exp term - see comment in CMFFEOP
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
               WLUEXP = WLU 
            ELSEIF (IVERS_FE_EXPFAC .EQ. 1) THEN
               WLUEXP = WLU * EXP(C1*(WAV0-WAVK)/MAX(T(L),TEFF))
            ELSEIF (IVERS_FE_EXPFAC .EQ. 2) THEN
               WLUEXP = WLU * EXP(C1*(WAV0-WAVK)/T(L))
            ELSE
              STOP '*** FATAL INTERNAL ERROR 1 in subr. FREQUINT'
            ENDIF

C*** wrh version XLAMKCM --> XLAM0CM in the next 2 eqns. (DRATLU, DRATUL)
            DRATLU =  
     >          PI8/C2 * XJLPLAIN * XLAMKCM * SIGMA * FWEIGHTL 
            DRATUL =  
     >          WLUEXP * PI8/C2 * 
     >          (C2 * WAVK3 + XJLPLAIN) * XLAMKCM * SIGMA * FWEIGHTL
     
C***        Prevent using rates where DEPARTURE coefficients are wrong due to POPMIN            
c            IF (POPLOW < POPMIN) THEN
c              DRATLU = 0.
c              DRATUL = 0.
c            ELSEIF (POPNUP < POPMIN) THEN
C***          Do not add any contributions from levels which could be switched off
c              DRATUL = 0.
c            ENDIF
            
            FERATLU(IND,L) = FERATLU(IND,L) + DRATLU
            FERATUL(IND,L) = FERATUL(IND,L) + DRATUL

C***       CALCULATION OF ACTUAL DERIVATIVES OF S

C***       DERIV. OF ETAI
            IF ((OPAFEI(L,IND) .GT. 0.) .AND. (POPNUP > 0.)) THEN
               DETAIDNUP = ETAFEI(L,IND) / POPNUP
            ELSE
               DETAIDNUP = 0.
            ENDIF

C***       DERIV. OF OPAI
            IF ((OPAFEI(L,IND) .GT. 0.)
     >                 .AND. (POPLOW > 0. .OR. POPNUP > 0.)) THEN
               OPANULL  = OPAFEI(L,IND) / (POPLOW - WLUEXP*POPNUP)
            ELSE
               OPANULL = 0.
            ENDIF
            DOPAIDLOW = OPANULL
            DOPAIDNUP = -OPANULL*WLUEXP

C***       DERIV. OF S WITH RESPECT TO LOW AND NUP
            DSDFELOW(L,IND) = DSDOPA(L)*DOPAIDLOW
            DSDFENUP(L,IND) = DSDOPA(L)*DOPAIDNUP + DSDETA(L)*DETAIDNUP

C***       ADD DERIV. AT THE NEW FREQUENCY TO DERIV. OF J
            DJDLOW = DJDS(L) * DSDFELOW(L,IND)
            DJDNUP = DJDS(L) * DSDFENUP(L,IND)

C***       INTEGRATION OF DIAGONAL WEIGHTS
            WFELOW(L,IND) = WFELOW(L, IND) + DJDLOW * FWEIGHTL * SIGMA
            WFENUP(L,IND) = WFENUP(L, IND) + DJDNUP * FWEIGHTL * SIGMA

         ENDDO
      ENDDO


C********************************************************************
C***  Averaging the intensity from fine to coarse frequency grid
C********************************************************************

C***  Note: The current (fine) frequency point bears the index K; the 
C***        integration step trails and hits the frequency index K-1 = OLD
C***  KOLD defines the current coarse interval (KCONT, KCONT+1)

C***  Check whether the current frequency point is the first which 
C***  has entered a new coarse interval
      IF (K .EQ. 0) THEN 
         BNEWINTERVALOLD = .FALSE.
         KCONT = 1
         BEMPTY = .FALSE.
      ELSE
         BNEWINTERVALOLD =  BNEWINTERVAL
      ENDIF
      BNEWINTERVAL = XLAMK .GE. XLAMBDA(KCONT+1) 

C***  Frequency Integration Weight for old finegrid point
C***  Integration weight nue*dnue for test function 
C***   Note: First and last interval get additional terms - see below!
      IF (K .EQ. 0) THEN

         IF (BWCELAB) THEN
            FWTEST = 1.
         ELSE
            FWTEST = 0.
         ENDIF
         DO L=1, ND
            SUMJ     (L) = .0
            SUMJW    (L) = .0
            SUMDJDSC (L) = .0
            SUMDJDSCW(L) = .0
            SUMOPA   (L) = .0
            SUMOPAW  (L) = .0
            SUMETA   (L) = .0
            SUMETAW  (L) = .0
         ENDDO
         SUMEMF = 0.
         SUMEMFW= 0.
         SUMFW  =0.
         SUMFWW =0.
         XNUEK = CLIGHT / XLAMK

         IF (OPC .EQ. 'DIAGMIN') THEN
           DO KK=1, NF
             DO L=1, ND
               WJC_MIN(L,KK) = 1.
             ENDDO
           ENDDO
         ENDIF

         DLEFT  = 0.
         DRIGHT = 0.

cccccccccccccccccccccc LTE-Test Facility cccccccccccccccc
ccc      DO L=1, ND
ccc         xjl(l) = xjlsave(l) 
ccc         xjlold(l) = xjloldsave(l) 
ccc      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         RETURN

      ELSE IF (K .EQ. 1) THEN
         XNUEKOLD  = XNUEK
         XNUEK     = CLIGHT / XLAMK
         XNUE2 = CLIGHT / XLAMBDA(1)
         XNUE3 = CLIGHT / XLAMBDA(2)

         DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
         DRIGHT = 1. - DLEFT
         FWKOLD      = 0.5 * (XNUEKOLD - XNUEK)
         FWNUEKOLD   = (2. * XNUEKOLD + XNUEK) *
     >                 (XNUEKOLD - XNUEK) / 6.        

      ELSE
         XNUEKOLDOLD = XNUEKOLD
         XNUEKOLD  = XNUEK
         XNUEK     = CLIGHT / XLAMK
         IF (BNEWINTERVALOLD) THEN
            XNUELEFT = CLIGHT / XLAMBDA(KCONT)
         ELSE
            XNUELEFT = XNUEKOLDOLD
         ENDIF
         IF (BNEWINTERVAL) THEN
            XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
         ELSE
            XNUERIGHT = XNUEK
         ENDIF

         DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
         DRIGHT = 1. - DLEFT
         FWKOLD      = 0.5 * (XNUELEFT-XNUERIGHT)
         FWNUEKOLD   = (XNUELEFT + XNUEKOLD + XNUERIGHT) *
     >                 (XNUELEFT-XNUERIGHT) / 6.        

      ENDIF

C***  Determine minimum of 1.-EXP(-TAU), if requested
      IF (OPC .EQ. 'DIAGMIN') THEN
        DLEFT_MIN  = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
        DRIGHT_MIN = 1. - DLEFT_MIN
C        IF (DLEFT_MIN .GT. DRIGHT_MIN) THEN
        IF (DLEFT_MIN .LT. DRIGHT_MIN) THEN
          DO L=1, ND
            WJC_MIN(L,KCONT) = AMIN1(WJC_MIN(L,KCONT), DJDS_OLD(L))
          ENDDO
        ELSE
          DO L=1, ND
            WJC_MIN(L,KCONT+1) = AMIN1(WJC_MIN(L,KCONT+1), DJDS_OLD(L))
          ENDDO
        ENDIF
c        write (0,'(a,4(1x,e12.5), 2(1x,f7.3), 3(1x,e12.5))')
c     >    '===', XLAMKOLD, xnuekold, xnue2, xnue3, dleft, dright, 
c     >    WJC_MIN(50,KCONT), WJC_MIN(50,KCONT+1), DJDS_OLD(50)
      ENDIF

C***  Here is the entry for an empty interval:
   10 CONTINUE

      IF (KCONT .LT. 1) THEN
         WRITE (0,*) '== WARNING from FREQUINT: KCONT .LT. 1'
         KCONT = 1
      ENDIF
      IF (KCONT .GT. NF-1) THEN
         WRITE (0,*) '== WARNING from FREQUINT: KCONT .GT. NF-1'
         KCONT = NF-1
      ENDIF

      IF (.NOT. BEMPTY) THEN
         DO L=1, ND
C***        Average the Jnue's from the fine frequency step for the coarse
C***        continuum frequency mesh used by STEAL
            SUMJ (L)  = SUMJ (L) + XJLOLD(L) * FWKOLD
C***        The same is done after multiplication with the testfunction NUE
            SUMJW(L)  = SUMJW(L) + XJLOLD(L) * FWNUEKOLD    
            !the line above can sometimes lead to a (floating invalid) COLI crash (why?)

C***        Average opacity and emissivity for the coarse grid            
            SUMOPA(L) = SUMOPA(L) + OPAKOLD(L) * FWKOLD
            SUMOPAW(L) = SUMOPAW(L) + OPAKOLD(L) * FWNUEKOLD
            SUMETA(L) = SUMETA(L) + ETAKOLD(L) * FWKOLD
            SUMETAW(L) = SUMETAW(L) + ETAKOLD(L) * FWNUEKOLD
            
C***        Old Dioagonal Weight for Continuum is saved
            DJDSCOLD(L) = DJDSC(L)

C***        Dioagonal Weight for Continuum is calculated:
            DJDSC(L) = 0.

C***        Deriv. of J with resp. to Sc at the current Frequency
            IF (OPAKNOTH(L) .GT. 0.) THEN
              OPAC  = OPA(L)*(1. - THOMSON(L))
              DSDSC(L) = OPAC/OPAKNOTH(L)
            ELSE
              DSDSC(L) = 0.
            ENDIF

            IF (OPC .EQ. 'DIAGTAU') THEN
C***          In this version DIAGTAU, the Scharmer-weight factor DJDSC is first 
C***            convertet into a kinf of TAU before averaged over the fine grid
C***            The back conversion is made in SUBR. FREQUNORM
              DJFO = DJDS(L)*DSDSC(L)
C!!!              Operator without Opacity Correction
C!!!              DJFO = DJDS(L)

C***          Special WCHARM plot
              IF (L .EQ. LPLOT_WCHARM) THEN
                WRITE (105, '(2(E20.10,1X))') XLAMK, DJFO
              ENDIF
              IF (DJFO .GE. 1.) THEN
                DJDSC(L) = 0.
              ELSE
                DJDSC(L) = -ALOG(1. - DJFO)
              ENDIF

C***        In this version the averaging is made directly with DJDSC
C***        WARNING: DERIVATIVE TERM UNCLEAR !!!
            ELSE IF (OPC .EQ. 'DIAG') THEN
C***          Scale DJDSC with ratio of Source funktions
              DJDSC(L) = DJDS(L)*DSDSC(L)
C!              IF (DJDSC(L) .LT. 0.) DJDSC(L)=0.

            ELSE IF (OPC .EQ. 'DIAGNC') THEN
C***          B) Add Deriv. at the actual Frequency
              DJDSC(L) = DJDSC(L) + DJDS(L)*DSDSC(L)
              IF (DJDSC(L) .LT. 0.) DJDSC(L)=0.
            ELSE IF (OPC .EQ. 'DIAGMIN') THEN
              DJFO = DJDS(L)
C***          Special WCHARM plot
              IF (L .EQ. LPLOT_WCHARM) THEN
                WRITE (105, '(2(E20.10,1X))') XLAMK, DJFO
              ENDIF
            ELSE IF (OPC .EQ. 'NONE') THEN
              DJDSC(L) = 0.
            ELSE
              WRITE (0,*) 'Unknown Version OPC = ', OPC
            ENDIF

C***        DJDSCOLD is the diagonal element for the continuum 
C***        for optional use in STEAL
            IF (BWCELAB) THEN
C***        Average DJDSCOLD from the fine frequency step for the coarse
C***        continuum frequency mesh used by STEAL
               SUMDJDSC (L)  = SUMDJDSC (L) + DJDSCOLD(L) * FWKOLD
C***        Same is done after multiplication with the testfunction NUE
               SUMDJDSCW(L)  = SUMDJDSCW(L) + DJDSCOLD(L) * FWNUEKOLD
            ELSE
C***        Integration of WJC ('non elaborated', trapezregel)
               WJC(L, KCONT)  = WJC(L, KCONT)
     >                         + DRIGHT * DJDSCOLD(L) * FWKOLD
               WJC(L, KCONT+1)= WJC(L, KCONT+1)
     >                         + DLEFT  * DJDSCOLD(L) * FWKOLD
            ENDIF

         ENDDO
C***     Integration of Emergent Flux
         SUMEMF  = SUMEMF + XHLOLD(1) * FWKOLD
         SUMEMFW = SUMEMFW+ XHLOLD(1) * FWNUEKOLD
         SUMFW   = SUMFW  + FWKOLD
         SUMFWW  = SUMFWW + FWNUEKOLD
C***     Frequency Weight for Normalization
         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT)  = FWTEST(KCONT)  + DRIGHT * FWKOLD
            FWTEST(KCONT+1)= FWTEST(KCONT+1)+ DLEFT  * FWKOLD
         ENDIF

      ENDIF

      IF (BNEWINTERVAL) THEN 
C***     Now the (half-)interval (KCONT, KCONT+1) is completed:
C***     XLAMKOLD to XLAMBDA(KCONT+1), i.e. XNUEKOLD ... XNUE2
C***     and subtract (!) the excess part of the last integration step
         XNUE1 = CLIGHT / XLAMBDA(KCONT  )
         XNUE2 = CLIGHT / XLAMBDA(KCONT+1)
         XNUE3 = CLIGHT / XLAMBDA(KCONT+2)

         QSPECIAL = (XNUEKOLD - XNUE2) / (XNUEKOLD - XNUEK)
         QSPECIAL2 = 1. - QSPECIAL

         IF (BEMPTY) THEN
            XNUELEFT = CLIGHT / XLAMBDA(KCONT)
         ELSE
            XNUELEFT = XNUEKOLD
         ENDIF

C***     Complete the integration with the part of the step from
C***     MAX(XLAMKOLD,XLAMBDA(KCONT) ... XLAMBDA(KCONT+1),
C***     i.e. XNUELEFT ... XNUE2
         FWSPECIAL    = 0.5 * (XNUELEFT - XNUE2)
         FWNUESPECIAL = (XNUELEFT - XNUE2) * (2.*XNUE2 + XNUELEFT) / 6.

         DO L=1, ND
            XJSPECIAL = QSPECIAL2 * XJLOLD(L) + QSPECIAL * XJL(L)
            SUMJ (L)  = SUMJ (L) + XJSPECIAL * FWSPECIAL
            SUMJW(L)  = SUMJW(L) + XJSPECIAL * FWNUESPECIAL
            DJDSCSPECIAL = QSPECIAL2 * DJDSCOLD(L) 
     >                    + QSPECIAL * DJDSC(L)

            IF (BWCELAB) THEN
               SUMDJDSC (L)  = SUMDJDSC (L) + DJDSCSPECIAL * FWSPECIAL
               SUMDJDSCW(L)  = SUMDJDSCW(L)
     >                         + DJDSCSPECIAL * FWNUESPECIAL
            ELSE
C***        Integration with DLEFT=0., DRIGHT=1.
               WJC(L, KCONT+1)= WJC(L, KCONT+1)
     >                          + DJDSCSPECIAL * FWSPECIAL
            ENDIF

            OPASPECIAL = QSPECIAL2 * OPAKOLD(L) + QSPECIAL * OPAK(L)
            SUMOPA (L)  = SUMOPA (L) + OPASPECIAL * FWSPECIAL
            SUMOPAW(L)  = SUMOPAW(L) + OPASPECIAL * FWNUESPECIAL
            
            ETASPECIAL = QSPECIAL2 * ETAKOLD(L) + QSPECIAL * ETAK(L)
            SUMETA (L)  = SUMETA (L) + ETASPECIAL * FWSPECIAL
            SUMETAW(L)  = SUMETAW(L) + ETASPECIAL * FWNUESPECIAL
            
         ENDDO

         XHSPECIAL = QSPECIAL2 * XHLOLD(1) + QSPECIAL * XHL(1)
         SUMEMF    = SUMEMF + XHSPECIAL * FWSPECIAL
         SUMEMFW   = SUMEMFW+ XHSPECIAL * FWNUESPECIAL
         SUMFW     = SUMFW  + FWSPECIAL
         SUMFWW    = SUMFWW + FWNUESPECIAL

         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT+1) = FWTEST(KCONT+1) + FWSPECIAL
         ENDIF

C***     Generate intergration weights at coarse frequency points
C***     which incorporate the test function, i.e. nue*dnue
         FWNUE1 = (2.*XNUE1 + XNUE2) * (XNUE1-XNUE2) / 6.
         FWNUE2 = (2.*XNUE2 + XNUE1) * (XNUE1-XNUE2) / 6.

         FW1 = 0.5 * (XNUE1 - XNUE2)         
         FW2 = FW1
         IF (KCONT == 1) THEN 
            FP = 1.
            DO L=1, ND
               XJCINT(L,1) = 0.
               OPAKINT(L,1) = 0.
               ETAKINT(L,1) = 0.
            ENDDO 
            EMCOLI(1) = 0.
         ELSE
            XNUE0 = CLIGHT / XLAMBDA(KCONT-1) 
            FP = (XNUE1 - XNUE2) / (XNUE0 - XNUE2) 
         ENDIF
         FQ = 1. - FP

         DO L=1, ND
C***       Note: At the left point, J has already the contribution 
C***                from the right end of the last interval
C***       Note: The contributions are now weighted 
            XJCINT(L,KCONT+1) =
     >               ( FWNUE1 * SUMJ(L) - FW1 * SUMJW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
            XJCINT(L,KCONT) = FQ * XJCINT(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMJ(L) - FW2 * SUMJW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)

C***       **** ELABORATED Integration of WJC ******************
            IF (BWCELAB) THEN
               WJC(L,KCONT+1) =
     >               ( FWNUE1 * SUMDJDSC(L) - FW1 * SUMDJDSCW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
               WJC(L,KCONT) = FQ * WJC(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMDJDSC(L) - FW2 * SUMDJDSCW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)
C***       ******************************************************
            ENDIF

            OPAKINT(L,KONT+1) = 
     >               ( FWNUE1 * SUMOPA(L) - FW1 * SUMOPAW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
            OPAKINT(L,KCONT) = FQ * OPAKINT(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMOPA(L) - FW2 * SUMOPAW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)

            ETAKINT(L,KONT+1) = 
     >               ( FWNUE1 * SUMETA(L) - FW1 * SUMETAW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
            ETAKINT(L,KCONT) = FQ * ETAKINT(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMETA(L) - FW2 * SUMETAW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)
     
         ENDDO
         
         EMCOLI(KCONT+1) = 4.*
     >               ( FWNUE1 * SUMEMF - FW1 * SUMEMFW ) /
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
         EMCOLI(KCONT)   = FQ * EMCOLI(KCONT) + 4.*FP *
     >               ( FWNUE2 * SUMEMF - FW2 * SUMEMFW)/
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)

C***     Now increment KCONT and check for empty coarse interval 
C***     (i.e. which contains no fine frequency point).
C***     In that case the fine Interval ends at XLAMBDA(KCONT+1)
         KCONT = KCONT + 1
         IF (XLAMK .GT. XLAMBDA(KCONT+1)) THEN
            BEMPTY = .TRUE.
            XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
         ELSE
            BEMPTY = .FALSE.
            XNUERIGHT = XNUEK
         ENDIF

C***     Start the next integration (KCONT, KCONT+1) with the old
C***     (half-)interval:
C***     XLAMBDA(KCONT) to MIN(XLAMK,XLAMBDA(KCONT+1),
C***     i.e. XNUE2 ... XNUERIGHT
         FWSPECIAL    = 0.5 * (XNUE2 - XNUERIGHT)
         FWNUESPECIAL = (XNUE2 - XNUERIGHT)
     >                  * (2.*XNUE2 + XNUERIGHT) / 6.

         DO L=1, ND
            XJSPECIAL = QSPECIAL2 * XJLOLD(L) + QSPECIAL * XJL(L)
            SUMJ (L)  = XJSPECIAL * FWSPECIAL
            SUMJW(L)  = XJSPECIAL * FWNUESPECIAL
            DJDSCSPECIAL = QSPECIAL2 * DJDSCOLD(L) 
     >                    + QSPECIAL * DJDSC(L)

            IF (BWCELAB) THEN
               SUMDJDSC (L)  = DJDSCSPECIAL * FWSPECIAL
               SUMDJDSCW(L)  = DJDSCSPECIAL * FWNUESPECIAL
            ELSE
C***           Integration with DLEFT=1., DRIGHT=0.
               WJC(L, KCONT) = WJC(L, KCONT) + DJDSCSPECIAL * FWSPECIAL
            ENDIF

         ENDDO
         XHSPECIAL = QSPECIAL2 * XHLOLD(1) + QSPECIAL * XHL(1)
         SUMEMF  = XHSPECIAL * FWSPECIAL
         SUMEMFW = XHSPECIAL * FWNUESPECIAL
         SUMFW   = FWSPECIAL
         SUMFWW  = FWNUESPECIAL
C***     Normalization
         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT) = FWTEST(KCONT) + FWSPECIAL
         ENDIF

C***     In case of empty Interval go again to 'New Interval':
         IF (BEMPTY) GOTO 10

      ENDIF


C**********************************************************************
C***  INTEGRATE TOTAL FLUX AND RADIATIVE FORCE
C**********************************************************************
C      IF (K .EQ. 0) THEN
      IF (K == 1) THEN
         FWMID = 0.5 * (XNUEKOLD - XNUEK)
      ELSE
         FWMID = 0.5 * (XNUEKOLDOLD - XNUEK)
      ENDIF
      IF (FWMID < 0) THEN
        WRITE (0,*) "K= ",K," FWMID = ", FWMID
      ENDIF
      DND = 0.5 * (XNUEKOLD + XNUEK) * VDOPUNIT / CLIGHTKM
      
cc      WRITE (0,*) 'FW(K)', K, FWEIGHTL, FWMID
C***  TEST: Use COLI FWEIGHTL instead of FWMID
      FWMID = FWEIGHTL

      DO L=1, ND
         IF (XJLOLD(L) .GT. 1.E-15) THEN
            XJTOTL(L) = XJTOTL(L) + XJLOLD(L)*FWMID
            bIgnoreK = .FALSE.
         ELSE
            bIgnoreK = .TRUE.
            iIgnoreK(L) = iIgnoreK(L) + 1
         ENDIF
         IF (BCOLIRAY) THEN
              XKTOTL(L) = XKTOTL(L) + XKLOLD(L)*FWMID
         ELSEIF (XJLOLD(L) .GT. 1.E-15) THEN
            !TEST: use only from 1.E-15
              XKTOTL(L) = XKTOTL(L) + XJLOLD(L)*EDDIFO(L)*FWMID
         ENDIF
         
         IF (L == 1) THEN
           HTOTOUTMINUS = HTOTOUTMINUS 
     >           + (EDDIHOUTOLD - EDDIHOUTOP) * XJLOLD(L) * FWMID
         ENDIF
         
         
         IF (L == ND) THEN
           IF (K > 0) THEN
C***         Calculate Eddington Flux at inner boundary
C***          via frequency integration of
C***         HTOTND: H_nu,ND = H_nu,diff + B_nu/2 - h_nu,in * J_nu,ND
C***         Also: Special H- for inner boundary, maybe from Kudritzki (1973, PhD Thesis)    
             HTOTMINUSND = HTOTMINUSND + EDDIHINMOLD*XJLOLD(L) * FWMID
             HTOTNDS = HTOTNDS + EDDIHINTOLD*XJLOLD(L) * FWMID
             HTOTND = HTOTND + ( XHID  + 
     >               BNUE(XLAMKOLD,T(ND))/2.- EDDIHINOLD*XJLOLD(L)) * FWMID
             HTOTNDCOR = HTOTNDCOR + 
     >            (BNUE(XLAMKOLD,T(ND))/2.- EDDIHINOLD*XJLOLD(L)) * FWMID
           ENDIF
           CYCLE
         ENDIF

C***     Option HTOTCUT = tkae only positive part of XHLOLD(L)
         IF (iHTOTCUT > 0) THEN
           XHLUSE = MAX(XHLOLD(L), 0.)
         ELSE
           XHLUSE = XHLOLD(L)
         ENDIF

C***    Radiation Pressure (total and per element)
         ARAD(L) = ARAD(L)
     >               + 0.5*(OPAKOLD(L)+OPAKOLD(L+1))*XHLUSE*FWMID

C+C***    Contributions of the Continuum
         OPTH = 0.5 * (OPAO(L)*THOMSONO(L) + OPAO(L+1)*THOMSONO(L+1))
         OPAC = 0.5 * (OPAO(L) + OPAO(L+1))
         DO NA=1, NATOM
           ARADELEM(NA,L) = ARADELEM(NA,L) + 0.5 *
     >        (OPAKOLDELEM(NA,L)+OPAKOLDELEM(NA,L+1))*XHLUSE*FWMID      
           ACONTELEM(NA,L) = ACONTELEM(NA,L) + 0.5 *
     >        (OPACOLDELEM(NA,L)+OPACOLDELEM(NA,L+1))*XHLUSE*FWMID      
           DO ION=1, MAXION
             OPAFULLION = 
     >          0.5 * ( OPAKOLDION(L,NA,ION) + OPAKOLDION(L+1,NA,ION) )
             ARADION(L,NA,ION) = ARADION(L,NA,ION) 
     >                                + OPAFULLION * XHLUSE * FWMID
             OPACONTION =
     >          0.5 * ( OPACOLDION(L,NA,ION) + OPACOLDION(L+1,NA,ION) )
             ACONTION(L,NA,ION) = ACONTION(L,NA,ION)
     >                                + OPACONTION * XHLUSE * FWMID        
           ENDDO
         ENDDO
         ACONT(L) = ACONT(L) + OPAC*XHLUSE*FWMID
         ATHOM(L) = ATHOM(L) + OPTH*XHLUSE*FWMID

         IF (iHTOTCUT <= 1) THEN
           XHLUSE = XHLOLD(L)
         ENDIF

         HTOTL(L) = HTOTL(L) + XHLUSE*FWMID
         HTOTLTEST(L) = HTOTLTEST(L) + XHLOLDTEST(L)*FWMID
         IF (BCOLIRAY) THEN
              XNTOTL(L) = XNTOTL(L) + XNLOLD(L)*FWMID
         ELSE
C***       fixed calculation bug with eddimix: add J-contrib for eps > 0 
C***       To date, the stored NTOT in model is not used further in the code -- ansander, Aug 2023
              XJLUSE = 0.5 * (XJLOLD(L) + XJLOLD(L+1))
              XNTOTL(L) = XNTOTL(L) + XHLUSE*EDDIGO(L)*FWMID
     >                     + EDDIGO(L)*EPSGO(L)*XJLUSE*FWMID 
         ENDIF

      END DO

C***  Calculation of OPAROSS at L=ND
      DBDT   = DBNUEDT(XLAMKOLD, T(ND))
      DBDTINT    = DBDTINT    + DBDT*FWMID
      DBDTOPAINT = DBDTOPAINT + DBDT*FWMID/OPAKOLD(ND)
C***  ND-1/2:
      TMID   = 0.5 * (T(ND) + T(ND-1))
      OPAKMID = 0.5 * (OPAKOLD(ND) + OPAKOLD(ND-1))
      DBDT_M   = DBNUEDT(XLAMKOLD, TMID)
      DBDTINT_M    = DBDTINT    + DBDT*FWMID
      DBDTOPAINT_M = DBDTOPAINT + DBDT*FWMID/OPAKMID
C- ***********************


C***  Unsoeld-Lucy Procedure for COLI
C***  Note: Integrands are only incremented here for the OLD frequency point!
C***        Final calculation is performed in SUBR. FREQUNORM

C***  NOTE: cut off long wavelength -- de-activated 22-Apr-2002 17:00:37 wrh
ccc      IF (XLAMKOLD .GT. 20000.) GOTO 70

C***  Cut-Off of lauge optical depths governed by UNLU_TAUMAX und UNLU_TAUMAX2
C***     is switched off if ZERO ( = infinity)
      BCUTOFF  = UNLU_TAUMAX  .NE. 0. 
      BCUTOFF2 = UNLU_TAUMAX2 .NE. 0. 

C***  Special Eddi for outer boundary
      EDDIHOUTJMEAN  = EDDIHOUTJMEAN + EDDIHOUTOLD * XJLOLD(1) * FWMID 


C***  Begin of Loop ******************************************
      DO L=1, ND

C***   Rosseland Opacity
       DBDT       = DBNUEDT(XLAMKOLD, T(L))
       OPAROSS(L) = OPAROSS(L) + DBDT * FWMID / OPAKOLD(L)
       OPAROSSCONT(L) = OPAROSSCONT(L) + DBDT * FWMID / OPAO(L)
       DO NA=1, NATOM
         IF (OPAKOLDELEM(NA,L) /= 0.) THEN
           OPAROSSELEM(NA, L) = OPAROSSELEM(NA, L) + DBDT * FWMID / 
     >                                OPAKOLDELEM(NA, L)      
         ENDIF
         OPATOTELEM(NA,L) = OPATOTELEM(NA,L) + OPAKOLDELEM(NA,L)*FWMID
       ENDDO
       
C***   simple frequency-integrated opacity       
       OPACTOT(L) = OPACTOT(L) + OPAO(L) * FWMID
       OPALTOT(L) = OPALTOT(L) + OPAKOLD(L) * FWMID

C***   Determine current pure line opacity       
       OPALOL = OPAKOLD(L) - OPAFE(L) - OPAO(L)
       
       IF (BCUTOFF2 .AND.
     >    OPAKNOTHO(L) * (RADIUS(L)-1.) .GT. UNLU_TAUMAX2) GOTO 45
       IF (XJLOLD(L) .LT. 1.E-50) GOTO 45
       FTCOLI(L) = FTCOLI(L)
     >  + (ETAKNOTHO(L) - OPAKNOTHO(L)*XJLOLD(L)/RADIUS2(L)) * FWMID
       FTCOLIB(L) = FTCOLIB(L) + OPAKNOTHO(L)*
     >  ( XJLOLD(L)/RADIUS2(L) - BNUE(XLAMKOLD,T(L)) )*FWMID
 45    continue


C***    For temperature corrections: Cut strong lines from flux integral:
        TCUT = MAX(TEFF, T(L))
c        IF (XHLOLD(L) > 0.5 * XJLOLD(L)) THEN
        IF (T(L) < TEFF) THEN
          XHLCUT = MIN(XHLOLD(L), 0.5 * BNUE(XLAMKOLD, TCUT))
        ELSE 
          XHLCUT = XHLOLD(L)
        ENDIF
        HTOTCUT(L) = HTOTCUT(L) + XHLCUT * FWMID

 
 
c       bOSKIP = (BCUTOFF .AND.
c     >    (OPAKNOTHO(L) * (RADIUS(L)-1.) .GT. UNLU_TAUMAX))
       bOSKIP = (BCUTOFF .AND.
     >    (OPALOL * (RADIUS(L)-1.) .GT. UNLU_TAUMAX))
       IF (bOSKIP) XLAMLASTOSKIP(L) = XLAMK
 
C***  CUT OFF VERY LARGE OPTICAL DEPTHS ("detailed balance") 
c       IF (BCUTOFF .AND.
c     >    (OPAKNOTHO(L) * (RADIUS(L)-1.) .GT. UNLU_TAUMAX)) GOTO 44
       IF (BCUTOFF .AND.
     >    (OPALOL * (RADIUS(L)-1.) .GT. UNLU_TAUMAX)) GOTO 44
C***  Special treatment for negative 'true' Opacity
C***   (Ist dieser Fallback schlau oder eher problematisch?)
       bKSKIPMEAN = .FALSE.
       IF (OPAKNOTHO(L) <= 1.E-30) THEN
         bKSKIPMEAN = .TRUE.
       ELSE
         SMEAN(L)    = SMEAN(L)    + ETAKOLD(L)/OPAKOLD(L)   * FWMID
         IF (ABS(ETAKOLD(L) > 0.)) THEN
C***       CUTOPAMEANTHRES is set via CARDS option OPAMEANTHRES
C***        recommended values between 1.E-5 and 1.E-3
C***        Setting too high values can lead to a crash in STEAL
C***        due to OPASMEAN being zero at one or more depth points
           IF (ABS(OPAKOLD(L)/ETAKOLD(L)*XJLOLD(L) - 1.)
     >           < CUTOPAMEANTHRES) 
     >       bKSKIPMEAN = .TRUE.             
         ENDIF
       ENDIF

       IF (.NOT. bKSKIPMEAN) THEN
         OPAJMEAN(L) = OPAJMEAN(L) + OPAKOLD(L)*XJLOLD(L)    * FWMID
         OPASMEAN(L) = OPASMEAN(L) + ETAKOLD(L)              * FWMID
         OPAPMEAN(L) = OPAPMEAN(L) 
     >                      + OPAKOLD(L)*BNUE(XLAMKOLD,T(L)) * FWMID
       ENDIF

ccc OLD BRANCH:
c       IF (bKSKIPMEAN) THEN
c          SMEAN(L)    = SMEAN(L)    + ETAKOLD(L)/OPAKOLD(L)     * FWMID
c          OPAJMEAN(L) = OPAJMEAN(L) + OPAKOLD(L)*XJLOLD(L)      * FWMID
c          OPASMEAN(L) = OPASMEAN(L) + ETAKOLD(L)                * FWMID
c          OPAPMEAN(L) = OPAPMEAN(L) 
c      >                        + OPAKOLD(L)*BNUE(XLAMKOLD,T(L)) * FWMID
c       ELSE
c          SMEAN(L)    = SMEAN(L)    + ETAKNOTHO(L)/OPAKNOTHO(L) * FWMID
c          OPAJMEAN(L) = OPAJMEAN(L) + OPAKNOTHO(L)*XJLOLD(L)    * FWMID
c          OPASMEAN(L) = OPASMEAN(L) + ETAKNOTHO(L)              * FWMID
c          OPAPMEAN(L) = OPAPMEAN(L) 
c     >                       + OPAKNOTHO(L)*BNUE(XLAMKOLD,T(L)) * FWMID
c       ENDIF

C***    Skipping huge lines test       
        IF (.NOT. bOSKIP .AND. XLAMLASTOSKIP(L) > 0.) THEN
          bOSKIP = (XLAMK < XLAMLASTOSKIP(L)
     >                      * (1.+2.*VELO(1)*VDOPUNIT/CLIGHTKM))
        ENDIF

       
        
C***    Calculation of flux intgegral without strong lines 
C***    (can be used in STEAL->TEMPCORR)
        IF (.NOT. bOSKIP) THEN
          IF (bOSKIPLAST(L)) THEN
C***        The last point was skipped, but now we are back to normal
C***         => account for missing integral part by linear
C***             interpolation over the parts where lines have been cut out
c            HTOTCUT(L) = HTOTCUT(L) 
c     >           + 0.5 * (XHLOLD(L) + XHBFSKIP(L)) * (FWSKIP + FWMID)
          ELSE
C***        Neither the current nor the last point was skipped          
c            HTOTCUT(L) = HTOTCUT(L) + XHLOLD(L) * FWMID
          ENDIF
C***      Reset all variables that store stuff for future skips:          
          XHBFSKIP(L) = XHLOLD(L)
          FWSKIP = 0.
        ELSE
C***      The current point        
C***      Add up frequency weights of all skipped points:        
          FWSKIP = FWSKIP + FWMID
        ENDIF               
C***    Store current "skip" status for next frequency point        
        bOSKIPLAST(L) = bOSKIP
       
C***    Account for PRINT DELTAH (debug/demonstration) option:
        IF (L > 0 .AND. L < ND .AND. L == LPRDH .AND. .NOT. bOSKIP) THEN
           VDR = VELO(L+1) / RADIUS(L+1)
           DELTAHNUE = 0.25 * BNUE(XLAMKOLD, TEFF) 
           DHMECH = 0.
           IF (L < ND-1) THEN
             DHMECH =  - ((GRADI(L+1)-VDR)*XKLOLD(L+1)+VDR*XJLOLD(L+1)) 
     >                    * VDOPUNIT / CLIGHTKM      
     >                    *  0.5 * (RADIUS(L)-RADIUS(L+2))
             DELTAHNUE = DELTAHNUE + DHMECH
           ENDIF
           DELTAHNUE = DELTAHNUE - XHLOLD(L)         
           
           WRITE (140, '(G20.10,9(2X,G20.10))')
     >                            XLAMKOLD, DELTAHNUE, DHMECH, XNUEK,
     >                            XJLOLD(L), XHLOLD(L), 
     >                      OPAKOLD(L), OPAO(L), 
     >                      BNUE(XLAMKOLD, TEFF), BNUE(XLAMKOLD, T(L))
        ENDIF
       
       
C***    Option for ALI-like amplification of the temperature correction 
C***      of the "first term" in the Unsoeld-Lucy procedure
C***      De-activated if GAMMAT=0.
        IF (GAMMAT .NE. .0) THEN

C***  A more sophisticated estimate of TAU = Delta Tau to nearest border:
C***  Assume that the opacity dilutes with r**(-2). Integration yields:
            DROUT = RADIUS(L) * ( 1. - RADIUS(L) / RADIUS(1))
            DRIN  = RADIUS(L) * ( RADIUS(L) - 1.)
            DR = AMIN1 ( DROUT, DRIN )

            TAU = OPAKOLD(L) * DR / GAMMAT
            IF (TAU .GT. 1.E-6) THEN
               ALONUE = 1. - (1. - EXP(-TAU)) / TAU
            ELSE
               ALONUE = 0.
            ENDIF
ccc  TEST!!!
c            ALONUE = DJDS(L)
            ALONUE = DJDSraw(L)
C***        Damp with GAMMAT instead of GAMMACOLI            
            ALONUE = -LOG(1.-ALONUE)
            ALONUE = ALONUE / GAMMAT
            ALONUE = 1.-EXP(-ALONUE)

            OPALAMBDAMEAN(L) = OPALAMBDAMEAN(L) + 
     >                        ETAKNOTHO(L) * ALONUE * FWMID
        ENDIF

   44  CONTINUE

C***    The "TC" quantities are for TEMPCORR (Unsoeld-Lucy)
        OPASMEANTC(L) = OPASMEAN(L)
        OPAJMEANTC(L) = OPAJMEAN(L)

        QFJMEAN(L)  = QFJMEAN(L)  + 
     >       QLFOLD(L)*EDDIFO(L)*XJLOLD(L) * FWMID


C***    Flux quantities defined at interstices - skip for L=ND !
        IF (L .EQ. ND) CYCLE

CCC     "EFFECTIVE OPACITY" ACCOUNTING FOR THE INFLUENCE OF  
CCC        THOMSON SCATTERING ON THE THEMALIZATION DEPTH 
C***    ADDITIONALLY, CUT OFF VERY LARGE OPTICAL DEPTHS ("detailed balance") 
        OPATRUE    =  0.5 * (OPAKNOTHO(L) + OPAKNOTHO(L+1)) 
        IF (BCUTOFF) THEN
           IF (OPATRUE * (RADIUS(L)-1.) .GT. UNLU_TAUMAX2) GOTO 46
        ELSE
           IF (XHLOLD(L) .LT. 1.E-50) GOTO 46
        ENDIF

ccc  merkwuerdige Artistik, war wohl auch anders gemeint. 
ccc testweise abgeschaltet wrh 28-Nov-2003 14:35:58
c        OPASCATTER = OPAKHOLD(L) - OPATRUE 
c        OPAOPA     = OPATRUE * OPASCATTER
c        IF (OPAOPA .GT. .0) THEN
c           OPAEFFECTIVE = SQRT (OPAOPA)
c           QOPAHMEAN(L) = QOPAHMEAN(L) + 
c     >          QLHOLD(L) * OPAEFFECTIVE * XHLOLD(L) * FWMID
c        ENDIF
ccc dafuer angeschaltet:

ccc    noch eine Variante: nach Ansicht von Goetz kommt hier: 
ccc    in der 1. Momentengleichung, aus der in TEMPCORR der "2. Term" 
ccc    ausgerechnet wird, die *volle* Opazitaet (incl. Thomson) hin:
ccc        wrh  7-Dec-2005 11:56:43
ccc    Also seine Empfehlung: Ersetze das folgende Statement
cc        QOPAHMEAN(L) = QOPAHMEAN(L) +
cc     >       QLHOLD(L) * OPATRUE * XHLOLD(L) * FWMID
ccc    durch:
        OPAFULL    =  0.5 * (OPAKOLD(L) + OPAKOLD(L+1)) 
        IF (.NOT. bKSKIPMEAN) THEN
          QOPAHMEAN(L) = QOPAHMEAN(L) +
     >       QLHOLD(L) * OPAFULL * XHLOLD(L) * FWMID

C***      HMEAN = HTOTL without optically very think frequences
C***           (same as excluded from mean opacities)
          HMEAN(L)= HMEAN(L) + XHLOLD(L) * FWMID
        ENDIF
        

 46     CONTINUE

      END DO
C***  End of L-Loop

C***  Jump label: integrations skipped for very long wavelengths 
   70 CONTINUE



cccccccccccccccccccccc LTE-Test Facility cccccccccccccccc
ccc      DO L=1, ND
ccc         xjl(l) = xjlsave(l) 
ccc         xjlold(l) = xjloldsave(l) 
ccc      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      RETURN
      END
