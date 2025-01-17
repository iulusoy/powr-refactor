C***  MAIN PROGRAM COLI ********************************************************
      SUBROUTINE COLI

C*******************************************************************************
C***  RADIATION TRANSFER IN THE COMOVING FRAME: CONTINUA AND LINES
C***   FORMAL SOLUTION FROM GIVEN POPULATION NUMBERS
C*******************************************************************************

      IMPLICIT NONE

C***  SET ARRAY DIMENSIONS  ********************************************
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM   =             30 
      INTEGER, PARAMETER :: NDIM      =           2560 
      INTEGER, PARAMETER :: NFDIM     =   2*NDIM + 400 
      INTEGER, PARAMETER :: MAXAUTO   =           3200 
      INTEGER, PARAMETER :: MAXIND    =          45000      !in COLI, MAXIND does NOT contain Fe transitions
      INTEGER, PARAMETER :: MAXKONT   =       NFDIM/2 
      INTEGER, PARAMETER :: NDDIM     =            89 
      INTEGER, PARAMETER :: NPDIM     =            94 
      INTEGER, PARAMETER :: NFLDIM    =           500
C***  FOR ARRAYS WHICH ARE STORED OVER THE WHOLE REDISTRIBUTION WIDTH
      INTEGER, PARAMETER :: NFRI      =          1000 
      INTEGER, PARAMETER :: NFRO      =           100 
      INTEGER, PARAMETER :: MAXHIST   =          4000 
      INTEGER, PARAMETER :: MAXPLOT   =            25 
ccc      INTEGER, PARAMETER :: MAXLIN    =            40 
      INTEGER, PARAMETER :: MAXLIN    =           100 
      INTEGER, PARAMETER :: MAXEXT    =            10 
      INTEGER, PARAMETER :: NDIMEDDIA = 2*NDDIM + 3 + 2 + 4
      INTEGER, PARAMETER :: NREDMAX   =           500 
      INTEGER, PARAMETER :: NIT       =             3 
      INTEGER, PARAMETER :: MAXXDAT =              10 

C***  MAX. NUMBER OF IRON SUPERLINES
      INTEGER, PARAMETER :: MAXFEIND  =           2500 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR 
      INTEGER, PARAMETER :: MAXION = 27 

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
c      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !vd20
      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !custom hydro
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

C***  Fine Integration in STEALCL
C      INTEGER, PARAMETER :: IFF_MAX =   80000      !std
c      INTEGER, PARAMETER :: IFF_MAX =  200000      !vd20
      INTEGER, PARAMETER :: IFF_MAX =  300000      !xxl
      
      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

C***  IFF_MAX_MS must be 1/8 * IFF_MAX !!!
C      INTEGER, PARAMETER :: IFF_MAX_MS =   IFF_MAX / 8 

      INTEGER :: ND, J, N, NF, NP, JOBNUM, NATOM, NAUTO, NCOLIP, NF2,
     >           LASTIND, LASTINDAUTO, LASTKON, LASTFE, NLINE,
     >           LASTINDALL, N_WITH_DRLEVELS, 
     >           LSOPA, LSINT, LEVDEBUG, LASERV, NEXTHYDRO, NEWWRC,
     >           L, NL, IPLOT, LPLOT, IVERS, LPLOT_WCHARM,  
     >           KWORDS, KBLOCKS, ISTATS, IDUMMY, IERR, IFRO, NXJO, NXK,
     >           NDEDDIA, ITMAX, ITSTART, IT, NZE1, NZE2,
     >           IW_COLIMO_F, IW_COLIMO_G, IW_COLIMO_G2, IW_COLIMO_J,
     >           IW_COLIRAY_IPLUS, IW_COLIRAY_U, IWARNJ0, NLASER,
     >           NEWBGC, LBMAX, NBLENDS, KSPACE, KOUT, K, KLAST, NK,
     >           KONCHECK, KONCHECK2, KONTACT, KONTAUP, NLACT, NDK,
     >           LINECHECK, ILINECHECK, KCL, KCU, KCDMAX, KCDELTA,
     >           KCCHECK, IFF_N, IF1, IF2, IF3, IF4, MAXFEACT, KCONT,
     >           ITACT, JP, IRAY, NA, LRTCMAX, IOPAMAX1_K, IWARNJMNEG,
     >           IVERS_FE_EXPFAC, INDDR, N_RUDAUTO, LINEINDEX, ILINE,
     >           iHTOTCUT

      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, IONGRND, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      CHARACTER(10), DIMENSION(MAXAUTO) :: LEVUPAUTO, LEVAUTO

      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      REAL, DIMENSION(MAXATOM) :: STAGE, ATMASS
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, EN
      INTEGER, DIMENSION(MAXKONT) :: KONTHLP
      CHARACTER(8), DIMENSION(MAXKONT) :: IGAUNT, KEYCBF      
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO, SIGMA1I, 
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(4, MAXIND) :: COCO
      
      CHARACTER(MAXHIST*8) :: MODHIST

      REAL, DIMENSION(NFDIM) :: FWEIGHT, HWRCEDDI, XLAMBDA, XLAMBDA2
      REAL, DIMENSION(NDDIM) :: RADIUS, RADIUS2, RADIUSH, RADIUSH2, 
     >                          OPA, ETA, THOMSON, ETANOTH, VMIC,
     >                          VELO, GRADI, ENTOT, T, RNE, PP, W0,
     >                          TA, TB, TC, UB, GA, H, QQ, S, VA, VB, 
     >                          OPAK, ETAK, ETAKNOTH, OPAKNOTH, 
     >                          XHBFSKIP, XLAMLASTOSKIP
      REAL, DIMENSION(NPDIM) :: P
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM
      REAL, DIMENSION(NDDIM, MAXLIN) :: OPAL, ETAL, XJLMEAN, 
     >                                  XLAMAPPMEAN,
     >                                  XLAMAPPUMEAN, XLAMAPPLMEAN
      REAL, DIMENSION(NDDIM, NFDIM) :: XJC
      REAL, DIMENSION(NFLDIM, MAXLIN) :: PHI
      REAL, DIMENSION(NDDIM, MAXLIN) :: PWEIGHT
      REAL, DIMENSION(NDDIM, NPDIM) :: V, U, Z, PPP
      REAL, DIMENSION(NFLDIM) :: XJNUED
      REAL, DIMENSION(NDDIM) :: XJCIND
     
      REAL, DIMENSION(3, NDDIM) :: EDDI
      INTEGER, DIMENSION(MAXIND) :: LINE, INDNUP, INDLOW
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
     
      LOGICAL, DIMENSION(NDIM, NDDIM) :: ZERO_RATES
      LOGICAL, DIMENSION(NDDIM) :: bFELASER, bOSKIPLAST
     
      REAL :: OPAMAX1_LAMBDA, OPAMAX1, HTOTOUTMINUS, EDDIHOUTJMEAN,
     >        DBDTOPAINT_M, DBDTINT_M, DBDTOPAINT, DBDTINT, CMFBAND,
     >        DXMAX, VDOPFE, DXFE, XLAM0FE, TEFF, RSTAR, VDOP, OPARND,
     >        RANGE1, RANGE2, EXLAM1, EXLAM2, DUNLU, DUNLUR, XLP1, XLP2,
     >        EMIXSTART, EMIXFIX, XLAM_FINE_START, XLAM_FINE_END,
     >        GAMMACOLI, GAMMAT, BWES, BWESEX, BWEXBLU, BWEXRED, RIPL2,
     >        UNLU_TAUMAX, UNLU_TAUMAX2, DFEINDR, DTDR, XLAM0, XLAM0LN,
     >        ALN, CMFBANDR, XKCMAX, BANDP, BANDM, DKR, XK, XLAMK, XNK,
     >        XLAMKOLD, DFKONT, EDDIHOUT, EDDIHIN, EDDINOUT, EDDININ,
     >        XHI, XHO, XHOM, XNOM, EDDIHOUTP, EDDINOUTP, BCORE, DELTAX, 
     >        FWEIGHTL, DBDR, XHID, XNO, XNI, XHOP, XNOP, SK1, SK2, 
     >        XNL_MID, XHL_MID, XJL_MID, XNENN, EDDIHOUTOLD, EDDIHINOLD,
     >        XNUEK, XNUEKOLD, XNUEKOLDOLD, SL, SLNOTH, ATMEAN, DELTA,
     >        TOTOUT, POPMIN, HTOTMINUSND, HTOTND, HTOTNDCOR,
     >        EDDIHINM, EDDIHINMOLD, EDDIHINT, EDDIHINTOLD,
     >        HTOTNDS, DELTAR, ARADTOT, ARADTESTTOT,
     >        DTDRIN, VDOPUNIT, CUTOPAMEANTHRES, EMIXMAX
     
CC***  NO MORE CONTINUATION LINES ALLOWED IN BERLIN
 
      !INITIALIZE GAUNTFF common block (just in case)
      INTEGER :: NTUP, NTLOW, NFUP, NFLOW, NHELP
      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP      

      REAL, DIMENSION(NDDIM) :: VMACH, VELORAW

      REAL, DIMENSION(NDDIM) :: WNUE, HTOTL, HTOTLlast, HTOTLTEST, 
     >                          UCONT, VCONT, HTOTCUT
      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM) :: ABXYZ
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(NPDIM) :: WP1, WP1LAST

      CHARACTER(255) :: HISTENTRY
      CHARACTER(120) :: ALPOPT, DRLINES_CARD
      CHARACTER(100) :: MODHEAD
      CHARACTER(10) :: BLENDMSG
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(8), DIMENSION(NDDIM) :: IWARN
      CHARACTER(8) :: NAME, CATMOD
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(2) :: CPPLABEL
      LOGICAL :: LASER, NOLAP, NEWV, DRNORUD, bALOTri,
     >           BLLIST, BFORMAL, BLFERR, bPLOTRTAU1, 
     >           bNoIronLaser, bDDVDOP

C***  ARRAYS USED FOR THE COMPLETE REDISTRIBUTION OF LINE PHOTONS
      REAL, DIMENSION(NDDIM+2, NFRO) :: XJLO
      REAL, DIMENSION(NDDIM+2, NFRI) :: XJLI
      REAL, DIMENSION(NDDIM) :: ETAC
      CHARACTER(1) :: CMODE
      CHARACTER(8) :: LACTION
      CHARACTER(4) :: REDISMODE, LASTREDISMODE, REDISMODED
      LOGICAL :: XJLERR

C***  ARRAYS USED FOR PLOT FACILITIES
      REAL, DIMENSION(NFLDIM) :: XPLOT, YPLOT
      INTEGER, DIMENSION(MAXPLOT) :: LPOPAB, LPOPABD, LPJNUE, LPJNUED,
     >                               LPSNUE, LPSNUED
      REAL, DIMENSION(MAXIND) :: XLAMSOR, XLAMMAX, XLAMMIN
      INTEGER, DIMENSION(MAXIND) :: LINFO
      
C***  Arrays for storing the relative blend complex:
      REAL, DIMENSION(NDDIM, MAXLIN) :: WS
      REAL, DIMENSION(MAXLIN) :: XLAM
      INTEGER, DIMENSION(MAXLIN) :: NUP, LBLEND, LOW, IND

C***  To Store Feautrier Matrices ST(94*95,89)
      REAL, DIMENSION((NPDIM+1)*NPDIM,NDDIM) :: ST
      LOGICAL :: BELIFI

C***  New Dimension for COLI
      REAL, DIMENSION(MAXIND) :: XKMIN, XKMAX, XKMID, XKRED
      REAL, DIMENSION(NFDIM) :: XKC, XKC2
      INTEGER, DIMENSION(MAXLIN) :: LIND, LINDS
      REAL, DIMENSION(NDDIM,NPDIM) :: CWM0, CWM2, CWM1, CWM3
      REAL, DIMENSION(NPDIM) :: CWM1O, CWM1I, CWM3O, CWM3I
      REAL, DIMENSION(NDDIM) :: XJL, XHL, XKL, XNL, XJLMO, XJLMOR2,
     >                          XHLMO, XJTOTL, XKTOTL, XNTOTL,
     >                          EDDIF, DLF, GLF, VLF, GLF2, VLF2, QLF,
     >                          OPAKH, COLIA, COLIB, COLIC, COLIW
      REAL, DIMENSION(NDDIM-1) :: EDDIG, DLH, GLH, VLH, GLH2, VLH2,
     >                            QLH, ALH, BLH, CLH
      REAL, DIMENSION(NDDIM,NIT) :: XJLMO_OLD, XHLMO_OLD, EDDIFO
      REAL, DIMENSION(NDDIM-1,NIT) :: EDDIGO
      REAL, DIMENSION(NDIMEDDIA, NFRO) :: EDDIA
c      dimension exlam1(maxlin), exlam2(maxlin)
      REAL, DIMENSION(0:NREDMAX, NDDIM) :: WREDI
      INTEGER, DIMENSION(NDDIM) :: NREDI, iIgnoreK
      REAL, DIMENSION(NDDIM) :: WRED1, WRED2, WREDISUM, FULFIL0, 
     >                          FULFIL1, OPAKNOTHO, ETAKNOTHO, 
     >                          CLMOETA, CLMOOPA
      REAL, DIMENSION(NDDIM,NPDIM,NIT) :: U_OLD, V_OLD
      REAL, DIMENSION(NIT) :: EDDIHOUTO, EDDIHINO, XHIO,
     >                        EDDINOUTO, EDDININO,
     >                        XHOMO, XNOMO, EDDIHOUTOP, EDDINOUTOP
      CHARACTER(8) :: OPC

C***  Interpolation of continuum opacities
      REAL, DIMENSION(NDDIM) :: OPACL, OPACU, ETACL, ETACU
      REAL, DIMENSION(MAXATOM) :: OPALEL, ETALEL
      REAL, DIMENSION(NDDIM, MAXATOM) :: OPALION, ETALION
      REAL, DIMENSION(MAXATOM, NDDIM) :: OPACLEL, OPACUEL, 
     >                                   ETACLEL, ETACUEL
      REAL, DIMENSION(NDDIM, MAXATOM, MAXION) :: OPACLION, OPACUION, 
     >                                           ETACLION, ETACUION
     
C***  Depth- and element-dependent VDOP     
      REAL, DIMENSION(NDDIM, MAXATOM) :: VDOPDD

C***  Unsoeld-Lucy
      REAL, DIMENSION(NDDIM) :: OPASMEAN, SMEAN, QFJMEAN, OPAJMEAN,
     >                          QOPAHMEAN, HMEAN, 
     >                          OPASMEANTC, OPAPMEAN,
     >                          OPAJMEANTC, OPAROSS, OPALAMBDAMEAN,
     >                          QLFOLD, QLHOLD, OPAKHOLD, 
     >                          TAUROSS, TAULAST

C***  New Dimensions for Short Characteristics
      REAL, DIMENSION(NDDIM,NPDIM) :: XIPLUS, XIMINUS
      REAL, DIMENSION(NDDIM,NPDIM,NIT) :: XIPLUS_OLD, XIMINUS_OLD
      REAL, DIMENSION(NDDIM,NIT) :: S_OLD, OPAK_OLD, EPSG, EPSGO
      REAL, DIMENSION(NDDIM) :: EPSGMAX, GEPSB, GEPSBO
      REAL, DIMENSION(NIT) :: EPSGITMAX
      
      LOGICAL :: BKONCHECK, BKONCHECK2, BPLOT, BPLOT2, BPDONE, BPSTART,
     >           BCOLIP, BCLEERR, BSTATIC, BUSEMO, BCOLIRAY, BSMALLJ,
     >           CLHLP, BXJCE, BSHORT, BITCONT, 
     >           BKUDRITZKI, BCOLIPP, bDEBUG,
     >           BSHORT_CHAR, BEPSGMAXERR, BEMIX, BEMIXFIX,
     >           bWARNK, bMAXEMIX, bNoNEGEDDIG

C***  Array to handle Laser Condition output 
      LOGICAL, DIMENSION(MAXIND) :: BLASERL
      
C***  Variables to handle Fine Integration in Steal
C***  FF_INFO : Information on the Frequency grid
C***         (1) : XLAM0
C***         (2) : ALN
C***         (3) : XLAM_FINE_START
C***         (4) : XLAM_FINE_END
C***               (3) + (4) as demanded by the Input Cards
C***         (5) : Start index of Integration
C***         (6) : End index of Integration
C***         (7) : IFF_N Total number of Fine Frequencies
C***         (8) : XLAM_FINE_START
C***         (9) : XLAM_FINE_END
C***               (8) + (9) wavelength of first and last fine point
C***        (10) : KSPACE
C***  BFF_ACT : True, if actual K is in Fine Interval
C***  IFF_DK stores in a compressed way the Delta-K values
C***    IFF_DK can have values between -100 and 100
      REAL, DIMENSION(10) :: FF_INFO
      LOGICAL :: BFF_ACT
C      INTEGER, DIMENSION(IFF_MAX) :: IFF_DK
C      INTEGER, DIMENSION(IFF_MAX,NDDIM) :: IFF_WCHARM
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX,NDDIM) :: IFF_WCHARM


      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW, INDFEACT
      REAL, DIMENSION(MAXFEIND) :: SIGMAACT, SIGMAINT
      REAL, DIMENSION(MAXFEIND, NDDIM) :: FERATLU, FERATUL   !<-- dies ist im Steal anders dimensioniert

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(NDDIM) :: OPAFE, ETAFE, OPAROSSCONT,
     >                          OPAKFE, OPAKFEOLD, 
     >                          OPAFEFT, ETAFEFT, OPAKFEFT, ETAKFEFT, 
     >                          OPAKFEFTOLD, ETAKFEFTOLD
      REAL, DIMENSION(NDDIM, MAXION) :: OPAFEION, ETAFEION      
C      REAL, DIMENSION(NDDIM, MAXIND) :: WFELOW, WFENUP      <--- ist das vielleicht richtig?
C!!! Dimensionierungsfehler, entdeckt am 16-Sep-2002 14:07:15 goetz/wrh
ccc     >              WFELOW(NDDIM,MAXFEIND), WFENUP(NDDIM,MAXFEIND)
      REAL, DIMENSION(NDDIM, MAXFEIND) :: WFELOW, WFENUP, XJFEMEAN,
     >                                    OPAFEI, ETAFEI, FTFE,
     >                                    DSDFELOW, DSDFENUP

      REAL, DIMENSION(INDEXMAX) :: SIGMAFE

      LOGICAL :: BFECHECK, BFEWING, BFEMODEL, bTDIFFUS

C***  Calculation of Radiative Accelleration
      REAL, DIMENSION(NDDIM,NIT) :: XJLOLD, XHLOLD
      REAL, DIMENSION(NDDIM) :: OPAKOLD, ETAKOLD, 
     >                          OPAKNOFENOTH, OPAKNOFENOTHO, 
     >                          ETAKNOFENOTH, ETAKNOFENOTHO
      REAL, DIMENSION(NDDIM) :: FTCOLI, FTCOLIB, 
     >                          ARAD, ACONT, ATHOM
C***  Calculation of corrected Flux
      REAL, DIMENSION(NDDIM,NIT) :: XKLOLD, XNLOLD
C***  Calculation of XJC's and coarse-grid opacities within COLI
      REAL, DIMENSION(NDDIM,NFDIM) :: XJCINT, OPAKINT, ETAKINT
      REAL, DIMENSION(NFDIM) :: FWTEST
C***  Derivatives of J with respect to S and JOLD
      REAL, DIMENSION(NDDIM) :: DJDSMOD, DJDSMOD_OLD, DJDSMOU, DJDSMOL
C***  For Test Output
      REAL, DIMENSION(NDDIM) :: DJDS, DJDS_OLD, DJDSraw
C***  Diagonal Op. for Continua
      REAL, DIMENSION(NDDIM,NFDIM) :: WJC
      REAL, DIMENSION(NDDIM) :: DSDSC
C***  New for determining the Minimum of 1-EXP(-Tau)
      REAL, DIMENSION(NDDIM,NFDIM) :: WJC_MIN
      REAL, DIMENSION(NDDIM) :: DJDSC,DJDSCOLD,ETANOTHO,OPAO,THOMSONO
C***  Temporary storage vectors for FREQUINT
      REAL, DIMENSION(NDDIM) :: SUMJ, SUMJW, SUMDJDSC, SUMDJDSCW,
     >                          SUMOPA, SUMOPAW, SUMETA, SUMETAW, 
     >                          SC, SCO

C***  Depth dependent Clumping
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C***  Emergent Flux calculated in COLI
      REAL, DIMENSION(NFDIM) :: EMCOLI

C***  Force Multipliers for STEAL->HYDROSOLVE
      LOGICAL :: bForceCOLIP, bForceCLEERR, BPLOTALPHA, bAdjVINF
      LOGICAL :: bKALPHA, bHYDROSOLVE, bKATEST, bKADONE, bNewLoop,
     >           bPrintAlpha
      REAL, DIMENSION(NDDIM) :: HTOTTEST, ARADTEST, ATHOMTEST,
     >                          ALPHAF, ALPHAALT, RHO, XMU, 
     >                          GAMMATHIN, CK, ARMDRESP,
     >                          ACONTTEST, ALINETEST,
     >                          ALINEMDRESP, ACONTMDRESP
      REAL, DIMENSION(MAXATOM, NDDIM) :: OPAKELEM, OPAKOLDELEM, 
     >                                   OPACELEM, OPACOLDELEM, 
     >                                   ETAKELEM, ETAKOLDELEM, 
     >                                   ETACELEM, ETACOLDELEM, 
     >                                   ARADELEM, ACONTELEM,
     >                                   OPAROSSELEM, OPATOTELEM
      REAL, DIMENSION(NDDIM, MAXATOM, MAXION) :: OPAKION, OPAKOLDION,
     >                                           OPACION, OPACOLDION,
     >                                           ETAKION, ETAKOLDION,
     >                                           ETACION, ETACOLDION,
     >                                           ARADION, ACONTION

      INTEGER :: iTypeAlpha, LPRDH
      REAL :: PARALAS, VINF, VELOMODFAK, DOTM4P, 
     >        XMG, XMSTAR, GAMMAPOWR, ACONTTOT, ACONTTESTTOT

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS       

      CHARACTER(10) :: TIM1

      REAL, PARAMETER :: WPI = 1.772453851  !sqrt(Pi)
      REAL, PARAMETER :: PI4 = 12.5663706   !4 * Pi
c      DATA WREDISUM /4.000420524/
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5  !Stephan-Boltzmann constant (CGS) / Pi
      REAL, PARAMETER :: CLIGHT = 2.99792458E5  !CLIGHT = VELOCITY OF LIGHT IN KM/S
      REAL, PARAMETER :: CC  = 2.9979E10    !C IN CGS UNITS
      REAL, PARAMETER :: AMU = 1.6605E-24   !ATOMIC MASS UNIT (~Hydrogen MASS) in g
      REAL, PARAMETER :: GCONST = 6.6727E-8 !Gravitational Constant (CGS)
      REAL, PARAMETER :: RGAS = 8.3145E7    !Gas Constant (CGS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !Solar mass (CGS = in g)

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hPLOT = 2       !write to plot file
      INTEGER, PARAMETER :: hMODEL = 3      !integer handle of the MODEL file (fort.3)
      INTEGER, PARAMETER :: hDMFILE = 16    !integer handle of the DMFILE file
      INTEGER, PARAMETER :: hEDDI = 17      !integer handle of the EDDI file (fort.17)
      INTEGER, PARAMETER :: hFEDAT = 21     !integer handle of the FEDAT file (fort.21)
      INTEGER, PARAMETER :: hHIST = 22      !write to MODHIST file
      INTEGER, PARAMETER :: hALO = 23       !integer handle of the ALO file -- added July 2023
      INTEGER, PARAMETER :: hFF = 27        !integer handle for FFASSET

      REAL, EXTERNAL :: BNUE
      INTEGER, EXTERNAL :: IDX

C***  Write Link Data (Program Version) tp CPR file
      WRITE (hCPR,'(2A)') 
     >      '>>> COLI started: Program Version from ', LINK_DATE
      WRITE (hCPR,'(4A)') 
     >      '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

C***  Preparations for Multi-COLI Option
      bForceCOLIP = .FALSE. !Force Coli+ instead of Coli
      bForceCLEERR = .FALSE.!Force creation of new EDDI file if true
      bKATEST = .TRUE.      !Test run (0.9*v, 0.9*v') first
      bKADONE = .FALSE.     !Test modifications done
      bNewLoop = .FALSE.    !default value for multiple coli runs
      bWARNK = .FALSE.

C***  If no Generic Element (Fe) is accounted for:
      LASTFE = 0 
      
C***  Entry Point for multiple COLI
 1    CONTINUE
      bFELASER = .FALSE.

      CALL INSTALL

ccc   obsolete and not gfortran-compatible; therefore commented:
c      IF (OPSYS .EQ. 'CRAY') THEN
c        CALL CLOCK(TIM1)
c      ELSE
c        CALL TIME(TIM1)
c      ENDIF

      CALL COLI_SETZERO(ND, NDDIM, NIT, NPDIM, NFDIM, MAXATOM, MAXFEIND,
     >                  MAXLIN, MAXIND, MAXION,
     >                  DBDTINT, DBDTOPAINT, EDDIHOUTJMEAN, 
     >             HTOTOUTMINUS, DBDTINT_M, DBDTOPAINT_M,
C***  with ND
     >             OPA, ETA,
     >             XJTOTL, HTOTL, XKTOTL, XNTOTL, ARAD, ACONT, ATHOM,
     >             FTCOLI, FTCOLIB, OPAKOLD, ETAKOLD, OPAKNOTHO, 
     >             ETAKNOTHO, OPAO, THOMSONO, ETANOTHO, 
     >             DJDSMOD_OLD, 
     >             OPASMEAN, QFJMEAN, OPAJMEAN, OPASMEANTC, OPAJMEANTC,
     >             OPAPMEAN, SMEAN, QLFOLD, EPSGMAX, OPAROSS, 
     >             OPALAMBDAMEAN,
     >             OPAROSSELEM, OPAROSSCONT, 
     >             OPAKFEOLD, iIgnoreK, OPAKFEFTOLD, ETAKFEFTOLD,
     >             OPAKNOFENOTHO, ETAKNOFENOTHO, 
     >             bOSKIPLAST, XLAMLASTOSKIP, XJLMO, XJLMOR2,
C***  with ND-1
     >             QOPAHMEAN, HMEAN, QLHOLD, OPAKHOLD, 
     >             HTOTLTEST, HTOTCUT,
C***  with NDDIM,NIT
     >             XJLOLD, XJLMO_OLD, EDDIFO, S_OLD, OPAK_OLD, EPSG,
C***  with NDDIM,NPDIM
     >             CWM0, CWM2, CWM1, CWM3,
C***  with NDDIM,NFDIM
     >             WJC,
C***  with NDDIM,NPDIM,NIT
     >             XIPLUS_OLD, XIMINUS_OLD,
C***  with NDDIM-1, NIT; In COLI it is also NDDIM,NIT
     >             XHLOLD, XHLMO_OLD, EDDIGO,
C***  with MAXFEIND, NDDIM
     >             XJFEMEAN, FERATLU, FERATUL, FTFE, WFELOW, WFENUP,
C***  with MAXATOM,ND
     >             OPAKELEM, OPAKOLDELEM, OPACELEM, OPACOLDELEM,
     >             ETAKELEM, ETAKOLDELEM, ETACELEM, ETACOLDELEM,
     >             OPATOTELEM,
C***  with MAXATOM,ND-1
     >             ARADELEM, ACONTELEM,
C***  with ND, MAXATOM, MAXION
     >             OPAKION, OPAKOLDION, OPACION, OPACOLDION,
     >             ETAKION, ETAKOLDION, ETACION, ETACOLDION,
C***  with ND-1, MAXATOM, MAXION
     >             ARADION, ACONTION,
C***  with MAXLIN
     >             LIND, LINDS, 
C***  with NDDIM, MAXLIN     
     >             WS, 
C***  with MAXIND  
     >             BLASERL, 
C***  with NFDIM 
     >             EMCOLI, 
C***  no Arrays
     >             HTOTMINUSND, HTOTND, HTOTNDS, HTOTNDCOR,
     >             OPAMAX1, OPAMAX1_LAMBDA, IOPAMAX1_K)

      
C***  -------------------------------------------------------------
C***  MAXIMUM HALF-BANDWIDTH FOR CMF LINE TRANSFER IN DOPPLER UNITS
      CMFBAND = 4.5 
C***  -------------------------------------------------------------
C***  MAXIMUM SPACING OF LINE-FREQUENCY POINTS
      DXMAX = 0.3
C***  -------------------------------------------------------------
C***  DUMMY FOR BACKJCU-CALL
      BLFERR = .FALSE.

C***  Initialize BELIFI; To Store Feautrier Matrices in Memory is now default
      BELIFI = .FALSE.

C***  Initialise BSTATIC; To calculate COLIRAYwith frequency deriv.
      BSTATIC = .FALSE.
c      bstatic = .true.

C***  Switch for the Kudritzki boundary condition
      BKUDRITZKI = .FALSE.

C!!!  test: call with routine name 'STEAL' instead of 'COLI' --> crash
ccc   wrh 31-Mar-2023

C***  READ ATOMIC DATA
       CALL DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $             EINST, ALPHA, SEXPO,
     $             ADDCON1, ADDCON2, ADDCON3, 
     $             IGAUNT, COCO, KEYCBB, ALTESUM,
     $             INDNUP, INDLOW, LASTIND, MAXIND, MAXATOM, NATOM,
     $             ELEMENT, SYMBOL, NOM, KODAT, ATMASS, STAGE,
     $             SIGMATHK, SEXPOK, EDGEK, NFIRST,
     $             NLAST, NAUTO, MAXAUTO, LOWAUTO, WAUTO, EAUTO, AAUTO,
     $             IONAUTO, KRUDAUT, KONTNUP, KONTLOW, LASTKON, MAXKONT,
     $             IONGRND, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >             'COLI', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

      LASTINDAUTO = LASTIND + NAUTO
      LASTINDALL  = LASTIND + NAUTO + LASTFE

C***  READ THE MODEL FILE
      CALL       RMODCOLI (RADIUS,ENTOT,RNE,T,VELO,GRADI,XLAMBDA, 
     >                     FWEIGHT, POPNUM, POPMIN,
     >                     RSTAR, VDOP, MODHEAD, JOBNUM, XJC,
     >                     P,ND,NDDIM,NF,NFDIM,N,NDIM,NP,NPDIM,Z,
     >                     TEFF,HWRCEDDI,EDDI,MODHIST,MAXHIST, 
     >                     DENSCON, FILLFAC, ABXYZ, NATOM, 
     >                     LASTREDISMODE, NCOLIP, NF2, XLAMBDA2, 
     >                     OPARND, EPSGMAX, BEPSGMAXERR, MAXXDAT, 
     >                     XDATA,XMSTAR, ZERO_RATES, HTOTMINUSND, 
     >                     HTOTNDCOR, HTOTLlast, TAULAST, 
     >                     NEXTHYDRO, DTDRIN, VMIC)
      WRITE(hCPR,'(A,I7)') '>>> This is job number ', JOBNUM
      XMG = GCONST * XMSTAR * XMSUN

C***  ALOs are only used in the next STEAL, so always create new ALO file
      CALL OPENEXMS (hALO, IDUMMY, IDUMMY, 'ALO', 'REPLACE', IERR)

C***  Set popnumbers at POPMIN to zero to avoid any effect on the
C***  radiative transfer calculations:
      CALL POPMIN_NULLING(ZERO_RATES, POPNUM, POPMIN, ND, N)
      
      !Determine if force multipliers are calculated
      ! this is either done if the next STEAL job 
      ! performs a hydro calculation (NEXTHYDRO = 1) or if
      ! forced by the special CARDS option FORCEMULTIPLIERS
      bKALPHA = .FALSE.
      IF (NEXTHYDRO == -13) THEN
        bKALPHA = .TRUE.        !alpha calulation for converged model
      ENDIF

C***  DECODE INPUT CARDS
      CALL DECCOLI (LSOPA, LSINT, MODHIST,
     >       NOLAP, LPOPAB, LPOPABD, LEVDEBUG, LPJNUE,
     >       LPJNUED, LASERV, PARALAS, LPSNUE, LPSNUED,
     $       MAXPLOT, RANGE1, RANGE2, EXLAM1, EXLAM2, MAXEXT,
     $       BLLIST,
     >       REDISMODE, NEWWRC, BCOLIRAY, BKUDRITZKI,
     >       CLHLP, BITCONT, BPLOT,  
     >       IPLOT, LPLOT, ND, OPC,
     >       IVERS, BEMIX, EMIXSTART, BEMIXFIX, EMIXFIX, 
     >       IVERS_FE_EXPFAC, BPLOTALPHA, ALPOPT, POPMIN,
     >       XLAM_FINE_START, XLAM_FINE_END, LPLOT_WCHARM, XLP1,
     >       XLP2, GAMMACOLI, GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2, 
     >       bKALPHA, bHYDROSOLVE, VELOMODFAK, iTypeAlpha, 
     >       bForceCOLIP, bTDIFFUS, bPLOTRTAU1, bNoIronLaser, 
     >       iHTOTCUT, bMAXEMIX, bALOTri,
     >       bNoNEGEDDIG, bDDVDOP, LPRDH,
     >       bDEBUG, CUTOPAMEANTHRES, DRLINES_CARD, EMIXMAX)

      DO LINEINDEX = 1, LASTIND+NAUTO
         LINE(LINEINDEX) = LINEINDEX
      ENDDO
     
C***  Setup VDOPDD array (for line profiles) and VDOPUNIT (grid spacing) variable      
      VDOPUNIT = VDOP
      CALL VDOPDD_SETUP(bDDVDOP, VDOPDD, VDOPUNIT, DXMAX,
     >                  T, VMIC, ATMASS, ND, NATOM)    
     
      IF ((.NOT. bKALPHA) 
     >      .AND. (bHYDROSOLVE) .AND. (NEXTHYDRO == 1)) THEN
        bKALPHA = .TRUE.
      ENDIF
      bAdjVINF = .FALSE.
      IF (VELOMODFAK < 0.) THEN
        VELOMODFAK = -1. * VELOMODFAK
        bAdjVINF = .TRUE.
      ENDIF
      IF (iTypeAlpha == 0) THEN
        !Switch off alpha/response calculation if not used in HYDRO step 
        bKALPHA = .FALSE.
      ELSE
        WRITE (hCPR,*) ' Requested Alpha-Type', iTypeAlpha
      ENDIF
      

C***  Open special Plot file wcharm.dat
      IF (LPLOT_WCHARM .GT. 0) THEN
        OPEN (UNIT=105, FILE='wcharm.dat', STATUS='UNKNOWN')
        WRITE (105,'(A,I3)') 
     >    'KASDEF LUN XMAX YMAX 0.5 0.0 0.3 L = ', LPLOT_WCHARM
        WRITE (105,'(A)') '* Fine Information: Lambda, WJC'
        WRITE (105,'(A)') 'N=?'
      ENDIF      


      EXLAM1 = XLAMBDA(1)
      EXLAM2 = XLAMBDA(NF)

      REDISMODED = REDISMODE

C***  prepare the rudimental settings for DRTRANSIT lines
      IF (NAUTO > 0) THEN
         CALL PREP_DRLINES (DRLINES_CARD, NAUTO, KRUDAUT, EAUTO)

         N_RUDAUTO = 0
         DO INDDR=1, NAUTO
            N_RUDAUTO = N_RUDAUTO + KRUDAUT(INDDR)
         ENDDO
         WRITE (0,*) 'DRLINES_CARD = '//DRLINES_CARD(:IDX(DRLINES_CARD))
         WRITE (0,'(A, I5)') 'Number of DRTRANSITs :', NAUTO
         WRITE (0,'(A, I5)') 'of those "rudimental":', N_RUDAUTO

C***  INSERT UPPER DR-LEVELS
        CALL DRLEVEL (N, NDIM, MAXIND, MAXAUTO, NAUTO, KRUDAUT,
     $                LOWAUTO, IONGRND, AAUTO, ELEVEL, 
     $                LEVEL, EINST, EION, WEIGHT, INDLOW, 
     $                INDNUP, LASTIND,  
     $                ND, T, ENTOT, RNE, POPNUM, DENSCON, 
     >                LEVUPAUTO, LEVAUTO, WAUTO, N_WITH_DRLEVELS)
      ENDIF

ccc   test output
c         do iline= 1, LASTIND+LASTFE+NAUTO, 100
c           write (0,'(2I6,x,A,x,I6,x,A,x,3F10.2)') 
c     >      iline, indlow(iline), LEVEL(indlow(iline)), 
c     >             indnup(iline), LEVEL(indnup(iline)),
c     >             elevel(indlow(iline)), elevel(indnup(iline)),
c     >             1.E8/(elevel(indnup(iline))-elevel(indlow(iline)))
c         enddo

C***  HALF BANDWIDTH OF ELECTRON REDISTRIBUTION INTEGRAL
C***     IN UNITS OF THE ELECTRON DOPPLER VELOCITY
C***     IN FORMCMF BWESEX CAN BE SET BY AN INPUT CARD, HERE IT
C***     IS SET TO 1.
      BWESEX = 1.
      BWES = 2. * BWESEX

C***  THE FREQUENCY RANGES (BOTH, CMF AND OBS. FRAME) ARE EXTENDED
C***      BY ONE HALF BANDWIDTH, MULTIPLIED WITH THE FOLLOWING FACTORS:
      BWEXRED = 2.
      BWEXBLU = 1.5

      !Fixed value for VINF (even if velocity field is modified later)      
      VINF = VELO(1)/VDOPUNIT

C***  IRON: WIDTH OF RED LINE WING IN IRON-DELTAX-UNITS
C      IF (BFEMODEL) THEN
C         DFEINDR = 2.*VELO(1)/VDOPFE/DXFE
C      ENDIF

      ATMEAN = 0.
      DO NA=1, NATOM
        ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO

C***  INTRODUCING DIMENSIONLESS VELOCITY UNITS : 
      DO L=1,ND
        VELO(L)=VELO(L)/VDOPUNIT
        GRADI(L)=GRADI(L)/VDOPUNIT
        RHO(L) = AMU*ATMEAN*ENTOT(L)
      ENDDO
 
C***  MASS STORAGE ON FILE 7 FOR THE FEAUTRIER MATRICES (WRCONT)
C***  ARRAY SIZE: SUM OVER 'LL' (SEE SUBR. ELIMIN) FROM L=1 TO ND-1
C***  NUMBER OF WORDS:
      KWORDS=(ND-1)*((NP+2)*(6*NP+6)-ND*(6*NP+9)+ND*(2*ND-1))/6
C***  NUMBER OF BLOCKS (1 BLOCK = 512 WORDS) IN MEMORY
      KBLOCKS=(KWORDS+511)/512+1
C***  SET MEMORY SIZE (ISTATS=1: COLLECT STATISTICS)
      ISTATS=1
      IF (BELIFI) THEN
        CALL OPENMS (7, IDUMMY, IDUMMY, 0, IERR)
      ENDIF

      CMODE = 'N'
C***  Setup all Variables for treatment od Redistribution
      IFRO = 0
      NXJO = 0
      NXK = 0
      NXJO = 0
C***  INITIALIZE XJLLOAD
      LACTION = 'LOAD'

C***  Dimension for EDDIA
C***    Note : XHI and XHO are also stored 
C***           XHI for the inner boundary in COLIMO
C***           XHO for Unsoeld-Lucy
C***           XHOM and EDDIHOUTP for special treatment of the outer boundary
C***           EDDIHINM for calculation of HTOTMINUSND (needed for innermost temperature correction)
      NDEDDIA = 2*ND + 3 + 2 + 4 + 2

C***  Counter for Number of EDDI Resets in COLIMO
      IW_COLIMO_F      = 0
      IW_COLIMO_G      = 0
      IW_COLIMO_G2     = 0
      IW_COLIMO_J      = 0
      IW_COLIRAY_IPLUS = 0
      IW_COLIRAY_U     = 0

C***  Decide on wether an COLI (no ray-by-ray) or an 
C***                      COLI+ (with ray-by-ray) is performed
C***  Default : COLI, i.e. no calculation of new EDDIEs
      BCOLIP = .FALSE.
      IF (bForceCOLIP) BCOLIP = .TRUE.

C***  Short Characteristics?
C***  Note: Other branch (COLIRAY) is not longer valid
      BSHORT_CHAR = .TRUE.

C***  Extended spacing at Continuum points which are not due to edges?
      BXJCE = .FALSE.

      IF (JOBNUM .GT. 3 .AND. JOBNUM .LT. 10 .AND. 
     >    .NOT. CLHLP .AND. BITCONT) THEN
        BSHORT = .TRUE. 
      ELSE
        BSHORT = .FALSE.
      ENDIF
      IF (BSHORT) THEN
        WRITE (hCPR,*) 'This COLI was a SHORT one!'
        GOTO 100
      ENDIF

C***  Check the jobnumber of last COLI+ (NOT DONE!!!!)
      IF (.NOT. bNewLoop) THEN
        !Important: Do not call CLOPENE again in multiple COLI runs
        ! and keep BCLEERR status from first call
        ! CLOPENE also updates MODEL string in EDDI file, but the MODEL
        ! is always the same in multiple COLI runs, so this should not be a problem
        BCLEERR = .FALSE.
        CALL CLOPENE (hEDDI, MODHEAD, JOBNUM, 
     >                NCOLIP, NDEDDIA, BCLEERR)
        WRITE (hCPR,*) 'Opening EDDI file...'
C***    Open ALO file
      ENDIF
      
C***  Determine type of COLI (normal COLI, COLI+, or COLI++)
C***  by evaluating the current status of the EDDIs and CARDS options
      IF (BCLEERR .OR. bForceCLEERR) THEN
C***    EDDI file needs to be created or renewed => COLI++     
        BCOLIP = .TRUE.
        BCOLIPP = .TRUE.
        IF (bForceCLEERR) THEN
          WRITE (hCPR,'(A,A4)') 
     >      ' Creation of new EDDI file forced: COLI++ required'
        ELSE
          WRITE (hCPR,'(A,A4)') 
     >      ' No valid EDDI file found: COLI++ required'
        ENDIF
      ELSE
C***    Check for at least one of the criteria for a COLI+      
        BCOLIPP = .FALSE.
        IF (BCOLIP) THEN
          WRITE (hCPR,'(A,I2,A)') 
     >      ' COLI+ forced'
        ELSEIF (NCOLIP == 0) THEN
          WRITE (hCPR,'(A,I2,A)') 
     >      ' Previous WRCONT detected: performing COLI+'
          BCOLIP = .TRUE.
        ELSEIF (NCOLIP >= NEWWRC) THEN
            WRITE (hCPR,'(A)') 
     >      ' EDDIs too old (NCOLIP >= NEWWRC): COLI+ required'
          BCOLIP = .TRUE.
          NCOLIP = 0
        ENDIF
      ENDIF


C***  VELOCITY SCALING FOR THE CALCULATION OF FORCE MULTIPLIERS
      IF (bKALPHA .AND. bKATEST) THEN
         IF (iTypeAlpha == 3) THEN
           CATMOD = ' DENSITY'
         ELSEIF (iTypeAlpha /= 11) THEN
           CATMOD = 'GRADIENT'
         ELSE
           CATMOD = 'VELOCITY'
         ENDIF
         WRITE (hCPR,'(A)') 
     >       ' *********************************************'
         WRITE (hCPR,'(A)')
     >       ' * CALCULATION OF FORCE MULTIPLIERS:         *'
         WRITE (hCPR,'(A,2X,F14.11,X,A)') 
     >       ' * '// CATMOD //' IS MODIFIED BY ',  VELOMODFAK, ' *'
         WRITE (hCPR,'(A)')
     >       ' *********************************************'
         IF ((iTypeAlpha /= 3) .AND. 
     >         (VELOMODFAK > 1. .OR. VELOMODFAK <= 0.)) THEN
c           WRITE (hCPR,'(2A)') 'FATAL ERROR: Velocity modifier must ',
c     >              'be greater than zero and lower than 1!'
c           STOP 'Fatal ERROR in Subr. COLI'
         ELSEIF (VELOMODFAK == 1.) THEN
           WRITE (hCPR,'(2A)') '**WARNING** GRADIENT is not modified!'
         ENDIF
         DO L=1, ND
            IF (iTypeAlpha /= 3) THEN
              IF (iTypeAlpha /= 12) THEN
                VELO(L)  = VELO(L)  * VELOMODFAK
              ENDIF
              IF (iTypeAlpha /= 11) THEN
                GRADI(L) = GRADI(L) * VELOMODFAK
              ENDIF
            ENDIF
            IF (iTypeAlpha == 3 .OR. iTypeAlpha == 31) THEN
              ENTOT(L) = ENTOT(L) / VELOMODFAK
            ENDIF
C            GRADI(L) = GRADI(L)
         ENDDO
         IF (bAdjVINF) THEN
           VINF = VELO(1)
         ENDIF
         bKADONE = .TRUE.
      ELSE
         bKATEST = .FALSE.
      ENDIF

C***  Reorder lines in sequence of increasing wavelengths
      NLINE = LASTIND + NAUTO
      CALL SEQLINECL(NLINE, LINE, EINST, INDLOW, 
     >       INDNUP, XLAMSOR, ELEVEL, NDIM, VDOPUNIT,
     >       CMFBAND, CLIGHT, XLAMMIN, XLAMMAX,
     >       VINF*VDOPUNIT, EXLAM1, EXLAM2, MAXEXT, 
     >       LASTIND, NAUTO, KRUDAUT )

C***  IRON: WIDTH OF RED LINE WING IN IRON-DELTAX-UNITS
      IF (BFEMODEL) THEN
         DFEINDR = 2.*VINF*VDOPUNIT/VDOPFE/DXFE
      ENDIF
     
      IF (BCOLIRAY) THEN
        WRITE (hCPR,*) 'Pure Coliray Branch is not longer valid'
        STOP 'Fatal ERROR in Subr. COLI'
      ENDIF

C***  Setup COLI iteration scheme
C      ITSTART  ITMAX
C         2       3        COLI++  (RAY-MO-RAY-MO)
C         1       2        COLI+   (MO-RAY-MO)
C         1       1        COLI    (MO)
      ITSTART = 1
      IF (BCOLIP) THEN
        IF (BCOLIPP) THEN
          ITSTART = 2
          ITMAX = 3
        ELSE
          ITMAX = 2
        ENDIF
      ELSE
        ITMAX = 1
      ENDIF

C***  Failsafe check: ITMAX must not be larger than dimensioned arrays      
      IF (ITMAX > NIT) THEN
        WRITE (hCPR,*) 'ITMAX > NIT : ', ITMAX, NIT
        STOP 'ERROR in Subr. COLI'
      ENDIF
      
C***  Increase COLI+ counter (stored in MODEL file)      
      IF (.NOT. bKATEST) THEN
        NCOLIP = NCOLIP + 1
C***    For COLI++ ensure that next COLI is a COLI+        
        IF (BCOLIPP) NCOLIP = 0
      ENDIF


C***  EDDIMIX factor for new EDDIG (if COLI+) is specified from MODEL-File, 
C***     or set ZERO, EMIXSTART or EMIXFIX, as requested 
C***     Note: Modification by wrh 17-Aug-2000 12:12:25 
C***           EMIXFIX is considered as minimum value, but EPSGMAX is applied  
C***           if neccessary  
      EPSGITMAX(1) = 0.

      DO IT=2, ITMAX
        DO L=1, ND-1     
          IF (.NOT. BEMIX) THEN
            EPSG(L,IT) = 0.   
          ELSE
            IF (BEMIXFIX) THEN
              EPSG(L,IT) = AMAX1 (EMIXFIX, EPSGMAX(L)) 
              IF (EMIXMAX .GT. .0) EPSG(L,IT)=AMIN1(EPSG(L,IT),EMIXMAX)
            ELSE
              IF (BEPSGMAXERR) THEN 
                EPSG(L,IT) = EMIXSTART
              ELSE
                EPSG(L,IT) = EPSGMAX(L)
              ENDIF
            ENDIF
          ENDIF          
          EPSGITMAX(IT) = MAX(EPSGITMAX(IT), EPSG(L,IT))
        ENDDO
c        IF (bMAXEMIX) THEN
cC***      New test option: Apply the maximum required EDDIMIX to ALL depth points        
c          DO L=1, ND-1
c            EPSG(L,IT) = EPSGITMAX(IT)
c          ENDDO
c        ENDIF
      ENDDO

      EPSGMAX = 0.  !needed for multiple call of COLI (from Goetz)  

C***  Initialize Counters for reading and writing the EDDIEs to 
C***    File fort.<hEDDI>
C***  EDDIA : Array where all eddis are stored
C***          NFRO frequencies are stored in one entry in the file
C***  NZE1 : is part of the Name of the arrays in fort.<hEDDI>
C***  NZE2 : is the actual Index in the array EDDIA
      NZE1 = 0
      NZE2 = 0


C***  CALCULATION OF "DTDR" FOR USE IN SUBR. CLDIFFUS
      IF (DTDRIN >= 0. .AND. (.NOT. bTDIFFUS)) THEN
C***    Use precise calculation from STEAL->CALCDRDRIN if available
C***     and TDIFFUS artistic is not switched on
        DTDR = DTDRIN
      ELSE
        CALL DIFDTDR (DTDR,TEFF,XJC,HWRCEDDI,T(ND),RADIUS(ND),
     >                ND,EN,POPNUM,POPMIN,RNE(ND),ENTOT(ND),RSTAR,
     >                NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     >                ALPHA,SEXPO,
     >                ADDCON1, ADDCON2, ADDCON3, 
     >                IGAUNT,NOM,NF,XLAMBDA,FWEIGHT, 
     >                MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >                KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC, 
     >                BKUDRITZKI, OPARND) 
      ENDIF

C***  NEWBGC = COUNTER OF LINES FOR WHICH NEW BACKGROUND CONT. RAD. 
C***           FIELD IS CALCULATED
      NEWBGC = 0

C***  IWARNJ0 = COUNTER OF ZERO BACKGROUND CONTINUUM INTENSITIES 
C***            ENCOUNTERED BY SUBR. BACKJCU
      IWARNJ0 = 0
 
C***  NLASER = COUNTER OF LINES FOR WHICH THE LASER CONDITION WAS 
C***           ENCOUNTERED
      NLASER = 0

      LBMAX=0
      NBLENDS=0
      
C**** IWARNJMNEG = Counter for negative XJL values obtained in COLIMO
C****  (Such values are now [Apr 2015] set to zero! )
      IWARNJMNEG = 0

C******************************************************************
C***  Start of new COLI code
C******************************************************************

C***  Testplots?
C***  BPLOT  = .TRUE. : Plot from Subroutines 
C***                      COLIMOP, PREREDIS, COLIMO
C***  BPLOT2 = .TRUE. : Plot from Subroutines 
C***                      COLIRAY, CMFCOOP
C***  Plot of the Lambda-spacing
C***    The Plotfile is coli_k.plot
C***    Geometry                  :  coli_geo.plot
C***    Lambda-spacing            :  coli_k.plot
C***    J,H,K... depth            :  coli_opa_l.plot
C***    J,H,K... frequency        :  coli_opa.plot
C***    Integrated J,H,K,L depth  :  coli_flux.plot

C***  BPLOT is set in DECCOLI; BPLOT2 is a hard switch
      BPLOT2 = .FALSE.
      IF (CLHLP)  BPLOT = .TRUE.
      IF (BPLOT2) BPLOT = .TRUE.
      IF (BPLOT) THEN
C***  Depth point for Plot over wavelength
c        IPLOT = ND / 2
c        IPLOT = 67
        RIPL2 = RADIUS(IPLOT)*RADIUS(IPLOT)
C***  Index for plot in depth. The next following after LPLOT is taken
c        LPLOT = 27190
        BPDONE = .FALSE.
        IF (BPLOT) THEN
          WRITE (hCPR,'(A,I3,A)') 'Plot in depth ', IPLOT, ' prepared'
          WRITE (hCPR,'(A,I6,A)') 
     >                          'Plot at frequency ',LPLOT,' prepared'
          WRITE (hCPR,'(A,2(1X,F12.5))') 'Plot-Range [A]: ',XLP1,XLP2
          DO NL=1, MAXLIN
            IF (NL .LT. 10) THEN
              WRITE (NAME,'(A3,I1,A4)') 'cl_', nl, '    '
            ELSE
              WRITE (NAME,'(A3,I2,A3)') 'cl_', nl, '   '
            ENDIF
            OPEN(UNIT=40+NL, FILE=NAME)
            WRITE (40+NL,'(A3)') 'N=?'
          ENDDO
          OPEN (UNIT=40, FILE='coli.dat')
          OPEN (UNIT=39, FILE='coli_cont.dat')
          OPEN (UNIT=38, FILE='coli.kasdef')
          OPEN (UNIT=37, FILE='coli_opa.dat')
          OPEN (UNIT=36, FILE='coli_opa_l.dat')
          OPEN (UNIT=35, FILE='coli_flux_l.dat')
          IF (BPLOT2) THEN
            IF (MAXLIN .GT. 49) THEN
              WRITE (hCPR,*) 'WARNING : File Number Mismatch in COLI'
            ENDIF
            OPEN (UNIT=90, FILE='colimo.dat')
            OPEN (UNIT=91, FILE='coliray.dat')
            OPEN (UNIT=92, FILE='cmfcoop.dat')
          ENDIF
          WRITE (37,'(A6,27A16)')
     >    ' * IND', 'LAMBDA', 'OPA-KONT', 
     >    'ETA-K', 'ETANOTH-K.', 'S-KONT', 
     >    'S-KONT(NOTH)', 'OPA-L', 'ETA-L', 'S-L', 'J-KONT', 'J', 
     >    'J-SCHLANGE', 'H-S', 'K-S', 'N-S', 'EDDIF', 'EDDIG', 
     >    'XJLMO', 'XJLMOR2', 'XHLMO', 
     >    'FULFIL0', 'FULFIL1', 'DJDSMO', 'DJDSMO', 'DJDS', 'SLNOTH'
          WRITE (36,'(A6,26A16)')
     >    ' * L', 'LAMBDA', 'OPA-KONT', 
     >    'ETA-K', 'ETANOTH-K.', 'S-KONT', 
     >    'S-KONT(NOTH)', 'OPA-L', 'ETA-L', 'S-L', 'J-KONT', 'J', 
     >    'J-SCHLANGE', 'H-S', 'K-S', 'N-S', 
     >    'EDDIF', 'EDDIG', 'XJLMO', 'XJLMOR2', 'XHLMO', 
     >    'FULFIL0', 'FULFIL1', 'DJDSMO', 'RADIUS'
          WRITE (35,'(A4,5A16)')
     >    ' * L', 'JTOT', 'HTOT', 'KTOT', 'NTOT', 'RADIUS'
          WRITE (40,'(A3)') 'N=?'
          WRITE (39,'(A3)') 'N=?'
          WRITE (38,'(A3)') 'N=?'
          WRITE (37,'(A3)') 'N=?'
          WRITE (36,'(A3)') 'N=?'
          WRITE (35,'(A3)') 'N=?'
          IF (BPLOT2) THEN
            WRITE (60,'(A3)') 'N=?'
            WRITE (61,'(A3)') 'N=?'
            WRITE (62,'(A3)') 'N=?'
          ENDIF
        ENDIF
      ENDIF

      
C***  Preparation of flux integration weights (test purpose only)
      CALL GENWP1 (NP, P, WP1, WP1LAST)
      CALL COLIWM(Z, P, RADIUS, ND, NP, 
     >            CWM0, CWM1, CWM2, CWM3, 
     >            BSHORT_CHAR)

C***  The weights could be checked 
C      IF (.FALSE.) THEN
C        CALL CHECKWM (CWM0, CWM1, CWM2, CWM3, 
C     >                CWM1O, CWM1I, CWM3O, CWM3I, 
C     >                WP1, WP1LAST, W0, RADIUS, ND, NP, Z, P)
C      ENDIF

C***  Precalculation of the geometry for COLIMO
      CALL COLIMOP(ND, RADIUS, GRADI, VELO, 
     >                   DLF, DLH, GLF2, GLH2, VLF2, VLH2, BPLOT,
     >                   RADIUS2, RADIUSH, RADIUSH2)

C     OLD CALL COLI_SETZERO position

      CALL PREPK(
C***  Input
     >    XLAMBDA, NF, XLAMBDA2, NF2, VDOPUNIT, CLIGHT, DXMAX,
     >    VINF, ND, CMFBAND, NLINE, MAXIND, XLAMSOR,
     >    BPLOT, BXJCE, LINE, 
C***  Output 
     >    XLAM0, XLAM0LN, ALN, CMFBANDR,
     >    XKMID, XKMIN, XKRED, XKMAX, XKC, XKC2, XKCMAX,
     >    KSPACE, BANDP, BANDM, KOUT, K, KLAST, NK,
     >    KONCHECK, KONCHECK2, KONTACT, KONTAUP, NLACT,
     >    LINECHECK, ILINECHECK,
     >    KCL, KCU, KCDMAX, KCDELTA, KCCHECK, BFF_ACT, IFF_N, 
     >    IFF_DK, IFF_MAX)
      
      
      IF (BCOLIPP) THEN
        CPPLABEL = '++'
      ELSE
        CPPLABEL = '  '
      ENDIF
      WRITE (hCPR,*)
      WRITE (hCPR,'(13A8)')
     >  'COLI+', 'BCLEERR', 'EDDIMIX', 
     >  'TDIFFUS', 'OB-VERS', 'KSPACE', 'BXJCE',
     >  'BPLOT', 'BPLOT2', 'BSHORT', 'CLHLP', 'VERS_FE'
      WRITE (hCPR,'(L6,A2,L5,3L8,2I8,5L8,I7)')
     >  BCOLIP, CPPLABEL, BCLEERR, BEMIX, 
     >  bTDIFFUS, IVERS, KSPACE, BXJCE,
     >  BPLOT, BPLOT2, BSHORT, CLHLP, IVERS_FE_EXPFAC
      WRITE (hCPR,*)


      WRITE (hCPR,'(A,I3,1X,I3)') 
     >  'Iteration Scheme : ITSTART, ITMAX = ', ITSTART, ITMAX

cc      WRITE (hCPR,*)
cc      WRITE (hCPR,'(5A12)') 
cc     >  'Operator:   ', 'Continuum', 'Lines', 'R-Lines', 'Fe-Lines'
cc      WRITE (hCPR,'(2A12,3I12)')
cc     >  ' ', OPC, NLINE, NUMRUD, LASTFE
cc      WRITE (hCPR,*)

C***  Prepare the Indices IF1..4 which are controlling the Storage of 
C***    the old Js from the last Iteration
C***    IF1 : First Index
C***    IF2 : Number of Elements
C***    IF3 : First Index relevant for the redistribution
C***    IF4 : Number of Elements relevant for the Redistribution
      IF1 = 1
      IF2 = 0
      IF3 = 1
      IF4 = 0
C***  Number of frequency-indicees needed for the redistribution
C***    For CONTINUUM or COHERENCE the default is 1
      DKR = 1.

C***  IRON: SET LOWEST LINE TO '1', NO ACTIVE FE-LINES 
      INDFEACT(1) = 1
      MAXFEACT    = 0
      BFECHECK    = .FALSE.
      BFEWING     = .FALSE.

      HTOTNDCOR = 0.

C*****************************************************
C***  Main-loop over all frequency-indices
C*****************************************************

      freqloop: DO !- - - - - - - - - - - - - - - - - - - - - - - - - -

C***  Preset PPP
      IF (K .LT. 2) THEN
        CALL PREP_PPP(ND, NP, NDDIM, NPDIM,
     >                PPP, Z, RADIUS, GRADI, VELO,
     >                K, BSTATIC)
      ENDIF

        XK = FLOAT(K)
        IF (K .EQ. 0) THEN
          XLAMK = XLAM0
C          XLAMKOLD = EXP(XLAM0LN - ALN)        !this line has been taken from Goetz
        ELSE
          XLAMK = EXP(XLAM0LN + XK*ALN)
        ENDIF
        XNK = FLOAT(NK)

C***  IRON: Test for active FE-Transitions
ccc      INDBAK = INDFEACT(1)
      IF (BFEMODEL)
     >      CALL FECHECK (XLAMK, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEINDR)


      CALL CHECK_CONT(BKONCHECK, BKONCHECK2,
     >             KONCHECK, KONCHECK2, NF, NF2, K, XK,
     >             XLAMKOLD, XKC, XKC2, BANDM, BANDP,
C***  Parameters for Continuum Interpolation 
     >             CMFBAND, CMFBANDR, DXMAX,
     >             KCCHECK, KCL, KCU, KCONT, KCDMAX, KCDELTA, 
     >             XLAMBDA, KONTACT, KONTAUP, DFKONT, BPLOT)
                 
C***  Explaining Variables :
C***    XLAMSOR      : Array of sorted wavelengths (1..NLINE)
C***    XKMIN, XKMAX : Indices where the starts and ends
C***    LINECHECK    : Actual line which is checked
C***    ILINECHECK   : Index of actual line in original list
C***    LINDS        : Array of linecheck (1..MAXLIN)
C***    LIND         : Array of ILINECHECK (1..MAXLIN)
C***    BLASERL      : Is true for Laser Lines

      CALL CHECK_LINES(XK, XKMIN, XKMID, XKMAX,
     >             LINECHECK, ILINECHECK, NLINE, MAXLIN, LEVEL, 
     >             LIND, LINDS, WS, BPLOT, RADIUS, NLACT, LINE, 
C***  Parameters for PRELINECL
     >             NUP, LOW, N, XLAM, NDIM, ND, XJLMEAN, XLAMAPPMEAN,
     >             bALOTri, XLAMAPPUMEAN, XLAMAPPLMEAN,
     >             ELEVEL, INDNUP, INDLOW, NDDIM,
C***  Parameters and also for LIOP 
     >             EINST, WEIGHT, XLAMSOR, ENTOT, POPNUM, RSTAR,
     >             OPAL, ETAL, VDOPUNIT, hMODEL, hALO)

CC***  Output of K
C        IF (K .GT. KOUT) THEN
C          WRITE (hCPR,'(A,I7,1X,I7,1X,F15.3)') 'K=', K, NK, XLAMK
C          KOUT = KOUT + 20000
C       ENDIF

C***  Force the writing of XJ and ED Arrays
      IF (FLOAT(K+1) .GT. XKCMAX) THEN
        CMODE = 'F'
ccc        write (hCPR,*) 'FLOAT(K) .EQ. XKCMAX -> CMODE=', CMODE
      ENDIF

C***  Check if actual frequency should be skipped
C***    Special treatment for Iron-Lines: TEST FOR ACTIVE WING
        IF (.NOT. BKONCHECK .AND. .NOT. BKONCHECK2 .AND.
     >      .NOT. BFEWING
     >      .AND. NLACT .EQ. 0 .AND. 
     >      K-KLAST .LT. KSPACE .AND. K .GT. 0 .AND. .NOT. 
     >      CMODE .EQ. 'F') THEN 
          K = K + 1
          IF (K .GT. KCU) THEN
C***        New Start Indices for Continuum Interpolation (KCL = KCU)
             KCL = K
             KCU = K
          ENDIF
          CYCLE
        ENDIF
        IF (BPLOT) THEN
          WRITE (40,'(I8,1X,F3.0)') K, 0.
        ENDIF

C***  Main loop is now restricted to necessary points
C***  Now the radiation transfer is calculated for index K
C***  Prepared quantities :
C***    K, XLAMK                   actual index and its wavelength
C***    BKONCHECK,  KONCHECK       near kontinuum and index of next kontinuum
C***    BKONCHECK2, KONCHECK2      near kontinuum and index of next true kontinuum
C***    NLACT, LIND(NL(1..MAXLIN)) Number of active lines and indices of the lines
C***    

C***  Check if actual K is in the Fine Integration Interval
      IF (XLAMK .GE. XLAM_FINE_START .AND.
     >    XLAMKOLD .LT. XLAM_FINE_START) THEN
        BFF_ACT = .TRUE.
        FF_INFO(5) = XK
        FF_INFO(8) = XLAMK
      ENDIF
      IF (XLAMK .GT. XLAM_FINE_END .AND.
     >    XLAMKOLD .LE. XLAM_FINE_END) THEN
        BFF_ACT = .FALSE.
        FF_INFO(6) = FLOAT(KLAST)
        FF_INFO(9) = XLAMKOLD
      ENDIF
      IF (BFF_ACT) THEN
        IFF_N = IFF_N + 1
        IF (IFF_N .GT. IFF_MAX) THEN
          WRITE (hCPR,*) 'Number of fine frequencies for ',
     >                'Integration in COLI exceeded'
          WRITE (hCPR,'(A,I8)') 'Present IFF_MAX = ', IFF_MAX
          STOP 'ERROR in Subr. COLI'
        ENDIF
      ENDIF

C***  Read old EDDIEs from file fort.<hEDDI>
      CALL CLLOADE (hEDDI, NZE1, NZE2, NFRO, EDDIA, NDEDDIA,
     >              EDDIF, EDDIG, ND,
     >              EDDIHOUT, EDDIHIN, EDDIHINM, EDDINOUT, EDDININ,
     >              BCLEERR, BCOLIP, XHI, XHO, EPSG(1,1), 
     >              XHOM, XNOM, EDDIHOUTP, EDDINOUTP, EDDIHINT)

C***  Overview on current EPSG
      IF (K .EQ. 0 .AND. BCOLIP) THEN
        IF (BEPSGMAXERR .AND. BEMIX) THEN 
          WRITE (hCPR, '(A,F5.2)') 
     >        'COLI: EDDIMIX started with EPSG=', EMIXSTART
        ENDIF
        CALL PRI_EPSG(BEMIX, BEMIXFIX, ITMAX, EPSG, NDDIM, ND, NIT,
     >                JOBNUM)
      ENDIF

        CALL BACKJC (XJC, ND, NF, XJCIND, XLAMBDA, XLAMK, RADIUS)
        CALL CMFCOOP (XLAMK,ND,T,RNE,POPNUM,ENTOT,RSTAR, LEVEL,
     >         OPA,ETANOTH,THOMSON,NDIM,N,NCHARG,WEIGHT,ELEVEL,EION,
     >         EINST,ALPHA,SEXPO,ADDCON1,ADDCON2,ADDCON3,IGAUNT,
     >         SIGMA1I,KONTLOW,KONTNUP,LASTKON,NATOM,KONTHLP,
     >         DENSCON,BPLOT,BPLOT2,IPLOT,K,KCL,KCU,KCDELTA,
     >         OPACL,OPACU,ETACL,ETACU,XLAM0LN,ALN, MAXXDAT, XDATA, 
     >         SIGMATHK,SEXPOK,EDGEK, MAXATOM, NOM, KODAT, RADIUS,
     >         MAINPRO, MAINLEV, IWARN, OPACLEL, OPACUEL, OPACELEM,
     >         ETACLEL, ETACUEL, ETACELEM, OPACLION, OPACUION,
     >         OPACION, ETACLION, ETACUION, ETACION,
     >         OPALEL, ETALEL, OPALION, ETALION)

C***  IRON: CALCULATE FE-OPACITY AND EMISSIVITY
        IF (BFECHECK) THEN
           CALL  CMFFEOP (XLAMK, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDEXMAX, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE, 
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT, 
     >                    SIGMAACT, OPAFEI, ETAFEI, 
     >                    OPAFEION, ETAFEION, T, IVERS_FE_EXPFAC, 
     >                    TEFF, NCHARG, MAXION, bFELASER, bNoIronLaser)
        ENDIF

C***  ---------------------------------------------------------
C***  ADD THE THOMSON EMISSIVITY, ACCORDING TO THE ACTIVE MODE: 
C***  Note: CMFCOOP has stored the true emissivity into ETANOTH 
C***  ---------------------------------------------------------
C***    Using the Continuum Radiation Field 
        IF (BCOLIPP) THEN
          DO L=1, ND
            ETA(L) = ETANOTH(L) + OPA(L) * THOMSON(L) * XJCIND(L)
          ENDDO
        ENDIF

        LASER = .FALSE.

C***  Calculate frequency Step in X
C        DELTAX = FLOAT(K-KLAST) * DXMAX
        NDK = K - KLAST
        IF (NDK .NE. 1) THEN
           DELTAX = (EXP(FLOAT(NDK)*ALN) - 1.)*CLIGHT/VDOPUNIT
        ELSE
           DELTAX = DXMAX
        ENDIF

c         open (unit=120, file='k_coli.dat', status='unknown')
c         write (120,'(a,i8,f20.10,3i8)') 
c     >     'k=', k, xlamk, klast, ndk, iff_n

C***  Store NDK = K - KLAST in the special array IFF_DK
        IF (BFF_ACT) THEN
          IFF_DK(IFF_N) = NDK - 100
        ENDIF

        FWEIGHTL = DELTAX * VDOPUNIT * 1.E13 / XLAMK
C***  Sum up Opacities
        CALL ADDOPA (ND, NDDIM, MAXLIN, MAXIND, LIND, LINDS, 
     >         XK, XKMID, XKRED, DELTAX, FWEIGHTL,
     >         PARALAS, LASER, LASERV, ALN, VDOPUNIT, VDOPDD,
     >         WS, ETAL, OPAL, ETA, ETANOTH, OPA, ETAK, ETAKNOTH, 
     >         OPAK, OPAKNOTH, THOMSON, PWEIGHT, NOM,
     >         OPAFE, ETAFE, BFECHECK, BLASERL, NUP, LOW, N, LEVEL, 
     >         OPAKFE, NATOM, MAXATOM, KODAT, OPACELEM, OPAKELEM,     
     >         ETACELEM, ETAKELEM, MAXION, OPACION, OPAKION,
     >         ETACION, ETAKION, NCHARG, OPAFEION, ETAFEION, 
     >         OPAKNOFENOTH, ETAKNOFENOTH, 
     >         OPAFEFT, ETAFEFT, OPAKFEFT, ETAKFEFT)

        IF (LASER) THEN
          NLASER=NLASER+1
C!!!          WRITE (*,'(A,I7)') 'Laser condition at Index k=', k
        ENDIF

C***  INNER BOUNDARY VALUES
        CALL CLDIFFUS (XLAMK,T,RADIUS,ND,BCORE,DBDR,DTDR,TEFF,
     >                 OPAK,XHID)

C***  Beginning of Iteration Loop
      DO IT=ITSTART, ITMAX

        IF ((IT .GT. ITSTART)) THEN
           DO L=1, ND
              IF (XJLMOR2(L) .GT. 1.E-100) THEN
                 ETAK(L) = ETAKNOTH(L) 
     >                     + (OPAK(L) - OPAKNOTH(L)) * XJLMOR2(L)
              ELSEIF (XJCIND(L) > 1.E-50) THEN
                 ETAK(L) = ETAKNOTH(L)
     >                     + (OPAK(L) - OPAKNOTH(L)) * XJCIND(L)  
              ELSE
                 ETAK(L) = ETAKNOTH(L)
              ENDIF
           ENDDO
        ELSE
           DO L=1, ND
              IF (XJCIND(L) .GT. 1.E-100) THEN
                 ETAK(L) = ETAKNOTH(L)
     >                + (OPAK(L) - OPAKNOTH(L)) * XJCIND(L)
              ELSE
                 ETAK(L) = ETAKNOTH(L)
              ENDIF
           ENDDO
        ENDIF

C***    Calculate Source Function if Short Characteristiks is applied
C***      In the other Version S is calculated in COLISET
        IF (BSHORT_CHAR) THEN
          CALL CALC_S(S, ETAK, OPAK, ND)
        ENDIF
        
        ITACT = IT

C***   Beginning of COLIRAY-Block
        IF (IT .GT. 1) THEN
          IF ( K .EQ. 0) THEN
            WRITE (hCPR,'(A,I2)') 'COLIRAY, IT=', IT
          ENDIF
C***     Preset Moments of the Intensity to be calculated in COLIRAY
          CALL SET_MOMZERO(ND, XJL, XHL, XKL, XNL,
     >                         XHO, XHI, XNO, XNI,
     >                         XHOM, XHOP, XNOM, XNOP)

C***  Begin of SHORTRAY Block
C***  Loop over the impact parameters
            DO JP=1, NP
              IRAY=(JP-1)*ND+1
              CALL SHORTRAY(K, XIPLUS(1,JP), XIPLUS_OLD(1,JP,ITACT), 
     >          RADIUS, ND, NP, JP,
     >          S, S_OLD(1,ITACT), 
     >          XIMINUS(1,JP), XIMINUS_OLD(1,JP,ITACT),
     >          BCORE, DBDR, XHID, XLAMK, ENTOT, DELTAX,
     >          OPAK, ETAK, OPAK_OLD(1,ITACT), 
     >          XJL, XHL, XKL, XNL,
     >          CWM0, CWM1, CWM2, CWM3, 
     >          XHO, XHI, XNO, XNI, 
     >          Z(IRAY,1), PPP(IRAY,1), 
     >          IPLOT, BPLOT2, IW_COLIRAY_IPLUS, IVERS, 
     >          XHOM, XHOP, XNOM, XNOP, OPA, bWARNK)
            ENDDO
C***  For special integration in coliray needed
C!!!          XHI = XHI + BCORE/2. + DBDR/OPAK(ND)/3.

            XIPLUS_OLD(1:ND,1:NP,IT) = XIPLUS(1:ND,1:NP) 
            XIMINUS_OLD(1:ND,1:NP,IT) = XIMINUS(1:ND,1:NP) 

C***  Notation   :
C***       XJL   : R^2 * J from ray-by-ray  Weight : CWM0
C***       XHL   :  "    H  "       "       Weight : CWM1 at interstices
C***       XKL   :  "    K  "       "       Weight : CWM2
C***       XNL   :  "    N  "       "       Weight : CWM3 at interstices
C***     XJLMO   : R^2 * J from the moment equation
C***     XJLMOR2 : XJLMO divided by RADIUS^2
C***     XHLMO   :  "    H  "    "    "       "     (by differenciation of J)
C***  At the boundaries H and N are calculated from U at the radius points
C***       XHO   : R^2 * H from ray-by-ray  Weight : CWM1O at outer boundary
C***       XHI   : R^2 * H from ray-by-ray  Weight : CWM1I at inner boundary
C***       XNO   : R^2 * N from ray-by-ray  Weight : CWM3O at outer boundary
C***       XNI   : R^2 * N from ray-by-ray  Weight : CWM3I at inner boundary

C***  Calcultate the Eddis (EDDIF and EDDIG) and EPGSMAX
            DO L=1, ND
C***      COLIRAY (the last, if COLI++) recommends EPSGMAX for the next COLI+
              IF (IT .EQ. ITMAX .AND. ABS(XJL(L)) .GT. 1.E-300) THEN
                EPSGMAX(L) = AMAX1(EPSGMAX(L),
     >                      -(GLF2(L)/VLF2(L)*XNL(L)+XHL(L))/XJL(L))
                EPSGMAX(L) = AMAX1(EPSGMAX(L), -XHL(L)/XJL(L))
              ENDIF
              IF (XJL(L) .GT. 0.) THEN
                EDDIF(L) = XKL(L) / XJL(L)

              ELSE
                EDDIF(L) = 1.
              ENDIF
              IF (L .EQ. ND) CYCLE
              XNL_MID = XNL(L) + XNL(L+1)
              XHL_MID = XHL(L) + XHL(L+1)
              XJL_MID = XJL(L) + XJL(L+1)
              XNENN = XHL_MID + EPSG(L,IT)*XJL_MID
              
              IF (ABS(XNENN) .GT. EXP(-499.D0)) THEN
                EDDIG(L) = XNL_MID / XNENN
                 
C***            CARDS option: Reset all negative EDDIG to
C***                          ensure proper momentum equation form
                IF (bNoNEGEDDIG .AND. EDDIG(L) < 0.) THEN
                  EDDIG(L) = 1. / (1. + EPSG(L,IT)) 
                ENDIF
                
              ELSE
c***            Changed the fallback to 1 / (1+eps) -- ansander, Aug 2023
c                EDDIG(L) = 1.
                EDDIG(L) = 1. / (1. + EPSG(L,IT))              
              ENDIF
            ENDDO
C***  Calcultate the Eddis (EDDIH and EDDIN)
C***    at outer and inner boundary
            IF (XJL(1) .GT. 0.) THEN
              EDDIHOUT  = XHO  / XJL(1)
              EDDINOUT  = XNO  / XJL(1)
              EDDIHOUTP = XHOP / XJL(1)
              EDDINOUTP = XNOP / XJL(1)
            ELSE
              EDDIHOUT  = 1.
              EDDINOUT  = 1.
              EDDIHOUTP = 1.
              EDDINOUTP = 1.
            ENDIF
            IF (XJL(ND) .GT. 0.) THEN
C***          Attention: Since 09 Feb 2016  XHI is the special definition 
C***                     of H_spec at the inner boundary, calculated in SHORTRAY
C***                     Before XHI was identical to XHL(ND) from SHORTRAY
              EDDIHIN  = XHI / XJL(ND)      !special Eddington factor for inner boundary
              EDDININ  = XNI / XJL(ND)      !to be updated (but currently unused)
            ELSE
              EDDIHIN  = 0.5                !Default value has to be 0.5 due to H_spec
              EDDININ  = 1.
            ENDIF
C***   End of SHORTRAY Block
          ENDIF

          IF ( K .EQ. 0) THEN
            WRITE (hCPR,'(A,I2)') 'COLIMO,  IT=', IT
          ENDIF
C*** ISU unused variable
C         CALL COLIMO(K, ND, RADIUS, OPAK, ETAK, ETAKNOTH, 
          CALL COLIMO(K, ND, RADIUS, OPAK, ETAKNOTH, 
     >              OPAKNOTH, 
C    >              S, XJLMO, XJLMOR2, XHLMO,
     >              XJLMO, XJLMOR2, XHLMO,
     >              XJLMO_OLD(1,ITACT), XHLMO_OLD(1,ITACT), 
     >              DLF, DLH, GLF, GLH, VLF, VLH,
     >              GLF2, GLH2, VLF2, VLH2, 
     >              QLF, QLH, OPAKH,
     >              EDDIF, EDDIFO(1,ITACT), EDDIG, EDDIGO(1,ITACT), 
C    >              EDDIHOUT, EDDIHIN, EDDIHOUTO(ITACT), 
     >              EDDIHOUT, EDDIHIN,
C    >              EDDINOUT, EDDININ, EDDINOUTO(ITACT), 
     >              ALH, BLH, CLH,
     >              COLIA, COLIB, COLIC, COLIW, DELTAX, 
C    >              BCORE, DBDR, XIMINUS, BPDONE, XLAMK,
     >              BCORE, XLAMK,
     >              DJDSMOD, bALOTri, DJDSMOU, DJDSMOL,
     >              FULFIL0, FULFIL1, BPLOT, BPLOT2, IPLOT, 
     >              IW_COLIMO_F, IW_COLIMO_G, IW_COLIMO_G2, BSTATIC, 
     >              CLMOETA, CLMOOPA, XHID, 
     >              RADIUS2, EPSG(1,ITACT), GEPSB, GEPSBO,  
C    >              XHOM, XHOMO(ITACT), XNOM, XNOMO(ITACT),
C    >              EDDIHOUTP, EDDINOUTP, 
C    >              EDDIHOUTOP(ITACT), EDDINOUTOP(ITACT), IWARNJMNEG)
     >              IWARNJMNEG)

***    Calculate Integrals over Frequency
          IF (IT .EQ. ITMAX) THEN                
          
            IF (bPLOTRTAU1) THEN
              CALL PLOTRTAU1COLI (hPLOT, K, XLAMK, XLAMKOLD, CMODE, 
     >                            ND, RADIUS, OPAK, OPA, FILLFAC, 
     >                            MAINPRO, MAINLEV, MODHEAD)
            ENDIF
          
            IF ( K .EQ. 0) THEN
              WRITE (hCPR,'(A,I2)') 
     >          'Frequency Integration always from COLIMO'
            ENDIF
            CALL FREQUINT (K, FWEIGHTL, XLAMK, XLAMKOLD,
     >               XLAMBDA, DXMAX, VELO, GRADI,
C***        Take J and H from COLIMO
     >               XJLMO, XJLMO_OLD(1,ITMAX), 
     >               XHLMO, XHLMO_OLD(1,ITMAX),
     >               OPAKOLD, NF, ND, NDDIM, NATOM,
     >               KONTACT, KONTAUP, DFKONT, LIND, PWEIGHT, MAXLIN,
     >               INDFEACT, SIGMAACT, MAXFEACT, LASTFE,
     >               OPAK, ETAK, OPAFEI, ETAFEI, IFENUP, IFELOW,
C***                Take Diagonal Weights from COLIMO
     >               DJDSMOD, DJDSMOD_OLD, DJDSraw,
     >               bALOTri, DJDSMOU, DJDSMOL,
     >               WEIGHT, N, POPNUM,
     >               XKLOLD(1,ITMAX), XNLOLD(1,ITMAX),
     >               OPAKNOTH, ETAKNOTH, OPAKNOTHO, ETAKNOTHO,
     >               EDDIFO(1,ITMAX), EDDIGO(1,ITMAX),
     >               EPSGO(1,ITMAX), RADIUS, BCOLIRAY,
     >               ETANOTH, OPA, THOMSON, T, ETAKOLD,
     >               QLFOLD, QLHOLD, OPAKHOLD, 
     >               EDDIHOUTOLD, EDDIHINOLD, XHID,
     >               FWEIGHT, OPAO, THOMSONO,
     >               OPAKOLDELEM, OPACOLDELEM, ETAKOLDELEM, 
     >               ETACOLDELEM, OPAKOLDION, OPACOLDION,
     >               ETAKOLDION, ETACOLDION, OPAFE,
     >               SUMJ, SUMJW, SUMDJDSC, SUMDJDSCW,
     >               SUMOPA, SUMOPAW, SUMETA, SUMETAW, SC, SCO,
     >               XJCINT, FWTEST, XJLMEAN, XJFEMEAN, 
     >               XLAMAPPMEAN, HTOTL,
     >               HTOTMINUSND, HTOTND, HTOTNDCOR,
     >               HTOTNDS, ARAD, ACONT, ATHOM, ARADELEM, ACONTELEM, 
     >               ARADION, ACONTION,
     >               iIgnoreK, XJTOTL, XKTOTL, XNTOTL, WFELOW, WFENUP, 
     >               FTCOLI, FTCOLIB, DSDFELOW, DSDFENUP, WJC, WJC_MIN, 
     >               DSDSC, DJDSC, DJDSCOLD, OPAKINT, ETAKINT,
     >               DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
     >               OPASMEAN, OPASMEANTC, SMEAN, QFJMEAN, OPAJMEAN, 
     >               OPAJMEANTC, OPAPMEAN, QOPAHMEAN, HMEAN, 
     >               EDDIHOUTJMEAN, XHOM, HTOTOUTMINUS,
     >               RADIUS2, OPC, 
     >               FERATUL, FERATLU, ELEVEL, EMCOLI,
     >               FTFE, IVERS_FE_EXPFAC, LPLOT_WCHARM, GAMMACOLI,
     >               OPAROSS, OPAROSSCONT, OPALAMBDAMEAN, 
     >               OPAROSSELEM, OPATOTELEM,
     >               GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2, TEFF,
C*** die folgenden SKALAREN Parameter werden ausgereicht, weil sie sonst
C*** ueberschrieben werden!!!!???
     >               XNUEK, XNUEKOLD, XNUEKOLDOLD, VDOPUNIT,
     >               IFRBSTA, IFRBEND, VDOPFE, DXFE, XLAM0FE, NCHARG,
     >               MAXION, POPMIN, TAULAST, 
     >               XHLOLD(1,ITMAX), HTOTLTEST, iHTOTCUT, HTOTCUT,
     >               bOSKIPLAST, XLAMLASTOSKIP, XHBFSKIP, LPRDH, bDEBUG,
     >               CUTOPAMEANTHRES)     
     
          ENDIF

C***   Save old Quantities (Dimensioned with NIT) for next Frequency
          XJLOLD(1:ND,IT)      = XJL(1:ND)
          XHLOLD(1:ND,IT)      = XHL(1:ND)
          XKLOLD(1:ND,IT)      = XKL(1:ND)
          XNLOLD(1:ND,IT)      = XNL(1:ND)
          S_OLD(1:ND,IT)       = S(1:ND)
          OPAK_OLD(1:ND,IT)    = OPAK(1:ND)
          EDDIFO(1:ND,IT)      = EDDIF(1:ND)
          XJLMO_OLD(1:ND,IT)   = XJLMO(1:ND)
          EDDIGO(1:ND-1,IT)    = EDDIG(1:ND-1)
          EPSGO(1:ND-1,IT)     = EPSG(1:ND-1,IT)
          XHLMO_OLD(1:ND-1,IT) = XHLMO(1:ND-1)
          EDDIHOUTO(IT)        = EDDIHOUT
          EDDIHINO(IT)         = EDDIHIN
          EDDINOUTO(IT)        = EDDINOUT
          EDDININO(IT)         = EDDININ
          XHIO(IT)             = XHI
          XHOMO(IT)            = XHOM
          XNOMO(IT)            = XNOM
          EDDIHOUTOP(IT)       = EDDIHOUTP
          EDDINOUTOP(IT)       = EDDINOUTP

C***  End of NIT Iteration Loop
        ENDDO

C***  Store old quantities which are dimensioned without NIT
        XLAMKOLD = XLAMK
        EDDIHOUTOLD = EDDIHOUT
        EDDIHINOLD = EDDIHIN
        EDDIHINMOLD = EDDIHINM
        EDDIHINTOLD = EDDIHINT
        OPAKOLD(1:ND)    = OPAK(1:ND)
        ETAKOLD(1:ND)    = ETAK(1:ND)
        OPAKNOTHO(1:ND)  = OPAKNOTH(1:ND)
        ETAKNOTHO(1:ND)  = ETAKNOTH(1:ND)
        QLFOLD(1:ND)     = QLF(1:ND)
        DJDSMOD_OLD(1:ND) = DJDSMOD(1:ND)
        ETANOTHO(1:ND)   = ETANOTH(1:ND)
        OPAO(1:ND)       = OPA(1:ND)
        THOMSONO(1:ND)   = THOMSON(1:ND)
        QLHOLD(1:ND-1)   = QLH(1:ND-1)
        OPAKHOLD(1:ND-1) = OPAKH(1:ND-1)
        OPAKFEOLD(1:ND)   = OPAKFE(1:ND)
        OPAKFEFTOLD(1:ND) = OPAKFEFT(1:ND)
        ETAKFEFTOLD(1:ND) = ETAKFEFT(1:ND)
        OPAKNOFENOTHO(1:ND) = OPAKNOFENOTH(1:ND)
        ETAKNOFENOTHO(1:ND) = ETAKNOFENOTH(1:ND)
        CALL FREQUBAK(NATOM, ND, OPAKELEM, OPAKOLDELEM)
        CALL FREQUBAK(NATOM, ND, OPACELEM, OPACOLDELEM)
        CALL FREQUBAK(NATOM, ND, ETAKELEM, ETAKOLDELEM)
        CALL FREQUBAK(NATOM, ND, ETACELEM, ETACOLDELEM)
        CALL FREQUBAKION(ND, NATOM, MAXION, OPAKION, OPAKOLDION)
        CALL FREQUBAKION(ND, NATOM, MAXION, OPACION, OPACOLDION)
        CALL FREQUBAKION(ND, NATOM, MAXION, ETAKION, ETAKOLDION)
        CALL FREQUBAKION(ND, NATOM, MAXION, ETACION, ETACOLDION)
        

C***  Test-output at given depth for all wavelengths between XLP1 and XLP2
        IF (BPLOT) THEN
          IF(XLAMK .GE. XLP1 .AND. XLAMK .LE. XLP2) THEN
            SK1 = ETA(IPLOT)/OPA(IPLOT)
            SK2 = ETANOTH(IPLOT)/OPA(IPLOT)
            SL = ETAK(IPLOT)/OPAK(IPLOT) * RADIUS2(IPLOT)
            SLNOTH = ETAKNOTH(IPLOT)/OPAKNOTH(IPLOT)
            WRITE (37,'(I8,27(1X,E15.7))') 
     >        K, XLAMK, OPA(IPLOT), ETA(IPLOT), ETANOTH(IPLOT), 
     >        SK1, SK2, 
     >        OPAK(IPLOT), ETAK(IPLOT), SL, 
     >        XJCIND(IPLOT)*RIPL2, XJLO(IPLOT, IFRO), 
     >        XJL(IPLOT), 
     >        XHL(IPLOT), XKL(IPLOT), XNL(IPLOT), 
     >        EDDIF(IPLOT), EDDIG(IPLOT), 
     >        XJLMO(IPLOT), XJLMOR2(IPLOT), 
     >        XHLMO(IPLOT), 
     >        FULFIL0(IPLOT), FULFIL1(IPLOT), DJDSMOD(IPLOT),
     >        DJDS(IPLOT), SLNOTH
C            XHI, XHID
c            CALL PLOT3D (XLAMK, ND, INDFEACT, MAXFEACT, BPSTART, 
c     >                   LEVEL, IFELOW, IFENUP, OPAFEI, OPAKH, 
c     >                   RSTAR, RHO, ENTOT, XHLMO, XJLMO, XNUEK, 
c     >                   OPAKNOTH, ETAKNOTH, RADIUS2)
            IF (BPSTART) BPSTART = .FALSE.
          ENDIF
C***  Test-output at given wavelentgh over depth
          IF (.NOT. BPDONE) THEN
            IF (K .GE. LPLOT) THEN
              BPDONE = .TRUE.
              DO L=1, ND
                SK1 = ETA(L)/OPA(L)
                SK2 = ETANOTH(L)/OPA(L)
                SL = ETAK(L)/OPAK(L) * RADIUS2(L)
                WRITE (36,'(I6,25(1X,E15.7))') 
     >            L, XLAMK, OPA(L), ETA(L), ETANOTH(L), SK1, SK2, 
     >            OPAK(L), ETAK(L), SL, 
     >            XJCIND(L), XJLO(L, IFRO), 
     >            XJL(L), XHL(L), XKL(L), XNL(L), 
     >            EDDIF(L), EDDIG(L), XJLMO(L), XJLMOR2(L), XHLMO(L), 
     >            FULFIL0(L), FULFIL1(L), DJDSMOD(L), RADIUS(L)
              ENDDO
            ENDIF
          ENDIF                
        ENDIF

C***  Store Wcharm-factors here (Variable = DJDSMOD) in array IFF_WCHARM
        IF (BFF_ACT .AND. (.NOT. bKATEST)) THEN
          CALL CLSAVEWC (IFF_WCHARM, IFF_MAX, IFF_N, ND, DJDSMOD, 
     >                    K, XLAMK)
        ENDIF

C***  Store Maximum Opacity
        IF (OPAK(1) .GT. OPAMAX1) THEN
          OPAMAX1        = OPAK(1)
          OPAMAX1_LAMBDA = XLAMK
          IOPAMAX1_K     = K
        ENDIF
 
C***  Output of K
        IF (K .GT. KOUT) THEN
          WRITE (hCPR,'(A,I7,1X,I7,1X,F15.3)') 'K=', K, NK, XLAMK
          KOUT = KOUT + 20000
        ENDIF
 
C***  Increase K
        KLAST = K
        K = K + 1
        NK = NK + 1

C***  The next line is to skip the first frequency points for test reasons
C!!!        if (k .eq. 1) k = 100000

C***  Save old EDDIEs to file fort.<hEDDI>
        IF (BCOLIP .AND. (.NOT. bKATEST)) THEN 
          CALL CLSAVEE (hEDDI, NZE1, NZE2, NFRO, EDDIA, NDEDDIA,
     >                EDDIF, EDDIG, ND,
     >                EDDIHOUT, EDDIHIN, EDDIHINM, EDDINOUT, EDDININ,
     >                BCLEERR, CMODE, XHI, XHO, EPSG(1,ITMAX), 
     >                XHOM, XNOM, 
     >                EDDIHOUTP, EDDINOUTP, EDDIHINT)
        ENDIF

C***  Store XJLMOR2 in XJLO and write XJLO to file (if necessary)
c        CALL CLSAVEJ (hEDDI, ND, NFRO, XK, IFRO, NXJO, NXK,
c     >                      XJLO, XJLMOR2, CMODE)

C***  Last frequency has been finished : EXIT Main Loop
        IF (CMODE .EQ. 'F') EXIT

C***  End of Main-Loop over all frequencies
      ENDDO freqloop !- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      WRITE (hCPR,'(A,I7,1X,I7,1X,F15.3)') 'K=', K, NK, XLAMK
      IF (BCOLIP .AND. (.NOT. bKATEST)) THEN
        WRITE (hCPR,*) 'EDDI file has been updated...'
      ENDIF

C***  Take EPSGMAX at Interstices as maximum from the adjacent full points.
C***  Multiply with a FACTOR for a save estimate. Write on model-file.
      IF (BCOLIP) THEN
        DO L=1, ND-1
          EPSGMAX(L) = 2. * AMAX1(EPSGMAX(L), EPSGMAX(L+1)) + 0.1
C***  Switching down of EPSG in small steps, not more than -20%
          EPSGMAX(L) = AMAX1 (EPSGMAX(L), 0.8 * EPSG(L,ITMAX))
        ENDDO
      ENDIF

      CALL FREQUNORM (ND, MAXION, OPASMEAN, OPASMEANTC, SMEAN, 
     >                QFJMEAN, XJTOTL, OPAJMEAN, OPAJMEANTC, OPAPMEAN, 
     >                QOPAHMEAN, HMEAN, EDDIHOUTJMEAN,
     >                RADIUS, RSTAR, DENSCON, 
     >                  FTCOLI, FTCOLIB, WJC, WJC_MIN, 
     >                  FWTEST, NF, OPC, FTFE, LASTFE, 
     >                  LPLOT_WCHARM, XLAMBDA, OPAROSS, OPALAMBDAMEAN, 
     >                  ARAD, ACONT, ATHOM, ENTOT, ABXYZ, ATMASS, NATOM,
     >                  T, GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2,TAUROSS,
     >                  ARADELEM, ACONTELEM, ARADION, ACONTION,
     >                  OPAROSSCONT, OPAROSSELEM, OPATOTELEM)

C***  Output of Maximum Opacity
      WRITE (hCPR,*)
      WRITE (hCPR,'(A, A,I6, A,F10.3, A,E10.3, A,F10.3)')
     >  'Maximum Opacity at Depth 1: ', 
     >  'K=', IOPAMAX1_K, 
     >  ';  Lambda=', OPAMAX1_LAMBDA, 
     >  ';  Opacity=', OPAMAX1, 
     >  ';  Delta-Tau (Opa*Delta_R)=', OPAMAX1 * (RADIUS(1)-RADIUS(2))
      WRITE (hCPR,*)
      WRITE (*,'(A, A,I6, A,F12.3, A,E12.5, A,F12.5)')
     >  'COLI> Maximum Opacity at Depth 1: ', 
     >  'K=', IOPAMAX1_K, 
     >  ';  Lambda=', OPAMAX1_LAMBDA, 
     >  ';  Opacity=', OPAMAX1, 
     >  ';  Delta-Tau (Opa*Delta_R)=', OPAMAX1 * (RADIUS(1)-RADIUS(2))

C***  Output of Warnings
      IF (IW_COLIMO_F > 0) THEN 
        WRITE (hCPR,'(A,I7,1X,I6)') 
     >    'EDDIF resets in COLIMO: JOBNUM, All: ', 
     >    JOBNUM, IW_COLIMO_F
      ENDIF
      IF (IW_COLIMO_G > 0) THEN
        WRITE (hCPR,'(A,I7,1X,I6,A,I6)') 
     >    'EDDIG resets in COLIMO: JOBNUM, All / Relevant: ', 
     >    JOBNUM, IW_COLIMO_G, ' / ', IW_COLIMO_G2

      ENDIF
      IF (IW_COLIMO_J .GT. 0) THEN
        WRITE (hCPR,'(A,I8)') 'XJLMO not used : ', IW_COLIMO_J
      ENDIF
      IF (IWARNJMNEG > 0) THEN
        WRITE (hCPR,'(A,I8)') 'XJLMO set to zero : ', IWARNJMNEG
      ENDIF
      IF (IW_COLIRAY_IPLUS .GT. 0) THEN
        WRITE (hCPR,'(A,I6,A)') 
     >    'XIPLUS set to zero for', IW_COLIRAY_IPLUS, ' Rays'
      ENDIF
      IF (IW_COLIRAY_U .GT. 0) THEN
        WRITE (hCPR,'(A,I6,A)') 
     >    'U not accounted for at ', IW_COLIRAY_U, ' Rays and Depths'
      ENDIF
      IF (IW_COLIRAY_IPLUS .GT. 0) THEN
        WRITE (hCPR,'(A,I6,A)') 
     >    'XIPLUS set to zero for', IW_COLIRAY_IPLUS, ' Rays'
      ENDIF

C**** Debug option: print out the number of frequency points with J < 0 
      IF (bDEBUG) THEN 
        DO L=1, ND
          WRITE (hCPR, '(A,I3,A,I8)') 'Frequencies ignored'
     >        // ' in FREQUINT at L = ', L,': ', iIgnoreK(L)
        ENDDO
      ENDIF 
       
      bNewLoop = .FALSE.

      !Preparations for the Calculation of Force-Multipliers
      ! (in this case COLI is passed two times, first run
      !    is a test run with modified velocity and v gradient)
      IF (bKALPHA .AND. bKATEST) THEN
         !If force multipliers are calculated (usually for STEAL->HYDROSOLVE)
         ! then save results of test run and force a second run
         IF (bKADONE) THEN
            bForceCOLIP = BCOLIP
            bForceCLEERR= BCLEERR
            bKATEST     = .FALSE.
            IF (iTypeAlpha == 3 .OR. iTypeAlpha == 31) THEN
              ARADTEST    = ARAD
            ELSE
              ARADTEST    = ARAD - ACONT
            ENDIF
            ALINETEST   = ARAD - ACONT
            ACONTTEST   = ACONT
            ATHOMTEST   = ATHOM
            HTOTTEST    = HTOTL
         ENDIF
         bNewLoop = .TRUE. 

C***  Calculation of Force-Multipliers:
      ! Following Eq. (4) of GH2005 and ALINES = ARAD - ACONT
      ! alpha can be calulated from ALINES = C/r^2 (dv/dt)^(ALPHAF)
      ! To avoid calculating C COLI is run two times, first with 
      ! 0.9 * v and 0.9 * v', then with usual v, v'. 
      ! From the first run, ARADTEST is calculated and then the
      ! ratio with the result from the second one is taken. As ALPHAF
      ! should be the same in both runs, we can directly obtain ALPHAF
      ! from the relation 0.9^(ALPHAF) = DELTA
      ! Instead of taking a simple ratio, ALINES is first divided by
      ! HTOT for flux normalization (GH2005)
C***  ALPHA TYPES: (can be set on the HYDRO cards)
C      0    - no calculation of ALPHAF or Response Factor (switches off double-COLI)
C      1,10 - calculation of ALPHAF via lamda (modifies v and dv/dr)
C       11  -   dito but only v is modified
C       12  -   dito but only dv/dr is modified
C      2    - calculation of ALPHAF via lamb**2 (modified v and dv/dr)
C      3    - calculation of Response Factor (ARMDRESP)
C      31   - calculation of Response Factor (ARMDRESP) incl. velo modificiation
      ELSE IF (bKALPHA .AND. bKADONE) THEN
         bPrintAlpha = .TRUE.
      
         !Debug output:
C         WRITE (hCPR,*)
C         WRITE (hCPR,*) "Debug output - Force multipliers"
         WRITE (hCPR,*) 
         WRITE (hCPR,*) ' Alpha-Type:', iTypeAlpha
         WRITE (hCPR,*) ' XMG = ', XMG
         DO L=1, ND
            RHO(L) = AMU*ATMEAN*ENTOT(L)
            XMU(L) = ATMEAN/(1.+RNE(L))
            VELORAW(L) = VELO(L) * VDOPUNIT
            VMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !v_mach in km/s
         ENDDO
         DOTM4P = RHO(1)*VELO(1)*VDOPUNIT*1.E5*(RADIUS(1)*RSTAR)**2
         ARADTOT = 0.
         ACONTTOT = 0.
         ARADTESTTOT = 0.
         ACONTTESTTOT = 0.
         DO L=1, ND-1
           DELTAR = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
           ARADTOT = ARADTOT + ARAD(L) * DELTAR
           ACONTTOT = ACONTTOT + ACONT(L) * DELTAR
           ARADTESTTOT = ARADTESTTOT + ARADTEST(L) * DELTAR
           ACONTTESTTOT = ACONTTESTTOT + ACONTTEST(L) * DELTAR
           DELTA = 
     >           (HTOTL(L)/HTOTTEST(L))*(ARADTEST(L)/(ARAD(L)-ACONT(L)))
            !Additional multiplication with VELOMODFAK from inferred density variation
C            DELTA = VELOMODFAK *
C     >           (HTOTL(L)/HTOTTEST(L))*(ARADTEST(L)/(ARAD(L)-ACONT(L)))
           ARMDRESP(L) = 1.
           
           IF (iTypeAlpha == 3) THEN
               ALPHAF(L) = 0.
               ACONTMDRESP(L) = ACONT(L)*VELOMODFAK/ACONTTEST(L)
               ALINEMDRESP(L) = (ARAD(L)-ACONT(L))
     >                                 * VELOMODFAK/ALINETEST(L)
               ARMDRESP(L) = ARAD(L)*VELOMODFAK/ARADTEST(L)
C               ARMDRESP(L) = ARAD(L)*VELOMODFAK/ARADTEST(L)
     >                           * HTOTTEST(L)/HTOTL(L)         
CC          scale to a 10 percent change
               ARMDRESP(L) = 1. - (1. - ARMDRESP(L)) 
     >                                    * 0.1 / (1./VELOMODFAK - 1.)
               WRITE (hCPR,'(A,I4,6(2X,G15.8))') ' Mdmod ', L, 
     >            ARAD(L)/XMG * (RADIUS(L)*RSTAR)**(2.),
     >            ARADTEST(L)/XMG * (RADIUS(L)*RSTAR)**(2.),
     >            ARADTEST(L)/ARAD(L), 
     >            ARMDRESP(L),      !Response of a_rad on Mdot change
     >            VELOMODFAK/ACONTMDRESP(L), VELOMODFAK/ALINEMDRESP(L)
           ELSEIF (DELTA > 0.) THEN               
             IF (iTypeAlpha == 2) THEN          
               ALPHAF(L) = ALOG10(DELTA*ATHOM(L)/ATHOMTEST(L))
     >                               /ALOG10(VELOMODFAK)/2.
               CK(L) = (ARAD(L)-ACONT(L) -ARADTEST(L))
     >            *(RADIUS(L)*RSTAR)**2 * RHO(L)/(1.-VELOMODFAK**2)
     >             / (GRADI(L)*VDOPUNIT*1.E5/RSTAR)
               ALPHAALT(L) = 
     >             LOG10(ARADTEST(L)/(ARAD(L)-ACONT(L))) 
     >                            / LOG10(VELOMODFAK)
             ELSE
               ALPHAF(L) = ALOG10(DELTA)/ALOG10(VELOMODFAK)
               CK(L) = (ARAD(L)-ACONT(L) -ARADTEST(L))
     >            *(RADIUS(L)*RSTAR)**2 * RHO(L) / (1 - VELOMODFAK)
     >             / (GRADI(L)*VDOPUNIT*1.E5/RSTAR)
             ENDIF
             GAMMATHIN(L) = 1./XMG*((RADIUS(L)*RSTAR)**(2.)*
     >                                             (ARAD(L)-ACONT(L))
     >          - CK(L)/RHO(L) * GRADI(L)*VDOPUNIT*1.E5/RSTAR)
             GAMMAPOWR = (RADIUS(L)*RSTAR)**(2.) /XMG * ACONT(L) 
     >                       + GAMMATHIN(L)
             IF (iTypeAlpha <= 2 .AND. bPrintAlpha) THEN
               ALPHAALT(L) = 
     >             LOG10(ARADTEST(L)/(ARAD(L)-ACONT(L))) 
     >                            / LOG10(VELOMODFAK)
               WRITE (hCPR,'(A,I4,10(2X,G15.8))') ' Thomson-corr: ',
     >           L , ATHOMTEST(L)/ATHOM(L), ALPHAF(L),  
     >           (CK(L)/(RADIUS(L)*RSTAR)**2 
     >                      /RHO(L)*GRADI(L)*VDOPUNIT*1.E5/RSTAR)/
     >           ( ARAD(L) - ACONT(L) ),
C     >           ( XMG / (RADIUS(L)*RSTAR)**2 *GAMMATHIN(L) ),
     >           CK(L), GAMMATHIN(L), GAMMAPOWR, 
     >            (RADIUS(L)*RSTAR)**2 /XMG * (ARAD(L) * (1.-ALPHAF(L))
     >             + ALPHAF(L) * ACONT(L)),
     >            CK(L)/DOTM4P, ALPHAF(L)*(ARAD(L)-ACONT(L))/
     >               (VELO(L)*GRADI(L)/RSTAR*(VDOPUNIT*1.E5)**2),
     >           1. - GAMMAPOWR
             ENDIF
           ENDIF

         ENDDO
         WRITE (hCPR,'(A,8(2X,G15.8))') 'TOT: ', ARADTESTTOT/ARADTOT, 
     >        ARADTOT/ARADTESTTOT, VELOMODFAK, 1./VELOMODFAK,
     >       ARADTESTTOT/ARADTOT/VELOMODFAK, 
     >       ARADTESTTOT/ARADTOT*VELOMODFAK,
     >       ACONTTESTTOT/ACONTTOT, 
     >       (ARADTESTTOT-ACONTTESTTOT)/(ARADTOT-ACONTTOT)
         WRITE (hCPR,*) 
         IF (BPLOTALPHA) THEN         
           CALL PLOTALPHA(ALPOPT, ND, RADIUS, ALPHAF,
     >                    VMACH, VELORAW, TAUROSS, ENTOT, ATMEAN, 
     >                    MODHEAD, JOBNUM, .TRUE.)
         ENDIF
      ENDIF


C***  Setup FF_INFO
      CALL SETUP_FF(FF_INFO, XLAM0, ALN,
     >              XLAM_FINE_START, XLAM_FINE_END, IFF_N, KSPACE, 
     >              IFF_DK, IFF_MAX)

c      open (unit=120, file='iff_dk_coli.dat', status='unknown')
c      write (120,'(i8)') (iff_dk(ii), ii=1, ff_info(7))

C***  Computation of total flux, TOTOUT, from the line-blanketed flux
C***  This overwrites on the MODEL file TOTOUT from WRCONT 
C***  -- new by wrh 28-Jan-2011 
      TOTOUT = .0
      DO K=1, NF
         TOTOUT = TOTOUT + EMCOLI(K) * FWEIGHT(K)
      ENDDO

      IF (.NOT. bNewLoop) THEN
        DO K=1, NF
          IF (LSOPA > 0) THEN
            IF (K == 1) THEN
              WRITE (hCPR,*) 'COLI> Opacities exported to OUT file'
            ENDIF
            CALL PRIOPACOLI (XLAMBDA(K), K, ND, LSOPA, RADIUS, 
     >                       OPAKINT, ETAKINT, THOMSON, FILLFAC, 
     >                       IWARN, MAINPRO, MAINLEV, JOBNUM, MODHEAD)
          ENDIF      
        ENDDO
            
        CALL WMODCOLI(XJCINT, FWTEST, ARAD, ACONT, ATHOM, ND, NF,
     >              RADIUS, ENTOT, RSTAR,
     >              XJFEMEAN, SIGMAINT, LASTFE, HTOTL, WFELOW,
     >              WFENUP, FTCOLI, FTCOLIB, REDISMODE, NCOLIP, 
     >              XJTOTL, XKTOTL, XNTOTL, WJC,
     >              DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
     >              OPASMEAN, OPASMEANTC, OPAPMEAN, 
     >              QFJMEAN, SMEAN,
     >              OPAJMEAN, OPAJMEANTC, QOPAHMEAN,
     >              EDDIHOUTJMEAN, HMEAN,
     >              HTOTOUTMINUS, HTOTMINUSND, HTOTND, 
     >              HTOTNDCOR, HTOTCUT,
     >              LASTIND, FERATUL, FERATLU, BCOLIP, EPSGMAX,
     >              FTFE, EMCOLI, FF_INFO, IFF_DK, IFF_MAX,
     >              IFF_WCHARM, OPALAMBDAMEAN, TOTOUT, bKALPHA,
     >              ALPHAF, ARMDRESP, GAMMATHIN, CK,
     >              RHO, XMU, TAUROSS, OPAROSS, OPAROSSELEM, 
     >              OPAROSSCONT, OPATOTELEM, NATOM, 
     >              ARADELEM, ACONTELEM, MAXION, ARADION, ACONTION,
     >              hMODEL, hALO)
        CALL STORE_FF(MODHEAD, JOBNUM, ND, NDDIM, 
     >                FF_INFO, IFF_DK, IFF_MAX, IFF_WCHARM)
      ENDIF

C***  Laserversion 2 : same as 1, but the laserlines are additionally
C***    listed in the file LASER_LINES. This causes in STEAL that those lines
C***    are not included in the Scharmer amplification.
      IF (LASERV .EQ. 2) THEN
C***    OPEN FILE "LASER_LINES" FOR TRANSFER OF INFORMATION ABOUT LASER LINES
C***    TO MAIN PROGRAM "STEAL"
        OPEN (11,FILE='LASER_LINES', STATUS='UNKNOWN')
        DO NL=1, LASTIND
          IF (BLASERL(NL)) WRITE (11,'(A4,I6,2X,3A)') 
     >        'LINE', NL,LEVEL(INDNUP(NL)), ' - ', LEVEL(INDLOW(NL)) 
        ENDDO
      CLOSE(11)
      ENDIF

      IF (NLASER .GT. 0) THEN
        WRITE (*,'(A,F4.2,A,I6,A)') 'COLI> Opacities set to ', PARALAS, 
     >    ' * Background Opacity at ', NLASER, ' frequencies'
      ENDIF




C***  Write integrated moments of the intensity into file
      IF (BPLOT) THEN
        DO L=1, ND
          WRITE (35,'(I4,5(1X,E15.7))')
     >      L, XJTOTL(L), HTOTL(L), XKTOTL(L), XNTOTL(L), RADIUS(L)
        ENDDO
      ENDIF

C***  Entry Point for Short Coli
C***    This is necessary when Model is started with Temperature 
C***    Corrections. The first COLIs and STEALs are omitted
  100 CONTINUE

C***  UPDATING THE MODEL HISTORY
      CALL COLIHIST (MODHIST, MAXHIST, JOBNUM, BCOLIP, 
     >               ITSTART, ITMAX, BSHORT, IVERS, 
     >               BEMIX, BEPSGMAXERR, EMIXSTART, 
     >               bNewLoop .AND. bKALPHA)
 
      IF (BELIFI) THEN
        CALL CLOSMS (7, IERR)
      ENDIF
      CALL CLOSMS (hMODEL, IERR)
      CALL CLOSMS (hALO, IERR)
      IF ((.NOT. BSHORT) .AND. (.NOT. bNewLoop)) THEN
        CALL CLOSMS (hEDDI, IERR)
        WRITE (hCPR,*) 'Closing EDDI file...'
      ENDIF

      IF (BPLOT) THEN
        CLOSE (40)
        CLOSE (39)
        CLOSE (38)
        CLOSE (37)
        CLOSE (36)
        CLOSE (35)
        CLOSE (90)
        CLOSE (91)
        CLOSE (92)
        CLOSE (121)
        CLOSE (122)
        CLOSE (123)
        CLOSE (131)
        CLOSE (132)
        CLOSE (133)
      ENDIF

      IF (LPRDH > 0) THEN
        CLOSE (140)
      ENDIF
      
C***  Close WCHARM Plot File
      IF (LPLOT_WCHARM .GT. 0) THEN
        CLOSE (105)
      ENDIF      

      !write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)

C***  Multiple Coli
      IF (bNewLoop) THEN
        WRITE (hCPR,'(A)') 
     >      'COLI> New loop forced: COLI is run another time'
        GOTO 1
      ENDIF

C!!!  END OF COLI
      IF (CLHLP) THEN
        WRITE (hCPR,*) 'WARNING : This was Output Only'
        WRITE (hCPR,*) 'Model was not stored'
      ELSE
        CALL JSYMSET ('G0', '0')
      ENDIF

      CALL STAMP (OPSYS, 'COLI', TIM1)

      STOP 'O.K.'

      RETURN
      END
