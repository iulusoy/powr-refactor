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
