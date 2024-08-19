C***  MAIN PROGRAM COMO  *******************************************************
      SUBROUTINE COMO
 
C*******************************************************************************
C***  CONTINUOUS RADIATION TRANSFER (MOMENT EQUATIONS) WITH GIVEN EDDI-FACTORS
C***  -- FORMAL SOLUTION FROM GIVEN POP NUMBERS
C***  NOTE: ETLA METHOD NOW DISABLED !
C*******************************************************************************
 
      IMPLICIT NONE
 
C***  SET ARRAY DIMENSION PARAMETERS
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM  =          26 
      INTEGER, PARAMETER :: NDIM     =        2560 
      INTEGER, PARAMETER :: NFDIM    = 2*NDIM + 400 
      INTEGER, PARAMETER :: MAXIND   =       45000 
      INTEGER, PARAMETER :: MAXFEIND =        1500 
      INTEGER, PARAMETER :: MAXKONT  =     NFDIM/2 
      INTEGER, PARAMETER :: NDDIM    =          89 
      INTEGER, PARAMETER :: MAXHIST  =        4000 
      INTEGER, PARAMETER :: MAXXDAT  =          10 
       
 
C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO  = 2850 

      REAL, DIMENSION(NFDIM) :: XLAMBDA
      INTEGER :: NF, ND, N, LASTKON, LASTIND, LAST, LASTFE

      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(NDDIM) :: IWARN
      INTEGER, DIMENSION(NFDIM) :: KEY
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, EN
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      REAL, DIMENSION(NDDIM) :: OPA, ETA, THOMSON, RADIUS, ENTOT, T,
     >                          RNE, HTOTC, HNUE, A, B, C, W
      REAL, DIMENSION(NFDIM) :: FWEIGHT, XDUMMY, HEDDI
      
      REAL, DIMENSION(NFDIM, MAXKONT) :: SIGMAKI
      REAL, DIMENSION(NDDIM, NFDIM) :: XJC
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM
      
      REAL, DIMENSION(3,NDDIM) :: EDDI
      REAL, DIMENSION(4,NDIM) :: ALTESUM
      REAL, DIMENSION(4,MAXIND) :: COCO
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      
      CHARACTER(8*MAXHIST) :: MODHIST

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      COMMON /IRON/ FEDUMMY, INDRB, INDRF, SIGMAFE,
     >              IFRBSTA, IFRBEND, IFENUP, IFELOW, SIGMAINT
      LOGICAL BFEMODEL, BPLOTRTAU1

      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK

      CHARACTER(255) :: HISTENTRY
      CHARACTER(100) :: MODHEAD
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8) :: NAME
      CHARACTER(24) :: BUFFER24
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

      INTEGER, PARAMETER :: NPLOTOPADIM  =  20
      CHARACTER(80), DIMENSION(NPLOTOPADIM) :: OPTIONPLOTOPA

      REAL :: RSTAR, RMAX, TAUBMAX, DTDR, OPARND, TEFF, KV, TAUB,
     >        XLAM0FE, DXFE, VDOPFE, OPAMAX, RMAXGOOD, DUMMY
      INTEGER L, K, KBWARN1, KBWARN2, KBWARN3, IWARNJN, 
     >        NPLOTOPA, LPLOHTOT, LSHTOT, LSINT, LSOPA,
     >        JOBNUM, NATOM, NAUTO, KBMAX, IDUMMY, IERR, 
     >        N_WITH_DRLEVELS
      
      LOGICAL :: NOTEMP, BUNLU, BKUDRITZKI 

C***  Depth dependent Clumping
c***  nach goetz
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

      INTEGER, EXTERNAL :: IDX

C***  Array indicating POPMIN levels (flagged by steal)
      LOGICAL, DIMENSION(NDIM,NDDIM) :: ZERO_RATES
      
C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS
      CHARACTER(10) :: TIM1

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !MODEL file (fort.3)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE 
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
C***  Write Link Data (Program Version) to CPR file
      WRITE (hCPR,'(2A)') '>>> COMO started: Program Version from '
     >                 ,LINK_DATE
      WRITE (hCPR,'(4A)') '>>> created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL

c      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
c        CALL CLOCK(TIM1)
c      ELSE
c        CALL TIME(TIM1)
c      ENDIF

C***  READING THE ATOMIC DATA FROM FILE DATOM
      CALL DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $            EINST, ALPHA, SEXPO, 
     $            ADDCON1, ADDCON2, ADDCON3, 
     $            IGAUNT, COCO, KEYCBB, ALTESUM,
     $            INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $            ELEMENT, SYMBOL, NOM, KODAT, ATMASS, STAGE, 
     $            SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $            NLAST, NAUTO, MAXAUTO, LOWAUTO, WAUTO, EAUTO, AAUTO,
     $            IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $            IONGRND, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'COMO', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS)

C***  READING OF THE MODEL FILE ****************************************
      CALL       REMOCO (ND, NDDIM, RADIUS, ENTOT, T, RNE, NF, NFDIM,
     $                  XLAMBDA, FWEIGHT, KEY, XJC, 
     $                  N, POPNUM, RSTAR,
     $                  MODHEAD, LAST, MAXHIST, MODHIST, JOBNUM,
     $                  NOTEMP, TEFF, HEDDI, EDDI, MAXXDAT, XDATA, 
     >                  DENSCON, FILLFAC, OPARND, ABXYZ, NATOM, 
     >                  ZERO_RATES, NDIM)
      WRITE(hCPR,'(A,I7)') '>>> This is job number ', JOBNUM

      CALL POPMIN_NULLING (ZERO_RATES, POPNUM, ND, N)

C***  DECODING INPUT OPTIONS *******************************************
      CALL DECOMO (LSOPA, LSINT, NOTEMP, LSHTOT, LPLOHTOT, MODHIST, 
     >       BUNLU, BPLOTRTAU1, NPLOTOPA, OPTIONPLOTOPA, NPLOTOPADIM)

C***  logical NOTEMP is obsolete, as it is always .FALSE.
      NOTEMP = .FALSE.

C***  PRECALCULATION OF THE BOUND-FREE CROSS SECTIONS SIGMAKI
      CALL       BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $                    XLAMBDA,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,
     $                    KONTNUP,KONTLOW,LASTKON)
 
C***  CALCULATE TEMPERATURE GRADIENT "DTDR" AT INNER BOUNDARY
      CALL DIFDTDR(DTDR,TEFF,XJC,HEDDI,T(ND),RADIUS(ND),ND,EN,POPNUM,
     $             RNE(ND),ENTOT(ND),RSTAR,NDIM,N,LEVEL,NCHARG,WEIGHT,
     $             ELEVEL,EION,EINST,ALPHA,SEXPO,
     $             ADDCON1, ADDCON2, ADDCON3, 
     $             IGAUNT,NOM,NF,
     $             XLAMBDA,FWEIGHT,
     $             MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     $             KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC, 
     >             BKUDRITZKI,OPARND)
 
      IWARNJN = 0
      KBWARN1 = 0
      KBWARN2 = 0
      KBWARN3 = 0
      TAUBMAX = .0
      RMAX = RADIUS(1)

      DO L=1, ND
        HTOTC(L)=0.
      ENDDO

C***  From here: ENTOT scaled to clum density 
      DO L=1, ND
         ENTOT(L) = ENTOT(L) * DENSCON(L)
      ENDDO

C***  LOOP OVER ALL FREQUENCY POINTS  ----------------------------------
C***  SOLUTION OF THE MOMENT EQUATION AT EACH FREQUENCY POINT (FORMAL SOLUTION)
      DO K=1, NF

        KV=1+ND*(K-1)
 
        CALL COOP (XLAMBDA(K),ND,T,RNE,POPNUM,ENTOT,RSTAR,
     >             OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     >             NOM,KODAT,NDIM,N,MAXATOM,
     >             LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     >             ALPHA, SEXPO,
     >             ADDCON1, ADDCON2, ADDCON3,
     >             IGAUNT, SIGMATHK, SEXPOK, EDGEK,
     >             K,NF,SIGMAKI,RADIUS,
     >             KONTNUP,KONTLOW,LASTKON,XDATA)

C***    Scale down OPA, ETA for filling factor
        DO L=1, ND
          OPA(L) = OPA(L) * FILLFAC(L)
          ETA(L) = ETA(L) * FILLFAC(L)
        ENDDO

        IF (K <= 999) THEN
          WRITE (UNIT=NAME, FMT='(A4,I3)') 'EDDI', K
        ELSE
          WRITE (UNIT=NAME, FMT='(A4,I4)') 'EDDI', K
        ENDIF                
        CALL READMS (hMODEL,EDDI,3*ND,NAME, IERR)
        IF (IERR == -10) THEN
          !No EDDI entries found in model 
          ! maybe this is a Goetz model (EDD instead of EDDI)
          WRITE (UNIT=NAME, FMT='(A3, I4, A1)') 'EDD', K, ' '
          CALL READMS (hMODEL,EDDI,3*ND,NAME, IERR)
          IF (IERR == -10) THEN
            !If no EDDI values are there MOMO will crash => STOP here
            WRITE (hCPR,*) '*** MODEL DOES NOT CONTAIN EDDI[K] VALUES'
            WRITE (hCPR,*) '*** Try starting with COLI or new WRSTART'
            STOP 'FATAL ERROR in COMO'
          ENDIF
        ENDIF
 
        CALL MOMO (OPA, ETA, THOMSON, EDDI, RADIUS, XJC(KV,1),A,B,C,W,
     >             ND, XLAMBDA(K), T, DTDR, TEFF, NOTEMP, HNUE)
 
        IF (LSOPA > 0) THEN
          CALL PRIOPA (XLAMBDA(K),K,ND,LSOPA,RADIUS,OPA, ETA,
     >                 THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)
        ENDIF
      
        IF (BPLOTRTAU1) THEN
          CALL PLOTRTAU1 (NF, XLAMBDA, K, ND, RADIUS, OPA, MAINPRO,
     >                    MAINLEV, MODHEAD)
        ENDIF

C***    UPDATE OF ARRAY "EDDI" BECAUSE OF IMPROVED H+, WHICH ENTERS
C***    THE BOUNDARY CONDITION (KUDRITZKI VERSION)
        CALL WRITMS (hMODEL,EDDI,3*ND,NAME,-1, IDUMMY, IERR)

C***    FREQUENCY INTEGRATION OF FLUX 
        DO L=2,ND
          HTOTC(L) = HTOTC(L) + HNUE(L) * FWEIGHT(K)
        ENDDO

C***    NEGATIVE XJC'S ARE SET TO ABSOLUTE VALUE
        DO L=0, ND-1
          IF (XJC(KV+L,1) .LT. 0.) THEN
            IWARNJN = IWARNJN + 1
            XJC(KV+L,1) = ABS(XJC(KV+L,1))
          ENDIF
        ENDDO

C***    WARNING FOR OPTICALLY THICK CONTINUUM AT BOUNDARY
        IF (XLAMBDA(K) .LT. 1.E5) THEN
          TAUB = OPA(1) * RMAX
          IF (TAUB .GT. TAUBMAX) THEN
              TAUBMAX = TAUB
              KBMAX = K
          ENDIF
          IF ((KBWARN1 == 0) .AND. (TAUB > 1.)) KBWARN1 = K 
          IF ((KBWARN1 /= 0) .AND. (KBWARN2 == 0) .AND. (TAUB <= 1.)) 
     >                                         KBWARN2 = K-1
          IF ((KBWARN1 /= 0) .AND. (KBWARN2 /= 0) .AND. (KBWARN3 == 0)
     >                     .AND. (TAUB > 1.))  KBWARN3 = K 
        ENDIF

C***    UPDATING THE CONTINUOUS RADIATION FIELD ON THE MODEL FILE
        WRITE (UNIT=NAME, FMT='(A3,I4)') 'XJC',K
        CALL WRITMS(hMODEL,XJC(KV,1),ND,NAME,-1, IDUMMY, IERR)

      ENDDO
C***  ENDLOOP  ---------------------------------------------------------

      CALL WRITMS (hMODEL, HTOTC(2), ND-1, 'HTOTC   ',-1,IDUMMY,IERR)

C***  UPDATING THE MODEL HISTORY
      WRITE(UNIT=BUFFER24, FMT=20) JOBNUM, 'TEMP    '
   20 FORMAT ('/',I7,'. COMO  ',A8)

      CALL ADDHISTENTRY(MODHIST, -1, MAXHIST, 24, BUFFER24)
      CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
 
      CALL CLOSMS (hMODEL, IERR)
      CALL JSYMSET ('G0','0')

C***  OUTPUT *****************
      IF (IWARNJN .GT. 0) 
     >  WRITE(hOUT,*) 'COMO> WARNING: ',IWARNJN,' negative Continuum ',
     >                'Intensities are set to absolute values'

C***  WARNING FOR OPTICALLY THICK CONTINUUM AT BOUNDARY
      IF (TAUBMAX .GT. 1.) THEN 
        RMAXGOOD = TAUBMAX * RMAX * 3. 
        OPAMAX   = TAUBMAX / RMAX

        WRITE (hOUT,6) XLAMBDA(KBWARN1),XLAMBDA(KBWARN2),KBWARN1,KBWARN2
    6   FORMAT 
     >  (' COMO> *** WARNING: CONTINUUM OPTICALLY THICK AT OUTER BOUNDARY'
     >  ' BETWEEN', F7.1, ' AND', F7.1, ' ANG  (K=', I4, ' ...', I4, ')')

        IF (KBWARN3 .NE. 0) WRITE (*,7) XLAMBDA(KBWARN3), KBWARN3
    7   FORMAT (20X, 'AND BEHIND ', E10.3, ' ANG  (K=', I4, ')')
        WRITE (*,8) OPAMAX, KBMAX, RMAXGOOD, RMAX
    8   FORMAT(20X, 'MAXIMUM OPACITY =', E10.3, ' AT K=', I4, 
     >  '   RECOMMENDED RMAX =', E10.3, '  (PRESENTLY: ', E10.3, ')' ) 
      ENDIF

      IF (TAUBMAX .LT. 0.1) THEN 
        RMAXGOOD = TAUBMAX * RMAX * 3. 
        OPAMAX   = TAUBMAX / RMAX
        WRITE (hOUT,'(A)') 
     >     ' COMO> *** WARNING: RMAX MUCH LARGER THAN NECESSARY'
        WRITE (hOUT,8) OPAMAX, KBMAX, RMAXGOOD, RMAX
      ENDIF

      IF (LSHTOT .GT. 0)   CALL PRIHTOT (HTOTC,ND,LSHTOT,JOBNUM,MODHEAD)
      IF (LPLOHTOT .GT. 0) CALL PLOHTOT (HTOTC,ND,JOBNUM,MODHEAD,TEFF)

      IF (LSINT.GT.0)
     $      CALL PRIMINT (XJC,ND,XLAMBDA,NF,LSINT,JOBNUM,MODHEAD)
 
C**  OPACITY PLOTS
      IF (NPLOTOPA .GT. 0) 
     >    CALL PLOTOPA (ND, MODHEAD, NATOM, ATMASS, ABXYZ,
     >           FILLFAC, OPTIONPLOTOPA, NPLOTOPA,
     >           T,RNE,POPNUM,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $           NOM,KODAT,NDIM,N,MAXATOM,
     $           LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $           SIGMATHK,SEXPOK,EDGEK,
     $           K,NF,SIGMAKI,RADIUS,
     $           KONTNUP,KONTLOW,LASTKON,XDATA)

      !write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)

      CALL STAMP (OPSYS, 'COMO', TIM1)

      STOP 'O.K.'
      END
