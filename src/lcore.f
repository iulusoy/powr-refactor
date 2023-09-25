      SUBROUTINE LCORE (XRED,XBLUE,GAMMAL,LASTIND,INDLOW,INDNUP,
     $      GAMMAR,IPRILC,MODHEAD,JOBNUM,OPALOLD,LASTINDAUTO,
     >      XLAMAPPMEAN, ALOMIN, bUSEALO, bLAMAPPCOLI,
     $      L,ND,VELO,GRADI,RADIUS,SLOLD,XMAX,ERXMIN,
     $      VDOPUNIT,RSTAR,ENTOT,EN,NF,XLAMBDA,ELEVEL,NOM,
     $      NDIM,N,NCHARG,WEIGHT,EINST,LINE,NSCHAR,BDIAG,GAMMAD)
C*******************************************************************************
C***  DETERMINATION OF THE CMF-FREQUENCIES XRED, XBLUE WHICH CONFINE THE
C***  OPTICALLY THICK LINE CORES
C*******************************************************************************
 
      INTEGER, INTENT(IN) :: LASTIND, LASTINDAUTO
      REAL, INTENT(IN) :: ALOMIN

      REAL, DIMENSION(ND, LASTINDAUTO) :: XLAMAPPMEAN

      COMMON / COMFUN / DELTAV,XMIN
      LOGICAL LINE(2)
      DIMENSION TAUMIN(0:2),GAMPRI(0:2)
      LOGICAL THIN
      DIMENSION VELO(ND),GRADI(ND),RADIUS(ND),ENTOT(ND)
      DIMENSION SLOLD(2),XRED(2),XBLUE(2)
      DIMENSION EINST(NDIM,NDIM),WEIGHT(NDIM),EN(NDIM),NCHARG(NDIM)
      DIMENSION XLAMBDA(NF)
      DIMENSION ELEVEL(NDIM)
      DIMENSION NOM(N)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND),OPALOLD(LASTIND)
      INTEGER NSCHAR

      REAL :: DELTAV, XMIN

      LOGICAL, DIMENSION(LASTINDAUTO), INTENT(INOUT) :: bUSEALO
      LOGICAL, DIMENSION(LASTIND) :: BDIAG

      LOGICAL :: bLAMAPPCOLI     !.true. if CARDS line XJLAPP COLI has been set

C***  BRANCH FOR GAMMAL=0 : LINE CORES ARE ZERO
      IF (GAMMAL.LE..0 .AND. GAMMAR.LE..0 .AND. GAMMAD.LE..0) THEN
        DO IND=1,LASTIND
          XRED(IND)=.0
          XBLUE(IND)=.0
          bUSEALO(IND) = .FALSE.
          SLOLD(IND)=TRANSFER('UNDEF', SLOLD(IND))
        ENDDO
        RETURN
      ENDIF
 
C***  BRANCH FOR GAMMAL .NE. .0
C***  VL = VELOCITY IN DOPPLER UNITS ,  GL = GRADIENT IN DOPPLER UNITS
      VL=VELO(L)/VDOPUNIT
      GL=GRADI(L)/VDOPUNIT
      RL=RADIUS(L)
      
      DO IND=1, LASTIND !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)
        XRED(IND)=0.
        XBLUE(IND)=0.
        
        SLOLD(IND)=TRANSFER('UNDEF', SLOLD(IND))
        OPALOLD(IND)=.0
C***    IRON: OMIT LINES TREATED WITH DIAGONAL OPERATOR, LASER CHECK
        IF (BDIAG(IND)) THEN
          XLAM=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
          CALL LIOP (EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     >               1,XLAM,ENTOT(L),EN,RSTAR,OPALOLD(IND),ETAL,
     >               VDOPUNIT)
          IF (OPALOLD(IND) .LE. .0) THEN
            SLOLD(IND)=TRANSFER('UNDEF', SLOLD(IND))
            GOTO 2
          ELSE
C***        Note: This setting of XBLUE is used as a switch for DERIV
C***              in order to enforce NRBs for the iron lines
cccccccccc            XBLUE(IND) = 1.
            SLOLD(IND) = ETAL/OPALOLD(IND)
            GOTO 2
          ENDIF
        ENDIF

C***    TAKE THE APPROPRIATE GAMMA : GAMMAR FOR RESONANCE LINE, GAMMAL ELSE
        IF (LOW .EQ. 1) THEN
          GAMMA=GAMMAR
        ELSE
          IF ((NOM(LOW) .NE. NOM(LOW-1)) .OR. (NOM(LOW).EQ.NOM(LOW-1)
     $        .AND. NCHARG(LOW).NE.NCHARG(LOW-1))) THEN
            GAMMA=GAMMAR
          ELSE
            GAMMA=GAMMAL
          ENDIF
        ENDIF
        IF (GAMMA .LE. .0) THEN
          IF (bLAMAPPCOLI) THEN
            XLAMAPPMEAN(L, IND) = 0.
            IF (.NOT. BDIAG(IND)) bUSEALO(IND) = .FALSE.
          ENDIF
          GOTO 2
        ELSEIF (bLAMAPPCOLI) THEN
C***    For Gamma damping, we first transformed into
C***      a TAU-like quantity, then reduce and transform back
C***      BOUNDARIES: No amplification
          IF (L == 1 .OR. L == ND) THEN
            XLAMAPPMEAN (L, IND) = .0
c            bUSEALO(IND) = .FALSE.
          ELSE
            TAU = -ALOG(1.-XLAMAPPMEAN(L,IND))
            XLAMAPPMEAN(L, IND) = 1. - EXP(-TAU/GAMMA)
            XLAMAPPMEAN(L, IND) = MAX(0., XLAMAPPMEAN(L, IND))
            IF (XLAMAPPMEAN(L, IND) < ALOMIN) XLAMAPPMEAN(L, IND) = .0
          ENDIF
      
        ENDIF

C***    RUDIMENTAL LINES : CORE IS SET TO ZERO
        IF (EINST(LOW,NUP) .EQ. -2.) GOTO 2
C***    BOUNDARY POINTS: OPTICALLY THIN (ZERO CORE)
        IF (L .EQ. ND) GOTO 2
        IF (L .EQ. 1) GOTO 2
C***    LINES WHICH WERE NOT TREATED IN THE RADIATION TRANSFER
C***    ARE ASSUMED TO HAVE ZERO CORE
        IF (.NOT.LINE(IND)) GOTO 2
        
C***    CALCULATE LINE OPACITY OPAL AT CURRENT DEPTH POINT AND LINE
        XLAM=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
        CALL LIOP (EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     >             1, XLAM, ENTOT(L),EN,RSTAR,OPALOLD(IND),ETAL,
     >             VDOPUNIT)
C***    LASER SECURITY : CORE OF LASER LINES ARE SET TO ZERO
        IF (OPALOLD(IND) .LE. .0) GOTO 2
      
C***  IF AMBIENT CONTINUUM IS OPTICALLY THICK : 
C***    RADIAL DIRECTION -----------------------------------------------
        TAU=OPALOLD(IND)/GL
        GDT=GAMMA /TAU
        DV1=VELO(1)/VDOPUNIT - VL
        DVND=VL - VELO(ND)/VDOPUNIT
        DELTAV=AMIN1(DV1,DVND)
        CALL COFREQ (XRR,XBR,XMAX,ERXMIN,GDT   ,THIN)
        IF (THIN) GOTO 2
        
C***    TRANSVERSAL DIRECTION  -----------------------------------------
        R1=RADIUS(1)
        TAUT=OPALOLD(IND)*RL/VL
        GDT=GAMMA /TAUT
        DELTAV=SQRT(1.-RL*RL/R1/R1)*VELO(1)/VDOPUNIT
        CALL COFREQ (XRT,XBT,XMAX,ERXMIN,   GDT,THIN)
        IF (THIN) GOTO 2
        
C***    FIND THE MINIMUM CORE FROM BOTH DIRECTIONS ---------------------
        XRED(IND)=AMAX1(XRR,XRT)
        XBLUE(IND)=AMIN1(XBR,XBT)
C***    STORE OLD LINE SOURCE FUNCTION
        SLOLD(IND)=ETAL/OPALOLD(IND)
        
    2   CONTINUE
C***    Non-core lines should not be considered for ALO (except IRON)
        IF (XRED(IND) >= XBLUE(IND) .AND. .NOT. BDIAG(IND)) THEN
c          IF (.NOT. bLAMAPPCOLI) bUSEALO(IND) = .FALSE.
        ENDIF

C***    STORE OPTICAL DEPTHS FOR PRINTOUT
        DO I=0,2
          IF (IND .EQ. IPRILC+I) THEN
            TAUMIN(I)=OPALOLD(IND)/AMAX1(GL,VL/RL)
            GAMPRI(I)=GAMMA
          ENDIF
        ENDDO
      ENDDO  !End of large IND loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      IF (IPRILC .GT. 0)
     $      CALL PRILC (IPRILC,LASTIND,XRED,XBLUE,TAUMIN,L,ND,ERXMIN,
     $                  MODHEAD,JOBNUM,GAMPRI)
 
      RETURN
      END
