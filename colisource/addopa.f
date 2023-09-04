      SUBROUTINE ADDOPA (ND, NDDIM, MAXLIN, MAXIND, LIND, LINDS, 
     >             XK, XKMID, XKRED, DELTAX, FWEIGHTL,
     >             PARALAS, LASER, LASERV, ALN, VDOPUNIT, VDOPDD,
     >             WS, ETAL, OPAL, ETA, ETANOTH, OPA, ETAK, ETAKNOTH, 
     >             OPAK, OPAKNOTH, THOMSON, PWEIGHT, NOM,
     >             OPAFE, ETAFE, BFECHECK, BLASERL, NUP, LOW, N, LEVEL, 
     >             OPAKFE, NATOM, MAXATOM, KODAT, OPACELEM, OPAKELEM,
     >             ETACELEM, ETAKELEM, MAXION, OPACION, OPAKION, 
     >             ETACION, ETAKION, NCHARG, OPAFEION, ETAFEION, 
     >             OPAKNOFENOTH, ETAKNOFENOTH,
     >             OPAFEFT, ETAFEFT, OPAKFEFT, ETAKFEFT)
C**********************************************************************
C***  TREATMENT OF LINE BLENDS:
C***  ADD OVERLAPPING OPACITIES AND EMISSIVITIES
C***  - ALSO GENERATES THE LASER WARNING (ONLY LASER-VERSION 1)
C***  - CALLED FROM: COLI, SETXJFINE ( <- with ND=1 )
C
C
C     The different OPA variables have the following meaning
C     (all of them are depth-dependent and only for the current fine freq.)
C       OPA       - sum of bf and ff opacity plus Thomson
C       OPAL      - bb opacity per line transition
C       OPAFE     - Fe bb opacity
C     By the end of the routine these values are filled as follows:
C       OPAK      - total opacity (bf, bb, ff, and Thomson)
C       OPAKNOTH  - total opacity w/o Thomson (bf, bb, ff)
C       THOMSON   - fraction of Thomson opacity relative to OPAK
C
C**********************************************************************

      INTEGER, INTENT(IN) :: ND, NDDIM, MAXLIN, MAXIND, NATOM, MAXATOM
      REAL, INTENT(IN) :: ALN, VDOPUNIT, FWEIGHTL

      DIMENSION LIND(MAXLIN),LINDS(MAXLIN),XKMID(MAXIND),XKRED(MAXIND)
      DIMENSION NUP(MAXLIN), LOW(MAXLIN)
      DIMENSION OPAL(NDDIM,MAXLIN), ETAL(NDDIM,MAXLIN)
      DIMENSION ETA(ND), ETANOTH(ND), OPA(ND)
      REAL, DIMENSION(NDDIM,MAXLIN) :: PWEIGHT, WS
      REAL, DIMENSION(NDDIM) :: OPAK, ETAK, OPAKNOTH, ETAKNOTH, THOMSON
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      REAL, DIMENSION(ND) :: OPAFE, ETAFE, OPAKFE, 
     >                       OPAKNOFENOTH, ETAKNOFENOTH,
     >                       OPAFEFT, ETAFEFT, OPAKFEFT, ETAKFEFT
      REAL, DIMENSION(NATOM, ND) :: OPAKELEM, OPACELEM,
     >                              ETACELEM, ETAKELEM
      REAL, DIMENSION(ND, NATOM, MAXION) :: OPAKION, OPACION,
     >                                      ETACION, ETAKION
      REAL, DIMENSION(ND, MAXION) :: OPAFEION, ETAFEION
      REAL, DIMENSION(ND, NATOM), INTENT(IN) :: VDOPDD
      INTEGER, DIMENSION(N) :: NOM, NCHARG

      CHARACTER(10), DIMENSION(N) :: LEVEL
      LOGICAL LASER, BFECHECK, BLASERL(MAXIND)
      INTEGER :: L, NL, NA, NAFE

      REAL, PARAMETER :: WPI      = 1.772454    !WPI = SQRT(PI)
      REAL, PARAMETER :: CLIGHTKM = 2.9979E5    !CLIGHT = SPEED OF LIGHT IN KILOMETER/SECOND
      
      
      LASER = .FALSE.
      NAFE = KODAT(26)      !Atomic number of iron

C***  REAL LINE OPACITIES FROM ALL INVENTED LINES
      DO L=1,ND
        ETAK(L)     = 0.0
        ETAKNOTH(L) = 0.0
        OPAK(L)     = 0.0
        OPAKNOTH(L) = 0.0
        OPAKFE(L)   = 0.0
        OPAKNOFENOTH(L) = 0.0
        ETAKNOFENOTH(L) = 0.0
        DO NA=1, NATOM
          OPAKELEM(NA,L) = 0.0
          ETAKELEM(NA,L) = 0.0
          DO ION=1, MAXION
            OPAKION(L,NA,ION) = 0.0
            ETAKION(L,NA,ION) = 0.0
          ENDDO
        ENDDO
        OPAKFEFT(L)   = 0.0
        ETAKFEFT(L)   = 0.0
      ENDDO

      DO NL=1, MAXLIN
        IF (LIND(NL) .EQ. 0) CYCLE
        IF (XK .GT. XKRED(LINDS(NL))) CYCLE
        
        NA = NOM(LOW(NL))
        ION = NCHARG(LOW(NL)) + 1

        DO L=1, ND
           DK = (XK - XKMID(LINDS(NL))) * DELTAX
     >                    * VDOPUNIT/VDOPDD(L,NA)           
           
           PWEIGHT(L,NL) = VDOPUNIT/VDOPDD(L,NA) * EXP(-DK*DK)
           WS(L,NL) = WS(L,NL) + PWEIGHT(L,NL) * FWEIGHTL
           PHI = PWEIGHT(L,NL) / WPI

           ETAK(L) = ETAK(L) + ETAL(L,NL)*PHI
           OPAK(L) = OPAK(L) + OPAL(L,NL)*PHI

           OPAKELEM(NA,L) = OPAKELEM(NA,L) + OPAL(L,NL)*PHI
           ETAKELEM(NA,L) = ETAKELEM(NA,L) + ETAL(L,NL)*PHI
           OPAKION(L,NA,ION) = OPAKION(L,NA,ION) + OPAL(L,NL)*PHI
           ETAKION(L,NA,ION) = ETAKION(L,NA,ION) + ETAL(L,NL)*PHI
           OPAG = (PARALAS-1.0) * MAX(OPA(L), 0.)        !Goetz version
           BLASERL(LIND(NL)) = 
     >          BLASERL(LIND(NL)) .OR.  OPAL(L,NL) .LT. OPAG
        ENDDO
      ENDDO


C***  IRON: ADD IRON LINE OPACITY AND EMISSIVITY IF NECESSARY
      IF (BFECHECK) THEN
         DO L=1, ND
            ETAK(L) = ETAK(L) + ETAFE(L)
            OPAK(L) = OPAK(L) + OPAFE(L)
            OPAKFE(L) = OPAKFE(L) + OPAFE(L)
            NA = KODAT(26)
            OPAKELEM(NA,L) = OPAKELEM(NA,L) + OPAFE(L)
            ETAKELEM(NA,L) = ETAKELEM(NA,L) + ETAFE(L)
            DO ION=1, MAXION
              OPAKION(L,NA,ION) = OPAKION(L,NA,ION) + OPAFEION(L,ION)
              ETAKION(L,NA,ION) = ETAKION(L,NA,ION) + ETAFEION(L,ION)
            ENDDO
         ENDDO
      ENDIF

C***  Laser treatment (all versions): 
C***       restrict Total Opacity to > PARALAS * CONTINUUM
      IF (LASERV > -1) THEN
        DO L=1, ND
c          OPAG = (PARALAS-1.0) * OPA(L)
          OPAG = (PARALAS-1.0) * MAX(OPA(L), 0.)
          LASER = LASER .OR. (OPAK(L).LT.OPAG)
          OPAK(L)=AMAX1(OPAG,OPAK(L))
        ENDDO
      ENDIF

C***  Add the Continuum Values to ETAK, ETAKNOTH and OPAK
      DO L=1, ND
        ETAKNOTH(L) = ETANOTH(L) + ETAK(L)
        ETAK(L)     = ETA(L) + ETAK(L)
        OPAG        = PARALAS * OPA(L)
        OPAKNOTH(L) = OPA(L) * (1.-THOMSON(L)) + OPAK(L)
        OPAK(L)     = OPA(L) + OPAK(L)
C***    Perform element- and ionspecific calulations 
C***     to enable radiative driving analysis (ACCELEM) plots 
        DO NA=1, NATOM
          IF (NA /= NAFE) THEN
            OPAKNOFENOTH(L) = OPAKNOFENOTH(L) 
     >                         + OPACELEM(NA,L)             
     >                         + OPAKELEM(NA,L)
            ETAKNOFENOTH(L) = ETAKNOFENOTH(L)
     >                         + ETACELEM(NA,L) + ETAKELEM(NA,L)
          ELSE
            OPAKFE(L)   = OPAKFE(L) + OPACELEM(NAFE,L) 
            OPAKFEFT(L) = OPAFEFT(L) + OPACELEM(NAFE,L)
            ETAKFEFT(L) = ETAFEFT(L) + ETACELEM(NAFE,L) 
          ENDIF
          OPAKELEM(NA,L) =  OPACELEM(NA,L) + OPAKELEM(NA,L)
          ETAKELEM(NA,L) =  ETACELEM(NA,L) + ETAKELEM(NA,L)
          DO ION=1, MAXION
            OPAKION(L,NA,ION) = OPACION(L,NA,ION) + OPAKION(L,NA,ION)
            ETAKION(L,NA,ION) = ETACION(L,NA,ION) + ETAKION(L,NA,ION)
          ENDDO
        ENDDO

      ENDDO      
      

      RETURN
      END
