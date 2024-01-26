      SUBROUTINE RMODCOLI(RADIUS,ENTOT,RNE,T,VELO,GRADI,XLAMBDA,FWEIGHT,
     $                   POPNUM,POPMIN,RSTAR,VDOP,MODHEAD,JOBNUM,XJC,
     $                   P,ND,NDDIM,NF,NFDIM,N,NDIM,NP,NPDIM,Z,
     $                   TEFF,HEDDI,EDDI,MODHIST,MAXHIST,
     >                   DENSCON, FILLFAC, ABXYZ, NATOM,
     >                   LASTREDISMODE, NCOLIP, NF2, XLAMBDA2,
     >                   OPARND, EPSGMAX, BEPSGMAXERR, MAXXDAT, XDATA,
     >                   XMSTAR, ZERO_RATES, HTOTMINUSND, HTOTNDCOR,
     >                   HTOTLlast, TAULAST, NEXTHYDRO, DTDRIN,
     >                   VMIC)
C*******************************************************************************
C***  READING OF THE MODEL FILE FOR MAIN PROGRAM "COLI" ************************
C*******************************************************************************

      !IMPLICIT NONE    !not yet possible - more code rework needed

      INTEGER, INTENT(IN) :: MAXHIST
      INTEGER, INTENT(INOUT) :: NF, NP, ND, NF2
      INTEGER, INTENT(OUT) :: NCOLIP, NEXTHYDRO

      REAL, DIMENSION(NF) :: HEDDI
      REAL, DIMENSION(2) :: XJC, EDDI
      REAL, DIMENSION(ND) :: EPSGMAX
      REAL, DIMENSION(NATOM) :: ABXYZ
      REAL, DIMENSION(NDDIM) :: ENTOT, TAULAST, VMIC
      REAL, DIMENSION(MAXXDAT) :: XDATA
      LOGICAL, DIMENSION(NDIM*NDDIM) :: ZERO_RATES
      REAL, DIMENSION(NDDIM*NDIM) :: POPNUM

      CHARACTER(MAXHIST*8) :: MODHIST

      REAL, DIMENSION(ND) :: RADIUS, T, VELO, GRADI, HTOTLlast

      REAL :: XMSTAR, POPMIN, HTOTMINUSND, HTOTNDCOR, DTDRIN, VTURBND
      LOGICAL :: BEPSGMAXERR


c*** tiefenabh. clumping nach goetz
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C *** NOTE: ASSURE 64-BIT TYPE FOR USE IN WRITMS
      CHARACTER(8) :: NAME
      CHARACTER(4) :: LASTREDISMODE

C***  Preset LASTREDISMODE and NCOLIP
      LASTREDISMODE = 'UNKN'
      NCOLIP = 0

      CALL OPENMS (3, IDUMMY, IDUMMY, 1, IERR)

C***  READ ALL RELEVANT DIMENSIONS
      CALL READMS (3, ND   , 1 , 'ND      ', IERR)
      CALL READMS (3, NP   , 1 , 'NP      ', IERR)
      CALL READMS (3, NF   , 1 , 'NF      ', IERR)
      CALL READMS (3, NF2  , 1 , 'NF2     ', IERR)
      CALL READMS (3, LAST , 1 , 'MODHIST ', IERR)
      OPARND = 0.
      CALL READMS (3, OPARND, 1, 'OPARND  ', IERR)
      HTOTMINUSND = 0.
      CALL READMS(3, HTOTMINUSND, 1, 'HTMND   ', IERR)
      HTOTNDCOR = 0.
      CALL READMS(3, HTOTNDCOR  , 1, 'HTNDCOR ', IERR)
C***  Read (theoretical) DTDR at inner boundary
C***  (calculated in STEAL->CALCDTDRIN)
      CALL READMS (3, DTDRIN, 1, 'DTDRIN  ', IERR)
      IF (IERR == -10. .OR. DTDRIN < 0.) THEN
C***    Fallback if not yet calculated or negative gradient
        DTDRIN = -99.
      ENDIF

C***  CHECK WHETHER DIMENSIONING IS SUFFICIENT
      IF (ND .GT.N DDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      IF (NP .GT. NPDIM) THEN
            CALL REMARK ('TOO MANY IMPACT-PARAMETER POINTS')
            STOP 'ERROR'
            ENDIF
      IF (NF .GT. NFDIM) THEN
            CALL REMARK ('TOO MANY FREQUENCY POINTS')
            STOP 'ERROR'
            ENDIF
      IF (LAST .GT. MAXHIST) THEN
            CALL REMARK ('OLD MODEL HISTORY TOO LONG')
            STOP 'ERROR'
            ENDIF

      CALL READMS (3, RADIUS ,  ND , 'R       ', IERR)
      CALL READMS (3, ENTOT  ,  ND , 'ENTOT   ', IERR)
      DENSCON(1:ND) = 0.0
      CALL READMS(3,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               WRITE(0,*)'Error in depth-dep. clumping during COLI!'
               STOP 'Error in depth-dep. clumping during COLI!'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
c      CALL READMS (3, DENSCON, ND,    'DENSCON ', IERR)
c      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
c      FILLFAC(1:ND) = 1. / DENSCON(1:ND)
      CALL READMS (3, RNE    ,  ND , 'RNE     ', IERR)
      CALL READMS (3, T      ,  ND , 'T       ', IERR)
      CALL READMS (3, VELO   ,  ND , 'VELO    ', IERR)
      CALL READMS (3, GRADI  ,  ND , 'GRADI   ', IERR)
      CALL READMS (3, P      ,  NP , 'P       ', IERR)
      CALL READMS (3, Z      ,ND*NP, 'Z       ', IERR)
      CALL READMS (3, XLAMBDA,  NF , 'XLAMBDA ', IERR)
      CALL READMS (3, XLAMBDA2, NF2, 'XLAMBD2 ', IERR)
      CALL READMS (3, FWEIGHT,  NF , 'FWEIGHT ', IERR)
      CALL READMS (3, HTOTLlast, ND, 'HTOTL   ', IERR)
      CALL READMS (3, POPNUM , ND*N, 'POPNUM  ', IERR)
      CALL READMS (3, POPMIN ,   1 , 'POPMIN  ', IERR)
      IF (IERR < 0) THEN
C***  POPMIN not on MODEL file => read from CARDS
        POPMIN = 1.E-100
      ENDIF
      CALL READMS (3, RSTAR  ,   1 , 'RSTAR   ', IERR)
      CALL READMS (3, MODHEAD,  13 , 'MODHEAD ', IERR)
      CALL READMS (3, VDOP   ,   1 , 'VDOP    ', IERR)
      CALL READMS (3, MODHIST, LAST, 'MODHIST ', IERR)
      CALL READMS (3, ABXYZ  ,NATOM, 'ABXYZ   ', IERR)
C***  Read EPSGMAX from last COLI+
      CALL READMS (3, EPSGMAX, ND-1,    'EPSGMAX ', IERR)
      BEPSGMAXERR = IERR .EQ. -10

      CALL READMS (3, TAULAST,  ND , 'TAUROSS ', IERR)

C***  READ TEFF FROM MODEL FILE
C***  -- CAUTION: TEFF RECORD MAY NOT EXIST IN VERY OLD BERLIN MODELS
      IERR=1
      CALL READMS (3,TEFF,1,'TEFF    ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK (' TEFF NOT ON MODEL FILE: USE OLDSTART')
         PRINT *,     ' TEFF NOT ON MODEL FILE: USE OLDSTART'
         STOP 'ERROR'
         ENDIF

C***  Flags for the POPMIN levels from STEAL
      CALL READMS (3,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
C*    Default if variable does not exist yet
      IF (IERR .EQ. -10) THEN
        DO I=1, N*ND
          ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF

C***  Read (micro-)turbulent velocity
      CALL READMS (3,VMIC, ND,  'VMIC    ', IERR)
      IF (IERR == -10) THEN
C***    old MODEL file => only one VTURB value
        CALL READMS(3,VTURBND,1, 'VTURB   ', IERR)
C***    update MODEL file: add new variable
        IF (IERR == -10) VTURBND = 0.
        DO L=1, ND
          VMIC(L) = VTURBND * SQRT(2.)
        ENDDO
        CALL WRITMS(3,VMIC, ND, 'VMIC    ',-1, IDUMMY, IERR)
      ENDIF


C***  READ ALL CONTINUUM INTENSITIES
      ND3=3*ND
      DO 15 K=1,NF
      WRITE (NAME,'(A3,I4,A1)') 'XJC',K, ' '
      CALL READMS(3,XJC(1+ND*(K-1)),ND,NAME, IERR)
      IF (K <= 999) THEN
        WRITE (NAME,'(A4,I3,A1)') 'EDDI',K, ' '
      ELSE
        WRITE (NAME,'(A4,I4)') 'EDDI', K
      ENDIF

      CALL READMS (3,EDDI,ND3,NAME, IERR)

      write(0,*) "array exceeded? inga"
      write(0,*) ND3
      write(0,*) K
      write(0,*) EDDI(ND3)
      HEDDI(K)=EDDI(ND3)
   15 CONTINUE

C***  Read REDISMODE and NCOLIP
      CALL READMS(3, LASTREDISMODE, 1, 'REDISMO ', IERR)
      CALL READMS(3, NCOLIP,        1,        'NCOLIP  ', IERR)

C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3,XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         XDATA(1) = 0.
      ENDIF

      !Read hydro stuff
      CALL READMS (3,XMSTAR,  1, 'XMSTAR  ', IERR)
      CALL READMS (3,NEXTHYDRO,1,'NXTHYDRO', IERR)
      !If NXTHYDRO does not exist in the MODEL (=old model)
      ! set value to -1 to deactivate any hydro-related stuff
      IF (IERR == -10) THEN
        NEXTHYDRO = -1
      ENDIF

C***  READ, INCREASE AND WRITE JOBNUMBER
      CALL READMS (3,JOBNUM,1,'JOBNUM  ', IERR)
      JOBNUM=JOBNUM+1
C      IF (JOBNUM .GE. 1000) JOBNUM=JOBNUM-1000
      CALL WRITMS (3,JOBNUM,1,'JOBNUM  ',-1, IDUMMY, IERR)

      RETURN
      END
