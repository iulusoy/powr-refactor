      SUBROUTINE POPSMALL(L, ND, N, POPMIN, BFIRSTITER, ZERO_RATES, 
     >                    IMAXPOP, NATOM, NFIRST, NLAST, EN, 
     >                    NRANK, NCHARG, RATCO, V1, ITNEL, TL, KODAT,
     >                    LEVEL, MAXATOM, ABXYZ, iZRType)
C***********************************************************************
C*** handling of small popnumbers (usually values lower than POPMIN)
C***  chooses between wrh, goetz and no flagging depending on CARDS
C***
C*** called by STEAL->LINPOP->COMA
C***********************************************************************
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: L, N, ND, NATOM, ITNEL, NRANK, 
     >                       MAXATOM, iZRType
      REAL, DIMENSION(NRANK) :: EN, V1
      REAL, DIMENSION(NRANK, NRANK) :: RATCO
      REAL, DIMENSION(NATOM) :: ABXYZ
      INTEGER, DIMENSION(N) :: NCHARG
      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST, IMAXPOP
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      CHARACTER(10), DIMENSION(N) :: LEVEL
      
      LOGICAL, DIMENSION(N, ND) :: ZERO_RATES
      LOGICAL :: BFIRSTITER, bLOW

      INTEGER :: NA, J, I, NCTEST, NFIRNA, NLANA, 
     >           NSMALL, NBOTTOM, NSTART, LL
      REAL :: POPMIN, POPLOW, TL
      
      INTEGER, SAVE, DIMENSION(100) :: NL_T, NL_TP, NL_B, NL_BM, NL_BP

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)                                                    
                        
      IF (iZRType == 0) THEN
C************ Method 0: Do not flag any levels at all ******************
        DO J=1, NRANK
          ZERO_RATES(J,L)= .FALSE.
        ENDDO

C***    COLUMN IMAXPOP (I.E. INDEX OF MAX. POPNUMBER PER ELEMENT):
C***         NUMBER CONSERVATION FOR EACH ELEMENT (NA)
C***    REMARK: TOTAL NUMBER CONSERVATION IS IMPLICITLY ENSURED
        DO NA=1, NATOM
           DO I=NFIRST(NA), NLAST(NA)
             RATCO(I,IMAXPOP(NA))=1.
           ENDDO
           V1(IMAXPOP(NA))=ABXYZ(NA)
        ENDDO        
        
      ELSEIF (iZRType == 1) THEN      
C************ Method 1: Classic preguessing wrh method *****************      
        DO NA=1,NATOM
          NFIRNA=NFIRST(NA)
          NLANA=NLAST(NA)

C***      Check if all Rate Coefficients in one column are non-zero
C***      (otherwise: the matrix is singular!)
C***      and store logical flag ZERO_RATES for later use 
          IF (ITNEL == 1) THEN  
            CALL FLAG_ZERORATES(NFIRNA, NLANA, RATCO, N, NRANK,
     >                          IMAXPOP(NA), EN, POPMIN, 
     >                          ZERO_RATES(1,L))
C***      new, wrh 25-Feb-2015:
C***      in case of iron, levels MUST be popmin if they are also flagged
C***      (update Mar 2015: system is now used for all elements)
C***      at the next-inner depth point 
cc            IF (NA .EQ. KODAT(26) .AND. L < ND) THEN
            IF (L < ND) THEN
               DO J=NFIRNA, NLANA
                  IF (ZERO_RATES(J,L+1)) THEN
                     DO LL=L, 1, -1
                        ZERO_RATES(J,LL)=.TRUE.
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
          ENDIF

C***      COLUMN IMAXPOP (I.E. INDEX OF MAX. POPNUMBER PER ELEMENT):
C***         NUMBER CONSERVATION FOR EACH ELEMENT (NA)
C***      REMARK: TOTAL NUMBER CONSERVATION IS IMPLICITLY ENSURED
          DO I=NFIRST(NA), NLAST(NA)
            RATCO(I,IMAXPOP(NA))=1.
          ENDDO
          V1(IMAXPOP(NA))=ABXYZ(NA)

C***      If ZERO_RATES: Replace diagonal element by 1.0
C***      and the rest of this column by 0.0
          DO J = NFIRNA, NLANA
            IF (.NOT. ZERO_RATES(J,L)) CYCLE
            DO I = NFIRNA, NLANA
               RATCO(I,J) = .0
               RATCO(J,I) = .0
            ENDDO
            RATCO(J,J) = 1.
            V1(J) = POPMIN
            EN(J) = POPMIN
          ENDDO
        ENDDO
        
      ELSEIF (iZRType == 2) THEN
C************ Method 2: Goetz' method: Top and bottom ladder ***********
C***    This method does not rely on any preguessing, but only on
C***    the current population. Level flags are only determined in the 
C***    first Scharmer iteration, but rates are switched off throughout.
      
        IF (NATOM > 100) THEN
          WRITE (hCPR,*) 'Error in POPSMALL> More than 100 Elements'
          STOP 'FATAL ERROR IN STEAL'
        ENDIF
        !Set limit slighly higher than POPMIN to account for renormalization effects
C        POPLOW = 1.005 * POPMIN
        POPLOW = POPMIN
        IF (BFIRSTITER) THEN
          DO J=1, N
            ZERO_RATES(J,L)= .FALSE.
          ENDDO
        ENDIF
      
        DO NA=1,NATOM  !~ ~ start of atom loop ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
          NFIRNA=NFIRST(NA)
          NLANA=NLAST(NA)
          
          IF (BFIRSTITER) THEN
C***        Flag levels with POPMIN: These are usually the lowest and/or
C***        highest levels of an atom. Therefore we start with the lowest
C***        level and set the switch "bLOW" if this level is not larger
C***        than POPMIN. This switch is changed to .FALSE. as soon as a
C***        level above POPMIN occurs.  Levels with POPMIN afterwards are
C***        interpreted as part of the highest (top) levels.
C***        Checks are only done for ion ground states, but higher levels
C***        are also switched off if the ground state is on POPMIN.
            NCTEST = -9999
            bLOW = (EN(NFIRNA) <= POPLOW)
            NL_T (NA) = 0   !first top level with/below POPLOW
            NL_TP(NA) = 0   !second top level with/below POPLOW 
            NL_B (NA) = 0   !last bottom level with/below POPLOW 
            NL_BM(NA) = 0   !last but one bottom level with/below POPLOW
            NL_BP(NA) = 0   !first bottom level above POPLOW
            DO J=NFIRNA, NLANA
              IF (NCHARG(J) /= NCTEST) THEN       !New Ion Ground State?
                NCTEST = NCHARG(J)
                IF (bLOW) THEN  !Are we in the 'Bottom' part of the model atom?
                  bLOW = (EN(J) <= POPLOW)
                  IF (bLOW) THEN  !If the Ion ground state is small -> set NL_B
                    NL_BM(NA) = NL_B(NA)
                    NL_B (NA) = J
                  ELSE
                    NL_BP(NA) = J
                  ENDIF
                ELSE 
C***              Else we are in the 'Top' part of the model atom             
C***              We only store the two ground levels which are at/below POPLOW and
C***              assume all higher ground states to be <= POPLOW as well
                  IF (EN(J) <= POPLOW) THEN !If the Ion gound state is small -> set NL_T
                    IF (NL_T(NA) == 0) THEN
                      NL_T (NA) = J
                    ELSE IF (NL_TP(NA) == 0) THEN
                      NL_TP(NA) = J
                    ENDIF
                  ENDIF
                ENDIF                
              ENDIF
            ENDDO
            IF (bLOW) THEN !If bLOW is still true, all ions are at POPLOW
              NL_T (NA) = 0
              NL_TP(NA) = 0
              NL_B (NA) = 0
              NL_BM(NA) = 0
              NL_BP(NA) = 0
            ENDIF
          ENDIF ! - - - end of IF BFIRSTITER
          
C***      Switch off all levels with POPLOW, except the upmost low one and the lowest top one.
C***      (This means we leave a control level, but switch off all higher levels of the same stage.)
C***      With this method, Levels can only be switched back on one ion after another
          NSMALL = NL_T(NA)
          IF (NSMALL > 0 .AND. NSMALL < NLANA) THEN
C***        Switch off remaining upper levels for the current element (one ground state level remains)
            DO I=NSMALL+1, NLANA
              IF (BFIRSTITER) THEN                
                EN(I) = POPMIN
C***            For top levels, we also force that levels once switched off cannot be 
C***            switched back on further outwards. This prevents artifacts caused by existence 
C***            of POPMIN and is especially stabilizing for iron superlevels.
                DO LL=L, 1, -1
                  ZERO_RATES(I,LL)=.TRUE.
                ENDDO
              ENDIF
            ENDDO
          ENDIF
C***      For the bottom, it is slightly more complicated as we have to switch
C***      off levels up to NSTART-1, but also want to keep a control level (NBOTTOM)
          NBOTTOM = NL_BM(NA)
          NSTART  = NL_B(NA)
C***      Old Goetz idea was NBOTTOM = NL_BP(NA) and NSTART = NL_BP(NA)
C***      but he claims in his code comments this is causing problems with Neon
          IF (NBOTTOM > 0.) THEN
            DO I=NFIRNA, NSTART-1
              IF (I==NBOTTOM) CYCLE
              IF (BFIRSTITER) THEN                
                EN(I) = POPMIN
                ZERO_RATES(I,L) = .TRUE.
              ENDIF
            ENDDO
          ENDIF

C***      If ZERO_RATES: Replace diagonal element by 1.0
C***      and the rest of this column by 0.0
          DO J = NFIRNA, NLANA
            IF (.NOT. ZERO_RATES(J,L)) CYCLE
            DO I = NFIRNA, NLANA
               RATCO(I,J) = .0
               RATCO(J,I) = .0
            ENDDO
            RATCO(J,J) = 1.
            V1(J) = POPMIN
          ENDDO
          
        ENDDO  !~ ~ end of atom loop ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

C***    COLUMN IMAXPOP (I.E. INDEX OF MAX. POPNUMBER PER ELEMENT):
C***         NUMBER CONSERVATION FOR EACH ELEMENT (NA)
C***    REMARK: TOTAL NUMBER CONSERVATION IS IMPLICITLY ENSURED
        DO NA=1, NATOM
           DO I=NFIRST(NA), NLAST(NA)
             RATCO(I,IMAXPOP(NA))=1.
           ENDDO
           V1(IMAXPOP(NA))=ABXYZ(NA)
        ENDDO
        
      ELSE
C***********************************************************************
        WRITE (hCPR,*) 'POPSMALL> Invalid method given for ZERO_RATES!'
        STOP 'FATAL ERROR IN STEAL'        
      ENDIF

      
      RETURN
      
      END
      
