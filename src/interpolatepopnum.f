      SUBROUTINE INTERPOLATEPOPNUM(POPNUM,  !old (and new) popnumber array
     >                             POPOLD,  !old POPNUM array
     >                             POPMIN,  !POPMIN value (replaces zeroes)
     >                             Rnew,    !new radius grid vector
     >                             Rold,    !old radius grid vector
     >                             ENTOTnew, !old total particle number
     >                             ENTOTold, !new total particle number
     >                             RNEnew,  !new relative electron number density
     >                             RNEold,
     >                             Tnew,    !new temperature structure
     >                             Told,    !old temperature structure
     >                             N,       !number of levels in DATOM
     >                             ND,      !number of depth points
     >                             ABXYZ,
     >                             NFIRST,  !Array with first level number of element blocks in level list
     >                             NLAST,   !similar to NFIRST, but with last level number
     >                             NATOM,   !Number of different elements in the model
     >                             NCHARG,  !charge per level
     >                             WEIGHT,  !Weight (g) factor per level
     >                             EION,    !Ionization energy per level 
     >                             ELEVEL,  !Level energy per level
     >                             NOM,     !corresponding element index for each level
     >                             bUseENTOT)     
C**********************************************************************
C***
C***    Interpolation of Popnumbers on new Radius-Grid
C***     crucial if radius grid has been updated in-between iterations
C***    called from: ENSURETAUMAX, HYDROSOLVE, HYDRO_REGRID
C***
C***    The interpolation is performed in  (log n_i) over (log n_tot)
C***    The older version (interpolation over radius) 
C***    can be activated by setting bUseENTOT to .FALSE.
C***
C**********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: N, ND, NATOM
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ
      REAL, DIMENSION(ND), INTENT(IN) :: Rnew, Rold, Tnew, Told,
     >                                   ENTOTold, ENTOTnew, RNEold
      REAL, DIMENSION(ND), INTENT(INOUT) :: RNEnew
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPOLD
      REAL, DIMENSION(ND, N), INTENT(INOUT) :: POPNUM
      INTEGER, DIMENSION(N), INTENT(IN) :: NCHARG, NOM
      REAL, DIMENSION(N), INTENT(IN) :: WEIGHT, EION, ELEVEL 
      REAL, INTENT(IN) :: POPMIN
      
      INTEGER, PARAMETER :: NIPMAX = 5000
      INTEGER, PARAMETER :: NDIPMAX = 500
      
      REAL, DIMENSION(NDIPMAX) :: ENTOToldLOG, RoldLOG, POPJLOG
      REAL, DIMENSION(NIPMAX) :: ENLTE, DEPARToldND
        
      INTEGER :: L, J, NA, NFIRNA, NLANA, NDr
      REAL :: SUMME, POPJLOGnewL, ENTOTnewLOGL, RnewLOGL, RNEL, TL, ENE

      LOGICAL, INTENT(IN) :: bUseENTOT
      
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      IF (N > NIPMAX) THEN
        WRITE (hCPR,'(A)') 'IPVECRLOG: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'IPVECRLOG: NIPMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'N = ', N, ', NIPMAX = ', NIPMAX
        STOP 'FATAL ERROR IN STEAL->INTERPOLATEPOPNUM'      
      ENDIF
      
      IF (ND > NDIPMAX) THEN
        WRITE (hCPR,'(A)') 'IPVECRLOG: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'IPVECRLOG: NDIPMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDIPMAX = ', NDIPMAX
        STOP 'FATAL ERROR IN STEAL->INTERPOLATEPOPNUM'      
      ENDIF
      
      
C***  Skip this routine if called for a not yet defined popnumber array
      IF (MAXVAL(POPOLD) <= 0.) THEN
        POPNUM = POPOLD
        RETURN
      ENDIF

C***  Calculate departure coefficients for old innermost depth point
C     (required for inner extrapolation)
      ENE = RNEold(ND) * ENTOTold(ND)
      TL = Told(ND)
      CALL LTEPOP (N, ENLTE, TL, ENE, 
     >             WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >             ABXYZ, NFIRST, NLAST, NATOM)
      DO J=1, N
        DEPARToldND(J) = POPOLD(ND,J) / ENLTE(J) 
      ENDDO

C***  Recommended branch: interpolation of (log n_i) over (log n_tot)
      IF (bUseENTOT) THEN
        DO J=1, N
C         prepare logarithmic  vectors
          DO L=1, ND
            ENTOToldLOG(L) = LOG10(ENTOTold(L))
C           Note: minimum popnumber ist POPMIN
            POPJLOG(L) = LOG10(MAX(POPOLD(L,J),POPMIN))
          ENDDO          
C         Perform interpolation
          dploop: DO L=1, ND
            ENTOTnewLOGL = LOG10(ENTOTnew(L))
            IF (ENTOTnewLOGL > ENTOToldLOG(ND)) THEN
               TL = Tnew(L)
C              !use old RNE(ND) for ENE here (can this be improved?)
               ENE = RNEold(ND) * ENTOTnew(L)
               CALL LTEPOP (N, ENLTE, TL, ENE, 
     >                      WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >                      ABXYZ, NFIRST, NLAST, NATOM)                            
               POPJLOGnewL = LOG10( ENLTE(J) * DEPARToldND(J) )
            ELSEIF (ENTOTnewLOGL < ENTOToldLOG(1)) THEN
              !less dense than old outermost value => take old outer boundary value
              POPJLOGnewL = POPJLOG(1)
            ELSE
              CALL SPLINPOX(POPJLOGnewL, ENTOTnewLOGL,
     >                     POPJLOG, ENTOToldLOG, ND)
            ENDIF
            POPNUM(L,J) = 10**(POPJLOGnewL)
          ENDDO dploop
        ENDDO

      ELSE
C***   double-logarithmic interpolation over radius
        DO J=1, N
          DO L=1, ND
            IF (L > 1) THEN
C***          Account for the possibility that the old radius might be shifted
C***          and thus we have to determine the number of points NDr in the "normal"
C***          radius Range 1...RMAX. 
C***          (One point for R < 1 is allowed to avoid inner cutoffs.)
              IF (Rold(L-1) < 1. .OR. Rold(L) <= 0.) EXIT
            ENDIF
            RoldLOG(L) = LOG10(Rold(L))
            NDr = L
C           Note: minimum popnumber ist POPMIN
            POPJLOG(L) = LOG10(MAX(POPOLD(L,J),POPMIN))
          ENDDO          
          DO L=1, ND            
            RnewLOGL = LOG10(Rnew(L))
            IF (RnewLOGL > RoldLOG(1)) THEN
C             If R outside the old grid, use old outermost value
              POPJLOGnewL = POPJLOG(1)
            ELSEIF (RnewLOGL < RoldLOG(NDr)) THEN
C             If R < 1 use innermost value (should never happen)
              POPJLOGnewL = POPJLOG(NDr)
            ELSE
              CALL SPLINPOX(POPJLOGnewL,RnewLOGL,POPJLOG,RoldLOG,NDr)
            ENDIF
            POPNUM(L,J) = 10**(POPJLOGnewL)
          ENDDO
        ENDDO
      ENDIF

C***  Renormalization      
      DO L=1, ND
        DO NA=1, NATOM
          SUMME=0.0
          NFIRNA = NFIRST(NA)
          NLANA = NLAST(NA)
          DO J = NFIRNA, NLANA
            SUMME = SUMME + POPNUM(L,J)
          ENDDO
          SUMME = SUMME / ABXYZ(NA)
          IF (SUMME /= 0.) THEN
            DO J=NFIRNA,NLANA
              POPNUM(L,J) = POPNUM(L,J) / SUMME
            ENDDO
          ENDIF
        ENDDO

        RNEL=0.0
        DO J=1, N
          RNEL = RNEL + NCHARG(J) * POPNUM(L,J)
        ENDDO
        RNEnew(L)=RNEL
      ENDDO

      RETURN
      END
