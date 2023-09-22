      SUBROUTINE ADAPOP (POPNUM, ND, N, POPOLD, NDOLD, NOLD, NCHARG,
     >          NATOM, ABXYZ, NFIRST, NLAST, RNE, NTRANS, POPLTE, 
     >          BDEPART, ADPWEIGHT, RADIUS, ROLD, POPHELP, TAURCONT, 
     >          TAURCONTOLD, POPLTE_OLD, ENTOT, ENTOTOLD, 
     >          BTAUR, bUseENTOT, POPMIN) 
C***********************************************************************
C***  TRANSFORMATION OF POPULATION NUMBERS FROM OLD TO NEW MODEL ATOM
C***  Radically simplified version: wrh 10-Aug-2007
C***  The assignment of levels is directed by the vector NTRANS(J)
C***    IF NTRANS(J) = -1 : no action, POPNUM stays from WRSTART 
C***    IF NTRANS(J) =  0 : POPNUM set to ZERO 
C***    IF NTRANS(J) >  0 : assign POPOLD with index NTRANS
C***********************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
 
      INTEGER, INTENT(IN) :: N, NOLD, ND, NDOLD, NATOM
      REAL, INTENT(IN) :: POPMIN

      INTEGER, DIMENSION(N) :: NCHARG, NTRANS
      REAL, DIMENSION(N) :: ADPWEIGHT(N)      

      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST
      REAL, DIMENSION(NATOM) :: ABXYZ
      
      REAL, DIMENSION(ND) :: RNE, RADIUS, TAURCONT, ENTOT
      REAL, DIMENSION(NDOLD) :: ROLD, TAURCONTOLD, ENTOTOLD, 
     >                          ENTOTOLDLOG, POPHELPJLOG

      REAL, DIMENSION(ND, N) :: POPNUM, POPLTE
      REAL, DIMENSION(ND, NOLD) :: POPOLD
      REAL, DIMENSION(NDOLD,NOLD) :: POPHELP, POPLTE_OLD

      INTEGER :: L, J, NFIRNA, NLANA, NA
      REAL :: SUM, POPJLOGL, ENTOTLOGL
      
      LOGICAL :: BDEPART, BTAUR, bUseENTOT
      
C***  The old popnumbers in POPHELP are interpolated with respect
C***  to the depth coordinate and stored in POPOLD (which then has
C***  still the old atomic levels, but the new radius grid)

C***  Using departure coefficients: replace old POPs by DEPARTs
      IF (BDEPART) THEN
         WRITE (0,*) 
     >      'Old DEPARTure coeficients used instead of POPNUMbers'
         DO L=1, NDOLD
            DO J=1, NOLD
               POPHELP(L,J) = POPHELP(L,J) / POPLTE_OLD(L,J)
            ENDDO
         ENDDO
         DO L=1, ND
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) / POPLTE(L,J)
            ENDDO
         ENDDO
      ENDIF


C***  Interpolation on Tau-Grid 
      IF (BTAUR) THEN
         WRITE (0,*) 'Interpolation of Popnumbers on Tau-Grid'
         DO L=1, ND
           DO J=1, NOLD
              IF (TAURCONT(L) .GT. TAURCONTOLD(NDOLD)) THEN
                 POPOLD(L,J) = POPHELP(NDOLD,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), TAURCONT(L), 
     >                    POPHELP(1,J), TAURCONTOLD, NDOLD)
              ENDIF 
           ENDDO
         ENDDO
      
      ELSEIF (bUseENTOT) THEN
        WRITE (0,*) 'Interpolation of Popnumbers over LOG density'
        DO J=1, NOLD
          !prepare necessary vectors
          DO L=1, NDOLD
            ENTOTOLDLOG(L) = LOG10(ENTOTOLD(L))
            POPHELPJLOG(L) = LOG10(MAX(POPHELP(L,J), POPMIN))
          ENDDO          
          !Perform interpolation on log(n_tot)
          dploop: DO L=1, ND
            ENTOTLOGL = LOG10(ENTOT(L))
            IF (ENTOTLOGL > ENTOTOLDLOG(NDOLD)) THEN
              !more dense than old innermost value => take old inner boundary value
              POPJLOGL = POPHELPJLOG(NDOLD)
            ELSEIF (ENTOTLOGL < ENTOTOLDLOG(1)) THEN
              !less dense than old outermost value => take old outer boundary value
              POPJLOGL = POPHELPJLOG(1)
            ELSE
              CALL SPLINPOX(POPJLOGL, ENTOTLOGL,
     >                     POPHELPJLOG, ENTOTOLDLOG, NDOLD)
            ENDIF
            POPOLD(L,J) = 10**(POPJLOGL)
          ENDDO dploop
        ENDDO              
      
      ELSE

C***  INTERPOLATION OF OLD POPNUMBERS TO THE NEW RADIUS GRID
        WRITE (0,*) 'Interpolation of Popnumbers on Radius-Grid'
        DO L=1, ND
           DO J=1, NOLD
              IF (RADIUS(L) .GT. ROLD(1)) THEN
                 POPOLD(L,J) = POPHELP(1,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), RADIUS(L), 
     >              POPHELP(1,J), ROLD, NDOLD)
              ENDIF
           ENDDO
         ENDDO

      ENDIF

C*****************************************************************
C***  Now the replacement of levels according to NTRANS
C*****************************************************************

C***  Loop over all depth points --------------------------------
      DO L=1, ND

C***  Copy old POPNUMs as assigned 
         DO J=1, N
            IF (NTRANS(J) .EQ. 0) THEN
               POPNUM(L,J) = POPMIN
            ELSE IF (NTRANS(J) .GT. 0 .AND. NTRANS(J) .LE. NOLD) THEN
               POPNUM(L,J) = POPOLD(L,NTRANS(J)) * ADPWEIGHT(J)
            ENDIF
         ENDDO

C***  If POPNUM actually contains DEPARTURE coeficients,
C***   convert them now back to POPNUMs
         IF (BDEPART) THEN
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) * POPLTE(L,J)
            ENDDO
         ENDIF

C***  Renormalization to the abundance of each element
 
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
         DO NA=1, NATOM
            SUM=0.0
            NFIRNA = NFIRST(NA)
            NLANA = NLAST(NA)
            DO J = NFIRNA, NLANA
               SUM = SUM + POPNUM(L,J)
            ENDDO
            SUM = SUM / ABXYZ(NA)
            IF (SUM /= 0.) THEN
               DO J=NFIRNA,NLANA
                  !POPMIN ensurance test (ansander, 2014)
                  IF (POPNUM(L,J) > POPMIN) THEN
                    POPNUM(L,J) = POPNUM(L,J) / SUM
                  ELSE
                    POPNUM(L,J) = POPMIN
                  ENDIF                   
               ENDDO
            ENDIF
         ENDDO
 
C***  Consistent electron density
         RNE(L) = .0
         DO J=1, N
           RNE(L) = RNE(L) + NCHARG(J) * POPNUM(L,J)
         ENDDO
 
      ENDDO   ! Deph points -------------------------------------! 

      RETURN
      END
