      SUBROUTINE POPMIN_NULLING (ZERO_RATES, POPNUM, POPMIN, ND, N)
C***********************************************************************
C***  Sets all POPNUMS that were flagged by ZERO_RATES to zero 
C***  Reason: these level populations have been set to POPMIN in steal 
C***      but should be set to 0.0 in all radiative-transfer programs
C***      in oder to avoid any artificial contributions to the 
C***      emissivities and opacities 
C***  Called from: WRCONT, COMO, COLI, FORMAL  
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND, N
      REAL, INTENT(IN) :: POPMIN

      REAL,    DIMENSION(ND, N) :: POPNUM
      LOGICAL, DIMENSION(N, ND) :: ZERO_RATES
    
      INTEGER :: L, J

      DO L=1, ND
        DO J=1, N
          IF (ZERO_RATES(J,L) .OR. POPNUM(L,J) < 1.1 * POPMIN) THEN
            POPNUM(L,J) = .0
          ENDIF
        ENDDO        
      ENDDO

      RETURN
      END
