      SUBROUTINE FREQUBAKION(ND, NATOM, NION, ARRKION, ARRKOLDION)
C***********************************************************************
C***  copies threedimensional arrays into "OLD"-variables before next
C***   fine frequency loop cycle (needed e.g. for FREQUINT)
C***  
C***  called by COLI
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NION, NATOM, ND

      REAL, DIMENSION(ND, NATOM, NION) :: ARRKION, ARRKOLDION

      INTEGER :: L, NA, ION

      DO L=1, ND
        DO NA=1, NATOM
          DO ION=1, NION
            ARRKOLDION(L, NA, ION) = ARRKION(L, NA, ION)
          ENDDO
        ENDDO
      ENDDO

      RETURN

      END
