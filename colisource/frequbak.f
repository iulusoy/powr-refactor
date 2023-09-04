      SUBROUTINE FREQUBAK(NATOM, ND, OPAKELEM, OPAKOLDELEM)
C***********************************************************************
C***  copies multidimensional arrays into "OLD"-variables before next
C***   fine frequency loop cycle (needed e.g. for FREQUINT)
C***  
C***  called by COLI
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOM, ND

      REAL, DIMENSION(NATOM, ND) :: OPAKELEM, OPAKOLDELEM

      INTEGER :: L, NA

      DO NA=1, NATOM
        DO L=1, ND
          OPAKOLDELEM(NA, L) = OPAKELEM(NA, L)
        ENDDO
      ENDDO

      RETURN

      END
