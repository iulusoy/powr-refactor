      SUBROUTINE POP_RENORM(POPNUM, ND, N, NATOM, NFIRST, NLAST, ABXYZ)
C***********************************************************************      
C***  Renormalization of POPNUMs 
C***
C***  called by LINPOP, ENSURETAUMAX, REGRIDOLDPOP
C***********************************************************************            
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, INTENT(IN) :: ND, N, NATOM
      
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ
      
      REAL, DIMENSION(ND, N), INTENT(INOUT) :: POPNUM
      
      REAL :: SUMME
      INTEGER :: L, J, NA, NFIRNA, NLANA
      
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
            DO J = NFIRNA,NLANA
              POPNUM(L,J) = POPNUM(L,J) / SUMME
            ENDDO
          ENDIF
        ENDDO

      ENDDO
      
      RETURN
      END
