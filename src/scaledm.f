      SUBROUTINE SCALEDM(DM, EN, N, bINV)

      IMPLICIT NONE
  
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: bINV

      REAL, DIMENSION(N, N) :: DM
      REAL, DIMENSION(N) :: EN

      INTEGER :: I, J

      DO I=1, N
        DO J=1, N
          IF (bINV) THEN
            DM(I,J) = DM(I,J) / EN(J)
          ELSE
            DM(I,J) = DM(I,J) * EN(J)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
    
