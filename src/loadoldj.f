      SUBROUTINE LOADOLDJ(XJCold, XLAMBDAold, NDold, NFold, 
     >                    NDDIM, NFDIM, hOldMODEL)
C***********************************************************************
C***  read out the XJC radiation field from an old model
C***
C***  called from WRSTART
C***********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDold, NDDIM, NFDIM, hOldMODEL
      INTEGER, INTENT(INOUT) :: NFold
      REAL, DIMENSION(NFDIM) :: XLAMBDAold
      REAL, DIMENSION(NDDIM * NFDIM) :: XJCold

      CHARACTER(8) :: NAME
      INTEGER :: K, IERR


      CALL READMS (hOldMODEL, NFold     ,     1, 'NF      ', IERR)
      CALL READMS (hOldMODEL, XLAMBDAold, NFold, 'XLAMBDA ', IERR)
      DO K=1, NFold
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL READMS(hOldMODEL, XJCold(1+NDold*(K-1)), NDold, NAME, IERR)
      ENDDO


      RETURN
 
      END

