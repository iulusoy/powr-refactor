      SUBROUTINE READMS(ICHANNEL, X, NDIM, NAME, IERR)
      IMPLICIT NONE
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************
      INTEGER ICHANNEL, NDIM, IERR, X, IDUMMY
      CHARACTER*8 NAME, NDUMMY

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, X, NDIM,
     >              'READ    ', IERR)

      RETURN
      END
