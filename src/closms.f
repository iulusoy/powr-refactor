      SUBROUTINE CLOSMS (ICHANNEL, IERR)
C************************************************************
C***  ROUTINE VON LARS KOESTERKE      8-Sep-1995 15:50:47
C************************************************************
      INTEGER ICHANNEL, IERR, IDUMMY
      CHARACTER*8 NDUMMY
      REAL X

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NDUMMY, NDUMMY, X,
     >              IDUMMY, 'CLOSE    ', IERR)

      RETURN
      END
