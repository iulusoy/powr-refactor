      SUBROUTINE OPENMS(ICHANNEL, IADR, MAXADR, IFNAME, IERR)
C************************************************************
C***  ROUTINE BY LARS KOESTERKE      8-Sep-1995 15:49:02
C************************************************************

      INTEGER :: IFNAME, ICHANNEL, IADR, MAXADR, IERR
      INTEGER :: DUMMY, IDUMMY
      CHARACTER(8) :: FNAME, CSTAT     !IMPORTANT: only the first 7 characters can be used!

      IF ((IFNAME /= 0) .AND. (IFNAME /= 1)) THEN
        WRITE(UNIT=FNAME, FMT='(A8)') IFNAME
      ELSE
        FNAME = '        '
      ENDIF
      CSTAT = 'AUTO    '
      DUMMY = 0
      IDUMMY = 0
      CALL CMSSTORE (ICHANNEL, IADR, MAXADR, CSTAT, FNAME,
     >              DUMMY, IDUMMY, 'OPEN    ', IERR)

      RETURN
      END
