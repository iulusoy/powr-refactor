      PROGRAM MAINmodify 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:26 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL modify 
      END
      SUBROUTINE ADDHISTENTRY(MODHIST,LAST,MAXHIST,ENTRYLEN,ENTRYSTR)
C***********************************************************************
C***  ADDS A STRING ENTRY TO THE MODEL HISTORY CHARACTER ARRAY
C     ENTRYSTR: new string to add
C     ENTRYLEN: length of ENTRYSTR
C     
C     Note: This routine updates MODHIST and the historical LAST parameter
C           to be compartible with older code that still uses ENDCODE/DECODE 
C           and declares MODHIST as an integer array which is filled with
C           Hollerith constants
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAXHIST, ENTRYLEN      
      CHARACTER(ENTRYLEN), INTENT(IN) :: ENTRYSTR

      INTEGER, INTENT(INOUT) :: LAST
      CHARACTER(8*MAXHIST), INTENT(INOUT) :: MODHIST

      INTEGER :: LASTCHAR, INFROM, INTO, BUFFERINT, CURLAST
      REAL :: ADDBYTES
      CHARACTER(8) :: BUFFER8


      !LAST can be set to -1 => read from MODHIST(1:8)
      IF (LAST < 0) THEN
        BUFFER8 = MODHIST(1:8)
        READ(UNIT=BUFFER8, FMT='(A8)') BUFFERINT
        CURLAST = BUFFERINT
      ELSE
        CURLAST = LAST
      ENDIF

      LASTCHAR = CURLAST * 8
      
      INFROM = LASTCHAR + 1
      INTO = LASTCHAR + ENTRYLEN
      MODHIST(INFROM:INTO) = ENTRYSTR

      ADDBYTES = REAL(ENTRYLEN) / 8.
      CURLAST = CURLAST + INT(ADDBYTES)
      IF (ADDBYTES - REAL(INT(ADDBYTES)) > 0) THEN
        !Entry has a length that is not a multiple of 8, one more byte needed)
        CURLAST = CURLAST + 1
      ENDIF

      !MODHIST(1:8) or MODHIST(1) in integer array definition contains the 
      !currently used length of MODHIST (in Bytes i.e. in CHARs / 8)
      ! The first bytes is for historical written as an integer, 
      !  NOT as a character containing an integer number
      !  THerefore (A8) is used als FORMAT instead of (I8) 
      WRITE(UNIT=BUFFER8, FMT='(A8)') CURLAST          
      MODHIST(1:8)=BUFFER8

      IF (LAST >= 0) THEN
        !Fill LAST only if not called with -1 (otherwise COLI will crash on subroutine call)
        LAST = CURLAST
      ENDIF

      RETURN
      END
      SUBROUTINE APPEND_AUTOLEVELS (N, N_WITH_DRLEVELS, NDIM, MAXIND, 
     $                  MAXAUTO, NAUTO, LOWAUTO, IONAUTO, EAUTO, ELEVEL, 
     $                  LEVEL, EION, WEIGHT, INDLOW, INDNUP, LASTIND,
     >                  LEVUPAUTO, LEVAUTO, WAUTO, 
     >                  NCHARG, IONGRND, NOM)
C******************************************************************************
C***  This subroutine appends autoionization levels 
C***  to the vectors LEVEL, WEIGHT and ELEVEL
C***  Index range: N+1 ... N_WITH_DRLEVELS
C***  THIS SUBROUTINE IS ESSENTIAL TO BE CALLED FROM COLI, 
C**   but also called from STEAL - PRIDAT in order to displaying 
C**   the autoionization levels ("DR-LEVELS") with PRINT DATOM 
C******************************************************************************

      DIMENSION WEIGHT(NDIM), ELEVEL(NDIM), EION(NDIM)
      DIMENSION NCHARG(NDIM), IONGRND(NDIM), NOM(NDIM)
      DIMENSION INDNUP(MAXIND),INDLOW(MAXIND)
      DIMENSION LOWAUTO(MAXAUTO),IONAUTO(MAXAUTO), WAUTO(MAXAUTO)
      DIMENSION EAUTO(MAXAUTO)
      CHARACTER*10 LEVEL(NDIM)
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

C***  Copy vector LEVUPAUTO to LEVAUTO, skipping multiple occurrences 
      NLEVEL_AUTO = 1
      LEVAUTO(1) = LEVUPAUTO(1)
      DO 5 ILEVEL=2, NAUTO
         DO ITEST=1, NLEVEL_AUTO 
            IF (LEVUPAUTO(ILEVEL) .EQ. LEVAUTO(ITEST)) GOTO 5
         ENDDO
         NLEVEL_AUTO = NLEVEL_AUTO + 1
         LEVAUTO(NLEVEL_AUTO) = LEVUPAUTO(ILEVEL)
    5 CONTINUE

C***  Number of levels is increased
      N_WITH_DRLEVELS = N + NLEVEL_AUTO
      IF (N_WITH_DRLEVELS .GT. NDIM) THEN
         WRITE (*,'(A, I5,A,I5)') 
     >    '*** NDIM=', NDIM, ' insufficient: needed', N_WITH_DRLEVELS
      ENDIF

C***  The new AUTO-Levels are appended to LEVEL list
      DO ILEVEL = 1, NLEVEL_AUTO
C***     ... the list of levels
         LEVEL(N+ILEVEL) = LEVAUTO(ILEVEL)
      ENDDO

C***  further level attributes are same as for lower level 
C***    of DRTRANSIT with same upper level name
      DO ILEVEL= N+1, N+NLEVEL_AUTO
         DO INDDR = 1, NAUTO
            IF (LEVEL(ILEVEL) .EQ. LEVUPAUTO(INDDR)) THEN
              NCHARG (ILEVEL) = NCHARG (LOWAUTO(INDDR))
              EION   (ILEVEL) = EION   (LOWAUTO(INDDR))
              IONGRND(ILEVEL) = IONGRND(LOWAUTO(INDDR))
              NOM    (ILEVEL) = NOM    (LOWAUTO(INDDR))

              WEIGHT (ILEVEL) = WAUTO (INDDR)
              ELEVEL (ILEVEL) = EAUTO (INDDR)

              GOTO 10
            ENDIF
         ENDDO
         STOP '*** FATAL INTERNAL ERROR IN subr. append_autolevels'
   10    CONTINUE 
      ENDDO

C***  NEW!! for stabilizing transitions, INDNUP points to the
C***      auto-ionizing level, and NOT to the next-higher ion
      DO 11 INDDR = 1, NAUTO
         DO ILEVEL= N+1, N+NLEVEL_AUTO
            IF (LEVEL(ILEVEL) .EQ. LEVUPAUTO(INDDR)) THEN
               IND = INDDR + LASTIND
               INDNUP (IND) = ILEVEL
               GOTO 11
            ENDIF
         ENDDO
         STOP '*** FATAL INTERNAL ERROR2 IN subr. append_autolevels'
   11 CONTINUE

C***  Energies of the appended auto-ionizing levels now relative 
C***   to the ground level, as for "normal" levels
      DO ILEVEL = N+1, N_WITH_DRLEVELS
         ELEVEL(ILEVEL) = ELEVEL(ILEVEL) + EION(ILEVEL)
      ENDDO

      RETURN
      END 
      SUBROUTINE CHANGE (NAMEOLD, NAMENEW, ICHANNEL)
C***  THIS SUBROUTINE REPLACES A NAME IN THE NAME-INDEX ARRAY IADR
C***  OF A MASS-STORAGE FILE

        CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAMEOLD, NAMENEW, 
     >                 DUMMY, NDUMMY, 'CHANGE', IERR)

      RETURN
      END
      SUBROUTINE CLOCK

      STOP 'CLOCK NOT IMPLEMENTED AT DEC/OSF'

      RETURN
      END
      SUBROUTINE CLOSMS (ICHANNEL, IERR)
C************************************************************
C***  ROUTINE VON LARS KOESTERKE      8-Sep-1995 15:50:47
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, IDUMMY, IDUMMY, DUMMY, 
     >              IDUMMY, 'CLOSE', IERR)

      RETURN
      END
      SUBROUTINE CMSSTORE (ICHANNEL, IADRDUMMY, MAXADRDUMMY, 
     >                     NAME, NAME2, X, NDIM, ACTION, IERR)
C*******************************************************************
C***
C***  Cray-MS-STORagE
C***
C***  ROUTINE IS THE ADAPTER BETWEEN THE ROUTINE STORAGE, WHICH
C***    EMULATES MASS-STORAGE, AND THE MS-CALL AT CRAY
C***  IT IS NEEDED BECAUSE AT THE CRAY THE INDEX ARRAYS ARE NOT
C***    GIVEN BY THE CALL OF READMS, WRITMS AND CLOSMS
C***  Version 1.0 :
C***                The actual file is immedeatly closed when a second
C***                file is used
C***  Version 2.0 :
C***                unchanged
C***  Version 3.0 :
C***                The Array IADR is now twodimensional. The actual
C***                file need not to be closed when a new file is used
C***  Version 3.1 :
C***                MAXADR increased from 25600 to 256000 (2000 Records)
C***                wrh 15-Mar-2005 16:33:59
C***  Version 3.2 :
C***                MSMAXCH increased from 3 to 5
C*******************************************************************

      PARAMETER (MSMAXCH = 5)
C***  NOTE: MAXADR MUST BE A MULTIPLE OF 128
C***        (I.E. IADRL IN ROUTINE STORAGE)
      PARAMETER (MAXADR = 2000 * 128)

      CHARACTER*8 STATUS, ACTION
      DIMENSION IFILE(MSMAXCH), IADR(MAXADR,MSMAXCH)

C      SAVE LASTCH, IFILE, NFILE, STATUS, ITRANS
      SAVE

C***  CHECK THE ICHANNEL NUMBER
      IF (ICHANNEL .LE. 0) THEN
        WRITE (0,*) ' NEGATIV ICHANNEL NUMBERS ARE NOT ALLOWED'
        STOP 'ERROR IN CMSSTORE'
      ENDIF

      IF (ICHANNEL .NE. LASTCH) THEN
C***  SEARCH FOR FILE NUMBER
        STATUS = 'UNKNOWN'
        DO I=1, NFILE
          IF (ICHANNEL .EQ. IFILE(I)) THEN
            ITRANS = I
            STATUS = 'KNOWN'
            GOTO 10
          ENDIF
        ENDDO
   10   CONTINUE
C***  CLOSE THE FILE USED LAST. This is now (Vers 3) done by SCLOSE
        IF (LASTCH .GT. 0) THEN
          CALL STORAGE (LASTCH, IADR(1,ITRANSLAST), 
     >                  MAXADR, 'DUMMY', 'DUMMY', X, NDIM, 
     >                  'SCLOSE', 'CRAY', IERR)
        ENDIF
C***  OPEN THE ACTUAL FILE. This is now (Vers 3) done by SOPEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), 
     >                  MAXADR, 'DUMMY', 'DUMMY', 
     >                  X, NDIM, 'SOPEN', 'CRAY', IERR)
        ENDIF
      ENDIF

      IF (ACTION(1:4) .EQ. 'OPEN') THEN
        IF (STATUS .NE. 'UNKNOWN') THEN
          WRITE (0,*) 'DO NOT OPEN AN OPEN FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF
        NFILE = NFILE + 1
        ITRANS = NFILE
        IF (NFILE .GT. MSMAXCH) THEN
          WRITE (0,*) 'DO NOT OPEN MORE FILES THAN MSMAXCH'
          STOP 'ERROR IN CMSSTORE'
        ENDIF
        CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                NAME, NAME2, X, NDIM, 
     >                'OPEN', 'CRAY', IERR)
        IFILE(NFILE) = ICHANNEL
        STATUS = 'KNOWN'

      ELSE IF (ACTION(1:4) .EQ. 'READ') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'READ', 'CRAY', IERR)
        ELSE
          WRITE (0,'(A,I3)') 
     >        'DO NOT READ IN A CLOSED FILE, CHANNEL=', ICHANNEL
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:6) .EQ. 'LENGTH') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'LENGTH', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT READ (LENGTH) IN A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:5) .EQ. 'WRITE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'WRITE', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT WRITE IN A CLOSED FILE'
          WRITE (0,*) 'ICHANNEL=',ICHANNEL
          WRITE (0,'(A5,A8)') 'NAME=',NAME
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:5) .EQ. 'CLOSE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'CLOSE', 'CRAY', IERR)
          STATUS = 'UNKNOWN'
          NFILE = NFILE - 1
          DO I=ITRANS, NFILE
            IFILE(I) = IFILE(I+1)
            DO J=1, MAXADR
              IADR(J,I) = IADR(J,I+1)
            ENDDO
          ENDDO
        ELSE
          WRITE (0,*) 'DO NOT CLOSE A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:6) .EQ. 'CHANGE') THEN
        IF (STATUS .EQ. 'KNOWN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  'CHANGE', 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT CHANGE IN A CLOSED FILE'
          WRITE (0,*) 'ICHANNEL=',ICHANNEL
          WRITE (0,'(A5,A8)') 'NAME=',NAME
          WRITE (0,'(A5,A8)') 'NAME2=',NAME2
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE IF (ACTION(1:4) .EQ. 'INFO') THEN
        IF (STATUS .EQ. 'OPEN') THEN
          CALL STORAGE (ICHANNEL, IADR(1,ITRANS), MAXADR, 
     >                  NAME, NAME2, X, NDIM, 
     >                  ACTION, 'CRAY', IERR)
        ELSE
          WRITE (0,*) 'DO NOT DO INFO ON A CLOSED FILE'
          STOP 'ERROR IN CMSSTORE'
        ENDIF

      ELSE
        WRITE (0,*) ' ACTION ', ACTION( :IDX(ACTION)), ' NOT KNOWN'
        STOP 'ERROR IN CMSSTORE'

      ENDIF

      IF (ACTION(1:5) .NE. 'CLOSE') THEN
        LASTCH = ICHANNEL
        ITRANSLAST = ITRANS
      ELSE
        LASTCH = 0
      ENDIF

      RETURN
      END
      SUBROUTINE COUNT(J, ISTR)
C **********************************************************************
C *** Writes a 3-digit integer number with leading zeros (001, 002 ...)
C ***  CALLED BY FEDAT
C **********************************************************************

      CHARACTER*3 ISTR

ccc   test: restore the error of Leindecker's version!
ccc      if (j .eq. 100) then
ccc         istr = '099'
ccc         return
ccc      endif 

C *** COUNTER FOR LEVEL NAMES
      IF (J .LT. 0 .OR. J .GT. 999) THEN
         STOP 'FATAL ERROR IN SUBR. COUNT'
      ELSE IF (J .LT. 10) THEN
         WRITE (ISTR,'(A2,I1)') '00', J
      ELSE IF (J .LT. 100) THEN
         WRITE (ISTR,'(A1,I2)') '0', J
      ELSE
         WRITE (ISTR,'(I3)') J
      ENDIF
    
      RETURN
    
      END
      SUBROUTINE DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,SIGMATHK,
     $                  SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND, KEYCBF,
     >                  ROUTINE, INDEXMAX, NFEREADMAX, MAXFEIND,
     >                  LASTFE, SIGMAFE, INDRB, INDRF,
     >                  IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY, 
     >                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >                  LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

c!!!!!! Folgende parameter wurden entfernt:
C!!!    CBFC, BOUND, EINSTINT, COCOFE, NCOMAX, NCO
c!!!    folgende Parameter sind neu: MAXFEIND, FEDUMMY
C!!!    umbenannt wurden: NMAX -> NFEREADMAX, DUMMY -> FEDUMMY


C*******************************************************************************
C***  READS ATOMIC DATA FROM TAPE4=DATOM  **************************************
C***  The decoded elements are indexed in sequence of their occurrence
C***   (J = 1 ... NATOM)
C***  The meaning of vector KODAT is weird:
C***  Each chemical element is assigned a position in this vector
C***    (originally by an arbitrary definition; 
C***    in the version from March 2007, the index is equal to 
C***    the corecharge (NZ, atomic number, Kernladungszahl)
C***  KODAT(NZ) contains the index J under which this element was found
C***    in the DATOM file; unused elements have KODAT(NZ)=0 
C*******************************************************************************
 
      INTEGER, INTENT(IN) :: NDIM, MAXIND, MAXION, MAXKONT, MAXAUTO

      INTEGER, DIMENSION(NDIM) :: NCHARG, IONGRND, MAINQN, NOM
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,ADDCON1,ADDCON2,ADDCON3
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      DIMENSION LOWAUTO(MAXAUTO),WAUTO(MAXAUTO),EAUTO(MAXAUTO)
     $         ,AAUTO(MAXAUTO),IONAUTO(MAXAUTO),KRUDAUT(MAXAUTO)
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(NDIM),LEVUP,LEVLOW, LEVION
      CHARACTER*10 ELEMENT(MAXATOM),NEWELE
      CHARACTER*8 IGLOW, ICBF, IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      CHARACTER*4 CEY,KEYCBB(MAXIND)
      CHARACTER*3 KRUDI,DRRUDI
      CHARACTER*2 SYMBOL(MAXATOM), KSHELLSYM
      CHARACTER(LEN=*), INTENT(IN) :: ROUTINE

      LOGICAL :: BFEMODEL
 
      DO 15 NA=1,MAXATOM
         DO ISTAGE=1, MAXION
           SIGMATHK(NA,ISTAGE)=.0
           SEXPOK  (NA,ISTAGE)=.0
           EDGEK   (NA,ISTAGE)=.0
         ENDDO
   15 KODAT(NA)=0
      DO 6 I=1,NDIM
      EION(I)=.0
      IONGRND(I)=-1
      ALTESUM(1,I)=-1.
C***  INITIALIZE TRANSTION MATRIX TO DETECT MISSING LINE TRANSITIONS
      DO 6 J=1,NDIM
    6 EINST(I,J)=-99.
    
      DO IND=1,MAXIND
       INDNUP(IND)=0
       INDLOW(IND)=0
       COCO(1,IND)=.0
       COCO(2,IND)=.0
       COCO(3,IND)=.0
       COCO(4,IND)=.0
       KEYCBB(IND)='    '
      ENDDO

      DO KON=1,MAXKONT
       KONTNUP(KON)=-1
       KONTLOW(KON)=0
       IGAUNT(KON)= ' '
      ENDDO

      DO 96 J=1,MAXAUTO
      WAUTO(J)=0.0
      EAUTO(J)=0.0
      AAUTO(J)=0.0
      IONAUTO(J)=-1
      KRUDAUT(J)=0
   96 LOWAUTO(J)=0
      NATOM=0
      N=0
      IND=0
      KONT=0
      NAUTO=0
      LEVSEQ=0

      BFEMODEL = .FALSE.
 
      OPEN (4, FILE='DATOM', STATUS='OLD')
    1 READ(4,2,END=3) KARTE
    2 FORMAT(A)

      IF (KARTE(:1) .EQ. '*' .OR. KARTE(:1) .EQ. ' ') GOTO 1
      IF (KARTE(:10) .EQ. 'ELEMENT   ' ) GOTO 5
      IF (KARTE(:10) .EQ. 'LEVEL     ' ) GOTO 10
      IF (KARTE(:10) .EQ. 'LINE      ' ) GOTO 20
      IF (KARTE(:10) .EQ. 'CONTINUUM ' ) GOTO 30
      IF (KARTE(:10) .EQ. 'K-SHELL   ' ) GOTO 40
      IF (KARTE(:10) .EQ. 'LTESUM    ' ) GOTO 50
      IF (KARTE(:10) .EQ. 'DRTRANSIT ' ) GOTO 60
      CALL REMARK ('UNRECOGNIZED DATA INPUT')
      GOTO 990
      
C***  ELEMENTS ---------------------------------------------------------
    5 CONTINUE
C***  DECODED ELEMENT IS ALREADY KNOWN
      NEWELE=KARTE(13:22)
      DO 19 NA=1,NATOM
      IF (NEWELE .EQ. ELEMENT(NA)) GOTO 1
   19 CONTINUE

C***  If DATOM was called with parameter ROUTINE = 'NOIRON, 
C****   the 'ELEMENT GENERIC' card is ignored. 
C***    This feature was introduced to facilitate 
C***    Sonja's NEWFORMAL_CARDS program
      IF (ROUTINE(1:6) == 'NOIRON' .AND. NEWELE == 'GENERIC') GOTO 1

C***  NEW ELEMENT DECODED:
      LEVSEQ=0
      NATOM=NATOM+1
      IF (NATOM .GT. MAXATOM) THEN
         CALL REMARK ('DATOM: MORE ELEMENTS THAN DIMENSIONED')
         GOTO 990
      ENDIF
      READ (KARTE,9) ELEMENT(NATOM),SYMBOL(NATOM),ATMASS(NATOM),
     $                    STAGE(NATOM)
    9 FORMAT (12X,A10,2X,A2,4X,F6.2,3X,F5.0)

      CALL FINDCHARGE (ELEMENT(NATOM), NZ)
      KODAT(NZ) = NATOM      

      IF (NZ .EQ. 0) THEN
         WRITE (0,*) 'UNKNOWN ELEMENT DECODED: ', ELEMENT(NATOM)
         GOTO 990
      ENDIF

C***  "GENERIC" MODEL ATOM OF IRON GROUP ELEMENTS DECODED
      IF  (NZ .EQ. 26) THEN
         BFEMODEL = .TRUE.
C***    DECODE INPUT CARD AGAIN
 109     FORMAT (12X,A10,2X,A2,4X,I6,2X,I6)
         READ (KARTE,109) ELEMENT(NATOM),SYMBOL(NATOM),IONLOW,IONTOP
C***    COMPLETE IRON-DATA IS READ IN FROM MASS STORAGE FILE 'FEDAT'
         NFIRSTFE = N + 1

C***     Iron line indices are arranged *behind* the DRTRANSITs
         LASTINDAUTO = LASTIND + NAUTO

         CALL FEDAT (ROUTINE, INDEXMAX, NFEREADMAX, IONLOW, IONTOP,
     &               MAXATOM, NDIM, MAXIND, MAXKONT, NATOM,
     &               N, LASTFE, LASTKON, LASTINDAUTO, MAXFEIND,
     &               EINST, SIGMAFE, INDRB, INDRF, IFENUP, 
     &               IFELOW, INDNUP, INDLOW, KONTNUP, KONTLOW,
     &               LEVEL, ELEMENT, SYMBOL, ATMASS, STAGE,
     &               ELEVEL, WEIGHT, EION, NCHARG, NOM, KODAT,
     &               NFIRST, NLAST, IFRBSTA, IFRBEND, FEDUMMY,
     &               VDOPFE, DXFE, XLAM0FE, SIGMAINT, KEYCBB)
C***  Fill vector with ionization energies
         DO I=NFIRSTFE+1, N-1
            IF (EION(I) .EQ. .0) EION(I) = EION(I-1)
         ENDDO
      ENDIF
      GOTO 1
 
C***  LEVELS -----------------------------------------------------------
   10 N=N+1
      IF (LEVSEQ .NE. 0) THEN
         CALL REMARK ('DATOM: LEVEL CARD OUT OF SEQUENCE')
         GOTO 990
      ENDIF
      IF(N .GT. NDIM) THEN
         CALL REMARK ('DATOM : MORE LEVELS THEN DIMENSIONED (NDIM)')
         GOTO 990
      ENDIF
      IF (NATOM .NE. 0) NOM(N)=NATOM

      READ (KARTE,11,ERR=985) LEVEL(N),NCHARG(N),NW,ELEVEL(N),E,MAINQN(N)
   11 FORMAT(12X,A10,1X,I2,1X,I4,2F10.0,1X,I2)

      WEIGHT(N)=FLOAT(NW)

C***  If EION is empty, it will be repeated from last level
C***  after checking that the latter belongs to the same element and ion
      IF (E .EQ. 0.0) THEN
         IF (N .GT. 1) THEN    
            IF (NOM   (N-1) .NE. NOM   (N) .OR. 
     >          NCHARG(N-1) .NE. NCHARG(N) ) THEN
C***            Setting missing ionization energy to zero;
C***              this may be OK, if it is the highest level of an element
C***              and will be checked later (in the continuum block)
                  EION(N) = .0
            ELSE
                EION(N) = EION(N-1)
            ENDIF
         ELSE
            CALL REMARK ('ERROR: FIRST LEVEL WITHOUT IONIZATION ENERGY')
            GOTO 990
         ENDIF 
      ELSE
         EION(N) = E
C***  If EION is given redundantly, it must be identical  
         IF (N .GT. 1) THEN    
            IF (NOM   (N-1) .EQ. NOM   (N) .AND. 
     >          NCHARG(N-1) .EQ. NCHARG(N) .AND.
     >          EION  (N-1) .NE. E ) THEN
                  CALL REMARK ('ERROR: DIFFERENT IONIZATION ENERGIES')
                  GOTO 990
            ENDIF
         ENDIF
      ENDIF

C***  If level energy entry is empty, it is calculated from MainQN 
C***  by RYDBERG's formula
C***  If both are not given, this may be a mistake. Unfortunately,
C***  old DATOM files use these empty files for He III, H II as meaning
C***  ELEVEL = 0.0 -- therefore no error check can be made here
      IF (KARTE(31:40) .EQ. '          ') THEN 
         IF (KARTE(52:53) .NE. '  ') THEN
            F=FLOAT(MAINQN(N))
            ELEVEL(N) = (1.-1./F/F)*EION(N)
         ELSE
            ELEVEL(N) = .0
         ENDIF
      ENDIF
      GOTO 1
 
C***  LINE TRANSITIONS  ------------------------------------------------
   20 READ (KARTE,21) LEVUP,LEVLOW,AUPLOW,KRUDI,CEY,CO1,CO2,CO3,CO4
   21 FORMAT(10X,A10,2X,A10,G10.0,2X,A3,A4,1X,4G7.0)
      LEVSEQ=1
C***  FIND UPPER INDEX
      DO 22 J=1,N
      NUP=J
      IF (LEVEL(J).EQ.LEVUP ) GOTO 23
   22 CONTINUE
      CALL REMARK ('UPPER LINE LEVEL NOT FOUND')
      GOTO 990
C***  FIND LOWER INDEX
   23 DO 24 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 25
   24 CONTINUE
      CALL REMARK ('LOWER LINE LEVEL NOT FOUND')
      WRITE (0,*) 'LOWER LEVEL = ',LEVLOW
      GOTO 990
   25 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK ('ERROR: LINE BETWEEN DIFFERENT ELEMENTS')
            GOTO 990
         ENDIF
         IF (EION(LOW) .EQ. .0) THEN
            CALL REMARK ('ERROR: MISSING IONIZATION ENERGY')
            GOTO 990
         ENDIF

      ENDIF
      IF (NCHARG(NUP) .NE. NCHARG(LOW)) THEN
         CALL REMARK ('LINE BETWEEN DIFFERENT IONIZATION STAGES')
         GOTO 990
      ENDIF
      IF (NUP.LE.LOW) THEN
         CALL REMARK ('LINE TRANSITION INDICES WRONG')
         GOTO 990
      ENDIF
      IF (ELEVEL(NUP) < ELEVEL(LOW)) THEN
         CALL REMARK ('LINE ERROR: UPPER LEVEL NOT IN FRONT')
         GOTO 990
      ENDIF
C***  CORRECT LINE TRANSITION DETECTED:
      IND=IND+1

C***  ERROR STOP
      IF (IND .GT. MAXIND) THEN
         CALL REMARK ('ERROR: IND .GT. MAXIND')
         GOTO 990
      ENDIF
      INDNUP(IND)=NUP
      INDLOW(IND)=LOW
      EINST(NUP,LOW) = AUPLOW
      KEYCBB(IND)=CEY
      COCO(1,IND)=CO1
      COCO(2,IND)=CO2
      COCO(3,IND)=CO3
      COCO(4,IND)=CO4
C***  RUDIMENTAL TRANSITIONS ARE MARKED BY -2. IN THE TRANSPOSED
C***  MATRIX ELEMENT  EINST(LOW,NUP)
      IF (KRUDI.NE.'   ') EINST(LOW,NUP)=-2.
      LASTIND=IND
      GOTO 1
 
C***  CONTINUUM TRANSITIONS    -----------------------------------------
   30 READ (KARTE,31,ERR=990) 
     >    LEVLOW, SIGMA, ALPLOW, SLOW, IGLOW, ICBF, LEVION
   31 FORMAT (10X,A10,3G10.0,1X,A8,1X,A8,1X,A10)
      LEVSEQ=1
C***  FIND LOWER INDEX
      DO 34 J=1,N
      LOW=J
      IF (LEVEL(J) .EQ. LEVLOW ) GOTO 35
   34 CONTINUE
      CALL REMARK ('LOWER CONTINUUM LEVEL NOT FOUND')
      GOTO 990
   35 CONTINUE
C***  FIND UPPER INDEX
      IF (LEVION .EQ. '          ') THEN
          EMIN=999999.
          NUP=0
          DO 39 J=1, N
          IF ((NATOM .GT. 1) .AND. (NOM(J) .NE. NOM(LOW))) GOTO 39
          IF ((NCHARG(J) .EQ. NCHARG(LOW)+1).AND.(ELEVEL(J) .LT. EMIN)) 
     $       THEN
             NUP=J
             EMIN=ELEVEL(J)
          ENDIF 
   39     CONTINUE
          IF (NUP .NE. 0) GOTO 33
      ELSE
          DO 32   J=1, N
          NUP=J
          IF ((LEVEL(J) .EQ. LEVION).AND.(NCHARG(J) .EQ. NCHARG(LOW)+1))
     $        GOTO 33
   32     CONTINUE
      ENDIF

      WRITE (0,*) 'ERROR: UPPER CONTINUUM LEVEL ', LEVION, 
     $         ' NOT FOUND'
      GOTO 990

   33 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK ('CONTINUUM BETWEEN DIFFERENT ELEMENTS')
            GOTO 990
         ENDIF
      ENDIF
C***  CORRECT CONTINUUM TRANSITION DETECTED:
      KONT=KONT+1
C***  ERROR STOP
      IF (KONT .GT. MAXKONT) THEN
         CALL REMARK ('ERROR: MORE CONTINUA THAN DIMENSIONED (MAXKONT)')
         GOTO 990
      ENDIF
      KONTNUP(KONT)=NUP
      KONTLOW(KONT)=LOW
      EINST(LOW,NUP)=SIGMA
      ALPHA(KONT)=ALPLOW
      SEXPO(KONT)=SLOW
      IGAUNT(KONT)=IGLOW
C***  KEYCBF(KONT)=ICBF
      IF (IGLOW .EQ. 'PIKB12' .OR. IGLOW .EQ. 'OPAPROIX') THEN
C*** READ FURTHER INFORMATION (ADDCON1-3) IN FOLLOWING LINE
  121   READ(4,122,END=3) KARTE
  122   FORMAT(A)
  130   READ (KARTE,131) ADD1, ADD2, ADD3
  131   FORMAT (21X, 3G10.0)
      ELSE
        ADD1 = 0.
        ADD2 = 0.
        ADD3 = 0.
      ENDIF
      ADDCON1(KONT) = ADD1
      ADDCON2(KONT) = ADD2
      ADDCON3(KONT) = ADD3
      LASTKON=KONT
      GOTO 1
 
C***  K-SHELL IONISATION
   40 CONTINUE
C***  DOES AKTUELL ELEMENT FIT TO K-SHELL-DATA ?
      IF (SYMBOL(NATOM) .NE. KARTE(16:17)) GOTO 982

C***  Are the K-shell data split for ionization stages?
      IF (KARTE(18:20) .NE. '   ') THEN
         READ (KARTE(18:20),'(I3)',ERR=980) ISTAGE 
         IF (ISTAGE .LT. 1) GOTO 981
         IF (ISTAGE .GT. NZ-2) GOTO 981
         READ (KARTE,41,ERR=983) SIGMATHK(NATOM,ISTAGE), 
     >                           SEXPOK  (NATOM,ISTAGE)
   41    FORMAT (20X,2G10.0)
C***     Special option: energy might be given in electron Volts
         IF (KARTE(52:53) .EQ. 'EV') THEN
            READ (KARTE(44:51), '(F8.0)', ERR=983) EVOLT
            EDGEK(NATOM,ISTAGE) = EVOLT * 1.E8 / 12397.7
         ELSE 
            READ (KARTE(44:53), '(F10.0)', ERR=983) EDGEK(NATOM,ISTAGE)
         ENDIF
      ELSE
         READ (KARTE,42,ERR=983) SIGMATHK(NATOM,1), 
     >                           SEXPOK  (NATOM,1),
     >                           EDGEK   (NATOM,1)
   42    FORMAT (20X,2G10.0,3X,F10.0)
         DO ISTAGE=2, NZ-2
            SIGMATHK(NATOM,ISTAGE) = SIGMATHK(NATOM,1)
            SEXPOK  (NATOM,ISTAGE) = SEXPOK  (NATOM,1)
            EDGEK   (NATOM,ISTAGE) = EDGEK   (NATOM,1)
         ENDDO
      ENDIF

      LEVSEQ=1
      GOTO 1

C***  SUM OF TRANSITIONS TO UPPER LEVELS WHICH ARE ASSUMED TO BE IN LTE
   50 CONTINUE
      WRITE (0,*) 'ERROR: the LTESUM branch is no longer supported' 
      WRITE (0,*) 'ERROR: the LTESUM branch is no longer supported' 
      STOP 'FATAL ERROR reported from subr. DATOM'

c      READ (KARTE,51) LEVLOW,IRANGE,ASUM,COEFF1,COEFF2
c   51 FORMAT (10X,A10,1X,A8,1X,G9.0,1X,F7.0,1X,F7.0)
c      LEVSEQ=1
cC***  FIND LOWER INDEX
c      DO 52 J=1,N
c      LOW=J
c      IF (LEVEL(J) .EQ. LEVLOW) GOTO 53
c   52 CONTINUE
c      CALL REMARK ('LOWER LTESUM LEVEL NOT FOUND')
c      GOTO 990
c   53 CONTINUE
c      ALTESUM(1,LOW)=ASUM
c      ALTESUM(2,LOW)=COEFF1
c      ALTESUM(3,LOW)=COEFF2
c      ENCODE (8,54,ALTESUM(4,LOW)) IRANGE
c      WRITE  (8,54,ALTESUM(4,LOW)) IRANGE
c   54 FORMAT (A8)
c      GOTO 1

C***  AUTOIONIZATION AND DIELECTRONIC RECOMBINATION  -------------------
C***  DRTRANSIT line encountered
   60 NAUTO=NAUTO+1
      IF (NAUTO .GT. MAXAUTO) THEN
         WRITE (0,*) '*** NAUTO .GT. MAXAUTO'
         GOTO 990
      ENDIF
C***  Iron line indices are arranged at the end of the range.
C***  Therefore, no DRTRANSIT lines may occur after the GENERIC element
      IF (LASTFE .GT. 0) THEN
         WRITE (0,*) '*** DRTRANSIT line found behind GENERIC element'
         GOTO 990
      ENDIF
      READ (KARTE,61,ERR=986) 
     >     LEVLOW,LEVUP,NW,EAUTO(NAUTO),AAUTO(NAUTO),LEVION,DRRUDI
   61 FORMAT(10X,A10,2X,A10,1X,I4,1X,F10.0,1X,G10.0,2X,A10,1X,A3)
      LEVUPAUTO(NAUTO) = LEVUP
C***  wrstart may perform some consistency checks
      IF (ROUTINE .EQ. 'WRSTART') THEN
         IF (NW .LE. 0) GOTO 987
         IF (KARTE(33:33) .NE. ' ') GOTO 988
         IF (KARTE(38:38) .NE. ' ') GOTO 988
      ENDIF
      WAUTO(NAUTO)=FLOAT(NW)
      LEVSEQ=1
C***  FIND LOWER INDEX
      DO 64 J=1,N
         LOW=J
         IF (LEVEL(J) .EQ. LEVLOW ) GOTO 65
   64 CONTINUE
      CALL REMARK ('LOWER LEVEL FOR DIEL. RECOMBINATION NOT FOUND')
      STOP 'DRLOWER'
   65 CONTINUE
C***  FIND INDEX OF PARENT ION
      IF (LEVION .EQ. '          ') GOTO 69
      DO 68 J=1, N
      NUP=J
      IF ((LEVEL(J) .EQ. LEVION) .AND. (NCHARG(LOW)+1 .EQ. NCHARG(J)))
     $    GOTO 63
   68 CONTINUE
      PRINT *, 'ERROR: PARENT ION FOR DR ', LEVION, 
     $         'NOT FOUND'
      CALL REMARK ('PARENT ION FOR DR NOT FOUND')
      GOTO 990
   63 IF (NATOM .GT. 1) THEN
         IF (NOM(NUP) .NE. NOM(LOW)) THEN
            CALL REMARK 
     $        ('DIELECTRONIC RECOMBINATION BETWEEN DIFFERENT ELEMENTS')
            STOP 'DRATOMS'
         ENDIF
      ENDIF
   69 CONTINUE
      IF (LEVION .EQ. '          ') THEN
         IONAUTO(NAUTO)=0
      ELSE
         IONAUTO(NAUTO)=NUP
      ENDIF
      LOWAUTO(NAUTO)=LOW

C***  RUDIMENTAL DR-TRANSITIONS 
C***  in the current DRTANSIT data there are no entries in that column 
C***  in the future, this could be used to override the global 
C***  setting via the DRLEVELS option in CARDS
      IF (DRRUDI .NE. '   ') KRUDAUT(NAUTO) = 1

      GOTO 1

C***  END OF INPUT DATA REACHED  ---------------------------------------
    3 CLOSE (4)
 
C***  OLD DATOM FILE (CONTAINING ONLY A HELIUM MODEL ATOM) RECOGNIZED
      IF (NATOM .EQ. 0) THEN
         NATOM=1
         ELEMENT(1)='HELIUM    '
         SYMBOL(1)='HE'
         ATMASS(1)=4.
         STAGE(1)=3.
         KODAT(1)=1
         DO I=1,N
            NOM(I)=1
         ENDDO
      ENDIF

C***  SOME CHECKS OF THE INPUT DATA (only in WRSTART, else skipped)) 
      If (ROUTINE .NE. 'WRSTART') GOTO 78

      IF (N .EQ. 0) THEN
            CALL REMARK ('NO ENERGY LEVELS RECOGNIZED')
            STOP '*** ERROR detected in subr. DATOM'
      ENDIF

C***  ALL ELEMENTS ARE CHECKED ONE BY ONE FOR EXISTENCE OF ANY LEVEL
      DO 84 NA=1,NATOM
         DO 86 I=1,N
            IF (NOM(I) .EQ. NA) GOTO 84
   86    CONTINUE
         CALL REMARK ('ERROR: ELEMENT WITHOUT ANY LEVEL DECODED')
         STOP '*** ERROR detected in subr. DATOM'
   84 CONTINUE

C***  LEVELS ARE CHECKED FOR CORRECT ELEMENT MEMBERSHIP
      DO 85 J=1,N
         IF (LEVEL(J)(:2) .NE. SYMBOL(NOM(J))) THEN
            IF (LEVEL(J)(:1) .NE. SYMBOL(NOM(J))(:1) 
     >         .OR. SYMBOL(NOM(J))(2:2) .NE. ' ') THEN 
               CALL REMARK ('WRONG ELEMENT MEMBERSHIP OF LEVELS')
               GOTO 990
            ENDIF
         ENDIF
   85 CONTINUE

C***  TRANSITIONS ARE CHECKED FOR COMPLETENESS
      DO 7 I=1,N
         DO 7 J=1,N
            IF (NOM(I) .NE. NOM(J)) GOTO 7
            IF (ELEMENT(NOM(I)) .EQ. 'GENERIC') GOTO 7
            IF (NCHARG(I) .NE. NCHARG(J)) GOTO 8
            IF (I.LE.J) GOTO 7
            IF (EINST(I,J) .LE. -99.  ) THEN
            CALL REMARK ('LINE TRANSITION MISSING OR f-VALUE IS GREATER
     >-OR-EQUAL TO 99.0')
            WRITE (0,*) 'LEVEL-NO.: ',I,J
            WRITE (0,*) 'LEVEL-NAMES: "',LEVEL(I),'"   "',LEVEL(J),'"'
     >                  ,' f-VALUE: ', EINST(I,J)
            STOP '*** ERROR detected in subr. DATOM'
            ENDIF
            GOTO 7

C***     Continuum transition
    8    IF (I .GE. J) GOTO 7
C***     CHARGES MUST DIFFER BY 1
         IF (NCHARG(I)+1 .NE. NCHARG(J)) GOTO 7
C***     THERE MUST BE AT LEAST ONE CONTINUUM TRANSITION FOR EACH LOWER LEVEL
         DO 77 KON=1, LASTKON
            IF (KONTLOW(KON) .EQ. I) GOTO 7
   77    CONTINUE
         IF (EINST(I,J) .LT. .0 ) THEN
            CALL REMARK ('CONTINUUM TRANSITION MISSING')
            WRITE (0,'(4A)') 'LOWER LEVEL: ', LEVEL(I), 
     >                     '  UPPER LEVEL: ', LEVEL(J)
            STOP '*** ERROR detected in subr. DATOM'
            ENDIF
    7 CONTINUE
C***  Checks for completeness finished

C***  Consistency check for DRTRANSIT lines:
      IF (NAUTO .GT. 0) THEN
        DO I=1, NAUTO
         DO J=1, NAUTO
            IF (LEVUPAUTO(I) .EQ. LEVUPAUTO(J)) THEN
               IF (WAUTO(I) .NE. WAUTO(J)) THEN
                  WRITE (0,'(A,2I4)') '*** ERROR in DRTRANSIT data: ' //
     >                        'different statistical weights for ' //
     >                        'level ' // LEVUPAUTO(I),  
     >                        INT(WAUTO(I)), INT(WAUTO(J))  
                  GOTO 989
               ENDIF
               IF (EAUTO(I) .NE. EAUTO(J)) THEN
                  WRITE (0,'(A,2X,2F10.2)') 
     >                '*** ERROR in DRTRANSIT data: different energies ' //
     >                'for level ' // LEVUPAUTO(I),  
     >                EAUTO(I), EAUTO(J)  
                  GOTO 989
               ENDIF
            ENDIF
         ENDDO
        ENDDO
      ENDIF

   78 CONTINUE
C***  End of data checks only performed in WRSTART *********************

C***  GENERATE VECTORS NFIRST, NLAST: FIRST AND LAST LEVEL OF EACH ELEMENT
      DO 90 NA=1,NATOM
      IF (NA .EQ. 1) THEN
          NFIRST(NA)=1
      ELSE
          NFIRST(NA)=NLAST(NA-1)+1
      ENDIF
      IF (NA .LT. NATOM) THEN
          NLAST(NA)= ISRCHEQ(N,NOM(1),1,NA+1) - 1
      ELSE
          NLAST(NA)=N
      ENDIF
   90 CONTINUE

C***  GENERATION OF VECTOR IONGRND: DEFAULT LEVEL FOR IONIZATION (LOWEST
C***  LEVEL OF PARENT ION)
      DO 92 J=1, N
      IONGRND(J)=0
      EMIN=999999.
      DO 92 I=1, N
      IF ((NOM(I) .EQ. NOM(J)) .AND. (NCHARG(I) .EQ. NCHARG(J)+1) .AND.
     $    (ELEVEL(I) .LT. EMIN)) THEN
         EMIN=ELEVEL(I)
         IONGRND(J)=I
      ENDIF
   92 CONTINUE

C***  Convert f-values into EINSTEIN coeficients A_up,low
C***  NEGATIVE LINE-CARD ENTRIES INDICATE OSCILLATOR STRENGTHS
      DO 66 IND=1,LASTIND
         NUP=INDNUP(IND)
         LOW=INDLOW(IND)
         AUPLOW=EINST(NUP,LOW)
         IF (AUPLOW .GE. 0.0) GOTO 66
         WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
         EINST(NUP,LOW)=-0.6669*WAVENUM*WAVENUM*AUPLOW*WEIGHT(LOW)/
     /               WEIGHT(NUP)
   66 CONTINUE

C****************
C***  DRTRANSITs 
C****************

C***  CHECK for max. number of line transitions incl. DRTRANSITS
      IF (LASTIND+NAUTO+LASTFE > 99999) THEN
        WRITE (0,*)
     >     '*** MORE THAN 99999 LINE TRANSITIONS ENCOUNTERED ***'
        WRITE (0,*) 'This is not compatible with the encoding of the'
        WRITE (0,'(A)') ' line index in the MODEL file variables'
     >     // ' XJLnnnnn.'
        STOP '*** FATAL ERROR IN DATOM'
      ENDIF

C***  ASSIGNMENT OF DEFAULT IONIZATION LEVEL (GROUND STATE OF PARENT ION)
C***  FOR DIELECTRONIC RECOMBINATION TRANSITIONS
C***  NOTE: ASSUMPTION IS THAT ALL DOUBLY EXCITED STATES AUTOIONIZE
C***        INTO THE GROUND STATE OF THE PARENT ION
      DO 97 I=1, NAUTO
         IF (IONAUTO(I) .EQ. 0) IONAUTO(I)=IONGRND(LOWAUTO(I))
         IF (IONAUTO(I) .NE. IONGRND(LOWAUTO(I))) STOP 'IONAUTO'
C***     Check that the stabilizing transitions have positive wavelength
         LOW=LOWAUTO(I)
C***     WAVENUMBER OF STABILIZING TRANSITION
         WSTABIL = EION(LOW) - ELEVEL(LOW) + EAUTO(I)
         IF (WSTABIL .LE. .0) THEN
           WRITE (0,*)  '*** INCONSISTENCY IN DRTRANSIT DATA (DATOM):'
           WRITE (0,*)  '*** STABILIZING LINE HAS NEGATIVE WAVELENGTH'
           WRITE (0,*)  '*** Transition ', LEVEL(LOWAUTO(I)), ' - ', 
     >               LEVEL(IONAUTO(I))
           STOP '*** FATAL ERROR DETECTED BY SUBR. DATOM'
         ENDIF
   97 CONTINUE

C***  Conversion of oscillator strength f (indicated by neg-sign)
C***   into Aup-low
C***   Bug removed (wrong: EAUTO(LOW) ) wrh 10-Apr-2003 18:47:27
      DO 67 J=1,NAUTO
      AAUTOJ=AAUTO(J)
      IF (AAUTOJ .LT. 0.0) THEN
         LOW=LOWAUTO(J)
         WAVENUM=EION(LOW)-ELEVEL(LOW)+EAUTO(J)
         AAUTO(J)=-0.6669*WAVENUM*WAVENUM*AAUTOJ*WEIGHT(LOW)/WAUTO(J)
      ENDIF
   67 CONTINUE

C***  The DRTRANSITs are arranged behind the regular line transitions in 
C***      the vectors INDLOW, INDNUP
C***      Here, INDNUP is the index of the next-higher ground level
C***      The upper level names are stored in LEVUPAUTO(NAUTO) 
      DO IND=1, NAUTO
        INDLOW(LASTIND+IND) = LOWAUTO(IND)
        INDNUP(LASTIND+IND) = IONAUTO(IND)
      ENDDO

C***  Append the auto-ionizing levels to those vectors that specify
C***    the energy levels (index range N+1 ... N_WITH_DRLEVELS)
      N_WITH_DRLEVELS = N 
      IF (NAUTO .GT. 0)
     >    CALL APPEND_AUTOLEVELS (N, N_WITH_DRLEVELS, NDIM, MAXIND,
     >            MAXAUTO, NAUTO, LOWAUTO, IONAUTO, EAUTO, ELEVEL,
     $            LEVEL, EION, WEIGHT, INDLOW, INDNUP, LASTIND,
     >            LEVUPAUTO, LEVAUTO, WAUTO, NCHARG, IONGRND, NOM)

      RETURN

C***  ERROR exits **************************

  980 WRITE (0,*) 'Ionization Stage must be an integer number'
      GOTO 990

  981 WRITE (0,*) 'Ionization Stage outside valid range'
      GOTO 990

  982 WRITE (0,*) 'K-SHELL-DATA DO NOT FIT TO CURRENT ELEMENT'
      GOTO 990

  983 WRITE (0,*) 'K-SHELL DATA COULD NOT BE DECODED AS NUMBERS'
      GOTO 990

  985 WRITE (0,*) 'ERROR WHEN DECODING LEVEL CARD'
      GOTO 990

  986 WRITE (0,*) 'ERROR WHEN DECODING DRTRANSIT CARD'
      GOTO 990

  987 WRITE (0,*) 'stat. weight < 0 read from DRTRANSIT card'
      GOTO 990

  988 WRITE (0,*) 'Non-blank entry falls into column gap'
      GOTO 990

C***  ERROR BRANCH ********************************************
  989 WRITE (0,'(A,2I4)') '*** You must provide DRTRANSIT data ' //
     >                    'ion the corrected version (after June 2023)'
      GOTO 990


  990 WRITE (0,*) 'The Error was detected when decoding DATOM line:'
      WRITE (0,*) KARTE
      STOP 'ERROR detected by Subr. DATOM'

      END
      SUBROUTINE FEDAT (ROUTINE, INDEXMAX, NFEREADMAX, IONLOW, IONTOP,
     &                  MAXATOM, NDIM, MAXIND, MAXKONT, NATOM,      
     &                  N, LASTFE, LASTKON, LASTINDAUTO, MAXFEIND, 
     &                  EINST, SIGMAFE, INDRB, INDRF, IFENUP, 
     &                  IFELOW, INDNUP, INDLOW, KONTNUP, KONTLOW,
     &                  LEVEL, ELEMENT, SYMBOL, ATMASS, STAGE,
     &                  ELEVEL, WEIGHT, EION, NCHARG, NOM, KODAT,
     &                  NFIRST, NLAST, IFRBSTA, IFRBEND, FEDUMMY,
     &                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, KEYCBB)

c!!!!!! Folgende parameter wurden entfernt: 
C!!!    CBFC, BOUND, EINSTINT, COCOFE, NCOMAX, NCO
c!!!    folgende Parameter sind neu: MAXFEIND, FEDUMMY
C!!!    umbenannt wurden: NMAX -> NFEREADMAX, DUMMY -> FEDUMMY

C **********************************************************************
C ***
C *** CALLED BY: SUBROUTINE DATOM
C ***    READS ALL RELEVANT ATOMIC DATA FOR IRON GROUP LINE BLANKETING
C ***    FROM A MASS-STORAGE FILE CREATED BY THE IRON-PACKAGE (TAPE 21)
C ***
C **********************************************************************

      IMPLICIT NONE

C***  Local dimensions:
C***  Maximum number of bound-bound transitions within one ion
c      PARAMETER ( NBBMAX  = 400 )
      INTEGER, PARAMETER :: NBBMAX = 999  !taken from CREATMS (from blanket program)
C***  NIONMAX must cover the number of Iron Ionization stages in the FEDAT file
      INTEGER, PARAMETER :: NIONMAX = 27
C***  NROMMAX: Roman Numbers encoded: cf. DATA-Statement for ROMNUM 
      INTEGER, PARAMETER :: NROMMAX = 27
C***  Local arrays
      INTEGER, DIMENSION(NIONMAX) :: NB, NTRA, NTRB
      INTEGER, DIMENSION(NBBMAX) :: NX_A, NX_B
      CHARACTER(2) :: CIONBUFFER
      CHARACTER(5), DIMENSION(NROMMAX) :: ROMNUM 
      CHARACTER(16), DIMENSION(NBBMAX) :: NAMARRAY
      CHARACTER(16) :: NAMBUFFER
      CHARACTER(LEN=8) :: NAME
      CHARACTER(LEN=24) :: GENER
      CHARACTER(LEN=40) :: IONNAME(NIONMAX)
      CHARACTER(LEN=3) :: ISTR
      
C***  Maximum number of superlevels within the same ion
      INTEGER, PARAMETER :: MAXLEVEL = 100
      CHARACTER(LEN=8), DIMENSION(MAXLEVEL) :: LEVNAMES

      LOGICAL :: BEXTEND

      INTEGER, INTENT(IN) :: IONLOW, IONTOP, NATOM, NDIM,
     >                       MAXIND, MAXATOM, MAXKONT, MAXFEIND, 
     >                       INDEXMAX, NFEREADMAX, LASTINDAUTO
      INTEGER, INTENT(INOUT) :: N, LASTFE, LASTKON
      
C***  Formal Parameters
      CHARACTER(LEN=*)  ROUTINE
      CHARACTER(LEN=2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(LEN=10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(LEN=10), DIMENSION(NDIM) :: LEVEL
      
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      INTEGER, DIMENSION(NDIM) :: NCHARG, NOM
      REAL, DIMENSION(NDIM) :: EION, ELEVEL, WEIGHT
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW

C***  The following arrays are ONLY used in STEAL, WRSTART 
C***  -> not necessary to define these arrays and MAXIND in other calls 
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      CHARACTER(LEN=4), DIMENSION(MAXIND) :: KEYCBB

C***  IRON-SPECIFIC ARRAYS; MAXFEIND = MAX. NUMBER OF IRON SUPERLINES
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFENUP, IFELOW,
     >                                IFRBSTA, IFRBEND
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      
      REAL :: SIGMAUL, SIGMALU, SIGMAINTCUR,
     >        WLOW, XF, XFOLD, DSIGMA, XLAM, XLOGSTEP, DXFE, 
     >        XLAM0FE, VDOPFE, XBAK, WNUP, XLAMCMIND
      INTEGER :: INDEX, INDEXS, LOWION, NUPION, I, J, K, NION,
     >           NDATA, LEVCOUNT, ILAM, NUP, LOW, INDSTA,
     >           INDEND, INDSELF, IND, IFREQRBSTA, IFREQRBEND,
     >           NOLD, KONT, IADR, MAXADR, IERR, NTEST, IBAK,
     >           NREADIN, IERRLEVNAM, KZEROS, KZEROE
      
      LOGICAL :: bFEINFO, bFEULSEP, bINTRA

C***  Constants:
      REAL, PARAMETER :: CLIGHT = 2.99792458E10   !C in cm/s
      REAL, PARAMETER :: PI8 = 25.1327412288      !PI8 = 8*PI
      REAL, PARAMETER :: FSIG = 2.6540E-2         !PI*e^2/m/c in CGS units
      
C***  Roman Numbers for Level names
      DATA ROMNUM / 'I....', 'II...', 'III..', 'IV...', 'V....', 
     >              'VI...', 'VII..', 'VIII.', 'IX...', 'X....',
     >              'XI...', 'XII..', 'XIII.', 'XIV..', 'XV...',
     >              'XVI..', 'XVII.', 'XVIII', 'XIX..', 'XX...',
     >              'XXI..', 'XXII.', 'XXIII', 'XXIV.', 'XXV..',
     >              'XXVI.', 'XXVII' /

      CALL OPENMS(21,IADR,MAXADR,1,IERR)

C *** READ GENERAL PARAMETERS OF GENERIC ION
      CALL READMS (21, GENER, 3, 'GENERIC ', IERR)
      READ (GENER(15:16), '(F2.0)') STAGE(NATOM)
      READ (GENER(17:24), '(F8.5)') ATMASS(NATOM)
C***  NION = Number of Iron Ionization stages in the FEDAT Data file 
      CALL READMS (21, NION,    1,      'NION    ', IERR)
      IF (NION .GT. NIONMAX .OR. IONTOP .GT. NIONMAX) THEN
         WRITE (0,'(A)') '*** Local dimension NIONMAX insufficient !' 
         WRITE (0,*) NION, IONTOP, NIONMAX
         STOP            '*** ERROR STOP IN FEDAT *********'
      ENDIF

C *** READ PROPERTIES AND NO. OF TRANSITIONS FOR ALL IONIZATION STAGES
      CALL READMS (21, IONNAME, NION*5, 'IONNAME ', IERR)

      CALL READMS (21, NTRA,    NION,   'NTRA_A  ', IERR)
      CALL READMS (21, NTRB,    NION,   'NTRA_B  ', IERR)

C *** READ PARAMETERS OF FREQUENCY GRID
      CALL READMS (21, VDOPFE,  1,      'VDOPP   ', IERR)
      CALL READMS (21, DXFE,    1,      'FSTEP   ', IERR)
      CALL READMS (21, XLAM0FE, 1,      'XLAMNULL', IERR)

C *** INITIALIZE COUNTERS
      INDEX    = 0
      LEVCOUNT = N
      IND      = 0
      KONT     = LASTKON


C *** READ NUMBER OF SUPERLEVELS PER IONIZATION STAGE
        CALL READMS (21, NB, NION, 'NLEV    ', IERR)

C **********************************************************************
C *** LOOP OVER IONIZATION STAGES
C **********************************************************************
      NOLD = N
      DO 10 I=IONLOW, IONTOP

C***    IONLOW == 0 is a placeholder for using the full levels for neutral Fe
        IF (I == 0) CYCLE 

        IF (IONTOP .GT. NION+1) THEN
           WRITE (0, '(A,I3)') 
     >         '*** IONTOP stage requested from DATOM file:', IONTOP 
           WRITE (0, '(A,I3)') 
     >         '*** Highest stage +1 available in FEDAT:', NION+1  
           STOP '*** ERROR DETECTED BY FEDAT' 
        ENDIF

C***    Note: A further ionization stage NION+1, which is not in the data, 
C***        is added if requested. This EXTEND stage has only one level
        BEXTEND = I .EQ. NION+1

C ***   REDUCTION TO 1 LEVEL FOR HIGHEST AND LOWEST IONISATIION STAGE
        IF ((I .EQ. IONLOW).OR.(I .EQ. IONTOP) .OR. BEXTEND) NB(I) = 1

        IF (N+NB(I) .GT. NDIM) THEN
           WRITE (0,'(A)') '*** Dimension NDIM insufficient !' 
           WRITE (0,'(A,I4)') '*** Present value NDIM = ', NDIM
           WRITE (0,'(A,I4)') '*** Required value = ', N+NB(I)
           WRITE (0,*)  I, N, NB(I)
           STOP            '*** ERROR STOP IN FEDAT *********'
        ENDIF

        IF (NB(I) .GT. MAXLEVEL) THEN
           WRITE (0,'(A)') '*** Dimension MAXLEVEL insufficient !' 
           WRITE (0,'(A)') '*** max. number of superlevels in one ion' 
           WRITE (0,'(A,I4)') '*** Present value = ', MAXLEVEL
           STOP            '*** ERROR STOP IN FEDAT *********'
        ENDIF

        IF (.NOT.BEXTEND) THEN

C***      Check dimension: Max number of bound-bound transitions within present ion
           IF (NTRA(I) .GT. NBBMAX .OR. NTRB(I) .GT. NBBMAX) THEN
              WRITE (0,'(A)') '*** Dimension NBBMAX insufficient !' 
              STOP '*** ERROR STOP IN FEDAT'
           ENDIF

C***       STORE MEAN ENERGIES AND STATISTICAL WEIGHTS        
           NAME = 'ELEV' // IONNAME(I)(2:4) 
           CALL READMS (21, ELEVEL(LEVCOUNT+1), NB(I), NAME, IERR)

           NAME = 'WEIG' // IONNAME(I)(2:4) 
           CALL READMS (21, WEIGHT(LEVCOUNT+1), NB(I), NAME, IERR)

C***       Levelnames with parity (if used) - new since 18-Jan-2016
           NAME = 'LEVN' // IONNAME(I)(2:4)
           CALL READMS (21, LEVNAMES, NB(I), NAME, IERRLEVNAM)

C***       READ NUMBER OF DATA-POINTS PER CROSS-SECTION (PRESENT ION I)
           NAME = 'N' // IONNAME(I)(2:4) // '_A  '
           CALL READMS (21, NX_A, NTRA(I), NAME, IERR)

           NAME = 'N' // IONNAME(I)(2:4) // '_B  '
           CALL READMS (21, NX_B, NTRB(I), NAME, IERR)

        ELSE
           ELEVEL(LEVCOUNT+1) = 0.
           WEIGHT(LEVCOUNT+1) = 1.
        ENDIF

        
C***  CREATE SUPERLEVEL NAMES, READ CHARGES AND IONIZATION ENERGIES      
        DO J=1,NB(I)
          N = LEVCOUNT+J
          IF (.NOT. BEXTEND) THEN
            READ(IONNAME(I)(2:3),'(I2)') NCHARG(N)        
            IF (J .EQ. 1) READ(IONNAME(I)(17:24),'(F8.0)') EION(N)
          ELSE
C***        For the extended Level no IONNAME exists
            READ(IONNAME(I-1)(2:3),'(I2)') NCHARG(N)
            NCHARG(N)=NCHARG(N)+1
            IF (J .EQ. 1) EION(N)=0.
          ENDIF
          NOM(N) = NATOM


          IF (IERRLEVNAM /= -10 .AND. .NOT. BEXTEND) THEN
C***        level names are already prepared in the FEDAT file
            LEVEL(N) = SYMBOL(NATOM) // LEVNAMES(J)
          ELSE
C***        use default level names if no predefined names are available
            IF (NCHARG(N)+1 .GT. NROMMAX) THEN
              WRITE (0,'(A)') '*** Roman Number for Ion. stage not known' 
              STOP            '*** ERROR STOP IN FEDAT *********'
            ENDIF
            LEVEL(N) = SYMBOL(NATOM) // ' ' // ROMNUM(NCHARG(N)+1) // '.'
            WRITE (LEVEL(N)(9:10),'(I2)') J
            IF (J .LE. 9) LEVEL(N)(9:9) = '.'
          ENDIF 
        ENDDO

        IF (N .GT. NOLD) THEN
          NFIRST(NATOM) = NOLD+1
          NLAST(NATOM) = N
        ENDIF
                
        NAME = 'A' // IONNAME(I)(2:4) // 'NAM ' 
        CALL READMS (21, NAMARRAY, 2*NTRA(I), NAME, IERR)
                
C**********************************************************************
C***    STORE RBB TRANSITION-DATA IN ONE-DIMENSIONAL ARRAY >>SIGMAFE<<
C***    LOOP OVER ALL BOUND-BOUND TRANSITIONS
        DO 20 J=1,NTRA(I)

C***      NO BB-TRANSITIONS FOR HIGHEST AND LOWEST IONISATION STAGES
          IF ((I .GE. IONTOP).OR.(I .LE. IONLOW)) GOTO 20

C***      READ LEVEL NUMBERS ASSOCIATED WITH TRANSITION          
          CIONBUFFER = NAMARRAY(J)(3:4)
          READ (UNIT=CIONBUFFER,FMT='(I2)') LOWION 
          CIONBUFFER = NAMARRAY(J)(6:7)
          READ (UNIT=CIONBUFFER,FMT='(I2)') NUPION
          LOW=LOWION+LEVCOUNT
          NUP=NUPION+LEVCOUNT

C***      OMIT TRANSITION IF LOW=NUP
          IF (LOW .EQ. NUP) GOTO 20

C***      BB-TRANSITION INDEX FOR IRON LINES (STARTING FROM 1)
          IND = IND + 1        
          IF (IND. GT. MAXFEIND) THEN
           WRITE (0,'(A)') '*** Dimension MAXFEIND insufficient !' 
           STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF     

C***      CREATE POINTER TO STARTING INDEX OF RBB TRANSITION-DATA         
          INDSTA = INDEX+1
          NDATA = NX_A(J)
          IF (NDATA+2 .GT. NFEREADMAX) THEN
              WRITE (0,'(A)') '*** Dim. NFEREADMAX insufficient for b-b!' 
              WRITE (0,'(A, I10)') 
     >          '*** dimensioned: NFEREADMAX = ', NFEREADMAX
              WRITE (0,'(A, I10)') 
     >          '*** required   : NFEREADMAX = ', NDATA+2
              STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          INDEND = INDSTA+NDATA
          IF (INDEND .GE. INDEXMAX) THEN 
              WRITE (0,'(A)') '*** Dimension INDEXMAX insufficient !' 
              WRITE (0,'(A, I10)') 
     >          '*** dimensioned: INDEXMAX =', INDEXMAX
              WRITE (0,'(A, I10)') 
     >          '*** required   : INDEXMAX =', INDEND
              STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          
C ***   STORE LEVEL-NUMBERS IN INDEX-ARRYS          
          IFENUP(IND)=NUP
          IFELOW(IND)=LOW
          
C ***   READ TRANSITION DATA          
          CALL COUNT(J, ISTR)
          NAME = 'A' // IONNAME(I)(2:4) // ISTR 
          CALL READMS (21, FEDUMMY, NDATA+2, NAME, IERR)

C ***   STORE FREQUENCY INDICES
          IFREQRBSTA = - INT(FEDUMMY(2))
          IFREQRBEND = - INT(FEDUMMY(1))

C ***   STORE CROSS-SECTIONS IN ARRAY >>SIGMAFE<<          
          XLOGSTEP = ALOG10(1. + VDOPFE*1.E5*DXFE/CLIGHT)
          SIGMAINTCUR = 0.
          KZEROS = 0
          KZEROE = 0
          DO K=1,NDATA
C***        Calculation of Lambda and Nu in cgs
             ILAM = IFREQRBSTA + K - 1
             XLAM = XLAM0FE*1.E-8 * 10.**(ILAM*XLOGSTEP)
             XF   = CLIGHT / XLAM
C ***       CROSS-SECTION
               SIGMAFE(INDEX+K) = FEDUMMY(NDATA-K+3) 
C***           Determine zero cross section regions at start and end
               IF (SIGMAFE(INDEX+K) <= 0.) THEN
                 IF (KZEROS >= 0) KZEROS = KZEROS + 1
                 KZEROE = KZEROE + 1
               ELSE 
C***             Non-zero cross section found:
                 IF (K == 1) THEN
C***               The first entry in the array is already non-zero.
C***               Therefore stop all further increasements of KZEROS by
C***               setting the counter to a negative value:
                   KZEROS = -1
                 ELSEIF (KZEROS > 0) THEN
C***               To stop further increasing of KZEROS after the first time
C***               this has occured, multiply the result with -1.
                   KZEROS = -1. * KZEROS
                 ENDIF
C***             Reset KZEROE since the formet part with zero cross-section
C***             was definately not at the end of the cross section array.
                 KZEROE = 0
               ENDIF
 
C ***          INTEGRATION OF EINSTEIN-COEFFICIENT AND NORM FOR SIGMA
               IF (K .GT. 1) THEN
                 DSIGMA = (SIGMAFE(INDEX+K)+SIGMAFE(INDEX+K-1))/2.
                 SIGMAINTCUR = SIGMAINTCUR + DSIGMA * (XFOLD - XF)
ccc the following statement is deactivated in libcr_cl version 16-Feb-1999
C***           NU^2 - Term is now accounted for
ccc             XNU = XLAMCMIND/XLAM
ccc                XNUMID = (XNU+XNUOLD)/2.
ccc                XNUMID2= XNUMID*XNUMID
ccc    IMPORTANT: --- consitent change required in CMFFEOP
ccc               --- new linking of both, COLI *and* STEALCL 
ccc                SIGMAINT(IND) = SIGMAINT(IND)
ccc     >                          + XNUMID2 * DSIGMA * (XFOLD - XF)
             ENDIF

             XFOLD = XF
ccc             XNUOLD= XNU
          ENDDO


C***      Band-Band transition 
          INDRB(IND) = INDSTA
            
          IFRBSTA(IND) = IFREQRBSTA
          IFRBEND(IND) = IFREQRBEND
C***      Remove empty cross sections regions from pointer range
C***      (leaves only at maximum one zero entry at beginning and end)
!             IF (ABS(KZEROS) > 0.) THEN
!               DO K=1, ABS(KZEROS)-1
!                 IF (SIGMAFE(INDEX+K) > 0.) STOP 'FATAL: KREDUCINGS FAILED!'
!               ENDDO
! c              WRITE (0,*) 'STA: ', IFREQRBSTA, KZEROS
!               IFRBSTA(IND) = IFREQRBSTA + ABS(KZEROS) - 1
!               INDRB(IND) = INDSTA + ABS(KZEROS) - 1
!             ENDIF
!             IF (ABS(KZEROE) > 0.) THEN
!               DO K=NDATA, NDATA-ABS(KZEROE)+1, -1
!                 IF (SIGMAFE(INDEX+K) > 0.) STOP 'FATAL: KREDUCINGE FAILED!'
!               ENDDO
! c              WRITE (0,*) 'END: ', IFREQRBEND, KZEROE
!               IFRBEND(IND) = IFREQRBEND - ABS(KZEROE) + 1
!             ENDIF
            
            
          SIGMAINT(IND) = SIGMAINTCUR                  
                    
          XLAMCMIND = 1./(ELEVEL(NUP) - ELEVEL(LOW))
          WLOW = WEIGHT(LOW)
          WNUP = WEIGHT(NUP)

          EINST(NUP,LOW) = SIGMAINTCUR * 
     >                       PI8*WLOW/WNUP/(XLAMCMIND*XLAMCMIND)
            
                    
C***      enhance index for next cross section reading     
          INDEX = INDEND

 20    CONTINUE
 
 21    CONTINUE

C **********************************************************************
C *** STORE RBF TRANSITION-DATA IN ARRAY >>EINST<< 
C *** LOOP OVER ALL BOUND-FREE TRANSITIONS
        NAME = 'B' // IONNAME(I)(2:4) // 'NAM'
        CALL READMS (21, NAMARRAY, 2*NTRB(I), NAME, IERR)

        DO 30 J=1, NTRB(I)

C ***   NO BF-TRANSITION FOR HIGHEST IONISATION STAGE
          IF (I .EQ. IONTOP) GOTO 30

C ***   READ LEVEL NUMBERS ASSOCIATED WITH TRANSITION          
          READ (NAMARRAY(J)(6:7),'(I2)') LOWION
          LOW = LOWION + LEVCOUNT
          NUP = LEVCOUNT + NB(I) + 1

C ***   MODEL ION REDUCED TO ONE LEVEL?
          IF (LOWION .GT. NB(I)) GOTO 30
          
C ***   INCREASE CONTINUUM-INDEX (ADDING UP TO "NORMAL" CONTINUA)
          KONT = KONT + 1
          IF (KONT. GT. MAXKONT) THEN
           WRITE (0,'(A)') '*** Dimension MAXKONT insufficient !' 
           STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF     

C ***   STORE LEVEL-NUMBERS IN INDEX-ARRAYS          
          KONTNUP(KONT) = NUP
          KONTLOW(KONT) = LOW
          
C ***   CREATE POINTER AT STARTING INDICES OF RBF TRANSITION-DATA         
          IF (J.EQ.1) THEN
            INDRF(KONT) = INDEX + 1
          ELSE
            INDRF(KONT) = NX_B(J-1) + 2 + INDRF(KONT-1)
          ENDIF
          
C ***   READ TRANSITION DATA          
          CALL COUNT(J, ISTR)
          NAME = 'B' // IONNAME(I)(2:4) // ISTR 
C ***   THE LAST INDEX IS THE COMPOSED CROSS SECTION AT THE EDGE 
          NREADIN = NX_B(J)+3
          IF (NREADIN .GT. NFEREADMAX) THEN
             WRITE (0,'(A)') '*** Dim. NFEREADMAX insufficient for b-f!' 
             STOP            '*** ERROR STOP IN FEDAT *********'
          ENDIF
          CALL READMS (21, FEDUMMY, NREADIN, NAME, IERR)

C ***     STORE COMPOUND THRESOLD CROSS-SECTIONS IN "EINST" [10**-18 CM**2]
          EINST(LOW, NUP) = FEDUMMY(NREADIN) * 1.E18

C ***     STORE CROSS-SECTIONS IN ARRAY >>SIGMAFE<<          
          DO K=1, NX_B(J)+2
            SIGMAFE(INDEX+K) = FEDUMMY(K)
          ENDDO
          INDEX = INDEX+NX_B(J) + 2

 30     CONTINUE
C***    END-OF LOOP OVER BOUND-FREE TRANSITIONS  **********

C***    ACTUALIZE LEVEL-COUNTER        
        LEVCOUNT = LEVCOUNT + NB(I)

 10   CONTINUE
C *** END OF LOOP OVER IONIZATION STAGES
C **********************************************************************

C ***  REORDER BB-TRANSITIONS TO INCREASING FREQUECY INDEX 'IFRBSTA'
 66   NTEST = 0
      DO K=2, IND
         IF (IFRBSTA(K) .LT. IFRBSTA(K-1)) THEN
            NTEST = 1
             
             IBAK = IFRBSTA(K)
             IFRBSTA(K) = IFRBSTA(K-1)
             IFRBSTA(K-1) = IBAK
             
             IBAK = IFRBEND(K)
             IFRBEND(K) = IFRBEND(K-1)
             IFRBEND(K-1) = IBAK
             
             IBAK = INDRB(K)
             INDRB(K) = INDRB(K-1)
             INDRB(K-1) = IBAK
             
             IBAK = IFENUP(K)
             IFENUP(K) = IFENUP(K-1)
             IFENUP(K-1) = IBAK

             IBAK = IFELOW(K)
             IFELOW(K) = IFELOW(K-1)
             IFELOW(K-1) = IBAK
             
             XBAK = SIGMAINT(K)
             SIGMAINT(K) = SIGMAINT(K-1)
             SIGMAINT(K-1) = XBAK
             
          ENDIF
       ENDDO 

       IF (NTEST .EQ. 1) GOTO 66
C ***  END OF BUBBLESORT

C *** SAVE ARRAY-LENGHTS
      LASTKON  = KONT
      LASTFE   = IND

C***  Dimension check
      IF (LASTINDAUTO+LASTFE .GT. MAXIND) THEN
         WRITE (0,'(A)') '*** Dimension MAXIND insufficient !' 
         WRITE (0,'(A,I5)') '*** Available: MAXIND = ', MAXIND
         WRITE (0,'(A,I5)') 
     >            '*** Required : MAXIND = ', LASTINDAUTO+LASTFE
         STOP            '*** ERROR STOP IN FEDAT *********'
      ENDIF
      
C***  Append the Superline indices to the line-transition vectors 
      DO K=1, LASTFE
         INDNUP(LASTINDAUTO+K) = IFENUP(K)
         INDLOW(LASTINDAUTO+K) = IFELOW(K)
C***     set valid keyword for bound-bound collision rates
         KEYCBB(LASTINDAUTO+K) = 'NULL'
      ENDDO

      CALL CLOSMS(21,IERR)
            
      RETURN
      END
      SUBROUTINE FINDCHARGE (THISNAME, NZ)
C*****************************************************************
C***  This subroutine searches for THISNAME in the list of 
C***  chemical elements and returns NZ = core charge 
C***  (Ordnungszahl im Periodischen System)
C***  If THISNAME is not found, NZ=0 is returned.
C***  Note that GENERIC is in the position of IRON 
C***  Called from: DATOM, DECSTAR
C*****************************************************************

      CHARACTER THISNAME*(*)
      PARAMETER (MAXELEM = 26)
      CHARACTER*10 ELEMNAME(MAXELEM)
      DATA ELEMNAME /'HYDROGEN  ', 'HELIUM    ', 'LITHIUM   ', 
     2               'BERYLLIUM ', 'BORON     ', 'CARBON    ', 
     3               'NITROGEN  ', 'OXYGEN    ', 'FLUORINE  ',
     4               'NEON      ', 'SODIUM    ', 'MAGNESIUM ', 
     5               'ALUMINIUM ', 'SILICON   ', 'PHOSPHORUS',
     6               'SULFUR    ', 'CHLORINE  ', 'ARGON     ', 
     7               'POTASSIUM ', 'CALCIUM   ', 'SCANDIUM  ', 
     8               'TITANIUM  ', 'VANADIUM  ', 'CHROMIUM  ', 
     9               'MANGANESE ', 'GENERIC   '/

      NZ = 0
      DO K = 1, MAXELEM
         IF (THISNAME .EQ. ELEMNAME(K)) THEN
           NZ = K
           EXIT
         ENDIF
      ENDDO

      END
      SUBROUTINE GETHISTENTRY(HISTENTRY,JOBNR,MODHIST,MAXHIST)
C***********************************************************************
C***  RETURNS AN ENTRY FOR A GIVEN JOB NUMER FROM THE MODEL HISTORY 
C       CHARACTER ARRAY
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JOBNR, MAXHIST
      CHARACTER(8*MAXHIST), INTENT(IN) :: MODHIST
      CHARACTER(*), INTENT(OUT) :: HISTENTRY

      INTEGER :: LASTCHAR, LAST, BUFFERINT, ISTART, IEND, ILEN,
     >           NRFOUND, IFROM, ITO, I
      CHARACTER(16) :: BUFFER8, JOBNRSTRING, JOBCANDSTR

      BUFFER8 = MODHIST(1:8)
      READ(UNIT=BUFFER8, FMT='(A8)') BUFFERINT
      LAST = BUFFERINT
      LASTCHAR = LAST * 8

      NRFOUND = -1
      ISTART = 0            !ISTART and IEND are used for reading out job numbers
      IEND = 0
      IFROM = 0             !IFROM and ITO are used to read out the history entry

      WRITE(UNIT=JOBNRSTRING, FMT='(I16)') JOBNR
      JOBNRSTRING = ADJUSTL(JOBNRSTRING)

      findjobnr: DO I=9, LASTCHAR
        SELECTCASE (MODHIST(I:I))
          CASE (" ")
            !Blank => Nothing happens
          CASE ("0", "1":"9", "A")
            !Nothing happens
            IF (NRFOUND >= 0) THEN
              NRFOUND = NRFOUND + 1
            ENDIF
          CASE (".")
            IF (NRFOUND > 0) THEN
              !success: number found and this is really a job number
              IF (IFROM == 0) THEN
                IEND = I - 1
                JOBCANDSTR = ADJUSTL(MODHIST(ISTART:IEND))
                IF (TRIM(JOBCANDSTR) == TRIM(JOBNRSTRING)) THEN
                  !This is the job we are searching for
                  IFROM = ISTART - 1
                ENDIF
                NRFOUND = -1
                ISTART = 0
                IEND = 0
              ELSE
                !This is already the jobnumber after the searched job
                ITO = ISTART - 2
                EXIT findjobnr
              ENDIF
            ELSE
              !Dot found at wrong positon => This is not a job number
              NRFOUND = -1
            ENDIF
          CASE ("/")
            IF (NRFOUND < 0) THEN
              ISTART = I + 1
              NRFOUND = 0
            ELSEIF (NRFOUND > 0) THEN
              !Second Slash found => This is not a job number
              ISTART = 0
              IEND = 0
              NRFOUND = -1
            ENDIF
          CASE DEFAULT
            NRFOUND = -1
        ENDSELECT        
        IF ((I == LASTCHAR) .AND. (IFROM /= 0)) THEN
          !Last job was the one that was searched for, end is
          !  determined by end of history instead of next entry
          ITO = LASTCHAR 
        ENDIF
      ENDDO findjobnr

      IF ((IFROM > 0) .AND. (ITO > 0)) THEN
        HISTENTRY = MODHIST(IFROM:ITO)
      ELSE
        HISTENTRY = 'FAILED'
      ENDIF

      RETURN
      END

       FUNCTION IDX(TEXT)
C***  Returns the number of non-blank characters of string TEXT
       CHARACTER*(*) TEXT
       DO IDX=LEN(TEXT),1,-1
          IF(TEXT(IDX:IDX).NE.' ') RETURN
       END DO
       IDX=0
       RETURN
       END
      SUBROUTINE INHIBIT (POPNUM, N, ND, NCHARG, RNE, 
     $                   NATOM, ABXYZ, NFIRST, NLAST, POPMIN)
C*******************************************************************************
C***  THE POPULATION NUMBERS ARE MODIFIED BY THIS SUBROUTINE,
C***     IN ORDER TO AVOID NUMBERS WHICH ARE TOO SMALL OR NEGATIVE
C***  CALLED FROM: STEAL, EXTRAP, MODIFY
C!!!  NOTE : NORMALLY RENORM IS A LOGICAL, BUT THE CRAY PREPROZESSOR (fpp) 
C!!!         HAS A BUG. NOW IRENORM IS USED WITH 
C!!!         RENORM = FALSE --> IRENORM = 0
C!!!         RENORM = TRUE  --> IRENORM = 1
C*******************************************************************************
 
      DIMENSION POPNUM(ND,N)
      DIMENSION NCHARG(N), RNE(ND)
      DIMENSION ABXYZ(NATOM), NFIRST(NATOM), NLAST(NATOM)
C!!!      LOGICAL RENORM

C**********************************************************************
C***  SET POPMIN : SMALLER POPNUMBERS ARE SET TO THIS VALUE
C***     New on 29-Sep-1999 13:51:23 (wrh)
C***     POPMIN is now set by the following CARDS option: 
C***     POPMIN = ...
C***     The Default value (1.E-100) is set in DECSTE
C**********************************************************************

C**********************************************************************
C***  INHIBIT NEGATIVE (AND/OR SMALL) POP.NUMBERS
C**********************************************************************
      NSMALL=0
      NEGWARN=0

      DO 19 L=1,ND
C!!!      RENORM = .FALSE.
      IRENORM = 0

      DO 16 NA=1,NATOM

      DO 16 J=NFIRST(NA),NLAST(NA)
      POPLJ=POPNUM(L,J)

C***  INHIBIT SMALL POP.NUMBERS 
      IF (POPLJ .LT. POPMIN .AND. POPLJ .GE. 0.0) THEN
C!!!            RENORM = .TRUE.
            IRENORM = 1
            NSMALL = NSMALL + 1
            POPNUM(L,J) = POPMIN
            ENDIF
C***  INHIBIT NEGATIVE POP. NUMBERS
      IF (POPLJ .LT. .0) THEN
C!!!            RENORM = .TRUE.
            IRENORM = 1
            NEGWARN = NEGWARN + 1
            POPNUM(L,J) = POPMIN
            ENDIF
   16 CONTINUE
 
C***  RENORMALIZATION OF THE POPULATION NUMBERS
C!!!      IF (RENORM) THEN
      IF (IRENORM .EQ. 1) THEN
         DO 15 NA=1,NATOM
         NFIRNA=NFIRST(NA)
         NLANA=NLAST(NA)
         SUM=0.0

         DO 18 J=NFIRNA,NLANA
   18    SUM = SUM + POPNUM(L,J)

         SUM = SUM / ABXYZ(NA)

         DO 15 J=NFIRNA,NLANA
         POPNUM(L,J) = POPNUM(L,J) / SUM
   15    CONTINUE

C***     RE-ADJUST THE ELECTRON DENSITY
         RNEL=0.0
         DO 20 J=1,N
         RNEL = RNEL + NCHARG(J) * POPNUM(L,J)
   20    CONTINUE
         RNE(L)=RNEL
      ENDIF
   19 CONTINUE
 

C**********************************************************************
C***  PRINTOUT OF WARNINGS AND ERROR MESSAGES  ************************
C**********************************************************************
      IF (NSMALL .GT. 0) PRINT 27, NSMALL, POPMIN
   27 FORMAT (10X, I4, 
     >  ' WARNINGS: POP. NUMBERS .LT.', 1PE8.1, ' SET TO THIS VALUE')
      IF (NEGWARN .GT. 0) PRINT 17, NEGWARN, POPMIN
   17 FORMAT (10X,I4,' WARNINGS: NEGATIVE POP. NUMBERS ARE SET TO ',
     $        1PE8.1 )
  
      RETURN
      END
      SUBROUTINE INTEPO (ND,N,RNE,NCHARG,POPNUM,ENTOT,NATOM,
     >                   NALIM, IONLIM, ABXYZ, NFIRST, NLAST,
     >                   IFRO, ITO, MODE, NOPOP)
C***********************************************************************
C***  LINEAR INTERPOLATION BETWEEN GIVEN DENSITY POINTS OR 
C***  EXTRAPOLATION FROM A GIVEN DENSITY POINT TO INNER OR OUTER BOUNDARY
C***
C***  THE ELECTRON DENSITY IS UPDATED ACCORDING TO THE NEW POPNUMBERS
C***
C***  note: limitation to ION blocks not yet implemented (only ATOM blocks)
C***********************************************************************
      INTEGER, INTENT(IN) :: ND, N, NATOM, NALIM, IONLIM
 
      DIMENSION RNE(ND),NCHARG(N),POPNUM(ND,N),ENTOT(ND)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)

      CHARACTER*7 MODE
      LOGICAL NOPOP

      IF (IFRO.LT.1.OR.IFRO.GT.ND) STOP 'IFRO'
      IF (ITO .LT.1.OR.ITO .GT.ND) STOP 'ITO'
      IF ((MODE .EQ. 'INTERPO').AND.(IFRO.GT.ITO)) THEN
         IHELP=IFRO
         IFRO=ITO
         ITO=IHELP
      ENDIF

      IF (NOPOP) RETURN
  
      IF (IFRO .NE. ITO) THEN
        ADENFR=LOG(ENTOT(IFRO))
        ADENTO=LOG(ENTOT(ITO))
        ADIFF=ADENFR-ADENTO
        IF (MODE .EQ. 'EXTR IN') ADIFF = -ADIFF
      ENDIF 

C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      RNE(L)=.0
      IF (IFRO .NE. ITO) THEN
        IF (MODE .EQ. 'EXTR IN') THEN
          ADEN=(ADENTO-LOG(ENTOT(L)))/ADIFF
        ELSE
          ADEN=(ADENFR-LOG(ENTOT(L)))/ADIFF
        ENDIF
      ENDIF

C***  LOOP FOR EACH ELEMENT  ------------------------------
      DO 1 NA=1,NATOM
      IF (NALIM > 0 .AND. NA /= NALIM) CYCLE        !ATOM limitation (useful with DM split)
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=.0
      DO 2 J=NFIRNA,NLANA

      IF (MODE .EQ. 'INTERPO') THEN
        IF (L.LE.IFRO.OR.L.GE.ITO) GOTO 4
        if (POPNUM(IFRO,J) .gt. 0.) then
        APOPJ=LOG(POPNUM(IFRO,J))
        else
        apopj = -1000
        endif
        if (POPNUM(ITO,J) .gt. 0.) then
        APOPD=APOPJ-LOG(POPNUM(ITO,J))
        else
        apopd = -1000
        endif
        POPNUM(L,J)=EXP(APOPJ-APOPD*ADEN)
      ENDIF

      IF (MODE .EQ. 'EXTR IN') THEN
        IF (L .LE. IFRO) GOTO 4
        IF (IFRO .EQ. ITO) THEN
          POPNUM(L,J) = POPNUM(IFRO,J)
        ELSE
          APOPJ=LOG(POPNUM(ITO,J))
          APOPD=APOPJ-LOG(POPNUM(IFRO,J))
C***      Extrapolated popnumber might not exceed the element abundance
          POPNUM(L,J)=AMIN1(EXP(APOPJ-APOPD*ADEN), ABXYZ(NA))
        ENDIF
      ENDIF

      IF (MODE .EQ. 'EXTROUT') THEN
        IF (L .GE. IFRO) GOTO 4
        IF (IFRO .EQ. ITO) THEN 
           POPNUM(L,J) = POPNUM(IFRO,J)
        ELSE
           APOPJ=LOG(POPNUM(IFRO,J))
           APOPD=APOPJ-LOG(POPNUM(ITO,J))
C***       Extrapolated popnumber might not exceed the element abundance
           POPNUM(L,J)=AMIN1(EXP(APOPJ-APOPD*ADEN), ABXYZ(NA))
c          ENDIF
        ENDIF
      ENDIF

    4 POPSUM=POPSUM+POPNUM(L,J)

    2   CONTINUE
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      POPSUM=POPSUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
      POPNUM(L,J)=POPNUM(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNUM(L,J)
    3 CONTINUE
 
    1 CONTINUE
C***  ENDLOOPS  --------------------------------------------------------
 
      RETURN
      END
      FUNCTION ISRCHEQ(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHEQ, ISRCHNE - Searches a vector for the first element equal or
C***       not equal to a target

C***  SYNOPSIS
C***       index = ISRCHEQ (n, x, incx, target)

C***       index = ISRCHNE (n, x, incx, target)

C***  IMPLEMENTATION
C***       Cray PVP systems

C***  DESCRIPTION
C***       ISRCHEQ searches a real or integer vector for the first element that
C***       is equal to a real or integer target.

C***       ISRCHNE searches a real or integer vector for the first element that
C***       is not equal to a real or integer target.

C***       These functions have the following arguments:

C***       index  Integer.  (output)
C***              Index of the first element equal or not equal to target.
C***              If target is not found, n+1 is returned.
C***              If n <= 0, 0 is returned.

C***       n      Integer.  (input)
C***              Number of elements to be searched.

C***       x      Real or integer array of dimension  (n-1)*|incx|+1.  (input)
C***              Array x contains the vector to be searched.

C***       incx   Integer.  (input)
C***              Increment between elements of the searched array.

C***       target Real or integer.  (input)
C***              Value for which to search in the array.

C***  The Fortran equivalent code for ISRCHEQ is as follows:

      INTEGER X(*), TARGET

      J=1
      ISRCHEQ=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).EQ.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHEQ=I

      RETURN
      END


      SUBROUTINE JSYMSET (SYMBOL,WORD)
C***********************************************************************
C***  USER -DEFINED UNICOS VERSION OF THE COS LIBRARY ROUTINE JSYMSET
C***  OPENS A FILE NAMED AS THE SPECIFIED COS SYMBOL AND
C***  WRITES THE CHARACTER STRING "WORD" INTO THIS FILE
C***********************************************************************
      CHARACTER SYMBOL*(*), WORD*(*)

      OPEN (66,FILE=SYMBOL, STATUS='UNKNOWN')
      WRITE (66,'(A)') WORD
      CLOSE (66)

      RETURN
      END
C***  MAIN PROGRAM MODIFY  *********************************************
      SUBROUTINE MODIFY
C***********************************************************************
C***  THIS PROGRAM INTERPOLATES THE LEVEL POPULATIONS BETWEEN TWO
C***  GIVEM POINTS
C***********************************************************************
 
      IMPLICIT NONE

C***  DEFINE ARRAY DIMENSIONS ******************************************
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM  =          26 
      INTEGER, PARAMETER :: NDIM     =        2560 
      INTEGER, PARAMETER :: NFDIM    = 2*NDIM + 400 
      INTEGER, PARAMETER :: MAXKONT  =     NFDIM/2 
      INTEGER, PARAMETER :: MAXIND   =       45000 
      INTEGER, PARAMETER :: MAXFEIND =        2500 
      INTEGER, PARAMETER :: NDDIM    =          89 
      INTEGER, PARAMETER :: MAXHIST  =        4000 

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 3200 
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      COMMON / COMAUTO / LOWAUTO, WAUTO, EAUTO, AAUTO, IONAUTO, KRUDAUT
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

      INTEGER, DIMENSION(NDIM) :: MAINQN, NOM, IONGRND, NCHARG
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      INTEGER, DIMENSION(MAXATOM) :: NFIRST, NLAST, KODAT
      REAL, DIMENSION(NDDIM) :: RNE, ENTOT, T
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM
      REAL, DIMENSION(MAXKONT) :: ADDCON1, ADDCON2, ADDCON3,
     >                            ALPHA, SEXPO
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXIND) :: INDLOW, INDNUP
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK
      CHARACTER(8*MAXHIST) :: MODHIST

      INTEGER :: I, N, L, IFRO, ITO, JOBNUM, LASTIND, NATOM, NAUTO, ND,
     >           LASTKON, LASTFE, JOBNUM_SAVE, IDUMMY, IERR, NALIM, 
     >           LAST, LSPOP, IARG, NARG, IONLIM, N_WITH_DRLEVELS
      REAL :: A, B, X, Q, FEDUMM, POPMIN, TMIN, TMPREAL,
     >        VDOPFE, DXFE, XLAM0FE, XL, FROM, FTO, ASECOND

      CHARACTER(255) :: HISTENTRY
      CHARACTER(100) :: MODHEAD
      CHARACTER(80) :: KARTE
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(7) :: MODE
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(10), DIMENSION(15) :: ARGUMENT
      CHARACTER(64) :: BUFFER64
      CHARACTER(72) :: BUFFER72
      LOGICAL NOTEMP, BAUTO, NOPOP
 
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl / hydro

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      COMMON /IRON/ FEDUMM, INDRB, INDRF, SIGMAFE, IFRBSTA, IFRBEND,
     >              IFENUP, IFELOW, SIGMAINT
      LOGICAL :: BFEMODEL

C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

      CHARACTER(10) :: TIM1, TIM2

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODINPUT = 1   !MODIFY_INPUT file
      INTEGER, PARAMETER :: hMODEL = 3      !MODEL file (fort.3)
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      INTEGER, EXTERNAL :: IDX

C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> MODIFY started: Program Version from '
     >      , LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

      CALL       DATOM (NDIM,N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'MODIFY', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

C***  READING OF THE MODEL FILE ----------------------------------------
      CALL OPENMS (hMODEL, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (hMODEL, ND, 1, 'ND      ', IERR)
      IF(ND > NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (hMODEL, MODHEAD, 13, 'MODHEAD ' , IERR)
      CALL READMS (hMODEL, JOBNUM ,  1, 'JOBNUM  ' , IERR)
C***  Save JOBNUM
      JOBNUM_SAVE = JOBNUM
      CALL READMS (hMODEL, POPNUM ,ND*N,'POPNUM  ' , IERR)
      CALL READMS (hMODEL, POPMIN ,  1, 'POPMIN  ' , IERR)
      IF (IERR == -10) THEN
        POPMIN = 1.E-99   !Note: STEAL has a default of 1.E-25
      ENDIF
      CALL READMS (hMODEL, RNE    , ND, 'RNE     ' , IERR)
      CALL READMS (hMODEL, ENTOT  , ND, 'ENTOT   ' , IERR)
      CALL READMS (hMODEL, T      , ND, 'T       ' , IERR)
C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (hMODEL,ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF
C***  NOT EXISTING RECORD 'ABXYZ': DEFAULT IS AN ATOMIC DATA FILE "DATOM"
C***  CONTAINING "HELIUM" AS THE ONLY ELEMENT
      IF (IERR .EQ. -10) THEN
         IF (NATOM .EQ. 1) THEN
            ABXYZ(1)=1.
         ELSE
            CALL REMARK ('NOT EXISTING RECORD ABXYZ')
            STOP 'ERROR'
         ENDIF
      ENDIF
      CALL READMS(hMODEL,LAST,1,'MODHIST ' , IERR)
      IF (LAST+9 >= MAXHIST) THEN
          !not really sufficient check because MODIFY can write more than one line
          ! so the limit might still be too small even if this check does not fail
          CALL REMARK ('MODHIST DIMENSION INSUFFICIENT')
          STOP 'ERROR'
      ENDIF
      CALL READMS (hMODEL,MODHIST,LAST,'MODHIST ', IERR)
C*********************************************************************** 

C***  DECODING INPUT CARDS  ********************************************
      OPEN (hMODINPUT, FILE='MODIFY_INPUT', STATUS='UNKNOWN')
      REWIND(hMODINPUT)
      LSPOP=-1
      MODE='NOTHING'
      NOTEMP = .FALSE.
      BAUTO = .FALSE.
      NOPOP = .FALSE.
      NALIM = 0
      IONLIM = 0

C***  Set minimum temperature for extrapolation
      TMIN = 3000.

C***  IF INPUT-FILE IS FINISHED (END=66) THEN UPDATE THE MODEL
    7 READ(hMODINPUT,'(A)',END=66) KARTE

      DO IARG=1, 15
       ARGUMENT(IARG) = ""
      ENDDO

      CALL SARGC (KARTE,NARG)
      DO IARG=1, NARG
       CALL SARGV (KARTE,IARG,ARGUMENT(IARG))
      ENDDO

      IF ( KARTE(:11) .EQ. 'AUTO MODIFY' ) THEN
C                           ===========
            BAUTO = .TRUE.

      ELSEIF ( KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                              =========
            READ(UNIT=KARTE, FMT='(9X,F10.0)') XL
            LSPOP=IFIX(XL)
            IF (LSPOP.EQ.0) LSPOP=1
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF ( ARGUMENT(1)(:9) .EQ. 'INTERPOLA' ) THEN
C                                    =========
            READ (ARGUMENT(4),'(F10.0)') FROM
            READ (ARGUMENT(7),'(F10.0)') FTO
            IFRO=IFIX(FROM)
            ITO=IFIX(FTO)
            MODE='INTERPO'
            IF ( NARG > 7 .AND. ARGUMENT(8) == 'ONLY' ) THEN
              READ (ARGUMENT(10),'(F10.0)') tmpREAL
              IF ( ARGUMENT(9) == 'ATOM' ) THEN
                NALIM = IFIX(tmpREAL)
              ELSEIF ( ARGUMENT(9) == 'ION' ) THEN
                IONLIM = IFIX(tmpREAL)
              ENDIF
            ENDIF
            JOBNUM = JOBNUM + 1
C***        DO INTERPOLATION
            GOTO 80

      ELSEIF ( ARGUMENT(1)(:5) .EQ. 'INNER' 
C                                    =====
     $    .OR. ARGUMENT(1)(:5) .EQ. 'OUTER') THEN
C                                    =====
            IF ( ARGUMENT(2)(:5) .EQ. 'EXTRA' ) THEN
C                                      =====
              READ (ARGUMENT(5),'(F10.0)') FROM
              READ (ARGUMENT(8),'(F10.0)') ASECOND
              IFRO=IFIX(FROM)
              ITO=IFIX(ASECOND)
              JOBNUM = JOBNUM + 1 
              IF ( ARGUMENT(1)(:5) .EQ. 'INNER' ) THEN
                IF (IFRO .LT. ITO) THEN
                  CALL REMARK('EXTRA IN: SECOND POINT WRONG')
                  STOP 'ERROR'
                ENDIF
                MODE='EXTR IN'
              ELSEIF (ARGUMENT(1)(:5) .EQ. 'OUTER') THEN  
                IF (IFRO .GT. ITO) THEN
                  CALL REMARK('EXTRA OUT: SECOND POINT WRONG')
                  STOP 'ERROR'
                ENDIF
                MODE='EXTROUT'
              ENDIF
C***          DO EXTRAPOLATION
              GOTO 80
         
            ELSE

              CALL REMARK('INNER / OUTER: SECOND ARGUMENT WRONG')
              STOP 'ERROR'              

            ENDIF

      ELSEIF (KARTE(:7) .EQ. 'NO TEMP') THEN
C                             =======
            NOTEMP = .TRUE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:4) .EQ. 'TEMP') THEN
C                             ====
            NOTEMP = .FALSE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:6) .EQ. 'NO POP') THEN
C                             ======
            NOPOP = .TRUE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:3) .EQ. 'POP') THEN 
C                             ===
            NOPOP = .FALSE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSE IF (ARGUMENT(1)(1:6) .EQ. 'POPMIN' ) THEN
C                                     =======
            READ (ARGUMENT(2),'(F10.0)') tmpREAL
            IF (POPMIN <= 1.E-98) THEN
              !use CARDS value only if not stored in MODEL file
              POPMIN = tmpREAL
            ENDIF
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (ARGUMENT(1) == 'TMIN') THEN 
C                             ====
            READ (ARGUMENT(2),'(F10.0)') TMIN
C***        READ NEXT INPUT CARD
            GOTO 7

      ENDIF

      GOTO 7
C***  LOOP OF DECODING-OF-INPUT-CARDS  **********************************

C***  INTER- OR EXTRAPOLATION OF THE POPULATIONS  **********************
C***  ENTRY-POINT FOR INTERPOLATION AND EXTRAPOLATION
   80 CONTINUE
          CALL INTEPO (ND,N,RNE,NCHARG,POPNUM,ENTOT,NATOM,NALIM,IONLIM,
     >                           ABXYZ,NFIRST,NLAST,IFRO,ITO,MODE,NOPOP)

C***  INTERPOLATION OF THE TEMPERATURE STRATIFICATION
      IF (.NOT. NOTEMP) THEN
          A = ALOG10(ENTOT(IFRO))
          B = ALOG10(ENTOT(ITO ))

          IF (MODE .EQ. 'INTERPO') THEN
              DO L=IFRO+1, ITO-1
               X = ALOG10(ENTOT(L))
               Q = (X-A) / (B-A)
               T(L) = (1.-Q) * T(IFRO) + Q * T(ITO)
              ENDDO
          ENDIF

          IF (MODE .EQ. 'EXTR IN') THEN
              IF (IFRO .EQ. ITO) THEN
                DO L=IFRO+1,ND
                 T(L) = T(IFRO)
                ENDDO
              ELSE
                DO L=IFRO+1,ND
                 X = ALOG10(ENTOT(L))
                 Q = (X-B) / (A-B)
                 T(L) = (1.-Q) * T(ITO) + Q * T(IFRO)
                ENDDO
              ENDIF
          ENDIF

          IF (MODE .EQ. 'EXTROUT') THEN
              IF ((IFRO .EQ. ITO)) THEN
                DO L=IFRO-1,1,-1
                 T(L) = T(IFRO)
                ENDDO
              ELSE
                DO L=IFRO-1,1,-1
                 X = ALOG10(ENTOT(L))
                 Q = (X-A) / (B-A)
                 T(L) = (1.-Q) * T(IFRO) + Q * T(ITO)
                ENDDO
              ENDIF
          ENDIF
      ENDIF

C***  TO AVOID NEGATIVE OR SMALL POPNUMBERS: CALL INHIBIT
      CALL INHIBIT (POPNUM, N, ND, NCHARG, RNE,
     >              NATOM, ABXYZ, NFIRST, NLAST, POPMIN)

C***  AVOID TEMPERATURES BELOW TMIN (especially negative T)
      DO L=1, ND
        T(L) = MAX(T(L), TMIN)
      ENDDO

C***  UPDATING THE MODEL HISTORY
C***  JOBNUM IS INCREASED FOR EVERY INTERPOLATION OR EXTRAPOLATION CARD
C***  JOBNUM is now set to JOBNUM_SAVE + 1 (if JOBNUM was increased)
      JOBNUM = MIN(JOBNUM, JOBNUM_SAVE+1)
C***  Now no Jobnumbers greater than 999 are written (now obsolete in i7 versions)
C      IF (JOBNUM <= 999) THEN
C        JOBNUM2 = JOBNUM
C      ELSE
C        JOBNUM2 = 999
C      ENDIF

      IF (MODE == 'INTERPO') THEN
C        ENCODE (55,50,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER64, FMT=50) JOBNUM, IFRO, ITO
   50   FORMAT ('/',I7,'. MODIFY  INTERPOLATION  FROM POINT ', 
     >          I3,' TO POINT ',I3,'    ')
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,64,BUFFER64)
C        LAST=LAST+7 !56
      ENDIF
      IF (MODE == 'EXTR IN') THEN
C        ENCODE (63,51,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER72, FMT=51) JOBNUM, IFRO, ITO
   51   FORMAT ('/',I7,'. MODIFY  INNER EXTRAPOLATION FROM POINT ',I3,
     >  ' SECOND POINT ',I3,'  ') 
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,72,BUFFER72)
C        LAST = LAST + 8 !64
      ENDIF
      IF (MODE == 'EXTROUT') THEN
C        ENCODE (63,52,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER72, FMT=52) JOBNUM, IFRO, ITO
   52   FORMAT ('/',I7,'. MODIFY  OUTER EXTRAPOLATION FROM POINT ',I3,
     >  ' SECOND POINT ',I3,'  ')
C        LAST = LAST + 8 !64
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,72,BUFFER72)
      ENDIF
      IF (NOTEMP) THEN
C        LAST = LAST + 1 !8
C        MODHIST(LAST) = ' NOTEMP'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
      ENDIF
      IF (BAUTO) THEN
C        LAST = LAST + 1 !8
C        MODHIST(LAST) = ' AUTO'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' AUTO   ')
      ENDIF
C      MODHIST(1)=LAST
      CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
C***  LOOK FOR NEXT INPUT-CARD
      GOTO 7
 
C***  ENTRY, IF INPUT FILE IS FINISHED
C***  UPDATING THE MODEL, IF MODE <> 'NOTHING'
   66 CONTINUE
      IF (MODE .EQ. 'NOTHING') THEN
         CALL REMARK ('NO INTERPOLATION OPTION DECODED')
         PRINT *, 'MODIFY: NO INTERPOLATION OPTION DECODED'
         GOTO 99
      ENDIF

C***  CYCLIC CHANGING OF THE ARRAY NAMES FOR THE LAST THREE ITERATIONS
      CALL CHANGE ('POP3    '  ,'HELP    '  , hMODEL)
      CALL CHANGE ('POP2    '  ,'POP3    '  , hMODEL)
      CALL CHANGE ('POP1    '  ,'POP2    '  , hMODEL)
      CALL CHANGE ('POPNUM  '  ,'POP1    '  , hMODEL)
      CALL CHANGE ('HELP    '  ,'POPNUM  '  , hMODEL)
      CALL CHANGE ('TOLD3   '  ,'THELP   '  , hMODEL)
      CALL CHANGE ('TOLD2   '  ,'TOLD3   '  , hMODEL)
      CALL CHANGE ('TOLD1   '  ,'TOLD2   '  , hMODEL)
      CALL CHANGE ('T       '  ,'TOLD1   '  , hMODEL)
      CALL CHANGE ('THELP   '  ,'T       '  , hMODEL)

      CALL WRITMS (hMODEL, POPNUM, ND*N, 'POPNUM   ', -1, IDUMMY, IERR) 
      CALL WRITMS (hMODEL, RNE   ,   ND, 'RNE      ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL, T     ,   ND, 'T        ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL, JOBNUM ,  1,  'JOBNUM   ', -1, IDUMMY, IERR)

      IF (LSPOP.GT.0)
     $   CALL    PRIMO (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $                  IFRO,ITO)
 
C***  CLOSE FILES
   99 CLOSE (hMODINPUT)
      CALL CLOSMS (hMODEL, IERR)
      CALL JSYMSET ('G0','0')

      !write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)

      CALL STAMP (OPSYS, 'MODIFY', TIM1)

      STOP 'O.K.'
      END
      SUBROUTINE OPENMS(ICHANNEL, IADR, MAXADR, IFNAME, IERR)
C************************************************************
C***  ROUTINE BY LARS KOESTERKE      8-Sep-1995 15:49:02
C************************************************************

      INTEGER :: IFNAME
      CHARACTER(8) :: FNAME, CSTAT     !IMPORTANT: only the first 7 characters can be used!

      IF ((IFNAME /= 0) .AND. (IFNAME /= 1)) THEN
        WRITE(UNIT=FNAME, FMT='(A8)') IFNAME
      ELSE
        FNAME = '        '
      ENDIF
      CSTAT = 'AUTO'

      CALL CMSSTORE (ICHANNEL, IADR, MAXADR, CSTAT, FNAME, 
     >              DUMMY, IDUMMY, 'OPEN', IERR)

      RETURN
      END
      SUBROUTINE PRIMO (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $                  IFRO,ITO)
C***********************************************************************
C***********************************************************************

      DIMENSION RNE(ND),POPNUM(ND,N)
      CHARACTER MODHEAD*100,LEVEL(N)*10

      PRINT 1,MODHEAD,JOBNUM,IFRO,ITO
    1 FORMAT (1X,  A,  20X,'JOB NO.',I3,
     $ //,20X,'RELATIVE NON-LTE POPULATION NUMBERS MODIFIED FROM',
     $ ' POINT',I3,' TO POINT',I3,/,20X,79('-'))
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2, (LEVEL(J),J=J1,J2)
    2 FORMAT (//,'  L EL.DENS.',10(2X,A10))
      IF (N.LT.J2) PRINT 11
   11 FORMAT (1X)
      PRINT 11
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      PRINT 9,L,RNE(L),(ALOG10(POPNUM(L,J)),J=J1,J2)
    9 FORMAT (I3,F7.3,2X,10F12.2)
    3 CONTINUE
      IF (J2.EQ.N) RETURN
      J1=J1+10
      GOTO 4
      END
      SUBROUTINE READMS(ICHANNEL, X, NDIM, NAME, IERR)
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, X, NDIM, 
     >              'READ', IERR)

      RETURN
      END
      SUBROUTINE REMARK (STRING)

      CHARACTER*(*) STRING

      WRITE (0,'(A)') STRING( :IDX(STRING))

      RETURN
      END
      SUBROUTINE SARGC(TEXT,N)

C**   Diese Funktion ist aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt die Anzahl der Argumente in einem String.
C**   Das Parsen des Strings uebernimmt die Routine sargp.
C**   Dort sind auch die syntaktischen Regeln beschrieben.

      CHARACTER*(*) TEXT
      INTEGER AS,AE,N,I

      I=0 				! Do not look for any Argument
      CALL SARGP(TEXT,N,I,AS,AE)

      RETURN

      END
      SUBROUTINE SARGP(TEXT,N,I,AS,AE)

C**	Die Subroutine sargp zerlegt einen String.
C**	Ermittelt werden:
C**	n : die Anzahl der Argumente und, 
C**     	falls i in [1..n],
C**	as und ae : Start- und Endindex des i-ten Arguments.
C**	Geparsed wird nach folgenden Regel:
C**	Leerzeichen werden grunsaetzlich nicht beachtet
C**	(Ausnahmen siehe unten).
C**	Argumente werden durch Leerzeichen oder Komma oder '=' oder ':'
C**     getrennt.
C**	Wird ein Argument durch Leerzeichen und Komma getrennt,
C**     so gilt dies als eine Trennung.
C**	Folgen zwei Kommata ohne Argument, also hoechstens durch
C**	Leerzeichen getrennt, gilt dies als Leerargument.
C**	Zeichen zwischen zwei doppelten Anfuehrungszeichen gelten
C**	als ein Argument, auch wenn Leerzeichen oder Kommata
C**	enthalten sind. In diesem Fall werden die Anfuehrungszeichen
C**	als nicht zum Argument gehoerig betrachtet.
C**	Interpretaion der Rueckgabewerte:
C**	n: immer die Anzahl der Argumente
C**	as=-1 -> Es wurde kein i-tes Argument gefunden (i nicht in [1..n])
C**	as=0  -> Das i-te Argument war leer (z.B. ",,")
C**	sonst : text(as:ae) = i-tes Argument

        CHARACTER*(*) TEXT
        INTEGER N                  ! Out: Anzahl der Argumente
        INTEGER I                  ! In : gesuchtes Argument
        INTEGER AS,AE              ! Out: erste u. letzte Position des 
                                   !         Arguments im String

	INTEGER TL,TI
        INTEGER STATE

	TL=LEN(TEXT)

	state=0
C	state=0 -> kein Argument aktiv
C	state=1 -> Argumentende gefunden, naechstes Komma
C			ist   k e i n   Leerargument
C	state=2 -> normales Argument aktiv
C	state=3 -> Argument in '"' aktiv


        N=0
        AS=-1
	AE=TL

	DO 100 TI=1,TL

C** Go here looking for next start of an argument
        IF (STATE .EQ. 0) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                            N=N+1
                            IF (N .EQ. I) THEN
                                AS=0
                            ENDIF
                            GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
			   IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 1) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			   STATE=0
			   GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
                           IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 2) THEN
           IF (TEXT(TI:TI) .EQ. ' ') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI
			  GOTO 100
	   ENDIF
	   IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			  STATE=0
                          IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
           GOTO 100
        ELSE 
C***  ! IF (STATE .EQ. 3)
           IF (TEXT(TI:TI) .EQ. '"') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
        ENDIF

100	CONTINUE

	RETURN
        END
      SUBROUTINE SARGV(TEXT,I,ARGTEXT)

C**   Diese Funktion ist fast aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt das i-te Argument in einem String.
C**   Die Regeln zur Argumenttrennung sind in sargp beschrieben.
C**   Ist das i-te Argument nicht vorhanden, wird argtext nicht veraendert.

      CHARACTER*(*) TEXT,ARGTEXT
      INTEGER I
      INTEGER N,AS,AE

      CALL SARGP(TEXT,N,I,AS,AE)

      IF (AS .EQ. -1) GOTO 10
      IF (AS .EQ. 0) THEN
             ARGTEXT=' '
      ELSE
             ARGTEXT=TEXT(AS:AE)
      ENDIF

10    RETURN

      END
      SUBROUTINE SECOND

      STOP 'SECOND NOT IMPLEMENTED AT DEC/OSF'

      RETURN
      END
      SUBROUTINE STAMP (OPSYS, PROGNAME, TIM1)
C**************************************************************
C***  Writes a time stamp in channel 0 (cpr-file)
C***    and the elapsed CPU time into standard-out
C**************************************************************

      CHARACTER OPSYS*(*), PROGNAME*(*), TIM1*(*), TIM2*10 
c      CHARACTER OPSYS*(*), PROGNAME*(*), TIM1*10, TIM2*10 
      REAL*4 DTIME, ETIME, TARRAY(2)



      IF (OPSYS .EQ. 'CRAY') THEN
        CALL CLOCK(TIM2)
        WRITE (0,'(A,A8,2X,3A,F8.1)')
     >      'My Wallclock: ', PROGNAME, TIM1, TIM2, 
     >      'CPU-sec.:', SECOND()
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', SECOND(), ' sec'
      ELSEIF (OPSYS .EQ. 'DEC/UNIX') THEN
           CALL DATE_AND_TIME (DATE=TIM1, TIME=TIM2)
ccc        CALL TIME(TIM2)
ccc        ET = ETIME(TARRAY)
ccc        DTIME gives decimals:
        ET = DTIME(TARRAY)
        WRITE (0,'(A,A8,2X,A,1X,A,1X,A,1X,F8.1)')
     >      'My Wallclock: ', PROGNAME, 
     >       TIM1(1:4) // '/' // TIM1(5:6) // '/' // TIM1(7:8),  
     >       TIM2(1:2) // ':' // TIM2(3:4) // ':' // TIM2(5:6),  
     >      'CPU-sec.:', ET
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', ET, ' sec'
      ELSEIF (OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM2)
        ET = ETIME(TARRAY)
        WRITE (0,'(A,A8,2X,2(A8,1X),A,3X,F8.1)')
     >      'My Wallclock: ', PROGNAME, TIM1, TIM2, 
     >      'CPU-sec.:', ET
        WRITE (*, '(1X,A,F8.1,A)') 
     >      PROGNAME(:IDX(PROGNAME))//'> CPU TIME:', ET, ' sec'
      ENDIF
 
      RETURN
      END
      SUBROUTINE STORAGE (ICHANNEL, IADR, MAXADR, 
     >                    CNAME,    !variable name (or file status if opening a file)
     >                    CNAME2,   !used for: 2nd var (CHANGE), format (MSINFO), kind (WRITE), file name (OPEN)
     >                    X, NDIM, ACTION, MODE, IERR)

C******************************************************************
C***  MASS-STORAGE EMULATOR FOR DEC/OSF1 BY LARS KOESTERKE
C***  VERSION 2.0     21-Jul-1997 13:41:16
C***    New Feautures : INFO-D
C***  VERSION 3.0     15-May-1998 14:16:39
C***    IADR is now stored for several MS-Files 
C***    SOPEN and SCLOSE are new Actions
C******************************************************************

C***  IRECL : RECORD LENGTH OF THE FILE OPENENED BY FORTRAN. THE VARIABLE 
C***                     LENGTH IS 4 BYTE
C***   => changed with compiler parameter "assume byterecl" to be compartible with gfortran
C***  IADRL = IRECL / 2, BECAUSE ALL VARIABLES, EVEN THE CHARACTERS, SHOULD 
C***                     HAVE 8 BYTE

C***  Variable-names have 8 Characters. Note that the Character ^ is not valid
C***    because it is used to transport blanks from the MSFILE plotutility
C***    mcplot.com to msinfo.com (INFO-D)

C***  For INTEL Compiler, 
C***    unless the compiler option "-assume byterecl" is set:
ccc      PARAMETER (IRECL = 256)
C***  For gfortran compiler:
      PARAMETER (IRECL = 1024)

      PARAMETER (IADRL = 128)

      DIMENSION IADR(MAXADR), ISCRATCH(IADRL), SCRATCH(IADRL)
      DIMENSION X(NDIM)
      INTEGER(KIND=1), DIMENSION(IADRL*8) :: SCRATCHBYTE   !for msinfo only - works just with 8 byte defaults

C***  MODE IS NOT USED SO FAR

      CHARACTER(1) :: CKIND
      CHARACTER(7) :: CSTATUS, FN
      CHARACTER(8) :: ACTION, MODE, KSTR, KINDSTR, CNAME, CNAME2
      CHARACTER(10) :: FMTSTR, CACTION

      REAL :: RDUMMY, RSCRATCH
      INTEGER :: IKIND, IDEFKIND, NDIMR, LFMT, IFMTLEN, INTSCRATCH

      LOGICAL BEXIST, BVARKN, BDIMEQ, BNEWKINDINFO, bWRINT

      INTEGER, EXTERNAL :: IDX  !function to obtain the (non-blank) length of a string

C      SAVE IRECL2, NIND2, NIND2U, NINDEX, NREC, NVAR, LASTCH
      SAVE

C***  IWARN = 0  :  No Warning
C***          1  :  Warnings when Reading unknown Variables
C***          2  :  Verbose output
      IWARN = 0

C***  Substitute ^ in blanks in CNAME and CNAME2
      IF (ACTION .EQ. 'INFO-D') THEN
        DO I=1,8
          IF (CNAME(I:I)  .EQ. '^') CNAME(I:I)  = ' '
          IF (CNAME2(I:I) .EQ. '^') CNAME2(I:I) = ' '
        ENDDO
C***  Check for the substring '^' in CNAME and CNAME2. This is not longer
C***  allowed
      ELSEIF (ACTION .EQ. 'READ' .OR. 
     >        ACTION .EQ. 'WRITE' .OR. 
     >        ACTION .EQ. 'LENGTH' .OR. 
     >        ACTION .EQ. 'CHANGE') THEN
        DO I=1,8
          IF (CNAME(I:I)  .EQ. '^' .OR. 
     >        CNAME2(I:I) .EQ. '^') THEN
            WRITE (0,*) 'The Substring ^ is not allowed in the Names', 
     >        CNAME, CNAME2
            STOP 'ERROR in Subr. STORAGE'
          ENDIF
        ENDDO
      ENDIF

      IF (MODE(1:4) .EQ. 'CRAY') THEN
        READ(UNIT=CNAME, FMT='(A8)') NAME
        READ(UNIT=CNAME2, FMT='(A8)') NAME2
c    1   FORMAT (A8)
      ELSE
        WRITE (0,*) 'MODE NICHT CRAY'
        WRITE (0,*) 'MODE=',MODE(1:4)
        STOP 'ERROR IN STORAGE'
      ENDIF

      IF (IWARN .GE. 2) THEN
        WRITE(0,'(A,A8,A,I2,2X,A,A)') 
     >        'storage: action=',action,' ICHANNEL=',ICHANNEL, 
     >        'NAME=',NAME
      ENDIF

C***  CHECK THE ICHANNEL NUMBER
      IF (ICHANNEL .LE. 0) THEN
        WRITE (0,*) ' NEGATIVE ICHANNEL NUMBERS ARE NOT ALLOWED'
        STOP 'ERROR IN STORAGE'
      ENDIF

C***  CHECK FOR MINIMUM LENGTH OF THE INDEX ARRAY IADR
      IF (MAXADR .LT. IRECL) THEN
        WRITE (0,*) ' DIMENSION (MAXADR) OF INDEX ARRAY IADR SMALLER',
     >              ' THAN RECORD LENGTH (IRECL) OF MASS-STORAGE FILE'
        STOP 'ERROR IN ROUTINE STORAGE'
      ENDIF

C***  NUMBER OF INDEX RECORDS (NIND) CLAIMED BY MAXADR
      NIND = (MAXADR - 1) / IADRL + 1
c      write (*,*) 'STORAGE : action, nind=', action,nind

C***                        ====
      IF (ACTION( :4) .EQ. 'OPEN') THEN
C***                        ====

C***  FILE-NAME ON DISK: fort.<ICHANNEL>
        IF (IDX(CNAME2) > 0) THEN
          FN = CNAME2(1:7)
        ELSEIF (ICHANNEL < 10) THEN
          WRITE(FN, '("fort.",I1,1X)') ICHANNEL
        ELSE
          WRITE(FN, '("fort.",I2   )') ICHANNEL
        ENDIF

C***  CHECK IF FILE DOES EXIST
        INQUIRE (FILE=FN, EXIST=BEXIST)
        CACTION = 'READWRITE'
        IF (CNAME == 'AUTO') THEN
          IF (BEXIST) THEN
            CSTATUS = 'OLD'
          ELSE
            CSTATUS = 'NEW'
          ENDIF
        ELSEIF (CNAME == 'READ') THEN
          CSTATUS = 'OLD'
          CACTION = 'READ'
        ELSE
          CSTATUS = CNAME(1:7)
          IF ((CSTATUS == 'REPLACE') .OR. (CSTATUS(1:3) == 'NEW')) THEN
            BEXIST = .FALSE.
          ENDIF
        ENDIF

        OPEN (UNIT=ICHANNEL,
     >        FILE=FN,
     >        ACCESS='DIRECT',
     >        FORM='UNFORMATTED',
     >        RECL=IRECL,
     >        STATUS=CSTATUS,
     >        ACTION=CACTION,
     >        ERR=90,
     >        IOSTAT=IOS)

C***  READ IADR IF FILE EXISTS
C***    READ FIRST RECORD
        IF (BEXIST) THEN
          READ (UNIT=ICHANNEL, REC=1, ERR=91, IOSTAT=IOS) ISCRATCH
          DO I=1, IADRL
            IADR(I) = ISCRATCH(I)
          ENDDO

C***  INTERPRET THE FIRST ELEMENTS
C***          IRECL2   : RECORD LENGTH OF THE EXISTING FILE
C***          NIND2    : NUMBER OF INDEX RECORDS OF THE EXISTING FILE
C***          NIND2U   : NUMBER OF INDEX RECORDS USED IN THE EXISTING FILE
C***          NINDEX   : NUMBER OF INDICES IN INDEX RECORD
C***          NREC     : TOTAL NUMBER OF RECORDS
C***          NVAR     : NUMBER OF VARIABLES
C***          THE INDEX ARRAY (IADR) IS USED FOR THE STORAGE OF THE 
C***            INFORMATION ABOUT THE ARRAYS STORED IN THE FILE
C***            FIVE ENTRYS ARE USED FOR EACH ARRAY
C***          IADR(11) : NAME OF THE ARRAY
C***          IADR(12) : NUMBER OF FIRST RECORD
C***          IADR(13) : NUMBER OF RECORDS
C***          IADR(14) : NUMBER OF VARIABLES
C***          IADR(15) : DATA TYPE and KIND (was UNUSED until May 2012)
C***          IADR(16) : etc. (next entry!)
          IRECL2 = IADR(1)
          NIND2  = IADR(2)
          NIND2U = IADR(3)
          NINDEX = IADR(4)
          NREC   = IADR(5)
          NVAR   = IADR(6)

C***  CONSISTENCY CHECKS

C***    COMPARE RECORD LENGTH
C***      the following is no longer fatal and therefore not reported
c          IF (IRECL .NE. IRECL2) THEN
c            WRITE (0,*) ' RECORD LENGTH OF FILE (IRECL2) AND',
c     >                  ' ROUTINE (IRECL) DO NOT MATCH'
c            WRITE (0,'(A7,I4,A8,I4)') ' IRECL=',IRECL, ' IRECL2=',IRECL2
c            WRITE (0,'(A)') '- If compiled with "ifort", make sure that'
c     >                   // ' the option "-assume byterecl" was used!'
cC!            STOP 'ERROR IN STORAGE'
c          ENDIF

C***    COMPARE NUMBER OF INDEX RECORDS
          IF (NIND .EQ. NIND2) THEN
          ELSE IF (NIND .GT. NIND2) THEN
cc  message suppressed! wrh 12-Aug-2008 13:58:44
cc            WRITE (0,*) 'INFO from STORAGE: ' // FN 
cc     >       // ' has less Index-Records than dimensioned > expanded'
          ELSE
            WRITE (0,*) 'Number of Index-Records in File is greater'
            WRITE (0,*) 'than in Index-Array'
C            WRITE (0,*) ' MIN IS TAKEN?'
            WRITE (0,'(A6,I4,A6,I4)') 'NIND2=', NIND2, ' NIND=', NIND
            WRITE (0,'(A,I3)') 'CHANNEL=',ICHANNEL
            STOP 'ERROR IN STORAGE'
          ENDIF

C***    READ THE REST OF THE INDEX ARRAY
          DO I=2, NIND2U
            READ (UNIT=ICHANNEL, REC=I) ISCRATCH
            DO J=1, IADRL
              IADR((I-1)*IADRL + J) = ISCRATCH(J)
            ENDDO
          ENDDO

        ELSE
C***  NEW FILE OPENED
          IRECL2 = IRECL
          NIND2  = (MAXADR - 1) / IADRL + 1
          NIND2U = 1
          NINDEX = 0
          NREC   = NIND2
          NVAR   = 0
        ENDIF

C***                        ====
      ELSE IF (ACTION( :5) .EQ. 'SOPEN') THEN
C***                        ====
          IRECL2 = IADR(1)
          NIND2  = IADR(2)
          NIND2U = IADR(3)
          NINDEX = IADR(4)
          NREC   = IADR(5)
          NVAR   = IADR(6)

C***                             =====
      ELSE IF (ACTION( :5) .EQ. 'WRITE') THEN
C***                             =====
C***  TRANSFORM SIZE IF KIND IS DIFFERENT 
        NDIMR = NDIM
        IF (IDX(CNAME2) > 0) THEN
          IF ( (CNAME2(1:1) == 'I') .OR.
     >         (CNAME2(1:1) == 'i') .OR.
     >         (CNAME2(1:1) == 'R') .OR.
     >         (CNAME2(1:1) == 'r') ) THEN
            !kind calculations only required for data types which can have different kinds (integer, real)
            KINDSTR = ''
            DO J=2, 8
              SELECTCASE(CNAME2(J:J))
                CASE ('1':'9', '0')
                  KINDSTR = TRIM(ADJUSTL(KINDSTR)) // CNAME2(J:J)
                CASE DEFAULT
                  EXIT
              ENDSELECT
            ENDDO
            READ(UNIT=KINDSTR, FMT='(I8)') IKIND       !transform number part from kind string into integer
            !Get default kind for used type
            IF ((CNAME2(1:1) == 'I') .OR. (CNAME2(1:1) == 'i')) THEN
              IDEFKIND = KIND(IDEFKIND)
            ELSE
              IDEFKIND = KIND(RDUMMY)
            ENDIF
            !rescale dimension with used kind
            NDIMR = NDIMR * IKIND / IDEFKIND
          ENDIF
        ENDIF

C***  CHECK IF VARIABLE IS KNOWN
        BVARKN = .FALSE.
        INUM = 0
        BNEWKINDINFO = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
c          WRITE (0,'(A9,A8,1x,a8,I6)') 'VARIABLE=',NAME,IADR(1+INDEX),
c     >    IADR(4+INDEX)
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            INUM   = IADR(4+INDEX)
            IF ((IDX(CNAME2) > 0) .AND. (NAME2 /= IADR(5+INDEX))) THEN
              !new identification string => can only be used if total byte size matches
              BNEWKINDINFO = .TRUE.
              IF (IADR(INDEX+5) == -1) THEN
                !no format set so far
                BDIMEQ = (NDIMR == INUM)                   
              ELSE
                WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)  !transform from integer into character
                KINDSTR = ''
                DO J=2, 8
                  SELECTCASE(KSTR(J:J))
                    CASE ('1':'9', '0')
                      KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                    CASE DEFAULT
                      EXIT
                  ENDSELECT
                ENDDO
                READ(UNIT=KINDSTR, FMT='(I8)') IKIND      !transform number part from kind string into integer
                IDEFKIND = KIND(RDUMMY)
                BDIMEQ = (NDIMR == NDIM * IKIND / IDEFKIND)
              ENDIF
            ELSEIF (INUM .EQ. NDIM) THEN
              BDIMEQ = .TRUE.
            ELSE
              BDIMEQ = .FALSE.
            ENDIF
            IF (BDIMEQ) THEN
              BVARKN = .TRUE.
              IFIRST = IADR(2+INDEX)
C***          EXIT THIS LOOP
              GOTO 40
            ELSE
              WRITE (0,*) ' INCONSISTENCE: VARIABLE FOUND BUT WITH',
     >                    ' DIFFERENT ARRAY LENGTH'
              WRITE (0,'(A9,A8)') 'VARIABLE=',NAME
              WRITE (0,'(A21,I6)') 'DIMENSION OF ARRAY : ', NDIM
              WRITE (0,'(A21,I6)') 'DIMENSION IN FILE  : ', INUM
              STOP 'ERROR IN STORAGE (ACTION : WRITE)'
            ENDIF
          ENDIF
        ENDDO
   40   CONTINUE

C***    NUMBER OF RECORDS USED FOR THE ARRAY WHICH WILL BE STORED
        NIND3 = (NDIMR - 1) / IADRL + 1
        IF (.NOT. BVARKN) THEN
C***    UPDATE INDEX AND APPEND NAME OF THE ARRAY
          NVAR = NVAR + 1
          IFIRST = NREC + 1
          NREC = NREC + NIND3
          NIND2U = (10 + NVAR*5 - 1) / IADRL + 1
C          IF (NIND2U .GT. NIND2) THEN
C            WRITE (0,*) ' INDEX ARRAY IS NOT LARGE ENOUGH TO',
C     >                  ' RECEPT A NEW ARRAY'
C            WRITE (0,*) 'NIND2=', NIND2
C            STOP 'ERROR IN STORAGE'
C          ENDIF
          NINDEX = NINDEX + 1
          INDEX = 10 + (NINDEX-1)*5
          IADR(INDEX+1) = NAME
          IADR(INDEX+2) = IFIRST
          IADR(INDEX+3) = NIND3
          IADR(INDEX+4) = NDIM
          IF (IDX(CNAME2) > 0) THEN
            IADR(INDEX+5) = NAME2
          ELSE
            IADR(INDEX+5) = -1
          ENDIF
        ELSEIF (BNEWKINDINFO) THEN
          IADR(INDEX+4) = NDIM
          IADR(INDEX+5) = NAME2
        ENDIF
        
C!!!  OLD VERSION (MIT UMKOPIEREN DES GESAMTEN ARRAYS)
C***    AS LONG AS THE ARRAY FILLS THE NEXT RECORD COPY IT TO SCRATCH
C***    THE LAST RECORD IS, IF NECESSARY, FILLD UP WITH 0.
          NREST = NDIMR
          DO I=1, NIND3
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)
              SCRATCH(J) = X(J+INDEX)
            ENDDO
            DO J=NREST+1, IADRL
              SCRATCH(J) = 0.
            ENDDO
            WRITE(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=93, IOSTAT=IOS) 
     >            SCRATCH
            NREST = NREST - IADRL
          ENDDO
C!!!          DO I=1, NIND3
C!!!            INDEX = (I-1) * IADRL
C!!!            WRITE(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=93, IOSTAT=IOS) 
C!!!     >            X(INDEX+1)
C!!!          ENDDO

C***                             ====
      ELSE IF (ACTION( :4) .EQ. 'READ') THEN
C***                             ====
C***  CHECK IF VARIABLE IS KNOWN
        NDIMR = NDIM
        BVARKN = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
c      write (0,'(A,2A8)') 'testname=',IADR(1+INDEX), name
c      write (0,'(A,2I)') 'testname=',IADR(1+INDEX), name
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            IF (IADR(5+INDEX) /= -1) THEN
              !Rescale NDIM (or copy NDIMR) if different KIND specified
              WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)  !transform from integer into character
              KINDSTR = ''
              IF ( (KSTR(1:1) == 'I') .OR.
     >             (KSTR(1:1) == 'i') .OR.
     >             (KSTR(1:1) == 'R') .OR.
     >             (KSTR(1:1) == 'r') ) THEN
                !kind calculations only required for data types which can have different kinds (integer, real)
                DO J=2, 8
                  SELECTCASE(KSTR(J:J))
                    CASE ('1':'9', '0')
                      KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                    CASE DEFAULT
                      EXIT
                  ENDSELECT
                ENDDO
                READ(UNIT=KINDSTR, FMT='(I8)') IKIND      !transform number part from kind string into integer
                !Get default kind for used type
                IF ((KSTR(1:1) == 'I') .OR. (KSTR(1:1) == 'i')) THEN
                  IDEFKIND = KIND(IDEFKIND)
                ELSE
                  IDEFKIND = KIND(RDUMMY)
                ENDIF
                !rescale dimension with used kind
                NDIMR = NDIMR * IKIND / IDEFKIND
              ENDIF
            ENDIF
            IF (IADR(4+INDEX) .EQ. NDIMR) THEN
              BVARKN = .TRUE.
            ELSEIF (IADR(4+INDEX) .GT. NDIMR) THEN
              BVARKN = .TRUE.
              IF (IWARN .GE. 2) THEN
                WRITE (0,*) ' WARNING from Subroutine STORAGE'
                WRITE (0,*) ' Inconsistence: Variable found but with',
     >                      ' different Array Length (Action: READ)'
                WRITE (0,'(A9,A8)') 'Variable=',NAME
                WRITE (0,'(A21,I6)') 'Dimension of Array : ', NDIMR
                INUM   = IADR(4+INDEX)
                WRITE (0,'(A21,I6)') 'Dimension in File  : ', INUM
              ENDIF
            ELSE
              BVARKN = .TRUE.
              INUM   = IADR(4+INDEX)
              IF (IWARN .GE. 1) THEN
                WRITE (0,*) ' WARNING from Subroutine STORAGE'
                WRITE (0,*) ' Inconsistence: Variable found but with',
     >                      ' different Array Length (Action: READ)'
                WRITE (0,'(A9,A8)') 'Variable=',NAME
                WRITE (0,'(A21,I6)') 'Dimension of Array : ', NDIMR
                WRITE (0,'(A21,I6)') 'Dimension in File  : ', INUM
C!!!                STOP 'ERROR IN STORAGE (ACTION: READ)'
              ENDIF
              NDIMR = MIN(NDIMR,INUM)
            ENDIF
C***        EXIT THIS LOOP
            GOTO 45
          ENDIF
        ENDDO
   45   CONTINUE
        IF (.NOT. BVARKN) THEN
          IF (IWARN .GE. 1) THEN
            WRITE (0,'(A,A8,A)') 
     >        ' WARNING from STORAGE: Variable ', NAME, 'not found'
          ENDIF
          IERR = -10
          RETURN
C!!!          STOP 'ERROR in STORAGE (ACTION: READ)'
        ELSE
          IERR = 0
          IFIRST = IADR(2+INDEX)
          INUM   = IADR(4+INDEX)
        ENDIF

C!!!  ALTE VERSION (MIT UMKOPIEREN DES GESAMTEN ARRAYS)
C***    AS LONG AS THE ARRAY FILLS THE NEXT RECORD COPY IT TO SCRATCH
C***    THE LAST RECORD IS, IF NECESSARY, FILLED UP WITH 0.
          NIND3 = (NDIMR-1) / IADRL + 1
          NREST = NDIMR
          DO I=1, NIND3
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >            SCRATCH
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)
              X(J+INDEX) = SCRATCH(J)
            ENDDO
            NREST = NREST - IADRL
          ENDDO

C!!!          NIND3 = (NDIM-1) / IADRL + 1
C!!!          DO I=1, NIND3-1
C!!!            INDEX = (I-1) * IADRL
C!!!            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
C!!!     >            X(INDEX+1)
C!!!          ENDDO
C!!!          INDEX = (NIND3-1) * IADRL
C!!!          READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
C!!!     >          SCRATCH
C!!!          NREST = NDIM - ((NIND3-1) * IADRL)
C!!!          DO J=1, NREST
C!!!            X(J+INDEX) = SCRATCH(J)
C!!!          ENDDO


C***  New option, returns in NDIM the array lenth of the variable
C***  wrh 10-Aug-2007 
C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'LENGTH') THEN
C***                             =====
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
          IF (IADR(1+INDEX) .EQ. NAME) THEN 
             NDIM = IADR(4+INDEX) 
             GOTO 14
          ENDIF
        ENDDO
        WRITE (0,*) 'ERROR: COULD NOT FIND LENGTH OF MS-VARIABLE ', 
     >              CNAME
        STOP 'FATAL ERROR in subr. STORAGE (action: LENGTH)'
   14   CONTINUE
C***                             =====
      ELSE IF (ACTION( :5) .EQ. 'CLOSE') THEN
C***                             =====

C***  STORE IADR(1..6)
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR
C***  WRITE INDEX RECORD TO FILE
        DO I=1, NIND2U
          DO J=1, IADRL
            ISCRATCH(J) = IADR((I-1)*IADRL + J)
          ENDDO
          WRITE (UNIT=ICHANNEL, REC=I, ERR=92, IOSTAT=IOS) ISCRATCH
        ENDDO
        CLOSE(ICHANNEL)

C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'SCLOSE') THEN
C***                             =====
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR

C***                             =====
      ELSE IF (ACTION( :6) .EQ. 'CHANGE') THEN
C***                             =====
        BVARKN = .FALSE.
        DO I=1, IADR(6)
          INDEX = 10 + (I-1)*5
c          write (0,'(a,i8,a,a8,a,a8)') 'index=',index, 
c     >                            'nameold=',IADR(1+INDEX),
c     >                            'testname=',cname
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            BVARKN = .TRUE.
            GOTO 10
          ENDIF
        ENDDO
   10   CONTINUE
        IF (BVARKN) THEN
          IADR(INDEX+1) = NAME2
        ELSE
          WRITE (0,*) 'WARNING: CHANGE: Variable not found'
c          WRITE (0,'(A9,A8,1x,a8)') 'VAR1,2',NAME,NAME2
        ENDIF

C***                             =====
      ELSE IF (ACTION( :4) .EQ. 'INFO' .AND. 
     >         ACTION( :6) .NE. 'INFO-D') THEN
C***                             ====
C        write (*,*) 'test---------------'
        IADR(1) = IRECL2
        IADR(2) = NIND2
        IADR(3) = NIND2U
        IADR(4) = NINDEX
        IADR(5) = NREC
        IADR(6) = NVAR
        WRITE(*,'(I9,4X,A)') IRECL2, ': RECORD LENGTH OF THE FILE'
        WRITE(*,'(I9,4X,A)') NIND2,  ': NUMBER OF INDEX RECORDS'
        WRITE(*,'(I9,4X,A)') NIND2U, ': NUMBER OF INDEX RECORDS USED'
        WRITE(*,'(I9,1X,A)') NINDEX, ': NUMBER OF INDICES IN INDEX REC'
        WRITE(*,'(I9,1X,A)') NREC,   ': TOTAL NUMBER OF RECORDS'
        WRITE(*,'(I9,4X,A)') NVAR,   ': NUMBER OF VARIABLES'

        IF (ACTION( :6) .EQ. 'INFO-L') THEN
C***                          ------
          DO I=1, NVAR
            INDEX = 10 + (I-1)*5
            IF (IADR(5+INDEX) == -1) THEN
              KSTR = 'default'
            ELSE
              WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)
            ENDIF
            WRITE (*,50) I, IADR(1+INDEX), 
     >                   IADR(2+INDEX), 
     >                   IADR(3+INDEX), IADR(4+INDEX), KSTR
   50       FORMAT (' * Variable Nr.', I6, 2X, 'Name=', A, 2X, 
     >              'First=', I6, 2X, 'Num-Rec=', I6, 2X, 
     >              'Num-Var=', I8, 2X, 'Vartype=', A)

          ENDDO
        ENDIF
      ELSE IF (ACTION( :6) .EQ. 'INFO-D') THEN
C***                             ------
C***  CHECK IF VARIABLE IS KNOWN
        BVARKN = .FALSE.
        DO I=1, NVAR
          INDEX = 10 + (I-1)*5
          ID = IDX(CNAME)
c        write (0,'(4(a,1x))') 'INFO-D :', IADR(1+INDEX), NAME, ':'
          IF (IADR(1+INDEX) .EQ. NAME) THEN
            BVARKN = .TRUE.
C***        EXIT THIS LOOP
            GOTO 44
          ENDIF
        ENDDO
   44   CONTINUE
        IF (.NOT. BVARKN) THEN
          WRITE (0,'(3A)') 'Variable ',CNAME, ' not known'
          STOP 'ERROR WHEN ACTION = INFO-D'
        ENDIF
        IERR = 0
        IFIRST = IADR(2+INDEX)
        INUM   = IADR(4+INDEX)

C***    Kurzuebersicht ueber die Variable mit der Nummer I
        INDEX = 10 + (I-1)*5
        IF (IADR(5+INDEX) == -1) THEN
          KSTR = 'default'
          CKIND = ' '
        ELSE
          WRITE(UNIT=KSTR, FMT='(A8)') IADR(5+INDEX)
          SELECTCASE(KSTR(1:1))
            CASE ('i','I') 
              CKIND = 'I'
            CASE ('r','R') 
              CKIND = 'R'
            CASE ('a','A')
              CKIND = 'A'
            CASE ('c','C') 
              CKIND = 'C'
            CASE DEFAULT
              CKIND = ' '
          ENDSELECT
          IF (CKIND /= ' ') THEN
            KINDSTR = ''
            DO J=2, 8
              SELECTCASE(KSTR(J:J))
                CASE ('1':'9', '0')
                  KINDSTR = TRIM(ADJUSTL(KINDSTR)) // KSTR(J:J)
                CASE DEFAULT
                  EXIT
              ENDSELECT
            ENDDO
            READ(UNIT=KINDSTR, FMT='(I8)') IKIND       !transform number part from kind string into integer
          ENDIF
        ENDIF
        WRITE (*,50) I, IADR(1+INDEX), 
     >               IADR(2+INDEX), 
     >               IADR(3+INDEX), IADR(4+INDEX), KSTR

C***    Ausgabe der Variable mit der Nummer I
        NIND3 = (INUM-1) / IADRL + 1
        NREST = INUM
        LFMT = IDX(CNAME2)
        bWRINT = .FALSE.
        IF ((CNAME2(1:4) == 'AUTO') .OR. (CNAME2(1:4) == 'auto')) THEN
          IF (CKIND == 'I') THEN
            IFMTLEN = INT(  LOG10(2. ** FLOAT(IKIND * 8)) + 2. )
            WRITE(UNIT=FMTSTR, FMT='(I7,A1)') IFMTLEN, ')'
            FMTSTR = '(I' // TRIM(ADJUSTL(FMTSTR))
          ELSEIF (CKIND == 'R') THEN
            FMTSTR = '(G16.6)'
          ELSEIF ((CKIND == 'A') .OR. (CKIND == 'C')) THEN
            FMTSTR = '(A)'
          ELSE
            WRITE (*,*) '* WARNING: No vartype specified,',
     >                    ' auto-format is not recommended'
            SELECTCASE (CNAME(1:1)) !guessing based on first letter (implicit fortran)
              CASE ('I':'N','i':'n')
                FMTSTR = '(I16)'
              CASE DEFAULT
                FMTSTR = '(G16.6)'
            ENDSELECT
          ENDIF  
        ELSEIF (       (CNAME2(1:1) /= '(') 
     >           .AND. (CNAME2(LFTM:LFTM) /= ')')       ) THEN
          FMTSTR = '(' // TRIM(ADJUSTL(CNAME2)) // ')'
        ELSE
          FMTSTR = CNAME2
        ENDIF
        FMTSTR = ADJUSTL(FMTSTR)
        IF (FMTSTR(1:2) == '(I') bWRINT = .TRUE.
        WRITE (*,*) 'N=?'
        DO I=1, NIND3
          IF (CKIND /= ' ') THEN
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >              SCRATCHBYTE
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST), IKIND
              WRITE (*,FMT=TRIM(FMTSTR)) (SCRATCHBYTE(JJ), JJ=J, J+IKIND-1)
            ENDDO
          ELSE
            READ(UNIT=ICHANNEL, REC=(IFIRST-1+I), ERR=94, IOSTAT=IOS)
     >              SCRATCH
            INDEX = (I-1) * IADRL
            DO J=1, MIN(IADRL, NREST)    
              IF (bWRINT) THEN
                INTSCRATCH = TRANSFER(SCRATCH(J), INTSCRATCH)
                WRITE (*,FMT=TRIM(FMTSTR)) INTSCRATCH
              ELSE 
                WRITE (*,FMT=TRIM(FMTSTR)) SCRATCH(J)
              ENDIF               
            ENDDO
          ENDIF
          NREST = NREST - IADRL
        ENDDO
        WRITE (*,*) 'FINISH'

      ELSE
        WRITE (0,*) ' ACTION ', ACTION( :IDX(ACTION)), ' NOT KNOWN'
        STOP 'ERROR IN STORAGE'

      ENDIF

      RETURN

C********** Error Stops ****************************************

   90 WRITE (0,*) ' ERROR WHEN OPENING MASS-STORAGE FILE'
      GOTO 99

   91 WRITE (0,*) ' ERROR WHEN READING MASS-STORAGE FILE (LABEL=91)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   92 WRITE (0,*) ' ERROR WHEN WRITING MASS-STORAGE FILE (LABEL=92)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   93 WRITE (0,*) ' ERROR WHEN WRITING MASS-STORAGE FILE (LABEL=93)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99

   94 WRITE (0,*) ' ERROR WHEN READING MASS-STORAGE FILE (LABEL=94)'
      WRITE (0,*) ' FILE NAME = ', FN
      GOTO 99


99    WRITE (0,'(A,I4)') 'Fortran Channel:', ICHANNEL
      WRITE (0,'(A,I4)') 'IOS=',IOS
      STOP 'ERROR IN SUBROUTINE STORAGE'

      END
      SUBROUTINE WRITMS(ICHANNEL, !file handle 
     >                         X, !(fortran) variable that should be saved
     >                      NDIM, !size of variable (length of array)
     >                      NAME, !index name in mass-storage file (must be <= 8 characters)
     >                  IKINDSTR, !kind string (use -1 for default kind)
     >                    IDUMMY, !unused integer dummy
     >                      IERR)
C************************************************************
C***  ROUTINE VON LASR KOESTERKE           8-Sep-1995 15:51:52
C************************************************************

      IMPLICIT NONE

      CHARACTER(8) :: BUFFER8, NAME
      INTEGER :: IKINDSTR, ICHANNEL, IDUMMY, IERR, NDIM
      REAL :: X

      !IKINDSTR contains a string containing type and format
      ! first character defines the type, numbers following specifiy the kind
      ! NOTE: for compartibility issues IKINDSTR must be read as integer by the subroutine
      IF (IKINDSTR /= -1) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') IKINDSTR    !transform from integer into character
      ELSE
        BUFFER8 = '        '
      ENDIF

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, BUFFER8, X, NDIM, 
     >              'WRITE', IERR)

      RETURN
      END
