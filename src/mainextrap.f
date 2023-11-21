      PROGRAM MAINextrap 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:18 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL extrap 
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
      SUBROUTINE AITKEN (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $                   ABXYZ,NFIRST,NLAST,TNEW,T,TOLD1,TOLD2,NOTEMP)
C***********************************************************************
C***  LOGARITHMIC AITKEN EXTRAPOLATION OF POPULATION NUMBERS
C***  THE ELECTRON DENSITY IS UPDATED ACCORDING TO THE NEW POPNUMBERS
C***********************************************************************
 
      DIMENSION RNE(ND),NCHARG(N)
      DIMENSION POPNEW(ND,N),POPNUM(ND,N),POP1(ND,N),POP2(ND,N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION TNEW(ND),T(ND),TOLD1(ND),TOLD2(ND)
      LOGICAL NOTEMP
      COMMON / COMNEGT / NEGT

C**********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION
      TMIN = 3000.
C**********************************************************************
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      RNE(L)=.0
 
C***  LOOP FOR EACH ELEMENT  -------------------------------
      DO 1 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=.0
      DO 2 J=NFIRNA,NLANA
C*** ZERO POPULATION NUMBERS ARE KEPT AT ZERO
      IF ( POPNUM(L,J).EQ.0. .OR. POP1(L,J).EQ.0. .OR. POP2(L,J)
     > .EQ.0. ) THEN
         POPNEW(L,J) = POPNUM(L,J)
         GOTO 8
         ENDIF
      A1   =ALOG10(POPNUM(L,J)/POP1(L,J))
C     THE AMPLIFICATION OF THE CORRECTIONS IS LIMITED TO A FACTOR OF 20.
      POPQ=POP1(L,J)/POP2(L,J)
      IF (POPQ .EQ. 1.) THEN 
         Q = 0.95
         ELSE
         A0   = ALOG10 (POPQ)
         Q = AMIN1(0.95,A1/A0 )
         ENDIF
      COR = A0/(1.-Q)
C***  THE CORRECTION ITSELF IS LIMITED TO A FACTOR OF 10.
      COR = AMIN1 (COR, 1. )
      COR = AMAX1 (COR,-1. )
      POPNEW(L,J)=POP2(L,J) * 10.**COR
    8 CONTINUE
      POPSUM=POPSUM+POPNEW(L,J)
    2 CONTINUE
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      POPSUM=POPSUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
      POPNEW(L,J)=POPNEW(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
    3 CONTINUE
 
C***  TEMPERATURE EXTRAPOLATION  -------------------------------
      IF ( .NOT. NOTEMP ) THEN
      IF ( T(L).EQ.0. .OR. TOLD1(L).EQ.0. .OR. TOLD2(L)
     > .EQ.0. ) THEN
         TNEWL = T(L)
         GOTO 18
         ENDIF
      A1   =ALOG10(T(L)/TOLD1(L))
C     THE AMPLIFICATION OF THE CORRECTIONS IS LIMITED TO A FACTOR OF 20.
      POPQ=TOLD1(L)/TOLD2(L)
      IF (POPQ .EQ. 1.) THEN 
         Q = 0.95
         ELSE
         A0   = ALOG10(POPQ)
         Q = AMIN1(0.95,A1/A0 )
         ENDIF
      COR = A0/(1.-Q)
C***  THE CORRECTION ITSELF IS LIMITED TO A FACTOR OF 10.
      COR = AMIN1 (COR, 1. )
      COR = AMAX1 (COR,-1. )
      TNEWL =TOLD2(L) * 10.**COR
C***  SET MINIMUM TEMPERATURE
      IF (TNEWL .LT. TMIN) THEN
         TNEW(L) = TMIN
         NEGT = NEGT + 1
         ELSE
         TNEW(L) = TNEWL
         ENDIF

   18 CONTINUE
      ENDIF
    1 CONTINUE
C***  ENDLOOPS  --------------------------------------------------------
 
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
      SUBROUTINE DECNOT (LAST,MODHIST,KARTE,NOTEMP,THISPRG)
C**********************************************************************
C***  IN CASE OF AN 'CONDITIONAL' NOTEMP-OPTION, THE OPTIONAL LIMIT
C***  FOR THE CORRECTION 'CORLIM' IS DECODED AND COMPARED TO THE
C***  LAST CORRECTION FOUND IN THE MODEL HISTORY (CORLAST)
C***  FORMAT OF THE INPUT OPTION:
C***  N0 TEMPERATURE CORRECTIONS WHILE COR. .GT.__________
C***  1234567890123456789012345678901234567890123456789012
C***  0000000001111111111222222222233333333334444444444555
C**********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LAST
      CHARACTER(LAST*8), INTENT(IN) :: MODHIST
      CHARACTER(80), INTENT(IN) :: KARTE
      CHARACTER(8), INTENT(IN) :: THISPRG
      LOGICAL, INTENT(OUT) :: NOTEMP
      
      CHARACTER(20) :: ACTPAR
      REAL :: CORLIM, CORLIM2, CORLAST
      INTEGER :: NPAR, I, J, IWRC, LASTCH, IEQUAL
      LOGICAL WRCONTFOUND, CORFOUND


C***  READ THE SPECIFIED LIMIT FROM THE INPUT CARD:
      CALL SARGC (KARTE, NPAR)
      CALL SARGV (KARTE, 7, ACTPAR)
      READ (ACTPAR, '(F20.0)') CORLIM
      IF (NPAR > 7) THEN
        CALL SARGV (KARTE, 8, ACTPAR)
        READ (ACTPAR, '(F20.0)') CORLIM2
        IF (CORLIM2 < CORLIM) THEN
          CORLIM2 = CORLIM
          WRITE (0,*) 'WARNING FROM SUBR. DECNOT'
          WRITE (0,*) 'SECOND LIMIT IS SMALLER THEN FIRST LMIT'
          WRITE (0,*) 'SECOND LIMIT IS IGNORED'
          WRITE (0,'(1X,2A)') 'KARTE=', KARTE( :70)
        ENDIF
      ELSE
        CORLIM2 = CORLIM
      ENDIF

c      READ (KARTE,3) CORLIM
c    3 FORMAT (42X,F10.0)

C***  LAST WRITTEN CHARACTER OF MODEL HISTORY
      LASTCH=LAST*8

C***  IF THIS PROGRAM IS NOT 'WRCONT': 
C***    GO BACK TO THE LAST WRCONT ENTRY
      WRCONTFOUND = .FALSE.

      IF (THISPRG( :6) /= 'WRCONT') THEN
         DO I=LASTCH,10,-1
           IWRC=I
           IF (MODHIST(I-5:I) == 'WRCONT') THEN
             WRCONTFOUND = .TRUE.
             LASTCH=IWRC
             EXIT
           ENDIF
         ENDDO
      ELSE
        WRCONTFOUND = .TRUE.
      ENDIF

      !exit subroutine if no WRCONT was found
      IF (.NOT. WRCONTFOUND) RETURN

C***  FIND LAST OCCURRENCE OF STRING 'COR.=' IN THE MODEL HISTORY
      CORFOUND = .FALSE.
      DO I=LASTCH,10,-1
        IEQUAL=I
        IF (MODHIST(I-4:I) == 'COR.=') THEN
          CORFOUND = .TRUE.
          EXIT
        ENDIF
      ENDDO

      !exit subroutine if no 'COR.=' was found or corrections are undefined
      IF ((.NOT. CORFOUND) .OR. 
     >     MODHIST(IEQUAL+1:IEQUAL+6) == 'UNDEF.') RETURN

C***  READ LAST CORRECTION:
      READ (MODHIST(IEQUAL+1:IEQUAL+8),4) CORLAST
    4 FORMAT (F8.0)

C***  NEW OPTION
C***  NOTEMP, IF CORLAST >  CORLIM2
C***  TEMP,   IF CORLAST <  CORLIM
C***  NO CHANGE, IF CORLIM <= CORLAST <= CORLIM2
      IF (CORLAST <= CORLIM) THEN
        NOTEMP = .FALSE.
      ELSEIF (CORLAST <= CORLIM2) THEN
C***  FIND LAST OCCURRENCE OF STRING 'COR.=' IN THE MODEL HISTORY
        LASTCH=LAST*8
        CORFOUND = .FALSE.
        DO I=LASTCH,10,-1
          IEQUAL=I
          IF (MODHIST(I-4:I) == 'COR.=') THEN
            CORFOUND = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT. CORFOUND) RETURN
C***  TRY TO FIND NEXT STRING 'TC='
        NOTEMP = .TRUE.
        DO J=I, LASTCH-3
          IF (MODHIST(J:J+2) == 'TC=') THEN
            NOTEMP = .FALSE.
            EXIT
          ENDIF
        ENDDO
      ELSE
        NOTEMP = .TRUE.
      ENDIF

C!!!  OLD PATH
C***  ALLOW TEMPERATURE CORRECTIONS, IF CORLAST .LT. CORLIM
C!!!      IF (CORLAST .LT. CORLIM) NOTEMP = .FALSE.

      RETURN
      END
C***  MAIN PROGRAM EXTRAP  *********************************************
      SUBROUTINE EXTRAP
C***********************************************************************
C***  THIS PROGRAM MAKES A (LOGARITHMIC) EXTRAPOLATION OF THE POPULATION
C***  NUMBERS FROM THE LAST THREE ITERATIONS
C***  AITKEN: EXTRAPOLATION OF POP.NUMBER FROM LAST 3 ITERATIONS AT
C***          SINGLE DEPTH POINT L
C***  NG: EXTRAPOLATION OF POP.NUMBER FROM LAST 3 ITERATIONS AT ALL
C***      DEPTH POINTS L (L=1...ND)
C***      IF (.NOT.NOTEMP): ALSO EXTRAPOLATION OF TEMPERATURE T(R)
C***  DEFAULT EXTRAPOLATION: ---  NG  ---
C***  THE NG-EXTRAPOLATION IS CONTROLLED FOR THE INCREASE OF POPNUMBERS
C***********************************************************************
 
C***  DEFINE ARRAY DIMENSIONS ******************************************
      INTEGER, PARAMETER :: MAXATOM =           26 
      INTEGER, PARAMETER :: NDIM    =         2560 
      INTEGER, PARAMETER :: NFDIM   = 2*NDIM + 400 
      INTEGER, PARAMETER :: MAXKONT =      NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR =         NDIM 
      INTEGER, PARAMETER :: MAXIND  =        45000 
      INTEGER, PARAMETER :: MAXFEIND  =       2500 
      INTEGER, PARAMETER :: NDDIM   =           89 
      INTEGER, PARAMETER :: MAXHIST =         4000 
C***  NUMBER OF ENTRYS STORED IN THE GAMMA HISTORY
      INTEGER, PARAMETER :: MAXGAHIST = 100

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 
      
C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 3200 
      COMMON / COMAUTO / LOWAUTO(MAXAUTO),WAUTO(MAXAUTO)
     $                  ,EAUTO(MAXAUTO),AAUTO(MAXAUTO),IONAUTO(MAXAUTO)
     $                  ,KRUDAUT(MAXAUTO)
      COMMON // NCHARG(NDIM),WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
     $ ,IONGRND(NDIM)
     $ ,MAINQN(NDIM),EINST(NDIM,NDIM),NOM(NDIM)
     $ ,ALTESUM(4,NDIM),COCO(4,MAXIND)
     $ ,ABXYZ(MAXATOM),KODAT(MAXATOM),ATMASS(MAXATOM),STAGE(MAXATOM)
     $ ,NFIRST(MAXATOM),NLAST(MAXATOM)
     $ ,RNE(NDDIM),RADIUS(NDDIM),W(NDDIM)
     $ ,TNEW(NDDIM),T(NDDIM),TOLD1(NDDIM),TOLD2(NDDIM), TOLD3(NDDIM)
     $ ,POPNEW(NDDIM,NDIM),POPNUM(NDDIM,NDIM),POP1(NDDIM,NDIM)
     $ ,POP2(NDDIM,NDIM), POP3(NDDIM,NDIM)
     $ ,ALPHA(MAXKONT),SEXPO(MAXKONT)
     $ ,ADDCON1(MAXKONT), ADDCON2(MAXKONT), ADDCON3(MAXKONT) 
     $ ,KONTNUP(MAXKONT),KONTLOW(MAXKONT)
     $ ,INDLOW(MAXIND),INDNUP(MAXIND)
     $ ,MDUMMY(MAXHIST)
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      CHARACTER(MAXHIST*8) :: MODHIST
      DIMENSION GAHIST(26,MAXGAHIST)
C***  PREVENT UNNECESSARY OUTPUT FROM SUBR. PRICORR:
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
      COMMON / COMNEGT / NEGT
      COMMON / COMITWA / ITWARN, ITMAX
      DATA NEGINTL, NEGT, ITWARN / 0, 0, 0 /
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      CHARACTER(255) :: HISTENTRY
      CHARACTER(40) :: BUFFER40
      CHARACTER(32) :: BUFFER32
      CHARACTER KARTE*80,MODHEAD*100
      CHARACTER LEVEL(NDIM)*10
      CHARACTER*10 ELEMENT(MAXATOM), ARG(5)
      CHARACTER*4 KEYCBB(MAXIND)
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(1), DIMENSION(NDDIM) :: CKONVER
      LOGICAL, DIMENSION(MAXATOM) :: TRACEELEM
      LOGICAL :: NGACC, NOTEMP, NOEXTEMP, NGWEIGHT, BEX4
 
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5     !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      LOGICAL :: BFEMODEL

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !write to MODEL file
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file
      
      
C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      CHARACTER TIM1*10, TIM2*10

c      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
c        CALL CLOCK(TIM1)
c      ELSE
c        CALL TIME(TIM1)
c      ENDIF

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
     >            'EXTRAP', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)
 
C***  READING OF THE MODEL FILE  ---------------------------------------
      CALL OPENMS(3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (3,JOBNUM,1,'JOBNUM  ', IERR)
      CALL READMS (3,MODHEAD,13,'MODHEAD ', IERR)
      CALL READMS(3,ND,1,'ND      ', IERR)
      IF (ND.GT.NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,RADIUS,ND,  'R       ', IERR)
      CALL READMS (3,POPNUM,ND*N,'POPNUM  ', IERR)
      CALL READMS (3,POP1  ,ND*N,'POP1    ', IERR)
      CALL READMS (3,POP2  ,ND*N,'POP2    ', IERR)
      CALL READMS (3,POP3  ,ND*N,'POP3    ', IERR)
C***  READING TEMPERATURE ARRAYS FOR EXTRAPOLATION OF T(R)
      CALL READMS (3,T    ,ND,   'T       ', IERR)
      IERR1=1
      CALL READMS (3,TOLD1,ND,   'TOLD1   ',IERR1)
      IERR2=1
      CALL READMS (3,TOLD2,ND,   'TOLD2   ',IERR2)
      IERR3=1
      CALL READMS (3,TOLD3,ND,   'TOLD3   ',IERR2)
      IF (((IERR1 .LT. 0) .AND. (IERR1 .NE. -10)) .OR. 
     $    ((IERR2 .LT. 0) .AND. (IERR2 .NE. -10)) .OR. 
     $    ((IERR3 .LT. 0) .AND. (IERR3 .NE. -10))) THEN
         CALL REMARK ('ERROR WHEN READING TOLD FROM MODEL FILE')
         STOP 'TOLD'
      ENDIF

C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3,ABXYZ,NATOM,'ABXYZ   ',IERR)
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
      CALL READMS(3,LAST,1,'MODHIST ', IERR)
      IF (LAST >= MAXHIST) THEN
            CALL REMARK ('MODHIST DIMENSION INSUFFICIENT')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,MODHIST,LAST,'MODHIST ', IERR)

C***  Read GAMMA History
      CALL READMS (3, GAHIST, 26*MAXGAHIST, 'GAHIST  ',IERR)
C***  Move Entries in GAMMA HISTORY
      DO I=MAXGAHIST-1, 1, -1
        DO J=1, 26
          GAHIST(J,I+1) = GAHIST(J,I)
        ENDDO
      ENDDO
C***  This Job is EXTRAP
      GAHIST(2,1) = 1.
C***  No Broyden Statistics and GAMMAS
        GAHIST( 3,1) = 0.
        GAHIST( 4,1) = 0.
        GAHIST( 5,1) = 0.
        GAHIST( 6,1) = 0.
        GAHIST(22,1) = 0.
        GAHIST(23,1) = 0.
        GAHIST(24,1) = 0.
        GAHIST(25,1) = 0.

C***  DEFAULTS (INPUT OPTIONS):  ***************************************
      LSPOP    = -1
      NGACC    = .TRUE.
      NOTEMP   = .FALSE.
      NOEXTEMP = .FALSE.
      BEX4     = .TRUE.
      NEWWRC   = 6

C***  DECODING INPUT CARDS  --------------------------------------------
      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
    7 READ(1,2, END=6) KARTE
    2 FORMAT (A)
      CALL SARGC (KARTE, NPAR)
      IF (NPAR .LT. 1) GOTO 7
      IF (NPAR .GT. 5) NPAR = 5
      DO 3 I=1, NPAR
    3 CALL SARGV (KARTE, I, ARG(I))

      IF ( KARTE(:8) .EQ. 'NO TEMPE' ) THEN
C                          ========
            NOTEMP=.TRUE.
            IF (KARTE(30:40) .NE. ' ')
     $          CALL DECNOT (LAST,MODHIST,KARTE,NOTEMP,'EXTRAP')
            GOTO 7
            ENDIF
      IF ( KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                          =========
            DECODE (80,4,KARTE) XL
    4       FORMAT (9X,F10.0)
            LSPOP=IFIX(XL)
            IF (LSPOP.EQ.0) LSPOP=1
            GOTO 7
            ENDIF
      IF (ARG(1) .EQ. 'EXTRAP' .AND. ARG(2) .EQ. 'NOTEMP') THEN
C                      ======                     ======
            NOEXTEMP = .TRUE.
            GOTO 7
            ENDIF
      IF (ARG(1) .EQ. 'EXTRAP' .AND. ARG(2) .EQ. 'NGWEIGHT') THEN
C                      ======                     ======
            NGWEIGHT = .TRUE.
            GOTO 7
            ENDIF

      IF ( KARTE(:6) .EQ. 'AITKEN' ) THEN
C                          ======
         NGACC = .FALSE.
         BEX4 = .FALSE.
         GOTO 7
         ENDIF
      IF (ARG(1) .EQ. 'NG3') THEN
C                      ===
            BEX4 = .FALSE.
            GOTO 7
            ENDIF

      IF (ARG(1) .EQ. 'NEWWRC') THEN
C                      ======
            READ (ARG(2),'(F10.0)') XL
            NEWWRC=IFIX(XL)
            GOTO 7
            ENDIF

      GOTO 7
    6 CLOSE (1)

      IF (NOEXTEMP) NOTEMP = .TRUE.

C***  NO TEMPERATURE EXTRAPOLATION IN CASE OF NON-EXISTING RECORDS
C***                               "TOLD."  (OLD MODEL FILE)
      IF ((IERR1 .EQ. -10) .OR. (IERR3 .EQ. -10)) NOTEMP=.TRUE.

      IF (BEX4) THEN
         CALL NG4 (ND, N, RNE, NCHARG, 
     >             POPNEW, POPNUM, POP1, POP2, POP3, 
     >             NATOM, ABXYZ, NFIRST, NLAST, 
     >             TNEW, T, TOLD1, TOLD2, TOLD3, NOTEMP, TRESH, NOUT)
      ELSE
C***  LOGARITHMIC EXTRAPOLATION  ***************************************
        IF (.NOT. NGACC) THEN
C***          THIS VERSION: NO EXTRAPOLATION OF T(R)
          DO 10 L=1,ND
   10     TNEW(L)=T(L)
          CALL AITKEN (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $                   ABXYZ,NFIRST,NLAST, TNEW, T, TOLD1, TOLD2, 
     $                   NOTEMP)
        ELSE
          CALL NG3 (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $      ABXYZ,NFIRST,NLAST,RADIUS,W,TNEW,T,TOLD1,TOLD2,
     $      NOTEMP, NGWEIGHT, NDONE, NTDONE)
        ENDIF
      ENDIF
 
C***  CYCLIC CHANGING OF THE ARRAY NAMES FOR THE LAST THREE ITERATIONS
      CALL CHANGE ('POP3    ','HELP    ', 3)
      CALL CHANGE ('POP2    ','POP3    ', 3)
      CALL CHANGE ('POP1    ','POP2    ', 3)
      CALL CHANGE ('POPNUM  ','POP1    ', 3)
      CALL CHANGE ('HELP    ','POPNUM  ', 3)
      CALL CHANGE ('TOLD3   ','THELP   ', 3)
      CALL CHANGE ('TOLD2   ','TOLD3   ', 3)
      CALL CHANGE ('TOLD1   ','TOLD2   ', 3)
      CALL CHANGE ('T       ','TOLD1   ', 3)
      CALL CHANGE ('THELP   ','T       ', 3)
 
C***  UPDATING THE MODEL FILE  -----------------------------------------
      JOBNUM=JOBNUM+1
C      IF (JOBNUM .GE. 1000) JOBNUM=JOBNUM-1000
      CALL WRITMS (3,JOBNUM,1,'JOBNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (3,POPNEW,ND*N,'POPNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (3,RNE,ND,'RNE     ',-1, IDUMMY, IERR)
      CALL WRITMS (3,TNEW,ND,'T       ',-1, IDUMMY, IERR)
 
C***  PRINTOUT  -------------------------------------------------------

C     Dummy setting for CKONVER array (mark all points as converged)
      DO L=1, ND
        CKONVER(L) = 'C'
      ENDDO
C     Dummy setting for TRACEELEM array 
      DO K=1, NATOM
        TRACEELEM(K) = .FALSE.
      ENDDO


      IF (LSPOP.GT.0)
     $  CALL PRIEX (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $              TNEW,NOTEMP)
      CALL INHIBIT (POPNEW, N, ND, NCHARG, RNE, 
     $              NATOM, ABXYZ, NFIRST, NLAST, 1.E-99)


C***  Commented Lines; full call as in STEAL
      CALL PRICORR (POPNEW, POPNUM, LEVEL, N, ND, MODHEAD, LSPOP,
C*   $              CORMAX, RTCMAX, JOBNUM, REDUCE, CKONVER,
     $              CORMAX, RTCMAX, JOBNUM, 1., CKONVER,     
C*   >              GAMMAC, GAMMAL, GAMMAR, GAMMAD,
     >               0.,      0.,      0.,     0., 
C*   $              T, TNEW, EPSILON, DELTAC, SMPOP,
     $              T, TNEW, 0.,       1.,    0.,  
C*   $              BUNLU,  DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, bTDIFFUS, TNDCORR,
     >              .NOT.NOTEMP, 0.,    0.,    0.,  0.,  .FALSE.,  0., 
C*   >              HNDCORFAC, GAHIST, MAXGAHIST, STHLP, ICMMODE,
     >              1., GAHIST, MAXGAHIST, .FALSE., 1, 
C*   >              IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >              0,                 0              ,
C*   >              TBTAU, TAUINT, NDOUT, NATOM, NOM, TRACEELEM)
     >                -1.,   -1., NOUT+1, NATOM, NOM, TRACEELEM)

      IF (NOTEMP) THEN
        GAHIST(26,1) = 0.
      ELSE
        GAHIST(26,1) = 1.
      ENDIF

      WRITE (*,'(1X,3A,1E8.2,A,I2)')
     >  'EXTRAP: POP. Numbers below treshhold and at outer points ', 
     >  'are not taken into ',
     >  'account : TRESH=',TRESH, '   NOUT=', NOUT

      IF (BEX4) THEN
        PRINT *,'This was a new NG4-extrapolation -------------------'
      ELSE
        IF (NGACC) THEN
          PRINT *,
     >      'This was an old NG3-extrapolation -------------------'
          PRINT *,'THE EXTRAPOLATION WAS DONE FOR ',NDONE,' OF ',ND,
     >            ' * ', N ,' POPNUMBERS AND ',NTDONE,' TEMPRERATURES'
        ELSE
          PRINT *,'This was an AITKEN-extrapolation  -------------'
        ENDIF
      ENDIF

C***  UPDATING THE MODEL HISTORY  --------------------------------------
      IF (CORMAX > 1.E-100) CORMAX=ALOG10(CORMAX)
      IF (BEX4) THEN
        WRITE(UNIT=BUFFER40, FMT=22) JOBNUM, CORMAX
C        ENCODE (32,22,MODHIST(LAST+1)) JOBNUM,CORMAX
C***  NOTE : The word Cor. must not be written in capitals in order to 
C***         distinguish from COR. of the Main Program STEAL
   22   FORMAT ('/',I7,'. EXTRAP4 Cor.=',F8.4)
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        IF (NOTEMP) THEN
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
        ENDIF
      ELSE
        IF (.NOT. NGACC) THEN
          WRITE(UNIT=BUFFER40, FMT=20) JOBNUM, CORMAX
   20     FORMAT ('/',I7,'. EXTRAP  COR.=',F8.4)
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        ELSE
          WRITE(UNIT=BUFFER40, FMT=30) JOBNUM, CORMAX
   30     FORMAT ('/',I7,'. EXTRAP(NG)  COR.=',F8.4)
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        ENDIF
        IF (NGACC) THEN
          WRITE(UNIT=BUFFER32, FMT=31) NDONE, ND, N
   31     FORMAT (I5,' DPS OF ',I2,'*',I3,' DONE')
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,32,BUFFER32)
          IF (.NOT. NOTEMP) THEN
            WRITE(UNIT=BUFFER32, FMT=32) RTCMAX, NTDONE
   32       FORMAT ('   TC=',F8.4, ' DONE FOR ',I2,' DPS')
            CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,32,BUFFER32)
          ENDIF
        ELSE
          IF (NOTEMP) THEN
            CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
          ENDIF
        ENDIF
      ENDIF
      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      CALL WRITMS (3, GAHIST, 26*MAXGAHIST, 'GAHIST  ',-1, IDUMMY, IERR)
 
C***  write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >              ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)      
 
   99 CALL CLOSMS (3, IERR)

C***  NEXT JOB: REPEAT-CYCLE if NEWWRC NE 1
      IF (NEWWRC .NE. 1) THEN
        CALL REMARK ('EXTRAP: NEXTJOB=REPEAT')
        PRINT *, 'EXTRAP: NEXTJOB=REPEAT'
        CALL JSYMSET ('G1','REPEAT')
      ELSE
C***  NEXT JOB: WRCONT-CYCLE if NEWWRC EQ 1
        CALL REMARK ('EXTRAP: NEXTJOB=WRCONT')
        PRINT *, 'EXTRAP: NEXTJOB=WRCONT'
        CALL JSYMSET ('G1','WRCONT')
      ENDIF
      CALL JSYMSET ('G3','MOREJOBS')

      CALL JSYMSET ('G0','0')

      CALL STAMP (OPSYS, 'EXTRAP', TIM1)

      STOP 'O.K.'
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
      SUBROUTINE GAUNTFF (GIII,NCHARGE,XLAM,TEMP)
C***********************************************************************
C***  THE MODULE COMPUTES THE G-III FACTOR FOR THERMAL BREMSSTRAHLUNG
C***  GIII = GAUNT FACTOR
C***  NCHARGE = REST CHARGE OF ION ( =1 FOR H+ AND HE+, =2 FOR HE++ )
C***  XLAM    = WAVELEGTH IN ANGSTROEM
C***  TEMP    = TEMPERATURE IN KELVIN
C***  THE GAUNT FACTORS ARE COMPILED FROM THE PUPLICATIONS OF BERGER (1956)
C***  APJ 124,P550  AND KARZAS AND LATTER (1961) APJ SUPPL 6,P167 (FIG 3,4,5)
C***  TO GET THE HERE TABULATED VALUES BICUBIC SPLINE INTERPOLATION WAS USED
C***  --------  TABULATED RANGE FOR NCHARGE = 1  -----------------------
C***  TEMPERATURE (KELVIN)             :  1577  ...  157 700
C***  WAVELENGTH (ANGSTROEM)           :   120  ...  1.2 E6
C***********************************************************************
 
      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP

      PARAMETER ( MAXION = 27 )
 
      REAL A(21,21)
      DIMENSION ZLOG(MAXION)
 
      DATA      ZLOG / 0.     , 0.30103, 0.47712, 0.60206, 0.69897,
     >                 0.77815, 0.84510, 0.90309, 0.95424, 1.0    ,
     >                 1.04139, 1.07918, 1.11394, 1.14613, 1.17609, 
     >                 1.20412, 1.23045, 1.25527, 1.27875, 1.30103,
     >                 1.32222, 1.34232, 1.36173, 1.38021, 1.39794,
     >                 1.41497, 1.43136 /
      DATA ((A(I,J),I=1,21),J=1,6) /
     $ 1.331, 1.274, 1.232, 1.200, 1.177, 1.158, 1.143, 1.130, 1.120,
     * 1.115, 1.112, 1.108, 1.103, 1.100, 1.100, 1.101, 1.098, 1.090,
     $ 1.080, 1.069, 1.057,
     $ 1.412, 1.347, 1.299, 1.258, 1.214, 1.173, 1.145, 1.132, 1.126,
     * 1.122, 1.118,1.114, 1.110, 1.106, 1.102, 1.099, 1.098, 1.096,
     * 1.090, 1.078, 1.060,
     $ 1.503, 1.418, 1.350, 1.289, 1.230, 1.178, 1.146, 1.133, 1.129,
     * 1.125, 1.121, 1.117, 1.113, 1.109, 1.104, 1.100, 1.100, 1.100,
     $1.098, 1.084, 1.062,
     $ 1.601, 1.489, 1.394, 1.311, 1.239, 1.182, 1.148, 1.135, 1.131,
     $ 1.128, 1.123, 1.119, 1.115, 1.111, 1.106, 1.102, 1.102, 1.105,
     $ 1.102, 1.087, 1.064,
     $ 1.704, 1.566, 1.443, 1.338, 1.254, 1.193, 1.157, 1.142, 1.136,
     $ 1.131, 1.126, 1.121, 1.118, 1.114, 1.109, 1.105, 1.105, 1.104,
     $ 1.103, 1.087, 1.064,
     $ 1.807, 1.652, 1.509, 1.386, 1.288, 1.218, 1.176, 1.156, 1.147,
     * 1.139, 1.131, 1.126, 1.122, 1.119, 1.114, 1.110, 1.109, 1.107,
     $ 1.101, 1.085, 1.063 /
      DATA ((A(I,J),I=1,21),J=7,12) /
     $ 1.911, 1.750, 1.600, 1.466, 1.352, 1.264, 1.207, 1.179, 1.165,
     $ 1.153, 1.142, 1.129, 1.126, 1.126, 1.121, 1.116, 1.111, 1.106,
     $ 1.096, 1.081, 1.060,
     $ 2.016, 1.849, 1.702, 1.562, 1.431, 1.322, 1.248, 1.209, 1.188,
     $ 1.172, 1.158, 1.142, 1.136, 1.135, 1.130, 1.122, 1.114, 1.104,
     $ 1.091, 1.076, 1.059,
     $ 2.121, 1.947, 1.790, 1.644, 1.500, 1.375, 1.288, 1.241, 1.216,
     $ 1.196, 1.179, 1.158, 1.149, 1.145, 1.137, 1.128, 1.116, 1.103,
     $ 1.090, 1.076, 1.061,
     $ 2.226, 2.046, 1.876, 1.713, 1.561, 1.429, 1.335, 1.281, 1.250,
     $ 1.225, 1.204, 1.177, 1.164, 1.157, 1.144, 1.132, 1.119, 1.105,
     $ 1.091, 1.078, 1.064,
     $ 2.331, 2.146, 1.970,1.802, 1.644,  1.507, 1.407, 1.343, 1.301,
     $ 1.265, 1.234, 1.201, 1.182, 1.172, 1.145, 1.138, 1.122, 1.105,
     $ 1.090, 1.076, 1.061,
     $ 2.500, 2.307, 2.123, 1.947, 1.766, 1.623, 1.513, 1.435, 1.374,
     $ 1.318, 1.268, 1.229, 1.204, 1.189, 1.168, 1.147, 1.126, 1.105,
     $ 1.085, 1.067, 1.046 /
      DATA ((A(I,J),I=1,21),J=13,18) /
     $ 2.659, 2.460, 2.268,2.084, 1.907, 1.755, 1.633, 1.539, 1.458,
     $ 1.379, 1.309, 1.262, 1.229, 1.208, 1.185, 1.158, 1.132, 1.106,
     $ 1.081, 1.055, 1.026,
     $ 2.809, 2.604, 2.407, 2.217, 2.035, 1.873, 1.740, 1.634, 1.538,
     $ 1.442, 1.356, 1.300, 1.259, 1.230, 1.202, 1.172, 1.141, 1.112,
     $ 1.082, 1.049, 1.012,
     $ 2.953, 2.743, 2.541, 2.345, 2.156, 1.973, 1.830, 1.714, 1.610,
     $ 1.506, 1.410, 1.344, 1.294, 1.255, 1.221, 1.187, 1.155, 1.124,
     $ 1.090, 1.050, 1.005,
     $ 3.091, 2.878, 2.670, 2.469, 2.276, 2.090, 1.923, 1.797, 1.683,
     $ 1.573, 1.471, 1.394, 1.333, 1.285, 1.242, 1.204, 1.170, 1.136,
     $ 1.099, 1.054, 1.001,
     $ 3.261, 3.042, 2.829, 2.622, 2.421, 2.227, 2.039, 1.895, 1.766,
     $ 1.646, 1.537, 1.451, 1.378, 1.319, 1.268, 1.224, 1.184, 1.146,
     $ 1.104, 1.055, 0.997,
     $ 3.432, 3.208, 2.989, 2.776, 2.569, 2.368, 2.175, 2.008, 1.858,
     $ 1.725, 1.610, 1.514, 1.429, 1.358, 1.299, 1.246, 1.199, 1.155,
     $ 1.108, 1.054, 0.993 /
      DATA ((A(I,J),I=1,21),J=19,21) /
     $ 3.593, 3.365, 3.142, 2.924, 2.711, 2.505, 2.305, 2.127, 1.956,
     $ 1.811, 1.689, 1.583, 1.485, 1.403, 1.335, 1.271, 1.216, 1.166,
     $ 1.114, 1.057, 0.994,
     $ 3.747, 3.515, 3.289, 3.067, 2.850, 2.638, 2.432, 2.233, 2.055,
     $ 1.903, 1.776, 1.660, 1.549, 1.455, 1.374, 1.301, 1.239, 1.183,
     $ 1.127, 1.066, 1.001,
     $ 3.895, 3.661, 3.431, 3.205, 2.984, 2.768, 2.558, 2.354, 2.152,
     $ 2.000, 1.873, 1.743, 1.620, 1.513, 1.415, 1.336, 1.269, 1.209,
     $ 1.151, 1.083, 1.018 /
 
      DATA NHELP / 0 /

C***  LOGARITHM OF TEMP AND XLAM
      TEMPLOG=ALOG10(TEMP)
      XLAMLOG=ALOG10(XLAM)
      GOTO 1
 
C***  ENTRY FOR CALLING WITH LOGARITHMIC ARGUMENTS
      ENTRY GFFLOG (GIII,NCHARGE,XLOG,TLOG)
      XLAMLOG=XLOG
      TEMPLOG=TLOG
    1 CONTINUE
 
      IF (NCHARGE .LE. 0 .OR. NCHARGE .GT. MAXION) THEN
         CALL REMARK ('NCHARGE OUTSIDE VALID RANGE')
         STOP 'ERROR in Subr. GAUNTFF'
         ENDIF
 
C***  AT THE FIRST CALL WITHIN A MAIN PROGRAMM, THE COUNTERS ARE SET ZERO.
C***  THESE COUNTERS INDICATE HOW OFTEN THIS ROUTINE HAS BEEN CALLED
C***  WITH PARAMETERS OUTSIDE OF THE TABULATED RANGE.
      IF (NHELP .EQ. 0) THEN
         NHELP=1
         NTUP=0
         NTLOW=0
         NFUP=0
         NFLOW=0
      ENDIF
 
C***  SCALE TEMP AND XLAM WITH THE REST IONIC CHARGE
      XZLOG=XLAMLOG+2.*ZLOG(NCHARGE)
      TZLOG=TEMPLOG-2.*ZLOG(NCHARGE)
 
C***  CALCULATE TABULAR INDICES (BROKEN NUMBERS)
      AT=10.*TZLOG    -30.98364
      AF=31.3944035-5.*XZLOG
 
      MT=AT
C***  LOWER TEMPERATURE BOUNDARY
      IF (MT.LT.1) THEN
         MT=1
         AT=1.
         NTLOW=NTLOW+1
         ENDIF
C***  UPPER TEMPERATURE BOUNDARY
      IF (MT.GT.20) THEN
         MT=20
         AT=20.
         NTUP=NTUP+1
         ENDIF
 
      MF=AF
C***  LOWER FREQUENCY BOUNDARY
      IF (MF.LT.1) THEN
         MF=1
         AF=1.
         NFLOW=NFLOW+1
         ENDIF
C***  UPPER FREQUENCY BOUNDARY
      IF (MF.GT.20) THEN
         GIII=1.
C***     THIS APPROXIMATION BECOMES WORSE IN THE X-RAY REGION.
C***     BEYOND NU=10**17 BETTER VALUES SHOULD BE CALCULATED
         IF (MF .GT. 24) NFUP=NFUP+1
         RETURN
         ENDIF
 
C***  INTERPOLATION
      XF=AF-MF
      GF1=A(MF,MT)+XF*(A(MF+1,MT)-A(MF,MT))
      GF2=A(MF,MT+1)+XF*(A(MF+1,MT+1)-A(MF,MT+1))
      GIII=GF1+(AT-MT)*(GF2-GF1)

      IF (GIII < 0.) THEN
        WRITE (0,*) "XLAM=", XLAM, " TEMP=", TEMP
        WRITE (0,*) "TZLOG=", TZLOG, " ZLOG(NCHARGE)=", ZLOG(NCHARGE)
        WRITE (0,*) "MF=", MT, " MT=", MT
        WRITE (0,*) "GF1=", GF1, " GF2=", GF2
        WRITE (0,*) "AT=", AT, " AT-MT=", AT-MT
        WRITE (0,*) "negativer Gauntfaktor: ", GIII
        STOP "------- FATAL ERROR: SOMETHING IS TOTALLY WRONG -------"
      ENDIF
 
c      giii = 1.

      RETURN
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
      SUBROUTINE NG3 (ND,N,RNE,NCHARG,POPNEW,POP,POP1,POP2,NATOM,ABXYZ,
     $         NFIRST,NLAST,RADIUS,W,TNEW,T,TOLD1,TOLD2,NOTEMP,NGWEIGHT,
     $         NDONE, NTDONE)
C***********************************************************************
C***  NG'S ACCELERATION METHOD (K.C.NG 1974, CHEM. PHYS. 61, 2680):
C***  EXTRAPOLATION BY CONSIDERING THE LAST THREE ITERATIONS
C***                      (NG'S PAPER: LAST FOUR ITERATIONS!!!)
C***  === VERSION ===: EXTRAPOLATION OF LOGARITHMIC POP. NUMBERS
C***                   IF (.NOT.NOTEMP): ALSO EXTRAPOLATION OF TEMPERATURE
C***  "INTEGRATION" WEIGHTS:  
C***    LOGICAL SWITCH: NGWEIGHT
C***      .TRUE. = WEIGHTED WITH 1/VARIABLE (OPTION: NGWEIGHT)
C***     .FALSE. = UNWEIGHTED               (DEFAULT)
C***
C***  TESTED : INCREASING CORRECTIONS IN GENERAL (A>1.)
C***             ==> NO EXTRAP DONE AT ALL
C***           INCREASING CORRECTIONS ONLY FOR SOME POPNUMBERS
C***             ==> NO EXTRAP DONE FOR THIS POPNUMBER
C***********************************************************************
 
      DIMENSION POPNEW(ND,N),POP(ND,N),POP1(ND,N),POP2(ND,N)
      DIMENSION RNE(ND),NCHARG(ND),RADIUS(ND),W(ND)
      DIMENSION TNEW(ND),T(ND),TOLD1(ND),TOLD2(ND)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      LOGICAL NOTEMP, NGWEIGHT
      COMMON / COMNEGT / NEGT

C***********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION OF T(R)
      TMIN=3000.
C***********************************************************************

C***  INITIALIZE COUNTER FOR NEGATIVE TEMPERATURE WARNING
      NEGT=0
      NDONE = 0
      NTDONE = 0

      DO 101 J=1,N
      DO 102 L=1,ND
        IF ( POP(L,J) .LT. 1.E-99 ) POP(L,J)=1.E-99
        IF ( POP1(L,J) .LT. 1.E-99 ) POP1(L,J)=1.E-99
        IF ( POP2(L,J) .LT. 1.E-99 ) POP2(L,J)=1.E-99
102   CONTINUE
101   CONTINUE

C***  INITIALIZATION OF "INTEGRATION" WEIGHTS (UNWEIGHTED)
      DO 2 L=1, ND
      W(L)=1.
    2 CONTINUE


***  LOOP FOR EACH LEVEL  ---------------------------------------------
      DO 10 J=1,N
      A1=0.0
      C1=0.0

C***  INITIALIZATION OF "INTEGRATION" WEIGHTS
      IF (NGWEIGHT) THEN
         DO 1 L=1,ND
         IF (POP(L,J) .GT. 0.) THEN
            W(L)=1./POP(L,J)
            ELSE
            W(L) = 1.
            ENDIF
    1    CONTINUE
         ENDIF

C***  CALCULATION OF THE EXTRAPOLATION CONSTANTS FOR NG'S ACCELERATION METHOD
      DO 5 L=1,ND
      WL=W(L)
      D01=ALOG10(POP(L,J)*POP2(L,J)/POP1(L,J)/POP1(L,J))
      DN=ALOG10(POP(L,J)/POP1(L,J))
      A1=A1+D01*D01*WL
      C1=C1+D01*DN*WL
    5 CONTINUE
      IF (A1 .EQ. 0.) THEN
         A=0.0
         ELSE
         A=C1/A1
         ENDIF
      DO 12 L=1,ND
        POPNEW(L,J) = POP(L,J)
   12 CONTINUE
C***  CONTROL IN GENERAL
      IF (A .LT. 1.) THEN
        DO 8 L=1,ND
C***      CONTROL FOR EACH POPNUMBER
          IF (ABS(POP2(L,J)-POP1(L,J)) .GE. ABS(POP1(L,J)-POP(L,J)))
     >    THEN
            NDONE = NDONE + 1
            POPNEW(L,J) = POP(L,J)*(POP1(L,J)/POP(L,J))**A
C***        CORRECTION NOT REATER THAN 10 TIMES OF LAST CORRECTION
            IF (ABS(POPNEW(L,J)-POP(L,J)) .GT. 
     >          10. *ABS(POP(L,J)-POP1(L,J))) THEN
              POPNEW(L,J) = POP(L,J) + 10.*(POP(L,J)-POP1(L,J))
            ENDIF
          ENDIF
    8   CONTINUE
      ENDIF
   10 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      DO 20 L=1,ND
      RNE(L)=0.0
      DO 20 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=0.0
      DO 15 J=NFIRNA,NLANA
      POPSUM=POPSUM+POPNEW(L,J)
   15 CONTINUE
      POPSUM=POPSUM/ABXYZ(NA)
      DO 17 J=NFIRNA,NLANA
      POPNEW(L,J)=POPNEW(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
   17 CONTINUE
   20 CONTINUE
 
C***  IN CASE OF TEMPERATURE CORRECTIONS: ALSO EXTRAPOLATION OF T(R)
      IF (.NOT. NOTEMP) THEN

C***     INITIALIZATION OF "INTEGRATION" WEIGHTS
         DO 11 L=1,ND
           IF (NGWEIGHT) THEN
             W(L)=1. / T(L)
           ELSE
             W(L) = 1.
           ENDIF
   11    CONTINUE
         A1=0.0
         C1=0.0
C***     CALCULATION OF THE EXTRAPOLATION CONSTANTS
         DO 25 L=1,ND
         WL=W(L)
         DN=T(L)-TOLD1(L)
         D01=DN-(TOLD1(L)-TOLD2(L))
         A1=A1+D01*D01*WL
         C1=C1+D01*DN*WL
   25    CONTINUE
      IF (A1 .EQ. 0.) THEN
         A=0.0
         ELSE
         A=C1/A1
         ENDIF
C!!!      WRITE(*,*) 'TESTOUTPUT TEMPERATUR-EXTRAP'
C!!!      WRITE(*,*) 'A=',A
C***  NO TEMPERATURE EXTRAPOLATION IF CORRECTION INCREASES IN GENERAL (A >= 1.)
         IF (A .LT. 1.) THEN
           DO 28 L=1, ND
C***  CONTROL FOR EACH DEPTH=POINT
C!!!      WRITE(*,*) 'L=',L,' T=',T(L),
C!!!     >           ' TOLD1=',TOLD1(L),' TOLD2=',TOLD2(L)
             IF (ABS(TOLD1(L)-TOLD2(L)) .GE. ABS(T(L)-TOLD1(L))) THEN
C!!!      WRITE(*,*) 'DONE'
               NTDONE = NTDONE + 1
               TNEW(L) = (1.-A) * T(L) + A * TOLD1(L)
C***  CORRECTION NOT GREATER THAN 10 TIMES OF THE LAST CORRECTION
               IF (ABS(TNEW(L) - T(L)) .GT. 10.*ABS(T(L) - TOLD1(L)))
     >           THEN
C!!!      WRITE(*,*) 'CORRECTION BRAKED'
                 TNEW(L) = T(L) + 10. * (T(L) - TOLD1(L))
               ENDIF
               IF (TNEW(L) .LT. TMIN) THEN
                 TNEW(L) = TMIN
                 NEGT = NEGT + 1
               ENDIF
             ELSE
               TNEW(L) = T(L)
             ENDIF
   28     CONTINUE
        ELSE
          DO 29 L=1, ND
            TNEW(L) = T(L)
   29     CONTINUE
        ENDIF
      ELSE
C***  NO TEMPERATURE EXTRAPOLATION
        DO 38 L=1, ND
          TNEW(L) = T(L)
   38   CONTINUE
      ENDIF

      RETURN
      END
      SUBROUTINE NG4 (ND, N, RNE, NCHARG, 
     >                POPNEW, POP, POP1, POP2, POP3, 
     >                NATOM, ABXYZ, NFIRST, NLAST,
     >                TNEW, T, TOLD1, TOLD2, TOLD3, NOTEMP, TRESH, NOUT)
C***********************************************************************
C***  NG'S Acceleration method (K.C.NG 1974, CHEM. PHYS. 61, 2680):
C***  Extrapolation by considering the last four iterations
C***  === Version ===: The acceleration is now the same as suggested 
C***                   by Ivan Hubeny :
C***                   a) All popnumbers at all depths are considered
C***                      to examine one acceleration factor. If Temperature 
C***                      corrections are taken into account, T is also
C***                      included
C***                   b) No logarithmic treatment
C***                   c) Popnumbers are always weigted 
C***********************************************************************
 
      DIMENSION POPNEW(ND,N), POP(ND,N), POP1(ND,N)
      DIMENSION POP2(ND,N), POP3(ND,N)
      DIMENSION RNE(ND),NCHARG(ND)
      DIMENSION TNEW(ND), T(ND), TOLD1(ND), TOLD2(ND), TOLD3(ND)
      DIMENSION ABXYZ(NATOM), NFIRST(NATOM), NLAST(NATOM)
      LOGICAL NOTEMP
      COMMON / COMNEGT / NEGT

C***********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION OF T(R)
      TMIN=3000.
C***********************************************************************

C***  INITIALIZE COUNTER FOR NEGATIVE TEMPERATURE WARNING
      NEGT=0

C***  Popnumbers smaller than TRESH are not accounted for
      TRESH = 1.E-15

C***  Do not EXTRAP the outer points
      N1 = 5
      IF (N1 .GT. ND/2) THEN
        N1 = ND/2
      ENDIF
      NOUT = N1

      DO J=1,N
        DO L=1,ND
          IF ( POP(L,J) .LT. 1.E-99 ) POP(L,J)=1.E-99
          IF ( POP1(L,J) .LT. 1.E-99 ) POP1(L,J)=1.E-99
          IF ( POP2(L,J) .LT. 1.E-99 ) POP2(L,J)=1.E-99
          IF ( POP3(L,J) .LT. 1.E-99 ) POP3(L,J)=1.E-99
        ENDDO
      ENDDO

***  LOOP FOR ALL LEVELS AND POPNUMBERS  ---------------------------------------------
      A1=0.0
      C1=0.0
      DO J=1,N
        DO L=1,ND
          IF (POP(L,J) .LT. TRESH .OR. L .LT. N1) CYCLE
          WT = 1. / POP(L,J)
          D0 = POP(L,J) - POP1(L,J)
          D1 = D0 - POP1(L,J) + POP2(L,J)
          D2 = D0 - POP2(L,J) + POP3(L,J)
          A1 = A1 + WT*D1*D1
          B1 = B1 + WT*D1*D2
          B2 = B2 + WT*D2*D2
          C1 = C1 + WT*D0*D1
          C2 = C2 + WT*D0*D2
        ENDDO
      ENDDO
C*** Temperature corrections
      DO L=N1, ND
        WT = 1. / T(L)
        D0 = T(L) - TOLD1(L)
        D1 = D0 - TOLD1(L) + TOLD2(L)
        D2 = D0 - TOLD2(L) + TOLD3(L)
        A1 = A1 + WT*D1*D1
        B1 = B1 + WT*D1*D2
        B2 = B2 + WT*D2*D2
        C1 = C1 + WT*D0*D1
        C2 = C2 + WT*D0*D2
      ENDDO

C***  ENDLOOP  ---------------------------------------------------------

      AB = B2*A1 - B1*B1
      IF (AB .NE. 0.) THEN
        A = (B2*C1 - B1*C2) / AB
        B = (A1*C2 - B1*C1) / AB
      ELSE
        WRITE (0,*) 'No EXTRAP (NG4) executed'
        WRITE (*,*) 'No EXTRAP (NG4) executed'
        A = 0.
        B = 0.
      ENDIF

      DO J=1,N
        DO L=1,ND
          IF (POP(L,J) .LT. TRESH .OR. L .LT. N1) THEN
            POPNEW(L,J) = POP(L,J)
          ELSE
            POPNEW(L,J) = (1.-A-B)*POP(L,J) + 
     >                    A*POP1(L,J) + B*POP2(L,J)
          ENDIF
        ENDDO
      ENDDO
C*** Temperature corrections
      IF (.NOT. NOTEMP) THEN
        DO L=1, ND
          IF (L .LT. N1) THEN
            TNEW(L) = T(L)
          ELSE
            TNEW(L) = (1.-A-B)*T(L) + 
     >                  A*TOLD1(L) + B*TOLD2(L)
          ENDIF
        ENDDO
      ELSE
        DO L=1, ND
          TNEW(L) = T(L)
        ENDDO
      ENDIF

C***  Popnumbers are scaled to ensure number conservation
C***  Electron density is updated
      DO L=1,ND
        RNE(L)=0.0
        DO NA=1,NATOM
          NFIRNA=NFIRST(NA)
          NLANA=NLAST(NA)
          POPSUM=0.0
          DO J=NFIRNA,NLANA
            POPSUM=POPSUM+POPNEW(L,J)
          ENDDO
          POPSUM=POPSUM/ABXYZ(NA)
          DO J=NFIRNA,NLANA
            POPNEW(L,J)=POPNEW(L,J)/POPSUM
            RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
          ENDDO
        ENDDO
      ENDDO

C*** Check for very low temperatures
      IF (.NOT. NOTEMP) THEN
        DO L=1, ND
          IF (TNEW(L) .LT. TMIN) THEN
            TNEW(L) = TMIN
            NEGT = NEGT + 1
          ENDIF
        ENDDO
      ENDIF

      RETURN
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
      SUBROUTINE PRICORR (POPNUM,POP1,LEVEL,N,ND,MODHEAD,LSPOP,
     >                    CORMAX,RTCMAX,JOBNUM,REDUCE,CKONVER,
     >                    GAMMAC,GAMMAL,GAMMAR, GAMMAD, 
     >                    TOLD, TNEW, EPSILON, DELTAC, SMALLPOP, BUNLU, 
     >                    DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB,
     >                    BTDIFFUS, TNDCORR, HNDCORFAC, GAHIST,  
     >                    MAXGAHIST, STHLP, ICMMODE, IWARN_NEG_XJCAPP, 
     >                    IWARN_NEG_XJLAPP, TBTAU, TAUINT, NDOUT, 
     >                    NATOM, NOM, TRACEELEM)
C*******************************************************************************
C***  PRINTOUT OF CORRECTION FACTORS OF POPULATION NUMBERS
C***  RELATIVE TO THE LAST ITERATION
C***  Called from: STEAL, EXTRAP
C*******************************************************************************
      INTEGER, INTENT(IN) :: ND, N, NDOUT, NATOM
 
      INTEGER, DIMENSION(N) :: NOM
      DIMENSION POPNUM(ND,N),  POP1(ND,N)
      DIMENSION TOLD(ND), TNEW(ND)
      DIMENSION GAHIST(26,MAXGAHIST)
      CHARACTER(1), DIMENSION(ND) :: CKONVER
      CHARACTER(10), DIMENSION(N) :: LEVEL
      CHARACTER MODHEAD*100
      CHARACTER PRILINE*130,NUMBER*12
      CHARACTER*10 LEVMAX,LEVMIN,LEVMA2,LEVMI2
      LOGICAL, DIMENSION(NATOM) :: TRACEELEM
      LOGICAL BUNLU, BTDIFFUS, STHLP, bConsiderDepth
      INTEGER :: NDMIN, ICMMODE
C***  COMMON /GIIIERR/ COUNTS THE ERROR MESSAGES FROM SUBR. GAUNTFF
      COMMON / GIIIERR /  NTUP,NTLOW,NFUP,NFLOW,NHELP
C***  NEGINTL COUNTS THE ERROR MESSAGES  FROM SUBR. LINPOP
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
C***  Warning counter for negative continuum opacities <-- LINPOP
      COMMON / COMNEGC / NEGINTC,INEGMINC,INEGMAXC,LNEGMINC,LNEGMAXC
C***  ITWARN COUNTS THE NOT CONVERGED NEWTON-RAPHSON ITERATIONS FROM
C***  SUBR. LINPOP
      COMMON / COMITWA / ITWARN, ITMAX

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

C*********************************************************
C***  SET NDMIN = FIRST DEPTH INDEX WHICH IS ACCOUNTED FOR
C***               IN THE CONVERGENCE CRITERION
C**********************************************************
      IF (NDOUT > 0) THEN
C***  CARDS option added which allows manual adjustment of NDMIN, by ansander at 02-Oct-2016
        NDMIN = NDOUT
      ELSE
C***      NDMIN = 9
C***  New Version, Lars 30-Jul-1998 10:15:26
        NDMIN = ND / 8 + 1
      ENDIF

C**********************************************************************
C***  SET SMALLPOP : SMALLER POPNUMBERS ARE NOT ACCOUNTED FOR
C***                 IN THE CONVERGENCE CRITERIUM!
C***  NOTE: THIS CONVERGENCE CRITERION SHOULD HARMONIZE WITH THE
C***        CONVERGENCE CRITERION IN SUBR. LINPOP !!
C**********************************************************************

C**********************************************************************
C***  TABLE OF RELATIVE CORRECTIONS (IF REQUESTED)
C**********************************************************************
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT ( /, 1X, A, 15X, 'JOB NO.', I7, //, 10X,
     $ 'RATIO OF POPULATION NUMBERS TO THOSE OF THE LAST ITERATION:')
      IF (LSPOP.LT.1) GOTO 8
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,' DEPTH',10(2X,A10))
      IF (N.LT.J2) PRINT 5
    5 FORMAT (1X)
      PRINT 5
 
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      WRITE (UNIT=PRILINE, FMT='(I6)') L
 
      DO 12 J=J1,J2
        IF (POP1(L,J) .NE. .0) THEN
          WRITE (UNIT=NUMBER, FMT='(F12.4)') POPNUM(L,J)/POP1(L,J)
        ELSE
          NUMBER='    INFINITE'
        ENDIF
        I=7+(J-J1)*12
        PRILINE(I:I+11)=NUMBER
   12 CONTINUE
 
      PRINT 13,PRILINE
   13 FORMAT (A)
    3 CONTINUE
 
      IF (J2.EQ.N) GOTO 8
      J1=J1+10
      GOTO 4
    8 CONTINUE
 
C***  PRINTOUT OF THE TEMPERATURE CORRECTIONS
      IF (BUNLU .AND. LSPOP .GT. 0) THEN
      PRINT 21
   21 FORMAT (////,40X,'RELATIVE TEMPERATURE CORRECTIONS',/,40X,32('='),
     $   //,40X,'DEPTH INDEX      REL.TEMP.CORR.',/)
      DO 22 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 22
      RTC=(TNEW(L)-TOLD(L))/TOLD(L)
      PRINT 23,L,RTC
   23 FORMAT (45X,I3,F15.4)
   22 CONTINUE
      ENDIF
 

C**********************************************************************
C***  CALCULATE MINIMUM AND MAXIMUM OF ALL CORRECTION FACTORS
C**********************************************************************
      QMAX=1.
      QMA2=1.
      QMIN=1.
      QMI2=1.
      LEVMAX='NONE'
      LEVMA2='NONE'
      LEVMIN='NONE'
      LEVMI2='NONE'
      LMAX=0
      LMA2=0
      LMIN=0
      LMI2=0
      NSMALL = 0
      NTRACE = 0

      DO 7 L=NDMIN,ND
      DO 7 J = 1, N
      POP1LJ=POP1(L,J)
C***  CHECK IF CURRENT LEVEL BELONGS TO TRACE ELEMENTS => not considered
      IF (TRACEELEM(NOM(J))) THEN
        NTRACE = NTRACE + 1
        CYCLE
      ENDIF
      
C***  TOO SMALL OR NEGATIVE POPNUMBERS ARE CONSIDERED AS BEEING = SMALLPOP  **
      IF (POP1LJ .LT. SMALLPOP) POP1LJ = SMALLPOP
      POPLJ=POPNUM(L,J)
C***  POPNUMBER LESS THAN SMALLPOP ARE NOT ACCOUNTED FOR THE CALCULATION
C***    OF THE MAXMIMUM CORRECTION
      IF (POPLJ .LT. SMALLPOP) THEN
        POPLJ = SMALLPOP
        NSMALL = NSMALL + 1
      ELSE         
        SELECTCASE (ICMMODE)
          CASE (0)
            bConsiderDepth = .TRUE.                    
          CASE (1)
            bConsiderDepth = (CKONVER(L) == 'C')
          CASE DEFAULT
            WRITE (hCPR,*) 'Invalid CORMAX-MODE, please change CARDS'
            STOP 'FATAL ERROR IN PRICORR'
        ENDSELECT
        IF (bConsiderDepth) THEN
         Q=POPLJ/POP1LJ
         IF (Q .GT. QMA2) THEN
           IF (Q .GT. QMAX) THEN
              QMA2=QMAX
              QMAX=Q
              LEVMA2=LEVMAX
              LEVMAX=LEVEL(J)
              LMA2=LMAX
              LMAX=L
              JMAX = J
           ELSE
              QMA2=Q
              LEVMA2=LEVEL(J)
              LMA2=L
              JMA2 = J
           ENDIF
         ENDIF
         IF (Q .LT. QMI2) THEN
           IF (Q .LT. QMIN) THEN
              QMI2=QMIN
              QMIN=Q
              LEVMI2=LEVMIN
              LEVMIN=LEVEL(J)
              LMI2=LMIN
              LMIN=L
              JMIN = J
           ELSE
              QMI2=Q
              LEVMI2=LEVEL(J)
              LMI2=L
              JMI2 = J
           ENDIF
         ENDIF
        ENDIF
      ENDIF
    7 CONTINUE

C***********************************************************************
      CORMAX=AMAX1(QMAX-1.,1./QMIN-1.)
C***********************************************************************

      IF (CORMAX .GT. .0) CORLOG=ALOG10(CORMAX)
      PRINT 9,QMAX,LEVMAX,LMAX,QMA2,LEVMA2,LMA2,QMIN,LEVMIN,LMIN,QMI2,
     $        LEVMI2,LMI2,NDMIN,ND,CORMAX,CORLOG
    9 FORMAT (/,15X,'MAX:',F8.4,'  (',A10,'  L=',I3,')',
     $           7X,'2ND:',F8.4,'  (',A10,'  L=',I3,')',
     $        /,15X,'MIN:',F8.4,'  (',A10,'  L=',I3,')',
     $           7X,'2ND:',F8.4,'  (',A10,'  L=',I3,')',
     $        /,15X,'DEPTH POINTS CONSIDERED:',I5,'  TO',I5,
     $          15X,'CORMAX=',F7.4,5X,'LOG=',F6.2,//)

      WRITE (0, '(A,I7,A,F6.2,A)') 'JOBNUM=', JOBNUM, 
     >       '   log max correction =', CORLOG, ' <<<<<<<'
 
C***  MAXIMUM TEMPERATURE CORRECTION
      IF (BUNLU) THEN
        RTCMAX=.0
        LRTCMAX=0
        DO L=1,ND
          RTC=(TNEW(L)-TOLD(L))/TOLD(L)
          IF (ABS(RTC) .GT. ABS(RTCMAX)) THEN
            RTCMAX=RTC
            LRTCMAX=L
          ENDIF
        ENDDO

        PRINT 35, RTCMAX, LRTCMAX, 
     >            DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB
   35   FORMAT (15X,'MAX.REL.TEMP.CORR.=',F8.4,'  AT L=',I3,
     >     '    (UNLU SCALINGS LOC, INT, RMAX, TB=', 4(F5.3,1X), ')')

        IF (DUNLU_TB .GT. .0) PRINT 36, TBTAU 
   36    FORMAT 
     >   (15X, 'Thermal Balance term suppressed for tau_Ross >', G12.5)

        IF (TAUINT > .0) PRINT 37, TAUINT
   37    FORMAT 
     >   (15X, 'INT and RMAX term damped depth-dependent, TAUINT=', F8.2)

        PRINT *     

      ENDIF
         
      IF (BTDIFFUS) 
     >   WRITE (*,'(15X,2A,E10.3,/)')
     >    'TEMPERATURE AT INNER BOUNDARY ADJUSTED ',
     >    'BY DIFFUSION APPROXIMATION : REL.CORR.= ', TNDCORR
 
      IF (.NOT. BTDIFFUS .AND. BUNLU) 
     >   WRITE (*, '(15X,2A,E10.3,/)') 
     >    'TEMPERATURE GRADIENT AT INNER BOUNDARY CORRECTED ',
     >    'FOR FLUX CONSISTENCY : REL.CORR.= ', HNDCORFAC-1.
     
      IF (GAMMAC.GT..0 .OR. GAMMAL.GT..0)
     $   PRINT 14, GAMMAC, GAMMAL, DELTAC, GAMMAR, GAMMAD
   14    FORMAT (10X,'SCHARMER-PARAMETERS:',  10X,
     $          'GAMMAC=',F9.1,'   GAMMAL=',F9.1,/,
     $     40X,'DELTAC=',F9.1,'   GAMMAR=',F9.1,'   GAMMAD=',F9.1,//)
 
      IF (REDUCE .NE. 1.) PRINT 10,REDUCE
   10 FORMAT (10X,'CORRECTIONS REDUCED BY FACTOR',F5.2,//)
 

C**********************************************************************
C***  PRINTOUT OF WARNINGS AND ERROR MESSAGES  ************************
C**********************************************************************
      IF (NSMALL .GT. 0) PRINT 27, NSMALL, SMALLPOP
   27 FORMAT ( 8X, I6, 
     >  ' WARNINGS: POP. NUMBERS .LT.', 1PE8.1, 
     >  ' NOT ACCOUNTED FOR IN THE MAX. CORRECTIONS')
      IF (NTRACE > 0) THEN
        WRITE (hOUT, FMT='(A,I6,A)') 'NOTE: ', NTRACE, 
     >   ' TRACE LEVELS NOT ACCOUNTED FOR IN THE MAX. CORRECTIONS'
      ENDIF

C***  PRINTOUT OF WARNINGS FROM SUBROUTINE LINPOP  *********************

      IF (NEGINTL .GT. 0)  
     $   PRINT 28, NEGINTL,INEGMIN, INEGMAX, LNEGMIN, LNEGMAX
   28 FORMAT ( 8X,I6,' WARNINGS: NEGATIVE LINE INTENSITIES REPLACED',
     $   ' BY ZERO  --  ', I4, ' < IND <', I4, ';', I5, ' < L <', I3)

      IF (NEGINTC .GT. 0)  
     $   PRINT 29, NEGINTC,INEGMINC, INEGMAXC, LNEGMINC, LNEGMAXC
   29 FORMAT ( 8X,I6,' WARNINGS: NEGATIVE CONT INTENSITIES REPLACED',
     $   ' BY ZERO  --  ', I4, ' <  K  <', I4, ';', I5, ' < L <', I3)

      IF (ITWARN .GT. 0) PRINT 33, ITWARN, ITMAX
   33 FORMAT ( 8X,I6,' WARNINGS: NEWTON-RAPHSON ITERATION ',
     $      'NOT CONVERGED (MAX. NUMBER OF ITERATIONS: ',I3,')')
 

C***  PRINTOUT OF WARNINGS FROM SUBROUTINE LINPOP  *********************

      IF (IWARN_NEG_XJCAPP .GT. 0) THEN
        write (0,*) 'Neg XJCAPP in XJAPP:', IWARN_NEG_XJCAPP
        PRINT 40, IWARN_NEG_XJCAPP
   40   FORMAT ( 8X,I6,' WARNINGS: Negative XJCAPP appeared in XJAPP')
      ENDIF

      IF (IWARN_NEG_XJLAPP .GT. 0) THEN
        write (0,*) 'Neg XJLAPP in XJAPP:', IWARN_NEG_XJLAPP
        PRINT 41, IWARN_NEG_XJLAPP
   41   FORMAT ( 8X,I6,' WARNINGS: Negative XJLAPP appeared in XJAPP')
      ENDIF

C***  PRINTOUT OF ERROR MESSAGES FROM SUBROUTINE GAUNTFF

      IF (CORMAX .LT. EPSILON) THEN
   26 FORMAT (7X,I7,' WARNINGS: CALLS OF GAUNTFF BEYOND ',A,' BOUND')
        IF (NTUP .GT.0) PRINT 26, NTUP , 'UPPER TEMPERATURE'
        IF (NTLOW.GT.0) PRINT 26, NTLOW, 'LOWER TEMPERATURE'
        IF (NFUP .GT.0) PRINT 26, NFUP , 'UPPER FREQUENCY'
        IF (NFLOW.GT.0) PRINT 26, NFLOW, 'LOWER FREQUENCY'
        ENDIF
 
C***  Update GAMMA HISTORY
      IF (.NOT. STHLP) THEN
        GAHIST(1,1)  = FLOAT(JOBNUM)
        GAHIST(3,1)  = GAMMAC
        GAHIST(4,1)  = GAMMAL
        GAHIST(5,1)  = GAMMAR
        GAHIST(6,1)  = GAMMAD
        GAHIST(7,1)  = FLOAT(LMAX)
        GAHIST(8,1)  = FLOAT(JMAX)
        GAHIST(9,1)  = QMAX
        GAHIST(10,1) = FLOAT(LMA2)
        GAHIST(11,1) = FLOAT(JMA2)
        GAHIST(12,1) = QMA2
        GAHIST(13,1) = FLOAT(LMIN)
        GAHIST(14,1) = FLOAT(JMIN)
        GAHIST(15,1) = QMIN
        GAHIST(16,1) = FLOAT(LMI2)
        GAHIST(17,1) = FLOAT(JMI2)
        GAHIST(18,1) = QMI2
        GAHIST(19,1) = CORMAX
        GAHIST(20,1) = LRTMAX
        GAHIST(21,1) = RTCMAX
      ENDIF

      RETURN
      END
      SUBROUTINE PRIEX (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $                  TNEW,NOTEMP)
C***********************************************************************
C***  OUTPUT OF EXTRAPOLATED POPNUMBERS (AND TEMPERATURE)
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, N, JOBNUM
      REAL, INTENT(IN) :: LSPOP
      REAL, DIMENSION(ND), INTENT(IN) :: RNE, TNEW
      REAL, DIMENSION(ND,N), INTENT(IN) :: POPNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      CHARACTER(10), DIMENSION(N), INTENT(IN) :: LEVEL
      LOGICAL, INTENT(IN) :: NOTEMP

      CHARACTER(12), DIMENSION(10) :: VALUE
      INTEGER :: L, J, J1, J2

      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A,  20X,'JOB NO.',I7,
     $ //,20X,'RELATIVE NON-LTE POPULATION NUMBERS EXTRAPOLATED ',
     $ 'FROM THE LAST THREE ITERATIONS',/,20X,79('-'))
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2, (LEVEL(J),J=J1,J2)
    2 FORMAT (//,'  L EL.DENS.',10(2X,A10))
      IF (N.LT.J2) PRINT 11
   11 FORMAT (1X)
      PRINT 11
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3

      DO 93 J=J1, J2
      IF (POPNUM(L,J) .GT. .0) THEN
         WRITE (VALUE(J),'(F12.2)') ALOG10(POPNUM(L,J))
         ELSE IF (POPNUM(L,J) .LT. .0) THEN
         VALUE(J) = '* NEGATIVE *'
         ELSE
         VALUE(J) = '*** ZERO ***'
         ENDIF

   93 CONTINUE

      PRINT 9, L,RNE(L), (VALUE(J), J=J1,J2)
    9 FORMAT (I3,F7.3,2X,10A12)
    3 CONTINUE
      IF (J2.EQ.N) GOTO 99
      J1=J1+10
      GOTO 4

C***  PRINTOUT OF THE TEMPERATURE STRATIFICATION
   99 CONTINUE
      IF (.NOT. NOTEMP) THEN
      PRINT 10
   10 FORMAT (///,40X,'TEMPERATURE STRATIFICATION',/,40X,26('='),
     $   //,40X,'DEPTH INDEX      T (KELVIN)',/)
      DO 12 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 12
      PRINT 13, L,TNEW(L)
   13 FORMAT (40X,I10,F15.0)
   12 CONTINUE
      ENDIF

      RETURN
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
