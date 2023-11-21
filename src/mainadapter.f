      PROGRAM MAINadapter 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:07 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL adapter 
      END
      SUBROUTINE ADAPOP (POPNUM, ND, N, POPOLD, NDOLD, NOLD, NCHARG,
     >          NATOM, ABXYZ, NFIRST, NLAST, RNE, NTRANS, POPLTE, 
     >          BDEPART, ADPWEIGHT, RADIUS, ROLD, POPHELP, TAURCONT, 
     >          TAURCONTOLD, POPLTE_OLD, ENTOT, ENTOTOLD, 
     >          BTAUR, bUseENTOT, POPMIN) 
C***********************************************************************
C***  TRANSFORMATION OF POPULATION NUMBERS FROM OLD TO NEW MODEL ATOM
C***  Radically simplified version: wrh 10-Aug-2007
C***  The assignment of levels is directed by the vector NTRANS(J)
C***    IF NTRANS(J) = -1 : no action, POPNUM stays from WRSTART 
C***    IF NTRANS(J) =  0 : POPNUM set to ZERO 
C***    IF NTRANS(J) >  0 : assign POPOLD with index NTRANS
C***********************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
 
      INTEGER, INTENT(IN) :: N, NOLD, ND, NDOLD, NATOM
      REAL, INTENT(IN) :: POPMIN

      INTEGER, DIMENSION(N) :: NCHARG, NTRANS
      REAL, DIMENSION(N) :: ADPWEIGHT(N)      

      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST
      REAL, DIMENSION(NATOM) :: ABXYZ
      
      REAL, DIMENSION(ND) :: RNE, RADIUS, TAURCONT, ENTOT
      REAL, DIMENSION(NDOLD) :: ROLD, TAURCONTOLD, ENTOTOLD, 
     >                          ENTOTOLDLOG, POPHELPJLOG

      REAL, DIMENSION(ND, N) :: POPNUM, POPLTE
      REAL, DIMENSION(ND, NOLD) :: POPOLD
      REAL, DIMENSION(NDOLD,NOLD) :: POPHELP, POPLTE_OLD

      INTEGER :: L, J, NFIRNA, NLANA, NA
      REAL :: SUM, POPJLOGL, ENTOTLOGL
      
      LOGICAL :: BDEPART, BTAUR, bUseENTOT
      
C***  The old popnumbers in POPHELP are interpolated with respect
C***  to the depth coordinate and stored in POPOLD (which then has
C***  still the old atomic levels, but the new radius grid)

C***  Using departure coefficients: replace old POPs by DEPARTs
      IF (BDEPART) THEN
         WRITE (0,*) 
     >      'Old DEPARTure coeficients used instead of POPNUMbers'
         DO L=1, NDOLD
            DO J=1, NOLD
               POPHELP(L,J) = POPHELP(L,J) / POPLTE_OLD(L,J)
            ENDDO
         ENDDO
         DO L=1, ND
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) / POPLTE(L,J)
            ENDDO
         ENDDO
      ENDIF


C***  Interpolation on Tau-Grid 
      IF (BTAUR) THEN
         WRITE (0,*) 'Interpolation of Popnumbers on Tau-Grid'
         DO L=1, ND
           DO J=1, NOLD
              IF (TAURCONT(L) .GT. TAURCONTOLD(NDOLD)) THEN
                 POPOLD(L,J) = POPHELP(NDOLD,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), TAURCONT(L), 
     >                    POPHELP(1,J), TAURCONTOLD, NDOLD)
              ENDIF 
           ENDDO
         ENDDO
      
      ELSEIF (bUseENTOT) THEN
        WRITE (0,*) 'Interpolation of Popnumbers over LOG density'
        DO J=1, NOLD
          !prepare necessary vectors
          DO L=1, NDOLD
            ENTOTOLDLOG(L) = LOG10(ENTOTOLD(L))
            POPHELPJLOG(L) = LOG10(MAX(POPHELP(L,J), POPMIN))
          ENDDO          
          !Perform interpolation on log(n_tot)
          dploop: DO L=1, ND
            ENTOTLOGL = LOG10(ENTOT(L))
            IF (ENTOTLOGL > ENTOTOLDLOG(NDOLD)) THEN
              !more dense than old innermost value => take old inner boundary value
              POPJLOGL = POPHELPJLOG(NDOLD)
            ELSEIF (ENTOTLOGL < ENTOTOLDLOG(1)) THEN
              !less dense than old outermost value => take old outer boundary value
              POPJLOGL = POPHELPJLOG(1)
            ELSE
              CALL SPLINPOX(POPJLOGL, ENTOTLOGL,
     >                     POPHELPJLOG, ENTOTOLDLOG, NDOLD)
            ENDIF
            POPOLD(L,J) = 10**(POPJLOGL)
          ENDDO dploop
        ENDDO              
      
      ELSE

C***  INTERPOLATION OF OLD POPNUMBERS TO THE NEW RADIUS GRID
        WRITE (0,*) 'Interpolation of Popnumbers on Radius-Grid'
        DO L=1, ND
           DO J=1, NOLD
              IF (RADIUS(L) .GT. ROLD(1)) THEN
                 POPOLD(L,J) = POPHELP(1,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), RADIUS(L), 
     >              POPHELP(1,J), ROLD, NDOLD)
              ENDIF
           ENDDO
         ENDDO

      ENDIF

C*****************************************************************
C***  Now the replacement of levels according to NTRANS
C*****************************************************************

C***  Loop over all depth points --------------------------------
      DO L=1, ND

C***  Copy old POPNUMs as assigned 
         DO J=1, N
            IF (NTRANS(J) .EQ. 0) THEN
               POPNUM(L,J) = POPMIN
            ELSE IF (NTRANS(J) .GT. 0 .AND. NTRANS(J) .LE. NOLD) THEN
               POPNUM(L,J) = POPOLD(L,NTRANS(J)) * ADPWEIGHT(J)
            ENDIF
         ENDDO

C***  If POPNUM actually contains DEPARTURE coeficients,
C***   convert them now back to POPNUMs
         IF (BDEPART) THEN
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) * POPLTE(L,J)
            ENDDO
         ENDIF

C***  Renormalization to the abundance of each element
 
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
         DO NA=1, NATOM
            SUM=0.0
            NFIRNA = NFIRST(NA)
            NLANA = NLAST(NA)
            DO J = NFIRNA, NLANA
               SUM = SUM + POPNUM(L,J)
            ENDDO
            SUM = SUM / ABXYZ(NA)
            IF (SUM /= 0.) THEN
               DO J=NFIRNA,NLANA
                  !POPMIN ensurance test (ansander, 2014)
                  IF (POPNUM(L,J) > POPMIN) THEN
                    POPNUM(L,J) = POPNUM(L,J) / SUM
                  ELSE
                    POPNUM(L,J) = POPMIN
                  ENDIF                   
               ENDDO
            ENDIF
         ENDDO
 
C***  Consistent electron density
         RNE(L) = .0
         DO J=1, N
           RNE(L) = RNE(L) + NCHARG(J) * POPNUM(L,J)
         ENDDO
 
      ENDDO   ! Deph points -------------------------------------! 

      RETURN
      END
C***  MAIN PROGRAM ADAPTER  ****************************************************
      SUBROUTINE ADAPTER
C***********************************************************************
C***  THIS PROGRAM IS AN ADAPTER BETWEEN MODEL FILES BASED ON DIFFERENT
C***  NUMBERS OF ATOMIC LEVELS
C***  IT TRANSFORMS POPULATION NUMBERS FROM OLD TO NEW ENERGY LEVELS AS
C***  SPECIFIED IN THE INPUT OPTIONS, THEREBY FILLING UP THE ADDITIONAL
C***  LEVELS WITH POP.NUMBERS CALCULATED FROM BOLTZMANN'S OR SAHA'S EQUATION
C***********************************************************************
      IMPLICIT NONE
 
C***  DEFINE ARRAY DIMENSIONS ******************************************
      INTEGER, PARAMETER :: MAXATOM =          26 
      INTEGER, PARAMETER :: NDIM    =        2560 
      INTEGER, PARAMETER :: NFDIM   = 2*NDIM + 400
      INTEGER, PARAMETER :: MAXKONT =     NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR =        NDIM 
      INTEGER, PARAMETER :: MAXIND  =       40000 
      INTEGER, PARAMETER :: MAXFEIND  =      2500 
      INTEGER, PARAMETER :: NDDIM   =          89 
      INTEGER, PARAMETER :: MAXHIST =        4000 
      INTEGER, PARAMETER :: MAXLEVELCARD =   1000       

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 
      
C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 3200 
      COMMON / COMAUTO / LOWAUTO(MAXAUTO),WAUTO(MAXAUTO)
     $                  ,EAUTO(MAXAUTO),AAUTO(MAXAUTO),IONAUTO(MAXAUTO)
     $                  ,KRUDAUT(MAXAUTO)

      INTEGER, DIMENSION(NDIM) :: NCHARG, IONGRND, MAINQN, NOM, NTRANS
      INTEGER, DIMENSION(MAXIND) :: INDLOW, INDNUP
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ADPWEIGHT
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO, 
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(4, NDIM) :: ALTESUM 
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      REAL, DIMENSION(NDDIM) :: RNE, RADIUS, ROLD
      REAL, DIMENSION(NDDIM,NDIM) :: POPNUM, POPOLD, POPHELP
      REAL, DIMENSION(MAXATOM,MAXION) :: EDGEK, SIGMATHK, SEXPOK

      CHARACTER(MAXHIST*8) :: MODHIST

      CHARACTER(64) :: BUFFER64
      CHARACTER(100) :: MODHEAD, OLDHEAD
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(4) :: KEYCBB(MAXIND)
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      LOGICAL :: OLDSTART, NEWATOM, BDEPART, BTAUR, bUseENTOT
      CHARACTER(80), DIMENSION(MAXLEVELCARD) :: LEVELCARD
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

      REAL, DIMENSION(NDDIM) :: TAURCONT, TAURCONTOLD, ENTOT, ENTOTOLD
      REAL, DIMENSION(NDDIM,NDIM) :: POPLTE, POPLTE_OLD

      REAL :: WAUTO, EAUTO, AAUTO, VDOPFE, DXFE, XLAM0FE, POPMIN

      INTEGER :: N, ND, NDOLD, NOLD, LAST, IERR, IDUMMY, LASTIND,
     >           NATOM, NAUTO, LOWAUTO, IONAUTO, LASTKON, KRUDAUT,
     >           LASTFE, NLEVELCARD, N_WITH_DRLEVELS
  
      INTEGER, EXTERNAL :: IDX

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl / hydro

      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      COMMON /IRON/ FEDUMMY, INDRB, INDRF, SIGMAFE, IFRBSTA, IFRBEND, 
     >              IFENUP, IFELOW, SIGMAINT
      LOGICAL :: BFEMODEL

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !write to MODEL file      
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
C***  Write Link Data (Program Version) to CPR file
      WRITE (hCPR,'(2A)') '>>> ADAPTER started: Program Version from '
     >                 ,LINK_DATE
      WRITE (hCPR,'(4A)') '>>> created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL

      IF (OPSYS == 'CRAY') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

C***  NEW ENERGY LEVELS FROM NEW MODEL ATOM: 
      CALL DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
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
     >            'ADAPTER', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)

C***  DECODING INPUT OPTIONS
      CALL DECADP (OLDSTART, BDEPART, BTAUR, bUseENTOT, POPMIN, 
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)

C***  NO OPTION "OLDSTART" DECODED: ADAPTER IS NOT NECESSARY
      IF (.NOT. OLDSTART) GOTO 999
 
C***  READ OLD AND NEW MODEL FILES
      CALL RMODADP (NDDIM, OLDHEAD, N, NOLD, NDIM, 
     >           NATOM, MODHEAD, ND, NDOLD, ABXYZ, LAST, MODHIST,
     >           RADIUS, ROLD, POPHELP, POPNUM, TAURCONT, TAURCONTOLD, 
     >           POPLTE, POPLTE_OLD, ENTOT, ENTOTOLD, BTAUR)

C***  Decoding the LEVEL cards
      CALL ADATRANS (NTRANS, ADPWEIGHT, N, NOLD, NEWATOM, POPMIN,
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)

C***  TRANSFORMING THE POP.NUMBERS  ************************************
      CALL ADAPOP (POPNUM, ND, N, POPOLD, NDOLD, NOLD, NCHARG,
     >          NATOM, ABXYZ, NFIRST, NLAST, RNE, NTRANS, POPLTE, 
     >          BDEPART, ADPWEIGHT, RADIUS, ROLD, POPHELP, TAURCONT,
     >          TAURCONTOLD, POPLTE_OLD, ENTOT, ENTOTOLD,
     >          BTAUR, bUseENTOT, POPMIN)

C***  UPDATING THE MODEL HISTORY  **************************************
      WRITE (UNIT=BUFFER64, FMT=60) OLDHEAD(15:32)
   60 FORMAT ('/     1A. ADAPTER:  POPNUMBERS FROM OLD MODEL ',A18)
      CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,64,BUFFER64)

      !write model history entry into explicit history file
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(BUFFER64))
      CLOSE(hHIST)


C***  PRINTOUT  *******************************************************
      CALL PRIADP (MODHEAD, OLDHEAD, NEWATOM, LEVEL, NTRANS, 
     >   ADPWEIGHT, N, NOLD)

C***  UPDATING THE MODEL FILE  *****************************************
      CALL WRITMS (3,POPNUM,ND*N,'POPNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (3,RNE,ND,'RNE     ',-1, IDUMMY, IERR)
      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      CALL CLOSMS (3, IERR)
 

  999 CALL JSYMSET ('G0','0')
 
      CALL STAMP (OPSYS, 'ADAPTER', TIM1)

      STOP 'O.K.'
      END
      SUBROUTINE ADATRANS (NTRANS, ADPWEIGHT, N, NOLD, NEWATOM, POPMIN,
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)

C*****************************************************************
C***  Interprets the LEVEL-Cards for the ADAPTER and writes 
C***  the vector NTRANS 
C***     index: level index in new model atom
C***     value: -1 --> level is not overwritten 
C***             0 --> level is set to POPMIN
C***            >0 --> corresponding index in old model atom 
C***
C***  LEVEL options modernized by wrh 29-Jul-2007 
C***  This subroutine was separated from DECADP, because only 
C***  after reading both MODEL files, NOLD becomes available
C***      wrh, 27-Nov-2010 
C*****************************************************************

      DIMENSION NTRANS(N), ADPWEIGHT(N)
      LOGICAL NEWATOM, CHRINSTR
      CHARACTER KARTE*80, REST*80, ACTPAR1*20, ACTPAR2*20
      CHARACTER LEVELCARD*80(MAXLEVELCARD)

C***  DEFAULT VALUES: 
      NEWATOM=.FALSE.
      DO I=1,N
        ADPWEIGHT(I) = 1.
        NTRANS(I)    = -1
      ENDDO

      NEWATOM = NLEVELCARD .GT. 0

C***  Simple OLDSTART  (NO "LEVEL" CARDS GIVEN)
      IF (.NOT. NEWATOM) THEN
         DO J=1, N
           NTRANS(J) = J
         ENDDO
C***    ...  but model atom became bigger? --> Set new levels ZERO
         IF (N .GT. NOLD) THEN
            NEWATOM = .TRUE.
            DO J=NOLD+1, N
               NTRANS(J) = 0
            ENDDO
            WRITE (0,'(A,1PG8.1)') 
     >            'NOTE (ADAPTER): NEW MODEL HAS MORE ' //
     >            'LEVELS THAN OLD --> SET to POPMIN:',POPMIN  
         ENDIF
      ENDIF


C***  Oldstart with LEVEL cards

C***  Loop over all LEVEL cards
      DO ILC =1, NLEVELCARD
         KARTE = LEVELCARD(ILC)
         REST = KARTE(6:)       
         CALL SARGC (REST, NPAR)
         IF (NPAR .LE. 0) GOTO 90        
         CALL SARGV (REST,1,ACTPAR1)
         LAST1 = IDX (ACTPAR1)

C***     If there are no more parameters, skip the next part!
         IF (NPAR .GE.2) THEN       
            CALL SARGREST (REST, NPAR, 2, IA, IE)
            REST = REST(IA:)
            CALL SARGV (REST,1,ACTPAR2)

C***        If last character of ACTPAR1 is "-" 
C***         -->  directly append next parameter 
            IF (ACTPAR1(LAST1:LAST1) .EQ. '-') THEN
               ACTPAR1(LAST1+1:) = ACTPAR2
               CALL SARGREST (REST, NPAR, 2, IA, IE)
               REST = REST(IA:)
            ENDIF

C***        If first character of next parameter is "-"
            IF (ACTPAR2(1:1) .EQ. '-') THEN
C***        ... and this is the whole parameter --> glue next parameter
               IF (IDX(ACTPAR2) .EQ. 1) THEN
                  CALL SARGV (REST,2,ACTPAR2)
                  ACTPAR1 = ACTPAR1(:LAST1) // '-' // ACTPAR2
                  CALL SARGC (REST, NPAR)
                  IF (NPAR .GE. 3) THEN 
                     CALL SARGREST (REST, NPAR, 3, IA, IE)
                     REST = REST(IA:)
                  ELSE
                     REST = ' '
                  ENDIF
               ELSE
                  ACTPAR1 = ACTPAR1(:LAST1) // ACTPAR2
                  IF (NPAR .GE. 2) THEN 
                     CALL SARGREST (REST, NPAR, 2, IA, IE)
                     REST = REST(IA:)
                  ELSE
                     REST = ' '
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            REST = ' '
         ENDIF

C***     Now ACTPAR1 should contain the compact string nnn-mmm
C***        --> decode nnn and mmm

         IF (CHRINSTR ('-', ACTPAR1)) THEN
            LAST1 = IDX(ACTPAR1)
            DO I=1, LAST1
               IF (ACTPAR1(I:I) .EQ. '-') THEN
                  MPOS = I
                  EXIT
               ENDIF
            ENDDO
            READ (ACTPAR1(:MPOS-1), '(I10)', ERR=91) NSTART
            READ (ACTPAR1(MPOS+1:), '(I10)', ERR=91) NSTOP
         ELSE
            READ (ACTPAR1, '(I10)', ERR=91) NSTART
            NSTOP = NSTART
         ENDIF

         NSTOPOLD = MIN(NSTOP, NOLD)  !this value is needed for SHIFT cards
         NSTOP = MIN0 (NSTOP, N)

         IF (NSTART .LE. 0) GOTO 96   ! Error exit

C***     Decode the further parameters

C***     No further parameters -- identify the levels
         IF (REST .EQ. ' ') THEN
            DO I = NSTART, NSTOP 
               IF (I .GT. NOLD) THEN
                  WRITE (0,'(A)') 'NOTE FROM ADAPTER:ADATRANS:'
                  WRITE (0,*) KARTE(:IDX(KARTE))
                  WRITE (0,'(A,I5)') ' LOOP TRUNCATED AT', I-1 
                  EXIT
               ENDIF
               NTRANS(I)=I
            ENDDO

         ELSE IF (REST(1:5) .EQ. 'SHIFT') THEN
C                                 =====
            REST = REST(6:)
            IF (REST .EQ. ' ') GOTO 92
            CALL SARGV (REST, 1, ACTPAR1)
            READ (ACTPAR1, '(I10)', ERR=92) ISHIFT
            NSTOPMAX = MIN0(NOLD, N-ISHIFT)
            NSTOP = MIN0(NSTOPOLD,NSTOPMAX)  ! corrected 1-Mar-2017 ansander
            DO I = NSTART, NSTOP 
               IF (I+ISHIFT .GT. N .OR. I .GT. NOLD) THEN
                  WRITE (0,'(A)') 'NOTE FROM ADAPTER:ADATRANS:'
                  WRITE (0,*) KARTE(:IDX(KARTE))
                  WRITE (0,'(A,I5)') ' LOOP TRUNCATED AT', I-1 
                  EXIT
               ENDIF
               NTRANS(I+ISHIFT)=I
            ENDDO


         ELSE IF (REST(1:4) .EQ. 'NULL') THEN
C                                 =====
C***       SET NEW LEVEL TO ZERO
            DO I = NSTART, NSTOP 
               NTRANS(I) = 0
            ENDDO

         ELSE IF (REST(1:6) .EQ. 'WEIGHT') THEN
            CALL SARGC (REST, NPAR)
            IF (NPAR .LT. 2) GOTO 94
            CALL SARGV (REST, 2, ACTPAR1)
            READ (ACTPAR1, '(F10.0)', ERR=94) ADPW
            DO I = NSTART, NSTOP 
               ADPWEIGHT(I) = ADPW
            ENDDO

         ELSE
C***        the only remaining possibility: next parameter is old level
            IF (NSTART .EQ. NSTOP) THEN
               CALL SARGV (REST, 1, ACTPAR1)
               READ (ACTPAR1, '(I10)', ERR=91) NOLDLEVEL
               NTRANS(NSTART) = NOLDLEVEL
            ELSE
C***           Unidentified further parameter
               GOTO 93
            ENDIF


         ENDIF
      ENDDO !--------- end of loop over all LEVEL cards


C***  Check assigned oldlevel for existence
      DO J=1, N
         IF (NTRANS(J) .GT. NOLD) THEN
            WRITE (0, '(A,/,A, I4, A, I4, A)')
     >          'ADAPTER: WARNING OF INTERNAL INCONSISTENCY',
     >          'non-existing old level:', NTRANS(J),
     >          ' asigned to new level:', J, ' - replaced by NULL'
            NTRANS(J) = 0
         ENDIF
      ENDDO

      RETURN

C***  ERROR EXITS *******************************

   90 WRITE (0,*) 'LEVEL option needs parameters!'
      GOTO 95

   91 WRITE (0,*) 'Error when decoding LEVEL index'
      GOTO 95 

   92 WRITE (0,*) 'Error when decoding SHIFT index'
      GOTO 95

   93 WRITE (0,*) 'Unidentified parameter:', REST(:IDX(REST))
      GOTO 95

   94 WRITE (0,*) 'A number must be given for the WEIGHT'
      GOTO 95

   96 WRITE (0,*) 'Invalid level index (.LE. 0) encountered'
      GOTO 95

   95 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) KARTE 
      STOP 'FATAL ERROR IN ADAPTER:ADPTRANS'

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
      LOGICAL FUNCTION CHRINSTR (TESTCHAR, STRING)
C********************************************************************
C***  Returns .TRUE. if the character TESTCHAR is member of STRING 
C********************************************************************
      CHARACTER TESTCHAR*1, STRING*(*)

      CHRINSTR = .FALSE.
      L = IDX(STRING)

      DO I=1, L
      CHRINSTR = CHRINSTR .OR. (TESTCHAR .EQ. STRING(I:I))
      ENDDO

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
      SUBROUTINE DECADP (OLDSTART, BDEPART, BTAUR, bUseENTOT, POPMIN, 
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)
C**********************************************************************
C***  DEECODE INPUT OPTIONS FOR PROGRAM 'ADAPTER'
C**********************************************************************

      LOGICAL OLDSTART, BDEPART, BTAUR, bUseENTOT
      CHARACTER KARTE*80, ACTPAR1*20
      CHARACTER LEVELCARD(MAXLEVELCARD)*80

C***  DEFAULT VALUES: 
      OLDSTART=.FALSE.
      BDEPART = .FALSE. 
      BTAUR   = .FALSE.
      bUseENTOT = .TRUE.
C***  POPMIN: NULL popnums are set to this value
      POPMIN = 1.E-25
      NLEVELCARD = 0

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1

    5 READ (1,9,END=99) KARTE
    9 FORMAT (A)

C***  LEVEL option, modernized by wrh 29-Jul-2007 

      IF (KARTE(:5) .EQ. 'LEVEL') THEN
C                         =====
         NLEVELCARD = NLEVELCARD + 1
         IF (NLEVELCARD .GT. MAXLEVELCARD) GOTO 93
         LEVELCARD(NLEVELCARD) = KARTE

      ELSE IF (KARTE(:6) .EQ. 'POPMIN' ) THEN
C                              ======
         CALL SARGV (KARTE, 2, ACTPAR1)
         READ (ACTPAR1,'(F10.0)', ERR=95) POPMIN
 
      ELSE IF (KARTE(:8) .EQ. 'OLDSTART') THEN
C                              ========
         OLDSTART=.TRUE.
         CALL SARGC (KARTE, NPAR)
         DO I=2, NPAR         
            CALL SARGV (KARTE, I, ACTPAR1)
            IF (ACTPAR1(:6) .EQ. 'DEPART') THEN
C                                 ======
               BDEPART = .TRUE.
            ELSE IF (ACTPAR1(:3) .EQ. 'TAU') THEN
C                                      ===
               BTAUR   = .TRUE.
            ELSE IF (ACTPAR1(:6) == 'RADIUS') THEN
C                                    ======
               bUseENTOT = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      GOTO 5

C***  End of CARDS reached
   99 CLOSE (1)

      RETURN

C***  ERROR EXITS *******************************

   93 WRITE (0,*) '*** FATAL ERROR:' 
      WRITE (0,*) '*** More LEVEL cards found than dimensioned'
      WRITE (0,*) '*** MAXLEVELCARD = ', MAXLEVELCARD
      GOTO 96

   95 WRITE (0,*) '*** FATAL ERROR:' 
      WRITE (0,*) '*** Parameter not readable as floating-point number'
      GOTO 96

   96 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) KARTE(:IDX(KARTE)) 
      STOP 'FATAL ERROR IN ADAPTER:DECADP'

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

       FUNCTION IDX(TEXT)
C***  Returns the number of non-blank characters of string TEXT
       CHARACTER*(*) TEXT
       DO IDX=LEN(TEXT),1,-1
          IF(TEXT(IDX:IDX).NE.' ') RETURN
       END DO
       IDX=0
       RETURN
       END
      SUBROUTINE INSTALL
C***********************************************************************
C***  Called from Main Programs
C***********************************************************************

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C*********************************************
C***  Hier Hauptschalter:
C***  INST = 1 : Cray
C***  INST = 2 : Potsdam DEC/UNIX
C***  INST = 3 : Potsdam SGI Origin 2000
C*********************************************

C     vvvvvvvv
      INST = 2
C     ^^^^^^^^

C******  Cray  **************************************
      IF (INST .EQ. 1) THEN

      OPSYS = 'CRAY'

C******  DEC  **************************************
      ELSE IF (INST .EQ. 2) THEN

      OPSYS = 'DEC/UNIX'

C******  SGI  **************************************
      ELSE IF (INST .EQ. 3) THEN

      OPSYS = 'SGI'

C****** ERROR  **************************************
      ELSE
      STOP 'ERROR IN INSTALL'
      ENDIF

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
      SUBROUTINE LENGTHMS(ICHANNEL, NDIM, NAME, IERR)
C************************************************************
C***  ROUTINE VON wrh
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, XDUMMY, 
     >               NDIM, 'LENGTH  ', IERR)

      RETURN
      END
      SUBROUTINE LIPO (F,X,FI,XI,N)
C***********************************************************************
C***  LINEAR INTERPOLATION: FOR GIVEN X, FIND F(X) FROM A TABLE FI(XI)
C ------ INDEXSUCHE DURCH BISECTION ------
C***********************************************************************

      DIMENSION FI(N),XI(N)

      NA=1
      A=XI(1)
      NB=N
      B=XI(N)
      IF((X-A)*(X-B).GT..0) STOP 'ERROR IN SUBR. LIPO'
    1 IF((NB-NA).EQ.1) GOTO 2
      NH=(NA+NB)/2
      H=XI(NH)
      IF((X-A)*(X-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
 
    2 P=(X-A)/(B-A)
      F=P*FI(NB)+(1.-P)*FI(NA)
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
      SUBROUTINE PRIADP (MODHEAD, OLDHEAD, NEWATOM, LEVEL, NTRANS,
     >                   ADPWEIGHT, N, NOLD)
C**********************************************************************
C***  PRINTOUT FOR PROGRAM 'ADAPTER'
C**********************************************************************

      CHARACTER*100 MODHEAD, OLDHEAD
      CHARACTER*10 LEVEL(N), MESSAGE(N)*4
      DIMENSION NTRANS(N), ADPWEIGHT(N)
      LOGICAL NEWATOM

C**********************************************************************
C***  PRINT HEADER                                                  ***
C**********************************************************************
      PRINT 11, MODHEAD, OLDHEAD
   11 FORMAT (1H1,/,1X,A,20X,'JOB NO. 1A',///,50X,
     $  'A D A P T E R',/,50X,13('='),//,
     $  10X,'POPULATION NUMBERS TAKEN FROM OLD MODEL:',/,10X,A,//)

C**********************************************************************
C***  PRINTOUT OF LEVEL CORRESPONDENCE (IN 5 COLUMNS), IF CHANGED   ***
C**********************************************************************
      IF (NEWATOM) THEN
         PRINT 13
   13    FORMAT (5(6X,'OLD  NEW      LEVEL'),/)
         DO 30 K=1, N
           IF (ADPWEIGHT(K) .NE. 1.) THEN
              WRITE (MESSAGE(K), '(1X,F3.1)') ADPWEIGHT(K)
           ELSE
              MESSAGE(K) = '    '
           ENDIF
   30    CONTINUE
C***     NUMBER OF ROWS REQUIRED:
         NROW= (N-1)/5 + 1
         DO 14 IROW=1,NROW
            ILAST=MIN0(N,IROW+NROW*4)
            PRINT 15, (NTRANS(I), I, LEVEL(I),MESSAGE(I),
     >              I=IROW,ILAST,NROW)
   15       FORMAT (4X,5(1X,I4,I5,1X,A10,A4))
   14    CONTINUE
      ELSE IF (N .LT. NOLD) THEN
         WRITE (*,*) '*************************************************'
         WRITE (0,'(A,I4,A,I4)')  'WARNING: N0', N, '   NOLD= ', NOLD
         WRITE (*,'(A,I4,A,I4)')  'WARNING: N0', N, '   NOLD= ', NOLD
         WRITE (*,*) 'NEW MODEL HAS LESS LEVELS THAN OLDSTART MODEL'
         WRITE (*,*) '*************************************************'
      ENDIF

      PRINT 97
   97 FORMAT (//)


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
      SUBROUTINE RMODADP (NDDIM, OLDHEAD, N, NOLD, NDIM, 
     >             NATOM, MODHEAD, ND, NDOLD, ABXYZ,LAST,MODHIST,
     >             RADIUS,ROLD,POPHELP, POPNUM, TAURCONT, TAURCONTOLD,
     >             POPLTE, POPLTE_OLD, ENTOT, ENTOTOLD, BTAUR)
C**********************************************************************
C***  ... READS OLD AND NEW MODEL DATA FOR PROGRAM 'ADAPTER'
C***  this routine also interpolates the population numbers 
C***    for the radius grid of the new model  
C**********************************************************************

      LOGICAL BTAUR
      REAL, DIMENSION(NDDIM) :: ENTOT, ENTOTOLD


C***  READING OF THE OLD MODEL FILE  -----------------------------------
      IERR=0
      CALL OPENMS (9,IDUMMY, IDUMMY,1,IERR)
      IF ((IERR .NE. 0) .AND. (IERR .NE. -5)) THEN
         CALL REMARK ('ADAPTER: ERROR WHEN OPENING OLD MODEL FILE')
         STOP 'ERROR'
      ENDIF

      CALL READMS (9,NDOLD,1,        'ND      ', IERR)
      IF (NDOLD .GT. NDDIM) THEN
         CALL REMARK ('TOO MANY DEPTH POINTS IN OLD MODEL FILE')
         STOP 'ERROR'
         ENDIF

C***  Use MS-routine STORAGE to find out the number of levels
      CALL LENGTHMS (9, NDNOLD, 'POPNUM  ', IERR)
      NOLD = NDNOLD / NDOLD 
c      write (0,*) 'NUMBER OF LEVELS IN OLD MODEL ATOM:', NOLD
      IF (NOLD .GT. NDIM) THEN
         WRITE (0,*) 'OLD MODEL HAS MORE LEVELS THAN DIMENSIONED'
         WRITE (0,'(A,I5)') 'AVAILABLE NDIM:', NDIM 
         WRITE (0,'(A,I5)') 'REQUIRED  NDIM:', NOLD
         STOP 'ERROR in ADAPTER:RMODADP'
      ENDIF 
      CALL READMS (9,OLDHEAD,13,     'MODHEAD ', IERR)
      CALL READMS (9,POPHELP,     NDNOLD,  'POPNUM  ', IERR)
      CALL READMS (9,POPLTE_OLD , NDNOLD,  'POPLTE  ', IERR)
      CALL READMS (9,ROLD   ,     NDOLD,   'R       ', IERR)

      IF (BTAUR) THEN
         CALL READMS (9,TAURCONTOLD, NDOLD, 'TAURCONT' , IERR)
         IF (IERR .EQ. -10) THEN
            WRITE (0,*) 
     >       '*** WARNING: TAURCONT not on MODEL file, take TAUROSS'
            CALL READMS (9,TAURCONTOLD, NDOLD, 'TAUROSS ' , IERR)
            IF (IERR .EQ. -10) THEN
               WRITE (0,*) 
     >          '*** WARNING: TAUROSS not on MODEL file,' 
     >          // ' TAU option in OLDSTART disabled'
               BTAUR = .FALSE.
            ENDIF        
         ENDIF
      ENDIF

      CALL READMS (9,ENTOTOLD,    NDOLD,   'ENTOT   ', IERR)

      CALL CLOSMS (9, IERR)

C***  READING OF THE NEW MODEL FILE  -----------------------------------
      CALL OPENMS (3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (3,ND,1,'ND      ', IERR)
      IF (ND .GT. NDDIM) THEN
         CALL REMARK ('INSUFFICIENT DIMENSION: NDDIM')
         STOP 'ERROR'
         ENDIF
      NDN = ND * N
      CALL READMS (3, MODHEAD, 13   , 'MODHEAD ', IERR)
      CALL READMS (3, POPNUM , NDN  , 'POPNUM  ', IERR)
      CALL READMS (3, POPLTE , NDN  , 'POPLTE  ', IERR)
      CALL READMS (3, ABXYZ  , NATOM, 'ABXYZ   ', IERR)
      CALL READMS (3, LAST   , 1    , 'MODHIST ', IERR)
      CALL READMS (3, MODHIST, LAST , 'MODHIST ', IERR)
      CALL READMS (3, RADIUS , ND   , 'R       ', IERR)
      CALL READMS (3, TAURCONT, ND  , 'TAURCONT', IERR)
      CALL READMS (3, ENTOT  , ND   , 'ENTOT   ', IERR)


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
      SUBROUTINE SARGREST (TEXT, N, I, IFIRST, LAST)
C**********************************************************************
C***  Ermittelt den Beginn des i-ten Arguments und das nichtleere Ende 
C***  gesamten restlichen Strings. 
C***  Balancierte "..." am Anfang und Ende werden entfernt, sonst 
C***  bleiben " erhalten. Faengt des Argument mit einem Trennzeichen
C***  an (=,:/), so muss (!!!, sonst Fehlerabbruch!) der Reststring in 
C***  "..." eingeschlossen werden, ebenso natuerlich wenn der 
C***  Reststring mit einem Blank beginnen soll. An spaeterer Position sind
C***  Trennzeichen ohne Wirkung. 
C**********************************************************************

      CHARACTER TEXT*(*)

      CALL SARGP (TEXT, N, I, IFIRST, ILAST)
      LAST = IDX(TEXT)

C***  Leerer Reststring
      IF (IFIRST .EQ. -1) THEN
         IFIRST = 1
         LAST   = 1

C***  Der Reststring  beginnt mit einem Trennzeichen ohne " davor
      ELSE IF (IFIRST .EQ. 0) THEN
         WRITE (6, '(A)') 'Error when parsing the following string:'
         WRITE (6, '(A)') TEXT
         WRITE (6, '(A, I3, A)') 'Argument ', I, 
     >                        ' begins with delimiter -> use "..."'
         STOP '>>> ERROR IN SUBROUTINE SARGREST <<<'
      ENDIF

      IF (IFIRST .GT. 1) THEN
C***  Reststring begins with "
         IF (TEXT(IFIRST-1:IFIRST-1) .EQ. '"') THEN 
C***     Remove balanced closing quote if present, else restore leading quote 
            IF (TEXT(LAST:LAST) .EQ. '"' .AND. LAST .GT. IFIRST) THEN
               LAST = LAST - 1
               ELSE
               IFIRST = IFIRST - 1
            ENDIF
         ENDIF
      ENDIF

c      print *, 'Test: ', text
c      print *, 'Test: ', ifirst 
c      print *, 'Test: ', text(ifirst:ifirst)
c      print *, 'Test: ', last
c      print *, 'Test: ', text(last:last)

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
      SUBROUTINE SPLINPOX(F, X, FI, XI, N, SAFE, LPARAM, DFDX, D2FD2X)
C***********************************************************************
C***  CUBIC SPLINE INTERPOLATION, READY-FOR-USE
C***  XI(I), FI(I)  TABLE WHICH DEFINES THE FUNCTION TO BE INTERPOLATED
C***  X             ARGUMENT FOR WHICH THE FUNCTION VALUE IS REQUIRED
C***  FX            RESULTING FUNCTION VALUE
C***  THE RESULTING FUNCTION IS A PIECEWISE CUBIC INTERPOLATION 
C***     POLYNOMIAL WITH CONTINUOUS DERIVATIVE. EXTREMA CAN ONLY OCCUR 
C***     AT GIVEN MESHPOINTS (MONOTONIC VERSION AFTER M. STEFFEN)
C
C     Unified version implementing all features from classic routines
C       SPLINPO, SPLINPO_FAST and SPLINP from Goetz branch (A. Sander, Jan 2012)
C
C     Due to the optional arguments (e.g. for the derivatives), all main routines
C     need to have the following interface block:
C
C      INTERFACE SPLINPO
C        SUBROUTINE SPLINPO(F, X, FI, XI, N, DFDX, D2FD2X)
C          INTEGER, INTENT(IN) :: N          
C          REAL, DIMENSION(N), INTENT(IN) :: XI, FI
C          REAL, INTENT(OUT) :: F
C          REAL, INTENT(IN) :: X
C          LOGICAL, INTENT(IN), OPTIONAL :: SAFE
C          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
C          REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X
C        END SUBROUTINE
C      END INTERFACE SPLINPO
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      
      REAL, DIMENSION(N), INTENT(IN) :: XI, FI
      REAL, INTENT(OUT) :: F
      REAL, INTENT(IN) :: X
      LOGICAL, INTENT(IN), OPTIONAL :: SAFE
      INTEGER, INTENT(IN), OPTIONAL :: LPARAM
      REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X

      REAL :: DN, DX, DXM, FS0,
     >        D1, D2, D3, D23, H11, H12, H13, H14,
     >        H21, H22, H23, H24, H31, H32, H33, H34,
     >        H41, H42, H43, H44,
     >        F1, F2, F3, F4, FSM, FSP, S3, S4, S5,
     >        P, P1, P2, P3, P4

      INTEGER :: L, I, LA, LB

      LOGICAL :: bFoundInterval, bSafeMode

      !File and channel handles
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqXX.cpr

      SAVE      !this is required due to the _SAME_X entry point

      IF (present(SAFE)) THEN
        bSafeMode = SAFE
      ELSE
        bSafeMode = .TRUE.              !safe mode is default
      ENDIF

C***  CHECK FOR STRICTLY MONOTONIC ORDER (safe mode only)
      IF (bSafeMode) THEN
        DN = XI(N) - XI(1)
        DO L=2, N
          DX = XI(L) - XI(L-1)
          IF (DX*DN <= .0) THEN
            WRITE (hCPR, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3       FORMAT (' *** BAD USE OF SUBROUTINE SPLINPOX:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
            CALL TRBK
            STOP 'ERROR'
          ENDIF
        ENDDO

C***    FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
        bFoundInterval = .FALSE.
        DO I=2, N
          L = I
          IF ( (X-XI(L-1)) * (X-XI(L)) <= .0) THEN
            bFoundInterval = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF (.NOT. bFoundInterval) THEN
          CALL REMARK('BAD USE OF SUBR. SPLINPOX - X OUTSIDE TABLE')
          WRITE (hCPR,'(A,G12.5,5X,A,G12.5,5X,A,I3,A,G12.5)') 
     >         ' X=', X,' X(1)=', XI(1),' X(', N, ')=', XI(N)
          CALL TRBK
          STOP 'ERROR'
        ENDIF
        IF (present(LPARAM)) THEN
          IF ((LPARAM > 0) .AND. (LPARAM /= L)) THEN
            !If LPARAM has been specified (and is > 0), it must have the correct value
            WRITE (hCPR,'(A, I5)') '**** LPARAM =', LPARAM
            WRITE (hCPR,'(A, I5)') '**** CORRECT L =', L
            STOP ' *** ERROR IN SUBR. SPLINPO - WRONG INDEX LPARAM'
          ENDIF
        ENDIF

      ENDIF

      !Note: Unlike in other languages, the following two IF statements
      !      cannot be combined, because the second check would be made
      !      even if present(LPARAM) already returns FALSE. Therere it
      !      would cause a crash if LPARAM is not set 
      IF (present(LPARAM)) THEN
        IF (LPARAM >= 0) THEN
          !Preset L if LPARAM has been specified and is inside 2..N
          L = LPARAM
        ENDIF
      ENDIF

C***  DETERMINATION OF THE COEFFICIENTS P1, P2, P3, P4 (CF. SUBR. CUBIC)
 
C***  SET UP THE COEFFICIENT MATRIX
      D1=1./(XI(L)-XI(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,N)

C      WRITE (hCPR,*) 'Debug: L-2=', L-2, '  LA=', LA

C***  Entry point for a subsequent call with same x point, 
C***     but different function 
C      ENTRY SPLINPOX_SAME_X (F, X, FI, XI, N, SAFE)

C***  FUNCTION TO BE INTERPOLATED: FI
      F1 = FI(L-1)
      F2 = FI(L)

      IF (.TRUE.) THEN
C***     Standard version: zentrierte Ableitungen an den Stuetzstellen. Das 
C***     ist bei nicht aequidistanten Stuetzstellen fragwuerdig, verringert 
C***     aber andererseits das Ueberschwingen insbesondere wenn man nicht 
C***     MONO verwendet 
         F3 = (FI(L) - FI(LA)) / (XI(L) - XI(LA))
         F4 = (FI(LB) - FI(L-1)) / (XI(LB) - XI(L-1))

      ELSE
C***  Alternative Version wrh 16-May-2002 13:39:29
C***  Statt der zentrierten Ableitung werden gewichtete Mittel der
C***  Ableitungen der angrenzenden Intervalle genommen. Das entspricht 
C***  der Ableitung eines Parabel-Fits durch drei Punkte.  
         FS0 = (FI(L) - FI(L-1)) / (XI(L) - XI(L-1))
         IF (L == 2) THEN 
            F3 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L) - XI(L-2))
            FSM = (FI(L-1) - FI(L-2)) / (XI(L-1) - XI(L-2))
            F3 = P * FSM + (1.-P) * FS0
         ENDIF
         IF (L == N) THEN 
            F4 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L+1) - XI(L-1))
            FSP = (FI(L+1) - FI(L)) / (XI(L+1) - XI(L))
            F4 = P * FSP + (1.-P) * FS0
         ENDIF
      ENDIF

C***  SET TRUE FOR MONO OPTION
      IF (.TRUE.) THEN

ccc   Diese bis heute (3-Sep-2002) verwendete Version erscheint mir 
ccc   merckwuerdig und an mehreren Stellen fehlerhaft! Ich lasse sie 
ccc   aus Dokumentationsgruenden hier stehen. Nachfolgend dann eine 
ccc   Version nach heutiger Erkenntnis. wrh  
c       S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
c       IF (LA .NE. L-2 .AND. LB .NE. L) THEN
c         S3 = S4
c         S5 = S4
c       ELSE
c         IF (LA .EQ. L-2) THEN
c           S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
c         ELSE
c           S3 = 1.4 * S4 - 0.5 * F3
c         ENDIF
c         IF (LB .EQ. L+1) THEN
c           S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
c         ELSE
c           S5 = 1.5 * S4 - 0.5 * F4
c         ENDIF
c       ENDIF

          S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
C***   We are not in the first interval:
          IF (LA /= L-2) THEN
            S3 = S4
          ELSE
            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
          ENDIF
C***   We are not in the last interval:
          IF (LB /= L+1) THEN
             S5 = S4
          ELSE
             S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
          ENDIF

       F3 = (SIGN(1.0,S3)+SIGN(1.0,S4))*MIN(ABS(S3),ABS(S4),0.5*ABS(F3))
       F4 = (SIGN(1.0,S4)+SIGN(1.0,S5))*MIN(ABS(S4),ABS(S5),0.5*ABS(F4))

      ENDIF
 
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
 

C***  EVALUATION OF THE INTERPOLATION POLYNOMIAL
      DXM = X - XI(L-1)
      DX  = XI(L) - X
      F = (P1 * DXM * DXM + P2 ) * DXM
     >  + (P3 * DX  * DX  + P4 ) * DX

C***  Calculation of derivatives (optional)
      !added on 06.10.2011 to provide the same functionality as goetz
      IF (present(DFDX)) THEN
        DFDX = 3. * P1 *  DXM * DXM + P2 - 3. * P3 * DX * DX - P4
      ENDIF
      IF (present(D2FD2X)) THEN
        D2FD2X = 6. * (P1 *  DXM + P3 * DX)
      ENDIF

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
      SUBROUTINE TRBK

      WRITE (*,*) 'TRACE BACK FACILITY (TRBK) IS NOT AVAILABLE AT',
     >            ' DEC/OSF1 AT THE MOMENT'

      CALL ABORT()
      RETURN
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
