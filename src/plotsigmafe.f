      SUBROUTINE PLOTSIGMAFE(PLOTOPT, MODHEAD, JOBNUM,
     >                       SIGMAFE, SIGMAFEUL, SIGNU3FE, 
     >                       INDRB, MAXFEIND, 
     >                       LASTINDALL, LASTFE, IFRBSTA, IFRBEND, 
     >                       LASTKON, INDRF, KONTLOW, KONTNUP,
     >                       LEVEL, N, INDNUP, INDLOW, NOM,
     >                       ELEVEL, EION,
     >                       INDEXMAX, NFEREADMAX, MAXATOM, KODAT,
     >                       VDOPFE, DXFE, XLAM0FE, bFEULSEP, bOwn)
C***  TODO: Update routine to plot also bound-free sections
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MAXFEIND, INDEXMAX, NFEREADMAX, MAXATOM,
     >                       LASTINDALL, LASTFE, LASTKON, N, JOBNUM
      REAL, INTENT(IN) :: XLAM0FE, DXFE, VDOPFE
      INTEGER, DIMENSION(LASTINDALL), INTENT(IN) :: INDNUP, INDLOW
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: KONTNUP, KONTLOW
      REAL, DIMENSION(INDEXMAX), INTENT(IN) :: SIGMAFE, SIGMAFEUL, 
     >                                         SIGNU3FE
      REAL, DIMENSION(N) :: ELEVEL, EION
      INTEGER, DIMENSION(N), INTENT(IN) :: NOM
      INTEGER, DIMENSION(MAXATOM), INTENT(IN) :: KODAT
      CHARACTER(10), DIMENSION(N) :: LEVEL
      CHARACTER PLOTOPT*(*)
      
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, IFRBEND, IFRBSTA, INDRF
      
      INTEGER, PARAMETER :: LAMDIMMAX =  500000         !should have the size of NFEREADMAX
      REAL, DIMENSION(LAMDIMMAX) :: X, Y, Y2

      CHARACTER(5) :: CIFE, CIND
      CHARACTER(10) :: LEV1, LEV2
      CHARACTER(100) :: MODHEAD
      CHARACTER(110) :: HEADLINE
      
      REAL :: XLOGSTEP, XLAMEDGE, EDGE, XLAMNULL
      INTEGER :: NPTS, IndexFE, IndexFELAM, I, J, ILAM, IND, KON,
     >           hPLOT, ILEV1, ILEV2, IERR, NPAR, LOW, NUP, 
     >           ILAMSTART, ILAMEND, IndexSIGMAINT, IndexFELAMstart
      
      LOGICAL :: bOwn, bFEULSEP
      
      INTEGER, EXTERNAL :: IDX

      !Physical constants
      REAL, PARAMETER :: CLIGHT = 2.99792458E10     !Speed of Light in cm/s

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      IF (LAMDIMMAX < NFEREADMAX) THEN
        WRITE (hCPR,*) '**** PLOTSIGMAFE: LAMDIMMAX too small '
     >    // '--> code might crash'
      ENDIF
      
      WRITE (0,*) 'PLOTSIGMAFE called'
C***  DECODE LEVELS (can be index number or name)      
      ILEV1 = 0
      ILEV2 = 0
      
C***  Debug output: list of levels      
c      DO I=1, N
c        WRITE (0,*) 'LEVEL: ', I, LEVEL(I)
c      ENDDO
      
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR == 4) THEN
        CALL SARGV (PLOTOPT, 3, LEV1)
        READ (LEV1(:IDX(LEV1)), '(I4)', IOSTAT=IERR) ILEV1
        IF (IERR > 0) THEN
          DO I=1, N
            IF (LEVEL(I) == LEV1) THEN
              ILEV1 = I
              EXIT
            ENDIF
          ENDDO      
        ELSE
          LEV1 = LEVEL(ILEV1)          
        ENDIF
        CALL SARGV (PLOTOPT, 4, LEV2)
        READ (LEV2(:IDX(LEV2)), '(I4)', IOSTAT=IERR) ILEV2
        IF (IERR > 0) THEN
          DO I=1, N
            IF (LEVEL(I) == LEV2) THEN
              ILEV2 = I
              EXIT
            ENDIF
          ENDDO      
        ELSE
          LEV2 = LEVEL(ILEV2)          
        ENDIF
      ELSE
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid parameters ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF            

      IF (ILEV1 == 0 .OR. ILEV2 == 0) THEN
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid levels ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF
      
C***  TODO: Update for SIGMAFEUL, SIGNU3FE case      
      IND = 0
      DO I=1, LASTINDALL
        IF (INDNUP(I) == ILEV1 .AND. INDLOW(I) == ILEV2) THEN
          LOW = ILEV2
          NUP = ILEV1
          IND = I
          EXIT
        ELSEIF (INDNUP(I) == ILEV2 .AND. INDLOW(I) == ILEV1) THEN
          LOW = ILEV1
          NUP = ILEV2
          IND = I
          EXIT
        ENDIF
      ENDDO

      IF (IND == 0) THEN
        KON = 0
        DO I=1, LASTKON
          IF (KONTNUP(I) == ILEV1 .AND. KONTLOW(I) == ILEV2) THEN
            LOW = ILEV2
            NUP = ILEV1
            KON = I
            EXIT
          ELSEIF (KONTNUP(I) == ILEV2 .AND. KONTLOW(I) == ILEV1) THEN
            LOW = ILEV1
            NUP = ILEV2
            KON = I
            EXIT
          ENDIF
        ENDDO
        IF (KON == 0) THEN
          WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid transition ************'
          WRITE (hCPR,*) '***** The following plot was aborted:'
          WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
          RETURN
        ELSEIF (NOM(LOW) /= KODAT(26)) THEN
          WRITE (hCPR,*) 'PLOTSIGMAFE: Non-iron transition ************'
          WRITE (hCPR,*) '***** The following plot was aborted:'
          WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
          RETURN
        ELSE
          EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
          XLAMEDGE = 1.E8/EDGE
          WRITE (hCPR,*) 'Transition found: ', KON, XLAMEDGE
        ENDIF        
      ENDIF
                  
      
C***  Check if transition is found and really an iron bound-bound transition
      IF (IND < (LASTINDALL - LASTFE) .AND. KON == 0) THEN
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid transition ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF
      
      
C***  PLOT FILE INITIALIZATION
      IF (bOwn) THEN
        hPLOT = 37       !write to seperate plot file
        OPEN (hPLOT, FILE='sigmafe.plot', STATUS='UNKNOWN')
        WRITE (hCPR,*) 'SIGMAFE PLOTTED FOR IND = ', IND
      ELSE
        hPLOT = 2       !write to file PLOT (becomes steal.plot)
        OPEN (hPLOT, FILE='PLOT', STATUS='UNKNOWN')
        CALL JSYMSET ('G2','TRANSFER')
        CALL REMARK ('SIGMAFE TO BE ROUTED')
      ENDIF
      
      HEADLINE = 'SIGMAFE:'//MODHEAD(13:)
C      HEADLINE = 'SIGMAFE: to be done'
      
      WRITE (hPLOT, '(A,I8)') '*NEXTPLOT: SIGMAFE ', IND
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT,'(A)') '\COLOR=1'
      WRITE (hPLOT,'(A)') '\PEN=1'
      
      XLOGSTEP = LOG10(1. + VDOPFE*1.E5*DXFE/CLIGHT)
      IF (IND > 0) THEN
        IndexFE = IND - (LASTINDALL - LASTFE)     !convert total IND to fe-only IND number
        IndexFELAM = INDRB(IndexFE)               !starting index inside SIGMAFE for current FeIND
        ILAMSTART = IFRBSTA(IndexFE)              !starting wavelength index
        ILAMEND   = IFRBEND(IndexFE)              !last wavelength index
        XLAMNULL = XLAM0FE
      ELSE
C***    Bound-free transitions are always created on a fixed grid 
C        with VDOP = 30 and FSTEP = 30
C        ( check ~wrh/Blanket/Data/Fgrids/coarsek_grid )
        IndexFE = KON
        IndexFELAM = INDRF(IndexFE) + 2
        ILAMSTART = - INT(SIGMAFE(IndexFELAM - 2))
        ILAMEND   = - INT(SIGMAFE(IndexFELAM - 1))
        WRITE (hCPR,'(A,4(2X,I8))') 'Indices found: ', 
     >              IndexFE, IndexFELAM, ILAMSTART, ILAMEND
        XLOGSTEP = LOG10(1. + 30*1.E5*30/CLIGHT)
        XLAMNULL = 1001.0993601
      ENDIF
      IndexFELAMstart = IndexFELAM
      DO 
        I = IndexFELAM - IndexFELAMstart + 1
        ILAM = ILAMSTART + I - 1
C        WRITE (hCPR, *) I, ILAM, XLAM0FE * 10.**(ILAM*XLOGSTEP)
        X(I) = XLAMNULL * 10.**(ILAM*XLOGSTEP)       !lambda from index
        Y(I) = SIGMAFE(IndexFELAM) / 1.E-15          !plot in 10^(-15) cm^2
        IF (KON > 0) THEN
          IF (X(I) > XLAMEDGE) THEN
            Y2(I) = 0.
          ELSE
            Y2(I) = (X(I)/XLAMEDGE)**3
          ENDIF
        ENDIF
C        Y(I) = SIGMAFE(IndexFELAM)
        IF (ILAM == ILAMEND) THEN
          IF (KON > 0) THEN
            IndexSIGMAINT = IndexFELAM + 1
          ENDIF
          EXIT
        ENDIF
        IndexFELAM = IndexFELAM + 1
      ENDDO
      NPTS = I      

      WRITE (UNIT=CIFE, FMT='(I5)') IndexFE
      WRITE (UNIT=CIND, FMT='(I5)') IND
      WRITE (hPLOT,'(5A)') 
     >  '\LUN XMAX YMAX R-0.5 U-0.5 0.2 Fe-Transition ',
     >    TRIM(ADJUSTL(CIFE)), ' (', TRIM(ADJUSTL(CIND)), ')'
      WRITE (hPLOT,'(2A)') 
     >  '\LUN XMAX YMAX R-0.5 U-1.0 0.25 Low: ', LEVEL(LOW)
      WRITE (hPLOT,'(2A)') 
     >  '\LUN XMAX YMAX R-0.5 U-1.5 0.25 Up: ', LEVEL(NUP)
      
      CALL PLOTANFS (hPLOT,HEADLINE, '&E'//HEADLINE,
     >        '\CENTER\#l# / \A',
     >        '\CENTER\#s#&TLU&M / 10&H-15&M\,cm&H2&M',
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     >        X, Y, NPTS, 'COLOR=4')
     
       IF (KON > 0) THEN
C***     Plot hydrogenic slope for comparison
         DO J=1, NPTS
           Y2(J) = Y2(J) * SIGMAFE(IndexSIGMAINT) / 1.E-15
         ENDDO
         CALL PLOTCONS (hPLOT, X, Y2, NPTS, 'COLOR=1')
       ENDIF
     
C      CALL PLOTCONS (hPLOT, X, TCORR,  ND,   'PEN=4 COLOR=2')

C      IF (bOwn) THEN
C        CLOSE (hPLOT)
C      ENDIF


      RETURN

      END
