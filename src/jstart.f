      SUBROUTINE JSTART (NF,XLAMBDA,KEY,ND,R,T,XJC,XJL,ELEVEL, N,EINST,
     $                   NDIM, INDNUP, INDLOW, LASTINDALL, LASTIND,
     >                   NAUTO, WSTABIL, R23, TAUROSS, 
     >                   LTESTART, BLACKEDGE, bOLDJ, XJCold, XLAMBDAold, Rold,
     >                   NFold, NDold, KRUDAUT)
C******************************************************************************
C***  START APPROXIMATION OF THE RADIATION FIELD
C***  XJC: CONTINUUM RADIATION FIELD
C***  XJL: LINE RADIATION FIELD
C******************************************************************************
 
      INTEGER, INTENT(IN) :: NFold, NDold, LASTIND, LASTINDALL, NAUTO, ND, NF
      REAL, DIMENSION(NF) :: XLAMBDA
      CHARACTER(8), DIMENSION(NF) :: KEY
      REAL, DIMENSION(ND) :: R, T, XJC, XJL, TAUROSS
      REAL, DIMENSION(NDIM) :: ELEVEL
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      INTEGER, DIMENSION(LASTINDALL) :: INDNUP, INDLOW
      INTEGER, DIMENSION(NAUTO) :: KRUDAUTO
      REAL, DIMENSION(NAUTO) :: WSTABIL
      REAL, DIMENSION(NDold) :: Rold
      REAL, DIMENSION(NFold) :: XLAMBDAold
      REAL, DIMENSION(NDold,NFold) :: XJCold
      CHARACTER(8) :: NAME
      REAL, DIMENSION(4) :: XJSP, TAUSP

      LOGICAL :: LTESTART, bOLDJ

C***  Find temperature at tau = 2/3
      CALL LIPO (T23, R23, T, R, ND)

C***  Find depth index inside tau = 1
      L1  = ISRCHFGT(ND,TAUROSS, 1, 1.5)
      IF (L1 .EQ. 0) L1 = ND-1

C***  Find depth index outside tau = 1/3
      L13 = ISRCHFGT(ND,TAUROSS, 1, 0.333333)
      IF (L13 .GT. 1) L13 = L13 - 1
      IF (L13 .EQ. 0) L13 = 1

C***  Construct vector with 4 elements for spline interpolation
C***   in order to smooth over the tau=1 discontinuity 
      IF (L13 .EQ. 1 .OR. L1 .GE. (ND-1)) THEN 
         NSP = 0
      ELSE
         NSP = 4
         ISP1 = L13 - 1       
         ISP2 = L13       
         ISP3 = L1       
         ISP4 = L1+1       
         TAUSP(1) = TAUROSS(ISP1)
         TAUSP(2) = TAUROSS(ISP2)     
         TAUSP(3) = TAUROSS(ISP3)     
         TAUSP(4) = TAUROSS(ISP4)     
      ENDIF

C***  CONTINUUM RADIATION FIELD XJC  *****************************************
C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS
C***  ( DEPTH VEKTOR AT EACH FREQUENCY POINT )
      DO 6 K=1,NF
      XLAM=XLAMBDA(K)
 
      IF (bOLDJ) THEN
C***     BRANCH FOR XJC FROM OLD MODEL
         DO L=1, ND
C***       Determine the corresponding depth point in the old model
C***       (This could be optimized with a more precise interpolation)
           lofind: DO LL=1, NDold
             IF (R(L) >= Rold(LL)) THEN
               Lold = LL
               EXIT lofind
             ENDIF
           ENDDO lofind
           WN = 1.E8/XLAM
           IF (XLAM < XLAMBDAold(1)) THEN
             XJC(L) = MIN(XJCold(L,1), BNUE(XLAM, T(L)))
           ELSEIF (XLAM > XLAMBDAold(NFold)) THEN 
             XJC(L) = MIN(XJCold(L,NFold), BNUE(XLAM, T(L))) 
           ELSE
             CALL XRUDI (XJC(L),WN,XJCold,XLAMBDAold,NDold,NFold,Lold)
           ENDIF
         ENDDO
      ELSEIF (LTESTART) THEN
C***     BRANCH FOR LTESTART
         DO L=1,ND
           XJC(L)=BNUE(XLAM,T(L))
           IF (R(L) > R23 .AND. XLAM < BLACKEDGE) THEN
             XJC(L) = XJC(L) * EXP(-20.*(2./3. - TAUROSS(L)))
           ENDIF
         ENDDO   
      ELSE
C***     BRANCH FOR GEOMETRICAL DILUTION OF BLACKBODY FIELD 
C***       AT PHOTOSPHERIC RADIUS
         BPHOT = BNUE(XLAM,T23)
cc         IF (XLAM .LT. BLACKEDGE) BPHOT = .0

         DO 20 L=1,ND
            IF (R23 .LT. R(L)) THEN
               W=0.5
               ARG=1. - R23*R23 / (R(L)*R(L))
               IF (ARG .GT. .0) W=0.5*(1.-SQRT(ARG))
               XJC(L) = BPHOT * W
cc naechste Zeile testweise
               IF (XLAM .LT. BLACKEDGE) XJC(L) = XJC(L) * 
     >             EXP(-20.*(2./3. - TAUROSS(L)))
            ELSE
               XJC(L)=BNUE(XLAM,T(L))
            ENDIF
   20    CONTINUE

C***  replaced by a new version with Spline interpolation, wrh  3-Mar-2006
         IF (NSP .EQ. 4) THEN
            XJSP(1) = XJC(ISP1)
            XJSP(2) = XJC(ISP2)
            XJSP(3) = XJC(ISP3)
            XJSP(4) = XJC(ISP4)
            DO L = L13+1, L1-1
               CALL SPLINPO (XJC(L), TAUROSS(L), XJSP, TAUSP, NSP)
            ENDDO
         ENDIF
      ENDIF
 
      WRITE (UNIT=NAME, FMT='(A3,I4,A1)') 'XJC', K, ' '
      CALL WRITMS(3,XJC,ND,NAME,-1,0, IERR)
    6 CONTINUE
C******************************************************************************
 
C***  LINE RADIATION FIELD XJL  ***********************************************
C***  ( DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND )
      DO 99 IND=1, LASTINDALL
C***    calc XLAM: Branch for normal lines and for superlevels
        IF (IND .LE. LASTIND .OR. IND .GT. LASTIND+NAUTO) THEN
           J=INDNUP(IND)
           I=INDLOW(IND)
           IF (EINST(I,J) .EQ. -2.) CYCLE
           XLAM=1.E8/(ELEVEL(J)-ELEVEL(I))
        ELSE
           INDDR = IND - LASTIND
           IF (KRUDAUT(INDDR) .EQ. 1) CYCLE
           XLAM = 1.E8 / WSTABIL(INDDR)
        ENDIF

C***    THIS VERSION: SAME APPROXIMATION AS FOR CONTINUUM
        IF (bOLDJ) THEN
C***      Branch for old radiation field: Just interpolate XJC from above
          DO L=1,ND
C***        Determine the corresponding depth point in the old model
C***        (This could be optimized with a more precise interpolation)
            lofind2: DO LL=1, NDold
              IF (R(L) >= Rold(LL)) THEN
                Lold = LL
                EXIT lofind2
              ENDIF
            ENDDO lofind2
            WN = 1.E8/XLAM
            CALL XRUDI (XJL(L),WN,XJCold,XLAMBDAold,NDold,NFold,Lold)
          ENDDO        
        ELSEIF (LTESTART) THEN
C***      BRANCH FOR LTESTART
          DO L=1,ND
            XJL(L)=BNUE(XLAM,T(L))
            IF (R(L) > R23 .AND. XLAM < BLACKEDGE) THEN
              XJL(L) = XJL(L) * EXP(-20.*(2./3. - TAUROSS(L)))
            ENDIF
          ENDDO
        ELSE
C***    BRANCH FOR GEOMETRICAL DILUTION OF BLACKBODY FIELD
          BPHOT = BNUE(XLAM,T23)

cc         IF (XLAM .LT. BLACKEDGE) BPHOT = .0

          DO 9 L=1,ND
            IF (R23 .LT. R(L)) THEN
               W=0.5
               ARG=1. - R23*R23 / (R(L)*R(L))
               IF (ARG .GT. .0) W=0.5*(1.-SQRT(ARG))
               XJL(L)=BPHOT*W
cc naechste Zeile testweise
               IF (XLAM .LT. BLACKEDGE) XJL(L) = XJL(L) * 
     >             EXP(-20.*(2./3. - TAUROSS(L)))
            ELSE
               XJL(L)=BNUE(XLAM,T(L))
            ENDIF
    9     CONTINUE

C***     Test output for the smoothing 
C***     For activation, set INDPLO=137 for HeII Lyman-alpha, for instance
         INDPLO = 0
         IF (IND .EQ. INDPLO) THEN
          OPEN (2, FILE='PLOT')
          NSPPLOT = L1 - L13 + 3
          CALL PLOTANFS (2,' ', ' ',
     $        'TAUROSS',
     $        'Radiation intensity',
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     >        TAUROSS(L13-1), XJL(L13-1), NSPPLOT, 
     >        'SYMBOL=2 SIZE=0.2 COLOR=4')
          CALL JSYMSET ('G2','TRANSFER')
         ENDIF

C***     replaced by a new version with Spline interpolation, wrh  3-Mar-2006
         IF (NSP .EQ. 4) THEN
            XJSP(1) = XJL(ISP1)
            XJSP(2) = XJL(ISP2)
            XJSP(3) = XJL(ISP3)
            XJSP(4) = XJL(ISP4)
            DO L = L13+1, L1-1
               CALL SPLINPO (XJL(L), TAUROSS(L), XJSP, TAUSP, NSP)
            ENDDO
         ENDIF

C***     Continue testplot (now after smoothing)
         IF (IND .EQ. INDPLO)
     >    CALL PLOTCONS (2, 
     >        TAUROSS(L13-1), XJL(L13-1), NSPPLOT, 
     >        'SYMBOL=8 SIZE=0.2 COLOR=2')

      ENDIF

C***  WRITE RADIATION FIELD TO THE MODEL FILE 
      IF (IND <= 9999) THEN
        WRITE (UNIT=NAME, FMT='(A3,I4,A1)') 'XJL',IND,' '
      ELSE
        WRITE (UNIT=NAME, FMT='(A3,I5)') 'XJL', IND
      ENDIF
      CALL WRITMS (3,XJL,ND,NAME,-1,0, IERR)
   99 CONTINUE
C***********************************************************************

      RETURN
      END
