      SUBROUTINE SETXJL (ITNEL, LASTIND, INDLOW, INDNUP, bUSEALO,
     >                   SLNEW, SLOLD, OPAL, OPALOLD, XJLAPP, BDIAG,
     >                   NBLENDS, IBLENDS, MAXLAP, LASTINDAUTO, NDIM,
     >                   ND, EINST, RUDLINE, ELEVEL, EN, EN1, WEIGHT, 
     >                   NATOM, XJL, XLAMAPPMEAN, XLAMZERO, RSTAR,
     >                   XJC, XJCAPP, XLAMBDA, NF, L, TL, TLOLD, 
     >                   NOTEMP, ENTOTL, VDOPUNIT, KRUDAUT, MAXAUTO,
     >                   WFELOW, WFENUP, SCNEW, SCOLIND, SCNEIND,
     >                   POPMIN)
C***********************************************************************
C***  CALCULATE LINE RADIATION FIELD WITH APPROXIMATE LAMBDA OPERATOR (ALO) TERMS
C***  USING THE PRECALCULATED LINE OPERATORS CALCULATES IN COLI->FREQUINT
C***
C***  ALO is taken from the diagonal elements of the LAMBDA operator,
C***  see COLI->COLIMO and its subroutine INVTRI
C***  
C***  THIS SUBROUTINE ALSO PROVIDES  (BY CALLING LIOP) THE LINE OPACITIES
C***  AND THE NEW LINE SOURCE FUNCTION
C***  ATTENTION: XJLAPP AND XJL ARE NOT DEFINED FOR RUDIMENTAL LINES!
C***          -  OPAL AND SLNEW ARE NOT DEFINED FOR RUDIMENTAL LINES
C*** 
C***    EN :  current population numbers for current depth point
C***    EN1:  old population numbers before current STEAL
C***
C***  CALLED FROM: SUBROUTINE COMA
C***********************************************************************
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: L, ND, NF, NDIM, NATOM, ITNEL,
     >                       LASTIND, LASTINDAUTO, MAXAUTO, MAXLAP
      REAL, INTENT(IN) :: POPMIN, ENTOTL, RSTAR, VDOPUNIT, TL, TLOLD
 
      REAL, DIMENSION(NDIM) :: ELEVEL, EN, EN1, WEIGHT
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      
      REAL, DIMENSION(ND, LASTINDAUTO) :: XJL, XLAMAPPMEAN
      REAL, DIMENSION(ND, LASTIND) :: WFELOW, WFENUP
      REAL, DIMENSION(LASTIND) :: XJLAPP, SLNEW, SLOLD, 
     >                            SCOLIND, SCNEIND,
     >                            OPAL, OPALOLD, XLAMZERO
      INTEGER, DIMENSION(LASTIND) :: INDNUP, INDLOW, NBLENDS
      INTEGER, DIMENSION(MAXLAP, LASTIND) :: IBLENDS
      REAL, DIMENSION(NF) :: XJCAPP, XLAMBDA, SCNEW
      REAL, DIMENSION(ND,NF) :: XJC 
            
      INTEGER, DIMENSION(MAXAUTO) :: KRUDAUT

      LOGICAL :: NOTEMP, OLDOPAL, BNOCORE
      LOGICAL, DIMENSION(LASTIND) :: BDIAG, bUSEALO, RUDLINE
      
      REAL, EXTERNAL :: BNUE
      
      REAL :: XLAM, WAVENUM, ETAL, DELTASL, DELTASC, DPNUP, DPLOW,
     >        XJCIND, XJCAPPIND, DELTAXJC, XSHIFT,
     >        SUMSLNEW, SUMSLOLD, SUMOPALNEW, SUMOPALOLD
      INTEGER :: IND, NUP, LOW, LB, INDLB, K

C***  Speed of light in km/s      
      REAL, PARAMETER :: CLIGHTKMS = 2.99792458E5

      
      
C***  First iteration for a depth point: SL_new = SL_old => no contribution from ALO 
      OLDOPAL = (L > 1) .AND. (ITNEL == 1)

C***  Prepare SLNEW(IND) and init XJLAPP for all line indices ----------
      DO IND=1,LASTIND
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)

C***  FOR RUDIMENTAL LINES: OPAL, SLNEW    NOT DEFINED
        IF (RUDLINE(IND)) THEN
          OPAL(IND)=5HUNDEF
          SLNEW(IND)=5HUNDEF
          CYCLE
        ENDIF

        XJLAPP(IND) = XJL(L, IND)

        XLAM = XLAMZERO(IND)

        IF (OLDOPAL) THEN
           OPAL(IND) = OPALOLD(IND)
        ELSE
C***       Call LIOP only for the current depth point
C***       (i.e. with ND = 1 and all depth-arrays reduced to one REAL variable)        
           CALL LIOP (EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     >         1,XLAM, ENTOTL, EN,RSTAR,OPAL(IND),ETAL,VDOPUNIT)
C***       Note: LIOP returns ETAL / DNUEDOP and OPAL / DNUEDOP, but since we
C***             always build the source function SL = ETAL / OPAL, the 
C***             VDOP-dependency cancels out (assuming coherent redis.)
        ENDIF

C***  LASER SECURITY
        IF (OPAL(IND) <= .0) THEN
            SLNEW(IND) = .0
            OPAL(IND)  = .0
            bUSEALO(IND) = .FALSE.
C***        If we have a LASER siuation, do not use ALO => cycle           
            CYCLE
        ENDIF

        IF (OLDOPAL) THEN
           SLNEW(IND) = SLOLD(IND)
C***       No ALO contribution => cycle           
           CYCLE
        ELSE
           SLNEW(IND)=ETAL/OPAL(IND)
        ENDIF 

      ENDDO
C***  End of preparation loop ------------------------------------------
        
C***  If the line opacities have not changed,
C***    we do not need the second loop since XJLAPP = XJL
      IF (OLDOPAL) GOTO 102      
        
C***  Second IND loop to fill XJL --------------------------------------
C***  (double structure required due to blends)        
      DO IND=1,LASTIND
      
        IF (BDIAG(IND)) THEN
C***      IRON LINES
          DPLOW = EN(LOW) - EN1(LOW)
          DPNUP = EN(NUP) - EN1(NUP)
          XJLAPP(IND) = XJL(L, IND) +
     >                    WFELOW(L, IND)*DPLOW + WFENUP(L, IND)*DPNUP
        ELSEIF (bUSEALO(IND)) THEN
C***      FOR NON-IRON LINES, USE ALO(IND) CALCULATED IN COLI->FREQUINT        

          DELTASL = 0.
          SUMSLNEW = 0.
          SUMSLOLD = 0.
          SUMOPALNEW = 0.
          SUMOPALOLD = 0.
          DO LB=1, NBLENDS(IND)
C***        Consider contributions to XJL from all blending lines          
            INDLB = IBLENDS(LB,IND)
c            IF (.NOT. bUSEALO(INDLB) .OR. OPALOLD(INDLB) <= 0.) CYCLE
            IF (OPALOLD(INDLB) <= 0.) CYCLE
              SUMOPALNEW = SUMOPALNEW + OPAL(INDLB)
              SUMOPALOLD = SUMOPALOLD + OPALOLD(INDLB)
              SUMSLOLD = SUMSLOLD + SLOLD(INDLB) * OPALOLD(INDLB)
              SUMSLNEW = SUMSLNEW + SLNEW(INDLB) * OPAL(INDLB)
ccc            DELTASL = DELTASL + SLNEW(INDLB) - SLOLD(INDLB)
c            IF (INDLB /= IND) THEN
C***          Account for the fact that other line blending into the 
C***            current IND transition do not have the full contribution
C***            but are weighted with the contribution of their profile 
C***            function at the rest wavelength of the considered transition.
c              XSHIFT = ( XLAMZERO(INDLB)/XLAMZERO(IND)-1. )
c     >                                         * CLIGHTKMS / VDOPUNIT
c              DELTASL = DELTASL * EXP(-XSHIFT*XSHIFT)
c            ENDIF    
c            XJLAPP(IND) = XJLAPP(IND) + DELTASL * XLAMAPPMEAN(L, INDLB)
          ENDDO
          IF (SUMOPALNEW > 0. .AND. SUMOPALOLD > 0.) THEN
            SUMSLNEW = SUMSLNEW / SUMOPALNEW
            SUMSLOLD = SUMSLOLD / SUMOPALOLD
            DELTASL = SUMSLNEW - SUMSLOLD
            XJLAPP(IND) = XJLAPP(IND) + DELTASL * XLAMAPPMEAN(L, IND)
          ENDIF
          
C***      Consider the change of the continuum contribution
C***      by interpolating the original and approximated XJC at the 
C***      wavelength of the current transition
          WAVENUM = 1.E8 / XLAM
          CALL XRUDI (XJCIND,    WAVENUM, XJC,    XLAMBDA, ND, NF, L)
          CALL XRUDI (XJCAPPIND, WAVENUM, XJCAPP, XLAMBDA,  1, NF, 1)
          DELTAXJC = XJCAPPIND - XJCIND
          XJLAPP(IND) = XJLAPP(IND) + DELTAXJC
C***      Alt. approach
          CALL XRUDI (SCNEIND(IND),WAVENUM,SCNEW,XLAMBDA,1,NF,1)
c          DELTASC = SCNEIND(IND) - SCOLIND(IND)
c          XJLAPP(IND) = XJLAPP(IND) + DELTASC * XLAMAPPMEAN(L, IND)
        ENDIF
        
C***    Failsafe to prevent any negative radiation field
        IF (XJLAPP(IND) < 0.) THEN
          XJLAPP(IND) = XJL(L, IND)
          bUSEALO(IND) = .FALSE.
          IF (BDIAG(IND)) THEN
C***        For iron lines, prevent the usage of the WFE values for the derivatives          
            WFELOW(L, IND) = 0.
            WFENUP(L, IND) = 0.
          ENDIF
        ENDIF
        
C***    ONLY IF TEMPERATURE CORRECTIONS ARE APPLIED:
C***    ADD SPECIAL TERM AT INNER BOUNDARY, WHICH ACCOUNTS FOR THE
C***      DIRECT TEMPERATURE-DEPENDENCE OF THE RADIATION FIELD
C***      VIA THE BOUNDARY CONDITION
c        IF (.NOT. NOTEMP  .AND.  L .EQ. ND) THEN 
c           XLAM=XLAMZERO(IND)
c           XJLAPP(IND)=XJLAPP(IND)+0.5*(BNUE(XLAM,TL)-BNUE(XLAM,TLOLD))
c        ENDIF
c  we have to remove this if we consider the effect of XJC directly since this is included there
        
      ENDDO
C***  End line loop ----------------------------------------------------

C***  XJL --> XJLAPP FOR DR-TRANSITIONS
  102 CONTINUE
      DO IND = LASTIND+1, LASTINDAUTO
         IF (KRUDAUT(IND-LASTIND) == 1) CYCLE
         XJLAPP(IND) = XJL(L,IND)
      ENDDO

      RETURN
      END
