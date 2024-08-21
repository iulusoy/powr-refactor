      SUBROUTINE PREP_GAMMARAD (iOLDRAD, bFULLHYDROSTAT, GEDD, 
     >      IGEddFix, GAMMARAD, NDold, ARAD, GLOG, RSTAR,
     >      RADIUSold, RSTARold, XMGold, TAURold, RCONold, GEFFLOG,
     >      STAPEL, ATMEAN, XLOGL, XMSTAR, RADGAMMASTART, 
     >      GEDDRAD, iOldStratification, bGAMMARADMEAN, bSaveGEFF)
C***********************************************************************
C***  This subroutine prepares GAMMARAD and related quentities for
C***  the hydrostatic equation
C***  Called from: WRSTART
C***********************************************************************

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NDold, IGEddFix, iOLDRAD, 
     >                       iOldStratification
      REAL, INTENT(IN) :: XMGold, XLOGL, STAPEL, ATMEAN,
     >                    RSTARold, RCONold, RADGAMMASTART
      
      INTEGER, PARAMETER :: NDMAX = 200
      
      REAL, DIMENSION(NDold), INTENT(IN) :: RADIUSold, TAURold, ARAD
      REAL, DIMENSION(NDold), INTENT(OUT) :: GAMMARAD
      REAL, DIMENSION(NDMAX) :: ARADL 

      REAL, INTENT(INOUT) :: GEDD, GEFFLOG, GLOG, XMSTAR, GEDDRAD

      LOGICAL, INTENT(IN) :: bFULLHYDROSTAT, bSaveGEFF, bGAMMARADMEAN

      INTEGER :: L
      REAL :: GEDDRADMEAN, GEDDRADMID, TAUMID, TAUNORM, WTAU, RSTAR

      !Physical constants
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (NDold > NDMAX) THEN
        WRITE (hCPR,'(A)') 'PREP_GAMMARAD: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'PREP_GAMMARAD: NDMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'NDold = ', NDold, ', NDMAX = ', NDMAX
        STOP 'FATAL ERROR IN WRSTART->PREP_GAMMARAD'
      ENDIF
      
      
      GEDDRAD = 0.      
      
c      WRITE (hCPR,*) 'Preparing subsonic stratification...'
      IF (iOLDRAD == 2 .AND. bFULLHYDROSTAT) THEN
        !calculate effective g via ARAD from OLD MODEL

        IF (GEDD > .0 .AND. IGEddFix == 2) THEN
          !Gamma value has been fixed via CARDS 
          DO L=1, NDold
            GAMMARAD(L) = GEDD
          ENDDO

        ELSE

C***      Interpolate ARAD from interstices to ARADL on full depth points  
          DO L=1, NDold
            IF (L==1) THEN
              ARADL(L) = ARAD(1)
            ELSEIF (L==NDold) THEN
              ARADL(L) = ARAD(NDold-1)
            ELSE
              ARADL(L) = 0.5 * (ARAD(L-1) + ARAD(L))
            ENDIF
            GAMMARAD(L) = ARADL(L)*(RADIUSold(L)*RSTARold)**2/XMGold
C***        restrict GAMMARAD 
            GAMMARAD(L) = MIN(GAMMARAD(L), 0.9)
          ENDDO

          !Calculate mean value, needed for g-g_eff-Relation
          ! if RADGAMMASTART: OLD is used (cannot be done in DECSTAR)
c          GEDDRADMEAN = GAMMARAD(NDold) * EXP(-TAURold(NDold))
          GEDDRADMEAN = GAMMARAD(NDold) 
          TAUNORM = 0.
          IF (RCONold >= RADIUSold(NDold-1)) THEN
            GEDDRADMEAN = 0.
            DO L=NDold-1, 1, -1
              IF (RADIUSold(L) > RCONold) EXIT  
              TAUMID = 0.5 * ( TAURold(L+1) + TAURold(L) )
              WTAU = (TAURold(L+1) - TAURold(L)) * EXP(-TAUMID)
              GEDDRADMID = 0.5 * ( GAMMARAD(L) +  GAMMARAD(L+1) )
              GEDDRADMEAN = GEDDRADMEAN + GEDDRADMID * WTAU
              TAUNORM = TAUNORM + WTAU
            ENDDO
            IF (TAUNORM > 0.) THEN
              GEDDRADMEAN = GEDDRADMEAN / TAUNORM
            ENDIF
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDDRADMEAN))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = LOG10((10.**GLOG)*(1.-GEDDRADMEAN))
          ENDIF
          
          IF (bGAMMARADMEAN) THEN
            !Use mean value of the inner part instead of individual values (CARDS option)
            DO L=1, NDold
              GAMMARAD(L) = GEDDRADMEAN
            ENDDO
          ENDIF
        ENDIF

      ELSEIF (GEDD > 0. .OR. RADGAMMASTART >= 0.) THEN
C***    for new models only g_THOMSON is available, 
C***    unless RADGAMMA-START is specified
        IF (IGEddFix == 0) THEN        
          IF (RADGAMMASTART > 0) THEN
            GEDD = RADGAMMASTART
          ELSE
            GEDD = 10.**(-4.51) * (STAPEL/ATMEAN) *(10.**XLOGL) /XMSTAR
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDD))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
          ENDIF
        ENDIF
        IF (RADGAMMASTART >= 0.) THEN
          DO L=1, NDold
            GAMMARAD(L) = RADGAMMASTART
          ENDDO
        ENDIF
      ENDIF

      IF (bFULLHYDROSTAT) THEN
        IF (iOLDRAD == 2) THEN
          GEDDRAD = GEDDRADMEAN
        ELSEIF (RADGAMMASTART >= 0.) THEN
          GEDDRAD = RADGAMMASTART
        ELSE
          GEDDRAD = GEDD      !poor start with Thomson only
        ENDIF        
      ENDIF
      
C***  Only if not copying the old stratification
      IF (iOldStratification == 0) THEN
        IF (bFULLHYDROSTAT) THEN
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Full Gamma (ND) = ', GEDDRAD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG

        ELSE
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Eddington Gamma = ', GEDD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG
        ENDIF
      ENDIF
      
      RETURN 
      END
