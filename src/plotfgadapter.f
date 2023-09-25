      SUBROUTINE PLOTFGADAPTER(PLOTOPT, XMSTAR, TAUROSS, ARAD, 
     >                         VELO, RADIUS, ND, ATMEAN, ENTOT, T, 
     >                         VTURB, XMU, RSTAR, Rcritical, WORKRATIO,
     >                         bNoARAD, MODHEAD, JOBNUM, KANAL)      
C***********************************************************************
C***
C***
C***  called from STEAL
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, KANAL, JOBNUM
      REAL, INTENT(IN) :: XMSTAR, ATMEAN, RSTAR, Rcritical, WORKRATIO
      LOGICAL, INTENT(IN) :: bNoARAD
      
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS, VELO, RADIUS, 
     >                                   ENTOT, T, VTURB, XMU
      REAL, DIMENSION(ND-1), INTENT(IN) :: ARAD
      
      REAL, DIMENSION(NDMAX) :: F_raw, G_raw, CKL, ALPHAL, GAMMARAD, RI,
     >                          GAMMARADL, AMACH, RCM, VCM
      
      CHARACTER(8) :: KEY
      CHARACTER(100) :: MODHEAD
      
      REAL :: XMG, DOTM4P
      INTEGER :: L, INDOUT
      
      CHARACTER PLOTOPT*(*)
      
C***  Physical constants      
      REAL, PARAMETER :: AMU = 1.66E-24             !Atomic mass unit (gramm)
      REAL, PARAMETER :: GCONST = 6.6727E-8         !Gravitational Constant (CGS)
      REAL, PARAMETER :: RGAS = 8.3145E7            !GAS CONSTANT in CGS 
      REAL, PARAMETER :: XMSUN = 1.989E33           !Solar mass (CGS = in g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOWNPLOT = 2    !current plot

      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTFGADAPTER: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTFGADAPTER: NDMAX DIM. INSUFFICIENT'
         WRITE (hCPR,'(A)') 'PLOTFGADAPTER: PLOTFGSTRAT SKIPPED'
         RETURN
      ENDIF
      
      IF (bNoARAD) THEN
         WRITE (hCPR,'(A)') 'PLOTFGADAPTER: Radiative acceleration ' //
     >      'not yet calculated ***'
         WRITE (hCPR,'(3A)') 'PLOTFGADAPTER: ',TRIM(PLOTOPT),' SKIPPED'
         RETURN
      ENDIF      
      
C***  Fill all quantities that are not needed outside hydro calculations 
C***  with zero and calculate ("true") sound speed (incl. turbulence)
      INDOUT = -1
      DO L=1, ND
        F_raw(L) = 0.
        G_raw(L) = 0.
        CKL(L) = 0.
        ALPHAL(L) = 0.
        AMACH(L) = SQRT(RGAS * T(L) / XMU(L) + (VTURB(L)*1.E5)**2) 
C***    Convert RADIUS and VELO into CGS units
        RCM(L) = RADIUS(L) * RSTAR
        VCM(L) = VELO(L) * 1.E5    
      ENDDO
      
      XMG = GCONST * XMSTAR * XMSUN
      
C***  Calculate GAMMARAD (ARAD is given on interstices)      
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
        GAMMARAD(L) = ARAD(L) / XMG * (RI(L)*RSTAR)**2
      ENDDO
      
C***  Convert GAMMARAD back to regular grid points      
      GAMMARADL(1) = GAMMARAD(1)
      DO L=2, ND-1
        CALL SPLINPOX(GAMMARADL(L), RADIUS(L), GAMMARAD, RI, ND-1)
      ENDDO
      GAMMARADL(ND) = GAMMARAD(L)
      
C***  Calculate current mass-loss rate via continuity equation      
      DOTM4P = ENTOT(1) * ATMEAN * AMU
     >               * VELO(1)*1.E5
     >               * RADIUS(1)*RADIUS(1)*RSTAR*RSTAR
      
C***  Call actual plotting subroutine
      KEY = 'TRANSFER'      
      CALL PLOTFGSTRAT(PLOTOPT, RCM, RSTAR, VCM, ND, DOTM4P, 
     >                 AMACH, XMG, CKL, ALPHAL,
     >                 GAMMARADL, GAMMARADL, F_raw, G_raw,
     >                 ENTOT, TAUROSS,
     >                 MODHEAD, INDOUT, JOBNUM, KANAL, KEY)


      RETURN 

      END
