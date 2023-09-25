      SUBROUTINE TDIFFUS(T, RADIUS, ND, TEFF, OPARND, TNDCORR,
     >                   bNoARAD, HTOTNDCOR, fTNDCOR)
C***********************************************************************
C*** Correction of the innermost temperature point via the (pure)
C*** diffusion approximation => THIS UPDATES T(ND)
C*** (requires OPARND to be given, typically from previous COLI run)
C*** 
C*** (until Feb 2016 this was performed inside the LINPOP routine)
C***      
C*** called from STEAL
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND   
      REAL, INTENT(IN) :: TEFF, OPARND, HTOTNDCOR, fTNDCOR
      REAL, INTENT(INOUT) :: TNDCORR
      LOGICAL, INTENT(IN) :: bNoARAD
      
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS
      REAL, DIMENSION(ND), INTENT(INOUT) :: T
      
      REAL :: TSAVE, TMID, TLAST, TL, DTDR, DTDRplus, DTDRminus
      INTEGER :: L, ICOUNT
      
      LOGICAL, PARAMETER :: bPRINT = .FALSE. !print debug information

      !Physical constants
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5      !Stephan-Boltzmann constant (CGS) / Pi
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
C***  Do not perform any corrections if Rosseland opacity at
C***  inner boundary is not available (i.e. no full COLI yet)
      IF (OPARND <= 0.) THEN
        TNDCORR = 0.
        RETURN
      ENDIF
      
C***  Store old boundary temperature      
      TSAVE = T(ND)
      TL = TSAVE
      
C***  Start value: always from ND-1
C***  (i.e. last point set by UNLU if TDIFFUS is active)
      T(ND) = T(ND-1)
      
C***  Iteration counter      
      ICOUNT = 0
      
      DO
        TMID = 0.5 * (T(ND-1) + TL)
        
        IF ((bNoARAD .OR. HTOTNDCOR == 0.) .OR. (fTNDCOR <= 0.)) THEN
c***      Classic version
          DTDR = 0.1875 * OPARND * TEFF*(TEFF/TMID)**3.
        ELSE
C***      Include inner boundary correction terms (seems to be unstable)
C***      ( switched on by fTNDCOR > 0 via TDIFFCOR CARDS line )
          DTDR = 0.1875 * OPARND * ( TEFF*(TEFF/TMID)**3. )
     >            - 0.75 * OPARND/STEBOLDPI 
     >                     * HTOTNDCOR / TMID**(3.)
        ENDIF
        
C***    Debug print
        IF (bPRINT .AND. (.NOT. bNoARAD) .AND. HTOTNDCOR /= 0.) THEN
           DTDRplus = 0.1875 * OPARND * ( TEFF*(TEFF/TMID)**3. )
           DTDRminus = - 0.75 * OPARND/STEBOLDPI 
     >                     * HTOTNDCOR / TMID**(3.)
           WRITE (hCPR,'(A,4(G15.5),2X,L1)') 'DTDR: ', 
     >        DTDRplus, DTDRminus, DTDR, 
     >        (T(L-1) + (RADIUS(L-1)-RADIUS(L))*DTDR ) / 1000., fTNDCOR
        ENDIF
        
        TLAST = TL
        
        TL = T(ND-1) + (RADIUS(ND-1)-RADIUS(ND))*DTDR
        ICOUNT = ICOUNT + 1

C***    Issue a warning if more than 10 iterations are needed
        IF (ICOUNT > 10) THEN
            WRITE (hCPR,'(3A,E12.5)') 
     >         'Troubles when determining Temperature at inner', 
     >         ' boundary by diffusion approximation', 
     >         ' ABS(TLAST/T(L)-1.) = ', ABS(TLAST/TL-1.)
        ENDIF
        
C***    New iteration if sufficient accuracy is not reached
C***     and we have not done more than 100 iterations
        IF (ABS(TLAST/TL-1.) < 1.E-6 .OR. ICOUNT > 100) EXIT
        
      ENDDO

C***  optional damping of the change in T(ND)
      IF (fTNDCOR > 0.) THEN
        TL = TSAVE + (TL - TSAVE) * fTNDCOR
      ENDIF
      
C***  Save new (converged) temperature
      T(ND) = TL
      
C***  Calculate relative correction to T(ND) for PRICORR routine
      TNDCORR = (TL-TSAVE) / TL
      
      RETURN
      END
