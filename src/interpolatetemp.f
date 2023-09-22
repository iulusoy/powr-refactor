      SUBROUTINE INTERPOLATETEMP(Tnew, Told, Rnew, Rold, 
     >                           ENTOTnew, ENTOTold, 
     >                           TAUold, TEFF, ND, bUseENTOT)
C***********************************************************************
C***  Interpolation of the electron temperature on a new grid
C***   note: TAUold is the FULL TAUROSS incl. lines (from COLI)
C***  
C***  called by ENSURETAUMAX
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: TEFF
      REAL, DIMENSION(ND), INTENT(INOUT) :: Tnew
      REAL, DIMENSION(ND), INTENT(IN) :: Rnew, Rold, Told,
     >                                   ENTOTold, ENTOTnew, TAUold

      LOGICAL, INTENT(IN) :: bUseENTOT
      
      INTEGER, PARAMETER :: NDextrap = 2
      REAL, DIMENSION(NDextrap) :: Textrap, TextrapENTOT
      REAL, DIMENSION(ND) :: ENTOToldLOG, ENTOTnewLOG, TAUnew, RoldLOG
      
      REAL :: TMold, TMnew, DTDR, Thelp, Rhelp, 
     >        RinOLD, TinOLD, RnewLOGL
      INTEGER :: L, NDr
      LOGICAL :: bTINEXTRAP, bTAUEXTRAP
      
      bTINEXTRAP = .FALSE.      !switch on/off inner temperature extrapolation
      bTAUEXTRAP = .FALSE.     
      
      IF (bUseENTOT) THEN

        DO L=1, ND
          ENTOToldLOG(L) = LOG10(ENTOTold(L))
          ENTOTnewLOG(L) = LOG10(ENTOTnew(L))
        ENDDO          
        
C       Estimate current Tau change         
        TAUnew(1) = 0.
        DO L=2, ND
          TAUnew(L) = TAUnew(L-1) + ENTOTnew(L)/ENTOTold(L)
     >         * (Rnew(L-1) - Rnew(L))/(Rold(L-1) - Rold(L))
     >         * (TAUold(L) - TAUold(L-1))          
        ENDDO
        TMnew = TAUnew(ND)
        
        IF (ENTOTnewLOG(ND) > ENTOToldLOG(ND)) THEN          
          CALL SPLINPOX(TMold, ENTOToldLOG(ND), TAUnew, ENTOTnewLOG, ND)
          TinOLD = Told(ND)          
        ENDIF
        
C       Perform interpolation
        dploop: DO L=1, ND
          IF (ENTOTnewLOG(L) > ENTOToldLOG(ND)) THEN
            IF (bTINEXTRAP) THEN
              Tnew(L) = TinOLD *
     >                ( (TAUnew(L) + 2./3.) / (TMold + 2./3.) )**(0.25)
            ELSE
              Tnew(L) = Told(ND)
            ENDIF
          ELSEIF (ENTOTnewLOG(L) < ENTOToldLOG(1)) THEN
            !less dense than old outermost value => take old outer boundary value
            Tnew(L) = Told(1)
          ELSE
            CALL SPLINPOX(Tnew(L), ENTOTnewLOG(L),
     >                    Told, ENTOToldLOG, ND)
          ENDIF
        ENDDO dploop

      ELSEIF (.TRUE.) THEN
C***    Interpolation over log radius  
        DO L=1, ND
          IF (L > 1) THEN
C***    Account for the possibility that the old radius might be shifted
C***    and thus we have to determine the number of points NDr in the "normal"
C***    radius Range 1...RMAX. 
C***    (One point for R < 1 is allowed to avoid inner cutoffs.)
            IF (Rold(L-1) < 1. .OR. Rold(L) <= 0.) EXIT
          ENDIF
          RoldLOG(L) = LOG10(Rold(L))
          NDr = L
        ENDDO                 
        DO L=1, ND
          RnewLOGL = LOG10(Rnew(L))
          IF (RnewLOGL > RoldLOG(1)) THEN
            Tnew(L) = Told(1)
          ELSEIF (RnewLOGL < RoldLOG(NDr)) THEN
            Tnew(L) = Told(NDr) 
          ELSE
            CALL SPLINPOX(Tnew(L), RnewLOGL, Told, RoldLOG, NDr)
          ENDIF
        ENDDO
      ELSE
C***    Interpolation over radius      
        DO L=1, ND
          IF (Rnew(L) > Rold(1)) THEN
            Tnew(L) = Told(1)
          ELSEIF (Rnew(L) < Rold(ND)) THEN
            Tnew(L) = Told(ND) 
          ELSE
            CALL SPLINPOX(Tnew(L), Rnew(L), Told, Rold, ND)
          ENDIF
        ENDDO
      ENDIF

      RETURN

      END
