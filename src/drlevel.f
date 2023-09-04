      SUBROUTINE DRLEVEL (N, NDIM, MAXIND, MAXAUTO, NAUTO, KRUDAUT, 
     $                    LOWAUTO, IONGRND, AAUTO, ELEVEL, 
     $                    LEVEL, EINST, EION, WEIGHT, INDLOW, 
     $                    INDNUP, LASTIND, 
     $                    ND, T, ENTOT, RNE, POPNUM, DENSCON, 
     >                    LEVUPAUTO, LEVAUTO, WAUTO, N_WITH_DRLEVELS)
C******************************************************************************
C***  POPNUM array is augmented by additional columns for the 
C***  auto-ionizing levels, assuming LTE relative to the next-upper
C***  ground level.
C***  Einstein ceofficients are inserted in the EINST array at the
C***  according positions (-2. for rudimental lines).
C***  After DRLEVEL, the DRTRANSIT lines are treated in the COLI
C**   completely analogues to "normal" lines. 
C***  The new number of levels (=N_WITH_DRLEVELS) never appears explicitely 
C***  in COLI elsewhere, but indirectly since the INDNUP entries can point
C***  to the appended level indices. 
C***  THIS SUBROUTINE MAY ONLY BE CALLED FROM COLI !!
C******************************************************************************

      DIMENSION WEIGHT(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),EINST(NDIM,NDIM)
      DIMENSION INDNUP(MAXIND),INDLOW(MAXIND), IONGRND(MAXIND)
      DIMENSION LOWAUTO(MAXAUTO), WAUTO(MAXAUTO)
      DIMENSION AAUTO(MAXAUTO), KRUDAUT(MAXAUTO)
      DIMENSION T(ND),ENTOT(ND),RNE(ND),POPNUM(ND,NDIM)
      CHARACTER*10 LEVEL(NDIM)
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)
      DIMENSION DENSCON(ND)

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

C***  Calculate and append LTE popnumbers for the AUTO levels
      DO L = 1, ND
         TL = T(L)
         SQRTL = SQRT(TL)
         TL32 = TL * SQRTL
C***     Density increased by DENSCON
         ENTOTL = ENTOT(L)*DENSCON(L)
         RNEL = RNE(L)
         DO ILEVEL = N+1, N_WITH_DRLEVELS
            INDION = IONGRND(ILEVEL)            
            POPION = POPNUM(L,INDION)
            BOLTZ = EXP(-C1 * (ELEVEL(ILEVEL)-EION(ILEVEL)) / TL)
            SAHAPHI = WEIGHT(ILEVEL) / WEIGHT(INDION) * CI / TL32 * BOLTZ 
            POPNUM(L,ILEVEL) = POPION * ENTOTL * RNEL * SAHAPHI

ccc   test output
cc           if (L .eq. 30) write (0,'(3i6,3G14.6)')
cc     >       ilevel, indion, int(weight(indion)), boltz
cc     >       , sahaphi, POPNUM(L,ILEVEL)/POPION 

         ENDDO
      ENDDO 

C***  Append  rows and columns to EINST array with A_up_low
C***  PRESET FOR NEW EINST-VALUES  = -11.   -- Why? 
***********************************************
      DO NUP = N+1, N_WITH_DRLEVELS
         DO LOW = N+1, N_WITH_DRLEVELS
            EINST(NUP,LOW) = -11.
         ENDDO
      ENDDO

      DO INDDR = 1, NAUTO
         IND = LASTIND + INDDR
         NUP = INDNUP(IND)
         LOW = INDLOW(IND)
         EINST(NUP,LOW) = AAUTO(INDDR)
         IF (KRUDAUT(INDDR) .EQ. 1) EINST(LOW,NUP) = -2.
      ENDDO

      RETURN

C***  ERROR BRANCH ********************************************
   99 WRITE (0,'(A,2I4)') '*** You must DRTRANSIT data ' //
     >                    'ion the corrected version (after June 2023)'
      STOP '*** ERROR STOP in subr. DRLEVEL'

      END 
