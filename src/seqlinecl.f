      SUBROUTINE SEQLINECL(NLINE, LINE, EINST, INDLOW, 
     >             INDNUP, XLAMSOR, ELEVEL,
     $             NDIM, VDOP, CMFBAND, CLIGHT, XLAMMIN, XLAMMAX,
     $             VELO, EXLAM1, EXLAM2, MAXEXT, LASTIND, NAUTO, KRUDAUT)
C**********************************************************************
C***  "SEQUENCE OF LINES"  
C***  THE LINES ARE SORTED IN THE SEQUENCE OF INCREASING WAVELENGTHS.
C***  RUDIMENTAL LINES ARE OMITTED.
C***  THE BOUNDARIES OF THE CMF BANDS ARE DETERMINED (XLAMMIN, XLAMMAX)
C***  - ATTENTION! CALLED FROM TWO DIFFERENT PLACES: 
C***    COLI 
C***    STEAL - LINPOP - COMA - SETXJFINE
C**********************************************************************

      DIMENSION EINST(NDIM,NDIM),INDLOW(NDIM),INDNUP(NDIM),ELEVEL(NDIM)
      DIMENSION XLAMSOR(NLINE),XLAMMIN(NLINE),XLAMMAX(NLINE),LINE(NLINE)
      DIMENSION KRUDAUT(NAUTO)
      LOGICAL BRUDTEST

C***  Count of  rudimental lines which are removed from the LINE vector
      NUMRUD=0

      DO NL = 1, LASTIND + NAUTO
         IND = LINE(NL)
         XLAM = 1.E8 / (ELEVEL(INDNUP(IND))-ELEVEL(INDLOW(IND)))
         IF (IND .LE. LASTIND) THEN
            BRUDTEST = (EINST(INDLOW(IND),INDNUP(IND)) .EQ. -2.) 
         ELSE
            BRUDTEST = (KRUDAUT(IND-LASTIND) .EQ. 1)
         ENDIF
         IF (BRUDTEST) THEN
            NUMRUD=NUMRUD+1
            XLAMSOR(NL)=-1.0
         ELSE
            XLAMSOR(NL) = XLAM
         ENDIF
      ENDDO

C**   Sort by Wavelength / Key to Sorted Field is Line(..)
      CALL RSORT(NLINE,XLAMSOR,LINE)

C**   Omit rudimentals from Line-list
      DO NL=1,NLINE-NUMRUD
         XLAMSOR(NL)=XLAMSOR(NL+NUMRUD)
         LINE (NL)=LINE (NL+NUMRUD)
      ENDDO

      NLINE = LASTIND + NAUTO - NUMRUD

C**   Calculate red and blue border for each Line
      DO 400 NL=1,NLINE

C***  "EXTEND" OPTION: 
C***  FOR LINES WHICH ARE FALLING IN THE SPECIFIED WAVELENGTH RANGE, 
C***  THE RED WING IS EXTENDED TO  2 * V-FINAL
      REDBAND = CMFBAND
      XLAMACT = XLAMSOR(NL)
      IF ( (XLAMACT - EXLAM1) * (XLAMACT - EXLAM2) .LT. .0 )
     >   REDBAND = CMFBAND + 2. * VELO / VDOP

      XLAMMIN(NL) = XLAMSOR(NL) * (1 - CMFBAND * VDOP / CLIGHT)
      XLAMMAX(NL) = XLAMSOR(NL) * (1 + REDBAND * VDOP / CLIGHT)
400   CONTINUE

      RETURN
      END
