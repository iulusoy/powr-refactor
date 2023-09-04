      SUBROUTINE CLSAVEE (NCHANE, NZE1, NZE2, NFRO, EDDIA, NDEDDIA, 
     >                    EDDIF, EDDIG, ND, EDDIHOUT, 
     >                    EDDIHIN, EDDIHINM, EDDINOUT, EDDININ, 
     >                    BCLEERR, CMODE, XHI, XHO, EPSG, 
     >                    XHOM, XNOM, EDDIHOUTP, EDDINOUTP, EDDIHINT)
C**********************************************************************
C***  Saves the old EDDIEs to file fort.<NCHANE>
C**********************************************************************

      DIMENSION EDDIA(NDEDDIA, NFRO)
      DIMENSION EDDIF(ND), EDDIG(ND-1)
      DIMENSION EPSG(ND)
      CHARACTER NAME*8, CMODE*1
      LOGICAL BCLEERR
      INTEGER :: IDUMMY
      REAL :: DBDTAUND, EDDIHINT

C***  Store the EDDIEs into EDDIA
C***  At this point the change of the EDDIEs could be calculated
C***  For that BCLEERR would be needed
      NDA = 2*ND - 1
      DO L=1, ND
        EDDIA(L, NZE2) = EDDIF(L)
      ENDDO
      DO L=1, ND-1
        EDDIA(ND+L, NZE2) = EDDIG(L)
      ENDDO
      EDDIA(NDA+1, NZE2) = EDDIHOUT
      EDDIA(NDA+2, NZE2) = EDDIHIN
      EDDIA(NDA+3, NZE2) = EDDINOUT
      EDDIA(NDA+4, NZE2) = EDDININ
      EDDIA(NDA+5, NZE2) = XHI
      EDDIA(NDA+6, NZE2) = XHO
      EDDIA(NDA+7, NZE2) = XHOM
      EDDIA(NDA+8, NZE2) = XNOM
      EDDIA(NDA+9, NZE2) = EDDIHOUTP
      EDDIA(NDA+10, NZE2) = EDDINOUTP
      EDDIA(NDA+11, NZE2) = EDDIHINM
      EDDIA(NDA+12, NZE2) = EDDIHINT

C***  Save array to file
      IF (NZE2 .EQ. NFRO .OR. CMODE .EQ. 'F') THEN
        WRITE (NAME,'(A2,I5,A1)') 'ED', NZE1, ' '
        CALL WRITMS (NCHANE, EDDIA, NDEDDIA*NFRO, NAME,-1,IDUMMY, IERR)
      ENDIF
      IF (CMODE .EQ. 'F') THEN
        CALL WRITMS (NCHANE, EPSG, ND-1, 'EPSG    ',  -1, IDUMMY, IERR)
      ENDIF      

      RETURN
      END
