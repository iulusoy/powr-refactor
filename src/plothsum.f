      SUBROUTINE PLOTHSUM (HTOTL, HTOTND, EDDIHOUTJMEAN, XJTOTL, HTOTM,
     >                     HTOTG, HTOTOBS, HTOTCMF0, HTOTCMF0ADV, 
     >                     FLUXEPS, ND, MODHEAD, JOBNUM, KANAL, TEFF, 
     >                     BHTOTERR, bINCADV)
C******************************************************************************
C***  DIRECT TRANSFER OF HSUM PLOT
C***  TOTAL (FREQUENCY-INTEGRATED) FLUX versus DEPTH INDEX
C***  Note: there are some more curves available, but de-activated
C***        wrh 14-Apr-2003 11:58:53
C******************************************************************************
    
      IMPLICIT NONE

      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, KANAL, JOBNUM
      REAL, INTENT(IN) :: TEFF

      CHARACTER(100) :: MODHEAD 
      CHARACTER(110) :: HEADLINE
      CHARACTER(8) :: CENTER

      REAL, DIMENSION(NDMAX) :: X, X1, Y, YEPS
      REAL, DIMENSION(ND) :: HTOTL, HTOTM, HTOTG, HTOTOBS, HTOTCMF0, 
     >                       HTOTCMF0ADV, XJTOTL
      REAL, INTENT(IN) :: HTOTND, EDDIHOUTJMEAN, FLUXEPS
      LOGICAL :: BHTOTERR, bINCADV

      INTEGER :: L
      REAL :: XMIN, XMAX, YMIN, YMAX, TRAD, HSUM, HTOTOUT

C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5


      IF (BHTOTERR) THEN
         WRITE (*, '(A)')
     >         'PLOTHSUM: WARNING: HTOTL  not present in MODEL file'
     >         //' --  PLOT HTOT disabled!'
         RETURN
      ENDIF 


      IF (ND .GT. NDMAX) THEN
        WRITE (0, '(A)') 'DIMENSION INSUFFICIENT - HTOT PLOT SUPPRESSED'       
        WRITE (0, '(A)') 'NON-FATAL ERROR IN SUBROUTINE PLOTHSUM'
        RETURN
      ENDIF
 
      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
      CALL JSYMSET ('G2','TRANSFER')
 
      HEADLINE = 'HSUM: M'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM

      WRITE (KANAL, '(A)') '*NEXTPLOT: PLOT HTOT'
      WRITE (KANAL, '(A)') 'PLOT: ' // HEADLINE
      WRITE (KANAL, '(A)') '\FONT=HELVET'      
      WRITE (KANAL, '(A)') '\DEFINECOLOR 9 0.6 0.60 0.6'
      WRITE (KANAL, '(A)') '\PEN=1'     !needed to get rid of definecolor pen bug 
      
C***  HTOTOBS: Conversion into radiation temperatures
      DO L=1, ND-1
         X(L) = FLOAT(L)
         HSUM = HTOTOBS(L)
         TRAD = (4. * ABS(HSUM) / STEBOLDPI)**0.25
         IF (HSUM .LT. .0) TRAD = -TRAD
         Y(L) = TRAD / 1000.
      ENDDO

cc      CALL PLOTANF (KANAL,HEADLINE, '&E'//HEADLINE,
cc     $        CENTER//'DEPTH INDEX L',
cc     $        CENTER//'T&Trad&M / KK',
cc     >             0., XMIN, XMAX, 5., 10., 0.,
cc     >             0., YMIN, YMAX, .2,  1., 0.,
cc     $        X,Y,ND-1,5)

C***  HTOT+HMECH+HGRAV: Conversion into radiation temperatures
cc      DO L=1, ND-1
cc         HSUM = HTOTOBS(L)+HTOTM(L)+HTOTG(L)
cc         TRAD = (4. * ABS(HSUM) / STEBOL)**0.25
cc         IF (HSUM .LT. .0) TRAD = -TRAD
cc         Y(L) = TRAD / 1000.
cc      ENDDO
cc      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 3') 

C***  HTOT+HGRAV: Conversion into radiation temperatures
cc      DO L=1, ND-1
cc         HSUM = HTOTOBS(L)+HTOTG(L)
cc         TRAD = (4. * ABS(HSUM) / STEBOL)**0.25
cc         IF (HSUM .LT. .0) TRAD = -TRAD
cc         Y(L) = TRAD / 1000.
cc      ENDDO
cc      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 3') 

C***  HTOTCMF0: Conversion into radiation temperatures
      DO L=1, ND
        IF (L /= ND) THEN
          X1(L) = FLOAT(L)+0.5
          Y(L) = HTOTCMF0(L) 
        ELSE
          X1(L) = FLOAT(ND)
          Y(L) = HTOTCMF0(ND-1) 
        ENDIF
        TRAD = (4. * ABS(Y(L)) / STEBOLDPI)**0.25
        IF (Y(L) .LE. .0) TRAD = -TRAD
        Y(L) = TRAD / 1000.
      ENDDO
cc      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 1 PEN = 3') 
      CALL PLOTANFS (KANAL,HEADLINE, '&E'//HEADLINE,
     $        CENTER//'DEPTH INDEX L',
     $        CENTER//'T&Trad&M / KK',
     >             0., 0., 0., 0., 0., 0.,
     >             0., 0., 0., 0., 0., 0.,
     $        X1,Y,ND, 'PEN=8')

C***  HTOTCMF0 +- FLUXEPS in TRAD      
      IF (FLUXEPS > 0.) THEN
        DO L=1, ND
          YEPS(L) = ABS(Y(L)) * (1. + FLUXEPS)**0.25
        ENDDO
        CALL PLOTCONS (KANAL,X1,YEPS,ND,'COLOR=9 PEN=1') 
        DO L=1, ND
          YEPS(L) = ABS(Y(L)) * (1. - FLUXEPS)**0.25
        ENDDO
        CALL PLOTCONS (KANAL,X1,YEPS,ND,'COLOR=9 PEN=1') 
      ENDIF

C***  HTOTCMF0 (ADV) in TRAD
cc      DO L=1, ND
cc        IF (L /= ND) THEN
cc          X1(L) = FLOAT(L)+0.5
cc          Y(L) = HTOTCMF0ADV(L) 
cc        ELSE
cc          X1(L) = FLOAT(ND)
cc          Y(L) = HTOTCMF0ADV(ND-1) 
cc        ENDIF
cc        TRAD = (4. * ABS(Y(L)) / STEBOLDPI)**0.25
cc        IF (Y(L) .LE. .0) TRAD = -TRAD
cc        Y(L) = TRAD / 1000.
cc      ENDDO

C***  HTOTCMF0ADV in TRAD      
cc      CALL PLOTCONS (KANAL,X1,Y,ND,'COLOR=4 PEN = 3') 
     
C***  HTOTL: Conversion into radiation temperatures
      HTOTOUT = EDDIHOUTJMEAN * XJTOTL(1)
c      WRITE (0,*) 'HTOTOUT = ', HTOTOUT, HTOTL(1)
      DO L=0, ND
        IF (L == 0) THEN
          X(1) = FLOAT(1)
          Y(1) = HTOTOUT
        ELSEIF (L == ND) THEN
          X(ND+1) = FLOAT(ND)
          Y(ND+1) = HTOTND
        ELSE
          X(L+1) = 0.5 * ( FLOAT(L) + FLOAT(L+1) )
          Y(L+1) = HTOTL(L) 
        ENDIF
        TRAD = (4. * ABS(Y(L+1)) / STEBOLDPI)**0.25
        IF (Y(L+1) .LE. .0) TRAD = -TRAD
        Y(L+1) = TRAD / 1000.
      ENDDO
            
C***  HTOTL in TRAD      
      CALL PLOTCONS (KANAL,X,Y,ND+1,'COLOR= 2 PEN = 3') 

C***  Continuum
      X(1) = 1.
      X(2) = FLOAT(ND-1)
      Y(1) = TEFF / 1000.
      Y(2) = TEFF / 1000. 
      CALL PLOTCONS (KANAL,X,Y,2, 'COLOR=3') 

C***  Add an invisible dataset to span TEFF +- 5000
      Y(1) = TEFF/1000. + 5
      Y(2) = TEFF/1000. - 5
      CALL PLOTCON (KANAL,X,Y,2,0)

      RETURN
      END
