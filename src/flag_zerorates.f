      SUBROUTINE FLAG_ZERORATES(NFIRNA, NLANA, RATCO, N, NRANK, 
     >                          IMAXPOP, EN, POPMIN, ZERO_RATES)
C******************************************************************
C***  Check if all Rate Coefficients in one column are non-zero
C***  (otherwise: the matrix is singular!)
C***  and store logical flag ZERO_RATES for later use
C***  Note: This subroutine is called once per each element
C***  Called from:
C***  STEALCL -  LINPOP - COMA - FLAG_ZERORATES
***   STEALCL -  NLTEPOP - FLAG_ZERORATES
C******************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N, NRANK,
     >                       NFIRNA, NLANA, IMAXPOP
      REAL, INTENT(IN) :: POPMIN
      
      REAL, DIMENSION(NRANK,NRANK) :: RATCO
      REAL, DIMENSION(NRANK) :: EN
      LOGICAL, DIMENSION(N) :: ZERO_RATES
      CHARACTER(10), DIMENSION(N) :: LEVEL
      
      REAL :: TEXC, POPLIMIT, GAINS, ENJ
      INTEGER :: I, J

C***  Threshold value: if estimated popnumber is lower, the level is flagged
      POPLIMIT = MAX(1.E-20 * POPMIN, 1.E-45)

C***  check colum J

C***  test modification, wrh 29-Aug-2008 16:23:38:
C***  Only the hightest levels can be switched off!
C***  I.e., for a level to be switched off, 
C***  the next-higher level must be already switched off 
C***  For this criterion, the loop over levels must run backwards:
      DO J = NLANA, NFIRNA, -1
        
C***     Flag if diagonal element too small, else divison below will fail!
         IF (ABS(RATCO(J,J)) .LT. 1.E-200) THEN
            ZERO_RATES(J) = .TRUE.
            CYCLE
         ENDIF

C***     Sum up the gains of level J
         GAINS = .0
         DO I = NFIRNA, NLANA
            IF (I .EQ. J) CYCLE
            IF (EN(I) > 1.1 * POPMIN) THEN
              GAINS = GAINS + EN(I) * RATCO(I,J)
            ENDIF
         ENDDO
         
C***     expected popnumber of level J
           ENJ = - GAINS / RATCO(J,J)
C***     mysteriously usefull??
           ENJ = ABS(ENJ)
         
C***     flag level, depending on predicted and on last popnumber
C***        the second condition provides some hysteresis
         ZERO_RATES(J) = ENJ .LT. POPLIMIT
     >       .OR. (EN(J) .LE. 1.1*POPMIN .AND. ENJ .LE. 100.*POPLIMIT)
        
C***     Never flag the IMAXPOP level!
C***     Never flag a level which was strongly populated last iteration
           IF (J == IMAXPOP .OR. EN(J) > 1000*POPMIN) THEN
              ZERO_RATES(J) = .FALSE.
              CYCLE
           ENDIF

C***  This version: Only the highest levels can be flagged
        IF (J .EQ. NLANA) CYCLE
         
        IF (.NOT. ZERO_RATES(J+1) ) THEN
          ZERO_RATES(J)=.FALSE. 
        ENDIF

      ENDDO

      RETURN
      END
