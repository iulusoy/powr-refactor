      FUNCTION ISAMAX (N, X, INC)
C***  Find index of vector element with highest absolute value
      DIMENSION X(N)

ccc   The following call works with our intel compiler
ccc   with gfortran one might need to add the library -lbas as option

      ISAMAX = IDAMAX (N, X, INC)

C***  If IDAMAX is still not found, one can comment the above call
C***  and activete the following hand-made version
c      XMAX = ABS(X(1))
c      IMAX = 1
c
c      DO I=2*INC, N, INC
c        IF (ABS(X(I)) .GT. XMAX) THEN
c          XMAX = ABS(X(I))
c          IMAX = I
c        ENDIF
c      ENDDO
c
c      ISAMAX = IMAX

      RETURN
      END
