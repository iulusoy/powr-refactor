      SUBROUTINE VSUB(A,B,N)
C***********************************************************************
C***  C***********************************************************************
      DIMENSION A(N),B(N)

      DO 1 I=1,N
    1 A(I)=A(I)-B(I)

      RETURN
      END