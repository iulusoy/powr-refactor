      SUBROUTINE LINTRIDIAGSOL (A,B,C,Q,DIAG,N,bALOTri,DIU,DIL)
C***********************************************************************
C***  TRIDIAGONALE MATRIX   -A(L)   B(L)  -C(L)
C***  RECHTE SEITE Q(L)
C***  i.e.  -A(L) J(L-1) + B(L) J(L) - C(L) J(L+1) = Q(L)
C***  LOESUNG AUF VEKTOR Q(L)
C***  DIAGONALE der INVERSEN der TRIDIAGONALEN MATRIX in vector DIAG
C***  Optional: Nebendiagonalen in DIU (upper) and DIL (lower)
C***
C***
C***  see Rybicki & Hummer (1991), Appendices A & B
C***
C***  called from COLIMO
C***********************************************************************
      IMPLICIT NONE

      INTEGER, PARAMETER :: NMAX = 200
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: bALOTri

      REAL, DIMENSION(N) :: A, B, C, Q, DIAG, DIU, DIL
      REAL, DIMENSION(NMAX) :: D, Z, E, F, H
      
      REAL :: DZNENNINV, FCI
      INTEGER :: I

C***  File and channel handles      
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)      
      
C***  Hard-coded switch for branching between optimized and non-optimized version      
c      LOGICAL, PARAMETER :: bOPTIMIZE = .TRUE.
      LOGICAL, PARAMETER :: bOPTIMIZE = .FALSE.

      IF (N > NMAX) THEN
        WRITE (hCPR,'(A)') 'LINTRIDIAGSOL: DIMENSION NMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I7))') ' NMAX= ', NMAX, ' N= ', N
        STOP 'FATAL ERROR in LINTRIDIAGSOL'
      ENDIF

      
C***  Branch between optimized and non-optimized version
C***  optimized version is numerically stable even for very small 
C***  optical depth increasements (see Rybicki & Hummer 1991, Appendix A)
      
C***  ------------------------------------------------------------------        
      IF (bOPTIMIZE) THEN
C***  optimized version for better delta-tau handling
C***  TODO: This version seems to have issues with div/0 (where?)
        H(1) = B(1) - C(1)
        F(1) = H(1) / C(1)
        Z(1) = Q(1) / B(1)
        DIAG(1) = 1./ B(1)
        DO I=2, N
          H(I) = - A(I) + B(I) - C(I)
          FCI  = (H(I) + A(I)*F(I-1)/(1.+F(I-1)))
          Z(I) = (Q(I) + A(I)*Z(I-1)) / (1.+FCI)        
          DIAG(I) = 1./(H(I)+A(I)+C(I)-A(I)/(1.+F(I-1)))
          IF (I < N) THEN
C***        Since C(N) == 0 we have to prevent division by zero          
C***        F(N) is never needed anyhow
            F(I) = FCI / C(I)
          ENDIF
        ENDDO
        E(N)=A(N)/B(N)
C***    Reuse Q for solution vector        
        Q(N)=Z(N)
C***    Backward elimination 
        DO I=N-1, 1, -1
          E(I) = A(I)/(B(I)-C(I)*E(I+1))
          DIAG(I) = DIAG(I)/(1.-E(I+1)/(1.+F(I)))
C***      Reuse Q for solution vector        
          Q(I)=Q(I+1)/(1.+F(I)) + Z(I)
        ENDDO      

C***  ------------------------------------------------------------------        
      ELSE
C***  standard version (non-optimized)            
        D(1)=C(1)/B(1)
        Z(1)=Q(1)/B(1)
        DIAG(1) = 1./B(1)
C***    Forward elimination      
        DO I=2, N
          DZNENNINV=1./(B(I)-A(I)*D(I-1))
          D(I) = C(I)*DZNENNINV
          Z(I) = (Q(I)+A(I)*Z(I-1))*DZNENNINV
          DIAG(I) = DZNENNINV
        ENDDO
        E(N)=A(N)/B(N)
C***    Reuse Q for solution vector        
        Q(N)=Z(N)
C***    Backward elimination 
        DO I=N-1, 1, -1
          E(I) = A(I)/(B(I)-C(I)*E(I+1))
          DIAG(I) = DIAG(I)/(1.-D(I)*E(I+1))
          IF (bALOTri) THEN
C***        Upper next-diagnonal entries
            DIU(I) = D(I)*DIAG(I+1)
C***        Lower next-diagonal entries
            DIL(I+1) = E(I+1)*DIAG(I)
          ENDIF
C***      Reuse Q for solution vector        
          Q(I)=D(I)*Q(I+1) + Z(I)
        ENDDO
      ENDIF
C***  ------------------------------------------------------------------        
      
      RETURN

      END
