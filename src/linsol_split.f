      SUBROUTINE LINSOL_SPLIT (X, A, B, N, NDIM,
     >                         NFIRST, NLAST, NATOM, NBLOCK, AOFF,
     >                         SCRATCH, ATEST, BTEST, DTEST, VERSION, 
     >                         NBSTART, NBEND, NBATOM, NBFirstIon, 
     >                         IMAXPOP, ZERO_RATES, LEVEL, 
     >                         iBLOCKINVERSION)
C**********************************************************************
C***  SOLVES THE LINEAR SYSTEM:  >>  X * A = B
C***    by using the atomic block structure of A
C***    (this means of course that A must be block-diagonal)
C***  NDIM  = ROW DIMENSION  OF A
C***  N     = RANK OF THE SYSTEM
C***  A     = COEFFICIENT MATRIX 
C***  B     = RIGHT=HAND SIDE VECTOR
C***  X     = UNKNOWN VECTOR  --> SOLUTION
C***  Call from Subroutine LINPOP
C**********************************************************************
C     
C     Notation info:
C     i = row index, j = column index
C     When using blocks, the indexes i and j ares used for the index inside the blocks only.
C     indBlockI and indBlockJ are used for referencing a block inside the large matrix
C      
C**********************************************************************
C     This version is based on the ideas from the subroutine 
C      linsol_split in the Goetz branch. However, there are certain
C      things in Goetz code that are not in the wrh version. These are:

C     BFAST - Logical to switch between two methods, today always
C        set to BFAST = .TRUE. because BFAST = .FALSE. simply means
C        that no block inversion is used
C
C     B_SKIP_POP(N) - Array of Logicals
C     NSTART(100), NSTOP(100)
C        Allows to skip certain levels, therefore NSTART and NSTOP
C        had to be introduced to store the "reduced" block start
C        This version works like all entries of B_SKIP_POP are .FALSE.
C        which is equal to 
C      NSTART(NA)= NFIRST(NA)
C      NSTOP(NA) = NLAST(NA)
C        for all entries in the Goetz branch. However, when it comes
C        to the matrix inversion part, NSTART and NSTOP are not used,
C        so this stuff might have been discarded

C**********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, N, NATOM, NBLOCK, iBLOCKINVERSION
      REAL, DIMENSION(NDIM,NDIM), INTENT(INOUT) :: A, AOFF
      REAL, DIMENSION(NDIM), INTENT(IN) :: B
      REAL, DIMENSION(NDIM), INTENT(OUT) :: X
      
C*** As we do not know the size of the blocks we init SCRATCH as a vector
C     (Note: SCRATCH is used for one block in the matrix, this can be
C            up to the dimension of the original matrix)
      REAL, DIMENSION(NDIM*NDIM), INTENT(INOUT) :: SCRATCH

      REAL, DIMENSION(NDIM,NDIM), INTENT(OUT) :: ATEST
      REAL, DIMENSION(NDIM), INTENT(OUT) :: BTEST, DTEST
      REAL, DIMENSION(NDIM) :: ROWINVIONIMAX, JBSTART, JBEND 
      INTEGER, DIMENSION(NDIM), INTENT(IN) :: NBSTART, NBEND, NBATOM
      INTEGER, DIMENSION(NATOM) :: NBPIVOT
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST,
     >                                         IMAXPOP, NBFirstIon
      CHARACTER(4), INTENT(INOUT) :: VERSION
      CHARACTER(10), DIMENSION(N), INTENT(IN) :: LEVEL
      LOGICAL, DIMENSION(NDIM), INTENT(IN) :: ZERO_RATES

      INTEGER :: IMPMAX, na, nn, i, j, nBlockDim, nPivotDim, iBLOCK,
     >           indBlockI, indBlockJ, nFirstLev, nLastLev,
     >           nVectorIndex, IMPRO, iPrimat
      REAL :: DIFF, SUMINVi

      LOGICAL :: bForceBlock, bPrintRightVector, bPIVOT

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      
      
C***  set iPrimat to zero to print the inverted matrix structure
      iPrimat = 1

      bForceBlock = .TRUE.          !Force block matrix structure 
      bPrintRightVector = .FALSE.   !Print B, BTEST and difference (debug option)


C***  DEFINE THE INVERSION SUBROUTINE: 'FREE' OR 'OWN' OR 'OWNL'
C***  The difference between OWN and OWNL is, that in case of a singularity
C***    the Subroutine OWNINV stops for VERSION = OWN while it returns
C***    SING in the other case. LINSOL can then skip this depth point and 
C***    continues with a MODIFY.
      VERSION = 'OWNL'

C***  NUMBER OF IMPROVEMENT-ITERATIONS
      IMPMAX = 1

      IF (N > NDIM) THEN
         PRINT *,'ERROR IN SUBROUTINE LINSOL: N .GT. NDIM'
         STOP 'ERROR'
      ENDIF
      
C***  SAVE MATRIX A BEFORE INVERSION 
C***  (OPTIONAL: ENSURE BLOCK STRUCTURE BY forceblock = .TRUE.)
C***  08.06.2010 forceblock MUST be used to avoid SUM-Error at some depth points
C     Side note: Only if GAMMA=0 has been set in the CARDS file
C      the matrix DM is fully block-diagonal. In all other cases
C      there are (usually only few) non-zero elements outside of
C      the blocks. There forceblock = .TRUE. should always be used
C      as otherwise these elements would be simply copied (without
C      any inversion) to the inverted matrix and could cause
C      strange errors
      IF (bForceBlock) THEN
C***     init ATEST as identity matrix
         ATEST = 0.
C         DO nn=1, NLAST(NATOM)
         DO nn=1, NDIM
            ATEST(nn, nn) = 1.
         ENDDO
C***     copy atom/ion blocks
         DO iBLOCK=1, NBLOCK
            nFirstLev = NBSTART(iBLOCK)
            nLastLev = NBEND(iBLOCK)
            nBlockDim = nLastLev - nFirstLev + 1
            DO j=1, nBlockDim
               indBlockJ = nFirstLev + j - 1
               DO i=1, nBlockDim
                  indBlockI = nFirstLev + i - 1
                  ATEST(indBlockI, indBlockJ) = A(indBlockI, indBlockJ)
               ENDDO
            ENDDO
         ENDDO
C***     copy additional diagonal elements (if there are any)
         IF (NDIM .GT. NLAST(NATOM)) THEN
            DO nn=NLAST(NATOM)+1, NDIM
               ATEST(nn, nn) = A(nn, nn)
            ENDDO
         ENDIF

         AOFF = 0.
         DO i=1, NDIM
            DO j=1, NDIM
               AOFF(i,j) = A(i,j) - ATEST(i,j)    !non-block elements are stored in AOFF
               A(i,j) = ATEST(i,j)
            ENDDO
         ENDDO
         IF (NDIM > NLAST(NATOM)) THEN
           !remove all non-block lines fom offcenter matrix
           DO J=NLAST(NATOM)+1, NDIM
             DO I=1, NDIM     
               IF (I /= J) THEN
                 AOFF(I, J) = 0.
                 AOFF(J, I) = 0.
               ENDIF
             ENDDO
           ENDDO
         ENDIF

C***     Restore number conservation column for ion split
         IF (iBLOCKINVERSION == 2) THEN          
           DO na=1, NATOM
              DO i=NFIRST(na), NLAST(na)
C                IF (AOFF(i, IMAXPOP(na)) == 1.) THEN
                IF (AOFF(i, IMAXPOP(na)) == 1. .AND.
     >                 .NOT. ZERO_RATES(i)) THEN
                  !The innermost IF is necessary to prevent overwriting level switchoff effects
                  ATEST(i, IMAXPOP(na)) = 1.
                  AOFF(i, IMAXPOP(na)) = 0.
                ENDIF
              ENDDO
           ENDDO
         ENDIF

         
C***  PRINT Matrix debug option
         IF (iPrimat .eq. 0) THEN
            write (*,*) 'LINSOL_SPLIT> PRIMAT'
            CALL PRIMAT(ATEST, N, NDIM, 'Matrix A')
            write (*,*) 'LINSOL_SPLIT> end of PRIMAT'
            write (*,*) 'LINSOL_SPLIT> PRIMAT'
            CALL PRIMAT(AOFF, N, NDIM, 'Matrix AOFF')
            write (*,*) 'LINSOL_SPLIT> end of PRIMAT'
            iPrimat = 1
         ENDIF
C***  end of block structure test
         
      ELSE      
         DO i=1, N
            DO j=1, N
               ATEST(i,j) = A(i,j)
            ENDDO
         ENDDO
      ENDIF
      
      DO iBLOCK=1, NBLOCK
        JBSTART(iBLOCK) = NBSTART(iBLOCK)
        JBEND(iBLOCK) = NBEND(iBLOCK)
C        WRITE (6,*) 'START, END: ', i, JBSTART(i), JBEND(i)
      ENDDO
      NBPIVOT = 0.
      IF (iBLOCKINVERSION == 2) THEN
        !Re-sort blocks to have pivot block (= block with IMAXPOP) first for each atom
        DO iBLOCK=1, NBLOCK
          NA = NBATOM(iBLOCK)
          IF (IMAXPOP(NA) >= NBSTART(iBLOCK) 
     >               .AND. IMAXPOP(NA) <= NBEND(iBLOCK)) THEN
            !"PIVOT" block has the number conservation directly in the block
            ! => This block is the one that has to be inverted first
            NBPIVOT(NA) = iBLOCK
          ENDIF 
        ENDDO
        DO na=1, NATOM
          WRITE (6,*) 'PIVOT, FI: ', na, NBPIVOT(na), NBFirstIon(na)
          IF (NBPIVOT(na) /= NBFirstIon(na)) THEN
            !If first block does not contain number conservation 
            !   => switch with PIVOT block, so that PIVOT block is inverted first
            !get block start level from block index from NBFirstIon
            JBSTART(NBFirstIon(na)) = NBSTART(NBPIVOT(na))
            JBEND(NBFirstIon(na)) = NBEND(NBPIVOT(na))
            JBSTART(NBPIVOT(na)) = NBSTART(NBFirstIon(na))
            JBEND(NBPIVOT(na)) = NBEND(NBFirstIon(na))
          ENDIF
        ENDDO
      ENDIF

C***  MATRIX INVERSION -------------------------------
      DO iBLOCK=1, NBLOCK         
         SCRATCH = 0.
C***     extract blocks
         nFirstLev = JBSTART(iBLOCK)
         nLastLev = JBEND(iBLOCK)
         nBlockDim = nLastLev - nFirstLev + 1
         !First ion block in an atom is pivot block (ensured above)
         NA = NBATOM(iBLOCK)
         IF (iBLOCKINVERSION == 2) THEN
            bPIVOT = (NBFirstIon(NA) == iBLOCK)
            WRITE (6,'(A,2(I3,2X),L1,2(2X,I5))') 'BLOCK: ', iBLOCK, NA,
     >                bPIVOT, nFirstLev, nLastLev
            IF (bPIVOT) THEN
              ROWINVIONIMAX = 0.
              nPivotDim = nBlockDim     !save dimension of pivot block
            ENDIF
         ELSE
           bPIVOT = .FALSE.
         ENDIF
         DO j=1, nBlockDim
            indBlockJ = nFirstLev + j - 1
            DO i=1, nBlockDim
               indBlockI = nFirstLev + i - 1
C***           store column vectors of current block in one giant vector
C***             because the dimension of the blocks is unknown
C***             when SCRATCH has to be initialized
C***           subroutine INV will treat SCRATCH as a Matrix
               nVectorIndex = i + (j - 1) * nBlockDim
               SCRATCH(nVectorIndex) = A(indBlockI, indBlockJ)
            ENDDO
         ENDDO            

C***     call INV and cancel everything if singularity encounters
C***       => returns inverted block in SCRATCH
         CALL INV (nBlockDim, nBlockDim, SCRATCH, VERSION)

C***     replace original block with inverted block in the matrix
         DO j=1, nBlockDim
            indBlockJ = nFirstLev + j - 1
            DO i=1, nBlockDim
               indBlockI = nFirstLev + i - 1
               nVectorIndex = i + (j - 1) * nBlockDim
               A(indBlockI, indBlockJ) = SCRATCH(nVectorIndex)
               IF (bPIVOT .AND. indBlockI == IMAXPOP(NA)) THEN
                 ROWINVIONIMAX(j) = A(indBlockI, indBlockJ)
               ENDIF
            ENDDO
         ENDDO
         
C***     Off-Diagonal-Elements (only for ION SPLIT)
         IF (iBLOCKINVERSION == 2 .AND. (.NOT. bPIVOT)) THEN
C***       off-diagonal block is rectangular: 
C***            rows = rows of current block
C***         columns = columns of IIONMAX (PIVOT) 
C***       The offcenter blocks are calculated as  -1 * MATMUL( A^(-1), C ) 
C***            with C consisting of identical rows. The value for this row
C***            is stored in ROWINVIONIMAX when the IONMAX (=PIVOT) block is inverted
           DO i=1, nBlockDim
             SUMINVi = 0.
             indBlockI = nFirstLev + i - 1
             IF (ZERO_RATES(indBlockI)) CYCLE
             DO j=1, nBlockDim
               !First build the sum of A^(-1) of row i in the regular block
               indBlockJ = nFirstLev + j - 1
               SUMINVi = SUMINVi + A(indBlockI, indBlockJ)
             ENDDO
             DO j=1, nPivotDim
               !Now fill the offcenter block with 
               indBlockJ = JBSTART(NBFirstIon(NA)) + j - 1
               A(indBlockI, indBlockJ) = - SUMINVi * ROWINVIONIMAX(j)
             ENDDO
           ENDDO
         ENDIF

         IF (VERSION .EQ. 'SING') THEN
            WRITE (hCPR,'(A,I5,3A)') 
     >          'SING in block with first index', nFirstLev, 
     >          ' (',LEVEL(i),')'
            RETURN
         ENDIF
      ENDDO

C***  invert additional diagonal elements
      IF (NDIM .GT. NLAST(NATOM)) THEN
         DO nn=NLAST(NATOM)+1, NDIM
            A(nn, nn) = 1. / A(nn, nn)
         ENDDO
      ENDIF
C***  End of MATRIX INVERSION ---------------------------


C***  -------/ the following is unchanged from LINSOL /---------

C***  SOLUTION BY:   X = B * A-INV
      X = MATMUL(B, A)

C***  ITERATIVE IMPROVEMENT: ("Nachbrenner") ------------------------
      DO IMPRO=1, IMPMAX
        BTEST = MATMUL(X, ATEST)  !ACTUAL RIGHT-HAND SIDE VECTOR:  BTEST := X * A
        BTEST = BTEST - B         ! RIGHT-HAND SIDE ERROR:  BTEST := BTEST - B
        DTEST = MATMUL(BTEST, A)  !SOLUTION OF:   DTEST * A = BTEST    BY:  DTEST = BTEST * A-INV
        X = X - DTEST             !SUBTRACTION OF THE ERROR TERM: X(NEW) := X(OLD) - DTEST
      ENDDO
      
C***  TEST PRINTOUT: RIGHT-HAND SIDE VECTOR (OPTINAL!!!)  .............
      IF (bPrintRightVector) THEN
        !CALL VMF (BTEST, X, ATEST, N, NDIM)
        BTEST = MATMUL(X, ATEST)

        PRINT 11, VERSION, IMPMAX
   11   FORMAT (//, 10X, 'VERSION: ', A, 10X, 'ITERATIONS:',I2,
     $    //,10X,'INDEX', 11X, 'B(I)', 6X,'BTEST(I)',7X,'DEVIATION')


    3   FORMAT (10X, I5, 3(3X,G12.3))
        DO i=1,N
          DIFF = B(i) - BTEST(i)
          PRINT 3, i, B(i), BTEST(i), DIFF
        ENDDO
      ENDIF
C......................................................................

      RETURN
      END
