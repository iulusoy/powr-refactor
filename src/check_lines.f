      SUBROUTINE CHECK_LINES(XK, XKMIN, XKMID, XKMAX, 
     >             LINECHECK, ILINECHECK, NLINE, MAXLIN, LEVEL,
     >             LIND, LINDS, WS, BPLOT, RADIUS, NLACT, LINE, 
C***  for PRELINECL
     >             NUP, LOW, N, XLAM, NDIM, ND, XJLMEAN, XLAMAPPMEAN,
     >             bALOTri, XLAMAPPUMEAN, XLAMAPPLMEAN,
     >             ELEVEL, INDNUP, INDLOW, NDDIM, 
C***  and also for LIOP
     >             EINST, WEIGHT, XLAMSOR, ENTOT, POPNUM, RSTAR, 
     >             OPAL, ETAL, VDOPUNIT, hMODEL, hALO)
C****************************************************************
C***  Adds new Lines to the List of Active Lines
C***    The Line Opacities are prepared
C***  Deletes Lines which are finished from the List
C***    The J-bars are writted to the Model File
C***
C***  Active lines are 1 ... NLACT and identified via the 
C***  LIND array. LIND returns the IND (line index) for an 
C***  active transition. LIND == 0 means that this entry is
C***  unused. Note that due to additions and removals LIND is
C***  not necessarily filled as a block, but instead might have
C***  a mix or zero and non-zero entries between 1 and MAXLIN.
C***
C***    Called by COLI
C****************************************************************

      INTEGER, INTENT(IN) :: hMODEL, hALO
      LOGICAL, INTENT(IN) :: bALOTri

      DIMENSION XKMIN(NLINE), XKMID(NLINE), XKMAX(NLINE), LINE(NLINE)
      DIMENSION LIND(MAXLIN), LINDS(MAXLIN), WS(NDDIM, MAXLIN)
      DIMENSION NUP(MAXLIN), LOW(MAXLIN), XLAM(MAXLIN)
      DIMENSION XLAMSOR(NLINE)
      DIMENSION RADIUS(ND)
      DIMENSION EINST(NDIM,NDIM), WEIGHT(NDIM)
      DIMENSION XJLMEAN(NDDIM,MAXLIN), XLAMAPPMEAN(NDDIM, MAXLIN)
      DIMENSION XLAMAPPUMEAN(NDDIM,MAXLIN), XLAMAPPLMEAN(NDDIM, MAXLIN)
      DIMENSION OPAL(NDDIM,MAXLIN), ETAL(NDDIM,MAXLIN)

      LOGICAL BPLOT
      REAL :: POPMIN

      CHARACTER(8) :: NAME
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      
      REAL, PARAMETER :: ALOMAX = 0.9999999999

C***  At the first call of this routine we have
C***  LINECHECK = 1  and  ILINECHECK = LINE(LINECHECK)
C***  This has been prepared in the subroutine PREPK, 
C***  together with the XKMIN and XKMAX arrays.

C***  In subsequent call, LINECHECK always keeps its current
C***  value. This is possible since the lines are ordered by
C***  wavelengths and thus each line has only one wavelength
C***  regime where it is active. (I.e., once its active region
C***  has passed, it will never be active again!)
      

C****************************************************************
C***  Check the currently active lines if they can be checked out
C****************************************************************
      DO NL=1, MAXLIN
         LACT = LIND(NL)
         LACTS = LINDS(NL)
         IF (LACT .EQ. 0) CYCLE
C***     Line has been finished? -> Store XJL and check-out
C***      (for iron lines, storing is done in WMODCOLI)
         IF (XKMAX(LACTS) .LT. XK) THEN
            DO L=1, ND
              XJLMEAN(L,NL) = XJLMEAN(L,NL) / WS(L,NL)
     >                         /RADIUS(L)/RADIUS(L)
              XLAMAPPMEAN(L,NL) = XLAMAPPMEAN(L,NL) / WS(L,NL)
     >                         /RADIUS(L)/RADIUS(L)
C***          Limit to values between 0 and ALOMAX:
              XLAMAPPMEAN(L,NL) = MAX(0., XLAMAPPMEAN(L,NL))
              XLAMAPPMEAN(L,NL) = MIN(XLAMAPPMEAN(L,NL), ALOMAX)
              IF (bALOTri) THEN
                XLAMAPPUMEAN(L,NL) = XLAMAPPUMEAN(L,NL) / WS(L,NL)
     >                         /RADIUS(L)/RADIUS(L)
                XLAMAPPLMEAN(L,NL) = XLAMAPPLMEAN(L,NL) / WS(L,NL)
     >                         /RADIUS(L)/RADIUS(L)
              ENDIF
            ENDDO

            IF (LACT <= 9999) THEN
              WRITE (NAME,'(A3,I4,A1)') 'XJL',LACT,' '
            ELSE
              WRITE (NAME,'(A3,I5)') 'XJL',LACT
            ENDIF
C***        Write J_L_bar into MODEL file
            CALL WRITMS(hMODEL, 
     >        XJLMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
C***        empty entries are marked by LIND=0
C***        Write freq.-integrated ALO into ALO file
            WRITE (NAME,'(A3,I5)') 'ALO', LACT
            CALL WRITMS(hALO, 
     >        XLAMAPPMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
C***        Tri-diagonal ALO is only calculated and stored
C***          if requested via CARDS line (XJLAPP COLI TRI)
            IF (bALOTri) THEN
              WRITE (NAME,'(A3,I5)') 'ALU', LACT
              CALL WRITMS(hALO, 
     >          XLAMAPPUMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
              WRITE (NAME,'(A3,I5)') 'ALL', LACT
              CALL WRITMS(hALO, 
     >          XLAMAPPLMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
            ENDIF
            LIND(NL) = 0
            NLACT = NLACT - 1
cc         Plotting facility not documented
cc         ELSE
cc            IF (BPLOT) THEN
cc              WRITE (40+NL,'(I8,1X,F3.0)') NINT(XK), FLOAT(NL)
cc            ENDIF
         ENDIF
      ENDDO


C***  **************************************
C***  Which lines are newly becoming active?
C***  **************************************
      DO
         IF (XK .LT. XKMIN(LINECHECK) .OR. LINECHECK .GT. NLINE) EXIT

         IF (XK .LE. XKMAX(LINECHECK)) THEN
C***     Search first free index
            DO NLNEW=1, MAXLIN
              IF (LIND(NLNEW) .EQ. 0) THEN
C***          Free table entry found at NLNEW
                LIND(NLNEW) = ILINECHECK
                LINDS(NLNEW) = LINECHECK
                DO L=1, ND
                  WS(L, NLNEW) = 0.
                  XLAMAPPMEAN(L, NLNEW) = 0.
                  IF (bALOTri) THEN
                    XLAMAPPUMEAN(L, NLNEW) = 0.
                    XLAMAPPLMEAN(L, NLNEW) = 0.
                  ENDIF
                ENDDO
C***            Prepare line quantities
                CALL PRELINE (NUP(NLNEW), LOW(NLNEW), ILINECHECK, N,
     >                 XLAM(NLNEW), ND, XJLMEAN(1,NLNEW),
     >                 ELEVEL, NDIM, INDNUP, INDLOW)
                CALL LIOP (EINST(NUP(NLNEW),LOW(NLNEW)),
     >                 WEIGHT(LOW(NLNEW)),
     >                 WEIGHT(NUP(NLNEW)), LOW(NLNEW), NUP(NLNEW), ND,
     >                 XLAMSOR(LINECHECK), ENTOT, POPNUM, RSTAR,
     >                 OPAL(1,NLNEW), ETAL(1,NLNEW), VDOPUNIT)
c                IF (BPLOT) THEN
c                  WRITE (38,'(A,F9.2,1X,F3.0,1X,A,I4,1X,I4,F11.3)')
c     >              'KASDEF LUN ', XK, FLOAT(NLNEW), '0. 0.1 0.2 &E',
c     >              LINECHECK, LINE(LINECHECK),
c     >              XLAMSOR(LINECHECK)
c                  WRITE (38,'(A,F9.2,A)')
c     >              'KASDEF SYM ', XKMID(LINECHECK),
c     >              ' 0.0 0. 0.2 -0.2 8'
c                  WRITE (38,'(A,F9.2,1X,F5.1,A)')
c     >              'KASDEF SYM ', XKMID(LINECHECK),
c     >              FLOAT(NLNEW), ' 0. -0.2 -0.2 8'
c                ENDIF
                NLACT = NLACT + 1
                GOTO 10
              ENDIF
            ENDDO

C***        If no free entry was found: ERROR STOP:
            WRITE (0,*) '*** ERROR: Capacity of MAXLIN exceeded'
            STOP '*** FATAL ERROR in Subr. COLI -> CHECK_LINES'

   10       CONTINUE

            LINECHECK = LINECHECK + 1
            ILINECHECK = LINE(LINECHECK)
         ENDIF
      ENDDO

      RETURN
      END
