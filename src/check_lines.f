      SUBROUTINE CHECK_LINES(XK, XKMIN, XKMID, XKMAX, 
     >             LINECHECK, ILINECHECK, NLINE, MAXLIN, LEVEL,
     >             LIND, LINDS, WS, BPLOT, RADIUS, NLACT, LINE, 
C***  for PRELINECL
     >             NUP, LOW, N, XLAM, NDIM, ND, XJLMEAN, ELEVEL, 
     >             INDNUP, INDLOW, NDDIM, 
C***  and also for LIOP
     >             EINST, WEIGHT, XLAMSOR, ENTOT, POPNUM, RSTAR, 
     >             OPAL, ETAL, VDOP)
C****************************************************************
C***  Check-out from the list those lines which are finished 
C***    --> J-bars are writted to the Model File
C***  Add new Lines to the List of Active Lines
C***    The Line Opacities are prepared
C***
C***    Called by COLI
C****************************************************************

      DIMENSION XKMIN(NLINE), XKMID(NLINE), XKMAX(NLINE), LINE(NLINE)
      DIMENSION LIND(MAXLIN), LINDS(MAXLIN), WS(MAXLIN)
      DIMENSION NUP(MAXLIN), LOW(MAXLIN), XLAM(MAXLIN)
      DIMENSION XLAMSOR(NLINE)
      DIMENSION RADIUS(ND)
      DIMENSION EINST(NDIM,NDIM), WEIGHT(NDIM)
      DIMENSION XJLMEAN(NDDIM,MAXLIN)
      DIMENSION OPAL(NDDIM,MAXLIN), ETAL(NDDIM,MAXLIN)

      LOGICAL BPLOT

      CHARACTER NAME*8
      CHARACTER LEVEL(NDDIM)*10

C****************************************************************
C***  Check the currently active lines if they can be checked out
C****************************************************************
      DO NL=1, MAXLIN
         LACT = LIND(NL)
         LACTS = LINDS(NL)
         IF (LACT .EQ. 0) CYCLE
C***     Line has been finished? -> Store XJL and check-out
         IF (XKMAX(LACTS) .LT. XK) THEN
            DO L=1, ND
              XJLMEAN(L,NL) = XJLMEAN(L,NL) / WS(NL)
     >                         /RADIUS(L)/RADIUS(L)
            ENDDO
            IF (LACT <= 9999) THEN
              WRITE (NAME,'(A3,I4,A1)') 'XJL',LACT,' '
            ELSE
              WRITE (NAME,'(A3,I5)') 'XJL',LACT
            ENDIF
            CALL WRITMS 
     >        (3, XJLMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
C***        empty entries are marked by LIND=0
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
                WS(NLNEW) = 0.
C***            Prepare line quantities
                CALL PRELINE (NUP(NLNEW), LOW(NLNEW), ILINECHECK, N, 
     >                 XLAM(NLNEW), ND, XJLMEAN(1,NLNEW), 
     >                 ELEVEL, NDIM, INDNUP, INDLOW)
                CALL LIOP (EINST(NUP(NLNEW),LOW(NLNEW)), 
     >                 WEIGHT(LOW(NLNEW)), 
     >                 WEIGHT(NUP(NLNEW)), LOW(NLNEW), NUP(NLNEW), ND, 
     >                 XLAMSOR(LINECHECK), ENTOT, POPNUM, RSTAR, 
     >                 OPAL(1,NLNEW), ETAL(1,NLNEW), VDOP)
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
