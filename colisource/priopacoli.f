      SUBROUTINE PRIOPACOLI (XLAM,K,ND,LSTEP,R,
     >                       OPA,ETA,THOMSON,FILLFAC,IWARN,
     >                       MAINPRO,MAINLEV,JOBNUM,MODHEAD)
C***********************************************************************
C***  PRINTOUT OF THE CONTINUUM OPACITIES ETC.
C***********************************************************************

      INTEGER, INTENT(IN) :: K, ND, LSTEP, JOBNUM
      REAL, INTENT(IN) :: XLAM
      
      REAL, DIMENSION(ND), INTENT(IN) :: OPA, ETA, THOMSON, R, FILLFAC
      CHARACTER(8), DIMENSION(ND) :: IWARN
      CHARACTER(10), DIMENSION(ND) :: MAINPRO, MAINLEV
      CHARACTER(100) :: MODHEAD
      
      REAL, EXTERNAL :: TRADFUN
      
      REAL :: TAU, TAUOLD, Q, S, RTAU1, TRAD, OPAL, OPALP, ETAL
      INTEGER :: L, IHELP, IHELPINIT
 
      IHELPINIT = 5HIHELP
      IF (IHELP.EQ.IHELPINIT) GOTO 1
      IHELP=5HIHELP
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1X,A,20X,'JOB NO.',I3,//,10X,
     $ 'CONTINUUM AND LINE OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,50('-'),//,
     $ ' FREQUENCY  DEPTH      OPACITY  THOMSON   OPTICAL    R(TAU=1)  '
     $ 'MAIN CONTRIBUTION      LASER    EMISSIVITY   SOURCE FUNCTION',/,
     $ '   INDEX    INDEX   (PER RSTAR) FRACTION   DEPTH               '
     $ 'PROCESS     LEVEL      WARNING    (...)        TRAD/KELVIN',/)
 
    1 PRINT 3
    3 FORMAT (' ')
      TAU=.0
      RTAU1=.0
 
      DO L=1, ND
        OPAL = OPA(L) * FILLFAC(L)
        ETAL = ETA(L) * FILLFAC(L)
        IF (L /= ND) THEN
          OPALP = OPA(L+1) * FILLFAC(L+1)
          TAUOLD=TAU
          TAU=TAU+0.5*(OPAL+OPALP)*(R(L)-R(L+1))
          IF (TAUOLD < 1. .AND. TAU >= 1.) THEN
            Q=(1.-TAUOLD)/(TAU-TAUOLD)
            RTAU1=(1.-Q)*R(L)+Q*R(L+1)
          ENDIF
        ENDIF               
       
        IF (((L-1)/LSTEP)*LSTEP /= (L-1) .AND. L /= ND) CYCLE
    
        IF (OPAL .LE..0) THEN
          TRAD=.0
        ELSE
          S=ETAL/OPAL/(1.-THOMSON(L))
          TRAD=TRADFUN (XLAM,S)
        ENDIF
        
        IF (L == ND) THEN
          PRINT 9,K,L,OPAL,THOMSON(L),TAU,RTAU1,
     $                  MAINPRO(L),MAINLEV(L),IWARN(L)(1:1),ETAL,TRAD
        ELSE
          PRINT 5,K,L,OPAL,THOMSON(L),TAUOLD,
     $                  MAINPRO(L),MAINLEV(L),IWARN(L)(1:1),ETAL,TRAD
        ENDIF
      ENDDO
    
      RETURN
 
    5 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2, 10X ,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
    9 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2,0P,F10.3,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
 
      END
