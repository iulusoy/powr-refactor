      SUBROUTINE WMODCOLI(XJCINT, FWTEST, ARAD, ACONT, ATHOM,
     >                    ND, NF, RADIUS, ENTOT, RSTAR,
     >                    XJFEMEAN, SIGMAINT, LASTFE, HTOTL, WFELOW, 
     >                    WFENUP, FTCOLI, FTCOLIB, REDISMODE, NCOLIP, 
     >                    XJTOTL, XKTOTL, XNTOTL, WJC,
     >                    DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
     >                    OPASMEAN, OPASMEANTC, OPAPMEAN, 
     >                    QFJMEAN, SMEAN, OPAJMEAN, OPAJMEANTC, 
     >                    QOPAHMEAN, EDDIHOUTJMEAN, HMEAN,
     >                    HTOTOUTMINUS, HTOTMINUSND, HTOTND, 
     >                    HTOTNDCOR, HTOTCUT,
     >                    LASTIND,  FERATUL, FERATLU, BCOLIP, EPSGMAX, 
     >                    FTFE, EMCOLI, FF_INFO, IFF_DK, IFF_MAX_MS,
     >                    IFF_WCHARM, OPALAMBDAMEAN, TOTOUT, bKALPHA,
     >                    ALPHAF, ARMDRESP, GAMMATHIN, CK, RHO, XMU, 
     >                    TAUROSS, OPAROSS, OPAROSSELEM, OPAROSSCONT, 
     >                    OPATOTELEM, NATOM, 
     >                    ARADELEM, ACONTELEM, MAXION,
     >                    ARADION, ACONTION, hMODEL, hALO)
C****************************************************************
C*** Handles all MODEL write (MSWRIT) statements from COLI
C***  only called by COLI main program
C****************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER, INTENT(IN) :: ND, NF, NATOM, MAXION, LASTFE, IFF_MAX_MS,
     >                       hMODEL, hALO

      REAL, DIMENSION(ND) :: ARAD, ACONT, ATHOM,
     >                       HTOTL, HTOTCUT, FTCOLI, FTCOLIB,
     >                       XJTOTL, XKTOTL, XNTOTL,
     >                       EPSGMAX, OPALAMBDAMEAN, SMEAN,
     >                       ALPHAF, ARMDRESP, GAMMATHIN, CK,
     >                       RHO, XMU, OPAROSSCONT,          
     >                       QOPAHMEAN, QOPAHMEANL, HMEAN,
     >                       OPAJMEANTC, OPAJMEAN, OPASMEAN,
     >                       OPASMEANTC, OPAPMEAN,
     >                       QFJMEAN

      REAL, DIMENSION(LASTFE) :: SIGMAINT
      REAL, DIMENSION(ND,NF) :: XJCINT, WJC
      REAL, DIMENSION(LASTFE,ND) :: FERATLU, FERATUL
      REAL, DIMENSION(ND,LASTFE) :: XJFEMEAN, FTFE, WFELOW, WFENUP
      REAL, DIMENSION(NF) :: EMCOLI, FWTEST
      REAL, DIMENSION(10) :: FF_INFO
C      INTEGER, DIMENSION(IFF_MAX_MS) ::  IFF_DK
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX_MS) :: IFF_DK   !erst, wenn MS-Storage aktualisiert
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX_MS,ND) :: IFF_WCHARM
      
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS, OPAROSS,
     >                                   RADIUS, ENTOT
      REAL, DIMENSION(NATOM, ND), INTENT(IN) :: OPAROSSELEM, OPATOTELEM
      REAL, DIMENSION(NATOM, ND-1), INTENT(IN) :: ARADELEM, ACONTELEM
      REAL, DIMENSION(ND-1, NATOM, MAXION), INTENT(IN) :: ARADION, 
     >                                                    ACONTION

      INTEGER :: L, INDF, NCOLIP, LASTIND, IND, INDFE, IDUMMY, IERR,
     >           IFF_N_MS, LARM
      REAL :: EDDIHOUTJMEAN, DBDTINT_M, DBDTOPAINT_M, OPARNDM,
     >        DBDTOPAINT, DBDTINT, OPARND, RSTAR, HTOTOUTMINUS,
     >        TOTOUT, HTOTMINUSND, HTOTND, HTOTNDCOR

      LOGICAL :: BCOLIP, bKALPHA
      
      CHARACTER(8) :: NAME, CKIND
      CHARACTER(4) :: REDISMODE


C***  NORMALIZE AND SAVE XJC'S
      DO INDF=1, NF
         DO L=1, ND
            IF (FWTEST(INDF) .NE. 0.) THEN
               XJCINT(L, INDF) = XJCINT(L, INDF) / FWTEST(INDF)
     >                                 /RADIUS(L)/RADIUS(L)
C!!!  The Division by FWTEST is now done in FREQUNORM
C!!!               WJC(L,INDF) = WJC(L,INDF)/FWTEST(INDF)
            ELSE
               XJCINT(L, INDF) = 0.
               WJC(L,INDF) = 0.
            ENDIF
         ENDDO
         WRITE (NAME,'(A3,I4)') 'XJC', INDF
         CALL WRITMS 
     >        (hMODEL, XJCINT(1, INDF), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WJC', INDF
         CALL WRITMS 
     >        (hALO,  WJC(1, INDF), ND, NAME, -1, IDUMMY, IERR)
      ENDDO

      DO INDF=1, NF
         EMCOLI(INDF)=EMCOLI(INDF)/FWTEST(INDF)
      ENDDO
      CALL WRITMS (3,EMCOLI,NF,'EMCOLI  ',-1, IDUMMY, IERR)


      CALL WRITMS(3, TAUROSS,ND, 'TAUROSS ', -1, IDUMMY, IERR)
      CALL WRITMS(3, OPAROSS,ND, 'OPAROSS ', -1, IDUMMY, IERR)
      LARM = NATOM * ND
      CALL WRITMS(3,OPAROSSELEM,LARM,'OPARELEM', -1, IDUMMY, IERR)
      CALL WRITMS(3,OPAROSSCONT,ND,  'OPARCONT', -1, IDUMMY, IERR)
      CALL WRITMS(3,OPATOTELEM,LARM,'OPALELEM', -1, IDUMMY, IERR)

      CALL WRITMS(3, ARAD, ND-1, 'ARAD    ', -1, IDUMMY, IERR)
      CALL WRITMS(3, ACONT, ND-1, 'ACONT   ', -1, IDUMMY, IERR)
      CALL WRITMS(3, ATHOM, ND-1, 'ATHOM   ', -1, IDUMMY, IERR)
      LARM = NATOM * (ND-1) 
      CALL WRITMS(3, ARADELEM, LARM, 'ARADELEM', -1, IDUMMY, IERR)
      CALL WRITMS(3, ACONTELEM,LARM, 'ACNTELEM', -1, IDUMMY, IERR)
      LARM = (ND-1) *  NATOM * MAXION
      CALL WRITMS(3, ARADION, LARM, 'ARADION ', -1, IDUMMY, IERR)
      CALL WRITMS(3, ACONTION,LARM, 'ACONTION', -1, IDUMMY, IERR)
      
C      CALL WRITMS(3, RHO,    ND, 'RHO     ', -1, IDUMMY, IERR)
C      CALL WRITMS(3, XMU,    ND, 'XMU     ', -1, IDUMMY, IERR)
C***  save radiative acceleration and fractions in internal units
C      (this will be used in STEAL->HYDROSOLVE)
      IF (bKALPHA) THEN
        CALL WRITMS(3, ALPHAF, ND, 'ALPHAF  ', -1, IDUMMY, IERR)
      ENDIF



C***  CALCULATE AND SAVE LOCAL RADIATIVE ENERGY LOSS IN CGS
      DO L=1, ND
         FTCOLI(L) = FTCOLI(L)/RSTAR
         FTCOLIB(L) = FTCOLIB(L)/RSTAR
      ENDDO
      CALL WRITMS(3, FTCOLI, ND, 'FTCOLI  ', -1, IDUMMY, IERR)
      CALL WRITMS(3, FTCOLIB,ND, 'FTCOLIB ', -1, IDUMMY, IERR)

C***  IRON: NORMALIZE AND SAVE MEAN INTENSITIES 'XJFEMEAN'
C           AND DIAGONAL WEIGHTS 'WFELOW' AND 'WFENUP'
      DO INDFE = 1, LASTFE
         DO L=1, ND
            XJFEMEAN(L, INDFE) = XJFEMEAN(L, INDFE) / SIGMAINT(INDFE)
     >                           /RADIUS(L)/RADIUS(L)
            WFELOW(L, INDFE)   = WFELOW(L, INDFE) / SIGMAINT(INDFE)
            WFENUP(L, INDFE)   = WFENUP(L, INDFE) / SIGMAINT(INDFE)
            FTFE(L,INDFE)      = FTFE(L,INDFE)/RSTAR
         ENDDO
         IND = LASTIND + INDFE
C***     XJL for iron lines is stored here, for all other lines
C***     this has already been done in the subroutine CHECK_LINES
         IF (IND <= 9999) THEN
            WRITE (NAME,'(A3,I4,A1)') 'XJL', IND, ' '
         ELSE
            WRITE (NAME,'(A3,I5)') 'XJL', IND
         ENDIF
         CALL WRITMS 
     >        (3, XJFEMEAN(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WFL', INDFE
         CALL WRITMS 
     >        (hALO, WFELOW(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WFU', INDFE
         CALL WRITMS 
     >        (hALO, WFENUP(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
C***     FTFE is for tests only and not used in other programs
C***     therefore it is not written to the MODEL anymore - wrh 17-May-2023
ccc         WRITE (NAME,'(A3,I4)') 'FTF', INDFE
ccc         CALL WRITMS
ccc     >        (3, FTFE(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
      ENDDO

C***  Radiative Rates for Iron superlevels
      IF (LASTFE .GT. 0) THEN
         DO L=1, ND
            WRITE (NAME, '(A5, I3)') 'FERLU', L
            CALL WRITMS 
     >           (3, FERATLU(1,L), LASTFE, NAME, -1, IDUMMY, IERR)
            WRITE (NAME, '(A5, I3)') 'FERUL', L
            CALL WRITMS 
     >           (3, FERATUL(1,L), LASTFE, NAME, -1, IDUMMY, IERR)
         ENDDO
      ENDIF


C***  Save the integrated Flux of all Lines
C***  Note that HTOTL and NTOTL are on interstices (i.e. only valid ND-1 entries)
C!!!  Unsinnige Laenge beibehalten wegen Inkompatibilitaet mit aelteren Modellen
      CALL WRITMS(3, HTOTL, ND, 'HTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XJTOTL, ND, 'JTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XKTOTL, ND, 'KTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XNTOTL, ND, 'NTOTL   ', -1, IDUMMY, IERR)
C***  non-negative HTOTL (can be used in TEMPCORR)
      CALL WRITMS(3, HMEAN, ND-1, 'HMEAN   ', -1, IDUMMY, IERR)
C***  HTOTL with too strong lines cut (can be used in TEMPCORR)      
      CALL WRITMS(3, HTOTCUT, ND-1,  'HTOTCUT ', -1, IDUMMY, IERR)

C***  Calculate and Save Rosseland Opacity at the inner Boundary
      OPARND = DBDTINT / DBDTOPAINT
      CALL WRITMS(3, OPARND, 1, 'OPARND  ', -1, IDUMMY, IERR)
      OPARNDM = DBDTINT_M / DBDTOPAINT_M
      CALL WRITMS(3, OPARNDM, 1, 'OPARNDM ', -1, IDUMMY, IERR)
      
C***  Save special Hminus at inner boundary (used in STEAL for T(ND)-correction)      
      CALL WRITMS(3, HTOTMINUSND, 1, 'HTMND   ', -1, IDUMMY, IERR)
      CALL WRITMS(3, HTOTND     , 1, 'HTOTND  ', -1, IDUMMY, IERR)
      CALL WRITMS(3, HTOTNDCOR  , 1, 'HTNDCOR ', -1, IDUMMY, IERR)
C***  Save special Hminus at outer boundary (used in STEAL for T-correction)      
      CALL WRITMS(3, HTOTOUTMINUS, 1,'HTOTOUTM', -1, IDUMMY, IERR)

C***  Write REDISMODE and NCOLIP
      CALL WRITMS(3, REDISMODE, 1, 'REDISMO ', -1, IDUMMY, IERR)
      CALL WRITMS(3, NCOLIP,    1, 'NCOLIP  ', -1, IDUMMY, IERR)

C***  Unsoeld-Lucy Terms for TEMPEQ
      CALL WRITMS (3,OPASMEAN,     ND, 'OPASMEAN', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAJMEAN,     ND, 'OPAJMEAN', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAPMEAN,     ND, 'OPAPMEAN', -1, IDUMMY, IERR)
      CALL WRITMS (3,QFJMEAN,      ND, 'QFJMEAN ', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPASMEANTC,   ND, 'OPASMTC ', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAJMEANTC,   ND, 'OPAJMTC ', -1, IDUMMY, IERR)
      CALL WRITMS (3,SMEAN,        ND, 'SMEAN   ', -1, IDUMMY, IERR)
      CALL WRITMS (3,QOPAHMEAN,  ND-1, 'QOPAHMEA', -1, IDUMMY, IERR)
      CALL WRITMS (3,EDDIHOUTJMEAN, 1, 'EDDIHOJM', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPALAMBDAMEAN,ND, 'OPALMEAN', -1, IDUMMY, IERR)
      
      IF (BCOLIP) THEN
        CALL WRITMS(3, EPSGMAX, ND-1, 'EPSGMAX ', -1, IDUMMY, IERR)
      ENDIF

      CALL WRITMS (3,TOTOUT,1, 'TOTOUT  ',-1, IDUMMY, IERR)

      RETURN
      END
