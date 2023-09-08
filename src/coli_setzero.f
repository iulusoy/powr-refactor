      SUBROUTINE COLI_SETZERO(ND, NDDIM, NIT, NPDIM, NFDIM, MAXATOM,
     >             MAXFEIND, MAXLIN, MAXIND, MAXION,
     >             DBDTINT, DBDTOPAINT, EDDIHOUTJMEAN, 
     >             HTOTOUTMINUS, DBDTINT_M, DBDTOPAINT_M,
C***  with ND (NDDIM)
     >             OPA, ETA,
     >             XJTOTL, HTOTL, XKTOTL, XNTOTL, ARAD, ACONT, ATHOM, 
     >             FTCOLI, FTCOLIB, OPAKOLD, ETAKOLD, OPAKNOTHO, 
     >             ETAKNOTHO, OPAO, THOMSONO, ETANOTHO,
     >             DJDSMOD_OLD, 
     >             OPASMEAN, QFJMEAN, OPAJMEAN, OPASMEANTC, OPAJMEANTC, 
     >             OPAPMEAN, SMEAN, QLFOLD, EPSGMAX, OPAROSS, 
     >             OPALAMBDAMEAN,
     >             OPAROSSELEM, OPAROSSCONT,
     >             OPAKFEOLD, iIgnoreK, OPAKFEFTOLD, ETAKFEFTOLD, 
     >             OPAKNOFENOTHO, ETAKNOFENOTHO,
     >             bOSKIPLAST, XLAMLASTOSKIP, XJLMO, XJLMOR2,
C***  with ND-1 (NDDIM-1)
     >             QOPAHMEAN, HMEAN, QLHOLD, OPAKHOLD,
     >             HTOTLTEST, HTOTCUT,
C***  with NDDIM,NIT
     >             XJLOLD, XJLMO_OLD, EDDIFO, S_OLD, OPAK_OLD, EPSG,
C***  with NDDIM,NPDIM
     >             CWM0, CWM2, CWM1, CWM3,
C***  with NDDIM,NFDIM,
     >             WJC,
C***  with NDDIM,NPDIM,NIT
     >             XIPLUS_OLD, XIMINUS_OLD,
C***  with NDDIM-1, NIT; In COLI it is also NDDIM,NIT
     >             XHLOLD, XHLMO_OLD, EDDIGO, 
C***  with MAXFEIND, NDDIM
     >             XJFEMEAN, FERATLU, FERATUL, FTFE, WFELOW, WFENUP,
C***  with MAXATOM,ND
     >             OPAKELEM, OPAKOLDELEM, OPACELEM, OPACOLDELEM,
     >             ETAKELEM, ETAKOLDELEM, ETACELEM, ETACOLDELEM,
     >             OPATOTELEM,
C***  with MAXATOM,ND-1
     >             ARADELEM, ACONTELEM,
C***  with ND, MAXATOM, MAXION
     >             OPAKION, OPAKOLDION, OPACION, OPACOLDION,
     >             ETAKION, ETAKOLDION, ETACION, ETACOLDION,
C***  with ND-1, MAXATOM, MAXION
     >             ARADION, ACONTION,
C***  with MAXLIN
     >             LIND, LINDS,
C***  with NDDIM, MAXLIN
     >             WS,
C***  with MAXIND
     >             BLASERL,
C***  with NFDIM
     >             EMCOLI,
C***  no Arrays
     >             HTOTMINUSND, HTOTND, HTOTNDS, HTOTNDCOR,
     >             OPAMAX1, OPAMAX1_LAMBDA, IOPAMAX1_K)

C****************************************************************
C***  Presets all given variables to zero
C***    Called by COLI
C***
C***  Note: For the line opacities only the "old"-variables are
C***        nulled here, the variables for the current K index
C***        are nulled in ADDOPA
C****************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NDDIM, NIT, NPDIM, NFDIM, MAXATOM,
     >                       MAXFEIND, MAXLIN, MAXIND, MAXION
      REAL, DIMENSION(NDDIM) :: XJTOTL, HTOTL, XKTOTL, XNTOTL, 
     >                       ARAD, ACONT, ATHOM, FTCOLI, FTCOLIB,
     >                       OPAKOLD, ETAKOLD, OPAKNOTHO, ETAKNOTHO, 
     >                       OPAO, THOMSONO, ETANOTHO, 
     >                       DJDSMOD_OLD, OPASMEAN, 
     >                       QFJMEAN, OPAJMEAN, OPASMEANTC, 
     >                       OPAJMEANTC, OPAPMEAN, SMEAN,
     >                       QLFOLD, 
     >                       EPSGMAX, OPAROSS,
     >                       OPAKNOFENOTHO, ETAKNOFENOTHO,
     >                       OPAKFEFTOLD, ETAKFEFTOLD, HTOTLTEST, 
     >                       HTOTCUT, XLAMLASTOSKIP, XJLMO, XJLMOR2
      LOGICAL, DIMENSION(NDDIM) :: bOSKIPLAST
      INTEGER, DIMENSION(NDDIM) :: iIgnoreK
      REAL, DIMENSION(NDDIM) :: OPA, ETA, 
     >                          OPAKFEOLD, OPAROSSCONT, 
     >                          OPALAMBDAMEAN
      REAL, DIMENSION(NDDIM,NIT) :: S_OLD, OPAK_OLD, EPSG
      REAL, DIMENSION(NDDIM,NPDIM) :: CWM0, CWM2, CWM1, CWM3
      REAL, DIMENSION(NDDIM, MAXFEIND) :: WFELOW, WFENUP, FTFE, 
     >                                    XJFEMEAN
      REAL, DIMENSION(NDDIM,NFDIM) :: WJC
      
      REAL, DIMENSION(NDDIM-1) :: QOPAHMEAN, HMEAN, QLHOLD, OPAKHOLD
      
      REAL, DIMENSION(NDDIM,NIT) :: XJLOLD, XJLMO_OLD, EDDIFO, 
     >                              XHLOLD, XHLMO_OLD, EDDIGO

      REAL, DIMENSION(NDDIM,NPDIM,NIT) :: XIPLUS_OLD, XIMINUS_OLD

      REAL, DIMENSION(NFDIM) :: EMCOLI
      
      REAL, DIMENSION(MAXFEIND,NDDIM) :: FERATLU, FERATUL

      REAL, DIMENSION(MAXATOM, NDDIM) :: OPAKELEM, OPAKOLDELEM, 
     >                                   OPACELEM, OPACOLDELEM,
     >                                   ETAKELEM, ETAKOLDELEM, 
     >                                   ETACELEM, ETACOLDELEM, 
     >                                   OPAROSSELEM, OPATOTELEM
      REAL, DIMENSION(MAXATOM, NDDIM-1) :: ARADELEM, ACONTELEM

      REAL, DIMENSION(NDDIM, MAXATOM, MAXION) :: OPAKION, OPAKOLDION, 
     >                                           OPACION, OPACOLDION,
     >                                           ETAKION, ETAKOLDION, 
     >                                           ETACION, ETACOLDION
      REAL, DIMENSION(NDDIM, MAXATOM, MAXION) :: ARADION, ACONTION
      
      INTEGER, DIMENSION(MAXLIN) :: LIND, LINDS
      REAL, DIMENSION(NDDIM, MAXLIN) :: WS
      
      LOGICAL, DIMENSION(MAXIND) :: BLASERL
      
      INTEGER :: IOPAMAX1_K, IND, L
      REAL :: OPAMAX1_LAMBDA, OPAMAX1, HTOTOUTMINUS, EDDIHOUTJMEAN,
     >        DBDTOPAINT_M, DBDTINT_M, DBDTOPAINT, DBDTINT, 
     >        HTOTMINUSND, HTOTND, HTOTNDCOR, HTOTNDS
     
C***  Set Zero
      DBDTINT        = 0.
      DBDTOPAINT     = 0.
      DBDTINT_M      = 0.
      DBDTOPAINT_M   = 0.
      EDDIHOUTJMEAN  = 0.
      HTOTOUTMINUS   = 0.
      HTOTMINUSND    = 0.
      HTOTND         = 0.
      HTOTNDS        = 0.
      HTOTNDCOR      = 0.
      HTOTLTEST      = 0.
      HTOTCUT        = 0.
 
      OPA            = 0.
      ETA            = 0.
      CWM0           = 0.
      CWM1           = 0.
      CWM2           = 0.
      CWM3           = 0.

      XJTOTL         = 0.
      HTOTL          = 0.
      XKTOTL         = 0.
      XNTOTL         = 0.
      ARAD           = 0.
      ACONT          = 0.
      ATHOM          = 0.
      FTCOLI         = 0.
      FTCOLIB        = 0.
      OPAKOLD        = 0.
      ETAKOLD        = 0.
      OPAKNOTHO      = 0.
      ETAKNOTHO      = 0.
      OPAO           = 0.
      THOMSONO       = 0.
      ETANOTHO       = 0.
      DJDSMOD_OLD    = 0.
      OPASMEAN       = 0.
      QFJMEAN        = 0.
      OPAJMEAN       = 0.
      OPASMEANTC     = 0.
      OPAJMEANTC     = 0.
      OPAPMEAN       = 0.
      OPAROSS        = 0.
      OPALAMBDAMEAN  = 0.
      SMEAN          = 0.
      QLFOLD         = 0.
      EPSGMAX        = 0.
      iIgnoreK       = 0
      OPAKNOFENOTHO  = 0.
      OPAKNOFENOTHO  = 0.
 
      QOPAHMEAN      = 0.
      HMEAN          = 0.
      QLHOLD         = 0.
      OPAKHOLD       = 0.

      XJLMO          = 0.
      XJLMOR2        = 0. 
      
      XLAMLASTOSKIP  = 0.
      DO L=1, NDDIM
        bOSKIPLAST(L) = .FALSE.
      ENDDO
 
      XJLOLD         = 0.
      XJLMO_OLD      = 0.
      EDDIFO         = 0.
      S_OLD          = 0.
      OPAK_OLD       = 0.
      EPSG           = 0.
C***  BLUE BOUNDARY CONDITION FOR SHORTRAY 
      XIPLUS_OLD     = 0.
      XIMINUS_OLD    = 0.

      XHLOLD         = 0.
      XHLMO_OLD      = 0.
      EDDIGO         = 0.
 
      XJFEMEAN       = 0.
      FERATLU        = 0.
      FERATUL        = 0.
      FTFE           = 0.
      WFELOW         = 0.
      WFENUP         = 0.

      LIND           = 0
      LINDS          = 0
      WS             = 0.
      DO IND=1, MAXIND
        BLASERL(IND) = .FALSE.
      ENDDO

      OPAMAX1        = 0.
      OPAMAX1_LAMBDA = 0.
      IOPAMAX1_K     = 0

      EMCOLI         = 0.
      OPAROSSELEM      = 0.
      OPAROSSCONT    = 0.
      OPAKFEOLD      = 0.
      OPAKFEFTOLD    = 0.
      ETAKFEFTOLD    = 0.

      OPACELEM       = 0.
      OPACOLDELEM    = 0.
      OPAKELEM       = 0.
      OPAKOLDELEM    = 0.
      ETACELEM       = 0.
      ETACOLDELEM    = 0.
      ETAKELEM       = 0.
      ETAKOLDELEM    = 0.
      
      OPATOTELEM     = 0.
      
      ARADELEM       = 0.
      ACONTELEM      = 0.

      OPACION       = 0.
      OPACOLDION    = 0.
      OPAKION       = 0.
      OPAKOLDION    = 0.
      ETACION       = 0.
      ETACOLDION    = 0.
      ETAKION       = 0.
      ETAKOLDION    = 0.

      ARADION       = 0.
      ACONTION      = 0.
      
      RETURN
      END
