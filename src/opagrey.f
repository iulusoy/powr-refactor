      SUBROUTINE OPAGREY (OPARL,EN,TL,RNEL,ENTOTL,RSTAR,N,
     $                   NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,
     $                   NOM,EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,
     $                   MAXION,SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $                   KONTNUP,KONTLOW,LASTKON,RL,XDATA)
C***********************************************************************
C***  --- ONLY CALLED FROM SUBR. GREY ---
C***  FAST VERSION (VECTORIZED "COOPFRQ"!!!) OF ORIGINAL SUBR. OPAROSS:
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT DEPTH POINT L
C***  FOR GIVEN POPNUMBERS
C***********************************************************************
 
      DIMENSION EN(N)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF),EXPFAC(NF),OPAC(NF),ETAC(NF)
      DIMENSION SIGMAFF(NF,0:MAXION)
 
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /

C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN P.100 )
      DATA CFF / 1.370E-23 /

C***  SIGMAE = ELECTRON SCATTERING CROSS SECTION ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /

C***  PRE-CALCULATE EXPONENTIAL FACTORS AND FREE-FREE CROSS SECTIONS 
C***  FOR THE TEMPERATURE OF THE CURRENT DEPTH POINT
      TLOG=ALOG10(TL)
      ROOTTL=SQRT(TL)
      DO 1 K=1,NF
      XLAM=XLAMBDA(K)
      W=1.E8/XLAM
      EXPFAC(K)=EXP(-C1*W/TL)
      W3=W*W*W
      PRESIG=CFF/W3/ROOTTL
      XLAMLOG=ALOG10(XLAM)
      SIGMAFF(K,0)=0.0
      DO 1 ION=1,MAXION
      CALL GFFLOG (GIII,ION,XLAMLOG,TLOG)
      SIGMAFF(K,ION)=PRESIG*FLOAT(ION*ION)*GIII
    1 CONTINUE

      CALL COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $             WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $             RNEL,TL,SIGMAFF,MAXION,RL,XDATA,
     $             SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $             KONTNUP,KONTLOW,LASTKON,OPATHOM)

C***  FREQUENCY INTEGRATION
C***  NOTE: FOR NUMERICAL REASONS, THE NORMALIZATION CONSTANT Q IS
C**         ALSO INTEGRATED NUMERICALLY!
      SUM=.0
      Q=.0
C***  THOMSON SCATTERING OPACITY
      OPATH=RNEL*SIGMAE*ENTOTL*RSTAR
      DO 2 K=1,NF
C***  DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
      DBDT=DBNUEDT(XLAMBDA(K),TL)
C***  CONSIDERING THE THOMSON OPACITY
      OPAC(K)=OPAC(K)+OPATH
      SUM=SUM+DBDT*FWEIGHT(K)/OPAC(K)
      Q=Q+DBDT*FWEIGHT(K)
    2 CONTINUE
 
      OPARL=Q/SUM      

      RETURN
      END