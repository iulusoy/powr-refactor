!> This module contains arrays for passing between steal and linpop.
!! Used by steal and linpop for now.
!! @param NDIM
!! @param NFDIM 
!! @param MAXIND
!! @param MAXFEIND
module linpop_arrays
use params
       REAL, DIMENSION(NDDIM) :: TNEW, RNE, ENTOT, TOLD, RADIUS, &
                                 OPA

       REAL, DIMENSION(NDIM) :: ENLTE, WEIGHT, NCHARG, EION, ELEVEL

       REAL, DIMENSION(NDIMP2) :: EN, V1, V2

       REAL, DIMENSION(NFDIM) :: FWEIGHT

       REAL, DIMENSION(MAXINDE) :: XJLAPP

       INTEGER, DIMENSION(NDDIM) :: ITNE

       CHARACTER(10), DIMENSION(NDIM) :: LEVEL

       LOGICAL, DIMENSION(MAXINDE) :: bUSEALO

       REAL, DIMENSION(NDDIM, NDIM) :: POPNUM, DEPART, POPLTE, POP1

       REAL, DIMENSION(NDIM, NDIM) :: EINST

       REAL, DIMENSION(NDDIM,NFDIM) :: XJC, WCHARM

       REAL, DIMENSION(NFDIM,NDDIM) :: SCOLD

       REAL, DIMENSION(NDDIM,MAXINDE) :: XJL, XLAMAPPMEAN, &
                                         XLAMAPPUMEAN, XLAMAPPLMEAN

       REAL, DIMENSION(NDIMP2,NDIMP2) :: DM

end module linpop_arrays