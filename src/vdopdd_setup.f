      SUBROUTINE VDOPDD_SETUP(bDDVDOP, VDOPDD, VDOPUNIT, DXMAX,
     >                         T, VMIC, ATMASS, ND, NATOM)
C***********************************************************************      
C**** fills the depth-dependent VDOP array VDOPDD 
C**** used in the line profile calculations for the radiative transfer
C**** In case of depth-dependent VDOP, the frequency-spacing parameter
C****  DXMAX is also adjusted
C****
C**** called from COLI
C***********************************************************************      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND, NATOM
      REAL, INTENT(IN) :: VDOPUNIT
      REAL, INTENT(INOUT) :: DXMAX
      LOGICAL, INTENT(IN) :: bDDVDOP
      
      REAL, DIMENSION(ND), INTENT(IN) :: T, VMIC
      REAL, DIMENSION(NATOM), INTENT(IN) :: ATMASS
      REAL, DIMENSION(ND, NATOM), INTENT(INOUT) :: VDOPDD

      REAL :: MASSNA, VMICCM2, VDOPMIN, VDOPMAX
      INTEGER :: L, NA

C***  Physical constants      
      REAL, PARAMETER :: BOLTZK = 1.3807E-16    !BOLTZMANN CONSTANT on cgs units
      REAL, PARAMETER :: AMU = 1.6605E-24       !Atomic mass unit (gramm) = m_H
      
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (bDDVDOP) THEN
        WRITE (hCPR,*) 'COLI uses depth-dependent Doppler profiles!'
      ENDIF
      VDOPMIN = VDOPUNIT
      VDOPMAX = VDOPUNIT
      DO L=1, ND
        DO NA=1, NATOM
          IF (bDDVDOP) THEN
C***        Modern option: depth-dependent VDOP for line profiles from
C***                       thermal and (micro-)turbulent broadening
            MASSNA = AMU * ATMASS(NA)
            VMICCM2 = (VMIC(L)*1.E5)**2
            VDOPDD(L, NA) = SQRT(2 * BOLTZK * T(L) / MASSNA + VMICCM2)
            VDOPDD(L, NA) = VDOPDD(L, NA) / 1.E5
            VDOPMIN = MIN(VDOPMIN, VDOPDD(L, NA))
            VDOPMAX = MAX(VDOPMAX, VDOPDD(L, NA))
          ELSE 
C***        Classic option: prescribed and depth-independent VDOP with value
C***                        similar to the spacing of the fine frequency grid
            VDOPDD(L, NA) = VDOPUNIT
          ENDIF
        ENDDO
      ENDDO
      
C***  Prevent undersampling due to lower VDOPUNIT resolution
      IF (VDOPMIN < VDOPUNIT) THEN
        WRITE (hCPR,*) '*** WARNING: MIN. VDOP SMALLER THAN VDOPUNIT!'
        WRITE (hCPR,'(2(A,F6.2))') '*** APPLIED VDOP WILL BE INCREASED '
     >     // 'FROM ', VDOPMIN, ' TO ', VDOPUNIT
        DO L=1, ND
          DO NA=1, NATOM
            VDOPDD(L, NA) = MAX(VDOPUNIT, VDOPDD(L, NA))
          ENDDO
        ENDDO
      ENDIF
      
      
      RETURN
      
      END
