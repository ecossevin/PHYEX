SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT (D, PALP, PBETA, PGAM, PLTT, PC, PT, PEXN, PCP, PLOCPEXN, PAMOIST, PATHETA, &
  & PDRVSATDT, PRVSAT, CST, PPABST)
    !     ########################################################################
    !!
    !!****  *COMPUTE_FUNCTION_THERMO* routine to compute several thermo functions
    !
    !!    AUTHOR
    !!    ------
    !!
    !!     JP Pinty      *LA*
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original   24/02/03
    !!     Modified: Wim de Rooy 06-02-2019
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    !
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: CST_t
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t

    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D  ! PHYEX variables dimensions structure
    REAL, INTENT(IN) :: PALP, PBETA, PGAM, PLTT, PC
    REAL, INTENT(IN), DIMENSION(D%NIJT, D%NKT) :: PT, PEXN, PCP
    !
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLOCPEXN
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PAMOIST, PATHETA
    REAL, INTENT(INOUT) :: PDRVSATDT(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PRVSAT(D%NIJT, D%NKT)
    TYPE(CST_t), INTENT(IN) :: CST
    REAL, INTENT(IN) :: PPABST(D%NIJT, D%NKT)

    REAL :: ZEPS
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2
    INTEGER :: JK, JIJ, IIJB, IIJE, IKT
    !
    !-------------------------------------------------------------------------------
    !
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT', 0, ZHOOK_HANDLE2)

    IKT = D%NKT
    IIJE = D%NIJE
    IIJB = D%NIJB
    ZEPS = CST%XMV / CST%XMD
    !
    !*       1.1 Lv/Cph at  t
    !
!$acc kernels
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
    PLOCPEXN(IIJB:IIJE, 1:IKT) = (PLTT + (CST%XCPV - PC)*(PT(IIJB:IIJE, 1:IKT) - CST%XTT)) / PCP(IIJB:IIJE, 1:IKT)
    !
    !*      1.2 Saturation vapor pressure at t
    !
    PRVSAT(IIJB:IIJE, 1:IKT) = EXP(PALP - PBETA / PT(IIJB:IIJE, 1:IKT) - PGAM*LOG(PT(IIJB:IIJE, 1:IKT)))
    !
    !*      1.3 saturation  mixing ratio at t
    !
    PRVSAT(IIJB:IIJE, 1:IKT) = PRVSAT(IIJB:IIJE, 1:IKT)*ZEPS / (PPABST(IIJB:IIJE, 1:IKT) - PRVSAT(IIJB:IIJE, 1:IKT))
    !
    !*      1.4 compute the saturation mixing ratio derivative (rvs')
    !
    PDRVSATDT(IIJB:IIJE, 1:IKT) = (PBETA / PT(IIJB:IIJE, 1:IKT) - PGAM) / PT(IIJB:IIJE, 1:IKT)*PRVSAT(IIJB:IIJE, 1:IKT)*(1. +  &
    & PRVSAT(IIJB:IIJE, 1:IKT) / ZEPS)
    !
    !*      1.5 compute Amoist
    !
    PAMOIST(IIJB:IIJE, 1:IKT) = 1.0 / (1.0 + PDRVSATDT(IIJB:IIJE, 1:IKT)*PLOCPEXN(IIJB:IIJE, 1:IKT))
    !
    !*      1.6 compute Atheta
    !
    PATHETA(IIJB:IIJE, 1:IKT) = PAMOIST(IIJB:IIJE, 1:IKT)*PEXN(IIJB:IIJE, 1:IKT)*PDRVSATDT(IIJB:IIJE, 1:IKT)
    !
    !*      1.7 Lv/Cph/Exner at t-1
    !
    PLOCPEXN(IIJB:IIJE, 1:IKT) = PLOCPEXN(IIJB:IIJE, 1:IKT) / PEXN(IIJB:IIJE, 1:IKT)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT', 1, ZHOOK_HANDLE2)
  END SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT


