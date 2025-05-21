  SUBROUTINE CLOUD_MODIF_LM (ZLM_CLOUD, D, TPFILE, PRT, PTKET, PDZZ, TZFIELD, ZLM, PZZ, ZETHETA, ZSHEAR, KRRI, IIJE, ZD,  &
  & ZEMOIST, CST, CSTURB, ZCOEF_AMPL, PTHVREF, ZDRTDZ, ZDTHLDZ, OOCEAN, IKT, IKU, ZPENTE, PTHLT, ZWORK2, GOCEAN, ZTHLM, ZRM, JK,  &
  & IKTE, ZLOCPEXNM, OCOMPUTE_SRC, PSRCT, PCOEF_AMPL_SAT, IKE, JIJ, ZAMOIST, IKA, ZCOEF_AMPL_CEI_NUL, ZALPHA, PDIRCOSZW, ZWORK1,  &
  & IKTB, PCEI, IKL, TURBN, PCEI_MIN, PDXX, ZVAR, IIJB, O2D, HTURBLEN_CL, PDYY, ZHOOK_HANDLE2, KRR, ZWORK2D, IKB, PCEI_MAX,  &
  & ZATHETA)
    !     #########################
    !!
    !!*****CLOUD_MODIF_LM routine to:
    !!       1/ change the mixing length in the clouds
    !!       2/ emphasize the mixing length in the cloud
    !!           by the coefficient ZCOEF_AMPL calculated here
    !!             when the CEI index is above ZCEI_MIN.
    !!
    !!
    !!      ZCOEF_AMPL ^
    !!                 |
    !!                 |
    !!  ZCOEF_AMPL_SAT -                       ---------- Saturation
    !!    (XDUMMY1)    |                      -
    !!                 |                     -
    !!                 |                    -
    !!                 |                   -
    !!                 |                  - Amplification
    !!                 |                 - straight
    !!                 |                - line
    !!                 |               -
    !!                 |              -
    !!                 |             -
    !!                 |            -
    !!                 |           -
    !!               1 ------------
    !!                 |
    !!                 |
    !!               0 -----------|------------|----------> PCEI
    !!                 0      ZCEI_MIN     ZCEI_MAX
    !!                        (XDUMMY2)    (XDUMMY3)
    !!
    !!
    !!
    !!    AUTHOR
    !!    ------
    !!     M. Tomasini   *CNRM METEO-FRANCE
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!     Original   09/07/04
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: CST_t
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    !
    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments
    !
    REAL, INTENT(INOUT) :: ZLM_CLOUD(D%NIJT, D%NKT)
    TYPE(DIMPHYEX_t), INTENT(IN) :: D
    TYPE(TFILEDATA), INTENT(INOUT) :: TPFILE
    REAL, INTENT(INOUT) :: PRT(D%NIJT, D%NKT, KRR)
    REAL, INTENT(IN) :: PTKET(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDZZ(D%NIJT, D%NKT)
    TYPE(TFIELDMETADATA), INTENT(INOUT) :: TZFIELD
    REAL, INTENT(INOUT) :: ZLM(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PZZ(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZETHETA(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZSHEAR(D%NIJT, D%NKT)
    INTEGER, INTENT(IN) :: KRRI
    INTEGER, INTENT(INOUT) :: IIJE
    REAL, INTENT(INOUT) :: ZD
    REAL, INTENT(INOUT) :: ZEMOIST(D%NIJT, D%NKT)
    TYPE(CST_t), INTENT(IN) :: CST
    TYPE(CSTURB_t), INTENT(IN) :: CSTURB
    REAL, INTENT(INOUT) :: ZCOEF_AMPL(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PTHVREF(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZDRTDZ(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZDTHLDZ(D%NIJT, D%NKT)
    LOGICAL, INTENT(IN) :: OOCEAN
    INTEGER, INTENT(INOUT) :: IKT
    INTEGER, INTENT(INOUT) :: IKU
    REAL, INTENT(INOUT) :: ZPENTE
    REAL, INTENT(INOUT) :: PTHLT(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZWORK2(D%NIJT, D%NKT)
    LOGICAL, INTENT(INOUT) :: GOCEAN
    REAL, INTENT(INOUT) :: ZTHLM(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZRM(D%NIJT, D%NKT, KRR)
    INTEGER, INTENT(INOUT) :: JK
    INTEGER, INTENT(INOUT) :: IKTE
    REAL, INTENT(INOUT) :: ZLOCPEXNM(D%NIJT, D%NKT)
    LOGICAL, INTENT(IN) :: OCOMPUTE_SRC
    REAL, INTENT(IN) :: PSRCT(MERGE(D%NIJT, 0, OCOMPUTE_SRC), MERGE(D%NKT, 0, OCOMPUTE_SRC))
    REAL, INTENT(IN) :: PCOEF_AMPL_SAT
    INTEGER, INTENT(INOUT) :: IKE
    INTEGER, INTENT(INOUT) :: JIJ
    REAL, INTENT(INOUT) :: ZAMOIST(D%NIJT, D%NKT)
    INTEGER, INTENT(INOUT) :: IKA
    REAL, INTENT(INOUT) :: ZCOEF_AMPL_CEI_NUL
    REAL, INTENT(INOUT) :: ZALPHA
    REAL, INTENT(IN) :: PDIRCOSZW(D%NIJT)
    REAL, INTENT(INOUT) :: ZWORK1(D%NIJT, D%NKT)
    INTEGER, INTENT(INOUT) :: IKTB
    REAL, INTENT(IN) :: PCEI(MERGE(D%NIJT, 0, OCLOUDMODIFLM), MERGE(D%NKT, 0, OCLOUDMODIFLM))
    INTEGER, INTENT(INOUT) :: IKL
    TYPE(TURB_t), INTENT(IN) :: TURBN
    REAL, INTENT(IN) :: PCEI_MIN
    REAL, INTENT(IN) :: PDXX(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: ZVAR
    INTEGER, INTENT(INOUT) :: IIJB
    LOGICAL, INTENT(IN) :: O2D
    CHARACTER(LEN=4), INTENT(IN) :: HTURBLEN_CL
    REAL, INTENT(IN) :: PDYY(D%NIJT, D%NKT)
    REAL(KIND=JPHOOK), INTENT(INOUT) :: ZHOOK_HANDLE2
    INTEGER, INTENT(IN) :: KRR
    REAL, INTENT(INOUT) :: ZWORK2D(D%NIJT)
    INTEGER, INTENT(INOUT) :: IKB
    REAL, INTENT(IN) :: PCEI_MAX
    REAL, INTENT(INOUT) :: ZATHETA(D%NIJT, D%NKT)
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    INITIALISATION
    !              --------------
    !
    IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM', 0, ZHOOK_HANDLE2)
    ZPENTE = (PCOEF_AMPL_SAT - 1.) / (PCEI_MAX - PCEI_MIN)
    ZCOEF_AMPL_CEI_NUL = 1. - ZPENTE*PCEI_MIN
    !
!$acc kernels
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
    ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = 1.
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !*       2.    CALCULATION OF THE AMPLIFICATION COEFFICIENT
    !              --------------------------------------------
    !
    ! Saturation
    !
!$acc kernels
!$mnh_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) >= PCEI_MAX)
      ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = PCOEF_AMPL_SAT
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    ! Between the min and max limits of CEI index, linear variation of the
    ! amplification coefficient ZCOEF_AMPL as a function of CEI
    !
!$acc kernels
!$mnh_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) < PCEI_MAX .and. PCEI(IIJB:IIJE, 1:IKT) > PCEI_MIN)
      ZCOEF_AMPL(IIJB:IIJE, 1:IKT) = ZPENTE*PCEI(IIJB:IIJE, 1:IKT) + ZCOEF_AMPL_CEI_NUL
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !
    !*       3.    CALCULATION OF THE MIXING LENGTH IN CLOUDS
    !              ------------------------------------------
    !
    IF (HTURBLEN_CL == TURBN%CTURBLEN) THEN
!$acc kernels
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
      ZLM_CLOUD(:, :) = ZLM(:, :)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    ELSE
      !
      !*         3.1 BL89 mixing length
      !           ------------------
      SELECT CASE (HTURBLEN_CL)
      CASE ('BL89', 'RM17', 'HM21')
!$acc kernels
!$mnh_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
        ZSHEAR(:, :) = 0.
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
        CALL BL89(D, CST, CSTURB, TURBN, PZZ, PDZZ, PTHVREF, ZTHLM, KRR, ZRM, PTKET, ZSHEAR, ZLM_CLOUD, OOCEAN)
        !
        !*         3.2 Delta mixing length
        !           -------------------
      CASE ('DELT')
        CALL DELT(PLM(D%NIJT, D%NKT), ODZ, ZWORK2, D, ZWORK1, IIJB, IKE, JIJ, O2D, PZZ(JIJ, IKTE + 1), PDYY, IKA, IKT, ZALPHA,  &
        & IKU, PDIRCOSZW(JIJ), IKTB, ZHOOK_HANDLE2, IIJE, ZD, GOCEAN, IKL, IKB, TURBN, JK, PDXX, IKTE, ODZ=.true.)
        !
        !*         3.3 Deardorff mixing length
        !           -----------------------
      CASE ('DEAR')
        CALL DEAR(PLM(D%NIJT, D%NKT), D, PRT(JIJ, JK + IKL, 1), PDZZ(JIJ, JK + IKL), PZZ(JIJ, IKTE + 1), PTKET(IIJB:IIJE, IKB),  &
        & ZETHETA(JIJ, JK), KRRI, IIJE, ZD, ZEMOIST(IIJB:IIJE, IKB), CST, PTHVREF(JIJ, JK), ZDRTDZ(IIJB:IIJE, IKB),  &
        & ZDTHLDZ(JIJ, JK), IKT, IKU, PTHLT(JIJ, JK + IKL), ZWORK2(IIJB:IIJE, 1:IKT), GOCEAN, JK, IKTE, ZLOCPEXNM, OCOMPUTE_SRC,  &
        & PSRCT, IKE, JIJ, ZAMOIST, IKA, ZALPHA, PDIRCOSZW(JIJ), ZWORK1(IIJB:IIJE, 1:IKT), IKTB, IKL, TURBN, PDXX, ZVAR, IIJB,  &
        & O2D, PDYY, ZHOOK_HANDLE2, KRR, ZWORK2D(IIJB:IIJE), IKB, ZATHETA)
        !
      END SELECT
    END IF
    !
    !*       4.    MODIFICATION OF THE MIXING LENGTH IN THE CLOUDS
    !              -----------------------------------------------
    !
    ! Impression before modification of the mixing length
    IF (TURBN%LTURB_DIAG .and. TPFILE%LOPENED) THEN
      TZFIELD = TFIELDMETADATA(CMNHNAME='LM_CLEAR_SKY', CSTDNAME='', CLONGNAME='LM_CLEAR_SKY', CUNITS='m', CDIR='XY', CCOMMENT= &
      & 'X_Y_Z_LM CLEAR SKY', NGRID=1, NTYPE=TYPEREAL, NDIMS=3, LTIMEDEP=.true.)
!$acc update self( ZLM )
      CALL IO_FIELD_WRITE_PHY(D, TPFILE, TZFIELD, ZLM)
    END IF
    !
    ! Amplification of the mixing length when the criteria are verified
    !
!$acc kernels
!$mnh_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (ZCOEF_AMPL(IIJB:IIJE, 1:IKT) /= 1.)
      ZLM(IIJB:IIJE, 1:IKT) = ZCOEF_AMPL(IIJB:IIJE, 1:IKT)*ZLM_CLOUD(IIJB:IIJE, 1:IKT)
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    ! Cloud mixing length in the clouds at the points which do not verified the CEI
    !
!$acc kernels
!$mnh_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
    WHERE (PCEI(IIJB:IIJE, 1:IKT) == -1.)
      ZLM(IIJB:IIJE, 1:IKT) = ZLM_CLOUD(IIJB:IIJE, 1:IKT)
    END WHERE
!$mnh_end_expand_where ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !
    !*       5.    IMPRESSION
    !              ----------
    !
    IF (TURBN%LTURB_DIAG .and. TPFILE%LOPENED) THEN
      TZFIELD = TFIELDMETADATA(CMNHNAME='COEF_AMPL', CSTDNAME='', CLONGNAME='COEF_AMPL', CUNITS='1', CDIR='XY', CCOMMENT= &
      & 'X_Y_Z_COEF AMPL', NGRID=1, NTYPE=TYPEREAL, NDIMS=3, LTIMEDEP=.true.)
!$acc update self( ZCOEF_AMPL )
      CALL IO_FIELD_WRITE_PHY(D, TPFILE, TZFIELD, ZCOEF_AMPL)
      !
      TZFIELD = TFIELDMETADATA(CMNHNAME='LM_CLOUD', CSTDNAME='', CLONGNAME='LM_CLOUD', CUNITS='m', CDIR='XY', CCOMMENT= &
      & 'X_Y_Z_LM CLOUD', NGRID=1, NTYPE=TYPEREAL, NDIMS=3, LTIMEDEP=.true.)
!$acc update self( ZLM_CLOUD )
      CALL IO_FIELD_WRITE_PHY(D, TPFILE, TZFIELD, ZLM_CLOUD)
      !
    END IF
    !
    IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM', 1, ZHOOK_HANDLE2)
  END SUBROUTINE CLOUD_MODIF_LM
  !
