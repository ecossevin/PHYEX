!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION_SPLIT
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                                   &OELEC, OSEDIM_BEARD, PTHVREFZIKB, &
                                   &PTSTEP, KRR, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                                   &PRS, PRT, &
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS,&
                                   &PEFIELDW, &
                                   &PSEA, PTOWN,  &
                                   &PINPRH, PFPR, &
                                   &PQHT, PQHS)
!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!                and modified for optimisation
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!  C. Barthe      03/2023: Add sedimentation of electric charges
!
!
!*      0. DECLARATIONS
!          ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_CST, ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_ELEC_DESCR,     ONLY: ELEC_DESCR_t
USE MODD_ELEC_PARAM,     ONLY: ELEC_PARAM_t
USE MODD_FIELDS_ADDRESS
!
!USE MODI_GAMMA, ONLY: GAMMA
#ifndef MNH_COMPILER_CCE
USE MODI_GAMMA
#endif
#if defined(TARGET_NV70)
USE MODI_GAMMA
#endif
!
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_ELEC_BEARD_EFFECT, ONLY: ELEC_BEARD_EFFECT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),              INTENT(IN)              :: D       !array dimensions
TYPE(CST_t),                   INTENT(IN)              :: CST
TYPE(RAIN_ICE_PARAM_t),        INTENT(IN)              :: ICEP
TYPE(RAIN_ICE_DESCR_t),        INTENT(IN)              :: ICED
TYPE(PARAM_ICE_t),             INTENT(IN)              :: PARAMI
TYPE(ELEC_PARAM_t),            INTENT(IN)              :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),            INTENT(IN)              :: ELECD   ! electrical descriptive csts
LOGICAL,                       INTENT(IN)              :: OELEC   ! if true, cloud electricity is activated
LOGICAL,                       INTENT(IN)              :: OSEDIM_BEARD ! if true, effect of electrical forces on sedim.
REAL,                          INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                       INTENT(IN)              :: KRR     ! Number of moist variable
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PT      ! Temperature at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT)       :: PRS     ! m.r. source
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)          :: PRT     ! m.r. at t
REAL, DIMENSION(D%NIJT),       INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(D%NIJT),       INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(D%NIJT),       INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(D%NIJT),       INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(D%NIJT),       INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
! variables for cloud electricity
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PQCT   ! Cloud water electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PQRT   ! Rain water electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PQIT   ! Pristine ice electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PQST   ! Snow/aggregate electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PQGT   ! Graupel electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(IN)    :: PQHT   ! Hail electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQCS   ! Cloud water electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQRS   ! Rain water electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQIS   ! Pristine ice electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQSS   ! Snow/aggregate electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQGS   ! Graupel electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHS   ! Hail electric charge source
REAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)), INTENT(IN) :: PEFIELDW ! Vertical E field 
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!*       0.2  declaration of local variables
!
!
INTEGER :: JIJ, JK
INTEGER :: IKTB, IKTE, IKB, IKL, IIJE, IIJB
INTEGER :: IKRR !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL :: GSEDIC !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL :: GPRESENT_PFPR, GPRESENT_PSEA
REAL    :: ZINVTSTEP
REAL, DIMENSION(D%NIJT)               :: ZCONC_TMP    ! Weighted concentration
REAL, DIMENSION(D%NIJT,D%NKTB:D%NKTE) :: ZW ! work array
REAL, DIMENSION(D%NIJT, D%NKT)        :: ZCONC3D, & !  droplet condensation
                                       & ZRAY,   & ! Cloud Mean radius
                                       & ZLBC,   & ! XLBC weighted by sea fraction
                                       & ZFSEDC
REAL, DIMENSION(D%NIJT, D%NKT, KRR) :: ZRS, &   ! Mixing ratios created during the time step
                                     & ZRT
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)) :: &
                   ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQHT, &      ! electric charge a t
                   ZPQCS, ZPQRS, ZPQIS, ZPQSS, ZPQGS, ZPQHS   ! electric charge created during the time step
CHARACTER (LEN=4) :: HCLOUD  ! Kind of microphysical scheme
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JRR
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT', 0, ZHOOK_HANDLE)

!-------------------------------------------------------------------------------
!
!
GSEDIC = PARAMI%LSEDIC
IKRR    = KRR
!
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
!
IF (PRESENT(PFPR)) THEN
  GPRESENT_PFPR = .TRUE.
ELSE
  GPRESENT_PFPR = .FALSE.
END IF
!
IF (PRESENT(PSEA)) THEN
  GPRESENT_PSEA = .TRUE.
ELSE
  GPRESENT_PSEA = .FALSE.
END IF
!
IF(IKRR==7) THEN
  HCLOUD='ICE4'
ELSE
  HCLOUD='ICE3'
ENDIF
!
!        O. Initialization of for sedimentation
!
!$acc kernels
ZINVTSTEP=1./PTSTEP
!$acc end kernels
IF (GPRESENT_PFPR) THEN
!$acc kernels
  PFPR(:,:,:) = 0.
!$acc end kernels
END IF
!
!*       1. Parameters for cloud sedimentation
!
IF (GSEDIC) THEN
!$acc kernels
  ZRAY(:,:)   = 0.
  ZLBC(:,:)   = ICED%XLBC(1)
  ZFSEDC(:,:) = ICEP%XFSEDC(1)
  ZCONC3D(:,:)= ICED%XCONC_LAND
  ZCONC_TMP(:)= ICED%XCONC_LAND
!$acc end kernels
  IF (GPRESENT_PSEA) THEN
!$acc kernels
!$acc loop independent
    DO JIJ = IIJB, IIJE
      ZCONC_TMP(JIJ)=PSEA(JIJ)*ICED%XCONC_SEA+(1.-PSEA(JIJ))*ICED%XCONC_LAND
    ENDDO
!$acc end kernels

!$acc kernels
!$acc loop independent collapse(2)
    DO JK=IKTB, IKTE
      DO JIJ = IIJB, IIJE
          ZLBC(JIJ,JK)   = PSEA(JIJ)*ICED%XLBC(2)+(1.-PSEA(JIJ))*ICED%XLBC(1)
          ZFSEDC(JIJ,JK) = (PSEA(JIJ)*ICEP%XFSEDC(2)+(1.-PSEA(JIJ))*ICEP%XFSEDC(1))
          ZFSEDC(JIJ,JK) = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)),ZFSEDC(JIJ,JK))
          ZCONC3D(JIJ,JK)= (1.-PTOWN(JIJ))*ZCONC_TMP(JIJ)+PTOWN(JIJ)*ICED%XCONC_URBAN
          ZRAY(JIJ,JK)   = 0.5*((1.-PSEA(JIJ))*GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)) + &
                         & PSEA(JIJ)*GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)/(GAMMA(ICED%XNUC2)))
      ENDDO
    END DO
!$acc end kernels
  ELSE
!$acc kernels
    ZCONC3D(:,:) = ICED%XCONC_LAND
    ZRAY(:,:)  = 0.5*(GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)))
!$acc end kernels
  END IF
!$acc kernels
!$acc loop independent collapse(2)
  DO JK=IKTB, IKTE
    DO JIJ = IIJB, IIJE
      ZRAY(JIJ,JK)      = MAX(1.,ZRAY(JIJ,JK))
      ZLBC(JIJ,JK)      = MAX(MIN(ICED%XLBC(1),ICED%XLBC(2)),ZLBC(JIJ,JK))
    ENDDO
  ENDDO
!$acc end kernels
ENDIF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately
!
!$acc kernels
!$acc loop independent collapse(2)
DO JK=IKTB, IKTE
  DO JIJ = IIJB, IIJE
    ! External tendecies
    IF (GSEDIC) THEN
      ZRS(JIJ,JK,IRC) = PRS(JIJ,JK,IRC)-PRT(JIJ,JK,IRC)*ZINVTSTEP
    ENDIF
    DO JRR=IRR, IKRR
      ZRS(JIJ,JK,JRR) = PRS(JIJ,JK,JRR)-PRT(JIJ,JK,JRR)*ZINVTSTEP
    ENDDO
    !
    ! mr values inside the time-splitting loop
    DO JRR=IRC, IKRR
      ZRT(JIJ,JK,JRR) = PRT(JIJ,JK,JRR)
    ENDDO
    !
    ZW(JIJ,JK) =1./(PRHODREF(JIJ,JK)* PDZZ(JIJ,JK))
    !
    ! Cloud electricity
    IF (OELEC) THEN
      IF (GSEDIC) ZPQCS(JIJ,JK) = PQCS(JIJ,JK) - PQCT(JIJ,JK) * ZINVTSTEP
      ZPQRS(JIJ,JK) = PQRS(JIJ,JK) - PQRT(JIJ,JK) * ZINVTSTEP
      ZPQIS(JIJ,JK) = PQIS(JIJ,JK) - PQIT(JIJ,JK) * ZINVTSTEP
      ZPQSS(JIJ,JK) = PQSS(JIJ,JK) - PQST(JIJ,JK) * ZINVTSTEP
      ZPQGS(JIJ,JK) = PQGS(JIJ,JK) - PQGT(JIJ,JK) * ZINVTSTEP
      IF (IKRR==7) ZPQHS(JIJ,JK) = PQHS(JIJ,JK) - PQHT(JIJ,JK) * ZINVTSTEP
      !
      ZQCT(JIJ,JK) = PQCT(JIJ,JK)
      ZQRT(JIJ,JK) = PQST(JIJ,JK)
      ZQIT(JIJ,JK) = PQIT(JIJ,JK)
      ZQST(JIJ,JK) = PQST(JIJ,JK)
      ZQGT(JIJ,JK) = PQGT(JIJ,JK)
      IF (IKRR==7) ZQHT(JIJ,JK) = PQHT(JIJ,JK)
    ENDIF
  ENDDO
ENDDO
!$acc end kernels
!
!
!*       2.1   for cloud
!
IF (GSEDIC) THEN
    CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                          &IKRR, OELEC, OSEDIM_BEARD, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                          &IRC, &
                          &ZRT(:,:,IRC), PRS(:,:,IRC), PINPRC, ZRS(:,:,IRC), &
                          &ZQCT, PQCS, ZPQCS, PEFIELDW, &
                          &PRAY=ZRAY, PLBC=ZLBC, PFSEDC=ZFSEDC, PCONC3D=ZCONC3D, &
                          &PFPR=PFPR)
ENDIF
!
!*       2.2   for rain
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                          &IKRR, OELEC, OSEDIM_BEARD, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                          &IRR, &
                          &ZRT(:,:,IRR), PRS(:,:,IRR), PINPRR, ZRS(:,:,IRR), &
                          &ZQRT, PQRS, ZPQRS, PEFIELDW, &
                          &PFPR=PFPR)
!
!*       2.3   for pristine ice
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI,  ELECP, ELECD, &
                          &IKRR, OELEC, OSEDIM_BEARD, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                          &IRI, &
                          &ZRT(:,:,IRI), PRS(:,:,IRI), PINPRI, ZRS(:,:,IRI), &
                          &ZQIT, PQIS, ZPQIS, PEFIELDW, &
                          &PFPR=PFPR)
!
!*       2.4   for aggregates/snow
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI,  ELECP, ELECD, &
                          &IKRR, OELEC, OSEDIM_BEARD, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                          &IRS, &
                          &ZRT(:,:,IRS), PRS(:,:,IRS), PINPRS, ZRS(:,:,IRS), &
                          &ZQST, PQSS, ZPQSS, PEFIELDW, &
                          &PFPR=PFPR)
!
!*       2.5   for graupeln
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI,  ELECP, ELECD, &
                          &IKRR, OELEC, OSEDIM_BEARD, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                          &IRG, &
                          &ZRT(:,:,IRG), PRS(:,:,IRG), PINPRG, ZRS(:,:,IRG), &
                          &ZQGT, PQGS, ZPQGS, PEFIELDW, &
                          &PFPR=PFPR)
!
!*       2.6   for hail
!
IF (IKRR==7) THEN
    CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                            &IKRR, OELEC, OSEDIM_BEARD, &
                            &PRHODREF, ZW, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                            &IRH, &
                            &ZRT(:,:,IRH), PRS(:,:,IRH), PINPRH, ZRS(:,:,IRH), &
                            &ZQHT, PQHS, ZPQHS, PEFIELDW, &
                            &PFPR=PFPR)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT', 1, ZHOOK_HANDLE)
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                              &KRR, OELEC, OSEDIM_BEARD, &
                              &PRHODREF, POORHODZ, PDZZ, PPABST, PTHT, PT, PTSTEP, &
                              &KSPE, &
                              &PRXT, PRXS, PINPRX, PPRXS, &
                              &PQXT, PQXS, PPQXS, PEFIELDW, &
                              &PRAY, PLBC, PFSEDC, PCONC3D, PFPR)
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
!
! parameters for electricity
USE MODD_ELEC_PARAM,     ONLY: ELEC_PARAM_t
USE MODD_ELEC_DESCR,     ONLY: ELEC_DESCR_t
!
USE MODI_MOMG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),             INTENT(IN)              :: D
TYPE(CST_t),                  INTENT(IN)              :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)              :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)              :: ICED
TYPE(PARAM_ICE_t),            INTENT(IN)              :: PARAMI
TYPE(ELEC_PARAM_t),           INTENT(IN)              :: ELECP        ! electrical parameters
TYPE(ELEC_DESCR_t),           INTENT(IN)              :: ELECD        ! electrical descriptive csts
INTEGER,                      INTENT(IN)              :: KRR
LOGICAL,                      INTENT(IN)              :: OELEC        ! if true, sedimentation of elec. charges
LOGICAL,                      INTENT(IN)              :: OSEDIM_BEARD ! if true, effect of electric forces on sedim.
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT,D%NKTB:D%NKTE), INTENT(IN)     :: POORHODZ ! One Over (Rhodref times delta Z)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PDZZ ! layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PPABST
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PTHT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PT
REAL,                          INTENT(IN)             :: PTSTEP  ! total timestep
INTEGER,                       INTENT(IN)             :: KSPE ! 1 for rc, 2 for rr...
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT)          :: PRXT ! mr of specy X
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT)          :: PRXS !Tendency of the specy KSPE
REAL, DIMENSION(D%NIJT),       INTENT(OUT)            :: PINPRX ! instant precip
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)             :: PPRXS ! external tendencie
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQXT  ! electric charge at t for specy KSPE
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQXS  ! tendency of the electric charge 
                                                                                    ! for specy KSPE
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN)    :: PPQXS ! external tendency
REAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)), INTENT(IN) :: PEFIELDW ! Vertical E field
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN), OPTIONAL    :: PRAY, PLBC, PFSEDC, PCONC3D
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT), OPTIONAL :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
INTEGER                         :: JI, JIJ, JK
LOGICAL                         :: GPRESENT_PFPR
REAL                            :: ZINVTSTEP
REAL                            :: ZZWLBDC, ZZRAY, ZZT, ZZWLBDA, ZZCC
REAL                            :: ZLBDA
REAL                            :: ZFSED, ZEXSED
REAL                            :: ZMRCHANGE
REAL, DIMENSION(D%NIJT)       :: ZMAX_TSTEP1D ! Maximum CFL in column
REAL, DIMENSION(D%NIJT,D%NKT)       :: ZMAX_TSTEP2D ! Maximum CFL in column
REAL, DIMENSION(SIZE(ICED%XRTMIN))   :: ZRSMIN
REAL, DIMENSION(D%NIJT)       :: ZREMAINT   ! Remaining time until the timestep end
LOGICAL :: ZANYREMAINT
REAL, DIMENSION(D%NIJT, 0:D%NKT+1) :: ZWSED   ! Sedimentation fluxes
INTEGER :: IKTB, IKTE, IKB, IKL, IIJE, IIJB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! local variables for cloud electricity
REAL :: ZEXT   ! e_x coefficient of the q(D) relation
REAL :: ZNCI   ! N_ci for ice crystal sedimentation
!REAL,    DIMENSION(D%NIJT,0:D%NKT+1) :: ZWSEDQ ! Sedimentation fluxes for electric charges
!REAL,    DIMENSION(D%NIJT,0:D%NKT+1) :: ZBEARDCOEFF ! effect of electric forces on sedimentation
!REAL,    DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(0:D%NKT+1,0,OELEC)) :: &
!++cb-- est-ce que cette declaration est correcte par rapport a ce qui est fait pour zwsed ?
REAL,    DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT+2,0,OELEC)) :: &
                                     ZWSEDQ      ! Sedimentation fluxes for electric charges
!REAL,    DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(0:D%NKT+1,0,OSEDIM_BEARD)) :: &
REAL,    DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)) :: &
                                     ZBEARDCOEFF ! effect of electric forces on sedimentation
!REAL,    DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(0:D%NKT+1,0,OSEDIM_BEARD)) :: &
REAL,    DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)) :: &
                                     ZLBDA3      ! slope parameter of the distribution
!LOGICAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(0:D%NKT+1,0,OSEDIM_BEARD)) :: GMASK
LOGICAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)) :: GMASK
REAL :: ZQCHANGE
REAL :: ZFQSED, ZEXQSED
REAL :: ZEXMIN, ZEXMAX
REAL :: ZLBX, ZLBEXX
REAL :: ZFQUPDX
REAL :: ZCXX
REAL :: ZFX
! end - local variables for cloud electricity
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 0, ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKL=D%NKL
IIJB=D%NIJB
IIJE=D%NIJE
!-------------------------------------------------------------------------------
IF (KSPE<2 .OR. KSPE>7) CALL PRINT_MSG(NVERB_FATAL,'GEN','INTERNAL_SEDIM_SPLIT','invalid species (KSPE variable)')
!
IF (PRESENT(PFPR)) THEN
  GPRESENT_PFPR = .TRUE.
ELSE
  GPRESENT_PFPR = .FALSE.
END IF
!
!$acc kernels
PINPRX(:) = 0.
ZINVTSTEP=1./PTSTEP
DO JI = 1, SIZE(ICED%XRTMIN)
  ZRSMIN(JI) = ICED%XRTMIN(JI) * ZINVTSTEP
END DO
ZREMAINT(:) = 0.
ZREMAINT(IIJB:IIJE) = PTSTEP
!$acc end kernels
!$acc update self(ZREMAINT)
!
ZANYREMAINT = .TRUE.
DO WHILE (ZANYREMAINT)
  !
  ! Effect of electrical forces on sedimentation
  IF (OELEC .AND. OSEDIM_BEARD) THEN
    DO JK = IKTB, IKTE
      DO JIJ = IIJB, IIJE
        IF (PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZREMAINT(JIJ)>0.) THEN
          GMASK(JIJ,JK) = .TRUE.
        ELSE
          GMASK(JIJ,JK) = .FALSE.
        END IF
      END DO
    END DO
  END IF
  !
  !*       1. Parameters for cloud sedimentation
  !
  !
  !*       2.    compute the fluxes
  !
  !
  IF(KSPE==IRC) THEN
    !******* for cloud
!$acc kernels
    ZWSED(:,:) = 0.
    IF (OELEC) THEN
      ZWSEDQ(:,:) = 0.
      ZLBDA3(:,:) = 0.
    END IF
!$acc end kernels
!$acc kernels
!$acc loop independent collapse(2) private(ZEXT)
    DO JK = IKTB,IKTE
      DO JIJ = IIJB,IIJE
        IF(PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZREMAINT(JIJ)>0.) THEN
          ZZWLBDC = PLBC(JIJ,JK) * PCONC3D(JIJ,JK) / &
                   &(PRHODREF(JIJ,JK) * PRXT(JIJ,JK))
          ZZWLBDC = ZZWLBDC**ICED%XLBEXC
          ZZRAY = PRAY(JIJ,JK) / ZZWLBDC
          ZZT = PTHT(JIJ,JK) * (PPABST(JIJ,JK)/CST%XP00)**(CST%XRD/CST%XCPD)
          ZZWLBDA = 6.6E-8*(101325./PPABST(JIJ,JK))*(ZZT/293.15)
          ZZCC = ICED%XCC*(1.+1.26*ZZWLBDA/ZZRAY)
          ZWSED(JIJ, JK) = PRHODREF(JIJ,JK)**(-ICED%XCEXVT +1 ) *   &
                             &ZZWLBDC**(-ICED%XDC)*ZZCC*PFSEDC(JIJ,JK) * PRXT(JIJ,JK)
!++cb++ nouveau : traitement de la sedimentation des charges portees par les gouttelettes
! A TESTER
          IF (OELEC) THEN
            ZEXT = PQXT(JIJ,JK) / ELECP%XFQUPDC * PRHODREF(JIJ,JK)
            IF (ABS(ZEXT) .GT. ELECP%XECMIN) THEN
              ZWSEDQ(JIJ,JK) = ELECP%XFQSEDC * ZEXT * PCONC3D(JIJ,JK) * &
                               PRHODREF(JIJ,JK)**(-ICED%XCEXVT) *       &
                               ZZCC * ZZWLBDC**(-ELECP%XEXQSEDC)
              IF (OSEDIM_BEARD) ZLBDA3(JIJ,JK) = ZZWLBDC
            ENDIF
          ENDIF
        ENDIF
!--cb--
      ENDDO
    ENDDO
!$acc end kernels
    IF (OELEC .AND. OSEDIM_BEARD) THEN
      CALL ELEC_BEARD_EFFECT(D, CST, HCLOUD, KSPE, GMASK, PT, PRHODREF, PTHVREFZIKB, &
                             PRXT, PQXT, PEFIELDW, ZLBDA3, ZBEARDCOEFF, ICED=ICED)
!$acc kernels
!$acc loop independent collapse(2)
      DO JK = IKTB,IKTE
        DO JIJ = IIJB,IIJE
          ZWSED(JIJ,JK)  = ZWSED(JIJ,JK)  * ZBEARDCOEFF(JIJ,JK)
          ZWSEDQ(JIJ,JK) = ZWSEDQ(JIJ,JK) * ZBEARDCOEFF(JIJ,JK)
        END DO
      END DO
!$acc end kernels
    END IF
  ELSEIF(KSPE==IRI) THEN
    ! ******* for pristine ice
!$acc kernels
    ZWSED(:,:) = 0.
    IF (OELEC) THEN
      ZWSEDQ(:,:) = 0.
      ZLBDA3(:,:) = 0.
    END IF
!$acc end kernels
!$acc kernels
!$acc loop independent collapse(2) private(ZEXT)
    DO JK = IKTB,IKTE
      DO JIJ = IIJB,IIJE
        IF(PRXT(JIJ, JK) .GT. MAX(ICED%XRTMIN(4), 1.0E-7) .AND. ZREMAINT(JIJ)>0.) THEN
          ZWSED(JIJ, JK) =  ICEP%XFSEDI * PRXT(JIJ, JK) *  &
                              & PRHODREF(JIJ,JK)**(1.-ICED%XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      LOG(PRHODREF(JIJ,JK)*PRXT(JIJ,JK)) )**ICEP%XEXCSEDI
          IF (OELEC) THEN
            ! N_ci from McF&H
            ZNCI = ELECP%XFCI * PRHODREF(JIJ,JK) * PRXT(JIJ,JK) * &
                          MAX(0.05E6,-0.15319E6-0.021454E6*LOG(PRHODREF(JIJ,JK)*PRXT(JIJ,JK)))**3.
            ! compute e_i of the q - D relationship
            ZEXT = PQXT(JIJ,JK) / ELECP%XFQUPDI *                         &
                  (PRHODREF(JIJ,JK) * PRXT(JIJ,JK))**(-ELECP%XEXFQUPDI) * &
                   ZNCI**(ELECP%XEXFQUPDI-1.)
            IF (ABS(ZEXT) .GT. ELECP%XEIMIN) THEN
              ZWSEDQ(JIJ,JK) = ELECP%XFQSEDI * ZEXT * PRXT(JIJ,JK) * &
                               PRHODREF(JIJ,JK)**(1.-ICED%XCEXVT) *  &
                               MAX( 0.05E6,-0.15319E6-0.021454E6*    & !    McF&H
                                    LOG(PRHODREF(JIJ,JK)*PRXT(JIJ,JK)) )**(3.*(1-ELECP%XEXQSEDI))
              IF (OSEDIM_BEARD) ZLBDA3(JIJ,JK) = (2.14E-3 * MOMG(ICED%XALPHAI,ICED%XNUI,1.7) *     &
                                                  ZNCI / (PRHODREF(JIJ,JK) * PRXT(JIJ,JK)))**0.588235
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
!$acc end kernels
    IF (OELEC .AND. OSEDIM_BEARD) THEN
      CALL ELEC_BEARD_EFFECT(D, CST, HCLOUD, KSPE, GMASK, PT, PRHODREF, PTHVREFZIKB, &
                             PRXT, PQXT, PEFIELDW, ZLBDA3, ZBEARDCOEFF, ICED=ICED)
!$acc kernels
!$acc loop independent collapse(2)
      DO JK = IKTB,IKTE
        DO JIJ = IIJB,IIJE
          ZWSED(JIJ,JK)  = ZWSED(JIJ,JK)  * ZBEARDCOEFF(JIJ,JK)
          ZWSEDQ(JIJ,JK) = ZWSEDQ(JIJ,JK) * ZBEARDCOEFF(JIJ,JK)
        END DO
      END DO
!$acc end kernels
    END IF
  ELSEIF(KSPE==IRS) THEN
    ! ******* for snow
!$acc kernels
    ZWSED(:,:) = 0.
    IF (OELEC) THEN
      ZWSEDQ(:,:) = 0.
      ZLBDA3(:,:) = 0.
    END IF
!$acc end kernels
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      !The following lines must be kept equal to the computation in the general case ("for other species" case below)
      ZFSED=ICEP%XFSEDS
      ZEXSED=ICEP%XEXSEDS
!$acc kernels
!$acc loop independent collapse(2)
      DO JK = IKTB,IKTE
        DO JIJ = IIJB,IIJE
          IF(PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZREMAINT(JIJ)>0.) THEN
            ZWSED(JIJ, JK) = ZFSED  * PRXT(JIJ, JK)**ZEXSED            &
                           &        * PRHODREF(JIJ, JK)**(ZEXSED-ICED%XCEXVT)
          ENDIF
        ENDDO
      ENDDO
!$acc end kernels
    ELSE
!$acc kernels
!$acc loop independent collapse(2) private(ZEXT)
      DO JK = IKTB,IKTE
        DO JIJ = IIJB,IIJE
          IF(PRXT(JIJ,JK)> ICED%XRTMIN(KSPE) .AND. ZREMAINT(JIJ)>0.) THEN
            IF (PARAMI%LSNOW_T .AND. PT(JIJ,JK)>263.15) THEN
              ZLBDA = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*PT(JIJ,JK))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
            ELSE IF (PARAMI%LSNOW_T) THEN
              ZLBDA = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226 -0.0106*PT(JIJ,JK))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
            ELSE
              ZLBDA=MAX(MIN(ICED%XLBDAS_MAX, ICED%XLBS * ( PRHODREF(JIJ,JK) * PRXT(JIJ,JK) )**ICED%XLBEXS),ICED%XLBDAS_MIN)
            END IF
            ZWSED(JIJ, JK) = ICEP%XFSEDS *  &
                 & PRXT(JIJ,JK)* &
                 & PRHODREF(JIJ,JK)**(1-ICED%XCEXVT) * &
                 & (1 + (ICED%XFVELOS/ZLBDA)**ICED%XALPHAS)** (-ICED%XNUS+ICEP%XEXSEDS/ICED%XALPHAS) * &
                 & ZLBDA ** (ICED%XBS+ICEP%XEXSEDS)
            IF (OELEC .AND. ZLBDA > 0.) THEN
              ! compute the e_x coefficient of the q - D relationship
              ZEXT = PRHODREF(JIJ,JK) * PQXT(JIJ,JK) / (ELECP%XFQUPDS * ZLBDA**(ICED%XCXS-ELECD%XFS))
              ZEXT = SIGN( MIN(ABS(ZEXT), ELECP%XESMAX), ZEXT)
              IF (ABS(ZEXT) > ELECP%XESMIN) THEN
                ZWSEDQ(JIJ,JK) = ELECP%XFQSEDS * ZEXT * PRXT(JIJ,JK)**ELECP%XEXQSEDS &
                                         * PRHODREF(JIJ,JK)**(ELECP%XEXQSEDS-ICED%XCEXVT)
                IF (OSEDIM_BEARD) ZLBDA3(JIJ,JK) = ZLBDA
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!$acc end kernels
      IF (OELEC .AND. OSEDIM_BEARD) THEN
        CALL ELEC_BEARD_EFFECT(D, CST, HCLOUD, KSPE, GMASK, PT, PRHODREF, PTHVREFZIKB,&
                               PRXT, PQXT, PEFIELDW, ZLBDA3, ZBEARDCOEFF, ICED=ICED)
!$acc kernels
!$acc loop independent collapse(2)
        DO JK = IKTB,IKTE
          DO JIJ = IIJB,IIJE
            ZWSED(JIJ,JK)  = ZWSED(JIJ,JK)  * ZBEARDCOEFF(JIJ,JK)
            ZWSEDQ(JIJ,JK) = ZWSEDQ(JIJ,JK) * ZBEARDCOEFF(JIJ,JK)
          END DO
        END DO
!$acc end kernels
      END IF
    ENDIF
  ELSE
    ! ******* for other species
    SELECT CASE(KSPE)
      CASE(IRR)
        ZFSED=ICEP%XFSEDR
        ZEXSED=ICEP%XEXSEDR
      CASE(IRG)
        ZFSED=ICEP%XFSEDG
        ZEXSED=ICEP%XEXSEDG
      CASE(IRH)
        ZFSED=ICEP%XFSEDH
        ZEXSED=ICEP%XEXSEDH
      CASE DEFAULT
        CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT', 'no sedimentation parameter for KSPE')
    END SELECT
    !
    IF (OELEC) THEN
      SELECT CASE(KSPE)
        CASE(IRR)
          ZFQSED = ELECP%XFQSEDR
          ZEXQSED = ELECP%XEXQSEDR
          ZEXMIN = ELECP%XERMIN
          ZEXMAX = ELECP%XERMAX
          ZLBX = ICED%XLBR
          ZLBEXX = ICED%XLBEXR
          ZFQUPDX = ELECP%XFQUPDR
          ZCXX = ELECD%XCXR
          ZFX = ELECD%XFR
        CASE(IRG)
          ZFQSED = ELECP%XFQSEDG
          ZEXQSED = ELECP%XEXQSEDG
          ZEXMIN = ELECP%XEGMIN
          ZEXMAX = ELECP%XEGMAX
          ZLBX = ICED%XLBG
          ZLBEXX = ICED%XLBEXG
          ZFQUPDX = ELECP%XFQUPDG
          ZCXX = ICED%XCXG
          ZFX = ELECD%XFG
        CASE(IRH)
          ZFQSED = ELECP%XFQSEDH
          ZEXQSED = ELECP%XEXQSEDH
          ZEXMIN = ELECP%XEHMIN
          ZEXMAX = ELECP%XEHMAX
          ZLBX = ICED%XLBH
          ZLBEXX = ICED%XLBEXH
          ZFQUPDX = ELECP%XFQUPDH
          ZCXX = ICED%XCXH
          ZFX = ELECD%XFH
        CASE DEFAULT
          CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT', 'no sedimentation parameter for KSPE')
      END SELECT
    END IF
    !
!$acc kernels
    ZWSED(:,:) = 0.
    IF (OELEC) THEN
      ZWSEDQ(:,:) = 0.
      ZLBDA3(:,:) = 0.
    END IF
!$acc end kernels
!$acc kernels
!$acc loop independent collapse(2) private(ZEXT)
    DO JK = IKTB,IKTE
      DO JIJ = IIJB,IIJE
        IF(PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZREMAINT(JIJ)>0.) THEN
          ZWSED(JIJ, JK) = ZFSED  * PRXT(JIJ, JK)**ZEXSED            &
                         &        * PRHODREF(JIJ, JK)**(ZEXSED-ICED%XCEXVT)
          IF (OELEC) THEN
            ! need lambda_x to compute e_x
            ZLBDA = ZLBX * (PRHODREF(JIJ,JK) * MAX(PRXT(JIJ,JK), ICED%XRTMIN(KSPE)))**ZLBEXX
            IF (ZLBDA > 0.) THEN
              ! compute the e_x coefficient of the q - D relationship
              ZEXT = PRHODREF(JIJ,JK) * PQXT(JIJ,JK) / (ZFQUPDX * ZLBDA**(ZCXX-ZFX))
              ZEXT = SIGN( MIN(ABS(ZEXT), ZEXMAX), ZEXT)
            END IF
            IF (ABS(ZEXT) > ZEXMIN) THEN
              ZWSEDQ(JIJ,JK) = ZFQSED * ZEXT * PRXT(JIJ,JK)**ZEXQSED &
                                      * PRHODREF(JIJ,JK)**(ZEXQSED-ICED%XCEXVT)
              IF (OSEDIM_BEARD) ZLBDA3(JIJ,JK) = ZLBDA
            END IF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
!$acc end kernels
    IF (OELEC .AND. OSEDIM_BEARD) THEN
      CALL ELEC_BEARD_EFFECT(D, CST, HCLOUD, KSPE, GMASK, PT, PRHODREF, PTHVREFZIKB, &
                             PRXT, PQXT, PEFIELDW, ZLBDA3, ZBEARDCOEFF, ICED=ICED)
!$acc kernels
!$acc loop independent collapse(2)
      DO JK = IKTB,IKTE
        DO JIJ = IIJB,IIJE
          ZWSED(JIJ,JK)  = ZWSED(JIJ,JK)  * ZBEARDCOEFF(JIJ,JK)
          ZWSEDQ(JIJ,JK) = ZWSEDQ(JIJ,JK) * ZBEARDCOEFF(JIJ,JK)
        END DO
      END DO
!$acc end kernels
    END IF
  ENDIF
!$acc kernels
  ZMAX_TSTEP1D(:) = ZREMAINT(:)

!$acc loop independent collapse(2)
  DO JK = IKTB,IKTE
    DO JIJ = IIJB,IIJE
      ZMAX_TSTEP2D(JIJ,JK) = ZREMAINT(JIJ)
    END DO
  END DO

!$acc loop independent collapse(2)
  DO JK = IKTB,IKTE
    DO JIJ = IIJB,IIJE
      IF(PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZWSED(JIJ, JK)>1.E-20 .AND. ZREMAINT(JIJ)>0.) THEN
        ZMAX_TSTEP2D(JIJ,JK) = PARAMI%XSPLIT_MAXCFL * PRHODREF(JIJ, JK) * &
                        & PRXT(JIJ, JK) * PDZZ(JIJ, JK) / ZWSED(JIJ, JK)
      ENDIF
    ENDDO
  ENDDO

!$acc loop independent collapse(2)
  DO JK = IKTB,IKTE
    DO JIJ = IIJB,IIJE
      IF(PRXT(JIJ,JK)>ICED%XRTMIN(KSPE) .AND. ZWSED(JIJ, JK)>1.E-20 .AND. ZREMAINT(JIJ)>0.) THEN
        ZMAX_TSTEP1D(JIJ) = MIN(ZMAX_TSTEP1D(JIJ), ZMAX_TSTEP2D(JIJ,JK))
      ENDIF
    ENDDO
  ENDDO
!$acc end kernels
!$acc kernels
!$acc loop independent
  DO JIJ = IIJB, IIJE
      ZREMAINT(JIJ) = ZREMAINT(JIJ) - ZMAX_TSTEP1D(JIJ)
      PINPRX(JIJ) = PINPRX(JIJ) + ZWSED(JIJ,IKB) / CST%XRHOLW * (ZMAX_TSTEP1D(JIJ) * ZINVTSTEP)
  ENDDO
!$acc end kernels
!$acc kernels
  DO JK = IKTB , IKTE
!$acc loop independent
    DO JIJ = IIJB, IIJE
      ZMRCHANGE = ZMAX_TSTEP1D(JIJ) * POORHODZ(JIJ,JK)*(ZWSED(JIJ,JK+IKL)-ZWSED(JIJ,JK))
      PRXT(JIJ,JK) = PRXT(JIJ,JK) + ZMRCHANGE + PPRXS(JIJ,JK) * ZMAX_TSTEP1D(JIJ)
      PRXS(JIJ,JK) = PRXS(JIJ,JK) + ZMRCHANGE * ZINVTSTEP
      IF (GPRESENT_PFPR) THEN
        PFPR(JIJ,JK,KSPE) = PFPR(JIJ,JK,KSPE) + ZWSED(JIJ,JK) * (ZMAX_TSTEP1D(JIJ) * ZINVTSTEP)
      ENDIF
      IF (OELEC) THEN
        ZQCHANGE = ZMAX_TSTEP1D(JIJ) * POORHODZ(JIJ,JK) * (ZWSEDQ(JIJ,JK+IKL) - ZWSEDQ(JIJ,JK))
        PQXT(JIJ,JK) = PQXT(JIJ,JK) + ZQCHANGE + PPQXS(JIJ,JK) * ZMAX_TSTEP1D(JIJ)
        PQXS(JIJ,JK) = PQXS(JIJ,JK) + ZQCHANGE * ZINVTSTEP
      ENDIF
    ENDDO
  ENDDO
!$acc end kernels
  !
!$acc kernels
  ZANYREMAINT = .FALSE.
!$acc loop independent
  DO JIJ=IIJB,IIJE
    IF(ZREMAINT(JIJ)>0.) THEN
      ZANYREMAINT = .TRUE.
    END IF
  END DO
!$acc end kernels
END DO
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 1, ZHOOK_HANDLE)
END SUBROUTINE INTERNAL_SEDIM_SPLI
!
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT
END MODULE MODE_ICE4_SEDIMENTATION_SPLIT
