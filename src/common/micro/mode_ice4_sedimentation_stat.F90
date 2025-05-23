!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION_STAT
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION_STAT(D, CST, ICEP, ICED, PARAMI, &
                                  &PTSTEP, KRR, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                                  &PRS, PRT, &
                                  &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN, &
                                  &PINPRH, PFPR)

!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!!      Ryad El Khatib 09-Oct-2019 Substantial re-write for optimization
!!       (outerunrolling, vectorization, memory cache saving, unrolling)
!  P. Wautelet 21/01/2021: initialize untouched part of PFPR
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
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
USE MODD_FIELDS_ADDRESS
USE MODI_GAMMA, ONLY: GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),             INTENT(IN)              :: D       !array dimensions
TYPE(CST_t),                  INTENT(IN)              :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)              :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)              :: ICED
TYPE(PARAM_ICE_t),            INTENT(IN)              :: PARAMI
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PT      ! Temperature
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT)       :: PRS     ! m.r. source
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)          :: PRT     ! m.r. at t
REAL, DIMENSION(D%NIJT),     INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(D%NIJT),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
LOGICAL :: LLSEA_AND_TOWN
INTEGER :: JRR, JIJ, JK, IKB, IKE,IKL, IIJB, IIJE, IKTB, IKTE
INTEGER :: ISHIFT, IK, IKPLUS
REAL :: ZINVTSTEP, ZGAC, ZGC, ZGAC2, ZGC2, ZRAYDEFO
REAL, DIMENSION(D%NIJT) :: ZTSORHODZ        ! TimeStep Over (Rhodref times delta Z)
REAL, DIMENSION(D%NIJT,0:1,2:KRR) :: ZSED   ! sedimentation flux array for each species and for above and current levels
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IIJB=D%NIJB
IIJE=D%NIJE
IKTB=D%NKTB
IKTE=D%NKTE
!
IF ( PRESENT( PFPR ) ) THEN
 !Set to 0. to avoid undefined values (in files)
 PFPR(:, : IKTB, :) = 0.
 PFPR(:, IKTE :, :) = 0.
END IF

!-------------------------------------------------------------------------------
!
!*       1.    compute the fluxes
!
ZINVTSTEP = 1./PTSTEP
ZGAC=GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)
ZGC=GAMMA(ICED%XNUC)
ZGAC2=GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)
ZGC2=GAMMA(ICED%XNUC2)
ZRAYDEFO=MAX(1.,0.5*(ZGAC/ZGC))
LLSEA_AND_TOWN=PRESENT(PSEA).AND.PRESENT(PTOWN)

!
!*       2.    compute the fluxes
!
! Start shift mechanism:
ISHIFT=0
CALL SHIFT

! Initialize vertical loop
DO JRR=IRC,KRR
  ZSED(:,IKPLUS,JRR) = 0.
ENDDO

! calculation sedimentation flux
DO JK = IKE , IKB, -1*IKL

  DO JIJ = IIJB, IIJE
    ZTSORHODZ(JIJ) =PTSTEP/(PRHODREF(JIJ,JK)*PDZZ(JIJ,JK))
  ENDDO
!
  DO JRR=IRC,KRR

    IF (JRR==IRC) THEN

      !******* for cloud
      IF (PARAMI%LSEDIC) THEN
        CALL CLOUD(PRT(:,JK,IRC))
      ELSE
        ZSED(:,IK,JRR)=0.
      ENDIF

    ELSEIF (JRR==IRR) THEN

      !*       2.2   for rain
      CALL OTHER_SPECIES(ICEP%XFSEDR,ICEP%XEXSEDR,PRT(:,JK,IRR))

    ELSEIF (JRR==IRI) THEN

      CALL PRISTINE_ICE(PRT(:,JK,IRI))

    ELSEIF (JRR==IRS) THEN

      !*       2.4   for aggregates/snow
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        CALL OTHER_SPECIES(ICEP%XFSEDS,ICEP%XEXSEDS,PRT(:,JK,IRS))
      ELSE
        CALL SNOW(PRT(:,JK,IRS))
      ENDIF

    ELSEIF (JRR==IRG) THEN

      !*       2.5   for graupeln
      CALL OTHER_SPECIES(ICEP%XFSEDG,ICEP%XEXSEDG,PRT(:,JK,IRG))

    ELSEIF (JRR==IRH) THEN

      !*       2.6   for hail
       CALL OTHER_SPECIES(ICEP%XFSEDH,ICEP%XEXSEDH,PRT(:,JK,IRH))

    ENDIF

  ENDDO ! JRR

  ! Wrap-up

  IF(PRESENT(PFPR)) THEN
    DO JRR=IRC,KRR
      PFPR(:,JK,JRR)=ZSED(:,IK,JRR)
    ENDDO
  ENDIF

  DO JRR=IRC, KRR
    DO JIJ = IIJB, IIJE
      PRS(JIJ,JK,JRR) = PRS(JIJ,JK,JRR)+ZTSORHODZ(JIJ)*(ZSED(JIJ,IKPLUS,JRR)-ZSED(JIJ,IK,JRR))*ZINVTSTEP
    ENDDO
  ENDDO

  IF (JK==IKB) THEN
    DO JIJ = IIJB, IIJE
      IF(PARAMI%LSEDIC) PINPRC(JIJ) = ZSED(JIJ,IK,2)/CST%XRHOLW
      PINPRR(JIJ) = ZSED(JIJ,IK,3)/CST%XRHOLW
      PINPRI(JIJ) = ZSED(JIJ,IK,4)/CST%XRHOLW
      PINPRS(JIJ) = ZSED(JIJ,IK,5)/CST%XRHOLW
      PINPRG(JIJ) = ZSED(JIJ,IK,6)/CST%XRHOLW
      IF (PRESENT(PINPRH) .AND. KRR==7) THEN
        PINPRH(JIJ) = ZSED(JIJ,IK,7)/CST%XRHOLW
      ENDIF
    ENDDO
  ENDIF

  ! shift mechanism : current level now takes the place of previous one
  CALL SHIFT

ENDDO ! JK

IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',1,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE CLOUD(PRXT)

    REAL, INTENT(IN)    :: PRXT(D%NIJT) ! mr of specy X

    REAL :: ZLBC    ! XLBC weighted by sea fraction
    REAL :: ZFSEDC
    REAL :: ZCONC3D ! droplet condensation
    REAL :: ZRAY    ! Cloud Mean radius
    REAL :: ZZWLBDA, ZZWLBDC, ZZCC
    INTEGER :: JIJ
    REAL :: ZQP
    REAL :: ZWSEDW1, ZWSEDW2

    !!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CLOUD',0,ZHOOK_HANDLE)

    DO JIJ = IIJB, IIJE
      !estimation of q' taking into account incoming ZWSED from previous vertical level
      ZQP=ZSED(JIJ,IKPLUS,JRR)*ZTSORHODZ(JIJ)
      IF ((PRXT(JIJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
        IF (LLSEA_AND_TOWN) THEN
          ZRAY   = MAX(1.,0.5*((1.-PSEA(JIJ))*ZGAC/ZGC+PSEA(JIJ)*ZGAC2/ZGC2))
          ZLBC   = MAX(MIN(ICED%XLBC(1),ICED%XLBC(2)),(PSEA(JIJ)*ICED%XLBC(2)+(1.-PSEA(JIJ))*ICED%XLBC(1)) )
          ZFSEDC = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)), (PSEA(JIJ)*ICEP%XFSEDC(2)+(1.-PSEA(JIJ))*ICEP%XFSEDC(1)) )
          ZCONC3D= (1.-PTOWN(JIJ))*(PSEA(JIJ)*ICED%XCONC_SEA+(1.-PSEA(JIJ))*ICED%XCONC_LAND) + &
                    PTOWN(JIJ)  *ICED%XCONC_URBAN
        ELSE
          ZRAY   = ZRAYDEFO
          ZLBC   = ICED%XLBC(1)
          ZFSEDC = ICEP%XFSEDC(1)
          ZCONC3D= ICED%XCONC_LAND
        ENDIF
        !calculation of w
        IF(PRXT(JIJ) > ICED%XRTMIN(JRR)) THEN
          ZZWLBDA=6.6E-8*(101325./PPABST(JIJ,JK))*(PTHT(JIJ,JK)/293.15)
          ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JIJ,JK)*PRXT(JIJ)))**ICED%XLBEXC
          ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
          ZWSEDW1=PRHODREF(JIJ,JK)**(-ICED%XCEXVT ) * ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC
        ELSE
          ZWSEDW1=0.
        ENDIF
        IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
          ZZWLBDA=6.6E-8*(101325./PPABST(JIJ,JK))*(PTHT(JIJ,JK)/293.15)
          ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JIJ,JK)*ZQP))**ICED%XLBEXC
          ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
          ZWSEDW2=PRHODREF(JIJ,JK)**(-ICED%XCEXVT ) * ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC
        ELSE
          ZWSEDW2=0.
        ENDIF
      ELSE
        ZWSEDW1=0.
        ZWSEDW2=0.
      ENDIF
!- duplicated code -------------------------------------------------------------------------
      IF (ZWSEDW2 /= 0.) THEN
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP) &
         & + FWSED2(ZWSEDW2,PTSTEP,PDZZ(JIJ,JK),ZSED(JIJ,IKPLUS,JRR))
      ELSE
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP)
      ENDIF
!-------------------------------------------------------------------------------------------
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CLOUD',1,ZHOOK_HANDLE)

  END SUBROUTINE CLOUD

  SUBROUTINE PRISTINE_ICE(PRXT)

    REAL, INTENT(IN)    :: PRXT(D%NIJT) ! mr of specy X
    INTEGER :: JIJ
    REAL :: ZQP
    REAL :: ZWSEDW1, ZWSEDW2

    !!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:PRISTINE_ICE',0,ZHOOK_HANDLE)

    ! ******* for pristine ice
    DO JIJ = IIJB, IIJE
      ZQP=ZSED(JIJ,IKPLUS,JRR)*ZTSORHODZ(JIJ)
      IF ((PRXT(JIJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
        !calculation of w
        IF ( PRXT(JIJ) > MAX(ICED%XRTMIN(JRR),1.0E-7 ) ) THEN
          ZWSEDW1= ICEP%XFSEDI *  &
                            & PRHODREF(JIJ,JK)**(-ICED%XCEXVT) * & !    McF&H
                            & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                            &      LOG(PRHODREF(JIJ,JK)*PRXT(JIJ)) )**ICEP%XEXCSEDI
        ELSE
          ZWSEDW1=0.
        ENDIF
        IF ( ZQP > MAX(ICED%XRTMIN(JRR),1.0E-7 ) ) THEN
          ZWSEDW2= ICEP%XFSEDI *  &
                            & PRHODREF(JIJ,JK)**(-ICED%XCEXVT) * & !    McF&H
                            & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                            &      LOG(PRHODREF(JIJ,JK)*ZQP) )**ICEP%XEXCSEDI
        ELSE
          ZWSEDW2=0.
        ENDIF
      ELSE
        ZWSEDW1=0.
        ZWSEDW2=0.
      ENDIF
!- duplicated code -------------------------------------------------------------------------
      IF (ZWSEDW2 /= 0.) THEN
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP) &
         & + FWSED2(ZWSEDW2,PTSTEP,PDZZ(JIJ,JK),ZSED(JIJ,IKPLUS,JRR))
      ELSE
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP)
      ENDIF
!-------------------------------------------------------------------------------------------
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:PRISTINE_ICE',1,ZHOOK_HANDLE)

  END SUBROUTINE PRISTINE_ICE

  SUBROUTINE SNOW(PRXT)

    REAL, INTENT(IN)    :: PRXT(D%NIJT) ! mr of specy X
    INTEGER :: JIJ
    REAL :: ZQP, ZLBDAS
    REAL :: ZWSEDW1, ZWSEDW2

    !!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:SNOW',0,ZHOOK_HANDLE)

    ! ******* for snow
    DO JIJ = IIJB, IIJE
      ZQP=ZSED(JIJ,IKPLUS,JRR)*ZTSORHODZ(JIJ)
      IF ((PRXT(JIJ) > ICED%XRTMIN(JRR)) ) THEN
        !Compute lambda_snow parameter
        IF (PARAMI%LSNOW_T) THEN 
          IF(PT(JIJ,JK)>CST%XTT-10.0) THEN
            ZLBDAS = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*PT(JIJ,JK))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
          ELSE
            ZLBDAS = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226-0.0106*PT(JIJ,JK))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
          END IF
        ELSE
          ZLBDAS  = MAX(MIN(ICED%XLBDAS_MAX,ICED%XLBS*(PRHODREF(JIJ,JK)*PRXT(JIJ))**ICED%XLBEXS),ICED%XLBDAS_MIN)
        END IF
        !calculation of w
        IF ( PRXT(JIJ) > ICED%XRTMIN(JRR) ) THEN
          ZWSEDW1= ICEP%XFSEDS *  &
                        & PRHODREF(JIJ,JK)**(-ICED%XCEXVT) * &
                        & (1+(ICED%XFVELOS/ZLBDAS)**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEXSEDS/ICED%XALPHAS)* &
                        & ZLBDAS**(ICED%XBS+ICEP%XEXSEDS) 
        ELSE
          ZWSEDW1=0.
        ENDIF
        IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
          ZWSEDW2= ICEP%XFSEDS *  &
                        & PRHODREF(JIJ,JK)**(-ICED%XCEXVT) * &
                        & (1+(ICED%XFVELOS/ZLBDAS)**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEXSEDS/ICED%XALPHAS)* &
                        & ZLBDAS**(ICED%XBS+ICEP%XEXSEDS) 
        ELSE
          ZWSEDW2=0.
        ENDIF
      ELSE
        ZWSEDW1=0.
        ZWSEDW2=0.
      ENDIF
!- duplicated code -------------------------------------------------------------------------
      IF (ZWSEDW2 /= 0.) THEN
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP) &
         & + FWSED2(ZWSEDW2,PTSTEP,PDZZ(JIJ,JK),ZSED(JIJ,IKPLUS,JRR))
      ELSE
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP)
      ENDIF
!-------------------------------------------------------------------------------------------
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:SNOW',1,ZHOOK_HANDLE)

  END SUBROUTINE SNOW

  SUBROUTINE OTHER_SPECIES(PFSED,PEXSED,PRXT)

    REAL, INTENT(IN)    :: PFSED
    REAL, INTENT(IN)    :: PEXSED
    REAL, INTENT(IN)    :: PRXT(D%NIJT) ! mr of specy X
    INTEGER :: JIJ
    REAL :: ZQP
    REAL :: ZWSEDW1, ZWSEDW2

    !!REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:OTHER_SPECIES',0,ZHOOK_HANDLE)

    ! for all but cloud and pristine ice :
    DO JIJ = IIJB, IIJE
      ZQP=ZSED(JIJ,IKPLUS,JRR)*ZTSORHODZ(JIJ)
      IF ((PRXT(JIJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
        !calculation of w
        IF ( PRXT(JIJ) > ICED%XRTMIN(JRR) ) THEN
          ZWSEDW1= PFSED *PRXT(JIJ)**(PEXSED-1)*PRHODREF(JIJ,JK)**(PEXSED-ICED%XCEXVT-1)
        ELSE
          ZWSEDW1=0.
        ENDIF
        IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
          ZWSEDW2= PFSED *ZQP**(PEXSED-1)*PRHODREF(JIJ,JK)**(PEXSED-ICED%XCEXVT-1)
        ELSE
          ZWSEDW2=0.
        ENDIF
      ELSE
        ZWSEDW1=0.
        ZWSEDW2=0.
      ENDIF
!- duplicated code -------------------------------------------------------------------------
      IF (ZWSEDW2 /= 0.) THEN
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP) &
         & + FWSED2(ZWSEDW2,PTSTEP,PDZZ(JIJ,JK),ZSED(JIJ,IKPLUS,JRR))
      ELSE
        ZSED(JIJ,IK,JRR)=FWSED1(ZWSEDW1,PTSTEP,PDZZ(JIJ,JK),PRHODREF(JIJ,JK),PRXT(JIJ),ZINVTSTEP)
      ENDIF
!-------------------------------------------------------------------------------------------
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:OTHER_SPECIES',1,ZHOOK_HANDLE)

  END SUBROUTINE OTHER_SPECIES

  SUBROUTINE SHIFT

    IKPLUS=ISHIFT
    IK=1-ISHIFT
    ISHIFT=1-ISHIFT

  END SUBROUTINE SHIFT
!
!
ELEMENTAL FUNCTION FWSED1(PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP) RESULT(PVAR)
  REAL, INTENT(IN) :: PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP
  REAL :: PVAR
! 5 multiplications only => cost = 5X
  PVAR = MIN(PRHODREF1*PDZZ1*PRXT1*PINVTSTEP,PWSEDW*PRHODREF1*PRXT1)
END FUNCTION FWSED1
!
ELEMENTAL FUNCTION FWSED2(PWSEDW,PTSTEP1,PDZZ1,PWSEDWSUP) RESULT(PVAR)
  REAL, INTENT(IN) :: PWSEDW,PTSTEP1,PDZZ1,PWSEDWSUP
  REAL :: PVAR
  PVAR = MAX(0.,1.-PDZZ1/(PTSTEP1*PWSEDW))*PWSEDWSUP
END FUNCTION FWSED2

END SUBROUTINE ICE4_SEDIMENTATION_STAT

END MODULE MODE_ICE4_SEDIMENTATION_STAT
