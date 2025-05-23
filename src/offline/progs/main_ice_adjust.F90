PROGRAM MAIN_ICE_ADJUST

USE XRD_GETOPTIONS,  ONLY: INITOPTIONS, GETOPTION, CHECKOPTIONS
USE GETDATA_ICE_ADJUST_MOD, ONLY: GETDATA_ICE_ADJUST
USE COMPUTE_DIFF,    ONLY: DIFF
USE MODI_ICE_ADJUST
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PHYEX,      ONLY: PHYEX_t
USE STACK_MOD,       ONLY: STACK
USE MODE_MNH_ZWORK,  ONLY: ZMNH_STACK, IMNH_BLOCK, YMNH_STACK, INUMPIN
USE OMP_LIB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
#ifdef _OPENACC
USE MODD_UTIL_PHYEX_T, ONLY: COPY_PHYEX_T, WIPE_PHYEX_T
#endif

IMPLICIT NONE

INTEGER      :: KLEV
INTEGER      :: KRR

REAL,    ALLOCATABLE   :: PRHODJ         (:,:,:)   
REAL,    ALLOCATABLE   :: PEXNREF        (:,:,:)   
REAL,    ALLOCATABLE   :: PRHODREF       (:,:,:)   
REAL,    ALLOCATABLE   :: PPABSM         (:,:,:)   
REAL,    ALLOCATABLE   :: PTHT           (:,:,:)   
REAL,    ALLOCATABLE   :: PSIGS          (:,:,:)   
REAL,    ALLOCATABLE   :: PMFCONV        (:,:,:)   
REAL,    ALLOCATABLE   :: PRC_MF         (:,:,:)   
REAL,    ALLOCATABLE   :: PRI_MF         (:,:,:)   
REAL,    ALLOCATABLE   :: PCF_MF         (:,:,:)   
REAL,    ALLOCATABLE   :: PTHS           (:,:,:)   
REAL,    ALLOCATABLE   :: PRS            (:,:,:,:) 
REAL,    ALLOCATABLE   :: PSRCS          (:,:,:)   
REAL,    ALLOCATABLE   :: PCLDFR         (:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HRC       (:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HCF       (:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HRI       (:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HCF       (:,:,:)   
REAL,    ALLOCATABLE   :: ZRS            (:,:,:,:) 
REAL,    ALLOCATABLE   :: ZZZ            (:,:,:)   
REAL,    ALLOCATABLE   :: ZSIGQSAT       (:,:)   
REAL,    ALLOCATABLE   :: ZICE_CLD_WGT   (:,:)   
REAL,    ALLOCATABLE   :: PWEIGHT_MF_CLOUD(:,:,:)
REAL,    ALLOCATABLE   :: ZDUM1          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM2          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM3          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM4          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM5          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM6          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM7          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM8          (:,:,:)
REAL,    ALLOCATABLE   :: ZDUM9          (:,:,:)

REAL,    ALLOCATABLE   :: PRS_OUT        (:,:,:,:) 
REAL,    ALLOCATABLE   :: PSRCS_OUT      (:,:,:)   
REAL,    ALLOCATABLE   :: PCLDFR_OUT     (:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HRC_OUT   (:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HCF_OUT   (:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HRI_OUT   (:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HCF_OUT   (:,:,:)   

INTEGER :: NPROMA, NGPBLKS, NFLEVG
INTEGER :: JLON, JLEV
INTEGER, TARGET :: IBL

TYPE(DIMPHYEX_t)         :: D, D0
TYPE(PHYEX_t)            :: PHYEX
LOGICAL                  :: LLCHECK
LOGICAL                  :: LLCHECKDIFF
LOGICAL                  :: LLDIFF
INTEGER                  :: IBLOCK1, IBLOCK2
INTEGER                  :: ISTSZ(2), JBLK1, JBLK2
INTEGER                  :: NTID, ITID

REAL, ALLOCATABLE, TARGET :: PSTACK8(:,:)
REAL(KIND=4), ALLOCATABLE, TARGET :: PSTACK4(:,:)
TYPE(STACK), TARGET :: YLSTACK

REAL(KIND=8) :: TS,TE
REAL(KIND=8) :: TSC, TEC, TSD, TED, ZTC, ZTD 
INTEGER :: ITIME, NTIME
INTEGER :: IRANK, ISIZE
LOGICAL :: LLVERBOSE, LLSTAT, LLBIND
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

CALL INITOPTIONS ()
NGPBLKS = 296
CALL GETOPTION ("--blocks", NGPBLKS)
NPROMA = 32
CALL GETOPTION ("--nproma", NPROMA)
NFLEVG = -1
CALL GETOPTION ("--nflevg", NFLEVG)
CALL GETOPTION ("--check",  LLCHECK)
CALL GETOPTION ("--checkdiff",  LLCHECKDIFF)
IBLOCK1 = 1
CALL GETOPTION ("--check-block-1", IBLOCK1)
IBLOCK2 = NGPBLKS
CALL GETOPTION ("--check-block-2", IBLOCK2)
CALL GETOPTION ("--stat", LLSTAT)
NTIME = 1
CALL GETOPTION ("--times", NTIME)
CALL GETOPTION ("--verbose", LLVERBOSE)
CALL GETOPTION ("--bind", LLBIND)
CALL CHECKOPTIONS ()

LLDIFF = .FALSE.

IRANK = 0
ISIZE = 1
IF (LLBIND) THEN
  CALL LINUX_BIND      (IRANK, ISIZE)
  CALL LINUX_BIND_DUMP (IRANK, ISIZE)
ENDIF

CALL GETDATA_ICE_ADJUST (NPROMA, NGPBLKS, NFLEVG, PRHODJ, PEXNREF, PRHODREF, PPABSM, PTHT, ZICE_CLD_WGT,     &
& ZSIGQSAT, PSIGS, PMFCONV, PRC_MF, PRI_MF, PCF_MF, PWEIGHT_MF_CLOUD, ZDUM1, ZDUM2, ZDUM3, ZDUM4, ZDUM5, &
& ZDUM6, ZDUM7, ZDUM8, ZDUM9, PTHS, PRS, PSRCS, PCLDFR, PHLC_HRC, PHLC_HCF, &
& PHLI_HRI, PHLI_HCF, ZRS, ZZZ, PRS_OUT, PSRCS_OUT, PCLDFR_OUT, PHLC_HRC_OUT, PHLC_HCF_OUT,       &
& PHLI_HRI_OUT, PHLI_HCF_OUT, LLVERBOSE)

KLEV = SIZE (PRS, 2)
KRR  = SIZE (PRS, 3)

IF (LLVERBOSE) PRINT *, " KLEV = ", KLEV, " KRR = ", KRR

PRINT *, " NPROMA = ", NPROMA, " KLEV = ", KLEV, " NGPBLKS = ", NGPBLKS

CALL INIT_PHYEX(KRR, PHYEX)

D0%NIT  = NPROMA
D0%NIB  = 1
D0%NIE  = NPROMA
D0%NJT  = 1
D0%NJB  = 1
D0%NJE  = 1
D0%NIJT = D0%NIT * D0%NJT
D0%NIJB = 1
D0%NIJE = NPROMA
D0%NKL  = -1
D0%NKT  = KLEV
D0%NKA  = KLEV
D0%NKU  = 1
D0%NKB  = KLEV 
D0%NKE  = 1
D0%NKTB = 1
D0%NKTE = KLEV

ISTSZ=0
ISTSZ(KIND(KLEV)/4) = NPROMA * KLEV
ISTSZ(KIND(PRHODJ)/4) = ISTSZ(KIND(PRHODJ)/4) + NPROMA * 80 * KLEV
#ifndef USE_STACK
ISTSZ(2) = ISTSZ(2) + CEILING(ISTSZ(1) / 2.)
ISTSZ(1) = 0
#endif
ALLOCATE (PSTACK4 (ISTSZ(1), NGPBLKS))
ALLOCATE (PSTACK8 (ISTSZ(2), NGPBLKS))
ZMNH_STACK => PSTACK8

TS = OMP_GET_WTIME ()

ZTD = 0.
ZTC = 0.

IF (LHOOK) CALL DR_HOOK ('MAIN',0,ZHOOK_HANDLE)

DO ITIME = 1, NTIME

  TSD = OMP_GET_WTIME ()

#ifdef _OPENACC
  CALL COPY_PHYEX_T(PHYEX)
#endif

!$acc data &
!$acc      & copyin  (D0, &
!$acc      &          ZSIGQSAT, PRHODJ, PEXNREF, PRHODREF, PSIGS, PMFCONV, PPABSM, ZZZ, PCF_MF, PRC_MF, PRI_MF, &
!$acc      &          PWEIGHT_MF_CLOUD, ZDUM1, ZDUM2, ZDUM3, ZDUM4, ZDUM5, ZDUM6, ZDUM7, ZDUM8, ZDUM9, ZRS, ZICE_CLD_WGT) &
!$acc      & copy    (PRS, PTHS) &
!$acc      & copyout (PSRCS, PCLDFR, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF) &
!$acc      & create  (PSTACK4, PSTACK8) 

  TSC = OMP_GET_WTIME ()

#ifdef USE_OPENMP
!$OMP PARALLEL PRIVATE (D, YLSTACK, ITID, JBLK1, JBLK2)

  NTID = OMP_GET_MAX_THREADS ()
  ITID = OMP_GET_THREAD_NUM ()
  JBLK1 = 1 +  (NGPBLKS * (ITID+0)) / NTID
  JBLK2 =      (NGPBLKS * (ITID+1)) / NTID
#else
  JBLK1 = 1
  JBLK2 = NGPBLKS
#endif

  D = D0

!$acc parallel loop gang vector present (PHYEX) private (YLSTACK, IBL, JLON, D) collapse (2)

  DO IBL = JBLK1, JBLK2

#ifdef USE_COLCALL
  DO JLON = 1, NPROMA
    D = D0
    D%NIB = JLON
    D%NIE = JLON
    D%NIJB = JLON
    D%NIJE = JLON
#endif

#ifdef USE_STACK
    !Using cray pointers, AROME mechanism
    YLSTACK%L(1) = LOC (PSTACK4 (1, IBL))
    YLSTACK%U(1) = YLSTACK%L(1) + ISTSZ(1) * KIND (PSTACK4)
    YLSTACK%L(2) = LOC (PSTACK8 (1, IBL))
    YLSTACK%U(2) = YLSTACK%L(2) + ISTSZ(2) * KIND (PSTACK8)
#else
    !Using fortran indexing, Meso-NH mechanism
    YLSTACK%L(2) = 1
    YLSTACK%U(2) = ISTSZ(2)
    IMNH_BLOCK => IBL
    YMNH_STACK => YLSTACK
    INUMPIN = 0
#endif

    !Workaround for NEC: PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF must be provided to prevent
    !                    crash during execution, even if they are unused
    CALL ICE_ADJUST (D, PHYEX%CST, PHYEX%RAIN_ICE_PARAMN, PHYEX%NEBN, PHYEX%TURBN, PHYEX%PARAM_ICEN, &
    & PHYEX%MISC%TBUCONF, PHYEX%MISC%KRR, PHYEX%MISC%HBUNAME,     &
    & PHYEX%MISC%PTSTEP, ZSIGQSAT (:, IBL), PRHODJ=PRHODJ (:, :, IBL), &
    & PEXNREF=PEXNREF (:, :, IBL),                                                                                           &
    & PRHODREF=PRHODREF (:, :, IBL), PSIGS=PSIGS (:, :, IBL), LMFCONV=PHYEX%MISC%LMFCONV, PMFCONV=PMFCONV (:, :, IBL), &
    & PPABST=PPABSM (:, :, IBL), PZZ=ZZZ (:, :, IBL), PEXN=PEXNREF (:, :, IBL), PCF_MF=PCF_MF (:, :, IBL),          &
    & PRC_MF=PRC_MF (:, :, IBL), PRI_MF=PRI_MF  (:, :, IBL), PWEIGHT_MF_CLOUD=PWEIGHT_MF_CLOUD (:, :, IBL),               &
    & PICLDFR=ZDUM1(:, :, IBL), PWCLDFR=ZDUM2(:, :, IBL), PSSIO=ZDUM3(:, :, IBL),                                      &
    & PSSIU=ZDUM4(:, :, IBL), PIFR=ZDUM5(:, :, IBL),                                                                      &
    & PRV=ZRS(:, :, 1, IBL), PRC=ZRS(:, :, 2, IBL),                                                                       &
    & PRVS=PRS(:, :, 1, IBL), PRCS=PRS(:, :, 2, IBL), PTH=ZRS(:, :, 0, IBL), PTHS=PTHS (:, :, IBL),                 &
    & OCOMPUTE_SRC=PHYEX%MISC%OCOMPUTE_SRC,                                                                                     &
    & PSRCS=PSRCS (:, :, IBL), PCLDFR=PCLDFR (:, :, IBL), PRR=ZRS(:, :, 3, IBL), PRI=ZRS(:, :, 4, IBL),             &
    & PRIS=PRS(:, :, 4, IBL), PRS=ZRS(:, :, 5, IBL), PRG=ZRS(:, :, 6, IBL), &
    & TBUDGETS=PHYEX%MISC%YLBUDGET, KBUDGETS=PHYEX%MISC%NBUDGET,    &
    & PICE_CLD_WGT=ZICE_CLD_WGT(:, IBL),                                                                                     &
    & PHLC_HRC=PHLC_HRC(:, :, IBL), PHLC_HCF=PHLC_HCF(:, :, IBL),                                                         &
    & PHLI_HRI=PHLI_HRI(:, :, IBL), PHLI_HCF=PHLI_HCF(:, :, IBL),                                                         &
    & PHLC_HRC_MF=ZDUM6(:, :, IBL), PHLC_HCF_MF=ZDUM7(:, :, IBL), &
    & PHLI_HRI_MF=ZDUM8(:, :, IBL), PHLI_HCF_MF=ZDUM9(:, :, IBL) &
#ifdef USE_STACK
    & , YDSTACK=YLSTACK &
#endif
    & )

#ifdef USE_COLCALL
    ENDDO
#endif

  ENDDO

#ifdef USE_OPENMP
!$OMP END PARALLEL
#endif

!$acc end parallel loop

  TEC = OMP_GET_WTIME ()

!$acc end data

#ifdef _OPENACC
  CALL WIPE_PHYEX_T(PHYEX)
#endif

  TED = OMP_GET_WTIME ()

  ZTC = ZTC + (TEC - TSC)
  ZTD = ZTD + (TED - TSD)

ENDDO

IF (LHOOK) CALL DR_HOOK ('MAIN',1,ZHOOK_HANDLE)

TE = OMP_GET_WTIME()

WRITE (*,'(A,F8.2,A)') 'elapsed time : ',TE-TS,' s'
WRITE (*,'(A,F8.4,A)') '          i.e. ',1000.*(TE-TS)/(NPROMA*NGPBLKS)/NTIME,' ms/gp'

PRINT *, " ZTD = ", ZTD, ZTD / REAL (NPROMA*NGPBLKS*NTIME)
PRINT *, " ZTC = ", ZTC, ZTC / REAL (NPROMA*NGPBLKS*NTIME)

IF (LLCHECK .OR. LLSTAT .OR. LLCHECKDIFF) THEN
  DO IBL = IBLOCK1, IBLOCK2
    PRINT *, " IBL = ", IBL
    CALL DIFF ("PSRCS",    PSRCS_OUT    (:,:,IBL), PSRCS    (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PCLDFR",   PCLDFR_OUT   (:,:,IBL), PCLDFR   (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLC_HRC", PHLC_HRC_OUT (:,:,IBL), PHLC_HRC (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLC_HCF", PHLC_HCF_OUT (:,:,IBL), PHLC_HCF (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLI_HRI", PHLI_HRI_OUT (:,:,IBL), PHLI_HRI (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLI_HCF", PHLI_HCF_OUT (:,:,IBL), PHLI_HCF (:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
  ENDDO
ENDIF

IF (LLCHECKDIFF) THEN
  IF (LLDIFF) THEN
    PRINT*, "THERE ARE DIFF SOMEWHERE"
  ELSE
    PRINT*, "THERE IS NO DIFF AT ALL"
  ENDIF
ENDIF

STOP

CONTAINS

SUBROUTINE INIT_PHYEX(KRR, PHYEX)

USE MODD_BUDGET, ONLY: TBUCONF_ASSOCIATE, NBUDGET_RI, TBUCONF, LBU_ENABLE, LBUDGET_U, LBUDGET_V, LBUDGET_W, LBUDGET_TH, &
                       LBUDGET_TKE, LBUDGET_RV, LBUDGET_RC, LBUDGET_RR, LBUDGET_RI, LBUDGET_RS, LBUDGET_RG, LBUDGET_RH, LBUDGET_SV
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODI_INI_PHYEX, ONLY: INI_PHYEX

IMPLICIT NONE

! -----------------------------------------------------------------------
!     DUMMY VARIABLES
INTEGER,          INTENT(IN)  :: KRR
TYPE(PHYEX_t),    INTENT(OUT) :: PHYEX

!-----------------------------------------------------------------------
!    LOCAL VARIABLES
INTEGER :: IULOUT, JRR
REAL :: ZDZMIN
CHARACTER(LEN=6) :: CPROGRAM
CHARACTER(LEN=4) :: CMICRO, CSCONV, CTURB
REAL             :: PTSTEP
! -----------------------------------------------------------------------

IULOUT=20
CPROGRAM='AROME'
ZDZMIN=20.
CMICRO='ICE3'
CSCONV='NONE'
CTURB='TKEL'
PTSTEP = 50.000000000000000
PHYEX%MISC%TPFILE%NLU=0

!Default values
CALL INI_PHYEX(CPROGRAM, PHYEX%MISC%TPFILE, .TRUE., IULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, CSCONV, CTURB, &
              &LDDEFAULTVAL=.TRUE., LDREADNAM=.FALSE., LDCHECK=.FALSE., KPRINT=0, LDINIT=.FALSE., &
              &PHYEX_OUT=PHYEX)

!Control parameters
PHYEX%MISC%PTSTEP       = PTSTEP
PHYEX%MISC%HBUNAME      = 'DEPI'
PHYEX%MISC%LMFCONV      = .TRUE.
PHYEX%MISC%KRR          = KRR
PHYEX%MISC%OCOMPUTE_SRC = .TRUE.

!Emulate the namelist reading
PHYEX%PARAM_ICEN%LCRIAUTI=.TRUE.
PHYEX%PARAM_ICEN%XCRIAUTI_NAM=0.2E-3
PHYEX%PARAM_ICEN%XT0CRIAUTI_NAM=-5.
PHYEX%PARAM_ICEN%XCRIAUTC_NAM=0.1E-2
PHYEX%PARAM_ICEN%LOCND2=.FALSE.
PHYEX%PARAM_ICEN%CSEDIM='STAT'
PHYEX%PARAM_ICEN%LWARM=.TRUE.
PHYEX%PARAM_ICEN%LSEDIC=.TRUE.
PHYEX%PARAM_ICEN%CSNOWRIMING='M90 '
PHYEX%PARAM_ICEN%XFRACM90=0.1 ! Fraction used for the Murakami 1990 formulation
PHYEX%PARAM_ICEN%LCONVHG=.TRUE. ! TRUE to allow the conversion from hail to graupel
PHYEX%PARAM_ICEN%LCRFLIMIT=.TRUE. !True to limit rain contact freezing to possible heat exchange
PHYEX%PARAM_ICEN%LFEEDBACKT=.TRUE. ! When .TRUE. feed back on temperature is taken into account
PHYEX%PARAM_ICEN%LEVLIMIT=.TRUE.   ! When .TRUE. water vapour pressure is limited by saturation
PHYEX%PARAM_ICEN%LNULLWETG=.TRUE.  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
PHYEX%PARAM_ICEN%LWETGPOST=.TRUE.  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
PHYEX%PARAM_ICEN%LNULLWETH=.TRUE.  ! Same as LNULLWETG but for hail
PHYEX%PARAM_ICEN%LWETHPOST=.TRUE.  ! Same as LWETGPOST but for hail
PHYEX%PARAM_ICEN%LSEDIM_AFTER=.FALSE. ! Sedimentation done after microphysics
PHYEX%PARAM_ICEN%XSPLIT_MAXCFL=0.8
PHYEX%PARAM_ICEN%LDEPOSC=.FALSE.  ! water deposition on vegetation
PHYEX%PARAM_ICEN%XVDEPOSC=0.02    ! deposition speed (2 cm.s-1)
PHYEX%PARAM_ICEN%CSUBG_RC_RR_ACCR='NONE'
PHYEX%PARAM_ICEN%CSUBG_RR_EVAP='NONE'
PHYEX%PARAM_ICEN%CSUBG_PR_PDF='SIGM'
PHYEX%NEBN%LSUBG_COND   = .TRUE.
PHYEX%NEBN%LSIGMAS = .TRUE.
PHYEX%NEBN%CFRAC_ICE_ADJUST='S' ! Ice/liquid partition rule to use in adjustment
PHYEX%NEBN%CFRAC_ICE_SHALLOW_MF='S' ! Ice/liquid partition rule to use in shallow_mf

!Param initialisation
CALL INI_PHYEX(CPROGRAM, PHYEX%MISC%TPFILE, .TRUE., IULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, CSCONV, CTURB, &
              &LDDEFAULTVAL=.FALSE., LDREADNAM=.FALSE., LDCHECK=.TRUE., KPRINT=2, LDINIT=.TRUE., &
              &PHYEX_IN=PHYEX, PHYEX_OUT=PHYEX)

!Budgets
CALL TBUCONF_ASSOCIATE
PHYEX%MISC%NBUDGET=NBUDGET_RI
DO JRR=1, PHYEX%MISC%NBUDGET
  PHYEX%MISC%YLBUDGET(JRR)%NBUDGET=JRR
ENDDO
LBU_ENABLE=.FALSE.                                                                                                       
LBUDGET_U=.FALSE.
LBUDGET_V=.FALSE.
LBUDGET_W=.FALSE.
LBUDGET_TH=.FALSE.
LBUDGET_TKE=.FALSE.
LBUDGET_RV=.FALSE.
LBUDGET_RC=.FALSE.
LBUDGET_RR=.FALSE.
LBUDGET_RI=.FALSE.
LBUDGET_RS=.FALSE.
LBUDGET_RG=.FALSE.
LBUDGET_RH=.FALSE.
LBUDGET_SV=.FALSE.
PHYEX%MISC%TBUCONF=TBUCONF

END SUBROUTINE INIT_PHYEX

END PROGRAM

