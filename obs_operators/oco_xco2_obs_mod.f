!-----------------------------------------------------------------------------
! OCO-2/OCO-3 Lite XCO2 observation operator for GEOS-Chem adjoint v09-02
!
! Supported input products:
!   OCO-2 Lite v11.2 (B11.2)
!   OCO-3 Lite v11   (B11)
!
! Required daily-file aliases:
!   ./data/OCO2/oco2_LtCO2_YYYYMMDD.nc4
!   ./data/OCO3/oco3_LtCO2_YYYYMMDD.nc4
!
! The aliases may be symbolic links to the original NASA filenames.  Both
! products use the same root variables read below.  Only quality flag zero is
! retained.  The operator acts on total CO2 tracer ID 1.
!-----------------------------------------------------------------------------

      MODULE OCO_XCO2_OBS_MOD

        IMPLICIT NONE
  
#     include "CMN_SIZE"
  
        PRIVATE
  
        INTEGER, PARAMETER :: MAXLEV_OCO = 20
        INTEGER, PARAMETER :: IDTCO2     = 1
        INTEGER, PARAMETER :: SAT_OCO2   = 2
        INTEGER, PARAMETER :: SAT_OCO3   = 3
  
        REAL*8, PARAMETER  :: TCVV_CO2 = 28.97d0 / 44.0d0
        REAL*8, PARAMETER  :: BAD_LIMIT = -9.0d5
        REAL*8, PARAMETER  :: TIME_EPS  = 1.0d-6
  
        CHARACTER(LEN=255), PARAMETER :: OCO2_TEMPLATE =
     &     './data/OCO2/oco2_LtCO2_YYYYMMDD.nc4'
        CHARACTER(LEN=255), PARAMETER :: OCO3_TEMPLATE =
     &     './data/OCO3/oco3_LtCO2_YYYYMMDD.nc4'
  
        TYPE OCO_XCO2_OBS
           INTEGER :: LOBS
           INTEGER :: SATELLITE
           INTEGER :: IGC
           INTEGER :: JGC
           INTEGER :: EDGE_FLAG
           INTEGER :: NSOUNDINGS
           REAL*8  :: LAT
           REAL*8  :: LON
           REAL*8  :: SEC_OF_DAY
           REAL*8  :: XCO2
           REAL*8  :: XCO2_ERR
           REAL*8  :: XCO2_APRIORI
           REAL*8  :: PRES(MAXLEV_OCO)
           REAL*8  :: PRES_WF(MAXLEV_OCO)
           REAL*8  :: AVG_KERNEL(MAXLEV_OCO)
           REAL*8  :: PRIOR_PROF(MAXLEV_OCO)
        END TYPE OCO_XCO2_OBS
  
        TYPE(OCO_XCO2_OBS), ALLOCATABLE, SAVE :: OCO(:)
  
        INTEGER, SAVE :: LOADED_DATE       = -1
        INTEGER, SAVE :: NOBS_LOADED       = 0
        INTEGER, SAVE :: OCO2_NOBS_LOADED  = 0
        INTEGER, SAVE :: OCO3_NOBS_LOADED  = 0
        INTEGER, SAVE :: OCO_NOBS_LAST     = 0
        INTEGER, SAVE :: OCO_NSOUND_LAST   = 0
        REAL*8,  SAVE :: OCO_COST_LAST     = 0d0
  
        ! These two values support later observation-error sensitivity tests.
        ! Effective sigma = MAX(OCO_ERR_FLOOR,
        !                       OCO_ERR_SCALE * xco2_uncertainty).
        REAL*8,  SAVE :: OCO_ERR_SCALE     = 1d0
        REAL*8,  SAVE :: OCO_ERR_FLOOR     = 0d0
        LOGICAL, SAVE :: USE_X2019_SCALE   = .FALSE.
  
        ! Observation-selection switches.  In superob mode the input file
        ! must contain superobs_count.  Edge filtering is recomputed from
        ! GEOS-Chem I/J indices, so it always follows the compiled grid.
        LOGICAL, SAVE :: USE_SUPEROBS_INPUT = .TRUE.
        LOGICAL, SAVE :: EXCLUDE_EDGE_OBS   = .TRUE.
        INTEGER, SAVE :: EDGE_WIDTH_CELLS   = 3
  
        ! Gridded diagnostics for the current observation window.
        REAL*8, SAVE :: CURR_GC_XCO2(IIPAR,JJPAR)  = 0d0
        REAL*8, SAVE :: CURR_OBS_XCO2(IIPAR,JJPAR) = 0d0
        REAL*8, SAVE :: CURR_DIFF(IIPAR,JJPAR)      = 0d0
        REAL*8, SAVE :: CURR_FORCING(IIPAR,JJPAR)   = 0d0
        REAL*8, SAVE :: CURR_COST(IIPAR,JJPAR)      = 0d0
        INTEGER, SAVE :: CURR_COUNT(IIPAR,JJPAR)    = 0
  
        ! Diagnostics accumulated over one optimizer evaluation.
        REAL*8, SAVE :: ALL_GC_XCO2(IIPAR,JJPAR)   = 0d0
        REAL*8, SAVE :: ALL_OBS_XCO2(IIPAR,JJPAR)  = 0d0
        REAL*8, SAVE :: ALL_DIFF(IIPAR,JJPAR)       = 0d0
        REAL*8, SAVE :: ALL_FORCING(IIPAR,JJPAR)    = 0d0
        REAL*8, SAVE :: ALL_COST(IIPAR,JJPAR)       = 0d0
        INTEGER, SAVE :: ALL_COUNT(IIPAR,JJPAR)     = 0
        INTEGER, SAVE :: DIAG_N_CALC                = -1

        ! Sampled operator values, indexed by retained OCO record.
        ! Reset on reread or optimizer evaluation; no duplicate rows.
        INTEGER, SAVE :: TRUTH_N_CALC = -9999
        LOGICAL, ALLOCATABLE, SAVE :: TRUTH_SEEN(:)
        INTEGER, ALLOCATABLE, SAVE :: TRUTH_GC_TIME(:)
        REAL*8, ALLOCATABLE, SAVE :: TRUTH_VALUE(:), TRUTH_SIGMA(:)
  
        PUBLIC :: CALC_OCO_XCO2_FORCE
        PUBLIC :: READ_OCO_XCO2_OBS
        PUBLIC :: SET_OCO_XCO2_OPTIONS
        PUBLIC :: OCO_COST_LAST, OCO_NOBS_LAST, OCO_NSOUND_LAST
        PUBLIC :: OCO2_NOBS_LOADED, OCO3_NOBS_LOADED
  
        CONTAINS
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE SET_OCO_XCO2_OPTIONS( ERR_SCALE, ERR_FLOOR,
     &                                 USE_X2019, SUPEROBS_MODE,
     &                                 DROP_EDGE, EDGE_CELLS )
  
        REAL*8,  INTENT(IN), OPTIONAL :: ERR_SCALE
        REAL*8,  INTENT(IN), OPTIONAL :: ERR_FLOOR
        LOGICAL, INTENT(IN), OPTIONAL :: USE_X2019
        LOGICAL, INTENT(IN), OPTIONAL :: SUPEROBS_MODE
        LOGICAL, INTENT(IN), OPTIONAL :: DROP_EDGE
        INTEGER, INTENT(IN), OPTIONAL :: EDGE_CELLS
  
        IF ( PRESENT(ERR_SCALE) ) THEN
           IF ( ERR_SCALE > 0d0 ) OCO_ERR_SCALE = ERR_SCALE
        ENDIF
        IF ( PRESENT(ERR_FLOOR) ) THEN
           IF ( ERR_FLOOR >= 0d0 ) OCO_ERR_FLOOR = ERR_FLOOR
        ENDIF
        IF ( PRESENT(USE_X2019) ) USE_X2019_SCALE = USE_X2019
        IF ( PRESENT(SUPEROBS_MODE) ) THEN
           USE_SUPEROBS_INPUT = SUPEROBS_MODE
        ENDIF
        IF ( PRESENT(DROP_EDGE) ) EXCLUDE_EDGE_OBS = DROP_EDGE
        IF ( PRESENT(EDGE_CELLS) ) THEN
           IF ( EDGE_CELLS >= 0 ) EDGE_WIDTH_CELLS = EDGE_CELLS
        ENDIF
  
        ! Force observations to be reread if options change between tests.
        LOADED_DATE = -1
  
        END SUBROUTINE SET_OCO_XCO2_OPTIONS
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE READ_OCO_XCO2_OBS( YYYYMMDD, NOBS )
  
        USE TIME_MOD, ONLY : EXPAND_DATE
  
        INTEGER, INTENT(IN)  :: YYYYMMDD
        INTEGER, INTENT(OUT) :: NOBS
  
        TYPE(OCO_XCO2_OBS), ALLOCATABLE :: OBS2(:), OBS3(:)
        CHARACTER(LEN=255)              :: FILE2, FILE3
        INTEGER                         :: N2, N3
  
        FILE2 = OCO2_TEMPLATE
        FILE3 = OCO3_TEMPLATE
        CALL EXPAND_DATE( FILE2, YYYYMMDD, 9999 )
        CALL EXPAND_DATE( FILE3, YYYYMMDD, 9999 )
  
        CALL READ_ONE_OCO_FILE( FILE2, YYYYMMDD, SAT_OCO2,
     &                        OBS2, N2 )
        CALL READ_ONE_OCO_FILE( FILE3, YYYYMMDD, SAT_OCO3,
     &                        OBS3, N3 )
  
        IF ( ALLOCATED(OCO) ) DEALLOCATE(OCO)
        IF ( ALLOCATED(TRUTH_SEEN) ) THEN
           DEALLOCATE(TRUTH_SEEN,TRUTH_GC_TIME,
     &                 TRUTH_VALUE,TRUTH_SIGMA)
        ENDIF
        NOBS = N2 + N3
        ALLOCATE( OCO(NOBS) )
        ALLOCATE(TRUTH_SEEN(NOBS),TRUTH_GC_TIME(NOBS),
     &            TRUTH_VALUE(NOBS),TRUTH_SIGMA(NOBS))
        TRUTH_SEEN = .FALSE.
        TRUTH_GC_TIME = 0
        TRUTH_VALUE = 0d0
        TRUTH_SIGMA = 0d0
  
        IF ( N2 > 0 ) OCO(1:N2)       = OBS2(1:N2)
        IF ( N3 > 0 ) OCO(N2+1:NOBS)  = OBS3(1:N3)
  
        IF ( ALLOCATED(OBS2) ) DEALLOCATE(OBS2)
        IF ( ALLOCATED(OBS3) ) DEALLOCATE(OBS3)
  
        LOADED_DATE      = YYYYMMDD
        NOBS_LOADED      = NOBS
        OCO2_NOBS_LOADED = N2
        OCO3_NOBS_LOADED = N3
  
        WRITE(6,*) ' OCO XCO2 retained: OCO-2=', N2,
     &           ' OCO-3=', N3, ' total=', NOBS
  
        END SUBROUTINE READ_OCO_XCO2_OBS
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE READ_ONE_OCO_FILE( FILENAME, YYYYMMDD, SATELLITE,
     &                              RECORDS, NKEEP )
  
        USE GRID_MOD,  ONLY : GET_IJ
        USE NETCDF
  
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER, INTENT(IN)          :: YYYYMMDD, SATELLITE
        TYPE(OCO_XCO2_OBS), ALLOCATABLE, INTENT(OUT) :: RECORDS(:)
        INTEGER, INTENT(OUT)         :: NKEEP
  
        INTEGER :: FID, DID_OBS, DID_LEV, NFILE, NLEV
        INTEGER :: VID_LAT, VID_LON, VID_DATE, VID_QF
        INTEGER :: VID_XCO2, VID_ERR, VID_XAP
        INTEGER :: VID_PRES, VID_PWF, VID_AK, VID_PRIOR
        INTEGER :: VID_NSUPER
        INTEGER :: STATUS, N, K, IIJJ(2), OBS_DATE
        LOGICAL :: EXISTS, VALID, HAS_NSUPER
  
        REAL*4, ALLOCATABLE :: TMP_LAT(:), TMP_LON(:)
        REAL*4, ALLOCATABLE :: TMP_XCO2(:), TMP_ERR(:), TMP_XAP(:)
        REAL*4, ALLOCATABLE :: TMP_PRES(:,:), TMP_PWF(:,:)
        REAL*4, ALLOCATABLE :: TMP_AK(:,:), TMP_PRIOR(:,:)
        INTEGER, ALLOCATABLE :: TMP_QF(:), TMP_DATE(:,:)
        INTEGER, ALLOCATABLE :: TMP_NSUPER(:)
  
        NKEEP = 0
        INQUIRE( FILE=TRIM(FILENAME), EXIST=EXISTS )
        IF ( .NOT. EXISTS ) THEN
           WRITE(6,*) ' OCO file not found; skipping: ', TRIM(FILENAME)
           RETURN
        ENDIF
  
        WRITE(6,*) ' Reading OCO Lite file: ', TRIM(FILENAME)
        CALL NC_CHECK( NF90_OPEN(TRIM(FILENAME),NF90_NOWRITE,FID),
     &                1000 + SATELLITE )
  
        CALL NC_CHECK( NF90_INQ_DIMID(FID,'sounding_id',DID_OBS),
     &               1010 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_DIMID(FID,'levels',DID_LEV),
     &               1020 + SATELLITE )
        CALL NC_CHECK( NF90_INQUIRE_DIMENSION(FID,DID_OBS,
     &               LEN=NFILE), 1030 + SATELLITE )
        CALL NC_CHECK( NF90_INQUIRE_DIMENSION(FID,DID_LEV,
     &               LEN=NLEV), 1040 + SATELLITE )
  
        IF ( NLEV /= MAXLEV_OCO ) THEN
           WRITE(6,*) ' Unexpected OCO level count: ', NLEV
           CALL NC_CHECK( NF90_CLOSE(FID), 1050 + SATELLITE )
           CALL OCO_ERROR_STOP( 'OCO Lite file must contain 20 levels')
        ENDIF
  
        CALL NC_CHECK( NF90_INQ_VARID(FID,'latitude',VID_LAT),
     &               1100 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'longitude',VID_LON),
     &               1110 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'date',VID_DATE),
     &               1120 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'xco2_quality_flag',VID_QF),
     &               1130 + SATELLITE )
  
        IF ( USE_X2019_SCALE ) THEN
           STATUS = NF90_INQ_VARID(FID,'xco2_x2019',VID_XCO2)
        ELSE
           STATUS = NF90_INQ_VARID(FID,'xco2',VID_XCO2)
        ENDIF
        CALL NC_CHECK( STATUS, 1140 + SATELLITE )
  
        CALL NC_CHECK( NF90_INQ_VARID(FID,'xco2_uncertainty',
     &               VID_ERR), 1150 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'xco2_apriori',VID_XAP),
     &               1160 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'pressure_levels',VID_PRES),
     &               1170 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'pressure_weight',VID_PWF),
     &               1180 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'xco2_averaging_kernel',
     &               VID_AK), 1190 + SATELLITE )
        CALL NC_CHECK( NF90_INQ_VARID(FID,'co2_profile_apriori',
     &               VID_PRIOR), 1200 + SATELLITE )
  
        STATUS = NF90_INQ_VARID(FID,'superobs_count',VID_NSUPER)
        HAS_NSUPER = STATUS == NF90_NOERR
        IF ( USE_SUPEROBS_INPUT .AND. .NOT. HAS_NSUPER ) THEN
           CALL NC_CHECK( NF90_CLOSE(FID), 1210 + SATELLITE )
           CALL OCO_ERROR_STOP(
     &      'Superob mode requires variable superobs_count' )
        ENDIF
  
        ALLOCATE( TMP_LAT(NFILE), TMP_LON(NFILE) )
        ALLOCATE( TMP_XCO2(NFILE), TMP_ERR(NFILE), TMP_XAP(NFILE) )
        ALLOCATE( TMP_QF(NFILE), TMP_DATE(7,NFILE) )
        ALLOCATE( TMP_NSUPER(NFILE) )
        ALLOCATE( TMP_PRES(MAXLEV_OCO,NFILE) )
        ALLOCATE( TMP_PWF(MAXLEV_OCO,NFILE) )
        ALLOCATE( TMP_AK(MAXLEV_OCO,NFILE) )
        ALLOCATE( TMP_PRIOR(MAXLEV_OCO,NFILE) )
  
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_LAT,TMP_LAT),
     &               1300 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_LON,TMP_LON),
     &               1310 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_DATE,TMP_DATE),
     &               1320 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_QF,TMP_QF),
     &               1330 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_XCO2,TMP_XCO2),
     &               1340 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_ERR,TMP_ERR),
     &               1350 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_XAP,TMP_XAP),
     &               1360 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_PRES,TMP_PRES),
     &               1370 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_PWF,TMP_PWF),
     &               1380 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_AK,TMP_AK),
     &               1390 + SATELLITE )
        CALL NC_CHECK( NF90_GET_VAR(FID,VID_PRIOR,TMP_PRIOR),
     &               1400 + SATELLITE )
        TMP_NSUPER(:) = 1
        IF ( HAS_NSUPER ) THEN
           CALL NC_CHECK( NF90_GET_VAR(FID,VID_NSUPER,TMP_NSUPER),
     &                  1405 + SATELLITE )
        ENDIF
        CALL NC_CHECK( NF90_CLOSE(FID), 1410 + SATELLITE )
  
        ! First pass: count quality-controlled records inside the model domain.
        DO N = 1, NFILE
           VALID = TMP_QF(N) == 0
           OBS_DATE = TMP_DATE(1,N) * 10000
     &            + TMP_DATE(2,N) * 100 + TMP_DATE(3,N)
           VALID = VALID .AND. OBS_DATE == YYYYMMDD
           VALID = VALID .AND. TMP_LAT(N) >= -90.0
           VALID = VALID .AND. TMP_LAT(N) <=  90.0
           VALID = VALID .AND. TMP_LON(N) >= -180.0
           VALID = VALID .AND. TMP_LON(N) <=  180.0
           VALID = VALID .AND. TMP_XCO2(N) > BAD_LIMIT
           VALID = VALID .AND. TMP_XAP(N)  > BAD_LIMIT
           VALID = VALID .AND. TMP_ERR(N)  > 0.0
           VALID = VALID .AND. MINVAL(TMP_PRES(:,N)) > 0.0
           VALID = VALID .AND. MINVAL(TMP_PWF(:,N)) > BAD_LIMIT
           VALID = VALID .AND. MINVAL(TMP_AK(:,N)) > BAD_LIMIT
           VALID = VALID .AND. MINVAL(TMP_PRIOR(:,N)) > BAD_LIMIT
           VALID = VALID .AND. TMP_NSUPER(N) > 0
           IF ( .NOT. VALID ) CYCLE
  
           IIJJ = GET_IJ( TMP_LON(N), TMP_LAT(N) )
           IF ( IIJJ(1) < 1 .OR. IIJJ(1) > IIPAR ) CYCLE
           IF ( IIJJ(2) < 1 .OR. IIJJ(2) > JJPAR ) CYCLE
           IF ( EXCLUDE_EDGE_OBS ) THEN
              IF ( IS_EDGE_CELL(IIJJ(1),IIJJ(2)) ) CYCLE
           ENDIF
           NKEEP = NKEEP + 1
        ENDDO
  
        ALLOCATE( RECORDS(NKEEP) )
        K = 0
  
        ! Second pass: copy retained observations into the compact array.
        DO N = 1, NFILE
           VALID = TMP_QF(N) == 0
           OBS_DATE = TMP_DATE(1,N) * 10000
     &            + TMP_DATE(2,N) * 100 + TMP_DATE(3,N)
           VALID = VALID .AND. OBS_DATE == YYYYMMDD
           VALID = VALID .AND. TMP_LAT(N) >= -90.0
           VALID = VALID .AND. TMP_LAT(N) <=  90.0
           VALID = VALID .AND. TMP_LON(N) >= -180.0
           VALID = VALID .AND. TMP_LON(N) <=  180.0
           VALID = VALID .AND. TMP_XCO2(N) > BAD_LIMIT
           VALID = VALID .AND. TMP_XAP(N)  > BAD_LIMIT
           VALID = VALID .AND. TMP_ERR(N)  > 0.0
           VALID = VALID .AND. MINVAL(TMP_PRES(:,N)) > 0.0
           VALID = VALID .AND. MINVAL(TMP_PWF(:,N)) > BAD_LIMIT
           VALID = VALID .AND. MINVAL(TMP_AK(:,N)) > BAD_LIMIT
           VALID = VALID .AND. MINVAL(TMP_PRIOR(:,N)) > BAD_LIMIT
           VALID = VALID .AND. TMP_NSUPER(N) > 0
           IF ( .NOT. VALID ) CYCLE
  
           IIJJ = GET_IJ( TMP_LON(N), TMP_LAT(N) )
           IF ( IIJJ(1) < 1 .OR. IIJJ(1) > IIPAR ) CYCLE
           IF ( IIJJ(2) < 1 .OR. IIJJ(2) > JJPAR ) CYCLE
           IF ( EXCLUDE_EDGE_OBS ) THEN
              IF ( IS_EDGE_CELL(IIJJ(1),IIJJ(2)) ) CYCLE
           ENDIF
  
           K = K + 1
           RECORDS(K)%LOBS          = N
           RECORDS(K)%SATELLITE     = SATELLITE
           RECORDS(K)%IGC           = IIJJ(1)
           RECORDS(K)%JGC           = IIJJ(2)
           RECORDS(K)%EDGE_FLAG     = 0
           IF ( IS_EDGE_CELL(IIJJ(1),IIJJ(2)) ) THEN
              RECORDS(K)%EDGE_FLAG  = 1
           ENDIF
           RECORDS(K)%NSOUNDINGS    = TMP_NSUPER(N)
           RECORDS(K)%LAT           = DBLE(TMP_LAT(N))
           RECORDS(K)%LON           = DBLE(TMP_LON(N))
           RECORDS(K)%SEC_OF_DAY    = DBLE(TMP_DATE(4,N))*3600d0
     &                            + DBLE(TMP_DATE(5,N))*60d0
     &                            + DBLE(TMP_DATE(6,N))
     &                            + DBLE(TMP_DATE(7,N))/1000d0
           RECORDS(K)%XCO2          = DBLE(TMP_XCO2(N))
           RECORDS(K)%XCO2_ERR      = DBLE(TMP_ERR(N))
           RECORDS(K)%XCO2_APRIORI  = DBLE(TMP_XAP(N))
           RECORDS(K)%PRES(:)       = DBLE(TMP_PRES(:,N))
           RECORDS(K)%PRES_WF(:)    = DBLE(TMP_PWF(:,N))
           RECORDS(K)%AVG_KERNEL(:) = DBLE(TMP_AK(:,N))
           RECORDS(K)%PRIOR_PROF(:) = DBLE(TMP_PRIOR(:,N))
        ENDDO
  
        DEALLOCATE( TMP_LAT, TMP_LON, TMP_XCO2, TMP_ERR, TMP_XAP )
        DEALLOCATE( TMP_QF, TMP_DATE, TMP_NSUPER )
        DEALLOCATE( TMP_PRES, TMP_PWF )
        DEALLOCATE( TMP_AK, TMP_PRIOR )
  
        END SUBROUTINE READ_ONE_OCO_FILE
  
  !-----------------------------------------------------------------------------
  
        LOGICAL FUNCTION IS_EDGE_CELL( I, J )
  
        INTEGER, INTENT(IN) :: I, J
  
        IS_EDGE_CELL = I <= EDGE_WIDTH_CELLS
     &          .OR. I > IIPAR - EDGE_WIDTH_CELLS
     &          .OR. J <= EDGE_WIDTH_CELLS
     &          .OR. J > JJPAR - EDGE_WIDTH_CELLS
  
        END FUNCTION IS_EDGE_CELL
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE CALC_OCO_XCO2_FORCE( COST_FUNC, COST_THIS_CALL )
  
        USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ, OBS_FREQ, N_CALC
        USE CHECKPT_MOD,    ONLY : CHK_STT
        USE DAO_MOD,        ONLY : AD
        USE PRESSURE_MOD,   ONLY : GET_PCENTER
        USE TIME_MOD,       ONLY : GET_NYMD, GET_NHMS
        USE TIME_MOD,       ONLY : YMD_EXTRACT, GET_TAUB, GET_TAU
  
        REAL*8, INTENT(INOUT)          :: COST_FUNC
        REAL*8, INTENT(OUT), OPTIONAL  :: COST_THIS_CALL
  
        INTEGER :: N, I, J, L, K, HH, MM, SS, NOBS, NUSED
        INTEGER :: NSOUND_USED
        REAL*8  :: T0, T1, SIGMA, DIFF, FORCE, JCOST
        REAL*8  :: XCO2_MODEL
        REAL*8  :: GC_PRES(LLPAR), GC_CO2(LLPAR)
        REAL*8  :: MAP(LLPAR,MAXLEV_OCO)
        REAL*8  :: OCO_CO2(MAXLEV_OCO)
        REAL*8  :: GC_CO2_ADJ(LLPAR)
  
        IF ( LOADED_DATE /= GET_NYMD() ) THEN
           CALL READ_OCO_XCO2_OBS( GET_NYMD(), NOBS )
        ENDIF
  
        CALL YMD_EXTRACT( GET_NHMS(), HH, MM, SS )
        IF ( TRUTH_N_CALC /= N_CALC ) THEN
           TRUTH_SEEN = .FALSE.
           TRUTH_N_CALC = N_CALC
        ENDIF
        T0 = DBLE(HH)*3600d0 + DBLE(MM)*60d0 + DBLE(SS)
        ! CALC_ADJ_FORCE_FOR_OBS is called every OBS_FREQ minutes.  Use
        ! that interval here; GET_TS_CHEM may be shorter than OBS_FREQ.
        T1 = T0 + OBS_FREQ*60d0
  
        JCOST = 0d0
        NUSED = 0
        NSOUND_USED = 0
  
        ! ALL_* must describe one cost-function evaluation, not a mixture
        ! of successive optimization iterations.
        IF ( DIAG_N_CALC /= N_CALC ) THEN
           ALL_GC_XCO2  = 0d0
           ALL_OBS_XCO2 = 0d0
           ALL_DIFF     = 0d0
           ALL_FORCING  = 0d0
           ALL_COST     = 0d0
           ALL_COUNT    = 0
           DIAG_N_CALC  = N_CALC
        ENDIF
  
        CURR_GC_XCO2  = 0d0
        CURR_OBS_XCO2 = 0d0
        CURR_DIFF     = 0d0
        CURR_FORCING  = 0d0
        CURR_COST     = 0d0
        CURR_COUNT    = 0
  
        DO N = 1, NOBS_LOADED
           IF ( OCO(N)%SEC_OF_DAY < T0 - TIME_EPS ) CYCLE
           IF ( OCO(N)%SEC_OF_DAY >= T1 - TIME_EPS ) CYCLE
           IF ( EXCLUDE_EDGE_OBS .AND. OCO(N)%EDGE_FLAG == 1 ) CYCLE
  
           I = OCO(N)%IGC
           J = OCO(N)%JGC
           IF ( MINVAL(AD(I,J,1:LLPAR)) <= 0d0 ) CYCLE
  
           DO L = 1, LLPAR
              GC_PRES(L) = GET_PCENTER(I,J,L)
              GC_CO2(L)  = CHK_STT(I,J,L,IDTCO2) * TCVV_CO2
     &                 * 1d6 / AD(I,J,L)
           ENDDO
  
           CALL BUILD_LOGP_MAP( GC_PRES, OCO(N)%PRES, MAP )
  
           OCO_CO2(:) = 0d0
           DO K = 1, MAXLEV_OCO
              DO L = 1, LLPAR
                 OCO_CO2(K) = OCO_CO2(K) + MAP(L,K)*GC_CO2(L)
              ENDDO
           ENDDO
  
           XCO2_MODEL = OCO(N)%XCO2_APRIORI
           DO K = 1, MAXLEV_OCO
              XCO2_MODEL = XCO2_MODEL
     &                 + OCO(N)%PRES_WF(K)*OCO(N)%AVG_KERNEL(K)
     &                 * ( OCO_CO2(K) - OCO(N)%PRIOR_PROF(K) )
           ENDDO
  
           SIGMA = MAX( OCO_ERR_FLOOR,
     &                OCO_ERR_SCALE*OCO(N)%XCO2_ERR )
           IF ( SIGMA <= 0d0 ) CYCLE

           TRUTH_SEEN(N) = .TRUE.
           TRUTH_GC_TIME(N) = GET_NHMS()
           TRUTH_VALUE(N) = XCO2_MODEL
           TRUTH_SIGMA(N) = SIGMA
  
           DIFF  = XCO2_MODEL - OCO(N)%XCO2
           FORCE = DIFF / (SIGMA*SIGMA)
           JCOST = JCOST + 0.5d0*DIFF*DIFF/(SIGMA*SIGMA)
           NUSED = NUSED + 1
           NSOUND_USED = NSOUND_USED + OCO(N)%NSOUNDINGS
  
           CURR_COUNT(I,J)    = CURR_COUNT(I,J) + 1
           CURR_GC_XCO2(I,J)  = CURR_GC_XCO2(I,J) + XCO2_MODEL
           CURR_OBS_XCO2(I,J) = CURR_OBS_XCO2(I,J) + OCO(N)%XCO2
           CURR_DIFF(I,J)     = CURR_DIFF(I,J) + DIFF
           CURR_FORCING(I,J)  = CURR_FORCING(I,J) + FORCE
           CURR_COST(I,J)     = CURR_COST(I,J)
     &                      + 0.5d0*DIFF*DIFF/(SIGMA*SIGMA)
  
           ALL_COUNT(I,J)     = ALL_COUNT(I,J) + 1
           ALL_GC_XCO2(I,J)   = ALL_GC_XCO2(I,J) + XCO2_MODEL
           ALL_OBS_XCO2(I,J)  = ALL_OBS_XCO2(I,J) + OCO(N)%XCO2
           ALL_DIFF(I,J)      = ALL_DIFF(I,J) + DIFF
           ALL_FORCING(I,J)   = ALL_FORCING(I,J) + FORCE
           ALL_COST(I,J)      = ALL_COST(I,J)
     &                      + 0.5d0*DIFF*DIFF/(SIGMA*SIGMA)
  
           GC_CO2_ADJ(:) = 0d0
           DO L = 1, LLPAR
              DO K = 1, MAXLEV_OCO
                 GC_CO2_ADJ(L) = GC_CO2_ADJ(L)
     &              + FORCE*OCO(N)%PRES_WF(K)
     &              * OCO(N)%AVG_KERNEL(K)*MAP(L,K)
              ENDDO
              STT_ADJ(I,J,L,IDTCO2) = STT_ADJ(I,J,L,IDTCO2)
     &           + GC_CO2_ADJ(L)*TCVV_CO2*1d6/AD(I,J,L)
           ENDDO
        ENDDO
  
        COST_FUNC     = COST_FUNC + JCOST
        OCO_COST_LAST = JCOST
        OCO_NOBS_LAST = NUSED
        OCO_NSOUND_LAST = NSOUND_USED
        IF ( PRESENT(COST_THIS_CALL) ) COST_THIS_CALL = JCOST
  
        WRITE(6,*) ' OCO XCO2 cost this step = ', JCOST,
     &           ' superobs = ', NUSED,
     &           ' source soundings = ', NSOUND_USED
  
        IF ( NUSED > 0 ) CALL MAKE_CURRENT_OCO_XCO2
        IF ( NUSED > 0 ) CALL WRITE_OCO_OSSE_TRUTH
        IF ( ABS(GET_TAUB()-GET_TAU()) < 1d-6 ) THEN
           CALL MAKE_AVERAGE_OCO_XCO2
        ENDIF
  
        END SUBROUTINE CALC_OCO_XCO2_FORCE
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE BUILD_LOGP_MAP( GC_PRES, RET_PRES, MAP )
  
        REAL*8, INTENT(IN)  :: GC_PRES(LLPAR)
        REAL*8, INTENT(IN)  :: RET_PRES(MAXLEV_OCO)
        REAL*8, INTENT(OUT) :: MAP(LLPAR,MAXLEV_OCO)
  
        INTEGER :: L, K, IMIN, IMAX
        REAL*8  :: PMIN, PMAX, DENOM, W2
  
        MAP  = 0d0
        IMIN = 1
        IMAX = 1
        PMIN = GC_PRES(1)
        PMAX = GC_PRES(1)
  
        DO L = 2, LLPAR
           IF ( GC_PRES(L) < PMIN ) THEN
              PMIN = GC_PRES(L)
              IMIN = L
           ENDIF
           IF ( GC_PRES(L) > PMAX ) THEN
              PMAX = GC_PRES(L)
              IMAX = L
           ENDIF
        ENDDO
  
        DO K = 1, MAXLEV_OCO
           IF ( RET_PRES(K) <= PMIN ) THEN
              MAP(IMIN,K) = 1d0
           ELSEIF ( RET_PRES(K) >= PMAX ) THEN
              MAP(IMAX,K) = 1d0
           ELSE
              DO L = 1, LLPAR-1
                 IF ( (RET_PRES(K)-GC_PRES(L))
     &              *(RET_PRES(K)-GC_PRES(L+1)) <= 0d0 ) THEN
                    DENOM = LOG(GC_PRES(L+1)) - LOG(GC_PRES(L))
                    IF ( ABS(DENOM) > 1d-14 ) THEN
                       W2 = ( LOG(RET_PRES(K))-LOG(GC_PRES(L)) )
     &                  / DENOM
                    ELSE
                       W2 = 0d0
                    ENDIF
                    MAP(L,K)   = 1d0 - W2
                    MAP(L+1,K) = W2
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
  
        END SUBROUTINE BUILD_LOGP_MAP
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE MAKE_CURRENT_OCO_XCO2
  
        CALL WRITE_OCO_XCO2_DIAG( .TRUE. )
  
        END SUBROUTINE MAKE_CURRENT_OCO_XCO2
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE MAKE_AVERAGE_OCO_XCO2
  
        CALL WRITE_OCO_XCO2_DIAG( .FALSE. )
  
        END SUBROUTINE MAKE_AVERAGE_OCO_XCO2
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE WRITE_OCO_XCO2_DIAG( CURRENT_WINDOW )
  
        USE ADJ_ARRAYS_MOD,    ONLY : N_CALC, EXPAND_NAME
        USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
        USE TIME_MOD,          ONLY : EXPAND_DATE, GET_NYMD, GET_NHMS
        USE NETCDF
  
        LOGICAL, INTENT(IN) :: CURRENT_WINDOW
  
        CHARACTER(LEN=255) :: FILENAME, OUTPUT_FILE
        INTEGER :: I, J, FID, LON_DIM_ID, LAT_DIM_ID
        INTEGER :: NPOINT
        REAL*4  :: TMP_DATA(IIPAR,JJPAR,6)
  
        IF ( CURRENT_WINDOW ) THEN
           FILENAME = 'gctm.oco_xco2.YYYYMMDD.hhmm.NN'
           CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
        ELSE
           FILENAME = 'gctm.oco_xco2.YYYYMMDD.NN'
           CALL EXPAND_DATE( FILENAME, GET_NYMD(), 9999 )
        ENDIF
        CALL EXPAND_NAME( FILENAME, N_CALC )
        OUTPUT_FILE = TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
  
        CALL NC_CHECK( NF90_CREATE(OUTPUT_FILE,NF90_CLOBBER,FID),
     &               8000 )
        CALL NC_CHECK( NF90_DEF_DIM(FID,'lon',IIPAR,LON_DIM_ID),
     &               8010 )
        CALL NC_CHECK( NF90_DEF_DIM(FID,'lat',JJPAR,LAT_DIM_ID),
     &               8020 )
        CALL NC_CHECK( NF90_ENDDEF(FID), 8030 )
  
        TMP_DATA = 0.0
        DO J = 1, JJPAR
           DO I = 1, IIPAR
              IF ( CURRENT_WINDOW ) THEN
                 NPOINT = CURR_COUNT(I,J)
                 IF ( NPOINT > 0 ) THEN
                    TMP_DATA(I,J,1) = REAL(CURR_GC_XCO2(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,2) = REAL(CURR_OBS_XCO2(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,3) = REAL(CURR_DIFF(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,4) = REAL(CURR_FORCING(I,J),4)
                    TMP_DATA(I,J,5) = REAL(CURR_COST(I,J),4)
                    TMP_DATA(I,J,6) = REAL(NPOINT,4)
                 ENDIF
              ELSE
                 NPOINT = ALL_COUNT(I,J)
                 IF ( NPOINT > 0 ) THEN
                    TMP_DATA(I,J,1) = REAL(ALL_GC_XCO2(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,2) = REAL(ALL_OBS_XCO2(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,3) = REAL(ALL_DIFF(I,J)
     &                                  /DBLE(NPOINT),4)
                    TMP_DATA(I,J,4) = REAL(ALL_FORCING(I,J),4)
                    TMP_DATA(I,J,5) = REAL(ALL_COST(I,J),4)
                    TMP_DATA(I,J,6) = REAL(NPOINT,4)
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
  
        CALL WRITE_NC_2D_FLOAT( FID, 'gc_xco2',
     &     'Mean GEOS-Chem XCO2 after OCO operator', 'ppm',
     &     TMP_DATA(:,:,1), LON_DIM_ID, LAT_DIM_ID )
        CALL WRITE_NC_2D_FLOAT( FID, 'obs_xco2',
     &     'Mean assimilated OCO XCO2', 'ppm', TMP_DATA(:,:,2),
     &     LON_DIM_ID, LAT_DIM_ID )
        CALL WRITE_NC_2D_FLOAT( FID, 'diff',
     &     'Mean model minus observation XCO2', 'ppm',
     &     TMP_DATA(:,:,3), LON_DIM_ID, LAT_DIM_ID )
        CALL WRITE_NC_2D_FLOAT( FID, 'forcing',
     &     'Sum of OCO scalar adjoint forcing', '1/ppm',
     &     TMP_DATA(:,:,4), LON_DIM_ID, LAT_DIM_ID )
        CALL WRITE_NC_2D_FLOAT( FID, 'cost',
     &     'Sum of OCO cost-function contribution', '1',
     &     TMP_DATA(:,:,5), LON_DIM_ID, LAT_DIM_ID )
        CALL WRITE_NC_2D_FLOAT( FID, 'count',
     &     'Number of assimilated OCO soundings', '1',
     &     TMP_DATA(:,:,6), LON_DIM_ID, LAT_DIM_ID )
  
        CALL NC_CHECK( NF90_CLOSE(FID), 8090 )
        WRITE(6,*) ' Wrote OCO XCO2 diagnostics: ', TRIM(OUTPUT_FILE)
  
        END SUBROUTINE WRITE_OCO_XCO2_DIAG
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE WRITE_NC_2D_FLOAT( FID, NAME, LONG_NAME, UNITS,
     &                              VAR, DIM1, DIM2 )
  
        USE NETCDF
  
        INTEGER, INTENT(IN)      :: FID, DIM1, DIM2
        CHARACTER(*), INTENT(IN) :: NAME, LONG_NAME, UNITS
        REAL*4, INTENT(IN)       :: VAR(IIPAR,JJPAR)
        INTEGER                  :: VID
  
        CALL NC_CHECK( NF90_REDEF(FID), 8100 )
        CALL NC_CHECK( NF90_DEF_VAR(FID,TRIM(NAME),NF90_FLOAT,
     &               (/DIM1,DIM2/),VID), 8110 )
        CALL NC_CHECK( NF90_PUT_ATT(FID,VID,'long_name',
     &               TRIM(LONG_NAME)), 8120 )
        CALL NC_CHECK( NF90_PUT_ATT(FID,VID,'units',
     &               TRIM(UNITS)), 8130 )
        CALL NC_CHECK( NF90_ENDDEF(FID), 8140 )
        CALL NC_CHECK( NF90_PUT_VAR(FID,VID,VAR), 8150 )
  
        END SUBROUTINE WRITE_NC_2D_FLOAT
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE NC_CHECK( STATUS, LOCATION )
  
        USE ERROR_MOD, ONLY : ERROR_STOP
        USE NETCDF
  
        INTEGER, INTENT(IN) :: STATUS, LOCATION
  
        IF ( STATUS /= NF90_NOERR ) THEN
           WRITE(6,*) TRIM(NF90_STRERROR(STATUS))
           WRITE(6,*) ' OCO netCDF error at location ', LOCATION
           CALL ERROR_STOP( 'netCDF error', 'oco_xco2_obs_mod' )
        ENDIF
  
        END SUBROUTINE NC_CHECK
  
  !-----------------------------------------------------------------------------
  
        SUBROUTINE OCO_ERROR_STOP( MESSAGE )
  
        USE ERROR_MOD, ONLY : ERROR_STOP
  
        CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
  
        CALL ERROR_STOP( MESSAGE, 'oco_xco2_obs_mod' )
  
        END SUBROUTINE OCO_ERROR_STOP
  
  !-----------------------------------------------------------------------------
  

      SUBROUTINE WRITE_OCO_OSSE_TRUTH
        USE NETCDF
        USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
        USE ADJ_ARRAYS_MOD, ONLY : EXPAND_NAME, OBS_FREQ
        USE TIME_MOD, ONLY : EXPAND_DATE

        INTEGER :: FID, DID, M, N, K, H
        INTEGER, ALLOCATABLE :: IDX(:)
        CHARACTER(LEN=255) :: FILENAME, FILE2, FILE3
        CHARACTER(LEN=3) :: SUFFIX

        M = COUNT(TRUTH_SEEN)
        IF ( M == 0 ) RETURN
        ALLOCATE(IDX(M))
        H = 0
        DO N = 1, NOBS_LOADED
            IF ( .NOT. TRUTH_SEEN(N) ) CYCLE
            H = H + 1
            IDX(H) = N
        ENDDO

        FILENAME = 'oco_osse_truth.YYYYMMDD.NN'
        CALL EXPAND_DATE(FILENAME,LOADED_DATE,9999)
        CALL EXPAND_NAME(FILENAME,TRUTH_N_CALC)
        FILENAME = TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
        FILE2 = OCO2_TEMPLATE
        FILE3 = OCO3_TEMPLATE
        CALL EXPAND_DATE(FILE2,LOADED_DATE,9999)
        CALL EXPAND_DATE(FILE3,LOADED_DATE,9999)

        CALL NC_CHECK(NF90_CREATE(FILENAME,NF90_CLOBBER,FID),9000)
        CALL NC_CHECK(NF90_DEF_DIM(FID,'nSamples',M,DID),9001)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'description',
     & 'Per-record H_OCO(x); truth only for a truth simulation'),9002)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'source_oco2',
     & TRIM(FILE2)),9003)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'source_oco3',
     & TRIM(FILE3)),9004)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'drop_edge',
     & MERGE(1,0,EXCLUDE_EDGE_OBS)),9005)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'edge_width',
     & EDGE_WIDTH_CELLS),9006)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'superobs_mode',
     & MERGE(1,0,USE_SUPEROBS_INPUT)),9007)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'use_x2019',
     & MERGE(1,0,USE_X2019_SCALE)),9008)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'obs_freq_minutes',
     & OBS_FREQ),9009)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'sampling',
     & 'Same time window, QC and edge selection as assimilation'),9010)
        CALL NC_CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'coverage',
     & 'Only records processed so far; not a completion marker'),9011)
        CALL NC_CHECK(NF90_ENDDEF(FID),9012)

        CALL OCO_TRUTH_INT(FID,DID,'observation_index',
     & 'Source file row, one based',
     & (/(OCO(IDX(N))%LOBS,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'satellite',
     & 'Satellite: 2=OCO-2, 3=OCO-3',
     & (/(OCO(IDX(N))%SATELLITE,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'date',
     & 'UTC date YYYYMMDD',
     & (/(LOADED_DATE,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'gc_time',
     & 'Model sampling time HHMMSS UTC',
     & (/(TRUTH_GC_TIME(IDX(N)),N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'model_i',
     & 'Model longitude index, one based',
     & (/(OCO(IDX(N))%IGC,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'model_j',
     & 'Model latitude index, one based',
     & (/(OCO(IDX(N))%JGC,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'edge_flag',
     & 'One for edge cell',
     & (/(OCO(IDX(N))%EDGE_FLAG,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'superobs_count',
     & 'Number of source soundings',
     & (/(OCO(IDX(N))%NSOUNDINGS,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'quality_flag',
     & 'Input quality flag, retained records are zero',
     & (/(0,N=1,M) /))

        CALL OCO_TRUTH_INT(FID,DID,'n_calc',
     & 'Optimizer evaluation',
     & (/(TRUTH_N_CALC,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'sec_of_day',
     & 'Observation UTC seconds since daily midnight','s',
     & (/(OCO(IDX(N))%SEC_OF_DAY,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'latitude',
     & 'Observation latitude','degrees_north',
     & (/(OCO(IDX(N))%LAT,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'longitude',
     & 'Observation longitude','degrees_east',
     & (/(OCO(IDX(N))%LON,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'xco2_true',
     & 'Unaveraged modeled XCO2 after OCO operator','ppm',
     & (/(TRUTH_VALUE(IDX(N)),N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'input_xco2',
     & 'Original input XCO2, not OSSE truth','ppm',
     & (/(OCO(IDX(N))%XCO2,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'xco2_uncertainty',
     & 'Original input uncertainty','ppm',
     & (/(OCO(IDX(N))%XCO2_ERR,N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'sigma_used',
     & 'Effective uncertainty after scale and floor','ppm',
     & (/(TRUTH_SIGMA(IDX(N)),N=1,M) /))

        CALL OCO_TRUTH_REAL(FID,DID,'xco2_apriori',
     & 'Retrieval a priori column','ppm',
     & (/(OCO(IDX(N))%XCO2_APRIORI,N=1,M) /))

       DO K = 1, MAXLEV_OCO
        WRITE(SUFFIX,'("_",I2.2)') K

        CALL OCO_TRUTH_REAL(FID,DID,
     &      'pressure_levels'//SUFFIX,
     &      'Retrieval pressure level','hPa',
     &      (/(OCO(IDX(N))%PRES(K),N=1,M)/))

        CALL OCO_TRUTH_REAL(FID,DID,
     &      'pressure_weight'//SUFFIX,
     &      'Retrieval pressure weighting','1',
     &      (/(OCO(IDX(N))%PRES_WF(K),N=1,M)/))

        CALL OCO_TRUTH_REAL(FID,DID,
     &      'xco2_averaging_kernel'//SUFFIX,
     &      'Retrieval averaging kernel','1',
     &      (/(OCO(IDX(N))%AVG_KERNEL(K),N=1,M)/))

        CALL OCO_TRUTH_REAL(FID,DID,
     &      'co2_profile_apriori'//SUFFIX,
     &      'Retrieval prior profile','ppm',
     &      (/(OCO(IDX(N))%PRIOR_PROF(K),N=1,M)/))

       ENDDO
        CALL NC_CHECK(NF90_CLOSE(FID),9013)
        DEALLOCATE(IDX)
        WRITE(6,*) ' Wrote OCO truth: ',TRIM(FILENAME),' rows=',M
        END SUBROUTINE WRITE_OCO_OSSE_TRUTH

        SUBROUTINE OCO_TRUTH_INT(FID,DID,NAME,LABEL,VAR)
        USE NETCDF
        INTEGER, INTENT(IN) :: FID,DID
        CHARACTER(*), INTENT(IN) :: NAME,LABEL
        INTEGER, INTENT(IN) :: VAR(:)
        INTEGER :: VID
        CALL NC_CHECK(NF90_REDEF(FID),9100)
        CALL NC_CHECK(NF90_DEF_VAR(FID,NAME,NF90_INT,
     & (/DID/),VID),9101)
        CALL NC_CHECK(NF90_PUT_ATT(FID,VID,'long_name',LABEL),9102)
        CALL NC_CHECK(NF90_ENDDEF(FID),9104)
        CALL NC_CHECK(NF90_PUT_VAR(FID,VID,VAR),9105)
      END SUBROUTINE OCO_TRUTH_INT

      SUBROUTINE OCO_TRUTH_REAL(FID,DID,NAME,LABEL,UNITS,VAR)
      USE NETCDF
        INTEGER, INTENT(IN) :: FID,DID
        CHARACTER(*), INTENT(IN) :: NAME,LABEL,UNITS
        REAL*8, INTENT(IN) :: VAR(:)
        INTEGER :: VID
        CALL NC_CHECK(NF90_REDEF(FID),9100)
        CALL NC_CHECK(NF90_DEF_VAR(FID,NAME,NF90_DOUBLE,
     & (/DID/),VID),9101)
        CALL NC_CHECK(NF90_PUT_ATT(FID,VID,'long_name',LABEL),9102)
        CALL NC_CHECK(NF90_PUT_ATT(FID,VID,'units',UNITS),9103)
        CALL NC_CHECK(NF90_ENDDEF(FID),9104)
        CALL NC_CHECK(NF90_PUT_VAR(FID,VID,VAR),9105)
        END SUBROUTINE OCO_TRUTH_REAL

      END MODULE OCO_XCO2_OBS_MOD
  
