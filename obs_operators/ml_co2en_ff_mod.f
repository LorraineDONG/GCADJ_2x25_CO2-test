!$Id: xco2en_obs_mod.f,v 1.0 2026/05/20
      MODULE XCO2EN_OBS_MOD

       IMPLICIT NONE

       !=================================================================
       ! MODULE VARIABLES
       !=================================================================

#     include "CMN_SIZE"

       ! Record to store data from each ML XCO2en obs
       TYPE ML_XCO2EN_OBS
              INTEGER     :: LOBS      
              REAL*8      :: LAT      
              REAL*8      :: LON      
              REAL*8      :: TIME     
              INTEGER     :: QF        
              REAL*8      :: XCO2EN      ! 机器学习重构的增量
              REAL*8      :: ERR         ! 观测不确定性
       END TYPE ML_XCO2EN_OBS

       TYPE(ML_XCO2EN_OBS), ALLOCATABLE :: OBS(:)

       ! 总 CO2 的示踪物 ID
       INTEGER, PARAMETER :: IDTCO2 = 1
       ! 定义化石燃料 CO2 的示踪物 ID
       INTEGER, PARAMETER :: IDTCO2_FF = 2
       
       ! 空气与 CO2 的分子量比
       REAL*8,  PARAMETER :: TCVV_CO2 = 28.97d0  / 44.0d0
       
       !-----------------------------------------------------------------
       ! Arrays for Diagnostics (Current Window & All-Time)
       !-----------------------------------------------------------------
       REAL*8 :: CURR_GC_XCO2EN(IIPAR, JJPAR)  = 0d0
       REAL*8 :: CURR_OBS_XCO2EN(IIPAR, JJPAR) = 0d0
       REAL*8 :: CURR_DIFF(IIPAR, JJPAR)       = 0d0
       REAL*8 :: CURR_FORCING(IIPAR, JJPAR)    = 0d0
       REAL*8 :: CURR_COST(IIPAR, JJPAR)       = 0d0
       INTEGER:: CURR_COUNT(IIPAR, JJPAR)      = 0
       
       REAL*8 :: ALL_GC_XCO2EN(IIPAR, JJPAR)   = 0d0
       REAL*8 :: ALL_OBS_XCO2EN(IIPAR, JJPAR)  = 0d0
       REAL*8 :: ALL_DIFF(IIPAR, JJPAR)        = 0d0
       REAL*8 :: ALL_FORCING(IIPAR, JJPAR)     = 0d0
       REAL*8 :: ALL_COST(IIPAR, JJPAR)        = 0d0
       INTEGER:: ALL_COUNT(IIPAR, JJPAR)       = 0

       ! Keep diagnostics and truth records separate for each N_CALC.
       INTEGER, SAVE :: DIAG_N_CALC = -9999
       INTEGER, SAVE :: CURRENT_NOBS = 0

       ! Compact per-observation truth buffer for the current input day.
       INTEGER, SAVE :: TRUTH_DATE    = -1
       INTEGER, SAVE :: TRUTH_N_CALC  = -1
       INTEGER, SAVE :: TRUTH_NMATCH  = 0
       CHARACTER(LEN=255), SAVE :: CURRENT_INPUT_FILE = ''

       INTEGER, ALLOCATABLE :: TRUTH_LOBS(:)
       INTEGER, ALLOCATABLE :: TRUTH_YYYYMMDD(:)
       INTEGER, ALLOCATABLE :: TRUTH_GC_TIME(:)
       INTEGER, ALLOCATABLE :: TRUTH_I(:), TRUTH_J(:)
       INTEGER, ALLOCATABLE :: TRUTH_QF(:), TRUTH_CALC(:)
       REAL*8,  ALLOCATABLE :: TRUTH_OBS_TIME(:)
       REAL*8,  ALLOCATABLE :: TRUTH_LAT(:), TRUTH_LON(:)
       REAL*8,  ALLOCATABLE :: TRUTH_XCO2EN(:)
       REAL*8,  ALLOCATABLE :: TRUTH_INPUT_XCO2EN(:)
       REAL*8,  ALLOCATABLE :: TRUTH_UNCERTAINTY(:)
       REAL*8,  ALLOCATABLE :: TRUTH_AIR_MASS(:)
       REAL*8,  ALLOCATABLE :: TRUTH_FF_MASS(:)

       CONTAINS
!-----------------------------------------------------------------------

       SUBROUTINE READ_XCO2EN_OBS( YYYYMMDD, NOBS )
              ! 读取机器学习重构的 XCO2en 数据
              USE DIRECTORY_MOD,          ONLY : DATA_DIR
              USE NETCDF
              USE TIME_MOD,               ONLY : EXPAND_DATE
              USE ERROR_MOD,              ONLY : ERROR_STOP
       
              IMPLICIT NONE 
       
              INTEGER, INTENT(IN)  :: YYYYMMDD
              INTEGER, INTENT(OUT) :: NOBS
              INTEGER, PARAMETER   :: DP = KIND(1.0D0)

              ! --- NetCDF ID Variables ---
              INTEGER :: FID, NT_ID, QF_ID, ERR_ID, TM_ID
              INTEGER :: XCO2EN_ID, LA_ID, LO_ID
       
              ! --- Temporary Arrays ---
              REAL(DP), ALLOCATABLE :: TMP_LAT(:), TMP_LON(:)
              REAL(DP), ALLOCATABLE :: TMP_TIME(:), TMP_XCO2EN(:)
              REAL(DP), ALLOCATABLE :: TMP_ERR(:)      
              INTEGER,  ALLOCATABLE :: TMP_QF(:)
       
              INTEGER :: N, STATUS
              LOGICAL :: FILE_EXISTS
              CHARACTER(LEN=255) :: READ_FILENAME


       !=================================================================
       ! READ_GOS_CO2_OBS begins here!
       !=================================================================

       ! 【dwh】根据实际文件名修改前缀
       READ_FILENAME = TRIM( 'ml-xco2en-YYYYMMDD.nc' )
       CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )
       READ_FILENAME = TRIM( './data/ML_XCO2EN/' ) //
     &                      TRIM( READ_FILENAME )

       ! CHANGE 1: Skip a missing daily observation file.
       ! A day without an observation file contributes zero to the
       ! cost function and adjoint forcing.  Discard yesterday's data.
       NOBS = 0
       IF ( ALLOCATED(OBS) ) DEALLOCATE( OBS )
       CURRENT_INPUT_FILE = ''

       INQUIRE( FILE=TRIM(READ_FILENAME), EXIST=FILE_EXISTS )
       IF ( .NOT. FILE_EXISTS ) THEN
              WRITE(6,*) '    - READ_XCO2EN_OBS: no file for ',
     &                   YYYYMMDD, '; skipping ', TRIM(READ_FILENAME)
              RETURN
       ENDIF

       WRITE(6,*) '    - READ_XCO2EN_OBS: reading: ', READ_FILENAME
       STATUS = NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID )
       IF ( STATUS /= NF90_NOERR ) THEN
              ! Cover a file removed between INQUIRE and NF90_OPEN.
              INQUIRE( FILE=TRIM(READ_FILENAME), EXIST=FILE_EXISTS )
              IF ( .NOT. FILE_EXISTS ) THEN
                     WRITE(6,*) '    - READ_XCO2EN_OBS: vanished; ',
     &                          'skipping ', TRIM(READ_FILENAME)
                     RETURN
              ENDIF
              CALL CHECK( STATUS, 0 )
       ENDIF
       CURRENT_INPUT_FILE = TRIM( READ_FILENAME )

       ! 【dwh】获取变量 ID (请确保这些字符串与你生成的nc文件内部变量名一致)
       CALL CHECK( NF90_INQ_DIMID( FID, "nSamples", NT_ID),  101)
       CALL CHECK( NF90_INQ_VARID( FID, "xco2_enhancement", 
     &                                  XCO2EN_ID ), 102)
       CALL CHECK( NF90_INQ_VARID( FID, "latitude",  LA_ID), 103)
       CALL CHECK( NF90_INQ_VARID( FID, "longitude", LO_ID), 104)
       CALL CHECK( NF90_INQ_VARID( FID, "time",      TM_ID), 105)
       CALL CHECK( NF90_INQ_VARID( FID, "qa_flag",QF_ID),106)
       CALL CHECK( NF90_INQ_VARID( FID, "uncertainty", ERR_ID),107)

       CALL CHECK(NF90_INQUIRE_DIMENSION(FID,NT_ID,LEN=NOBS), 201)

       WRITE(6,*) ' Number of observations (NOBS): ', NOBS
       IF (ALLOCATED(OBS)) DEALLOCATE(OBS)
       ALLOCATE( OBS(NOBS) )

       ALLOCATE( TMP_LAT(NOBS), TMP_LON(NOBS), TMP_TIME(NOBS) )
       ALLOCATE( TMP_XCO2EN(NOBS), TMP_ERR(NOBS), TMP_QF(NOBS) )

       CALL CHECK( NF90_GET_VAR( FID, LA_ID,     TMP_LAT   ), 301 )
       CALL CHECK( NF90_GET_VAR( FID, LO_ID,     TMP_LON   ), 302 )
       CALL CHECK( NF90_GET_VAR( FID, TM_ID,     TMP_TIME  ), 303 )
       CALL CHECK( NF90_GET_VAR( FID, XCO2EN_ID, TMP_XCO2EN), 304 )
       CALL CHECK( NF90_GET_VAR( FID, QF_ID,     TMP_QF    ), 305 )
       CALL CHECK( NF90_GET_VAR( FID, ERR_ID,    TMP_ERR   ), 306 )

       CALL CHECK( NF90_CLOSE( FID ), 999 )

       DO N = 1, NOBS
              OBS(N)%LOBS   = N
              OBS(N)%LAT    = TMP_LAT(N)
              OBS(N)%LON    = TMP_LON(N)
              OBS(N)%TIME   = TMP_TIME(N)
              OBS(N)%QF     = TMP_QF(N)
              OBS(N)%XCO2EN = TMP_XCO2EN(N)
              OBS(N)%ERR    = TMP_ERR(N)
       END DO

       DEALLOCATE( TMP_LAT, TMP_LON, TMP_TIME )
       DEALLOCATE( TMP_XCO2EN, TMP_ERR, TMP_QF )

       END SUBROUTINE READ_XCO2EN_OBS
!-----------------------------------------------------------------------

       SUBROUTINE CHECK( STATUS, LOCATION )
       ! Reference to f90 modules
       USE ERROR_MOD,    ONLY  : ERROR_STOP
       USE NETCDF

       ! Arguments
       INTEGER, INTENT(IN)    :: STATUS
       INTEGER, INTENT(IN)    :: LOCATION

       !=================================================================
       ! CHECK begins here!
       !=================================================================

       IF ( STATUS /= NF90_NOERR ) THEN
              WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
              WRITE(6,*) 'At location = ', LOCATION
              CALL ERROR_STOP('netCDF error', 'xco2en_obs_mod')
       ENDIF

       ! Return to calling program
       END SUBROUTINE CHECK

!-----------------------------------------------------------------------
       SUBROUTINE CLEAR_OSSE_TRUTH

       IF ( ALLOCATED(TRUTH_LOBS) ) THEN
              DEALLOCATE( TRUTH_LOBS, TRUTH_YYYYMMDD, TRUTH_GC_TIME,
     &                  TRUTH_I, TRUTH_J, TRUTH_QF, TRUTH_CALC,
     &                  TRUTH_OBS_TIME, TRUTH_LAT, TRUTH_LON,
     &                  TRUTH_XCO2EN, TRUTH_INPUT_XCO2EN,
     &                  TRUTH_UNCERTAINTY, TRUTH_AIR_MASS,
     &                  TRUTH_FF_MASS )
       ENDIF

       TRUTH_DATE   = -1
       TRUTH_N_CALC = -1
       TRUTH_NMATCH = 0

       END SUBROUTINE CLEAR_OSSE_TRUTH

!-----------------------------------------------------------------------
       SUBROUTINE INIT_OSSE_TRUTH( NMAX, YYYYMMDD, ICALC )

       INTEGER, INTENT(IN) :: NMAX, YYYYMMDD, ICALC

       CALL CLEAR_OSSE_TRUTH
       IF ( NMAX <= 0 ) RETURN

       ALLOCATE( TRUTH_LOBS(NMAX), TRUTH_YYYYMMDD(NMAX),
     &          TRUTH_GC_TIME(NMAX), TRUTH_I(NMAX), TRUTH_J(NMAX),
     &          TRUTH_QF(NMAX), TRUTH_CALC(NMAX),
     &          TRUTH_OBS_TIME(NMAX), TRUTH_LAT(NMAX),
     &          TRUTH_LON(NMAX), TRUTH_XCO2EN(NMAX),
     &          TRUTH_INPUT_XCO2EN(NMAX),
     &          TRUTH_UNCERTAINTY(NMAX), TRUTH_AIR_MASS(NMAX),
     &          TRUTH_FF_MASS(NMAX) )

       TRUTH_DATE   = YYYYMMDD
       TRUTH_N_CALC = ICALC
       TRUTH_NMATCH = 0

       END SUBROUTINE INIT_OSSE_TRUTH

!-----------------------------------------------------------------------
       SUBROUTINE STORE_OSSE_TRUTH( NT, I, J, GC_TIME,
     &                             SIM_XCO2EN, AIR_MASS, FF_MASS )
       USE ERROR_MOD, ONLY : ERROR_STOP

       INTEGER, INTENT(IN) :: NT, I, J, GC_TIME
       REAL*8,  INTENT(IN) :: SIM_XCO2EN, AIR_MASS, FF_MASS
       INTEGER             :: K

       IF ( .NOT. ALLOCATED(TRUTH_LOBS) ) RETURN

       K = TRUTH_NMATCH + 1
       IF ( K > SIZE(TRUTH_LOBS) ) THEN
              CALL ERROR_STOP( 'OSSE truth buffer is too small',
     &                       'STORE_OSSE_TRUTH' )
       ENDIF

       TRUTH_NMATCH          = K
       TRUTH_LOBS(K)         = OBS(NT)%LOBS
       TRUTH_YYYYMMDD(K)     = TRUTH_DATE
       TRUTH_OBS_TIME(K)     = OBS(NT)%TIME
       TRUTH_GC_TIME(K)      = GC_TIME
       TRUTH_LAT(K)          = OBS(NT)%LAT
       TRUTH_LON(K)          = OBS(NT)%LON
       TRUTH_I(K)            = I
       TRUTH_J(K)            = J
       TRUTH_QF(K)           = OBS(NT)%QF
       TRUTH_CALC(K)         = TRUTH_N_CALC
       TRUTH_XCO2EN(K)       = SIM_XCO2EN
       TRUTH_INPUT_XCO2EN(K) = OBS(NT)%XCO2EN
       TRUTH_UNCERTAINTY(K)  = OBS(NT)%ERR
       TRUTH_AIR_MASS(K)     = AIR_MASS
       TRUTH_FF_MASS(K)      = FF_MASS

       END SUBROUTINE STORE_OSSE_TRUTH

!-----------------------------------------------------------------------
       SUBROUTINE CALC_XCO2EN_FORCE( COST_FUNC )
              USE ADJ_ARRAYS_MOD,   ONLY : STT_ADJ, N_CALC
              USE CHECKPT_MOD,      ONLY : CHK_STT
              USE DAO_MOD,          ONLY : AD
              USE GRID_MOD,         ONLY : GET_IJ
              USE TIME_MOD,         ONLY : GET_NYMD, GET_NHMS
              USE TIME_MOD,         ONLY : GET_TS_CHEM,GET_TAUB,GET_TAU
       
              REAL*8, INTENT(INOUT)       :: COST_FUNC
              INTEGER                     :: NTSTART, NTSTOP, NT
              INTEGER                     :: IIJJ(2), I, J, L, NOBS
              REAL*8                      :: TOTAL_AIR_MASS
              REAL*8                      :: TOTAL_FF_MASS
              REAL*8                      :: SIM_XCO2EN
              REAL*8                      :: DIFF_SCALAR
              REAL*8                      :: OBS_ERR_SQ
              REAL*8                      :: FORCE_SCALAR
              REAL*8                      :: TMP_COST
              REAL*8                      :: ADJ_GRADIENT
              REAL*8                      :: OLD_COST
              LOGICAL, SAVE               :: FIRST = .TRUE.
       
              print*, '     - CALC_XCO2EN_FORCE '
              OLD_COST = COST_FUNC

              ! Never mix diagnostics from different optimizer calls.
              IF ( DIAG_N_CALC /= N_CALC ) THEN
                     ALL_GC_XCO2EN  = 0d0
                     ALL_OBS_XCO2EN = 0d0
                     ALL_DIFF       = 0d0
                     ALL_FORCING    = 0d0
                     ALL_COST       = 0d0
                     ALL_COUNT      = 0
                     CURRENT_NOBS   = 0
                     CALL CLEAR_OSSE_TRUTH
                     DIAG_N_CALC    = N_CALC
              ENDIF
       
              IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN
                     CALL READ_XCO2EN_OBS( GET_NYMD(), CURRENT_NOBS )
                     CALL INIT_OSSE_TRUTH( CURRENT_NOBS, GET_NYMD(),
     &                                  N_CALC )
              ENDIF

              NOBS = CURRENT_NOBS
       
              CALL GET_NT_RANGE( NOBS, GET_NHMS(), NTSTART, NTSTOP )
       
              IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN
                     IF ( ABS(GET_TAUB() - GET_TAU()) < 1D-6 ) THEN
                     CALL MAKE_AVERAGE_XCO2EN
                     CALL WRITE_OSSE_TRUTH
                     ENDIF
                     RETURN
              ENDIF
       
              ! 初始化当前时段的诊断数组
              CURR_GC_XCO2EN  = 0d0
              CURR_OBS_XCO2EN = 0d0
              CURR_DIFF       = 0d0
              CURR_FORCING    = 0d0
              CURR_COST       = 0d0
              CURR_COUNT      = 0

              DO NT  = NTSTART, NTSTOP, -1
                     ! 质量控制
                     IF ( OBS(NT)%QF > 1 ) CYCLE
       
                     IIJJ  = GET_IJ( REAL(OBS(NT)%LON,4),
     &                    REAL(OBS(NT)%LAT,4))
                     I     = IIJJ(1)
                     J     = IIJJ(2)
                     
                     ! --- 正向模型计算 ---
                    TOTAL_AIR_MASS = SUM( AD(I,J,1:LLPAR) )
                    TOTAL_FF_MASS =SUM( CHK_STT(I,J,1:LLPAR,IDTCO2_FF))
                    SIM_XCO2EN  = ( TOTAL_FF_MASS / TOTAL_AIR_MASS ) 
     &                              * TCVV_CO2 * 1.0d6

                     ! Per-record H_delta(x_truth), before any averaging.
                     CALL STORE_OSSE_TRUTH( NT, I, J, GET_NHMS(),
     &                 SIM_XCO2EN, TOTAL_AIR_MASS, TOTAL_FF_MASS )

                     ! --- 代价函数与强迫计算 ---
                     DIFF_SCALAR  = SIM_XCO2EN - OBS(NT)%XCO2EN
                     OBS_ERR_SQ   = OBS(NT)%ERR ** 2
                     FORCE_SCALAR = DIFF_SCALAR / OBS_ERR_SQ
                     TMP_COST  = 0.5d0 * (DIFF_SCALAR**2) / OBS_ERR_SQ   
                     COST_FUNC    = COST_FUNC + TMP_COST

                     ! --- 累加诊断变量 ---
                     CURR_COUNT(I,J)     = CURR_COUNT(I,J) + 1
                     CURR_GC_XCO2EN(I,J)=CURR_GC_XCO2EN(I,J)+SIM_XCO2EN
                     CURR_OBS_XCO2EN(I,J)= CURR_OBS_XCO2EN(I,J) +
     &                                                    OBS(NT)%XCO2EN
                     CURR_DIFF(I,J)    = CURR_DIFF(I,J)   + DIFF_SCALAR
                     CURR_FORCING(I,J) = CURR_FORCING(I,J)+FORCE_SCALAR
                     CURR_COST(I,J)    = CURR_COST(I,J)     + TMP_COST

                     ALL_COUNT(I,J)      = ALL_COUNT(I,J) + 1
                     ALL_GC_XCO2EN(I,J)= ALL_GC_XCO2EN(I,J)+ SIM_XCO2EN
                     ALL_OBS_XCO2EN(I,J) = ALL_OBS_XCO2EN(I,J) +
     &                                                    OBS(NT)%XCO2EN
                     ALL_DIFF(I,J)       = ALL_DIFF(I,J)   + DIFF_SCALAR
                     ALL_FORCING(I,J) = ALL_FORCING(I,J) + FORCE_SCALAR
                     ALL_COST(I,J)       = ALL_COST(I,J)     + TMP_COST

                     ! --- 伴随反向传播 ---
                     ADJ_GRADIENT = FORCE_SCALAR * TCVV_CO2 * 1.0d6 
     &                      / TOTAL_AIR_MASS

                     DO L = 1, LLPAR
                     STT_ADJ(I,J,L,IDTCO2_FF) =STT_ADJ(I,J,L,IDTCO2_FF)
     &                                     + ADJ_GRADIENT
                     ENDDO
              ENDDO  ! NT
       
              IF ( FIRST ) FIRST = .FALSE.
              print*, ' Updated value of COST_FUNC = ', COST_FUNC
              print*, ' XCO2EN contribution = ', COST_FUNC - OLD_COST
       
              CALL MAKE_CURRENT_XCO2EN
              ! Keep a recoverable daily 1-D truth file as the backward
              ! simulation advances through each matched time window.
              CALL WRITE_OSSE_TRUTH

              IF ( ABS(GET_TAUB() - GET_TAU()) < 1D-6 ) THEN
                     CALL MAKE_AVERAGE_XCO2EN
                     CALL WRITE_OSSE_TRUTH
              ENDIF

       END SUBROUTINE CALC_XCO2EN_FORCE

!-----------------------------------------------------------------------

       SUBROUTINE GET_NT_RANGE( NTOT, HHMMSS, NTSTART, NTSTOP)
              USE TIME_MOD,    ONLY : YMD_EXTRACT
              INTEGER, INTENT(IN)  :: NTOT, HHMMSS
              INTEGER, INTENT(OUT) :: NTSTART, NTSTOP
              ! CHANGE 2: Initialize and clear the daily record cursor.
              INTEGER, SAVE        :: NTSAVE = 0
              LOGICAL              :: FOUND_ALL_RECORDS
              INTEGER              :: NTEST, HH, MM, SS
              REAL*8           :: GC_HH_FRAC, H1_FRAC, OBS_FRAC
              REAL*8           :: OBS_TIME_VAL,HH_OBS,MM_OBS,SS_OBS
       
              FOUND_ALL_RECORDS  = .FALSE.
              NTSTART = 0
              NTSTOP  = 0
              IF ( NTOT <= 0 ) THEN
                     NTSAVE = 0
                     RETURN
              ENDIF
              IF ( HHMMSS == 230000 ) NTSAVE = NTOT
              
              CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )
              GC_HH_FRAC = REAL(HH,8) / 24d0
              H1_FRAC    = 1d0 / 24d0
              IF ( NTSAVE == 0 ) RETURN
       
              OBS_TIME_VAL = OBS(NTSAVE)%TIME
              HH_OBS = INT( OBS_TIME_VAL / 10000d0 )
              MM_OBS = INT( (OBS_TIME_VAL - HH_OBS*10000d0) / 100d0 )
              SS_OBS = OBS_TIME_VAL - HH_OBS*10000d0 - MM_OBS*100d0
              OBS_FRAC = (HH_OBS + MM_OBS/60d0 + SS_OBS/3600d0) / 24d0
       
              IF ( OBS_FRAC + H1_FRAC < GC_HH_FRAC ) THEN
                     RETURN
              ELSEIF ( OBS_FRAC + H1_FRAC >= GC_HH_FRAC ) THEN
                     NTSTART = NTSAVE
                     NTEST   = NTSTART
                     DO WHILE ( FOUND_ALL_RECORDS == .FALSE. )
                            NTEST = NTEST - 1
                            IF ( NTEST == 0 ) THEN
                            NTSTOP = 1
                            FOUND_ALL_RECORDS = .TRUE.
                            NTSAVE = 0
                            ELSE
                     OBS_TIME_VAL = OBS(NTEST)%TIME
                     HH_OBS = INT( OBS_TIME_VAL / 10000d0 )
                     MM_OBS=INT((OBS_TIME_VAL-HH_OBS*10000d0)/100d0)
                     SS_OBS=OBS_TIME_VAL-HH_OBS*10000d0-MM_OBS*100d0
                     OBS_FRAC=(HH_OBS+MM_OBS/60d0+SS_OBS/3600d0)/24d0
                            IF ( OBS_FRAC + H1_FRAC < GC_HH_FRAC ) THEN
                                   NTSTOP = NTEST + 1
                                   FOUND_ALL_RECORDS = .TRUE.
                                   NTSAVE = NTEST
                            ENDIF
                            ENDIF
                     ENDDO
              ENDIF
       END SUBROUTINE GET_NT_RANGE
       
!-----------------------------------------------------------------------
       SUBROUTINE MAKE_CURRENT_XCO2EN
              USE NETCDF
              USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
              USE ADJ_ARRAYS_MOD,   ONLY: N_CALC, EXPAND_NAME
              USE TIME_MOD,         ONLY: EXPAND_DATE,GET_NYMD,GET_NHMS
              
              CHARACTER(LEN=255)     :: FILENAME, OUTPUT_FILE
              INTEGER    :: I, J, FID, LON_DIM_ID, LAT_DIM_ID
              REAL*4                 :: TMP_DATA(IIPAR,JJPAR,6) 
              
              FILENAME='gctm.xco2en.YYYYMMDD.hhmm.NN'
              CALL EXPAND_DATE(FILENAME, GET_NYMD(), GET_NHMS())
              CALL EXPAND_NAME(FILENAME, N_CALC)
              OUTPUT_FILE=TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
              
              CALL CHECK(NF90_CREATE(OUTPUT_FILE,NF90_CLOBBER,FID),0)
              CALL CHECK(NF90_DEF_DIM(FID,'lon',IIPAR,LON_DIM_ID),100)
              CALL CHECK(NF90_DEF_DIM(FID,'lat',JJPAR,LAT_DIM_ID),101)
              CALL CHECK(NF90_ENDDEF(FID), 8888)
              
              TMP_DATA = 0.0
              
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
              DO J = 1, JJPAR
              DO I = 1, IIPAR
                     IF (CURR_COUNT(I,J) > 0) THEN
                     TMP_DATA(I,J,1) = REAL(CURR_GC_XCO2EN(I,J),4)  / 
     &                              REAL(CURR_COUNT(I,J),4)
                     TMP_DATA(I,J,2) = REAL(CURR_OBS_XCO2EN(I,J),4) / 
     &                              REAL(CURR_COUNT(I,J),4)
                     TMP_DATA(I,J,3) = REAL(CURR_DIFF(I,J),4)       / 
     &                              REAL(CURR_COUNT(I,J),4)
                     TMP_DATA(I,J,4) = REAL(CURR_FORCING(I,J),4)
                     TMP_DATA(I,J,5) = REAL(CURR_COST(I,J),4)
                     TMP_DATA(I,J,6) = REAL(CURR_COUNT(I,J),4)
                     ENDIF
              ENDDO
              ENDDO
              !$OMP END PARALLEL DO
              
              CALL WRITE_NC_2D_FLOAT(FID,'gc_xco2en','Model XCO2_En', 
     & 'ppm', TMP_DATA(:,:,1), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'obs_xco2en','Obs XCO2_En', 
     & 'ppm', TMP_DATA(:,:,2), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'diff','Model - Obs', 
     & 'ppm', TMP_DATA(:,:,3), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'forcing','Adjoint Forcing', 
     & '1/ppm', TMP_DATA(:,:,4),IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'cost','Cost Function', 
     & 'unitless',TMP_DATA(:,:,5), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'count','# of obs', 
     & 'unitless',TMP_DATA(:,:,6), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              
              CALL CHECK(NF90_CLOSE(FID), 9999)
              END SUBROUTINE MAKE_CURRENT_XCO2EN

       !-----------------------------------------------------------------
              SUBROUTINE MAKE_AVERAGE_XCO2EN
              USE NETCDF
              USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
              USE ADJ_ARRAYS_MOD,    ONLY: N_CALC, EXPAND_NAME
              USE TIME_MOD,          ONLY: EXPAND_DATE, GET_NYMD
              
              CHARACTER(LEN=255)     :: FILENAME, OUTPUT_FILE
              INTEGER       :: I, J, FID, LON_DIM_ID, LAT_DIM_ID
              REAL*4                 :: TMP_DATA(IIPAR,JJPAR,6)
              
              FILENAME='gctm.xco2en.YYYYMMDD.NN'
              CALL EXPAND_DATE(FILENAME, GET_NYMD(), 9999)
              CALL EXPAND_NAME(FILENAME, N_CALC)
              OUTPUT_FILE=TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
              
              CALL CHECK(NF90_CREATE(OUTPUT_FILE, NF90_CLOBBER,FID),0)
              CALL CHECK(NF90_DEF_DIM(FID,'lon',IIPAR,LON_DIM_ID), 100)
              CALL CHECK(NF90_DEF_DIM(FID,'lat',JJPAR,LAT_DIM_ID), 101)
              CALL CHECK(NF90_ENDDEF(FID), 8888)
              
              TMP_DATA = 0.0
              
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
              DO J = 1, JJPAR
              DO I = 1, IIPAR
                     IF (ALL_COUNT(I,J) > 0) THEN
                     TMP_DATA(I,J,1) = REAL(ALL_GC_XCO2EN(I,J),4)  / 
     &                              REAL(ALL_COUNT(I,J),4)
                     TMP_DATA(I,J,2) = REAL(ALL_OBS_XCO2EN(I,J),4) / 
     &                              REAL(ALL_COUNT(I,J),4)
                     TMP_DATA(I,J,3) = REAL(ALL_DIFF(I,J),4)       / 
     &                              REAL(ALL_COUNT(I,J),4)
                     TMP_DATA(I,J,4) = REAL(ALL_FORCING(I,J),4)
                     TMP_DATA(I,J,5) = REAL(ALL_COST(I,J),4)
                     TMP_DATA(I,J,6) = REAL(ALL_COUNT(I,J),4)
                     ENDIF
              ENDDO
              ENDDO
              !$OMP END PARALLEL DO
              
              CALL WRITE_NC_2D_FLOAT(FID,'gc_xco2en','Model XCO2_En',
     & 'ppm', TMP_DATA(:,:,1), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'obs_xco2en','Obs XCO2_En',
     & 'ppm', TMP_DATA(:,:,2), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'diff','Model - Obs', 
     & 'ppm', TMP_DATA(:,:,3), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'forcing','Adjoint Forcing', 
     & '1/ppm', TMP_DATA(:,:,4), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'cost','Cost Function', 
     & 'unitless',TMP_DATA(:,:,5),IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              CALL WRITE_NC_2D_FLOAT(FID,'count','# of obs', 
     & 'unitless',TMP_DATA(:,:,6),IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
              
              CALL CHECK(NF90_CLOSE(FID), 9999)
              END SUBROUTINE MAKE_AVERAGE_XCO2EN

       !-----------------------------------------------------------------
              SUBROUTINE WRITE_OSSE_TRUTH
              USE NETCDF
              USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
              USE ADJ_ARRAYS_MOD,    ONLY: EXPAND_NAME
              USE TIME_MOD,          ONLY: EXPAND_DATE

              CHARACTER(LEN=255) :: FILENAME, OUTPUT_FILE
              INTEGER            :: FID, OBS_DIM_ID

              IF ( TRUTH_NMATCH <= 0 ) RETURN
              IF ( TRUTH_DATE < 0 .OR. TRUTH_N_CALC < 0 ) RETURN

              FILENAME='xco2en_osse_truth.YYYYMMDD.NN'
              CALL EXPAND_DATE(FILENAME, TRUTH_DATE, 9999)
              CALL EXPAND_NAME(FILENAME, TRUTH_N_CALC)
              OUTPUT_FILE=TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'

              CALL CHECK(NF90_CREATE(OUTPUT_FILE,NF90_CLOBBER,FID),
     &                 10000)
              CALL CHECK(NF90_DEF_DIM(FID,'nSamples',TRUTH_NMATCH,
     &                              OBS_DIM_ID),10001)
              CALL CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,'title',
     &           'Per-observation ML XCO2 enhancement OSSE truth'),
     &           10002)
              CALL CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,
     &           'truth_definition',
     &           'H_delta(x)=column mean of IDTCO2_FF tracer'),10003)
              CALL CHECK(NF90_PUT_ATT(FID,NF90_GLOBAL,
     &           'source_observation_file',
     &           TRIM(CURRENT_INPUT_FILE)),10004)
              CALL CHECK(NF90_ENDDEF(FID),10005)

              CALL WRITE_NC_1D_INT(FID,'observation_index',
     &           'Index in source observation file','1',
     &           TRUTH_LOBS(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'date','UTC date','YYYYMMDD',
     &           TRUTH_YYYYMMDD(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'obs_time',
     &           'Observation time from source file','HHMMSS UTC',
     &           TRUTH_OBS_TIME(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'gc_time',
     &           'GEOS-Chem time step used for sampling','HHMMSS UTC',
     &           TRUTH_GC_TIME(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'latitude','Latitude',
     &           'degrees_north',TRUTH_LAT(1:TRUTH_NMATCH),
     &           TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'longitude','Longitude',
     &           'degrees_east',TRUTH_LON(1:TRUTH_NMATCH),
     &           TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'model_i',
     &           'GEOS-Chem longitude index','1',
     &           TRUTH_I(1:TRUTH_NMATCH),TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'model_j',
     &           'GEOS-Chem latitude index','1',
     &           TRUTH_J(1:TRUTH_NMATCH),TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'quality_flag',
     &           'Quality flag from source file','1',
     &           TRUTH_QF(1:TRUTH_NMATCH),TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_INT(FID,'n_calc',
     &           'Optimizer evaluation number','1',
     &           TRUTH_CALC(1:TRUTH_NMATCH),TRUTH_NMATCH,OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'xco2en_true',
     &           'Modeled per-observation OSSE truth','ppm',
     &           TRUTH_XCO2EN(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'input_xco2en',
     &           'XCO2 enhancement in source observation file','ppm',
     &           TRUTH_INPUT_XCO2EN(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'uncertainty',
     &           'Observation uncertainty from source file','ppm',
     &           TRUTH_UNCERTAINTY(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'total_air_mass',
     &           'Air column mass used by observation operator','kg',
     &           TRUTH_AIR_MASS(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)
              CALL WRITE_NC_1D_DOUBLE(FID,'total_ff_co2_mass',
     &           'Fossil-fuel CO2 column mass used by operator','kg',
     &           TRUTH_FF_MASS(1:TRUTH_NMATCH),TRUTH_NMATCH,
     &           OBS_DIM_ID)

              CALL CHECK(NF90_CLOSE(FID),10999)
              WRITE(6,*) '    - WRITE_OSSE_TRUTH: wrote ',
     &                 TRIM(OUTPUT_FILE), ' records=', TRUTH_NMATCH

              END SUBROUTINE WRITE_OSSE_TRUTH

       !-----------------------------------------------------------------
              SUBROUTINE WRITE_NC_1D_INT(FID,NAMES,LONGNAMES,UNITS,
     &                                 VAR,N,DIM1)
              USE NETCDF
              INTEGER,      INTENT(IN) :: FID, N, DIM1
              CHARACTER(*), INTENT(IN) :: NAMES, LONGNAMES, UNITS
              INTEGER,      INTENT(IN) :: VAR(N)
              INTEGER                  :: VID

              CALL CHECK(NF90_REDEF(FID),11000)
              CALL CHECK(NF90_DEF_VAR(FID,TRIM(NAMES),NF90_INT,
     &                 (/DIM1/),VID),11001)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'long_name',
     &                 TRIM(LONGNAMES)),11002)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'units',TRIM(UNITS)),
     &                 11003)
              CALL CHECK(NF90_ENDDEF(FID),11004)
              CALL CHECK(NF90_PUT_VAR(FID,VID,VAR),11005)

              END SUBROUTINE WRITE_NC_1D_INT

       !-----------------------------------------------------------------
              SUBROUTINE WRITE_NC_1D_DOUBLE(FID,NAMES,LONGNAMES,UNITS,
     &                                    VAR,N,DIM1)
              USE NETCDF
              INTEGER,      INTENT(IN) :: FID, N, DIM1
              CHARACTER(*), INTENT(IN) :: NAMES, LONGNAMES, UNITS
              REAL*8,       INTENT(IN) :: VAR(N)
              INTEGER                  :: VID

              CALL CHECK(NF90_REDEF(FID),12000)
              CALL CHECK(NF90_DEF_VAR(FID,TRIM(NAMES),NF90_DOUBLE,
     &                 (/DIM1/),VID),12001)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'long_name',
     &                 TRIM(LONGNAMES)),12002)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'units',TRIM(UNITS)),
     &                 12003)
              CALL CHECK(NF90_ENDDEF(FID),12004)
              CALL CHECK(NF90_PUT_VAR(FID,VID,VAR),12005)

              END SUBROUTINE WRITE_NC_1D_DOUBLE

       !-----------------------------------------------------------------
              SUBROUTINE WRITE_NC_2D_FLOAT(FID,NAMES,LONGNAMES,
     &                              UNITS,VAR,IIN,JJN,DIM1,DIM2)
              USE NETCDF
              INTEGER,      INTENT(IN) :: FID, IIN, JJN, DIM1, DIM2
              CHARACTER(*), INTENT(IN) :: NAMES, LONGNAMES, UNITS
              REAL*4,       INTENT(IN) :: VAR(IIN,JJN)
              INTEGER                  :: VID
              
              CALL CHECK(NF90_REDEF(FID), 9000)
              CALL CHECK(NF90_DEF_VAR(FID,TRIM(NAMES),
     &                 NF90_FLOAT,(/DIM1,DIM2/),VID),9001)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'longname',
     &                 TRIM(LONGNAMES)),9002)
              CALL CHECK(NF90_PUT_ATT(FID,VID,'unit',TRIM(UNITS)),9003)
              CALL CHECK(NF90_ENDDEF(FID), 9004)
              CALL CHECK(NF90_PUT_VAR(FID,VID,VAR), 9005)
              END SUBROUTINE WRITE_NC_2D_FLOAT
       
       END MODULE XCO2EN_OBS_MOD
      
