!$Id: gosat_co2_mod.f,v 1.2 2011/02/23 00:08:48 daven Exp $
      MODULE GOSAT_CO2_MOD

            IMPLICIT NONE
      
            !=================================================================
            ! MODULE VARIABLES
            !=================================================================

#     include "CMN_SIZE"  

            ! Parameters
            INTEGER, PARAMETER           :: MAXLEV = 15
      
            ! Record to store data from each GOS obs
            TYPE GOS_CO2_OBS
                  INTEGER     :: LGOS      
                  REAL        :: LAT      
                  REAL        :: LON      
                  REAL        :: TIME     
                  INTEGER     :: QF        
                  REAL*8      :: CO2
                  REAL*8      :: PRES(MAXLEV)
                  REAL*8      :: PRES_WF(MAXLEV)
                  REAL*8      :: PRIOR
                  REAL*8      :: AVG_KERNEL(MAXLEV)
                  REAL*8      :: PRIOR_PROF(MAXLEV)
                  REAL*8      :: CO2_ERR
            END TYPE GOS_CO2_OBS
      
            TYPE(GOS_CO2_OBS), ALLOCATABLE :: GOS(:)
      
            ! IDTCO2 isn't defined in tracerid_mod because people just assume
            ! it is one. Define it here for now as a temporary patch.
            INTEGER, PARAMETER :: IDTCO2   = 1
            ! Same thing for TCVV(IDTCO2)
            REAL*8,  PARAMETER :: TCVV_CO2 = 28.97d0  / 44d0
            
            !-----------------------------------------------------------------
            ! New Arrays for Diagnostics (Current Window)
            !-----------------------------------------------------------------
            REAL*8 :: CURR_GC_CO2(IIPAR, JJPAR)   = 0d0
            REAL*8 :: CURR_OBS_CO2(IIPAR, JJPAR)  = 0d0
            REAL*8 :: CURR_DIFF(IIPAR, JJPAR)     = 0d0
            REAL*8 :: CURR_FORCING(IIPAR, JJPAR)  = 0d0
            REAL*8 :: CURR_COST(IIPAR, JJPAR)     = 0d0
            INTEGER:: CURR_COUNT(IIPAR, JJPAR)    = 0
            
            !-----------------------------------------------------------------
            ! New Arrays for Diagnostics (Average/Total)
            !-----------------------------------------------------------------
            REAL*8 :: ALL_GC_CO2(IIPAR, JJPAR)    = 0d0
            REAL*8 :: ALL_OBS_CO2(IIPAR, JJPAR)   = 0d0
            REAL*8 :: ALL_DIFF(IIPAR, JJPAR)      = 0d0
            REAL*8 :: ALL_FORCING(IIPAR, JJPAR)   = 0d0
            REAL*8 :: ALL_COST(IIPAR, JJPAR)      = 0d0
            INTEGER:: ALL_COUNT(IIPAR, JJPAR)     = 0

            CONTAINS
      !------------------------------------------------------------------------------
      
            SUBROUTINE READ_GOS_CO2_OBS( YYYYMMDD, NGOS )
      !
      !******************************************************************************
      !  Subroutine READ_GOS_CO2_OBS reads the file and passes back info contained
      !  therein. (dkh, 10/12/10)
      !******************************************************************************
      !
            ! Reference to f90 modules
            USE DIRECTORY_MOD,          ONLY : DATA_DIR
            USE NETCDF
            USE TIME_MOD,               ONLY : EXPAND_DATE
      
            IMPLICIT NONE 
      
            INTEGER, INTENT(IN)  :: YYYYMMDD
            INTEGER, INTENT(OUT) :: NGOS
            INTEGER, PARAMETER :: DP = KIND(1.0D0)
            REAL(DP), PARAMETER :: FILL_VAL = -999.0_DP
      
            ! --- NetCDF ID Variables ---
            INTEGER :: FID, NT_ID, AP_ID, QF_ID, ERR_ID, PF_ID, TM_ID
            INTEGER :: CO2_ID, PS_ID, LA_ID, LO_ID, AK_ID, PWF_ID
      
            ! --- Temporary Arrays ---
            REAL(DP), ALLOCATABLE :: TMP_LAT(:)
            REAL(DP), ALLOCATABLE :: TMP_LON(:)
            REAL(DP), ALLOCATABLE :: TMP_TIME(:)
            REAL(DP), ALLOCATABLE :: TMP_CO2(:)
            REAL(DP), ALLOCATABLE :: TMP_PRES(:,:)
            REAL(DP), ALLOCATABLE :: TMP_PWF(:,:)
            REAL(DP), ALLOCATABLE :: TMP_PF(:,:)
            REAL(DP), ALLOCATABLE :: TMP_AK(:,:) 
            REAL(DP), ALLOCATABLE :: TMP_PRIOR(:)    
            REAL(DP), ALLOCATABLE :: TMP_ERR(:)      
            INTEGER,  ALLOCATABLE :: TMP_QF(:)
      
            ! --- Loop Variables ---
            INTEGER :: N
            CHARACTER(LEN=255) :: READ_FILENAME
      
            !=================================================================
            ! READ_GOS_CO2_OBS begins here!
            !=================================================================
      
            ! filename root
            READ_FILENAME = TRIM( 'acos-v27-YYYYMMDD.nc' )
      
            ! Expand date tokens in filename
            CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )
      
            ! Construct complete filename
            READ_FILENAME = TRIM( './data/GOSAT-CO2/' ) //
     &                TRIM( READ_FILENAME )
      
            WRITE(6,*) '    - READ_GOSAT_CO2_OBS: reading file: ',
     &   READ_FILENAME
       
           ! Open file and assign file id (FID)
           CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE,FID),0)

           !--------------------------------
           ! Get data record IDs
           !--------------------------------
           CALL CHECK( NF90_INQ_DIMID( FID, "nSamples", NT_ID),  101)
           CALL CHECK( NF90_INQ_VARID( FID, "xco2",    CO2_ID ), 102)
           CALL CHECK( NF90_INQ_VARID( FID, "xco2_averaging_kernel", 
     &                                                  AK_ID ), 103)
           CALL CHECK( NF90_INQ_VARID( FID, "pressure", PS_ID ), 104)
           CALL CHECK( NF90_INQ_VARID( FID,"xco2_apriori",AP_ID),105)
           CALL CHECK( NF90_INQ_VARID( FID, "latitude",   LA_ID),106)
           CALL CHECK( NF90_INQ_VARID( FID, "longitude",  LO_ID),107)
           CALL CHECK( NF90_INQ_VARID( FID, "time",     TM_ID ), 108)
           CALL CHECK( NF90_INQ_VARID( FID, "xco2_quality_flag",
     &                                                   QF_ID), 109)
           CALL CHECK( NF90_INQ_VARID( FID, "xco2_uncertainty",
     &                                                  ERR_ID), 110)
           CALL CHECK( NF90_INQ_VARID( FID,"co2_profile_apriori",
     &                                                    PF_ID),111)
           CALL CHECK( NF90_INQ_VARID( FID,"pres_weight_func",
     &                                                   PWF_ID),112)
          
           ! READ number of retrievals, NGOS
           CALL CHECK(NF90_INQUIRE_DIMENSION(FID,NT_ID,LEN=NGOS),201)
     
           WRITE(6,*) ' Number of retrievals (NGOS): ', NGOS
           IF (ALLOCATED(GOS)) DEALLOCATE(GOS)
           ALLOCATE( GOS(NGOS) )
     
           !--------------------------------
           ! Read Data
           !--------------------------------
           ALLOCATE( TMP_LAT(NGOS), TMP_LON(NGOS), TMP_TIME(NGOS) )
           ALLOCATE( TMP_CO2(NGOS), TMP_PRIOR(NGOS), TMP_ERR(NGOS))
           ALLOCATE( TMP_QF(NGOS) )
           ALLOCATE( TMP_PRES(MAXLEV, NGOS) )
           ALLOCATE( TMP_PWF(MAXLEV, NGOS) )
           ALLOCATE( TMP_PF(MAXLEV, NGOS) )
           ALLOCATE( TMP_AK(MAXLEV, NGOS) )
     
           CALL CHECK( NF90_GET_VAR( FID, LA_ID,  TMP_LAT  ), 301 )
           CALL CHECK( NF90_GET_VAR( FID, LO_ID,  TMP_LON  ), 302 )
           CALL CHECK( NF90_GET_VAR( FID, TM_ID,  TMP_TIME ), 303 )
           CALL CHECK( NF90_GET_VAR( FID, CO2_ID, TMP_CO2  ), 304 )
           
           CALL CHECK( NF90_GET_VAR( FID, PS_ID,  TMP_PRES ), 305 )
           CALL CHECK( NF90_GET_VAR( FID, AK_ID,  TMP_AK   ), 306 )
           CALL CHECK( NF90_GET_VAR( FID, PF_ID,  TMP_PF   ), 307 )
     
           CALL CHECK( NF90_GET_VAR( FID, AP_ID,  TMP_PRIOR), 308 )
           CALL CHECK( NF90_GET_VAR( FID, QF_ID,  TMP_QF   ), 309 )
           CALL CHECK( NF90_GET_VAR( FID, ERR_ID, TMP_ERR  ), 310 )
           CALL CHECK( NF90_GET_VAR( FID, PWF_ID, TMP_PWF  ), 311 )
     
           CALL CHECK( NF90_CLOSE( FID ), 999 )
     
           DO N = 1, NGOS
                 GOS(N)%LGOS = N
                 
                 GOS(N)%LAT  = REAL( TMP_LAT(N),  KIND=4 )
                 GOS(N)%LON  = REAL( TMP_LON(N),  KIND=4 )
                 GOS(N)%TIME = REAL( TMP_TIME(N), KIND=4 )
                 GOS(N)%QF   = TMP_QF(N)
        
                 GOS(N)%CO2     = TMP_CO2(N)
                 GOS(N)%PRIOR   = TMP_PRIOR(N)
                 GOS(N)%CO2_ERR = TMP_ERR(N)
        
                 GOS(N)%PRES(:)       = TMP_PRES(:, N)
                 GOS(N)%PRES_WF(:)    = TMP_PWF(:, N)
                 GOS(N)%AVG_KERNEL(:) = TMP_AK(:, N)
                 GOS(N)%PRIOR_PROF(:) = TMP_PF(:, N)
           END DO
     
           DEALLOCATE( TMP_LAT, TMP_LON, TMP_TIME, TMP_CO2 )
           DEALLOCATE( TMP_PRIOR, TMP_ERR, TMP_QF, TMP_PWF )
           DEALLOCATE( TMP_PRES, TMP_AK, TMP_PF)
     
           END SUBROUTINE READ_GOS_CO2_OBS
      !------------------------------------------------------------------------------
      
            SUBROUTINE CHECK( STATUS, LOCATION )
      !
      !******************************************************************************
      ! Subroutine CHECK checks the status of calls to netCDF libraries routines
      ! (dkh, 02/15/09)
      !******************************************************************************
      !
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
              CALL ERROR_STOP('netCDF error', 'gosat_co2_mod')
            ENDIF
      
            ! Return to calling program
            END SUBROUTINE CHECK
      
      !------------------------------------------------------------------------------
      
            SUBROUTINE CALC_GOS_CO2_FORCE( COST_FUNC )
      !
      !******************************************************************************
      !  Subroutine CALC_GOS_CO2_FORCE calculates the adjoint forcing from the GOSAT
      !  CO2 observations and updates the cost function. (dkh, 10/12/10)
      !******************************************************************************
      !
            ! Reference to f90 modules
            USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
            USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
            USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
            USE CHECKPT_MOD,        ONLY : CHK_STT
            USE COMODE_MOD,         ONLY : CSPEC, JLOP
            USE DAO_MOD,            ONLY : AD
            USE DAO_MOD,            ONLY : AIRDEN
            USE DAO_MOD,            ONLY : BXHEIGHT
            USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
            USE GRID_MOD,           ONLY : GET_IJ
            USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
            USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
            USE TIME_MOD,           ONLY : GET_TS_CHEM
            USE TIME_MOD,           ONLY : GET_TAUB, GET_TAU
            USE TRACER_MOD,         ONLY : XNUMOLAIR
            USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      
            ! Size params (Already included via header in module, but kept for safety)
            ! # include      "CMN_SIZE"      
      
            ! Arguments
            REAL*8, INTENT(INOUT)       :: COST_FUNC
      
            ! Local variables
            INTEGER                     :: NTSTART, NTSTOP, NT
            INTEGER                     :: IIJJ(2), I,      J
            INTEGER                     :: L,       LL,     LGOS
            INTEGER                     :: JLOOP
            REAL*8                      :: GC_PRES(LLPAR)
            REAL*8                      :: GC_CO2_NATIVE(LLPAR)
            REAL*8                      :: GC_CO2(MAXLEV)
            REAL*8                      :: GC_PSURF
            REAL*8                      :: MAP(LLPAR,MAXLEV)
            REAL*8                      :: CO2_HAT(MAXLEV)
            REAL*8                      :: CO2_PERT(MAXLEV)
            REAL*8                      :: FORCE(MAXLEV)
            REAL*8                      :: DIFF(MAXLEV)
            REAL*8                      :: XCO2_HAT
            REAL*8, allocatable         :: NEW_COST(:)
            REAL*8                      :: OLD_COST, TMP_COST
            INTEGER,SAVE                :: NGOS
            REAL*8                      :: DIFF_SCALAR   
            REAL*8                      :: FORCE_SCALAR   
            REAL*8                      :: OBS_ERROR_SQ   
      
            REAL*8                      :: GC_CO2_NATIVE_ADJ(LLPAR)
            REAL*8                      :: CO2_HAT_ADJ(MAXLEV)
            REAL*8                      :: CO2_PERT_ADJ(MAXLEV)
            REAL*8                      :: GC_CO2_ADJ(MAXLEV)
            REAL*8                      :: DIFF_ADJ(MAXLEV)
      
            LOGICAL, SAVE               :: FIRST = .TRUE.
            INTEGER                     :: IOS
            CHARACTER(LEN=255)          :: FILENAME
      
      
            !=================================================================
            ! CALC_GOS_CO2_FORCE begins here!
            !=================================================================
      
            print*, '     - CALC_GOS_CO2_FORCE '
      
            ! Save a value of the cost function first
            OLD_COST = COST_FUNC
      
            ! Check if it is the last hour of a day
            IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN
               ! Read the GOS CO2 file for this day
               CALL READ_GOS_CO2_OBS( GET_NYMD(), NGOS )
            ENDIF
      
            ! Get the range of GOS retrievals for the current hour
            CALL GET_NT_RANGE( NGOS, GET_NHMS(), NTSTART, NTSTOP )
      
            IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN
               print*, ' No matching GOS CO2 obs for this hour'
               ! Even if no obs, we might need to write empty diagnostics or averages
               IF ( ABS(GET_TAUB() - GET_TAU()) < 1D-6 ) THEN
                  CALL MAKE_AVERAGE_GOS_CO2
               ENDIF
               RETURN
            ENDIF
      
            IF ( ALLOCATED(NEW_COST) ) DEALLOCATE(NEW_COST)
            ALLOCATE( NEW_COST(NGOS) )
            NEW_COST = 0D0
            
            !------------------------------------------------------------
            ! Initialize Current Diag Arrays
            !------------------------------------------------------------
            CURR_GC_CO2  = 0d0
            CURR_OBS_CO2 = 0d0
            CURR_DIFF    = 0d0
            CURR_FORCING = 0d0
            CURR_COST    = 0d0
            CURR_COUNT   = 0

            DO NT  = NTSTART, NTSTOP, -1
               ! quality screening
               IF ( GOS(NT)%QF > 1 ) THEN
                  print*, ' BAD QF, skipping record ', NT
                  CYCLE
               ENDIF
      
               ! skip antarctica
               IF ( GOS(NT)%LAT < -60d0 ) THEN
                  print*, ' Skipp data with latitude < 60 S ', NT
                  CYCLE
               ENDIF
      
               ! For safety, initialize these up to LLGOS
               GC_CO2(:)       = 0d0
               MAP(:,:)        = 0d0
               CO2_HAT_ADJ(:)  = 0d0
               FORCE(:)        = 0d0
      
               ! Copy LGOS to make coding a bit cleaner
               LGOS = MAXLEV
      
               ! Get grid box of current record
               IIJJ  = GET_IJ( REAL(GOS(NT)%LON,4), REAL(GOS(NT)%LAT,4))
               I     = IIJJ(1)
               J     = IIJJ(2)
      
               ! Get GC pressure levels (mbar)
               DO L = 1, LLPAR
                  GC_PRES(L) = GET_PCENTER(I,J,L)
               ENDDO
      
               ! Get GC surface pressure (mbar)
               GC_PSURF = GET_PEDGE(I,J,1)
      
               ! Calculate the interpolation weight matrix
               MAP(1:LLPAR,1:LGOS)
     &         = GET_INTMAP( LLPAR, GC_PRES(:),           GC_PSURF,
     &                    LGOS,  GOS(NT)%PRES(1:LGOS), GC_PSURF  )
      
               ! Get CO2 values at native model resolution
               GC_CO2_NATIVE(:) = CHK_STT(I,J,:,IDTCO2)
      
               ! Convert from kg/box to ppm
               GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2 
     &                         * 1.0d6 / AD(I,J,:)
      
               ! Interpolate GC CO2 column to TES grid
               DO LL = 1, LGOS
                  GC_CO2(LL) = 0d0
                  DO L = 1, LLPAR
                     GC_CO2(LL) = GC_CO2(LL)
     &                    + MAP(L,LL) * GC_CO2_NATIVE(L)
                  ENDDO
               ENDDO
      
               !--------------------------------------------------------------
               ! Apply GOS observation operator
               !--------------------------------------------------------------
               XCO2_HAT = GOS(NT)%PRIOR
               DO L = 1, LGOS
                  CO2_PERT(L) = GC_CO2(L) - GOS(NT)%PRIOR_PROF(L)
                  XCO2_HAT = XCO2_HAT+GOS(NT)%AVG_KERNEL(L) * 
     &                       GOS(NT)%PRES_WF(L) * CO2_PERT(L)
               ENDDO
               
               !--------------------------------------------------------------
               ! Calculate cost function, given S is error in vmr
               ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
               !--------------------------------------------------------------
               
               DIFF_SCALAR = XCO2_HAT - GOS(NT)%CO2
               OBS_ERROR_SQ = GOS(NT)%CO2_ERR ** 2
      
               !--------------------------------------------------------------
               ! Begin adjoint calculations
               !--------------------------------------------------------------
               FORCE_SCALAR = DIFF_SCALAR / OBS_ERROR_SQ
               TMP_COST = 0.5d0 * (DIFF_SCALAR**2) / OBS_ERROR_SQ   
               NEW_COST(NT)=NEW_COST(NT)+TMP_COST
               
               !--------------------------------------------------------------
               ! Accumulate Diagnostics
               !--------------------------------------------------------------
               CURR_COUNT(I,J)   = CURR_COUNT(I,J) + 1
               CURR_GC_CO2(I,J)  = CURR_GC_CO2(I,J)  + XCO2_HAT
               CURR_OBS_CO2(I,J) = CURR_OBS_CO2(I,J) + GOS(NT)%CO2
               CURR_DIFF(I,J)    = CURR_DIFF(I,J)    + DIFF_SCALAR
               CURR_FORCING(I,J) = CURR_FORCING(I,J) + FORCE_SCALAR
               CURR_COST(I,J)    = CURR_COST(I,J)    + TMP_COST

               ! Accumulate to All-Time Arrays
               ALL_COUNT(I,J)    = ALL_COUNT(I,J) + 1
               ALL_GC_CO2(I,J)   = ALL_GC_CO2(I,J)  + XCO2_HAT
               ALL_OBS_CO2(I,J)  = ALL_OBS_CO2(I,J) + GOS(NT)%CO2
               ALL_DIFF(I,J)     = ALL_DIFF(I,J)    + DIFF_SCALAR
               ALL_FORCING(I,J)  = ALL_FORCING(I,J) + FORCE_SCALAR
               ALL_COST(I,J)     = ALL_COST(I,J)    + TMP_COST
               
               !--------------------------------------------------------------
               ! Adjoint propagation
               !--------------------------------------------------------------
               GC_CO2_ADJ(:) = 0d0
      
               DO L = 1, LGOS          
                  GC_CO2_ADJ(L) = FORCE_SCALAR * GOS(NT)%AVG_KERNEL(L) 
     &                         * 1.0d6 * GOS(NT)%PRES_WF(L)            
               ENDDO
      
               GC_CO2_NATIVE_ADJ(:) = 0d0
      
               DO L = 1, LLPAR
                  DO LL = 1, LGOS
                     GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L) 
     &               + MAP(L,LL) * GC_CO2_ADJ(LL)
                  ENDDO
               ENDDO
      
               GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) 
     &    * TCVV_CO2 / AD(I,J,:)
               STT_ADJ(I,J,:,IDTCO2) = STT_ADJ(I,J,:,IDTCO2) 
     &    + GC_CO2_NATIVE_ADJ(:)
      
            ENDDO  ! NT
      
            ! Update cost function
            COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))
      
            IF ( FIRST ) FIRST = .FALSE.
            print*, ' Updated value of COST_FUNC = ', COST_FUNC
            print*, ' GOS contribution = ', COST_FUNC - OLD_COST
      
            !--------------------------------------------------------------
            ! Output Current Diagnostics
            !--------------------------------------------------------------
            CALL MAKE_CURRENT_GOS_CO2

            !--------------------------------------------------------------
            ! Check for End of Simulation to Output Averages
            !--------------------------------------------------------------
            IF ( ABS(GET_TAUB() - GET_TAU()) < 1D-6 ) THEN
                  CALL MAKE_AVERAGE_GOS_CO2
            ENDIF

            ! Return to calling program
            END SUBROUTINE CALC_GOS_CO2_FORCE
      
      !------------------------------------------------------------------------------
      
            SUBROUTINE GET_NT_RANGE( NTES, HHMMSS, NTSTART, NTSTOP)
      !
      !******************************************************************************
      !  Subroutine GET_NT_RANGE returns range of records for current hour.
      !  Modified to access GOS directly from module variables.
      !  Assumes GOS(N)%TIME is HHMMSS (e.g. 130000).
      !******************************************************************************
            USE ERROR_MOD,    ONLY : ERROR_STOP
            USE TIME_MOD,     ONLY : YMD_EXTRACT
      
            ! Arguments
            INTEGER, INTENT(IN)   :: NTES      ! Total records (NGOS)
            INTEGER, INTENT(IN)   :: HHMMSS    ! Model Time (e.g. 130000)
            INTEGER, INTENT(OUT)  :: NTSTART, NTSTOP
      
            ! Local variables
            INTEGER, SAVE         :: NTSAVE
            LOGICAL               :: FOUND_ALL_RECORDS
            INTEGER               :: NTEST
            INTEGER               :: HH, MM, SS
            REAL*8                :: GC_HH_FRAC, H1_FRAC
            REAL*8                :: OBS_FRAC
            REAL*8                :: OBS_TIME_VAL
            REAL*8                :: HH_OBS, MM_OBS, SS_OBS
      
            ! Initialize
            FOUND_ALL_RECORDS  = .FALSE.
            NTSTART            = 0
            NTSTOP             = 0
      
            ! Reset NTSAVE at start of day (Assuming backward run, 230000 is first)
            IF ( HHMMSS == 230000 ) NTSAVE = NTES
      
            ! Model time fraction
            CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )
            GC_HH_FRAC = REAL(HH,8) / 24d0
            H1_FRAC    = 1d0 / 24d0
      
            IF ( NTSAVE == 0 ) THEN
                  print*, 'All records processed.'
                  RETURN
            ENDIF
      
            ! Helper internal function to calculate fraction from HHMMSS
            ! Using current NTSAVE
            OBS_TIME_VAL = GOS(NTSAVE)%TIME
            HH_OBS = INT( OBS_TIME_VAL / 10000d0 )
            MM_OBS = INT( (OBS_TIME_VAL - HH_OBS*10000d0) / 100d0 )
            SS_OBS = OBS_TIME_VAL - HH_OBS*10000d0 - MM_OBS*100d0
            OBS_FRAC = (HH_OBS + MM_OBS/60d0 + SS_OBS/3600d0) / 24d0
      
            ! Check if current NTSAVE is too early (smaller fraction) than model time window
            ! Model window is [GC_HH_FRAC - 1h, GC_HH_FRAC]
            ! If ObsTime + 1h < ModelTime, then ObsTime < ModelTime - 1h.
            IF ( OBS_FRAC + H1_FRAC < GC_HH_FRAC ) THEN
                  print*, 'No records reached yet'
                  RETURN
            ELSEIF ( OBS_FRAC + H1_FRAC >= GC_HH_FRAC ) THEN
                  ! Found valid record (Start from high index)
                  NTSTART = NTSAVE
                  NTEST   = NTSTART
                  
                  DO WHILE ( FOUND_ALL_RECORDS == .FALSE. )
                  NTEST = NTEST - 1
                  
                  IF ( NTEST == 0 ) THEN
                        NTSTOP = 1
                        FOUND_ALL_RECORDS = .TRUE.
                        NTSAVE = 0 ! All done
                  ELSE
                        ! Check next record time
                        OBS_TIME_VAL = GOS(NTEST)%TIME
                        HH_OBS = INT( OBS_TIME_VAL / 10000d0 )
                        MM_OBS=INT((OBS_TIME_VAL-HH_OBS*10000d0)/100d0)
                        SS_OBS=OBS_TIME_VAL-HH_OBS*10000d0-MM_OBS*100d0
                        OBS_FRAC=(HH_OBS+MM_OBS/60d0+SS_OBS/3600d0)/24d0
                        
                        ! If this record is outside window (too early)
                        IF ( OBS_FRAC + H1_FRAC < GC_HH_FRAC ) THEN
                        NTSTOP = NTEST + 1
                        FOUND_ALL_RECORDS = .TRUE.
                        NTSAVE = NTEST
                        ENDIF
                  ENDIF
                  ENDDO
            ENDIF
      
            END SUBROUTINE GET_NT_RANGE
      
      !------------------------------------------------------------------------------
      
            FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )
     &         RESULT      ( HINTERPZ )
      !
      !******************************************************************************
      !  Function GET_INTMAP linearly interpolates column quatities
      !  based upon the centered (average) pressue levels.
      !******************************************************************************
      !
            ! Reference to f90 modules
            USE ERROR_MOD,     ONLY : ERROR_STOP
            USE PRESSURE_MOD,  ONLY : GET_BP
      
            ! Arguments
            INTEGER            :: LGC_TOP, LTM_TOP
            REAL*8             :: GC_PRESC(LGC_TOP)
            REAL*8             :: TM_PRESC(LTM_TOP)
            REAL*8             :: GC_SURFP
            REAL*8             :: TM_SURFP
      
            ! Return value
            REAL*8             :: HINTERPZ(LGC_TOP, LTM_TOP)
      
            ! Local variables
            INTEGER  :: LGC, LTM
            REAL*8   :: DIFF, DELTA_SURFP
            REAL*8   :: LOW, HI
      
            !=================================================================
            ! GET_HINTERPZ_2 begins here!
            !=================================================================
      
            HINTERPZ(:,:) = 0D0
      
            ! Rescale GC grid according to TM surface pressure
            DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
      
            ! Loop over each pressure level of TM grid
            DO LTM = 1, LTM_TOP
      
               ! Find the levels from GC that bracket level LTM
               DO LGC = 1, LGC_TOP - 1
      
                  LOW = GC_PRESC(LGC+1)
                  HI  = GC_PRESC(LGC)
                  IF (LGC == 0) HI = TM_SURFP
      
                  ! Linearly interpolate value on the LTM grid
                  IF ( TM_PRESC(LTM) <= HI .and.
     &           TM_PRESC(LTM)  > LOW) THEN
      
                     DIFF                = HI - LOW
                     HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM) ) / DIFF
                     HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW) / DIFF
      
      
                  ENDIF
      
                ENDDO
             ENDDO
      
             ! Bug fix:  a more general version allows for multiples TES pressure
             ! levels to exist below the lowest GC pressure.  (dm, dkh, 09/30/10)
             ! New code:
             ! Loop over each pressure level of TM grid
             DO LTM = 1, LTM_TOP
                IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
                   HINTERPZ(1,LTM)         = 1D0
                   HINTERPZ(2:LGC_TOP,LTM) = 0D0
                ENDIF
             ENDDO
      
            ! Return to calling program
            END FUNCTION GET_INTMAP

      !------------------------------------------------------------------------------
      ! NEW SUBROUTINES FOR NETCDF OUTPUT
      !------------------------------------------------------------------------------

            SUBROUTINE MAKE_CURRENT_GOS_CO2
            !******************************************************************
            ! This subroutine outputs diagnostic data for current assimilation 
            ! time window. Adapted from tropomi_no2_obs_mod.
            !******************************************************************
            USE NETCDF
            USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
            USE ADJ_ARRAYS_MOD,    ONLY: N_CALC
            USE ADJ_ARRAYS_MOD,    ONLY: EXPAND_NAME
            USE TIME_MOD,          ONLY: EXPAND_DATE
            USE TIME_MOD,          ONLY: GET_NYMD, GET_NHMS
            
            ! Local variables
            CHARACTER(LEN=255)     :: FILENAME, OUTPUT_FILE
            INTEGER                :: I, J
            INTEGER                :: FID, LON_DIM_ID, LAT_DIM_ID
            REAL*4                 :: TMP_DATA(IIPAR,JJPAR,6) ! 6 variables
            
            FILENAME='gctm.gosat.co2.YYYYMMDD.hhmm.NN'
            FILENAME=TRIM(FILENAME)
            
            ! Expand date tokens
            CALL EXPAND_DATE(FILENAME, GET_NYMD(), GET_NHMS())
            
            ! Expand iteration number
            CALL EXPAND_NAME(FILENAME, N_CALC)
            
            OUTPUT_FILE=TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
            
            WRITE(6,*) '- MAKE_CURRENT_GOS_CO2: writing ', 
     &       TRIM(OUTPUT_FILE)
            
            ! Open file
            CALL CHECK(NF90_CREATE(OUTPUT_FILE,NF90_CLOBBER,FID),0)
            
            ! Define dimensions
            CALL CHECK(NF90_DEF_DIM(FID,'lon',IIPAR,LON_DIM_ID),100)
            CALL CHECK(NF90_DEF_DIM(FID,'lat',JJPAR,LAT_DIM_ID),101)
            
            ! End define mode
            CALL CHECK(NF90_ENDDEF(FID), 8888)
            
            TMP_DATA = 0.0
            
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
            DO J = 1, JJPAR
              DO I = 1, IIPAR
                IF (CURR_COUNT(I,J) > 0) THEN
                  ! Average for concentrations
                  TMP_DATA(I,J,1) = REAL(CURR_GC_CO2(I,J),4)  / 
     &             REAL(CURR_COUNT(I,J),4)
                  TMP_DATA(I,J,2) = REAL(CURR_OBS_CO2(I,J),4) / 
     &             REAL(CURR_COUNT(I,J),4)
                  TMP_DATA(I,J,3) = REAL(CURR_DIFF(I,J),4)    / 
     &             REAL(CURR_COUNT(I,J),4)
                  ! Sum for forcing and cost (or as provided in TROPOMI example, raw accumulation)
                  ! TROPOMI code stores forcing and cost directly without division in tmp_no2
                  TMP_DATA(I,J,4) = REAL(CURR_FORCING(I,J),4)
                  TMP_DATA(I,J,5) = REAL(CURR_COST(I,J),4)
                  TMP_DATA(I,J,6) = REAL(CURR_COUNT(I,J),4)
                ENDIF
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
            
            ! Save variables
            CALL WRITE_NC_2D_FLOAT(FID,'gc_xco2',
     &       'GEOS-Chem XCO2 (with AK)', 
     & 'ppm', TMP_DATA(:,:,1), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'gosat_xco2','GOSAT XCO2', 
     & 'ppm', TMP_DATA(:,:,2), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'diff','Model - Obs', 
     & 'ppm', TMP_DATA(:,:,3), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'forcing',
     &       'Adjoint Forcing (diff/err^2)', '1/ppm', 
     &         TMP_DATA(:,:,4),IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,
     & 'cost','Cost Function', 'unitless', 
     &        TMP_DATA(:,:,5), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,
     & 'count','# of obs', 'unitless', 
     &       TMP_DATA(:,:,6), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
            
            CALL CHECK(NF90_CLOSE(FID), 9999)
            
            END SUBROUTINE MAKE_CURRENT_GOS_CO2

            SUBROUTINE MAKE_AVERAGE_GOS_CO2
            !******************************************************************
            ! This subroutine outputs diagnostic data for the whole simulation.
            !******************************************************************
            USE NETCDF
            USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
            USE ADJ_ARRAYS_MOD,    ONLY: N_CALC
            USE ADJ_ARRAYS_MOD,    ONLY: EXPAND_NAME
            USE TIME_MOD,          ONLY: EXPAND_DATE
            USE TIME_MOD,          ONLY: GET_NYMD
            
            ! Local variables
            CHARACTER(LEN=255)     :: FILENAME, OUTPUT_FILE
            INTEGER                :: I, J
            INTEGER                :: FID, LON_DIM_ID, LAT_DIM_ID
            REAL*4                 :: TMP_DATA(IIPAR,JJPAR,6)
            
            FILENAME='gctm.gosat.co2.YYYYMMDD.NN'
            FILENAME=TRIM(FILENAME)
            
            ! Expand date tokens
            CALL EXPAND_DATE(FILENAME, GET_NYMD(), 9999)
            
            ! Expand iteration number
            CALL EXPAND_NAME(FILENAME, N_CALC)
            
            OUTPUT_FILE=TRIM(DIAGADJ_DIR)//TRIM(FILENAME)//'.nc'
            
            WRITE(6,*) '- MAKE_AVERAGE_GOS_CO2:', TRIM(OUTPUT_FILE)
            
            ! Open file
            CALL CHECK(NF90_CREATE(OUTPUT_FILE, NF90_CLOBBER,FID),0)
            
            ! Define dimensions
            CALL CHECK(NF90_DEF_DIM(FID,'lon',IIPAR,LON_DIM_ID), 100)
            CALL CHECK(NF90_DEF_DIM(FID,'lat',JJPAR,LAT_DIM_ID), 101)
            
            ! End define mode
            CALL CHECK(NF90_ENDDEF(FID), 8888)
            
            TMP_DATA = 0.0
            
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
            DO J = 1, JJPAR
              DO I = 1, IIPAR
                IF (ALL_COUNT(I,J) > 0) THEN
                  TMP_DATA(I,J,1) = REAL(ALL_GC_CO2(I,J),4)  / 
     &             REAL(ALL_COUNT(I,J),4)
                  TMP_DATA(I,J,2) = REAL(ALL_OBS_CO2(I,J),4) / 
     &             REAL(ALL_COUNT(I,J),4)
                  TMP_DATA(I,J,3) = REAL(ALL_DIFF(I,J),4)    / 
     &             REAL(ALL_COUNT(I,J),4)
                  TMP_DATA(I,J,4) = REAL(ALL_FORCING(I,J),4)
                  TMP_DATA(I,J,5) = REAL(ALL_COST(I,J),4)
                  TMP_DATA(I,J,6) = REAL(ALL_COUNT(I,J),4)
                ENDIF
              ENDDO
            ENDDO
            !$OMP END PARALLEL DO
            
            ! Save variables
            CALL WRITE_NC_2D_FLOAT(FID,'gc_xco2',
     &       'GEOS-Chem XCO2 (with AK)', 
     & 'ppm', TMP_DATA(:,:,1), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'gosat_xco2','GOSAT XCO2',
     & 'ppm', TMP_DATA(:,:,2), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'diff','Model - Obs', 
     & 'ppm', TMP_DATA(:,:,3), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,'forcing',
     &       'Adjoint Forcing (diff/err^2)', 
     & '1/ppm', TMP_DATA(:,:,4), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,
     &      'cost','Cost Function', 'unitless', 
     & TMP_DATA(:,:,5), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
                 
            CALL WRITE_NC_2D_FLOAT(FID,
     &      'count','# of obs', 'unitless', 
     & TMP_DATA(:,:,6), IIPAR, JJPAR, LON_DIM_ID, LAT_DIM_ID)
            
            CALL CHECK(NF90_CLOSE(FID), 9999)
            
            END SUBROUTINE MAKE_AVERAGE_GOS_CO2

            SUBROUTINE WRITE_NC_2D_FLOAT(FID,NAMES,LONGNAMES,UNITS,VAR,
     &   IIN,JJN,DIM1,DIM2)
            !******************************************************************
            ! Helper to write 2D float variables to NetCDF
            !******************************************************************
            USE NETCDF
            
            INTEGER,      INTENT(IN) :: FID
            INTEGER,      INTENT(IN) :: IIN, JJN
            INTEGER,      INTENT(IN) :: DIM1, DIM2
            CHARACTER(*), INTENT(IN) :: NAMES, LONGNAMES, UNITS
            REAL*4,       INTENT(IN) :: VAR(IIN,JJN)
            
            INTEGER                  :: VID
            
            ! Open define mode
            CALL CHECK(NF90_REDEF(FID), 9000)
            
            ! Define variable
            CALL CHECK(NF90_DEF_VAR(FID,TRIM(NAMES),
     &       NF90_FLOAT,(/DIM1,DIM2/),VID),9001)
            
            ! Define attributes
            CALL CHECK(NF90_PUT_ATT(FID,VID,'longname',
     &       TRIM(LONGNAMES)),9002)
            CALL CHECK(NF90_PUT_ATT(FID,VID,'unit',TRIM(UNITS)),9003)
            
            ! End define mode
            CALL CHECK(NF90_ENDDEF(FID), 9004)
            
            ! Save variable
            CALL CHECK(NF90_PUT_VAR(FID,VID,VAR), 9005)
            
            END SUBROUTINE WRITE_NC_2D_FLOAT
      
      !------------------------------------------------------------------------------
      
      END MODULE GOSAT_CO2_MOD