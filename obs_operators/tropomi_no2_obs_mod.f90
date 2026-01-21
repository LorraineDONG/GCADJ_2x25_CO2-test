! $Id: tropomi_no2_obs_mod.f90
module TROPOMI_NO2_OBS_MOD

!********************************************************************
! This module is used to read TROPOMI L2 NO2 observations and compute
! adjoint forcing and cost function. The TROPOMI NO2 data is preprocessed
! by Python code.
!
! zhendlu (zhendong-lu@uiowa.edu)
! 01/23/2021
!********************************************************************

implicit none

  ! header files
# include "define.h"
# include "CMN_SIZE"
    
  ! keep internal variables and subroutines from being seen outside
  ! 'tropomi_no2_obs_mod.f90', except calc_tropomi_no2_force.

  private

  public::calc_tropomi_no2_force

  !==================================================================
  ! module variables
  !==================================================================
  
  ! constants
  integer, parameter  :: time_window=60 ! unit:min
  integer, parameter  :: levels=47 ! number of vertical levels in TROPOMI POMINO

  ! variables
  type tropomi_no2_obs
    real   :: lon     ! longitude, degree_east
    real   :: lat     ! latitude, degree_north
    real*8 :: time    ! time since 2010-01-01, second
    real*8 :: NO2     ! TROPOMI tropospheric NO2 VCD, molec/cm2
    real*8 :: err     ! NO2 VCD error, molec/cm2
    real   :: amf     ! TROPOMI air mass factor, unitless
    real   :: ak(levels)       ! averaging kernel, unitless
    real   :: surf_P(levels)   ! surface pressure, hPa

  end type tropomi_no2_obs

  type(tropomi_no2_obs), allocatable::tropomi_no2(:)

  integer  :: N_NO2  ! number of observations for current day
  integer  :: N_CURR ! number of observations for current time window

  ! the flags array marks whether or not a specific observation is in the
  ! current time window and simulation domain
  logical, allocatable :: flags(:)

  ! GC NO2 vertical column for current time window, molec/cm2
  real*8   :: curr_gc_no2(IIPAR,JJPAR)
  ! GC NO2 slant column for current time window, molec/cm2
  real*8   :: curr_gc_s_no2(IIPAR,JJPAR)
  ! GC NO2 tropospheric air mass factor for current time window, unitless
  real*8   :: curr_gc_amf(IIPAR,JJPAR)
  ! TROPOMI NO2 vertical column for current time window, molec/cm2
  real*8   :: curr_no2(IIPAR,JJPAR)
  ! TROPOMI NO2 slant column for current time window, molec/cm2
  real*8   :: curr_s_no2(IIPAR,JJPAR)
  ! TROPOMI NO2 VCD applying averaging kernel and GC NO2 profile, molec/cm2
  real*8   :: curr_vak_no2(IIPAR,JJPAR)
  ! TROPOMI NO2 tropospheric air mass factor for current time window, unitless
  real*8   :: curr_amf(IIPAR,JJPAR)
  ! difference between current GC and TROPOMI NO2 column, molec/cm2
  real*8   :: curr_diff_no2(IIPAR,JJPAR)
  ! adjoint forcing in current time window, cm2/molec
  real*8   :: curr_forcing(IIPAR,JJPAR)
  ! cost function in current time window, unitless
  real*8   :: curr_cost(IIPAR,JJPAR)
  ! number of observations in current time window
  integer  :: curr_count(IIPAR,JJPAR)

  ! same as above but for the whole simulation period
  real*8   :: all_gc_no2(IIPAR,JJPAR)   = 0D0
  real*8   :: all_gc_s_no2(IIPAR,JJPAR) = 0D0
  real*8   :: all_gc_amf(IIPAR,JJPAR)   = 0D0
  real*8   :: all_no2(IIPAR,JJPAR)      = 0D0
  real*8   :: all_s_no2(IIPAR,JJPAR)    = 0D0
  real*8   :: all_vak_no2(IIPAR,JJPAR)  = 0D0
  real*8   :: all_amf(IIPAR,JJPAR)      = 0D0
  real*8   :: all_diff_no2(IIPAR,JJPAR) = 0D0
  real*8   :: all_forcing(IIPAR,JJPAR)  = 0D0
  real*8   :: all_cost(IIPAR,JJPAR)     = 0D0
  integer  :: all_count(IIPAR,JJPAR)    = 0

  !==================================================================
  ! module subroutines
  !==================================================================

  contains

  !******************************************************************

  subroutine init_tropomi_no2_obs

  ! This subroutine is used to initialize the TROPOMI NO2 observation, i.e. to
  ! check if NO2 observation is specified by the adjoint input files.

    use tracerid_mod,   only: IDNO2
    use adj_arrays_mod, only: obs_this_species, cname, id2c
    use error_mod,      only: error_stop

!#   include "comode.h"  ! import namegas

    if(obs_this_species(id2c(IDNO2))) then
      write (6,100) IDNO2, cname(id2c(IDNO2))
      write (6,*) 'IDNO2:',IDNO2
      write (6,*) 'ID2C(IDNO2):', id2c(IDNO2)
    else
      call error_stop( 'ERROR: NO2 observation is not specified in OBSERVATION &
    MENU in input.gcadj.', 'INIT_TROPOMI_NO2_OBS(tropomi_no2_obs_mod.f90)' )
    end if

    100 format(1x,'Tracer ID: ',i4,3x,'Use TROPOMI L2 ', a6)

  end subroutine init_tropomi_no2_obs

  !******************************************************************
  
  subroutine read_tropomi_no2_obs(YYYYMMDD)

  ! This subroutine is used to read TROPOMI L2 NO2 data.
  !
  ! Arguments as input:
  ! (1) YYYYMMDD (integer): current year-month-day
  ! 
  ! Module variables as output:
  ! (1) N_NO2    (integer): number of TROPOMI NO2 observations for current day
  ! (2) tropomi_no2 (tropomi_no2_obs): TROPOMI NO2 observation

    use netcdf
    use time_mod, only: expand_date, get_nymd

    ! arguments
    integer, intent(in) :: YYYYMMDD

    ! local variables
    integer             :: i,j
    integer             :: FID, N_ID
    integer             :: lon_ID, lat_ID, time_ID
    integer             :: NO2_ID, err_ID, amf_ID
    integer             :: surf_P_ID, ak_ID
    character(len=255)  :: filename
    character(len=255)  :: dir
    character(len=4)    :: tmp
    real*4, allocatable :: tmp4(:), tmpak(:,:), tmpsp(:,:)
    real*8, allocatable :: tmp8(:)
    logical             :: lf

    ! tropomi data filename root
    filename='POMINO_GC_adj_NO2_YYYYMMDD.nc'

    ! expand date tokens in the filename
    call expand_date(filename,YYYYMMDD,9999)

    ! add file directory
    dir='./data/POMINO-NO2/'
    filename=trim(dir)//trim(filename)

    ! check whether the file exist
    inquire(file=trim(filename), exist=lf)
    if (.not. lf) then
      N_NO2=0
      write(6,120) get_nymd()
      return
    end if

    120 format(1x, '- read_tropomi_no2_obs (warning): No data file on', i10)

    write(6,100) trim(filename)
    100 format(1x, '- read_tropomi_no2_obs: reading file ', a)

    ! open the file and assign a file ID (FID) to it
    call check(nf90_open(filename, nf90_nowrite, FID),0)

    ! get data record IDs
    call check(nf90_inq_dimid(FID, 'n_obs',     N_ID),    100)
    call check(nf90_inq_varid(FID, 'lon',     lon_ID),    101)
    call check(nf90_inq_varid(FID, 'lat',     lat_ID),    102)
    call check(nf90_inq_varid(FID, 'TAI10',  time_ID),    103)
    call check(nf90_inq_varid(FID, 'no2',     NO2_ID),    104)
    call check(nf90_inq_varid(FID, 'ERR',     err_ID),    105)
    call check(nf90_inq_varid(FID, 'AK',       ak_ID),    106)
    call check(nf90_inq_varid(FID,'surf_press',surf_P_ID),107)
    call check(nf90_inq_varid(FID, 'AMF',     amf_ID),    110)

    ! read number of observations, N_NO2
    call check(nf90_inquire_dimension(FID,N_ID,tmp,N_NO2),200)

    write(6,110) N_NO2, get_nymd()
    110 format(1x,'Number of TROPOMI NO2 observations is',i10,' on',i10)

    if (N_NO2==0) then
      ! close the file
      call check(nf90_close(FID),9999)
      return
    end if

    ! read 1-D data
    ! allocate temporary array
    allocate(tmp4(N_NO2))
    allocate(tmpak(levels, N_NO2))
    allocate(tmpsp(levels, N_NO2))
    allocate(tmp8(N_NO2))
    tmp4    = 0E0
    tmpak   = 0E0
    tmpsp   = 0E0
    tmp8    = 0D0
    
    ! allocate TROPOMI NO2 observations array
    if (allocated(tropomi_no2)) deallocate(tropomi_no2)
    allocate(tropomi_no2(N_NO2))

    ! allocate the flags array
    if (allocated(flags)) deallocate(flags)
    allocate(flags(N_NO2))

    ! read longitude
    call check(nf90_get_var(FID,lon_ID,tmp4),301)
    tropomi_no2(1:N_NO2)%lon = tmp4(1:N_NO2)

    ! read latitude
    call check(nf90_get_var(FID,lat_ID,tmp4),302)
    tropomi_no2(1:N_NO2)%lat = tmp4(1:N_NO2)

    ! read time
    call check(nf90_get_var(FID,time_ID,tmp8),303)
    tropomi_no2(1:N_NO2)%time = tmp8(1:N_NO2)

    ! read tropospheric NO2 VCD
    call check(nf90_get_var(FID,NO2_ID,tmp8),304)
    tropomi_no2(1:N_NO2)%NO2 = tmp8(1:N_NO2)

    ! read NO2 VCD error
    call check(nf90_get_var(FID,err_ID,tmp8),305)
    tropomi_no2(1:N_NO2)%err = tmp8(1:N_NO2)

    ! read averaging kernel
    call check(nf90_get_var(FID,ak_ID,tmpak),306)
    do i = 1,levels
      tropomi_no2(1:N_NO2)%ak(i) = tmpak(i,1:N_NO2)
    end do

    ! read surface pressure
    call check(nf90_get_var(FID,surf_P_ID,tmpsp),307)
    do i = 1,levels
      tropomi_no2(1:N_NO2)%surf_P(i) = tmpsp(i,1:N_NO2)
    end do

    ! read air mass factor
    call check(nf90_get_var(FID,amf_ID,tmp4),310)
    tropomi_no2(1:N_NO2)%amf = tmp4(1:N_NO2)

    ! close the file
    call check(nf90_close(FID),9999)

    deallocate(tmp4)
    deallocate(tmpak)
    deallocate(tmpsp)
    deallocate(tmp8)

  end subroutine read_tropomi_no2_obs

  !******************************************************************

  subroutine calc_tropomi_no2_force(cost_func)

  ! This subroutine is used to calculate the adjoint forcing from TROPOMI L2 NO2
  ! observation and update the cost function.
  !
  ! Arguments as input/output:
  ! (1) cost_func (real, double precision): cost function, unitless

    use adj_arrays_mod, only: obs_freq, id2c
    use comode_mod,     only: jlop
    use comode_mod,     only: cspec_after_chem, cspec_after_chem_adj
    use dao_mod,        only: bxheight,airden
    use grid_mod,       only: get_ij
    use time_mod,       only: get_nymd, get_nhms
    use time_mod,       only: get_taub, get_tau, get_taue
    use tracerid_mod,   only: IDNO2
    use error_mod,      only: error_stop
    use tropopause_mod, only: its_in_the_trop
    use pressure_mod,   only: get_pcenter, get_pedge ! hPa
    use tracer_mod,     only: xnumolair

    ! arguments
    real*8              :: cost_func

    ! local variables
    integer             :: nc, nt, k
    integer             :: iijj(2), i, j, l, jloop

    real*8              :: old_cost, tmp_cost
    real*8, allocatable :: new_cost(:)
    real*8              :: diff
    real*8              :: forcing

    real*8              :: gc_no2_conc(LLPAR)  ! GC NO2 concentraction, molec/cm3
    real*8              :: gc_no2_layer(LLPAR) ! GC NO2 layer column, molec/cm2
    real*8              :: gc_no2_column       ! GC NO2 total column, molec/cm2
    ! TROPOMI NO2 VCD after applying averaging kernel and GC NO2 profile
    real*8              :: vak, amf_gc
    real*8              :: sw(levels)          ! TROPOMI scattering weight
    real*8              :: swp(levels)         ! TM5 pressure levels
    real*8              :: sw_gc(LLPAR), dp(LLPAR)       

    logical             :: first = .True.

    print *, '- calc_tropomi_no2_force: TROPOMI NO2 forcing...'

    ! initialize
    if (first) then
      call init_tropomi_no2_obs
      first = .False.
    end if

    ! save the value of the old cost function
    old_cost=cost_func

    ! check if it is the last time of a day
    if (obs_freq>60) then
      call error_stop('ERROR: 236000 - obs_freq * 100 is not valid!',&
                      'CALC_TROPOMI_NO2_FORCE(tropomi_no2_obs_mod.f90)' )
    end if

    if (get_nhms()==236000-obs_freq*100) then
      call read_tropomi_no2_obs(get_nymd())
    end if

    ! for the conditions that there are no observations for current day
    if (N_NO2==0) then
      print *, '- calc_tropomi_no2_force: No TROPOMI NO2 observations for the &
                current day'
      ! output some diagnostic data for the whole simulation period
      if (abs(get_taub()-get_tau())<1e-6) then
        call make_average_tropomi_no2
      endif
      return
    endif
 
    ! get observations for the current time window
    call get_obs

    ! for the conditions that there are no observations for current time window
    if (N_CURR==0) then
      print *, '- calc_tropomi_no2_force: No TROPOMI NO2 observations for the &
                current time window'
      ! output some diagnostic data for the whole simulation period
      if (abs(get_taub()-get_tau())<1e-6) then
        call make_average_tropomi_no2
      endif
      return
    endif

    ! reset variables for current time window
    curr_gc_no2   = 0D0
    curr_gc_s_no2 = 0D0
    curr_gc_amf   = 0D0
    curr_no2      = 0D0
    curr_s_no2    = 0D0
    curr_vak_no2  = 0D0
    curr_amf      = 0D0
    curr_diff_no2 = 0D0
    curr_forcing  = 0D0
    curr_cost     = 0D0
    curr_count    = 0

    if(allocated(new_cost)) deallocate(new_cost)
    allocate(new_cost(N_CURR))
    new_cost=0D0

    nc=0 ! index for observations in current time window

    ! loop over all observations in current day
    do nt=1,N_NO2

      ! only consider observations in current time window and simulation domain
      if (flags(nt)) then

        ! gridbox indices in GC for the observation
        iijj=get_ij(real(tropomi_no2(nt)%lon, 4), real(tropomi_no2(nt)%lat, 4))
        i=iijj(1)
        j=iijj(2)

        ! NO2 outside the troposphere is set to 0
        gc_no2_conc  = 0D0
        gc_no2_layer = 0D0
        sw_gc        = 0D0
        dp           = 0D0

        do l=1,LLPAR
          if(its_in_the_trop(i,j,l)) then
            jloop=jlop(i,j,l)
            ! units of gc_no2_conc : molec/cm3;
            !          bxheight    : m;
            !          gc_no2_layer: molec/cm2.
            gc_no2_conc(l)=cspec_after_chem(jloop,id2c(IDNO2))
            gc_no2_layer(l)=gc_no2_conc(l)*bxheight(i,j,l)*100.0
          endif
        enddo

        gc_no2_column=sum(gc_no2_layer)
        curr_count(i,j) = curr_count(i,j)+1
        ! calculate sum here, mean value will be calculated right before output
        curr_gc_no2(i,j)= curr_gc_no2(i,j)+gc_no2_column
        curr_no2(i,j)=curr_no2(i,j)+tropomi_no2(nt)%NO2

        ! calculate TM5 pressure levels and scattering weights
        do k=1, levels
          sw(k)=tropomi_no2(nt)%amf*tropomi_no2(nt)%ak(k)
          swp(k)=tropomi_no2(nt)%surf_P(k)
        enddo

        ! interpolate the TM5 scattering weights to GC vertical layers
        do l=1, LLPAR
          dp(l)=get_pedge(i,j,l)-get_pedge(i,j,l+1)
          ! unit conversion from molec/cm3 to mixing ratio (mol/mol)
          ! units of gc_no2_conc: molec/cm3;
          !          airden     : kg/m3;
          !          xnumolair  : molec/kg.
          gc_no2_conc(l)=gc_no2_conc(l)*1D6/(airden(l,i,j)*xnumolair)

          ! for GC pressures that even higher than the highest TM5 pressure
          if (get_pcenter(i,j,l)>swp(1)) then
            sw_gc(l)=sw(2)+(sw(1)-sw(2))*(get_pcenter(i,j,l)-swp(2))/(swp(1)-swp(2))
            cycle
          endif

          do k=2, levels
            if((get_pcenter(i,j,l)<swp(k-1)) .and. (get_pcenter(i,j,l)>=swp(k))) then
              sw_gc(l)=sw(k)+(sw(k-1)-sw(k))*(get_pcenter(i,j,l)-swp(k))/&
                       (swp(k-1)-swp(k))
              exit
            endif
          enddo
        enddo

        ! apply averaging kernel and GC NO2 profiles through GC air mass factor
        amf_gc=sum(gc_no2_conc*dp*sw_gc)/sum(gc_no2_conc*dp)
        vak=tropomi_no2(nt)%NO2*tropomi_no2(nt)%amf/amf_gc

        curr_gc_s_no2(i,j)=curr_gc_s_no2(i,j)+gc_no2_column*amf_gc
        curr_gc_amf(i,j)=curr_gc_amf(i,j)+amf_gc
        curr_s_no2(i,j)=curr_s_no2(i,j)+tropomi_no2(nt)%NO2*tropomi_no2(nt)%amf
        curr_vak_no2(i,j)=curr_vak_no2(i,j)+vak
        curr_amf(i,j)=curr_amf(i,j)+tropomi_no2(nt)%amf

        ! calculate adjoint forcing and cost function
        diff=gc_no2_column*amf_gc-tropomi_no2(nt)%NO2*tropomi_no2(nt)%amf
        forcing=diff/((tropomi_no2(nt)%err*tropomi_no2(nt)%amf)**2)
        tmp_cost=0.5*forcing*diff

        nc=nc+1
        curr_diff_no2(i,j)=curr_diff_no2(i,j)+diff
        curr_forcing(i,j)=curr_forcing(i,j)+forcing
        curr_cost(i,j)=curr_cost(i,j)+tmp_cost

        all_forcing(i,j)=all_forcing(i,j)+forcing
        all_cost(i,j)=all_cost(i,j)+tmp_cost
        new_cost(nc)=new_cost(nc)+tmp_cost

        ! pass the adjoint forcing back to adjoint tracer array
        do l=1, LLPAR
          if(its_in_the_trop(i,j,l)) then
            jloop=jlop(i,j,l)
            cspec_after_chem_adj(jloop,id2c(IDNO2))     = &
                cspec_after_chem_adj(jloop,id2c(IDNO2)) + &
                forcing*bxheight(i,j,l)*100.0*amf_gc
          endif
        enddo

      endif
    enddo
 
    ! update cost function
    cost_func=cost_func+sum(new_cost)
    print *,'Current COST_FUNC:', cost_func
    print *,'TROPOMI NO2 contribution to COST_FUNC:', cost_func-old_cost
    print *,'MIN/MAX cspec_after_chem_adj:', minval(cspec_after_chem_adj),&
                                             maxval(cspec_after_chem_adj)
    print *,'MIN/MAX cspec_after_chem_adj at:',minloc(cspec_after_chem_adj),&
                                               maxloc(cspec_after_chem_adj)
    print *,'MIN/MAX new_cost:', minval(new_cost), maxval(new_cost)
    print *,'MIN/MAX new_cost at:', minloc(new_cost), maxloc(new_cost)

    call make_current_tropomi_no2

    if(abs(get_taub()-get_tau())<1e-6) then
      call make_average_tropomi_no2
    endif

  end subroutine calc_tropomi_no2_force

  !******************************************************************

  subroutine make_current_tropomi_no2

  ! This subroutine is used to output some diagnostic data for current
  ! assimilation time window.

  use netcdf
  use directory_adj_mod, only: diagadj_dir
  use adj_arrays_mod,    only: n_calc
  use adj_arrays_mod,    only: expand_name
  use time_mod,          only: expand_date
  use time_mod,          only: get_nymd, get_nhms

  ! local variables
  character(len=255)     :: filename, output_file
  integer                :: i, j
  integer                :: fid, lon_dim_id, lat_dim_id
  real*4                 :: tmp_no2(IIPAR,JJPAR,11)

  filename='gctm.tropomi.no2.YYYYMMDD.hhmm.NN'
  filename=trim(filename)

  ! expand the date tokens
  call expand_date(filename,get_nymd(),get_nhms())

  ! expand the iteration number token
  call expand_name(filename,n_calc)

  output_file=trim(diagadj_dir)//trim(filename)//'.nc'

  write(6,100) trim(output_file)
  100 format(1x,'- make_current_tropomi_no2: writing ',a)

  ! open file
  call check(nf90_create(output_file,nf90_clobber,fid),0)

  ! define dimensions
  call check(nf90_def_dim(fid,'lon',IIPAR,lon_dim_id),100)
  call check(nf90_def_dim(fid,'lat',JJPAR,lat_dim_id),101)

  ! end define mode
  call check(nf90_enddef(fid),8888)

  tmp_no2=0.0

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
  do j=1, JJPAR
    do i=1, IIPAR
      if(curr_count(i,j)>0) then
        tmp_no2(i,j,1)  = real(curr_gc_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,2)  = real(curr_gc_s_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,3)  = real(curr_gc_amf(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,4)  = real(curr_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,5)  = real(curr_s_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,6)  = real(curr_vak_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,7)  = real(curr_amf(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,8)  = real(curr_diff_no2(i,j),4) / real(curr_count(i,j),4)
        tmp_no2(i,j,9)  = real(curr_forcing(i,j),4)
        tmp_no2(i,j,10) = real(curr_cost(i,j),4)
        tmp_no2(i,j,11) = real(curr_count(i,j),4)

        all_gc_no2(i,j)=all_gc_no2(i,j)+tmp_no2(i,j,1)
        all_gc_s_no2(i,j)=all_gc_s_no2(i,j)+tmp_no2(i,j,2)
        all_gc_amf(i,j)=all_gc_amf(i,j)+tmp_no2(i,j,3)
        all_no2(i,j)=all_no2(i,j)+tmp_no2(i,j,4)
        all_s_no2(i,j)=all_s_no2(i,j)+tmp_no2(i,j,5)
        all_vak_no2(i,j)=all_vak_no2(i,j)+tmp_no2(i,j,6)
        all_amf(i,j)=all_amf(i,j)+tmp_no2(i,j,7)
        all_diff_no2(i,j)=all_diff_no2(i,j)+tmp_no2(i,j,8)
        all_count(i,j)=all_count(i,j)+1
      endif
    enddo
  enddo
!!$OMP END PARALLEL DO

  ! create and save variables in netcdf file

  ! curr_gc_no2
  call write_nc_2d_float(fid,'gc_no2','GEOS-Chem vertical column NO2', &
       'molec/cm2', tmp_no2(:,:,1), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_gc_s_no2
  call write_nc_2d_float(fid,'gc_s_no2','GEOS-Chem slant column NO2', &
       'molec/cm2', tmp_no2(:,:,2), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_gc_amf
  call write_nc_2d_float(fid,'gc_amf','GEOS-Chem air mass factor', &
       'unitless', tmp_no2(:,:,3), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_no2
  call write_nc_2d_float(fid,'tropomi_no2','TROPOMI vertical column NO2', &
       'molec/cm2', tmp_no2(:,:,4), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_s_no2
  call write_nc_2d_float(fid,'tropomi_s_no2','TROPOMI slant column NO2', &
       'molec/cm2', tmp_no2(:,:,5), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_vak_no2
  call write_nc_2d_float(fid,'tropomi_vak_no2', &
       'TROPOMI vertical column NO2 corrected by GEOS-Chem air mass factor', &
       'molec/cm2', tmp_no2(:,:,6), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_amf
  call write_nc_2d_float(fid,'tropomi_amf','TROPOMI air mass factor', &
       'unitless', tmp_no2(:,:,7), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_diff
  call write_nc_2d_float(fid,'diff','gc_no2 - tropomi_vak_no2', &
       'molec/cm2', tmp_no2(:,:,8), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_forcing
  call write_nc_2d_float(fid,'forcing','adjoint forcing (diff/err^2)', &
       'cm2/molec', tmp_no2(:,:,9), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_cost
  call write_nc_2d_float(fid,'cost','cost function', &
       'unitless', tmp_no2(:,:,10), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! curr_count
  call write_nc_2d_float(fid,'count','# of observations in each gridbox', &
       'unitless', tmp_no2(:,:,11), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! close file
  call check(nf90_close(fid),9999)

  end subroutine make_current_tropomi_no2

  !******************************************************************

  subroutine make_average_tropomi_no2

  ! This subroutine is used to output some diagnostic data for the whole
  ! simulation time period.

  use netcdf
  use directory_adj_mod, only: diagadj_dir
  use adj_arrays_mod,    only: n_calc
  use adj_arrays_mod,    only: expand_name
  use time_mod,          only: expand_date
  use time_mod,          only: get_nymd

  ! local variables
  character(len=255)     :: filename, output_file
  integer                :: i, j
  integer                :: fid, lon_dim_id, lat_dim_id
  real*4                 :: tmp_no2(IIPAR,JJPAR,11)

  filename='gctm.tropomi.no2.YYYYMMDD.NN'
  filename=trim(filename)

  ! expand the date tokens
  call expand_date(filename,get_nymd(),9999)

  ! expand the iteration number token
  call expand_name(filename,n_calc)

  output_file=trim(diagadj_dir)//trim(filename)//'.nc'

  write(6,100) trim(output_file)
  100 format(1x,'- make_average_tropomi_no2: writing ',a)

  ! open file
  call check(nf90_create(output_file,nf90_clobber,fid),0)

  ! define dimensions
  call check(nf90_def_dim(fid,'lon',IIPAR,lon_dim_id),100)
  call check(nf90_def_dim(fid,'lat',JJPAR,lat_dim_id),101)

  ! end define mode
  call check(nf90_enddef(fid),8888)

  tmp_no2=0.0

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
  do j=1,JJPAR
    do i=1,IIPAR
      if(all_count(i,j)>0) then
        tmp_no2(i,j,1)=real(all_gc_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,2)=real(all_gc_s_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,3)=real(all_gc_amf(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,4)=real(all_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,5)=real(all_s_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,6)=real(all_vak_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,7)=real(all_amf(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,8)=real(all_diff_no2(i,j),4)/real(all_count(i,j),4)
        tmp_no2(i,j,9)=real(all_forcing(i,j),4)
        tmp_no2(i,j,10)=real(all_cost(i,j),4)
        tmp_no2(i,j,11)=real(all_count(i,j),4)
      endif
    enddo
  enddo
!!$OMP END PARALLEL DO

  ! create and save variables in netcdf file

  ! all_gc_no2
  call write_nc_2d_float(fid,'gc_no2','GEOS-Chem vertical column NO2', &
       'molec/cm2', tmp_no2(:,:,1), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_gc_s_no2
  call write_nc_2d_float(fid,'gc_s_no2','GEOS-Chem slant column NO2', &
       'molec/cm2', tmp_no2(:,:,2), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_gc_amf
  call write_nc_2d_float(fid,'gc_amf','GEOS-Chem air mass factor', &
       'unitless', tmp_no2(:,:,3), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_no2
  call write_nc_2d_float(fid,'tropomi_no2','TROPOMI vertical column NO2', &
       'molec/cm2', tmp_no2(:,:,4), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_s_no2
  call write_nc_2d_float(fid,'tropomi_s_no2','TROPOMI slant column NO2', &
       'molec/cm2', tmp_no2(:,:,5), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_vak_no2
  call write_nc_2d_float(fid,'tropomi_vak_no2', &
       'TROPOMI vertical column NO2 corrected by GEOS-Chem air mass factor', &
       'molec/cm2', tmp_no2(:,:,6), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_amf
  call write_nc_2d_float(fid,'tropomi_amf','TROPOMI air mass factor', &
       'unitless', tmp_no2(:,:,7), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_diff
  call write_nc_2d_float(fid,'diff','gc_no2 - tropomi_vak_no2', &
       'molec/cm2', tmp_no2(:,:,8), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_forcing
  call write_nc_2d_float(fid,'forcing','adjoint forcing (diff/err^2)', &
       'cm2/molec', tmp_no2(:,:,9), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_cost
  call write_nc_2d_float(fid,'cost','cost function', &
       'unitless', tmp_no2(:,:,10), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! all_count
  call write_nc_2d_float(fid,'count','# of observations in each gridbox', &
       'unitless', tmp_no2(:,:,11), IIPAR,JJPAR, lon_dim_id, lat_dim_id)

  ! close file
  call check(nf90_close(fid),9999)

  end subroutine make_average_tropomi_no2

  !******************************************************************

  subroutine write_nc_2d_float(fid,names,longnames,units,var,iin,jjn,dim1,dim2)

  ! This subroutine is used to create and save 2D variable in the netcdf file.
  !
  ! Auguments as input:
  ! (1) fid         (integer): file ID of the netcdf file we are writing
  ! (2) names     (character): name of the variable
  ! (3) longnames (character): long name of the variable
  ! (4) units     (character): unit of the variable
  ! (5) var    (real, kind=4): data array of the variable
  ! (6) iin         (integer): length of 1st dimension of the variable
  ! (7) jjn         (integer): length of 2nd dimension of the variable
  ! (8) dim1        (integer): dimension ID 1 of the variable
  ! (9) dim2        (integer): dimension ID 2 of the variable

  use netcdf

  integer,      intent(in) :: fid
  integer,      intent(in) :: iin, jjn
  integer,      intent(in) :: dim1, dim2
  character(*), intent(in) :: names, longnames, units
  real*4,       intent(in) :: var(iin,jjn)

  integer                  :: vid

  ! open define mode
  call check(nf90_redef(fid),9000)

  ! define the variable
  call check(nf90_def_var(fid,trim(names),nf90_float,(/dim1,dim2/),vid),9001)

  ! define attributes
  call check(nf90_put_att(fid,vid,'longname',trim(longnames)),9002)
  call check(nf90_put_att(fid,vid,'unit',trim(units)),9003)

  ! end define mode
  call check(nf90_enddef(fid),9004)

  ! save variable
  call check(nf90_put_var(fid,vid,var),9005)

  end subroutine write_nc_2d_float

  !******************************************************************

  subroutine get_obs

  ! This subroutine is used to find all the TROPOMI NO2 observations in the
  ! current time window and simulation domain.
  !
  ! Module variables as input:
  ! (1) N_NO2 (integer)  : number of observations in the current day
  !
  ! Module variables as output:
  ! (1) flags (logical)  : whether or not a specific observation is in current
  !                        time window and simulation domain
  ! (2) N_CURR (integer) : number of observations in the current time window

    use grid_mod, only: get_xedge, get_yedge
    use time_mod, only: get_jd, get_tau, get_nymd, get_nhms

    ! local variables
    real*8               :: half_time_window
    real*8               :: window_begin, window_end
    real*8               :: jd85, jd10, jd10_85
    real*8               :: curr_tau
    integer              :: nt

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
    real*8, save         :: xedge_min, xedge_max
    real*8, save         :: yedge_min, yedge_max
#endif

    ! get the difference between jd10 and jd85 (jd: Astronomical Julian Date)
    ! In GEOS-Chem, reference time is 1/1/1985; while in TROPOMI, it is 1/1/2010
    jd85    = get_jd(19850000,000000)
    jd10    = get_jd(20100000,000000)
    jd10_85 = (jd10-jd85)*24D0 ! days to hours
    !print *,'jd10_85 in hrs', jd10_85

    ! get current TAU from GEOS-Chem
    curr_tau = get_tau()
    !print *, 'current TAU_85(hr)', curr_tau

    ! transfer current GEOS-Chem TAU into TROPOMI TAU
    curr_tau = curr_tau-jd10_85
    !print *, 'current TAU_10(hr)', curr_tau

    ! transfer TAU unit from hours to seconds
    curr_tau = curr_tau*3600D0
    !print *, 'current TAU_10(s)', curr_tau

    ! get half time window
    half_time_window = time_window/2D0
    half_time_window = half_time_window*60D0 ! minutes to seconds
    !print *, 'half time window', half_time_window

    ! get current time window
    window_begin = curr_tau-half_time_window
    window_end   = curr_tau+half_time_window
    !print *, 'window begin', window_begin
    !print *, 'window end', window_end
    !print *, 'TROPOMI first 10 time', tropomi_no2(1:10)%time

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
    xedge_min = get_xedge(1)
    xedge_max = get_xedge(IIPAR+1)
    yedge_min = get_yedge(1)
    yedge_max = get_yedge(JJPAR+1)
    print *, 'Nested domain'
    print *, 'XEDGE_MIN:', xedge_min
    print *, 'XEDGE_MAX:', xedge_max
    print *, 'YEDGE_MIN:', yedge_min
    print *, 'YEDGE_MAX:', yedge_max
#endif

    N_CURR=0

    ! find observations in current time window and simulation domain
    do nt=1,N_NO2

      if (  (tropomi_no2(nt)%time>=window_begin)&
      .and. (tropomi_no2(nt)%time<window_end)&
#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
      .and. (tropomi_no2(nt)%lon>=xedge_min)&
      .and. (tropomi_no2(nt)%lon<=xedge_max)&
      .and. (tropomi_no2(nt)%lat>=yedge_min)&
      .and. (tropomi_no2(nt)%lat<=yedge_max)&
#endif
         ) then

        flags(nt)=.true.
        N_CURR=N_CURR+1

      else
        flags(nt)=.false.
      endif
    enddo

    write(6,100) N_CURR, get_nhms()
    100 format(1x,'Number of TROPOMI NO2 observations is',i10,' at',i10.6)
    
  end subroutine get_obs

  !******************************************************************

  subroutine check(status, location)

  ! This subroutine is used to check the status of calls to netCDF libraries.
  ! (dkh, 02/15/09)
  !
  ! Arguments as input:
  ! (1) status (integer): completion status of netCDF library call
  ! (2) location (integer): location where netCDF library call was made

    use error_mod, only: error_stop
    use netcdf

    ! arguments
    integer, intent(in) :: status
    integer, intent(in) :: location

    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      print *, 'At location = ', location
      call error_stop('netCDF error','tropomi_no2_obs_mod')
    end if

  end subroutine check

end module TROPOMI_NO2_OBS_MOD
