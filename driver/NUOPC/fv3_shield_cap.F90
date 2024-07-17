module fv3_shield_cap

  !-----------------------------------------------------------------------------
  ! Basic NUOPC Model cap
  ! Interfaces the FV3-SHiELD coupled atmos-land-oml model with other
  ! NUOPC compliant model components (e.g. WW3 waves and MOM6 ocean)
  !
  ! Following template:
  ! https://github.com/esmf-org/nuopc-app-prototypes/tree/main/AtmOcnMedProto/
  ! Reference fv3 cap for ufs:
  ! https://github.com/NOAA-EMC/fv3atm/blob/9743346431c46642958712690e2c2733763ce5de/fv3_cap.F90
  ! Reference for input files:
  ! https://ufs-weather-model.readthedocs.io/en/ufs-v1.0.0/InputsOutputs.html
  ! Reference MOM6 cap for interfacing with GFDL's FMS:
  ! https://ncar.github.io/MOM6/APIs/mom__cap_8F90_source.html
  !
  !-----------------------------------------------------------------------------
  ! Authors:
  !
  ! (1) Stephen G. Penny (Aug-2022)
  !     steve.penny@sofarocean.com
  !
  !-----------------------------------------------------------------------------


  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS    => SetServices
    
  !-----------------------------------------------------------------------------
  ! add use statements for your model's initialization
  ! and run subroutines
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Used by: 
  ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/main/simple_coupler/coupler_main.F90
  ! in docker container: fv3_gfsphysics/simple_coupler/coupler_main.F90
  ! program coupler_main
  !-----------------------------------------------------------------------
  !
  !   program that couples component models for the atmosphere,
  !   ocean (amip), land, and sea-ice using the exchange module. 
  !
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Insert use modules and leading instantiations from coupler_main
  !-----------------------------------------------------------------------
  ! 11/20/23:
  ! Updating with:
  ! https://github.com/NOAA-GFDL/FMScoupler/blob/2023.03/SHiELD/coupler_main.F90
  ! -and- later missed updates here:
  ! https://github.com/NOAA-GFDL/FMScoupler/pull/87/files

  ! https://github.com/NOAA-GFDL/FMS/blob/2023.02/libFMS.F90
  ! NOTE (11/20/23): starting with FMS 2023.02, shield team removed all explicit callouts to FMS modules and replaced with this:
  use FMS
  use FMSconstants,      only: fmsconstants_init

  use  atmos_model_mod,  only: atmos_model_init, atmos_model_end,  &
                               update_atmos_model_dynamics,        &
                               update_atmos_radiation_physics,     &
                               update_atmos_model_state,           &
                               atmos_data_type, atmos_model_restart

  !--- FMS old io
  #ifdef use_deprecated_io
  use fms_io_mod, only: fms_io_exit!< This can't be removed until fms_io is not used at all
  #endif

  ! After upgrading to FMS 2023.02, the rest here can be removed:
! https://github.com/NOAA-GFDL/FMS/blob/7898b51d901c9c804d56b5c8d76bb0169c1c65fc/libFMS.F90#L768
  use time_manager_mod,  only: FmsTime_type => time_type, & 
                               fms_time_manager_set_calendar_type => set_calendar_type, &
                               fms_time_manager_set_time => set_time,    &
                               fms_time_manager_set_date => set_date, &
                               fms_time_manager_days_in_month => days_in_month, &
                               fms_time_manager_month_name => month_name, &
                               fms_time_manager_date_to_string => date_to_string, &
                               fms_time_manager_get_date => get_date, &
                               operator(+), operator (<), operator (>),   &
                               operator (/=), operator (/), operator (==),&
                               operator (*), THIRTY_DAY_MONTHS, JULIAN,   &
                               GREGORIAN, NOLEAP, NO_CALENDAR

! Other subroutines used by:
! https://github.com/NOAA-EMC/fv3atm/blob/86ba9011c891d11665dfafbf0e37c973f835f587/module_fcst_grid_comp.F90#LL31C3-L39C80
!                              get_atmos_model_ungridded_dim,             &
!                              atmos_model_exchange_phase_1,              &
!                              atmos_model_exchange_phase_2,              &
!                              addLsmask2grid, atmos_model_get_nth_domain_info

  use mpp_mod,           only: fms_mpp_input_nml_file => input_nml_file, &
                               fms_mpp_stdlog => stdlog
  use fms_affinity_mod,  only: fms_affinity_init, &
                               fms_affinity_set

  use fms_mod,           only: fms_init, fms_end,   &
                               fms_check_nml_error => check_nml_error, &
                               fms_error_mesg => error_mesg, &
                               fms_write_version_number => write_version_number, &
                               fms_mpp_uppercase => uppercase
  use fms2_io_mod,       only: fms2_io_ascii_read =>   ascii_read, &
                               fms2_io_file_exists => file_exists
  use mpp_mod,           only: fms_mpp_init => mpp_init, &
                               fms_mpp_pe => mpp_pe, &
                               fms_mpp_root_pe => mpp_root_pe, &
                               fms_mpp_npes => mpp_npes, & 
                               fms_mpp_get_current_pelist => mpp_get_current_pelist, &
                               fms_mpp_set_current_pelist => mpp_set_current_pelist, &
                               fms_mpp_error => mpp_error, &
                               NOTE, FATAL, WARNING
  use mpp_mod,           only: fms_mpp_clock_id => mpp_clock_id, &
                               fms_mpp_clock_begin => mpp_clock_begin, &
                               fms_mpp_clock_end => mpp_clock_end, &
                               fms_mpp_sync => mpp_sync

  use mpp_domains_mod,   only: fms_mpp_domains_get_global_domain => mpp_get_global_domain, &
                               fms_mpp_domains_global_field => mpp_global_field, &
                               CORNER
  use memutils_mod,      only: fms_memutils_print_memuse_stats => print_memuse_stats
  use sat_vapor_pres_mod,only: fms_sat_vapor_pres_init => sat_vapor_pres_init

  use diag_manager_mod,  only: fms_diag_init => diag_manager_init, &
                               fms_diag_end => diag_manager_end, &
                               fms_diag_get_base_date => get_base_date, &
                               fms_diag_set_time_end => diag_manager_set_time_end

  use data_override_mod, only: fms_data_override_init => data_override_init

  implicit none

  character(len=128) :: version = 'fv3_shield_cap.F90 - 2023/11/20 steve.penny@sofarocean.com'
  character(len=128) :: tag = 'FMSCoupler_SHiELD'


  !---- model defined-types ----
  type(atmos_data_type), save :: Atm

  ! ----- coupled model time -----
  type (FmsTime_type) :: Time_atmos, Time_init, Time_end,  &
                      Time_step_atmos, Time_step_ocean, &
                      Time_restart, Time_step_restart,  &
                      Time_start_restart, Time_restart_aux, &
                      Time_step_restart_aux, Time_start_restart_aux, &
                      Time_duration_restart_aux, Time_restart_end_aux

  integer :: num_cpld_calls, num_atmos_calls, nc, na, ret

  ! ----- coupled model initial date -----
  integer :: date_init(6)
  integer :: calendar_type = -99

  ! ----- timing flags -----
  integer :: initClock, mainClock, termClock
  integer, parameter :: timing_level = 1

  ! ----- namelist -----
  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /) !< The date that the current integration starts with
  character(len=17) :: calendar = '                 '   !< The calendar type used by the current integration.  Valid values are
                                                        !! consistent with the time_manager module: 'gregorian', 'julian',
                                                        !! 'noleap', or 'thirty_day'. The value 'no_calendar' cannot be used
                                                        !! because the time_manager's date !! functions are used.
                                                        !! All values must be lower case.
  logical :: force_date_from_namelist = .false.   !< Flag that determines whether the namelist variable current_date should override
                                                  !! the date in the restart file `INPUT/coupler.res`.  If the restart file does not
                                                  !! exist then force_date_from_namelist has no effect, the value of current_date
                                                  !! will be used.
  integer :: years=0                              !< Number of years the current integration will be run
  integer :: months=0                             !< Number of months the current integration will be run
  integer :: days=0                               !< Number of days the current integration will be run
  integer :: hours=0                              !< Number of hours the current integration will be run
  integer :: minutes=0                            !< Number of minutes the current integration will be run
  integer :: seconds=0                            !< Number of seconds the current integration will be run
  integer :: dt_atmos = 0                         !< Atmospheric model time step in seconds
  integer :: dt_ocean = 0                         !< Ocean model time step in seconds - NOT USED IN THIS MODEL
  integer :: restart_days = 0                     !< Time interval in days to write out intermediate restart files
  integer :: restart_secs = 0                     !< Time interval in seconds to write out intermediate restart files
  integer :: restart_start_days = 0               !< Start time in days to write out intermediate restart files
  integer :: restart_start_secs = 0               !< Start time in seconds to write out intermediate restart files
  integer :: restart_days_aux = 0                 !< Time interval in days for auxiliary restart files
  integer :: restart_secs_aux = 0                 !< Time interval in seconds for auxiliary restart files
  integer :: restart_start_days_aux = 0           !< Start time in days for auxiliary restart files
  integer :: restart_start_secs_aux = 0           !< Start time in days for auxiliary restart files
  integer :: restart_duration_days_aux = 0        !< Duration in days for auxiliary restart files
  integer :: restart_duration_secs_aux = 0        !< Duration in seconds for auxiliary restart files
  integer :: atmos_nthreads = 1                   !< Number of OpenMP threads to use in the atmosphere
  logical :: use_hyper_thread = .false.           !< If .TRUE>, affinity placement (if activated) will consider virtual cores
                                                  !! in the placement algorithm
  integer :: iau_offset = 0

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, &
                          years, months, days, hours, minutes, seconds, &
                          iau_offset, dt_atmos, dt_ocean, atmos_nthreads, &
                          use_hyper_thread, restart_secs, restart_days, &
                          restart_start_secs, restart_start_days, &
                          restart_secs_aux, restart_days_aux, &
                          restart_start_secs_aux, restart_start_days_aux, &
                          restart_duration_secs_aux, restart_duration_days_aux

  ! ----- local variables -----
  character(len=32) :: timestamp
  logical :: intrm_rst, intrm_rst_1step
  character(ESMF_MAXSTR) :: msgString, input_string, output_string
  
  ! End insert
  !-----------------------------------------------------------------------------
  
  private
  
  ! For debugging at runtime  (defaults):
  logical :: use_mlm = .false.
  logical :: use_gridcreate_addedges = .false.
  logical :: use_sfcpropwinds = .false.  ! Whether to use winds from intdiag data structure (default) or Sfcprop data structure

  ! For running with atmosphere-ocean fluxes
  logical :: use_aofluxes = .false.

  ! For setting up internal grid object
  logical :: use_mosaic=.false.
  character(ESMF_MAXSTR) :: mosaic_filename, mosaic_tileFilePath
  integer :: layout_x=-1, layout_y=-1
  integer :: tilesize=-1   ! The number of elements on a dimension of each tile (e.g. C96 -> 96)

  ! Provides ability to change coupling variables at runtime:
  character(32), dimension(:), allocatable :: Sa_import_names
  character(32), dimension(:), allocatable :: Sa_export_names

  integer :: Sa_import_len=0, Sa_export_len=0

  logical :: dodebug_Advertise = .false.
  logical :: dodebug_Realize = .false.
  logical :: dodebug_DataInitialize = .false.
  logical :: dodebug_Advance = .false.
  logical :: dodebug_Finalize = .false.
  logical :: dodebug_bypassImport = .false.
  logical :: dodebug_bypassExport = .false.
  logical :: dodebug_getImport = .false.
  logical :: dodebug_getImport_writenetcdf = .false.
  logical :: dodebug_setExport = .false.
  logical :: dodebug_setExport_writenetcdf = .false.
  logical :: dodebug_setExport_usesrf = .false.
  logical :: dodebug_setExport_usecoupling = .false.
  logical :: dodebug_setExport_useneutral = .false.
  logical :: dodebug_setExport_multistep = .false.
  logical :: dodebug_setExport_datainit_uv = .false.
  logical :: dodebug_couplerinit =.false.
  logical :: dodebug_gridcreate = .false.

  real(ESMF_KIND_R8), dimension(:,:), pointer :: farray_u, farray_v, farray_t2m, &  ! Atmosphere exports
                                                 farray_un, farray_vn,           &
                                                 farray_rhoa, farray_tsfc,       &
                                                 farray_slmsk, farray_oceanfrac, &  ! Atmosphere masks
                                                 farray_landfrac,                &
                                                 farray_lakefrac,                & 
                                                 farray_fice,                    & 
                                                 farray_psurf,                   &  ! 
                                                 farray_q2m, farray_tisfc,       &
                                                 farray_uustar,                  &
                                                 farray_ffmm, farray_f10m,       &
                                                 farray_hpbl,                    &
                                                 farray_tml, farray_mld,         &  ! Atmosphere's mixed layer model (MLM)
                                                 farray_huml, farray_hvml,       &  ! See: https://github.com/NOAA-GFDL/SHiELD_physics/blob/69d13c245348264f1bf12fd871261d5e25e36695/GFS_layer/GFS_typedefs.F90#L191
                                                 farray_ts_som,                  & 
                                                 farray_c, farray_z,             &  ! Atmosphere imports from wave model
                                                 farray_sst,                     &  ! Atmosphere imports from ocean model
                                                 farray_ssu, farray_ssv,         &
                                                 farray_astdiff,                 &  ! Post-processed air-sea temperature difference on atmosphere grid
                                                 farray_taux,                    &  ! Atmosphere-ocean fluxes computed in atmosphere model: (Added 7/16/24)
                                                 farray_tauy,                    &  ! See supported MOM6 imports here:
                                                 farray_rain,                    &  ! https://github.com/NOAA-GFDL/MOM6/blob/2c1a9d32fce72828b7091e6f623c0ec20069e637/config_src/drivers/nuopc_cap/mom_cap_methods.F90#L105
                                                 farray_lwnet,                   &
                                                 farray_sen,                     &
                                                 farray_evap,                    &
                                                 farray_swndr,                   &
                                                 farray_swndf,                   &
                                                 farray_swvdr,                   &
                                                 farray_swvdf

  
! ! For using CMEPS with cpl_scalars:  (Added 7/1/24)
  ! set from config
  integer                   :: verbosity = 0
  character(len=80), public :: flds_scalar_name = ''
  integer, public           :: flds_scalar_num = 0
  integer, public           :: flds_scalar_index_nx = 0
  integer, public           :: flds_scalar_index_ny = 0
  integer, public           :: flds_scalar_index_ntile = 0
  logical                   :: isPresent, isSet
  character(ESMF_MAXSTR)    :: cvalue
  
  public :: SetServices
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine SetServices(model, rc)
  !-----------------------------------------------------------------------------
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config)           :: config
    
    rc = ESMF_SUCCESS
    
    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    
    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, specRoutine=Advertise, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, specRoutine=Realize, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_CompSpecialize(model, specLabel=Label_DataInitialize, specRoutine=DataInitialize, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! This adds checking of the import fields to make sure they are appropriately tagged
    ! with the correct current time as understood by the component model
    call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, specRoutine=CheckImport, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_CompSpecialize(model, specLabel=label_Advance, specRoutine=Advance, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, specRoutine=Finalize, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

  end subroutine SetServices

  
  !-----------------------------------------------------------------------------
  subroutine Advertise(model, rc) 
  ! See for reference:
  ! https://ncar.github.io/MOM6/APIs/mom__cap_8F90_source.html
  !-----------------------------------------------------------------------------

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State) :: importState, exportState

    ! local variables (for mpi with fms)
    type(esmf_vm)               :: vm
    integer                     :: mpi_comm_esmf, mpi_comm_fv3, localPet
    integer                     :: ierr

    type(ESMF_Clock)            :: clock, dclock, mclock
    type(ESMF_Config)           :: config

    ! Test:
    character(128) :: pelist_name
    integer, allocatable, dimension(:) :: pelist
    integer :: pelist_commID
    logical :: dodebug

    integer :: i,j
    character(len=*),parameter         :: subname='(fv3_shield_cap:Advertise)'
    
    rc = ESMF_SUCCESS 
    
    !--------------------------------------------------------------------------
    !STEVE: add runtime debugging
    ! query Driver for localPet
    call ESMF_GridCompGet(model, config=config, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Use the fv3-shield's active coupled mixed layer model
    call support_get_logical(model,input_string='use_mlm',input_var=use_mlm,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Access atmoshere model fluxes for export to CMEPS mediator
    call support_get_logical(model,input_string='use_aofluxes',input_var=use_aofluxes,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Add edges to the grid definition (NOT SUPPORTED, as of ESMF 8.5.0)
    call support_get_logical(model,input_string='use_gridcreate_addedges',input_var=use_gridcreate_addedges,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Whether to use winds from intdiag data structure (default) or Sfcprop data structure
    call support_get_logical(model,input_string='use_sfcpropwinds',input_var=use_sfcpropwinds,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Whether to read in grid mosiac or generate internally:
    call support_get_logical(model,input_string='use_mosaic',input_var=use_mosaic,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_string(model,input_string='mosaic_filename',input_var=mosaic_filename,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_string(model,input_string='mosaic_tileFilePath',input_var=mosaic_tileFilePath,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='layout_x',input_var=layout_x,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='layout_y',input_var=layout_y,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='tilesize',input_var=tilesize,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ---
    ! Get list of import variables. default="Sw_z0rlen,Sw_charno"
    ! ---
    call support_get_string(model,input_string='Sa_import_names',input_var=input_string,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Handle empty case
    if (trim(input_string).eq.'none' .or. trim(input_string).eq.'' .or. trim(input_string).eq.'n/a') then
        input_string = ''
    endif

    call string_split(trim(input_string),Sa_import_names)

    if (allocated(Sa_import_names)) then
        Sa_import_len = size(Sa_import_names)

        write (msgString,*) "ATM::Advertise:: Sa_import_len = ", Sa_import_len
        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        msgString="ATM::Advertise:: Sa_import_names = "
        do i = 1,Sa_import_len
            msgString=trim(msgString)//' '//trim(Sa_import_names(i))
        enddo
        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif
    ! ---

    ! 
    call support_get_logical(model,input_string='dodebug_Advertise',input_var=dodebug_Advertise,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_Realize',input_var=dodebug_Realize,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_DataInitialize',input_var=dodebug_DataInitialize,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_Advance',input_var=dodebug_Advance,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_Finalize',input_var=dodebug_Finalize,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_bypassImport',input_var=dodebug_bypassImport,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_bypassExport',input_var=dodebug_bypassExport,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_getImport',input_var=dodebug_getImport,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_getImport_writenetcdf',input_var=dodebug_getImport_writenetcdf,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport',input_var=dodebug_setExport,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_writenetcdf',input_var=dodebug_setExport_writenetcdf,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_usesrf',input_var=dodebug_setExport_usesrf,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_usecoupling',input_var=dodebug_setExport_usecoupling,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_useneutral',input_var=dodebug_setExport_useneutral,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_multistep',input_var=dodebug_setExport_multistep,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_setExport_datainit_uv',input_var=dodebug_setExport_datainit_uv,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_couplerinit',input_var=dodebug_couplerinit,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_logical(model,input_string='dodebug_gridcreate',input_var=dodebug_gridcreate,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return


    ! set cpl_scalars from config. Default to null values for standalone
    ! This was added 7/1/24 to support CMEPS history files, following:
    ! https://github.com/NOAA-EMC/fv3atm/pull/794
    ! "This will enable FV3 to sen[d] to CMEPS the 3 dimensions of either a CSG [cubed sphere grid] or Regional domain. 
    !  For the CSG, the indices of the scalar field will contain the nx,ny and number of tiles."
    ! Code from ufs:
    ! https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/fv3_cap.F90#L278C5-L323C10

    call NUOPC_CompAttributeGet(model, name="Verbosity", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    if (isPresent .and. isSet) then
        verbosity = ESMF_UtilString2Int( cvalue, &
                                         specialStringList=(/'none', 'low', 'high','max'/), &
                                         specialValueList =(/0, 1,  99,  255/), rc=rc )
        call ESMF_LogWrite(trim(subname)//' verbosity = '//trim(cvalue), ESMF_LOGMSG_INFO)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    call support_get_string(model,input_string='ScalarFieldName',input_var=flds_scalar_name,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='ScalarFieldCount',input_var=flds_scalar_num,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='ScalarFieldIdxGridNX',input_var=flds_scalar_index_nx,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='ScalarFieldIdxGridNY',input_var=flds_scalar_index_ny,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_get_integer(model,input_string='ScalarFieldIdxGridNTile',input_var=flds_scalar_index_ntile,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !--------------------------------------------------------------------------

    ! Get the esmf/nuopc mpi information to pass to fms
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call ESMF_vmget(vm, localPet=localPet, mpicommunicator=mpi_comm_esmf, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
!5/24/23! The returned MPI communicator spans the same MPI processes that the VM is defined on.
    
    call MPI_Comm_dup(mpi_comm_esmf, mpi_comm_fv3, ierr)
!5/24/23! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user code.
    ! https://mpitutorial.com/tutorials/introduction-to-groups-and-communicators/
    
    !-------------------------------------------------
    ! FV3-SHiELD model initialization routines:
    !-------------------------------------------------


    ! TRACE start: fms_init ----------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::fms_init", rc=rc)

    call fms_init(mpi_comm_fv3)

    initClock = fms_mpp_clock_id( '-Initialization' )
    call fms_mpp_clock_begin (initClock) !nesting problem

    call fms_sat_vapor_pres_init()
    call fmsconstants_init()

    call ESMF_TraceRegionExit("fv3_shield_cap::fms_init", rc=rc)
    ! TRACE end: fms_init ------------------------------------------------------


    ! Test:
!   write (output_string, "(I0)") mpi_comm_esmf
!   call ESMF_LogWrite("fv3_shield_cap:: mpi_comm_esmf = "//output_string, ESMF_LOGMSG_INFO)
!   write (output_string, "(I0)") mpi_comm_fv3
!   call ESMF_LogWrite("fv3_shield_cap:: mpi_comm_fv3 = "//output_string, ESMF_LOGMSG_INFO)
!   allocate( pelist(fms_mpp_npes()) )
!   call fms_mpp_get_current_pelist(pelist, name=pelist_name, commID=pelist_commID)
!   write (output_string, "(I0)") pelist_commID
!   call ESMF_LogWrite("fv3_shield_cap:: pelist_commID = "//output_string, ESMF_LOGMSG_INFO)


    ! Get driver clock from model to pass to coupler_init
    call NUOPC_ModelGet(model, driverClock=dclock, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! TRACE start: coupler_init ------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init", rc=rc)
    
    if (dodebug_Advertise) print *, "fv3_shield_cap:: calling coupler_init..."
    call coupler_init(model=model,clock=dclock,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call fms_memutils_print_memuse_stats('after coupler init')

    call fms_mpp_set_current_pelist()
    
    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init", rc=rc)
    ! TRACE end: coupler_init --------------------------------------------------


    !-------------------------------------------------
    ! ESMF-NUOPC commands:
    ! Advertise the model's import and export fields
    ! note: this is simply a way for models to 
    !       communicate the standard names of fields
    !-------------------------------------------------

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !----------
    ! Import
    ! see: https://github.com/NOAA-EMC/fv3gfs/blob/main/sorc/fv3gfs.fd/NEMS/src/module_EARTH_GRID_COMP.F90
    !----------

    do i=1,Sa_import_len
        if (trim(Sa_import_names(i)).eq.'Sw_charno') then
            call advertise_support(importState,standard_name="Sw_charno",LongName="charnock parameter",ShortName="charno",rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        elseif (trim(Sa_import_names(i)).eq.'Sw_z0rlen') then
            ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L272
            call advertise_support(importState,standard_name="Sw_z0rlen",LongName="roughness length",ShortName="z0rlen",rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        elseif (trim(Sa_import_names(i)).eq.'So_sst') then
            call advertise_support(importState,standard_name="So_sst",LongName="sea surface foundation temperature",ShortName="sst",rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
!       elseif (trim(Sa_import_names(i)).eq.'So_ssu') then
!           call advertise_support(importState,standard_name="So_ssu",LongName="sea surface zonal current",ShortName="ucur",rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return
!       elseif (trim(Sa_import_names(i)).eq.'So_ssv') then
!           call advertise_support(importState,standard_name="So_ssv",LongName="sea surface meridional current",ShortName="vcur",rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return
        else
            output_string = "Sa_import_names(i) = "//trim(Sa_import_names(i))//" is not a valid import variable. Choose from: Sw_z0rlen, Sw_charno, So_sst" !, So_ssu, So_ssv"
            call ESMF_LogWrite("fv3_shield_cap::Advertise:: "//trim(output_string), ESMF_LOGMSG_INFO)
        endif
    enddo


    !----------
    ! Export
    ! NOTE: These are easy to add, and don't need to be included in the mediator.
    !       Allocation of data storage/pointers occurs next in 'Realize' subroutine
    !
    !----------

    call advertise_support(exportState,standard_name="Sa_u10m",LongName="zonal wind at 10m height",ShortName="u10m",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call advertise_support(exportState,standard_name="Sa_v10m",LongName="meridional wind at 10m height",ShortName="v10m",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call advertise_support(exportState,standard_name="Sa_u10n",LongName="zonal neutral wind at 10m height",ShortName="u10n",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call advertise_support(exportState,standard_name="Sa_v10n",LongName="meridional neutral wind at 10m height",ShortName="v10n",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L313
    call advertise_support(exportState,standard_name="Sa_t2m",LongName="atmosphere temperature at 2m height",ShortName="t2m",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call advertise_support(exportState,standard_name="Sa_rhoa",LongName="surface atmospheric density",ShortName="rhoa",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L269
    call advertise_support(exportState,standard_name="Sa_tsfc",LongName="surface temperature seen by atmosphere",ShortName="tsfc",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Add post-processing fields
    call advertise_support(exportState,standard_name="Sa_astdiff",LongName="air-sea temperature difference",ShortName="astdiff",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call advertise_support(exportState,standard_name="Sa_oceanfrac",LongName="ocean fraction [0:1]",ShortName="oceanfrac",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L267
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/GFS_layer/GFS_typedefs.F90#L112
    call advertise_support(exportState,standard_name="Sa_psurf",LongName="surface pressure (Pa)",ShortName="psurf",rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Scalar data used by CMEPS describing grid
    if (flds_scalar_num > 0) then
        call NUOPC_Advertise(exportState, trim(flds_scalar_name), trim(flds_scalar_name), rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! --------------------------------------------------------------------------------
    ! The remaining fields are only used if dodebug_setExport_writenetcdf=.true.
    ! --------------------------------------------------------------------------------
    if (dodebug_setExport_writenetcdf) then

!       ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L268
        call advertise_support(exportState,standard_name="Sa_slmsk",LongName="sea/land mask array (sea:0,land:1,seaice:2)",ShortName="slmsk",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_landfrac",LongName="land fraction [0:1]",ShortName="landfrac",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_lakefrac",LongName="lake fraction [0:1]",ShortName="lakefrac",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/69d13c245348264f1bf12fd871261d5e25e36695/GFS_layer/GFS_typedefs.F90#L206
        call advertise_support(exportState,standard_name="Sa_fice",LongName="ice fraction over open water grid",ShortName="fice",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L314
        call advertise_support(exportState,standard_name="Sa_q2m",LongName="atmosphere humidity at 2m height",ShortName="q2m",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/69d13c245348264f1bf12fd871261d5e25e36695/GFS_layer/GFS_typedefs.F90#L200C66-L200C103
        call advertise_support(exportState,standard_name="Sa_tisfc",LongName="surface temperature over ice fraction",ShortName="tisfc",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L290
        call advertise_support(exportState,standard_name="Sa_uustar",LongName="friction velocity",ShortName="uustar",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L296
        call advertise_support(exportState,standard_name="Sa_ffmm",LongName="neutral wind computation component 1",ShortName="ffmm",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

!       ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L298
        call advertise_support(exportState,standard_name="Sa_f10m",LongName="neutral wind computation component 3",ShortName="f10m",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/69d13c245348264f1bf12fd871261d5e25e36695/GFS_layer/GFS_typedefs.F90#L1231C39-L1231C43
        call advertise_support(exportState,standard_name="Sa_hpbl",LongName="pbl height (m)",ShortName="hpbl",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    endif

    ! If using fv3-shield mixed layer model:
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/69d13c245348264f1bf12fd871261d5e25e36695/GFS_layer/GFS_typedefs.F90#L191
    if (use_mlm) then
        call advertise_support(exportState,standard_name="Sa_ts_som",LongName="predicted SST in SOM or MLM",ShortName="ts_som",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_mld",LongName="ocean mixed layer depth (MLD)",ShortName="mld",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_tml",LongName="ocean mixed layer temp",ShortName="tml",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_huml",LongName="ocean zonal current * MLD",ShortName="huml",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Sa_hvml",LongName="ocean meridional current * MLD",ShortName="hvml",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! Add FLUXES - 7/16/24
    ! 
    if (use_aofluxes) then
        call advertise_support(exportState,standard_name="Faxa_taux",LongName="inst_zonal_moment_flx",ShortName="taux",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_tauy",LongName="inst_merid_moment_flx_atm",ShortName="tauy",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_rain",LongName="inst_prec_rate",ShortName="prate",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_lwnet",LongName="inst_net_lw_flx",ShortName="lwnet",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_sen",LongName="inst_sensi_heat_flx",ShortName="sen",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_evap",LongName="inst_evap_rate",ShortName="evap",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_swndr",LongName="inst_down_sw_ir_dir_flx",ShortName="swndr",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_swndf",LongName="inst_down_sw_ir_dif_flx",ShortName="swndf",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_swvdr",LongName="inst_down_sw_vis_dir_flx",ShortName="swvdr",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call advertise_support(exportState,standard_name="Faxa_swvdf",LongName="inst_down_sw_vis_dif_flx",ShortName="swvdf",rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif



    ! -----------------------------------------------
    ! Set internal fv3 clock timing information
    ! -----------------------------------------------
    call fms_mpp_clock_end (initClock) !end initialization
    mainClock = fms_mpp_clock_id( '-Main Loop' )
    call fms_mpp_clock_begin(mainClock) !begin main loop


  contains

    ! ----------------------------------
    subroutine support_get_string(model,input_string,input_var,rc)
      type(ESMF_GridComp), intent(in) :: model
      character(*), intent(in)        :: input_string
      character(*), intent(inout)     :: input_var
      integer,intent(out)             :: rc
      ! local
      character(ESMF_MAXSTR)          :: valueString
      logical                         :: isPresent, isSet
      character(len=*),parameter      :: subname='(fv3_shield_cap:Advertise:support_get_string)'

      call NUOPC_CompAttributeGet(model, name=trim(input_string), value=valueString, isPresent=isPresent, isSet=isSet, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      if (isPresent .and. isSet) then
          input_var = trim(valueString)
      else
          valueString = trim(input_var)//' (default)'
      endif

      call ESMF_LogWrite(trim(subname)//': '//trim(input_string)//' = '//trim(valueString), ESMF_LOGMSG_INFO, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return
    end subroutine support_get_string
    ! ----------------------------------

    ! ----------------------------------
    subroutine support_get_integer(model,input_string,input_var,rc)
      type(ESMF_GridComp), intent(in) :: model
      character(*), intent(in)        :: input_string
      integer, intent(inout)          :: input_var
      integer,intent(out)             :: rc
      ! local
      character(ESMF_MAXSTR)          :: valueString
      logical                         :: isPresent, isSet
      character(len=*),parameter      :: subname='(fv3_shield_cap:Advertise:support_get_integer)'

      call NUOPC_CompAttributeGet(model, name=trim(input_string), value=valueString, isPresent=isPresent, isSet=isSet, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      if (isPresent .and. isSet) then
          read(valueString,*) input_var
          input_var = ESMF_UtilString2Int(valueString, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
      else
          valueString = ESMF_UtilStringInt2String(input_var, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
          valueString = valueString//' (default)'
      endif

      call ESMF_LogWrite(trim(subname)//': '//trim(input_string)//' = '//trim(valueString), ESMF_LOGMSG_INFO, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return
    end subroutine support_get_integer
    ! ----------------------------------

    ! ----------------------------------
    subroutine support_get_logical(model,input_string,input_var,rc)
      type(ESMF_GridComp), intent(in) :: model
      character(*), intent(in)        :: input_string
      logical, intent(inout)          :: input_var
      integer,intent(out)             :: rc
      ! local
      character(ESMF_MAXSTR)          :: valueString
      logical                         :: isPresent, isSet
      character(len=*),parameter      :: subname='(fv3_shield_cap:Advertise:support_get_logical)'

      call NUOPC_CompAttributeGet(model, name=trim(input_string), value=valueString, isPresent=isPresent, isSet=isSet, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      if (isPresent .and. isSet) then
          read(valueString,*) input_var
      else
          valueString = '(default)'
      endif

      if (input_var) then
          call ESMF_LogWrite(trim(subname)//': '//trim(input_string)//' = .true. '//trim(valueString), ESMF_LOGMSG_INFO, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
      else
          call ESMF_LogWrite(trim(subname)//': '//trim(input_string)//' = .false. '//trim(valueString), ESMF_LOGMSG_INFO, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
      endif

    end subroutine support_get_logical
    ! ----------------------------------

    subroutine advertise_support(state,standard_name,LongName,ShortName,rc)
      ! https://earthsystemmodeling.org/docs/release/ESMF_8_4_2/NUOPC_refdoc.pdf
      ! (page 98-99)
      type(ESMF_State), intent(inout) :: state
      character(*)                    :: standard_name, LongName, ShortName
      integer, intent(out)            :: rc


      call NUOPC_Advertise(state, &
                           TransferOfferGeomObject="will provide", &
                           StandardName=trim(standard_name), &
                           LongName=trim(LongName), &
                           ShortName=trim(ShortName), &
                           name=trim(standard_name), &
                           SharePolicyField="share", &
                           SharePolicyGeomObject="share", &
                           rc=rc)

      if (CheckError(rc,__LINE__,__FILE__)) return

    end subroutine advertise_support

    subroutine string_split(string_input,string_array,sep_in)
        ! https://fortran-lang.org/ja/learn/best_practices/allocatable_arrays/
        character(len=*), intent(in) :: string_input
        character(len=32), dimension(:), allocatable, intent(out) :: string_array
        character(len=*), intent(in), optional :: sep_in
        character(1) :: sep = ','
        integer :: i,n

        ! First check if the string is blank
        if (trim(string_input)=='') then
            return
        endif

        ! Overwrite default delimiter
        if (present(sep_in)) sep = sep_in

        ! Find the number of strings
        n = 1
        do i=1, len_trim(string_input)
            if (string_input(i:i) == sep) n = n + 1
        enddo

        ! Read in the list of strings
        allocate (string_array(n))
        read (unit=string_input,fmt=*) string_array

      end subroutine string_split

  end subroutine Advertise
  

  !-----------------------------------------------------------------------------
  subroutine Realize(model, rc) 
  !-----------------------------------------------------------------------------
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)         :: importState, exportState
    type(ESMF_Field)         :: field
    type(ESMF_Grid)          :: gridIn
    type(ESMF_Grid)          :: gridOut

    type(ESMF_Config)        :: config

    logical :: dodebug
    integer :: i,j


    dodebug = dodebug_Realize
    
    rc = ESMF_SUCCESS  

    !-----------------------------------------
    ! query for importState and exportState
    !-----------------------------------------
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------------
    ! read runtime parameters from config
    ! ----------------------------------------------
    call ESMF_GridCompGet(model, config=config, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !-----------------------------------------
    ! Create grid for ESMF fields
    !-----------------------------------------
    if (use_mosaic) then !(default)
        call grid_create_mosaic(model, grid=gridOut, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    else  !(EMC version)
        !NOTE: model is defined as an inout argument here, but in-only above
        call grid_create_internal(model, grid=gridOut, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! for now use the same grid for imports and exports, 
    ! assuming conversions will be done by mediator
    gridIn = gridOut 

    !-----------------------------------------
    ! Allocate pointers for import and export fields
    !-----------------------------------------
    !STEVE: review how it is done with WW3 here:
    !       https://github.com/wavespotter/EarthSystemModel/blob/aaaa9d379ea351794243d654ff5f94ea58d60edf/wmesmfmd.F90#L454


    !--------------------------------------------------------------------------
    ! Exportable fields
    !--------------------------------------------------------------------------

    !------------------
    ! exportable field: eastward_wind_at_10m_height
    !------------------
!   field = ESMF_FieldCreate(name="Sa_u10m", grid=gridOut, typekind=ESMF_TYPEKIND_R8, rc=rc)
!   if (CheckError(rc,__LINE__,__FILE__)) return


    !STEVE:ISSUE: should winds / wave info be on center of corner? (ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER,ESMF_STAGGERLOC_EDGE1)
    !NOTE: " 2) the bounds and counts retrieved from GridGet are DE specific or equivalently PET specific, which means that the Fortran array shape could be different from one PET to another."
    ! from : https://earthsystemmodeling.org/docs/release/ESMF_8_1_0/ESMF_refdoc.pdf
    !        (page 345)

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_u10m",farray_u,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: northward_wind_at_10m_height
    !------------------
!   field = ESMF_FieldCreate(name="Sa_v10m", grid=gridOut, typekind=ESMF_TYPEKIND_R8, rc=rc)
!   if (CheckError(rc,__LINE__,__FILE__)) return

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_v10m",farray_v,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: eastward_neutral_wind_at_10m_height
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_u10n",farray_un,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: northward_neutral_wind_at_10m_height
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_v10n",farray_vn,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    
    !------------------
    ! exportable field: atmosphere temperature at 2m height
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_t2m",farray_t2m,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: atmospheric surface density
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_rhoa",farray_rhoa,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: underlying surface temperature
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_tsfc",farray_tsfc,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! exportable field: astdiff (post-processed)
    !------------------

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_astdiff",farray_astdiff,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Allocate space for this field's data
    call fpointer_allocate(field,gridOut,"Sa_oceanfrac",farray_oceanfrac,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    
    !------------------
    ! exportable field: cpl_scalar grid information for CMEPS
    !------------------

    ! Create a grid for the scalar data transfer to CMEPS
    ! from: https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/cpl/module_cplscalars.F90#L62
    if (flds_scalar_num>0) then
        call support_Realize_scalars(exportState, rc=rc)
    endif

    ! --------------------------------------------------------------------------------
    ! The remaining fields are only used if dodebug_setExport_writenetcdf=.true.
    ! --------------------------------------------------------------------------------
    if (dodebug_setExport_writenetcdf) then

        !------------------
        ! exportable field: sea/land/ice mask
        !------------------

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_slmsk",farray_slmsk,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_landfrac",farray_landfrac,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_lakefrac",farray_lakefrac,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_fice",farray_fice,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
    
        !------------------
        ! exportable field: surface pressure
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_psurf",farray_psurf,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: atmospheric humidity at 2m
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_q2m",farray_q2m,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: surface ice temperature
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_tisfc",farray_tisfc,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: friction velocity
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_uustar",farray_uustar,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: ffmm
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_ffmm",farray_ffmm,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: f10m
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_f10m",farray_f10m,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        !------------------
        ! exportable field: pbl height (m)
        !------------------
    
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_hpbl",farray_hpbl,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    endif

    !------------------
    ! exportable field: mixed layer model (MLM) fields
    !------------------

    if (use_mlm) then
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_ts_som",farray_ts_som,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_tml",farray_tml,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_mld",farray_mld,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_huml",farray_huml,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Sa_hvml",farray_hvml,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
    endif


    ! Add FLUXES - 7/16/24

    if (use_aofluxes) then
        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_taux",farray_taux,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_tauy",farray_tauy,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_rain",farray_rain,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_lwnet",farray_lwnet,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_sen",farray_sen,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_evap",farray_evap,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_swndr",farray_swndr,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_swndf ",farray_swndf ,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_swvdr",farray_swvdr,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

        ! Allocate space for this field's data
        call fpointer_allocate(field,gridOut,"Faxa_swvdf",farray_swvdf,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 

    endif
    
    !--------------------------------------------------------------------------
    ! Importable fields
    !--------------------------------------------------------------------------

    do i=1,Sa_import_len
        if (trim(Sa_import_names(i)).eq.'Sw_charno') then
            ! Allocate space for this field's data
            call fpointer_allocate(field,gridIn,"Sw_charno",farray_c,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

            call NUOPC_Realize(importState, field=field, rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

        elseif (trim(Sa_import_names(i)).eq.'Sw_z0rlen') then
            ! Allocate space for this field's data
            call fpointer_allocate(field,gridIn,"Sw_z0rlen",farray_z,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

            call NUOPC_Realize(importState, field=field, rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

        elseif (trim(Sa_import_names(i)).eq.'So_sst') then
            ! Allocate space for this field's data
            call fpointer_allocate(field,gridIn,"So_sst",farray_sst,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

            call NUOPC_Realize(importState, field=field, rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return

!       elseif (trim(Sa_import_names(i)).eq.'So_ssu') then
!           ! Allocate space for this field's data
!           call fpointer_allocate(field,gridIn,"So_ssu",farray_ssu,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return

!           call NUOPC_Realize(importState, field=field, rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return

!       elseif (trim(Sa_import_names(i)).eq.'So_ssv') then
!           ! Allocate space for this field's data
!           call fpointer_allocate(field,gridIn,"So_ssv",farray_ssv,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return

!           call NUOPC_Realize(importState, field=field, rc=rc)
!           if (CheckError(rc,__LINE__,__FILE__)) return

        else
            output_string = "Sa_import_names(i) = "//trim(Sa_import_names(i))//" is not a valid import variable. Choose from: Sw_z0rlen, Sw_charno, So_sst" !, So_ssu, So_ssv"
            call ESMF_LogWrite("fv3_shield_cap::Realize:: "//trim(output_string), ESMF_LOGMSG_INFO)
        endif
    enddo


    contains


    ! -------------- - - - -- - - -- ------- - - - - - - - - -
    subroutine fpointer_allocate(field,grid,fieldname,farray,staggerloc,rc)
    !STEVE: Make sure the field arrays are allocated and ready to recieve data (8/30/23)
    ! https://earthsystemmodeling.org/docs/release/ESMF_8_1_0/ESMF_refdoc.pdf
    ! (page 333)
      type(ESMF_Field), intent(inout)  :: field
      type(ESMF_Grid),  intent(in)     :: grid
      character(*), intent(in)         :: fieldname
      real(ESMF_KIND_R8), dimension(:,:), pointer, intent(inout) :: farray ! for grid
      type (ESMF_StaggerLoc),intent(in):: staggerloc
      integer, intent(out)             :: rc
      type(ESMF_DistGrid)              :: distgrid
      type(ESMF_Array)                 :: array
      integer                          :: xdim, ydim, zdim
      real(ESMF_KIND_R8), parameter    :: l_fill_value = 0._ESMF_KIND_R8
      integer, dimension(2)            :: fa_shape, gec, &
                                          compLBnd, compUBnd, exclLBnd, exclUBnd, totalLBnd, totalUBnd, &
                                          comp_count, excl_count, total_count
      real(ESMF_KIND_R8), parameter    :: PI = 3.14159265
      integer :: i,j
    
      ! create a 2D data Field from a Grid and Array.
      call ESMF_GridGet(grid=grid, staggerloc=staggerloc, distgrid=distgrid, rc=rc) !exclusiveCount=gec, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      call ESMF_GridGetFieldBounds(grid=grid, localDe=0, staggerloc=staggerloc, totalCount=fa_shape, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      ! Allocate an array to store the field data
      allocate( farray(fa_shape(1), fa_shape(2)) ) !, fa_shape(3)) )

      ! create an Array
      array = ESMF_ArrayCreate(distgrid, farray, indexflag=ESMF_INDEX_DELOCAL, rc=rc)

      ! create a Field
      field = ESMF_FieldCreate(name=trim(fieldname), grid=grid, array=array, rc=rc) !indexflag=ESMF_INDEX_DELOCAL, rc=rc)

      ! Fill in some arbitrary initialization values

      ! retrieve the Fortran data pointer from the Field
      call ESMF_FieldGet(field=field, localDe=0, farrayPtr=farray, rc=rc)

      ! retrieve the Fortran data pointer from the Field and bounds
      call ESMF_FieldGet(field=field, localDe=0, &
                         farrayPtr=farray, &
                         computationalLBound=compLBnd, computationalUBound=compUBnd, &
                         exclusiveLBound=exclLBnd, exclusiveUBound=exclUBnd, &
                         totalLBound=totalLBnd, totalUBound=totalUBnd, &
                         computationalCount=comp_count, &
                         exclusiveCount=excl_count, &
                         totalCount=total_count, &
                         rc=rc)

      ! iterate through the total bounds of the field data pointer
      ! to create an arbitrary wavelike initial condition
      do j = totalLBnd(2), totalUBnd(2)
          do i = totalLBnd(1), totalUBnd(1)
              farray(i, j) = sin(2*i/total_count(1)*PI) + &
                             sin(4*j/total_count(2)*PI)
          enddo
      enddo

    end subroutine fpointer_allocate
    ! -------------- - - - -- - - -- ------- - - - - - - - - -


    ! See:
    ! https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/cpl/module_cplfields.F90#L495
    ! https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/cpl/module_cplscalars.F90#L61
    subroutine support_Realize_scalars(exportState, rc)
        type(ESMF_State), intent(inout) :: exportState
        integer, intent(out)            :: rc

        ! For CMEPS cpl_scalar:
        type(ESMF_Distgrid) :: distgrid
        type(ESMF_Grid)     :: grid
        real(ESMF_KIND_R8), parameter :: d_fill_value = 0._ESMF_KIND_R8
        real(ESMF_KIND_R8) :: l_fill_value

        ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
        distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        grid = ESMF_GridCreate(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        field = ESMF_FieldCreate(name=trim(flds_scalar_name),  &
                                 grid=grid,  &
                                 typekind=ESMF_TYPEKIND_R8,  &
                                 ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/),  &
                                 gridToFieldMap=(/2/),  &
                                 rc=rc) ! num of scalar values
        if (CheckError(rc,__LINE__,__FILE__)) return

        call NUOPC_Realize(exportState, field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! -- initialize field value
        l_fill_value = d_fill_value
        call ESMF_FieldFill(field, dataFillScheme="const", const1=l_fill_value, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    end subroutine support_Realize_scalars


  end subroutine Realize


  !-----------------------------------------------------------------------------
  subroutine DataInitialize(model, rc)
  !-----------------------------------------------------------------------------
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)                  :: exportState
    type(ESMF_Field)                  :: field
    integer                           :: i, j

    type(ESMF_VM)                     :: vm
    integer                           :: localDeCount, localPet
    integer                           :: fieldCount, n
    character(len=64),allocatable     :: fieldNameList(:)

    type(ESMF_Config)           :: config

    character(128) :: varname = 'noname'

    ! local variables
    character(len=*),parameter         :: subname='(fv3_shield_cap:DataInitialize)'
    character(ESMF_MAXSTR) :: cname
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TRACE start: DataInitialize ----------------------------------------------------
    call ESMF_TraceRegionEnter(trim(subname), rc=rc)

    ! ----------------------------------------
    ! Prep
    ! ----------------------------------------
    call ESMF_GridCompGet(model, name=cname, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    if (verbosity.gt.0) call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! Access the atmosphere component export state 
    ! so we can initialize the exchanged fields
    call NUOPC_ModelGet(model, exportState=exportState, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------
    ! Access internal atmospheric model state 
    ! and assign it to 'export' here
    ! ----------------------------------------
    call setExport(model, exportState, rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------
    ! Do some checks 
    ! ----------------------------------------

    ! Query for the number of fields in the atmospheric component
    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Get the list of names for each of those fields
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    do n=1, fieldCount
        call ESMF_LogWrite("ATM_DataInitialize:: updating field: "//trim(fieldNameList(n)), ESMF_LOGMSG_INFO)

        call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        if (verbosity.gt.0) call ESMF_LogWrite("ATM_DataInitialize:: updated. "//trim(fieldNameList(n)), ESMF_LOGMSG_INFO)
    enddo
    deallocate(fieldNameList)

    ! check whether all Fields in the exportState are "Updated"
    if (verbosity.gt.0) call ESMF_LogWrite("FV3-SHiELD - check for InitializeDataComplete...", ESMF_LOGMSG_INFO)
    if (NUOPC_IsUpdated(exportState)) then
        call NUOPC_CompAttributeSet(model, name="InitializeDataComplete", value="true", rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_LogWrite("FV3-SHiELD - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! indicate that data initialization is complete (breaking out of init-loop)
    call NUOPC_CompAttributeSet(model, name="InitializeDataComplete", value="true", rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Write out the initialized fields
    if (dodebug_DataInitialize) then

        ! -----------------------
        ! Write all fields in the export state
        ! -----------------------
        call NUOPC_Write(exportState, fileNamePrefix="ATM_DataInitialize_exportState_tile*_", overwrite=.true., rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    endif

    if (verbosity.gt.0) call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

    call ESMF_TraceRegionExit(trim(subname), rc=rc)
    ! TRACE end: DataInitialize ------------------------------------------------------

  end subroutine DataInitialize


  !-----------------------------------------------------------------------------
  subroutine CheckImport(model, rc)
  !-----------------------------------------------------------------------------
  !
  ! Checks model component clock and compares to internal model clock
  !
  ! See EMC fv3atm example:
  ! https://github.com/NOAA-EMC/fv3atm/blob/a13a239a746cb95c74bbe5841f300cb75f8b80d9/fv3_cap.F90#L1282
  !-----------------------------------------------------------------------------
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='fv3_shield_cap::CheckImport'
    integer                     :: n, nf
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: currTime, invalidTime
    type(ESMF_State)            :: importState
    logical                     :: isValid
    type(ESMF_Field), pointer   :: fieldList(:)
    character(len=128)          :: fldname
    integer                     :: date_esmf(6), date_fv3(6)

    rc = ESMF_SUCCESS

    ! query the Component for its clock
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    date_esmf(1:6) = 0
    call ESMF_TimeGet(time=currTime,yy=date_esmf(1),mm=date_esmf(2),dd=date_esmf(3),h=date_esmf(4), m=date_esmf(5),s=date_esmf(6),rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! set up invalid time (by convention)
!   call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, h=00, m=00, s=00, rc=rc)
!   if (CheckError(rc,__LINE__,__FILE__)) return

    ! Get the current time from the model

    !----- compute current date ------
    ! https://github.com/NOAA-GFDL/FMS/blob/06b94a7f574e7794684b8584391744ded68e2989/time_manager/time_manager.F90#L1145
    ! NOTE: the first argument is input and the rest are returned
    call fms_time_manager_get_date (Atm%Time, date_fv3(1), date_fv3(2), date_fv3(3), date_fv3(4), date_fv3(5), date_fv3(6))

    ! Compare the ESMF model clock and the FV3 internal clock
    if (.not. all(date_esmf.EQ.date_fv3)) then
        print *, trim(subname)//": ERROR - the ESMF internal clock and fv3 model clock do not agree:"
        print *, "date_esmf = ", date_esmf(1:6)
        print *, "date_fv3  = ", date_fv3(1:6)
        print *, "EXITING..."
        stop(1)
    endif

  end subroutine CheckImport


  !-----------------------------------------------------------------------------
  subroutine Advance(model, rc)
  !-----------------------------------------------------------------------------
  !
  ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
  !
  ! Because of the way that the internal Clock was set in SetClock(),
  ! its timeStep is likely smaller than the parent timeStep. As a consequence
  ! the time interval covered by a single parent timeStep will result in
  ! multiple calls to the Advance() routine. Every time the currTime
  ! will come in by one internal timeStep advanced. This goes until the
  ! stopTime of the internal Clock has been reached.    
  !
  !-----------------------------------------------------------------------------
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep

    ! FMS time type
    type (FmsTime_type)            :: Time
    integer                     :: date_esmf(6)
    
    ! Just for printing:
    type(ESMF_VM)               :: vm
    integer :: localPet

    character(16)               :: timestamp_str

    rc = ESMF_SUCCESS

    !-----------------------------------------
    ! query the Component for its clock, importState and exportState
    !-----------------------------------------
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! Print the current time
    !------------------
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="ATM::Advance::------>Advancing FV3-SHiELD ATM/LND/OML from: ", unit=msgString, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! Write to log
    !------------------
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !------------------
    ! Get the current time and time step
    !------------------
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    date_esmf=0
    call ESMF_TimeGet (currTime,                           &
                       YY=date_esmf(1), MM=date_esmf(2), DD=date_esmf(3), &
                       H=date_esmf(4),  M =date_esmf(5), S =date_esmf(6), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return


    !-----------------------------------------
    ! Get the vm/localPet to print the current time on root pe
    !-----------------------------------------
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
 
    call ESMF_vmget(vm, localPet=localPet, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
 
    if (localPet == 0) write(*,'(A,6I5)') 'fv3_shield_cap.F90:Advance:: ESMF CurrTime =', date_esmf

    Time_atmos = fms_time_manager_set_date (date_esmf(1), date_esmf(2), date_esmf(3), date_esmf(4), date_esmf(5), date_esmf(6))


    !-----------------------------------------
    !STEVE: extract importState here...
    !-----------------------------------------

    ! TRACE start: getImport ---------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::getImport", rc=rc)

    if (dodebug_bypassImport) then
        call ESMF_LogWrite("fv3_shield_cap:: [DEBUG] bypassing getImport...", ESMF_LOGMSG_INFO, rc=rc)
    else
        call getImport(model,importState,rc)
    endif
    if (CheckError(rc,__LINE__,__FILE__)) return

    call ESMF_TraceRegionExit("fv3_shield_cap::getImport", rc=rc)
    ! TRACE end: getImport -----------------------------------------------------

    
    !----------------------------------------------------------------------------------
    ! FV3-SHiELD model advance routines:
    ! Step the model dynamics/physics/state
    !----------------------------------------------------------------------------------

    ! TRACE start: update ------------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::update", rc=rc)

    if (dodebug_Advance) print *, "fv3_shield_cap:: calling update_atmos_model_dynamics..."
    call ESMF_TraceRegionEnter("fv3_shield_cap::update::dynamics", rc=rc)
    call update_atmos_model_dynamics (Atm)
    call ESMF_TraceRegionExit("fv3_shield_cap::update::dynamics", rc=rc)

    if (dodebug_Advance) print *, "fv3_shield_cap:: calling update_atmos_radiation_physics..."
    call ESMF_TraceRegionEnter("fv3_shield_cap::update::physics", rc=rc)
    call update_atmos_radiation_physics (Atm)
    call ESMF_TraceRegionExit("fv3_shield_cap::update::physics", rc=rc)

    if (dodebug_Advance) print *, "fv3_shield_cap:: calling update_atmos_model_state..."
    call ESMF_TraceRegionEnter("fv3_shield_cap::update::state", rc=rc)
    call update_atmos_model_state (Atm)
    call ESMF_TraceRegionExit("fv3_shield_cap::update::state", rc=rc)

    call ESMF_TraceRegionExit("fv3_shield_cap::update", rc=rc)
    ! TRACE end: update --------------------------------------------------------


    !-----------------------------------------
    !STEVE: extract exportState from Atm here...
    !-----------------------------------------
    
    ! TRACE start: setExport ---------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::setExport", rc=rc)

    if (dodebug_bypassExport) then
        call ESMF_LogWrite("fv3_shield_cap:: [DEBUG] bypassing getExport...", ESMF_LOGMSG_INFO, rc=rc)
    else
        call setExport(model, exportState, rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif
    
    call ESMF_TraceRegionExit("fv3_shield_cap::setExport", rc=rc)
    ! TRACE end: setExport -----------------------------------------------------


    !-----------------------------------------
    !--- store intermediate restart files
    !-----------------------------------------
    ! Start insert
    if (intrm_rst) then
        if (nc /= num_cpld_calls) then
        if (intrm_rst_1step .and. nc == 1) then
            timestamp = fms_time_manager_date_to_string (Time_atmos)
            call atmos_model_restart(Atm, timestamp)
            call coupler_restart(timestamp)
        endif
        if (Time_atmos == Time_restart .or. Time_atmos == Time_restart_aux) then
            if (Time_atmos == Time_restart) then
            timestamp = fms_time_manager_date_to_string (Time_restart)
            else
            timestamp = fms_time_manager_date_to_string (Time_restart_aux)
            endif
            call atmos_model_restart(Atm, timestamp)
            call coupler_restart(timestamp)
            if (Time_atmos == Time_restart) &
                Time_restart = Time_restart + Time_step_restart
            if ((restart_secs_aux > 0 .or. restart_days_aux > 0) .and. &
                Time_atmos == Time_restart_aux .and. &
                Time_restart_aux < Time_restart_end_aux) then
            Time_restart_aux = Time_restart_aux + Time_step_restart_aux
            endif
        endif
        endif
    endif

    call fms_memutils_print_memuse_stats('after full step')
    ! End insert

    !------------------
    ! Print the new time after model timestep
    !------------------
    call ESMF_TimePrint(currTime + timeStep, &
      preString="ATM::Advance::---------------------> to: ", unit=msgString, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! NOTE: this is just to keep track of the time step in the module-global
    !       variable Time_atmos, used by coupler_end, to maintain consistency
    !       with the ESMF automatic time advance
    if (dodebug_Advance .and. localPet == 0) print *, "fv3_shield_cap:: stepping clock for next time step..."
    Time_atmos = Time_atmos + Time_step_atmos  !STEVE: replace with NUOPC clock

    !------------------
    ! Write to log
    !------------------
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

  end subroutine Advance


  !-----------------------------------------------------------------------------
  subroutine getImport(model, importState, rc)
  !-----------------------------------------------------------------------------

    use atmos_model_mod, only: Atm_block, IPD_Data
!   use atmos_model_mod, only: IPD_Control, IPD_Data, IPD_Diag, IPD_Restart

    type(ESMF_GridComp) :: model
    type(ESMF_State),intent(inout) :: importState
    integer,intent(out) :: rc

    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: dataPtr_ch(:,:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_z0(:,:)   => null()

    integer :: i,j, isc,iec,jsc,jec, nb,ix
    character(16)               :: timestamp_str

    !  1. Purpose :
    !
    !     Get import fields and put in internal data structures

    !STEVE: update to just loop through all fields in the importState - 11/10/23
    do i=1,Sa_import_len
        if (trim(Sa_import_names(i)).eq.'Sw_charno') then
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (before) support_getImport Sw_charno"
            call support_getImport(importState,itemName='Sw_charno',rc=rc)    
            if (CheckError(rc,__LINE__,__FILE__)) return
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (after) support_getImport Sw_charno"

        elseif (trim(Sa_import_names(i)).eq.'Sw_z0rlen') then
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (before) support_getImport Sw_z0rlen"
            call support_getImport(importState,itemName='Sw_z0rlen',rc=rc)    
            if (CheckError(rc,__LINE__,__FILE__)) return
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (after) support_getImport Sw_z0rlen"

        elseif (trim(Sa_import_names(i)).eq.'So_sst') then
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (before) support_getImport So_sst"
            call support_getImport(importState,itemName='So_sst',rc=rc)    
            if (CheckError(rc,__LINE__,__FILE__)) return
            if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (after) support_getImport So_sst"
 
!       elseif (trim(Sa_import_names(i)).eq.'So_ssu') then
!           if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (before) support_getImport So_ssu"
!           call support_getImport(importState,itemName='So_ssu',rc=rc)    
!           if (CheckError(rc,__LINE__,__FILE__)) return
!           if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (after) support_getImport So_ssu"

!       elseif (trim(Sa_import_names(i)).eq.'So_ssv') then
!           if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (before) support_getImport So_ssv"
!           call support_getImport(importState,itemName='So_ssv',rc=rc)    
!           if (CheckError(rc,__LINE__,__FILE__)) return
!           if (dodebug_getImport) print *, "fv3_shield_cap.F90::getImport:: (after) support_getImport So_ssv"

        else
            output_string = "Sa_import_names(i) = "//trim(Sa_import_names(i))//" is not a valid import variable. Choose from: Sw_z0rlen, Sw_charno, So_sst" !, So_ssu, So_ssv"
            call ESMF_LogWrite("fv3_shield_cap::getImport:: "//trim(output_string), ESMF_LOGMSG_INFO)
        endif
    enddo


    ! Write out the data fields in netcdf form
    if (dodebug_getImport_writenetcdf) then

        call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::getImport::dodebug_getImport_writenetcdf: WRITING..."), ESMF_LOGMSG_INFO, rc=rc)

        ! -----------------------
        ! Write all fields in the import state
        ! -----------------------
        if (dodebug_setExport_multistep) then

            timestamp_str=fms_time_manager_date_to_string(Time_atmos)  !STEVE: using built-in date-to-string converter from fv3
            call NUOPC_Write(importState, fileNamePrefix="ATM_getImport_importState_"//trim(timestamp_str)//"_tile*_", overwrite=.true., rc=rc)

            if (CheckError(rc,__LINE__,__FILE__)) return
        else
            call NUOPC_Write(importState, fileNamePrefix="ATM_getImport_importState_tile*_", overwrite=.true., rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        endif

    else
        call ESMF_LogWrite(trim("[NO DEBUG] fv3_shield_cap::getImport::dodebug_getImport_writenetcdf: SKIPPING."), ESMF_LOGMSG_INFO, rc=rc)
    endif


  contains 
          
    subroutine support_getImport(importState,itemName,rc)
      use atmos_model_mod, only: Atm_block, IPD_Data

      type(ESMF_State),intent(inout) :: importState
      character(*), intent(in)       :: itemName
      integer,intent(out)            :: rc
      type(ESMF_Field)               :: field
      real(ESMF_KIND_R8), pointer    :: dataPtr(:,:)   => null()

      integer :: i,j, isc,iec,jsc,jec, nb,ix

      call ESMF_StateGet(importState, itemName=itemName, field=field, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      !STEVE: now, set the values in the fv3-shield model
      !       IPD interface data layer
      isc = Atm_block%isc
      iec = Atm_block%iec
      jsc = Atm_block%jsc
      jec = Atm_block%jec

      do j=jsc,jec
        do i=isc,iec

          if (dodebug_getImport) print *, "i,j = ", i,j

          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)

          if (dodebug_getImport) then
              print *, "i,j,nb,ix = ", i,j,nb,ix
              print *, "associated(IPD_Data(nb)%Sfcprop%zorl     = ", associated(IPD_Data(nb)%Sfcprop%zorl)
              print *, "associated(IPD_Data(nb)%Sfcprop%charnock = ", associated(IPD_Data(nb)%Sfcprop%charnock)
              print *, "associated(IPD_Data(nb)%Sfcprop%tsfc     = ", associated(IPD_Data(nb)%Sfcprop%tsfc)
              print *, "isc,iec, jsc,jec = ", isc,iec,jsc,jec
          endif

          ! Only update ocean points
          ! Note: the land/sea mask in the fv3 is actually sea/land/ice, and represented by 0/1/2
          if (IPD_Data(nb)%Sfcprop%slmsk(ix) == 0) then
            if (dodebug_getImport) print *, "ocean point"

            if (trim(itemName)==trim('Sw_z0rlen')) then
                IPD_Data(nb)%Sfcprop%zorl(ix)  = dataPtr(i-isc+1,j-jsc+1) * 100.0 !! convert m --> cm expected by atmosphere
!               IPD_Data(nb)%Sfcprop%zorlo(ix) = dataPtr(i-isc+1,j-jsc+1)  ! Ocean roughness length
                !STEVE: "Error: 'zorlo' at (1) is not a member of the 'gfs_sfcprop_type' structure; did you mean 'zorl'?"
            elseif (trim(itemName)==trim('Sw_charno')) then
                IPD_Data(nb)%Sfcprop%charnock(ix)  = dataPtr(i-isc+1,j-jsc+1)
!           elseif (trim(itemName)==trim('So_sst')) then
!               IPD_Data(nb)%Sfcprop%tsfc(ix)  = dataPtr(i-isc+1,j-jsc+1)
            else
                print *, "fv3_shield_cap::setImport::support_setImport:: variable name {"//trim(itemName)//"} not recognized. EXITING..."
                stop(1)
            endif

          else
            if (dodebug_getImport) print *, "non-ocean point"
          endif

        enddo
      enddo


   end subroutine support_getImport

  end subroutine getImport


  !-----------------------------------------------------------------------------
  subroutine setExport(model, exportState, rc)
  ! See:
  ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/ff84897f63000f085a5d45c1deaa826ad51989a9/GFS_layer/GFS_typedefs.F90
  !-----------------------------------------------------------------------------

    use atmos_model_mod, only: Atm_block, IPD_Data
!   use atmos_model_mod, only: IPD_Control, IPD_Data, IPD_Diag, IPD_Restart
    use atmosphere_mod, only: Atm, mygrid
    use module_block_data,  only: block_atmos_copy, block_data_copy

    ! Driver version used here:
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/FV3-202210-public/atmos_drivers/coupled/atmos_model.F90
    ! newer version, moved out of the shield repo:
    ! https://github.com/NOAA-GFDL/atmos_drivers/blob/main/SHiELD/atmos_model.F90

    !STEVE: in atmosphere.F90::atmosphere_mod:: type(fv_atmos_type), allocatable, target :: Atm(:)
    !       This contains the model state information as advanced in subroutine atmosphere_dynamics:
    !       https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/e45141a6369b5a6c4676a8f12fc5d1509b51c6ba/driver/SHiELD/atmosphere.F90
    !
    !       Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz,
    !       Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,
    !
    !       Additional details defined here:
    !       https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/fv_arrays.F90
    !
    !       u_srf, v_srf,
    !       
    ! In atmos_model.F90::atmos_model_mod:: type (atmos_data_type) :: Atm
    ! This contains grid information, e.g.
    ! pelist, mygrid, pe (current pe), Time, Time_step, Time_init, grid
    !
    ! NOTE:
    ! IPD_* aliases GFD_ types
    ! 
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/GFS_layer/GFS_typedefs.F90#L113
    !!----------------------------------------------------------------
    !! GFS_statein_type
    !!   prognostic state variables with layer and level specific data
    !!----------------------------------------------------------------
    !  type GFS_statein_type
    !!--- prognostic variables
    !   real (kind=kind_phys), pointer :: pgr  (:)     => null()  !< surface pressure (Pa) real
    !   real (kind=kind_phys), pointer :: ugrs (:,:)   => null()  !< u component of layer wind
    !   real (kind=kind_phys), pointer :: vgrs (:,:)   => null()  !< v component of layer wind
    !
    ! Also:
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/GFS_layer/GFS_typedefs.F90#L337
    !!---------------------------------------------------------------------
    !! GFS_coupling_type
    !!   fields to/from other coupled components (e.g. land/ice/ocean/etc.)
    !!---------------------------------------------------------------------
    !  type GFS_coupling_type
    !   real (kind=kind_phys), pointer :: u10mi_cpl  (:) => null()   !< instantaneous U10m
    !   real (kind=kind_phys), pointer :: v10mi_cpl  (:) => null()   !< instantaneous V10m
    !
    ! get the values from the fv3-shield model and place in export fields...
    ! Expected shape:
    ! (From: https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/fv_arrays.F90)
    ! allocate (    Atm%u(isd:ied  ,jsd:jed+1,npz) )
    ! allocate (    Atm%v(isd:ied+1,jsd:jed  ,npz) )
    ! allocate (    Atm%u_srf(is:ie,js:je) )
    ! allocate (    Atm%v_srf(is:ie,js:je) )
    ! allocate (    Atm%ua(isd:ied  ,jsd:jed  ,npz) )
    ! allocate (    Atm%va(isd:ied  ,jsd:jed  ,npz) )
    ! allocate (    Atm%uc(isd:ied+1,jsd:jed  ,npz) )
    ! allocate (    Atm%vc(isd:ied  ,jsd:jed+1,npz) )
    ! allocate (    Atm%w(isd:ied, jsd:jed  ,npz  ) )
    !
    ! !---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
    ! real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt, qv_dt
    !
    ! Note that "Atm%gridstruct" has details of the grid cells and offsets, if needed by ESMF

    type(ESMF_GridComp) :: model
    type(ESMF_State),intent(inout) :: exportState
    integer,intent(out) :: rc

    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), dimension(:,:), pointer :: dataPtr_u_r8
    real(ESMF_KIND_R8), dimension(:,:), pointer :: dataPtr_v_r8
    real(ESMF_KIND_R8), dimension(:,:), pointer :: dataPtr_r8

    integer :: ni,nj,nk, i,j, ib,jb,nb,ix, blen
    integer :: isc,iec,jsc,jec
    integer :: isd,ied,jsd,jed
    integer :: ierr

    type(esmf_vm)                          :: vm
    integer                                :: mpi_comm_esmf, localPet

    integer                                     :: n,rank
    logical                                     :: isFound
    type(ESMF_TypeKind_Flag)                    :: datatype
    character(len=ESMF_MAXSTR)                  :: fieldName

    character(128) :: varname = 'noname'
    character(16)                               :: timestamp_str

    type(ESMF_Clock)          :: clock
    type(ESMF_Time)           :: start_time, current_time

    rc = ESMF_SUCCESS

    !  1. Purpose :
    !
    !     Set export fields taken from internal atmosphere model data structures

    ! ----------------------------------------------------------------------
    ! Use a block data copy to update the fields with appropriate MPI breakdown
    ! ----------------------------------------------------------------------
    call support_setExport_blockdatacopy(exportState,itemName='Sa_u10m',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_v10m',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_u10n',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_v10n',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_t2m',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_rhoa',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call support_setExport_blockdatacopy(exportState,itemName='Sa_tsfc',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! For sending cpl_scalar grid information to CMEPS
    ! from: https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/cpl/module_cplscalars.F90#L113
    if (flds_scalar_num>0) then
        call support_setExport_scalar(exportState, itemName=trim(flds_scalar_name), rc=rc)
    endif

    if (dodebug_setExport_writenetcdf) then

        call support_setExport_blockdatacopy(exportState,itemName='Sa_slmsk',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_oceanfrac',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_landfrac',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_lakefrac',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_fice',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_psurf',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_q2m',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_tisfc',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_uustar',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_ffmm',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_f10m',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    
        call support_setExport_blockdatacopy(exportState,itemName='Sa_hpbl',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    endif
    
    ! ----------------------------------------------------------------------
    ! Post-process
    ! Air-sea temperature difference needs to be computed (though this
    !     could happen in the mediator instead)
    ! Surface winds need to be populated in the data initialize step
    !     to avoid all-zero values upon initialization due to the
    !     use of the Diag data object in fv3-shield
    ! ----------------------------------------------------------------------
    call support_setExport_postprocess(exportState,itemName='Sa_astdiff',rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! query for clock, importState and exportState
!   call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!   if (CheckError(rc,__LINE__,__FILE__)) return

    call NUOPC_ModelGet(model, driverClock=clock, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, startTime=start_time, currTime=current_time, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if (current_time == start_time .and. dodebug_setExport_datainit_uv) then
        call support_setExport_postprocess(exportState,itemName='Sa_u10n',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_postprocess(exportState,itemName='Sa_v10n',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_postprocess(exportState,itemName='Sa_u10m',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_postprocess(exportState,itemName='Sa_v10m',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! ----------------------------------------------------------------------
    ! Mixed layer model (MLM) fields
    ! ----------------------------------------------------------------------
    if (use_mlm) then
        
        call support_setExport_blockdatacopy(exportState,itemName='Sa_ts_som',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_blockdatacopy(exportState,itemName='Sa_mld',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_blockdatacopy(exportState,itemName='Sa_tml',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_blockdatacopy(exportState,itemName='Sa_huml',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call support_setExport_blockdatacopy(exportState,itemName='Sa_hvml',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

    endif

    ! Add FLUXES - 7/16/24

    if (use_aofluxes) then
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_taux',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_tauy',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_rain',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_lwnet',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_sen',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_evap',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_swndr',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_swndf',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_swvdr',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
        call support_setExport_blockdatacopy(exportState,itemName='Faxa_swvdf',rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return 
        
    endif


    ! DEBUG by looking at different representations of surface winds in the atmosphere model
    ! -----------------------
    !STEVE:NOTE:: https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/atmos_model.F90#L2829
    if (dodebug_setExport) then
        isc = Atm_block%isc
        iec = Atm_block%iec
        jsc = Atm_block%jsc
        jec = Atm_block%jec
!       nk  = Atm_block%npz

        do nb = 1,Atm_block%nblks
            blen = Atm_block%blksz(nb)

        ! ---------------------------------------------------------------------
        ! Write debug output to ESMF PET log files
        ! ---------------------------------------------------------------------
        !STEVE: For debugging, write out to the PET-specific output file what each of the near-surface fields looks like,
        !       and then output what was used in the coupling to be passed to the mediator.
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: in dodebug loop..."), ESMF_LOGMSG_INFO, rc=rc)
            do ix = 1, blen
                ib = Atm_block%index(nb)%ii(ix)
                jb = Atm_block%index(nb)%jj(ix)
                i = ib - isc + 1
                j = jb - jsc + 1 

                write (output_string, "(I0)") ix
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ix = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                write (output_string, "(I0)") i
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: i = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                write (output_string, "(I0)") j
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: j = "//output_string), ESMF_LOGMSG_INFO, rc=rc)

                ! Backup - standard 10m winds
                write (output_string, "(16F8.3)") IPD_Data(nb)%intdiag%u10m(ix)
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%intdiag%u10m(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                write (output_string, "(16F8.3)") IPD_Data(nb)%intdiag%v10m(ix)
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%intdiag%v10m(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)

                ! Default - use 10m neutral winds
                write (output_string, "(16F8.3)") IPD_Data(nb)%intdiag%u10n(ix)
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%intdiag%u10n(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                write (output_string, "(16F8.3)") IPD_Data(nb)%intdiag%v10n(ix)
                call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%intdiag%v10n(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)

                ! Attempt to test neutral winds from Sfcprop computation
                !if (dodebug_setExport_useneutral) then
                !   write (output_string, "(16F8.3)") IPD_Data(nb)%Sfcprop%u10n(ix)
                !   call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%Sfcprop%u10n(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                !   write (output_string, "(16F8.3)") IPD_Data(nb)%Sfcprop%v10n(ix)
                !   call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%Sfcprop%v10n(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                !endif

                ! Attempt to test coupling data object
                if (dodebug_setExport_usecoupling) then
                    write (output_string, "(16F8.3)") IPD_Data(nb)%coupling%u10mi_cpl(ix)
                    call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%coupling%u10mi_cpl(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                    write (output_string, "(16F8.3)") IPD_Data(nb)%coupling%v10mi_cpl(ix)
                    call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: IPD_Data(nb)%coupling%v10mi_cpl(ix) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                endif

                ! Attempt to test surface data object
                !if (dodebug_setExport_usesrf) then
                    write (output_string, "(16F8.3)") Atm(mygrid)%u_srf(i,j)
                    call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: Atm(mygrid)%u_srf(i,j) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                    write (output_string, "(16F8.3)") Atm(mygrid)%v_srf(i,j)
                    call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: Atm(mygrid)%v_srf(i,j) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
                !endif

!               ! Write out what was actually stored to the dataPtr
!               write (output_string, "(16F8.3)") dataPtr_u_r8(i,j)
!               call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: dataPtr_u_r8(i,j) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
!               write (output_string, "(16F8.3)") dataPtr_v_r8(i,j)
!               call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: dataPtr_v_r8(i,j) = "//output_string), ESMF_LOGMSG_INFO, rc=rc)

            enddo
        enddo
    endif


    ! -------------------------------------------------------------------------
    ! Write out the data fields in netcdf form for debugging
    ! -------------------------------------------------------------------------
    if (dodebug_setExport_writenetcdf) then

        call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport::dodebug_setExport_writenetcdf: WRITING..."), ESMF_LOGMSG_INFO, rc=rc)

        ! -----------------------
        ! Write all fields in the export state
        ! -----------------------
        if (dodebug_setExport_multistep) then

            timestamp_str=fms_time_manager_date_to_string(Time_atmos)  !STEVE: using built-in date-to-string converter from fv3
            call NUOPC_Write(exportState, fileNamePrefix="ATM_setExport_exportState_"//trim(timestamp_str)//"_tile*_", overwrite=.true., rc=rc)

            if (CheckError(rc,__LINE__,__FILE__)) return
        else
            call NUOPC_Write(exportState, fileNamePrefix="ATM_setExport_exportState_tile*_", overwrite=.true., rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        endif

    else
        call ESMF_LogWrite(trim("[NO DEBUG] fv3_shield_cap::setExport::dodebug_setExport_writenetcdf: SKIPPING."), ESMF_LOGMSG_INFO, rc=rc)
    endif


  contains

    subroutine support_setExport_blockdatacopy(exportState,itemName,rc)
    ! From Lucas Harris (8/24/2022)
    ! To answer your questions:
    ! 1 & 3) Conveniently, the whole physics state is exposed to the rest of the model. In atmos_model.F90 and atmosphere.F90 this is called IPD_Data% and can be read or modified either in the driver or in the dynamical core. I think what you are most interested in is IPD_Data%Intdiag%u10m and IPD_Data%Intdiag%v10m for winds, and for surface roughness IPD_Data%Sfcprop%zorl and IPD_Data%Sfcprop%zorl%zorlo. The variables %zorl is used pretty heavily in FV3GFS_io.F90 and in GFS_Physics_Driver.F90 but I believe it is not overwritten anywhere; I am not sure whether %zorlo is used for anything after initialization. There is also a pretty extensive coupling datatype defined in GFS_typedefs.F90; it has u10m and v10m, but not any surface roughness information. Using IPD_Data%coupling may be the best choice since I believe the UFS Coupling efforts use it quite heavily:
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/main/GFS_layer/GFS_typedefs.F90
    ! Note that in some cases IPD_Data is decomposed further into blocks for better vectorization; you can see how the blocks are used in atmosphere.F90::atmos_phys_driver_statein().
    ! -----------------------

      use atmos_model_mod,    only: Atm_block, IPD_Data
      use module_block_data,  only: block_atmos_copy, block_data_copy

      type(ESMF_State),intent(inout) :: exportState
      character(*), intent(in) :: itemName
      integer,intent(out) :: rc

      type(ESMF_Field)                            :: field
      real(ESMF_KIND_R8), dimension(:,:), pointer :: dataPtr_r8 !=> null()
      logical                                     :: isFound

      integer :: ni,nj,nk, i,j, ib,jb,nb,ix, blen
      integer :: isc,iec,jsc,jec
      integer :: isd,ied,jsd,jed
      integer :: ierr
      integer                                     :: n,rank
      type(ESMF_TypeKind_Flag)                    :: datatype

      character(128) :: varname = 'noname'

      call ESMF_StateGet(exportState, itemName=trim(itemName), field=field, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      isFound = ESMF_FieldIsCreated(field, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      ! Get some info
      call ESMF_FieldGet(field, name=fieldname, rank=rank, typekind=datatype, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      ! Check that datatypes are what is expected before continuing
      if (datatype == ESMF_TYPEKIND_R8) then
          output_string = "[DEBUG] fv3_shield_cap::setExport:: "//trim(itemName)//" datatype is r8"
          if (dodebug_setExport) call ESMF_LogWrite(trim(output_string), ESMF_LOGMSG_INFO, rc=rc)
      elseif (datatype == ESMF_TYPEKIND_R4) then
          print *, "fv3_shield_cap::setExport:: u10m datatype is r4 but it is expected to be r8. EXITING..."
          stop(1)
      endif

      ! Get the pointer for the specified variable
      ! https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/atmos_model.F90#L2886C72-L2886C72
      call ESMF_FieldGet(field, farrayPtr=dataPtr_r8, localDE=0, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      ! -----------------------
      ! Get array indices
      ! -----------------------
      isc = Atm_block%isc
      iec = Atm_block%iec
      jsc = Atm_block%jsc
      jec = Atm_block%jec
!     nk  = Atm_block%npz

      ! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/c6edc12ff1e14b4a77792145c06bd225ec13eb50/driver/SHiELD/atmosphere.F90#L1417C4-L1424C9
      if (dodebug_setExport) then
          !STEVE:NOTE:: https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/module_fcst_grid_comp.F90#L821
          write (output_string, "(I0)") mygrid
          call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: mygrid = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
          write (output_string, "(I0)") isc
          call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: isc = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
          write (output_string, "(I0)") iec
          call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: iec = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
          write (output_string, "(I0)") jsc
          call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: jsc = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
          write (output_string, "(I0)") jec
          call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: jec = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
      endif

      ! -----------------------
      ! Assign the atmos fields to the ESMF data pointers
      ! -----------------------
      !STEVE:NOTE:: https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/atmos_model.F90#L2829
      do nb = 1,Atm_block%nblks
          blen = Atm_block%blksz(nb)

          if (dodebug_setExport) then
              write (output_string, "(I0)") nb
              call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: nb = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
              write (output_string, "(I0)") blen
              call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: blen = "//output_string), ESMF_LOGMSG_INFO, rc=rc)
          endif

          if (.not. associated(dataPtr_r8)) then
              call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - dataPtr_r8 for "//trim(itemName)//" is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
              stop(1)
          endif


          ! https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/atmos_model.F90#L2920C15-L2920C105
          if     (trim(itemName)=='Sa_u10m') then

              if (.not. use_sfcpropwinds) then
                  if (.not. associated(IPD_Data(nb)%intdiag%u10m)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%u10m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%u10m, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              else
                  if (.not. associated(IPD_Data(nb)%Sfcprop%u10m)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%u10m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%u10m, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              endif

          elseif (trim(itemName)=='Sa_v10m') then

              if (.not. use_sfcpropwinds) then
                  if (.not. associated(IPD_Data(nb)%intdiag%v10m)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%v10m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%v10m, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              else
                  if (.not. associated(IPD_Data(nb)%Sfcprop%v10m)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%v10m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%v10m, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              endif

          elseif (trim(itemName)=='Sa_u10n') then

              if (.not. use_sfcpropwinds) then
                  if (.not. associated(IPD_Data(nb)%intdiag%u10n)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%u10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%u10n, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              else
                  if (.not. associated(IPD_Data(nb)%Sfcprop%u10n)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%u10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%u10n, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              endif

          elseif (trim(itemName)=='Sa_v10n') then

              if (.not. use_sfcpropwinds) then
                  if (.not. associated(IPD_Data(nb)%intdiag%v10n)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%v10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%v10n, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              else
                  if (.not. associated(IPD_Data(nb)%Sfcprop%v10n)) then
                      call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%v10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                      stop(1)
                  endif

                  call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%v10n, Atm_block, nb, rc=rc)
                  if (CheckError(rc,__LINE__,__FILE__)) return

              endif

          elseif (trim(itemName)=='Sa_hpbl') then

              if (.not. associated(IPD_Data(nb)%intdiag%hpbl)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%hpbl is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%hpbl, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_t2m') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%t2m)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%t2m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/0ff11e85972c48057f8d4d5f08322d9fdc0fb7f2/FV3GFS/FV3GFS_io.F90#L6561C34-L6561C41
              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%t2m, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_rhoa') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%rhoa)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%rhoa is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%rhoa, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_tsfc') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%tsfc)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%tsfc is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%tsfc, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_slmsk') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%slmsk)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%slmsk is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%slmsk, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_oceanfrac') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%oceanfrac)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%oceanfrac is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%oceanfrac, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_landfrac') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%landfrac)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%landfrac is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%landfrac, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_lakefrac') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%lakefrac)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%lakefrac is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%lakefrac, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_fice') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%fice)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%fice is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%fice, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_psurf') then

              if (.not. associated(IPD_Data(nb)%intdiag%psurf)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%psurf is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%psurf, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_q2m') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%q2m)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%q2m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%q2m, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_tisfc') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%tisfc)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%tisfc is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%tisfc, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_uustar') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%uustar)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%uustar is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%uustar, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_ffmm') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%ffmm)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%ffmm is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%ffmm, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_f10m') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%f10m)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%f10m is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%f10m, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_ts_som') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%ts_som)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%ts_som is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%ts_som, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_mld') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%mld)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%mld is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%mld, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_tml') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%tml)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%tml is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%tml, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_huml') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%huml)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%huml is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%huml, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Sa_hvml') then

              if (.not. associated(IPD_Data(nb)%Sfcprop%hvml)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%hvml is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%hvml, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

      ! Add FLUXES - 7/16/24
          elseif (trim(itemName)=='Faxa_taux') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1331
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L404

              if (.not. associated(IPD_Data(nb)%intdiag%dusfci)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%dusfci is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%dusfci, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_tauy') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1332
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L405

              if (.not. associated(IPD_Data(nb)%intdiag%dvsfci)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%dvsfci is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%dvsfci, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_rain') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1357C39-L1357C42
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L384
              ! WARNING - this says 'total precip', but CMEPS seems to expect instantaneous precip rate

              if (.not. associated(IPD_Data(nb)%intdiag%pfr)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%pfr is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%pfr Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_lwnet') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1324
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L408
              
              if (.not. associated(IPD_Data(nb)%intdiag%dlwsfci)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%dlwsfci is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%dlwsfci, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_sen') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1333
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L406
              
              if (.not. associated(IPD_Data(nb)%intdiag%dtsfci)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%dtsfci is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%dtsfci, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_evap') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L1336
              
              if (.not. associated(IPD_Data(nb)%intdiag%epi)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%intdiag%epi is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%intdiag%epi, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_swndr') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L341
              
              if (.not. associated(IPD_Data(nb)%coupling%nirbmdi)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%coupling%nirbmdi is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%coupling%nirbmdi, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_swndf') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L342
              
              if (.not. associated(IPD_Data(nb)%coupling%nirdfdi)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%coupling%nirdfdi is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%coupling%nirdfdi, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_swvdr') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L343
              
              if (.not. associated(IPD_Data(nb)%coupling%visbmdi)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%coupling%visbmdi is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%coupling%visbmdi, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return

          elseif (trim(itemName)=='Faxa_swvdf') then
              ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/c8c5d3061b317266cec688acec3e7e33ad8a78a0/GFS_layer/GFS_typedefs.F90#L344
              
              if (.not. associated(IPD_Data(nb)%coupling%visdfdi)) then
                  call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%coupling%visdfdi is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
                  stop(1)
              endif

              call block_data_copy(dataPtr_r8, IPD_Data(nb)%coupling%visdfdi, Atm_block, nb, rc=rc)
              if (CheckError(rc,__LINE__,__FILE__)) return


          else
              print *, "fv3_shield_cap::setExport::support_setExport:: variable name {"//trim(itemName)//"} not recognized. EXITING..."
              stop(1)
          endif

      enddo

    end subroutine support_setExport_blockdatacopy


    subroutine support_setExport_postprocess(exportState,itemName,rc)
      use atmos_model_mod,    only: Atm_block, IPD_Data
      use module_block_data,  only: block_atmos_copy, block_data_copy

      type(ESMF_State),intent(inout)  :: exportState
      character(*), intent(in)        :: itemName
      integer,intent(out)             :: rc

      type(ESMF_Field)                            :: field
      real(ESMF_KIND_R8), dimension(:,:), pointer :: dataPtr_r8, dataPtr_r8_a, dataPtr_r8_b !=> null()
      integer :: nb


      ! Get a pointer to the input object
      call ESMF_StateGet(exportState, itemName=itemName, field=field, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return

      call ESMF_FieldGet(field, farrayPtr=dataPtr_r8, localDE=0, rc=rc)
      if (CheckError(rc,__LINE__,__FILE__)) return


      ! Do variable-specific operations
      if (trim(itemName)=='Sa_astdiff') then

          ! Compute the air-sea temperature difference

          call ESMF_StateGet(exportState, itemName='Sa_t2m', field=field, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
          call ESMF_FieldGet(field, farrayPtr=dataPtr_r8_a, localDE=0, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
          
          call ESMF_StateGet(exportState, itemName='Sa_tsfc', field=field, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
          call ESMF_FieldGet(field, farrayPtr=dataPtr_r8_b, localDE=0, rc=rc)
          if (CheckError(rc,__LINE__,__FILE__)) return
          
          ! Compute the air-sea temperature difference (in K)
          dataPtr_r8 = dataPtr_r8_a - dataPtr_r8_b

      elseif (trim(itemName)=='Sa_u10m' .or. trim(itemName)=='Sa_u10n' .or. trim(itemName)=='Sa_v10m' .or. trim(itemName)=='Sa_v10n') then

          ! If so, replace all-zeros with model's current surface winds estimate
          if (trim(itemName)=='Sa_u10m' .or. trim(itemName)=='Sa_u10n') then
              ! Get the surface winds and use them here as a replacement for this time step (should only be at initialization)
              call get_srfwinds(dataPtr_r8, 'u')

!             do nb = 1,Atm_block%nblks
!                 if (.not. associated(IPD_Data(nb)%statein%ugrs)) then
!                     call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%statein%ugrs is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
!                     stop(1)
!                 endif
!                 call block_data_copy(dataPtr_r8, IPD_Data(nb)%statein%ugrs, Atm_block, nb, rc=rc)
!                 if (CheckError(rc,__LINE__,__FILE__)) return
!             enddo

!             do nb = 1,Atm_block%nblks
!                 if (.not. associated(IPD_Data(nb)%Sfcprop%u10n)) then
!                     call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%u10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
!                     stop(1)
!                 endif
!                 call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%u10n, Atm_block, nb, rc=rc)
!                 if (CheckError(rc,__LINE__,__FILE__)) return
!             enddo

          elseif (trim(itemName)=='Sa_v10m' .or. trim(itemName)=='Sa_v10n') then
              ! Get the surface winds and use them here as a replacement for this time step (should only be at initialization)
              call get_srfwinds(dataPtr_r8, 'v')

!             do nb = 1,Atm_block%nblks
!                 if (.not. associated(IPD_Data(nb)%statein%vgrs)) then
!                     call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%statein%vgrs is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
!                     stop(1)
!                 endif
!                 call block_data_copy(dataPtr_r8, IPD_Data(nb)%statein%vgrs, Atm_block, nb, rc=rc)
!                 if (CheckError(rc,__LINE__,__FILE__)) return
!             enddo

!             do nb = 1,Atm_block%nblks
!                 if (.not. associated(IPD_Data(nb)%Sfcprop%v10n)) then
!                     call ESMF_LogWrite(trim("[DEBUG] fv3_shield_cap::setExport:: ERROR - IPD_Data(nb)%Sfcprop%v10n is not associated. EXITING..."), ESMF_LOGMSG_INFO, rc=rc)
!                     stop(1)
!                 endif
!                 call block_data_copy(dataPtr_r8, IPD_Data(nb)%Sfcprop%v10n, Atm_block, nb, rc=rc)
!                 if (CheckError(rc,__LINE__,__FILE__)) return
!             enddo
          endif


      endif

    end subroutine support_setExport_postprocess


    subroutine get_srfwinds(dataPtr_r8, varchar)
!     use GFS_typedefs,       only: kind_phys
      use atmos_model_mod, only: Atm_block
      use atmosphere_mod, only: Atm, mygrid, get_bottom_wind

      real(ESMF_KIND_R8), dimension(:,:), pointer, intent(inout) :: dataPtr_r8
      character(1), intent(in) :: varchar  ! either 'u' or 'v'
      integer :: isc,iec,jsc,jec
      integer :: i,j
!     real(ESMF_KIND_R4), dimension(:,:), allocatable :: u_bot, v_bot

      logical :: dodebug = .false.

      ! -------------------------------------------------------------------------
      ! Check fv3 dynamical core arrays to make sure they are available
      ! -------------------------------------------------------------------------
      if (.not. allocated(Atm(mygrid)%u_srf) ) then
          print *, "fv3_shield_cap::setExport::get_srfwinds: ERROR - Atm(mygrid)%u_srf is not allocated." 
          print *, "EXITING..."
          stop(1)
      elseif (.not. allocated(Atm(mygrid)%v_srf) ) then
          print *, "fv3_shield_cap::setExport::get_srfwinds: ERROR - Atm(mygrid)%v_srf is not allocated."
          print *, "EXITING..."
          stop(1)
      endif

      ! -------------------------------------------------------------------------
      ! Get array indices
      ! -------------------------------------------------------------------------
      isc = Atm_block%isc
      iec = Atm_block%iec
      jsc = Atm_block%jsc
      jec = Atm_block%jec
  !   nk  = Atm_block%npz

      ! -------------------------------------------------------------------------
      ! Assign the atmos fields to the ESMF data pointers
      ! -----------------------
      !STEVE:NOTE:: https://github.com/NOAA-EMC/fv3atm/blob/deeac5f0acb875f8a282200ebdc69d6157232163/atmos_model.F90#L2829
      ! -----------------------
      ! Set the u/v winds fields 
      ! using Atm(mygrid)%u_srf and Atm(mygrid)%v_srf
      !
      ! Note: the ATM data structure comes from the fv3 dynamical core:
      ! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/efcb02f06dec753518e7fb51c90e300d36a7455f/driver/SHiELD/atmosphere.F90#L1112C13-L1112C45
      ! -------------------------------------------------------------------------
!     if (.not.allocated(u_bot)) allocate(u_bot(isc:iec,jsc:jec))
!     if (.not.allocated(v_bot)) allocate(v_bot(isc:iec,jsc:jec))

!     call get_bottom_wind ( u_bot, v_bot )
      if (varchar=='u') then
!         dataPtr_r8(isc:iec,jsc:jec) = real(u_bot,kind=8)
          do j=jsc,jec
              do i=isc,iec
                  dataPtr_r8(i,j) = real(Atm(mygrid)%u_srf(i,j),kind=8)
              enddo
          enddo
      elseif (varchar=='v') then
!         dataPtr_r8(isc:iec,jsc:jec) = real(v_bot,kind=8)
          do j=jsc,jec
              do i=isc,iec
                  dataPtr_r8(i,j) = real(Atm(mygrid)%v_srf(i,j),kind=8)
              enddo
          enddo
      endif

!     if (allocated(u_bot)) deallocate(u_bot)
!     if (allocated(v_bot)) deallocate(v_bot)

    end subroutine get_srfwinds


    ! Added for interface with CMEPS, 7/1/24
    subroutine support_setExport_scalar(exportState,itemName,rc)
    ! from: https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/cpl/module_cplscalars.F90#L113

!       use atmosphere_mod, only: Atm
        use atmosphere_mod, only: atmosphere_resolution !, atmosphere_domain

        type(ESMF_State),intent(inout)  :: exportState
        character(*), intent(in)        :: itemName
        integer,intent(out)             :: rc
        integer                         :: mlon, mlat

        character(len=*),parameter      :: subname='(fv3_shield_cap:support_setExport_scalar)'
        real(ESMF_KIND_R8)              :: scalardim(3)
        character(ESMF_MAXSTR)          :: msgString

        ! Get local PET
        call ESMF_VMGetCurrent(vm, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        ! from: https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/atmos_model.F90#L572
        call atmosphere_resolution (mlon, mlat, global=.true.)

        scalardim = 0.0
        ! cpl_scalars for export state
        if (dodebug_setExport .and. localPet==0) then
            print *, "ATM::support_setExport_scalar:: mlon = ", mlon
            print *, "ATM::support_setExport_scalar:: mlat = ", mlat
        endif
        scalardim(1) = real(mlon,8)
        scalardim(2) = real(mlat,8)
!       scalardim(3) = 1.0
!       if (.not. Atm%regional)scalardim(3) = 6.0
        scalardim(3) = 6.0

        call ESMF_StateGet(exportState, itemName=trim(flds_scalar_name), field=field, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        if (localPet == 0) then
            call ESMF_FieldGet(field, farrayPtr=dataPtr_r8, rc=rc) !localDE=0, rc=rc)
            if (flds_scalar_index_nx < 0 .or. flds_scalar_index_nx > flds_scalar_num) then
                write (msgString,*) trim(subname)//": ERROR in flds_scalar_index_nx = ", flds_scalar_index_nx
                call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
                rc = ESMF_FAILURE
                stop(1)
            endif
            if (flds_scalar_index_ny < 0 .or. flds_scalar_index_ny > flds_scalar_num) then
                write (msgString,*) trim(subname)//": ERROR in flds_scalar_index_ny = ", flds_scalar_index_ny
                call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
                rc = ESMF_FAILURE
                stop(2)
            endif
            if (flds_scalar_index_ntile < 0 .or. flds_scalar_index_ntile > flds_scalar_num) then
                call ESMF_LogWrite(trim(subname)//": ERROR in flds_scalar_index_ntile", ESMF_LOGMSG_INFO)
                write (msgString,*) trim(subname)//": ERROR in flds_scalar_index_ntile = ", flds_scalar_index_ntile
                call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
                rc = ESMF_FAILURE
                stop(3)
            endif
            ! From: https://github.com/NOAA-EMC/fv3atm/blob/10cd0231282388da16d22a0aae22a1722b773720/module_fcst_grid_comp.F90#L514
            dataPtr_r8(flds_scalar_index_nx,1) = scalardim(1)
            dataPtr_r8(flds_scalar_index_ny,1) = scalardim(2)
            dataPtr_r8(flds_scalar_index_ntile,1) = scalardim(3)
        endif

    end subroutine support_setExport_scalar

  end subroutine setExport


  !-----------------------------------------------------------------------------
  subroutine Finalize(model, rc)
  !-----------------------------------------------------------------------------

    ! input arguments
    type(ESMF_GridComp)        :: model
    integer, intent(out)       :: rc

    ! local variables
    character(len=*),parameter :: subname='fv3_shield_cap::Finalize'
    integer                    :: i, urc
    type(ESMF_VM)              :: vm
    
    logical :: dodebug
    
    dodebug = dodebug_Finalize

    rc = ESMF_SUCCESS
 
    ! TRACE start: Finalize ----------------------------------------------------
    call ESMF_TraceRegionEnter(trim(subname), rc=rc)

    call fms_mpp_set_current_pelist()
    call fms_mpp_clock_end(mainClock)

    termClock = fms_mpp_clock_id( '-Termination' )
    call fms_mpp_clock_begin(termClock)

    call coupler_end

    call fms_mpp_set_current_pelist()
    call fms_mpp_clock_end(termClock)

    call fms_end

    !STEVE: clean up after module-wide allocation in the Realize routine
    !Note run-time error: forrtl: severe (173): A pointer passed to DEALLOCATE points to an object that cannot be deallocated
    ! This is apparently a bug in intel:
    ! https://community.intel.com/t5/Intel-Fortran-Compiler/A-pointer-passed-to-DEALLOCATE-points-to-an-object-that-cannot/td-p/939766
    ! and appears to have been fixed with intel oneAPI HPC toolkit 2023.2:
    ! https://community.intel.com/t5/Intel-Fortran-Compiler/A-pointer-passed-to-DEALLOCATE-points-to-an-object-that-cannot/td-p/939766/page/2
    !STEVE: more work is needed to identify the best way to deal with this.
!   if (associated(farray_u))         deallocate(farray_u)
!   if (associated(farray_v))         deallocate(farray_v)
!   if (associated(farray_t2m))       deallocate(farray_t2m)
!   if (associated(farray_un))        deallocate(farray_un)
!   if (associated(farray_vn))        deallocate(farray_vn)
!   if (associated(farray_rhoa))      deallocate(farray_rhoa)
!   if (associated(farray_tsfc))      deallocate(farray_tsfc)
!   if (associated(farray_slmsk))     deallocate(farray_slmsk)
!   if (associated(farray_oceanfrac)) deallocate(farray_oceanfrac)
!   if (associated(farray_landfrac))  deallocate(farray_landfrac)
!   if (associated(farray_lakefrac))  deallocate(farray_lakefrac)
!   if (associated(farray_fice))      deallocate(farray_fice)
!   if (associated(farray_psurf))     deallocate(farray_psurf)
!   if (associated(farray_q2m))       deallocate(farray_q2m)
!   if (associated(farray_tisfc))     deallocate(farray_tisfc)
!   if (associated(farray_uustar))    deallocate(farray_uustar)
!   if (associated(farray_ffmm))      deallocate(farray_ffmm)
!   if (associated(farray_f10m))      deallocate(farray_f10m)
!   if (associated(farray_hpbl))      deallocate(farray_hpbl)
!   if (associated(farray_tml))       deallocate(farray_tml)
!   if (associated(farray_mld))       deallocate(farray_mld)
!   if (associated(farray_huml))      deallocate(farray_huml)
!   if (associated(farray_hvml))      deallocate(farray_hvml)
!   if (associated(farray_ts_som))    deallocate(farray_ts_som)
!   if (associated(farray_c))         deallocate(farray_c)
!   if (associated(farray_z))         deallocate(farray_z)
!   if (associated(farray_sst))       deallocate(farray_sst)
!   if (associated(farray_astdiff))   deallocate(farray_astdiff)


    call ESMF_TraceRegionExit(trim(subname), rc=rc)
    ! TRACE end: Finalize -----------------------------------------------------


  end subroutine Finalize
  

  !-----------------------------------------------------------------------------
  subroutine coupler_init(model,clock,rc)
  !-----------------------------------------------------------------------------

    ! Added ESMF tools
    type(ESMF_GridComp),intent(in):: model
    type(ESMF_Clock), intent(in)  :: clock  
    integer, intent(out)          :: rc
    type(ESMF_Time)               :: CurrTime, StartTime, StopTime
    type(ESMF_TimeInterval)       :: timeStep, runDuration
    type(esmf_vm)                 :: vm
    integer                       :: localPet, petCount

    ! For checking coupler.res file
    integer                       :: io_unit, calendar_type_res, date_res(6), date_init_res(6), date_end(6)

    character(128) :: pelist_name
    integer :: pelist_commID
    ! For testing:
    integer :: mpi_comm_esmf
    logical :: isroot = .false.
    
    ! STEVE: old fields (11/20/23)
    logical :: use_namelist
    logical, allocatable, dimension(:,:) :: mask
    real,    allocatable, dimension(:,:) :: glon_bnd, glat_bnd

    ! Original coupler_init heading:
    !-----------------------------------------------------------------------
    !   initialize all defined exchange grids and all boundary maps
    !   Note: This is the FMS coupler for the fv3 / land / OML components
    !-----------------------------------------------------------------------
    integer :: total_days, total_seconds, unit, ierr, io
    integer :: n 
    integer :: date_esmf(6), flags
    type (FmsTime_type) :: Run_length
    character(len=9) :: month

    character(len=:), dimension(:), allocatable :: restart_file !< Restart file saved as a string
    integer :: time_stamp_unit !< Unif of the time_stamp file
    integer :: ascii_unit  !< Unit of a dummy ascii file

    rc = ESMF_SUCCESS


    ! TRACE start: vmget -------------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::vmget", rc=rc)

    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call ESMF_vmget(vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if ( fms_mpp_pe() == fms_mpp_root_pe() ) then
        isroot = .true.
    else
        isroot = .false.
    end if
    
    if (dodebug_couplerinit .and. isroot) write(*,*) 'fv3_shield_cap.F90::coupler_init:: petCount=', petCount

    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::vmget", rc=rc)
    ! TRACE end: vmget ---------------------------------------------------------


    !-----------------------------------------------------------------------
    !----- initialization timing identifiers ----
    !-----------------------------------------------------------------------

    !----- read namelist -------
    !----- for backwards compatibilty read from file coupler.nml -----

    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: reading fms_mpp_input_nml_file..."
    read(fms_mpp_input_nml_file, nml=coupler_nml, iostat=io)
    ierr = fms_check_nml_error(io, 'coupler_nml')

    !----- write namelist to logfile -----
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_write_version_number..."
    call fms_write_version_number (version, tag)
    if (isroot) write(fms_mpp_stdlog(),nml=coupler_nml)


    ! TRACE start: Atm -------------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::Atm", rc=rc)

    !----- allocate and set the pelist (to the global pelist) -----
    allocate( Atm%pelist(fms_mpp_npes()) )
    call fms_mpp_get_current_pelist(Atm%pelist)
!   call fms_mpp_get_current_pelist(Atm%pelist, name=pelist_name, commID=pelist_commID)

    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::Atm", rc=rc)
    ! TRACE end: Atm ---------------------------------------------------------


    ! TRACE start: time_steps --------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::time_steps", rc=rc)

    !STEVE: the calendar type should really be set by ESMF, not the input namelist
    !       if anything, we should check here to make sure it is not different...

    !----- read restart file -----
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: read INPUT/coupler.res (if it exists)..."
    if (fms2_io_file_exists('INPUT/coupler.res')) then
        call fms2_io_ascii_read('INPUT/coupler.res', restart_file)
        read(restart_file(1), *) calendar_type
        read(restart_file(2), *) date_init
        read(restart_file(3), *) date_esmf
        deallocate(restart_file)
    else
        force_date_from_namelist = .true.
    endif

    !----- use namelist value (either no restart or override flag on) ---
    if ( force_date_from_namelist ) then

        if ( sum(current_date) <= 0 ) then
            call fms_error_mesg ('program coupler', 'no namelist value for current_date', FATAL)
        else
            date_esmf = current_date
        endif

        !----- override calendar type with namelist value -----

        select case( fms_mpp_uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'GREGORIAN' )
            calendar_type = GREGORIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        case default
            call fms_mpp_error ( FATAL, 'fv3_shield_cap.F90:: coupler_nml entry calendar must '// &
                                    'be one of JULIAN|GREGORIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
        end select

    endif

    call fms_time_manager_set_calendar_type (calendar_type)

    !STEVE: ISSUE: check here to make sure calendar type is not different than what was set by ESMF
    !       then set it here for fv3 using fms_time_manager_set_calendar_type...


    !STEVE: Replace the input date/time with the NUOPC-ESMF date/time here:
    ! date = year, month, day, hour, min, sec
    ! date_init = 
    !
    !STEVE: 
    ! From: subroutine fcst_initialize(fcst_comp, importState, exportState, clock, rc)
    ! In: https://github.com/NOAA-EMC/fv3atm/blob/b9d61f2f0c51f834d5be26d5fbbd581f6d4e63d2/module_fcst_grid_comp.F90

    !-----------------------------------------------------------------------
    !***  set atmos time - using the ESMF clock
    !-----------------------------------------------------------------------
    call ESMF_ClockGet(clock, CurrTime=CurrTime, StartTime=StartTime, StopTime=StopTime, runDuration=runDuration, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    date_init = 0
    call ESMF_TimeGet (StartTime,                      &
                       YY=date_init(1), MM=date_init(2), DD=date_init(3), &
                       H=date_init(4),  M =date_init(5), S =date_init(6), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if ( isroot ) write(*,'(A,6I5)') 'fv3_shield_cap.F90:coupler_init:: ESMF StartTime=', date_init

    date_esmf=0
    call ESMF_TimeGet (CurrTime,                           &
                       YY=date_esmf(1), MM=date_esmf(2), DD=date_esmf(3), &
                       H=date_esmf(4),  M =date_esmf(5), S =date_esmf(6), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if ( isroot ) write(*,'(A,6I5)') 'fv3_shield_cap.F90:coupler_init:: ESMF CurrTime =', date_esmf

    date_end=0
    call ESMF_TimeGet (StopTime,                                       &
                       YY=date_end(1), MM=date_end(2), DD=date_end(3), &
                       H=date_end(4),  M =date_end(5), S =date_end(6), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if ( isroot ) write(*,'(A,6I5)') 'fv3_shield_cap.F90:coupler_init:: ESMF StopTime =', date_end

!STEVE: get interval here to avoid computation below
    call ESMF_TimeIntervalGet (runDuration, s=total_seconds, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if ( isroot ) write(*,'(A,6I5)') 'fv3_shield_cap.F90:coupler_init:: ESMF TimeInterval (in seconds) =', total_seconds

    !------------------------------------------------------------------------
    !STEVE: back to orig coupler_init:
    !------------------------------------------------------------------------

    !----- write current/initial date actually used to logfile file -----
    if ( isroot ) then
        write (fms_mpp_stdlog(),16) date_esmf(1),trim(fms_time_manager_month_name(date_esmf(2))),date_esmf(3:6)
    endif
 16 format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')

    !------ setting affinity ------
    !$  call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads)
    !$  call omp_set_num_threads(atmos_nthreads)

    !-----------------------------------------------------------------------
    !------ initialize diagnostics manager ------
    !! @details Open and read diag_table. Select fields and files for diagnostic output.
    ! https://github.com/NOAA-GFDL/FMS/blob/06b94a7f574e7794684b8584391744ded68e2989/diag_manager/diag_manager.F90#L3770
    ! NOTE: TIME_INIT is an input-only variable. It does this:
    ! "set the diag_init_time if time_init present.  Otherwise, set it to base_time"
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_diag_init..."
    call fms_diag_init (TIME_INIT=date_esmf)

    !----- always override initial/base date with diag_manager value -----
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_diag_get_base_date..."
    call fms_diag_get_base_date ( date_init(1), date_init(2), date_init(3), date_init(4), date_init(5), date_init(6)  )

    !----- use current date if no base date ------
    if ( date_init(1) == 0 ) date_init = date_esmf
!   if ( isroot ) write(*,'(A,6I5)') 'fv3_shield_cap.F90:coupler_init:: Actual StartTime=', date_init

    !----- set initial and current time types ------
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_time_manager_set_date (Time_init)..."
    Time_init  = fms_time_manager_set_date (date_init(1), date_init(2), date_init(3), date_init(4), date_init(5), date_init(6))

    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_time_manager_set_date (Time_atmos)..."
    Time_atmos = fms_time_manager_set_date (date_esmf(1), date_esmf(2), date_esmf(3), date_esmf(4), date_esmf(5), date_esmf(6))

    !-----------------------------------------------------------------------
    !----- compute the ending time (compute days in each month first) -----
    !-----------------------------------------------------------------------

    !STEVE: set the rest
    Time_end   = fms_time_manager_set_date (date_end(1), date_end(2), date_end(3), date_end(4), date_end(5), date_end(6))

    !STEVE: added to use ESMF total run duration
    Run_length = fms_time_manager_set_time (total_seconds,0)

    !Need to pass Time_end into diag_manager for multiple thread case.
    if (dodebug_couplerinit) print *, "fv3_shield_cap::coupler_init:: calling fms_diag_set_time_end..."
    call fms_diag_set_time_end(Time_end)

    !-----------------------------------------------------------------------
    !----- write time stamps (for start time and end time) ------
    !-----------------------------------------------------------------------

    if ( isroot ) then
        open(newunit = time_stamp_unit, file='time_stamp.out', status='replace', form='formatted')

        month = fms_time_manager_month_name(date_esmf(2))
        write (time_stamp_unit,20) date_esmf, month(1:3)

        month = fms_time_manager_month_name(date_end(2))
        write (time_stamp_unit,20) date_end, month(1:3)

        close(time_stamp_unit)
    endif

20  format (6i4,2x,a3)

    !-----------------------------------------------------------------------
    !----- compute the time steps ------
    !-----------------------------------------------------------------------
    
    if (dodebug_couplerinit) print *, "fv3_shield_cap.F90::coupler_init:: calling fms_time_manager_set_time(s)..."

    Time_step_atmos = fms_time_manager_set_time (dt_atmos,0)
    Time_step_ocean = fms_time_manager_set_time (dt_ocean,0)
    num_cpld_calls  = Run_length / Time_step_ocean
    num_atmos_calls = Time_step_ocean / Time_step_atmos

    ! For restarts
    Time_step_restart = fms_time_manager_set_time (restart_secs, restart_days)
    if (restart_start_secs > 0 .or. restart_start_days > 0) then
        Time_start_restart = fms_time_manager_set_time (restart_start_secs, restart_start_days)
        Time_restart = Time_atmos + Time_start_restart
    else
        Time_restart = Time_atmos + Time_step_restart
    end if
    
    ! For restart_aux
    Time_step_restart_aux = fms_time_manager_set_time (restart_secs_aux, restart_days_aux)
    Time_duration_restart_aux = fms_time_manager_set_time (restart_duration_secs_aux, restart_duration_days_aux)
    Time_start_restart_aux = fms_time_manager_set_time (restart_start_secs_aux, restart_start_days_aux)
    Time_restart_aux = Time_atmos + Time_start_restart_aux
    Time_restart_end_aux = Time_restart_aux + Time_duration_restart_aux

    intrm_rst = .false.
    if (restart_days > 0 .or. restart_secs > 0) intrm_rst = .true.
    intrm_rst_1step = .false.
    if (intrm_rst .and. restart_start_secs == 0 .and. restart_start_days == 0) intrm_rst_1step = .true.

    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::time_steps", rc=rc)
    ! TRACE end: time_steps ----------------------------------------------------


    !-----------------------------------------------------------------------
    !------------------- some error checks ---------------------------------
    !-----------------------------------------------------------------------

    ! TRACE start: error_checks ------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::error_checks", rc=rc)

    !----- initial time cannot be greater than current time -------
    if ( Time_init > Time_atmos ) &
         call fms_error_mesg ('program coupler',  'initial time is greater than current time', FATAL)

    !----- make sure run length is a multiple of ocean time step ------
    if ( num_cpld_calls * Time_step_ocean /= Run_length )  &
         call fms_error_mesg ('program coupler', 'run length must be multiple of ocean time step', FATAL)

    ! ---- make sure cpld time step is a multiple of atmos time step ----
    if ( num_atmos_calls * Time_step_atmos /= Time_step_ocean )  &
         call fms_error_mesg ('program coupler', 'atmos time step is not a multiple of the ocean time step', FATAL)

    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::error_checks", rc=rc)
    ! TRACE end: error_checks --------------------------------------------------


    !-----------------------------------------------------------------------
    !------ initialize component models ------
    !-----------------------------------------------------------------------

    ! TRACE start: atmos_init --------------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::atmos_init", rc=rc)

    if (dodebug_couplerinit) print *, "fv3_shield_cap.F90::coupler_init:: calling atmos_model_init..."
    call  atmos_model_init (Atm,  Time_init, Time_atmos, Time_step_atmos, iau_offset) !STEVE: iau_offset needed in more recent versions of fv3 (main post 202204)

    ! STEVE:TRACE
    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::atmos_init", rc=rc)

    if (dodebug_couplerinit) print *, "fv3_shield_cap.F90::coupler_init:: calling fms_memutils_print_memuse_stats..."
    call fms_memutils_print_memuse_stats('after atmos model init')

    !------ initialize data_override -----
    if (.NOT.Atm%bounded_domain) call fms_data_override_init (Atm_domain_in  = Atm%domain)
                             ! Atm_domain_in  = Atm%domain, &
                             ! Ice_domain_in  = Ice%domain, &
                             ! Land_domain_in = Land%domain )
    
    ! TRACE start: RESTART_check -----------------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::coupler_init::RESTART_check", rc=rc)
    ! TRACE end: atmos_init ----------------------------------------------------

    !-----------------------------------------------------------------------
    !---- open and close dummy file in restart dir to check if dir exists --
    if ( isroot ) then !one pe should do this check only in case of a nest
        open(newunit = ascii_unit, file='RESTART/file', status='replace', form='formatted')
        close(ascii_unit,status="delete")
    endif
    !-----------------------------------------------------------------------

    call ESMF_TraceRegionExit("fv3_shield_cap::coupler_init::RESTART_check", rc=rc)
    ! TRACE end: RESTART_check -------------------------------------------------


  end subroutine coupler_init


  !-----------------------------------------------------------------------------
  subroutine coupler_restart(time_stamp)
  !-----------------------------------------------------------------------------
    character(len=32), intent(in), optional :: time_stamp

    integer :: date_esmf(6)
    integer :: restart_unit !< Unit for the coupler restart file
    character(len=128)                      :: file_res

    !----- compute current date ------
    call fms_time_manager_get_date (Time_atmos, date_esmf(1), date_esmf(2), date_esmf(3),  &
                               date_esmf(4), date_esmf(5), date_esmf(6))

    !----- write restart file ------
    file_res = 'RESTART/coupler.res'
    if (present(time_stamp)) then
      file_res = 'RESTART/'//trim(time_stamp)//'.coupler.res'
    endif
    call fms_mpp_set_current_pelist()
    if (fms_mpp_pe() == fms_mpp_root_pe())then
        open(newunit = restart_unit, file=trim(file_res), status='replace', form='formatted')
        write(restart_unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write(restart_unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write(restart_unit, '(6i6,8x,a)' )date_esmf, &
             'Current model time: year, month, day, hour, minute, second'
        close(restart_unit)
    endif

  end subroutine coupler_restart


  !-----------------------------------------------------------------------------
  subroutine coupler_end
  !-----------------------------------------------------------------------------

    integer :: date_fv3(6)
    integer :: restart_unit !< Unit for the coupler restart file

    call atmos_model_end (Atm)

    !----- retrieve current date ------
!   call fms_time_manager_get_date (Atm%Time, date_fv3(1), date_fv3(2), date_fv3(3),  &
!                              date_fv3(4), date_fv3(5), date_fv3(6))

    !----- check time versus expected ending time ----
    if (Atm%Time /= Time_end) call fms_error_mesg ('program coupler',  &
              'final time does not match expected ending time', WARNING)

    !----- write restart file ------
    call coupler_restart()

    !----- write restart file ------
    call fms_mpp_set_current_pelist()
    if (fms_mpp_pe() == fms_mpp_root_pe())then
        open(newunit = restart_unit, file='RESTART/coupler.res', status='replace', form='formatted')
        write(restart_unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write(restart_unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write(restart_unit, '(6i6,8x,a)' )date_fv3, &
             'Current model time: year, month, day, hour, minute, second'
        close(restart_unit)
    endif

    !----- final output of diagnostic fields ----
    call fms_diag_end (Atm%Time)

    !----- to be removed once fms_io is fully deprecated -----
    #ifdef use_deprecated_io
    call fms_io_exit()
    #endif

  end subroutine coupler_end

  
  subroutine grid_create_mosaic(model, grid, rc)
!   use GFS_typedefs,       only: kind_phys
    use IPD_typedefs,       only: kind_phys !, kind_dbl_prec

    type(ESMF_GridComp), intent(in)   :: model
    type(ESMF_Grid), intent(out)      :: grid
    integer, intent(out)              :: rc

    character(len=80)     :: name
    type(ESMF_Info)       :: info
    integer,dimension(2,6):: decomptile                  !define delayout for the 6 cubed-sphere tiles
    integer,dimension(2)  :: regdecomp                   !define delayout for the nest grid
    type(ESMF_Decomp_Flag):: decompflagPTile(2,6)
    type(ESMF_TypeKind_Flag) :: grid_typekind
    character(3)          :: myGridStr
    type(ESMF_DistGrid)   :: distgrid
    type(ESMF_Array)      :: array
    type(ESMF_StaggerLoc) :: staggers(4)
    integer               :: num_staggers, ti

    type(ESMF_Config)      :: config
    integer                :: layout(2)

    ! From: https://github.com/NOAA-EMC/noahmp/blob/569e354ababbde7a7cd68647533769a5c966468d/drivers/ccpp/machine.F#L9
    integer, parameter :: kind_sngl_prec = 4 
    integer, parameter :: kind_dbl_prec = 8

    rc = ESMF_SUCCESS


    ! ----------------------------------------------
    ! read runtime parameters from config
    ! ----------------------------------------------
    call ESMF_GridCompGet(model, config=config, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    !-----------------------------------------
    ! create a Grid object for Fields
    ! https://earthsystemmodeling.org/docs/release/ESMF_8_3_0/ESMC_crefdoc/node5.html
    ! https://earthsystemmodeling.org/docs/release/ESMF_8_3_0/ESMC_crefdoc/node5.html#SECTION05051000000000000000
    !-----------------------------------------

    ! https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html#SECTION050862700000000000000
    ! Create a six-tile ESMF_Grid for a Cubed Sphere grid using regular decomposition. 
    ! Each tile can have different decomposition. The grid coordinates are generated 
    ! based on the algorithm used by GEOS-5, The tile resolution is defined by tileSize.

!   gridIn = ESMF_GridCreateCubedSphereReg(tileSize=6,           &
!              regDecompPTile=1, decompflagPTile=,                        &
!              coordSys, coordTypeKind,                                &
!              deLabelList, staggerLocList,                            &
!              delayout, indexflag, name, transformArgs, rc)

    ! Set up decomposition for each tile:
    ! http://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html#SECTION050862900000000000000
    !
    ! "List of DE counts for each dimension. The second index steps through the tiles. The total 
    ! deCount is determined as the sum over the products of regDecompPTile elements for each tile. 
    ! By default every tile is decomposed in the same way. If the total PET count is less than 6, 
    ! one tile will be assigned to one DE and the DEs will be assigned to PETs sequentially, 
    ! therefore, some PETs may have more than one DEs. If the total PET count is greater than 6, 
    ! the total number of DEs will be multiple of 6 and less than or equal to the total PET count. 
    ! For instance, if the total PET count is 16, the total DE count will be 12 with each tile 
    ! decomposed into 1x2 blocks. The 12 DEs are mapped to the first 12 PETs and the remainding 
    ! 4 PETs have no DEs locally, unless an optional delayout is provided."

    do ti=1,6
        decomptile(:,ti)=(/layout_x,layout_y/) ! Tile i out of 6
        decompflagPTile(:,ti) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)

        ! https://earthsystemmodeling.org/docs/release/ESMF_8_5_0/ESMF_refdoc/node9.html#const:decompflag
        ! ESMF_DECOMP_SYMMEDGEMAX
        ! Decompose elements across the DEs in a symmetric fashion. Start with the maximum number 
        ! of elements at the two edge DEs, and assign a decending number of elements to the 
        ! DEs as the center of the decomposition is approached from both sides.
    enddo

    ! Use the existing grid_spec.nc mosaic file to get the grid definition
    ! https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html#SECTION050862900000000000000
    ! http://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html#SECTION050831200000000000000


    if (kind_phys == kind_dbl_prec) then
        grid_typekind = ESMF_TYPEKIND_R8
        call ESMF_LogWrite("ATM::grid_create_internal:: grid_typekind = ESMF_TYPEKIND_R8", ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    else
        grid_typekind = ESMF_TYPEKIND_R4
        call ESMF_LogWrite("ATM::grid_create_internal:: grid_typekind = ESMF_TYPEKIND_R4", ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    ! ----------------------------------------------
    ! Read "filename" and "tileFilePath" inputs from a config file
    ! (note: this must be changed for different resolution grids)
    ! ----------------------------------------------
    

    ! TRACE start: ESMF_GridCreateMosaic ---------------------------------------
    call ESMF_TraceRegionEnter("fv3_shield_cap::ESMF_GridCreateMosaic", rc=rc)

    staggers(1) = ESMF_STAGGERLOC_CENTER
    staggers(2) = ESMF_STAGGERLOC_CORNER
    num_staggers = 2
    if (use_gridcreate_addedges) then
        ! WARNING: not supported for cubed sphere grids as of ESMF 8.5.0, 11/1/23
        ! [staggerLocList]
        ! The list of stagger locations to fill with coordinates. Only ESMF_STAGGERLOC_CENTER and ESMF_STAGGERLOC_CORNER are supported. If not present, no coordinates will be added or filled.
        ! https://earthsystemmodeling.org/docs/release/ESMF_8_5_0/ESMF_refdoc/node5.html#SECTION050862700000000000000
        staggers(3) = ESMF_STAGGERLOC_EDGE2
        staggers(4) = ESMF_STAGGERLOC_EDGE1
        num_staggers = 4
    endif

    ! See:
    ! https://earthsystemmodeling.org/docs/release/ESMF_8_5_0/ESMF_refdoc/node5.html#SECTION050831200000000000000
    grid = ESMF_GridCreateMosaic(filename=mosaic_filename, &
                                    name='fv3-mosaic-grid',  &
!                                   staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
!                                   staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), & !SGP: 10/27/23 adding edge as well
                                    staggerLocList=staggers(1:num_staggers), & !SGP: 10/31/23 adding edges via config file option
                                    regDecompPTile=decomptile, &
                                    decompflagPTile=decompflagPTile, &
                                    coordTypeKind=grid_typekind, &
                                    tileFilePath=mosaic_tileFilePath, &
                                    rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call ESMF_TraceRegionExit("fv3_shield_cap::ESMF_GridCreateMosaic", rc=rc)
    ! TRACE end: ESMF_GridCreateMosaic -----------------------------------------

    ! Add the sea-land mask to the grid data structure
    ! STEVE: this leverages the EMC cap function addLsmask2grid,
    !        consider developing a more precise method in the future
    !
    ! https://earthsystemmodeling.org/docs/release/ESMF_8_5_0/ESMF_refdoc/node5.html#SECTION050831700000000000000
    ! Like coordinates items are also created on stagger locations. When adding or accessing item data, the 
    ! stagger location is specified to tell the Grid method where in the cell to get the data. The different 
    ! stagger locations may also have slightly different index ranges and sizes.
    call addLsmask2grid(grid, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if (dodebug_gridcreate) then
        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_CENTER distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CORNER, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_CORNER distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_EDGE2 distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_EDGE1 distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    contains

  end subroutine grid_create_mosaic


  !-----------------------------------------------------------------------------
  subroutine grid_create_internal(model, grid, rc)
  !
  ! https://github.com/NOAA-EMC/fv3atm/blob/7a6751d1d438c21e47043e6c34b77d585ecbaef7/module_fcst_grid_comp.F90#L157
  !-----------------------------------------------------------------------------
!   use GFS_typedefs,       only: kind_phys !, kind_sngl_prec
    use IPD_typedefs,       only: kind_phys !, kind_dbl_prec
  
    type(ESMF_GridComp), intent(inout) :: model
    type(ESMF_Grid), intent(out)       :: grid
    integer, intent(out)               :: rc

    character(len=80)     :: name
    type(ESMF_Info)       :: info
    integer,dimension(2,6):: decomptile                  !define delayout for the 6 cubed-sphere tiles
    integer,dimension(2)  :: regdecomp                   !define delayout for the nest grid
    type(ESMF_Decomp_Flag):: decompflagPTile(2,6)
    type(ESMF_TypeKind_Flag) :: grid_typekind
    character(3)          :: myGridStr
    type(ESMF_DistGrid)   :: distgrid
    type(ESMF_Array)      :: array

    type(ESMF_Config)      :: config
    logical                :: do_write_grid
    integer                :: ngrids, mygrid, ti

    ! From: https://github.com/NOAA-EMC/noahmp/blob/569e354ababbde7a7cd68647533769a5c966468d/drivers/ccpp/machine.F#L9
    integer, parameter :: kind_sngl_prec = 4 
    integer, parameter :: kind_dbl_prec = 8

    ! See here:
    ! https://github.com/NOAA-GFDL/SHiELD_physics/blob/2882fdeb429abc2349a8e881803ac67b154532c3/atmos_drivers/coupled/atmos_model.F90#L124C6-L127C124
    ! And compare to here:
    ! https://github.com/NOAA-EMC/fv3atm/blob/1250b416f526f102c021bf1ab62f583bcba2d249/atmos_model.F90#L152
    integer, dimension(2) :: nxy, nxyp
    real(kind=kind_phys), pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
    real(kind=kind_phys), pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
    real(kind=kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
    real(kind=kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(model, name=name, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------------
    ! read runtime parameters from config
    ! ----------------------------------------------
    call ESMF_GridCompGet(model, config=config, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! Check to make sure that necessary inputs were supplied
    if (layout_x < 0 .or. layout_y < 0 .or. tilesize < 0) then
        print *, "fv3_shield_cap::grid_create_internal: layout_x, layout_y, and tilesize (e.g. 96, 192, 384, or 768) must be provided in the fv3_shield_cap.config file."
        print *, "layout_x = ", layout_x
        print *, "layout_y = ", layout_y
        print *, "tilesize = ", tilesize
        print *, "EXITING..."
        rc = ESMF_FAILURE
        return
    endif

    ! Try to resolve the type kind used by the model and use the same
    if (kind_phys == kind_dbl_prec) then
        grid_typekind = ESMF_TYPEKIND_R8
        call ESMF_LogWrite("ATM::grid_create_internal:: grid_typekind = ESMF_TYPEKIND_R8", ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    else
        grid_typekind = ESMF_TYPEKIND_R4
        call ESMF_LogWrite("ATM::grid_create_internal:: grid_typekind = ESMF_TYPEKIND_R4", ESMF_LOGMSG_INFO, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    do ti=1,6
        decomptile(:,ti)=(/layout_x,layout_y/) ! Tile i out of 6
        decompflagPTile(:,ti) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
    enddo

    grid = ESMF_GridCreateCubedSphere(tileSize=tilesize, &
                                      coordSys=ESMF_COORDSYS_SPH_RAD, &
                                      coordTypeKind=grid_typekind, &
                                      regDecompPTile=decomptile, &
                                      decompflagPTile=decompflagPTile, &
                                      name="fv3-internal-grid", &
                                      rc=rc)
                                 !    staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
    if (CheckError(rc,__LINE__,__FILE__)) return


    ! ----------------------------------------------
    ! - Create coordinate arrays around allocations held within Atm data structure and set in Grid
    ! ----------------------------------------------
    if (dodebug_gridcreate) then
        print *, "ATM::grid_create_internal:: shape(Atm%lon) = ", shape(Atm%lon)
        print *, "ATM::grid_create_internal:: shape(Atm%lat) = ", shape(Atm%lat)
    endif
    nxy = shape(Atm%lon)
    allocate (lon(nxy(1),nxy(2)), lat(nxy(1),nxy(2)))
    lon = real(Atm%lon,kind=kind_phys)
    lat = real(Atm%lat,kind=kind_phys)

    if (dodebug_gridcreate) then
        print *, "ATM::grid_create_internal:: shape(Atm%lon_bnd) = ", shape(Atm%lon_bnd)
        print *, "ATM::grid_create_internal:: shape(Atm%lat_bnd) = ", shape(Atm%lat_bnd)
    endif
    nxyp = shape(Atm%lon_bnd)
    allocate (lon_bnd(nxyp(1),nxyp(2)), lat_bnd(nxyp(1),nxyp(2)))
    lon_bnd = real(Atm%lon_bnd,kind=kind_phys)
    lat_bnd = real(Atm%lat_bnd,kind=kind_phys)


    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, distgrid=distgrid, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

!   array = ESMF_ArrayCreate(distgrid, farray=Atm%lon, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    array = ESMF_ArrayCreate(distgrid, farray=lon, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

!   array = ESMF_ArrayCreate(distgrid, farray=Atm%lat, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    array = ESMF_ArrayCreate(distgrid, farray=lat, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if (dodebug_gridcreate) then
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_CENTER distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CORNER, distgrid=distgrid, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

!   array = ESMF_ArrayCreate(distgrid, farray=Atm%lon_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    array = ESMF_ArrayCreate(distgrid, farray=lon_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

!   array = ESMF_ArrayCreate(distgrid, farray=Atm%lat_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    array = ESMF_ArrayCreate(distgrid, farray=lat_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    if (dodebug_gridcreate) then
        print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_CORNER distgrid output:"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_DistGridPrint(distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
    endif

    if (use_gridcreate_addedges) then
        ! WARNING: not supported as of ESMF 8.5.0, 11/1/23
        ! [staggerLocList]
        ! The list of stagger locations to fill with coordinates. Only ESMF_STAGGERLOC_CENTER and ESMF_STAGGERLOC_CORNER are supported. If not present, no coordinates will be added or filled.
        ! https://earthsystemmodeling.org/docs/release/ESMF_8_5_0/ESMF_refdoc/node5.html#SECTION050862700000000000000

        ! TESTING - 10/31/23
        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

!       array = ESMF_ArrayCreate(distgrid, farray=lon, indexflag=ESMF_INDEX_DELOCAL, rc=rc)    ! 10/31/23: ERROR: "Value unrecognized or out of range  - Internal subroutine call returned Error"
        array = ESMF_ArrayCreate(distgrid, farray=lon_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)    ! 10/31/23: ERROR: "Value unrecognized or out of range  - Internal subroutine call returned Error"
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_EDGE2, array=array, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        array = ESMF_ArrayCreate(distgrid, farray=lat_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_EDGE2, array=array, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        if (dodebug_gridcreate) then
            print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_EDGE2 distgrid output:"
            if (CheckError(rc,__LINE__,__FILE__)) return
            call ESMF_DistGridPrint(distgrid, rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        endif

        ! TESTING - 10/31/23
        call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, distgrid=distgrid, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        array = ESMF_ArrayCreate(distgrid, farray=lon_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_EDGE1, array=array, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        array = ESMF_ArrayCreate(distgrid, farray=lat_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return
        call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_EDGE1, array=array, rc=rc)
        if (CheckError(rc,__LINE__,__FILE__)) return

        if (dodebug_gridcreate) then
            print *, "ATM::grid_create_mosaic:: ESMF_STAGGERLOC_EDGE1 distgrid output:"
            if (CheckError(rc,__LINE__,__FILE__)) return
            call ESMF_DistGridPrint(distgrid, rc=rc)
            if (CheckError(rc,__LINE__,__FILE__)) return
        endif

    endif ! (add_edges)

    ! ----------------------------------------------
    !TODO: Consider aligning mask treatment with coordinates... 
    ! ----------------------------------------------
    call addLsmask2grid(grid, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------------
    ! - Add Attributes used by output
    ! ----------------------------------------------
    call ESMF_AttributeAdd(grid, convention="NetCDF", purpose="FV3", attrList=(/"ESMF:gridded_dim_labels"/), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return
    call ESMF_AttributeSet(grid, convention="NetCDF", purpose="FV3", name="ESMF:gridded_dim_labels", valueList=(/"grid_xt", "grid_yt"/), rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    ! ----------------------------------------------
    ! - Write grid to netcdf file
    ! ----------------------------------------------
!   call ESMF_ConfigGetAttribute(config, tilesize, label="do_write_grid:", default=.false., rc=rc)
!   if (CheckError(rc,__LINE__,__FILE__)) return
!
!   if( do_write_grid ) then
!     mygrid = Atm%mygrid
!     write (myGridStr,"(I0)") mygrid
!     call wrt_fcst_grid(grid, "diagnostic_FV3_fcstGrid."//trim(mygridStr)//".nc", regridArea=.TRUE., rc=rc)
!     if (CheckError(rc,__LINE__,__FILE__)) return
!   endif

    ! ----------------------------------------------
    ! - Hold on to the grid by GridComp
    ! ----------------------------------------------
    call ESMF_GridCompSet(model, grid=grid, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

  end subroutine grid_create_internal


  !-----------------------------------------------------------------------------
  subroutine addLsmask2grid(grid, rc)
  !-----------------------------------------------------------------------------
    use atmos_model_mod, only: IPD_Control  ! To get the local(mpi) coordinates
    use atmos_model_mod, only: IPD_Data     ! To get the land-sea mask
    use atmos_model_mod, only: Atm_block    ! To get the block/index numbers
    use GFS_typedefs,    only: GFS_init_type, GFS_kind_phys => kind_phys

    type(ESMF_Grid), intent(inout) :: grid
    integer, optional, intent(out) :: rc

    real(kind=GFS_kind_phys), parameter :: zero    = 0.0_GFS_kind_phys,     &
                                           one     = 1.0_GFS_kind_phys,     &
                                           epsln   = 1.0e-10_GFS_kind_phys, &
                                           zorlmin = 1.0e-7_GFS_kind_phys

!  local vars
    integer isc, iec, jsc, jec
    integer i, j, nb, ix
    integer, allocatable  :: lsmask(:,:)
    integer(kind=ESMF_KIND_I4), pointer  :: maskPtr(:,:)

    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1
    allocate(lsmask(isc:iec,jsc:jec))
!
!$omp parallel do default(shared) private(i,j,nb,ix)
    do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
! use land sea mask: land:1, ocean:0
        lsmask(i,j) = floor(one + epsln - IPD_data(nb)%SfcProp%oceanfrac(ix))
      enddo
    enddo

!
! Get mask
    call ESMF_GridAddItem(grid, itemflag=ESMF_GRIDITEM_MASK, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=maskPtr, rc=rc)
    if (CheckError(rc,__LINE__,__FILE__)) return

!
!$omp parallel do default(shared) private(i,j)
    do j=jsc,jec
      do i=isc,iec
        maskPtr(i-isc+1,j-jsc+1) = lsmask(i,j)
      enddo
    enddo
    deallocate(lsmask)

  end subroutine addLsmask2grid


  !-----------------------------------------------------------------------------
  !> Returns true if ESMF_LogFoundError() determines that rc is an error code. Otherwise false.
  !-----------------------------------------------------------------------------
  logical function CheckError(rc, line, file)
  !-----------------------------------------------------------------------------
    integer, intent(in) :: rc            !< return code to check
    integer, intent(in) :: line          !< Integer source line number
    character(len=*), intent(in) :: file !< User-provided source file name
    integer :: lrc
    CheckError = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
      CheckError = .true.
    endif
  end function CheckError


end module fv3_shield_cap
