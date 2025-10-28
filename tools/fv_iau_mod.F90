!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!-------------------------------------------------------------------------------
!> @brief incremental analysis update module
!> @author Xi.Chen - author of fv_treat_da_inc.F90
!> @author Philip Pegion <philip.pegion@noaa.gov>
!> @date 09/13/2017
!
!>  REVISION HISTORY:
!>  09/13/2017 - Initial Version based on fv_treat_da_inc.F90
!-------------------------------------------------------------------------------

#ifdef OVERLOAD_R4
#define _GET_VAR1 get_var1_real
#else
#define _GET_VAR1 get_var1_double
#endif

module fv_iau_mod

  use fms2_io_mod,         only: file_exists,          &
                                 FmsNetcdfDomainFile_t,&
                                 FmsNetcdfFile_t,      &
                                 close_file,           &
                                 open_file
  use mpp_mod,             only: mpp_error,           &
                                 FATAL,               &
                                 NOTE,                &
                                 mpp_pe,              &
                                 mpp_npes,            &
                                 mpp_get_current_pelist
  use mpp_domains_mod,     only: domain2d, mpp_get_ntile_count

  use constants_mod,       only: pi=>pi_8
  use fv_arrays_mod,       only: fv_atmos_type,       &
                                 fv_grid_type,        &
                                 fv_grid_bounds_type, &
                                 R_GRID
  use fv_mp_mod,           only: is_master
  use sim_nc_mod,          only: open_ncfile,         &
                                 close_ncfile,        &
                                 get_ncdim1,          &
                                 get_var1_double,     &
                                 get_var3_r4,         &
                                 get_var1_real, check_var_exists
#ifdef GFS_PHYS
  use IPD_typedefs,        only: IPD_init_type, IPD_control_type, &
                                 kind_phys
#endif
  use block_control_mod,   only: block_control_type
  use fv_treat_da_inc_mod, only: remap_coef
  use tracer_manager_mod,  only: get_tracer_names,get_tracer_index, get_number_tracers
  use field_manager_mod,   only: MODEL_ATMOS
  implicit none

  private

#ifndef GFS_PHYS
    integer, parameter :: kind_phys = 8
#endif

  real,allocatable::s2c(:,:,:)
!  real:: s2c(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je,4)
!  integer, dimension(Atm(1)%bd%is:Atm(1)%bd%ie,Atm(1)%bd%js:Atm(1)%bd%je):: &
!      id1, id2, jdc
  integer,allocatable,dimension(:,:) :: id1,id2,jdc

  real :: deg2rad,dt,rdt
  integer :: im,jm,km,nfiles,ncid
  integer :: is,  ie,  js,  je
  integer :: npz,ntracers
  character(len=32), allocatable :: tracer_names(:)
  integer, allocatable :: tracer_indicies(:)
  real, allocatable :: ak(:), bk(:)
  type(domain2d) :: fv_domain

  real(kind=4), allocatable:: wk3(:,:,:)
  type iau_internal_data_type
    real,allocatable :: ua_inc(:,:,:)
    real,allocatable :: va_inc(:,:,:)
    real,allocatable :: temp_inc(:,:,:)
    real,allocatable :: delp_inc(:,:,:)
    real,allocatable :: delz_inc(:,:,:)
    real,allocatable :: tracer_inc(:,:,:,:)
  end type iau_internal_data_type
  type iau_external_data_type
    real,allocatable :: ua_inc(:,:,:)
    real,allocatable :: va_inc(:,:,:)
    real,allocatable :: temp_inc(:,:,:)
    real,allocatable :: delp_inc(:,:,:)
    real,allocatable :: delz_inc(:,:,:)
    real,allocatable :: tracer_inc(:,:,:,:)
    logical          :: in_interval = .false.
    logical          :: drymassfixer = .false.
  end type iau_external_data_type
  type iau_state_type
      type(iau_internal_data_type):: inc1
      type(iau_internal_data_type):: inc2
      real(kind=kind_phys)        :: hr1
      real(kind=kind_phys)        :: hr2
      real(kind=kind_phys)        :: wt
      real(kind=kind_phys)        :: wt_normfact
  end type iau_state_type
  type(iau_state_type) :: IAU_state

  public iau_external_data_type

#ifdef GFS_PHYS
  public IAU_initialize, getiauforcing

contains
subroutine IAU_initialize (IPD_Control, IAU_Data,Init_parm, domain_read)
    type (IPD_control_type), intent(in) :: IPD_Control
    type (IAU_external_data_type), intent(inout) :: IAU_Data
    type (IPD_init_type),    intent(in) :: Init_parm
    type(domain2d), intent(in) :: domain_read
    ! local

    character(len=128) :: fname
    real, dimension(:,:,:), allocatable:: u_inc, v_inc
    real, allocatable:: lat(:), lon(:),agrid(:,:,:)
    real(kind=kind_phys) sx,wx,wt,normfact,dtp

    integer:: i, j, k, nstep, kstep
    integer:: i1, i2, j1
    integer:: jbeg, jend

    logical:: found
    integer nfilesall
    integer, allocatable :: idt(:)

    fv_domain = domain_read

    is  = IPD_Control%isc
    ie  = is + IPD_Control%nx-1
    js  = IPD_Control%jsc
    je  = js + IPD_Control%ny-1
    call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    allocate (tracer_names(ntracers))
    allocate (tracer_indicies(ntracers))
    do i = 1, ntracers
       call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
       tracer_indicies(i)  = get_tracer_index(MODEL_ATMOS,tracer_names(i))
    enddo
  
! determine number of increment files to read, and the valid forecast hours

   nfilesall = size(IPD_Control%iau_inc_files)
   nfiles = 0
   if (is_master()) print*,'in iau_init',trim(IPD_Control%iau_inc_files(1)),IPD_Control%iaufhrs(1)
   do k=1,nfilesall
      if (trim(IPD_Control%iau_inc_files(k)) .eq. '' .or. IPD_Control%iaufhrs(k) .lt. 0) exit
      if (is_master()) then
         print *,k,trim(adjustl(IPD_Control%iau_inc_files(k)))
      endif
      nfiles = nfiles + 1
   enddo
   if (is_master()) print *,'nfiles = ',nfiles
   if (nfiles < 1) then
      return
   endif
   if (nfiles > 1) then
      allocate(idt(nfiles-1))
      idt = IPD_Control%iaufhrs(2:nfiles)-IPD_Control%iaufhrs(1:nfiles-1)
      do k=1,nfiles-1
         if (idt(k) .ne. IPD_Control%iaufhrs(2)-IPD_Control%iaufhrs(1)) then
           print *,'forecast intervals in iaufhrs must be constant'
           call mpp_error (FATAL,' forecast intervals in iaufhrs must be constant')
         endif
      enddo
      deallocate(idt)
   endif
   if (is_master()) print *,'iau interval = ',IPD_Control%iau_delthrs,' hours'
   dt = (IPD_Control%iau_delthrs*3600.)
   rdt = 1.0/dt

   npz = IPD_Control%levs

   allocate(ak(npz+1), bk(npz+1))
   ak = Init_parm%ak
   bk = Init_parm%bk

   if (.not.IPD_Control%iau_on_cubed_sphere) then
   !  set up interpolation weights to go from GSI's gaussian grid to cubed sphere
      deg2rad = pi/180.

      fname = 'INPUT/'//trim(IPD_Control%iau_inc_files(1))

      if( file_exists(fname) ) then
         call open_ncfile( fname, ncid )        ! open the file
         call get_ncdim1( ncid, 'lon',   im)
         call get_ncdim1( ncid, 'lat',   jm)
         call get_ncdim1( ncid, 'lev',   km)

         if (km.ne.npz) then
         if (is_master()) print *, 'km = ', km
         call mpp_error(FATAL, &
               '==> Error in IAU_initialize: km is not equal to npz')
         endif

         if(is_master())  write(*,*) fname, ' DA increment dimensions:', im,jm,km

         allocate (  lon(im) )
         allocate (  lat(jm) )

         call _GET_VAR1 (ncid, 'lon', im, lon )
         call _GET_VAR1 (ncid, 'lat', jm, lat )
         call close_ncfile(ncid)

         ! Convert to radians
         do i=1,im
            lon(i) = lon(i) * deg2rad
         enddo
         do j=1,jm
            lat(j) = lat(j) * deg2rad
         enddo

      else
         call mpp_error(FATAL,'==> Error in IAU_initialize: Expected file '&
            //trim(fname)//' for DA increment does not exist')
      endif

      ! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      ! populate agrid
   !    print*,'is,ie,js,je=',is,ie,js,ie
   !    print*,'size xlon=',size(Init_parm%xlon(:,1)),size(Init_parm%xlon(1,:))
   !    print*,'size agrid=',size(agrid(:,1,1)),size(agrid(1,:,1)),size(agrid(1,1,:))
      allocate(s2c(is:ie,js:je,4))
      allocate(id1(is:ie,js:je))
      allocate(id2(is:ie,js:je))
      allocate(jdc(is:ie,js:je))
      allocate(agrid(is:ie,js:je,2))
      do j = 1,size(Init_parm%xlon,2)
         do i = 1,size(Init_parm%xlon,1)
   !         print*,i,j,is-1+j,js-1+j
            agrid(is-1+i,js-1+j,1)=Init_parm%xlon(i,j)
            agrid(is-1+i,js-1+j,2)=Init_parm%xlat(i,j)
         enddo
      enddo
      call remap_coef( is, ie, js, je, is, ie, js, je, &
         im, jm, lon, lat, id1, id2, jdc, s2c, &
         agrid)
      deallocate ( lon, lat,agrid )
   end if

    allocate(IAU_Data%ua_inc(is:ie, js:je, npz))
    allocate(IAU_Data%va_inc(is:ie, js:je, npz))
    allocate(IAU_Data%temp_inc(is:ie, js:je, npz))
    allocate(IAU_Data%delp_inc(is:ie, js:je, npz))
    allocate(IAU_Data%delz_inc(is:ie, js:je, npz))
    allocate(IAU_Data%tracer_inc(is:ie, js:je, npz,ntracers))
! allocate arrays that will hold iau state
    allocate (iau_state%inc1%ua_inc(is:ie, js:je, npz))
    allocate (iau_state%inc1%va_inc(is:ie, js:je, npz))
    allocate (iau_state%inc1%temp_inc (is:ie, js:je, npz))
    allocate (iau_state%inc1%delp_inc (is:ie, js:je, npz))
    allocate (iau_state%inc1%delz_inc (is:ie, js:je, npz))
    allocate (iau_state%inc1%tracer_inc(is:ie, js:je, npz,ntracers))

    iau_state%hr1=IPD_Control%iaufhrs(1)
    iau_state%wt = 1.0 ! IAU increment filter weights (default 1.0)
    iau_state%wt_normfact = 1.0
    if (IPD_Control%iau_filter_increments) then
       ! compute increment filter weights, sum to obtain normalization factor
       dtp=IPD_control%dtp
       nstep = 0.5*IPD_Control%iau_delthrs*3600/dtp
       ! compute normalization factor for filter weights
       normfact = 0.
       do k=1,2*nstep+1
          kstep = k-1-nstep
          sx     = acos(-1.)*kstep/nstep
          wx     = acos(-1.)*kstep/(nstep+1)
          if (kstep .ne. 0) then
             wt = sin(wx)/wx*sin(sx)/sx
          else
             wt = 1.0
          endif
          normfact = normfact + wt
          if (is_master()) print *,'filter wts',k,kstep,wt
       enddo
       iau_state%wt_normfact = (2*nstep+1)/normfact
    endif
    call read_iau_forcing(IPD_Control,iau_state%inc1,'INPUT/'//trim(IPD_Control%iau_inc_files(1)))
    if (nfiles.EQ.1) then  ! only need to get incrments once since constant forcing over window
       call setiauforcing(IPD_Control,IAU_Data,iau_state%wt)
    endif

    if (nfiles.GT.1) then  !have multiple files, but only read in 2 at a time and interpoalte between them
       allocate (iau_state%inc2%ua_inc(is:ie, js:je, npz))
       allocate (iau_state%inc2%va_inc(is:ie, js:je, npz))
       allocate (iau_state%inc2%temp_inc (is:ie, js:je, npz))
       allocate (iau_state%inc2%delp_inc (is:ie, js:je, npz))
       allocate (iau_state%inc2%delz_inc (is:ie, js:je, npz))
       allocate (iau_state%inc2%tracer_inc(is:ie, js:je, npz,ntracers))
       iau_state%hr2=IPD_Control%iaufhrs(2)
       call read_iau_forcing(IPD_Control,iau_state%inc2,'INPUT/'//trim(IPD_Control%iau_inc_files(2)))
    endif
!   print*,'in IAU init',dt,rdt
    IAU_data%drymassfixer = IPD_control%iau_drymassfixer

end subroutine IAU_initialize

subroutine getiauforcing(IPD_Control,IAU_Data)

   implicit none
   type (IPD_control_type), intent(in) :: IPD_Control
   type(IAU_external_data_type),  intent(inout) :: IAU_Data
   real(kind=kind_phys) t1,t2,sx,wx,wt,dtp
   integer n,i,j,k,sphum,kstep,nstep,itnext

   IAU_Data%in_interval=.false.
   if (nfiles.LE.0) then
       return
   endif

   if (nfiles .eq. 1) then
       t1 = IPD_Control%iaufhrs(1)-0.5*IPD_Control%iau_delthrs
       t2 = IPD_Control%iaufhrs(1)+0.5*IPD_Control%iau_delthrs
   else
       t1 = IPD_Control%iaufhrs(1)
       t2 = IPD_Control%iaufhrs(nfiles)
   endif
   if (IPD_Control%iau_filter_increments) then
      ! compute increment filter weight
      ! t1 is beginning of window, t2 end of window
      ! IPD_Control%fhour current time
      ! in window kstep=-nstep,nstep (2*nstep+1 total)
      ! time step IPD_control%dtp
      dtp=IPD_control%dtp
      nstep = 0.5*IPD_Control%iau_delthrs*3600/dtp
      ! compute normalized filter weight
      kstep = ((IPD_Control%fhour-t1) - 0.5*IPD_Control%iau_delthrs)*3600./dtp
      if (IPD_Control%fhour >= t1 .and. IPD_Control%fhour < t2) then
         sx     = acos(-1.)*kstep/nstep
         wx     = acos(-1.)*kstep/(nstep+1)
         if (kstep .ne. 0) then
            wt = (sin(wx)/wx*sin(sx)/sx)
         else
            wt = 1.
         endif
         iau_state%wt = iau_state%wt_normfact*wt
         !if (is_master()) print *,'kstep,t1,t,t2,filter wt=',kstep,t1,IPD_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
      else
         iau_state%wt = 0.
      endif
   endif

   if (nfiles.EQ.1) then
!  on check to see if we are in the IAU window,  no need to update the
!  tendencies since they are fixed over the window
      if ( IPD_Control%fhour < t1 .or. IPD_Control%fhour >= t2 ) then
!         if (is_master()) print *,'no iau forcing',t1,IPD_Control%fhour,t2
         IAU_Data%in_interval=.false.
      else
         if (IPD_Control%iau_filter_increments) call setiauforcing(IPD_Control,IAU_Data,iau_state%wt)
         if (is_master()) print *,'apply iau forcing t1,t,t2,filter wt=',t1,IPD_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
         IAU_Data%in_interval=.true.
      endif
      return
   endif

   if (nfiles > 1) then
      itnext=2
      if (IPD_Control%fhour < t1 .or. IPD_Control%fhour >= t2) then
!         if (is_master()) print *,'no iau forcing',IPD_Control%iaufhrs(1),IPD_Control%fhour,IPD_Control%iaufhrs(nfiles)
         IAU_Data%in_interval=.false.
      else
         if (is_master()) print *,'apply iau forcing t1,t,t2,filter wt=',t1,IPD_Control%fhour,t2,iau_state%wt/iau_state%wt_normfact
         IAU_Data%in_interval=.true.
         do k=nfiles,1,-1
            if (IPD_Control%iaufhrs(k) > IPD_Control%fhour) then
               itnext=k
            endif
         enddo
!         if (is_master()) print *,'itnext=',itnext
         if (IPD_Control%fhour >= iau_state%hr2) then ! need to read in next increment file
            iau_state%hr1=iau_state%hr2
            iau_state%hr2=IPD_Control%iaufhrs(itnext)
            iau_state%inc1=iau_state%inc2
            if (is_master()) print *,'reading next increment file',trim(IPD_Control%iau_inc_files(itnext))
            call read_iau_forcing(IPD_Control,iau_state%inc2,'INPUT/'//trim(IPD_Control%iau_inc_files(itnext)))
         endif
         call updateiauforcing(IPD_Control,IAU_Data,iau_state%wt)
      endif
   endif
   sphum=get_tracer_index(MODEL_ATMOS,'sphum')
 end subroutine getiauforcing

subroutine updateiauforcing(IPD_Control,IAU_Data,wt)

   implicit none
   type (IPD_control_type),        intent(in) :: IPD_Control
   type(IAU_external_data_type),  intent(inout) :: IAU_Data
   real(kind_phys) delt,wt
   integer i,j,k,l

!   if (is_master()) print *,'in updateiauforcing',nfiles,IPD_Control%iaufhrs(1:nfiles)
   delt = (iau_state%hr2-(IPD_Control%fhour))/(IAU_state%hr2-IAU_state%hr1)
   do j = js,je
      do i = is,ie
         do k = 1,npz
            IAU_Data%ua_inc(i,j,k)    =(delt*IAU_state%inc1%ua_inc(i,j,k)    + (1.-delt)* IAU_state%inc2%ua_inc(i,j,k))*rdt*wt
            IAU_Data%va_inc(i,j,k)    =(delt*IAU_state%inc1%va_inc(i,j,k)    + (1.-delt)* IAU_state%inc2%va_inc(i,j,k))*rdt*wt
            IAU_Data%temp_inc(i,j,k)  =(delt*IAU_state%inc1%temp_inc(i,j,k)  + (1.-delt)* IAU_state%inc2%temp_inc(i,j,k))*rdt*wt
            IAU_Data%delp_inc(i,j,k)  =(delt*IAU_state%inc1%delp_inc(i,j,k)  + (1.-delt)* IAU_state%inc2%delp_inc(i,j,k))*rdt*wt
            IAU_Data%delz_inc(i,j,k)  =(delt*IAU_state%inc1%delz_inc(i,j,k)  + (1.-delt)* IAU_state%inc2%delz_inc(i,j,k))*rdt*wt
            do l=1,ntracers
               IAU_Data%tracer_inc(i,j,k,l) =(delt*IAU_state%inc1%tracer_inc(i,j,k,l) + (1.-delt)* IAU_state%inc2%tracer_inc(i,j,k,l))*rdt*wt
            enddo
         enddo
       enddo
   enddo
 end subroutine updateiauforcing


 subroutine setiauforcing(IPD_Control,IAU_Data,wt)

 implicit none
 type (IPD_control_type),        intent(in) :: IPD_Control
 type(IAU_external_data_type),  intent(inout) :: IAU_Data
 real(kind_phys) delt, dt,wt
 integer i,j,k,l,sphum
!  this is only called if using 1 increment file
 if (is_master()) print *,'in setiauforcing',rdt
 do j = js,je
    do i = is,ie
       do k = 1,npz
          IAU_Data%ua_inc(i,j,k)    =wt*IAU_state%inc1%ua_inc(i,j,k)*rdt
          IAU_Data%va_inc(i,j,k)    =wt*IAU_state%inc1%va_inc(i,j,k)*rdt
          IAU_Data%temp_inc(i,j,k)  =wt*IAU_state%inc1%temp_inc(i,j,k)*rdt
          IAU_Data%delp_inc(i,j,k) =wt*IAU_state%inc1%delp_inc(i,j,k)*rdt
          IAU_Data%delz_inc(i,j,k) =wt*IAU_state%inc1%delz_inc(i,j,k)*rdt
          do l = 1,ntracers
             IAU_Data%tracer_inc(i,j,k,l) =wt*IAU_state%inc1%tracer_inc(i,j,k,l)*rdt
          enddo
       enddo
    enddo
 enddo
 sphum=get_tracer_index(MODEL_ATMOS,'sphum')
 end subroutine setiauforcing


subroutine check_ak_bk_consistency(ak_f, bk_f)
   ! Check that ak, bk in IAU file match those in Atm
   implicit none
   real, intent(in), dimension(npz+1) :: ak_f, bk_f
   integer :: k

   if (size(ak_f) /= size(ak)) then
      call mpp_error(FATAL, '==> Error in check_ak_bk_consistency: ak dimension mismatch')
   endif
   if (size(bk_f) /= size(bk)) then
      call mpp_error(FATAL, '==> Error in check_ak_bk_consistency: bk dimension mismatch')
   endif

   do k = 1, size(ak_f)
      if (abs(ak_f(k) - ak(k)) > 1.0e-6) then
         call mpp_error(FATAL, '==> Error in check_ak_bk_consistency: ak values do not match')
      endif
   enddo

   do k = 1, size(bk_f)
      if (abs(bk_f(k) - bk(k)) > 1.0e-6) then
         call mpp_error(FATAL, '==> Error in check_ak_bk_consistency: bk values do not match')
      endif
   enddo
end subroutine check_ak_bk_consistency

subroutine read_iau_forcing_cubed_sphere(increments, fname_time)
   ! Read the IAU on the cubed sphere from the specified file
   ! and populate the increments structure accordingly
   implicit none

   type(iau_internal_data_type), intent(inout):: increments
   character(len=100),  intent(in) :: fname_time

   !locals
   type(FmsNetcdfDomainFile_t) :: FV_tile_IAU, Tra_IAU
   type(FmsNetcdfFile_t)       :: Fv_IAU
   real, allocatable:: ak_f(:), bk_f(:)
   integer :: l, ntiles
   integer, allocatable, dimension(:) :: pes !< Array of the pes in the current pelist
   character(len=6) :: stile_name

   allocate ( ak_f(npz+1) )
   allocate ( bk_f(npz+1) )

   allocate(pes(mpp_npes()))
   fname = 'INPUT/'//trim(fname_time)//'_fv_iau_core.res.nc'
   call mpp_get_current_pelist(pes)
   if (open_file(Fv_IAU,fname,"read", is_restart=.false., pelist=pes)) then
      call read_data(Fv_IAU, 'ak', ak_f(:))
      call read_data(Fv_IAU, 'bk', bk_f(:))
      call close_file(Fv_IAU)
   else
      call mpp_error(NOTE,'==> Warning from read_iau_forcing_cubed_sphere: Expected file '//trim(fname)//' does not exist')
   endif
   deallocate(pes)

   call check_ak_bk_consistency(ak_f, bk_f)
   deallocate(ak_f)
   deallocate(bk_f)
   
   ntiles = mpp_get_ntile_count(fv_domain)
   if(ntiles == 1) then !
   ! In remap_restart theis conditional also checks for .and. .not. Atm(1)%neststruct%nested) then
   ! TO DO: ensure that nested grids are supported in the correct way
      stile_name = '.tile1'
   else
      stile_name = ''
   endif

   fname = 'INPUT/'//trim(fname_time)//'_fv_iau_core.res'//trim(stile_name)//'.nc'
   if (open_file(Fv_tile_IAU, fname, "read", fv_domain, is_restart=.false.)) then
      call read_data(Fv_tile_IAU, 'u_inc', increments%ua_inc)
      call read_data(Fv_tile_IAU, 'v_inc', increments%va_inc)
      call read_data(Fv_tile_IAU, 'T_inc', increments%temp_inc)
      call read_data(Fv_tile_IAU, 'delp_inc', increments%delp_inc)
      call read_data(Fv_tile_IAU, 'delz_inc', increments%delz_inc)
      call close_file(Fv_tile_IAU)
   else
      call mpp_error(NOTE,'==> Warning from read_iau_forcing_cubed_sphere: Expected file '//trim(fname)//' does not exist')
   endif

   fname = 'INPUT/'//trim(fname_time)//'_fv_iau_tracer.res'//trim(stile_name)//'.nc'
   if (open_file(Tra_IAU, fname, "read", fv_domain, is_restart=.false.)) then
      do l=1,ntracers
         call read_data(Tra_IAU, trim(tracer_names(l))//'_inc', increments%tracer_inc(:,:,:,l))
      enddo
      call close_file(Tra_IAU)
   else
      call mpp_error(NOTE,'==> Warning from read_iau_forcing_cubed_sphere: Expected file '//trim(fname)//' does not exist')
   endif

end subroutine read_iau_forcing_cubed_sphere

subroutine read_iau_forcing_gaussian_grid(IPD_Control, increments, fname)
   ! Read the IAU on the Gaussian grid from the specified file
   ! and populate the increments structure accordingly
   implicit none
   
   type (IPD_control_type), intent(in) :: IPD_Control
   type(iau_internal_data_type), intent(inout):: increments
   character(len=*),  intent(in) :: fname

   !locals
   integer:: i, j, l, km
   integer:: j1
   integer:: jbeg, jend
   integer :: is,  ie,  js,  je

   is  = IPD_Control%isc
   ie  = is + IPD_Control%nx-1
   js  = IPD_Control%jsc
   je  = js + IPD_Control%ny-1

   if( file_exists(fname) ) then
      call open_ncfile( fname, ncid )        ! open the file
   else
      call mpp_error(FATAL,'==> Error in read_iau_forcing: Expected file '&
         //trim(fname)//' for DA increment does not exist')
   endif

   km = IPD_Control%levs

   ! Find bounding latitudes:
   jbeg = jm-1;         jend = 2
   do j=js,je
   do i=is,ie
      j1 = jdc(i,j)
      jbeg = min(jbeg, j1)
      jend = max(jend, j1+1)
   enddo
   enddo

   allocate ( wk3(1:im,jbeg:jend, 1:km) )
   ! read in 1 time level
   call interp_inc('T_inc',increments%temp_inc(:,:,:),jbeg,jend)
   call interp_inc('delp_inc',increments%delp_inc(:,:,:),jbeg,jend)
   call interp_inc('delz_inc',increments%delz_inc(:,:,:),jbeg,jend)
   call interp_inc('u_inc',increments%ua_inc(:,:,:),jbeg,jend)   ! can these be treated as scalars?
   call interp_inc('v_inc',increments%va_inc(:,:,:),jbeg,jend)
   do l=1,ntracers
      call interp_inc(trim(tracer_names(l))//'_inc',increments%tracer_inc(:,:,:,l),jbeg,jend)
   enddo
   call close_ncfile(ncid)
   deallocate (wk3)

end subroutine read_iau_forcing_gaussian_grid

subroutine read_iau_forcing(IPD_Control,increments,fname)
    type (IPD_control_type), intent(in) :: IPD_Control
    type(iau_internal_data_type), intent(inout):: increments
    character(len=*),  intent(in) :: fname

    if (IPD_control%iau_on_cubed_sphere) then
      call read_iau_forcing_cubed_sphere(increments, fname)
    else
      call read_iau_forcing_gaussian_grid(IPD_Control, increments, fname)
    end if

end subroutine read_iau_forcing

subroutine interp_inc(field_name,var,jbeg,jend)
! interpolate increment from GSI gaussian grid to cubed sphere
! everying is on the A-grid, earth relative
 character(len=*), intent(in) :: field_name
 real, dimension(is:ie,js:je,1:km), intent(inout) :: var
 integer, intent(in) :: jbeg,jend
 integer:: i1, i2, j1, k,j,i,ierr
 call check_var_exists(ncid, field_name, ierr)
 if (ierr == 0) then
    call get_var3_r4( ncid, field_name, 1,im, jbeg,jend, 1,km, wk3 )
 else
    if (is_master()) print *,'warning: no increment for ',trim(field_name),' found, assuming zero'
    wk3 = 0.
 endif
 do k=1,km
    do j=js,je
       do i=is,ie
          i1 = id1(i,j)
          i2 = id2(i,j)
          j1 = jdc(i,j)
          var(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k)+&
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
       enddo
    enddo
 enddo
end subroutine interp_inc

#endif

end module fv_iau_mod


