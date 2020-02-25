!--------------- LANL CICE NUOPC CAP -----------------
! This is the LANL CICE model cap component that's NUOPC compiant.
!
! Author:  Fei.Liu@gmail.com
!
! 5/10/13
! This is now acting as a cap/connector between NUOPC driver and LANL CICE code.
!
! 8/27/18: Denise Worthen (denise.worthen@noaa.gov)
!  * `GridAttachArea` - when set to "true", this option indicates that CICE grid attaches cell area
!   using internal values computed in CICE. The default value is "false", so grid cell area will
!   be computed in ESMF.
! 9/24/18: Denise Worthen (denise.worthen@noaa.gov)
! * corrected calculation of sea surface slope gradients across PET boundaries; corrected vector halo
!   filling for move to CICE U-grid
! 10/3/18: Denise Worthen (denise.worthen@noaa.gov)
! * calculation of slope of sea surface set non-op; slopes are obtained by import of fields from the
!   ocean component
! 06/19/19: Denise Worthen (denise.worthen@noaa.gov)
! * removal of unused code and variables; basic tidying up of code in prep for unification with NCAR

module cice_cap_mod

  use ice_blocks, only: nx_block, ny_block, nblocks_tot, block, get_block, &
                        get_block_parameter
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_domain, only: nblocks, blocks_ice, halo_info, distrb_info
  use ice_distribution, only: ice_distributiongetblockloc
  use ice_constants, only: Tffresh, rad_to_deg
  use ice_calendar,  only: dt
  use ice_flux
  use ice_grid, only: TLAT, TLON, ULAT, ULON, hm, tarea, ANGLET, ANGLE, &
                      dxt, dyt, t2ugrid_vector
  use ice_constants, only: field_loc_center, field_loc_NEcorner, field_type_scalar, field_type_vector
  use ice_boundary, only: ice_HaloUpdate

  use ice_state
  use CICE_RunMod
  use CICE_InitMod
  use CICE_FinalMod 

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  implicit none
  private
  public SetServices

  type cice_internalstate_type
  end type

  type cice_internalstate_wrapper
    type(cice_internalstate_type), pointer :: ptr
  end type

  integer   :: import_slice = 0
  integer   :: export_slice = 0

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
#ifdef CMEPS
    real(ESMF_KIND_R8), dimension(:,:), pointer :: farrayPtr
#else
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
#endif
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToIce_num = 0
  type (fld_list_type) :: fldsToIce(fldsMax)
  integer :: fldsFrIce_num = 0
  type (fld_list_type) :: fldsFrIce(fldsMax)

  integer :: lsize    ! local number of gridcells for coupling
  character(len=256) :: tmpstr
  character(len=2048):: info
  logical :: isPresent
  integer :: dbrc     ! temporary debug rc value

  type(ESMF_Grid), save :: ice_grid_i
  logical :: write_diagnostics = .false.
  logical :: overwrite_timeslice = .false.
  logical :: profile_memory = .false.
  logical :: grid_attach_area = .false.
  ! local helper flag for halo debugging
  logical :: HaloDebug = .false.

#ifdef CMEPS
  character(ESMF_MAXSTR) :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0
  integer                :: flds_scalar_index_nx = 0
  integer                :: flds_scalar_index_ny = 0
  integer                :: flds_scalar_index_nextsw_cday = 0
#endif

  contains
  !-----------------------------------------------------------------------
  !------------------- CICE code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(cice_cap:SetServices)'

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    ! No need to change clock settings
    call ESMF_MethodAdd(gcomp, label=model_label_SetClock, &
      userRoutine=SetClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_MethodAdd(gcomp, label=model_label_Advance, &
      userRoutine=ModelAdvance_slow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=cice_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_FieldsSetup()

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    character(len=10)     :: value
    type(ESMF_VM)         :: vm
    integer               :: lpet

    character(240)        :: msgString
    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
     write(msgString,'(a12,i8)')'CICE lpet = ',lpet
     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write_diagnostics=(trim(value)=="true")
    write(msgString,'(A,l6)')'CICE_CAP: Dumpfields = ',write_diagnostics
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_AttributeGet(gcomp, name="OverwriteSlice", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    overwrite_timeslice=(trim(value)/="false")
    write(msgString,'(A,l6)')'CICE_CAP: OverwriteSlice = ',overwrite_timeslice
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    profile_memory=(trim(value)/="false")
    write(msgString,'(A,l6)')'CICE_CAP: Profile_memory = ',profile_memory
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_AttributeGet(gcomp, name="GridAttachArea", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    grid_attach_area=(trim(value)=="true")
    write(msgString,'(A,l6)')'CICE_CAP: GridAttachArea = ',grid_attach_area
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeP0
  
  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    integer                                :: mpi_comm
#ifdef CMEPS
    character(ESMF_MAXSTR)                 :: cvalue
    character(ESMF_MAXSTR)                 :: logmsg
    logical                                :: isPresent, isSet
#endif
    character(len=*),parameter  :: subname='(cice_cap:InitializeAdvertise)'

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Initialize(mpi_comm)
    
#ifdef CMEPS
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       call fld_list_add(fldsToIce_num, fldsToIce, trim(flds_scalar_name), "will provide")
       call fld_list_add(fldsFrIce_num, fldsFrIce, trim(flds_scalar_name), "will provide")
    else
       call ESMF_LogWrite(trim(subname)//' Need to set attribute ScalarFieldName', ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    else
       call ESMF_LogWrite(trim(subname)//' Need to set attribute ScalarFieldCount', ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    else
       call ESMF_LogWrite(trim(subname)//' Need to set attribute ScalarFieldIdxGridNX', ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    else
       call ESMF_LogWrite(trim(subname)//' Need to set attribute ScalarFieldIdxGridNY', ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    else
       call ESMF_LogWrite(trim(subname)//' Need to set attribute ScalarFieldIdxNextSwCday', ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    endif
#endif

    call CICE_AdvertiseFields(importState, fldsToIce_num, fldsToIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_AdvertiseFields(exportState, fldsFrIce_num, fldsFrIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) trim(subname),' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut
    type(ESMF_DistGrid)                    :: distgrid
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    integer                                :: npet
    integer                                :: i,j,iblk, n, i1,j1, DE
    integer                                :: ilo,ihi,jlo,jhi
    integer                                :: ig,jg,cnt
    integer                                :: peID,locID
    integer, pointer                       :: indexList(:)
    integer, pointer                       :: deLabelList(:)
    integer, pointer                       :: deBlockList(:,:,:)
    integer, pointer                       :: petMap(:)
    integer, pointer                       :: i_glob(:),j_glob(:)
    integer                                :: lbnd(2),ubnd(2)
    type(block)                            :: this_block
    type(ESMF_DELayout)                    :: delayout
    real(ESMF_KIND_R8), pointer            :: tarray(:,:)     
    real(ESMF_KIND_R8), pointer :: coordXcenter(:,:)
    real(ESMF_KIND_R8), pointer :: coordYcenter(:,:)
    real(ESMF_KIND_R8), pointer :: coordXcorner(:,:)
    real(ESMF_KIND_R8), pointer :: coordYcorner(:,:)
    integer(ESMF_KIND_I4), pointer :: gridmask(:,:)
    real(ESMF_KIND_R8), pointer :: gridarea(:,:)
    character(len=*),parameter  :: subname='(cice_cap:InitializeRealize)'

    rc = ESMF_SUCCESS

    ! We can check if npet is 4 or some other value to make sure
    ! CICE is configured to run on the correct number of processors.

    ! create a Grid object for Fields
    ! we are going to create a single tile displaced pole grid from a gridspec
    ! file. We also use the exact decomposition in CICE so that the Fields
    ! created can wrap on the data pointers in internal part of CICE

    write(tmpstr,'(a,2i8)') trim(subname)//' ice nx,ny = ',nx_global,ny_global
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    allocate(deBlockList(2,2,nblocks_tot))
    allocate(petMap(nblocks_tot))
    allocate(deLabelList(nblocks_tot))

    write(tmpstr,'(a,1i8)') trim(subname)//' nblocks = ',nblocks_tot
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do n = 1, nblocks_tot
       deLabelList(n) = n
       call get_block_parameter(n,ilo=ilo,ihi=ihi,jlo=jlo,jhi=jhi, &
          i_glob=i_glob,j_glob=j_glob)
       deBlockList(1,1,n) = i_glob(ilo)
       deBlockList(1,2,n) = i_glob(ihi)
       deBlockList(2,1,n) = j_glob(jlo)
       deBlockList(2,2,n) = j_glob(jhi)
       call ice_distributionGetBlockLoc(distrb_info,n,peID,locID)
       petMap(n) = peID - 1
       write(tmpstr,'(a,2i8)') trim(subname)//' IDs  = ',n,peID
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(tmpstr,'(a,3i8)') trim(subname)//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(tmpstr,'(a,3i8)') trim(subname)//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       write(tmpstr,'(a,3i8)') trim(subname)//' petMap = ',n,petMap(n),nblocks_tot
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo

    delayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    allocate(connectionList(2))
    ! bipolar boundary condition at top row: nyg
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nx_global+1, 2*ny_global+1/), &
      orientationVector=(/-1, -2/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! periodic boundary condition along first dimension
    call ESMF_DistGridConnectionSet(connectionList(2), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nx_global, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nx_global,ny_global/), &
        deBlockList=deBlockList, &
        delayout=delayout, &
        connectionList=connectionList, &
        rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    deallocate(deLabelList)
    deallocate(deBlockList)
    deallocate(petMap)
    deallocate(connectionList)

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(indexList(cnt))
    write(tmpstr,'(a,i8)') trim(subname)//' distgrid cnt= ',cnt
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(tmpstr,'(a,4i8)') trim(subname)//' distgrid list= ',indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    deallocate(IndexList)

    gridIn = ESMF_GridCreate(distgrid=distgrid, &
       coordSys = ESMF_COORDSYS_SPH_DEG, &
       gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
       rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Attach area to the Grid optionally. By default the cell areas are computed.
    if(grid_attach_area) then
      call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
         staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    do iblk = 1,nblocks
       DE = iblk-1
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           computationalLBound=lbnd, computationalUBound=ubnd, &
           farrayPtr=coordXcenter, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=coordYcenter, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       write(tmpstr,'(a,5i8)') trim(subname)//' iblk center bnds ',iblk,lbnd,ubnd
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       if (lbnd(1) /= 1 .or. lbnd(2) /= 1 .or. ubnd(1) /= ihi-ilo+1 .or. ubnd(2) /= jhi-jlo+1) then
          write(tmpstr,'(a,5i8)') trim(subname)//' iblk bnds ERROR '
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
          rc = ESMF_FAILURE
          return
       endif

       call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=gridmask, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if(grid_attach_area) then
       call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, localDE=DE, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=gridarea, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       do j1 = lbnd(2),ubnd(2)
       do i1 = lbnd(1),ubnd(1)
          i = i1 + ilo - lbnd(1)
          j = j1 + jlo - lbnd(2)
          gridarea(i1,j1) = tarea(i,j,iblk)
       enddo
       enddo
       write(tmpstr,'(a,5i8)') trim(subname)//' setting ESMF_GRIDITEM_AREA using tarea '
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
      endif

       do j1 = lbnd(2),ubnd(2)
       do i1 = lbnd(1),ubnd(1)
          i = i1 + ilo - lbnd(1)
          j = j1 + jlo - lbnd(2)
          coordXcenter(i1,j1) = TLON(i,j,iblk) * rad_to_deg
          coordYcenter(i1,j1) = TLAT(i,j,iblk) * rad_to_deg
          gridmask(i1,j1) = nint(hm(i,j,iblk))
       enddo
       enddo

       call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           computationalLBound=lbnd, computationalUBound=ubnd, &
           farrayPtr=coordXcorner, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           farrayPtr=coordYcorner, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       write(tmpstr,'(a,5i8)') trim(subname)//' iblk corner bnds ',iblk,lbnd,ubnd
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       ! ULON and ULAT are upper right hand corner from TLON and TLAT
       ! corners in ESMF need to be defined lon lower left corner from center
       ! ULON and ULAT have ghost cells, leverage that to fill corner arrays
       do j1 = lbnd(2),ubnd(2)
       do i1 = lbnd(1),ubnd(1)
          i = i1 + ilo - lbnd(1)
          j = j1 + jlo - lbnd(2)
          coordXcorner(i1,j1) = ULON(i-1,j-1,iblk) * rad_to_deg
          coordYcorner(i1,j1) = ULAT(i-1,j-1,iblk) * rad_to_deg
       enddo
       enddo

    enddo

    call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') trim(subname)//' gridIn center1 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') trim(subname)//' gridIn center2 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') trim(subname)//' gridIn corner1 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') trim(subname)//' gridIn corner2 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    gridOut = gridIn ! for now out same as in
    ice_grid_i = gridIn

    call CICE_RealizeFields(importState, gridIn , fldsToIce_num, fldsToIce, "Ice import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_RealizeFields(exportState, gridOut, fldsFrIce_num, fldsFrIce, "Ice export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

! Have to be careful with reset since states are pointing directly into cice arrays
!    call state_reset(ImportState, value=-99._ESMF_KIND_R8, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
#ifndef CMEPS
    call state_reset(ExportState, value=-99._ESMF_KIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

#ifdef CMEPS
    call ice_export(exportState)

    call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    write(tmpstr,'(a,3i8)') trim(subname)//' nx_block, ny_block, nblocks = ',nx_block,ny_block,nblocks
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(info,*) trim(subname),' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)

  end subroutine InitializeRealize
  
  !-----------------------------------------------------------------------------

  ! CICE model uses same clock as parent gridComp
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep, timestep
    character(len=*),parameter  :: subname='(cice_cap:SetClock)'

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! tcraig: dt is the cice thermodynamic timestep in seconds
    call ESMF_TimeIntervalSet(timestep, s=nint(dt), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    call ESMF_TimeIntervalSet(stabilityTimeStep, s=nint(dt), rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine SetClock

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance_slow(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
#ifdef CMEPS
    type(ESMF_Field)                       :: lfield
#else
    type(ESMF_Field)                       :: lfield,lfield2d
#endif
    type(ESMF_Grid)                        :: grid
#ifndef CMEPS
    real(ESMF_KIND_R8), pointer            :: fldptr(:,:,:)
    real(ESMF_KIND_R8), pointer            :: fldptr2d(:,:)
#endif
    type(block)                            :: this_block
    character(len=64)                      :: fldname
    integer                                :: i,j,iblk,n,i1,i2,j1,j2
    integer                                :: ilo,ihi,jlo,jhi
    real(ESMF_KIND_R8)                     :: ue, vn, ui, vj
    real(ESMF_KIND_R8)                     :: sigma_r, sigma_l, sigma_c
    type(ESMF_StateItem_Flag)              :: itemType
    ! imports
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer :: dataPtr_mdlwfx(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swir(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swif(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lprec(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fprec(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sst(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sss(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssz(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssm(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncz(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncm(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fmpot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rhoabot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_Tbot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_pbot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_qbot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_zlvl(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ubot(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vbot(:,:)
#else
    real(ESMF_KIND_R8), pointer :: dataPtr_mdlwfx(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swir(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swif(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lprec(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fprec(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sst(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sss(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssz(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssm(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncz(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncm(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fmpot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_rhoabot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_Tbot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_pbot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_qbot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_zlvl(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ubot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vbot(:,:,:)
#endif
    ! exports
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ifrac(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_itemp(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairxT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairyT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnxT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnyT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthru(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flwout(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsens(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flat(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fhocn(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fresh(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsalt(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vice(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vsno(:,:)
#else
    real(ESMF_KIND_R8), pointer :: dataPtr_mask(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ifrac(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_itemp(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairxT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairyT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnxT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnyT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthru(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flwout(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsens(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flat(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fhocn(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fresh(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsalt(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vice(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vsno(:,:,:)
#endif
    character(240)              :: msgString
    character(len=*),parameter  :: subname='(cice_cap:ModelAdvance_slow)'

    rc = ESMF_SUCCESS
    if(profile_memory) call ESMF_VMLogMemInfo("Entering CICE Model_ADVANCE: ")
    write(info,*) trim(subname),' --- run phase 1 called --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing CICE from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", &
      unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  if(write_diagnostics) then
    import_slice = import_slice + 1

    call state_diagnose(importState, 'cice_import', rc)
    do i = 1,fldsToice_num
      fldname = fldsToice(i)%shortname
      call ESMF_StateGet(importState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(importState, itemName=trim(fldname), field=lfield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(lfield,grid=grid,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

#ifdef CMEPS
        call ESMF_FieldWrite(lfield, fileName='field_ice_import_'//trim(fldname)//'.nc', &
          timeslice=import_slice, overwrite=overwrite_timeslice, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#else
        ! create a copy of the 3d data in lfield but in a 2d array, lfield2d
        lfield2d = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=trim(fldname), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_FieldGet(lfield  , farrayPtr=fldptr  , rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(lfield2d, farrayPtr=fldptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        fldptr2d(:,:) = fldptr(:,:,1)

        call ESMF_FieldWrite(lfield2d, fileName='field_ice_import_'//trim(fldname)//'.nc', &
          timeslice=import_slice, overwrite=overwrite_timeslice, rc=rc) 
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_FieldDestroy(lfield2d, noGarbage=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#endif
      endif
    enddo
  endif  ! write_diagnostics 

    call State_getFldPtr(importState,'inst_temp_height_lowest',dataPtr_Tbot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'inst_spec_humid_height_lowest',dataPtr_qbot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'inst_zonal_wind_height_lowest',dataPtr_ubot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'inst_merid_wind_height_lowest',dataPtr_vbot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'inst_pres_height_lowest',dataPtr_pbot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_down_lw_flx',dataPtr_mdlwfx,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_down_sw_vis_dir_flx',dataPtr_swvr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_down_sw_vis_dif_flx',dataPtr_swvf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_down_sw_ir_dir_flx',dataPtr_swir,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_down_sw_ir_dif_flx',dataPtr_swif,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_prec_rate',dataPtr_lprec,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_fprec_rate',dataPtr_fprec,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sea_surface_temperature',dataPtr_sst,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'s_surf',dataPtr_sss,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sea_surface_slope_zonal',dataPtr_sssz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sea_surface_slope_merid',dataPtr_sssm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ocn_current_zonal',dataPtr_ocncz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ocn_current_merid',dataPtr_ocncm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'freezing_melting_potential',dataPtr_fmpot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'inst_height_lowest',dataPtr_zlvl,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'air_density_height_lowest',dataPtr_rhoabot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    do iblk = 1,nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
       do i = ilo,ihi
          ! i1=1:120,j1=1:540
          i1 = i - ilo + 1
          j1 = j - jlo + 1
#ifdef CMEPS
          rhoa   (i,j,iblk) = dataPtr_rhoabot(i1,j1)  ! import directly from mediator  
          if(dataPtr_pbot(i1,j1) .gt. 0.0) &
          potT   (i,j,iblk) = dataPtr_Tbot   (i1,j1) * (100000./dataPtr_pbot(i1,j1))**0.286 ! Potential temperature (K)
          Tair   (i,j,iblk) = dataPtr_Tbot   (i1,j1)  ! near surface temp, maybe lowest level (K)
          Qa     (i,j,iblk) = dataPtr_qbot   (i1,j1)  ! near surface humidity, maybe lowest level (kg/kg)
          zlvl   (i,j,iblk) = dataPtr_zlvl   (i1,j1)  ! height of the lowest level (m) 
          flw    (i,j,iblk) = dataPtr_mdlwfx (i1,j1)  ! downwelling longwave flux
          swvdr  (i,j,iblk) = dataPtr_swvr   (i1,j1)  ! downwelling shortwave flux, vis dir
          swvdf  (i,j,iblk) = dataPtr_swvf   (i1,j1)  ! downwelling shortwave flux, vis dif
          swidr  (i,j,iblk) = dataPtr_swir   (i1,j1)  ! downwelling shortwave flux, nir dir
          swidf  (i,j,iblk) = dataPtr_swif   (i1,j1)  ! downwelling shortwave flux, nir dif
          fsw(i,j,iblk) = swvdr(i,j,iblk)+swvdf(i,j,iblk)+swidr(i,j,iblk)+swidf(i,j,iblk)
          frain  (i,j,iblk) = dataPtr_lprec  (i1,j1)  ! flux of rain (liquid only)
          fsnow  (i,j,iblk) = dataPtr_fprec  (i1,j1)  ! flux of frozen precip
          sss    (i,j,iblk) = dataPtr_sss    (i1,j1)  ! sea surface salinity (maybe for mushy layer)
          sst    (i,j,iblk) = dataPtr_sst    (i1,j1) - 273.15  ! sea surface temp (may not be needed?)
          frzmlt (i,j,iblk) = dataPtr_fmpot  (i1,j1)
!          ! --- rotate these vectors from east/north to i/j ---
          uocn   (i,j,iblk) = dataPtr_ocncz  (i1,j1)
          vocn   (i,j,iblk) = dataPtr_ocncm  (i1,j1)
          uatm   (i,j,iblk) = dataPtr_ubot   (i1,j1)
          vatm   (i,j,iblk) = dataPtr_vbot   (i1,j1)
          ss_tltx(i,j,iblk) = dataPtr_sssz   (i1,j1)
          ss_tlty(i,j,iblk) = dataPtr_sssm   (i1,j1)
#else
          rhoa   (i,j,iblk) = dataPtr_rhoabot(i1,j1,iblk)  ! import directly from mediator  
          if(dataPtr_pbot(i1,j1,iblk) .gt. 0.0) &
          potT   (i,j,iblk) = dataPtr_Tbot   (i1,j1,iblk) * (100000./dataPtr_pbot(i1,j1,iblk))**0.286 ! Potential temperature (K)
          Tair   (i,j,iblk) = dataPtr_Tbot   (i1,j1,iblk)  ! near surface temp, maybe lowest level (K)
          Qa     (i,j,iblk) = dataPtr_qbot   (i1,j1,iblk)  ! near surface humidity, maybe lowest level (kg/kg)
          zlvl   (i,j,iblk) = dataPtr_zlvl   (i1,j1,iblk)  ! height of the lowest level (m) 
          flw    (i,j,iblk) = dataPtr_mdlwfx (i1,j1,iblk)  ! downwelling longwave flux
          swvdr  (i,j,iblk) = dataPtr_swvr   (i1,j1,iblk)  ! downwelling shortwave flux, vis dir
          swvdf  (i,j,iblk) = dataPtr_swvf   (i1,j1,iblk)  ! downwelling shortwave flux, vis dif
          swidr  (i,j,iblk) = dataPtr_swir   (i1,j1,iblk)  ! downwelling shortwave flux, nir dir
          swidf  (i,j,iblk) = dataPtr_swif   (i1,j1,iblk)  ! downwelling shortwave flux, nir dif
          fsw(i,j,iblk) = swvdr(i,j,iblk)+swvdf(i,j,iblk)+swidr(i,j,iblk)+swidf(i,j,iblk)
          frain  (i,j,iblk) = dataPtr_lprec  (i1,j1,iblk)  ! flux of rain (liquid only)
          fsnow  (i,j,iblk) = dataPtr_fprec  (i1,j1,iblk)  ! flux of frozen precip
          sss    (i,j,iblk) = dataPtr_sss    (i1,j1,iblk)  ! sea surface salinity (maybe for mushy layer)
          sst    (i,j,iblk) = dataPtr_sst    (i1,j1,iblk) - 273.15  ! sea surface temp (may not be needed?)
          frzmlt (i,j,iblk) = dataPtr_fmpot  (i1,j1,iblk)
!          ! --- rotate these vectors from east/north to i/j ---
          uocn   (i,j,iblk) = dataPtr_ocncz  (i1,j1,iblk)
          vocn   (i,j,iblk) = dataPtr_ocncm  (i1,j1,iblk)
          uatm   (i,j,iblk) = dataPtr_ubot   (i1,j1,iblk)
          vatm   (i,j,iblk) = dataPtr_vbot   (i1,j1,iblk)
          ss_tltx(i,j,iblk) = dataPtr_sssz   (i1,j1,iblk)
          ss_tlty(i,j,iblk) = dataPtr_sssm   (i1,j1,iblk)
#endif
       enddo
       enddo
    enddo

    do iblk = 1, nblocks

       do j = 1,ny_block
          do i = 1,nx_block
          ! ocean
          ue = uocn(i,j,iblk)
          vn = vocn(i,j,iblk)
          uocn(i,j,iblk) =  ue*cos(ANGLET(i,j,iblk)) + vn*sin(ANGLET(i,j,iblk))  ! x ocean current
          vocn(i,j,iblk) = -ue*sin(ANGLET(i,j,iblk)) + vn*cos(ANGLET(i,j,iblk))  ! y ocean current

          ue = ss_tltx(i,j,iblk)
          vn = ss_tlty(i,j,iblk)
          ss_tltx(i,j,iblk) =  ue*cos(ANGLET(i,j,iblk)) + vn*sin(ANGLET(i,j,iblk))  ! x ocean surface slope
          ss_tlty(i,j,iblk) = -ue*sin(ANGLET(i,j,iblk)) + vn*cos(ANGLET(i,j,iblk))  ! y ocean surface slope

          ! atm
          ue = uatm(i,j,iblk)
          vn = vatm(i,j,iblk)
          uatm(i,j,iblk) =  ue*cos(ANGLET(i,j,iblk)) + vn*sin(ANGLET(i,j,iblk))  ! x wind
          vatm(i,j,iblk) = -ue*sin(ANGLET(i,j,iblk)) + vn*cos(ANGLET(i,j,iblk))  ! y wind
          wind(i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
          enddo !i
       enddo !j
    enddo !iblk

    ! Interpolate ocean dynamics variables from T-cell centers to
    ! U-cell centers.
    ! Atmosphere variables are needed in T cell centers in
    ! subroutine stability and are interpolated to the U grid
    ! later as necessary.
    ! note: t2ugrid call includes HaloUpdate at location center
    ! followed by call to move the vectors
    ! halos are returned as zeros
       call t2ugrid_vector(uocn)
       call t2ugrid_vector(vocn)
       call t2ugrid_vector(ss_tltx)
       call t2ugrid_vector(ss_tlty)

    write(info,*) trim(subname),' --- run phase 2 called --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)
    if(profile_memory) call ESMF_VMLogMemInfo("Before CICE_Run")
    call CICE_Run
    if(profile_memory) call ESMF_VMLogMemInfo("Afterr CICE_Run")
    write(info,*) trim(subname),' --- run phase 3 called --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

    !---- local modifications to coupling fields -----

    call State_getFldPtr(exportState,'ice_mask',dataPtr_mask,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'ice_fraction',dataPtr_ifrac,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'sea_ice_surface_temperature',dataPtr_itemp,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_vis_dir_albedo',dataPtr_alvdr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_vis_dif_albedo',dataPtr_alvdf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_ir_dir_albedo',dataPtr_alidr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_ir_dif_albedo',dataPtr_alidf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_air_ice_zonal',dataPtr_strairxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_air_ice_merid',dataPtr_strairyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_ocn_ice_zonal',dataPtr_strocnxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_ocn_ice_merid',dataPtr_strocnyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'net_heat_flx_to_ocn',dataPtr_fhocn,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_fresh_water_to_ocean_rate',dataPtr_fresh,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_salt_rate',dataPtr_fsalt,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_ice_volume',dataPtr_vice,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_snow_volume',dataPtr_vsno,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn',dataPtr_fswthru,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_vis_dir_flx',dataPtr_fswthruvdr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_vis_dif_flx',dataPtr_fswthruvdf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_ir_dir_flx',dataPtr_fswthruidr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_ir_dif_flx',dataPtr_fswthruidf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_up_lw_flx_ice',dataPtr_flwout,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sensi_heat_flx_atm_into_ice',dataPtr_fsens,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_laten_heat_flx_atm_into_ice',dataPtr_flat,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_evap_rate_atm_into_ice',dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    dataPtr_ifrac = 0._ESMF_KIND_R8
    dataPtr_itemp = 0._ESMF_KIND_R8
    dataPtr_mask = 0._ESMF_KIND_R8
    do iblk = 1,nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
       do i = ilo,ihi
          i1 = i - ilo + 1
          j1 = j - jlo + 1
#ifdef CMEPS
          if (hm(i,j,iblk) > 0.5) dataPtr_mask(i1,j1) = 1._ESMF_KIND_R8
          dataPtr_ifrac   (i1,j1) = aice(i,j,iblk)   ! ice fraction (0-1)
          if (dataPtr_ifrac(i1,j1) > 0._ESMF_KIND_R8) &
             dataPtr_itemp   (i1,j1) = Tffresh + trcr(i,j,1,iblk)  ! surface temperature of ice covered portion (degK)
          dataPtr_alvdr   (i1,j1) = alvdr(i,j,iblk)  ! albedo vis dir
          dataPtr_alidr   (i1,j1) = alidr(i,j,iblk)  ! albedo nir dir
          dataPtr_alvdf   (i1,j1) = alvdf(i,j,iblk)  ! albedo vis dif
          dataPtr_alidf   (i1,j1) = alidf(i,j,iblk)  ! albedo nir dif
          dataPtr_fswthru (i1,j1) = fswthru(i,j,iblk) ! flux of shortwave through ice to ocean
          dataPtr_fswthruvdr (i1,j1) = fswthruvdr(i,j,iblk) ! flux of vis dir shortwave through ice to ocean
          dataPtr_fswthruvdf (i1,j1) = fswthruvdf(i,j,iblk) ! flux of vis dif shortwave through ice to ocean
          dataPtr_fswthruidr (i1,j1) = fswthruidr(i,j,iblk) ! flux of ir dir shortwave through ice to ocean
          dataPtr_fswthruidf (i1,j1) = fswthruidf(i,j,iblk) ! flux of ir dif shortwave through ice to ocean
          dataPtr_flwout  (i1,j1) = flwout(i,j,iblk)   ! longwave outgoing (upward), average over ice fraction only
          dataPtr_fsens   (i1,j1) =  fsens(i,j,iblk)   ! sensible
          dataPtr_flat    (i1,j1) =   flat(i,j,iblk)   ! latent
          dataPtr_evap    (i1,j1) =   evap(i,j,iblk)   ! evaporation (not ~latent, need separate field)
          dataPtr_fhocn   (i1,j1) =  fhocn(i,j,iblk)   ! heat exchange with ocean 
          dataPtr_fresh   (i1,j1) =  fresh(i,j,iblk)   ! fresh water to ocean
          dataPtr_fsalt   (i1,j1) =  fsalt(i,j,iblk)   ! salt to ocean
          dataPtr_vice    (i1,j1) =   vice(i,j,iblk)   ! sea ice volume
          dataPtr_vsno    (i1,j1) =   vsno(i,j,iblk)   ! snow volume
          ! --- rotate these vectors from i/j to east/north ---
          ui = strairxT(i,j,iblk)
          vj = strairyT(i,j,iblk)
          dataPtr_strairxT(i1,j1) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! air ice stress
          dataPtr_strairyT(i1,j1) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! air ice stress
          ui = -strocnxT(i,j,iblk)
          vj = -strocnyT(i,j,iblk)
          dataPtr_strocnxT(i1,j1) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! ice ocean stress
          dataPtr_strocnyT(i1,j1) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! ice ocean stress
#else
          if (hm(i,j,iblk) > 0.5) dataPtr_mask(i1,j1,iblk) = 1._ESMF_KIND_R8
          dataPtr_ifrac   (i1,j1,iblk) = aice(i,j,iblk)   ! ice fraction (0-1)
          if (dataPtr_ifrac(i1,j1,iblk) > 0._ESMF_KIND_R8) &
             dataPtr_itemp   (i1,j1,iblk) = Tffresh + trcr(i,j,1,iblk)  ! surface temperature of ice covered portion (degK)
          dataPtr_alvdr   (i1,j1,iblk) = alvdr(i,j,iblk)  ! albedo vis dir
          dataPtr_alidr   (i1,j1,iblk) = alidr(i,j,iblk)  ! albedo nir dir
          dataPtr_alvdf   (i1,j1,iblk) = alvdf(i,j,iblk)  ! albedo vis dif
          dataPtr_alidf   (i1,j1,iblk) = alidf(i,j,iblk)  ! albedo nir dif
          dataPtr_fswthru (i1,j1,iblk) = fswthru(i,j,iblk) ! flux of shortwave through ice to ocean
          dataPtr_fswthruvdr (i1,j1,iblk) = fswthruvdr(i,j,iblk) ! flux of vis dir shortwave through ice to ocean
          dataPtr_fswthruvdf (i1,j1,iblk) = fswthruvdf(i,j,iblk) ! flux of vis dif shortwave through ice to ocean
          dataPtr_fswthruidr (i1,j1,iblk) = fswthruidr(i,j,iblk) ! flux of ir dir shortwave through ice to ocean
          dataPtr_fswthruidf (i1,j1,iblk) = fswthruidf(i,j,iblk) ! flux of ir dif shortwave through ice to ocean
          dataPtr_flwout  (i1,j1,iblk) = flwout(i,j,iblk)   ! longwave outgoing (upward), average over ice fraction only
          dataPtr_fsens   (i1,j1,iblk) =  fsens(i,j,iblk)   ! sensible
          dataPtr_flat    (i1,j1,iblk) =   flat(i,j,iblk)   ! latent
          dataPtr_evap    (i1,j1,iblk) =   evap(i,j,iblk)   ! evaporation (not ~latent, need separate field)
          dataPtr_fhocn   (i1,j1,iblk) =  fhocn(i,j,iblk)   ! heat exchange with ocean 
          dataPtr_fresh   (i1,j1,iblk) =  fresh(i,j,iblk)   ! fresh water to ocean
          dataPtr_fsalt   (i1,j1,iblk) =  fsalt(i,j,iblk)   ! salt to ocean
          dataPtr_vice    (i1,j1,iblk) =   vice(i,j,iblk)   ! sea ice volume
          dataPtr_vsno    (i1,j1,iblk) =   vsno(i,j,iblk)   ! snow volume
          ! --- rotate these vectors from i/j to east/north ---
          ui = strairxT(i,j,iblk)
          vj = strairyT(i,j,iblk)
          dataPtr_strairxT(i1,j1,iblk) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! air ice stress
          dataPtr_strairyT(i1,j1,iblk) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! air ice stress
          ui = -strocnxT(i,j,iblk)
          vj = -strocnyT(i,j,iblk)
          dataPtr_strocnxT(i1,j1,iblk) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! ice ocean stress
          dataPtr_strocnyT(i1,j1,iblk) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! ice ocean stress
#endif
       enddo
       enddo
    enddo

    !-------------------------------------------------

  if(write_diagnostics) then
    call state_diagnose(exportState, 'cice_export', rc)

    export_slice = export_slice + 1

    do i = 1,fldsFrIce_num
      fldname = fldsFrIce(i)%shortname
      call ESMF_StateGet(exportState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(exportState, itemName=trim(fldname), field=lfield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(lfield,grid=grid,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

#ifdef CMEPS
        call ESMF_FieldWrite(lfield, fileName='field_ice_export_'//trim(fldname)//'.nc', &
          timeslice=export_slice, overwrite=overwrite_timeslice, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#else
        ! create a copy of the 3d data in lfield but in a 2d array, lfield2d
        lfield2d = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=trim(fldname), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_FieldGet(lfield  , farrayPtr=fldptr  , rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(lfield2d, farrayPtr=fldptr2d, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        fldptr2d(:,:) = fldptr(:,:,1)

        call ESMF_FieldWrite(lfield2d, fileName='field_ice_export_'//trim(fldname)//'.nc', &
          timeslice=export_slice, overwrite=overwrite_timeslice,rc=rc) 
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMF_FieldDestroy(lfield2d, noGarbage=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#endif
      endif
    enddo
  endif  ! write_diagnostics 
    write(info,*) trim(subname),' --- run phase 4 called --- ',rc
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

! Dump out all the cice internal fields to cross-examine with those connected with mediator
! This will help to determine roughly which fields can be hooked into cice

   call dumpCICEInternal(ice_grid_i, import_slice, "inst_zonal_wind_height10m", "will provide", strax)
   call dumpCICEInternal(ice_grid_i, import_slice, "inst_merid_wind_height10m", "will provide", stray)
   call dumpCICEInternal(ice_grid_i, import_slice, "inst_pres_height_surface" , "will provide", zlvl)
   call dumpCICEInternal(ice_grid_i, import_slice, "ocn_current_zonal", "will provide", uocn)
   call dumpCICEInternal(ice_grid_i, import_slice, "ocn_current_merid", "will provide", vocn)
   call dumpCICEInternal(ice_grid_i, import_slice, "sea_surface_slope_zonal", "will provide", ss_tltx)
   call dumpCICEInternal(ice_grid_i, import_slice, "sea_surface_slope_merid", "will provide", ss_tlty)
   call dumpCICEInternal(ice_grid_i, import_slice, "sea_surface_salinity", "will provide", sss)
   call dumpCICEInternal(ice_grid_i, import_slice, "sea_surface_temperature", "will provide", sst)

!--------- export fields from Sea Ice -------------

   call dumpCICEInternal(ice_grid_i, export_slice, "ice_fraction"                    , "will provide", aice)
   call dumpCICEInternal(ice_grid_i, export_slice, "stress_on_air_ice_zonal"         , "will provide", strairxT)
   call dumpCICEInternal(ice_grid_i, export_slice, "stress_on_air_ice_merid"         , "will provide", strairyT)
   call dumpCICEInternal(ice_grid_i, export_slice, "stress_on_ocn_ice_zonal"         , "will provide", strocnxT)
   call dumpCICEInternal(ice_grid_i, export_slice, "stress_on_ocn_ice_merid"         , "will provide", strocnyT)
   call dumpCICEInternal(ice_grid_i, export_slice, "mean_sw_pen_to_ocn"              , "will provide", fswthru)
   if(profile_memory) call ESMF_VMLogMemInfo("Leaving CICE Model_ADVANCE: ")

  end subroutine ModelAdvance_slow 


#ifdef CMEPS
  subroutine ice_export(exportState)
    type(ESMF_State) :: exportState

    ! local variables
    integer :: rc
    type(block)                            :: this_block
    character(len=64)                      :: fldname
    integer                                :: i,j,iblk,n,i1,i2,j1,j2
    integer                                :: ilo,ihi,jlo,jhi
    real(ESMF_KIND_R8)                     :: ue, vn, ui, vj
   
    real(ESMF_KIND_R8), pointer :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ifrac(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_itemp(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairxT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairyT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnxT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnyT(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthru(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruvdf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidr(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthruidf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flwout(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsens(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flat(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fhocn(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fresh(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsalt(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vice(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_vsno(:,:)

    call State_getFldPtr(exportState,'ice_mask',dataPtr_mask,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'ice_fraction',dataPtr_ifrac,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'sea_ice_surface_temperature',dataPtr_itemp,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_vis_dir_albedo',dataPtr_alvdr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_vis_dif_albedo',dataPtr_alvdf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_ir_dir_albedo',dataPtr_alidr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'inst_ice_ir_dif_albedo',dataPtr_alidf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_air_ice_zonal',dataPtr_strairxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_air_ice_merid',dataPtr_strairyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_ocn_ice_zonal',dataPtr_strocnxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'stress_on_ocn_ice_merid',dataPtr_strocnyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'net_heat_flx_to_ocn',dataPtr_fhocn,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_fresh_water_to_ocean_rate',dataPtr_fresh,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_salt_rate',dataPtr_fsalt,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_ice_volume',dataPtr_vice,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_snow_volume',dataPtr_vsno,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn',dataPtr_fswthru,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_vis_dir_flx',dataPtr_fswthruvdr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_vis_dif_flx',dataPtr_fswthruvdf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_ir_dir_flx',dataPtr_fswthruidr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sw_pen_to_ocn_ir_dif_flx',dataPtr_fswthruidf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_up_lw_flx_ice',dataPtr_flwout,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_sensi_heat_flx_atm_into_ice',dataPtr_fsens,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_laten_heat_flx_atm_into_ice',dataPtr_flat,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mean_evap_rate_atm_into_ice',dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    dataPtr_ifrac = 0._ESMF_KIND_R8
    dataPtr_itemp = 0._ESMF_KIND_R8
    dataPtr_mask = 0._ESMF_KIND_R8
    do iblk = 1,nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
       do i = ilo,ihi
          i1 = i - ilo + 1
          j1 = j - jlo + 1
#ifdef CMEPS
          if (hm(i,j,iblk) > 0.5) dataPtr_mask(i1,j1) = 1._ESMF_KIND_R8
          dataPtr_ifrac   (i1,j1) = aice(i,j,iblk)   ! ice fraction (0-1)
          if (dataPtr_ifrac(i1,j1) > 0._ESMF_KIND_R8) &
             dataPtr_itemp   (i1,j1) = Tffresh + trcr(i,j,1,iblk)  ! surface temperature of ice covered portion (degK)
          dataPtr_alvdr   (i1,j1) = alvdr(i,j,iblk)  ! albedo vis dir
          dataPtr_alidr   (i1,j1) = alidr(i,j,iblk)  ! albedo nir dir
          dataPtr_alvdf   (i1,j1) = alvdf(i,j,iblk)  ! albedo vis dif
          dataPtr_alidf   (i1,j1) = alidf(i,j,iblk)  ! albedo nir dif
          dataPtr_fswthru (i1,j1) = fswthru(i,j,iblk) ! flux of shortwave through ice to ocean
          dataPtr_fswthruvdr (i1,j1) = fswthruvdr(i,j,iblk) ! flux of vis dir shortwave through ice to ocean
          dataPtr_fswthruvdf (i1,j1) = fswthruvdf(i,j,iblk) ! flux of vis dif shortwave through ice to ocean
          dataPtr_fswthruidr (i1,j1) = fswthruidr(i,j,iblk) ! flux of ir dir shortwave through ice to ocean
          dataPtr_fswthruidf (i1,j1) = fswthruidf(i,j,iblk) ! flux of ir dif shortwave through ice to ocean
          dataPtr_flwout  (i1,j1) = flwout(i,j,iblk)   ! longwave outgoing (upward), average over ice fraction only
          dataPtr_fsens   (i1,j1) =  fsens(i,j,iblk)   ! sensible
          dataPtr_flat    (i1,j1) =   flat(i,j,iblk)   ! latent
          dataPtr_evap    (i1,j1) =   evap(i,j,iblk)   ! evaporation (not ~latent, need separate field)
          dataPtr_fhocn   (i1,j1) =  fhocn(i,j,iblk)   ! heat exchange with ocean 
          dataPtr_fresh   (i1,j1) =  fresh(i,j,iblk)   ! fresh water to ocean
          dataPtr_fsalt   (i1,j1) =  fsalt(i,j,iblk)   ! salt to ocean
          dataPtr_vice    (i1,j1) =   vice(i,j,iblk)   ! sea ice volume
          dataPtr_vsno    (i1,j1) =   vsno(i,j,iblk)   ! snow volume
          ! --- rotate these vectors from i/j to east/north ---
          ui = strairxT(i,j,iblk)
          vj = strairyT(i,j,iblk)
          dataPtr_strairxT(i1,j1) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! air ice stress
          dataPtr_strairyT(i1,j1) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! air ice stress
          ui = -strocnxT(i,j,iblk)
          vj = -strocnyT(i,j,iblk)
          dataPtr_strocnxT(i1,j1) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! ice ocean stress
          dataPtr_strocnyT(i1,j1) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! ice ocean stress
#else
          if (hm(i,j,iblk) > 0.5) dataPtr_mask(i1,j1,iblk) = 1._ESMF_KIND_R8
          dataPtr_ifrac   (i1,j1,iblk) = aice(i,j,iblk)   ! ice fraction (0-1)
          if (dataPtr_ifrac(i1,j1,iblk) > 0._ESMF_KIND_R8) &
             dataPtr_itemp   (i1,j1,iblk) = Tffresh + trcr(i,j,1,iblk)  ! surface temperature of ice covered portion (degK)
          dataPtr_alvdr   (i1,j1,iblk) = alvdr(i,j,iblk)  ! albedo vis dir
          dataPtr_alidr   (i1,j1,iblk) = alidr(i,j,iblk)  ! albedo nir dir
          dataPtr_alvdf   (i1,j1,iblk) = alvdf(i,j,iblk)  ! albedo vis dif
          dataPtr_alidf   (i1,j1,iblk) = alidf(i,j,iblk)  ! albedo nir dif
          dataPtr_fswthru (i1,j1,iblk) = fswthru(i,j,iblk) ! flux of shortwave through ice to ocean
          dataPtr_fswthruvdr (i1,j1,iblk) = fswthruvdr(i,j,iblk) ! flux of vis dir shortwave through ice to ocean
          dataPtr_fswthruvdf (i1,j1,iblk) = fswthruvdf(i,j,iblk) ! flux of vis dif shortwave through ice to ocean
          dataPtr_fswthruidr (i1,j1,iblk) = fswthruidr(i,j,iblk) ! flux of ir dir shortwave through ice to ocean
          dataPtr_fswthruidf (i1,j1,iblk) = fswthruidf(i,j,iblk) ! flux of ir dif shortwave through ice to ocean
          dataPtr_flwout  (i1,j1,iblk) = flwout(i,j,iblk)   ! longwave outgoing (upward), average over ice fraction only
          dataPtr_fsens   (i1,j1,iblk) =  fsens(i,j,iblk)   ! sensible
          dataPtr_flat    (i1,j1,iblk) =   flat(i,j,iblk)   ! latent
          dataPtr_evap    (i1,j1,iblk) =   evap(i,j,iblk)   ! evaporation (not ~latent, need separate field)
          dataPtr_fhocn   (i1,j1,iblk) =  fhocn(i,j,iblk)   ! heat exchange with ocean 
          dataPtr_fresh   (i1,j1,iblk) =  fresh(i,j,iblk)   ! fresh water to ocean
          dataPtr_fsalt   (i1,j1,iblk) =  fsalt(i,j,iblk)   ! salt to ocean
          dataPtr_vice    (i1,j1,iblk) =   vice(i,j,iblk)   ! sea ice volume
          dataPtr_vsno    (i1,j1,iblk) =   vsno(i,j,iblk)   ! snow volume
          ! --- rotate these vectors from i/j to east/north ---
          ui = strairxT(i,j,iblk)
          vj = strairyT(i,j,iblk)
          dataPtr_strairxT(i1,j1,iblk) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! air ice stress
          dataPtr_strairyT(i1,j1,iblk) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! air ice stress
          ui = -strocnxT(i,j,iblk)
          vj = -strocnyT(i,j,iblk)
          dataPtr_strocnxT(i1,j1,iblk) = ui*cos(ANGLET(i,j,iblk)) - vj*sin(ANGLET(i,j,iblk))  ! ice ocean stress
          dataPtr_strocnyT(i1,j1,iblk) = ui*sin(ANGLET(i,j,iblk)) + vj*cos(ANGLET(i,j,iblk))  ! ice ocean stress
#endif
       enddo
       enddo
    enddo

    !-------------------------------------------------
  end subroutine ice_export
#endif



  subroutine cice_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)                        :: currTime
    character(len=*),parameter  :: subname='(cice_cap:cice_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) trim(subname),' --- finalize called --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Finalize

    write(info,*) trim(subname),' --- finalize completed --- '
    call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine cice_model_finalize

  subroutine CICE_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(cice_cap:CICE_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine CICE_AdvertiseFields

  subroutine CICE_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(cice_cap:CICE_RealizeFields)'
 
    rc = ESMF_SUCCESS

    do i = 1, nfields

      if (field_defs(i)%assoc) then
#ifdef CMEPS
        write(info, *) trim(subname), tag, ' Field ', trim(field_defs(i)%shortname), ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2)
        call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)

        if (trim(field_defs(i)%shortname) == trim(flds_scalar_name)) then
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(field_defs(i)%shortname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
          ! Create the scalar field
          call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          field = ESMF_FieldCreate(grid=grid, &
            farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
            name=field_defs(i)%shortname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if
#else
        write(info, *) trim(subname), tag, ' Field ', trim(field_defs(i)%shortname), ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO, rc=dbrc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
!          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
!          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#endif
      else
#ifdef CMEPS

        if (trim(field_defs(i)%shortname) == trim(flds_scalar_name)) then
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(field_defs(i)%shortname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
          ! Create the scalar field
          call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
            name=field_defs(i)%shortname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if
#else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
!          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#endif
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(trim(subname) // tag // " Field "// trim(field_defs(i)%stdname) // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(trim(subname) // tag // " Field "// trim(field_defs(i)%stdname) // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

  end subroutine CICE_RealizeFields

  !-----------------------------------------------------------------------------

  subroutine state_diagnose(State, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of state
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    character(len=*), intent(in), optional :: string
    integer, intent(out), optional  :: rc

    ! local variables
    integer                     :: i,j,n
    integer                     :: fieldCount
    character(len=64) ,pointer  :: fieldNameList(:)
    character(len=64)           :: lstring
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:)
#else
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:,:)
#endif
    integer                     :: lrc
    character(len=*),parameter  :: subname='(cice_cap:state_diagnose)'

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldNameList, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    do n = 1, fieldCount
      call State_GetFldPtr(State, fieldNameList(n), dataPtr, rc=lrc)
      if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(fieldNameList(n))//'  ', &
        minval(dataPtr),maxval(dataPtr),sum(dataPtr)
!      write(tmpstr,'(A)') trim(subname)//' '//trim(lstring)//':'//trim(fieldNameList(n))
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    deallocate(fieldNameList)

    if (present(rc)) rc = lrc

  end subroutine state_diagnose

  !-----------------------------------------------------------------------------

  subroutine state_reset(State, value, rc)
    ! ----------------------------------------------
    ! Set all fields to value in State
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    real(ESMF_KIND_R8), intent(in), optional :: value
    integer, intent(out), optional  :: rc

    ! local variables
    integer                     :: i,j,k,n
    integer                     :: fieldCount
    character(len=64) ,pointer  :: fieldNameList(:)
    real(ESMF_KIND_R8)          :: lvalue
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:)
#else
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:,:)
#endif
    character(len=*),parameter :: subname='(cice_cap:state_reset)'

    if (present(rc)) rc = ESMF_SUCCESS

    lvalue = 0._ESMF_KIND_R8
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    do n = 1, fieldCount
      call State_GetFldPtr(State, fieldNameList(n), dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

#ifdef CMEPS
      do j=lbound(dataPtr,2),ubound(dataPtr,2)
      do i=lbound(dataPtr,1),ubound(dataPtr,1)
         dataPtr(i,j) = lvalue
      enddo
      enddo
#else
      do k=lbound(dataPtr,3),ubound(dataPtr,3)
      do j=lbound(dataPtr,2),ubound(dataPtr,2)
      do i=lbound(dataPtr,1),ubound(dataPtr,1)
         dataPtr(i,j,k) = lvalue
      enddo
      enddo
      enddo
#endif

    enddo
    deallocate(fieldNameList)

  end subroutine state_reset

  !-----------------------------------------------------------------------------

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
#else
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:,:)
#endif
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(cice_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr

  !-----------------------------------------------------------------------------
  logical function FieldBundle_FldChk(FB, fldname, rc)
    type(ESMF_FieldBundle), intent(in) :: FB
    character(len=*)      ,intent(in) :: fldname
    integer, intent(out), optional :: rc

    ! local variables
    integer :: lrc
    character(len=*),parameter :: subname='(module_MEDIATOR:FieldBundle_FldChk)'

    if (present(rc)) rc = ESMF_SUCCESS

    FieldBundle_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=lrc)
    if (present(rc)) rc = lrc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
       FieldBundle_FldChk = .true.
    endif

  end function FieldBundle_FldChk

  !-----------------------------------------------------------------------------

  subroutine FieldBundle_GetFldPtr(FB, fldname, fldptr, rc)
    type(ESMF_FieldBundle), intent(in) :: FB
    character(len=*)      , intent(in) :: fldname
#ifdef CMEPS
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
#else
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:,:)
#endif
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(module_MEDIATOR:FieldBundle_GetFldPtr)'

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine FieldBundle_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine CICE_FieldsSetup
    character(len=*),parameter  :: subname='(cice_cap:CICE_FieldsSetup)'

!--------- import fields to Sea Ice -------------

! tcraig, don't point directly into cice data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_height_lowest"            , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height_lowest"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_spec_humid_height_lowest" , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_zonal_wind_height_lowest" , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_merid_wind_height_lowest" , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_pres_height_lowest"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_lw_flx"              , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dir_flx"      , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dif_flx"      , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dir_flx"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dif_flx"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_prec_rate"                , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_fprec_rate"               , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_temperature"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "s_surf"                        , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_zonal"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_merid"       , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_zonal"             , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_merid"             , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "freezing_melting_potential"    , "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "air_density_height_lowest"     , "will provide")
! this field is not used; however something about it being the 2nd field listed as 'toice' in 
! nems mediator requires it to be present
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height2m"            , "will provide")

!--------- export fields from Sea Ice -------------

    call fld_list_add(fldsFrIce_num, fldsFrIce, "sea_ice_surface_temperature"     , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dir_albedo"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dir_albedo"          , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dif_albedo"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dif_albedo"          , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "ice_mask"                        , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "ice_fraction"                    , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_zonal"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_merid"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_zonal"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_merid"         , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn"              , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn_vis_dir_flx"  , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn_vis_dif_flx"  , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn_ir_dir_flx"   , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn_ir_dif_flx"   , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_up_lw_flx_ice"              , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sensi_heat_flx_atm_into_ice", "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_laten_heat_flx_atm_into_ice", "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_evap_rate_atm_into_ice"     , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_fresh_water_to_ocean_rate"  , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_salt_rate"                  , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "net_heat_flx_to_ocn"             , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_ice_volume"                 , "will provide")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_snow_volume"                , "will provide")

  end subroutine CICE_FieldsSetup

  !-----------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
#ifdef CMEPS
    real(ESMF_KIND_R8), dimension(:,:), optional, target :: data
#else
    real(ESMF_KIND_R8), dimension(:,:,:), optional, target :: data
#endif
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(cice_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  subroutine dumpCICEInternal(grid, slice, stdname, nop, farray)

    type(ESMF_Grid)          :: grid
    integer, intent(in)      :: slice
    character(len=*)         :: stdname
    character(len=*)         :: nop
    real(ESMF_KIND_R8), dimension(:,:,:), target :: farray

    type(ESMF_Field)         :: field
    real(ESMF_KIND_R8), dimension(:,:), pointer  :: f2d
    integer                  :: i,j,rc

    if(.not. write_diagnostics) return ! remove this line to debug field connection

    field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
      indexflag=ESMF_INDEX_DELOCAL, &
      name=stdname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldGet(field, farrayPtr=f2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do j = lbound(f2d,2),ubound(f2d,2)
     do i = lbound(f2d,1),ubound(f2d,1)
      f2d(i,j) = farray(i+1,j+1,1)
     enddo
    enddo

#ifdef CMEPS
    call ESMF_FieldWrite(field, fileName='field_ice_internal_'//trim(stdname)//'.nc', &
      timeslice=slice, overwrite=overwrite_timeslice, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#else
    call ESMF_FieldWrite(field, fileName='field_ice_internal_'//trim(stdname)//'.nc', &
      timeslice=slice, overwrite=overwrite_timeslice, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    call ESMF_FieldDestroy(field, noGarbage=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

#ifdef CMEPS
  subroutine State_SetScalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output arguments
    real(ESMF_KIND_R8), intent(in)     :: scalar_value
    integer,            intent(in)     :: scalar_id
    type(ESMF_State),   intent(inout)  :: State
    character(len=*),   intent(in)     :: flds_scalar_name
    integer,            intent(in)     :: flds_scalar_num
    integer,            intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: lfield
    type(ESMF_VM)     :: vm
    real(ESMF_KIND_R8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(state_setscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (mytask == 0) then
       call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       endif
       farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine State_SetScalar

  subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
    ! ----------------------------------------------
    ! create a field with scalar data on the root pe
    ! ----------------------------------------------
    use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
    use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

    type(ESMF_Field) , intent(inout) :: field
    character(len=*) , intent(in)    :: flds_scalar_name
    integer          , intent(in)    :: flds_scalar_num
    integer          , intent(inout) :: rc

    ! local variables
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid
    character(len=*), parameter :: subname='(ice_import_export:SetScalarField)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
         ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    end subroutine SetScalarField
#endif
  !-----------------------------------------------------------------------------
end module cice_cap_mod
