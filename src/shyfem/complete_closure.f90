!=======================================================================
! complete_closure.f90
!
! Layer-dependent, bottom-mounted closure for SHYFEM.
!
! Conventions:
!   file flag = 0 : gate open
!   file flag = 1 : gate closed
!   geyer     = 0 : gate closed
!   geyer     = 1 : gate open
!
! The gate height is measured vertically upwards from the bed.  Layers
! above the crest remain open.  A layer intersected by the crest receives
! a proportional reduction.
!
! iff_ts_init is called with nintp=2 so that SHYFEM retains the two time
! records required by its time-series interpolation machinery. The returned
! value is subsequently rounded and clipped to 0/1 by coclose_read_flag.
!=======================================================================

module coclose

  implicit none

  type :: coclose_entry
    integer :: isc = 0
    integer, allocatable :: kboc(:)       !external/global closure nodes
    integer, allocatable :: kboc_local(:) !rank-local nodes; 0 if not present
    integer, allocatable :: ieboc(:)      !affected elements
    real :: height = 3.0                  !height above bed [m]
    real :: geyer = 1.0                   !0 closed, 1 open
    integer :: last_flag = 0              !file flag: 0 open, 1 closed
    character(len=80) :: cfile = ''
    integer :: idfile = 0
  end type coclose_entry

  integer, save :: icoclose = 0
  integer, save :: ncoclose = 0
  integer, save, private :: ndim_coclose = 0

  type(coclose_entry), allocatable, save :: coclose_entries(:)

  ! Layer-dependent multiplier used by hydrodynamic.f90.
  ! First index: vertical layer; second index: element.
  real, allocatable, save :: rclosurev(:,:)

contains

  subroutine coclose_alloc_entry

    type(coclose_entry), allocatable :: aux(:)
    integer :: olddim

    if (ndim_coclose == 0) then
      ndim_coclose = 4
      allocate(coclose_entries(ndim_coclose))
    else
      olddim = ndim_coclose
      ndim_coclose = 2*ndim_coclose
      allocate(aux(ndim_coclose))
      aux(1:olddim) = coclose_entries(1:olddim)
      call move_alloc(aux,coclose_entries)
    end if

  end subroutine coclose_alloc_entry

!-----------------------------------------------------------------------

  subroutine coclose_new_entry(id)

    integer, intent(out) :: id

    ncoclose = ncoclose + 1
    if (ncoclose > ndim_coclose) call coclose_alloc_entry

    id = ncoclose
    coclose_entries(id)%isc = 0
    coclose_entries(id)%height = 3.0
    coclose_entries(id)%geyer = 1.0
    coclose_entries(id)%last_flag = 0
    coclose_entries(id)%cfile = ''
    coclose_entries(id)%idfile = 0

  end subroutine coclose_new_entry

!-----------------------------------------------------------------------

  subroutine coclose_read_flag(id,dtime,flag)

    integer, intent(in) :: id
    double precision, intent(in) :: dtime
    integer, intent(out) :: flag

    real :: file_value

    call iff_ts_intp1(coclose_entries(id)%idfile,dtime,file_value)

    flag = nint(file_value)
    flag = max(0,min(1,flag))

  end subroutine coclose_read_flag

!-----------------------------------------------------------------------

  subroutine coclose_set_vertical_mask(id)

    use basin, only : nel
    use levels, only : ilhv, jlhv
    use mod_layer_thickness, only : hdeov
    use shympi

    integer, intent(in) :: id

    integer :: i,ie,l,ilevel,jlevel
    real :: gate_height,height_left,thickness
    real :: blocked_fraction,layer_factor,geyer

    gate_height = max(0.0,coclose_entries(id)%height)
    geyer = coclose_entries(id)%geyer
    geyer = max(0.0,min(1.0,geyer))

    do i=1,size(coclose_entries(id)%ieboc)
      ie = coclose_entries(id)%ieboc(i)
      if (ie <= 0 .or. ie > nel) cycle

      ilevel = ilhv(ie)                  !deepest active layer
      jlevel = jlhv(ie)                  !uppermost active layer
      height_left = gate_height

      ! SHYFEM layers are traversed here from the bed upwards.
      do l=ilevel,jlevel,-1
        thickness = hdeov(l,ie)
        if (thickness <= 0.0) cycle

        blocked_fraction = min(1.0,height_left/thickness)
        blocked_fraction = max(0.0,blocked_fraction)

        ! Fully above crest: blocked_fraction=0 -> layer_factor=1.
        ! Fully below crest: blocked_fraction=1 -> layer_factor=geyer.
        layer_factor = 1.0 - blocked_fraction*(1.0-geyer)

        ! min() allows two independently defined closures to overlap safely.
        rclosurev(l,ie) = min(rclosurev(l,ie),layer_factor)

        height_left = max(0.0,height_left-thickness)
      end do
    end do

  end subroutine coclose_set_vertical_mask

end module coclose

!=======================================================================
! Read one $coclo section. The short label is intentional: legacy SHYFEM
! section/parameter names are commonly limited to six characters.
!
! Proposed input:
!   $coclo1 Gnocca
!       kboc  = 1137 1236 857
!       height = 3.0
!       cfile = 'gnocca_closure.dat'
!   $end
!=======================================================================

subroutine rdcoclos(isc)

  use coclose
  use arrays

  implicit none

  integer, intent(in) :: isc

  character(len=6) :: name
  character(len=80) :: text
  double precision :: dvalue
  integer :: id,iweich,ikboc
  integer :: nrdnxt

  call coclose_new_entry(id)
  coclose_entries(id)%isc = isc

  ikboc = 0
  call init_array(coclose_entries(id)%kboc)

  do
    iweich = nrdnxt(name,dvalue,text)
    call to_lower(name)

    if (iweich == 0) exit
    if (iweich < 0 .or. iweich > 3) goto 98

    if (iweich == 2 .and. name /= 'kboc') goto 93
    if (iweich == 3 .and. name /= 'cfile') goto 93

    if (name == 'kboc') then
      call append_to_array(coclose_entries(id)%kboc,ikboc,nint(dvalue))
    else if (name == 'height') then
      coclose_entries(id)%height = real(dvalue)
    else if (name == 'cfile') then
      coclose_entries(id)%cfile = text
    else
      goto 96
    end if
  end do

  call trim_array(coclose_entries(id)%kboc,ikboc)
  return

93 continue
  write(6,*) 'Variable not allowed in coclose context: ',name
  stop 'error stop: rdcoclos'
96 continue
  write(6,*) 'Unrecognised coclose variable: ',name
  stop 'error stop: rdcoclos'
98 continue
  write(6,*) 'iweich = ',iweich
  write(6,*) 'Read error in complete closure section: ',isc
  stop 'error stop: rdcoclos'

end subroutine rdcoclos

!=======================================================================
! Validate input and create the rank-local node mapping.
!
! kboc is deliberately preserved in external/global numbering. In an MPI
! run, ipint() returns 0 on ranks that do not contain a given global node;
! this is normal and must not be treated as a global mesh error. The matching
! local index is stored in kboc_local, or 0 when absent on the current rank.
! Call this after all STR sections have been read and before coclose_init.
!=======================================================================

subroutine ckcoclos

  use coclose
  use shympi

  implicit none

  integer :: id,i,knode,kint,nlocal
  integer :: ipint
  logical :: bstop

  bstop = .false.

  do id=1,ncoclose
    if (.not. allocated(coclose_entries(id)%kboc)) then
      write(6,*) 'coclose section has no kboc: ',id
      bstop = .true.
      cycle
    end if

    if (size(coclose_entries(id)%kboc) <= 0) then
      write(6,*) 'coclose section has empty kboc: ',id
      bstop = .true.
    end if

    if (coclose_entries(id)%height <= 0.0) then
      write(6,*) 'coclose height must be positive: ',id
      bstop = .true.
    end if

    if (len_trim(coclose_entries(id)%cfile) == 0) then
      write(6,*) 'coclose section has no cfile: ',id
      bstop = .true.
    end if

    if (allocated(coclose_entries(id)%kboc_local)) then
      deallocate(coclose_entries(id)%kboc_local)
    end if
    allocate(coclose_entries(id)%kboc_local( &
             size(coclose_entries(id)%kboc)))
    coclose_entries(id)%kboc_local = 0
    nlocal = 0

    do i=1,size(coclose_entries(id)%kboc)
      knode = coclose_entries(id)%kboc(i)
      if (knode <= 0) then
        write(6,*) 'invalid coclose node number: ',knode
        bstop = .true.
      else
        kint = ipint(knode)
        if (kint > 0) then
          coclose_entries(id)%kboc_local(i) = kint
          nlocal = nlocal + 1
        end if
      end if
    end do

    write(6,*) 'coclose local node mapping: section=',id, &
               ' local=',nlocal, &
               ' global=',size(coclose_entries(id)%kboc)
  end do

  if (bstop) stop 'error stop: ckcoclos'

end subroutine ckcoclos

!=======================================================================
! Initialise affected elements, the time-series files and rclosurev.
!=======================================================================

subroutine coclose_init

  use coclose
  use basin
  use levels, only : nlvdi
  use arrays
  use shympi

  implicit none

  logical, save :: binit = .false.
  logical :: belem
  integer :: id,i,ie,k,nvert,idfile
  integer :: nieboc
  integer :: kvert(10)
  integer, allocatable :: index(:)
  double precision :: dtime
  real :: getpar

  if (binit) return
  binit = .true.

  ! Six-character name, compatible with the legacy SHYFEM parameter reader.
  icoclose = nint(getpar('icoclo'))
  if (icoclose <= 0) return

  if (ncoclose <= 0) then
    write(6,*) 'icoclo enabled but no $coclo sections were read'
    stop 'error stop: coclose_init'
  end if

  allocate(rclosurev(nlvdi,nel))
  rclosurev = 1.0
  allocate(index(nkn))

  call get_act_dtime(dtime)

  do id=1,ncoclose
    index = 0
    if (.not. allocated(coclose_entries(id)%kboc_local)) then
      write(6,*) 'coclose local nodes not initialized: ',id
      stop 'error stop: coclose_init'
    end if

    do i=1,size(coclose_entries(id)%kboc_local)
      k = coclose_entries(id)%kboc_local(i)
      if (k > 0 .and. k <= nkn) index(k) = 1
    end do

    call init_array(coclose_entries(id)%ieboc)
    nieboc = 0

    do ie=1,nel
      call nindex(ie,nvert,kvert)
      if (nvert > size(kvert)) then
        stop 'error stop: coclose_init nvert'
      end if

      belem = .false.
      do i=1,nvert
        if (index(kvert(i)) > 0) belem = .true.
      end do

      if (belem) then
        call append_to_array_i(coclose_entries(id)%ieboc,nieboc,ie)
      end if
    end do

    call trim_array(coclose_entries(id)%ieboc,nieboc)

    ! nintp=2: two-record interpolation window; nvar=1: one gate flag.
    ! coclose_read_flag rounds and clips the interpolated value to 0/1.
    call iff_ts_init(dtime,trim(coclose_entries(id)%cfile),2,1,idfile)
    coclose_entries(id)%idfile = idfile
    coclose_entries(id)%last_flag = 0
    coclose_entries(id)%geyer = 1.0

    write(6,*) 'complete closure initialized: ',id, &
               ' global nodes=',size(coclose_entries(id)%kboc), &
               ' local nodes=',count(coclose_entries(id)%kboc_local > 0), &
               ' elements=',size(coclose_entries(id)%ieboc), &
               ' height=',coclose_entries(id)%height
  end do

  deallocate(index)

end subroutine coclose_init

!=======================================================================
! Update all closures. Call once before the element momentum calculations
! in every time step.
!=======================================================================

subroutine coclose_handle(ic)

  use coclose
  use shympi

  implicit none

  integer, intent(out) :: ic
  integer :: id,flag
  double precision :: dtime

  ic = 0
  if (icoclose <= 0) return
  if (.not. allocated(rclosurev)) return

  call get_act_dtime(dtime)

  ! Rebuild the complete vertical mask because layer thickness can change
  ! with water level in z-level/z-star configurations.
  rclosurev = 1.0

  do id=1,ncoclose
    call coclose_read_flag(id,dtime,flag)

    ! File flag: 0=open, 1=closed.  Geyer uses the opposite convention.
    coclose_entries(id)%geyer = 1.0-real(flag)

    if (flag /= coclose_entries(id)%last_flag) then
      ic = 1
      write(6,*) 'complete closure changed: id=',id, &
                 ' flag=',flag,' geyer=',coclose_entries(id)%geyer
      coclose_entries(id)%last_flag = flag
    end if

    call coclose_set_vertical_mask(id)
  end do

end subroutine coclose_handle

!=======================================================================

subroutine prcoclos

  use coclose
  use shympi

  implicit none

  integer :: id

  if (ncoclose <= 0) return

  write(6,*)
  write(6,*) '====== complete closure sections ======'
  write(6,*) 'number of sections: ',ncoclose

  do id=1,ncoclose
    write(6,*) 'section: ',coclose_entries(id)%isc
    write(6,*) '  nodes:  ',size(coclose_entries(id)%kboc)
    write(6,*) '  height: ',coclose_entries(id)%height
    write(6,*) '  cfile:  ',trim(coclose_entries(id)%cfile)
  end do

end subroutine prcoclos

subroutine tscoclos
  implicit none
  call prcoclos
end subroutine tscoclos
