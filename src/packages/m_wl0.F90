module m_wl0

  ! ModAEM 1.8
  ! Copyright(c) 1995-2008 WHPA Inc. and Vic Kelson
  !
  ! This program is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License
  ! as published by the Free Software Foundation; either version 2
  ! of the License, or(at your option) any later version.
  !
  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
  ! You should have received a copy of the GNU General Public License
  ! along with this program; if not, write to the Free Software
  ! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
  !
  ! Contact the author by e-mail at: vic@wittmanhydro.com
  ! Or by regular mail at:
  ! WHPA, Inc
  ! 320 W 8th St
  ! Bloomington, IN 47401

  !! module m_wl0
  !!
  !! Element module for 2-D discharge specified wells
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! well elements. Each element requires a location, radius and discharge.

  use u_constants
  use u_io
  use f_well
  use f_bwl
  use m_aqu
  use u_math

  implicit none

  public

  type :: WL0_WELL
    !! type WL0_WELL
    !!
    !! Type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: rDischarge
    !!     Discharge of the well; positive indicates extraction
    !!     and negative indicates injection.
    !!   integer :: iID
    !!     Identification label for the well(for interaction with e.g. GUIs)
    !!   integer :: iFWLIndex
    !!     Index for the well entry in the FWL module
    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rDischarge
    real(kind=AE_REAL) :: rRadius
    real(kind=AE_REAL) :: rCheckHead
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iFWLIndex
    real(kind=AE_REAL) :: rUnstressedHead
    real(kind=AE_REAL) :: rEfficiency
    logical :: lAdjustDischarge
    real(kind=AE_REAL) :: rDesiredHead
    real(kind=AE_REAL) :: rDesHeadRelaxation
    real(kind=AE_REAL) :: rLastDischarge
    ! Bessel well components for partially-penetrating wells
    logical :: lPpWell
    real(kind=AE_REAL) :: rScrBot
    real(kind=AE_REAL) :: rScrTop
    real(kind=AE_REAL) :: rKhKv
    real(kind=AE_REAL) :: rWtblHead
    real(kind=AE_REAL) :: rScrHead
    type(BWL_WELL), pointer :: bwl
    ! Options for the drawdown option
    logical :: lDdn
    real(kind=AE_REAL) :: rDdnDischarge
    real(kind=AE_REAL) :: rDdnDhd
  end type WL0_WELL

  type :: WL0_COLLECTION
    !! type WL0_COLLECTION
    !!
    !! Type that holds all wells in a layer
    !!
    !! Members:
    !!   type(WL0_WELL), dimension(:), pointer :: Wells
    !!     Array of WL0_WELL objects for the layer; dimensioned for the maximum
    !!     number of wells according to the input file(see WL0_Read)
    !!   integer :: iCount
    !!     The actual number of wells in use in the layer
    !!
    type(WL0_WELL), dimension(:), pointer :: Wells
    integer(kind=AE_INT) :: iCount
    ! Iterator information
    integer(kind=AE_INT) :: iIterWell
    ! True if the drawdown simulation is active
    logical :: lDdnEnabled
  end type WL0_COLLECTION


contains


  function WL0_Create(io) result(wl0)
    !! function WL0_Create
    !!
    !! Creates a new WL0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    wl0 => WL0_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(WL0_COLLECTION), pointer :: wl0
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(wl0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL0_Create: allocation failed")
    nullify(wl0%Wells)
    wl0%iCount = 0
    wl0%lDdnEnabled = .false.

    return
  end function WL0_Create


  subroutine WL0_Alloc(io, wl0)
    !! Subroutine WL0_Alloc
    !!
    !! Allocates wells for the WL0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call WL0_Alloc(lw0, iNWL)
    !!
    !! Arguments:
    !!    (in)    type(WL0_COLLECTION), pointer :: wl0
    !!              The WL0_COLLECTION object to be used
    !!    (in)    type(IO_STATUS), pointer :: io
    !!              Dimension information will be in the input buffer
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNWL
    integer(kind=AE_INT) :: iStat

    iNWL = iIO_GetInteger(io, 'iNWL', minimum = 0)
    allocate(wl0%Wells(iNWL), stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL0_Alloc: allocation failed")

    return
  end subroutine WL0_Alloc


  subroutine WL0_Destroy(io, wl0)
    !! subroutine WL0_Destroy
    !!
    !! Frees memory allocated for an WL0 Wells and WL0 Collection object
    !!
    !! Calling Sequence:
    !!     call WL0_Destroy(wl0)
    !!
    !! Arguments:
    !!  type(WL0_COLLECTION), pointer :: wl0
    !!              Pointer to the WL0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_Destroy: WL0_Create has not been called")
    end if

    if (associated(wl0%Wells)) then
      deallocate(wl0%Wells, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "WL0_Destroy: deallocation of Wells failed")
    end if
    deallocate(wl0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL0_Destroy: deallocation failed")

    return
  end subroutine WL0_Destroy


  subroutine WL0_New(io, wl0, Well)
    !! function WL0_New
    !!
    !! Adds a new WL0_WELL object to the WL0_COLLECTION 'wl0'
    !!
    !! Calling Sequence:
    !!    call WL0_New(io, wl0, Well)
    !!
    !! Arguments:
    !!    (in)    type(WL0_COLLECTION), pointer :: wl0
    !!              The WL0_COLLECTION object to be used
    !!    (in)    type(WL0_WELL), pointer :: Well
    !!              Vector that defines the points along the barrier
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(WL0_WELL) :: Well
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_New: WL0_Create has not been called")
    end if

    call IO_Assert(io, (wl0%iCount < size(wl0%Wells)), &
         "WL0_New: Space exhausted")

    wl0%iCount = wl0%iCount + 1
    wl0%Wells(wl0%iCount) = Well

    return
  end subroutine WL0_New


  function iWL0_GetID(io, wl0, iIndex) result(iID)
    !! Returns the ID number for the well at index 'iIndex'
    type(WL0_COLLECTION), pointer :: wl0
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= wl0%iCount), "Internal error -- no such index")
    iID = wl0%Wells(iIndex)%iID

    return
  end function iWL0_GetID


  subroutine WL0_PreSolve(io, wl0)
    !! subroutine WL0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call WL0_PreSolve(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io

    return
  end subroutine WL0_PreSolve


  function iWL0_GetInfo(io, wl0, iOption, iIteration) result(iValue)
    !! function WL0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iWL0_GetInfo(io, wl0, iOption)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_GetInfo: WL0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = wl0%iCount
      case default
        iValue = 0
    end select

    return
  end function iWL0_GetInfo


  subroutine WL0_SetupFunctions(io, wl0, fwl, aqu)
    !! subroutine WL0_Setup
    !!
    !! This routine sets up the functions in f_well for the well elements
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call WL0_Setup(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             FWL_COLLECTION to be used
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(FWL_COLLECTION), pointer :: fwl
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWel
    type(WL0_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_Setup: WL0_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "WL0_Setup: Illegal FWL_COLLECTION object")
    end if

    do iWel = 1, wl0%iCount
      ! Create a well function in FWL for each well
      wel => wl0%Wells(iWel)
      call FWL_New(io, fwl, wel%cZ, wel%rDischarge, wel%rRadius, ELEM_WL0, iWel, -1, -1, wel%iFWLIndex)
      if (wel%lPpWell) then
        wel%bwl => BWL_New(io, wel%cZ, wel%rRadius, wel%rScrBot, wel%rScrTop, &
                           rAQU_Base(io, aqu, wel%cZ), &
                           rAQU_Base(io, aqu, wel%cZ) + &
                             rAQU_SatdThickness(io, aqu, wel%cZ, real(cAQU_Potential(io, aqu, wel%cZ))), &
                           rAQU_HydCond(io, aqu, wel%cZ), &
                           wel%rKhKv)
      end if
    end do

    return
  end subroutine WL0_SetupFunctions


  subroutine WL0_Update(io, wl0, fwl)
    !! subroutine WL0_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    WL0_Update(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(FWL_COLLECTION), pointer :: fwl
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWell
    type(WL0_WELL), pointer :: wel


    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_Update: WL0_Create has not been called")
    end if

    do iWell = 1, wl0%iCount
      wel => wl0%Wells(iWell)
      call FWL_Update(io, fwl, wel%iFWLIndex, wel%rDischarge)
    end do

    return
  end subroutine WL0_Update


  function WL0_EnableDrawdown(io, wl0) result(iChanges)
    !! Enables the elements that are subject to the "drawdown" flag
    !! This allows the modeler to compute "drawdown" simulations based on the difference between
    !! an "unstressed" condition and a "stressed" condition. Each element module that implements
    !! the DDN operator provides an XXX_EnableDrawdown() subroutine that makes the necessary
    !! adjustments.
    !!
    !! Returns the number of changes made when enabling the drawdown simulation
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL0_COLLECTION), pointer :: wl0
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iwel
    type(WL0_WELL), pointer :: wel

    iChanges = 0
    do iwel = 1, wl0%iCount
      wel => wl0%Wells(iwel)
      wel%rUnstressedHead = wel%rCheckHead
      if (wel%lDdn) then
        iChanges = iChanges + 1
        wel%rDischarge = wel%rDdnDischarge
        wl0%lDdnEnabled = .true.
        if (wel%rDdnDhd /= rHUGE) then
          wel%rDesiredHead = wel%rDdnDhd
          wel%lAdjustDischarge = .true.
        end if
      end if
    end do

    return
  end function WL0_EnableDrawdown


  subroutine WL0_SolvePartialPenetration(io, wl0, aqu)
    !! Solves all the partially-penetrating wells
    !! Traverses all the wells in wl0 and uses WL0_SolvePPWell to determine the additional
    !! drawdown due to partial penetration based on the work of Bakker (2001)
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL0_COLLECTION), pointer :: wl0
    type(AQU_COLLECTION), pointer :: aqu
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iwel
    type(WL0_WELL), pointer :: wel
    real(kind=AE_REAL) :: rAqTop
    logical :: lChange

    call IO_MessageText(io, "Solving partially-penetrating wells")
    do iwel = 1, wl0%iCount
      wel => wl0%Wells(iwel)
      if (wel%lPpWell) then
        rAqTop = rAQU_Base(io, aqu, wel%cZ) + &
                 rAQU_SatdThickness(io, aqu, wel%cZ, real(cAQU_Potential(io, aqu, wel%cZ)))
        call BWL_Solve(io, wel%bwl, rAqTop, lChange)
      end if
    end do

    return
  end subroutine WL0_SolvePartialPenetration


  function WL0_AdjustDischarges(io, wl0, aqu) result(iChange)
    !! Uses the "desired head" setting to iteratively adjust the pumping rate of wells in the model
    !! If an adjustable well is also a partially-penetrating well, the partial-penetration correction
    !! is also to be used. In all cases, the well efficiency will be included in the calculation.
    !! Returns iChange = the number of wells that are changed.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL0_COLLECTION), pointer :: wl0
    type(AQU_COLLECTION), pointer :: aqu
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChange
    ! [ LOCALS ]
    real(kind=AE_REAL) :: wt_head, scr_head
    integer(kind=AE_INT) :: i
    real(kind=AE_REAL) :: rAqTop, spec_cap
    logical :: lChange
    type(WL0_WELL), pointer :: wel

    iChange = 0
    do i = 1, wl0%iCount
      wel => wl0%Wells(i)
      if (wel%lAdjustDischarge) then
        ! Is it a partially-penetrating well?
        if (wel%lPpWell) then
          rAqTop = rAQU_Base(io, aqu, wel%cZ) + &
                   rAQU_SatdThickness(io, aqu, wel%cZ, real(cAQU_Potential(io, aqu, wel%cZ)))
          call BWL_Solve(io, wel%bwl, rAqTop, lChange)
        end if
        ! Get the head at the well screen
        call WL0_GetHeadAtWell(io, wel, wt_head, scr_head)
        spec_cap = wel%rDischarge / (wel%rUnstressedHead - scr_head)
        wel%rDischarge = wel%rDischarge + wel%rDesHeadRelaxation * spec_cap * (scr_head - wel%rDesiredHead)
        write (unit=IO_MessageBuffer, &
              fmt='("  Adjusted well ", i6, " Q:", g12.5, " DH:", g12.5,' // &
                  '" MH:", g12.5, "E:", g12.5)') &
              wel%iID, wel%rDischarge, wel%rDesiredHead, scr_head, wel%rDesiredHead - scr_head
        call IO_MessageText(io)
        iChange = iChange+1
      end if
    end do

    return
  end function WL0_AdjustDischarges

  subroutine WL0_FindWell(io, wl0, iWellID, cZWell, rDischarge, rRadius, iFWLIndex, lFound)
    !! subroutine WL0_FindWell
    !!
    !! Finds the well specified by the Well ID and returns its parameters
    !!
    !! Calling Sequence:
    !!    lFound = lWL0_FindWell(io, wl0, iWellID, cZWell, rDischarge, rRadius, iFWLIndex)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    integer :: iWellID
    !!             The well ID number
    !!   (out)   complex :: cZWell
    !!             Location of the well
    !!   (out)   real :: rDischarge
    !!             The discharge of the well
    !!   (out)   real :: rRadius
    !!             The radius of the well
    !!   (out)   integer :: iFWLIndex
    !!             Index of the well in f_well
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    integer(kind=AE_INT), intent(in) :: iWellID
    complex(kind=AE_REAL), intent(out) :: cZWell
    real(kind=AE_REAL), intent(out) :: rRadius, rDischarge
    integer(kind=AE_INT), intent(out) :: iFWLIndex
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(WL0_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_FindWell: WL0_Create has not been called")
    end if

    lFound = .false.
    do i = 1, wl0%iCount
      wel => wl0%Wells(i)
      if (wel%iID == iWellID) then
        lFound = .true.
        cZWell = wel%cZ
        rDischarge = wel%rDischarge
        rRadius = wel%rRadius
        iFWLIndex = wel%iFWLIndex
      end if
    end do

    return
  end subroutine WL0_FindWell


  subroutine WL0_FindWellPointer(io, wl0, iWellID, well, lFound)
    !! subroutine WL0_FindWell
    !!
    !! Finds the well specified by the Well ID and returns a pointer to it
    !!
    !! Calling Sequence:
    !!    lFound = lWL0_FindWellPointer(io, wl0, iWellID, well, lfound)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    integer :: iWellID
    !!             The well ID number
    !!   (out)   type(WL0_WELL :: well
    !!             Pointer to the well
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!   (in)    type(IO_STATUS), pointer :: io
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    integer(kind=AE_INT), intent(in) :: iWellID
    type(WL0_WELL), pointer :: well
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_FindWell: WL0_Create has not been called")
    end if

    lFound = .false.
    do i = 1, wl0%iCount
      well => wl0%Wells(i)
      if (well%iID == iWellID) then
        lFound = .true.
        return
      end if
    end do

    return
  end subroutine WL0_FindWellPointer


  subroutine WL0_ResetIterator(io, wl0)
    !! subroutine WL0_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call WL0_ResetIterator(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: wl0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_ResetIterator: WL0_Create has not been called")
    end if

    wl0%iIterWell = 0

    return
  end subroutine WL0_ResetIterator


  function WL0_NextIterator(io, wl0) result(itr)
    !! function WL0_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call WL0_NextIterator(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: wl0
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_NextIterator: WL0_Create has not been called")
    end if

    wl0%iIterWell = wl0%iIterWell + 1
    if (wl0%iIterWell > wl0%iCount) then
      nullify(itr)
    else
      allocate(itr)
      itr%iValueSelector = VALUE_POTENTIAL
      itr%iElementString = wl0%iIterWell
      allocate(itr%cZ(1))
      itr%cZ(1) = wl0%Wells(wl0%iIterWell)%cZ
    end if

    return
  end function WL0_NextIterator


  subroutine WL0_SetIterator(io, wl0, aqu, fwl, itr, cValue)
    !! function WL0_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call WL0_SetIterator(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: wl0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(AQU_COLLECTION), pointer :: aqu
    type(FWL_COLLECTION), pointer :: fwl
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(WL0_WELL), pointer :: wel
    complex(kind=AE_REAL) :: cPotWithoutWell

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_NextIterator: WL0_Create has not been called")
      call IO_Assert(io, (wl0%iIterWell <= wl0%iCount), &
           "WL0_SetIterator: Iterator out of range")
    end if
    wel => wl0%Wells(itr%iElementString)
    wel%rCheckHead = rAQU_PotentialToHead(io, aqu, real(cValue, AE_REAL), wel%cZ)
    ! If we're not currently in a drawdown simulation, the "desired head option" needs an estimate
    ! of the drawdown for estimation purposes
    if (.not. (wel%lDdn .and. wl0%lDdnEnabled)) then
      cPotWithoutWell = cValue - cFWL_Potential(io, fwl, wel%cZ+wel%rRadius, wel%iFWLIndex, 1)
      wel%rUnstressedHead = rAQU_PotentialToHead(io, aqu, real(cPotWithoutWell), wel%cZ)
    end if

    return
  end subroutine WL0_SetIterator


  subroutine WL0_GetHeadAtWell(io, wel, wt_head, scr_head)
    !! Computes the estimated water-table head and screen head for well 'i' in the WL0_COLLECTION.
    !! Includes the effects of partial penetration and well efficiency. Returns the estimated
    !! water-table elevation (more precisely, the head at the top of the saturated thickness) and
    !! the head at the well screen.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL0_WELL), pointer :: wel
    real(kind=AE_REAL), intent(out) :: wt_head
    real(kind=AE_REAL), intent(out) :: scr_head

    ! Look up PP well values if necessary
    if (wel%lPpWell) then
      wt_head = wel%rCheckHead - wel%rDischarge * wel%bwl%rInfl(1)
      scr_head = wel%rCheckHead - wel%rDischarge * wel%bwl%rInfl(wel%bwl%iLayer)
    else
      wt_head = wel%rCheckHead
      scr_head = wel%rCheckHead
    end if

    ! Adjust for efficiency
    scr_head = wel%rUnstressedHead - (wel%rUnstressedHead - wel%rCheckHead) / wel%rEfficiency

    return
  end subroutine WL0_GetHeadAtWell

  subroutine WL0_Read(io, wl0)
    !! subroutine WL0_Read
    !!
    !! Reads the wells for the specified WL0_COLLECTION from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call WL0_Read(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be populated
    !!
    !! The format of the WL0 section of the input file appears as follows:
    !! WL0
    !! DIM NWells
    !!     x y q r id
    !!     ... Up to NWells
    !!
    !! NOTE: It is assumed that the WL0 line was found by the caller
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(5), parameter :: dirDirectives = (/ dirEND, dirPPW, dirEFF, dirDHD, dirDDN /)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rDischarge, rRad, rScrBot, rScrTop, rKhKv, rDesiredHead, rEfficiency, &
                          rDesHeadRelaxation, rDdnDischarge
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxWel
    logical :: lFlag
    type(WL0_WELL), pointer :: wel

    call IO_MessageText(io, "  Reading WL0 module input")

    call IO_Assert(io, (associated(wl0)), "WL0_Read: WL0_Create has not been called")

    ! Process input
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "WL0_Read: I/O Error")
          exit
        case (kOpFileEOF)
          ! EOF is unexpected for all ModWL0 "ifXXXRead" routines.
          ! Report the condition, but proceed as if EOD was found.
          call IO_Assert(io, .false., "WL0_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(wl0%Wells)), "WL0_Read: No space allocated")
          call IO_Assert(io, (wl0%iCount < size(wl0%Wells)), "WL0_Read: Space exhausted")
          cZ = cIO_GetCoordinate(io, 'cZ', extents=.true.)
          rDischarge = rIO_GetReal(io, 'rDischarge')
          rRad = rIO_GetReal(io, 'rRad', minimum = rTINY)
          iID = iIO_GetInteger(io, 'iID')
          wl0%iCount = wl0%iCount+1
          wel => wl0%Wells(wl0%iCount)
          wel%cZ = cZ
          wel%rDischarge = rDischarge
          wel%rRadius = rRad
          wel%iID = iID
          wel%lPpWell = .false.
          wel%rScrBot = -rHUGE
          wel%rScrTop = rHUGE
          wel%rKhKv = rONE
          wel%rEfficiency = rONE
          wel%rDesiredHead = -rHUGE
          wel%lAdjustDischarge = .false.
          wel%lDdn = .false.
          ! No FWL is declared here; see WL0_Setup
          wel%iFWLIndex = -1
        case (kOpPPW)
          ! Partially-penetrating well info follows
          call IO_Assert(io, (wl0%iCount>0), "No current well to modify")
          wel => wl0%Wells(wl0%iCount)
          rScrBot = rIO_GetReal(io, 'rScrBot')
          rScrTop = rIO_GetReal(io, 'rScrTop', minimum=rScrBot+rTINY)
          rKhKv = rIO_GetReal(io, 'rKhKv', def=rONE)
          wel%lPpWell = .true.
          wel%rScrBot = rScrBot
          wel%rScrTop = rScrTop
          wel%rKhKv = rKhKv
        case (kOpEFF)
          ! Well efficiency  follows
          call IO_Assert(io, (wl0%iCount>0), "No current well to modify")
          wel => wl0%Wells(wl0%iCount)
          rEfficiency = rIO_GetReal(io, 'rEfficiency', def=rONE)
          wel%rEfficiency = rEfficiency
        case (kOpDHD)
          ! Desired head follows
          call IO_Assert(io, (wl0%iCount>0), "No current well to modify")
          wel => wl0%Wells(wl0%iCount)
          rDesiredHead = rIO_GetReal(io, 'rDesiredHead', def=rONE)
          rDesHeadRelaxation = rIO_GetReal(io, 'rDesHeadRelaxation', def=rONE_TENTH)
          wel%lAdjustDischarge = .true.
          wel%rDesiredHead = rDesiredHead
          wel%rDesHeadRelaxation = rDesHeadRelaxation
        case (kOpDDN)
          ! Pumping rate for drawdown simulation follows
          call IO_Assert(io, (wl0%iCount>0), "No current well to modify")
          wel => wl0%Wells(wl0%iCount)
          rDdnDischarge = rIO_GetReal(io, 'rDdnDischarge', def=rONE)
          rDesiredHead = rIO_GetReal(io, 'rDesiredHead', def=rHUGE)
          rDesHeadRelaxation = rIO_GetReal(io, 'rDesHeadRelaxation', def=rONE_TENTH)
          wel%lAdjustDischarge = .false.
          wel%rDdnDischarge = rDdnDischarge
          wel%rDdnDhd = rDesiredHead
          wel%rDesHeadRelaxation = rDesHeadRelaxation
          wel%lDdn = .true.
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
      end select
    end do

    call IO_MessageText(io, "  Leaving WL0 module")

    return
  end subroutine WL0_Read


  subroutine WL0_Inquiry(io, wl0, aqu, iLU)
    !! subroutine WL0_Inquiry
    !!
    !! Writes an inquiry report for all wells to iLU
    !!
    !! Calling Sequence:
    !!    call WL0_Inquiry(io, wl0, iLU)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(WL0_WELL), pointer :: wel
    real(kind=AE_REAL) :: wt_head, scr_head, rDdn

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_Inquiry: WL0_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#WL0, ID, X, Y, DISCHARGE, RADIUS, CHECK_HEAD, UN_HEAD, WT_HEAD, SCR_HEAD, DDN, DRY"")")
    do i = 1, wl0%iCount
      wel => wl0%Wells(i)
      call WL0_GetHeadAtWell(io, wel, wt_head, scr_head)
      if (wl0%lDdnEnabled) then
        rDdn = wel%rUnstressedHead - scr_head
      else
        rDdn = rHUGE
      end if
      write (unit=iLU, &
             fmt="(""WL0"", 1("", "", i9), 9("", "", e16.8), "", "", l1)" &
             ) wel%iID, &
             cIO_WorldCoords(io, wel%cZ), &
             wel%rDischarge, &
             wel%rRadius, &
             wel%rCheckHead, &
             wel%rUnstressedHead, &
             wt_head, &
             scr_head, &
             rDdn, &
             wel%rCheckHead <= rAQU_Base(io, aqu, wel%cZ)
    end do

    return
  end subroutine WL0_Inquiry


  subroutine WL0_Report(io, wl0)
    !! subroutine WL0_Report
    !!
    !! Writes a debugging report for all wells to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call WL0_Report(wl0)
    !!
    !! Arguments:
    !!   (in)    type(WL0_COLLECTION), pointer :: wl0
    !!             WL0_COLLECTION to be used
    !!
    ! [ ARGUMENTS ]
    type(WL0_COLLECTION), pointer :: wl0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    real(kind=AE_REAL) :: wt_head, scr_head, rDdn
    type(WL0_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl0)), &
           "WL0_Inquiry: WL0_Create has not been called")
    end if

    call HTML_Header('Module WL0', 1)
    call HTML_Header('Discharge-specified well information', 2)

    if (associated(wl0%Wells)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of wells', wl0%iCount)
      call HTML_AttrInteger('Number of FWL functions', iWL0_GetInfo(io, wl0, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iWL0_GetInfo(io, wl0, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iWL0_GetInfo(io, wl0, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iWL0_GetInfo(io, wl0, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iWL0_GetInfo(io, wl0, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      call HTML_Header('Wells', 4)

      call HTML_StartTable()
      call HTML_TableHeader((/'Well    ', 'ID      ', 'FWL #   ', 'X       ', 'Y       ', 'Q       ', 'R       ', &
                              'UN Head ', 'DF Head ', 'WT Head ', 'SCR Head' /))
      do i = 1, wl0%iCount
        wel => wl0%Wells(i)
        call WL0_GetHeadAtWell(io, wel, wt_head, scr_head)
        if (wl0%lDdnEnabled) then
          rDdn = wel%rUnstressedHead - scr_head
        else
          rDdn = rHUGE
        end if

        call HTML_StartRow()
        call HTML_ColumnInteger((/i, wel%iID, wel%iFWLIndex/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, wel%cZ)/))
        call HTML_ColumnReal((/ wel%rDischarge, wel%rRadius, wel%rCheckHead, wel%rUnstressedHead, &
                                wt_head, scr_head, rDdn /))
        call HTML_EndRow()
      end do
      call HTML_EndTable()

    else
      call HTML_Header('No wells defined', 3)
    end if

    return
  end subroutine WL0_Report

end module m_wl0
