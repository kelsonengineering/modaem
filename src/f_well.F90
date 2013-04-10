module f_well

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


  !! module f_well(fwl)
  !!
  !! Written by Victor A.Kelson
  !!
  !! Revision History:
  !!   1.0.0   22 March 1999
  !!           First "source code release" version, adapted from
  !!           previous testing version.
  !!
  !! Module of data structures and functions for 2-D, steady-state
  !! wells.
  !!
  !! This module encapsulates containers of well functions
  !! in layers.  The various data structures are organized for
  !! efficient computation on vector or parallel machinery.
  !!
  !! Limitations:
  !!   At present, this module has no hooks for series expansion
  !!   block implementation.  No known computational bugs exist.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   io         --  Universal ModAEM I/O functions and constants
  !!   i_well     --  Influence function module for wells
  !!
  use u_constants
  use u_io
  use i_well

  implicit none

  public


  type, public :: FWL_WELL
    !! type FWL_WELL
    !!
    !! PUBLIC type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: rDischarge
    !!     Discharge of the well; positive indicates a well withdrawal
    !!     and negative and injection rate.
    !!
    complex(kind=AE_REAL) :: cZC
    real(kind=AE_REAL) :: rRadius
    real(kind=AE_REAL) :: rDischarge
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
  end type FWL_WELL


  type, public :: FWL_COLLECTION
    !! type FWL_COLLECTION
    !!
    !! PUBLIC type that holds the well entries for a layer
    !!
    !! Members:
    !!   type(FWL_WELL), pointer :: Wells(:)
    !!     FWL_WELL structures that hold well information
    !!   integer :: iCount
    !!     Number of FWL_WELL structures currently in use
    !!
    type(FWL_WELL), dimension(:), pointer :: Wells
    integer(kind=AE_INT) :: iCount
  end type


contains


  function FWL_Create(io, iMax) result(fwl)
    !! function FWL_Create
    !!
    !! Creates a new FWL_COLLECTION object
    !!
    !! Calling Sequence:
    !!    fwl => FWL_Create(iMax)
    !!
    !! Arguments:
    !!    (in)      integer :: iMax
    !!                The maximum number of wells to be stored
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FWL_COLLECTION), pointer :: fwl
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(fwl, stat = iStat)
    call IO_Assert(io, (iStat == 0), "FWL_Create: allocation failed")

    allocate(fwl%Wells(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FWL_Create: allocation failed")
    fwl%iCount = 0

    return
  end function FWL_Create


  subroutine FWL_Alloc(io, fwl, iMax)
    !! subroutine FWL_Alloc
    !!
    !! Dimensions the internal buffers for iaMax wells in layer iL
    !!
    !! Calling Sequence:
    !!    call FWL_Alloc(io, fwl, iMax)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             The FWL_COLLECTION object to be used
    !!   (in)    integer :: iMax
    !!             The maximum number of wells in fwl
    !!
    !! Note: If allocation fails, causes a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Alloc: FWL_New has not been called")
      call IO_Assert(io, (.not. associated(fwl%Wells)), "FWL_Alloc: Wells already allocated")
    end if

    ! Now, allocate space for the specified layer and initialize
    allocate(fwl%Wells(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FWL_Alloc: Allocation failed")
    fwl%Wells = FWL_WELL(cZERO, rZERO, rZERO, -1, -1, -1, -1)
    fwl%iCount = 0

    return
  end subroutine FWL_Alloc


  subroutine FWL_New(io, fwl, cZc, rDischarge, rRadius, iElementType, iElementString, iElementVertex, iElementFlag, iRV)
    !! subroutine FWL_New
    !!
    !! Makes a new well entry. On call, the geometry and discharge of
    !! the well are provided. The internal well structures are then
    !! set up. Returns iRV = the index into the well table on success or
    !! iRV < 0 on failure.
    !!
    !! Calling Sequence:
    !!    call FWL_New(io, fwl, cZC, rDischarge, iRV)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    complex :: cZC
    !!             Complex coordinates of the center of the well
    !!   (in)    real :: rDischarge
    !!             The discharge of the well
    !!   (out)   integer :: iRV
    !!             On success, the index in fwl%Wells used
    !!
    !! Note: On failure, forces a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), intent(in) :: cZc
    real(kind=AE_REAL), intent(in) :: rRadius
    real(kind=AE_REAL), intent(in) :: rDischarge
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(out) :: iRV
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), &
           "FWL_New: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), &
           "FWL_New: FWL_Alloc has not been called")
      call IO_Assert(io, (fwl%iCount < size(fwl%Wells)), &
           "FWL_New: Space exhausted")
    end if

    fwl%iCount = fwl%iCount + 1
    wel => fwl%Wells(fwl%iCount)
    wel = FWL_WELL(cZc, rRadius, rDischarge, iElementType, iElementString, iElementVertex, iElementFlag)
    iRV = fwl%iCount

    return
  end subroutine FWL_New


  subroutine FWL_Update(io, fwl, iWL, rDischarge)
    !! subroutine FWL_Update
    !!
    !! Updates the discharge for a well entry.
    !!
    !! Calling Sequence:
    !!    call FWL_Update(io, fwl, iWL, rDischarge)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    integer :: iWL
    !!             The index for the well in fwl
    !!   (in)    real :: rDischarge
    !!             The discharge of the well
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    integer(kind=AE_INT), intent(in) :: iWL
    real(kind=AE_REAL), intent(in) :: rDischarge
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
      call IO_Assert(io, (iWL <= size(fwl%Wells)), "FWL_Update: Space exhausted")
    end if

    wel => fwl%Wells(iWL)
    wel%rDischarge = rDischarge

    return
  end subroutine FWL_Update


  subroutine FWL_GetInfluence(io, fwl, iWhich, iWL1, iNWL, cPathZ, cOrientation, cF)
    !! subroutine FWL_GetInfluence
    !!
    !! Retrieves arrays of influence functions for use in matrix generation
    !!
    !! Calling Sequence:
    !!    call FWL_GetInfluence(io, fwl, iWhich, iWL1, iNWL, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    integer :: iWhich
    !!             The influence function to be computed;  iWhich values are
    !!                INFLUENCE_P   - Complex potential
    !!                kInfluenceQ   - Complex discharge
    !!                kInfluenceF   - Integrated flux
    !!                kInfluenceG   - Areal infiltration
    !!                kInfluenceQ   - Extraction rate
    !!                kInfluenceD   - Difference in potential
    !!                kInfluenceZ   - All zeroes
    !!   (in)    integer :: iWL1
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL
    !!             The number of consecutive wells to be computed
    !!   (in)    complex :: cPathZ(:)
    !!             Complex coordinates of the control point(s) to be used. For
    !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
    !!             cPathZ(1) is used; for iWhich = INFLUENCE_F, the influence
    !!             function is computed along the path cPathZ(:)
    !!   (in)    complex :: cOrientation
    !!             Orientation normal vector for iWhich = INFLUENCE_W
    !!   (out)   complex :: cF(1:iNWL, 3)
    !!             The returned influence functions.  Indexes 1:iNWL relate
    !!             to well indices iWL1:iWL1+iNWL-1, respectively.
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    integer(kind=AE_INT), intent(in) :: iWhich, iWL1, iNWL
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWL2, i, j
    real(kind=AE_REAL), dimension(1, 1) :: rF
    complex(kind=AE_REAL), dimension(1, 1) :: cW, cUnit
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_GetInfluence: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_GetInfluence: FWL_Alloc has not been called")
      call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
           "FWL_GetInfluence: Bad index range")
      call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
           "FWL_GetInfluence: Bad index range")
      call IO_Assert(io, (size(cF, 1) >= iNWL .and. size(cF, 2) >= 1), "FWL_GetInfluence: Invalid result vector")
    end if

    select case (iWhich)
      case (INFLUENCE_P)
        do i = 1, iNWL
          wel => fwl%Wells(i)
          cF(i, :, :) = cIWL_InfluenceP(cPathZ(1), wel%cZC)
        end do
      case (INFLUENCE_D)
        do i = 1, iNWL
          wel => fwl%Wells(i)
          cF(i, :, :) = cIWL_InfluenceP(cPathZ(1), wel%cZC) - cIWL_InfluenceP(cPathZ(2), wel%cZC)
        end do
      case (INFLUENCE_W)
        do i = 1, iNWL
          wel => fwl%Wells(i)
          cUnit = cOrientation/abs(cOrientation)
          cW = cIWL_InfluenceW(cPathZ(1), wel%cZC)
          cF(i, :, :) = cUnit * conjg(cW)
        end do
      case (INFLUENCE_F)
        do i = 1, iNWL
          wel => fwl%Wells(i)
          rF = rZERO
          do j = 1, ubound(cPathZ, 1)-1
            rF = rF + cIWL_InfluenceF(cPathZ(j), cPathZ(j+1), wel%cZC)
          end do
          cF(i, :, :) = -cmplx(rZERO, rF, AE_REAL)
        end do
      case (INFLUENCE_G)
        cF(i, :, :) = rZERO
      case (INFLUENCE_Q)
        do i = 1, iNWL
          cF(i, :, :) = rONE
        end do
      case (INFLUENCE_Z)
        cF = rZERO
    end select

    return
  end subroutine FWL_GetInfluence


  function cFWL_Potential(io, fwl, cZ, iWL1, iNWL) result(cOmega)
    !! complex function cFWL_Potential
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    cOmega = cFWL_Potential(fwl, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    complex(kind=AE_REAL) :: cOmega
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cP
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
             "FWL_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
               "FWL_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fwl%iCount
    end if

    ! Sum up the contributions of all wells
    cOmega = cZERO
    do i = iStart, iEnd
      wel => fwl%Wells(i)
      cP = cIWL_InfluenceP(cZ, wel%cZC)
      cOmega = cOmega + wel%rDischarge*cP(1, 1)
    end do

    return
  end function cFWL_Potential


  function cFWL_Discharge(io, fwl, cZ, iWL1, iNWL) result(cQ)
    !! complex function cFWL_Discharge
    !!
    !! Computes the complex discharge due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    cQ = cFWL_Discharge(fwl, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    complex(kind=AE_REAL) :: cQ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cW
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
             "FWL_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
               "FWL_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fwl%iCount
    end if

    ! Sum up the contributions of all wells
    cQ = cZERO
    do i = iStart, iEnd
      wel => fwl%Wells(i)
      cW = cIWL_InfluenceW(cZ, wel%cZC)
      cQ = cQ + conjg(wel%rDischarge * cW(1, 1))
    end do

    return
  end function cFWL_Discharge


  function rFWL_Flow(io, fwl, cPathZ, iWL1, iNWL) result(rFlow)
    !! complex function cFWL_Discharge
    !!
    !! Computes the complex discharge due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rFlow = fFWL_Flow(fwl, cPathZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rFlow
    type(IO_STATUS), pointer :: io
    ! Locals
    integer(kind=AE_INT) :: i, j, iStart, iEnd
    real(kind=AE_REAL), dimension(1, 1) :: rS
    complex(kind=AE_REAL) :: cZC, cZL, cZMap
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
             "FWL_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
               "FWL_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fwl%iCount
    end if

    rFlow = rZERO
    ! Sum up the contributions of all wells
    do i = iStart, iEnd
      wel => fwl%Wells(i)
      do j = 1, ubound(cPathZ, 1)-1
        ! SPECIAL CASE: Does the line segment pass within the well bore?
        cZC = rHALF * (cPathZ(j) + cPathZ(j+1))
        cZL = rHALF * (cPathZ(j) - cPathZ(j+1))
        cZMap = (wel%cZC - cZC) / cZL
        if (wel%rRadius >= abs(aimag(cZMap))*abs(cZL)) then
          rS = rZERO
        else
          rS = cIWL_InfluenceF(cPathZ(j), cPathZ(j+1), wel%cZC)
        end if
        rFlow = rFlow + wel%rDischarge*rS(1, 1)
      end do
    end do

    return
  end function rFWL_Flow


  function rFWL_Recharge(io, fwl, cZ, iWL1, iNWL) result(rG)
    !! complex function rFWL_Recharge
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rG = rFWL_Recharge(fwl, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rG
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cG
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
             "FWL_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
               "FWL_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fwl%iCount
    end if

    ! Sum up the contributions of all wells
    rG = cZERO
    do i = iStart, iEnd
      wel => fwl%Wells(i)
      cG = cIWL_InfluenceG(cZ, wel%cZC)
      rG = rG + wel%rDischarge*cG(1, 1)
    end do

    return
  end function rFWL_Recharge


  function rFWL_Extraction(io, fwl, iWL1, iNWL) result(rQ)
    !! complex function rFWL_Extraction
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rQ = rFWL_Extraction(fwl, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rQ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cQ
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(fwl%Wells, 1) .and. iWL1 <= ubound(fwl%Wells, 1)), &
             "FWL_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(fwl%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(fwl%Wells, 1)), &
               "FWL_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fwl%iCount
    end if

    ! Sum up the contributions of all wells
    rQ = cZERO
    do i = iStart, iEnd
      wel => fwl%Wells(i)
      cQ = cIWL_InfluenceQ()
      rQ = rQ + wel%rDischarge*cQ(1, 1)
    end do

    return
  end function rFWL_Extraction


  function lFWL_CheckPoint(io, fwl, cZ, rTol, cZFix, rStrength, iElementType, iElementString, iElementVertex, &
             iElementFlag) result(lFound)
    !! function lFWL_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance.
    !! If rTol > 0, the tolerance is used in the test, otherwise the global constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rTol
    complex(kind=AE_REAL), intent(out) :: cZFix
    real(kind=AE_REAL), intent(out) :: rStrength
    integer(kind=AE_INT), intent(out) :: iElementType
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lFound
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(FWL_WELL), pointer :: wel
    real(kind=AE_REAL), dimension(:), allocatable :: rDist2
    real(kind=AE_REAL) :: rCheckTol

    if (fwl%iCount > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a well
      allocate(rDist2(fwl%iCount))
      rDist2 = (real(cZ)-real(fwl%Wells(1:fwl%iCount)%cZC, AE_REAL))**2 + &
               (aimag(cZ)-aimag(fwl%Wells(1:fwl%iCount)%cZC))**2
      i = minloc(rDist2, DIM = 1)
      wel => fwl%Wells(i)
      if (rDist2(i) < rCheckTol**2 .or. rDist2(i) < wel%rRadius**2) then
        lFound = .true.
        if (rDist2(i) > wel%rRadius**2) then
          cZFix = wel%cZC + rCheckTol * rMOVE_FACTOR
        else
          cZFix = wel%cZC + wel%rRadius * 1.001_AE_REAL
        end if
        rStrength = wel%rDischarge
        iElementType = wel%iElementType
        iElementString = wel%iElementString
        iElementVertex = wel%iElementVertex
        iElementFlag = wel%iElementFlag
      else
        lFound = .false.
        cZFix = cZERO
      end if
      deallocate(rDist2)
    else
      lFound = .false.
      cZFix = cZERO
    end if
    return
  end function lFWL_CheckPoint


  subroutine FWL_Report(io, fwl)
    !! subroutine FWL_Report
    !!
    !! Writes a report of all well information to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call FWL_Report(io, fwl)
    !!
    !! Arguments:
    !!    (in)   type(FWL_COLLECTION), pointer :: fwl
    !!             The FWL_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(FWL_COLLECTION), pointer :: fwl
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(FWL_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(fwl)), "FWL_Update: FWL_Create has not been called")
      call IO_Assert(io, (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called")
    end if

    call HTML_Header('Function module FWL', 1)
    call HTML_Header('Information about well functions', 2)

    if (fwl%iCount > 0) then
      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'XC    ', 'YC    ', 'Radius', 'Q     '/))
      do i = 1, fwl%iCount
        wel => fwl%Wells(i)
        call HTML_StartRow()
        call HTML_ColumnInteger((/i/))
        call HTML_ColumnComplex((/wel%cZC/))
        call HTML_ColumnReal((/wel%rRadius, wel%rDischarge/), 'e13.6')
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No well functions defined', 3)
    end if

    return
  end subroutine FWL_Report

end module f_well
