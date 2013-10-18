module f_reference

  ! ModAEM 2.0
  ! Copyright(c) 1995-2013 Vic Kelson and Layne Christensen Company
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
  ! Contact the author by e-mail at:
  !     victor.kelson@layne.com or vic.kelson@gmail.com
  !
  ! Or by regular mail at:
  ! 	Vic Kelson
  ! 	Chief Modeler
  ! 	Layne Christensen
  ! 	320 W 8th St
  ! 	Bloomington, IN 47401

  !! module f_reference(frf)
  !!
  !! Written by Victor A.Kelson
  !!
  !! Revision History:
  !!   2.0.0   17 October 2013
  !!           New module added for ModAEM-2.0
  !!
  !! Module of data structures and functions for the reference flow field,
  !! that is, the constant of integration and an optional uniform flow discharge.
  !!
  !! Module use:
  !!   constants   --  Universal ModAEM constant declarations
  !!   io          --  Universal ModAEM I/O functions and constants
  !!   i_reference --  Influence function module for the reference flow field
  !!
  use u_constants
  use u_io

  implicit none

  public

	!> A "collection" of one item that provides the reference flow field for a
	!> model layer.
  type, public :: FRF_COLLECTION
    type(IO_STATUS), pointer :: io !< Provides common I/O functions
    complex(kind=AE_REAL) :: cRefPoint !< Location of the reference point
    real(kind=AE_REAL) :: rRefHead !< Specified head at the reference point
    real(kind=AE_REAL) :: rRefPot !< Specified potential at the reference point
    complex(kind=AEM_REAL) :: cUniformFlow !< Uniform flow discharge rate
    real(kind=AE_REAL) :: rStrength !< Computed "strength" of the reference point
    integer(kind=AE_INT) :: iElementType !< Element type for the solver
    integer(kind=AE_INT) :: iElementString !< Element string ID for the solver
    integer(kind=AE_INT) :: iElementVertex !< Element vertex ID for the solver
    integer(kind=AE_INT) :: iElementFlag !< Element flag for the solver
  end type


contains

	!> Creates an FRF_COLLECTION object.
	!> @param[in] io An IO_STATUS object to be used for e.g. exception handling
	!> @return frf An FRF_COLLECTION object
	!> This routine creates the new collection object and zeroes out all the
	!> fields. They are set by other functions within the module.
  function FRF_Create(io) result(frf)
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FRF_COLLECTION), pointer :: frf
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(frf, stat = iStat)
    call IO_Assert(io, (iStat == 0), "FRF_Create: allocation failed")

		frf%cRefPoint = cZERO
		frf%rRefHead = rZERO
		frf%rRefPot = rZERO
		frf%cUniformFlow = cZERO
		frf%rStrength = rZERO
		frf%iElementType = 0


		return
  end function FRF_Create


  subroutine FRF_Alloc(io, frf, iMax)
    !! subroutine FRF_Alloc
    !!
    !! Dimensions the internal buffers for iaMax wells in layer iL
    !!
    !! Calling Sequence:
    !!    call FRF_Alloc(io, frf, iMax)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer
    !!             The FRF_COLLECTION object to be used
    !!   (in)    integer :: iMax
    !!             The maximum number of wells in frf
    !!
    !! Note: If allocation fails, causes a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Alloc: FRF_New has not been called")
      call IO_Assert(io, (.not. associated(frf%Wells)), "FRF_Alloc: Wells already allocated")
    end if

    ! Now, allocate space for the specified layer and initialize
    allocate(frf%Wells(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FRF_Alloc: Allocation failed")
    frf%Wells = FRF_WELL(cZERO, rZERO, rZERO, -1, -1, -1, -1)
    frf%iCount = 0

    return
  end subroutine FRF_Alloc


  subroutine FRF_New(io, frf, cZc, rDischarge, rRadius, iElementType, iElementString, iElementVertex, iElementFlag, iRV)
    !! subroutine FRF_New
    !!
    !! Makes a new well entry. On call, the geometry and discharge of
    !! the well are provided. The internal well structures are then
    !! set up. Returns iRV = the index into the well table on success or
    !! iRV < 0 on failure.
    !!
    !! Calling Sequence:
    !!    call FRF_New(io, frf, cZC, rDischarge, iRV)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
    !!   (in)    complex :: cZC
    !!             Complex coordinates of the center of the well
    !!   (in)    real :: rDischarge
    !!             The discharge of the well
    !!   (out)   integer :: iRV
    !!             On success, the index in frf%Wells used
    !!
    !! Note: On failure, forces a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
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
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), &
           "FRF_New: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), &
           "FRF_New: FRF_Alloc has not been called")
      call IO_Assert(io, (frf%iCount < size(frf%Wells)), &
           "FRF_New: Space exhausted")
    end if

    frf%iCount = frf%iCount + 1
    wel => frf%Wells(frf%iCount)
    wel = FRF_WELL(cZc, rRadius, rDischarge, iElementType, iElementString, iElementVertex, iElementFlag)
    iRV = frf%iCount

    return
  end subroutine FRF_New


  subroutine FRF_Update(io, frf, iWL, rDischarge)
    !! subroutine FRF_Update
    !!
    !! Updates the discharge for a well entry.
    !!
    !! Calling Sequence:
    !!    call FRF_Update(io, frf, iWL, rDischarge)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), frf
    !!             The FRF_COLLECTION object to be used
    !!   (in)    integer :: iWL
    !!             The index for the well in frf
    !!   (in)    real :: rDischarge
    !!             The discharge of the well
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    integer(kind=AE_INT), intent(in) :: iWL
    real(kind=AE_REAL), intent(in) :: rDischarge
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
      call IO_Assert(io, (iWL <= size(frf%Wells)), "FRF_Update: Space exhausted")
    end if

    wel => frf%Wells(iWL)
    wel%rDischarge = rDischarge

    return
  end subroutine FRF_Update


  subroutine FRF_GetInfluence(io, frf, iWhich, iWL1, iNWL, cPathZ, cOrientation, cF)
    !! subroutine FRF_GetInfluence
    !!
    !! Retrieves arrays of influence functions for use in matrix generation
    !!
    !! Calling Sequence:
    !!    call FRF_GetInfluence(io, frf, iWhich, iWL1, iNWL, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
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
    type(FRF_COLLECTION), pointer :: frf
    integer(kind=AE_INT), intent(in) :: iWhich, iWL1, iNWL
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWL2, i, j
    real(kind=AE_REAL), dimension(1, 1) :: rF
    complex(kind=AE_REAL), dimension(1, 1) :: cW, cUnit
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_GetInfluence: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_GetInfluence: FRF_Alloc has not been called")
      call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
           "FRF_GetInfluence: Bad index range")
      call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
           "FRF_GetInfluence: Bad index range")
      call IO_Assert(io, (size(cF, 1) >= iNWL .and. size(cF, 2) >= 1), "FRF_GetInfluence: Invalid result vector")
    end if

    select case (iWhich)
      case (INFLUENCE_P)
        do i = 1, iNWL
          wel => frf%Wells(i)
          cF(i, :, :) = cIWL_InfluenceP(cPathZ(1), wel%cZC)
        end do
      case (INFLUENCE_D)
        do i = 1, iNWL
          wel => frf%Wells(i)
          cF(i, :, :) = cIWL_InfluenceP(cPathZ(1), wel%cZC) - cIWL_InfluenceP(cPathZ(2), wel%cZC)
        end do
      case (INFLUENCE_W)
        do i = 1, iNWL
          wel => frf%Wells(i)
          cUnit = cOrientation/abs(cOrientation)
          cW = cIWL_InfluenceW(cPathZ(1), wel%cZC)
          cF(i, :, :) = cUnit * conjg(cW)
        end do
      case (INFLUENCE_F)
        do i = 1, iNWL
          wel => frf%Wells(i)
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
  end subroutine FRF_GetInfluence


  function cFRF_Potential(io, frf, cZ, iWL1, iNWL) result(cOmega)
    !! complex function cFRF_Potential
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    cOmega = cFRF_Potential(frf, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
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
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    complex(kind=AE_REAL) :: cOmega
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cP
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
             "FRF_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
               "FRF_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = frf%iCount
    end if

    ! Sum up the contributions of all wells
    cOmega = cZERO
    do i = iStart, iEnd
      wel => frf%Wells(i)
      cP = cIWL_InfluenceP(cZ, wel%cZC)
      cOmega = cOmega + wel%rDischarge*cP(1, 1)
    end do

    return
  end function cFRF_Potential


  function cFRF_Discharge(io, frf, cZ, iWL1, iNWL) result(cQ)
    !! complex function cFRF_Discharge
    !!
    !! Computes the complex discharge due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    cQ = cFRF_Discharge(frf, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
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
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    complex(kind=AE_REAL) :: cQ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cW
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
             "FRF_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
               "FRF_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = frf%iCount
    end if

    ! Sum up the contributions of all wells
    cQ = cZERO
    do i = iStart, iEnd
      wel => frf%Wells(i)
      cW = cIWL_InfluenceW(cZ, wel%cZC)
      cQ = cQ + conjg(wel%rDischarge * cW(1, 1))
    end do

    return
  end function cFRF_Discharge


  function rFRF_Flow(io, frf, cPathZ, iWL1, iNWL) result(rFlow)
    !! complex function cFRF_Discharge
    !!
    !! Computes the complex discharge due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rFlow = fFRF_Flow(frf, cPathZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
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
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rFlow
    type(IO_STATUS), pointer :: io
    ! Locals
    integer(kind=AE_INT) :: i, j, iStart, iEnd
    real(kind=AE_REAL), dimension(1, 1) :: rS
    complex(kind=AE_REAL) :: cZC, cZL, cZMap
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
             "FRF_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
               "FRF_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = frf%iCount
    end if

    rFlow = rZERO
    ! Sum up the contributions of all wells
    do i = iStart, iEnd
      wel => frf%Wells(i)
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
  end function rFRF_Flow


  function rFRF_Recharge(io, frf, cZ, iWL1, iNWL) result(rG)
    !! complex function rFRF_Recharge
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rG = rFRF_Recharge(frf, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
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
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rG
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cG
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
             "FRF_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
               "FRF_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = frf%iCount
    end if

    ! Sum up the contributions of all wells
    rG = cZERO
    do i = iStart, iEnd
      wel => frf%Wells(i)
      cG = cIWL_InfluenceG(cZ, wel%cZC)
      rG = rG + wel%rDischarge*cG(1, 1)
    end do

    return
  end function rFRF_Recharge


  function rFRF_Extraction(io, frf, iWL1, iNWL) result(rQ)
    !! complex function rFRF_Extraction
    !!
    !! Computes the complex potential due to the current set of wells.
    !!
    !! Calling Sequence:
    !!    rQ = rFRF_Extraction(frf, cZ, iWL1, iNWL)
    !!
    !! Arguments:
    !!   (in)    type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
    !!   (in)    integer :: iWL1 [OPTIONAL]
    !!             The index for the first well to be used
    !!   (in)    integer :: iNWL [OPTIONAL]
    !!             The number of consecutive wells to be computed
    !! Note:
    !!   If iWL1 is not provided, all wells will be used
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    integer(kind=AE_INT), intent(in), optional :: iWL1, iNWL
    real(kind=AE_REAL) :: rQ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cQ
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iWL1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iWL1 >= lbound(frf%Wells, 1) .and. iWL1 <= ubound(frf%Wells, 1)), &
             "FRF_GetInfluence_ILS: Bad index range")
      end if
      iStart = iWL1
      if (present(iNWL)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iWL1+iNWL-1) >= lbound(frf%Wells, 1) .and. (iWL1+iNWL-1) <= ubound(frf%Wells, 1)), &
               "FRF_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNWL-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = frf%iCount
    end if

    ! Sum up the contributions of all wells
    rQ = cZERO
    do i = iStart, iEnd
      wel => frf%Wells(i)
      cQ = cIWL_InfluenceQ()
      rQ = rQ + wel%rDischarge*cQ(1, 1)
    end do

    return
  end function rFRF_Extraction


  function lFRF_CheckPoint(io, frf, cZ, rTol, cZFix, rStrength, iElementType, iElementString, iElementVertex, &
             iElementFlag) result(lFound)
    !! function lFRF_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance.
    !! If rTol > 0, the tolerance is used in the test, otherwise the global constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
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
    type(FRF_WELL), pointer :: wel
    real(kind=AE_REAL), dimension(:), allocatable :: rDist2
    real(kind=AE_REAL) :: rCheckTol

    if (frf%iCount > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a well
      allocate(rDist2(frf%iCount))
      rDist2 = (real(cZ)-real(frf%Wells(1:frf%iCount)%cZC, AE_REAL))**2 + &
               (aimag(cZ)-aimag(frf%Wells(1:frf%iCount)%cZC))**2
      i = minloc(rDist2, DIM = 1)
      wel => frf%Wells(i)
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
  end function lFRF_CheckPoint


  subroutine FRF_Report(io, frf)
    !! subroutine FRF_Report
    !!
    !! Writes a report of all well information to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call FRF_Report(io, frf)
    !!
    !! Arguments:
    !!    (in)   type(FRF_COLLECTION), pointer :: frf
    !!             The FRF_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(FRF_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(frf)), "FRF_Update: FRF_Create has not been called")
      call IO_Assert(io, (associated(frf%Wells)), "FRF_Update: FRF_Alloc has not been called")
    end if

    call HTML_Header('Function module FRF', 1)
    call HTML_Header('Information about well functions', 2)

    if (frf%iCount > 0) then
      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'XC    ', 'YC    ', 'Radius', 'Q     '/))
      do i = 1, frf%iCount
        wel => frf%Wells(i)
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
  end subroutine FRF_Report

end module f_well
