module m_wl1

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

  !! module m_wl1
  !!
  !! Element module for 2-D head specified wells
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!

  use u_constants
  use u_io
  use f_well
  use f_bwl
  use u_matrix
  use m_aqu

  implicit none

  public

  type :: WL1_WELL
    !! type WL1_WELL
    !!
    !! Type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: cZHead
    !!     The point at which the head is specified
    !!   real :: rHead
    !!     The given head
    !!   real :: rHeadCorrection
    !!     The correction on the head as computed by another module(for derived packages)
    !!   integer :: iID
    !!     Identification label for the well(for interaction with e.g. GUIs)
    !!   integer :: iFWLIndex
    !!     Index for the well entry in the FWL module
    !!   real :: rStrength
    !!     Computed strength of the element    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rRadius
    complex(kind=AE_REAL) :: cZHead
    real(kind=AE_REAL) :: rHead
    real(kind=AE_REAL) :: rHeadCorrection
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iFWLIndex
    real(kind=AE_REAL) :: rStrength
    real(kind=AE_REAL) :: rCheckPot
    real(kind=AE_REAL) :: rDFHead
    real(kind=AE_REAL) :: rError
    ! Bessel well components for partially-penetrating wells
    logical :: lPpWell
    real(kind=AE_REAL) :: rScrBot
    real(kind=AE_REAL) :: rScrTop
    real(kind=AE_REAL) :: rKhKv
    real(kind=AE_REAL) :: rWtblHead
    real(kind=AE_REAL) :: rScrHead
    type(BWL_WELL), pointer :: bwl
  end type WL1_WELL

  type :: WL1_COLLECTION
    !! type WL1_COLLECTION
    !!
    !! Type that holds head-specified wells in a layer
    !!
    !! Members:
    !!   type(WL1_WELL), dimension(:), pointer :: Wells
    !!     Array of WL1_WELL objects for the layer; dimensioned for the maximum
    !!     number of wells according to the input file(see WL1_Read)
    !!   integer :: iCount
    !!     The number of strings actually in use
    !!
    type(WL1_WELL), dimension(:), pointer :: Wells
    integer(kind=AE_INT) :: iCount
    ! Iterator information
    integer(kind=AE_INT) :: iIterWell
    ! Debug flag
    logical :: lDebug
  end type WL1_COLLECTION

  ! Module flags for matrix generator routines
  integer(kind=AE_INT), private, parameter :: WL1Vertex = 1


contains


  function WL1_Create(io) result(wl1)
    !! function WL1_Create
    !!
    !! Creates a new WL1_COLLECTION object
    !!
    !! Calling Sequence:
    !!    WL1 => WL1_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(wl1, stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL1_Create: allocation failed")
    nullify(wl1%Wells)
    wl1%iCount = 0

    return
  end function WL1_Create


  subroutine WL1_Alloc(io, wl1)
    !! Subroutine WL1_Alloc
    !!
    !! Allocates Wells for the WL1_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call WL1_Alloc(io, wl1, iNStr)
    !!
    !! Arguments:
    !!    (in)    type(WL1_COLLECTION), pointer :: wl1
    !!              The WL1_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNWL
    integer(kind=AE_INT) :: iStat

    iNWL = iIO_GetInteger(io, 'iNWL', minimum = 0)
    allocate(wl1%Wells(iNWL), stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL1_Alloc: allocation failed")

    return
  end subroutine WL1_Alloc



  subroutine WL1_Destroy(io, wl1)
    !! subroutine WL1_Destroy
    !!
    !! Frees memory allocated for wl1 Wells and the WL1 Collection object
    !!
    !! Calling Sequence:
    !!     call WL1_Destroy(wl1)
    !!
    !! Arguments:
    !!  type(WL1_COLLECTION), pointer :: wl1
    !!              Pointer to the WL1_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat


    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Destroy: WL1_Create has not been called")
    end if

    if (associated(wl1%Wells)) then
      deallocate(wl1%Wells, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "WL1_Destroy: deallocation of Wells failed")
    end if
    deallocate(wl1, stat = iStat)
    call IO_Assert(io, (iStat == 0), "WL1_Destroy: deallocation failed")

    return
  end subroutine WL1_Destroy



  subroutine WL1_New(io, wl1, Well)
    !! function WL1_New
    !!
    !! Adds a new WL1_WELL object to the WL1_COLLECTION 'wl1'
    !!
    !! Calling Sequence:
    !!    call WL1_New(io, wl1, Well)
    !!
    !! Arguments:
    !!    (in)    type(WL1_COLLECTION), pointer :: wl1
    !!              The WL1_COLLECTION object to be used
    !!    (in)    type(WL1_WELL), pointer :: Well
    !!              The new well object
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(WL1_WELL) :: Well
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_New: WL1_Create has not been called")
    end if

    call IO_Assert(io, (wl1%iCount < size(wl1%Wells)), &
         "WL1_New: Space exhausted")

    wl1%iCount = wl1%iCount + 1
    wl1%Wells(wl1%iCount) = Well

    return
  end subroutine WL1_New


  function iWL1_GetID(io, wl1, iIndex) result(iID)
    !! Returns the ID number for the well at index 'iIndex'
    type(WL1_COLLECTION), pointer :: wl1
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= wl1%iCount), "Internal error -- no such index")
    iID = wl1%Wells(iIndex)%iID

    return
  end function iWL1_GetID


  subroutine WL1_PreSolve(io, wl1)
    !! subroutine WL1_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call WL1_PreSolve(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io

    return
  end subroutine WL1_PreSolve


  function iWL1_GetInfo(io, wl1, iOption, iIteration) result(iValue)
    !! function WL1_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iWL1_GetInfo(io, wl1, iOption)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_GetInfo: WL1_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = wl1%iCount
      case (SIZE_EQUATIONS)
        iValue = wl1%iCount
      case (SIZE_UNKNOWNS)
        iValue = wl1%iCount
      case (INFO_REGENERATE)
        if (iIteration < 2) then
          iValue = 1
        else
          iValue = 0
        end if
      case default
        iValue = wl1%iCount
    end select

    return
  end function iWL1_GetInfo


  subroutine WL1_SetupFunctions(io, wl1, fwl, aqu)
    !! subroutine WL1_SetupFunctions
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call WL1_Setup(io, wl1, fwl, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(FWL_COLLECTION), pointer :: fwl
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iWel, iWL
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(WL1_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Setup: WL1_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "WL1_Setup: Illegal FWL_COLLECTION object")
    end if

    do iWel = 1, wl1%iCount
      wel => wl1%Wells(iWel)

      call FWL_New(io, fwl, wel%cZ, rZERO, wel%rRadius, ELEM_WL1, iWel, -1, -1, iWL)
      wel%iFWLIndex = iWL
      if (wel%lPpWell) then
        wel%bwl => BWL_New(io, wel%cZ, wel%rRadius, wel%rScrBot, wel%rScrTop, &
                           rAQU_Base(io, aqu, wel%cZ), &
                           rAQU_Base(io, aqu, wel%cZ) + &
                             rAQU_SatdThickness(io, aqu, wel%cZ, real(cAQU_Potential(io, aqu, wel%cZ))), &
                           rAQU_HydCond(io, aqu, wel%cZ), &
                           wel%rKhKv)
        wel%cZHead = wel%cZ + wel%rRadius
      end if
    end do

    return
  end subroutine WL1_SetupFunctions


  subroutine WL1_SetupMatrix(io, wl1, aqu, mat)
    !! subroutine WL1_SetupMatrix
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call WL1_Setup(io, wl1, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iWell, iWL, iEQ
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(WL1_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Setup: WL1_Create has not been called")
      call IO_Assert(io, (associated(aqu)), &
           "WL1_Setup: Illegal AQU_COLLECTION object")
      call IO_Assert(io, (associated(mat)), &
           "WL1_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix generator entries for all segments
    do iWell = 1, wl1%iCount
      ! Set up the unknown variables
      wel => wl1%Wells(iWell)
      call MAT_CreateVariable(io, mat, ELEM_WL1, iWell, 0, WL1Vertex)

      iEQ = MAT_CreateEquation(io, mat, (/wel%cZHead/), EQN_HEAD, ELEM_WL1, 0, iWell, 0, cZERO, rZERO)
    end do

    return
  end subroutine WL1_SetupMatrix


  function iWL1_Prepare(io, wl1, aqu, iIteration) result(iChanges)
    !! subroutine WL1_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_wl1
    !!
    !! Calling Sequence:
    !!    call WL1_Setup(io, wl1, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(WL1_WELL), pointer :: wel
    integer(kind=AE_INT) :: iChanges, iwel
    logical :: lChange
    real(kind=AE_REAL) :: rAqTop

    iChanges = 0
    do iwel = 1, wl1%iCount
      wel => wl1%Wells(iwel)
      if (wel%lPpWell) then
        rAqTop = rAQU_Base(io, aqu, wel%cZ) + &
                 rAQU_SatdThickness(io, aqu, wel%cZ, real(cAQU_Potential(io, aqu, wel%cZ)))
        call BWL_Solve(io, wel%bwl, rAqTop, lChange)
        if (lChange) iChanges = iChanges+1
      end if
    end do

    return
  end function iWL1_Prepare


  function rWL1_GetCoefficientMultiplier(io, wl1, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier
    type(WL1_WELL), pointer :: wel

    ! Make the correction for the partially-penetrating well
    wel => wl1%Wells(iElementVertex)
    if (wel%lPpWell ) then
      rMultiplier = BWL_GetCoefficientMultiplier(io, wel%bwl)
    else
      rMultiplier = rONE
    end if

    return
  end function rWL1_GetCoefficientMultiplier


  subroutine WL1_ComputeCoefficients(io, wl1, fwl, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine WL1_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the WL1
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call WL1_ComputeCoefficients(io, wl1, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    complex :: cPathZ(:)
    !!             The control point(or control path) to be used
    !!   (in)    integer :: iEqType
    !!             The equation type
    !!   (in)    integer :: iElementType
    !!             The element type that created this equation
    !!   (in)    integer :: iElementString
    !!             The element string corresponding to this equation
    !!   (in)    integer :: iElementVertex
    !!             The element vertex corresponding to this equation
    !!   (in)    integer :: iElementFlag
    !!             The element flag(if any) for this equation
    !!   (in)    complex :: cOrientation
    !!             Orientation unit vector(for discharge-based equations)
    !!   (out)   real :: rARow(:)
    !!             The output row of coefficients(to be concatenated with
    !!             row portions for other element modules.
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(FWL_COLLECTION), pointer :: fwl
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    real(kind=AE_REAL), intent(in) :: rGhbResistance
    integer(kind=AE_INT), intent(in) :: iIteration
    real(kind=AE_REAL), intent(in) :: rMultiplier
    real(kind=AE_REAL), dimension(:), intent(out) :: rARow
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iCol, iVtx, iWL1, iNWL, iWhich
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cWLF, cWLW
    type(WL1_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_ComputeCoefficients: WL1_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "WL1_Setup: Illegal FWL_COLLECTION object")
    end if

    if (wl1%iCount > 0) then

      iCol = 0
      rARow = rZERO
      ! ASSUMES that WL1_Setup routine created consecutive well entries
      wel => wl1%Wells(1)
      iWL1 = wel%iFWLIndex
      iNWL = wl1%iCount
      allocate(cWLF(1:iNWL, 1, 1), cWLW(1:iNWL, 1, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "WL1_ComputeCoefficients: Allocation failed")

      ! Get the appropriate incluence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FWL_GetInfluence(io, fwl, INFLUENCE_P, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_BDYGHB)
          call FWL_GetInfluence(io, fwl, INFLUENCE_P, iWL1, iNWL, (/rHALF*sum(cPathZ)/), cOrientation, cWLF(1:iNWL, :, :))
          call FWL_GetInfluence(io, fwl, INFLUENCE_F, iWL1, iNWL, cPathZ, cOrientation, cWLW(1:iNWL, :, :))
          cWLF = cWLF - rGhbResistance * cWLW
        case (EQN_FLOW)
          call FWL_GetInfluence(io, fwl, INFLUENCE_F, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_INHO)
          call FWL_GetInfluence(io, fwl, INFLUENCE_P, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_DISCHARGE)
          call FWL_GetInfluence(io, fwl, INFLUENCE_W, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_RECHARGE)
          call FWL_GetInfluence(io, fwl, INFLUENCE_G, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_CONTINUITY)
          call FWL_GetInfluence(io, fwl, INFLUENCE_Q, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_POTENTIALDIFF)
          call FWL_GetInfluence(io, fwl, INFLUENCE_D, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
        case (EQN_TOTALFLOW)
          call FWL_GetInfluence(io, fwl, INFLUENCE_Z, iWL1, iNWL, cPathZ, cOrientation, cWLF(1:iNWL, :, :))
      end select

      ! Compute the matrix coefficients
      do iVtx = 1, iNWL
        iCol = iCol+1
        rARow(iCol) = real(cWLF(iVtx, 1, 1))
        ! Make the correction for the partially-penetrating well
        wel => wl1%Wells(iVtx)
        if (iElementType == ELEM_WL1 .and. iElementVertex == iVtx .and. wel%lPpWell ) then
          rARow(iCol) = rARow(iCol) + BWL_ComputeCoefficient(io, wel%bwl)
        end if
      end do

      deallocate(cWLF, cWLW)
    end if

    rARow = rARow * rMultiplier

    return
  end subroutine WL1_ComputeCoefficients


  function rWL1_ComputeRHS(io, wl1, aqu, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rWL1_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rWL1_ComputeRHS(io, wl1, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag, rSpecValue)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION object to be used
    !!   (in)    integer :: iElementType
    !!             Element type(either ELAM_AQU or ELEM_IN0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!   (in)    real :: rSpecValue
    !!             The new result value from the solution vector
    !!
    !! Return Value:
    !!   real :: rRHS
    !!     The RHS value for the module
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRHS
    ! [ LOCALS ]
    type(WL1_WELL), pointer :: wel

    ! Set up the unknown variables
    wel => wl1%Wells(iElementVertex)
    if (lDirect) then
      rRHS = rAQU_HeadToPotential(io, aqu, wel%rHead, wel%cZHead)
    else
      rRHS = rAQU_HeadToPotential(io, aqu, wel%rHead, wel%cZHead) - wel%rCheckPot
    end if

    ! Make the correction for the partially-penetrating well?
    if (wel%lPpWell ) then
      rRHS = (rRHS - BWL_ComputeRHS(io, wel%bwl, wel%rStrength)) * BWL_GetCoefficientMultiplier(io, wel%bwl)
    end if

    return
  end function rWL1_ComputeRHS


  subroutine WL1_StoreResult(io, wl1, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine WL1_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the WL1 module.
    !!
    !! Calling Sequence:
    !!    WL1_StoreResult(io, wl1, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(always ELEM_WL1)
    !!   (in)    integer :: iElementWell
    !!             Element well number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(WL1_WELL), pointer :: wel

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_StoreResult: WL1_Create has not been called")
      call IO_Assert(io, (iElementString >= 1 .and. iElementString <= wl1%iCount), &
           "WL1_StoreResult: Bad element well ID")
    end if

    wel => wl1%Wells(iElementString)
    if (lDirect) then
      wel%rStrength = rValue
    else
      wel%rStrength = wel%rStrength + rValue
    end if

    return
  end subroutine WL1_StoreResult


  subroutine WL1_Update(io, wl1, fwl)
    !! subroutine WL1_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    WL1_Update(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(FWL_COLLECTION), pointer :: fwl
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWell
    type(WL1_WELL), pointer :: wel


    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Update: WL1_Create has not been called")
    end if

    do iWell = 1, wl1%iCount
      wel => wl1%Wells(iWell)
      call FWL_Update(io, fwl, wel%iFWLIndex, wel%rStrength)
    end do

    return
  end subroutine WL1_Update


  subroutine WL1_ResetIterator(io, wl1)
    !! subroutine WL1_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call WL1_ResetIterator(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: wl1
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_ResetIterator: WL1_Create has not been called")
    end if

    wl1%iIterWell = 0

    return
  end subroutine WL1_ResetIterator


  function WL1_NextIterator(io, wl1) result(itr)
    !! function WL1_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call WL1_NextIterator(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: wl1
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_NextIterator: WL1_Create has not been called")
    end if

    wl1%iIterWell = wl1%iIterWell + 1
    if (wl1%iIterWell > wl1%iCount) then
      nullify(itr)
    else
      allocate(itr)
      itr%iValueSelector = VALUE_POTENTIAL
      itr%iElementString = wl1%iIterWell
      allocate(itr%cZ(1))
      itr%cZ(1) = wl1%Wells(wl1%iIterWell)%cZHead
    end if

    return
  end function WL1_NextIterator


  subroutine WL1_SetIterator(io, wl1, aqu, itr, cValue)
    !! function WL1_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call WL1_SetIterator(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION to be used
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: wl1
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(AQU_COLLECTION), pointer :: aqu
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    real(kind=AE_REAL) :: rPot
    type(WL1_WELL), pointer :: wel
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_NextIterator: WL1_Create has not been called")
      call IO_Assert(io, (wl1%iIterWell <= wl1%iCount), &
           "WL1_SetIterator: Iterator out of range")
    end if

    rPot = real(cValue, AE_REAL)
    wel => wl1%Wells(itr%iElementString)
    wel%rCheckPot = rPot
    wel%rDFHead = rAQU_PotentialToHead(io, aqu, rPot, wel%cZHead)
    if ( wel%lPpWell ) then
        wel%rError = wel%rDFHead - &
                     wel%rStrength * wel%bwl%rInfl(wel%bwl%iLayer) - &
                     wel%rHead
    else
        wel%rError = wel%rDFHead - wel%rHead
    end if

    return
  end subroutine WL1_SetIterator


  subroutine WL1_Read(io, wl1)
    !! subroutine WL1_Read
    !!
    !! Populates an WL1_COLLECTION using data from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call WL1_Read(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer :: wl1
    !!             WL1_COLLECTION to be populated
    !!
    !! The format of the WL1 section of the input file appears as follows:
    !! WL1 NWells
    !!     x y r x y head id
    !!     ... Up to NWells
    !! END
    !!
    !! NOTE: It is assumed that the WL1 line was found by the caller

    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(2), parameter :: dirDirectives = (/ dirEND, dirPPW /)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rRad, rHead, rScrBot, rScrTop, rKhKv
    complex(kind=AE_REAL) :: cZ, cZHead
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iID
    logical :: lFlag
    type(WL1_WELL), pointer :: wel

    call IO_MessageText(io, "  Reading WL1 module input")

    call IO_Assert(io, (associated(wl1)), "WL1_Read: WL1_Create has not been called")

    ! Process input
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "WL1_Read: I/O Error")
          exit
        case (kOpFileEOF)
          ! EOF is unexpected for all ModWL1 "ifXXXRead" routines.
          ! Report the condition, but proceed as if EOD was found.
          call IO_Assert(io, .false., "WL1_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(wl1%Wells)), "WL1_Read: No space allocated")
          call IO_Assert(io, (wl1%iCount < size(wl1%Wells)), "WL1_Read: Space exhausted")
          ! Oh, all right...  Retrieve the data and put it in there.
          cZ = cIO_GetCoordinate(io, "cZ", extents=.true.)
          rRad = rIO_GetReal(io, "rRad", minimum = rTINY)
          cZHead = cIO_GetCoordinate(io, "cZHead", extents=.true.)
          rHead = rIO_GetReal(io, "rHead")
          iID = iIO_GetInteger(io, "iID")
          wl1%iCount = wl1%iCount+1
          wel => wl1%Wells(wl1%iCount)
          wel%cZ = cZ
          wel%rRadius = rRad
          wel%cZHead = cZHead
          wel%rHead = rHead
          wel%rHeadCorrection = rZERO
          wel%iID = iID
          ! No FWL is declared here; see WL1_Setup
          wel%iFWLIndex = -1
          wel%rStrength = rZERO
        case (kOpPPW)
          ! Partially-penetrating well info follows
          call IO_Assert(io, (wl1%iCount>0), "No current well to modify")
          wel => wl1%Wells(wl1%iCount)
          rScrBot = rIO_GetReal(io, 'rScrBot')
          rScrTop = rIO_GetReal(io, 'rScrTop', minimum=rScrBot+rTINY)
          rKhKv = rIO_GetReal(io, 'rKhKv', def=rONE)
          wel%lPpWell = .true.
          wel%rScrBot = rScrBot
          wel%rScrTop = rScrTop
          wel%rKhKv = rKhKv
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
      end select
    end do

    call IO_MessageText(io, "  Leaving WL1 module")

  end subroutine WL1_Read


  subroutine WL1_Inquiry(io, wl1, aqu, iLU)
    !! subroutine WL1_Inquiry
    !!
    !! Writes an inquiry report for all line-sinks to iLU
    !!
    !! Calling Sequence:
    !!    call WL1_Inquiry(io, wl1, iLU)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(WL1_WELL), pointer :: wel
    real(kind=AE_REAL) :: pp_wtbl, pp_head

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Inquiry: WL1_Create has not been called")
    end if


    write (unit=iLU, &
           fmt="(""#WL1, ID, X, Y, RADIUS, CP_X, CP_Y, SPEC_HEAD, DISCHARGE, MOD_HEAD, DF_HEAD, PP_WTBL, PP_HEAD, ERROR"")")
    do i = 1, wl1%iCount
      wel => wl1%Wells(i)
      if (wel%lPpWell) then
        pp_wtbl = wel%rDFHead - wel%rStrength * wel%bwl%rInfl(1)
        pp_head = wel%rDFHead - wel%rStrength * wel%bwl%rInfl(wel%bwl%iLayer)
      else
        pp_wtbl = rHUGE
        pp_head = rHUGE
      end if
      write (unit=iLU, &
             fmt="(""WL1"", 1("", "", i9), 8("", "", e16.8))" &
             ) wel%iID, &
             cIO_WorldCoords(io, wel%cZ), &
             wel%rRadius, &
             cIO_WorldCoords(io, wel%cZHead), &
             wel%rHead, &
             wel%rStrength, &
             rAQU_PotentialToHead(io, aqu, wel%rCheckPot, wel%cZHead), &
             wel%rDFHead, &
             pp_wtbl, &
             pp_head, &
             wel%rError
    end do

    return
  end subroutine WL1_Inquiry


  subroutine WL1_Report(io, wl1)
    !! subroutine WL1_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call WL1_Report(wl1)
    !!
    !! Arguments:
    !!   (in)    type(WL1_COLLECTION), pointer
    !!             WL1_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(WL1_COLLECTION), pointer :: wl1
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(WL1_WELL), pointer :: wel
    integer(kind=AE_INT) :: i
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    real(kind=AE_REAL) :: pp_wtbl, pp_head

    if (io%lDebug) then
      call IO_Assert(io, (associated(wl1)), &
           "WL1_Report: WL1_Create has not been called")
    end if

    call HTML_Header('Module WL1', 1)
    call HTML_Header('Head-specified well information', 2)

    if (associated(wl1%Wells)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of wells', wl1%iCount)
      call HTML_AttrInteger('Number of FWL functions', iWL1_GetInfo(io, wl1, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iWL1_GetInfo(io, wl1, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iWL1_GetInfo(io, wl1, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iWL1_GetInfo(io, wl1, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iWL1_GetInfo(io, wl1, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      call HTML_Header('Wells', 4)

      call HTML_StartTable()
      call HTML_TableHeader((/'Well   ', 'ID     ', 'FWL #  ', 'X      ', 'Y      ', 'Head   ', 'R      ', &
           'Str    ', 'DF Head', 'PP Wtbl', 'PP Head', 'Error  '/))
      do i = 1, wl1%iCount
        wel => wl1%Wells(i)
        if (wel%lPpWell) then
          pp_wtbl = wel%rDFHead - wel%rStrength * wel%bwl%rInfl(1)
          pp_head = wel%rDFHead - wel%rStrength * wel%bwl%rInfl(wel%bwl%iLayer)
        else
          pp_wtbl = rHUGE
          pp_head = rHUGE
        end if
        call HTML_StartRow()
        call HTML_ColumnInteger((/i, wel%iID, wel%iFWLIndex/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, wel%cZ)/))
        call HTML_ColumnReal((/wel%rHead, wel%rRadius, wel%rStrength, wel%rDFHead, pp_wtbl, pp_head, wel%rError/))
        call HTML_EndRow()
      end do
      call HTML_EndTable()

    else
      call HTML_Header('No wells defined', 3)
    end if

    return
  end subroutine WL1_Report


  subroutine WL1_Save(io, wl1, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL1_COLLECTION), pointer :: wl1
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iwel
    type(WL1_WELL), pointer :: wel

    ! Output records will be of the form ELEM_WL1, IWEL, IRAD, IVTX, 0, SIGMA
    do iwel = 1, wl1%iCount
      wel => wl1%Wells(iwel)
      if (mode == IO_MODE_BINARY) then
        write (unit=LU_SCRATCH) ELEM_WL1, iwel, 1, 1, wel%rStrength
      else
        write (unit=LU_SCRATCH, fmt=*) "WL1", iwel, 1, 1, wel%rStrength
      end if
    end do

    return
  end subroutine WL1_Save


  subroutine WL1_Load(io, wl1, fwl, mode)
    !! Loads the WL1 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(WL1_COLLECTION), pointer :: wl1
    type(FWL_COLLECTION), pointer :: fwl
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, iwel, ivtx, iflg, istat
    real(kind=AE_REAL) :: rwelength
    character(len=3) :: smodule
    type(WL1_WELL), pointer :: wel

    ! Scans the entire precondition file for the WL1 data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, iwel, ivtx, iflg, rwelength
        if (imodule /= ELEM_WL1) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, iwel, ivtx, iflg, rwelength
        if (uppercase (trim(smodule)) /= "WL1") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      call IO_Assert(io, iwel > 0 .and. iwel <= wl1%iCount, "WL1 well not found")
      wel => wl1%Wells(iwel)
      call IO_Assert(io, ivtx == 1, "WL1 vertex not found")
      call IO_Assert(io, iflg == 1, "WL1 welength index not found")
      wel%rStrength = rwelength
    end do

    ! Now, populate the internal data weluctures
    call WL1_Update(io, wl1, fwl)

    return
  end subroutine WL1_Load

end module m_wl1
