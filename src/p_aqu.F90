module p_aqu

  ! ModAEM 1.9
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

  !! module m_aqu
  !!
  !! Element module for 2-D uniform(optionally bounded) aquifers
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_well     --  Function module for collections of wells
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for the overlying
  !! 2-D aquifer that works together with the other element modules
  !!
  !! The AQU module is a very special case in ModAEM.  It has three different
  !! "personalities", all in one module:
  !!     (a)  An ELEMENT MODULE that sets up equations to be solved
  !!     (b)  A FUNCTION MODULE that provides the low-level function s(io, like "f_aqu")
  !!     (c)  A UTILITY MODULE that allows for convenient AEM calculations e.g.
  !!          conversion of heads and potentials, looking up aquifer properties, etc.
  !! In addition, the AQU module must allow the other m_xxx element modules to access
  !! its "global" conversion functions.
  !!
  !! This module also incorporates the functionality previously in m_in0,
  !! providing inhomogeneity domain support directly.

  use u_constants
  use u_io
  use u_matrix
  use u_polygon
  use u_domain

  use f_reference
  use f_dipole
  use f_linesink
  use f_aem

  use i_linesink

  implicit none

  public

  type, public :: AQU_VERTEX
    !! type AQU_VERTEX
    !!
    !! Type that holds information for one vertex along the edge of an inhomogeneity
    !!
    !! Members:
    !!   complex :: cZ
    !!     The complex coordinate of the vertex
    !!   real :: rStrength(3)
    !!     The doublet strength at the vertex
    !!   type(FDP_DIPOLE), pointer :: pFDP
    !!     Pointer to the vertex entry in the FDP module.
    !!
    complex(kind=AE_REAL) :: cZ
    type(FDP_DIPOLE), pointer :: pFDP
    real(kind=AE_REAL), dimension(3) :: rStrength
    real(kind=AE_REAL), dimension(3) :: rCheckPot
    real(kind=AE_REAL), dimension(3) :: rLeftH
    real(kind=AE_REAL), dimension(3) :: rRightH
    real(kind=AE_REAL), dimension(3) :: rLeftT
    real(kind=AE_REAL), dimension(3) :: rRightT
    real(kind=AE_REAL), dimension(3) :: rError
    complex(kind=AE_REAL), dimension(3) :: cCPZ
    complex(kind=AE_REAL), dimension(3) :: cCPZOpp
  end type AQU_VERTEX

  type, public :: AQU_STRING
    !! type AQU_STRING
    !!
    !! Type that holds information for one inhomogeneity string.
    !! Strings are defined separately from the domains.
    !!
    !! Members:
    !!   type(AQU_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of AQU_VERTEX objects
    !!   integer :: iLeftID
    !!     ID of the AQU_DOMAIN that is to the left of the string
    !!   integer :: iRightID
    !!     ID of the AQU_DOMAIN that is to the right of the string
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   logical :: lClosed
    !!     .true. if the string closes on itself
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(AQU_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iLeftID
    integer(kind=AE_INT) :: iRightID
    real(kind=AE_REAL) :: rLeftB
    real(kind=AE_REAL) :: rRightB
    integer(kind=AE_INT) :: iNPts
    logical :: lClosed
    integer(kind=AE_INT) :: iID
  end type AQU_STRING

  type, public :: AQU_BDYELEMENT
    !! type AQU_BDYELEMENT
    !!
    !! Type that holds a vertex for the aquifer perimeter
    !!
    !! Members:
    !!   complex :: cZ1
    !!     The complex coordinate of the first vertex(x, y)
    !!   complex :: cZ2
    !!     The complex coordinate of the second vertex(x, y)
    !!   real :: rSpecHead
    !!     The specified head for the segment
    !!   real :: rSpecFlux
    !!     The specified flux for the segment
    !!   logical :: iBdyFlag
    !!     Flag : if .false., base computations on the specified flux
    !!            if .true., base computations on the specified head
    !!   integer :: iFDPIndex
    !!     Pointer in FDP to the dipole
    !!   real :: rLength
    !!     Length of the line segment
    !!   real :: rStrength
    !!     Sink density of the line segment
    !!
    complex(kind=AE_REAL) :: cZ1
    complex(kind=AE_REAL) :: cZ2
    complex(kind=AE_REAL) :: cCPZ1
    complex(kind=AE_REAL) :: cCPZ2
    complex(kind=AE_REAL) :: cZC
    complex(kind=AE_REAL) :: cZL
    real(kind=AE_REAL) :: rSpecHead
    real(kind=AE_REAL) :: rSpecFlux
    real(kind=AE_REAL) :: rGhbDistance
    real(kind=AE_REAL) :: rLength
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iBdyFlag
    real(kind=AE_REAL) :: rCheckPot
    real(kind=AE_REAL) :: rCheckFlux
    logical :: lMoveCPZ1                  ! Flag -- should CPZ1 be moved?
    logical :: lMoveCPZ2                  ! Flag -- should CPZ2 be moved?
    logical :: lFSHeadSpec                ! Flag -- is the free surface segment head-specified?
    complex(kind=AE_REAL) :: cFSZ1        ! Location of a "control point" at the first end
    complex(kind=AE_REAL) :: cFSZ2        ! Location of a "control point" at the second end
    integer(kind=AE_INT) :: iEqIndex
    real(kind=AE_REAL) :: rSaveSpecHead   ! Saves the user-specified head
    type(FLS_LINESINK), pointer :: pFLS   ! Pointer into aem%fls for this element
  end type AQU_BDYELEMENT


  type, public :: AQU_COLLECTION
    !! type AQU_COLLECTION
    !!
    !! Type that holds information for a layer
    !!
    !! Members:
    !!   type(FRF_COLLECTION), pointer :: frf
    !!     The singleton reference flow field (constant of integration + uniform flow)
    !!   type(AQU_BDYELEMENT) :: Boundary(:)
    !!     Vector of vertices which make up the perimeter of the aquifer
    !!   integer :: iNBdy
    !!     Number of perimeter points currently in use
    !!   type(AQU_DOMAIN) :: Domains(:)
    !!     Vector of inhomogeneity domains
    !!   integer :: iNDom
    !!     Number of domains in use
    !!   type(AQU_STRING) :: Strings(:)
    !!     Vector of inhomogeneity strings
    !!   integer :: iNStr
    !!     Number of strings in use
    !!
    type(FRF_COLLECTION), pointer :: frf
    type(AQU_BDYELEMENT), dimension(:), pointer :: BdyElements
    integer(kind=AE_INT) :: iNBdy
    ! Reference to the domain collection (owned by AEM_DOMAIN, assigned at creation)
    type(DOM_COLLECTION), pointer :: dom
    ! Inhomogeneity strings
    type(AQU_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    ! Perimeter of the aquifer
    complex(kind=AE_REAL), dimension(:), pointer :: cPerimeter
    integer(kind=AE_INT) :: iNPerimeter
    ! Regen on free surface
    logical :: lFSRegen
    ! Boundary condition locations have been previously set (and possibly moved)
    logical :: lBoundaryInPlace
    ! Turns on the new boundary control-point option
    logical :: lNewBdy
    ! Debug flag
    logical :: lDebug
  end type AQU_COLLECTION

  ! Locally-required constants
  ! Matrix generator flag for reference point equation
  integer(kind=AE_INT), private, parameter :: kAQUReference = 1
  integer(kind=AE_INT), private, parameter :: kAQUBoundary = 2
  ! Fake a small well radius to ensure that the proper fluxes are computed
  ! in FWL_Flow
  real(kind=AE_REAL), private, parameter :: PERIMETER_TOLERANCE = 1.0e-2_AE_REAL
  real(kind=AE_REAL), private, parameter :: BDY_MOVE_FACTOR = 1.0e-4_AE_REAL
  ! Matrix generator element flags (for inhomogeneities)
  integer(kind=AE_INT), private, parameter :: kAQU_Vertex = 1
  integer(kind=AE_INT), private, parameter :: kAQU_Center = 2
  integer(kind=AE_INT), private, parameter :: kAQU_Vertex2 = 3
  ! How far do I move the control points off the line segments?
  real(kind=AE_REAL), private, parameter :: AQU_NORMAL_OFFSET = -1.0e-7_AE_REAL
  real(kind=AE_REAL), private, parameter :: AQU_OPPOSITE_OFFSET = 1.0e-7_AE_REAL

  real(kind=AE_REAL) :: MOVEFACTOR = 1.0001_AE_REAL



contains

  !! ELEMENT MODULE ROUTINES
  !! These routines allow AQU to behave as a ModAEM Element Module


  function AQU_Create(io, dom_coll, frf_coll, iNStr) result(aqu)
    !! function AQU_Create
    !!
    !! Creates a new AQU_COLLECTION object and attaches it to existing domain and reference collections.
    !!
    !! Calling Sequence:
    !!    aqu => AQU_Create(io, dom_coll, frf_coll, iNStr)
    !!
    !! Arguments:
    !!   type(IO_STATUS), pointer :: io
    !!   type(DOM_COLLECTION), pointer :: dom_coll
    !!     The domain collection (owned by AEM_DOMAIN) to attach to this aquifer
    !!   type(FRF_COLLECTION), pointer :: frf_coll
    !!     The reference flow field collection (owned by AEM_DOMAIN) to attach to this aquifer
    !!   integer :: iNStr
    !!     Number of inhomogeneity strings to pre-allocate
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    type(FRF_COLLECTION), pointer :: frf_coll
    integer(kind=AE_INT), intent(in) :: iNStr

    ! [ RETURN VALUE ]
    type(AQU_COLLECTION), pointer :: aqu
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(aqu, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AQU_Create: allocation failed")

    ! Attach the domain and reference flow field collections; allocate string space
    aqu%dom => dom_coll
    aqu%frf => frf_coll
    allocate(aqu%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "AQU_Create: allocation failed")
    aqu%Strings%iID = -1

    aqu%iNStr = 0

    nullify(aqu%BdyElements)
    aqu%iNBdy = 0
    nullify(aqu%cPerimeter)
    aqu%iNPerimeter = 0
    ! First time up, no regen
    aqu%lFSRegen = .false.
    ! Boundary isn't set yet
    aqu%lBoundaryInPlace = .false.
    ! Default to the old control-point strategy
    aqu%lNewBdy = .false.

    return
  end function AQU_Create


  subroutine AQU_SetPrecondition(io, aqu, lPre)
    ! Sets the precondition flag. If it is set, it is presumed that this is the first iteration, and the
    ! saturated thickness / transmissivity are based on the "avg-head" setting.
    type(AQU_COLLECTION), pointer :: aqu
    logical, intent(in) :: lPre
    type(IO_STATUS), pointer :: io

    aqu%dom%lPrecondition = lPre

    return
  end subroutine AQU_SetPrecondition


  subroutine AQU_SetReference(io, aqu, cRefPoint, rRefHead, cRefUniformFlow)
    !! subroutine AQU_SetReference
    !!
    !! Sets up the reference point and uniform flow
    !!
    !! Calling Sequence:
    !!    call AQU_SetReference(io, aqu, cRefPoint, rRefHead, cRefUniformFlow)
    !!
    !! Arguments:
    !!    (in)    type(AQU_COLLECTION), pointer :: aqu
    !!              The AQU_COLLECTION object to be populated
    !!    (in)    complex :: cRefPoint
    !!              Location of the reference point
    !!    (in)    real :: rRefHead
    !!              The head at the reference point
    !!    (in)    complex :: cRefUniformFlow
    !!              Far-field uniform flow
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cRefPoint
    real(kind=AE_REAL), intent(in) :: rRefHead
    complex(kind=AE_REAL), intent(in) :: cRefUniformFlow
    type(IO_STATUS), pointer :: io

    call IO_Assert(io, associated(aqu), "AQU_SetReference: No AQU_COLLECTION object")

    call FRF_SetReference(io, aqu%frf, cRefPoint, rRefHead, cRefUniformFlow)

    return
  end subroutine AQU_SetReference


  subroutine AQU_PreSolve(io, aqu)
    !! subroutine AQU_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call AQU_PreSolve(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, istat, iDom, iStr, iVtx
    character(len=255) :: sBuf
    type(DOM_DOMAIN), pointer :: dom
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx
    complex(kind=AE_REAL), dimension(:), allocatable :: cTempZ
    real(kind=AE_REAL) :: rH, rT

    ! Fix up the perimeter stuff as necessary
    if (aqu%iNBdy > 0) then
      allocate(aqu%cPerimeter(aqu%iNBdy), stat = istat)
      call IO_Assert(io, (istat == 0), "Could not allocate aquifer perimeter")
      call PGN_Build(aqu%BdyElements(1:aqu%iNBdy)%cZ1, &
           aqu%BdyElements(1:aqu%iNBdy)%cZ2, &
           aqu%BdyElements(1:aqu%iNBdy)%iBdyFlag==BDY_HEAD, &
           PERIMETER_TOLERANCE, &
           aqu%cPerimeter, &
           aqu%iNPerimeter)
      call IO_Assert(io, (aqu%iNPerimeter > 0), "Could not build the aquifer perimeter")
    else
      allocate(aqu%cPerimeter(0), stat = istat)
      aqu%iNPerimeter = 0
    end if

    ! Orient domains and resolve nesting topology via u_domain
    call DOM_PreSolve(io, aqu%dom)

    ! If no strings were explicitly defined, create strings for the domains
    if (aqu%iNStr == 0) then
      print *, 'Copying AQU domains into strings...'
      call IO_Assert(io, (size(aqu%Strings) >= aqu%dom%iNDom-1), &
           "AQU_PreSolve: Insufficient number of strings is available")
      do iDom = 2, aqu%dom%iNDom
        aqu%iNStr = aqu%iNStr+1
        dom => aqu%dom%Domains(iDom)
        str => aqu%Strings(aqu%iNStr)
        allocate(str%Vertices(dom%iNPts), stat = iStat)
        call IO_Assert(io, (iStat == 0), &
             "AQU_PreSolve: Space exhausted")
        str%iNPts = dom%iNPts
        do iVtx = 1, dom%iNPts
          vtx => str%Vertices(iVtx)
          vtx%cZ = dom%cZ(iVtx)
          vtx%rStrength = rZERO
          vtx%rCheckPot = rZERO
          nullify(vtx%pFDP)
        end do
        str%iLeftID = dom%iInsideDomain
        str%iRightID = dom%iOutsideDomain
        str%lClosed = .true.
        str%iID = dom%iID
      end do
    end if

    ! Set the RightB and LeftB base elevations for all strings
    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      dom => DOM_FindDomainID(io, aqu%dom, str%iLeftID)
      str%rLeftB = dom%rBase
      if (dom%rAvgHead > (dom%rBase+dom%rThickness)) then
        rH = dom%rThickness
      else
        rH = dom%rAvgHead - dom%rBase
      end if
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        vtx%rLeftH = rH
        vtx%rLeftT = rH * dom%rHydCond
      end do
      dom => DOM_FindDomainID(io, aqu%dom, str%iRightID)
      str%rRightB = dom%rBase
      if (dom%rAvgHead > (dom%rBase+dom%rThickness)) then
        rH = dom%rThickness
      else
        rH = dom%rAvgHead - dom%rBase
      end if
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        vtx%rRightH = rH
        vtx%rRightT = rH * dom%rHydCond
      end do
    end do

    return
  end subroutine AQU_PreSolve


  function iAQU_GetInfo(io, aqu, iOption, iIteration) result(iValue)
    !! function AQU_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iAQU_GetInfo(io, aqu, iOption)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr
    type(AQU_STRING), pointer :: str

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_GetInfo: AQU_Create has not been called")
    end if

    ! First, count the inhomogeneity contributions (was iIN0_GetInfo)
    iValue = 0
    select case (iOption)
      case (SIZE_FLS)
        iValue = aqu%iNBdy
      case (SIZE_FDP)
        do iStr = 1, aqu%iNStr
          str => aqu%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + str%iNPts
          else
            iValue = iValue + str%iNPts-1
          end if
        end do
      case (SIZE_EQUATIONS)
        do iStr = 1, aqu%iNStr
          str => aqu%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + 2*str%iNPts
          else
            iValue = iValue + 2*str%iNPts-1
          end if
        end do
        iValue = iValue + aqu%iNBdy
        if (aqu%frf%lReference .or. aqu%iNBdy == 0) iValue = iValue + 1
      case (SIZE_UNKNOWNS)
        do iStr = 1, aqu%iNStr
          str => aqu%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + 2*str%iNPts
          else
            iValue = iValue + 2*str%iNPts-1
          end if
        end do
        iValue = iValue + aqu%iNBdy
        if (aqu%frf%lReference .or. aqu%iNBdy == 0) iValue = iValue + 1
      case (INFO_REGENERATE)
        ! Force a regen on the first two iterations
        if (iIteration < 2) then
          iValue = 1
        else
          iValue = aqu%dom%iRegenerate
        end if
        if (iIteration < 2 .or. aqu%lFSRegen) then
          iValue = 1
        end if
      case default
        iValue = 0
    end select

    return
  end function iAQU_GetInfo


  subroutine AQU_SetupFunctions(io, aqu, fdp, fls)
    !! subroutine AQU_SetupFunctions
    !!
    !! This routine sets up the functions in f_well and f_dipole for the perimeter
    !! and sets up the equation for the reference point.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    type(FLS_COLLECTION), pointer :: fls
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iStr, iVtx, iBdy
    real(kind=AE_REAL) :: rSum
    complex(kind=AE_REAL) :: cZ1, cZ2
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "AQU_Setup: AQU_Create has not been called")

    ! Set the internal values for the BDY elements and register each one in fls
    if (aqu%iNBdy > 0) then
      aqu%BdyElements(1:aqu%iNBdy)%cZC = rHALF * (aqu%BdyElements(1:aqu%iNBdy)%cZ2 + aqu%BdyElements(1:aqu%iNBdy)%cZ1)
      aqu%BdyElements(1:aqu%iNBdy)%cZL = rHALF * (aqu%BdyElements(1:aqu%iNBdy)%cZ2 - aqu%BdyElements(1:aqu%iNBdy)%cZ1)
      do iBdy = 1, aqu%iNBdy
        bdy => aqu%BdyElements(iBdy)
        bdy%pFLS => FLS_New(io, fls, bdy%cZ1, bdy%cZ2, cZERO, ELEM_AQU, kAQUBoundary, iBdy, 0)
      end do
    end if

    ! Set the confined potential (was IN0_SetupFunctions)
    ! Build dipoles for all segments
    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), &
           "AQU_SetupFunctions: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        cZ1 = vtx%cZ
        if (iVtx < str%iNPts) then
          cZ2 = str%Vertices(iVtx+1)%cZ
          vtx%pFDP => FDP_New(io, fdp, cZ1, cZ2, (/cZERO, cZERO, cZERO/), ELEM_IN0, iStr, iVtx, -1)
        else
          if (str%lClosed) then
            cZ2 = str%Vertices(1)%cZ
            vtx%pFDP => FDP_New(io, fdp, cZ1, cZ2, (/cZERO, cZERO, cZERO/), ELEM_IN0, iStr, iVtx, -1)
          end if
        end if
      end do
    end do

    return
  end subroutine AQU_SetupFunctions


  subroutine AQU_SetupMatrix(io, aqu, mat)
    !! subroutine AQU_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the module
    !! and sets up the equation for the reference point.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cCPZ0, cCPZ1, cCPZ2, cNormal, cN1, cN2, cUnitOffset
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3
    complex(kind=AE_REAL), dimension(6) :: cCPResult1, cCPResult2
    complex(kind=AE_REAL), dimension(3) :: cCPVtx
    real(kind=AE_REAL) :: rElev, rOffsetDistance
    integer(kind=AE_INT) :: iBdy, iStr, iVtx
    integer(kind=AE_INT) :: iNWL, iNPD, iNDP, iNEQ, iNUN, iEQ
    integer(kind=AE_INT) :: irv
    type(AQU_BDYELEMENT), pointer :: this, prev, next
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx
    character(len=255) :: sBuf

    call IO_Assert(io, (associated(aqu)), "AQU_Setup: AQU_Create has not been called")

    ! Build matrix entries for inhomogeneity strings (was IN0_SetupMatrix)
    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "AQU_SetupMatrix: Illegal MAT_MATRIX object")
    end if

    ! Build matrix entries for all segments
    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      ! Set up the unknown variables
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        ! Make a doublet vertex variable
        call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kAQU_Vertex)
        call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kAQU_Center)
        if (.not. str%lClosed .and. iVtx == str%iNPts-1) then
          ! Create the 'tip' variable, only if not closed, then exit the loop!
          call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kAQU_Vertex2)
          exit
        end if
      end do

      ! Set up control point sets and equations
      if (str%lClosed) then
        ! Closed boundary
        do iVtx = 1, str%iNPts
          vtx => str%Vertices(iVtx)
          ! Now, create the equation entries...
          cZ1 = vtx%cZ
          if (iVtx < str%iNPts) then
            cZ2 = str%Vertices(iVtx+1)%cZ
          else
            cZ2 = str%Vertices(1)%cZ
          end if
          ! Compute control point locations
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, AQU_NORMAL_OFFSET, &
               (/0.001_AE_REAL, 0.501_AE_REAL/))
          ! Vertex entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(1)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kAQU_Vertex, cZ2-cZ1, rZERO)
          ! Center entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(2)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kAQU_Center, cZ2-cZ1, rZERO)
          vtx%cCPZ(1:2) = cCPResult1(1:2)
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, AQU_NORMAL_OFFSET, &
               (/0.001_AE_REAL, 0.501_AE_REAL/))
          vtx%cCPZOpp(1:2) = cCPResult1(1:2)
        end do
      else
        ! Open string
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          ! Now, create the equation entries...
          cZ1 = vtx%cZ
          cZ2 = str%Vertices(iVtx+1)%cZ
          ! Compute control point locations
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, AQU_NORMAL_OFFSET, &
               (/0.001_AE_REAL, rHALF, rONE-0.001_AE_REAL/))
          ! Vertex entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(1)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kAQU_Vertex, cZ2-cZ1, rZERO)
          vtx%cCPZ(kAQU_Vertex) = cCPResult1(1)
          ! Center entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(2)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kAQU_Center, cZ2-cZ1, rZERO)
          vtx%cCPZ(kAQU_Center) = cCPResult1(2)
          if (iVtx == str%iNPts-1) then
            ! End of an open string
            iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(3)/), EQN_INHO, ELEM_IN0, &
                                     iStr, iVtx, kAQU_Vertex2, cZ2-cZ1, rZERO)
            vtx%cCPZ(kAQU_Vertex2) = cCPResult1(3)
          end if
        end do
      end if

    end do

    ! Build the matrix entries for the perimeter(if any)
    if (aqu%iNBdy > 0) then
      prev => aqu%BdyElements(aqu%iNBdy)
      rOffsetDistance = BDY_MOVE_FACTOR * sum(aqu%BdyElements(1:aqu%iNBdy)%rLength) / aqu%iNBdy
      do iBdy = 1, aqu%iNBdy
        this => aqu%BdyElements(iBdy)
        call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUBoundary, iBdy, 0)

        ! The locations for flux inquiries are moved "inside" the element
        cNormal = cI * (this%cZ2 - this%cZ1)
        if (.not. aqu%lBoundaryInPlace) then
          if (.not. aqu%lNewBdy) then
            this%cCPZ1 = this%cZ1 + (this%cZ2-this%cZ1)*BDY_MOVE_FACTOR + &
                             cNormal*BDY_MOVE_FACTOR*abs(this%cZ2-this%cZ1)
            this%cCPZ2 = this%cZ2 - (this%cZ2-this%cZ1)*BDY_MOVE_FACTOR +  &
                             cNormal*BDY_MOVE_FACTOR*abs(this%cZ2-this%cZ1)
          else
            ! Find the unit-normals to this element and the next...
            call AQU_GetNeighborBdy(io, aqu, this, prev, next)
            cN1 = cI * (this%cZ2 - this%cZ1) / abs(this%cZ2 - this%cZ1)
            cN2 = cI * (next%cZ2 - next%cZ1) / abs(next%cZ2 - next%cZ1)
            ! The unit vector pointing from this element's end 2 to the control point
            ! bisects the direction of cN1 and cN2
            cUnitOffset = (cN1 + cN2) / abs(cN1 + cN2)
            this%cCPZ2 = this%cZ2 + rOffsetDistance * cUnitOffset
            next%cCPZ1 = next%cZ1 + rOffsetDistance * cUnitOffset
          end if
          ! Save the starting locations for the vertices
          this%cFSZ1 = this%cCPZ1
          this%cFSZ2 = this%cCPZ2
        end if

        select case (this%iBdyFlag)
          case (BDY_HEAD)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/rHALF*(this%cCPZ1+this%cCPZ2)/), &
                                                   EQN_HEAD, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
          case (BDY_FLUX)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                   EQN_FLOW, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
          case (BDY_GHB)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                   EQN_BDYGHB, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, &
                 this%rGhbDistance/this%rLength)
          case (BDY_FREESURF)
            if (.not. aqu%lBoundaryInPlace) then
              rElev = max(aimag(this%cFSZ1), aimag(this%cFSZ2))
              if (this%rSpecHead >= rElev) then
                this%lFSHeadSpec = .true.
              else
                this%lFSHeadSpec = .false.
              end if
            end if

            if (this%lFSHeadSpec) then
              ! It's a head
              this%iEqIndex = MAT_CreateEquation(io, mat, (/rHALF*(this%cCPZ1+this%cCPZ2)/), &
                                                     EQN_HEAD, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
            else
              ! It's a moveable boundary
              this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                     EQN_FLOW, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
            end if
        end select
        prev => aqu%BdyElements(iBdy)
      end do
      aqu%lBoundaryInPlace = .true.

    end if
    ! Build a matrix entry for the reference point(if any)
    if (aqu%frf%lReference) then
      call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUReference, 0, 0)
      iEQ = MAT_CreateEquation(io, mat, (/aqu%frf%cRefPoint/), EQN_HEAD, ELEM_AQU, kAQUReference, 0, 0, cZERO, rZERO)
    else if (aqu%iNBdy == 0) then
      call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUReference, 0, 0)
      iEQ = MAT_CreateEquation(io, mat, (/cZERO/), EQN_CONTINUITY, ELEM_AQU, kAQUReference, 0, 0, cZERO, rZERO)
    end if

    return
  end subroutine AQU_SetupMatrix


  function iAQU_Prepare(io, aqu, iIteration) result(iChanges)
    !! subroutine AQU_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_aqu
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iAQU_Prepare


  subroutine AQU_GetNeighborBdy(io, aqu, this, prev, next)
    !! For profile models, returns 'prev' and 'next' as the previous and
    !! next (corresponding to cZ1 and cZ2, respectively) BDY elements adjoining
    !! the element 'this'.
    !! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    type(AQU_BDYELEMENT), pointer :: this
    type(AQU_BDYELEMENT), pointer :: prev
    type(AQU_BDYELEMENT), pointer :: next
    !! [ LOCALS ]
    integer(kind=AE_INT) :: ibdy
    type(AQU_BDYELEMENT), pointer :: bdy

    ! Find the neighbors
    prev => NULL()
    next => NULL()
    do ibdy=1,aqu%iNBdy
      bdy => aqu%BdyElements(ibdy)
      if (abs(this%cZ1 - bdy%cZ2) < PERIMETER_TOLERANCE) prev => aqu%BdyElements(ibdy)
      if (abs(this%cZ2 - bdy%cZ1) < PERIMETER_TOLERANCE) next => aqu%BdyElements(ibdy)
    end do

    call IO_Assert(io, (associated(prev) .and. associated(next)), "Failed to find neighboring BDY elements")

    return
  end subroutine AQU_GetNeighborBdy


  function rAQU_GetCoefficientMultiplier(io, aqu, iElementType, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier
    ! [ LOCALS ]
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx

    select case (iElementType)
      case (ELEM_IN0)
        str => aqu%Strings(iElementString)
        vtx => str%Vertices(iElementVertex)
        rMultiplier = vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)
      case default
        rMultiplier = rONE
    end select

    return
  end function rAQU_GetCoefficientMultiplier


  subroutine AQU_ComputeCoefficients(io, aqu, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor, &
               iIteration, rMultiplier, rARow)
    !! subroutine AQU_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the AQU
    !! elements.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    complex(kind=AE_REAL), intent(in) :: cOrientation
    real(kind=AE_REAL), intent(in) :: rGhbFactor
    integer(kind=AE_INT), intent(in) :: iIteration
    real(kind=AE_REAL), intent(in) :: rMultiplier
    real(kind=AE_REAL), dimension(:), intent(out) :: rARow
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStat, iCol, iBdy, iAQUNUN
    integer(kind=AE_INT) :: iSkipCol, iVtx, iNDP, iWhich, irv, iBaseCol, iStr
    complex(kind=AE_REAL) :: cZC, cZL, cBigZ1, cBigZ2, cUnit
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDPW
    complex(kind=AE_REAL), dimension(1, 3, 1) :: cDPJ
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: this_vtx, next_vtx

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_ComputeCoefficients: No AQU_COLLECTION object")
    end if
    ! Set it up
    rARow = rZERO

    ! Build entries for the inhomogeneities (was IN0_ComputeCoefficients)
    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), &
           "AQU_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iSkipCol = 0
    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      ! Assume: the AQU_SetupFunctions routine creates consecutive dipole entries
      if (str%lClosed) then
        iNDP = str%iNPts
      else
        iNDP = str%iNPts-1
      end if
      allocate(cDPF(0:iNDP+1, 3, 1), cDPW(0:iNDP+1, 3, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "AQU_ComputeCoefficients: Allocation failed")

      ! Now, compute the matrix coefficients
      cDPF = cZERO
      ! Get the appropriate influence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_BDYGHB)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, str%Vertices(1)%pFDP, iNDP, (/rHALF*sum(cPathZ)/), cOrientation, cDPF(1:iNDP, :, :))
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPW(1:iNDP, :, :))
          cDPF = cDPF + rGhbFactor*cDPW
        case (EQN_FLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_INHO)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_DISCHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_W, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_RECHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_G, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_CONTINUITY)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Q, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_POTENTIALDIFF)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_D, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_TOTALFLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Z, str%Vertices(1)%pFDP, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
      end select

      do iVtx = 1, iNDP
        this_vtx => str%Vertices(iVtx)
        if (iVtx < iNDP) then
          next_vtx => str%Vertices(iVtx+1)
        else
          next_vtx => str%Vertices(1)
        end if

        iBaseCol = 2*iVtx-1 + iSkipCol

        ! Compute the contributions for all line-doublets
        ! Vertex 1 contribution
        rARow(iBaseCol) = rARow(iBaseCol)   - rMultiplier * aimag(cDPF(iVtx, 1, 1))
        ! Center contribution
        rARow(iBaseCol+1) = rARow(iBaseCol+1) - rMultiplier * aimag(cDPF(iVtx, 2, 1))
        ! Vertex 2 contribution
        if (str%lClosed) then
          ! For closed strings...
          if (iVtx < iNDP) then
            rARow(iBaseCol+2) = rARow(iBaseCol+2) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
          else
            rARow(1+iSkipCol) = rARow(1+iSkipCol) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
          end if
        else
          ! For open strings, store the last term
          rARow(iBaseCol+2) = rARow(iBaseCol+2) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
        end if

        ! Do I compute the additional term for the inhomogeneity?
        if ((iEqType == EQN_INHO) .and. &
            (iElementType == ELEM_IN0) .and. &
            (iElementString == iStr) .and. &
            (iElementVertex == iVtx)) then
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_J, str%Vertices(iVtx)%pFDP, 1, cPathZ(1:1), cOrientation, cDPJ)
          ! Vertex 1 contribution
          rARow(iBaseCol) = rARow(iBaseCol)   - this_vtx%rRightT(1)*aimag(cDPJ(1, 1, 1))
          ! Center contribution
          rARow(iBaseCol+1) = rARow(iBaseCol+1) - this_vtx%rRightT(2)*aimag(cDPJ(1, 2, 1))
          ! Vertex 2 contribution
          if (str%lClosed) then
            ! For closed strings...
            if (iVtx < iNDP) then
              rARow(iBaseCol+2) = rARow(iBaseCol+2) - next_vtx%rRightT(1)*aimag(cDPJ(1, 3, 1))
            else
              rARow(1+iSkipCol) = rARow(1+iSkipCol) - next_vtx%rRightT(1)*aimag(cDPJ(1, 3, 1))
            end if
          else
            ! For open strings, update the last term
            rARow(iBaseCol+2) = rARow(iBaseCol+2) - this_vtx%rRightT(3)*aimag(cDPJ(1, 3, 1))
          end if
        end if
      end do

      ! Now, advance through the buffer
      if (str%lClosed) then
        iSkipCol = iSkipCol + 2*str%iNPts
      else
        iSkipCol = iSkipCol + 2*str%iNPts - 1
      end if

      ! No memory leaks, please!
      deallocate(cDPF, cDPW)
    end do

    iCol = iSkipCol + 1

    ! Build entries for the boundary elements
    if (associated(aqu%BdyElements)) then
      ! Get the appropriate influence functions for the boundary condition type
      allocate(cBigZ(size(cPathZ)), stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AQU_ComputeCoefficients: Allocation failed for big-Z')
      do iBdy = 1, aqu%iNBdy
        bdy => aqu%BdyElements(iBdy)
        cZC = rHALF*(bdy%cZ2+bdy%cZ1)
        cZL = rHALF*(bdy%cZ2-bdy%cZ1)
        cBigZ = (cPathZ - cZC) / cZL
        select case (iEqType)
          case (EQN_HEAD)
            cTemp = real(cILS_InfluenceP(cBigZ(1), cZL))
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_BDYGHB)
            cTemp = real(cILS_InfluenceP(rHALF*(cBigZ(1)+cBigZ(2)), cZL))
            cTemp = cTemp + rGhbFactor * real(cILS_InfluenceF(cBigZ(1), cBigZ(2), cZL))
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_FLOW)
            do i = 1, size(cBigZ)-1
              cTemp = cILS_InfluenceF(cBigZ(i), cBigZ(i+1), cZL)
              rARow(iCol) = rMultiplier * real(cTemp(1, 1))
            end do
          case (EQN_INHO)
            cTemp = cILS_InfluenceP(cBigZ(1), cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_DISCHARGE)
            cUnit = cOrientation / abs(cOrientation)
            cTemp = cILS_InfluenceW(cBigZ(1), cZL)
            rARow(iCol) = rMultiplier * real(cUnit * cTemp(1, 1))
          case (EQN_POTENTIALDIFF)
            cTemp = cILS_InfluenceP(cBigZ(1), cZL)-cILS_InfluenceP(cBigZ(2), cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_CONTINUITY)
            cBigZ1 = (bdy%cCPZ1 - cZC) / cZL
            cBigZ2 = (bdy%cCPZ2 - cZC) / cZL
            cTemp = cILS_InfluenceF(cBigZ1, cBigZ2, cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_TOTALFLOW)
            rARow(iCol) = rZERO
        end select

        iCol = iCol+1
      end do
      deallocate(cBigZ)
    end if

    ! Set up the closure condition
    if (aqu%frf%lReference .or. aqu%iNBdy == 0) then
      select case (iEqType)
        case (EQN_HEAD)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_BDYGHB)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_FLOW)
          rARow(iCol) = rZERO
        case (EQN_INHO)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_DISCHARGE)
          rARow(iCol) = rZERO
        case (EQN_RECHARGE)
          rARow(iCol) = rZERO
        case (EQN_CONTINUITY)
          rARow(iCol) = rZERO
      end select
      iCol = iCol+1
    end if

    return
  end subroutine AQU_ComputeCoefficients


  function rAQU_ComputeRHS(io, aqu, fdp, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rAQU_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
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
    real(kind=AE_REAL) :: rCorr
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx, next_vtx

    select case (iElementType)
      case (ELEM_AQU)
        select case (iElementString)
          case (kAQUBoundary)
            bdy => aqu%BdyElements(iElementVertex)
            select case (bdy%iBdyFlag)
              case (BDY_HEAD)
                if (lDirect) then
                  rRHS = rDOM_HeadToPotential(io, aqu%dom, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2))
                else
                  rRHS = rDOM_HeadToPotential(io, aqu%dom, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - bdy%rCheckPot
                end if
              case (BDY_FLUX)
                if (lDirect) then
                  rRHS = bdy%rSpecFlux
                else
                  rRHS = bdy%rSpecFlux - bdy%rCheckFlux
                end if
              case (BDY_GHB)
                rRHS = rDOM_HeadToPotential(io, aqu%dom, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - &
                       bdy%rCheckPot - bdy%rCheckFlux*bdy%rGhbDistance/bdy%rLength
              case (BDY_FREESURF)
                if (bdy%lFSHeadSpec) then
                  if (lDirect) then
                    rRHS = rDOM_HeadToPotential(io, aqu%dom, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2))
                  else
                    rRHS = rDOM_HeadToPotential(io, aqu%dom, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - bdy%rCheckPot
                  end if
                else
                  if (lDirect) then
                    rRHS = bdy%rSpecFlux
                  else
                    rRHS = bdy%rSpecFlux - bdy%rCheckFlux
                  end if
                end if
            end select
          case (kAQUReference)
            if (aqu%frf%lReference) then
              if (lDirect) then
                rRHS = rDOM_HeadToPotential(io, aqu%dom, aqu%frf%rRefHead, aqu%frf%cRefPoint)
              else
                rRHS = rDOM_HeadToPotential(io, aqu%dom, aqu%frf%rRefHead, aqu%frf%cRefPoint) - aqu%frf%rCheck
              end if
            else
              rRHS = rZERO
            end if
        end select
      case (ELEM_IN0)
        ! Compute the inhomogeneity RHS (was rIN0_ComputeRHS)
        ! Assertions...
        if (io%lDebug) then
          call IO_Assert(io, (iElementString <= aqu%iNStr), &
               "rAQU_ComputeRHS: Illegal string number")
        end if
        ! Compute the correction needed in the jump...
        str => aqu%Strings(iElementString)
        if (io%lDebug) then
          call IO_Assert(io, (iElementVertex <= str%iNPts), &
               "rAQU_ComputeRHS: Illegal vertex number")
          if (iElementVertex == str%iNPts .and. .not. str%lClosed) then
            call IO_Assert(io, (iElementFlag == kAQU_Vertex2), &
                 "rAQU_ComputeRHS: No center strength at last vertex of open string")
          end if
        end if
        vtx => str%Vertices(iElementVertex)

        ! The sign on the "Jump" is negative due to the i*i factor in the computations
        if (lDirect) then
          rRHS = - (vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)) * vtx%rCheckPot(iElementFlag) + &
                 vtx%rLeftT(iElementFlag) * vtx%rRightT(iElementFlag) * &
                 (str%rLeftB-str%rRightB + rHALF*(vtx%rLeftH(iElementFlag)-vtx%rRightH(iElementFlag)))
        else
          rRHS = - (vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)) * vtx%rCheckPot(iElementFlag) - &
                 vtx%rRightT(iElementFlag) * rFDP_PotentialJump(io, vtx%pFDP, vtx%cCPZ(iElementFlag)) + &
                 vtx%rLeftT(iElementFlag) * vtx%rRightT(iElementFlag) * &
                 (str%rLeftB-str%rRightB + rHALF*(vtx%rLeftH(iElementFlag)-vtx%rRightH(iElementFlag)))
        end if
    end select

    return
  end function rAQU_ComputeRHS


  subroutine AQU_StoreResult(io, aqu, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine AQU_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the AQU module.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(AQU_BDYELEMENT), pointer :: this
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), "AQU_StoreResult: AQU_Create has not been called")
    end if

    select case (iElementType)
      case (ELEM_AQU)
        select case (iElementString)
          case (kAQUReference)
            if (lDirect) then
              aqu%frf%rSolConst = rValue
            else
              aqu%frf%rSolConst = aqu%frf%rSolConst + rValue
            end if
          case (kAQUBoundary)
            this => aqu%BdyElements(iElementVertex)
            if (lDirect) then
              this%rStrength = rValue
            else
              this%rStrength = this%rStrength + rValue
            end if
            if (associated(this%pFLS)) &
              this%pFLS%cSigma = cmplx(this%rStrength, rZERO, AE_REAL)
        end select
      case (ELEM_IN0)
        ! Store inhomogeneity result (was IN0_StoreResult)
        if (io%lDebug) then
          call IO_Assert(io, (iElementString >= 1 .and. iElementString <= aqu%iNStr), &
               "AQU_StoreResult: Bad element string ID")
        end if

        str => aqu%Strings(iElementString)

        if (io%lDebug) then
          call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= str%iNPts), &
               "AQU_StoreResult: Bad element vertex ID")
        end if

        ! All is well.  Store the result...
        vtx => str%Vertices(iElementVertex)
        if (lDirect) then
          vtx%rStrength = rZERO
        end if

        select case (iElementFlag)
          case (kAQU_Vertex)
            vtx%rStrength(1) = vtx%rStrength(1) + rValue
          case (kAQU_Center)
            vtx%rStrength(2) = vtx%rStrength(2) + rValue
          case (kAQU_Vertex2)
            vtx%rStrength(3) = vtx%rStrength(3) + rValue
        end select
    end select

    return
  end subroutine AQU_StoreResult


  subroutine AQU_Update(io, aqu, fdp)
    !! subroutine AQU_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iBdy, iVtx, iRV
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(AQU_BDYELEMENT), pointer :: this
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: this_vtx, next_vtx

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), "AQU_Update: AQU_Create has not been called")
    end if

    ! Update the inhomogeneities (was IN0_Update)
    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), &
           "AQU_Update: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      if (str%lClosed) then
        do iVtx = 1, str%iNPts
          this_vtx => str%Vertices(iVtx)
          cRho1 = cmplx(rZERO, this_vtx%rStrength(1), AE_REAL)
          cRho2 = cmplx(rZERO, this_vtx%rStrength(2), AE_REAL)
          if (iVtx < str%iNPts) then
            next_vtx => str%Vertices(iVtx+1)
          else
            next_vtx => str%Vertices(1)
          end if
          cRho3 = cmplx(rZERO, next_vtx%rStrength(1), AE_REAL)
          this_vtx%pFDP%cRho = (/cRho1, cRho2, cRho3/)
        end do
      else
        do iVtx = 1, str%iNPts-1
          this_vtx => str%Vertices(iVtx)
          cRho1 = cmplx(rZERO, this_vtx%rStrength(1), AE_REAL)
          cRho2 = cmplx(rZERO, this_vtx%rStrength(2), AE_REAL)
          if (iVtx < str%iNPts-1) then
            next_vtx => str%Vertices(iVtx+1)
            cRho3 = cmplx(rZERO, next_vtx%rStrength(1), AE_REAL)
          else
            cRho3 = cmplx(rZERO, this_vtx%rStrength(3), AE_REAL)
          end if
          this_vtx%pFDP%cRho = (/cRho1, cRho2, cRho3/)
        end do
      end if
    end do

    ! If we've done this once, it's not the initial iteration anymore
    aqu%dom%lInitialIteration = .false.

    return
  end subroutine AQU_Update


  subroutine AQU_ComputeCheck(io, aqu, aem, lLinearize)
    !! Updates check potentials, fluxes, and linearization for all AQU elements.
    type(AQU_COLLECTION), pointer :: aqu
    type(AEM_DOMAIN), pointer :: aem
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx
    type(DOM_DOMAIN), pointer :: left, right
    complex(kind=AE_REAL), dimension(2) :: cZFlow
    real(kind=AE_REAL) :: H1, H2, rHead
    integer(kind=AE_INT) :: iBdy, iStr, iVtx, iFlag, iMaxFlag, iNVtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(aqu)), "AQU_ComputeCheck: AQU_Create has not been called")
    end if

    if (aqu%frf%lReference) then
      aqu%frf%rCheck = real(cAEM_Potential(io, aem, aqu%frf%cRefPoint, .false.), AE_REAL)
    else
      aqu%frf%rCheck = rAEM_Extraction(io, aem, .false.)
    end if

    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      bdy%rCheckPot = real(cAEM_Potential(io, aem, rHALF*(bdy%cCPZ1+bdy%cCPZ2), .false.), AE_REAL)
      cZFlow = (/bdy%cCPZ1, bdy%cCPZ2/)
      bdy%rCheckFlux = rAEM_Flow(io, aem, cZFlow, .false.)
    end do

    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      if (str%lClosed) then
        iNVtx = str%iNPts
      else
        iNVtx = str%iNPts - 1
      end if
      do iVtx = 1, iNVtx
        vtx => str%Vertices(iVtx)
        iMaxFlag = kAQU_Center
        if (.not. str%lClosed .and. iVtx == str%iNPts - 1) iMaxFlag = kAQU_Vertex2
        do iFlag = kAQU_Vertex, iMaxFlag
          vtx%rCheckPot(iFlag) = real(cAEM_Potential(io, aem, vtx%cCPZ(iFlag), .false.), AE_REAL)
          aqu%dom%iRegenerate = 0
          if (lLinearize) then
            left => DOM_FindDomainID(io, aqu%dom, str%iLeftID)
            right => DOM_FindDomainID(io, aqu%dom, str%iRightID)
            rHead = rHALF * (rDOM_PotentialToHead(io, aqu%dom, vtx%rCheckPot(iFlag), vtx%cCPZ(iFlag)) + &
                    rDOM_PotentialToHead(io, aqu%dom, vtx%rCheckPot(iFlag), vtx%cCPZ(iFlag)))
            H1 = rDOM_DomainThickness(io, aqu%dom, left, rHead)
            H2 = rDOM_DomainThickness(io, aqu%dom, right, rHead)
            if (H1 /= vtx%rLeftH(iFlag) .or. H2 /= vtx%rRightH(iFlag)) then
              aqu%dom%iRegenerate = 1
              vtx%rLeftH(iFlag) = H1
              vtx%rRightH(iFlag) = H2
              vtx%rLeftT(iFlag) = H1 * left%rHydCond
              vtx%rRightT(iFlag) = H2 * right%rHydCond
            end if
          end if
          vtx%rError(iFlag) = -(vtx%rLeftT(iFlag) - vtx%rRightT(iFlag)) * vtx%rCheckPot(iFlag) - &
              vtx%rRightT(iFlag) * rFDP_PotentialJump(io, vtx%pFDP, vtx%cCPZ(iFlag)) + &
              vtx%rLeftT(iFlag) * vtx%rRightT(iFlag) * &
              (str%rLeftB - str%rRightB + rHALF*(vtx%rLeftH(iFlag) - vtx%rRightH(iFlag)))
        end do
      end do
    end do

    return
  end subroutine AQU_ComputeCheck




  subroutine AQU_Read(io, aqu)
    !! subroutine AQU_Read
    !!
    !! Reads the aquifer information for layer iL from the input LU
    !!
    !! The format of the aqu section of the input file appears as follows:
    !!
    !! AQU ninho base thickness hyd-cond porosity
    !! REF(z0) refhead(Qx0, Qy0)    [ optional ]
    !! BDY nvertices                 [ optional ]
    !!     (x1, y1) (x2, y2) head flux ghb-distance flag
    !!     ...
    !! END
    !! IN0 [optional]
    !! DOM nvertices base thickness hyd-cond porosity
    !!     (x, y)
    !!     ...
    !! END
    !!
    !! NOTE: It is assumed that the AQU line was found by the caller
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! Locals -- for Directive parsing
    type(DIRECTIVE), dimension(7), parameter :: dirDirectives = &
                       (/dirEND, dirREF, dirSWI, dirBDY, dirOPT, dirDOM, dirSTR/)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMax
    real(kind=AE_REAL) :: rBase
    real(kind=AE_REAL) :: rThickness
    real(kind=AE_REAL) :: rHydCond
    real(kind=AE_REAL) :: rPorosity
    real(kind=AE_REAL) :: rRefHead
    complex(kind=AE_REAL) :: cRefPoint
    complex(kind=AE_REAL) :: cRefUniformFlow
    complex(kind=AE_REAL) :: cZ, cZ1, cZ2
    real(kind=AE_REAL) :: rSpecHead
    real(kind=AE_REAL) :: rSpecFlux
    real(kind=AE_REAL) :: rGhbDistance
    real(kind=AE_REAL) :: rInactiveValue
    integer(kind=AE_INT) :: iBdyFlag
    integer(kind=AE_INT) :: iNBdy
    integer(kind=AE_INT) :: iNIsl
    integer(kind=AE_INT) :: iID, iLeftID, iRightID
    logical :: lClosed
    type(AQU_BDYELEMENT), pointer :: vtx
    ! State variables for the parser
    integer(kind=AE_INT) :: iParseMode
    integer(kind=AE_INT), parameter :: PARSE_NONE = 0
    integer(kind=AE_INT), parameter :: PARSE_BDY = 1
    integer(kind=AE_INT), parameter :: PARSE_DOM = 2
    integer(kind=AE_INT), parameter :: PARSE_STR = 3
    type(DOM_DOMAIN), pointer :: dom
    character(len=132) :: sOption
    character(len=32) :: sTag

    call IO_MessageText(io, "  Reading AQU module input")
    call IO_Assert(io, (associated(aqu)), "AQU_Read: AQU_Create has not been called")

    iParseMode = PARSE_NONE
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "AQU_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all Mod "ifXXXRead" routines.
          call IO_Assert(io, .false., "AQU_Read: Unexpected EOF")
        case (kOpData)
          ! A data line was found. If we have a specified perimeter, add the point
          ! to the perimeter.
          select case (iParseMode)
            case (PARSE_BDY)
              call IO_Assert(io, (aqu%iNBdy < size(aqu%BdyElements)), &
                   "AQU_Read: Space exhausted")
              cZ1 = cIO_GetCoordinate(io, 'cZ1', extents=.true.)
              cZ2 = cIO_GetCoordinate(io, 'cZ2', extents=.true.)
              rSpecHead = rIO_GetReal(io, 'rSpecHead')
              rSpecFlux = rIO_GetReal(io, 'rSpecFlux')
              rGhbDistance = rIO_GetReal(io, 'rGhbDistance', minimum = rZERO)
              iBdyFlag = iIO_GetInteger(io, 'iBdyFlag')
              call IO_Assert(io, (iBdyFlag == BDY_HEAD .or. iBdyFlag == BDY_FLUX .or. &
                                  iBdyFlag == BDY_GHB .or. iBdyFlag == BDY_FREESURF), &
                   'Illegal boundary type')
              aqu%iNBdy = aqu%iNBdy+1
              vtx => aqu%BdyElements(aqu%iNBdy)
              vtx%cZ1 = cZ1
              vtx%cZ2 = cZ2
              vtx%rLength = abs(cZ2-cZ1)
              vtx%rSpecHead = rSpecHead
              vtx%rSaveSpecHead = rSpecHead
              vtx%rSpecFlux = rSpecFlux
              vtx%rGhbDistance = rGhbDistance
              vtx%iBdyFlag = iBdyFlag
            case (PARSE_DOM)
              dom => aqu%dom%Domains(aqu%dom%iNDom)
              call IO_Assert(io, (associated(dom%cZ)), "AQU_Read: No DOM directive")
              call IO_Assert(io, (dom%iNPts < size(dom%cZ)), "AQU_Read: DOM space exhausted")
              write (unit=sTag, fmt=*) 'cZ dom', aqu%dom%iNDom, dom%iNPts+1
              if (dom%iNPts > 0) then
                cZ = cIO_GetCoordinate(io, sTag, extents=.true., check_points=dom%cZ(1:dom%iNPts))
              else
                cZ = cIO_GetCoordinate(io, sTag, extents=.true.)
              end if
              dom%iNPts = dom%iNPts + 1
              dom%cZ(dom%iNPts) = cZ
            case (PARSE_STR)
              write (unit=sTag, fmt=*) 'cZ str', aqu%iNStr, aqu%Strings(aqu%iNStr)%iNPts+1
              associate (str => aqu%Strings(aqu%iNStr))
                if (str%iNPts > 0) then
                  cZ = cIO_GetCoordinate(io, sTag, extents=.true., &
                       check_points=str%Vertices(1:str%iNPts)%cZ)
                else
                  cZ = cIO_GetCoordinate(io, sTag, extents=.true.)
                end if
              end associate
              call AQU_AppendStringVertex(io, aqu, cZ)
            case default
              call IO_Assert(io, .false., "AQU_Read: Unexpected data record")
          end select
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpREF)
          ! Reference record was found. Read the reference point and uniform flow rate
          cRefPoint = cIO_GetCoordinate(io, 'cRefPoint', extents=.true.)
          rRefHead = rIO_GetReal(io, 'rRefHead')
          cRefUniformFlow = cIO_GetComplex(io, 'cUniformFlow')
          dom => DOM_FindDomain(io, aqu%dom, cRefPoint)
          call IO_Assert(io, (rRefHead > dom%rBase), "AQU_Read: Base elevation below aquifer top")
          call FRF_SetReference(io, aqu%frf, cRefPoint, rRefHead, cRefUniformFlow)
        case (kOpBDY)
          ! Have we already specified boundary elements?
          call IO_Assert(io, (.not. associated(aqu%BdyElements)), &
               "AQU_Read: BDY statement has already been encountered")
          iMax = iIO_GetInteger(io, 'iMax')
          call IO_Assert(io, (iMax > 0), "AQU_Read: Illegal dimension")
          allocate(aqu%BdyElements(iMax), stat = iStat)
          call IO_Assert(io, (iStat == 0), "AQU_Read: Allocation failed")
          do iNBdy = 1, iMax
            associate (b => aqu%BdyElements(iNBdy))
              b%cZ1 = cZERO;  b%cZ2 = cZERO
              b%cCPZ1 = cZERO; b%cCPZ2 = cZERO
              b%cZC = cZERO;  b%cZL = cZERO
              b%rSpecHead = rZERO; b%rSpecFlux = rZERO
              b%rGhbDistance = rZERO; b%rLength = rZERO
              b%rStrength = rZERO; b%iBdyFlag = BDY_FLUX
              b%rCheckPot = rZERO; b%rCheckFlux = rZERO
              b%lMoveCPZ1 = .false.; b%lMoveCPZ2 = .false.
              b%lFSHeadSpec = .false.
              b%cFSZ1 = cZERO; b%cFSZ2 = cZERO
              b%iEqIndex = 0; b%rSaveSpecHead = rZERO
              nullify(b%pFLS)
            end associate
          end do
          aqu%iNBdy = 0
          iParseMode = PARSE_BDY
        case (kOpOPT)
          sOption = sIO_GetField(io, "sOption", allowed=(/"NEWBOUNDARY"/), force_uppercase=.true.)
          if (sOption == "NEWBOUNDARY") then
            aqu%lNewBdy = .true.
            call IO_MessageText(io, "New boundary condition logic is enabled")
          end if
        case (kOpDOM)
          call IO_Assert(io, associated(aqu%dom), "AQU_Read: No DOM_COLLECTION (AQU not yet created)")
          call DOM_Read(io, aqu%dom)
          iParseMode = PARSE_DOM
        case (kOpSTR)
          call IO_Assert(io, associated(aqu%dom), "AQU_Read: No DOM_COLLECTION (AQU not yet created)")
          iMax = iIO_GetInteger(io, 'iMax')
          iLeftID = iIO_GetInteger(io, 'iLeftID')
          iRightID = iIO_GetInteger(io, 'iRightID')
          lClosed = lIO_GetLogical(io, 'lClosed')
          iID = iIO_GetInteger(io, 'iID', forbidden=aqu%Strings%iID)
          call AQU_BeginString(io, aqu, iMax, iLeftID, iRightID, lClosed, iID)
          iParseMode = PARSE_STR
        case default
      end select
    end do
    call IO_MessageText(io, "  Leaving AQU module")

    return
  end subroutine AQU_Read


  subroutine AQU_BeginString(io, aqu, iMax, iLeftID, iRightID, lClosed, iID)
    !! subroutine AQU_BeginString
    !!
    !! Creates a new inhomogeneity string entry in aqu%Strings.
    !! Called by AEM_ReadInho after reading the STR directive header.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iMax, iLeftID, iRightID, iID
    logical, intent(in) :: lClosed
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iVtx
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx

    call IO_Assert(io, (associated(aqu%Strings)), "AQU_BeginString: No strings have been allocated")
    call IO_Assert(io, (aqu%iNStr < size(aqu%Strings)), "AQU_BeginString: Space exhausted")

    aqu%iNStr = aqu%iNStr + 1
    str => aqu%Strings(aqu%iNStr)
    str%iNPts = 0
    str%iLeftID = iLeftID
    str%iRightID = iRightID
    str%lClosed = lClosed
    str%iID = iID
    allocate(str%Vertices(iMax), stat=iStat)
    call IO_Assert(io, (iStat == 0), "AQU_BeginString: Allocation failed")
    do iVtx = 1, size(str%Vertices)
      vtx => str%Vertices(iVtx)
      vtx%cZ = cZERO
      vtx%rStrength = rZERO
      vtx%rCheckPot(:) = rZERO
      nullify(vtx%pFDP)
    end do

    return
  end subroutine AQU_BeginString


  subroutine AQU_AppendStringVertex(io, aqu, cZ)
    !! subroutine AQU_AppendStringVertex
    !!
    !! Appends one vertex coordinate to the current (last) inhomogeneity string.
    !! Called by AEM_ReadInho for each data line following a STR directive.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(AQU_STRING), pointer :: str

    call IO_Assert(io, (aqu%iNStr > 0), "AQU_AppendStringVertex: No string in progress")
    str => aqu%Strings(aqu%iNStr)
    call IO_Assert(io, (associated(str%Vertices)), "AQU_AppendStringVertex: No STR directive")
    call IO_Assert(io, (str%iNPts < size(str%Vertices)), "AQU_AppendStringVertex: Space exhausted")

    str%iNPts = str%iNPts + 1
    str%Vertices(str%iNPts)%cZ = cZ

    return
  end subroutine AQU_AppendStringVertex


  subroutine AQU_Report(io, aqu)
    ! Reports aquifer information
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN, iBdy
    integer(kind=AE_INT) :: iDom, iStr, iVtx, i
    type(AQU_BDYELEMENT), pointer :: bdy
    type(DOM_DOMAIN), pointer :: dom
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx
    real(kind=AE_REAL) :: rCheckHead
    real(kind=AE_REAL), dimension(2) :: rError

    ! Aquifer data from this module
    call HTML_Header('Module AQU', 1)
    call HTML_Header('Basic aquifer information', 2)

    if (aqu%frf%lReference) then
      call HTML_Header('Reference aquifer', 3)
      call HTML_StartTable()
      call HTML_AttrComplex('Reference point', aqu%frf%cRefPoint)
      call HTML_AttrReal('Reference head', aqu%frf%rRefHead)
      call HTML_AttrComplex('Uniform flow', aqu%frf%cRefUniformFlow)
      call HTML_AttrReal('Solution constant', aqu%frf%rSolConst)
      call HTML_EndTable()
    else
      call HTML_Header('No reference point(continuity of flow condition)', 3)
      call HTML_StartTable()
      call HTML_AttrComplex('Reference point', cIO_WorldCoords(io, aqu%frf%cRefPoint))
      call HTML_AttrReal('Reference head', aqu%frf%rRefHead)
      call HTML_AttrComplex('Uniform flow', aqu%frf%cRefUniformFlow)
      call HTML_AttrReal('Solution constant', aqu%frf%rSolConst)
      call HTML_EndTable()
    end if

    call HTML_Header('Usage', 3)
    call HTML_StartTable()
    call HTML_AttrInteger('Number of boundary elements', aqu%iNBdy)
    call HTML_AttrInteger('Number of FWL functions', iAQU_GetInfo(io, aqu, SIZE_FWL, 0))
    call HTML_AttrInteger('Number of FPD functions', iAQU_GetInfo(io, aqu, SIZE_FPD, 0))
    call HTML_AttrInteger('Number of FDP functions', iAQU_GetInfo(io, aqu, SIZE_FDP, 0))
    call HTML_AttrInteger('Number of equations', iAQU_GetInfo(io, aqu, SIZE_EQUATIONS, 0))
    call HTML_AttrInteger('Number of unknowns', iAQU_GetInfo(io, aqu, SIZE_UNKNOWNS, 0))
    call HTML_EndTable()

    ! Report the boundary elements
    if (aqu%iNBdy > 0) then
      call HTML_Header('Boundary elements', 3)
      call HTML_StartTable()
      call HTML_TableHeader((/'Index   ', 'X1      ', 'Y1      ', 'X2      ', 'Y2      ', 'Flag    ', &
           'Spec Hd ', 'Spec Fl ', 'GHB Dis ', 'Strength', 'M Head  ', 'M Flux  ', &
           'Head Err', 'Flux Err'/))
      do iBdy = 1, aqu%iNBdy
        call HTML_StartRow()
        bdy => aqu%BdyElements(iBdy)
        rCheckHead = rDOM_PotentialToHead(io, aqu%dom, bdy%rCheckPot, rHALF*(bdy%cZ1+bdy%cZ2))
        call HTML_ColumnInteger((/iBdy/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, bdy%cZ1), cIO_WorldCoords(io, bdy%cZ2)/))
        call HTML_ColumnInteger((/bdy%iBdyFlag/))
        call HTML_ColumnReal((/ bdy%rSpecHead, &
                                rIO_WorldLength(io, bdy%rSpecFlux, bdy%cZ2-bdy%cZ1), &
                                bdy%rGhbDistance, bdy%rStrength, &
                                rCheckHead, &
                                rIO_WorldLength(io, bdy%rCheckFlux, bdy%cZ2-bdy%cZ1) /))
        if (bdy%iBdyFlag == BDY_HEAD) then
          call HTML_ColumnReal((/rCheckHead - bdy%rSpecHead/))
          call HTML_ColumnText((/'--'/))
        else if (bdy%iBdyFlag == BDY_FLUX) then
          call HTML_ColumnText((/'--'/))
          call HTML_ColumnReal((/bdy%rCheckFlux - bdy%rSpecFlux/))
        else if (bdy%iBdyFlag == BDY_GHB) then
          call HTML_ColumnText((/'--'/))
          call HTML_ColumnText((/'--'/))
        else if (bdy%iBdyFlag == BDY_FREESURF) then
          if (bdy%lFSHeadSpec) then
            call HTML_ColumnReal((/rCheckHead - bdy%rSpecHead/))
            call HTML_ColumnText((/'--'/))
          else
            call HTML_ColumnText((/'--'/))
            call HTML_ColumnReal((/bdy%rCheckFlux - bdy%rSpecFlux/))
          end if
        end if
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No boundary elements', 3)
    end if

    ! Write out the inhomogeneities (was IN0_Report)
    call HTML_Header('Module IN0 (Inhomogeneities)', 1)
    call HTML_Header('Inhomogeneity information', 2)

    if (.not. associated(aqu%dom%Domains)) then
      call HTML_Header('No domains allocated', 3)
    else
      call HTML_StartTable()
      call HTML_AttrInteger('Number of domains', aqu%dom%iNDom)
      call HTML_AttrInteger('Number of strings', aqu%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iAQU_GetInfo(io, aqu, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iAQU_GetInfo(io, aqu, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iAQU_GetInfo(io, aqu, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iAQU_GetInfo(io, aqu, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iAQU_GetInfo(io, aqu, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iDom = 1, aqu%dom%iNDom
        call HTML_Header('Domain Information', 3)
        call HTML_StartTable()
        dom => aqu%dom%Domains(iDom)
        call HTML_AttrInteger('Domain number ', iDom)
        call HTML_AttrInteger('ID', dom%iID)
        call HTML_AttrReal('Base elevation', dom%rBase)
        call HTML_AttrReal('Thickness', dom%rThickness)
        call HTML_AttrReal('Hydraulic conductivity', dom%rHydCond)
        call HTML_AttrReal('Porosity', dom%rPorosity)
        call HTML_EndTable()

        if (.not. associated(dom%cZ)) then
          call HTML_Header('No vertices -- infinite aquifer', 4)
        else
          call HTML_Header('Vertices', 4)
          call HTML_StartTable()
          call HTML_TableHeader((/'X', 'Y'/))
          do iVtx = 1, dom%iNPts
            call HTML_StartRow()
            call HTML_ColumnComplex((/cIO_WorldCoords(io, dom%cZ(iVtx))/))
            call HTML_EndRow()
          end do
          call HTML_EndTable()
        end if
      end do

      do iStr = 1, aqu%iNStr
        call HTML_Header('String Information', 3)
        call HTML_StartTable()
        str => aqu%Strings(iStr)
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrInteger('Left domain', str%iLeftID)
        call HTML_AttrInteger('Right domain', str%iRightID)
        call HTML_AttrLogical('Closed?', str%lClosed)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)
        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'X     ', 'Y     ', 'FDP # ', 'V1 Str', 'C Str ', 'V1 Err', 'C Err '/))
        do iVtx = 1, str%iNPts
          vtx => str%Vertices(iVtx)
          do i = 1, 2
            if (vtx%rStrength(i) /= rZERO) then
              rError(i) = vtx%rError(i) / vtx%rStrength(i)
            else
              rError(i) = rZERO
            end if
          end do
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx/))
          call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
          if (associated(vtx%pFDP)) then
            call HTML_ColumnInteger((/vtx%pFDP%iIndex/))
          else
            call HTML_ColumnText((/'--'/))
          end if
          call HTML_ColumnReal(vtx%rStrength(1:2), 'e13.6')
          call HTML_ColumnReal(rError, 'e13.6')
          call HTML_EndRow()
        end do
        call HTML_EndTable()
      end do

    end if

    return
  end subroutine AQU_Report


  subroutine AQU_Inquiry(io, aqu, iLU, lCSV)
    !! subroutine AQU_Inquiry
    !!
    !! Writes an inquiry report for all barriers to iLU
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iLU
    logical, intent(in), optional :: lCSV
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iBdy, iStr, iVtx, i
    real(kind=AE_REAL) :: rLength, rARecip
    real(kind=AE_REAL), dimension(2) :: rError
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx, next
    logical :: lDoCSV
    lDoCSV = .false.
    if (present(lCSV)) lDoCSV = lCSV

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_Inquiry: AQU_Create has not been called")
    end if

    if (aqu%frf%lReference) then
      if (lDoCSV) then
        write (unit=iLU, &
               fmt="(""tag, ibdy, istr, ref_x, ref_y, const, spec_head, model_head, error"")")
      else
        write (unit=iLU, &
               fmt="(""#[REFERENCE]AQU, 0, 0, REF_X, REF_Y, CONST, SPEC_HEAD, MODEL_HEAD, ERROR"")")
      end if
      write (unit=iLU, &
             fmt="(""AQU"", 2("", "", i9), 5("", "", e16.8))" &
             ) 0, &
             0, &
             cIO_WorldCoords(io, aqu%frf%cRefPoint), &
             aqu%frf%rSolConst, &
             aqu%frf%rRefHead, &
             rDOM_PotentialToHead(io, aqu%dom, aqu%frf%rCheck, aqu%frf%cRefPoint)
    end if

    if (aqu%iNBdy > 0) then
      if (lDoCSV) then
        write (unit=iLU, fmt="()")
        write (unit=iLU, &
               fmt="(""tag, ibdy, flag, x1, y1, x2, y2, length, strength, spec_head, spec_flux, " // &
               "ghb_distance, model_head, model_flux"")")
      else
        write (unit=iLU, &
               fmt="(""#[BOUNDARY]AQU, IBDY, FLAG, X1, Y1, X2, Y2, LENGTH, STRENGTH, SPEC_HEAD, SPEC_FLUX, " // &
               "GHB_DISTANCE, MODEL_HEAD, MODEL_FLUX"")")
      end if
      do iBdy = 1, aqu%iNBdy
        bdy => aqu%BdyElements(iBdy)
        write (unit=iLU, &
               fmt="(""BDY"", 2("", "", i9), 11("", "", e16.8))" &
               ) iBdy, &
               bdy%iBdyFlag, &
               cIO_WorldCoords(io, bdy%cZ1), &
               cIO_WorldCoords(io, bdy%cZ2), &
               bdy%rLength, bdy%rStrength, &
               bdy%rSpecHead, &
               bdy%rSpecFlux, &
               bdy%rGhbDistance, &
               rDOM_PotentialToHead(io, aqu%dom, bdy%rCheckPot, rHALF*(bdy%cZ1+bdy%cZ2)), &
               rIO_WorldLength(io, bdy%rCheckFlux, bdy%cZ2-bdy%cZ1)
      end do
    end if

    ! Inhomogeneity inquiry (was IN0_Inquiry)
    if (io%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_Inquiry: AQU_Create has not been called")
    end if

    if (lDoCSV) then
      write (unit=iLU, fmt="()")
      write (unit=iLU, &
             fmt="(""tag, id, vtx, x1, y1, x2, y2, length, strength1, strength2, pot1, pot2, error1, error2"")")
    else
      write (unit=iLU, &
             fmt="(""#IN0, ID, VTX, X1, Y1, X2, Y2, LENGTH, STRENGTH1, STRENGTH2, POT1, POT1, ERROR1, ERROR2"")")
    end if
    do iStr = 1, aqu%iNStr
      str => aqu%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        vtx => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        if (iVtx < str%iNPts) then
          rLength = abs(str%Vertices(iVtx+1)%cZ - vtx%cZ)
        else
          rLength = rZERO
        end if
        rARecip = rAQU_ARecip(io, aqu, iStr)
        do i = 1, 2
          if (vtx%rStrength(i) /= rZERO) then
            rError(i) = vtx%rError(i) / vtx%rStrength(i)
          else
            rError(i) = rZERO
          end if
        end do
        write (unit=iLU, &
               fmt="('IN0', 2(', ', i9), 11(', ', e16.8))" &
               ) str%iID, iVtx, cIO_WorldCoords(io, vtx%cZ), cIO_WorldCoords(io, next%cZ), &
                 rLength, vtx%rStrength(1), vtx%rStrength(2), vtx%rCheckPot(1), vtx%rCheckPot(2), &
                 rError(1), rError(2)
      end do
    end do

    return
  end subroutine AQU_Inquiry


  !! FUNCTION MODULE ROUTINES
  !! These routines allow AQU to behave as a ModAEM Function Module


  function cAQU_Potential(io, aqu, cZ) result(cOmega)
    !! complex function cAQU_Potential
    !!
    !! Computes the complex potential due to the reference aquifer at cZ
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    cOmega = cFRF_Potential(io, aqu%frf, cZ)

    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cZ-bdy%cZC) / bdy%cZL
      cTemp = cILS_InfluenceP(cBigZ, bdy%cZL)
      cOmega = cOmega + cTemp(1, 1) * bdy%rStrength
    end do

  end function cAQU_Potential


  function cAQU_Discharge(io, aqu, cZ) result(cQ)
    !! complex function cAQU_Discharge
    !!
    !! Computes the complex discharge due to the reference aquifer at cZ
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    cQ = cFRF_Discharge(io, aqu%frf, cZ)

    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cZ-bdy%cZC) / bdy%cZL
      cTemp = cILS_InfluenceW(cBigZ, bdy%cZL)
      cQ = cQ + conjg((cTemp(1, 1) * bdy%rStrength)/bdy%cZL)
    end do

    return
  end function cAQU_Discharge


  function rAQU_Flow(io, aqu, cPathZ) result(rFlow)
    !! real function rAQU_Flow
    !!
    !! Computes the integrated flow due to the current set of dipoles.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy, iStat, i
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    ! The reference flow field satisfies Laplace's equation, so only endpoint
    ! values of the potential are needed for its flow contribution.
    rFlow = rFRF_Flow(io, aqu%frf, cPathZ)

    allocate(cBigZ(size(cPathZ)), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'rAQU_Flow: Allocation failed')
    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cPathZ-bdy%cZC) / bdy%cZL
      do i = 1, size(cBigZ)-1
        cTemp = cILS_InfluenceF(cBigZ(i), cBigZ(i+1), bdy%cZL)
        rFlow = rFlow + real(cTemp(1, 1) * bdy%rStrength)
      end do
    end do
    deallocate(cBigZ)

    return
  end function rAQU_Flow


  function rAQU_Recharge(io, aqu, cZ) result(rG)
    !! complex function rAQU_Recharge
    !!
    !! Computes the recharge rate due to the aquifer
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG
    ! [ LOCALS ]

    rG = rZERO

    return
  end function rAQU_Recharge


  function rAQU_Extraction(io, aqu) result(rQ)
    !! complex function rAQU_Extraction
    !!
    !! Computes the extraction due to the aquifer
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (.not. aqu%frf%lReference), "AQU_Extraction: Invalid for problems with uniform flow")
    rQ = sum(aqu%BdyElements(1:aqu%iNBdy)%rCheckFlux)

    return
  end function rAQU_Extraction


  function rAQU_BranchCut(io, aqu, cPathZ) result(rBC)
    !! real function rAQU_BranchCut
    !!
    !! Computes the net branch cut.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBC
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy, iStat, i
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    rBC = rZERO

    allocate(cBigZ(size(cPathZ)), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'rAQU_Flow: Allocation failed')
    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cPathZ-bdy%cZC) / bdy%cZL
      do i = 1, size(cBigZ)-1
        cTemp = cILS_InfluenceB(cBigZ(i), cBigZ(i+1), bdy%cZL)
        rBC = rBC + aimag(cTemp(1, 1) * bdy%rStrength)
      end do
    end do
    deallocate(cBigZ)

    return
  end function rAQU_BranchCut


  function lAQU_CheckPoint(io, aqu, cZ, rTol, cZFix, rStrength, iElementType, iElementString, &
             iElementVertex, iElementFlag) result(lFound)
    !! function lAQU_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance of the
    !! perimeter of a pond.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
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
    type(AQU_BDYELEMENT), pointer :: bdy
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    real(kind=AE_REAL) :: rCheckTol
    complex(kind=AE_REAL) :: cBZFix
    ! [ LOCAL PARAMETERS ]
    real(kind=AE_REAL), parameter :: rEND_FIX_POSITION = rTWO

    if (aqu%iNBdy > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a perimeter line sink
      allocate(cMapZ(aqu%iNBdy))
      cMapZ(:) = (cZ-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL

      lFound = .false.
      cZFix = cZERO
      do i = 1, aqu%iNBdy
        bdy => aqu%BdyElements(i)
        if (abs(aimag(cMapZ(i))) < rCheckTol) then
          ! We're near a segment. Are we at an end?
          if (abs(real(cMapZ(i))-cONE) < rCheckTol) then
            ! We're at the upper tip
            lFound = .true.
            cBZFix = cmplx(rONE-rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * bdy%cZL + bdy%cZC
            rStrength = bdy%rStrength
            iElementType = ELEM_AQU
            iElementString = kAQUBoundary
            iElementVertex = i
            iElementFlag = -1
            exit
          else if (abs(real(cMapZ(i))+cONE) < rCheckTol) then
            ! We're at the lower tip
            lFound = .true.
            cBZFix = cmplx(-rONE+rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * bdy%cZL + bdy%cZC
            rStrength = bdy%rStrength
            iElementType = ELEM_AQU
            iElementString = kAQUBoundary
            iElementVertex = i
            iElementFlag = -1
            exit
          end if
        end if
      end do
      deallocate(cMapZ)
    else
      lFound = .false.
      cZFix = cZERO
    end if

    return
  end function lAQU_CheckPoint


  function lAQU_CheckIntersection(io, aqu, cZ1, cZ2, cZInt, cZBefore, cZAfter, &
             iElementType, iElementString, iElementVertex, iElementFlag, cENormal) result(lFound)
    !! function lAQU_CheckIntersection
    !!
    !! Checks the specified line segment Z1-Z2 and returns .true. if there is an intersection.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), intent(out) :: cZInt, cZBefore, cZAfter
    integer(kind=AE_INT), intent(out) :: iElementType
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    complex(kind=AE_REAL), intent(out) :: cENormal
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lFound
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iLoc
    type(AQU_BDYELEMENT), pointer :: bdy
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1, cMapZ2, cMapZInt
    real(kind=AE_REAL), dimension(:), allocatable :: rDist
    complex(kind=AE_REAL) :: cMyZ1, cMyZ2, cZFix
    integer(kind=AE_INT) :: iEType, iEStr, iEVtx, iEFlg
    real(kind=AE_REAL) :: rStr
    real(kind=AE_REAL) :: rCheckTol
    ! [ LOCAL PARAMETERS ]
    real(kind=AE_REAL), parameter :: rPOINT_TOLERANCE_FACTOR = rTWO
    real(kind=AE_REAL), parameter :: rINTERSECTION_TOLERANCE = 1.0e-4_AE_REAL

    ! First, make sure we're not in tolerance
    if (lAQU_CheckPoint(io, aqu, cZ1, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ1 = cZFix
    else
      cMyZ1 = cZ1
    end if
    if (lAQU_CheckPoint(io, aqu, cZ2, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ2 = cZFix
    else
      cMyZ2 = cZ2
    end if

    if (aqu%iNBdy > 0) then
      ! Find the smallest distance to a dipole
      allocate(cMapZ1(aqu%iNBdy), cMapZ2(aqu%iNBdy), cMapZInt(aqu%iNBdy), rDist(aqu%iNBdy))
      cMapZ1(:) = (cMyZ1-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL
      cMapZ2(:) = (cMyZ2-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL
      ! Find the points of intersection with the line extensions
      cMapZInt = cmplx(rHUGE, rHUGE, AE_REAL)
      where (sign(rONE, aimag(cMapZ1)) /= sign(rONE, aimag(cMapZ2)))
        cMapZInt = cMapZ1 + abs(aimag(cMapZ1))/(abs(aimag(cMapZ1))+abs(aimag(cMapZ2))) * (cMapZ2-cMapZ1)
      end where
      ! Find the distances to line(s) that intersect
      rDist = rHUGE
      where (abs(real(cMapZInt)) <= rONE)
        rDist = abs(cMapZInt-cMapZ1)
      end where
      ! Now, find the CLOSEST intersection point
      iLoc = minloc(rDist, DIM = 1)
      if (rDist(iLoc) < rHUGE) then
        bdy => aqu%BdyElements(iLoc)
        ! Reverse-map to get the intersection
        cZInt = cMapZInt(iLoc) * bdy%cZL + bdy%cZC
        cZBefore = cZInt - rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        cZAfter = cZInt + rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        if (abs(cZ1-cZInt) < abs(cZBefore-cZInt)) cZBefore = cZ1
        if (abs(cZ2-cZInt) < abs(cZAfter-cZInt)) cZAfter = cZ2
        lFound = .true.
        iElementType = ELEM_AQU
        iElementString = 1
        iElementVertex = iLoc
        iElementFlag = -1
        cENormal = cI * bdy%cZL / abs(bdy%cZL)
      else
        lFound = .false.
        cZBefore = cZERO
        cZAfter = cZERO
        iElementType = -1
        iElementString = -1
        iElementVertex = -1
        iElementFlag = -1
        cENormal = cZERO
      end if
      deallocate(cMapZ1, cMapZ2, cMapZInt, rDist)
    else
      lFound = .false.
      cZBefore = cZERO
      cZAfter = cZERO
      iElementType = -1
      iElementString = -1
      iElementVertex = -1
      iElementFlag = -1
      cENormal = cZERO
    end if

    return
  end function lAQU_CheckIntersection


  !! UTILITY ROUTINES
  !! These routines provide computational aids for AEM models


  function rAQU_ActiveArea(io, aqu, poly) result(rArea)
    !! real function rAQU_ActiveArea
    !!
    !! Computes the area of 'poly' that lies within the active area of the aquifer.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rArea

    rArea = PGN_Area(poly)

    return
  end function rAQU_ActiveArea


  function rAQU_InterfaceElevation(io, aqu, cZ, rPot) result(rIfcElev)
    !! real function rAQU_InterfaceElevation
    !!
    !! Interface calculations (not yet available)
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rIfcElev

    call IO_MessageText(io, "Interface calculations are not yet available")
    rIfcElev = rDOM_Base(io, aqu%dom, cZ)

    return
  end function rAQU_InterfaceElevation


  function lAQU_CheckActive(io, aqu, cZ) result(lActive)
    !! function lAQU_CheckActive
    !!
    !! Returns .true. if cZ is in the active area of the model
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lActive

    if (aqu%iNPerimeter==0) then
      lActive = .true.
    else
      lActive = (PGN_CheckPoint(aqu%cPerimeter(1:aqu%iNPerimeter), cZ) > 0)
    end if

    return
  end function lAQU_CheckActive


  function rAQU_ARecip(io, aqu, iString) result(rARecip)
    !! function rAQU_ARecip
    !!
    !! Returns the factor A*, defined as the reciprocal of(k+ - k-) / k-
    !! where k+ is the inside conductivity and k- is the outside conductivity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iString
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rARecip
    ! [ LOCALS ]
    type(AQU_STRING), pointer :: str
    type(DOM_DOMAIN), pointer :: left, right

    if (io%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "rAQU_ARecip: AQU_Create has not been called")
      call IO_Assert(io, (iString >= 1 .and. iString <= aqu%dom%iNDom), &
           "rAQU_ARecip: Bad domain index")
    end if

    str => aqu%Strings(iString)
    left => DOM_FindDomainID(io, aqu%dom, str%iLeftID)
    right => DOM_FindDomainID(io, aqu%dom, str%iRightID)

    rARecip = right%rHydCond / (left%rHydCond - right%rHydCond)

    return
  end function rAQU_ARecip


  subroutine AQU_Save(io, aqu, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: ibdy, istr, ivtx, iflg
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx

    ! Write the solution constant
    if (mode == IO_MODE_BINARY) then
      write (unit=LU_SCRATCH) ELEM_AQU, 0, 1, 1, aqu%frf%rSolConst
    else
      write (unit=LU_SCRATCH, fmt=*) "AQU", 0, 1, 1, aqu%frf%rSolConst
    end if
    ! Output records will be of the form ELEM_AQU, IBDY, 1, 1, SIGMA
    do ibdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(ibdy)
      if (mode == IO_MODE_BINARY) then
        write (unit=LU_SCRATCH) ELEM_AQU, ibdy, 1, 1, bdy%rStrength
      else
        write (unit=LU_SCRATCH, fmt=*) "AQU", ibdy, 1, 1, bdy%rStrength
      end if
    end do
    ! Save the inhomogeneities
    do istr = 1, aqu%iNStr
      str => aqu%Strings(istr)
      do ivtx = 1, str%iNPts
        vtx => str%Vertices(ivtx)
        do iflg = 1, ubound(vtx%rStrength, 1)
          if (mode == IO_MODE_BINARY) then
            write (unit=LU_SCRATCH) ELEM_IN0, istr, ivtx, iflg, vtx%rStrength(iflg)
          else
            write (unit=LU_SCRATCH, fmt=*) "IN0", istr, ivtx, iflg, vtx%rStrength(iflg)
          end if
        end do
      end do
    end do

    return
  end subroutine AQU_Save


  subroutine AQU_Load(io, aqu, fdp, mode)
    !! Loads the AQU and IN0 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, ibdy, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(AQU_BDYELEMENT), pointer :: bdy
    type(AQU_STRING), pointer :: str
    type(AQU_VERTEX), pointer :: vtx

    ! Scans the entire precondition file for the AQU data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, ibdy, ivtx, iflg, rstrength
        if (istat < 0) exit
        call IO_Assert(io, istat == 0, "I/O error on precondition file")
        if (imodule == ELEM_AQU) then
          if (ibdy == 0) then
            aqu%frf%rSolConst = rstrength
          else
            call IO_Assert(io, ibdy > 0 .and. ibdy <= aqu%iNBdy, "AQU bdyl not found")
            bdy => aqu%BdyElements(ibdy)
            bdy%rStrength = rstrength
          end if
        else if (imodule == ELEM_IN0) then
          call IO_Assert(io, ibdy > 0 .and. ibdy <= aqu%iNStr, "IN0 string not found")
          str => aqu%Strings(ibdy)
          call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "IN0 vertex not found")
          vtx => str%Vertices(ivtx)
          call IO_Assert(io, iflg > 0 .and. iflg <= ubound(vtx%rStrength, 1), "IN0 strength index not found")
          vtx%rStrength(iflg) = rstrength
        end if
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, ibdy, ivtx, iflg, rstrength
        if (istat < 0) exit
        call IO_Assert(io, istat == 0, "I/O error on precondition file")
        if (uppercase(trim(smodule)) == "AQU") then
          if (ibdy == 0) then
            aqu%frf%rSolConst = rstrength
          else
            call IO_Assert(io, ibdy > 0 .and. ibdy <= aqu%iNBdy, "AQU bdyl not found")
            bdy => aqu%BdyElements(ibdy)
            bdy%rStrength = rstrength
          end if
        else if (uppercase(trim(smodule)) == "IN0") then
          call IO_Assert(io, ibdy > 0 .and. ibdy <= aqu%iNStr, "IN0 string not found")
          str => aqu%Strings(ibdy)
          call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "IN0 vertex not found")
          vtx => str%Vertices(ivtx)
          call IO_Assert(io, iflg > 0 .and. iflg <= ubound(vtx%rStrength, 1), "IN0 strength index not found")
          vtx%rStrength(iflg) = rstrength
        end if
      end if
    end do

    ! Populate the internal data structures
    call AQU_Update(io, aqu, fdp)

    return
  end subroutine AQU_Load



end module p_aqu
