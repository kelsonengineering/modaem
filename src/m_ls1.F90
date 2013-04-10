module m_ls1

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

  !! module m_ls1
  !!
  !! Element module for 2-D head specified line-sinks
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! line-sink elements. Elements are defined as strings of points, with the
  !! discharge rate(per length) specified for each segment.
  !!
  !! Note: ModAEM uses the "traditional" line-sink function (io, Strack, 1989)
  !! for matrix generation, but uses strings of dipoles terminated by a well
  !! for computational performance, once a solution is achieved.

  use u_constants
  use u_io
  use f_well
  use f_dipole
  use u_matrix
  use m_aqu

  implicit none

  public

  type :: LS1_VERTEX
    !! type LS1_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rHead
    !!     The specified head at the vertex
    !!   real :: rDPStrength
    !!     The dipole strength at the vertex
    !!   real :: rLength
    !!     The segment length
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rHead
    real(kind=AE_REAL) :: rDPStrength
    real(kind=AE_REAL) :: rStrength
    real(kind=AE_REAL) :: rLength
    integer(kind=AE_INT) :: iFDPIndex
    complex(kind=AE_REAL), dimension(1) :: cCPZ
    real(kind=AE_REAL) :: rCPHead
    real(kind=AE_REAL) :: rCheckPot
    real(kind=AE_REAL) :: rCheckHead
  end type LS1_VERTEX

  type :: LS1_STRING
    !! type LS1_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   type(LS1_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of LS1_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(LS1_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iFWLIndex
    integer(kind=AE_INT) :: iID
  end type LS1_STRING

  type :: LS1_COLLECTION
    !! type LS1_COLLECTION
    !!
    !! Type that holds information for all LS1 elements in a layer
    !!
    !! Members:
    !!   type(LS1_STRING), dimension(:), pointer :: Strings
    !!     A vector of LS1_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !!
    type(LS1_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    ! Iterator Information
    integer(kind=AE_INT) :: iIterStr
    integer(kind=AE_INT) :: iIterVtx
  end type LS1_COLLECTION

  ! Module flags for matrix generator routines
  integer(kind=AE_INT), private, parameter :: LS1Vertex = 1


contains

  !**pd Modified to allocate only the collection object

  function LS1_Create(io) result(ls1)
    !! function LS1_Create
    !!
    !! Creates a new LS1_COLLECTION object
    !!
    !! Calling Sequence:
    !!    ls1 => LS1_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(ls1, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS1_Create: allocation failed")
    nullify(ls1%Strings)
    ls1%iNStr = 0

    return
  end function LS1_Create


  subroutine LS1_Alloc(io, ls1)
    !! Subroutine LS1_Alloc
    !!
    !! Allocates Strings for the LS1_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call LS1_Alloc(io, ls1, iNStr)
    !!
    !! Arguments:
    !!    (in)    type(LS1_COLLECTION), pointer :: ls0
    !!              The LS1_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iStat

    iNStr = iIO_GetInteger(io, 'iNStr', minimum = 0)
    allocate(ls1%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS1_Alloc: allocation failed")

  end subroutine LS1_Alloc


  subroutine LS1_Destroy(io, ls1)
    !! subroutine LS1_Destroy
    !!
    !! Frees memory allocated for ls1 Linesinks and strings of vertices
    !! and the ls1 Collection object
    !!
    !! Calling Sequence:
    !!     call LS1_Destroy(ls1)
    !!
    !! Arguments:
    !!  type(LS1_COLLECTION), pointer :: ls1
    !!              Pointer to the LS1_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: i
    type(LS1_STRING), pointer :: str


    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Destroy: LS1_Create has not been called")
    end if

    ! First deallocate each string of vertices
    do i = 1, ls1%iNStr
      str => ls1%Strings(i)
      deallocate(str%Vertices, stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS1_Destroy: deallocation of Vertices failed")
    end do
    ! Then deallocate the strings table
    if (associated(ls1%Strings)) then
      deallocate(ls1%Strings, stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS1_Destroy: deallocation of Strings failed")
    end if
    ! Now deallocate the collection
    deallocate(ls1, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS1_Destroy: deallocation failed")

    return
  end subroutine LS1_Destroy



  subroutine LS1_New(io, ls1, Vertices, iNPts)
    !! function LS1_New
    !!
    !! Adds a new LS1_STRING object to the LS1_COLLECTION 'ls1'
    !!
    !! Calling Sequence:
    !!    call LS1_New(io, ls1, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(LS1_COLLECTION), pointer :: ls1
    !!              The LS1_COLLECTION object to be used
    !!    (in)    type(LS1_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(LS1_VERTEX), dimension(:) :: Vertices
    integer(kind=AE_INT), intent(in) :: iNPts
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(LS1_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_New: LS1_Create has not been called")
    end if

    call IO_Assert(io, (ls1%iNStr < size(ls1%Strings)), &
         "LS1_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(Vertices)), &
         "LS1_New: Size of provided vertices is inconsistent")

    ls1%iNStr = ls1%iNStr + 1
    str => ls1%Strings(ls1%iNStr)
    allocate(str%Vertices(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS1_New: Allocation failed")
    str%Vertices = Vertices(1:iNPts)
    str%iNPts = iNPts

    return
  end subroutine LS1_New


  function iLS1_GetID(io, ls1, iIndex) result(iID)
    !! Returns the ID number for the well at index 'iIndex'
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= ls1%iNStr), "Internal error -- no such index")
    iID = ls1%Strings(iIndex)%iID

    return
  end function iLS1_GetID


  subroutine LS1_PreSolve(io, ls1)
    !! subroutine LS1_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call LS1_PreSolve(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io

    return
  end subroutine LS1_PreSolve


  function iLS1_GetInfo(io, ls1, iOption, iIteration) result(iValue)
    !! function LS1_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iLS1_GetInfo(io, ls1, iOption)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr
    type(LS1_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_GetInfo: LS1_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = ls1%iNStr
      case (SIZE_FDP)
        do iStr = 1, ls1%iNStr
          str => ls1%Strings(iStr)
          iValue = iValue + str%iNPts-1
        end do
      case (SIZE_EQUATIONS)
        do iStr = 1, ls1%iNStr
          str => ls1%Strings(iStr)
          iValue = iValue + str%iNPts-1
        end do
      case (SIZE_UNKNOWNS)
        do iStr = 1, ls1%iNStr
          str => ls1%Strings(iStr)
          iValue = iValue + str%iNPts-1
        end do
      case (INFO_REGENERATE)
        if (iIteration < 2) then
          iValue = 1
        else
          iValue = 0
        end if
      case default
        iValue = 0
    end select

    return
  end function iLS1_GetInfo


  subroutine LS1_SetupFunctions(io, ls1, fwl, fdp)
    !! subroutine LS1_Setup
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call LS1_Setup(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, i, iDP, iWL
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, next_vtx, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Setup: LS1_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "LS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "LS1_Setup: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      ! Build dipoles for all segments
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        this_vtx%rLength = abs(next_vtx%cZ - this_vtx%cZ)
        call FDP_New(io, fdp, this_vtx%cZ, next_vtx%cZ, (/cZERO, cZERO, cZERO/), ELEM_LS1, iStr, iVtx, -1, this_vtx%iFDPIndex)
      end do

      ! Put a well at the end of the string
      last_vtx => str%Vertices(str%iNPts)
      call FWL_New(io, fwl, last_vtx%cZ, rZERO, rZERO, ELEM_LS1, iStr, -1, -1, iWL)
      str%iFWLIndex = iWL
    end do

    return
  end subroutine LS1_SetupFunctions


  subroutine LS1_SetupMatrix(io, ls1, aqu, mat)
    !! subroutine LS1_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the linesinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call LS1_Setup(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, i, iDP, iWL, iEQ
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, next_vtx, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Setup: LS1_Create has not been called")
      call IO_Assert(io, (associated(aqu)), &
           "LS1_Setup: Illegal AQU_COLLECTION object")
      call IO_Assert(io, (associated(mat)), &
           "LS1_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix generator entries for all segments
    do iStr = 1, ls1%iNStr
      ! Set up the unknown variables
      ! Vertex entry -- all vertices
      str => ls1%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        call MAT_CreateVariable(io, mat, ELEM_LS1, iStr, iVtx, LS1Vertex)
      end do

      ! Set up control points and equations -- One equation per segment
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        rHead1 = this_vtx%rHead
        rHead2 = next_vtx%rHead
        ! Compute control points for the segment.  There is one unknown
        ! per segment.
        call MAT_ComputeControlPoints(io, this_vtx%cZ, next_vtx%cZ, 1, cCPResult, rZERO)
        ! Now, create the equation entry...
        this_vtx%cCPZ(1) = cCPResult(2)
        this_vtx%rCPHead = rHead1 + (rHead2-rHead1)*(this_vtx%cCPZ(1)-this_vtx%cZ)/(next_vtx%cZ-this_vtx%cZ)
        iEQ = MAT_CreateEquation(io, mat, this_vtx%cCPZ, EQN_HEAD, ELEM_LS1, iStr, iVtx, 0, cZERO, rZERO)
      end do
    end do

    return
  end subroutine LS1_SetupMatrix


  function iLS1_Prepare(io, ls1, aqu, iIteration) result(iChanges)
    !! subroutine LS1_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_ls1
    !!
    !! Calling Sequence:
    !!    call LS1_Setup(io, ls1, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iLS1_Prepare


  function rLS1_GetCoefficientMultiplier(io, ls1, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    rMultiplier = rONE

    return
  end function rLS1_GetCoefficientMultiplier


  subroutine LS1_ComputeCoefficients(io, ls1, fwl, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine LS1_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the LS1
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call LS1_ComputeCoefficients(io, ls1, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
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
    type(LS1_COLLECTION), pointer :: ls1
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
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
    integer(kind=AE_INT) :: iStat, iCol, iStr, iVtx, iDP1, iNDP
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDPW
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: first_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_ComputeCoefficients: LS1_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "LS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "LS1_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iCol = 0
    rARow = rZERO
    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      first_vtx => str%Vertices(1)
      ! ASSUMES that LS1_Setup routine created consecutive dipole entries
      iDP1 = first_vtx%iFDPIndex
      iNDP = str%iNPts-1
      allocate(cDPF(0:iNDP, 1, 1), cDPW(0:iNDP, 1, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS1_ComputeCoefficients: Allocation failed")

      ! Get the appropriate incluence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_BDYGHB)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, (/rHALF*sum(cPathZ)/), cOrientation, cDPF(1:iNDP, :, :))
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPW(1:iNDP, :, :))
          cDPF = cDPF + rGhbResistance*cDPW
        case (EQN_FLOW)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_INHO)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_DISCHARGE)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_W, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_RECHARGE)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_G, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_CONTINUITY)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_Q, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_POTENTIALDIFF)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_D, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_TOTALFLOW)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_Z, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
      end select

      do iVtx = 1, iNDP
        iCol = iCol+1
        rARow(iCol) = real(cDPF(iVtx, 1, 1))
      end do

      deallocate(cDPF, cDPW)
    end do

    rARow = rARow * rMultiplier

    return
  end subroutine LS1_ComputeCoefficients


  function rLS1_ComputeRHS(io, ls1, aqu, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rLS1_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rLS1_ComputeRHS(io, ls1, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag, rSpecValue)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
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
    type(LS1_COLLECTION), pointer :: ls1
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
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: vtx

    call IO_Assert(io, (iElementString > 0 .and. iElementString <= ls1%iNStr), &
         'Illegal string index')
    str => ls1%Strings(iElementString)
    call IO_Assert(io, (iElementVertex > 0 .and. iElementVertex <= str%iNPts), &
         'Illegal vertex index')
    vtx => str%Vertices(iElementVertex)

    ! For LS1, compute the RHS by subtracting the previous result from the
    ! desired potential at the control-point
    if (lDirect) then
      rRHS = rAQU_HeadToPotential(io, aqu, vtx%rCPHead, vtx%cCPZ(1))
    else
      rRHS = rAQU_HeadToPotential(io, aqu, vtx%rCPHead, vtx%cCPZ(1)) - vtx%rCheckPot
    end if
    return
  end function rLS1_ComputeRHS


  subroutine LS1_StoreResult(io, ls1, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine LS1_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the LS1 module.
    !!
    !! Calling Sequence:
    !!    LS1_StoreResult(io, ls1, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(always ELEM_LS1)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_StoreResult: LS1_Create has not been called")
      call IO_Assert(io, (iElementString >= 1 .and. iElementString <= ls1%iNStr), &
           "LS1_StoreResult: Bad element string ID")
    end if

    str => ls1%Strings(iElementString)

    if (io%lDebug) then
      call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= str%iNPts-1), &
           "LS1_StoreResult: Bad element vertex ID")
    end if

    this_vtx => str%Vertices(iElementVertex)
    next_vtx => str%Vertices(iElementVertex+1)
    if (lDirect) then
      this_vtx%rStrength = rValue
    else
      this_vtx%rStrength = this_vtx%rStrength + rValue
    end if
    next_vtx%rDPStrength = this_vtx%rDPStrength + this_vtx%rStrength * this_vtx%rLength

    return
  end subroutine LS1_StoreResult


  subroutine LS1_Update(io, ls1, fwl, fdp)
    !! subroutine LS1_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    LS1_Update(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, next_vtx


    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Update: LS1_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "LS1_Update: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        cRho1 = cmplx(this_vtx%rDPStrength, rZERO, AE_REAL)
        cRho3 = cmplx(next_vtx%rDPStrength, rZERO, AE_REAL)
        cRho2 = rHALF * (cRho1+cRho3)
        call FDP_Update(io, fdp, this_vtx%iFDPIndex, (/cRho1, cRho2, cRho3/))
      end do
      ! Put a well at the end of the string
      call FWL_Update(io, fwl, str%iFWLIndex, real(cRho3))
    end do

    return
  end subroutine LS1_Update


  subroutine LS1_FindStringPointer(io, ls1, iLSID, LSString, lFound)
    !! subroutine LS1_FindStringPointer
    !!
    !! Finds the linesink string specified by the ID and returns a pointer to it
    !!
    !! Calling Sequence:
    !!    call LS1_FindStringPointer(io, ls1, iLSID, LSString, lfound)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (in)    integer :: iLSID
    !!             The linesink string ID number
    !!   (out)   type(LS0_STRING) :: LSString
    !!             Pointer to the linesink string
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: iLSID
    type(LS1_STRING), pointer :: LSString
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS0_FindStringPointer: LS1_Create has not been called")
    end if

    lFound = .false.
    do i = 1, ls1%iNStr
      LSString => ls1%Strings(i)
      if (LSString%iID == iLSID) then
        lFound = .true.
        return
      end if
    end do

    return
  end subroutine LS1_FindStringPointer


  subroutine LS1_ResetIterator(io, ls1)
    !! subroutine LS1_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call LS1_ResetIterator(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls1
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_ResetIterator: LS1_Create has not been called")
    end if

    ls1%iIterStr = 1
    ls1%iIterVtx = 0

    return
  end subroutine LS1_ResetIterator


  function LS1_NextIterator(io, ls1) result(itr)
    !! function LS1_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS1_NextIterator(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls1
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_NextIterator: LS1_Create has not been called")
    end if

    if (ls1%iIterStr > ls1%iNStr) then
      nullify(itr)
      return
    end if

    ls1%iIterVtx = ls1%iIterVtx + 1
    if (ls1%iIterVtx > ls1%Strings(ls1%iIterStr)%iNPts-1) then
      ls1%iIterStr = ls1%iIterStr+1
      ls1%iIterVtx = 1
      if (ls1%iIterStr > ls1%iNStr) then
        nullify(itr)
        return
      end if
    end if

    allocate(itr)
    itr%iElementType = ELEM_AQU
    itr%iElementString = ls1%iIterStr
    itr%iElementVertex = ls1%iIterVtx
    itr%iValueSelector = VALUE_POTENTIAL
    allocate(itr%cZ(1))
    itr%cZ(1) = ls1%Strings(ls1%iIterStr)%Vertices(ls1%iIterVtx)%cCPZ(1)

    return
  end function LS1_NextIterator


  subroutine LS1_SetIterator(io, ls1, aqu, itr, cValue)
    !! function LS1_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS1_SetIterator(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be used
    !!   (in)    type(ITERATOR_RESULT), pointer :: itr
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: ls1
    !!             Tracks error conditions
    !!     Pointer to the information for data retrieval
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(AQU_COLLECTION), pointer :: aqu
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    type(IO_STATUS), pointer :: io
    type(LS1_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_NextIterator: LS1_Create has not been called")
      call IO_Assert(io, (ls1%iIterStr <= ls1%iNStr), &
           "LS1_SetIterator: Iterator out of range")
    end if

    vtx => ls1%Strings(itr%iElementString)%Vertices(itr%iElementVertex)
    vtx%rCheckPot = real(cValue, AE_REAL)
    vtx%rCheckHead = rAQU_PotentialToHead(io, aqu, vtx%rCheckPot, vtx%cCPZ(1))

    return
  end subroutine LS1_SetIterator


  subroutine LS1_Read(io, ls1)
    !! subroutine LS1_Read
    !!
    !! Populates an LS1_COLLECTION using data from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call LS1_Read(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer :: ls1
    !!             LS1_COLLECTION to be populated
    !!
    !! The format of the LS1 section of the input file appears as follows:
    !! LS1
    !! STR NVertices
    !! x y head
    !! ... Up to NVertices
    !! STR NVertices
    !! x y head
    !! ... Up to NVertices
    !! ... Up to NStrings
    !! END
    !!
    !! NOTE: It is assumed that the LS1 line was found by the caller

    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(2), parameter :: dirDirectives = (/dirEND, dirSTR/)
    ! [ LOCALS ]
    character(len=132) :: sTag
    real(kind=AE_REAL) :: rHead
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxStr, iMaxVtx
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iStr
    logical :: lFlag
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, last_vtx

    call IO_MessageText(io, "  Reading LS1 module input")

    call IO_Assert(io, (associated(ls1)), "LS1_Read: LS1_Create has not been called")

    ! Use IO_InputRecord to process the model input file.
    nullify(str, this_vtx, last_vtx)
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "LS1_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModLS1 "ifXXXRead" routines.
          call IO_Assert(io, .false., "LS1_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(str)), "LS1_Read: No current string")
          call IO_Assert(io, (str%iNPts < size(str%Vertices)), "LS1_Read: Space exhausted")
          write (unit=sTag, fmt=*) 'cZ ', ls1%iNStr, str%iNPts+1
          cZ = cIO_GetCoordinate(io, sTag, extents=.true.)
          rHead = rIO_GetReal(io, 'rHead')

          str%iNPts = str%iNPts+1
          this_vtx => str%Vertices(str%iNPts)
          this_vtx%cZ = cZ
          this_vtx%rHead = rHead
          this_vtx%rDPStrength = rZERO
          this_vtx%rStrength = rZERO
          if (str%iNPts > 1) then
            last_vtx%rLength = abs(this_vtx%cZ - last_vtx%cZ)
          end if
          this_vtx%iFDPIndex = -1
          last_vtx => this_vtx
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpSTR)
          !****************************************************************************
          ! Here for the STR command -- create a new string of line-sinks
          ! the maximum number of vertices is in the input record
          !****************************************************************************
          call IO_Assert(io, (ls1%iNStr < size(ls1%Strings)), "LS1_Read: Space exhausted")
          iMaxVtx = iIO_GetInteger(io, 'iMaxVtx', minimum = 2)
          iID = iIO_GetInteger(io, 'iID')
          ! OKAY! Allocate the vertices...
          ls1%iNStr = ls1%iNStr+1
          str => ls1%Strings(ls1%iNStr)
          allocate(str%Vertices(iMaxVtx), stat = iStat)
          call IO_Assert(io, (iStat == 0), "LS1_Read: Allocation failed")
          ! Made it!
          str%iFWLIndex = -1         ! No FWL function yet!
          str%iID = iID
          str%iNPts = 0      ! Initialize the vertex counter
      end select
    end do

    call IO_MessageText(io, "  Leaving LS1 module")

  end subroutine LS1_Read


  subroutine LS1_Inquiry(io, ls1, iLU)
    !! subroutine LS1_Inquiry
    !!
    !! Writes an inquiry report for all line-sinks to iLU
    !!
    !! Calling Sequence:
    !!    call LS1_Inquiry(io, ls1, iLU)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rLength, rSigma
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Inquiry: LS1_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#LS1, ID, VTX, X1, Y1, X2, Y2, LENGTH, SPEC_HEAD, STRENGTH, MOD_HEAD, ERROR"")")
    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        rLength = rIO_WorldLength(io, this%rLength, next%cZ-this%cZ)
        rSigma = rIO_WorldLength(io, this%rStrength, next%cZ-this%cZ)
        write (unit=iLU, &
               fmt="(""LS1"", 2("", "", i9), 9("", "", e16.8))" &
               ) str%iID, &
               iVtx, &
               cIO_WorldCoords(io, this%cZ), &
               cIO_WorldCoords(io, next%cZ), &
               rLength, &
               this%rCPHead, &
               rSigma, &
               this%rCheckHead, &
               this%rCheckHead-this%rCPHead
      end do
    end do

    return
  end subroutine LS1_Inquiry


  subroutine LS1_Report(io, ls1)
    !! subroutine LS1_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call LS1_Report(ls1)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: vtx, next
    real(kind=AE_REAL) :: rSigma

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Report: LS1_Create has not been called")
    end if

    call HTML_Header('Module LS1', 1)
    call HTML_Header('Head-specified area sink information', 2)

    if (associated(ls1%Strings)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of strings', ls1%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iLS1_GetInfo(io, ls1, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iLS1_GetInfo(io, ls1, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iLS1_GetInfo(io, ls1, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iLS1_GetInfo(io, ls1, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iLS1_GetInfo(io, ls1, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iStr = 1, ls1%iNStr
        str => ls1%Strings(iStr)
        call HTML_Header('Line-sink string definition', 3)
        call HTML_StartTable()
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrInteger('FWL index', str%iFWLIndex)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'FDP # ', 'X     ', 'Y     ', 'Head  ', 'Sigma ', 'Check ', 'M Head', 'Error '/))
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          next => str%Vertices(iVtx+1)
          rSigma = rIO_WorldLength(io, vtx%rStrength, next%cZ-vtx%cZ)
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, vtx%iFDPIndex/))
          call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
          call HTML_ColumnReal((/vtx%rCPHead, rSigma, vtx%rCheckPot, vtx%rCheckHead, vtx%rCheckHead-vtx%rCPHead/))
          call HTML_EndRow()
        end do
        vtx => str%Vertices(str%iNPts)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iVtx/))
        call HTML_ColumnText((/'--'/))
        call HTML_ColumnComplex((/vtx%cZ/))
        call HTML_ColumnText((/'--', '--', '--', '--', '--'/))
        call HTML_EndRow()
        call HTML_EndTable()
      end do
    else
      call HTML_Header('No line-sinks defined', 3)
    end if

    return
  end subroutine LS1_Report


  subroutine LS1_Save(io, ls1, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(LS1_COLLECTION), pointer :: ls1
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istr, ivtx, iflg
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: vtx

    ! Output records will be of the form ELEM_LS1, IWEL, IRAD, IVTX, 0, SIGMA
    do istr = 1, ls1%iNStr
      str => ls1%Strings(istr)
      do ivtx = 1, str%iNPts
        vtx => str%Vertices(ivtx)
        if (mode == IO_MODE_BINARY) then
          write (unit=LU_SCRATCH) ELEM_LS1, istr, ivtx, 1, vtx%rStrength
        else
          write (unit=LU_SCRATCH, fmt=*) "LS1", istr, ivtx, 1, vtx%rStrength
        end if
      end do
    end do

    return
  end subroutine LS1_Save


  subroutine LS1_Load(io, ls1, fwl, fdp, mode)
    !! Loads the LS1 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(LS1_COLLECTION), pointer :: ls1
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, istr, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: vtx
    type(LS1_VERTEX), pointer :: flg

    ! Scans the entire precondition file for the LS1 data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, istr, ivtx, iflg, rstrength
        if (imodule /= ELEM_LS1) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, istr, ivtx, iflg, rstrength
        if (uppercase (trim(smodule)) /= "LS1") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      call IO_Assert(io, istr > 0 .and. istr <= ls1%iNStr, "LS1 string not found")
      str => ls1%Strings(istr)
      call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "LS1 vertex not found")
      vtx => str%Vertices(ivtx)
      call IO_Assert(io, iflg == 1, "LS1 strength index not found")
      vtx%rStrength = rstrength
    end do

    ! Now, populate the internal data structures
    call LS1_Update(io, ls1, fwl, fdp)

    return
  end subroutine LS1_Load

end module m_ls1
