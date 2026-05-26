module p_ls1

  ! ModAEM 2.0
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
  !!   f_linesink   --  Function module for collections of line-sinks
  !!
  !! This module provides the necessary functionality for head-specified
  !! line-sink elements. Elements are defined as strings of points, with the
  !! specified head at each vertex.

  use u_constants
  use u_io
  use f_linesink
  use f_aem
  use u_matrix

  implicit none

  public

  type :: LS1_VERTEX
    !! type LS1_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZ
    !!     The complex coordinate of the vertex
    !!   real :: rHead
    !!     The specified head at the vertex
    !!   real :: rStrength
    !!     The sink density (sigma) for the segment starting at this vertex
    !!   real :: rLength
    !!     The segment length
    !!   type(FLS_LINESINK), pointer :: pFLS
    !!     Pointer to the vertex entry in the FLS module. Note: nullified for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rHead
    real(kind=AE_REAL) :: rStrength
    real(kind=AE_REAL) :: rLength
    type(FLS_LINESINK), pointer :: pFLS
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
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(LS1_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iNPts
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
    !! subroutine LS1_New
    !!
    !! Adds a new LS1_STRING object to the LS1_COLLECTION 'ls1'
    !!
    !! Calling Sequence:
    !!    call LS1_New(io, ls1, Vertices, iNPts)
    !!
    !! Arguments:
    !!    (in)    type(LS1_COLLECTION), pointer :: ls1
    !!              The LS1_COLLECTION object to be used
    !!    (in)    type(LS1_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the string
    !!    (in)    integer :: iNPts
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
    !! Returns the ID number for the line-sink string at index 'iIndex'
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
      case (SIZE_FLS)
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


  subroutine LS1_SetupFunctions(io, ls1, fls)
    !! subroutine LS1_SetupFunctions
    !!
    !! This routine sets up the functions in f_linesink for the line-sinks.
    !! All linesinks are created with zero strength; strengths are determined
    !! by the matrix solution and applied via LS1_Update.
    !!
    !! Calling Sequence:
    !!    call LS1_SetupFunctions(io, ls1, fls)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FLS_COLLECTION), pointer
    !!             FLS_COLLECTION function object
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(FLS_COLLECTION), pointer :: fls
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_SetupFunctions: LS1_Create has not been called")
      call IO_Assert(io, (associated(fls)), &
           "LS1_SetupFunctions: Illegal FLS_COLLECTION object")
    end if

    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        this_vtx%rLength = abs(next_vtx%cZ - this_vtx%cZ)
        this_vtx%pFLS => FLS_New(io, fls, this_vtx%cZ, next_vtx%cZ, cZERO, ELEM_LS1, iStr, iVtx, -1)
      end do
    end do

    return
  end subroutine LS1_SetupFunctions


  subroutine LS1_SetupMatrix(io, ls1, mat)
    !! subroutine LS1_SetupMatrix
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
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


  function iLS1_Prepare(io, ls1, iIteration) result(iChanges)
    !! subroutine LS1_Prepare
    !!
    !! Prepares the module for a new iteration -- do-nothing for LS1
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
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


  subroutine LS1_ComputeCoefficients(io, ls1, fls, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine LS1_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the LS1
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call LS1_ComputeCoefficients(io, ls1, fls, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FLS_COLLECTION), pointer
    !!             FLS_COLLECTION object to be used
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
    type(FLS_COLLECTION), pointer :: fls
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
    integer(kind=AE_INT) :: iStat, iCol, iStr, iVtx, iNLS
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cLSF, cLSW
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: first_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_ComputeCoefficients: LS1_Create has not been called")
      call IO_Assert(io, (associated(fls)), &
           "LS1_ComputeCoefficients: Illegal FLS_COLLECTION object")
    end if

    iCol = 0
    rARow = rZERO
    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      first_vtx => str%Vertices(1)
      ! ASSUMES that LS1_SetupFunctions created consecutive linesink entries
      iNLS = str%iNPts-1
      allocate(cLSF(0:iNLS, 1, 1), cLSW(0:iNLS, 1, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS1_ComputeCoefficients: Allocation failed")

      ! Get the appropriate influence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FLS_GetInfluence(io, fls, INFLUENCE_P, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_BDYGHB)
          call FLS_GetInfluence(io, fls, INFLUENCE_P, first_vtx%pFLS, iNLS, (/rHALF*sum(cPathZ)/), cOrientation, cLSF(1:iNLS, :, :))
          call FLS_GetInfluence(io, fls, INFLUENCE_F, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSW(1:iNLS, :, :))
          cLSF = cLSF + rGhbResistance*cLSW
        case (EQN_FLOW)
          call FLS_GetInfluence(io, fls, INFLUENCE_F, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_INHO)
          call FLS_GetInfluence(io, fls, INFLUENCE_P, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_DISCHARGE)
          call FLS_GetInfluence(io, fls, INFLUENCE_W, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_RECHARGE)
          call FLS_GetInfluence(io, fls, INFLUENCE_G, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_CONTINUITY)
          call FLS_GetInfluence(io, fls, INFLUENCE_Q, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_POTENTIALDIFF)
          call FLS_GetInfluence(io, fls, INFLUENCE_D, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
        case (EQN_TOTALFLOW)
          call FLS_GetInfluence(io, fls, INFLUENCE_Z, first_vtx%pFLS, iNLS, cPathZ, cOrientation, cLSF(1:iNLS, :, :))
      end select

      do iVtx = 1, iNLS
        iCol = iCol+1
        rARow(iCol) = real(cLSF(iVtx, 1, 1))
      end do

      deallocate(cLSF, cLSW)
    end do

    rARow = rARow * rMultiplier

    return
  end subroutine LS1_ComputeCoefficients


  function rLS1_ComputeRHS(io, ls1, aem, iEqType, iElementType, iElementString, iElementVertex, &
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
    type(AEM_DOMAIN), pointer :: aem
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

    if (lDirect) then
      rRHS = rDOM_HeadToPotential(io, aem%dom, vtx%rCPHead, vtx%cCPZ(1))
    else
      rRHS = rDOM_HeadToPotential(io, aem%dom, vtx%rCPHead, vtx%cCPZ(1)) - vtx%rCheckPot
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
    type(LS1_VERTEX), pointer :: this_vtx

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
    if (lDirect) then
      this_vtx%rStrength = rValue
    else
      this_vtx%rStrength = this_vtx%rStrength + rValue
    end if

    return
  end subroutine LS1_StoreResult


  subroutine LS1_Update(io, ls1, fls)
    !! subroutine LS1_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    LS1_Update(io, ls1, fls)
    !!
    !! Arguments:
    !!   (in)    type(LS1_COLLECTION), pointer
    !!             LS1_COLLECTION object to be used
    !!   (in)    type(FLS_COLLECTION), pointer
    !!             FLS_COLLECTION function object
    !!
    ! [ ARGUMENTS ]
    type(LS1_COLLECTION), pointer :: ls1
    type(FLS_COLLECTION), pointer :: fls
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Update: LS1_Create has not been called")
      call IO_Assert(io, (associated(fls)), &
           "LS1_Update: Illegal FLS_COLLECTION object")
    end if

    do iStr = 1, ls1%iNStr
      str => ls1%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        vtx => str%Vertices(iVtx)
        vtx%pFLS%cSigma = cmplx(vtx%rStrength, rZERO, AE_REAL)
      end do
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


  subroutine LS1_ComputeCheck(io, ls1, aem)
    !! Updates check potential and head for all LS1 linesink segments.
    type(LS1_COLLECTION), pointer :: ls1
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    type(LS1_VERTEX), pointer :: vtx
    integer(kind=AE_INT) :: iStr, iVtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), "LS1_ComputeCheck: LS1_Create has not been called")
    end if

    do iStr = 1, ls1%iNStr
      do iVtx = 1, ls1%Strings(iStr)%iNPts - 1
        vtx => ls1%Strings(iStr)%Vertices(iVtx)
        vtx%rCheckPot = real(cAEM_Potential(io, aem, vtx%cCPZ(1), .false.), AE_REAL)
        vtx%rCheckHead = rDOM_PotentialToHead(io, aem%dom, vtx%rCheckPot, vtx%cCPZ(1))
      end do
    end do

    return
  end subroutine LS1_ComputeCheck


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
          this_vtx%rStrength = rZERO
          if (str%iNPts > 1) then
            last_vtx%rLength = abs(this_vtx%cZ - last_vtx%cZ)
          end if
          nullify(this_vtx%pFLS)
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
          str%iID = iID
          str%iNPts = 0      ! Initialize the vertex counter
      end select
    end do

    call IO_MessageText(io, "  Leaving LS1 module")

  end subroutine LS1_Read


  subroutine LS1_Inquiry(io, ls1, iLU, lCSV)
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
    logical, intent(in), optional :: lCSV
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rLength, rSigma
    type(LS1_STRING), pointer :: str
    type(LS1_VERTEX), pointer :: this, next
    logical :: lDoCSV
    lDoCSV = .false.
    if (present(lCSV)) lDoCSV = lCSV

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls1)), &
           "LS1_Inquiry: LS1_Create has not been called")
    end if

    if (lDoCSV) then
      write (unit=iLU, fmt="(""tag, id, vtx, x1, y1, x2, y2, length, spec_head, strength, mod_head, error"")")
    else
      write (unit=iLU, &
             fmt="(""#LS1, ID, VTX, X1, Y1, X2, Y2, LENGTH, SPEC_HEAD, STRENGTH, MOD_HEAD, ERROR"")")
    end if
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
      call HTML_AttrInteger('Number of FLS functions', iLS1_GetInfo(io, ls1, SIZE_FLS, 0))
      call HTML_AttrInteger('Number of equations', iLS1_GetInfo(io, ls1, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iLS1_GetInfo(io, ls1, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iStr = 1, ls1%iNStr
        str => ls1%Strings(iStr)
        call HTML_Header('Line-sink string definition', 3)
        call HTML_StartTable()
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'FLS # ', 'X     ', 'Y     ', 'Head  ', 'Sigma ', 'Check ', 'M Head', 'Error '/))
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          next => str%Vertices(iVtx+1)
          rSigma = rIO_WorldLength(io, vtx%rStrength, next%cZ-vtx%cZ)
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, vtx%pFLS%iIndex/))
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


  subroutine LS1_Load(io, ls1, fls, mode)
    !! Loads the LS1 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(LS1_COLLECTION), pointer :: ls1
    type(FLS_COLLECTION), pointer :: fls
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
    call LS1_Update(io, ls1, fls)

    return
  end subroutine LS1_Load

end module p_ls1
