module m_hb0

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

  !! module m_hb0
  !!
  !! Element module for 2-D horizontal no-flow boundaries
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for no-flow strings
  !! (not closed domains at this time). The strings are constructed using
  !! line-doublets(implemented as line-dipoles of imaginary strength).

  use u_constants
  use u_io
  use u_matrix
  use f_dipole

  implicit none

  public

  type :: HB0_VERTEX
    !! type HB0_VERTEX
    !!
    !! Type that holds information for one vertex along a no-flow string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rStrength(1)
    !!     The doublet strength at the vertex
    !!   real :: rStrength(2)
    !!     The doublet strength at the center of the next segment
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iFDPIndex
    real(kind=AE_REAL), dimension(2) :: rStrength
    real(kind=AE_REAL), dimension(2) :: rCheck
    complex(kind=AE_REAL), dimension(:), pointer :: cVertexCPZ
    complex(kind=AE_REAL), dimension(:), pointer :: cCenterCPZ
  end type HB0_VERTEX

  type :: HB0_STRING
    !! type HB0_STRING
    !!
    !! Type that holds information for one no-flow string
    !!
    !! Members:
    !!   type(HB0_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of HB0_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(HB0_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iID
  end type HB0_STRING

  type :: HB0_COLLECTION
    !! type HB0_COLLECTION
    !!
    !! Type that holds information for all HB0 elements in a layer
    !!
    !! Members:
    !!   type(HB0_STRING), dimension(:), pointer :: Strings
    !!     A vector of HB0_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !!
    type(HB0_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    ! Iterator information
    integer(kind=AE_INT) :: iIterStr
    integer(kind=AE_INT) :: iIterVtx
    integer(kind=AE_INT) :: iIterFlag
  end type HB0_COLLECTION

  ! Matrix generator element flags
  integer(kind=AE_INT), private, parameter :: kHB0_Vertex = 1
  integer(kind=AE_INT), private, parameter :: kHB0_Center = 2

  real(kind=AE_REAL), parameter :: MOVEFACTOR = 1.0001_AE_REAL
  real(kind=AE_REAL), parameter :: NORMAL_OFFSET = rZERO


contains


  function HB0_Create(io) result(hb0)
    !! function HB0_Create
    !!
    !! Creates a new HB0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    hb0 => HB0_Create()
    !!
    !! Arguments:
    !!
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    type(HB0_COLLECTION), pointer :: hb0
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(hb0,stat = iStat)
    call IO_Assert(io, (iStat == 0), "HB0_Create: allocation failed")
    nullify(hb0%Strings)
    hb0%iNStr = 0

    return
  end function HB0_Create


  subroutine HB0_Alloc(io, hb0)
    !! Subroutine HB0_Alloc
    !!
    !! Allocates Strings for the HB0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call HB0_Alloc(io, hb0,iNStr)
    !!
    !! Arguments:
    !!    (in)    type(HB0_COLLECTION), pointer :: ls0
    !!              The HB0_COLLECTION object to be used
    !!    (in)    integer :: iNStr
    !!              The number of strings to make space for
    !!
    !!
    !! Return Value:
    !!
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iNStr

    iNStr = iIO_GetInteger(io, "iNStr", minimum=0)
    allocate(hb0%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "HB0_Alloc: allocation failed, ")
    hb0%Strings(:)%iID = -1

    return
  end subroutine HB0_Alloc


  subroutine HB0_Destroy(io, hb0)
    !! subroutine HB0_Destroy
    !!
    !! Frees memory allocated for HB0 horizontal no-flow boundaries
    !! and strings of vertices and the HB0 Collection object
    !!
    !! Calling Sequence:
    !!     call HB0_Destroy(ls0)
    !!
    !! Arguments:
    !!  type(HB0_COLLECTION), pointer :: hb0
    !!              Pointer to the HB0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: i
    type(HB0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Destroy: HB0_Create has not been called")
    end if

    ! First deallocate each string of vertices
    do i = 1, hb0%iNStr
      str => hb0%Strings(i)
      deallocate(str%Vertices, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "hb0_Destroy: deallocation of Vertices failed")
    end do
    ! Then deallocate the strings table
    if (associated(hb0%Strings)) then
      deallocate(hb0%Strings, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "HB0_Destroy: deallocation of Strings failed")
    end if
    ! Now deallocate the collection
    deallocate(hb0,stat = iStat)
    call IO_Assert(io, (iStat == 0), "HB0_Destroy: deallocation failed")

    return
  end subroutine HB0_Destroy


  subroutine HB0_New(io, hb0,Vertices, iNPts)
    !! function HB0_New
    !!
    !! Adds a new HB0_STRING object to the HB0_COLLECTION 'hb0'
    !!
    !! Calling Sequence:
    !!    call HB0_New(io, hb0,Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(HB0_COLLECTION), pointer :: hb0
    !!              The HB0_COLLECTION object to be used
    !!    (in)    type(HB0_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(HB0_VERTEX), dimension(:) :: Vertices
    integer(kind=AE_INT), intent(in) :: iNPts
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(HB0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_New: HB0_Create has not been called")
    end if

    call IO_Assert(io, (hb0%iNStr < size(hb0%Strings)), &
         "HB0_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(Vertices)), &
         "HB0_New: Size of provided vertices is inconsistent")

    hb0%iNStr = hb0%iNStr + 1
    str => hb0%Strings(hb0%iNStr)
    allocate(str%Vertices(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "HB0_New: Allocation failed")
    str%Vertices = Vertices(1:iNPts)
    str%iNPts = iNPts

    return
  end subroutine HB0_New


  subroutine HB0_PreSolve(io, hb0)
    !! subroutine HB0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call HB0_PreSolve(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io

    return
  end subroutine HB0_PreSolve


  function iHB0_GetInfo(io, hb0,iOption, iIteration) result(iValue)
    !! function HB0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iHB0_GetInfo(io, hb0,iOption)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr
    type(HB0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_GetInfo: HB0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FDP)
        do iStr = 1, hb0%iNStr
          str => hb0%Strings(iStr)
          iValue = iValue + str%iNPts-1
        end do
      case (SIZE_EQUATIONS)
        do iStr = 1, hb0%iNStr
          str => hb0%Strings(iStr)
          iValue = iValue + 2*str%iNPts-1
        end do
      case (SIZE_UNKNOWNS)
        do iStr = 1, hb0%iNStr
          str => hb0%Strings(iStr)
          iValue = iValue + 2*str%iNPts-1
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
  end function iHB0_GetInfo


  subroutine HB0_SetupFunctions(io, hb0,fdp)
    !! subroutine HB0_SetupFunctions
    !!
    !! This routine sets up the functions in f_well and f_dipole for the no-flows
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call HB0_Setup(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3
    complex(kind=AE_REAL), dimension(6) :: cCPResult1, cCPResult2
    complex(kind=AE_REAL), dimension(3) :: cCPVtx
    character(len=255) :: sBuf
    integer(kind=AE_INT) :: irv
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Setup: HB0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "HB0_Setup: Illegal FDP_COLLECTION object")
    end if

    ! Build dipoles for all segments
    do iStr = 1, hb0%iNStr
      str => hb0%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        vtx => str%Vertices(iVtx)
        cZ1 = vtx%cZ
        cZ2 = str%Vertices(iVtx+1)%cZ
        call FDP_New(io, fdp, cZ1, cZ2, (/cZERO, cZERO, cZERO/), ELEM_HB0, iStr, iVtx, -1, vtx%iFDPIndex)
      end do
    end do

    return
  end subroutine HB0_SetupFunctions


  subroutine HB0_SetupMatrix(io, hb0,mat)
    !! subroutine HB0_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the module
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call HB0_Setup(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer :: mat
    !!             MAT_MATRIX object to be used
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3
    complex(kind=AE_REAL), dimension(6) :: cCPResult1, cCPResult2
    complex(kind=AE_REAL), dimension(3) :: cCPVtx
    character(len=255) :: sBuf
    integer(kind=AE_INT) :: irv
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx
    integer(kind=AE_INT) :: iStat, iEQ

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Setup: HB0_Create has not been called")
      call IO_Assert(io, (associated(mat)), &
           "HB0_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix entries for all segments
    do iStr = 1, hb0%iNStr
      str => hb0%Strings(iStr)
      ! Set up the unknown variables
      ! Vertex entry -- all vertices except the first and last
      ! No dipole strength at either end of the string
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        ! Make a doublet vertex variable
        call MAT_CreateVariable(io, mat, ELEM_HB0, iStr, iVtx, kHB0_Vertex)
        ! Make a doublet center variable for all but the last entry...
        if (iVtx < str%iNPts) then
          call MAT_CreateVariable(io, mat, ELEM_HB0, iStr, iVtx, kHB0_Center)
        end if
      end do

      ! Set up control point sets and equations
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        ! Now, create the equation entries...
        ! Two equations only for vertices 1 - #Vertices
        if (iVtx == 1) then
          cZ1 = vtx%cZ
          cZ2 = str%Vertices(iVtx+1)%cZ
          ! Compute control intervals
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 2, cCPResult1, NORMAL_OFFSET, &
               (/1.0e-7_AE_REAL, 0.1_AE_REAL, &
               0.6_AE_REAL, rONE-1.0e-7_AE_REAL/))

          ! Vertex entry
          allocate(vtx%cVertexCPZ(2), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_SetupMatrix: Space Exhausted")
          vtx%cVertexCPZ = (/cCPResult1(1), cCPResult1(2)/)
          iEQ = MAT_CreateEquation(io, mat, vtx%cVertexCPZ, EQN_FLOW, ELEM_HB0, &
                                   iStr, iVtx, kHB0_Vertex, cZ2-cZ1, rZERO)

          ! Doublet center control interval
          allocate(vtx%cCenterCPZ(2), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_SetupMatrix: Space Exhausted")
          vtx%cCenterCPZ = (/cCPResult1(2), cCPResult1(3)/)
          iEQ = MAT_CreateEquation(io, mat, vtx%cCenterCPZ, EQN_FLOW, ELEM_HB0, &
                                   iStr, iVtx, kHB0_Center, cZ2-cZ1, rZERO)

        else if (iVtx == str%iNPts) then
          ! ONLY create a vertex control interval for the last vertex
          cZ1 = str%Vertices(iVtx-1)%cZ
          cZ2 = vtx%cZ

          ! Compute control intervals
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 2, cCPResult1, NORMAL_OFFSET, &
               (/1.0e-7_AE_REAL, 0.4_AE_REAL, &
               0.9_AE_REAL, rONE-1.0e-7_AE_REAL/))

          ! Vertex entry
          allocate(vtx%cVertexCPZ(2), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_SetupMatrix: Space Exhausted")
          vtx%cVertexCPZ = (/cCPResult1(3), cCPResult1(4)/)
          iEQ = MAT_CreateEquation(io, mat, vtx%cVertexCPZ, EQN_FLOW, ELEM_HB0, &
                                   iStr, iVtx, kHB0_Vertex, cZ2-cZ1, rZERO)

        else
          ! Compute control points for the segment.  There are two unknowns
          ! per segment, and we need two control points per segment(no overspecification
          ! in this version)
          cZ1 = str%Vertices(iVtx-1)%cZ
          cZ2 = vtx%cZ
          cZ3 = str%Vertices(iVtx+1)%cZ
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 2, cCPResult1, NORMAL_OFFSET, &
               (/1.0e-7_AE_REAL, 0.25_AE_REAL, &
               0.75_AE_REAL, rONE-1.0e-7_AE_REAL/))
          call MAT_ComputeControlPoints(io, cZ2, cZ3, 2, cCPResult2, NORMAL_OFFSET, &
               (/1.0e-7_AE_REAL, 0.25_AE_REAL, &
               0.75_AE_REAL, rONE-1.0e-7_AE_REAL/))

          ! Vertex strength control interval is based on neighboring doublets...
          allocate(vtx%cVertexCPZ(4), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_SetupMatrix: Space Exhausted")
          vtx%cVertexCPZ = (/cCPResult1(3), cCPResult1(4), cCPResult2(1), cCPResult2(2)/)
          iEQ = MAT_CreateEquation(io, mat, vtx%cVertexCPZ, EQN_FLOW, ELEM_HB0, iStr, iVtx, &
                                   kHB0_Vertex, cZ2-cZ1, rZERO)

          ! Center equations
          allocate(vtx%cCenterCPZ(2), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_SetupMatrix: Space Exhausted")
          vtx%cCenterCPZ = (/cCPResult2(2), cCPResult2(3)/)
          iEQ = MAT_CreateEquation(io, mat, vtx%cCenterCPZ, EQN_FLOW, ELEM_HB0, iStr, iVtx, &
                                   kHB0_Center, cZ2-cZ1, rZERO)

        end if
      end do
    end do
  end subroutine HB0_SetupMatrix


  function iHB0_Prepare(io, hb0,iIteration) result(iChanges)
    !! subroutine HB0_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_hb0
    !!
    !! Calling Sequence:
    !!    call HB0_Setup(io, hb0,aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer
    !!             HB0_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iHB0_Prepare


  function rHB0_GetCoefficientMultiplier(io, hb0,iElementString, iElementVertex, &
             iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    rMultiplier = rONE

    return
  end function rHB0_GetCoefficientMultiplier


  subroutine HB0_ComputeCoefficients(io, hb0,fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine HB0_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the HB0
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call HB0_ComputeCoefficients(io, hb0,cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
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
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
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
    integer(kind=AE_INT) :: iStat, iCol, iStr, iVtx, iDP1, iNDP, iWhich, irv
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDPW
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_ComputeCoefficients: HB0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "HB0_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iCol = 0
    rARow = rZERO
    do iStr = 1, hb0%iNStr
      str => hb0%Strings(iStr)
      ! Assume: the HB0_Setup routine creates consecutive dipole entries
      iDP1 = str%Vertices(1)%iFDPIndex
      iNDP = str%iNPts-1
      allocate(cDPF(0:iNDP+1, 3, 1), cDPW(0:iNDP+1, 3, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "HB0_ComputeCoefficients: Allocation failed")

      ! Get the appropriate influence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_BDYGHB)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, (/rHALF*sum(cPathZ)/), &
               cOrientation, cDPF(1:iNDP, :, :))
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, &
               cDPW(1:iNDP, :, :))
          cDPF = cDPF + rGhbResistance*cDPW
        case (EQN_FLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_INHO)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_DISCHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_W, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_RECHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_G, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_CONTINUITY)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Q, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_POTENTIALDIFF)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_D, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
        case (EQN_TOTALFLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Z, iDP1, iNDP, cPathZ, cOrientation, &
               cDPF(1:iNDP, :, :))
      end select

      ! Now, compute the matrix coefficients
      cDPF(0, :, :) = cZERO
      cDPF(iNDP+1, :, :) = cZERO
      do iVtx = 1, iNDP+1
        ! Vertex coefficient
        iCol = iCol+1
        rARow(iCol) = -aimag(cDPF(iVtx-1, 3, 1) + cDPF(iVtx, 1, 1))
        ! Center coefficient for all but the last
        if (iVtx <= iNDP) then
          iCol = iCol+1
          rARow(iCol) = -aimag(cDPF(iVtx, 2, 1))
        end if
      end do

      ! No memory leaks, please!
      deallocate(cDPF, cDPW)
    end do

    rARow = rARow * rMultiplier

    return
  end subroutine HB0_ComputeCoefficients


  function rHB0_ComputeRHS(io, hb0,iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rHB0_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rHB0_ComputeRHS(io, hb0,iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    integer :: iElementType
    !!             Element type(either ELAM_AQU or ELEM_HB0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    !! Return Value:
    !!   real :: rRHS
    !!     The RHS value for the module
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
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
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    call IO_Assert(io, (iElementString > 0 .and. iElementString <= hb0%iNStr), &
         "HB0_ComputeRHS: Illegal string index")
    str => hb0%Strings(iElementString)

    call IO_Assert(io, (iElementVertex > 0 .and. iElementVertex <= str%iNPts), &
         "HB0_ComputeRHS: Illegal vertex index")
    vtx => str%Vertices(iElementVertex)

    if (lDirect) then
      rRHS = rZERO
    else
      if (iElementFlag == kHB0_Vertex) then
        rRHS = -vtx%rCheck(1)
      else if (iElementFlag == kHB0_Center) then
        rRHS = -vtx%rCheck(2)
      end if
    end if
    return
  end function rHB0_ComputeRHS


  subroutine HB0_ComputeCheck(io, hb0,iElementString, iElementVertex, iElementFlag, cCPZ, rFlow)
    !! function rHB0_ComputeCheck
    !!
    !! Returns the check value for the specified domain and vertex
    !!
    !! Calling Sequence:
    !!    rCheck = rHB0_ComputeCheck(io, hb0,iElementString, iElementVertex, rFlow)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    complex, dimension(:) :: cCPZ
    !!             Control point(s) to be used in coefficient calculations
    !!   (in)    integer :: iElementString
    !!             The domain number to be examined
    !!   (in)    integer :: iElementVertex
    !!             The vertex number to be examined
    !!   (in)    integer :: iElementFlag
    !!             The equation flag for the equation
    !!   (in)    real :: rFlow
    !!             The modeled integrated flow for the vertex
    !!   (in)    type(IO_status), pointer :: io
    !!             pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    complex(kind=AE_REAL), dimension(:), intent(in) :: cCPZ
    real(kind=AE_REAL), intent(in) :: rFlow
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx
    integer(kind=AE_INT) :: idp

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_ComputeCheck: HB0_Create has not been called")
      call IO_Assert(io, (iElementString > 0 .and. iElementString <= hb0%iNStr), &
           "HB0_ComputeCheck: Illegal string number")
    end if
    str => hb0%Strings(iElementString)

    vtx => str%Vertices(iElementVertex)
    if (iElementFlag == kHB0_Vertex) then
      vtx%rCheck(1) = rFlow
    else if (iElementFlag == kHB0_Center) then
      vtx%rCheck(2) = rFlow
    end if

    return
  end subroutine HB0_ComputeCheck


  subroutine HB0_StoreResult(io, hb0,rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine HB0_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the HB0 module.
    !!
    !! Calling Sequence:
    !!    HB0_StoreResult(io, hb0,cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(always ELEM_HB0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!             For HB0, the constants kHB0_Vertex and kHBO_Center are used.
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_StoreResult: HB0_Create has not been called")
      call IO_Assert(io, (iElementString >= 1 .and. iElementString <= hb0%iNStr), &
           "HB0_StoreResult: Bad element string ID")
    end if

    str => hb0%Strings(iElementString)

    if (io%lDebug) then
      call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= str%iNPts), &
           "HB0_StoreResult: Bad element vertex ID")
    end if

    ! All is well.  Store the result...
    vtx => str%Vertices(iElementVertex)
    select case (iElementFlag)
      case (kHB0_Vertex)
        if (lDirect) then
          vtx%rStrength(1) = rValue
        else
          vtx%rStrength(1) = vtx%rStrength(1) + rValue
        end if
      case (kHB0_Center)
        if (lDirect) then
          vtx%rStrength(2) = rValue
        else
          vtx%rStrength(2) = vtx%rStrength(2) + rValue
        end if
    end select

    return
  end subroutine HB0_StoreResult


  subroutine HB0_Update(io, hb0, fdp)
    !! subroutine HB0_StoreResult
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    HB0_Update(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, irv
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Update: HB0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "HB0_Update: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, hb0%iNStr
      str => hb0%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        cRho1 = cmplx(rZERO, this_vtx%rStrength(1), AE_REAL)
        cRho2 = cmplx(rZERO, this_vtx%rStrength(2), AE_REAL)
        if (iVtx == str%iNPts) then
          cRho3 = cZERO
        else
          cRho3 = cmplx(rZERO, next_vtx%rStrength(1), AE_REAL)
        end if
        call FDP_Update(io, fdp, this_vtx%iFDPIndex, (/cRho1, cRho2, cRho3/))
      end do
    end do
  end subroutine HB0_Update


  subroutine HB0_ResetIterator(io, hb0)
    !! subroutine HB0_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call HB0_ResetIterator(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: hb0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_ResetIterator: HB0_Create has not been called")
    end if

    hb0%iIterStr = 1
    hb0%iIterVtx = 0
    hb0%iIterFlag = kHB0_Center

    return
  end subroutine HB0_ResetIterator


  function HB0_NextIterator(io, hb0) result(itr)
    !! function HB0_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call HB0_NextIterator(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: hb0
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_NextIterator: HB0_Create has not been called")
    end if

    if (hb0%iIterStr > hb0%iNStr) then
      nullify(itr)
      return
    end if

    ! Do I increment the vertex?
    if (hb0%iIterFlag == kHB0_Center) then
      hb0%iIterFlag = kHB0_Vertex
      hb0%iIterVtx = hb0%iIterVtx + 1
      ! Is that the end of the string?
      if (hb0%iIterVtx > hb0%Strings(hb0%iIterStr)%iNPts) then
        hb0%iIterStr = hb0%iIterStr+1
        hb0%iIterVtx = 1
        if (hb0%iIterStr > hb0%iNStr) then
          nullify(itr)
          return
        end if
      end if
      ! Handle the case of the end of the  string...
    else if (hb0%iIterVtx == hb0%Strings(hb0%iIterStr)%iNPts) then
      hb0%iIterStr = hb0%iIterStr+1
      hb0%iIterVtx = 1
      if (hb0%iIterStr > hb0%iNStr) then
        nullify(itr)
        return
      end if
      ! Otherwise, we just change from vertex to center...
    else
      hb0%iIterFlag = kHB0_Center
    end if

    str => hb0%Strings(hb0%iIterStr)
    vtx => str%Vertices(hb0%iIterVtx)
    allocate(itr)
    itr%iElementType = ELEM_HB0
    itr%iElementString = hb0%iIterStr
    itr%iElementVertex = hb0%iIterVtx
    itr%iElementFlag = hb0%iIterFlag
    itr%iValueSelector = VALUE_FLOW
    if (hb0%iIterFlag == kHB0_Vertex) then
      allocate(itr%cZ(size(vtx%cVertexCPZ)), stat = iStat)
      call IO_Assert(io, (iStat == 0), "HB0_NextIterator: Space exhausted")
      itr%cZ = vtx%cVertexCPZ
    else
      allocate(itr%cZ(size(vtx%cCenterCPZ)), stat = iStat)
      call IO_Assert(io, (iStat == 0), "HB0_NextIterator: Space exhausted")
      itr%cZ = vtx%cCenterCPZ
    end if

    return
  end function HB0_NextIterator


  subroutine HB0_SetIterator(io, hb0,itr, cValue)
    !! function HB0_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call HB0_SetIterator(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION to be used
    !!   (in0    type(ITERATOR_RESULT), pointer :: itr
    !!             Pointer to the information for data retrieval
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: hb0
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(HB0_COLLECTION), pointer :: hb0
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_NextIterator: HB0_Create has not been called")
      call IO_Assert(io, (hb0%iIterStr <= hb0%iNStr), &
           "HB0_SetIterator: Iterator out of range")
    end if

    hb0%Strings(itr%iElementString)%Vertices(itr%iElementVertex)%rCheck(itr%iElementFlag) = &
                                                                                            real(cValue, AE_REAL)

    return
  end subroutine HB0_SetIterator


  subroutine HB0_Read(io, hb0)
    !! subroutine HB0_Read
    !!
    !! Populates an HB0_COLLECTION with input from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call HB0_Read(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION) :: hb0
    !!             HB0_COLLECTION to be populated
    !!
    !! The format of the HB0 section of the input file appears as follows:
    !!
    !! HB0
    !! STR NVertices ID
    !! x y
    !! ... Up to NVertices
    !! STR NVertices ID
    !! x y s
    !! ... Up to NVertices
    !! ... Up to NStrings
    !!
    !! NOTE: It is assumed that the HB0 line was found by the caller
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io

    ! Locals -- for Directive parsing
    type(DIRECTIVE), dimension(2), parameter :: dirDirectives = (/dirEND, dirSTR/)
    ! Locals -- Input values
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxStr, iMaxVtx
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iStr
    logical :: lFlag
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    call IO_MessageText(io, "  Reading HB0 module input")

    call IO_Assert(io, (associated(hb0)), "HB0_Read: HB0_Create has not been called")

    ! Use ifIO_InputRecord to process the model input file.
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "HB0_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModHB0 "ifXXXRead" routines.
          call IO_Assert(io, .false., "HB0_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(str)), "HB0_Read: No current string")
          call IO_Assert(io, (str%iNPts < size(str%Vertices)), "HB0_Read: Space exhausted")
          cZ = cIO_GetCoordinate(io, "cZ", check_points=str%Vertices%cZ)
          call IO_Assert(io, (iStat == 0), "HB0_Read: I/O Error")
          str%iNPts = str%iNPts+1
          vtx => str%Vertices(str%iNPts)
          vtx%cZ = cZ
          vtx%iFDPIndex = -1
          vtx%rStrength(1) = rZERO
          vtx%rStrength(2) = rZERO
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpSTR)
          !****************************************************************************
          ! Here for the STR command -- create a new barrier string
          ! the maximum number of vertices is in the input record
          !****************************************************************************
          call IO_Assert(io, (associated(hb0%Strings)), "HB0_Read: No strings allocated")
          call IO_Assert(io, (hb0%iNStr < size(hb0%Strings)), "HB0_Read: Space exhausted")
          iMaxVtx = iIO_GetInteger(io, "iMaxVtx", minimum=3)
          iID = iIO_GetInteger(io, "iID", forbidden=hb0%Strings(:)%iID )
          hb0%iNStr = hb0%iNStr+1
          allocate(hb0%Strings(hb0%iNStr)%Vertices(iMaxVtx), stat = iStat)
          call IO_Assert(io, (iStat == 0), "HB0_Read: Allocation failed")
          ! Made it!
          hb0%Strings(hb0%iNStr)%iID = iID
          hb0%Strings(hb0%iNStr)%iNPts = 0
          str => hb0%Strings(hb0%iNStr)
      end select
    end do

    call IO_MessageText(io, "  Leaving HB0 module")

    return
  end subroutine HB0_Read


  subroutine HB0_Inquiry(io, hb0, iLU)
    !! subroutine HB0_Inquiry
    !!
    !! Writes an inquiry report for all barriers to iLU
    !!
    !! Calling Sequence:
    !!    call HB0_Inquiry(io, hb0,iLU)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rLength
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: this, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Inquiry: HB0_Create has not been called")
    end if


    write (unit=iLU, &
           fmt="(""#HB0, VTX, X1, Y1, X2, Y2, LENGTH, STRENGTH_1, STRENGTH_2, FLUX_1, FLUX_2"")")
    do iStr = 1, hb0%iNStr
      str => hb0%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        rLength = abs(str%Vertices(iVtx+1)%cZ - this%cZ)
        write (unit=iLU, &
               fmt="(""HB0"", 2("", "", i9), 9("", "", e16.8))") &
               str%iID, &
               iVtx, &
               cIO_WorldCoords(io, this%cZ), &
               cIO_WorldCoords(io, next%cZ), &
               rLength, &
               this%rStrength(1), &
               this%rStrength(2), &
               this%rCheck(1), &
               this%rCheck(2)
      end do
    end do

    return
  end subroutine HB0_Inquiry


  subroutine HB0_Report(io, hb0)
    !! subroutine HB0_Report
    !!
    !! Writes a debugging report for all no-flows to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call HB0_Report(hb0)
    !!
    !! Arguments:
    !!   (in)    type(HB0_COLLECTION), pointer :: hb0
    !!             HB0_COLLECTION object to be used
    !!
    ! [ io%lDebug ]
    type(HB0_COLLECTION), pointer :: hb0
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(hb0)), &
           "HB0_Report: HB0_Create has not been called")
    end if

    call HTML_Header('Module HB0', 1)
    call HTML_Header('No-flow boundary information', 2)

    if (associated(hb0%Strings)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of strings', hb0%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iHB0_GetInfo(io, hb0,SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iHB0_GetInfo(io, hb0,SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iHB0_GetInfo(io, hb0,SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iHB0_GetInfo(io, hb0,SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iHB0_GetInfo(io, hb0,SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iStr = 1, hb0%iNStr
        str => hb0%Strings(iStr)
        call HTML_Header('Line-sink string definition', 3)
        call HTML_StartTable()
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'FDP # ', 'X     ', 'Y     ', 'V Str ', 'C Str ', &
             'V Flux', 'C Flux'/))
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, vtx%iFDPIndex/))
          call HTML_ColumnComplex((/vtx%cZ/))
          call HTML_ColumnReal((/vtx%rStrength(1), vtx%rStrength(2), vtx%rCheck(1), &
               vtx%rCheck(2)/))
          call HTML_EndRow()
        end do
        vtx => str%Vertices(str%iNPts)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iVtx/))
        call HTML_ColumnText((/'--'/))
        call HTML_ColumnComplex((/vtx%cZ/))
        call HTML_ColumnText((/'--', '--', '--', '--'/))
        call HTML_EndRow()
        call HTML_EndTable()
      end do
    else
      call HTML_Header('No boundaries defined', 3)
    end if

    return
  end subroutine HB0_Report


  subroutine HB0_Save(io, hb0, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(HB0_COLLECTION), pointer :: hb0
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istr, ivtx, iflg
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx

    ! Output records will be of the form ELEM_HB0, IWEL, IRAD, IVTX, 0, SIGMA
    do istr = 1, hb0%iNStr
      str => hb0%Strings(istr)
      do ivtx = 1, str%iNPts
        vtx => str%Vertices(ivtx)
        do iflg = 1, ubound(vtx%rStrength,1)
          if (mode == IO_MODE_BINARY) then
            write (unit=LU_SCRATCH) ELEM_HB0, istr, ivtx, iflg, vtx%rStrength(iflg)
          else
            write (unit=LU_SCRATCH, fmt=*) "HB0", istr, ivtx, iflg, vtx%rStrength(iflg)
          end if
        end do
      end do
    end do

    return
  end subroutine HB0_Save


  subroutine HB0_Load(io, hb0, fdp, mode)
    !! Loads the HB0 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(HB0_COLLECTION), pointer :: hb0
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, istr, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(HB0_STRING), pointer :: str
    type(HB0_VERTEX), pointer :: vtx
    type(HB0_VERTEX), pointer :: flg

    ! Scans the entire precondition file for the HB0 data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, istr, ivtx, iflg, rstrength
        if (imodule /= ELEM_HB0) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, istr, ivtx, iflg, rstrength
        if (uppercase (trim(smodule)) /= "HB0") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      call IO_Assert(io, istr > 0 .and. istr <= hb0%iNStr, "HB0 string not found")
      str => hb0%Strings(istr)
      call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "HB0 vertex not found")
      vtx => str%Vertices(ivtx)
      call IO_Assert(io, iflg > 0 .and. iflg <= ubound(vtx%rStrength,1), "HB0 strength index not found")
      vtx%rStrength(iflg) = rstrength
    end do

    ! Now, populate the internal data structures
    call HB0_Update(io, hb0, fdp)

    return
  end subroutine HB0_Load

end module m_hb0
