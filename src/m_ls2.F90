module m_ls2

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

  !! module m_ls2
  !!
  !! Element module for 2-D head specified line-sinks with resistance
  !! in three modes: GHB, RIVER, DRAIN
  !!
  !!   GHB -- The boundary condition is always enforced
  !!   RIVER -- Flux drops to a fixed value if head is below the resistance layer
  !!   DRAIN -- The linesink can never supply water to the aquifer
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

  type :: LS2_VERTEX
    !! type LS2_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rDepth
    !!     The elevation of the bottom of the resistance layer
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
    ! INPUT DATA
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rHead
    real(kind=AE_REAL) :: rDepth
    ! COMPUTED BY THE MODULE
    real(kind=AE_REAL) :: rCPHead
    real(kind=AE_REAL) :: rDPStrength
    real(kind=AE_REAL) :: rStrength
    real(kind=AE_REAL) :: rLength
    integer(kind=AE_INT) :: iFDPIndex
    complex(kind=AE_REAL), dimension(1) :: cCPZ
    real(kind=AE_REAL) :: rCPDepth
    logical :: lEnabled
    logical :: lRouteEnabled
    real(kind=AE_REAL) :: rCheckPot     ! Check values for RHS calculation
    real(kind=AE_REAL) :: rCheckHead
    real(kind=AE_REAL) :: rSolutionPot      ! Check values for RHS calculation
    real(kind=AE_REAL) :: rSolutionHead     ! Check values for linearization
    ! For routing
    real(kind=AE_REAL) :: rBaseFlow
    real(kind=AE_REAL) :: rSegmentOverlandFlow
    real(kind=AE_REAL) :: rOverlandFlow
    real(kind=AE_REAL) :: rStreamFlow
  end type LS2_VERTEX

  type :: LS2_STRING
    !! type LS2_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   type(LS2_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of LS2_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !!   real :: rCW
    !!     The product of the entry resistance for the linesink and its width
    !!   integer :: iMode
    !!     Flag:
    !!            LS2_MODE_GHB   -- Treat as 'general head boundary'
    !!            LS2_MODE_RIVER -- Treat as 'river' boundary
    !!            LS2_MODE_DRAIN -- Treat as 'drain' boundary
    !!     if .false. treat as a stream segment
    !!   real :: rConductance
    !!     Element conductance, defined as w * k_c / t_c where w is the "width"
    !!     of the element, k_c is the hydraulic conductivity of the resistance
    !!     layer and t_c is the thickness of the resistance layer
    !!
    ! INPUT DATA
    type(LS2_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iFWLIndex
    integer(kind=AE_INT) :: iID
    real(kind=AE_REAL) :: rConductance
    integer(kind=AE_INT) :: iMode
    ! For routing
    logical :: lRoute
    integer(kind=AE_INT) :: iDownstreamID
    integer(kind=AE_INT) :: iRouteMode
    real(kind=AE_REAL) :: rInFlow
    real(kind=AE_REAL) :: rOverlandFlow
    real(kind=AE_REAL) :: rGage
    type(LS2_STRING), pointer, dimension(:) :: UpstreamStrings
  end type LS2_STRING

  type :: LS2_COLLECTION
    !! type LS2_COLLECTION
    !!
    !! Type that holds information for all LS2 elements in a layer
    !!
    !! Members:
    !!   type(LS2_STRING), dimension(:), pointer :: Strings
    !!     A vector of LS2_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !!
    type(LS2_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iRegenerate
    ! Iterator Information
    integer(kind=AE_INT) :: iIterStr
    integer(kind=AE_INT) :: iIterVtx
  end type LS2_COLLECTION

  ! Module flags for matrix generator routines
  integer(kind=AE_INT), private, parameter :: LS2Vertex = 1

  ! Flags for routing
  integer(kind=AE_INT), private, parameter :: ROUTE_END = -1
  integer(kind=AE_INT), private, parameter :: ROUTE_MODE_NORMAL = 0
  integer(kind=AE_INT), private, parameter :: ROUTE_MODE_BUDGET = 1


contains


  function LS2_Create(io) result(ls2)
    !! function LS2_Create
    !!
    !! Creates a new LS2_COLLECTION object
    !!
    !! Calling Sequence:
    !!    ls2 => LS2_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(ls2, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS2_Create: allocation failed")
    nullify(ls2%Strings)
    ls2%iNStr = 0
    ls2%iRegenerate = 0

    return
  end function LS2_Create


  subroutine LS2_Alloc(io, ls2)
    !! Subroutine LS2_Alloc
    !!
    !! Allocates Strings for the LS2_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call LS2_Alloc(io, ls2, iNStr)
    !!
    !! Arguments:
    !!    (in)    type(LS2_COLLECTION), pointer :: ls2
    !!              The LS2_COLLECTION object to be used
    !!    (in)    integer :: iNStr
    !!              The number of strings to make space for
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iStat

    iNStr = iIO_GetInteger(io, 'iNStr', minimum = 0)
    allocate(ls2%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS2_Alloc: allocation failed")
    ls2%Strings(:)%iID = -1

    return
  end subroutine LS2_Alloc


  subroutine LS2_Destroy(io, ls2)
    !! subroutine LS2_Destroy
    !!
    !! Frees memory allocated for ls2 Linesinks and strings of vertices
    !! and the ls2 Collection object
    !!
    !! Calling Sequence:
    !!     call LS2_Destroy(ls2)
    !!
    !! Arguments:
    !!  type(LS2_COLLECTION), pointer :: ls2
    !!              Pointer to the LS2_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: i
    type(LS2_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Destroy: LS2_Create has not been called")
    end if

    ! First deallocate each string of vertices
    do i = 1, ls2%iNStr
      str => ls2%Strings(i)
      deallocate(str%Vertices, stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS2_Destroy: deallocation of Vertices failed")
    end do
    ! Then deallocate the strings table
    if (associated(ls2%Strings)) then
      deallocate(ls2%Strings, stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS2_Destroy: deallocation of Strings failed")
    end if
    ! Now deallocate the collection
    deallocate(ls2, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS2_Destroy: deallocation failed")

    return
  end subroutine LS2_Destroy


  subroutine LS2_New(io, ls2, Vertices, iNPts)
    !! function LS2_New
    !!
    !! Adds a new LS2_STRING object to the LS2_COLLECTION 'ls2'
    !!
    !! Calling Sequence:
    !!    call LS2_New(io, ls2, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(LS2_COLLECTION), pointer :: ls2
    !!              The LS2_COLLECTION object to be used
    !!    (in)    type(LS2_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(LS2_VERTEX), dimension(:) :: Vertices
    integer(kind=AE_INT), intent(in) :: iNPts
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(LS2_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_New: LS2_Create has not been called")
    end if

    call IO_Assert(io, (ls2%iNStr < size(ls2%Strings)), &
         "LS2_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(Vertices)), &
         "LS2_New: Size of provided vertices is inconsistent")

    ls2%iNStr = ls2%iNStr + 1
    str => ls2%Strings(ls2%iNStr)
    allocate(str%Vertices(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS2_New: Allocation failed")
    str%Vertices = Vertices(1:iNPts)
    str%iNPts = iNPts
    str%lRoute = .false.
    str%iDownstreamID = ROUTE_END
    str%iRouteMode = ROUTE_MODE_NORMAL
    str%rInFlow = rZERO
    str%rOverlandFlow = rZERO
    str%rGage = rZERO

    return
  end subroutine LS2_New


  function iLS2_GetID(io, ls2, iIndex) result(iID)
    !! Returns the ID number for the well at index 'iIndex'
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= ls2%iNStr), "Internal error -- no such index")
    iID = ls2%Strings(iIndex)%iID

    return
  end function iLS2_GetID


  subroutine LS2_PreSolve(io, ls2)
    !! subroutine LS2_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call LS2_PreSolve(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rTotalLen
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx

    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)
      rTotalLen = rZERO
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        this_vtx%rLength = abs(next_vtx%cZ-this_vtx%cZ)
        this_vtx%rStrength = rZERO
        rTotalLen = rTotalLen + this_vtx%rLength
      end do
      ! Assign the overland flow
      if (str%lRoute) then
        do iVtx = 1, str%iNPts-1
          this_vtx => str%Vertices(iVtx)
          this_vtx%rSegmentOverlandFlow = str%rOverlandFlow * this_vtx%rLength / rTotalLen
        end do
      end if
    end do

    return
  end subroutine LS2_PreSolve


  function iLS2_GetInfo(io, ls2, iOption, iIteration) result(iValue)
    !! function LS2_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iLS2_GetInfo(io, ls2, iOption)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr, iVtx
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_GetInfo: LS2_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = ls2%iNStr
      case (SIZE_FDP)
        do iStr = 1, ls2%iNStr
          str => ls2%Strings(iStr)
          iValue = iValue + str%iNPts-1
        end do
      case (SIZE_EQUATIONS)
        do iStr = 1, ls2%iNStr
          str => ls2%Strings(iStr)
          do iVtx = 1, str%iNPts-1
            vtx => str%Vertices(iVtx)
            if (vtx%lEnabled) iValue = iValue + 1
          end do
        end do
      case (SIZE_UNKNOWNS)
        do iStr = 1, ls2%iNStr
          str => ls2%Strings(iStr)
          do iVtx = 1, str%iNPts-1
            vtx => str%Vertices(iVtx)
            if (vtx%lEnabled) iValue = iValue + 1
          end do
        end do
      case (INFO_REGENERATE)
        !print *,'INFO_REGEN',ls2%iRegenerate
        if (iIteration < 2) then
          iValue = 1
        else
          iValue = ls2%iRegenerate
        end if
      case default
        iValue = 0
    end select

    return
  end function iLS2_GetInfo


  subroutine LS2_SetupFunctions(io, ls2, fwl, fdp)
    !! subroutine LS2_SetupFunctions
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call LS2_Setup(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iStr2, iVtx, i, iDP, iWL, iUpstreamCount, iStat
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(LS2_STRING), pointer :: str, str2
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Setup: LS2_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "LS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "LS2_Setup: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)

      ! Build dipoles for all segments
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        this_vtx%rLength = abs(next_vtx%cZ - this_vtx%cZ)
        call FDP_New(io, fdp, this_vtx%cZ, next_vtx%cZ, (/cZERO, cZERO, cZERO/), ELEM_LS2, iStr, iVtx, -1, this_vtx%iFDPIndex)
      end do

      ! Put a well at the end of the string
      last_vtx => str%Vertices(str%iNPts)
      call FWL_New(io, fwl, last_vtx%cZ, rZERO, rZERO, ELEM_LS2, iStr, -1, -1, iWL)
      str%iFWLIndex = iWL
    end do

    return
  end subroutine LS2_SetupFunctions


  subroutine LS2_SetupMatrix(io, ls2, aqu, mat)
    !! subroutine LS2_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the module
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call LS2_Setup(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, i, iDP, iWL, iEQ
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead, rDepth1, rDepth2
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Setup: LS2_Create has not been called")
      call IO_Assert(io, (associated(aqu)), &
           "LS2_Setup: Illegal AQU_COLLECTION object")
      call IO_Assert(io, (associated(mat)), &
           "LS2_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix generator entries for all segments
    do iStr = 1, ls2%iNStr
      ! Set up the unknown variables
      ! Vertex entry -- all vertices
      str => ls2%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        if (str%Vertices(iVtx)%lEnabled) then
          call MAT_CreateVariable(io, mat, ELEM_LS2, iStr, iVtx, LS2Vertex)
        end if
      end do

      ! Set up control points and equations -- One equation per segment
      do iVtx = 1, str%iNPts-1
        if (str%Vertices(iVtx)%lEnabled) then
          this_vtx => str%Vertices(iVtx)
          next_vtx => str%Vertices(iVtx+1)
          rHead1 = this_vtx%rHead
          rHead2 = next_vtx%rHead
          rDepth1 = this_vtx%rDepth
          rDepth2 = next_vtx%rDepth
          ! Compute control points for the segment.  There is one unknown
          ! per segment.
          call MAT_ComputeControlPoints(io, this_vtx%cZ, next_vtx%cZ, 1, cCPResult, rZERO)
          ! Now, create the equation entry...
          this_vtx%cCPZ(1) = cCPResult(2)
          this_vtx%rCPHead = rHead1 + (rHead2-rHead1)*(cCPResult(2)-this_vtx%cZ)/(next_vtx%cZ-this_vtx%cZ)
          this_vtx%rCPDepth = rDepth1 + (rDepth2-rDepth1)*(cCPResult(2)-this_vtx%cZ)/(next_vtx%cZ-this_vtx%cZ)
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult(2)/), EQN_HEAD, ELEM_LS2, iStr, iVtx, 0, cZERO, rZERO)
        end if
      end do
    end do

    return
  end subroutine LS2_SetupMatrix


  function iLS2_Prepare(io, ls2, aqu, iIteration) result(iChanges)
    !! subroutine LS2_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! For LS2, this routine handles the turning on and off of elements due to
    !! percolating rivers, negative drains, and eventually streamflow routing.
    !!
    !! Calling Sequence:
    !!    call LS2_Setup(wl1, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer :: io
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx, this_vtx, next_vtx
    character(len=132) :: sBuf

    ! Clear the "regenerate" flag unless something interesting happens.
    ls2%iRegenerate = 0

    iChanges = 0
    ! Don't try to adjust the elements on the first iteration
    if (iIteration > 1) then
      do iStr = 1, ls2%iNStr
        str => ls2%Strings(iStr)
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          select case (str%iMode)
            case (LS2_MODE_GHB)
              ! Placeholder -- does nothing
              continue
            case (LS2_MODE_RIVER)
              ! Handle percolating rivers
              ! If the head is below the bottom of the resistance layer, turn it
              ! off, and set the strength for percolating conditions...
              if (vtx%rCheckHead < vtx%rCPDepth) then
                ! Turn off the vertex only
                if (vtx%lEnabled) then
                  write (unit=sBuf, &
                         fmt="(""  Disabling river element string "", i8, "" segment "", i8, 1x, f10.2, 1x, f10.2)" &
                         ) str%iID, iVtx, vtx%rCheckHead, vtx%rCPDepth
                  call IO_MessageText(io, sBuf)
                  iChanges = iChanges+1
                  vtx%rStrength = -(vtx%rCPHead - vtx%rCPDepth) * str%rConductance
                  vtx%lEnabled = .false.
                  ls2%iRegenerate = 1
                end if
                ! ...if not check to see if it needs to be turned back on.
              else
                if (.not. vtx%lEnabled .and. vtx%lRouteEnabled) then
                  write (unit=sBuf, &
                         fmt="(""  Enabling river element string  "", i8, "" segment "", i8, 1x, f10.2, 1x, f10.2)" &
                         ) str%iID, iVtx, vtx%rCheckHead, vtx%rCPDepth
                  call IO_MessageText(io, sBuf)
                  iChanges = iChanges+1
                  if (vtx%lRouteEnabled) vtx%lEnabled = .true.
                  ls2%iRegenerate = 1
                end if
              end if
            case (LS2_MODE_DRAIN)
              ! Handle losing drains
              if (vtx%lEnabled) then
                ! If the drain was on for the last iteration, check it...
                if (vtx%rStrength <= rZERO) then
                  ! Turn off the vertex
                  write (unit=sBuf, &
                         fmt="(""  Disabling drain element string "", i8, "" segment "", i8)" &
                         ) str%iID, iVtx
                  call IO_MessageText(io, sBuf)
                  vtx%rStrength = rZERO
                  vtx%lEnabled = .false.
                  ls2%iRegenerate = 1
                  iChanges = iChanges+1
                end if
              else
                ! ...if not check to see if it needs to be turned back on.
                if (rAQU_PotentialToHead(io, aqu, vtx%rCheckPot, vtx%cCPZ(1)) > vtx%rCPHead) then
                  write (unit=sBuf, &
                         fmt="(""  Enabling drain element string  "", i8, "" segment "", i8)" &
                         ) str%iID, iVtx
                  call IO_MessageText(io, sBuf)
                  vtx%lEnabled = .true.
                  ls2%iRegenerate = 1
                  iChanges = iChanges+1
                end if
              end if
          end select

          ! First, force a matrix regen if any of the vertices is unconfined(darn shame)
          if (vtx%lEnabled .and. &
              lAQU_IsConfined(io, aqu, vtx%cCPZ(1), vtx%rCheckPot) &
            ) ls2%iRegenerate = 1
        end do

        ! Update the dipole strengths also
        do iVtx = 1, str%iNPts-1
          this_vtx => str%Vertices(iVtx)
          next_vtx => str%Vertices(iVtx+1)
          next_vtx%rDPStrength = this_vtx%rDPStrength + this_vtx%rStrength * this_vtx%rLength
        end do
      end do
    end if

    print *, 'Made ', iChanges, ' changes to resistance line-sinks'

    return
  end function iLS2_Prepare


  function rLS2_GetCoefficientMultiplier(io, ls2, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    rMultiplier = rONE

    return
  end function rLS2_GetCoefficientMultiplier


  subroutine LS2_ComputeCoefficients(io, ls2, aqu, fwl, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine LS2_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the LS2
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call LS2_ComputeCoefficients(io, ls2, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
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
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
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
    integer(kind=AE_INT) :: iStat, iCol, iStr, iVtx, iDP1, iNDP, iWhich
    real(kind=AE_REAL) :: rCorr, rH, rPot
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDPW
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx, first_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_ComputeCoefficients: LS2_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "LS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "LS2_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iCol = 0
    rARow = rZERO
    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)
      first_vtx => str%Vertices(1)
      ! ASSUMES that LS2_Setup routine created consecutive dipole entries
      iDP1 = first_vtx%iFDPIndex
      iNDP = str%iNPts-1
      allocate(cDPF(0:iNDP, 1, 1), cDPW(0:iNDP, 1, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS2_ComputeCoefficients: Allocation failed")

      ! Get the appropriate influence functions for the boundary condition type
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
        vtx => str%Vertices(iVtx)
        if (vtx%lEnabled) then
          iCol = iCol+1
          rARow(iCol) = real(cDPF(iVtx, 1, 1))
          if (iElementType == ELEM_LS2 .and. &
              iElementString == iStr .and. &
              iElementVertex == iVtx) then
            ! Need to add the resistance correction term

            ! On the first iteration, use the CP head as a head estimate...
            if (iIteration == 1) then
              vtx%rSolutionPot = rAQU_HeadToPotential(io, aqu, vtx%rCPHead, vtx%cCPZ(1))
              vtx%rSolutionHead = vtx%rCPHead
            end if

            if (lAQU_IsConfined(io, aqu, vtx%cCPZ(1), vtx%rSolutionPot)) then
              rCorr = -rAQU_Transmissivity(io, aqu, vtx%cCPZ(1), vtx%rSolutionPot) / str%rConductance
            else
              if (iIteration == 1) then
                rCorr = - rAQU_HydCond(io, aqu, vtx%cCPZ(1)) * &
                        (vtx%rCPHead - rAQU_Base(io, aqu, vtx%cCPZ(1))) /  &
                        str%rConductance
              else
                rCorr = -rHALF * rAQU_HydCond(io, aqu, vtx%cCPZ(1)) * &
                        (vtx%rSolutionHead + vtx%rCPHead - rTWO*rAQU_Base(io, aqu, vtx%cCPZ(1))) /  &
                        str%rConductance
              end if
            end if
            rARow(iCol) = rARow(iCol) + rCorr
          end if
        end if
      end do

      deallocate(cDPF, cDPW)
    end do

    rARow = rARow * rMultiplier

    return
  end subroutine LS2_ComputeCoefficients


  function rLS2_ComputeRHS(io, ls2, aqu, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rLS2_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rLS2_ComputeRHS(io, ls2, aqu, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION object to be used
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
    !!
    !! Return Value:
    !!   real :: rRHS
    !!     The RHS value for the module
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
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
    real(kind=AE_REAL) :: rCorr, rPot
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx

    call IO_Assert(io, (iElementString > 0 .and. iElementString <= ls2%iNStr), &
         'Illegal string index')
    str => ls2%Strings(iElementString)
    call IO_Assert(io, (iElementVertex > 0 .and. iElementVertex <= str%iNPts), &
         'Illegal vertex index')
    vtx => str%Vertices(iElementVertex)

    ! For LS2, compute the RHS by subtracting the previous result from the
    ! desired potential at the control-point
    rRHS = rAQU_HeadToPotential(io, aqu, vtx%rCPHead, vtx%cCPZ(1)) - vtx%rCheckPot

    ! Need to add the resistance correction term
    if (lDirect) then
      rCorr = rZERO
    else
      if (lAQU_IsConfined(io, aqu, vtx%cCPZ(1), vtx%rSolutionPot)) then
        rCorr = vtx%rStrength * rAQU_Transmissivity(io, aqu, vtx%cCPZ(1), vtx%rSolutionPot) / str%rConductance
      else
        rCorr = vtx%rStrength * rHALF * rAQU_HydCond(io, aqu, vtx%cCPZ(1)) * &
                (vtx%rSolutionHead + vtx%rCPHead - rTWO*rAQU_Base(io, aqu, vtx%cCPZ(1))) / str%rConductance
      end if
    end if
    rRHS = rRHS + rCorr

    return
  end function rLS2_ComputeRHS


  subroutine LS2_StoreResult(io, ls2, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine LS2_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the LS2 module.
    !!
    !! Calling Sequence:
    !!    LS2_StoreResult(io, ls2, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(always ELEM_LS2)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_StoreResult: LS2_Create has not been called")
      call IO_Assert(io, (iElementString >= 1 .and. iElementString <= ls2%iNStr), &
           "LS2_StoreResult: Bad element string ID")
    end if

    str => ls2%Strings(iElementString)
    if (io%lDebug) then
      call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= str%iNPts-1), &
           "LS2_StoreResult: Bad element vertex ID")
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
  end subroutine LS2_StoreResult


  subroutine LS2_Update(io, ls2, fwl, fdp)
    !! subroutine LS2_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    LS2_Update(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Update: LS2_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "LS2_Update: Illegal FDP_COLLECTION object")
    end if

    ! Update dipoles and wells
    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this_vtx => str%Vertices(iVtx)
        next_vtx => str%Vertices(iVtx+1)
        cRho1 = cmplx(this_vtx%rDPStrength, rZERO, AE_REAL)
        cRho3 = cmplx(next_vtx%rDPStrength, rZERO, AE_REAL)
        cRho2 = rHALF * (cRho1 + cRho3)
        call FDP_Update(io, fdp, this_vtx%iFDPIndex, (/cRho1, cRho2, cRho3/))
      end do
      ! Put a well at the end of the string
      call FWL_Update(io, fwl, str%iFWLIndex, real(cRho3))
    end do

    return
  end subroutine LS2_Update


  subroutine LS2_FindStringPointer(io, ls2, iLSID, LSString, lFound)
    !! subroutine LS2_FindStringPointer
    !!
    !! Finds the linesink string specified by the ID and returns a pointer to it
    !!
    !! Calling Sequence:
    !!    call LS2_FindStringPointer(io, ls2, iLSID, LSString, lfound)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (in)    integer :: iLSID
    !!             The linesink string ID number
    !!   (out)   type(LS2_STRING) :: LSString
    !!             Pointer to the linesink string
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iLSID
    type(LS2_STRING), pointer :: LSString
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS0_FindStringPointer: LS2_Create has not been called")
    end if

    lFound = .false.
    do i = 1, ls2%iNStr
      LSString => ls2%Strings(i)
      if (LSString%iID == iLSID) then
        lFound = .true.
        return
      end if
    end do

    return
  end subroutine LS2_FindStringPointer


  subroutine LS2_ResetIterator(io, ls2)
    !! subroutine LS2_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call LS2_ResetIterator(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls2
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_ResetIterator: LS2_Create has not been called")
    end if

    ls2%iIterStr = 1
    ls2%iIterVtx = 0

    return
  end subroutine LS2_ResetIterator


  function LS2_NextIterator(io, ls2) result(itr)
    !! function LS2_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS2_NextIterator(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls2
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_NextIterator: LS2_Create has not been called")
    end if

    if (ls2%iIterStr > ls2%iNStr) then
      nullify(itr)
      return
    end if

    ls2%iIterVtx = ls2%iIterVtx + 1
    if (ls2%iIterVtx > ls2%Strings(ls2%iIterStr)%iNPts-1) then
      ls2%iIterStr = ls2%iIterStr+1
      ls2%iIterVtx = 1
      if (ls2%iIterStr > ls2%iNStr) then
        nullify(itr)
        return
      end if
    end if

    allocate(itr)
    itr%iElementType = ELEM_LS2
    itr%iElementString = ls2%iIterStr
    itr%iElementVertex = ls2%iIterVtx
    itr%iValueSelector = VALUE_POTENTIAL
    allocate(itr%cZ(1))
    itr%cZ(1) = ls2%Strings(ls2%iIterStr)%Vertices(ls2%iIterVtx)%cCPZ(1)

    return
  end function LS2_NextIterator


  subroutine LS2_SetIterator(io, ls2, aqu, itr, cValue, lLinearize)
    !! function LS2_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS2_SetIterator(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be used
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: ls2
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cValue
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io
    type(LS2_VERTEX), pointer :: vtx
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_NextIterator: LS2_Create has not been called")
      call IO_Assert(io, (ls2%iIterStr <= ls2%iNStr), &
           "LS2_SetIterator: Iterator out of range")
    end if

    vtx => ls2%Strings(itr%iElementString)%Vertices(itr%iElementVertex)
    vtx%rCheckPot = real(cValue, AE_REAL)
    vtx%rCheckHead = rAQU_PotentialToHead(io, aqu, vtx%rCheckPot, vtx%cCPZ(1))
    if (lLinearize) then
      vtx%rSolutionPot = vtx%rCheckPot
      vtx%rSolutionHead = vtx%rCheckHead
      !print *, 'cp, chk', ls2%Strings(itr%iElementString)%iID, itr%iElementVertex, vtx%rCPHead, vtx%rCheckHead, vtx%cCPZ(1)
      if (.not. lAQU_IsConfined(io, aqu, vtx%cCPZ(1), vtx%rCheckPot)) ls2%iRegenerate = 1
      !print *,'ls2 update',itr%iElementString,itr%iElementVertex,vtx%rSolutionPot,vtx%rSolutionHead,vtx%rCPHead,vtx%lEnabled
    end if

    return
  end subroutine LS2_SetIterator


  function iLS2_DoRouting(io, ls2, aqu, iCurrentIteration, lSwitch, lReportFlows) result(iChanges)
    !! subroutine LS2_DoRouting
    !!
    !! Performs routing calculations
    !!
    !! Calling Sequence:
    !!    call LS2_DoRouting(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be populated
    !!   (in)    type(IO_STATUS), pointer :: io
    !!             I/O acces information
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iCurrentIteration
    logical, intent(in) :: lSwitch
    logical, intent(in) :: lReportFlows
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges
    ! [ LOCALS ]
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, next_vtx
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rStreamFlow

    ! Route all the strings that have routing enabled, but no downstream routes
    iChanges = 0
    ! Don't do this on the first two iterations
    if (iCurrentIteration > 2) then
      do iStr = 1, ls2%iNStr
        str => ls2%Strings(iStr)
        if (str%lRoute .and. str%iDownstreamID < 0) then
          rStreamFlow = rZERO
          iChanges = iChanges + iLS2_RouteString(io, ls2, aqu, str%iID, rStreamFlow, lSwitch, 0)
          if (lReportFlows) then
            print *, 'Total routing for string ', str%iID, ':'
            print *, '  Stream flow      ', rStreamFlow
          end if
        end if
      end do

      ! Update the dipole strengths also
      do iStr = 1, ls2%iNStr
        str => ls2%Strings(iStr)
        do iVtx = 1, str%iNPts-1
          this_vtx => str%Vertices(iVtx)
          next_vtx => str%Vertices(iVtx+1)
          next_vtx%rDPStrength = this_vtx%rDPStrength + this_vtx%rStrength * this_vtx%rLength
        end do
      end do
    end if
    print *,'Number of routing changes: ',iChanges

    return
  end function iLS2_DoRouting

  recursive function iLS2_RouteString(io, ls2, aqu, iID, rStreamFlow, lSwitch, iChangesIn) result(iChanges)
    !! subroutine LS2_RouteString
    !!
    !! Performs routing calculations for a given string, storing routed values
    !! into the string's vertices. Returns the flow out the bottom of the string.
    !! NOTE: This routine has side-effects and is not parallel-safe.
    !!
    !! Calling Sequence:
    !!    call LS2_RouteString(io, ls2, iID, rBaseFlow, rOverlandFlow, lSwitch)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be populated
    !!   (in)    integer, intent(in) :: iID
    !!             ID number for the string to be routed
    !!   (out)   real :: rBaseFlow
    !!             The base flow out of the string
    !!   (out)   real :: rOverlandFlow
    !!             The overland flow out of the string
    !!   (in)    logical :: lSwitch
    !!             If .false., all line-sink toggling is disabled
    !!   (in)    type(IO_STATUS), pointer :: io
    !!             I/O acces information
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iID
    real(kind=AE_REAL), intent(out) :: rStreamFlow
    logical, intent(in) :: lSwitch
    integer(kind=AE_INT), intent(in) :: iChangesIn
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges
    ! [ LOCALS ]
    type(LS2_STRING), pointer :: this, str
    type(LS2_VERTEX), pointer :: vtx, next_vtx
    integer(kind=AE_INT) :: iVtx, iStr
    logical :: lFound
    real(kind=AE_REAL) :: rSQi

    ! Look up the string to be routed
    call LS2_FindStringPointer(io, ls2, iID, this, lFound)
    call IO_Assert(io, lFound, "LS2_RouteString: No string found.")


    ! Return zeroes if my routing is turned off
    if (.not. this%lRoute) then
      rStreamFlow = rZERO
      return
    end if

    ! First, find the amount of water flowing in from above the string
    rStreamFlow = rStreamFlow + this%rInFlow
    ! Tally up upstream routes
    iChanges = iChangesIn
    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)
      if (str%iID == this%iID) continue
      if (this%iID == str%iDownstreamID) then
        rSQi = rZERO
        iChanges = iLS2_RouteString(io, ls2, aqu, str%iID, rSQi, lSwitch, iChanges)
        !print *,'QB', str%iID,rBQi, rIQi, rOQi
        rStreamFlow = rStreamFlow + rSQi
      end if
    end do

    ! route water down this string...
    vtx => this%Vertices(1)
    vtx%rStreamFlow = rStreamFlow
    do iVtx = 1, this%iNPts-1
      vtx => this%Vertices(iVtx)
      next_vtx => this%Vertices(iVtx+1)
      rStreamFlow = rStreamFlow + vtx%rLength * vtx%rStrength + vtx%rSegmentOverlandFlow
      next_vtx%rStreamFlow = rStreamFlow
      !print *,'accum', str%iID, iVtx, vtx%rStrength, vtx%rSegmentOverlandFlow, ' -- ', rStreamFlow
      if (lSwitch) then
        select case (this%iRouteMode)
          case (ROUTE_MODE_BUDGET)
            continue
          case (ROUTE_MODE_NORMAL)
            if (vtx%rCheckHead > vtx%rCPHead .and. .not. vtx%lRouteEnabled) then
              ! A disabled reach where the aquifer head is higher than stream head was found
              ! Turn it back on
              write (unit=IO_MessageBuffer, fmt=*) "  Reenabled stream reach [head] ID =", this%iID, &
                                                   " Vertex =", iVtx, vtx%rCPHead, &
                                                   vtx%rCheckHead, rStreamFlow
              call IO_MessageText(io)
              vtx%lRouteEnabled = .true.
              vtx%lEnabled = .true.
              ls2%iRegenerate = 1
              iChanges = iChanges + 1
            else if (rStreamFlow > rZERO .and. .not. vtx%lRouteEnabled) then
              write (unit=IO_MessageBuffer, fmt=*) "  Reenabled stream reach [Q>0] ID =", this%iID, " Vertex =", iVtx, rStreamFlow
              call IO_MessageText(io)
              vtx%lRouteEnabled = .true.
              vtx%lEnabled = .true.
              ls2%iRegenerate = 1
              iChanges = iChanges + 1
            else if (rStreamFlow <= rZERO .and. vtx%lRouteEnabled) then
              ! A dry reach is found
              write (unit=IO_MessageBuffer, fmt=*) "  Negative streamflow ID =", this%iID, " Vertex =", iVtx, vtx%rStreamFlow, &
                     vtx%rStrength*vtx%rLength, next_vtx%rStreamFlow
              call IO_MessageText(io)
              next_vtx%rStreamFlow = rZERO
              !if (vtx%rStreamFlow > rZERO) then
              !  vtx%rStrength = vtx%rStreamFlow / vtx%rLength
              !else
              vtx%rStrength = rZERO
              !end if
              vtx%lEnabled = .false.
              vtx%lRouteEnabled = .false.
              ls2%iRegenerate = 1
              iChanges = iChanges + 1
              cycle
            end if
          case default
            call IO_Assert(io, .false., 'Illegal routing model specified')
        end select
      end if
    end do

    return
  end function iLS2_RouteString


  function rLS2_Gage(io, ls2, iID) result(rQ)
    !! subroutine LS2_Gage
    !!
    !! Returns the gage value at the end of the string with ID = iID.
    !!
    !! Calling Sequence:
    !!    q = rLS2_Gage(io, ls2, iID)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be populated
    !!   (in)    integer, intent(in) :: iID
    !!             ID number for the string to be routed
    !!   (in)    type(IO_STATUS), pointer :: io
    !!             I/O acces information
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iID
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    logical :: lFound
    type(LS2_STRING), pointer :: this
    type(LS2_VERTEX), pointer :: vtx

    ! Look up the string to be routed
    call LS2_FindStringPointer(io, ls2, iID, this, lFound)
    call IO_Assert(io, lFound, "LS2_RouteString: No string found.")

    ! The gage value is the total flow at the next-to-last vertex
    vtx => this%Vertices(this%iNPts)
    rQ = vtx%rStreamFlow

    return
  end function rLS2_Gage


  subroutine LS2_Read(io, ls2)
    !! subroutine LS2_Read
    !!
    !! Populates an LS2_COLLECTION using data from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call LS2_Read(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer :: ls2
    !!             LS2_COLLECTION to be populated
    !!
    !! The format of the LS2 section of the input file appears as follows:
    !! LS2
    !! STR NVertices
    !! x y head
    !! ... Up to NVertices
    !! STR NVertices
    !! x y head
    !! ... Up to NVertices
    !! ... Up to NStrings
    !!
    !! NOTE: It is assumed that the LS2 line was found by the caller

    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(3), parameter :: dirDirectives = (/dirEND, dirSTR, dirROU/)
    ! [ LOCALS ]
    character(len=132) :: sTag
    real(kind=AE_REAL) :: rHead
    real(kind=AE_REAL) :: rConductance
    real(kind=AE_REAL) :: rDepth
    real(kind=AE_REAL) :: rInFlow
    real(kind=AE_REAL) :: rOverlandFlow
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxStr, iMaxVtx
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iStr
    integer(kind=AE_INT) :: iMode
    integer(kind=AE_INT) :: iDownstreamID
    logical :: lFlag
    logical :: lDefaultEnabled
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this_vtx, last_vtx

    call IO_MessageText(io, "  Reading LS2 module input")
    call IO_Assert(io, (associated(ls2)), "LS2_Read: LS2_Create has not been called")

    nullify(str, this_vtx, last_vtx)
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpFileEOF)
          ! EOF is unexpected for all ModLS2 "ifXXXRead" routines.
          call IO_Assert(io, .false., "LS2_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(str)), "LS2_Read: No current string")
          call IO_Assert(io, (str%iNPts < size(str%Vertices)), "LS2_Read: Space exhausted")
          write (unit=sTag, fmt=*) 'cZ ', ls2%iNStr, str%iNPts+1
          cZ = cIO_GetCoordinate(io, sTag, extents=.true., check_points=str%Vertices(:)%cZ)
          rHead = rIO_GetReal(io, 'rHead')
          rDepth = rIO_GetReal(io, 'rDepth', maximum = rHead-rTINY)

          str%iNPts = str%iNPts+1
          this_vtx => str%Vertices(str%iNPts)
          this_vtx%cZ = cZ
          this_vtx%rHead = rHead
          this_vtx%rDepth = rDepth
          this_vtx%rDPStrength = rZERO
          this_vtx%lEnabled = lDefaultEnabled
          this_vtx%lRouteEnabled = .true.
          ! Start out percolating?
          if (lDefaultEnabled) then
            this_vtx%rStrength = rZERO
          else
            this_vtx%rStrength = str%rConductance * (this_vtx%rDepth - this_vtx%rHead)
          end if
          if (str%iNPts > 1) then
            last_vtx%rLength = abs(this_vtx%cZ - last_vtx%cZ)
          end if
          this_vtx%iFDPIndex = -1
          this_vtx%rBaseFlow = rZERO
          this_vtx%rOverlandFlow = rZERO
          this_vtx%rStreamFlow = rZERO
          ! Length and overland flow are calculated by LS2_PreSolve
          this_vtx%rLength = rZERO
          this_vtx%rCheckPot = rZERO
          last_vtx => this_vtx
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpSTR)
          !****************************************************************************
          ! Here for the STR command -- create a new string of line-sinks
          ! the maximum number of vertices is in the input record
          !****************************************************************************
          call IO_Assert(io, (ls2%iNStr < size(ls2%Strings)), "LS2_Read: Space exhausted")
          iMaxVtx = iIO_GetInteger(io, 'iMaxVtx', minimum = 2)
          iMode = iIO_GetInteger(io, 'iMode', &
                  allowed=(/LS2_MODE_GHB, LS2_MODE_RIVER, LS2_MODE_DRAIN/))
          rConductance = rIO_GetReal(io, 'rConductance')
          iID = iIO_GetInteger(io, 'iID')
          lDefaultEnabled = lIO_GetLogical(io, 'lDefaultEnabled', def = .true.)
          ! OKAY! Allocate the vertices...
          ls2%iNStr = ls2%iNStr+1
          str => ls2%Strings(ls2%iNStr)
          allocate(str%Vertices(iMaxVtx), stat = iStat)
          call IO_Assert(io, (iStat == 0), "LS2_Read: Allocation failed")
          str%iFWLIndex = -1
          str%iID = iID
          str%iNPts = 0
          str%iMode = iMode
          str%rConductance = rConductance
          str%iRouteMode = ROUTE_MODE_BUDGET
          str%iDownstreamId = -1
          str%Vertices(:)%cZ = cHUGE
        case (kOpROU)
          !****************************************************************************
          ! Here for the ROU command -- Add routing information for the current string
          !****************************************************************************
          call IO_Assert(io, (ls2%iNStr > 0), "LS2_Read: No current string for ROU directive")
          iMode = iIO_GetInteger(io, 'iMode')
          iDownstreamID = iIO_GetInteger(io, 'iDownstreamID')
          rInFlow = rIO_GetReal(io, 'rInFlow')
          rOverlandFlow = rIO_GetReal(io, 'rOverlandFlow')
          str => ls2%Strings(ls2%iNStr)
          str%lRoute = .true.
          str%iDownstreamID = iDownstreamID
          str%iRouteMode = iMode
          str%rInFlow = rInFlow
          str%rOverlandFlow = rOverlandFlow
          str%rGage = rZERO
      end select
    end do

    call IO_MessageText(io, "  Leaving LS2 module")

  end subroutine LS2_Read


  subroutine LS2_Inquiry(io, ls2, iLU)
    !! subroutine LS2_Inquiry
    !!
    !! Writes an inquiry report for all line-sinks to iLU
    !!
    !! Calling Sequence:
    !!    call LS2_Inquiry(io, ls2, iLU)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rLength
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: this, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Inquiry: LS2_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#LS2, ID, VTX, MODE, ENABLED, COND, X1, Y1, X2, Y2, LENGTH, SPEC_HEAD, DEPTH, STRENGTH, MOD_HEAD, " // &
           "FLUX_ERROR, ROUTE?, ROUTE_MODE, DOWNSTREAM_ID, BASEFLOW, OVFLOW, STRFLOW"")")
    do iStr = 1, ls2%iNStr
      str => ls2%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        write (unit=iLU, &
               fmt="(""LS2"", 3("", "", i8), 1("", "", l1), 11("", "", e16.8), "", "", l1, 2("", "", i8), 3("", "", e16.8))" &
               ) str%iID, iVtx, str%iMode, this%lEnabled, &
               str%rConductance, &
               cIO_WorldCoords(io, this%cZ), &
               cIO_WorldCoords(io, next%cZ), &
               this%rLength, &
               this%rCPHead, &
               this%rCPDepth, &
               this%rStrength, &
               this%rCheckHead, &
               this%rStrength-str%rConductance*(this%rCheckHead-this%rCPHead), &
               str%lRoute, &
               str%iRouteMode, &
               str%iDownstreamID, &
               this%rBaseFlow, &
               this%rOverlandFlow, &
               this%rStreamFlow
      end do
    end do

    return
  end subroutine LS2_Inquiry


  subroutine LS2_Report(io, ls2, aqu)
    !! subroutine LS2_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call LS2_Report(ls2)
    !!
    !! Arguments:
    !!   (in)    type(LS2_COLLECTION), pointer
    !!             LS2_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS2_COLLECTION), pointer :: ls2
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    real(kind=AE_REAL) :: rCheckFlux
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls2)), &
           "LS2_Report: LS2_Create has not been called")
    end if

    call HTML_Header('Module LS2', 1)
    call HTML_Header('Resistance line sink information', 2)

    if (associated(ls2%Strings)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of strings', ls2%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iLS2_GetInfo(io, ls2, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iLS2_GetInfo(io, ls2, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iLS2_GetInfo(io, ls2, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iLS2_GetInfo(io, ls2, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iLS2_GetInfo(io, ls2, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iStr = 1, ls2%iNStr
        str => ls2%Strings(iStr)
        call HTML_Header('Line-sink string definition', 3)
        call HTML_StartTable()
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrInteger('FWL index', str%iFWLIndex)
        call HTML_AttrReal('Conductance', str%rConductance)
        call HTML_AttrInteger('Mode', str%iMode)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/ 'Vertex', 'FDP # ', 'X     ', 'Y     ', 'Length', &
                                 'Head  ', 'Cond  ', 'Depth ', 'Sigma ', 'M Head', &
                                 'Enable', 'S Flow', 'Route ', 'Err % ', 'Err Q ' /))
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, vtx%iFDPIndex/))
          call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
          call HTML_ColumnReal((/vtx%rLength, vtx%rCPHead, str%rConductance, vtx%rCPDepth, &
          	   vtx%rStrength, rAQU_PotentialToHead(io, aqu, vtx%rCheckPot, vtx%cCPZ(1))/))
          call HTML_ColumnLogical((/vtx%lEnabled/))
          call HTML_ColumnReal((/vtx%rStreamFlow/))
          call HTML_ColumnLogical((/vtx%lRouteEnabled/))
          if (vtx%lEnabled .and. (vtx%rCheckHead /= vtx%rCPHead)) then
            rCheckFlux = str%rConductance*(vtx%rCheckHead-vtx%rCPHead)
            call HTML_ColumnReal((/rHUNDRED*(vtx%rStrength-rCheckFlux)/rCheckFlux, &
                 vtx%rLength*(vtx%rStrength-rCheckFlux)/))
          else
            call HTML_ColumnText((/'--', '--'/))
          end if
          call HTML_EndRow()
        end do
        vtx => str%Vertices(str%iNPts)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iVtx/))
        call HTML_ColumnText((/'--'/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
        call HTML_ColumnText((/'--', '--', '--', '--', '--', '--'/))
        call HTML_EndRow()
        call HTML_EndTable()
      end do
    else
      call HTML_Header('No line-sinks defined', 3)
    end if

    return
  end subroutine LS2_Report


  subroutine LS2_Save(io, ls2, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(LS2_COLLECTION), pointer :: ls2
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istr, ivtx, iflg
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx

    ! Output records will be of the form ELEM_LS2, IWEL, IRAD, IVTX, 0, SIGMA
    do istr = 1, ls2%iNStr
      str => ls2%Strings(istr)
      do ivtx = 1, str%iNPts
        vtx => str%Vertices(ivtx)
        if (mode == IO_MODE_BINARY) then
          write (unit=LU_SCRATCH) ELEM_LS2, istr, ivtx, 1, vtx%rStrength
        else
          write (unit=LU_SCRATCH, fmt=*) "LS2", istr, ivtx, 1, vtx%rStrength
        end if
      end do
    end do

    return
  end subroutine LS2_Save


  subroutine LS2_Load(io, ls2, fwl, fdp, mode)
    !! Loads the LS2 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(LS2_COLLECTION), pointer :: ls2
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, istr, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(LS2_STRING), pointer :: str
    type(LS2_VERTEX), pointer :: vtx
    type(LS2_VERTEX), pointer :: flg

    ! Scans the entire precondition file for the LS2 data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, istr, ivtx, iflg, rstrength
        if (imodule /= ELEM_LS2) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, istr, ivtx, iflg, rstrength
        if (uppercase (trim(smodule)) /= "LS2") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      call IO_Assert(io, istr > 0 .and. istr <= ls2%iNStr, "LS2 string not found")
      str => ls2%Strings(istr)
      call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "LS2 vertex not found")
      vtx => str%Vertices(ivtx)
      call IO_Assert(io, iflg == 1, "LS2 strength index not found")
      vtx%rStrength = rStrength
    end do

    ! Now, populate the internal data structures
    call LS2_Update(io, ls2, fwl, fdp)

    return
  end subroutine LS2_Load

end module m_ls2
