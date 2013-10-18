module m_ls0

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

  !! module m_ls0
  !!
  !! Element module for 2-D discharge specified line-sinks
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
  use u_matrix
  use f_well
  use f_dipole
  use m_aqu

  implicit none

  public

  type :: LS0_VERTEX
    !! type LS0_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rSigma
    !!     The sink density at the vertex
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rSigma
    integer(kind=AE_INT) :: iFDPIndex
    real(kind=AE_REAL) :: rCheckHead
  end type LS0_VERTEX

  type :: LS0_STRING
    !! type LS0_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   type(LS0_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of LS0_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(LS0_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iFWLIndex
    integer(kind=AE_INT) :: iID
  end type LS0_STRING

  type :: LS0_COLLECTION
    !! type LS0_COLLECTION
    !!
    !! Type that holds information for all LS0 elements in a layer
    !!
    !! Members:
    !!   type(LS0_STRING), dimension(:), pointer :: Strings
    !!     A vector of LS0_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !!
    type(LS0_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    ! Iterator Information
    integer(kind=AE_INT) :: iIterStr
    integer(kind=AE_INT) :: iIterVtx
  end type LS0_COLLECTION


contains


  function LS0_Create(io) result(ls0)
    !! function LS0_Create
    !!
    !! Creates a new LS0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    ls0 => LS0_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(ls0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS0_Create: allocation failed")
    nullify(ls0%Strings)
    ls0%iNStr = 0

    return
  end function LS0_Create


  subroutine LS0_Alloc(io, ls0)
    !! Subroutine LS0_Alloc
    !!
    !! Allocates Strings for the LS0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call LS0_Alloc(ls0)
    !!
    !! Arguments:
    !!    (in)    type(LS0_COLLECTION), pointer :: ls0
    !!              The LS0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iStat

    iNStr = iIO_GetInteger(io, 'iNStr', minimum = 0)
    allocate(ls0%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS0_Alloc: allocation failed")

  end subroutine LS0_Alloc


  subroutine LS0_Destroy(io, ls0)
    !! subroutine LS0_Destroy
    !!
    !! Frees memory allocated for LS0 Linesinks and strings of vertices
    !! and the LS0 Collection object
    !!
    !! Calling Sequence:
    !!     call LS0_Destroy(ls0)
    !!
    !! Arguments:
    !!  type(LS0_COLLECTION), pointer :: ls0
    !!              Pointer to the LS0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: i
    type(LS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_Destroy: LS0_Create has not been called")
    end if

    ! First deallocate each string of vertices
    do i = 1, ls0%iNStr
      str => ls0%Strings(i)
      deallocate(str%Vertices, stat = iStat)
      call IO_Assert(io, (iStat == 0), "LS0_Destroy: deallocation of Vertices failed")
    end do
    ! Then deallocate the strings table
    if (associated(ls0%Strings)) then
      deallocate(ls0%Strings, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "LS0_Destroy: deallocation of Strings failed")
    end if
    ! Now deallocate the collection
    deallocate(ls0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS0_Destroy: deallocation failed")

    return
  end subroutine LS0_Destroy


  subroutine LS0_New(io, ls0, Vertices, iNPts)
    !! function LS0_New
    !!
    !! Adds a new LS0_STRING object to the LS0_COLLECTION 'ls0'
    !!
    !! Calling Sequence:
    !!    call LS0_New(io, ls0, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(LS0_COLLECTION) :: ls0
    !!              The LS0_COLLECTION object to be used
    !!    (in)    type(LS0_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(LS0_VERTEX), dimension(:) :: Vertices
    integer(kind=AE_INT), intent(in) :: iNPts
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(LS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_New: LS0_Create has not been called")
    end if

    call IO_Assert(io, (ls0%iNStr < size(ls0%Strings)), &
         "LS0_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(Vertices)), &
         "LS0_New: Size of provided vertices is inconsistent")

    ls0%iNStr = ls0%iNStr + 1
    str => ls0%Strings(ls0%iNStr)
    allocate(str%Vertices(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "LS0_New: Allocation failed")
    str%Vertices = Vertices(1:iNPts)
    str%iNPts = iNPts

    return
  end subroutine LS0_New


  function iLS0_GetID(io, ls0, iIndex) result(iID)
    !! Returns the ID number for the string at index 'iIndex'
    type(LS0_COLLECTION), pointer :: ls0
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= ls0%iNStr), "Internal error -- no such index")
    iID = ls0%Strings(iIndex)%iID

    return
  end function iLS0_GetID


  subroutine LS0_PreSolve(io, ls0)
    !! subroutine LS0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call LS0_PreSolve(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io

    return
  end subroutine LS0_PreSolve


  function iLS0_GetInfo(io, ls0, iOption, iIteration) result(iValue)
    !! function LS0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iLS0_GetInfo(wl0, iOption)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr
    type(LS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_GetInfo: LS0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = ls0%iNStr
      case (SIZE_FDP)
        do iStr = 1, ls0%iNStr
          str => ls0%Strings(ls0%iNStr)
          iValue = iValue + str%iNPts-1
        end do
      case default
        iValue = 0
    end select

    return
  end function iLS0_GetInfo


  subroutine LS0_SetupFunctions(io, ls0, fwl, fdp)
    !! subroutine LS0_Setup
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call LS0_Setup(io, ls0, fwl, fdp)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer
    !!             LS0_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION function object
    !!   (in)    type(LS0_COLLECTION), pointer
    !!             FDP_COLLECTION function object
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rSigma, rDisch
    complex(kind=AE_REAL) :: cZ1, cZ2, cRho1, cRho2, cRho3
    type(LS0_STRING), pointer :: str
    type(LS0_VERTEX), pointer :: this_vtx, next, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_Setup: LS0_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "LS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "LS0_Setup: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, ls0%iNStr
      str => ls0%Strings(iStr)
      cRho1 = cZERO
      ! Build dipoles for all segments
      do iVtx = 1, str%iNPts-1  ! Set up nVertices-1 dipoles...
        this_vtx => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        cZ1 = this_vtx%cZ
        cZ2 = next%cZ

        ! Compute the dipole strengths in terms of the given sink density.  NOTE: these
        ! computations may be easily adjusted for linear-strength line-sinks.
        rSigma = rHALF * (this_vtx%rSigma + next%rSigma)
        rDisch = rSigma*abs(cZ2-cZ1)
        cRho3 = cRho1 + rDisch
        cRho2 = rHALF*(cRho1+cRho3)
        call FDP_New(io, fdp, cZ1, cZ2, (/cRho1, cRho2, cRho3/), ELEM_LS0, iStr, iVtx, -1, this_vtx%iFDPIndex)
        cRho1 = cRho3                               ! Move on to the next one with this Rho value
      end do

      ! Put a well at the end of the string
      last_vtx => str%Vertices(str%iNPts)
      cZ1 = last_vtx%cZ
      call FWL_New(io, fwl, cZ1, real(cRho3), rZERO, ELEM_LS0, iStr, -1, -1, str%iFWLIndex)
    end do

    return
  end subroutine LS0_SetupFunctions


  subroutine LS0_FindStringPointer(io, ls0, iLSID, LSString, lFound)
    !! subroutine LS0_FindStringPointer
    !!
    !! Finds the linesink string specified by the ID and returns a pointer to it
    !!
    !! Calling Sequence:
    !!    call LS0_FindStringPointer(io, ls0, iLSID, LSString, lfound)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (in)    integer :: iLSID
    !!             The linesink string ID number
    !!   (out)   type(LS0_STRING) :: LSString
    !!             Pointer to the linesink string
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    integer(kind=AE_INT), intent(in) :: iLSID
    type(LS0_STRING), pointer :: LSString
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_FindStringPointer: LS0_Create has not been called")
    end if

    lFound = .false.
    do i = 1, ls0%iNStr
      LSString => ls0%Strings(i)
      if (LSString%iID == iLSID) then
        lFound = .true.
        return
      end if
    end do

    return
  end subroutine LS0_FindStringPointer


  subroutine LS0_ResetIterator(io, ls0)
    !! subroutine LS0_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call LS0_ResetIterator(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_ResetIterator: LS0_Create has not been called")
    end if

    ls0%iIterStr = 1
    ls0%iIterVtx = 0

    return
  end subroutine LS0_ResetIterator


  function LS0_NextIterator(io, ls0) result(itr)
    !! function LS0_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS0_NextIterator(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: ls0
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_NextIterator: LS0_Create has not been called")
    end if

    if (ls0%iIterStr > ls0%iNStr) then
      nullify(itr)
      return
    end if

    ls0%iIterVtx = ls0%iIterVtx + 1
    if (ls0%iIterVtx > ls0%Strings(ls0%iIterStr)%iNPts-1) then
      ls0%iIterStr = ls0%iIterStr+1
      ls0%iIterVtx = 1
      if (ls0%iIterStr > ls0%iNStr) then
        nullify(itr)
        return
      end if
    end if

    allocate(itr)
    itr%iElementType = ELEM_AQU
    itr%iElementString = ls0%iIterStr
    itr%iElementVertex = ls0%iIterVtx
    itr%iValueSelector = VALUE_HEAD
    allocate(itr%cZ(1))
    itr%cZ(1) = rHALF * (ls0%Strings(ls0%iIterStr)%Vertices(ls0%iIterVtx)%cZ + &
                ls0%Strings(ls0%iIterStr)%Vertices(ls0%iIterVtx+1)%cZ)

    return
  end function LS0_NextIterator


  subroutine LS0_SetIterator(io, ls0, aqu, itr, cValue)
    !! function LS0_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call LS0_SetIterator(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be used
    !!   (in)    type(ITERATOR_RESULT), pointer :: itr
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: ls0
    !!             Tracks error conditions
    !!     Pointer to the information for data retrieval
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(AQU_COLLECTION), pointer :: aqu
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_NextIterator: LS0_Create has not been called")
      call IO_Assert(io, (ls0%iIterStr <= ls0%iNStr), &
           "LS0_SetIterator: Iterator out of range")
    end if

    ls0%Strings(itr%iElementString)%Vertices(itr%iElementVertex)%rCheckHead = real(cValue, AE_REAL)

    return
  end subroutine LS0_SetIterator


  subroutine LS0_Read(io, ls0)
    !! subroutine LS0_Read
    !!
    !! Reads the line-sinks for the specified LS0_COLLECTIOn from kIO_InputLU
    !!
    !! Calling Sequence:
    !!    call LS0_Read(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer :: ls0
    !!             LS0_COLLECTION to be populated
    !!
    !! The format of the LS0 section of the input file appears as follows:
    !! LS0
    !! DIM Nstrings
    !! STR NVertices
    !! x y strength
    !! ... Up to NVertices
    !! STRING NVertices
    !! x y strength
    !! ... Up to NVertices
    !! ... Up to NStrings
    !!
    !! NOTE: It is assumed that the LS0 line was found by the caller

    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(2), parameter :: dirDirectives = (/dirEND, dirSTR/)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rSigma
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxStr, iMaxVtx
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iStr
    logical :: lFlag
    type(LS0_STRING), pointer :: str
    type(LS0_VERTEX), pointer :: vtx, prev

    call IO_MessageText(io, "  Reading LS0 module input")

    call IO_Assert(io, (associated(ls0)), "LS0_Read: LS0_Create was not called")

    ! Use ifIO_InputRecord to process the model input file.
    nullify(str, vtx)
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "LS0_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModLS0 "ifXXXRead" routines.
          call IO_Assert(io, .false., "LS0_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(str)), "LS0_Read: No STR command was found")
          call IO_Assert(io, (str%iNPts < size(str%Vertices)), "LS0_Read: Space exhausted")
          cZ = cIO_GetCoordinate(io, 'cZ', extents=.true.)
          rSigma = rIO_GetReal(io, 'rSigma')

          str%iNPts = str%iNPts+1
          vtx => str%Vertices(str%iNPts)
          vtx%cZ = cZ
          if ( str%iNPts > 1) then
            prev => str%Vertices(str%iNPts-1)
            prev%rSigma = rIO_LocalLength(io, rSigma, vtx%cZ-prev%cZ)
          end if
          vtx%iFDPIndex = -1
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpSTR)
          !****************************************************************************
          ! Here for the STR command -- create a new string of line-sinks
          ! the maximum number of vertices is in the input record
          !****************************************************************************
          call IO_Assert(io, (associated(ls0%Strings)), "LS0_Read: No space is allocated")
          call IO_Assert(io, (ls0%iNStr < size(ls0%Strings)), "LS0_Read: Space exhausted")
          ! Retrive the number of vertices desired...
          iMaxVtx = iIO_GetInteger(io, 'iMaxVtx', minimum = 2)
          iID = iIO_GetInteger(io, 'iID')
          ! Allocate the vertices...
          ls0%iNStr = ls0%iNStr+1
          str => ls0%Strings(ls0%iNStr)
          allocate(str%Vertices(iMaxVtx), stat = iStat)
          call IO_Assert(io, (iStat == 0), "LS0_Read: Allocation failed")
          ! Made it!
          str%iFWLIndex = -1         ! No FWL function yet!
          str%iID = iID
          write (unit=IO_MessageBuffer, &
                 fmt="("" LS0_Read: "", i6, "" vertices allocated"")" &
                 ) iMaxVtx
          call IO_MessageText(io)
          str%iNPts = 0      ! Initialize the vertex counter
      end select
    end do

    call IO_MessageText(io, "  Leaving LS0 module")

  end subroutine LS0_Read


  subroutine LS0_Inquiry(io, ls0, iLU)
    !! subroutine LS0_Inquiry
    !!
    !! Writes an inquiry report for all line-sinks to iLU
    !!
    !! Calling Sequence:
    !!    call LS0_Inquiry(io, ls0, iLU)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer
    !!             LS0_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    real(kind=AE_REAL) :: rLength, rSigma
    type(LS0_STRING), pointer :: str
    type(LS0_VERTEX), pointer :: this, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_Inquiry: LS0_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#LS0, ID, VTX, X1, Y1, X2, Y2, LENGTH, STRENGTH, MOD_HEAD"")")
    do iStr = 1, ls0%iNStr
      str => ls0%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        this => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        rLength = rIO_WorldLength(io, abs(next%cZ-this%cZ), next%cZ-this%cZ)
        rSigma = rIO_WorldLength(io, this%rSigma, next%cZ-this%cZ)
        write (unit=iLU, &
               fmt="(""LS0"", 2("", "", i9), 7("", "", e16.8))" &
               ) str%iID, &
               iVtx, &
               cIO_WorldCoords(io, this%cZ), &
               cIO_WorldCoords(io, next%cZ), &
               rLength, &
               rSigma, &
               this%rCheckHead
      end do
    end do

    return
  end subroutine LS0_Inquiry


  subroutine LS0_Report(io, ls0)
    !! subroutine LS0_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call LS0_Report(ls0)
    !!
    !! Arguments:
    !!   (in)    type(LS0_COLLECTION), pointer
    !!             LS0_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(LS0_COLLECTION), pointer :: ls0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    real(kind=AE_REAL) :: rSigma
    type(LS0_STRING), pointer :: str
    type(LS0_VERTEX), pointer :: vtx, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(ls0)), &
           "LS0_Inquiry: LS0_Create has not been called")
    end if

    call HTML_Header('Module LS0', 1)
    call HTML_Header('Discharge-specified area sink information', 2)

    if (associated(ls0%Strings)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of line-sinks', ls0%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iLS0_GetInfo(io, ls0, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iLS0_GetInfo(io, ls0, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iLS0_GetInfo(io, ls0, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iLS0_GetInfo(io, ls0, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iLS0_GetInfo(io, ls0, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iStr = 1, ls0%iNStr
        str => ls0%Strings(iStr)
        call HTML_Header('Line-sink string definition', 3)
        call HTML_StartTable()
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrInteger('FWL index', str%iFWLIndex)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'FDP # ', 'X     ', 'Y     ', 'Sigma ', 'M Head'/))
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          next => str%Vertices(iVtx+1)
          rSigma = rIO_WorldLength(io, vtx%rSigma, next%cZ-vtx%cZ)
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, vtx%iFDPIndex/))
          call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
          call HTML_ColumnReal((/rSigma, vtx%rCheckHead/))
          call HTML_EndRow()
        end do
        vtx => str%Vertices(str%iNPts)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iVtx/))
        call HTML_ColumnText((/'--'/))
        call HTML_ColumnComplex((/vtx%cZ/))
        call HTML_ColumnText((/'--', '--'/))
        call HTML_EndRow()
        call HTML_EndTable()
      end do
    else
      call HTML_Header('No line-sinks defined', 3)
    end if

    return
  end subroutine LS0_Report

end module m_ls0
