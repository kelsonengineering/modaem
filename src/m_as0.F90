module m_as0

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

  !! module m_as0
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
  use u_polygon
  use u_matrix
  use f_well
  use f_dipole
  use i_areasink
  use m_aqu

  implicit none

  public

  type :: AS0_STRING
    !! type AS0_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   complex, pointer :: poly(:)
    !!     A vector of coordinates that define the polygon
    !!   integer, pointer :: iFDPIndex
    !!     A vector of indices into an FDP_COLLECTION object
    !!     (a one-to-one correspondance exists between elements
    !!     of poly and iFDPIndex)
    !!   real :: rN
    !!     The exfiltration rate
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !!
    complex(kind=AE_REAL), dimension(:), pointer :: poly
    integer(kind=AE_INT), dimension(:), pointer :: iFDPIndex
    real(kind=AE_REAL) :: rArea
    real(kind=AE_REAL) :: rN
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iFWLIndex
    integer(kind=AE_INT) :: iID
    ! Active area of the area sink
    real(kind=AE_REAL) :: rActiveArea
  end type AS0_STRING

  type :: AS0_COLLECTION
    !! type AS0_COLLECTION
    !!
    !! Type that holds information for all AS0 elements in a layer
    !!
    !! Members:
    !!   integer :: iAS0Flag
    !!     Flag -- is the collection at top or bottom? Only used for the inquiry routine.
    !!   type(AS0_STRING), dimension(:), pointer :: Strings
    !!     A vector of AS0_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !!
    integer(kind=AE_INT) :: iAS0Flag
    type(AS0_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
  end type AS0_COLLECTION


contains


  function AS0_Create(iAS0Flag) result(as0)
    !! function AS0_Create
    !!
    !! Creates a new AS0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    as0 => AS0_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in) :: iAS0Flag
    ! [ RETURN VALUE ]
    type(AS0_COLLECTION), pointer :: as0
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStat

    allocate(io, as0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AS0_Create: allocation failed")
    nullify(as0%Strings)
    as0%iAS0Flag = iAS0Flag
    as0%iNStr = 0

    return
  end function AS0_Create


  subroutine AS0_Alloc(io, as0)
    !! Subroutine AS0_Alloc
    !!
    !! Allocates Strings for the AS0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call AS0_Alloc(io, as0, iNStr)
    !!
    !! Arguments:
    !!    (in)    type(AS0_COLLECTION), pointer :: as0
    !!              The AS0_COLLECTION object to be used
    !!    (in)    integer :: iNStr
    !!              The number of strings to make space for
    !!
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iStat

    iNStr = iIO_GetInteger(io, "iNStr", minimum=0)
    allocate(as0%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "AS0_Alloc: allocation failed")
    as0%Strings(:)%iID = -1

    return
  end subroutine AS0_Alloc


  !**pd  New AS0_Destroy subroutine

  subroutine AS0_Destroy(io, as0)
    !! subroutine AS0_Destroy
    !!
    !! Frees memory allocated for AS0 Area-sinks and strings of vertices
    !! and the AS0 Collection object
    !!
    !! Calling Sequence:
    !!     call AS0_Destroy(as0)
    !!
    !! Arguments:
    !!  type(AS0_COLLECTION), pointer :: as0
    !!              Pointer to the AS0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: i
    type(AS0_STRING), pointer :: str


    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_Destroy: AS0_Create has not been called")
    end if

    ! First deallocate each string of vertices
    do i = 1, as0%iNStr
      str => as0%Strings(i)
      deallocate(str%poly, stat = iStat)
      call IO_Assert(io, (iStat == 0), "AS0_Destroy: deallocation of Vertices failed")
    end do
    ! Then deallocate the strings table
    if (associated(as0%Strings)) then
      deallocate(as0%Strings, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "AS0_Destroy: deallocation of Strings failed")
    end if
    ! Now deallocate the collection
    deallocate(as0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AS0_Destroy: deallocation failed")

    return
  end subroutine AS0_Destroy



  subroutine AS0_New(io, as0, cZ, iNPts, rN)
    !! function AS0_New
    !!
    !! Adds a new AS0_STRING object to the AS0_COLLECTION 'as0'
    !!
    !! Calling Sequence:
    !!    call AS0_New(io, as0, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(AS0_COLLECTION) :: as0
    !!              The AS0_COLLECTION object to be used
    !!    (in)    complex :: cZ(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    real :: rN
    !!              Exfiltration rates at the top and bottom of the element, respectively
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    complex(kind=AE_REAL), dimension(:) :: cZ
    real(kind=AE_REAL) :: rN
    integer(kind=AE_INT), intent(in) :: iNPts
    type(IO_STATUS), pointer :: io
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStat
    type(AS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_New: AS0_Create has not been called")
    end if

    call IO_Assert(io, (as0%iNStr < size(as0%Strings)), &
         "AS0_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(cZ)), &
         "AS0_New: Size of provided vertices is inconsistent")

    as0%iNStr = as0%iNStr + 1
    str => as0%Strings(as0%iNStr)
    allocate(str%poly(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "AS0_New: Allocation failed")
    str%poly = cZ(1:iNPts)
    str%iNPts = iNPts

    return
  end subroutine AS0_New


  subroutine AS0_PreSolve(io, as0, aqu)
    !! subroutine AS0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call AS0_PreSolve(as0)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer :: as0
    !!             AS0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iStat
    complex(kind=AE_REAL), dimension(:), allocatable :: poly
    type(AS0_STRING), pointer :: str

    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      ! First, fix up the polygon and reallocate the space in the string
      allocate(poly(str%iNPts), stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AS0_PreSolve: Allocation failed')
      poly = str%poly(1:str%iNPts)
      call PGN_OrientLeftmost(poly)
      deallocate(str%poly, stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AS0_PreSolve: Deallocation failed')
      allocate(str%poly(str%iNPts), str%iFDPIndex(str%iNPts), stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AS0_PreSolve: Allocation failed')
      str%poly = poly
      str%iFDPIndex = -1
      str%rArea = PGN_Area(str%poly)
      ! Compute the active area
      str%rActiveArea = rAQU_ActiveArea(io, aqu, str%poly(1:str%iNPts))
      !    print *, "AS0: ID = ", str%iID, "Area", str%rArea, "Active", str%rActiveArea
      deallocate(poly, stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AS0_PreSolve: Deallocation failed')
    end do

    return
  end subroutine AS0_PreSolve


  function iAS0_GetInfo(io, as0, iOption, iIteration) result(iValue)
    !! function AS0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iAS0_GetInfo(wl0, iOption)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer :: as0
    !!             AS0_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_GetInfo: AS0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = as0%iNStr
      case (SIZE_FDP)
        do iStr = 1, as0%iNStr
          str => as0%Strings(iStr)
          iValue = iValue + str%iNPts
        end do
      case default
        iValue = 0
    end select

    return
  end function iAS0_GetInfo


  subroutine AS0_SetupFunctions(io, as0, fwl, fdp)
    !! subroutine AS0_Setup
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call AS0_Setup(io, as0, fwl, fdp)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION function object
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             FDP_COLLECTION function object
    !!   (in)    type(AAU_COLLECTION), pointer
    !!             AQU_COLLECTION overlapping the model
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStr, iSeg
    real(kind=AE_REAL) :: rPotJump, rPsiJump
    real(kind=AE_REAL), dimension(1, 1) :: rInfl
    complex(kind=AE_REAL) :: cZ0, cZC, cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(2) :: seg
    type(AS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_Setup: AS0_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "AS0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "AS0_Setup: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      cZ0 = str%poly(1)
      cRho1 = cZERO
      ! Build dipoles for all segments
      do iSeg = 1, str%iNPts  ! Set up one dipole per segment(closed!)
        seg = PGN_Segment(str%poly, iSeg)
        cZC = rHALF * (seg(1)+seg(2))
        ! Compute the dipole strengths in terms of the given sink density.
        ! Center
        rInfl = rIAS_InfluenceP(cZC, cZ0)
        rPotJump = str%rN * rInfl(1, 1)
        rInfl = rIAS_InfluenceF(seg(1), cZC, cZ0)
        rPsiJump = str%rN * rInfl(1, 1)
        cRho2 = cmplx(real(cRho1)-rPsiJump, -rPotJump, AE_REAL)
        ! End 2
        rInfl = rIAS_InfluenceP(seg(2), cZ0)
        rPotJump = str%rN * rInfl(1, 1)
        rInfl = rIAS_InfluenceF(seg(1), seg(2), cZ0)
        rPsiJump = str%rN * rInfl(1, 1)
        cRho3 = cmplx(real(cRho1)-rPsiJump, -rPotJump, AE_REAL)
        ! Now, make the dipole!
        call FDP_New(io, fdp, seg(1), seg(2), (/cRho1, cRho2, cRho3/), ELEM_AS0, iStr, iSeg, -1, str%iFDPIndex(iSeg))
        cRho1 = cRho3                               ! Move on to the next one with this Rho value
      end do

      ! Put a well at the end of the string
      call FWL_New(io, fwl, cZ0, real(cRho3), rZERO, ELEM_AS0, iStr, -1, -1, str%iFWLIndex)

    end do

    return
  end subroutine AS0_SetupFunctions


  subroutine AS0_FindStringPointer(io, as0, iASID, ASString, lFound)
    !! subroutine AS0_FindStringPointer
    !!
    !! Finds the area-sink string specified by the ID and returns a pointer to it
    !!
    !! Calling Sequence:
    !!    call AS0_FindStringPointer(io, as0, iASID, ASString, lfound)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer :: as0
    !!             AS0_COLLECTION to be used
    !!   (in)    integer :: iASID
    !!             The area-sink string ID number
    !!   (out)   type(AS0_STRING) :: ASString
    !!             Pointer to the area-sink string
    !!   (out)   logical :: lFound
    !!             .true. if the well was found
    !!             .false. if the well was not found
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    integer(kind=AE_INT), intent(in) :: iASID
    type(AS0_STRING), pointer :: ASString
    logical :: lFound
    type(IO_STATUS), pointer :: io
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_FindStringPointer: AS0_Create has not been called")
    end if

    lFound = .false.
    do i = 1, as0%iNStr
      ASString => as0%Strings(i)
      if (ASString%iID == iASID) then
        lFound = .true.
        return
      end if
    end do

    return
  end subroutine AS0_FindStringPointer


  subroutine AS0_Read(io, as0)
    !! subroutine AS0_Read
    !!
    !! Reads the line-sinks for the specified AS0_COLLECTIOn from kIO_InputLU
    !!
    !! Calling Sequence:
    !!    call AS0_Read(as0)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer :: as0
    !!             AS0_COLLECTION to be populated
    !!
    !! The format of the AS0 section of the input file appears as follows:
    !! AS0
    !! STR NVertices N_top N_bottom ID
    !!  (x, y)
    !! ... Up to NVertices
    !! STR NVertices N_top N_bottom ID
    !!  (x, y)
    !! ... Up to NVertices
    !! ... Up to NStrings
    !!
    !! NOTE: It is assumed that the AS0 line was found by the caller

    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(2), parameter :: dirDirectives = (/dirEND, dirSTR/)
    ! [ LOCAAS ]
    real(kind=AE_REAL) :: rN
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxStr, iMaxVtx
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iStr
    logical :: lFlag
    type(AS0_STRING), pointer :: str

    call IO_MessageText(io, "  Reading AS0 module input")

    call IO_Assert(io, (associated(as0)), "AS0_Read: AS0_Create was not called")

    ! Use ifIO_InputRecord to process the model input file.
    nullify(str)
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "AS0_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModAS0 "ifXXXRead" routines.
          call IO_Assert(io, .false., "AS0_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(str)), "AS0_Read: No STR command was found")
          call IO_Assert(io, (str%iNPts < size(str%poly)), "AS0_Read: Space exhausted")
          cZ = cIO_GetCoordinate(io, "cZ", extents=.true.)
          str%iNPts = str%iNPts+1
          str%poly(str%iNPts) = cZ
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpSTR)
          !****************************************************************************
          ! Here for the STR command -- create a new string of line-sinks
          ! the maximum number of vertices is in the input record
          !****************************************************************************
          call IO_Assert(io, (associated(as0%Strings)), "AS0_Read: No space is allocated")
          call IO_Assert(io, (as0%iNStr < size(as0%Strings)), "AS0_Read: Space exhausted")
          iMaxVtx = iIO_GetInteger(io, "iMaxVtx", minimum=3)
          rN = rIO_GetReal(io, "rN")
          iID = iIO_GetInteger(io, "iID", forbidden=as0%Strings(:)%iID)
          as0%iNStr = as0%iNStr+1
          str => as0%Strings(as0%iNStr)
          allocate(str%poly(iMaxVtx), stat = iStat)
          call IO_Assert(io, (iStat == 0), "AS0_Read: Allocation failed")
          ! Made it!
          str%iFWLIndex = -1         ! No FWL function yet!
          str%rN = rN
          str%iID = iID
          str%iNPts = 0      ! Initialize the vertex counter
      end select
    end do

    call IO_MessageText(io, "  Leaving AS0 module")

  end subroutine AS0_Read


  subroutine AS0_Inquiry(io, as0, iLU)
    !! subroutine AS0_Inquiry
    !!
    !! Writes an inquiry report for all area-sinks to iLU
    !!
    !! Calling Sequence:
    !!    call AS0_Inquiry(io, as0, iLU)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_Inquiry: AS0_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#AS0, ID, FLAG, AREA, STRENGTH, ACTIVE_A"")")
    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      write (unit=iLU, &
             fmt="(""AS0"", 2("", "", i9), 3("", "", e16.8))" &
             ) str%iID, as0%iAS0Flag, str%rArea, str%rN, str%rActiveArea
    end do

    return
  end subroutine AS0_Inquiry


  subroutine AS0_Report(io, as0, sLabel)
    !! subroutine AS0_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call AS0_Report(as0)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(IO_STATUS), pointer :: io
    character(len=*), optional :: sLabel
    ! [ LOCAAS ]
    integer(kind=AE_INT) :: iStr, iVtx
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    type(AS0_STRING), pointer :: str
    character(len=32) :: sMyLabel

    if (present(sLabel)) then
      sMyLabel = sLabel
    else
      sLabel = ''
    end if

    if (io%lDebug) then
      call IO_Assert(io, (associated(as0)), &
           "AS0_Inquiry: AS0_Create has not been called")
    end if

    call HTML_Header('Module AS0 ' // trim(sMyLabel), 1)
    call HTML_Header('Discharge-specified area sink information', 2)

    if (associated(as0%Strings)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of area-sinks', as0%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iAS0_GetInfo(io, as0, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iAS0_GetInfo(io, as0, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iAS0_GetInfo(io, as0, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iAS0_GetInfo(io, as0, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iAS0_GetInfo(io, as0, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      call HTML_Header('Area-sink definition', 3)

      do iStr = 1, as0%iNStr
        str => as0%Strings(iStr)
        call HTML_StartTable()
        call HTML_AttrInteger('Area-sink #', str%iID)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrReal('Gamma', str%rN)
        call HTML_AttrReal('Area', str%rArea)
        call HTML_AttrReal('Active Area', str%rActiveArea)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)

        call HTML_StartTable()
        call HTML_TableHeader((/'      ', 'FDP # ', 'X     ', 'Y     '/))
        do iVtx = 1, str%iNPts
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx, str%iFDPIndex(iVtx)/))
          call HTML_ColumnComplex((/str%poly(iVtx)/))
          call HTML_EndRow()
        end do
        call HTML_EndTable()
      end do

    else
      call HTML_Header('No area-sinks defined', 3)
    end if

    return
  end subroutine AS0_Report

  ! Module AS0 requires some additional functions for the computation of potentials. These
  ! are the "InsidePotential", "InsideDischarge", "InsideRecharge", and "InsideInfiltration"
  ! functions for the module


  function rAS0_InsidePotential(io, as0, cZ) result(rPot)
    !! function AS0_InsidePotential
    !!
    !! Returns the potential contributions of all the inside functions at the point cZ
    !!
    !! Calling Sequence:
    !!    call AS0_InsidePotential(io, as0, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in question
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPot
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    rPot = rZERO
    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      if (PGN_Contains(str%poly, cZ)) then
        rP = rIAS_InfluenceP(cZ, str%poly(1))
        rPot = rPot + str%rN*rP(1, 1)
      end if
    end do

    return
  end function rAS0_InsidePotential


  function cAS0_InsideDischarge(io, as0, cZ) result(cW)
    !! function AS0_InsidePotential
    !!
    !! Returns the discharge contributions of all the inside functions at the point cZ
    !!
    !! Calling Sequence:
    !!    call AS0_InsideDischarge(io, as0, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in question
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cW
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cP(1, 1)
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    cW = rZERO
    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      if (PGN_Contains(str%poly, cZ)) then
        cP = conjg(cIAS_InfluenceW(cZ, str%poly(1)))
        cW = cW + str%rN*cP(1, 1)
      end if
    end do

    return
  end function cAS0_InsideDischarge


  function rAS0_InsideRecharge(io, as0, cZ) result(rG)
    !! function AS0_InsideRecharge
    !!
    !! Returns the recharge contributions of all the inside functions at the point cZ
    !!
    !! Calling Sequence:
    !!    call AS0_InsideRecharge(io, as0, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in question
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    rG = rZERO
    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      if (PGN_Contains(str%poly, cZ)) then
        rP = rIAS_InfluenceG(cZ, str%poly(1))
        rG = rG + str%rN*rP(1, 1)
      end if
    end do

    return
  end function rAS0_InsideRecharge


  function rAS0_InsideFlow(io, as0, cPathZ) result(rF)
    !! function AS0_InsideRecharge
    !!
    !! Returns the inside flow for all area-sinks across a path
    !!
    !! Calling Sequence:
    !!    call AS0_InsideRecharge(io, as0, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!   (in)    complex :: cZ(:)
    !!             The path to be examined
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rF
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iZ, iStr, iSeg
    complex(kind=AE_REAL), dimension(:), allocatable :: cSegChop
    complex(kind=AE_REAL), dimension(2) :: seg
    integer(kind=AE_INT) :: iChop
    type(AS0_STRING), pointer :: str

    rF = rZERO
    do iZ = 1, size(cPathZ)-1
      do iStr = 1, as0%iNStr
        str => as0%Strings(iStr)
        allocate(cSegChop(str%iNPts+1))
        call PGN_ChopSegment(str%poly, (/cPathZ(iZ), cPathZ(iZ+1)/), cSegChop, iChop)
        do iSeg = 1, iChop-1
          seg = PGN_Segment(cSegChop, iSeg)
          if (PGN_Contains(str%poly, rHALF*(seg(1)+seg(2)))) then
            rP = -rIAS_InfluenceF(seg(1), seg(2), str%poly(1))
            rF = rF + str%rN*rP(1, 1)
          end if
        end do
        deallocate(cSegChop)
      end do
    end do

    return
  end function rAS0_InsideFlow


  function rAS0_Extraction(io, as0) result(rQ)
    !! function AS0_Extraction
    !!
    !! Returns the extraction rate of all area-sinks
    !!
    !! Calling Sequence:
    !!    call AS0_Extraction(as0)
    !!
    !! Arguments:
    !!   (in)    type(AS0_COLLECTION), pointer
    !!             AS0_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(AS0_COLLECTION), pointer :: as0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr
    type(AS0_STRING), pointer :: str

    rQ = rZERO
    do iStr = 1, as0%iNStr
      str => as0%Strings(iStr)
      rQ = rQ - str%rN * str%rActiveArea
    end do
    rQ = rZERO

    return
  end function rAS0_Extraction

end module m_as0
