module u_io

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

  !! io.f90 -- Common I/O framework for ModAEM
  !!
  !! This is a Fortran MODULE for common Input/Output functions
  !!
  !! Revision History
  !!     v0.1  1/20/97  Vic Kelson
  !!                    Initial development version
  !!     v1.0  6/23/98  Vic Kelson
  !!                    Production(WhAEM) version
  !!     v1.1  4/02/01  Vic Kelson
  !!                    Code cleanup for LGPL release plans
  !!
  !! This module contains common functions for Input/Output operations,
  !! and provides an overall framework for ModAEM I/O.  All ModAEM I/O
  !! is handled by use of this framework. As a result, the code can
  !! be easily ported to new I/O frameworks, e.g. GUIs or networks.

  use u_constants

  implicit none

  public

  integer(kind=AE_INT), private, parameter :: IO_MAX_FIELDS = 50
  integer(kind=AE_INT), private, parameter :: IO_MAX_TUNING = 50

  type, public :: IO_STATUS
    logical :: lFilesOpen
    integer(kind=AE_INT) :: iRecord
    character(len=80) :: sErrorText
    ! Coordinate transformations
    complex(kind=AE_REAL) :: cZOrigin
    real(kind=AE_REAL) :: rScaling
    real(kind=AE_REAL) :: rRotation
    real(kind=AE_REAL), dimension(3, 3) :: rWtoL, rLtoW
    ! PROFILE mode features
    logical :: lProfile
    real(kind=AE_REAL) :: rKHoverKV
    real(kind=AE_REAL) :: rSqrtA
    ! Region-windowing operations
    real(kind=AE_REAL) :: rXMin, rXMax
    real(kind=AE_REAL) :: rYMin, rYMax
    ! Buffers for string start / end
    integer(kind=AE_INT) :: iFieldCount
    integer(kind=AE_INT) :: iThisField
    character(len=255), dimension(IO_MAX_FIELDS) :: sFields
    ! Exception count is increased by each exception
    integer(kind=AE_INT) :: iExceptionCount
    ! Optional tuning parameters (first letter of the name indicates the type)
    character (len=255), dimension(IO_MAX_TUNING) :: sTuning
    integer(kind=AE_INT) :: iTuningCount
    character (len=255), dimension(IO_MAX_TUNING) :: sTuningValues
    real(kind=AE_REAL), dimension(IO_MAX_TUNING) :: rTuningValues
    complex(kind=AE_REAL), dimension(IO_MAX_TUNING) :: cTuningValues
    logical, dimension(IO_MAX_TUNING) :: lTuningValues
    ! If true, debug reporting is enabled
    logical :: lDebug
    ! If true, exception crashing is enabled
    logical :: lCrash
  end type IO_STATUS

  ! I/O Logical Unit Assignments
  ! Only the output LUs are a 'PUBLIC' parameter. It is expected that the
  ! remaining LUs are programmatically accessed by support routines.
  integer(kind=AE_INT), private, parameter :: LU_INPUT = 20
  integer(kind=AE_INT), public, parameter :: LU_SCRATCH = 21
  integer(kind=AE_INT), private, parameter :: LU_ERROR = 22
  integer(kind=AE_INT), public, parameter :: LU_OUTPUT = 23
  integer(kind=AE_INT), public, parameter :: LU_GRID = 24
  integer(kind=AE_INT), public, parameter :: LU_INQ = 25
  integer(kind=AE_INT), public, parameter :: LU_CHECK = 26
  integer(kind=AE_INT), public, parameter :: LU_TRACE = 27
  integer(kind=AE_INT), public, parameter :: LU_EXT_INPUT = 28
  integer(kind=AE_INT), public, parameter :: LU_EXT_OUTPUT = 29

  ! Error status codes
  integer(kind=AE_INT), public, parameter :: IO_OK = 0
  integer(kind=AE_INT), public, parameter :: IO_BUFFEREMPTY = -1
  integer(kind=AE_INT), public, parameter :: IO_READERROR = -2
  ! Modes for save/load
  integer(kind=AE_INT), public, parameter :: IO_MODE_BINARY = 0
  integer(kind=AE_INT), public, parameter :: IO_MODE_ASCII = 1

  ! Place for writing messages
  character(len=255), public :: IO_MessageBuffer


contains


  function IO_Create(cZOrigin, rScaling, rRotation, lDebug, lCrash, lProfile, rKHoverKV) result(io)
    !! Creates an IO_STATUS object
    !!
    !! Arguments:
    !!    complex :: cZOrigin
    !!        The local origin for coordinate transformations
    !!    real :: rScale
    !!        Scaling factor(e.g. meters to ft)
    !!    real :: rRotation
    !!        Rotation of transformed domain(degrees CCW)
    !!    logical :: lDebug
    !!        Report debugging messages?
    !!    logical :: lCrash
    !!        Crash on error?
    !!    logical :: lProfile
    !!        True if model is in PROFILE mode
    !!    real :: rKHoverKV
    !!        VERTICAL anisotropy of the local model domain if in PROFILE mode
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZOrigin
    real(kind=AE_REAL), intent(in) :: rScaling
    real(kind=AE_REAL), intent(in) :: rRotation
    logical, intent(in) :: lDebug
    logical, intent(in) :: lCrash
    logical, intent(in) :: lProfile
    real(kind=AE_REAL), intent(in) :: rKHoverKV
    ! [ RETURN VALUE ]
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    real(kind=AE_REAL) :: rTheta
    real(kind=AE_REAL), dimension(3, 3) :: rT, rR, rS, rTi, rRi, rSi

    allocate(io)
    io%lFilesOpen = .false.
    io%iRecord = 0
    io%sErrorText = ""

    io%cZOrigin = cZOrigin
    io%rScaling = rScaling
    io%rRotation = rRotation
    rT = reshape((/rONE, rZERO, rZERO, &
                   rZERO, rONE, rZERO, &
                  -real(cZOrigin, AE_REAL), -aimag(cZOrigin), rONE/), &
                 (/3, 3/))
    rTi = reshape((/rONE, rZERO, rZERO, &
                    rZERO, rONE, rZERO, &
                    real(cZOrigin, AE_REAL), aimag(cZOrigin), rONE/), &
                  (/3, 3/))

    rS = reshape((/rScaling, rZERO, rZERO, &
                   rZERO, rScaling, rZERO, &
                   rZERO, rZERO, rONE/), &
                 (/3, 3/))
    rSi = reshape((/rONE/rScaling, rZERO, rZERO, &
                    rZERO, rONE/rScaling, rZERO, &
                    rZERO, rZERO, rONE/), &
                  (/3, 3/))

    rTheta = rRotation * rPIOVER180
    rR = reshape((/cos(rTheta), sin(rTheta), rZERO, &
                   -sin(rTheta), cos(rTheta), rZERO, &
                   rZERO, rZERO, rONE/), &
                 (/3, 3/))
    rRi = reshape((/cos(rTheta), -sin(rTheta), rZERO, &
                    sin(rTheta), cos(rTheta), rZERO, &
                    rZERO, rZERO, rONE/), &
                  (/3, 3/))

    io%rWtoL = matmul(matmul(rT, rR), rS)
    io%rLtoW = matmul(matmul(rSi, rRi), rTi)

    ! Profile entries
    io%lProfile = lProfile
    io%rKHoverKV = rKHoverKV
    io%rSqrtA = sqrt(rKHoverKV)

    ! Initialize the model element window
    io%rXMin = rHUGE
    io%rXMax = -rHUGE
    io%rYMin = rHUGE
    io%rYMax = -rHUGE

    io%iFieldCount = 0
    io%iThisField = 0
    io%sFields = repeat(" ", 255)

    io%lDebug = lDebug
    io%lCrash = lCrash

    return
  end function IO_Create


  subroutine IO_OpenAll(io, sBaseName)
    !! subroutine IO_OpenAll
    !!
    !! Opens the standard files for a ModAEM run. The standard I/O channels
    !! opened are: (1) the input LU [LU_INPUT]; (2) the error output LU
    !! [LU_ERROR]; (3) the report output LU [LU_OUTPUT]. All input is
    !! managed using the IO_InputRecord routine -- it is declared PRIVATE.
    !! General purpose reporting from add-on modules is allowed -- write to
    !! LU_OUTPUT.
    !!
    !! Arguments:
    !!    character(len=*) :: sBaseName
    !!        Contains the "base" file name for the standard LUs. The
    !!        extensions added are: .aem(input), .out(output) and
    !!        .err(error).

    ! [ ARGUMENTS ]
    character(len=*), intent(in) :: sBaseName
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sFName, sRecord, sName
    real(kind=AE_REAL) :: rValue
    integer(kind=AE_INT) :: iStatus
    logical :: lExists

    ! Initialize the input LU.  This will always be a disk file with
    ! the extension ".aem"
    sFName = trim(sBaseName) // ".aem"
    call IO_MessageText(io, "Processing file " // trim(sFName) // "...")
    open(unit=LU_INPUT, file = sFName, action = "READ", status="OLD", iostat=iStatus)
    call IO_Assert(io, (iStatus == 0), "IO_OpenAll: Could not open input file " // trim(sFName))

    ! Initialize the model output LU. This will always be a disk file with
    ! the extension ".out". This is the only LU which is directly written
    ! by the element modules output routines(via PUBLIC PARAMETER LU_OUTPUT
    sFName = trim(sBaseName) // ".out.html"
    open(unit=LU_OUTPUT, file = sFName, action = "WRITE", status="REPLACE", iostat=iStatus)
    call IO_Assert(io, (iStatus == 0), "IO_OpenAll: Could not open input file " // trim(sFName))

    ! Initialize the model output LU. This will always be a disk file with
    ! the extension ".out". This is the only LU which is directly written
    ! by the element modules" output routines(via PUBLIC PARAMETER LU_OUTPUT
    sFName = trim(sBaseName) // ".err"
    open(unit=LU_ERROR, file = sFName, action = "WRITE", status="REPLACE", iostat=iStatus)
    call IO_Assert(io, (iStatus == 0), "IO_OpenAll: Could not open input file " // trim(sFName))

    ! THIS CODE NEEDS TO BE COMPLETED    
    !! Load a tuning file, modaem.tun, if there is one. The tuning file
    !! has space-delimited name-value pairs. Use IO_GetTuningParameter() to
    !! retrieve a user-specified value
    !inquire (file="modaem.tun", exist=lExists)
    !if ( lExists ) then
    !    call IO_MessageText("[Reading tuning file modaem.tun]")
    !    open(unit=LU_SCRATCH, file="modaem.tun", action="READ")
    !    while (.true.)
    !        read(unit=LU_SCRATCH, )
    !        read (unit=LU_INPUT, fmt="(a255)", iostat=iStatus) sRecord
    !        if (iStatus < 0) then
    !            call IO_SplitRecord(io, sRecord)
    !            if ( io%iFieldCount == 0 ) cycle
    !            sName = sIO_GetField()

    !**********************************************************************
    ! Files have been successfully opened!
    ! Add any additional I/O initialization code here
    !**********************************************************************
    io%lFilesOpen = .true.
    io%iRecord = 0

    return
  end subroutine IO_OpenAll


  subroutine IO_CloseAll(io)
    !! subroutine IO_CloseAll
    !!
    !! Closes the standard I/O channels. NOTE: Any other open files MUST
    !! be closed by the analysis routines that opened them.
    !!
    !! Arguments: (none)
    type(IO_STATUS), pointer :: io

    if (io%lFilesOpen) then
      close(unit=LU_INPUT)
      close(unit=LU_OUTPUT)
      close(unit=LU_ERROR)
      io%lFilesOpen = .false.
    end if

  end subroutine IO_CloseAll


  subroutine IO_InputRecord(io, dirDirectives, iOpCode)
    ! iIO_InputRecord reads a record from the input file, then interprets
    ! it.  The valid directives for the current module are provided as the
    ! first parameter.  As the record is read, it is recorded(with the
    ! current line number).
    !
    ! Returns in iOpCode:
    !  kOpError  If a runtime error is detected
    !  kOpEof    If end-of-file is reached
    !  kOpData   If no directive is found
    !  (OpCode)  If read is successful and a directive is found, the OpCode for
    !            the directive is returned(for use in a SELECT CASE structure)
    !
    !  ALSO: For return values kOpPause, kOpData and directive OpCodes
    !        columns 6-end of record are returned in ca255Record
    !  ALSO: For "pause" lines, the function returns the implementation-
    !        specific keycode from ifIOMessageText()
    !
    ! Comment lines:
    !         #  If a # is found in column 1, the record is considered to be a comment
    !            and is skipped by this routine.
    !         &  If a & is found in column 1, the record is a message to the
    !            user.  It is reported using IO_MessageText()
    !
    type(DIRECTIVE), dimension(:), intent(in) :: dirDirectives
    integer(kind=AE_INT), intent(out) :: iOpCode
    type(IO_STATUS), pointer :: io
    ! Locals
    character(len=255) :: sRecord
    character(len=3) :: sTestDir
    integer(kind=AE_INT) :: iStatus
    integer(kind=AE_INT) :: ic

    ! Here we go.  Proceed only if the file is open!
    call IO_Assert(io, io%lFilesOpen, "INTERNAL ERROR -- NO FILES OPEN")
    ! The input loop reads records until valid program input is found,
    ! skipping comments and blank lines, and scanning for program directives.
    do
      ! Get a record from the input file, check for runtime error
      read (unit=LU_INPUT, fmt="(a255)", iostat=iStatus) sRecord
      if (io%lDebug) then
        ! Echo
        write (unit=IO_MessageBuffer, fmt=*) " < ", trim(sRecord), " > "
        call IO_MessageText(io)
      end if

      if (iStatus < 0) then
        iOpCode = kOpFileEOF
        return
      else if (iStatus > 0) then
        call IO_Assert(io, .false., "Read error on model input file")
      else
        ! Echo it to the run log
        io%iRecord = io%iRecord+1
        write (unit=IO_MessageBuffer, fmt="(1x, i5.5, "" >> "", a)") io%iRecord, trim(sRecord)
        call IO_ErrorText(io)
        ! Remove leading blanks and &^&^%#&!&@*&^&^!@ tabs
        do ic=1, len(sRecord)
          if (ichar(sRecord(ic:ic)) == 9) sRecord(ic:ic) = " "
        end do
        sRecord = adjustl(sRecord)
        if (sRecord(1:1) == "#" .or. len_trim(sRecord) == 0 ) cycle    ! Comment!
        ! Look for a directive...
        iOpCode = kOpData  ! Default is a data record
        do ic = 1, size(dirDirectives)
          sTestDir = uppercase (sRecord(1:3))
          if (dirDirectives(ic)%sText == sTestDir) then
            iOpCode = dirDirectives(ic)%iOpcode
            call IO_SplitRecord(io, sRecord(4:))
            return
          end if
        end do
        ! If I got here, it looks like a data record!
        call IO_SplitRecord(io, sRecord)
        return
      end if
    end do

    return
  end subroutine IO_InputRecord


  subroutine IO_SplitRecord(io, sRecord)
    !! Scans the record sRecord and stores the start/end positions of all substrings
    !! in the buffers in 'io'
    ! [ ARGUMENTS ]
    character(len=*), intent(in) :: sRecord
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iCount, iMode, iFieldStart, iFieldEnd
    integer(kind=AE_INT), parameter :: SEARCH_NONBLANK = 0
    integer(kind=AE_INT), parameter :: SEARCH_BLANK = 1
    integer(kind=AE_INT), parameter :: SEARCH_PAREN = 2
    integer(kind=AE_INT), parameter :: SEARCH_QUOTE = 3

    iMode = SEARCH_NONBLANK
    iCount = 0
    iFieldStart = 0
    iFieldEnd = 0
    io%iFieldCount = 0
    io%iThisField = 0
    do i = 1, len_trim(sRecord)
      select case (iMode)
        case (SEARCH_NONBLANK)
          if (sRecord(i:i) /= ' ') then
            if (sRecord(i:i) == '(') then
              iFieldStart = i
              iMode = SEARCH_PAREN
            else if (sRecord(i:i) == '"') then
              iFieldStart = i+1
              iMode = SEARCH_QUOTE
            else
              iFieldStart = i
              iMode = SEARCH_BLANK
            end if
          end if
        case (SEARCH_BLANK)
          if (sRecord(i:i) == ' ') then
            io%iFieldCount = io%iFieldCount + 1
            io%sFields(io%iFieldCount) = sRecord(iFieldStart:i)
            iMode = SEARCH_NONBLANK
          end if
        case (SEARCH_PAREN)
          if (sRecord(i:i) == ')') then
            io%iFieldCount = io%iFieldCount + 1
            io%sFields(io%iFieldCount) = sRecord(iFieldStart:i)
            iMode = SEARCH_NONBLANK
          end if
        case (SEARCH_QUOTE)
          if (sRecord(i:i) == '"') then
            io%iFieldCount = io%iFieldCount + 1
            io%sFields(io%iFieldCount) = sRecord(iFieldStart:i-1)
            iMode = SEARCH_NONBLANK
          end if
      end select
    end do
    call IO_Assert(io, iMode /= SEARCH_PAREN, "Unbalanced parentheses")
    call IO_Assert(io, iMode /= SEARCH_QUOTE, "Unbalanced quotation marks")
    if (iMode /= SEARCH_NONBLANK) then
      io%iFieldCount = io%iFieldCount + 1
      io%sFields(io%iFieldCount) = sRecord(iFieldStart:i-1)
    end if
    io%iThisField = 1

    return
  end subroutine IO_SplitRecord


  subroutine IO_SwitchFields(io)
    !! Switches the order of the first two I/O fields (for backwards compatibility/version control with GMS)
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ LOCAL ]
    character(len=255) :: sField

    call IO_Assert(io, io%iThisField < io%iFieldCount, "Attempt to switch less than two fields")
    sField = io%sFields(io%iThisField)
    io%sFields(io%iThisField) = io%sFields(io%iThisField+1)
    io%sFields(io%iThisField+1) = sField

    return
  end subroutine IO_SwitchFields


  function sIO_GetField(io, tag, def, allowed, force_uppercase) result(sField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a text string
    !! 'tag' is a text tag for the field(for testing and debugging)
    !! 'def' is a default value for the field; if omitted, field must be non-blank
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    character(len=*), intent(in), optional :: def
    character(len=*), intent(in), dimension(:), optional :: allowed
    logical, intent(in), optional :: force_uppercase
    ! [ RETURN VALUE ]
    character(len=255) :: sField

    if (io%iThisField > io%iFieldCount) then
      sField = repeat(' ', 255)
    else
      sField = io%sFields(io%iThisField)
      if (present(force_uppercase) ) then
        if (force_uppercase) sField = uppercase (sField)
      end if
      io%iThisField = io%iThisField+1
    end if

    if (len_trim(sField) == 0) then
      if (present(def)) then
        sField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    end if

    if (present(allowed)) then
      call IO_Assert(io, any(allowed(:)(1:len_trim(sField))==trim(sField)), "Invalid field entry " // trim(sField))
    end if

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "[", trim(tag), " = '", trim(sField), "']"
      call IO_MessageText(io)
    end if

    return
  end function sIO_GetField


  function cIO_GetComplex(io, tag, def) result(cField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a complex value
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    complex(kind=AE_REAL), intent(in), optional :: def
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cField
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (io%iThisField > io%iFieldCount .or. &
        len_trim(io%sFields(io%iThisField)) == 0) then
      if (present(def)) then
        cField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    else
      read (unit=io%sFields(io%iThisField), fmt=*, iostat=istat) cField
      call IO_Assert(io, istat == 0, "Complex read failed: " // io%sFields(io%iThisField))
    end if
    io%iThisField = io%iThisField+1

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "  [", trim(tag), " = ", cField, "]"
      call IO_MessageText(io)
    end if

    return
  end function cIO_GetComplex


  function cIO_GetCoordinate(io, tag, def, extents, check_points, check_tol) result(cField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a complex coordinate
    !! value, making use of the coordinate rotation/scaling options.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    complex(kind=AE_REAL), intent(in), optional :: def
    logical, intent(in), optional :: extents
    complex(kind=AE_REAL), dimension(:), intent(in), optional :: check_points
    real(kind=AE_REAL), intent(in), optional :: check_tol
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cField
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (io%iThisField > io%iFieldCount .or. &
        len_trim(io%sFields(io%iThisField)) == 0) then
      if (present(def)) then
        cField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    else
      read (unit=io%sFields(io%iThisField), fmt=*, iostat=istat) cField
      call IO_Assert(io, istat == 0, "Coordinate read failed: " // io%sFields(io%iThisField))
    end if
    io%iThisField = io%iThisField+1

    cField = cIO_LocalCoords(io, cField)
    if (present(extents)) then
      if (extents) then
        io%rXMin = min(io%rXMin, real(cField, AE_REAL))
        io%rXMax = max(io%rXMin, real(cField, AE_REAL))
        io%rYMin = min(io%rYMin, aimag(cField))
        io%rYMax = max(io%rYMin, aimag(cField))
      end if
    end if

    if (present(check_tol)) then
      if (present(check_points)) then
        call IO_Assert(io, .not. any(abs(cField-check_points) < check_tol), "Found coincident points")
      end if
    else
      if (present(check_points)) then
        call IO_Assert(io, .not. any(abs(cField-check_points) == rZERO), "Found coincident points")
      end if
    end if

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "  [", trim(tag), " = coord", cField, "]"
      call IO_MessageText(io)
    end if

    return
  end function cIO_GetCoordinate


  function rIO_GetReal(io, tag, def, minimum, maximum) result(rField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a real value
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    real(kind=AE_REAL), intent(in), optional :: def
    real(kind=AE_REAL), intent(in), optional :: minimum
    real(kind=AE_REAL), intent(in), optional :: maximum
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rField
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (io%iThisField > io%iFieldCount .or. &
        len_trim(io%sFields(io%iThisField)) == 0) then
      if (present(def)) then
        rField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    else
      read (unit=io%sFields(io%iThisField), fmt=*, iostat=istat) rField
      call IO_Assert(io, istat == 0, "Real read failed: " // io%sFields(io%iThisField))
    end if
    io%iThisField = io%iThisField+1

    if (present(minimum) .and. rField < minimum) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", rField, ") below minimum of ", minimum
      call IO_Assert(io, .false., "")
    end if

    if (present(maximum) .and. rField > maximum) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", rField, ") above maximum of ", maximum
      call IO_Assert(io, .false., "")
    end if

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "  [", trim(tag), " = ", rField, "]"
      call IO_MessageText(io)
    end if

    return
  end function rIO_GetReal


  function iIO_GetInteger(io, tag, def, minimum, maximum, allowed, forbidden) result(iField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a real value
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    integer(kind=AE_INT), intent(in), optional :: def
    integer(kind=AE_INT), intent(in), optional :: minimum
    integer(kind=AE_INT), intent(in), optional :: maximum
    integer(kind=AE_INT), intent(in), dimension(:), optional :: allowed
    integer(kind=AE_INT), intent(in), dimension(:), optional :: forbidden
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iField
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (io%iThisField > io%iFieldCount .or. &
        len_trim(io%sFields(io%iThisField)) == 0) then
      if (present(def)) then
        iField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    else
      read (unit=io%sFields(io%iThisField), fmt=*, iostat=istat) iField
      call IO_Assert(io, istat == 0, "Integer read failed: " // io%sFields(io%iThisField))
    end if
    io%iThisField = io%iThisField+1

    if (present(minimum) .and. iField < minimum) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", iField, ") below minimum of ", minimum
      call IO_Assert(io, .false., "")
    end if

    if (present(maximum) .and. iField > maximum) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", iField, ") above maximum of ", maximum
      call IO_Assert(io, .false., "")
    end if

    if (present(allowed)) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", iField, ") is not allowed"
      call IO_Assert(io, any(allowed==iField), "")
    end if

    if (present(forbidden)) then
      write (unit=IO_MessageBuffer, fmt=*) "Value for ", tag, " (", iField, ") is forbidden"
      call IO_Assert(io, .not. any(forbidden==iField), "")
    end if

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "  [", trim(tag), " = ", iField, "]"
      call IO_MessageText(io)
    end if

    return
  end function iIO_GetInteger


  function lIO_GetLogical(io, tag, def) result(lField)
    !! Returns the 'next' field, as parsed by IO_SplitRecord, as a real value
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in) :: tag
    logical, intent(in), optional :: def
    ! [ RETURN VALUE ]
    logical :: lField
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (io%iThisField > io%iFieldCount .or. &
        len_trim(io%sFields(io%iThisField)) == 0) then
      if (present(def)) then
        lField = def
      else
        call IO_Assert(io, .false., "Required field " // trim(tag) // " was not found")
      end if
    else
      read (unit=io%sFields(io%iThisField), fmt=*, iostat=istat) lField
      call IO_Assert(io, istat == 0, "Logical read failed: " // io%sFields(io%iThisField))
    end if
    io%iThisField = io%iThisField+1

    if (io%lDebug) then
      write (unit=IO_MessageBuffer, fmt=*) "  [", trim(tag), " = ", lField, "]"
      call IO_MessageText(io)
    end if

    return
  end function lIO_GetLogical


  function cIO_LocalCoords(io, cWZ) result(cLZ)
    !! Projects the world coordinate cWZ into the local coordinate system
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cWZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cLZ
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp

    ! First, translate and scale, then rotate
    rTemp = (/real(cWZ, AE_REAL), aimag(cWZ), rONE/)
    rTemp = matmul(io%rWtoL, rTemp)
    cLZ = cmplx(rTemp(1), rTemp(2), AE_REAL)

    ! Now, apply the anisotropy if in profile mode
    if ( io%lProfile ) then
      cLZ = cmplx(real(cLZ) / io%rSqrtA, aimag(cLZ), AE_REAL)
    end if

    return
  end function cIO_LocalCoords


  function cIO_WorldCoords(io, cLZ) result(cWZ)
    !! Projects the local coordinate cLZ into the world coordinate system
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cLZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cWZ
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp

    ! First, apply the anisotropy if in profile mode
    if ( io%lProfile ) then
      cWZ = cmplx(real(cLZ) * io%rSqrtA, aimag(cLZ), AE_REAL)
    else
      cWZ = cLZ
    end if

    ! Now, translate and scale, then rotate
    rTemp = (/real(cWZ, AE_REAL), aimag(cWZ), rONE/)
    rTemp = matmul(io%rLtoW, rTemp)
    cWZ = cmplx(rTemp(1), rTemp(2), AE_REAL)

    return
  end function cIO_WorldCoords


  function cIO_LocalDischarge(io, cWZ) result(cLZ)
    !! Projects the world discharge vector cWZ into the local coordinate system
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cWZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cLZ
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp

    ! First, translate and scale, then rotate
    rTemp = (/real(cWZ, AE_REAL), aimag(cWZ), rONE/)
    rTemp = matmul(io%rWtoL, rTemp)
    cLZ = cmplx(rTemp(1), rTemp(2), AE_REAL)

    ! Now, apply the anisotropy if in profile mode
    if ( io%lProfile ) then
      cLZ = cmplx(real(cLZ) * io%rSqrtA, aimag(cLZ), AE_REAL)
    end if

    return
  end function cIO_LocalDischarge


  function cIO_WorldDischarge(io, cLZ) result(cWZ)
    !! Projects the local discharge vector cLZ into the world coordinate system
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cLZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cWZ
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp

    ! First, apply the anisotropy if in profile mode
    if ( io%lProfile ) then
      cWZ = cmplx(real(cLZ) / io%rSqrtA, aimag(cLZ), AE_REAL)
    else
      cWZ = cLZ
    end if

    ! Now, translate and scale, then rotate
    rTemp = (/real(cWZ, AE_REAL), aimag(cWZ), rONE/)
    rTemp = matmul(io%rLtoW, rTemp)
    cWZ = cmplx(rTemp(1), rTemp(2), AE_REAL)

    return
  end function cIO_WorldDischarge


  function rIO_LocalLength(io, rWSigma, cOrientation) result(rLSigma)
    !! Converts given sink density for a line-sink or flux across a line-segment
    !! to the world coordinate system
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    real(kind=AE_REAL), intent(in) :: rWSigma
    complex(kind=AE_REAL), intent(in) :: cOrientation
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rLSigma
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp
    real(kind=AE_REAL) :: rOrientation

    if ( io%lProfile ) then
      ! Account for scaling and anisotropy
      rOrientation = atan2(real(cOrientation), aimag(cOrientation))
      rLSigma = io%rScaling * abs(cmplx(rWSigma*cos(rOrientation), &
                                        rWSigma*sin(rOrientation)*io%rSqrtA, AE_REAL))
    else
      ! Account for scaling only
      rLSigma = rWSigma / io%rScaling
    end if

    return
  end function rIO_LocalLength

  function rIO_WorldLength(io, rLSigma, cOrientation) result(rWSigma)
    !! Converts given sink density for a line-sink or flux across a line-segment
    !! to the world coordinate system
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    real(kind=AE_REAL), intent(in) :: rLSigma
    complex(kind=AE_REAL), intent(in) :: cOrientation
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rWSigma
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(3) :: rTemp
    real(kind=AE_REAL) :: rOrientation

    if ( io%lProfile ) then
      ! Account for scaling and anisotropy
      rOrientation = atan2(real(cOrientation), aimag(cOrientation))
      rWSigma = abs(cmplx(rLSigma*cos(rOrientation), &
                          rLSigma*sin(rOrientation)/io%rSqrtA, AE_REAL)) / io%rScaling
    else
      ! Account for scaling only
      rWSigma = rLSigma / io%rScaling
    end if

    return
  end function rIO_WorldLength


  subroutine IO_ErrorText(io, sMessage)
    ! Writes a line of text to the error LU
    character(len=*), intent(in), optional :: sMessage
    type(IO_STATUS), pointer :: io

    if (present(sMessage)) then
      if (io%lFilesOpen) then
        write (unit=LU_ERROR, fmt="(1x, a)") trim(sMessage)
      else
        write (unit=*, fmt="(1x, a)") trim(sMessage)
      end if
    else
      if (io%lFilesOpen) then
        write (unit=LU_ERROR, fmt="(1x, a)") trim(IO_MessageBuffer)
      else
        write (unit=*, fmt="(1x, a)") trim(IO_MessageBuffer)
      end if
    end if
    return
  end subroutine IO_ErrorText


  subroutine IO_MessageText(io, sMessage)
    ! Writes a run-time message to STDOUT. If the message is missing, uses IO_MessageBuffer
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    character(len=*), intent(in), optional :: sMessage

    ! Send the message...
    if (.not. present(sMessage)) then
      write (unit=*, fmt=*) trim(IO_MessageBuffer)
      if (io%lFilesOpen) then
        call IO_ErrorText(io, "      **  " // trim(IO_MessageBuffer))
      end if
    else if (trim(sMessage) == "\n") then
      write (unit=*, fmt=*)
    else
      write (unit=*, fmt=*) trim(sMessage)
      if (io%lFilesOpen) then
        call IO_ErrorText(io, "      **  " // trim(sMessage))
      end if
    end if

  end subroutine IO_MessageText


  subroutine IO_Assert(io, lTest, sText)
    ! Asserts that lTest is .true.  If .not. lTest, then print the message
    ! sText and halt execution.  Used for debug error testing.
    logical :: lTest
    character(len=*), intent(in) :: sText
    type(IO_STATUS), pointer :: io
    ! this is to force a traceback
    type(IO_STATUS), pointer :: dummy

    if (.not. lTest) then
      if (len_trim(sText) == 0) then
        call IO_ErrorText(io, "Assertion failed: " // IO_MessageBuffer)
      else
        call IO_ErrorText(io, "Assertion failed: " // sText)
      end if

      call IO_CloseAll(io)
      if (io%lCrash) then
        ! This forces a stack dump with most compilers -- aids in debugging
        print *, "Assertion failed: ", dummy
      end if
      call exit(1)
    end if

    return
  end subroutine IO_Assert

  ! The following are some convenience routines for writing HTML output on the output file
  ! This is a new feature for ModAEM 1.4+. All routines write directly to the output LU. This
  ! is by no means a complete implementation of HTML output. In fact, all output is written
  ! to tables, with the exception of the header information.


  subroutine HTML_Start()
    ! Writes the header to the HTML output file
    write (unit=LU_OUTPUT, &
           fmt='("<html><head></head><body align=""center"">")' &
           )
    return
  end subroutine HTML_Start


  subroutine HTML_End()
    ! Writes the trailer to the HTML output file
    call HTML_Header('End of output', 1)
    write (unit=LU_OUTPUT, &
           fmt='("</body></html>")' &
           )
    return
  end subroutine HTML_End


  subroutine HTML_HR()
    ! Writes a horizontal rule to the HTML output file
    write (unit=LU_OUTPUT, &
           fmt='("<hr>")' &
           )
    return
  end subroutine HTML_HR


  subroutine HTML_Header(sText, iLevel)
    ! Writes a header < H? > where ? is the optional parameter iLevel
    character(len=*), intent(in) :: sText
    integer(kind=AE_INT), intent(in), optional :: iLevel
    integer(kind=AE_INT) :: iMyLevel

    if (present(iLevel)) then
      iMyLevel = iLevel
    else
      iMyLevel = 1
    end if

    ! Put a rule if it's a level-1 header
    if (iMyLevel == 1) then
      call HTML_HR()
    end if

    write (unit=LU_OUTPUT, &
           fmt='("<h", i1.1, "> ", a, " </h", i1.1, ">")' &
           ) iMyLevel, sText, iMyLevel
    return
  end subroutine HTML_Header


  subroutine HTML_StartTable()
    ! Begins a table
    write (unit=LU_OUTPUT, &
           fmt='("<table border=""1"" cellspacing=""0"" cellpadding=""2"" bordercolor=""grey"" ' // &
           'bordercolorlight=""grey"" bordercolordark=""grey"" align=""center"">")' &
           )
    return
  end subroutine HTML_StartTable


  subroutine HTML_EndTable()
    ! Ends a table
    write (unit=LU_OUTPUT, &
           fmt='("</table>")' &
           )
    return
  end subroutine HTML_EndTable


  subroutine HTML_StartRow()
    ! Begins a table
    write (unit=LU_OUTPUT, &
           fmt='("<tr>")' &
           )
    return
  end subroutine HTML_StartRow


  subroutine HTML_EndRow()
    ! Ends a table
    write (unit=LU_OUTPUT, &
           fmt='("</tr>")' &
           )
    return
  end subroutine HTML_EndRow


  subroutine HTML_TableHeader(sValues, colspan)
    ! Writes text fields into a table(optionally with 'colspan')
    character(len=*), dimension(:), intent(in) :: sValues
    integer(kind=AE_INT), intent(in), optional :: colspan
    !
    integer(kind=AE_INT) :: i

    do i = 1, size(sValues)
      if (present(colspan)) then
        write (unit=LU_OUTPUT, &
               fmt='("<th colspan=""", i1.1, """>", a, "</th>")' &
               ) colspan, trim(sValues(i))
      else
        write (unit=LU_OUTPUT, &
               fmt='("<th>", a, "</th>")' &
               ) trim(sValues(i))
      end if
    end do

    return
  end subroutine HTML_TableHeader


  subroutine HTML_ColumnText(sValues, colspan)
    ! Writes text fields into a table(optionally with 'colspan')
    character(len=*), dimension(:), intent(in) :: sValues
    integer(kind=AE_INT), intent(in), optional :: colspan
    !
    integer(kind=AE_INT) :: i

    do i = 1, size(sValues)
      if (present(colspan)) then
        write (unit=LU_OUTPUT, &
               fmt='("<td colspan=""", i1.1, """>", a, "</td> ")' &
               ) colspan, trim(sValues(i))
      else
        write (unit=LU_OUTPUT, &
               fmt='("<td align=""center"">", a, "</td>")' &
               ) trim(sValues(i))
      end if
    end do

    return
  end subroutine HTML_ColumnText


  subroutine HTML_ColumnInteger(iValues, sFmt)
    ! Writes(optionally formatted) integers in columns
    integer(kind=AE_INT), dimension(:), intent(in) :: iValues
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt
    integer(kind=AE_INT) :: i

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'i12'
    end if

    do i = 1, size(iValues)
      write (unit=LU_OUTPUT, &
             fmt='("<td align=""center"">", ' // sMyFmt // ', "</td>")' &
             ) iValues(i)
    end do

    return
  end subroutine HTML_ColumnInteger


  subroutine HTML_ColumnLogical(lValues, sFmt)
    ! Writes(optionally formatted) integers in columns
    logical, dimension(:), intent(in) :: lValues
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt
    integer(kind=AE_INT) :: i

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'l1'
    end if

    do i = 1, size(lValues)
      write (unit=LU_OUTPUT, &
             fmt='("<td align=""center"">", ' // sMyFmt // ', "</td>")' &
             ) lValues(i)
    end do

    return
  end subroutine HTML_ColumnLogical


  subroutine HTML_ColumnReal(rValues, sFmt)
    ! Writes an(optionally formatted) real in a column
    real(kind=AE_REAL), dimension(:), intent(in) :: rValues
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt
    integer(kind=AE_INT) :: i

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'g12.5'
    end if

    do i = 1, size(rValues)
      write (unit=LU_OUTPUT, &
             fmt='("<td align=""center"">", ' // sMyFmt // ', "</td>")' &
             ) rValues(i)
    end do

    return
  end subroutine HTML_ColumnReal


  subroutine HTML_ColumnComplex(cValues, sFmt)
    ! Writes an(optionally formatted) complex in two consecutive columns
    complex(kind=AE_REAL), dimension(:), intent(in) :: cValues
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt
    integer(kind=AE_INT) :: i

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'g12.5'
    end if

    do i = 1, size(cValues)
      call HTML_ColumnReal((/real(cValues(i))/), sMyFmt)
      call HTML_ColumnReal((/aimag(cValues(i))/), sMyFmt)
    end do

    return
  end subroutine HTML_ColumnComplex


  subroutine HTML_AttrInteger(sTitle, iValue, sFmt)
    ! Writes an(optionally formatted) integer "attribute row" with title
    character(len=*), intent(in) :: sTitle
    integer(kind=AE_INT), intent(in) :: iValue
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'i12'
    end if

    call HTML_StartRow()
    call HTML_ColumnText((/sTitle/))
    call HTML_ColumnInteger((/iValue/), sMyFmt)
    call HTML_EndRow()

    return
  end subroutine HTML_AttrInteger


  subroutine HTML_AttrLogical(sTitle, lValue, sFmt)
    ! Writes an(optionally formatted) integer "attribute row" with title
    character(len=*), intent(in) :: sTitle
    logical, intent(in) :: lValue
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'l1'
    end if

    call HTML_StartRow()
    call HTML_ColumnText((/sTitle/))
    call HTML_ColumnLogical((/lValue/), sMyFmt)
    call HTML_EndRow()

    return
  end subroutine HTML_AttrLogical


  subroutine HTML_AttrReal(sTitle, rValue, sFmt)
    ! Writes an(optionally formatted) real "attribute row" with title
    character(len=*), intent(in) :: sTitle
    real(kind=AE_REAL), intent(in) :: rValue
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'g12.5'
    end if

    call HTML_StartRow()
    call HTML_ColumnText((/sTitle/))
    call HTML_ColumnReal((/rValue/), sMyFmt)
    call HTML_EndRow()

    return
  end subroutine HTML_AttrReal


  subroutine HTML_AttrComplex(sTitle, cValue, sFmt)
    ! Writes an(optionally formatted) complex "attribute row" with title
    character(len=*), intent(in) :: sTitle
    complex(kind=AE_REAL), intent(in) :: cValue
    character(len=*), intent(in), optional :: sFmt
    ! [ LOCALS ]
    character(len=32) :: sMyFmt

    if (present(sFmt)) then
      sMyFmt = sFmt
    else
      sMyFmt = 'g12.5'
    end if

    call HTML_StartRow()
    call HTML_ColumnText((/sTitle/))
    call HTML_ColumnReal((/real(cValue)/), sMyFmt)
    call HTML_ColumnReal((/aimag(cValue)/), sMyFmt)
    call HTML_EndRow()

    return
  end subroutine HTML_AttrComplex


  subroutine IO_GetWindow(io, cLL, cUR)
    !! Returns the corners of the saved element extents, subject to a minimum window size
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    complex(kind=AE_REAL), intent(out) :: cLL
    complex(kind=AE_REAL), intent(out) :: cUR

    cLL = cmplx(io%rXMin, io%rYMin)
    cUR = cmplx(io%rXMax, io%rYMax)

    return
  end subroutine IO_GetWindow

end module u_io
#ifdef __UNITTEST__

program unittest

  use u_constants
  use s_unittest
  use u_io

  implicit none

  ! Main program -- checks the character buffer functions in module u_io
  character(len=255) :: sBuf
  character(len=255), dimension(15) :: sTest
  logical :: lFail
  real(kind=AE_REAL) :: rV
  real(kind=AE_REAL), dimension(3, 3) :: rCheck
  complex(kind=AE_REAL) :: cV, cTest
  integer(kind=AE_INT) :: i, iV
  logical :: lV
  character(len=255) :: sV
  type(IO_STATUS), pointer :: io

  lFail = .false.
  sBuf = "TEST(100.0, -100.0) (1.0, 1.0) 1.0 -1.0 1.0 0.0 1 -1 1 0 T ""THIS IS A TEST"" HELLO"
  sTest = (/"TEST", &
          "(100.0, -100.0)", &
          "(1.0, 1.0)", &
          "1.0", &
          "-1.0", &
          "1.0", &
          "0.0", &
          "1", &
          "-1", &
          "1", &
          "0", &
          "1", &
          "T", &
          "THIS IS A TEST", &
          "HELLO"/)

  print *, "** PERFORMING UNIT TESTS FOR MODULE u_io **"
  print *

  print *, "Testing IO_Create"
  io => IO_Create(cmplx(-1.0, 1.0, AE_REAL), 0.1_AE_REAL, 30.0_AE_REAL, .true., .false.)
  call complexValueCheck(io, "Origin: ", io%cZOrigin, cmplx(-1.0, 1.0, AE_REAL), lFail)
  call realValueCheck(io, "Scaling: ", io%rScaling, 0.1_AE_REAL, lFail)
  call realValueCheck(io, "Rotation: ", io%rRotation, 30.0_AE_REAL, lFail)

  print *, "Testing IO_SplitRecord"
  call IO_SplitRecord(io, sBuf)
  call integerValueCheck(io, "Number of fields: ", io%iFieldCount, 8, lFail)
  call characterArrayCheck(io, "Field contents: ", io%sFields, sTest, 8, lFail)

  print *, "Testing individual fields"
  call characterValueCheck(io, "Field check: ", sIO_GetField(io, "CHARTEST"), sTest(1), lFail)
  cTest = cmplx(100.0_AE_REAL, -100.0_AE_REAL, AE_REAL)
  call complexValueCheck(io, "Complex check:", cIO_GetComplex(io, "CPLXTEST"), cTest, lFail)
  cTest = cmplx(1.03660254_AE_REAL, -0.86339746_AE_REAL, AE_REAL)
  cV = cIO_GetCoordinate(io, "COORDTEST")
  call complexValueCheck(io, "Coordinate check:", cV, cTest, lFail)
  cTest = cmplx(1.0, 1.0, AE_REAL)
  cV = cIO_WorldCoords(cV)
  call complexValueCheck(io, "Reverse coordinate check:", cV, cTest, lFail)

  call realValueCheck(io, "Real check:", rIO_GetReal(io, "REALTEST"), rONE, lFail)
  call realValueCheck(io, "Real check(low limit):", rIO_GetReal(io, "REALTEST", minimum = rZERO), rONE, lFail)
  call realValueCheck(io, "Real check(high limit):", rIO_GetReal(io, "REALTEST", maximum = rZERO), -rONE, lFail)
  call realValueCheck(io, "Real check(both limits):", rIO_GetReal(io, "REALTEST", minimum = -rONE, maximum = rONE), rZERO, lFail)

  call integerValueCheck(io, "Int check:", iIO_GetInteger(io, "INTTEST"), 1, lFail)
  call integerValueCheck(io, "Int check(low limit):", iIO_GetInteger(io, "INTTEST", minimum = 0), 1, lFail)
  call integerValueCheck(io, "Int check(high limit):", iIO_GetInteger(io, "INTTEST", maximum = 0), -1, lFail)
  call integerValueCheck(io, "Int check(both limits):", iIO_GetInteger(io, "INTTEST", minimum = -1, maximum = 1), 0, lFail)

  if (lFail) then
    print *, "One or more unit tests failed"
  else
    print *, "*** ALL TESTS SUCCESSFUL ***"
  end if

  ! Try to fire an exception!
  print *, "Testing missing field -- this should fail!"
  call characterValueCheck(io, "Field check(w/o default): ", sIO_GetField(io, "DEFAULTTEST"), "DEFAULT", lFail)

  stop
end program unittest
#endif
