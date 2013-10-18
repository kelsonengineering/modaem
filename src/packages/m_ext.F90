module m_ext

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

  ! Data extraction/inquiry/check module
  ! V.A.Kelson  10/23/98
  !
  ! This module allows for the extraction of results from a ModAEM solution
  !
  ! Commands:
  !
  !   FIL < filename > Output file. No extension will be appended.
  !
  !   HEA  (X, Y) L label            Extracts the head at the point(X, Y) in layer L
  !
  !   DIS  (X, Y) L label            Extracts the discharge at the point(X, Y) in layer L
  !
  !   VEL  (X, Y) L label            Extracts the velocity at the point(X, Y) in layer L
  !
  !   FLO  (X1, Y1) (X2, Y2) L label  Extracts the integrated flow across the specified path
  !
  !   WL0  L                        Writes check information for the specified module
  !   WL1  L
  !   LS0  L
  !   LS1  L
  !   HB0  L
  !
  ! For all extraction commands which requre a label, the label is returned in the extraction file.
  ! All data are written to space-delimeted files with constant column widths.  A header record is
  ! prepended to the file.

  use u_constants
  use u_io
  use m_wl0
  use m_wl1
  use m_ls0
  use m_ls1
  use m_hb0
  use m_aem

  implicit none

  public

  ! Global variables for this module


contains


  subroutine EXT_Read(io, aem)
    ! This reads and processes the input commands for the EXT package

    ! Argument list
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io

    ! Locals -- for Directive parsing
    type(DIRECTIVE), dimension(10), parameter :: dirDirectives = &
                       (/dirINP, dirOUT, dirOPT, dirEND, dirDIM, &
                       dirHEA, dirDIS, dirPOT, dirVEL, dirFLO/)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode         ! Opcode
    integer(kind=AE_INT) :: iStat           ! RunTime error code
    ! Placeholders for tracing parameters
    real(kind=AE_REAL) :: rX, rY, rZ, rX2, rY2
    integer(kind=AE_INT) :: iInputDimension
    logical :: lFlag, lInputMode, lOutputMode, lOutputOpen
    character(len=132) :: sFile, sInputFile, sMode

    call IO_MessageText(io, "Reading EXT module input")

    ! Clear the status flags
    lOutputOpen = .false.
    iInputDimension = 2
    call IO_MessageText(io, "Entering inquiry module EXT")
    call IO_Assert(io, (associated(aem)), "EXT_Read: No AEM_DOMAIN object")

    ! The remainder of this routine uses ifIOInputRecord to process
    ! the model input file.
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation.
          call IO_Assert(io, .false., "EXT_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModGRI "ifXXXRead" routines.
          call IO_Assert(io, .false., "EXT_Read: Unexpected EOF")
        case (kOpEND)
          ! END mark was found. Exit the file parser.
          if (lOutputOpen) then
            close(unit=LU_EXT_OUTPUT)
          end if
          exit
        case (kOpINP)
          !****************************************************************************
          ! Here for the INP command -- specify the input coordinate file
          !****************************************************************************
          sMode = sIO_GetField(io, "sMode", allowed=(/"ASC", "BIN"/), force_uppercase=.true.)
          sFile = sIO_GetField(io, "sFile")
          ! Save the input mode and file name for opening later
          if (sMode == 'ASC') then
            lInputMode = .false.
          else
            lInputMode = .true.
          end if
          sInputFile = sFile
        case (kOpOUT)
          !****************************************************************************
          ! Here for the OUT command -- specify the output file
          !****************************************************************************
          sMode = sIO_GetField(io, "sMode", allowed=(/"ASC", "BIN"/), force_uppercase=.true.)
          sFile = sIO_GetField(io, "sFile")
          if (lOutputOpen) then
            close(unit=LU_EXT_OUTPUT)
          end if
          ! Open the file in the appropriate mode
          if (uppercase (sMode) == 'ASC') then
            lOutputMode = .false.
            open(unit=LU_EXT_OUTPUT, file = sFile, iostat=iStat)
            call IO_Assert(io, (iStat == 0), 'Cannot open ASCII output file '//sFile)
            lOutputOpen = .true.
          else
            lOutputMode = .true.
            open(unit=LU_EXT_OUTPUT, file = sFile, access='STREAM', iostat=iStat)
            call IO_Assert(io, (iStat == 0), "Cannot open BINARY output file "//sFile)
            lOutputOpen = .true.
          end if
        case (kOpDIM)
          !****************************************************************************
          ! Here for the DIM command -- specify the input coordinate dimension
          !****************************************************************************
          iInputDimension = iIO_GetInteger(io, "iInputDimension", allowed=(/2, 3/), def=2)
        case (kOpHEA, kOpDIS, kOpPOT, kOpVEL, kOpFLO)
          print *, 'EXTRACTING ', iOpCode
          !****************************************************************************
          ! Here for the extraction commands
          !****************************************************************************
          ! Open the input file
          if (lInputMode) then
            open(unit=LU_EXT_INPUT, file = sInputFile, access='STREAM', iostat=iStat)
            call IO_Assert(io, (iStat == 0), "Cannot open input file " // sInputFile)
          else
            open(unit=LU_EXT_INPUT, file = sInputFile, iostat=iStat)
            call IO_Assert(io, (iStat == 0), "Cannot open input file " // sInputFile)
          end if
          ! Process the input file
          do
            ! Fetch the coordinate(s)
            if (lInputMode) then
              if (iInputDimension == 2) then
                read (LU_EXT_INPUT, iostat=iStat) rX, rY
                if (iStat /= 0) exit
                if (iOpCode == kOpFLO) then
                  read (LU_EXT_INPUT, iostat=iStat) rX2, rY2
                  if (iStat /= 0) exit
                end if
              else
                read (LU_EXT_INPUT, iostat=iStat) rX, rY, rZ
                if (iStat /= 0)  exit
              end if
            else
              if (iInputDimension == 2) then
                read (LU_EXT_INPUT, fmt=*, iostat=iStat) rX, rY
                if (iStat /= 0) exit
                if (iOpCode == kOpFLO) then
                  read (LU_EXT_INPUT, fmt=*, iostat=iStat) rX2, rY2
                  if (iStat /= 0) exit
                end if
              else
                read (LU_EXT_INPUT, fmt=*, iostat=iStat) rX, rY, rZ
                if (iStat /= 0) exit
              end if
            end if
            ! Write the output for this point
            if (iOpCode == kOpHEA .and. iInputDimension == 2) then
              if (lOutputMode) then
                write (unit=LU_EXT_OUTPUT) rAEM_Head(io, aem, cmplx(rX, rY, AE_REAL))
              else
                write (unit=LU_EXT_OUTPUT, fmt='(e14.7)') rAEM_Head(io, aem, cmplx(rX, rY, AE_REAL))
              end if
            else if (iOpCode == kOpPOT .and. iInputDimension == 2) then
              if (lOutputMode) then
                write (unit=LU_EXT_OUTPUT) cAEM_Potential(io, aem, cmplx(rX, rY, AE_REAL))
              else
                write (unit=LU_EXT_OUTPUT, fmt='(e14.7, 1x, e14.7)') cAEM_Potential(io, aem, cmplx(rX, rY, AE_REAL))
              end if
            else if (iOpCode == kOpDIS .and. iInputDimension == 2) then
              if (lOutputMode) then
                write (unit=LU_EXT_OUTPUT) cAEM_Discharge(io, aem, cmplx(rX, rY, AE_REAL))
              else
                write (unit=LU_EXT_OUTPUT, fmt='(e14.7, 1x, e14.7)') cAEM_Discharge(io, aem, cmplx(rX, rY, AE_REAL))
              end if
            else if (iOpCode == kOpVEL .and. iInputDimension == 2) then
              if (lOutputMode) then
                write (unit=LU_EXT_OUTPUT) cAEM_Velocity(io, aem, cmplx(rX, rY, AE_REAL))
              else
                write (unit=LU_EXT_OUTPUT, fmt='(e14.7, 1x, e14.7)') cAEM_Velocity(io, aem, cmplx(rX, rY, AE_REAL))
              end if
            else if (iOpCode == kOpFLO .and. iInputDimension == 2) then
              if (lOutputMode) then
                write (unit=LU_EXT_OUTPUT) rAEM_Flow(io, aem, (/cmplx(rX, rY, AE_REAL), cmplx(rX2, rY2, AE_REAL)/))
              else
                write (unit=LU_EXT_OUTPUT, fmt='(e14.7, 1x)') rAEM_Flow(io, aem, &
                                                              (/ cmplx(rX, rY, AE_REAL), &
                                                                 cmplx(rX2, rY2, AE_REAL) /))
              end if
            else
              call IO_Assert(io, .false., "Illegal extraction mode")
            end if
          end do
          close(unit=LU_EXT_INPUT)
      end select
    end do

    call IO_MessageText(io, "Leaving EXT module")

    return
  end subroutine EXT_Read

end module m_ext
