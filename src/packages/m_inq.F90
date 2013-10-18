module m_inq

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
  !   HEA  (X, Y) label              Extracts the head at the point(X, Y)
  !
  !   DIS  (X, Y) label            Extracts the discharge at the point(X, Y)
  !
  !   VEL  (X, Y) label            Extracts the velocity at the point(X, Y)
  !
  !   FLO  (X1, Y1) (X2, Y2) label  Extracts the integrated flow across the specified path
  !
  !   BDY  (X1, Y1) (X2, Y2) label  Extracts head and flux information across the
  !                               specified path
  !
  !   WL0                         Writes check information for the specified module
  !   WL1
  !   LS0
  !   LS1
  !   LS2
  !   HB0
  !   AS0
  !   AQU
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

  ! Flags
  logical, private, save :: fINQFileOpen = .false.


contains


  subroutine INQ_Read(io, aem)
    ! This reads and processes the input commands for the INQ package

    ! Argument list
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io

    ! Locals -- for Directive parsing

#ifndef __GPL__
    type(DIRECTIVE), dimension(20), parameter :: dirDirectives = &
#else
                       type(DIRECTIVE), dimension(19), parameter :: dirDirectives = &
#endif
                       (/dirFIL, &
                       dirHEA, dirDIS, dirPOT, dirVEL, dirFLO, dirRCH, &
                       dirSAT, dirBDY, dirGAG, &
                       dirWL0, dirWL1, dirLS0, &  !**pd
    dirLS1, dirLS2, dirHB0, &
              dirAS0, &
#ifndef __GPL__
              dirCW0, &
#endif
              dirAQU, dirEND/)

    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode         ! Opcode
    integer(kind=AE_INT) :: iStat           ! RunTime error code
    ! Placeholders for tracing parameters
    real(kind=AE_REAL) :: rValue, rHead, rFlux
    complex(kind=AE_REAL) :: cZ1, cZ2, cValue
    integer(kind=AE_INT) :: iLabel, iID, iGageId
    logical :: lFlag
    character(len=132) :: sFile

    call IO_MessageText(io, "Reading INQ module input")

    ! Clear the status flags
    fINQFileOpen = .false.
    call IO_MessageText(io, "Entering inquiry module INQ")
    call IO_Assert(io, (associated(aem)), "INQ_Read: No AEM_DOMAIN object")

    ! The remainder of this routine uses ifIOInputRecord to process
    ! the model input file.
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation.
          call IO_Assert(io, .false., "INQ_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModGRI "ifXXXRead" routines.
          call IO_Assert(io, .false., "INQ_Read: Unexpected EOF")
        case (kOpEND)
          ! END mark was found. Exit the file parser.
          if (fINQFileOpen) then
            close(unit=LU_INQ)
          end if
          exit
        case (kOpFIL)
          !****************************************************************************
          ! Here for the FIL command -- open the output file
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          ! If the output file is open, close it now.
          if (fINQFileOpen) then
            close(unit=LU_INQ)
            fINQFileOpen = .false.
          end if
          ! Open the new output file.  If the open fails, write a message.
          sFile = trim(sFile)
          open(unit=LU_INQ, file = sFile, action = "WRITE", status="REPLACE", iostat=iStat)
          call IO_Assert(io, (iStat == 0), "INQ_Read: Open failed")
          fINQFileOpen = .true.
        case (kOpHEA)
          !****************************************************************************
          ! Here for the HEA command -- extract the head
          !****************************************************************************
          ! Retrieve the location and label
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_Head(io, aem, cZ1)
          write (unit=LU_INQ, fmt="(""HEA"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
        case (kOpDIS)
          !****************************************************************************
          ! Here for the DIS command -- extract the discharge
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cIO_WorldDischarge(io, cAEM_Discharge(io, aem, cZ1))
          write (unit=LU_INQ, fmt="(""DIS"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
        case (kOpPOT)
          !****************************************************************************
          ! Here for the POT command -- extract the potential
          !****************************************************************************
          ! Retrieve the location and label
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cAEM_Potential(io, aem, cZ1)
          write (unit=LU_INQ, fmt="(""POT"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
        case (kOpVEL)
          !****************************************************************************
          ! Here for the VEL command -- extract the velocity
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cIO_WorldDischarge(io, cAEM_Velocity(io, aem, cZ1))
          write (unit=LU_INQ, fmt="(""VEL"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
        case (kOpFLO)
          !****************************************************************************
          ! Here for the FLO command -- extract the flow between points
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rIO_WorldLength(io, rAEM_Flow(io, aem, (/cZ1, cZ2/)), cZ2-cZ1)
          write (unit=LU_INQ, fmt="(""FLO"", "", "", i9, 1x, 5("", "", e18.11))") iLabel, cZ1, cZ2, rValue
        case (kOpBDY)
          !****************************************************************************
          ! Here for the BDY command -- extract data for an inset BDY element
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          iLabel = iIO_GetInteger(io, "iLabel")
          rFlux = rAEM_Flow(io, aem, (/cZ1, cZ2/))
          rHead = rAEM_Head(io, aem, rHALF*(cZ1+cZ2))
          write (unit=LU_INQ, fmt="(""BDY"", "", "", i9, 1x, 6("", "", e18.11))") iLabel, cZ1, cZ2, rHead, rFlux
        case (kOpRCH)
          !****************************************************************************
          ! Here for the RCH command -- extract the net recharge rate
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_Recharge(io, aem, cZ1)
          write (unit=LU_INQ, fmt="(""RCH"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
        case (kOpSAT)
          !****************************************************************************
          ! Here for the SAT command -- extract the saturated thickness
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_SatdThick(io, aem, cZ1)
          write (unit=LU_INQ, fmt="(""SAT"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
        case (kOpGAG)
          !****************************************************************************
          ! Here for the GAG command -- extract the gage value at the end of an LS2 string
          !****************************************************************************
          ! Retrieve the location and label
          ! Retrieve the maximum time
          call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
          iGageID = iIO_GetInteger(io, "iGageID")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rLS2_Gage(io, aem%ls2, iID)
          write (unit=LU_INQ, fmt='("GAG", ", ", i9, ", ", 1x, i9, 1x, ", ", e18.11)') iID, iGageId, rValue
        case (kOpWL0)
          !****************************************************************************
          ! Here for the WL0 command -- extract information about WL0 elements
          !****************************************************************************
          call WL0_Inquiry(io, aem%wl0, aem%aqu, LU_INQ)
        case (kOpWL1) !**pd
          !****************************************************************************
          ! Here for the WL1 command -- extract information about WL1 elements
          !****************************************************************************
          call WL1_Inquiry(io, aem%wl1, aem%aqu, LU_INQ) !**pd
        case (kOpLS0)
          !****************************************************************************
          ! Here for the LS0 command -- extract information about LS0 elements
          !****************************************************************************
          call LS0_Inquiry(io, aem%ls0, LU_INQ)
        case (kOpLS1)
          !****************************************************************************
          ! Here for the LS1 command -- extract information about LS1 elements
          !****************************************************************************
          call LS1_Inquiry(io, aem%ls1, LU_INQ)
        case (kOpLS2)
          !****************************************************************************
          ! Here for the LS2 command -- extract information about LS2 elements
          !****************************************************************************
          call LS2_Inquiry(io, aem%ls2, LU_INQ)
        case (kOpHB0)
          !****************************************************************************
          ! Here for the HB0 command -- extract information about HB0 elements
          !****************************************************************************
          call HB0_Inquiry(io, aem%hb0, LU_INQ)
        case (kOpAS0)
          !****************************************************************************
          ! Here for the AS0 command -- extract information about AS0 elements
          !****************************************************************************
          call AS0_Inquiry(io, aem%as0_top, LU_INQ)
          call AS0_Inquiry(io, aem%as0_bottom, LU_INQ)
        case (kOpAQU)
          !****************************************************************************
          ! Here for the AQU command -- extract information about AEM elements
          !****************************************************************************
          call AQU_Inquiry(io, aem%aqu, LU_INQ)
#ifndef __GPL__
        case (kOpCW0)
          !****************************************************************************
          ! Here for the CW0 command -- extract information about collector wells
          !****************************************************************************
          call CW0_Inquiry(io, aem%cw0, aem%aqu, LU_INQ)
#endif
        case default
          continue
      end select
    end do

    call IO_MessageText(io, "Leaving INQ module")

    return
  end subroutine INQ_Read

end module m_inq

