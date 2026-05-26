module a_inq

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

  ! Data extraction/inquiry/check module
  ! V.A.Kelson  10/23/98
  !
  ! Commands:
  !
  !   OPT < inq | csv >   Select output format. 'inq' (default) writes all
  !                        output to a single file named by FIL.  'csv' writes
  !                        one CSV file per directive type; files are named
  !                        <prefix>_<type>.csv (with FIL prefix) or
  !                        <type>_inq.csv (without FIL).
  !
  !   FIL < filename >     Traditional: output file (no extension appended).
  !                        CSV mode: prefix for all CSV file names.
  !
  !   HEA  (X, Y) label              Extracts the head at the point(X, Y)
  !   DIS  (X, Y) label              Extracts the discharge at the point(X, Y)
  !   VEL  (X, Y) label              Extracts the velocity at the point(X, Y)
  !   FLO  (X1, Y1) (X2, Y2) label   Extracts the integrated flow across the path
  !   BDY  (X1, Y1) (X2, Y2) label   Extracts head and flux across the path
  !   WL0 / WL1 / LS0 / LS1 / LS2 / HB0 / AS0 / AQU / CW0
  !                                   Writes element-check information

  use u_constants
  use u_io
  use f_aem
  use p_packages

  implicit none

  public

  logical, private, save :: fINQFileOpen = .false.

  ! Indices for the per-type spatial CSV tracking array
  integer(kind=AE_INT), private, parameter :: ICSV_HEA = 1
  integer(kind=AE_INT), private, parameter :: ICSV_DIS = 2
  integer(kind=AE_INT), private, parameter :: ICSV_POT = 3
  integer(kind=AE_INT), private, parameter :: ICSV_VEL = 4
  integer(kind=AE_INT), private, parameter :: ICSV_FLO = 5
  integer(kind=AE_INT), private, parameter :: ICSV_BDY = 6
  integer(kind=AE_INT), private, parameter :: ICSV_RCH = 7
  integer(kind=AE_INT), private, parameter :: ICSV_SAT = 8
  integer(kind=AE_INT), private, parameter :: ICSV_GAG = 9
  integer(kind=AE_INT), private, parameter :: ICSV_NSPATIAL = 9

contains


  function sCSVFileName(sBase, sModule) result(sFName)
    ! Returns <sBase>_<sModule>.csv if sBase is non-empty, else <sModule>_inq.csv
    character(len=*), intent(in) :: sBase, sModule
    character(len=256) :: sFName
    if (len_trim(sBase) > 0) then
      sFName = trim(sBase) // "_" // trim(sModule) // ".csv"
    else
      sFName = trim(sModule) // "_inq.csv"
    end if
  end function sCSVFileName


  subroutine vCSVSpatialOpen(sFile, iIdx, lFirst, iLU, lWroteHeader, iStat)
    ! Opens a spatial CSV file for the given directive index.
    ! On the first write in a session it creates/truncates the file (lWroteHeader returns .true.
    ! so the caller knows to write the header).  On subsequent writes it appends.
    character(len=*), intent(in) :: sFile
    integer(kind=AE_INT), intent(in) :: iIdx
    logical, dimension(ICSV_NSPATIAL), intent(inout) :: lFirst
    integer(kind=AE_INT), intent(out) :: iLU
    logical, intent(out) :: lWroteHeader
    integer(kind=AE_INT), intent(out) :: iStat

    lWroteHeader = lFirst(iIdx)
    if (lFirst(iIdx)) then
      open(newunit=iLU, file=trim(sFile), action="WRITE", &
           status="REPLACE", iostat=iStat)
      lFirst(iIdx) = .false.
    else
      open(newunit=iLU, file=trim(sFile), action="WRITE", &
           status="OLD", position="APPEND", iostat=iStat)
    end if
  end subroutine vCSVSpatialOpen


  subroutine INQ_Read(io, pkg)
    ! Reads and processes all input commands for the INQ module

    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io

    type(DIRECTIVE), dimension(21), parameter :: dirDirectives = &
                       (/dirFIL, dirOPT, &
                       dirHEA, dirDIS, dirPOT, dirVEL, dirFLO, dirRCH, &
                       dirSAT, dirBDY, dirGAG, &
                       dirWL0, dirWL1, dirLS0, &
                       dirLS1, dirLS2, dirHB0, &
                       dirAS0, &
                       dirCW0, &
                       dirAQU, dirEND/)

    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    real(kind=AE_REAL) :: rValue, rHead, rFlux
    complex(kind=AE_REAL) :: cZ1, cZ2, cValue
    integer(kind=AE_INT) :: iLabel, iID, iGageId
    logical :: lFlag
    character(len=132) :: sFile
    character(len=16) :: sOpt

    ! CSV mode state
    logical :: lCSV
    character(len=256) :: sCsvBase
    integer(kind=AE_INT) :: iCSVLU
    logical :: lWroteHeader
    logical, dimension(ICSV_NSPATIAL) :: lCSVSpatialFirst

    call IO_MessageText(io, "Reading INQ module input")

    fINQFileOpen = .false.
    lCSV = .false.
    sCsvBase = ""
    lCSVSpatialFirst = .true.
    call IO_MessageText(io, "Entering inquiry module INQ")
    call IO_Assert(io, (associated(pkg)), "INQ_Read: No PKG_DOMAIN object")

    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          call IO_Assert(io, .false., "INQ_Read: I/O Error")
        case (kOpFileEOF)
          call IO_Assert(io, .false., "INQ_Read: Unexpected EOF")
        case (kOpEND)
          if (fINQFileOpen) close(unit=LU_INQ)
          exit

        case (kOpOPT)
          !****************************************************************************
          ! Here for the OPT command -- select output format
          !****************************************************************************
          sOpt = sIO_GetField(io, "sOpt", allowed=(/"inq", "csv", "INQ", "CSV"/))
          select case (trim(sOpt))
            case ("csv", "CSV")
              lCSV = .true.
            case default
              lCSV = .false.
          end select

        case (kOpFIL)
          !****************************************************************************
          ! Here for the FIL command
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          if (lCSV) then
            sCsvBase = trim(sFile)
          else
            if (fINQFileOpen) then
              close(unit=LU_INQ)
              fINQFileOpen = .false.
            end if
            open(unit=LU_INQ, file=sFile, action="WRITE", status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: Open failed")
            fINQFileOpen = .true.
          end if

        case (kOpHEA)
          !****************************************************************************
          ! Here for the HEA command -- extract the head
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_Head(io, pkg%aem, cZ1)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "hea"), ICSV_HEA, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, head"")")
            write (unit=iCSVLU, fmt="(i9, 3("", "", e18.11))") iLabel, cZ1, rValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""HEA"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
          end if

        case (kOpDIS)
          !****************************************************************************
          ! Here for the DIS command -- extract the discharge
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cIO_WorldDischarge(io, cAEM_Discharge(io, pkg%aem, cZ1))
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "dis"), ICSV_DIS, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, discharge_x, discharge_y"")")
            write (unit=iCSVLU, fmt="(i9, 4("", "", e18.11))") iLabel, cZ1, cValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""DIS"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
          end if

        case (kOpPOT)
          !****************************************************************************
          ! Here for the POT command -- extract the potential
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cAEM_Potential(io, pkg%aem, cZ1)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "pot"), ICSV_POT, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, potential_r, potential_i"")")
            write (unit=iCSVLU, fmt="(i9, 4("", "", e18.11))") iLabel, cZ1, cValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""POT"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
          end if

        case (kOpVEL)
          !****************************************************************************
          ! Here for the VEL command -- extract the velocity
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          cValue = cIO_WorldDischarge(io, cAEM_Velocity(io, pkg%aem, cZ1))
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "vel"), ICSV_VEL, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, velocity_x, velocity_y"")")
            write (unit=iCSVLU, fmt="(i9, 4("", "", e18.11))") iLabel, cZ1, cValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""VEL"", "", "", i9, 1x, 4("", "", e18.11))") iLabel, cZ1, cValue
          end if

        case (kOpFLO)
          !****************************************************************************
          ! Here for the FLO command -- extract the flow between points
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rIO_WorldLength(io, rAEM_Flow(io, pkg%aem, (/cZ1, cZ2/)), cZ2-cZ1)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "flo"), ICSV_FLO, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x1, y1, x2, y2, flow"")")
            write (unit=iCSVLU, fmt="(i9, 5("", "", e18.11))") iLabel, cZ1, cZ2, rValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""FLO"", "", "", i9, 1x, 5("", "", e18.11))") iLabel, cZ1, cZ2, rValue
          end if

        case (kOpBDY)
          !****************************************************************************
          ! Here for the BDY command -- extract data for an inset BDY element
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          iLabel = iIO_GetInteger(io, "iLabel")
          rFlux = rAEM_Flow(io, pkg%aem, (/cZ1, cZ2/))
          rHead = rAEM_Head(io, pkg%aem, rHALF*(cZ1+cZ2))
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "bdy"), ICSV_BDY, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x1, y1, x2, y2, head, flux"")")
            write (unit=iCSVLU, fmt="(i9, 6("", "", e18.11))") iLabel, cZ1, cZ2, rHead, rFlux
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""BDY"", "", "", i9, 1x, 6("", "", e18.11))") iLabel, cZ1, cZ2, rHead, rFlux
          end if

        case (kOpRCH)
          !****************************************************************************
          ! Here for the RCH command -- extract the net recharge rate
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_Recharge(io, pkg%aem, cZ1)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "rch"), ICSV_RCH, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, recharge"")")
            write (unit=iCSVLU, fmt="(i9, 3("", "", e18.11))") iLabel, cZ1, rValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""RCH"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
          end if

        case (kOpSAT)
          !****************************************************************************
          ! Here for the SAT command -- extract the saturated thickness
          !****************************************************************************
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rAEM_SatdThick(io, pkg%aem, cZ1)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "sat"), ICSV_SAT, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""label, x, y, sat_thickness"")")
            write (unit=iCSVLU, fmt="(i9, 3("", "", e18.11))") iLabel, cZ1, rValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt="(""SAT"", "", "", i9, 1x, 3("", "", e18.11))") iLabel, cZ1, rValue
          end if

        case (kOpGAG)
          !****************************************************************************
          ! Here for the GAG command -- extract the gage value at the end of an LS2 string
          !****************************************************************************
          iGageID = iIO_GetInteger(io, "iGageID")
          iLabel = iIO_GetInteger(io, "iLabel")
          rValue = rLS2_Gage(io, pkg%ls2, iID)
          if (lCSV) then
            call vCSVSpatialOpen(sCSVFileName(sCsvBase, "gag"), ICSV_GAG, lCSVSpatialFirst, iCSVLU, lWroteHeader, iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            if (lWroteHeader) write (unit=iCSVLU, fmt="(""gage_id, label, flow"")")
            write (unit=iCSVLU, fmt="(i9, "", "", i9, "", "", e18.11)") iID, iGageId, rValue
            close(iCSVLU)
          else
            call IO_Assert(io, fINQFileOpen, "INQ_Read: No output file")
            write (unit=LU_INQ, fmt='("GAG", ", ", i9, ", ", 1x, i9, 1x, ", ", e18.11)') iID, iGageId, rValue
          end if

        case (kOpWL0)
          !****************************************************************************
          ! Here for the WL0 command -- extract information about WL0 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "wl0"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call WL0_Inquiry(io, pkg%wl0, pkg%aqu, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call WL0_Inquiry(io, pkg%wl0, pkg%aqu, LU_INQ)
          end if

        case (kOpWL1)
          !****************************************************************************
          ! Here for the WL1 command -- extract information about WL1 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "wl1"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call WL1_Inquiry(io, pkg%wl1, pkg%aqu, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call WL1_Inquiry(io, pkg%wl1, pkg%aqu, LU_INQ)
          end if

        case (kOpLS0)
          !****************************************************************************
          ! Here for the LS0 command -- extract information about LS0 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "ls0"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call LS0_Inquiry(io, pkg%ls0, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call LS0_Inquiry(io, pkg%ls0, LU_INQ)
          end if

        case (kOpLS1)
          !****************************************************************************
          ! Here for the LS1 command -- extract information about LS1 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "ls1"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call LS1_Inquiry(io, pkg%ls1, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call LS1_Inquiry(io, pkg%ls1, LU_INQ)
          end if

        case (kOpLS2)
          !****************************************************************************
          ! Here for the LS2 command -- extract information about LS2 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "ls2"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call LS2_Inquiry(io, pkg%ls2, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call LS2_Inquiry(io, pkg%ls2, LU_INQ)
          end if

        case (kOpHB0)
          !****************************************************************************
          ! Here for the HB0 command -- extract information about HB0 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "hb0"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call HB0_Inquiry(io, pkg%hb0, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call HB0_Inquiry(io, pkg%hb0, LU_INQ)
          end if

        case (kOpAS0)
          !****************************************************************************
          ! Here for the AS0 command -- extract information about AS0 elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "as0"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call AS0_Inquiry(io, pkg%as0_top, iCSVLU, lCSV=.true.)
            call AS0_Inquiry(io, pkg%as0_bottom, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call AS0_Inquiry(io, pkg%as0_top, LU_INQ)
            call AS0_Inquiry(io, pkg%as0_bottom, LU_INQ)
          end if

        case (kOpAQU)
          !****************************************************************************
          ! Here for the AQU command -- extract information about AEM elements
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "aqu"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call AQU_Inquiry(io, pkg%aqu, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call AQU_Inquiry(io, pkg%aqu, LU_INQ)
          end if

        case (kOpCW0)
          !****************************************************************************
          ! Here for the CW0 command -- extract information about collector wells
          !****************************************************************************
          if (lCSV) then
            open(newunit=iCSVLU, file=sCSVFileName(sCsvBase, "cw0"), action="WRITE", &
                 status="REPLACE", iostat=iStat)
            call IO_Assert(io, (iStat == 0), "INQ_Read: CSV open failed")
            call CW0_Inquiry(io, pkg%cw0, pkg%aem, iCSVLU, lCSV=.true.)
            close(iCSVLU)
          else
            call CW0_Inquiry(io, pkg%cw0, pkg%aem, LU_INQ)
          end if

        case default
          continue
      end select
    end do

    call IO_MessageText(io, "Leaving INQ module")

    return
  end subroutine INQ_Read

end module a_inq
