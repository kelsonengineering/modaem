module a_tr0

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

  ! Analysis module for streamline tracing in 3-D.
  ! V.A.Kelson  10/20/98
  !
  ! This module allows for the tracing of streamlines in two-dimensions, and supports
  ! a variety of options for output and selection of release locations
  !
  use u_constants
  use u_io
  use m_wl0
  use m_ls0
  use m_ls1
  use m_hb0
  use m_wl1
  use m_aem

  implicit none

  public

  ! Global variables for this module
  real(kind=AE_REAL), private, save :: rTR0WinX1, rTR0WinY1, rTR0WinX2, rTR0WinY2
  real(kind=AE_REAL), private, save :: rTR0MaxTime,rTR0WinDiagonal
  real(kind=AE_REAL), private, save :: rTR0Step, rTR0Prox, rTR0Frac, rTR0WellProx
  real(kind=AE_REAL), private, save :: rTR0Retard, rTR0Decay

  ! Local constants
  integer(kind=AE_INT), private, parameter :: kTR0Timeout = 0
  integer(kind=AE_INT), private, parameter :: kTR0LeftWindow = 1
  integer(kind=AE_INT), private, parameter :: kTR0HitElement = 2
  integer(kind=AE_INT), private, parameter :: kTR0AquiferDry = 3
  integer(kind=AE_INT), private, parameter :: kTR0Stagnation = 4
  integer(kind=AE_INT), private, parameter :: kTR0LeftTop = 5
  integer(kind=AE_INT), private, parameter :: kTR0LeftBottom = 6
  ! Check factor for termination at HB0 elements
  ! Default step factor -- fraction of tracing window size to be used for stepsize
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_STEP_FACTOR = 0.003_AE_REAL
  ! Default proximity factor -- number of step sizes away from elements for reduced step size
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_PROX = 10.0_AE_REAL
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_WELLPROX = 10.0_AE_REAL
  ! Termination factors -- multiply by step for the termination tolerance at BCs
  real(kind=AE_REAL), private, parameter :: TR0_WLX_TERMINATION_FACTOR = 0.2_AE_REAL
  real(kind=AE_REAL), private, parameter :: TR0_LS0_TERMINATION_FACTOR = 0.2_AE_REAL
  real(kind=AE_REAL), private, parameter :: TR0_LS1_TERMINATION_FACTOR = 0.2_AE_REAL
  real(kind=AE_REAL), private, parameter :: TR0_HB0_TERMINATION_FACTOR = 0.1_AE_REAL
  ! Default fraction of basic step size to use near elements
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_FRAC = 0.1_AE_REAL
  ! Default maximum tracing time
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_MAX_TIME = HUGE(AE_REAL)
  ! Default transport parameters(not yet implemented)
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_RETARDATION = rONE
  real(kind=AE_REAL), private, parameter :: TR0_DEFAULT_DECAY = rZERO
  ! Stagnation tolerance(number of step sizes the particle must move in 20 steps)
  real(kind=AE_REAL), private, parameter :: TR0_STAGNATION_TOLERANCE = 1.2_AE_REAL
  ! Distance factor for releasing particles from a well
  real(kind=AE_REAL), private, parameter :: TR0_WL0_DISTANCE_FACTOR = 4.0_AE_REAL


contains


  subroutine TR0_Read(io, aem)
    ! This reads and processes the input commands for the TR0 package
    ! Commands:
    !
    !   FIL < filename > -        Output file. The extension '.tr0' will be added to the filename.
    !                               For each pathline, the output file will contain:
    !                               A header line of "START X0, Y0, Z0, T0, C0, L, Dir", where:
    !                                 (X0, Y0, Z0) - Starting coordinates
    !                                 (T0) - Starting time(usually 0.0)
    !                                 (C0) - Starting concentration
    !                                 (L)  - Layer for tracing
    !                                 (Dir) - Direction(-1 backward) (1 forward)
    !                               For each point in each pathline, records of X, Y, Z, T, C follow, where
    !                                 (X, Y, Z is the coordinate of a point on the streamline)
    !                                 (T is the time of travel)
    !                                 (C is the concentration)
    !                               At the end of the pathline, "END flag, elem, vertex" is written where
    !                                 (flag) - Termination flag(0 = timeout, 1 = out of window, 2 = hit boundary,
    !                                                            3 = aquifer dry, 4 = stagnation point)
    !                                 (elem) - If at an element, the element type
    !                                 (vertex) - If at an element, the vertex ID
    !
    !                               NOTE: Since TR0 is a 2-D tracing package, the Z-coordinate remains
    !                               constant throughout all points
    !
    !   TRA  retard decay           Transport parameters along streamlines
    !                               retard is the retardation factor, and decay is the
    !                               first-order decay coefficient. Concentrations are computed
    !                               along pathlines using an analytical solution.
    !
    !   WIN  (X1, Y1) (X2, Y2)        Sets the tracing window. Particles that leave the window are terminated.
    !                               Default tuning parameters(see below) are derived from the window size.
    !
    !   TUN  step prox frac small   Tuning parameters for the tracing algorithm
    !                                 (step) - The base step size
    !                                 (prox) - The proximity(in terms of the current step size)
    !                                          to boundary conditions for reducing the step size
    !                                 (frac) - Factor for step size reductions
    !                                 (small) - Smallest allowable step size
    !
    !   TIM  max-time               Specifies the maximum time allowed for particle tracing
    !
    !   POI  (X0, Y0) Z0 T0 C0 Dir                   Releases a single point at the specified location/
    !                                                layer/direction
    !
    !   LIN  (X1, Y1) (X1, Y1) N Z0 T0 C0 dir        Releases N particles along the line from(X0, Y0, Z0) to
    !                                                (X1, Y1, Z0).
    !
    !   GRI  (X1, Y1) (X2, Y2) N Z0 T0 C0 dir        Releases a grid of particles in the subwindow
    !                                                (X1, Y1) - (X2, Y2) with resolution N
    !
    !   WL0  id N Z0 T0 C0          Releases N particles in reverse from the well bore of WL0 element #id
    !                               in layer L. The direction of travel is determined by the sign of the
    !                               well pumping rate(forward for injection wells, backward for pumping
    !                               wells).
    !   LS2  N T0 C0 Dir			  Releases N particles forward from "losing" resistance line sinks (Dir>0),
    !                               or backward from "gaining" resistance line sinks (Dir<0). Particles are
    !								  released from the bottom of the resistance layer of each line sink.
    !   CW0  id N Z0 T0 C0 R bufsz  Releases N particles a radius R from the center of collector well #id
    !                               in layer L. The direction of travel is determined by the sign of the
    !                               well pumping rate(forward for injection wells, backward for pumping
    !                               wells). The bufsz argument is the size of a temporary buffer for
    !                               holding vertices -- typically 200 or more is sufficient.


    ! Argument list
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! Locals -- for Directive parsing
#ifndef __GPL__
    type(DIRECTIVE), dimension(12), parameter :: dirDirectives = (/ dirCW0, &
#else
                       type(DIRECTIVE), dimension(11), parameter :: dirDirectives = (/&
#endif
    dirEND, &
              dirWIN, dirTRA, dirFIL, dirTUN, dirTIM, &
              dirPOI, dirLIN, dirGRI, dirWL0, dirLS2 &
              /)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    ! Placeholders for tracing parameters
    real(kind=AE_REAL) :: rX0, rY0, rX1, rX2, rY1, rY2, rZ0
    real(kind=AE_REAL) :: rStep, rProx, rFrac, rMaxTime, rWellProx, rDecay, rRetard
    real(kind=AE_REAL) :: rT0, rC0, rDX, rDY, rThisZ0, rThisT0, rThisC0
    complex(kind=AE_REAL) :: cZ0, cZ1, cZ2, cZWell, cWWell, cDelta
    real(kind=AE_REAL) :: rOrient, rAngle, rDischarge, rRadius
    integer(kind=AE_INT) :: i, j, iDir, iFlag, iElementType, iElementID, iElementVtx
    integer(kind=AE_INT) :: iWellID, iNParticles, iFWLIndex, iRes, iResX, iResY
    real(kind=AE_REAL) :: rMinX, rMaxX, rMinY, rMaxY, rDelta, rSizeX, rSizeY, rMinWidth
    character(len=132) :: sFile
    logical :: lFileOpen, lWindowSet, lWellFound, lFlag
    ! Buffers for CW0 tracing
    real(kind=AE_REAL) :: rCDistance, rZ0Save, rBottom
    ! Specials for LS2 tracing
    type(LS2_COLLECTION), pointer :: ls2
    type(LS2_STRING), pointer :: ls2_str
    type(LS2_VERTEX), pointer :: this_ls2_vtx, next_ls2_vtx
    integer (kind=AE_INT) :: istr, ivtx

    ! Scratch LUs for CW0 tracing
    integer(kind=AE_INT), parameter :: iFwdLU = 130

    ! Clear the status flags
    lFileOpen = .false.
    lWindowSet = .false.
    rTR0MaxTime = TR0_DEFAULT_MAX_TIME
    rTR0Retard = TR0_DEFAULT_RETARDATION
    rTR0Decay = TR0_DEFAULT_DECAY
    call IO_MessageText(io, "Entering 2-D trace module TR0")

    call IO_Assert(io, (associated(aem)), "TR0_Read: No AEM_DOMAIN object")

    ! The remainder of this routine uses IO_InputRecord to process
    ! the model input file.
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "TR0_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all xxx_Read routines.
          call IO_Assert(io, .false., "TR0_Read: Unexpected EOF")
          exit
        case (kOpEND)
          ! END mark was found. Exit the file parser.
          exit
        case (kOpFIL)
          !****************************************************************************
          ! Here for the FIL command -- open the output file
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          ! If the output file is open, close it now.
          if (lFileOpen) then
            close(unit=LU_TRACE)
            lFileOpen = .false.
          end if
          ! Open the new output file.  If the open fails, write a message.
          sFile = trim(sFile) // ".tr0"
          open(unit=LU_TRACE, file = sFile, action = "WRITE", status="REPLACE", iostat=iStat)
          call IO_Assert(io, (iStat == 0), "TR0_Read: Open failed")
          lFileOpen = .true.
        case (kOpWIN)
          !****************************************************************************
          ! Here for the WIN command -- Set the window
          !****************************************************************************
          ! Retrive the window coordinates
          lWindowSet = .false.
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          rTR0WinX1 = min(real(cZ1), real(cZ2))
          rTR0WinY1 = min(aimag(cZ1), aimag(cZ2))
          rTR0WinX2 = max(real(cZ1), real(cZ2))
          rTR0WinY2 = max(aimag(cZ1), aimag(cZ2))
          rTR0WinDiagonal = sqrt((rTR0WinX2-rTR0WinX1)**2 + (rTR0WinY2-rTR0WinY1)**2)
          ! Compute the default tuning parameters
          rTR0Step = TR0_DEFAULT_STEP_FACTOR * rTR0WinDiagonal
          rTR0Prox = TR0_DEFAULT_PROX
          rTR0Frac = TR0_DEFAULT_FRAC
          rTR0WellProx = TR0_DEFAULT_WELLPROX
          lWindowSet = .true.
          call IO_ErrorText(io, " >> Default tuning parameters are selected for the new window")
        case (kOpTUN)
          !****************************************************************************
          ! Here for the TUN command -- set the tracing algorithm tuning
          !****************************************************************************
          ! Retrieve the tuning parameters
          rTR0Step = rIO_GetReal(io, "rStep", def=TR0_DEFAULT_STEP_FACTOR * rTR0WinDiagonal)
          rTR0Prox = rIO_GetReal(io, "rProx", def=TR0_DEFAULT_PROX, minimum=0.01_AE_REAL)
          rTR0Frac = rIO_GetReal(io, "rFrac", def=TR0_DEFAULT_PROX, minimum=0.01_AE_REAL)
          rTR0WellProx = rIO_GetReal(io, "rWellProx", minimum=1.0E-4_AE_REAL)
          call IO_ErrorText(io, " >> Specified tuning parameters are in use")
        case (kOpTRA)
          !****************************************************************************
          ! Here for the TRA command -- set the transport parameters
          !****************************************************************************
          ! Retrieve the tuning parameters
          rTR0Retard = rIO_GetReal(io, "rTR0Retard", def=rZERO, minimum=rZERO)
          rTR0Decay = rIO_GetReal(io, "rTR0Decay", def=rZERO, minimum=rZERO)
        case (kOpTIM)
          !****************************************************************************
          ! Here for the TIM command -- set the maximum tracing time
          !****************************************************************************
          ! Retrieve the maximum time
          rTR0MaxTime = rIO_GetReal(io, "rMaxTime", minimum=rZERO)
          write (unit=IO_MessageBuffer, fmt=*) ">> Streamline termination time = ",rTR0MaxTime
          call IO_ErrorText(io)
        case (kOpPOI)
          !****************************************************************************
          ! Here for the POI command -- release a single point
          !****************************************************************************
          ! Retrieve the maximum time
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          cZ0 = cIO_GetCoordinate(io, "cZ0")
          rZ0 = rIO_GetReal(io, "rZ0", def=rHUGE)
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          iDir = iIO_GetInteger(io, "iDir", def=1, allowed=(/-1, 1/))
          ! For the purpose of the command interpreter, the returned flag information is ignored
          write (unit=IO_MessageBuffer, fmt="("" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)") cZ0, iDir
          call IO_MessageText(io)
          rX0 = real(cZ0)
          rY0 = aimag(cZ0)
          call TR0_Trace(io, aem, rX0, rY0, rZ0, rT0, rC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
        case (kOpLIN)
          !****************************************************************************
          ! Here for the LIN command -- release points along a line
          !****************************************************************************
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          rZ0 = rIO_GetReal(io, "rZ0", def=rHUGE)
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          iRes = iIO_GetInteger(io, "iRes", def=10, minimum=1)
          iDir = iIO_GetInteger(io, "iDir", def=1, allowed=(/-1, 1/))
          call IO_Assert(io, (iRes > 2), "TR0_Read: Bad line resolution")
          ! Fix iDir, then trace away for all the points in a grid
          iDir = iDir/abs(iDir)
          write (unit=LU_TRACE, fmt="('LINE ', 2i10)") iRes
          cDelta = (cZ2-cZ1) / (iRes-1)
          do i = 1, iResX
            write (unit=LU_TRACE, fmt="('POINT', 2i10)") i
            cZ0 = (i-1)*cDelta + cZ1
            rX0 = real(cZ0)
            rY0 = aimag(cZ0)
            rThisZ0 = rZ0
            rThisT0 = rT0
            rThisC0 = rC0
            write (unit=IO_MessageBuffer, &
                   fmt="("" Point: "", i5, "" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)" &
                   ) i, rX0, rY0, iDir
            call IO_MessageText(io)
            call TR0_Trace(io, aem, rX0, rY0, rThisZ0, rThisT0, rThisC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
          end do
        case (kOpGRI)
          !****************************************************************************
          ! Here for the GRI command -- release points in a grid
          !****************************************************************************
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          cZ2 = cIO_GetCoordinate(io, "cZ2")
          iRes = iIO_GetInteger(io, "iRes", def=10, minimum=2)
          rZ0 = rIO_GetReal(io, "rZ0", def=rHUGE)
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          iDir = iIO_GetInteger(io, "iDir", def=1, allowed=(/-1, 1/))
          ! Adjust the grid
          rX1 = real(cZ1)
          rY1 = aimag(cZ1)
          rX2 = real(cZ2)
          rY2 = aimag(cZ2)
          rSizeX = rX2-rX1
          rSizeY = rY2-rY1
          if (rSizeX >= rSizeY) then
            iResX = iRes
            rDelta = rSizeX/real(iRes-1, AE_REAL)
            iResY = int(rSizeY/rDelta)+1
          else
            iResY = iRes
            rDelta = rSizeY/real(iRes-1, AE_REAL)
            iResX = int(rSizeX/rDelta)+1
          end if
          ! Make the real window span the center of the specified window
          rMinX = rHALF * (rX1+rX2 - (iResX-1)*rDelta)
          rMaxX = rMinX + (iResX-1)*rDelta
          rMinY = rHALF * (rY1+rY2 - (iResY-1)*rDelta)
          rMaxY = rMinY + (iResY-1)*rDelta
          ! Fix iDir, then trace away for all the points in a grid
          iDir = iDir/abs(iDir)
          write (unit=LU_TRACE, fmt="('GRID ', 2i10)") iResX, iResY
          do j = 1, iResY
            do i = 1, iResX
              write (unit=LU_TRACE, fmt="('CELL ', 2i10)") j, i
              rX0 = (i-1)*rDelta + rMinX
              rY0 = rMaxY - (j-1)*rDelta
              rThisZ0 = rZ0
              rThisT0 = rT0
              rThisC0 = rC0
              write (unit=IO_MessageBuffer, &
                     fmt="("" Row: "", i5, "" Col: "", i5, "" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)" &
                     ) j, i, rX0, rY0, iDir
              call IO_MessageText(io)
              call TR0_Trace(io, aem, rX0, rY0, rThisZ0, rThisT0, rThisC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
            end do
          end do
        case (kOpWL0)
          !****************************************************************************
          ! Here for the WL0 command -- release points from a well
          !****************************************************************************
          ! This routine requires that the location of the well element be found(using ifWL0FindWell)
          ! and then using the discharge vector at the well(less the discharge due to the well) to
          ! orient the release points.  This routine does not check to see if the well has pumped the
          ! aquifer dry; TR0_Trace will report that in the END record.
          ! Retrieve the maximum time
          !
          ! Modified 6/25/2001 by VAK. Added rMinWidth value for Idaho delineations. If the
          ! computed "ultimate" capture zone width is less than rMinWidth, then ModAEM traces
          ! only one particle back from the wellscreen, in the direction of flow.
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          iWellID = iIO_GetInteger(io, "iWellID")
          iNParticles = iIO_GetInteger(io, "iNParticles")
          rZ0 = rIO_GetReal(io, "rZ0", def=rZERO)
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          rMinWidth = rIO_GetReal(io, "rMinWidth", def=rZERO)
          call WL0_FindWell(io, aem%wl0, iWellID, cZWell, rDischarge, rRadius, iFWLIndex, lWellFound)
          call IO_Assert(io, lWellFound, "TR0_Read: Specified well not found")
          ! Forward trace from injection wells; backward trace from pumping wells
          if (rDischarge < 0) then
            iDir = 1
            call IO_ErrorText(io, " >> Forward tracing is selected from injection well")
          else
            iDir = -1
            call IO_ErrorText(io, " >> Reverse tracing is selected from pumping well")
          end if
          ! Now, compute the orientation of the discharge vector, less the well
          cWWell = cAEM_DischargeAtWell(io, aem, cZWell, iFWLIndex)
          if (cWWell == cZERO) then
            rOrient = rZERO
          else
            rOrient = atan2(aimag(cWWell), real(cWWell))
          end if
          ! Distribute particles about the well bore, distributed evenly about the orientation vector
          if (abs(rDischarge)/abs(cWWell) <= rMinWidth) then
            write (unit=IO_MessageBuffer, &
                   fmt="('Found small well with ID ', i5, ' -- tracing one particle')" &
                   ) iWellID
            call IO_MessageText(io)
            iNParticles = 1
          end if
          write (unit=IO_MessageBuffer, fmt="(""Releasing "", i5, "" particles from well ID "", i5)") iNParticles, iWellID
          write (unit=LU_TRACE, fmt="('WL0   ', i10)") iWellID
          call IO_MessageText(io)
          do i = 1, iNParticles
            rAngle = rOrient + (float(i-1)+rHALF) * 2*rPI / float(iNParticles)
            rX1 = real(cZWell) + rRadius * TR0_WL0_DISTANCE_FACTOR * cos(rAngle)
            rY1 = aimag(cZWell) + rRadius * TR0_WL0_DISTANCE_FACTOR * sin(rAngle)
            write (unit=IO_MessageBuffer, fmt="("" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)") rX1, rY1, iDir
            call IO_MessageText(io)
            rThisZ0 = rZ0
            rThisT0 = rT0
            rThisC0 = rC0
            call TR0_Trace(io, aem, rX1, rY1, rThisZ0, rThisT0, rThisC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
          end do
        case (kOpLS2)
          !****************************************************************************
          ! Here for the LS2 command -- release points from resistance line-sinks
          !****************************************************************************
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          iNParticles = iIO_GetInteger(io, "iNParticles")
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          iDir = iIO_GetInteger(io, "iDir")
          ls2 => aem%ls2
          do istr = 1, ls2%iNStr
            ls2_str => ls2%Strings(istr)
            do ivtx = 1, ls2_str%iNPts-1
              this_ls2_vtx => ls2_str%Vertices(ivtx)
              next_ls2_vtx => ls2_str%Vertices(ivtx+1)
              ! Skip this segment?
              if (this_ls2_vtx%rStrength >= rZERO .and. iDir >= 0) cycle
              if (this_ls2_vtx%rStrength < rZERO .and. iDir <= 0) cycle
              ! Now, trace the particles...
              write (unit=IO_MessageBuffer, fmt="(""Releasing "", i5, "" particles from line sink ID "", i5, ""segment "", i5)") &
                     iNParticles, ls2_str%iID, ivtx
              write (unit=LU_TRACE, fmt="('LS2   ', i10, i10, e17.9)") ls2_str%iID, ivtx, &
                     this_ls2_vtx%rStrength*this_ls2_vtx%rLength/float(iNParticles)
              cDelta = (next_ls2_vtx%cZ - this_ls2_vtx%cZ) / float(iNParticles)
              do i = 1, iNParticles
                cZ0 = this_ls2_vtx%cZ + (float(i)-0.5) * cDelta
                rX1 = real(cZ0)
                rY1 = aimag(cZ0)
                write (unit=IO_MessageBuffer, fmt="("" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)") rX1, rY1, iDir
                call IO_MessageText(io)
                rThisZ0 = this_ls2_vtx%rCPDepth
                rThisT0 = rT0
                rThisC0 = rC0
                call TR0_Trace(io, aem, rX1, rY1, rThisZ0, rThisT0, rThisC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
              end do
            end do
          end do
#ifndef __GPL__
        case (kOpCW0)
          !****************************************************************************
          ! Here for the CW0 command -- release points from a collector well
          !****************************************************************************
          ! This routine requires that the location of the well element be found(using iCW0FindWell),
          ! then to release points some distance away from the well, both forwards and backwards.
          ! The forward trace provides the path into the well; it is written to a temporary file,
          ! then re-read to get its travel time and the points reversed for writing to the final         !
          ! This routine does not check to see if the well has pumped the
          ! aquifer dry; TR0_Trace will report that in the END record.
          ! Retrieve the maximum time
          call IO_Assert(io, (iStat == 0), "TR0_Read: I/O Error")
          call IO_Assert(io, lFileOpen, "TR0_Read: No output file")
          iWellID = iIO_GetInteger(io, "iWellID")
          iNParticles = iIO_GetInteger(io, "iNParticles", def=8, minimum=1)
          rZ0 = rIO_GetReal(io, "rZ0", def=rZERO)
          rT0 = rIO_GetReal(io, "rT0", def=rZERO)
          rC0 = rIO_GetReal(io, "rC0", def=rZERO)
          rCDistance = rIO_GetReal(io, "rC0", minimum=rONE)
          call IO_Assert(io, lWindowSet, "TR0_Read: No window specified")
          rZ0Save = rZ0

          ! Get the well location...
          call CW0_FindWell(io, aem%cw0, iWellID, cZWell, rDischarge, lWellFound)
          call IO_Assert(io, lWellFound, "TR0_Read: Specified well not found")

          ! Forward trace from injection wells; backward trace from pumping wells
          if (rDischarge < 0) then
            iDir = 1
            call IO_ErrorText(io, " >> Forward tracing is selected from injection well")
          else
            iDir = -1
            call IO_ErrorText(io, " >> Reverse tracing is selected from pumping well")
          end if
          rOrient = rZERO

          ! Put the particles around the well
          write (unit=IO_MessageBuffer, fmt="(""Releasing "", i5, "" particles from collector well ID "", i5)") iNParticles, iWellID
          write (unit=LU_TRACE, fmt="('CW0   ', i10)") iWellID
          call IO_MessageText(io)
          open(unit=iFwdLU, status='SCRATCH')
          do i = 1, iNParticles
            rAngle = rOrient + (float(i-1)+rHALF) * 2*rPI / float(iNParticles)
            rX1 = real(cZWell) + rCDistance * cos(rAngle)
            rY1 = aimag(cZWell) + rCDistance * sin(rAngle)
            rThisT0 = rT0
            rThisC0 = rC0
            write (unit=IO_MessageBuffer, fmt="("" Starting point: "", 2(e16.9, 1x), "" Direction: "", i2)") rX1, rY1, iDir
            call IO_MessageText(io)
            ! Here's where it gets interesting... Make a temporary file for the forward trace and trace into it.
            ! The forward trace happens at the aquifer bottom to ensure that all the water that _can_ make it
            ! to the collector get there.
            rBottom = rAQU_Base(io, aem%aqu, cmplx(rX1, rY1, AE_REAL))
            rZ0 = rZ0Save
            call TR0_Trace(io, aem, rX1, rY1, rBottom, rThisT0, rThisC0, -iDir, iFlag, &
                           iElementType, iElementID, iElementVtx, iFwdLU)
            print *, '   forward trace stopped at ', iFlag, iElementType, iElementID, iElementVtx
            ! Now, trace from where we hit the well, starting from time = 0
            rT0 = rZERO
            if (iElementType == ELEM_CW0 .and. iElementID == iWellID) then
              call TR0_Trace(io, aem, rX1, rY1, rZ0, rT0, rC0, iDir, iFlag, iElementType, iElementID, iElementVtx)
            end if
          end do
          close(unit=iFwdLU)
#endif
        case default
      end select
    end do
    ! Return the most recent error code found.
    return
  end subroutine TR0_Read


  subroutine TR0_Trace(io, aem, rX0, rY0, rZ0, rT0, rC0, iDir, iFlag, iElementType, iElementID, iElementVtx, iOptionalLU)
    ! Computes a single streamline starting at cZO, using the predictor-corrector method
    !
    ! All parameters receive the position at the end of the streamline.
    !
    type(AEM_DOMAIN), pointer :: aem
    real(kind=AE_REAL), intent(inout) :: rX0, rY0, rZ0, rT0, rC0
    integer(kind=AE_INT), intent(in) :: iDir
    integer(kind=AE_INT), intent(out) :: iFlag, iElementType, iElementID, iElementVtx
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iOptionalLU
    ! Locals
    integer(kind=AE_INT) :: iPrev, iNSteps, i1, i2, iE, iV
    integer(kind=AE_INT) :: iEType, iEString, iEVertex, iEFlag
    real(kind=AE_REAL) :: rBaseTol, rStagTol, rDZ20, rStep, rDeltaTime, rTime, rConc
    real(kind=AE_REAL) :: rBase, rThickness, rF, rZ, rAlphaT, rAlphaB, rStrength, rQ0, rQ1
    complex(kind=AE_REAL) :: cZO, cVO, cZN, cVN, cV, cZInt1, cZInt2, cZInt, cDis, cPrevDis
    complex(kind=AE_REAL) :: cZBefore, cZAfter, cENormal, cQBefore, cQAfter, cZFix
    complex(kind=AE_REAL), dimension(20) :: cPrevZ
    logical :: lTimeout, lFlag
    integer(kind=AE_INT) :: iMyLU
    !
    real(kind=AE_REAL), parameter :: rVERTEX_TOL = 1.0e-5_AE_REAL

    ! Select the LU number to be written...
    if (present(iOptionalLU)) then
      iMyLU = iOptionalLU
    else
      iMyLU = LU_TRACE
    end if

    ! Step size
    rStep = rTR0Step
    ! Tolerance used to determine whether an element has been reached
    rBaseTol = 1.1_AE_REAL * rTR0Step

    ! Make sure we're not on a vertex...
    cZN = cmplx(rX0, rY0, AE_REAL)
    if (lFWL_CheckPoint(io, aem%fwl, cZN, rVERTEX_TOL, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag)) cZN = cZFix
    if (lFDP_CheckPoint(io, aem%fdp, cZN, rVERTEX_TOL, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag)) cZN = cZFix
    if (lAQU_CheckPoint(io, aem%aqu, cZN, rVERTEX_TOL, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag)) cZN = cZFix

    ! Initialize
    cZO = cZN
    cVN = cAEM_Velocity(io, aem, cZN)
    cPrevDis = cAEM_Discharge(io, aem, cZN)
    rTime = rT0
    rConc = rC0
    iPrev = 1
    iNSteps = 0
    rStagTol = TR0_STAGNATION_TOLERANCE * rTR0Step
    lTimeout = .false.

    ! Is the particle in the aquifer?
    rBase = rAQU_Base(io, aem%aqu, cZO)
    rThickness = rAQU_SatdThickness(io, aem%aqu, cZO, real(cAEM_Potential(io, aem, cZO)))
    if (rZ0 > rBase + rThickness) then
      rZ = rBase + rThickness
    else if (rZ0 < rBase) then
      rZ = rBase
    else
      rZ = rZ0
    end if
    rF = (rZ-rBase)/rThickness

    ! Write the start record
    write (unit=iMyLU, fmt="(""START "", 5(e17.9, 1x), i10, 1x, e17.9)") rX0, rY0, rZ, rT0, rC0, iDir, abs(cPrevDis)

    do
      ! Write out the results of the previous step
      write (unit=iMyLU, fmt="(""      "", 5(e17.9, 1x))") cZN, rZ, rTime, rConc

      ! Prepare to take a step along the streamline
      cZO = cZN
      cVO = cVN
      cPrevZ(iPrev) = cZO
      iPrev = iPrev+1
      if (iPrev > 20) iPrev = mod(iPrev, 20)

      ! Abort if we are out of the window
      if (real(cZN) < rTR0WinX1 .or. real(cZN) > rTR0WinX2 .or. &
          aimag(cZN) < rTR0WinY1 .or. aimag(cZN) > rTR0WinY2) then
        iFlag = kTR0LeftWindow
        iElementType = 0
        iElementID = 0
        iElementVtx = 0
        exit
      end if

      ! Abort if this is a stagnation point
      if (abs(cVO) == rZERO) then
        iFlag = kTR0Stagnation
        iElementType = 0
        iElementID = 0
        iElementVtx = 0
        exit
      end if

      ! Abort if we have not moved in the past 20 steps
      if (iNSteps > 20) then
        i1 = iPrev
        i2 = iPrev-1
        if (i2 < 1) i2 = 20
        rDZ20 = abs(cPrevZ(i1) - cPrevZ(i2))
        if (rDZ20 < rStagTol) then
          iFlag = kTR0Stagnation
          iElementType = 0
          iElementID = 0
          iElementVtx = 0
          exit
        end if
      end if

      ! Abort if the maximum travel time is achieved
      if (lTimeout) then
        iFlag = kTR0Timeout
        iElementType = 0
        iElementID = 0
        iElementVtx = 0
        exit
      end if

      ! Abort if the head is below the aquifer bottom
      if (real(cAEM_Potential(io, aem, cZO)) <= rZERO) then
        iFlag = kTR0AquiferDry
        iElementType = 0
        iElementID = 0
        iElementVtx = 0
        exit
      end if

      ! Look to see if we're in a well...
      if (lFWL_CheckPoint(io, aem%fwl, cZN, rTR0Step*TR0_WLX_TERMINATION_FACTOR, cZFix, &
          rStrength, iEType, iEString, iEVertex, iEFlag)) then
        ! Check direction...
        if (sign(rONE, rStrength) == sign(rONE, real(iDir, AE_REAL))) then
          ! We are in a well!
          iFlag = kTR0HitElement
          iElementType = iEtype
          iElementID = iEString
          iElementVtx = 0
          exit
        end if
      end if

      ! PHEW. Now, proceed with the tracing...
      ! Adjust step size if necessary...  Note: all the xxxCheckPoint calls are made without regard to iDir
      if (lFWL_CheckPoint(io, aem%fwl, cZN, rTR0WellProx*rTR0Step, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag) .or. &
          lFDP_CheckPoint(io, aem%fdp, cZN, rTR0Prox*rTR0Step, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag) .or. &
          lAQU_CheckPoint(io, aem%aqu, cZN, rTR0Prox*rTR0Step, cZFix, rStrength, iEType, iEString, iEVertex, iEFlag) &
          ) then
        rStep = rTR0Frac * iDir * rTR0Step
      else
        rStep = iDir * rTR0Step
      end if

      ! *******************************************************************************************
      ! Determine the first guess for the complex coordinate of the next point on the streamline,
      ! Then check to see if we've hit something.
      cZN = cZO + rStep*cVO/abs(cVO)
      if (lFDP_CheckIntersection(io, aem%fdp, cZO, cZN, cZInt, cZBefore, cZAfter, &
          iEType, iEString, iEVertex, iEFlag, cENormal)) then
        ! We hit a dipole!
        cZN = cZBefore
      end if

      if (lAQU_CheckIntersection(io, aem%aqu, cZO, cZN, cZInt, cZBefore, cZAfter, &
          iEType, iEString, iEVertex, iEFlag, cENormal)) then
        ! We hit an edge element! Update the position and leave.
        cZN = cZBefore
      end if

      if (lFWL_CheckPoint(io, aem%fwl, cZN, rTR0Step*TR0_WLX_TERMINATION_FACTOR, cZFix, rStrength, &
          iEType, iEString, iEVertex, iEFlag)) then
        ! Check direction...
        if (sign(rONE, rStrength) == sign(rONE, real(iDir, AE_REAL))) then
          ! Move vertically with recharge until reaching the element
          if (TR0_RechargeMove(io, aem, cZO, cZN, rF, iDir, rZ, iFlag, iElementType, iElementID, iElementVtx)) then
            cV = cVO                        ! Needed for time computations(below)
            exit
          else
            ! We are in a well!
            cV = cVO                        ! Needed for time computations(below)
            iFlag = kTR0HitElement
            iElementType = iEType
            iElementID = iEString
            iElementVtx = iEVertex
            exit
          end if
        end if
      end if
      ! End of "first guess" checks
      ! *******************************************************************************************

      ! *******************************************************************************************
      ! Compute an improved location using a predictor-corrector strategy, then
      ! Compute the complex discharge Qx-iQy at the first guess
      cVN = cAEM_Velocity(io, aem, cZN)
      ! Estimate the average value of the complex discharge on the interval.
      cV = .5*(cVO+cVN)
      ! Obtain the second approximation of the point on the streamline.
      cZN = cZO + rStep*cV/abs(cV)

      if (lFDP_CheckIntersection(io, aem%fdp, cZO, cZN, cZInt, cZBefore, cZAfter, &
          iEType, iEString, iEVertex, iEFlag, cENormal)) then
        ! We hit a dipole!
        ! Move vertically with recharge until reaching the element
        cZN = cZBefore
        if (TR0_RechargeMove(io, aem, CZO, cZN, rF, iDir, rZ, iFlag, iElementType, iElementID, iElementVtx)) then
          cV = cVO
          exit
        end if
        ! Now, move vertically according to the amount of water removed by the element
        cQBefore = cAEM_Discharge(io, aem, cZBefore)
        cQAfter = cAEM_Discharge(io, aem, cZAfter)
        ! Does the flow reverse at the element(e.g. a line-sink)?
        if (iEType == ELEM_LS0 .or. iEType == ELEM_LS1 .or. iEType == ELEM_LS2) then
          if (sign(rONE, real(cENormal*cQBefore)) /= sign(rONE, real(cENormal*cQAfter))) then
            iFlag = kTR0HitElement
            iElementType = iEType
            iElementID = iEString
            iElementVtx = iEVertex
            rZ = rAQU_SatdThickness(io, aem%aqu, cZN, real(cAEM_Potential(io, aem, cZN), AE_REAL)) + &
                 rAQU_Base(io, aem%aqu, cZN)
            exit
          end if
          ! Now, jump the element to the point after
          cZN = cZAfter
          if (iDir > 0) then
            rF = rF * real(cENormal*cQBefore)/real(cENormal*cQAfter)
          else
            rF = rF * real(cENormal*cQAfter)/real(cENormal*cQBefore)
          end if
          if (rF > rONE) then
            iFlag = kTR0HitElement
            iElementType = iEType
            iElementID = iEString
            iElementVtx = iEVertex
            rZ = rAQU_SatdThickness(io, aem%aqu, cZN, real(cAEM_Potential(io, aem, cZN), AE_REAL)) + &
                 rAQU_Base(io, aem%aqu, cZN)
            exit
          end if
#ifndef __GPL__
        else if (iEType == ELEM_CW0) then
          iFlag = kTR0HitElement
          iElementType = iEType
          iElementID = iEString
          iElementVtx = iEVertex
          rZ = rAQU_SatdThickness(io, aem%aqu, cZN, real(cAEM_Potential(io, aem, cZN), AE_REAL)) + &
               rAQU_Base(io, aem%aqu, cZN)
          exit
#endif
        else
          cZN = cZAfter
        end if

      else if (lAQU_CheckIntersection(io, aem%aqu, cZO, cZN, cZInt, cZBefore, cZAfter, &
             iEType, iEString, iEVertex, iEFlag, cENormal)) then
        ! We hit an edge element! Update the position and leave.
        ! Move vertically with recharge until reaching the element
        cZN = cZBefore
        if (TR0_RechargeMove(io, aem, CZO, cZN, rF, iDir, rZ, iFlag, iElementType, iElementID, iElementVtx)) then
          cV = cVO
          exit
        end if
        ! This is the end of the pathline.
        iFlag = kTR0HitElement
        iElementType = iEType
        iElementID = iEString
        iElementVtx = iEVertex
        exit
      else
        ! Move vertically with recharge
        if (TR0_RechargeMove(io, aem, CZO, cZN, rF, iDir, rZ, iFlag, iElementType, iElementID, iElementVtx)) then
          cV = cVO
          exit
        end if
      end if
      ! End of predictor-corrector algorithm
      ! *******************************************************************************************

      ! Did we hit a well?
      if (lFWL_CheckPoint(io, aem%fwl, cZN, rTR0Step*TR0_WLX_TERMINATION_FACTOR, cZFix, rStrength, &
          iEType, iEString, iEVertex, iEFlag)) then
        ! Check direction...
        if (sign(rONE, rStrength) == sign(rONE, real(iDir, AE_REAL))) then
          ! We are in a well!
          iFlag = kTR0HitElement
          iElementType = iEType
          iElementID = iEString
          iElementVtx = iEVertex
          exit
        end if
      end if

      ! Compute the time-of-travel and report the results for this point
      rDeltaTime = abs(cZN-cZO)/abs(cV)
      !print *, rTime, abs(cV), rDeltaTime
      if (rTime+rDeltaTime > rTR0MaxTime) then
        if (abs(cZN-cZO) > rZERO) then
          cZN = cZO + (cZN-cZO)*(rTR0MaxTime-rTime)/rDeltaTime
          rTime = rTR0MaxTime
          lTimeout = .true.
        end if
      else
        rTime = rTime + rDeltaTime
      end if
      rConc = rC0 * exp(-rTR0Decay * (rTime-rT0)/rTR0Retard)

      ! Repeat the procedure for the next point.
      iNSteps = iNSteps + 1

    end do

    ! Write the end record, but first, look up the ID number for the element that was hit...
    select case (iElementType)
      case (ELEM_LS0)
        iElementID = iLS0_GetID(io, aem%ls0, iElementID)
      case (ELEM_LS1)
        iElementID = iLS1_GetID(io, aem%ls1, iElementID)
      case (ELEM_LS2)
        iElementID = iLS2_GetID(io, aem%ls2, iElementID)
      case (ELEM_WL0)
        iElementID = iWL0_GetID(io, aem%wl0, iElementID)
      case (ELEM_WL1)
        iElementID = iWL1_GetID(io, aem%wl1, iElementID)
#ifndef __GPL__
      case (ELEM_CW0)
        iElementID = iCW0_GetID(io, aem%cw0, iElementID)
#endif
    end select

    write (unit=iMyLU, fmt="(""END   "", 5(e17.9, 1x), 4i10)") cZN, rZ, rTime, rConc, &
           iFlag, iElementType, iElementID, iElementVtx
    ! Return the position and time at the endpoint
    rX0 = real(cZN)
    rY0 = aimag(cZN)
    rZ0 = rZ
    rT0 = rTime
    rC0 = rConc

    return
  end subroutine TR0_Trace


  function TR0_RechargeMove(io, aem, CZO, cZN, rF, iDir, rZ, iFlag, iElementType, iElementID, iElementVtx) result(lStop)
    ! Moves the fraction rF with recharge. If the particle leaves the aquifer along
    ! the path cZO-cZN, adjusts cZN and returns .true.
    ! [ PARAMETERS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZO
    complex(kind=AE_REAL), intent(inout) :: cZN
    real(kind=AE_REAL), intent(inout) :: rF, rZ
    integer(kind=AE_INT), intent(in) :: iDir
    integer(kind=AE_INT), intent(inout) :: iFlag, iElementType, iElementID, iElementVtx
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lStop
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rFN, rDS, rNT, rNB, rQ1

    rQ1 = abs(cAEM_Discharge(io, aem, cZN))
    rNT = rAEM_TopRecharge(io, aem, cZN)
    rNB = rAEM_BottomRecharge(io, aem, cZN)
    rDS = abs(cZN-cZO)
    if (iDir > 0) then
      rFN = (rF*rQ1 + rF*rNT*rDS + (rONE-rF)*rNB*rDS) / rQ1
    else
      rFN = (rF*rQ1 - rF*rNT*rDS - (rONE-rF)*rNB*rDS) / rQ1
    end if
    iElementType = 0
    iElementID = 0
    iElementVtx = 0

    if (rFN > rONE) then
      cZN = cZO + (cZN-cZO)*(rONE-rF)/(rFN-rF)
      iFlag = kTR0LeftTop
      rZ = rAQU_SatdThickness(io, aem%aqu, cZN, real(cAEM_Potential(io, aem, cZN), AE_REAL)) + &
           rAQU_Base(io, aem%aqu, cZN)
      rF = rONE
      lStop = .true.
    else if (rFN < rZERO) then
      cZN = cZO + (cZN-cZO)*(rZERO-rF)/(rFN-rF)
      iFlag = kTR0LeftTop
      rZ = rAQU_Base(io, aem%aqu, cZN)
      rF = rZERO
      lStop = .true.
    else
      rF = rFN
      rZ = rF * rAQU_SatdThickness(io, aem%aqu, cZN, real(cAEM_Potential(io, aem, cZN), AE_REAL)) + &
           rAQU_Base(io, aem%aqu, cZN)
      lStop = .false.
    end if

  end function TR0_RechargeMove

end module a_tr0

