module a_grid

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

  use u_constants
  use u_io
  use m_aem
  use m_aqu

  implicit none

  public

  !> Constants for grid generation
  !> These constants specify the type of data to be stored in a grid by GRI_MakeGrid()
  integer(kind=AE_INT), parameter :: GRID_HEAD = 1 
  integer(kind=AE_INT), parameter :: GRID_POTENTIAL = 2  
  integer(kind=AE_INT), parameter :: GRID_STREAMFUNCTION = 3  
  integer(kind=AE_INT), parameter :: GRID_QX = 4  
  integer(kind=AE_INT), parameter :: GRID_QY = 5  
  integer(kind=AE_INT), parameter :: GRID_VX = 6  
  integer(kind=AE_INT), parameter :: GRID_VY = 7  
  integer(kind=AE_INT), parameter :: GRID_FLOW = 8 
  integer(kind=AE_INT), parameter :: GRID_PRESSURE_HEAD = 9
  integer(kind=AE_INT), parameter :: GRID_VELOCITY_MAGNITUDE = 10

  !> Constants for file type specifications
  !> These determine the format of the grid file to be saved by GRI_MakeGrid()
  integer(kind=AE_INT), parameter :: FILE_SURFER = 1
  integer(kind=AE_INT), parameter :: FILE_MATLAB = 2
  integer(kind=AE_INT), parameter :: FILE_OPENDX = 3
  integer(kind=AE_INT), parameter :: FILE_XYZ = 4
  integer(kind=AE_INT), parameter :: FILE_ARC = 5

  real(kind=AE_REAL), parameter :: rSMALL_GRID = 1.0e-1

contains

  !> Processes model input for the GRI (grid) module
  !>
  !> Test doc
  subroutine GRI_Read(io, aem)

    ! Argument list
    type(AEM_DOMAIN), pointer :: aem 
    type(IO_STATUS), pointer :: io  
    ! Locals -- for Directive parsing
    integer(kind=AE_INT), parameter :: &
                          iEND=1001, iOPT=1002, iWIN=1003, iDIM=1004, &
                          iHEA=2000, iPOT=2001, iPSI=2002, iDIS=2003, iVEL=2004, &
                          iQ_X=2005, iQ_Y=2006, iFLO=2007, iPRH=2008 
    type(DIRECTIVE), dimension(13) :: dirDirectives = (/ &
                     DIRECTIVE(iEND, 'END'), &
                     DIRECTIVE(iOPT, 'OPT'), &
                     DIRECTIVE(iWIN, 'WIN'), &
                     DIRECTIVE(iDIM, 'DIM'), &
                     DIRECTIVE(iHEA, 'HEA'), &
                     DIRECTIVE(iPOT, 'POT'), &
                     DIRECTIVE(iPSI, 'PSI'), &
                     DIRECTIVE(iDIS, 'DIS'), &
                     DIRECTIVE(iVEL, 'VEL'), &
                     DIRECTIVE(iQ_X, 'Q_X'), &
                     DIRECTIVE(iQ_Y, 'Q_Y'), &
                     DIRECTIVE(iFLO, 'FLO'), &
                     DIRECTIVE(iPRH, 'PRH') /)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode         ! Opcode
    integer(kind=AE_INT) :: iStat           ! RunTime error code
    ! Placeholders for grid parameters
    complex(kind=AE_REAL) :: cLL, cUR, cLLi, cURi, cBcOrigin, cBcCheck
    integer(kind=AE_INT) :: iRes, iOption, ic, iBcQuadrant
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rMissingValue, rScale
    logical :: lFlag, lCheckActive
    character(len=132) :: sFile, sOption, sFixOption

    ! The remainder of this routine uses ifIO_InputRecord to process
    ! the model input file.
    iOption = FILE_SURFER
    rMissingValue = -901.0_AE_REAL
    rScale = rONE
    call IO_MessageText(io, "Entering GRI module")
    call IO_Assert(io, (associated(aem)), "GRI_Read: No AEM_DOMAIN object")
    ! Get the default settings
    call IO_GetWindow(io, cLL, cUR)
    iRes = 30
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "GRI_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all ModGRI "ifXXXRead" routines.
          ! Report the condition, but proceed as if EOD was found.
          call IO_Assert(io, .false., "GRI_Read: Unexpected EOF")
        case (iEND)
          ! END mark was found. Exit the file parser.
          call IO_MessageText(io, "Leaving GRI module")
          exit
        case (iDIM)
          !****************************************************************************
          ! Here for the DIM command -- dimension a grid
          !****************************************************************************
          iRes = iIO_GetInteger(io, "iRes", def=30)
          call IO_Assert(io, (iRes > 2), "GRI_Read: Illegal resolution")
        case (iWIN)
          !****************************************************************************
          ! Here for the WIN command -- Set the window
          !****************************************************************************
          cLLi = cIO_GetComplex(io, "cLL")
          cURi = cIO_GetComplex(io, "cRR")
          cLL = cmplx(min(real(cLLi), real(cURi)), min(aimag(cLLi),aimag(cURi)))
          cUR = cmplx(max(real(cLLi), real(cURi)), max(aimag(cLLi),aimag(cURi)))
        case (iHEA)
          !****************************************************************************
          ! Here for the HEA command -- generate and write a grid of heads
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_HEAD, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale)
        case (iPOT)
          !****************************************************************************
          ! Here for the POT command -- generate and write a grid of potentials
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_POTENTIAL, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale)
        case (iPSI)
          !****************************************************************************
          ! Here for the PSI command -- generate and write a grid of streamfunctions
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          cBcOrigin = cIO_GetCoordinate(io, "cBcOrigin", def=cHUGE)
          cBcCheck = cIO_GetCoordinate(io, "cBcCheck", def=cHUGE)
          iBcQuadrant = iIO_GetInteger(io, "iBcQuadrant", def=0)
          call GRI_MakeGrid(io, aem, GRID_STREAMFUNCTION, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale, cBcOrigin, cBcCheck, iBcQuadrant)
        case (iQ_X)
          !****************************************************************************
          ! Here for the QX command -- generate and write a grid of discharges
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_QX, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale)
        case (iQ_Y)
          !****************************************************************************
          ! Here for the QY command -- generate and write a grid of discharges
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_QY, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale)
        case (iFLO)
          !****************************************************************************
          ! Here for the FLO command -- generate and write a grid of integrated flows
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          cBcOrigin = cIO_GetCoordinate(io, "cBcOrigin", def=cZERO)
          call GRI_MakeGrid(io, aem, GRID_FLOW, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale, cBcOrigin)
        case (iPRH)
          !****************************************************************************
          ! Here for the PRH command -- generate and write a grid of pressure head
          !                             (only valid in PROFILE mode)
          !****************************************************************************
          call IO_Assert(io, io%lProfile, "Pressure-head grids are valid only in PROFILE models")
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_PRESSURE_HEAD, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale, cBcOrigin)
        case (iVEL)
          !****************************************************************************
          ! Here for the PRH command -- generate and write a grid of pressure head
          !                             (only valid in PROFILE mode)
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          call GRI_MakeGrid(io, aem, GRID_VELOCITY_MAGNITUDE, cLL, cUR, iRes, sFile, iOption, io%lProfile, &
               lCheckActive, rMissingValue, rScale, cBcOrigin)
        case (iOPT)
          !****************************************************************************
          ! Here for the QY command -- generate and write a grid of discharges
          !****************************************************************************
          sOption = sIO_GetField(io, "sOption", &
                    allowed=(/ "SURFER ", "MATLAB ", "OPENDX ", "XYZ    ", &
                               "ARC    ", "PROFILE", "MVALUE "/), &
                    def="SURFER", &
                    force_uppercase=.true.)
          lCheckActive = lIO_GetLogical(io, "lCheckActive", def=.true.)
          rMissingValue = rIO_GetReal(io, "rMissingValue", 1.70141E+38_AE_REAL)
          if (trim(sOption) == "SURFER") then
            iOption = FILE_SURFER
            call IO_MessageText(io, " Selecting SURFER output")
          else if (trim(sOption) == "MATLAB") then
            iOption = FILE_MATLAB
            call IO_MessageText(io, " Selecting Matlab output")
          else if (trim(sOption) == "OPENDX") then
            iOption = FILE_OPENDX
            call IO_MessageText(io, " Selecting OpenDX output")
          else if (trim(sOption) == "XYZ") then
            iOption = FILE_XYZ
            call IO_MessageText(io, " Selecting XYZ file output")
          else if (trim(sOption) == "ARC") then
            iOption = FILE_ARC
            call IO_MessageText(io, " Selecting ARC grid output")
            rScale = rIO_GetReal(io, "rScale", def=rONE)
            write (unit=IO_MessageBuffer, fmt=*) " Scale set to ", rScale
            call IO_MessageText(io)
          else
            call IO_Assert(io, (.false.), "GRI_Read: Illegal file format option")
          end if
      end select
    end do
    ! Return the most recent error code found.
    return
  end subroutine GRI_Read


  subroutine GRI_MakeGrid(io, aem, iFunc, cLL, cUR, iRes, sFile, iOption, lProfile, lCheckActive, &
               rMissingValue, rScale, cBcOrigin, cBcCheck, iBcQuadrant)
    !! subroutine MakeGrid
    !!
    !! Makes a grid of function 'func' for the AEM_DOMAIN 'aem', writing the
    !! results to sFile
    !!
    !! Arguments
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              AEM_DOMAIN to be used
    !!    (in)    integer :: iFunc
    !!              Function to be gridded(see constants defined above)
    !!                GRID_HEAD            Head(makes '.head.grd' file)
    !!                GRID_POTENTIAL       Potential(makes '.pot.grd' file)
    !!                GRID_STREAMFUNCTION  Streamfunction (io, makes '.psi.grd' file)
    !!                GRID_GRID_QX              X-discharge(makes '.GRID_QX.grd' file)
    !!                GRID_GRID_QY              Y-discharge(makes '.GRID_QY.grd' file)
    !!                GRID_GRID_VX              X-Velocity(makes 'GRID_VX.grd' file)
    !!                GRID_GRID_VY              Y-Velocity(makes 'GRID_VY.grd' file)
    !!                GRID_GRID_FLOW            Integrated flow(makes 'GRID_FLOW.grd' file)
    !!    (in)    complex :: cLL
    !!              Lower-left corner of the grid
    !!    (in)    complex :: cUR
    !!              Upper-rigth corner of the grid
    !!    (in)    integer :: iRes
    !!              Number of points along the longer axis
    !!    (in)    character :: sFile
    !!              File to be written, without extensions(they will be added
    !!              based on the function gridded; see above)
    !!    (in)    integer :: iOption
    !!              Option for file format
    !!    (in)    complex :: cBcOrigin
    !!              Origin point for branch-cut checks (see below), or origin point for GRD_FLOW
    !!    (in)    complex :: cBcCheck
    !!              If in PROFILE mode, check for branch cuts between grid points and cBcCheck
    !!    (in)    integer :: iBcQuadrant
    !!              In which quadrant (relative to cBcOrigin) do we check for branch cuts?
    !!              iBcQuadrant=0    Never
    !!              iBcQuadrant=1    Upper-right
    !!              iBcQuadrant=2    Upper-left
    !!              iBcQuadrant=3    Lower-left
    !!              iBcQuadrant=4    Lower-right
    !!
    !! Note: Uses the LU_GRID LU for output
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(in) :: iFunc
    complex(kind=AE_REAL), intent(in) :: cLL
    complex(kind=AE_REAL), intent(in) :: cUR
    integer(kind=AE_INT), intent(in) :: iRes
    character(len=*), intent(in) :: sFile
    integer(kind=AE_INT), intent(in) :: iOption
    logical, intent(in) :: lProfile
    logical, intent(in) :: lCheckActive
    real(kind=AE_REAL), intent(in) :: rMissingValue
    real(kind=AE_REAL), intent(in) :: rScale
    complex(kind=AE_REAL), intent(in), optional :: cBcOrigin
    complex(kind=AE_REAL), intent(in), optional :: cBcCheck
    integer(kind=AE_INT), intent(in), optional :: iBcQuadrant
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sFileName, sDataName
    character(len=4) :: sExt
    integer(kind=AE_INT) :: iStat, iResX, iResY, i, j, jj, iMyResX, iMyResY
    real(kind=AE_REAL) :: rSizeX, rSizeY, rDel, rDeltaX, rDeltaY, rX, rY, rBC
    real(kind=AE_REAL) :: rMinX, rMaxX, rMinY, rMaxY, rMyMinX, rMyMinY, rMyMaxX, rMyMaxY
    real(kind=AE_REAL), dimension(:, :), allocatable :: rGrid
    ! [ LOCAL CONSTANTS ]
    integer(kind=AE_INT), parameter :: BC_NEVER = 0
    integer(kind=AE_INT), parameter :: BC_UPPER_RIGHT = 1
    integer(kind=AE_INT), parameter :: BC_UPPER_LEFT = 2
    integer(kind=AE_INT), parameter :: BC_LOWER_LEFT = 3
    integer(kind=AE_INT), parameter :: BC_LOWER_RIGHT = 4

    call IO_Assert(io, (associated(aem)), 'MakeGrid: No AEM_DOMAIN object')
    call IO_Assert(io, (abs(cUR-cLL) > rSMALL_GRID), 'MakeGrid: Bad window')
    call IO_Assert(io, (iRes > 2), 'MakeGrid: Bad resolution')

    rSizeX = real(cUR-cLL)
    rSizeY = aimag(cUR-cLL)
    if ( .not. lProfile ) then
      if (rSizeX >= rSizeY) then
        iResX = iRes
        rDel = rSizeX/real(iRes-1, AE_REAL)
        iResY = int(rSizeY/rDel)+1
      else
        iResY = iRes
        rDel = rSizeY/real(iRes-1, AE_REAL)
        iResX = int(rSizeX/rDel)+1
      end if
      rDeltaX = rDel
      rDeltaY = rDel
    else
      call IO_Assert(io, iOption/=FILE_ARC, "Cannot use ARC grids in PROFILE mode")
      rDeltaX = rSizeX/real(iRes-1, AE_REAL)
      iResX = iRes
      rDeltaY = rSizeY/real(iRes-1, AE_REAL)
      iResY = iRes
    end if

    ! Make the real window span the center of the specified window
    rMinX = rHALF * (real(cLL+cUR) - (iResX-1)*rDeltaX)
    rMaxX = rMinX + (iResX-1)*rDeltaX
    rMinY = rHALF * (aimag(cLL+cUR) - (iResY-1)*rDeltaY)
    rMaxY = rMinY + (iResY-1)*rDeltaY

    ! Ugh. We need to change a few things for ARC grids...
    if (iOption == FILE_ARC .or. iOption == FILE_XYZ) then
      rMyMinX = rMinX+0.5*rDeltaX
      rMyMinY = rMinY+0.5*rDeltaY
      iMyResX = iResX-1
      iMyResY = iResY-1
      rMyMaxX = rMinX+iMyResX*rDeltaX
      rMyMaxY = rMinY+iMyResY*rDeltaY
    else
      rMyMinX = rMinX
      rMyMinY = rMinY
      iMyResX = iResX
      iMyResY = iResY
      rMyMaxX = rMinX+(iMyResX-1)*rDeltaX
      rMyMaxY = rMinY+(iMyResY-1)*rDeltaY
    end if
    ! Make space for the gridded data
    allocate(rGrid(iMyResY, iMyResX), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'MakeGrid: Allocation failed')

    ! Generate the grid
    select case (iFunc)
      case (GRID_HEAD)
        call IO_MessageText(io, "  Generating a grid of heads")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))
              if (lProfile) then
                if (rGrid(j, i) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_PRESSURE_HEAD)
        call IO_MessageText(io, "  Generating a grid of pressure-heads")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) - rY
            end if
          end do
        end do
      case (GRID_POTENTIAL)
        call IO_MessageText(io, "  Generating a grid of potentials")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = real(cAEM_Potential(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_STREAMFUNCTION)
        call IO_MessageText(io, "  Generating a grid of the streamfunction")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = aimag(cAEM_Potential(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) then
                  rGrid(j, i) = rMissingValue
                else
                  rGrid(j, i) = rGrid(j, i)
                  rBC = rZERO
                  select case (iBcQuadrant)
                    case (BC_UPPER_RIGHT)
                      if (rX > real(cBcOrigin) .and. rY > aimag(cBcOrigin)) then
                        rBC = rAQU_BranchCut(io, aem%aqu, (/ cmplx(rX, rY, AE_REAL), cBcCheck /))
                      end if
                    case (BC_UPPER_LEFT)
                      if (rX < real(cBcOrigin) .and. rY > aimag(cBcOrigin)) then
                        rBC = rAQU_BranchCut(io, aem%aqu, (/ cmplx(rX, rY, AE_REAL), cBcCheck /))
                      end if
                    case (BC_LOWER_LEFT)
                      if (rX < real(cBcOrigin) .and. rY < aimag(cBcOrigin)) then
                        rBC = rAQU_BranchCut(io, aem%aqu, (/ cmplx(rX, rY, AE_REAL), cBcCheck /))
                      end if
                    case (BC_LOWER_RIGHT)
                      if (rX > real(cBcOrigin) .and. rY < aimag(cBcOrigin)) then
                        rBC = rAQU_BranchCut(io, aem%aqu, (/ cmplx(rX, rY, AE_REAL), cBcCheck /))
                      end if
                  end select
                  rGrid(j, i) = rGrid(j, i) - rBC
                end if
              end if
            end if
          end do
        end do
      case (GRID_QX)
        call IO_MessageText(io, "  Generating a grid of Qx")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = real(cIO_WorldDischarge(io, &
                                 cAEM_Discharge(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_QY)
        call IO_MessageText(io, "  Generating a grid of Qy")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = aimag(cIO_WorldDischarge(io, &
                                 cAEM_Discharge(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_VX)
        call IO_MessageText(io, "  Generating a grid of Vx")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = real(cIO_WorldDischarge(io, &
                                 cAEM_Velocity(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_VY)
        call IO_MessageText(io, "  Generating a grid of Vy")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = aimag(cIO_WorldDischarge(io, &
                                  cAEM_Velocity(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_VELOCITY_MAGNITUDE)
        call IO_MessageText(io, "  Generating a grid of velocities")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
              rGrid(j, i) = abs(cIO_WorldDischarge(io, &
                                 cAEM_Velocity(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
      case (GRID_FLOW)
        call IO_MessageText(io, "  Generating a grid of integrated flows")
        do j = 1, iMyResY
          rY = rMyMinY + (j-1)*rDeltaY
          do i = 1, iMyResX
            rX = rMyMinX + (i-1)*rDeltaX
            if (lCheckActive .and. .not. lAQU_CheckActive(io, aem%aqu, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)))) then
              rGrid(j,i) = rMissingValue
            else
!              rGrid(j, i) = rIO_WorldLength(io, &
!                                  rAEM_Flow(io, aem, (/ cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)), cBcOrigin /)), &
!                                  cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))-cBcOrigin)
              rGrid(j, i) = rAEM_Flow(io, aem, (/ cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL)), cBcOrigin /))
              if (lProfile) then
                if (rAEM_Head(io, aem, cIO_LocalCoords(io,cmplx(rX, rY, AE_REAL))) < rY) rGrid(j, i) = rMissingValue
              end if
            end if
          end do
        end do
    end select

    ! Write the grid file
    if (iOption == FILE_SURFER) then
      sExt = ".grd"
    else if (iOption == FILE_MATLAB) then
      sExt = ".m"
    else if (iOption == FILE_OPENDX) then
      sExt = ".dx"
    else if (iOption == FILE_XYZ) then
      sExt = ".dat"
    else if (iOption == FILE_ARC) then
      sExt = ".grd"
    end if
    select case (iFunc)
      case (GRID_HEAD)
        sFileName = trim(sFile) // "_head" // sExt
        sDataName = 'head'
      case (GRID_PRESSURE_HEAD)
        sFileName = trim(sFile) // "_pressure_head" // sExt
        sDataName = 'pressure_head'
      case (GRID_POTENTIAL)
        sFileName = trim(sFile) // "_potential" // sExt
        sDataName = 'potential'
      case (GRID_STREAMFUNCTION)
        sFileName = trim(sFile) // "_psi" // sExt
        sDataName = 'streamfunction'
      case (GRID_QX)
        sFileName = trim(sFile) // "_qx" // sExt
        sDataName = 'qx'
      case (GRID_QY)
        sFileName = trim(sFile) // "_qy" // sExt
        sDataName = 'qy'
      case (GRID_VX)
        sFileName = trim(sFile) // "_vx" // sExt
        sDataName = 'vx'
      case (GRID_VY)
        sFileName = trim(sFile) // "_vy" // sExt
        sDataName = 'vy'
      case (GRID_VELOCITY_MAGNITUDE)
        sFileName = trim(sFile) // "_vel" // sExt
        sDataName = 'vy'
      case (GRID_FLOW)
        sFileName = trim(sFile) // "_flow" // sExt
        sDataName = 'flow'
    end select

    open(unit=LU_GRID, file = trim(sFileName), iostat=iStat)
    call IO_Assert(io, (iStat == 0), 'MakeGrid: Could not open file')
    call IO_MessageText(io, "  Writing grid to " // trim(sFileName))

    if (iOption == FILE_SURFER) then
      write (unit=LU_GRID, &
             fmt="('DSAA')" &
             )
      write (unit=LU_GRID, &
             fmt="(2(1x, i5))" &
             ) iMyResX, iMyResY
      write (unit=LU_GRID, &
             fmt="(2(1x, g15.8))" &
             ) rMyMinX, rMyMaxX
      write (unit=LU_GRID, &
             fmt="(2(1x, g15.8))" &
             ) rMyMinY, rMyMaxY
      write (unit=LU_GRID, &
             fmt="(2(1x, g15.8))" &
             ) minval(rGrid,rGrid/=rMissingValue), maxval(rGrid,rGrid/=rMissingValue)
      do j = 1, iMyResY
        write (unit=LU_GRID, &
               fmt=* &
               ) rGrid(j, :)
      end do
      ! That's it
      close(unit=LU_GRID)

    else if (iOption == FILE_MATLAB) then
      write (unit=LU_GRID, &
             fmt="('# Created by ModAEM 1.4')" &
             )
      write (unit=LU_GRID, &
             fmt="('# name: x')" &
             )
      write (unit=LU_GRID, &
             fmt="('# type: matrix')" &
             )
      write (unit=LU_GRID, &
             fmt="('# rows: ', i10)" &
             ) iMyResX
      write (unit=LU_GRID, &
             fmt="('# columns: 1')" &
             )
      do i = 1, iMyResX
        write (unit=LU_GRID, &
               fmt=* &
               ) rMyMinX + (i-1)*rDeltaX
      end do

      write (unit=LU_GRID, &
             fmt="('# Created by ModAEM 1.4')" &
             )
      write (unit=LU_GRID, &
             fmt="('# name: y')" &
             )
      write (unit=LU_GRID, &
             fmt="('# type: matrix')" &
             )
      write (unit=LU_GRID, &
             fmt="('# rows: ', i10)" &
             ) iMyResY
      write (unit=LU_GRID, &
             fmt="('# columns: 1')" &
             )
      do j = 1, iMyResY
        write (unit=LU_GRID, &
               fmt=* &
               ) rMyMinY + (j-1)*rDeltaY
      end do

      write (unit=LU_GRID, &
             fmt="('# Created by ModAEM 1.4')" &
             )
      write (unit=LU_GRID, &
             fmt="('# name: ', a20)" &
             ) sDataName
      write (unit=LU_GRID, &
             fmt="('# type: matrix')" &
             )
      write (unit=LU_GRID, &
             fmt="('# rows: ', i10)" &
             ) iMyResX
      write (unit=LU_GRID, &
             fmt="('# columns: ', i10)" &
             ) iMyResY
      do i = 1, iMyResX
        write (unit=LU_GRID, &
               fmt=* &
               ) rGrid(:, i)
      end do
      ! That's it
      close(unit=LU_GRID)

    else if (iOption == FILE_OPENDX) then
      write (unit=LU_GRID, &
             fmt="('grid_x = ', i5)" &
             ) iMyResX
      write (unit=LU_GRID, &
             fmt="('grid_y = ', i5)" &
             ) iMyResY
      write (unit=LU_GRID, &
             fmt="('origin_x = ', g12.5)" &
             ) rMyMinX
      write (unit=LU_GRID, &
             fmt="('origin_y = ', g12.5)" &
             ) rMyMinY
      write (unit=LU_GRID, &
             fmt="('step_x = ', g12.5)" &
             ) rDeltaX
      write (unit=LU_GRID, &
             fmt="('step_y = ', g12.5)" &
             ) rDeltaY
      do j = iMyResY, 1, -1
        write (unit=LU_GRID, &
               fmt=* &
               ) rGrid(j, :)
      end do
      ! That's it
      close(unit=LU_GRID)

    else if (iOption == FILE_XYZ) then
      do j = 1, iMyResY
        do i = 1, iMyResX
          write (unit=LU_GRID, &
                 fmt="(3(g12.5, 1x))" &
                 ) rMyMinX+(i-1)*rDeltaX, rMyMinY+(j-1)*rDeltaY, rGrid(j, i)
        end do
      end do
      ! That's it
      close(unit=LU_GRID)

    else if (iOption == FILE_ARC) then
      write (unit=LU_GRID, &
             fmt="('ncols', i5)" &
             ) iMyResX
      write (unit=LU_GRID, &
             fmt="('nrows', i5)" &
             ) iMyResY
      write (unit=LU_GRID, &
             fmt="('xllcorner', g16.9)" &
             ) rMinX
      write (unit=LU_GRID, &
             fmt="('yllcorner', g16.9)" &
             ) rMinY
      write (unit=LU_GRID, &
             fmt="('cellsize', g16.9)" &
             ) rDeltaX
      write (unit=LU_GRID, &
             fmt="('nodata_value', g16.9)" &
             ) rMissingValue
      do j = 1, iMyResY
        do i = 1, iMyResX
          write (unit=LU_GRID, &
                 fmt='(i10, 1x)', &
                 advance='NO' &
                 ) int(rGrid(iMyResY-j+1, i)*rScale)
        end do
      end do
      write (unit=LU_GRID, &
             fmt=*)
      ! That's it
      close(unit=LU_GRID)

    end if
    deallocate(rGrid, stat = iStat)
    call IO_Assert(io, (iStat == 0), 'MakeGrid: Deallocation failed')

    return
  end subroutine GRI_MakeGrid

end module a_grid
