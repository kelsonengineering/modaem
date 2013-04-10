program modaem

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

  !
  ! This is a Fortran MODULE for general-purpose AEM functions. The regional aquifer
  ! properties, conversion from head-to-potential, potential/velocity computations, and
  ! reporting functions are included here.
  !

  use u_constants
  use u_io
  use m_aqu
  use m_wl0
  use m_wl1
  use m_ls0
  use m_ls1
  use m_hb0
  use m_aem
  use m_inq
  use m_ext
  use a_grid
  use a_tr0
  use a_stdio
  use a_obs

  implicit none

  ! [ RETURN VALUE ]
  type(AEM_DOMAIN), pointer :: aem
  ! [ LOCALS ]
  integer(kind=AE_INT), parameter :: &
                          iEOD = 1000, iAEM = 1001, iRPT = 1002, iSOL = 1003, iGRI = 1004, iOBS = 1005, &
                          iTR0 = 1006, iEXT = 1007, iPRE = 1008, iSAV = 1009, iINQ = 1010, iDDN = 1011, &
                          iHEA = 2000, iPOT = 2001, iDIS = 2002, iFLO = 2003, iGRA = 2004, iLAP = 2005, &
                          iVEL = 2006, iRCH = 2007, iTHK = 2008, iIFC = 2009, iSTD=2010

  type(DIRECTIVE), dimension(23) :: dirDirectives = (/ &
                     DIRECTIVE(iEOD, 'EOD'), &
                     DIRECTIVE(iAEM, 'AEM'), &
                     DIRECTIVE(iRPT, 'RPT'), &
                     DIRECTIVE(iSOL, 'SOL'), &
                     DIRECTIVE(iDDN, 'DDN'), &
                     DIRECTIVE(iGRI, 'GRI'), &
                     DIRECTIVE(iOBS, 'OBS'), &
                     DIRECTIVE(iTR0, 'TR0'), &
                     DIRECTIVE(iEXT, 'EXT'), &
                     DIRECTIVE(iPRE, 'PRE'), &
                     DIRECTIVE(iSAV, 'SAV'), &
                     DIRECTIVE(iINQ, 'INQ'), &
                     DIRECTIVE(iHEA, 'HEA'), &
                     DIRECTIVE(iPOT, 'POT'), &
                     DIRECTIVE(iDIS, 'DIS'), &
                     DIRECTIVE(iVEL, 'VEL'), &
                     DIRECTIVE(iRCH, 'RCH'), &
                     DIRECTIVE(iTHK, 'THK'), &
                     DIRECTIVE(iFLO, 'FLO'), &
                     DIRECTIVE(iGRA, 'GRA'), &
                     DIRECTIVE(iLAP, 'LAP'), &
                     DIRECTIVE(iIFC, 'IFC'), &
                     DIRECTIVE(iSTD, 'STD') &
                     /)
  ! Directives
  character(len=255) :: sCommandLine
  character(len=132) :: sBaseName
  character(len=132) :: sFileName
  integer(kind=AE_INT) :: iOpCode
  integer(kind=AE_INT) :: iStat
  integer(kind=AE_INT) :: iretval
  integer(kind=AE_INT) :: iNIter, iNPolishIter
  integer(kind=AE_INT) :: i
  logical :: lDebug, lCrash, lProfile
  ! Placeholders for function module test calls
  complex(kind=AE_REAL) :: cZ1, cZ2, cZC, cZE1, cZE2, cDis, cPC, cPN, cPS, cPE, cPW, cVel
  real(kind=AE_REAL) :: rR, rTol, rDPhiDX, rDPhiDY, rLapl, rThk, rRch, rRelaxation
  complex(kind=AE_REAL) :: cZOrigin
  real(kind=AE_REAL) :: rScale, rRotation, rKHoverKV
  character(len=255) :: sFname, sMode
  type(IO_STATUS), pointer :: io

  ! First, open the names file(with directives), modaem.nam
  open(unit=LU_SCRATCH, file = 'modaem.nam', status='old', iostat=iStat)
  call IO_Assert(io, (iStat == 0), "ModAEM: Could not open names file modaem.nam for input")
  !** FUTURE: Put processing directives(e.g. SQL connections and what-not) here!
  ! Retrieve the AEM file name from modaem.nam
  read (unit=LU_SCRATCH, fmt=*, iostat=iStat) sBaseName
  call IO_Assert(io, (iStat == 0), "ModAEM: Could not read AEM base file name")
  ! Is there a coordinate conversion record?
  read (unit=LU_SCRATCH, fmt=*, iostat=iStat) cZOrigin, rScale, rRotation, lDebug, lCrash
  if (iStat /= 0) then
    print *, "No reprojection information was provided"
    cZOrigin = cZERO
    rScale = rONE
    rRotation = rZERO
    lDebug = .false.
    lCrash = .true.
  end if
  ! Is there a profile modeling record?
  read (unit=LU_SCRATCH, fmt=*, iostat=iStat) lProfile, rKHoverKV
  if (iStat /= 0) then
    print *, "No profile modeling information was provided"
    lProfile = .false.
    rKHoverKV = rONE
  end if
  close(unit=LU_SCRATCH)

  ! First Create theIO_STATUS Object
  io => IO_Create(cZOrigin, rScale, rRotation, lDebug, lCrash, lProfile, rKHoverKV)
  ! Turn off the debug flag
  lDebug = .false.
  ! Wire up the files...
  call IO_OpenAll(io, sBaseName)
  nullify(aem)

  ! Here we go!
  do
    call IO_InputRecord(io, dirDirectives, iOpCode)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation.
        call IO_Assert(io, .false., "AEM_Read: I/O Error")
      case (kOpFileEOF)
        ! EOF is unexpected for all ModAEM "ifXXXRead" routines.
        call IO_Assert(io, .false., "AEM_Read: Unexpected EOF")
      case (kOpData)
        ! A data line was found. The main ModAEM module has no default
        ! data lines. Report the condition.
        call IO_Assert(io, .false., "AEM_Read: Unexpected data record")
      case (iEOD)
        ! EOD mark was found. Exit the file parser.
        call IO_MessageText(io, "EOD Encountered -- Job complete")
        exit
      case (iAEM)
        ! Populate the AEM_DOMAIN object from the input file according to these sizes...
        call IO_Assert(io, (.not. associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has already been created")
        aem => AEM_Create(io)
        call AEM_Read(io, aem)
      case (iRPT)
        ! Generates a report
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        call HTML_Start()
        call HTML_Header('ModAEM 1.8.0-dev', 1)
        !call AEM_Report(io, aem)
        if (associated(aem%aqu)) call AQU_Report(io, aem%aqu)
        if (associated(aem%wl0)) call WL0_Report(io, aem%wl0)
        if (associated(aem%wl1)) call WL1_Report(io, aem%wl1)
        if (associated(aem%pd0)) call PD0_Report(io, aem%pd0)
        if (associated(aem%ls0)) call LS0_Report(io, aem%ls0)
        if (associated(aem%ls1)) call LS1_Report(io, aem%ls1)
        if (associated(aem%ls2)) call LS2_Report(io, aem%ls2, aem%aqu)
        if (associated(aem%hb0)) call HB0_Report(io, aem%hb0)
        if (associated(aem%as0_top)) call AS0_Report(io, aem%as0_top, '(aquifer top)')
        if (associated(aem%as0_bottom)) call AS0_Report(io, aem%as0_bottom, '(aquifer bottom)')
#ifndef __GPL__
        if (associated(aem%cw0)) call CW0_Report(io, aem%cw0, aem%aqu)
#endif
        if (io%lDebug) then
          if (associated(aem%mat)) call MAT_Report(io, aem%mat, 'mat', .false.)
          if (associated(aem%fwl)) call FWL_Report(io, aem%fwl)
          if (associated(aem%fpd)) call FPD_Report(io, aem%fpd)
          if (associated(aem%fdp)) call FDP_Report(io, aem%fdp)
        end if

        call HTML_End()
      case (iSAV)
        ! Save the model solution information to a .pre file
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        sFname = sIO_GetField(io, "sFname")
        sMode = sIO_GetField(io, "sMode", def="BINARY", allowed=(/"ASCII ", "BINARY"/))
        if (trim(sMode) == "BINARY") then
          call AEM_Save(io, aem, sFname, IO_MODE_ASCII)
        else
          call AEM_Save(io, aem, sFname, IO_MODE_ASCII)
        end if
      case (iPRE)
        ! Precondition the model solution information from a .pre file
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        sMode = sIO_GetField(io, "sMode", def="BINARY", allowed=(/"ASCII ", "BINARY"/))
        if (trim(sMode) == "BINARY") then
          call AEM_Load(io, aem, sFname, IO_MODE_ASCII)
        else
          call AEM_Load(io, aem, sFname, IO_MODE_ASCII)
        end if
      case (iSOL)
        ! Enter the SOLver module
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        iNIter = iIO_GetInteger(io, "iNIter", def=4, minimum=1)
        rRelaxation = rIO_GetReal(io, "rRelaxation", def=rONE, minimum=0.1_AE_REAL, maximum=rONE)
        iNPolishIter = iIO_GetInteger(io, "iNPolishIter", def=1, minimum=0)
        call AEM_Solve(io, aem, iNIter, iNPolishIter, rRelaxation)
      case (iDDN)
        ! Switch on the "drawdown" elements
        call AEM_EnableDrawdown(io, aem)
      case (iGRI)
        ! Enter the GRId module
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        call GRI_Read(io, aem)
      case (iINQ)
        ! Enter the inquiry module
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        call INQ_Read(io, aem)
      case (iOBS)
        ! Enter the observations module
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        call OBS_Read(io, aem)
      case (iEXT)
        ! Enter the extract module
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        call EXT_Read(io, aem)
      case (iTR0)
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        ! Enter the TR0 2-D tracing module
        call IO_Assert(io, (associated(aem)), "AEM_Read: No DIM statement")
        call TR0_Read(io, aem)
      case (iHEA)
        ! Report the head to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io, 'cZ1')
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Head       at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, rAEM_Head(io, aem, cZ1)
        call IO_ErrorText(io)
      case (iPOT)
        ! Report the complex potential to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io, 'cZ1')
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Potential  at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, cAEM_Potential(io, aem, cZ1)
        call IO_ErrorText(io)
      case (iGRA)
        ! Report the estimated gradient in Phi
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io, 'cZ1')
        rTol = rIO_GetReal(io, 'rTol', def=0.0001_AE_REAL)
        rDPhiDX = (real(cAEM_Potential(io, aem, cZ1-cmplx(rTol, rZERO, AE_REAL))) - &
                  real(cAEM_Potential(io, aem, cZ1+cmplx(rTol, rZERO, AE_REAL)))) / (rTWO*rTol)
        rDPhiDY = (real(cAEM_Potential(io, aem, cZ1-cmplx(rZERO, rTol, AE_REAL))) - &
                  real(cAEM_Potential(io, aem, cZ1+cmplx(rZERO, rTol, AE_REAL)))) / (rTWO*rTol)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Gradient   at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, cIO_WorldDischarge(io, cmplx(rDPhiDX, rDPhiDY, AE_REAL))
        call IO_ErrorText(io)
      case (iLAP)
        ! Report the estimated laplacian
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        rTol = rIO_GetReal(io,"rTol",def=0.0001_AE_REAL)
        cPE = cAEM_Potential(io, aem, cZ1+cmplx(rTol, rZERO, AE_REAL))
        cPW = cAEM_Potential(io, aem, cZ1-cmplx(rTol, rZERO, AE_REAL))
        cPN = cAEM_Potential(io, aem, cZ1+cmplx(rZERO, rTol, AE_REAL))
        cPS = cAEM_Potential(io, aem, cZ1-cmplx(rZERO, rTol, AE_REAL))
        cPC = cAEM_Potential(io, aem, cZ1)
        rLapl = (cPE+cPW+cPN+cPS - rFOUR*cPC) / (rTol*rTol)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Laplacian   at "", d13.6, 1x, d13.6, "" = "", d13.6)" &
               ) cZ1, rLapl
        call IO_ErrorText(io)
      case (iDIS)
        ! Report the complex discharge to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        cDis = cAEM_Discharge(io, aem, cZ1)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Discharge  at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, cIO_WorldDischarge(io, cDis)
        call IO_ErrorText(io)
      case (iVEL)
        ! Report the velocity to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        cVel = cAEM_Velocity(io, aem, cZ1)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Velocity   at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, cIO_WorldDischarge(io, cVel)
        call IO_ErrorText(io)
      case (iRCH)
        ! Report the recharge to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        rRch = rAEM_Recharge(io, aem, cZ1)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Recharge   at "", d13.6, 1x, d13.6, "" = "", d13.6)" &
               ) cZ1, rRch
        call IO_ErrorText(io)
      case (iTHK)
        ! Report the saturated thickness to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        rThk = rAEM_SatdThick(io, aem, cZ1)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Thickness at "", d13.6, 1x, d13.6, "" = "", d13.6, 1x, d13.6)" &
               ) cZ1, rThk
        call IO_ErrorText(io)
      case (iFLO)
        ! Report the total flow to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        cZ2 = cIO_GetCoordinate(io,"CZ2",def=cZ1)
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Flow       between "", d13.6, 1x, d13.6, "" and "", d13.6, 1x, d13.6, "" = "", d13.6)" &
               ) cZ1, cZ2, rIO_WorldLength(io, rAEM_Flow(io, aem, (/cZ1, cZ2/)), cZ2-cZ1)
        call IO_ErrorText(io)
      case (iIFC)
        ! Report the interface elevation to the error LU
        call IO_Assert(io, (associated(aem)), &
             "AEM_Read: the AEM_DOMAIN has not been created")
        cZ1 = cIO_GetCoordinate(io,"CZ1")
        write (unit=IO_MessageBuffer, &
               fmt="("" >> Interface elev  at "", d13.6, 1x, d13.6, "" = "", d13.6)" &
               ) cZ1, rAEM_InterfaceElevation(io, aem, cZ1)
        call IO_ErrorText(io)
      case (iSTD)
        ! Enter the STDIO module
        call STD_IO(io, aem)
      case default
    end select
  end do

  !**pd Destroy the aem object and all member objects
  if (associated(aem)) call AEM_Destroy(io, aem)

  call exit(0)
end program modaem
