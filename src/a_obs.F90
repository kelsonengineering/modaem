module a_obs

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

  ! Observation module
  ! V.A.Kelson  04/06/2006
  !
  ! This module allows for the computation of calibration observations
  !
  ! Input format:
  !
  !   # Enter the observations module, with up to ngrp groups
  !   OBS ngrp
  !
  !     GRP group-name group-id group-weight
  !
  !         Builds a group with the specified group name and weight, with
  !         ID number as specified
  !
  !     FIL file-name [delimiter]
  !
  !         Opens an output file. No extension is appended. Output will be
  !         written to that file until another FIL directive is encountered.
  !
  !         If provided, [delimiter] specifies the delimiter in the output file
  !         Options are: SPACE, TAB, COMMA
  !
  !     HEA tgt-name (X,Y) tgt-value obs-weight group-id
  !
  !         Reports the head at (x,y) as target-name.
  !
  !     DIS tgt-name (X,Y) tgt-value obs-weight group-id
  !
  !         Reports the magnitude of the discharge vector at (x,y) as target-name.
  !
  !     VEL tgt-name (X,Y) tgt-value obs-weight group-id
  !
  !         Reports the magnitude of the discharge vector at (x,y) as target-name.
  !
  !     FLO tgt_name tgt-value obs-weight group-id
  !         (x1,y1)
  !         (x2,y2)
  !         ...
  !         (xN,yN)
  !     END
  !
  !         Computes the integrated flux across the path as target-name.
  !
  !     GR2 tgt_name (x1,y1) (x2,y2) tgt-value obs-weight group-id
  !
  !         Reports the two-point gradient along (x1,y1)-(x2,y2) as target-name.
  !
  !     GR3 tgt_name (x1,y1) (x2,y2) (x3,y3) tgt-value obs-weight group-id
  !
  !         Reports the triangulated gradient among [(xi,yi)] as target-name.
  !
  !     GAG tgt-name ls2-str-id tgt-value obs-weight group-id
  !
  !         Reports the streamflow at the end point of LS2 string with the
  !         specified ID as target-name.
  !

  use u_constants
  use u_io
  use m_aem
  use m_ls2

  implicit none

  private

  public :: OBS_Read

  type :: OBS_GROUP
    ! Information about an observation group
    character(len=255) :: sName
    integer(kind=AE_INT) :: iID
    real(kind=AE_REAL) :: rWeight
  end type OBS_GROUP

  ! Global variables for this module

  ! True if the output file is open
  logical, private, save :: lOBSFileOpen = .false.


contains


  function OBS_GetGroup(io, grp_array, iNGrp, iID) result(grp)
    !! Returns the group with ID iGrpID from an array of OBS_GROUP objects,
    !! using elements 1..iNGrp
    ! [ PARAMETERS ]
    type(IO_STATUS), pointer :: io
    type(OBS_GROUP), dimension(:), pointer :: grp_array
    integer(kind=AE_INT), intent(in) :: iNGrp
    integer(kind=AE_INT), intent(in) :: iID
    ! [ RETURN VALUE ]
    type(OBS_GROUP), pointer :: grp
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    nullify(grp)
    do i = 1, iNGrp
      if ( grp_array(i)%iID == iID ) then
        grp => grp_array(i)
        return
      end if
    end do

    return
  end function OBS_GetGroup


  subroutine OBS_Read(io, aem)
    !! This reads and processes the input commands for the OBS package
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    integer(kind=AE_INT), parameter :: OP_GRP=1001, OP_FIL=1002, OP_HEA=1003, OP_DIS=1004, &
                            OP_FLO=1005, OP_GR2=1006, OP_GR3=1007, OP_GAG=1008, OP_VEL=1009, &
                            OP_END=1101
    type(DIRECTIVE), dimension(10), parameter :: dirDirectives = &
                       (/ DIRECTIVE(OP_END, "END"), &
                          DIRECTIVE(OP_GRP, "GRP"), &
                          DIRECTIVE(OP_FIL, "FIL"), &        ! Done
                          DIRECTIVE(OP_HEA, "HEA"), &        ! Done
                          DIRECTIVE(OP_DIS, "DIS"), &        ! Done
                          DIRECTIVE(OP_VEL, "VEL"), &        ! Done
                          DIRECTIVE(OP_FLO, "FLO"), &        ! Done
                          DIRECTIVE(OP_GR2, "GR2"), &
                          DIRECTIVE(OP_GR3, "GR3"), &
                          DIRECTIVE(OP_GAG, "GAG") &
                       /)
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iOpCode, iStat
    ! Placeholders for tracing parameters
    character(len=255) :: sName, sFile, sGrp
    character(len=10) :: sOpt
    character(len=1) :: sDelim
    integer(kind=AE_INT) :: iMaxGrp, iNGrp, iGrpID, iGageID, iGID
    real(kind=AE_REAL) :: rTgtValue, rTgtWeight, rValue, rWeight
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3, cValue
    type(OBS_GROUP), dimension(:), pointer :: ObsGroups
    type(OBS_GROUP), pointer :: grp

    ! Let's get ready to R U M B L E!
    call IO_MessageText(io, "Reading OBS module input")

    ! Clear the status flags
    lOBSFileOpen = .false.
    call IO_MessageText(io, "Entering observations module OBS")
    call IO_Assert(io, (associated(aem)), "OBS_Read: No AEM_DOMAIN object")

    ! How many groups? Note that we always create a group called "GLOBAL"
    ! with group weight 1.0 for observations that specify no group-id
    iMaxGrp = iIO_GetInteger(io, "iMaxGrp", def=1, minimum=1)
    allocate(ObsGroups(iMaxGrp), stat=iStat)
    call IO_Assert(io, iStat == 0, "Could not allocate observation groups")
    ! Note -- we initialize all of them to the default, but only the
    ! portion (1, iNGrp) is actually used.
    iNGrp = 1
    ObsGroups(iNGrp) = OBS_GROUP("GLOBAL", 0, 1.0)

    ! Uses ifIOInputRecord to process the model input file.
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpFileEOF)
          ! EOF is unexpected for all ModGRI "ifXXXRead" routines.
          call IO_Assert(io, .false., "OBS_Read: Unexpected EOF")
        case (OP_END)
          ! END mark was found. Exit the file parser.
          if (lOBSFileOpen) then
            close(unit=LU_SCRATCH)
            lOBSFileOpen = .false.
          end if
          exit
        case (OP_FIL)
          !****************************************************************************
          ! Here for the FIL command -- open the output file
          !****************************************************************************
          sFile = sIO_GetField(io, "sFile")
          sOpt = sIO_GetField(io, "sDelim", force_uppercase=.true., &
                 allowed=(/ "SPACE", "TAB  ", "COMMA" /), &
                 def="COMMA" )
          if (sOpt == "SPACE") then
            sDelim = " "
          else if (sDelim == "TAB") then
            sDelim = char(9)
          else
            sDelim = ","
          end if
          ! If the output file is open, close it now.
          if (lOBSFileOpen) then
            close(unit=LU_SCRATCH)
            lOBSFileOpen = .false.
          end if
          ! Open the new output file.  If the open fails, write a message.
          sFile = trim(sFile)
          open(unit=LU_SCRATCH, file = sFile, action = "WRITE", status="REPLACE", iostat=iStat, recl=200)
          call IO_Assert(io, (iStat == 0), "OBS_Read: Open failed")
          lOBSFileOpen = .true.
        case (OP_GRP)
          call IO_Assert(io, iNGrp < iMaxGrp, "No space for additional groups")
          iNGrp = iNGrp + 1
          sGrp = sIO_GetField(io, "sGrp")
          iGID = iIO_GetInteger(io, "iGID")
          rWeight = rIO_GetReal(io, "rWeight")
          ObsGroups(iNGrp) = OBS_GROUP(sGrp, iGID, rWeight)
        case (OP_HEA)
          !****************************************************************************
          ! Here for the HEA command -- head observation
          !     HEA tgt-name (X,Y) tgt-value obs-weight group-id
          !****************************************************************************
          call IO_Assert(io, lOBSFileOpen, "OBS_Read: No output file")
          sName = sIO_GetField(io, "sName")
          call IO_Assert(io, index(trim(sName), sDelim) == 0, "Target name contains the delimiter character")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          rTgtValue = rIO_GetReal(io, "rTgtValue")
          rTgtWeight = rIO_GetReal(io, "rTgtWeight", def=rONE, minimum=rZERO, maximum=rONE)
          iGrpID = iIO_GetInteger(io, "iGrpID", def=0)
          grp => OBS_GetGroup(io, ObsGroups, iNGrp, iGrpID)
          call IO_Assert(io, associated(grp), "Group ID does not exist")
          rValue = rAEM_Head(io, aem, cZ1)
          write (unit=LU_SCRATCH, fmt=*) trim(sName), sDelim, &
                 trim(grp%sName), sDelim, &
                 rTgtValue, sDelim, &
                 rValue, sDelim, &
                 rTgtValue-rValue, sDelim, &
                 rTgtWeight * grp%rWeight * (rTgtValue-rValue)
        case (OP_DIS)
          !****************************************************************************
          ! Here for the DIS command -- extract the magnitude of the discharge
          !     DIS tgt-name (X,Y) tgt-value obs-weight group-id
          !****************************************************************************
          call IO_Assert(io, lOBSFileOpen, "OBS_Read: No output file")
          sName = sIO_GetField(io, "sName")
          call IO_Assert(io, index(trim(sName), sDelim) == 0, "Target name contains the delimiter character")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          rTgtValue = rIO_GetReal(io, "rTgtValue")
          rTgtWeight = rIO_GetReal(io, "rTgtWeight", def=rONE, minimum=rZERO, maximum=rONE)
          iGrpID = iIO_GetInteger(io, "iGrpID", def=0)
          grp => OBS_GetGroup(io, ObsGroups, iNGrp, iGrpID)
          call IO_Assert(io, associated(grp), "Group ID does not exist")
          rValue = abs(cIO_WorldDischarge(io, cAEM_Discharge(io, aem, cZ1)))
          write (unit=LU_SCRATCH, fmt=*) trim(sName), sDelim, &
                 trim(grp%sName), sDelim, &
                 rTgtValue, sDelim, &
                 rValue, sDelim, &
                 rTgtValue-rValue, sDelim, &
                 rTgtWeight * grp%rWeight * (rTgtValue-rValue)
        case (OP_VEL)
          !****************************************************************************
          ! Here for the VEL command -- extract the magnitude of the velocity
          !     VEL tgt-name (X,Y) tgt-value obs-weight group-id
          !****************************************************************************
          call IO_Assert(io, lOBSFileOpen, "OBS_Read: No output file")
          sName = sIO_GetField(io, "sName")
          call IO_Assert(io, index(trim(sName), sDelim) == 0, "Target name contains the delimiter character")
          cZ1 = cIO_GetCoordinate(io, "cZ1")
          rTgtValue = rIO_GetReal(io, "rTgtValue")
          rTgtWeight = rIO_GetReal(io, "rTgtWeight", def=rONE, minimum=rZERO, maximum=rONE)
          iGrpID = iIO_GetInteger(io, "iGrpID", def=0)
          grp => OBS_GetGroup(io, ObsGroups, iNGrp, iGrpID)
          call IO_Assert(io, associated(grp), "Group ID does not exist")
          rValue = abs(cIO_WorldDischarge(io, cAEM_Velocity(io, aem, cZ1)))
          write (unit=LU_SCRATCH, fmt=*) trim(sName), sDelim, &
                 trim(grp%sName), sDelim, &
                 rTgtValue, sDelim, &
                 rValue, sDelim, &
                 rTgtValue-rValue, sDelim, &
                 rTgtWeight * grp%rWeight * (rTgtValue-rValue)
        case (OP_FLO)
          !****************************************************************************
          ! Here for the FLO command -- extract the flow between points
          !     FLO tgt_name tgt-value obs-weight group-id
          !         (x1,y1)
          !         (x2,y2)
          !         ...
          !         (xN,yN)
          !     END
          !****************************************************************************
          call IO_Assert(io, lOBSFileOpen, "OBS_Read: No output file")
          sName = sIO_GetField(io, "sName")
          call IO_Assert(io, index(trim(sName), sDelim) == 0, "Target name contains the delimiter character")
          rTgtValue = rIO_GetReal(io, "rTgtValue")
          rTgtWeight = rIO_GetReal(io, "rTgtWeight", def=rONE, minimum=rZERO, maximum=rONE)
          iGrpID = iIO_GetInteger(io, "iGrpID", def=0)
          grp => OBS_GetGroup(io, ObsGroups, iNGrp, iGrpID)
          call IO_Assert(io, associated(grp), "Group ID does not exist")
          cZ1 = cHUGE
          do
            call IO_InputRecord(io, (/DIRECTIVE(OP_END, "END")/), iOpCode)
            select case (iOpCode)
              case (kOpData)
                cZ2 = cIO_GetCoordinate(io, "CZ2")
                if (cZ1 == cHUGE) then
                  rValue = rZERO
                  cZ1 = cZ2
                else
                  rValue = rValue + rIO_WorldLength(io, rAEM_Flow(io, aem, (/cZ1, cZ2/)), cZ2-cZ1)
                end if
              case (OP_END)
                exit
            end select
          end do
          write (unit=LU_SCRATCH, fmt=*) trim(sName), sDelim, &
                 trim(grp%sName), sDelim, &
                 rTgtValue, sDelim, &
                 rValue, sDelim, &
                 rTgtValue-rValue, sDelim, &
                 rTgtWeight * grp%rWeight * (rTgtValue-rValue)
        case (OP_GR2)
            call IO_MessageText(io, "GR2 two-point gradient observations are not yet implemented")
        case (OP_GR3)
            call IO_MessageText(io, "GR3 three-point gradient observations are not yet implemented")
        case (OP_GAG)
            call IO_MessageText(io, "GAG stream gage observations are not yet implemented")
      end select
    end do

    call IO_MessageText(io, "Leaving OBS module")

    return
  end subroutine OBS_Read

end module a_obs
