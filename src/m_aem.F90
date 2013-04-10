module m_aem

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

  ! AEMver control module for ModAEM. Allows selection of AEMver types, handles
  ! matrix generation and AEMution.

  use u_constants
  use u_io
  use u_matrix
  use f_well
  use f_dipole
  use f_pond
  use m_aqu
  use m_wl0
  use m_wl1
  use m_pd0
  use m_ls0
  use m_ls1
  use m_ls2
  use m_ls3
  use m_hb0
  use m_as0
#ifndef __GPL__
  use m_cw0
#endif

  implicit none

  public

  type :: AEM_DOMAIN
    !! type AEM_DOMAIN
    !!
    !! Type that holds a complete AEM domain
    !!
    !! Members:
    !!   type(FWL_COLLECTION), pointer :: fwl
    !!     Collection of well functions
    !!   type(FPD_COLLECTION), pointer :: fpd
    !!     Collection of pond functions
    !!   type(FDP_COLLECTION), pointer :: fdp
    !!     Collection of dipole functions
    !!   type(WL0_COLLECTION), pointer :: wl0
    !!     WL0 elements -- discharge wells
    !!   type(WL1_COLLECTION), pointer :: wl1
    !!     WL1 elements -- head-specified wells
    !!   type(PD0_COLLECTION), pointer :: pd0
    !!     PD0 elements -- discharge wells
    !!   type(LS0_COLLECTION), pointer :: ls0
    !!     LS0 elements -- discharge linesinks
    !!   type(LS1_COLLECTION), pointer :: ls1
    !!     LS1 elements -- head linesinks
    !!   type(LS2_COLLECTION), pointer :: ls2
    !!     LS2 elements -- head linesinks with resistance
    !!   type(LS3_COLLECTION), pointer :: ls3
    !!     LS3 elements -- next-generation head linesinks with resistance
    !!   type(HB0_COLLECTION), pointer :: hb0
    !!     HB0 elements -- no-flow barriers
    !!   type(AS0_COLLECTION), pointer :: as0_top
    !!     AS0 elements -- discharge area-sinks(at top of the aquifer)
    !!   type(AS0_COLLECTION), pointer :: as0_bottom
    !!     AS0 elements -- discharge area-sinks(at bottom of the aquifer)
#ifndef __GPL__
    !!   type(CW0_COLLECTION), pointer :: cw0
    !!     CW0 elements -- collector wells
#endif
    !!   type(MAT_MATRIX), pointer :: mat
    !!     Solution matrix object
    !!   integer :: iCurrentIteration
    !!     Current iteration number
    !!   integer :: iAQUStart
    !!     Starting equation number for the AQU module
    !!   integer :: iAQUNUnk
    !!     Number of unknowns for the AQU module
    !!   integer :: iLS1Start
    !!     Starting equation number for the LS1 module
    !!   integer :: iLS1NUnk
    !!     Number of unknowns for the LS1 module
    !!   integer :: iHB0Start
    !!     Starting equation number for the HB0 module
    !!   integer :: iHB0NUnk
    !!     Number of unknowns for the HB0 module
    !!   integer :: iWL1Start
    !!     Starting equation number for the WL1 module
    !!   integer :: iWL1NUnk
    !!     Number of unknowns for the WL1 module
    !!
    integer(kind=AE_INT) :: iCurrentIteration
    integer(kind=AE_INT) :: iAQUStart
    integer(kind=AE_INT) :: iAQUNUnk
    integer(kind=AE_INT) :: iLS1Start
    integer(kind=AE_INT) :: iLS1NUnk
    integer(kind=AE_INT) :: iLS2Start
    integer(kind=AE_INT) :: iLS2NUnk
    integer(kind=AE_INT) :: iLS3Start
    integer(kind=AE_INT) :: iLS3NUnk
    integer(kind=AE_INT) :: iHB0Start
    integer(kind=AE_INT) :: iHB0NUnk
    integer(kind=AE_INT) :: iWL1Start
    integer(kind=AE_INT) :: iWL1NUnk
#ifndef __GPL__
    integer(kind=AE_INT) :: iCW0Start
    integer(kind=AE_INT) :: iCW0NUnk
#endif

    ! True if a solution is present
    logical :: lSolutionPresent

    ! True if the drawdown elements are enabled
    logical :: lDrawdown

    ! Pointers to other collection objects
    type(AQU_COLLECTION), pointer :: aqu
    type(FWL_COLLECTION), pointer :: fwl
    type(FPD_COLLECTION), pointer :: fpd
    type(FDP_COLLECTION), pointer :: fdp
    type(WL0_COLLECTION), pointer :: wl0
    type(PD0_COLLECTION), pointer :: pd0
    type(LS0_COLLECTION), pointer :: ls0
    type(LS1_COLLECTION), pointer :: ls1
    type(LS2_COLLECTION), pointer :: ls2
    type(LS3_COLLECTION), pointer :: ls3
    type(HB0_COLLECTION), pointer :: hb0
    type(WL1_COLLECTION), pointer :: wl1
    type(AS0_COLLECTION), pointer :: as0_top
    type(AS0_COLLECTION), pointer :: as0_bottom
#ifndef __GPL__
    type(CW0_COLLECTION), pointer :: cw0
#endif
    type(MAT_MATRIX), pointer :: mat
  end type AEM_DOMAIN

  real(kind=AE_REAL), private, parameter :: MOVEPOINT = 1.0e-3_AE_REAL

contains


  function AEM_Create(io) result(aem)
    !! function AEM_Create
    !!
    !! Creates a new AEM_DOMAIN object
    !!
    !! Calling Sequence:
    !!    AEM => AEM_Create(io, )
    !!
    !! Arguments:
    !!
    !! Return Value:
    !!   On success, AEM points to a new AEM_DOMAIN object
    !!   On failure(allocation error), fatal error
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(aem, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AEM_Create: allocation failed")

    ! Initialize
    aem%iCurrentIteration = 0
    aem%iAQUStart = 0
    aem%iAQUNUnk = 0
    aem%iLS1Start = 0
    aem%iLS1NUnk = 0
    aem%iLS2NUnk = 0
    aem%iLS2NUnk = 0
    aem%iLS3NUnk = 0
    aem%iLS3NUnk = 0
    aem%iHB0Start = 0
    aem%iHB0NUnk = 0
    aem%iWL1Start = 0
    aem%iWL1Start = 0
    aem%lDrawdown = .false.
    aem%lSolutionPresent = .false.

    ! Create the module objects
    nullify(aem%fwl)
    nullify(aem%fpd)
    nullify(aem%fdp)
    nullify(aem%aqu)
    nullify(aem%mat)
    ! Init pointers to NULL for now
    aem%wl0 => WL0_Create(io)
    aem%pd0 => PD0_Create(io)
    aem%ls0 => LS0_Create(io)
    aem%ls1 => LS1_Create(io)
    aem%ls2 => LS2_Create(io)
    aem%ls3 => LS3_Create(io)
    aem%hb0 => HB0_Create(io)
    aem%wl1 => WL1_Create(io)
    aem%as0_top => AS0_Create(AS0_TOP)
    aem%as0_bottom => AS0_Create(AS0_BOTTOM)
#ifndef __GPL__
    aem%cw0 => CW0_Create(io)
#endif
    ! Create the matrix
    aem%mat => MAT_Create(io)
    return
  end function AEM_Create



  subroutine AEM_Destroy(io, aem)
    !! subroutine AEM_Destroy
    !!
    !! Frees memory allocated for an AEM_DOMAIN object
    !!
    !! Calling Sequence:
    !!     call AEM_Destroy(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    call WL0_Destroy(io, aem%wl0)
    call PD0_Destroy(io, aem%pd0)
    call LS0_Destroy(io, aem%ls0)
    call LS1_Destroy(io, aem%ls1)
    call LS2_Destroy(io, aem%ls2)
    call LS3_Destroy(io, aem%ls3)
    call HB0_Destroy(io, aem%hb0)
    call WL1_Destroy(io, aem%wl1)
    call AS0_Destroy(io, aem%as0_top)
    call AS0_Destroy(io, aem%as0_bottom)
    call MAT_Destroy(io, aem%mat)
#ifndef __GPL__
    call CW0_Destroy(io, aem%cw0)
#endif
    deallocate(aem, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AEM_Destroy: deallocation failed")

    return
  end subroutine AEM_Destroy


  function cAEM_Potential(io, aem, cZ, lNoCheck) result(cOmega)
    !! function cAEM_Potential
    !!
    !! Returns the complex potential at cZ
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Potential(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg, cZFix
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    real(kind=AE_REAL) :: rStrength
    logical :: lCheck

    if (present(lNoCheck)) then
      lCheck = .not. lNoCheck
    else
      lCheck = .true.
    end if

    cZArg = cZ
    if (lCheck) then
      do
        ! Adjust the location if necessary
        if (lFDP_CheckPoint(io, aem%fdp, cZArg, rVERTEXTOL, cZFix, rStrength, &
            iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFWL_CheckPoint(io, aem%fwl, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFPD_CheckPoint(io, aem%fpd, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lAQU_CheckPoint(io, aem%aqu, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cOmega = cAQU_Potential(io, aem%aqu, cZArg) + &
             cFWL_Potential(io, aem%fwl, cZArg) + &
             cFPD_Potential(io, aem%fpd, cZArg) + &
             cFDP_Potential(io, aem%fdp, cZArg) + &
             rAS0_InsidePotential(io, aem%as0_top, cZArg) + &
             rAS0_InsidePotential(io, aem%as0_bottom, cZArg)

    return
  end function cAEM_Potential


  function cAEM_Discharge(io, aem, cZ, lNoCheck) result(cQ)
    !! function cAEM_Discharge
    !!
    !! Returns the discharge vector at cZ
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg, cZFix
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    logical :: lCheck

    if (present(lNoCheck)) then
      lCheck = .not. lNoCheck
    else
      lCheck = .true.
    end if

    cZArg = cZ
    if (lCheck) then
      do
        ! Adjust the location if necessary
        if (lFDP_CheckPoint(io, aem%fdp, cZArg, rVERTEXTOL, cZFix, rStrength, &
            iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFWL_CheckPoint(io, aem%fwl, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFPD_CheckPoint(io, aem%fpd, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lAQU_CheckPoint(io, aem%aqu, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cQ = cAQU_Discharge(io, aem%aqu, cZArg) + &
         cFWL_Discharge(io, aem%fwl, cZArg) + &
         cFPD_Discharge(io, aem%fpd, cZArg) + &
         cFDP_Discharge(io, aem%fdp, cZArg) + &
         cAS0_InsideDischarge(io, aem%as0_top, cZArg) + &
         cAS0_InsideDischarge(io, aem%as0_bottom, cZArg) + &
         cZERO

    return
  end function cAEM_Discharge


  function rAEM_Recharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! function rAEM_Recharge
    !!
    !! Returns the(net) recharge rate at cZ
    !!
    !! Calling Sequence:
    !!    rGamma = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rGamma
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ
    ! NOTE: No singularities in recharge(in ModAEM)

    rGamma = rAQU_Recharge(io, aem%aqu, cZArg) + &
             rFWL_Recharge(io, aem%fwl, cZArg) + &
             rFPD_Recharge(io, aem%fpd, cZArg) + &
             rFDP_Recharge(io, aem%fdp, cZArg) + &
             rAS0_InsideRecharge(io, aem%as0_top, cZArg) + &
             rAS0_InsideRecharge(io, aem%as0_bottom, cZArg)

    return
  end function rAEM_Recharge


  function rAEM_TopRecharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! function rAEM_Recharge
    !!
    !! Returns the(top) recharge rate at cZ
    !!
    !! Calling Sequence:
    !!    rGamma = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rGamma
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ
    ! NOTE: No singularities in recharge(in ModAEM)

    rGamma = rAQU_Recharge(io, aem%aqu, cZArg) + &
             rFWL_Recharge(io, aem%fwl, cZArg) + &
             rFPD_Recharge(io, aem%fpd, cZArg) + &
             rFDP_Recharge(io, aem%fdp, cZArg) + &
             rAS0_InsideRecharge(io, aem%as0_top, cZArg)

    return
  end function rAEM_TopRecharge


  function rAEM_BottomRecharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! function rAEM_Recharge
    !!
    !! Returns the(bottom) recharge rate at cZ
    !!
    !! Calling Sequence:
    !!    rGamma = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rGamma
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ
    ! NOTE: No singularities in recharge(in ModAEM)

    rGamma = rAS0_InsideRecharge(io, aem%as0_bottom, cZArg)

    return
  end function rAEM_BottomRecharge


  function rAEM_Extraction(io, aem, lNoCheck) result(rQ)
    !! function rAEM_Extraction
    !!
    !! Returns the net extraction rate of the model(should be zero!)
    !!
    !! Calling Sequence:
    !!    rQ = cAEM_Extraction(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]

    rQ = rAQU_Extraction(io, aem%aqu) + &
         rFWL_Extraction(io, aem%fwl) + &
         rFPD_Extraction(io, aem%fpd) + &
         rFDP_Extraction(io, aem%fdp) + &
         rAS0_Extraction(io, aem%as0_top) - &
         rAS0_Extraction(io, aem%as0_bottom)

    return
  end function rAEM_Extraction


  function rAEM_HeadAtWell(io, aem, cZ, iFWLIndex, lNoCheck) result(rH)
    !! function cAEM_DischargeAtWell
    !!
    !! Returns the discharge vector at cZ, excluding the well with index iFWLIndex
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iFWLIndex
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg, cZFix
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    logical :: lCheck

    if (present(lNoCheck)) then
      lCheck = .not. lNoCheck
    else
      lCheck = .true.
    end if

    cZArg = cZ
    if (lCheck) then
      do
        ! Adjust the location if necessary
        if (lFDP_CheckPoint(io, aem%fdp, cZArg, rVERTEXTOL, cZFix, rStrength, &
            iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFWL_CheckPoint(io, aem%fwl, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFPD_CheckPoint(io, aem%fpd, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lAQU_CheckPoint(io, aem%aqu, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    rH = rAQU_PotentialToHead(io, aem%aqu, &
                              real(cAEM_Potential(io, aem, cZArg) - &
                                   cFWL_Potential(io, aem%fwl, cZArg, iFWLIndex, 1)), &
                              cZArg)

    return
  end function rAEM_HeadAtWell


  function cAEM_DischargeAtWell(io, aem, cZ, iFWLIndex, lNoCheck) result(cQ)
    !! function cAEM_DischargeAtWell
    !!
    !! Returns the discharge vector at cZ, excluding the well with index iFWLIndex
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iFWLIndex
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZArg, cZFix
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    logical :: lCheck

    if (present(lNoCheck)) then
      lCheck = .not. lNoCheck
    else
      lCheck = .true.
    end if

    cZArg = cZ
    if (lCheck) then
      do
        ! Adjust the location if necessary
        if (lFDP_CheckPoint(io, aem%fdp, cZArg, rVERTEXTOL, cZFix, rStrength, &
            iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFWL_CheckPoint(io, aem%fwl, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lFPD_CheckPoint(io, aem%fpd, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else if (lAQU_CheckPoint(io, aem%aqu, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cQ = cAEM_Discharge(io, aem, cZArg) - &
         cFWL_Discharge(io, aem%fwl, cZArg, iFWLIndex, 1)

    return
  end function cAEM_DischargeAtWell


  function rAEM_Head(io, aem, cZ, lNoCheck) result(rHead)
    !! function rAEM_Head
    !!
    !! Returns the head at cZ
    !!
    !! Calling Sequence:
    !!    rHead = rAEM_Head(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHead
    ! [ LOCALS ]

    rHead = rAQU_PotentialToHead(io, aem%aqu, real(cAEM_Potential(io, aem, cZ, lNoCheck)), cZ)

    return
  end function rAEM_Head


  function rAEM_SatdThick(io, aem, cZ, lNoCheck) result(rSatdThick)
    !! function rAEM_SatdThick
    !!
    !! Returns the saturated thickness at cZ
    !!
    !! Calling Sequence:
    !!    rHead = rAEM_SatdThick(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rSatdThick
    ! [ LOCALS ]

    rSatdThick = rAQU_SatdThickness(io, aem%aqu, cZ, real(cAEM_Potential(io, aem, cZ, lNoCheck)))

    return
  end function rAEM_SatdThick


  function rAEM_InterfaceElevation(io, aem, cZ) result(rIfcElev)
    !! function rAEM_SatdThick
    !!
    !! Returns the sea-water interface elevation at cZ
    !!
    !! Calling Sequence:
    !!    rHead = rAEM_InterfaceElevation(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rIfcElev
    ! [ LOCALS ]

    rIfcElev = rAQU_InterfaceElevation(io, aem%aqu, cZ, real(cAEM_Potential(io, aem, cZ, .true.)))

    return
  end function rAEM_InterfaceElevation


  function rAEM_Flow(io, aem, cZ, lNoCheck) result(rFlow)
    !! function rAEM_Flow
    !!
    !! Returns the integrated flow across the path cZ
    !!
    !! Calling Sequence:
    !!    rFlow = rAEM_Flow(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), dimension(:), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: cZArg
    complex(kind=AE_REAL) :: cZFix
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iStat, i
    logical :: lCheck

    if (present(lNoCheck)) then
      lCheck = .not. lNoCheck
    else
      lCheck = .true.
    end if

    allocate(cZArg(size(cZ)), stat = iStat)
    call IO_Assert(io, (iStat == 0), "rAEM_Flow: Allocation failed")
    cZArg = cZ

    if (lCheck) then
      do i = 1, size(cZArg)
        do
          ! Adjust the location if necessary
          if (lFDP_CheckPoint(io, aem%fdp, cZArg(i), rVERTEXTOL, cZFix, rStrength, &
              iElementType, iElementString, iElementVertex, iElementFlag)) then
            cZArg(i) = cZFix
            cycle
          else if (lFWL_CheckPoint(io, aem%fwl, cZArg(i), rVERTEXTOL, cZFix, rStrength, &
                 iElementType, iElementString, iElementVertex, iElementFlag)) then
            cZArg(i) = cZFix
            cycle
          else if (lFPD_CheckPoint(io, aem%fpd, cZArg(i), rVERTEXTOL, cZFix, rStrength, &
                 iElementType, iElementString, iElementVertex, iElementFlag)) then
            cZArg(i) = cZFix
            cycle
          else if (lAQU_CheckPoint(io, aem%aqu, cZArg(i), rVERTEXTOL, cZFix, rStrength, &
                 iElementType, iElementString, iElementVertex, iElementFlag)) then
            cZArg(i) = cZFix
            cycle
          else
            exit
          end if
        end do
      end do
    end if

    rFlow = rAQU_Flow(io, aem%aqu, cZArg) + &
            rFWL_Flow(io, aem%fwl, cZArg) + &
            rFPD_Flow(io, aem%fpd, cZArg) + &
            rFDP_Flow(io, aem%fdp, cZArg) + &
            rAS0_InsideFlow(io, aem%as0_top, cZArg) + &
            rAS0_InsideFlow(io, aem%as0_bottom, cZArg)
    deallocate(cZArg)

    return
  end function rAEM_Flow


  function cAEM_Velocity(io, aem, cZ, lNoCheck) result(cV)
    !! function cAEM_DischargeAtWell
    !!
    !! Returns the head at cZ
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem, cZ)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    complex :: cZ
    !!              Complex coordinate in question
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cV
    ! [ LOCALS ]

    cV = cAEM_Discharge(io, aem, cZ, lNoCheck)
    cV = cAQU_DischargeToVelocity(io, aem%aqu, cAEM_Discharge(io, aem, cZ, lNoCheck), &
         cZ, real(cAEM_Potential(io, aem, cZ, lNoCheck)))

    return
  end function cAEM_Velocity


  subroutine AEM_AllocateFunctions(io, aem)
    !! subroutine AEM_Alloc
    !!
    !! Allocates the space in the MAT, FWL, FPD, and FDP modules
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNFWL, iNFPD, iNFDP
    integer(kind=AE_INT) :: i

    !******************************************************************************************
    ! Compute the required space for all function modules
    !******************************************************************************************
    iNFWL = 0
    iNFPD = 0
    iNFDP = 0

    !******************************************************************************************
    ! Given-strength(no equations necessary) go here
    !******************************************************************************************
    ! WL0 Module(given wells)
    iNFWL = iNFWL + iWL0_GetInfo(io, aem%wl0, SIZE_FWL, 0)
    iNFPD = iNFPD + iWL0_GetInfo(io, aem%wl0, SIZE_FPD, 0)
    iNFDP = iNFDP + iWL0_GetInfo(io, aem%wl0, SIZE_FDP, 0)
    ! PD0 Module(given ponds)
    iNFWL = iNFWL + iPD0_GetInfo(io, aem%pd0, SIZE_FWL, 0)
    iNFPD = iNFPD + iPD0_GetInfo(io, aem%pd0, SIZE_FPD, 0)
    iNFDP = iNFDP + iPD0_GetInfo(io, aem%pd0, SIZE_FDP, 0)
    ! LS0 Module(given linesinks)
    iNFWL = iNFWL + iLS0_GetInfo(io, aem%ls0, SIZE_FWL, 0)
    iNFPD = iNFPD + iLS0_GetInfo(io, aem%ls0, SIZE_FPD, 0)
    iNFDP = iNFDP + iLS0_GetInfo(io, aem%ls0, SIZE_FDP, 0)
    ! AS0 Module(given area-sinks) TOP
    iNFWL = iNFWL + iAS0_GetInfo(io, aem%as0_top, SIZE_FWL, 0)
    iNFPD = iNFPD + iAS0_GetInfo(io, aem%as0_top, SIZE_FPD, 0)
    iNFDP = iNFDP + iAS0_GetInfo(io, aem%as0_top, SIZE_FDP, 0)
    ! AS0 Module(given area-sinks) BOTTOM
    iNFWL = iNFWL + iAS0_GetInfo(io, aem%as0_bottom, SIZE_FWL, 0)
    iNFPD = iNFPD + iAS0_GetInfo(io, aem%as0_bottom, SIZE_FPD, 0)
    iNFDP = iNFDP + iAS0_GetInfo(io, aem%as0_bottom, SIZE_FDP, 0)
    !******************************************************************************************
    ! Modify for new element types:
    ! Add additional given-strength modules above...
    !******************************************************************************************

    !******************************************************************************************
    ! Unknown-strength(main matrix equations necessary) go here
    !******************************************************************************************
    ! AQU Module(reference point, etc.)
    iNFWL = iNFWL + iAQU_GetInfo(io, aem%aqu, SIZE_FWL, 0)
    iNFPD = iNFPD + iAQU_GetInfo(io, aem%aqu, SIZE_FPD, 0)
    iNFDP = iNFDP + iAQU_GetInfo(io, aem%aqu, SIZE_FDP, 0)
    ! LS1 Module(head-specified linesinks)
    iNFWL = iNFWL + iLS1_GetInfo(io, aem%ls1, SIZE_FWL, 0)
    iNFPD = iNFPD + iLS1_GetInfo(io, aem%ls1, SIZE_FPD, 0)
    iNFDP = iNFDP + iLS1_GetInfo(io, aem%ls1, SIZE_FDP, 0)
    ! LS2 Module(head-specified linesinks with resistance/routing)
    iNFWL = iNFWL + iLS2_GetInfo(io, aem%ls2, SIZE_FWL, 0)
    iNFPD = iNFPD + iLS2_GetInfo(io, aem%ls2, SIZE_FPD, 0)
    iNFDP = iNFDP + iLS2_GetInfo(io, aem%ls2, SIZE_FDP, 0)
    ! LS3 Module(head-specified linesinks with resistance/routing)
    iNFWL = iNFWL + iLS3_GetInfo(io, aem%ls3, SIZE_FWL, 0)
    iNFPD = iNFPD + iLS3_GetInfo(io, aem%ls3, SIZE_FPD, 0)
    iNFDP = iNFDP + iLS3_GetInfo(io, aem%ls3, SIZE_FDP, 0)
    ! HB0 Module(no-flow boundaries)
    iNFWL = iNFWL + iHB0_GetInfo(io, aem%hb0, SIZE_FWL, 0)
    iNFPD = iNFPD + iHB0_GetInfo(io, aem%hb0, SIZE_FPD, 0)
    iNFDP = iNFDP + iHB0_GetInfo(io, aem%hb0, SIZE_FDP, 0)
    ! WL1 Module(head-specified wells)
    iNFWL = iNFWL + iWL1_GetInfo(io, aem%wl1, SIZE_FWL, 0)
    iNFPD = iNFPD + iWL1_GetInfo(io, aem%wl1, SIZE_FPD, 0)
    iNFDP = iNFDP + iWL1_GetInfo(io, aem%wl1, SIZE_FDP, 0)
#ifndef __GPL__
    ! CW0 Module(collector wells)
    iNFWL = iNFWL + iCW0_GetInfo(io, aem%cw0, SIZE_FWL, 0)
    iNFPD = iNFPD + iCW0_GetInfo(io, aem%cw0, SIZE_FPD, 0)
    iNFDP = iNFDP + iCW0_GetInfo(io, aem%cw0, SIZE_FDP, 0)
#endif

    !******************************************************************************************
    ! Modify for new element types:
    ! Add additional unknown-strength modules above...
    !******************************************************************************************

    !******************************************************************************************
    ! Now, allocate the space for all the functions and the matrix generator
    !******************************************************************************************
    ! Well functions...
    aem%fwl => FWL_Create(io, iNFWL)
    !******************************************************************************************
    ! Pond functions...
    aem%fpd => FPD_Create(io, iNFPD)
    !******************************************************************************************
    ! dipole functions...
    aem%fdp => FDP_Create(io, iNFDP)
    !******************************************************************************************
    ! Modify for new function types:
    ! Add additional function modules above...
    !******************************************************************************************

    !******************************************************************************************
    ! If we made it here, allocation is complete. Now use the element module setup
    ! routines to set up the functions and matrix generator
    !******************************************************************************************
    call AQU_SetupFunctions(io, aem%aqu, aem%fdp)
    call WL0_SetupFunctions(io, aem%wl0, aem%fwl, aem%aqu)
    call PD0_SetupFunctions(io, aem%pd0, aem%fpd)
    call LS0_SetupFunctions(io, aem%ls0, aem%fwl, aem%fdp)
    call LS1_SetupFunctions(io, aem%ls1, aem%fwl, aem%fdp)
    call LS2_SetupFunctions(io, aem%ls2, aem%fwl, aem%fdp)
    call LS3_SetupFunctions(io, aem%ls3, aem%fwl, aem%fdp)
    call HB0_SetupFunctions(io, aem%hb0, aem%fdp)
    call WL1_SetupFunctions(io, aem%wl1, aem%fwl, aem%aqu)
    call AS0_SetupFunctions(io, aem%as0_top, aem%fwl, aem%fdp)
    call AS0_SetupFunctions(io, aem%as0_bottom, aem%fwl, aem%fdp)
#ifndef __GPL__
    call CW0_SetupFunctions(io, aem%cw0, aem%fwl, aem%fdp)
#endif

    return
  end subroutine AEM_AllocateFunctions


  subroutine AEM_AllocateMatrix(io, aem, iIteration)
    !! subroutine AEM_Alloc
    !!
    !! Allocates the space in the MAT module
    !!
    !! Calling Sequence:
    !!    cOmega = cAEM_Discharge(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNEQ, iNUN
    integer(kind=AE_INT) :: i

    !******************************************************************************************
    ! Compute the required space for all function modules
    !******************************************************************************************
    iNEQ = 0
    iNUN = 0

    !******************************************************************************************
    ! Unknown-strength(main matrix equations necessary) go here
    !******************************************************************************************
    ! AQU Module(reference point, etc.)
    iNEQ = iNEQ + iAQU_GetInfo(io, aem%aqu, SIZE_EQUATIONS, iIteration)
    aem%iAQUStart = 1
    aem%iAQUNUnk = iAQU_GetInfo(io, aem%aqu, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iAQUNUnk
    ! LS1 Module(head-specified linesinks)
    iNEQ = iNEQ + iLS1_GetInfo(io, aem%ls1, SIZE_EQUATIONS, iIteration)
    aem%iLS1Start = aem%iAQUNUnk + 1
    aem%iLS1NUnk = iLS1_GetInfo(io, aem%ls1, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iLS1NUnk
    ! LS2 Module(head-specified linesinks with resistance/routing)
    iNEQ = iNEQ + iLS2_GetInfo(io, aem%ls2, SIZE_EQUATIONS, iIteration)
    aem%iLS2Start = aem%iAQUNUnk + aem%iLS1NUnk + 1
    aem%iLS2NUnk = iLS2_GetInfo(io, aem%ls2, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iLS2NUnk
    ! LS3 Module(head-specified linesinks with resistance/routing)
    iNEQ = iNEQ + iLS3_GetInfo(io, aem%ls3, SIZE_EQUATIONS, iIteration)
    aem%iLS3Start = aem%iAQUNUnk + aem%iLS1NUnk + aem%iLS2NUnk + 1
    aem%iLS3NUnk = iLS3_GetInfo(io, aem%ls3, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iLS3NUnk
    ! HB0 Module(no-flow boundaries)
    iNEQ = iNEQ + iHB0_GetInfo(io, aem%hb0, SIZE_EQUATIONS, iIteration)
    aem%iHB0Start = aem%iAQUNUnk + aem%iLS1NUnk + aem%iLS2NUnk + aem%iLS3NUnk + 1
    aem%iHB0NUnk = iHB0_GetInfo(io, aem%hb0, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iHB0NUnk
    ! WL1 Module(head-specified wells)
    iNEQ = iNEQ + iWL1_GetInfo(io, aem%wl1, SIZE_EQUATIONS, iIteration)
    aem%iWL1Start = aem%iAQUNUnk + aem%iLS1NUnk + aem%iLS2NUnk + aem%iLS3NUnk + aem%iHB0NUnk + 1
    aem%iWL1NUnk = iWL1_GetInfo(io, aem%wl1, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iWL1NUnk
#ifndef __GPL__
    ! CW0 Module(head-specified wells)
    iNEQ = iNEQ + iCW0_GetInfo(io, aem%cw0, SIZE_EQUATIONS, iIteration)
    aem%iCW0Start = aem%iAQUNUnk + aem%iLS1NUnk + aem%iLS2NUnk + aem%iLS3NUnk + aem%iHB0NUnk + aem%iWL1NUnk + 1
    aem%iCW0NUnk = iCW0_GetInfo(io, aem%cw0, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + aem%iCW0NUnk
#endif

    !******************************************************************************************
    ! Modify for new element types:
    ! Add additional unknown-strength modules above...
    !******************************************************************************************

    !******************************************************************************************
    ! Now, allocate the matrix and control-point arrays.
    call MAT_Alloc(io, aem%mat, iNEQ, iNUN)

    !******************************************************************************************
    ! If we made it here, allocation is complete. Now use the element module setup
    ! routines to set up the functions and matrix generator
    !******************************************************************************************
    if (aem%iAQUNUnk /= 0) call AQU_SetupMatrix(io, aem%aqu, aem%mat)
    if (aem%iLS1NUnk /= 0) call LS1_SetupMatrix(io, aem%ls1, aem%aqu, aem%mat)
    if (aem%iLS2NUnk /= 0) call LS2_SetupMatrix(io, aem%ls2, aem%aqu, aem%mat)
    if (aem%iLS3NUnk /= 0) call LS3_SetupMatrix(io, aem%ls3, aem%aqu, aem%mat)
    if (aem%iHB0NUnk /= 0) call HB0_SetupMatrix(io, aem%hb0, aem%mat)
    if (aem%iWL1NUnk /= 0) call WL1_SetupMatrix(io, aem%wl1, aem%aqu, aem%mat)
#ifndef __GPL__
    if (aem%iCW0NUnk /= 0) call CW0_SetupMatrix(io, aem%cw0, aem%aqu, aem%mat)
#endif

    return
  end subroutine AEM_AllocateMatrix


  subroutine AEM_PreSolve(io, aem)
    !! Calls all the pre-solve steps
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! Run the pre-solve stuff
    call AQU_PreSolve(io, aem%aqu)
    call HB0_PreSolve(io, aem%hb0)
    call WL0_PreSolve(io, aem%wl0)
    call WL1_PreSolve(io, aem%wl1)
    call LS0_PreSolve(io, aem%ls0)
    call LS1_PreSolve(io, aem%ls1)
    call LS2_PreSolve(io, aem%ls2)
    call LS3_PreSolve(io, aem%ls3)
    call PD0_PreSolve(io, aem%pd0)
    call AS0_PreSolve(io, aem%as0_top, aem%aqu)
    call AS0_PreSolve(io, aem%as0_bottom, aem%aqu)
#ifndef __GPL__
    call CW0_PreSolve(io, aem%cw0)
#endif

    return
  end subroutine AEM_PreSolve


  subroutine AEM_Solve(io, aem, iNIter, iNPolishIter, rRelaxation)
    !! subroutine AEM_Solve
    !!
    !! Performs the solution process
    !!
    !! Calling Sequence:
    !!    call AEM_Solve(io, aem, iNIter)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!    (in)    integer :: iNIter
    !!              Number of iterations to be performed
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(inout) :: iNIter
    integer(kind=AE_INT), intent(in) :: iNPolishIter
    real(kind=AE_REAL), intent(in) :: rRelaxation
    type(IO_STATUS), pointer :: io
    ! Locals
    integer(kind=AE_INT) :: iIter, i, iChanges, iUpdate

    ! Heeeere we go!
    if (aem%lSolutionPresent) then
      call IO_MessageText(io, "Solution is present -- continuing solve")
    else
      call IO_MessageText(io, "Allocating space for functions and matrix")
      call AEM_PreSolve(io, aem)
      call AEM_AllocateFunctions(io, aem)
      aem%lSolutionPresent = .true.
    end if

    ! Set up the WL1 Bessel functions if necessary
    iChanges = iWL1_Prepare(io, aem%wl1, aem%aqu, aem%iCurrentIteration)

    ! Iteration loop...
    do iIter = 1, iNIter
      call IO_MessageText(io, '\n')
      aem%iCurrentIteration = aem%iCurrentIteration+1
      write (unit=IO_MessageBuffer, fmt="(""Iteration: "", i5)") aem%iCurrentIteration
      call IO_MessageText(io)

      ! Call the 'Modify' routines to polish up the boundary conditions...
      if (iIter>1) then
        iUpdate = iAQU_Prepare(io, aem%aqu, aem%iCurrentIteration)

        ! When drawdown elements are disabled, there are no well discharge adjustments
        iUpdate = iUpdate + WL0_AdjustDischarges(io, aem%wl0, aem%aqu)

        ! Check the free surface
        iUpdate = iUpdate + AEM_UpdateFreeSurface(io, aem)

        ! What other elements need adjustments?
        if (iIter>1 .and. iIter < iNIter-iNPolishIter+1 ) then
          iUpdate = iUpdate + &
                    iWL1_Prepare(io, aem%wl1, aem%aqu, aem%iCurrentIteration) + &
                    iLS1_Prepare(io, aem%ls1, aem%aqu, aem%iCurrentIteration) + &
                    iLS2_DoRouting(io, aem%ls2, aem%aqu, aem%iCurrentIteration, .true., .false.) + &
                    iLS2_Prepare(io, aem%ls2, aem%aqu, aem%iCurrentIteration) + &
                    iLS3_DoRouting(io, aem%ls3, aem%aqu, aem%iCurrentIteration, .true., .false.) + &
                    iLS3_Prepare(io, aem%ls3, aem%aqu, aem%iCurrentIteration) + &
                    iHB0_Prepare(io, aem%hb0, aem%iCurrentIteration)
#ifndef __GPL__
          iUpdate = iUpdate + &
                    iCW0_Prepare(io, aem%cw0, aem%aqu, aem%iCurrentIteration)
        end if
#endif
        if (iUpdate > 0) then
          call AEM_Update(io, aem)
          call AEM_ComputeCheck(io, aem, .false.)
        end if
      end if

      ! Shall I regenerate the matrix coefficients?
      !    print *,'AQU',iAQU_GetInfo(io, aem%aqu, INFO_REGENERATE, aem%iCurrentIteration)
      !    print *,'WL1',iWL1_GetInfo(io, aem%wl1, INFO_REGENERATE, aem%iCurrentIteration)
      !    print *,'LS1',iLS1_GetInfo(io, aem%ls1, INFO_REGENERATE, aem%iCurrentIteration)
      !    print *,'LS2',iLS2_GetInfo(io, aem%ls2, INFO_REGENERATE, aem%iCurrentIteration)
      !    print *,'LS3',iLS3_GetInfo(io, aem%ls3, INFO_REGENERATE, aem%iCurrentIteration)
      !    print *,'HB0',iHB0_GetInfo(io, aem%hb0, INFO_REGENERATE, aem%iCurrentIteration)
#ifndef __GPL__
      !    print *,'CW0',iCW0_GetInfo(io, aem%cw0, INFO_REGENERATE, aem%iCurrentIteration)
#endif

      if (iAQU_GetInfo(io, aem%aqu, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
          iWL1_GetInfo(io, aem%wl1, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
          iLS1_GetInfo(io, aem%ls1, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
          iLS2_GetInfo(io, aem%ls2, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
          iLS3_GetInfo(io, aem%ls3, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
          iHB0_GetInfo(io, aem%hb0, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
#ifndef __GPL__
          iCW0_GetInfo(io, aem%cw0, INFO_REGENERATE, aem%iCurrentIteration) /= 0 .or. &
#endif
          .false. &
          ) then
        call IO_MessageText(io, "Allocating matrix")
        call MAT_Clear(io, aem%mat)
        call AEM_AllocateMatrix(io, aem, aem%iCurrentIteration)
        call IO_MessageText(io, "Generating matrix")
        call AEM_GenerateMatrix(io, aem, aem%iCurrentIteration)
        call IO_MessageText(io, "Decomposing matrix")
        call MAT_Decompose(io, aem%mat)
        print *, 'Condition Number: ', aem%mat%rCond
      end if

      ! Initialize the CHECK values on the first go-round
      if (aem%iCurrentIteration == 1) then
        call AEM_ComputeCheck(io, aem, .true.)
      end if

      if (io%lDebug) then
        write (IO_MessageBuffer, '(i3)') aem%iCurrentIteration
      end if

      call IO_MessageText(io, "Generating solution...")
      ! Compute the right-hand side
      call AEM_GenerateRHS(io, aem, aem%iCurrentIteration, .false.)
      ! Solve for the new coefficients
      if (aem%iCurrentIteration == 1) then
        call AEM_SolveMatrix(io, aem, rONE, .false.)
      else
        call AEM_SolveMatrix(io, aem, rRelaxation, .false.)
      end if
      ! Update the elements
      call AEM_Update(io, aem)

      ! Compute the new check values
      call AEM_ComputeCheck(io, aem, .true.)

      ! Clear the AQU preconditioning
      call AQU_SetPrecondition(io, aem%aqu, .false.)

    end do

    ! Do routing along LS2 strings for final summary
    print *, 'Performing stream routing calculations'
    iChanges = iLS2_DoRouting(io, aem%ls2, aem%aqu, aem%iCurrentIteration, .false., .true.)

    ! Do routing along LS3 strings for final summary
    iChanges = iLS3_DoRouting(io, aem%ls3, aem%aqu, aem%iCurrentIteration, .false., .true.)

    ! Check the discharge-specified partially-penetrating wells
    call WL0_SolvePartialPenetration(io, aem%wl0, aem%aqu)

    call IO_MessageText(io)
    call IO_MessageText(io, "Solution complete")

    return
  end subroutine AEM_Solve


  function rAEM_GetCoefficientMultiplier(io, aem, iElementType, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Looks up the proper multiplier
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    select case (iElementType)
      case (ELEM_AQU)
        rMultiplier = rAQU_GetCoefficientMultiplier(io, aem%aqu, iElementString, iElementVertex, iElementFlag)
      case (ELEM_IN0)
        rMultiplier = rIN0_GetCoefficientMultiplier(io, aem%aqu%in0, iElementString, iElementVertex, iElementFlag)
      case (ELEM_LS1)
        rMultiplier = rLS1_GetCoefficientMultiplier(io, aem%ls1, iElementString, iElementVertex, iElementFlag)
      case (ELEM_LS2)
        rMultiplier = rLS2_GetCoefficientMultiplier(io, aem%ls2, iElementString, iElementVertex, iElementFlag)
      case (ELEM_LS3)
        rMultiplier = rLS3_GetCoefficientMultiplier(io, aem%ls3, iElementString, iElementVertex, iElementFlag)
      case (ELEM_HB0)
        rMultiplier = rHB0_GetCoefficientMultiplier(io, aem%hb0, iElementString, iElementVertex, iElementFlag)
      case (ELEM_WL1)
        rMultiplier = rWL1_GetCoefficientMultiplier(io, aem%wl1, iElementString, iElementVertex, iElementFlag)
#ifndef __GPL__
      case (ELEM_CW0)
        rMultiplier = rCW0_GetCoefficientMultiplier(io, aem%cw0, iElementString, iElementVertex, iElementFlag)
#endif
    end select

    return
  end function rAEM_GetCoefficientMultiplier


  subroutine AEM_GenerateMatrix(io, aem, iIteration)
    !! subroutine AEM_GenerateMatrix
    !!
    !! Generates the matrix coefficients
    !!
    !! Calling Sequence:
    !!    call AEM_GenerateMatrix(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iRow
    real(kind=AE_REAL) :: rMultiplier, rGhbDistance
    real(kind=AE_REAL), dimension(:), allocatable :: rARow
    integer(kind=AE_INT) :: iStat
    complex(kind=AE_REAL), dimension(10) :: cCPZ
    integer(kind=AE_INT) :: iNCP, iEqType, ic
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    complex(kind=AE_REAL) :: cOrientation

    ! First allocate the storage for the matrix row being generated
    allocate(rARow(aem%mat%iNVar), stat = iStat)
    call IO_Assert(io, (iStat == 0), "AEM_GenerateMatrix: Allocation failed")

    ! This is the main matrix generator loop. It generates the matrix row-by-row, and
    ! is parallelizable.
    print *, 'Solving a system of ', aem%mat%iNEqn, ' equations'
    do iRow = 1, aem%mat%iNEqn
      ! Get the matrix generator information for the row
      call MAT_GetEquation(io, aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, iElementString, &
           iElementVertex, iElementFlag, cOrientation, rGhbDistance)
      ! rMultiplier is a factor that an element module provides that is used to scale all
      ! the coefficients in a row, e.g. for Igor Jankovic's inhomogeneity formulation
      rMultiplier = rAEM_GetCoefficientMultiplier(io, aem, iElementType, iElementString, &
                    iElementVertex, iElementFlag)


      !*****************************************************************************************************
      ! Correct the matrix coefficients and compute the right-hand side constant for the equation; each
      ! element module with unknown strength coefficients should have a call here.
      !*****************************************************************************************************

      ! The array constructor below is used to ensure that the proper array size is received
      if (aem%iAQUNUnk > 0) then
        call AQU_ComputeCoefficients(io, aem%aqu, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iAQUStart:aem%iAQUStart+aem%iAQUNUnk-1))
      end if

      ! Generate the coefficients for the LS1 element module and store them in the matrix row
      if (aem%iLS1NUnk > 0) then
        call LS1_ComputeCoefficients(io, aem%ls1, aem%fwl, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iLS1Start:aem%iLS1Start+aem%iLS1NUnk-1))
      end if

      ! Generate the coefficients for the LS2 element module and store them in the matrix row
      if (aem%iLS2NUnk > 0) then
        call LS2_ComputeCoefficients(io, aem%ls2, aem%aqu, aem%fwl, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iLS2Start:aem%iLS2Start+aem%iLS2NUnk-1))
      end if

      ! Generate the coefficients for the LS3 element module and store them in the matrix row
      if (aem%iLS3NUnk > 0) then
        call LS3_ComputeCoefficients(io, aem%ls3, aem%aqu, aem%fwl, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iLS3Start:aem%iLS3Start+aem%iLS3NUnk-1))
      end if

      ! Generate the coefficients for the HB0 element module and store them in the matrix row
      if (aem%iHB0NUnk > 0) then
        call HB0_ComputeCoefficients(io, aem%hb0, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iHB0Start:aem%iHB0Start+aem%iHB0NUnk-1))
      end if

      ! Generate the coefficients for the LS1 element module and store them in the matrix row
      if (aem%iWL1NUnk > 0) then
        call WL1_ComputeCoefficients(io, aem%wl1, aem%fwl, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iWL1Start:aem%iWL1Start+aem%iWL1NUnk-1))
      end if

#ifndef __GPL__
      ! Generate the coefficients for the LS1 element module and store them in the matrix row
      if (aem%iCW0NUnk > 0) then
        call CW0_ComputeCoefficients(io, aem%cw0, aem%aqu, aem%fwl, aem%fdp, (/(cCPZ(ic), ic = 1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(aem%iCW0Start:aem%iCW0Start+aem%iCW0NUnk-1))
      end if
#endif

      !*****************************************************************************************************
      ! Add more element modules here
      !*****************************************************************************************************

      ! Store the row into the matrix
      call MAT_SetRow(io, aem%mat, iRow, rARow(1:aem%mat%iNVar), rZERO)

    end do

    deallocate(rARow)

    return
  end subroutine AEM_GenerateMatrix


  subroutine AEM_GenerateRHS(io, aem, iIteration, lDirect)
    !! subroutine AEM_GenerateRHS
    !!
    !! Generates the right-hand side vector
    !!
    !! Calling Sequence:
    !!    call AEM_GenerateMatrix(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iRow
    real(kind=AE_REAL) :: rRHS
    integer(kind=AE_INT) :: iStat
    complex(kind=AE_REAL), dimension(10) :: cCPZ
    integer(kind=AE_INT) :: iNCP, iEqType, iElementType, ic
    integer(kind=AE_INT) :: iElementString, iElementVertex, iElementFlag
    complex(kind=AE_REAL) :: cOrientation
    real(kind=AE_REAL) :: rGhbDistance

    ! This is the main matrix generator loop. It generates the matrix row-by-row, and
    ! is parallelizable.
    do iRow = 1, aem%mat%iNEqn
      ! Get the matrix generator information for the row
      call MAT_GetEquation(io, aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, iElementString, &
           iElementVertex, iElementFlag, cOrientation, rGhbDistance)

      ! Call the proper module to compute the right-hand side
      select case (iElementType)
        case (ELEM_AQU)
          rRHS = rAQU_ComputeRHS(io, aem%aqu, aem%fdp, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_IN0)
          rRHS = rAQU_ComputeRHS(io, aem%aqu, aem%fdp, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_LS1)
          rRHS = rLS1_ComputeRHS(io, aem%ls1, aem%aqu, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_LS2)
          rRHS = rLS2_ComputeRHS(io, aem%ls2, aem%aqu, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_LS3)
          rRHS = rLS3_ComputeRHS(io, aem%ls3, aem%aqu, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_HB0)
          rRHS = rHB0_ComputeRHS(io, aem%hb0, iEqType, iElementType, iElementString, iElementVertex, &
                 iElementFlag, iIteration, lDirect)
        case (ELEM_WL1)
          rRHS = rWL1_ComputeRHS(io, aem%wl1, aem%aqu, iEqType, iElementType, iElementString, &
                 iElementVertex, iElementFlag, iIteration, lDirect)
#ifndef __GPL__
        case (ELEM_CW0)
          rRHS = rCW0_ComputeRHS(io, aem%cw0, aem%aqu, iEqType, iElementType, iElementString, &
                 iElementVertex, iElementFlag, iIteration, lDirect)
#endif
      end select

      call MAT_SetRHS(io, aem%mat, iRow, rRHS)
    end do

    return
  end subroutine AEM_GenerateRHS


  subroutine AEM_SolveMatrix(io, aem, rRelaxation, lDirect)
    !! subroutine AEM_SolveMatrix
    !!
    !! Generates the right-hand side vector
    !!
    !! Calling Sequence:
    !!    call AEM_SolveMatrix(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    real(kind=AE_REAL), intent(in) :: rRelaxation
    logical :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    real(kind=AE_REAL) :: rValue
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag

    call MAT_Solve(io, aem%mat)

    ! Now, finish the job by setting the final results for all elements...
    do i = 1, aem%mat%iNEqn
      ! Get the result from the matrix AEMver
      call MAT_GetVariable(io, aem%mat, i, rValue, iElementType, &
           iElementString, iElementVertex, iElementFlag)

      ! Store it in the proper element data structures
      select case (iElementType)
        case (ELEM_AQU)
          call AQU_StoreResult(io, aem%aqu, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_IN0)
          call AQU_StoreResult(io, aem%aqu, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS1)
          call LS1_StoreResult(io, aem%ls1, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS2)
          call LS2_StoreResult(io, aem%ls2, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS3)
          call LS3_StoreResult(io, aem%ls3, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_HB0)
          call HB0_StoreResult(io, aem%hb0, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_WL1)
          call WL1_StoreResult(io, aem%wl1, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
#ifndef __GPL__
        case (ELEM_CW0)
          call CW0_StoreResult(io, aem%cw0, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
#endif
      end select
    end do

    return
  end subroutine AEM_SolveMatrix


  subroutine AEM_Update(io, aem)
    !! subroutine AEM_Update
    !!
    !! Update all element modules.  The update procedures should:
    !!   1)  Update the strength parameters for all functional modules
    !!   2)  Update any element module-dependent features(e.g. turning
    !!       off "percolating" line-sinks.
    !!   3)  Report any "check" information about the AEMution to this
    !!       point
    !!
    !! Calling Sequence:
    !!    call AEM_Update(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]

    ! Aquifer module
    call AQU_Update(io, aem%aqu, aem%fdp)
    ! LS1 module
    call LS1_Update(io, aem%ls1, aem%fwl, aem%fdp)
    ! LS2 module
    call LS2_Update(io, aem%ls2, aem%fwl, aem%fdp)
    ! LS3 module
    call LS3_Update(io, aem%ls3, aem%fwl, aem%fdp)
    ! HB0 module
    call HB0_Update(io, aem%hb0, aem%fdp)
    ! WL0 module (for adjustable-discharge wells)
    call WL0_Update(io, aem%wl0, aem%fwl)
    ! WL1 module
    call WL1_Update(io, aem%wl1, aem%fwl)
#ifndef __GPL__
    ! CW0 module
    call CW0_Update(io, aem%cw0, aem%fwl, aem%fdp)
#endif

    ! Add new elements here!

    return
  end subroutine AEM_Update


  function AEM_UpdateFreeSurface(io, aem) result(iUpdate)
    !! Adjusts the boundary condition elevations for the free surface BDY elements. This will
    !! potentially modify the elevation of the control intervals for a free surface (the intervals
    !! for the no-flow condition), and may change some boundary conditions from HEAD to FLUX
    !! specification.
    !!
    !! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem
    !! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iUpdate
    !! [ LOCALS ]
    integer(kind=AE_INT) :: ibdy, iNCP, iEqnType, iElementID, iElementString, iElementVertex, iElementFlag
    real(kind=AE_REAL) :: ph0, ph1, ph2, rzn, rGhbFactor, rMaxMove, rdz
    complex(kind=AE_REAL), dimension(2)  :: cCPZ
    complex(kind=AE_REAL) :: cOrientation
    type(AQU_COLLECTION), pointer :: aqu
    type(AQU_BDYELEMENT), pointer :: this, next, prev

    iUpdate = 1

    ! Only runs in PROFILE mode
    if ( .not. io%lProfile ) return

    ! Set a convenience pointer and make all the elements fair game for moving.
    aqu => aem%aqu
    ! We'll only force a matrix regen if a point gets moved
    aqu%lFSRegen = .false.
    ! Initially, use the saved specified head (in case we're iterating on a seepage face)
    aqu%BdyElements%rSpecHead = aqu%BdyElements%rSaveSpecHead
    ! Clear the head-specified flag (we'll recheck this below)
    aqu%BdyElements%lFSHeadSpec = .false.
    ! And for starters, assume we need to move _all_ the ends (we'll identify that also)
    aqu%BdyElements%lMoveCPZ1 = .true.
    aqu%BdyElements%lMoveCPZ2 = .true.
    ! Step 1 -- identify the points that might need to be moved...
    do ibdy=1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      call AQU_GetNeighborBdy(io, aqu, this, prev, next)
      if (this%iBdyFlag /= BDY_FREESURF) then
        this%lMoveCPZ1 = .false.
        prev%lMoveCPZ2 = .false.
        this%lMoveCPZ2 = .false.
        next%lMoveCPZ1 = .false.
      else if (prev%iBdyFlag /= BDY_FREESURF) then
        this%lMoveCPZ1 = .false.
        prev%lMoveCPZ2 = .false.
      else if (next%iBdyFlag /= BDY_FREESURF) then
        this%lMoveCPZ2 = .false.
        next%lMoveCPZ1 = .false.
      end if
    end do

    ! Here we traverse all BDY elements that _may_ require a move, and adjust
    ! the boundary condition elevations for each, if necessary. For efficiency,
    ! if this%cCPZ1 moves, so does prev%cCPZ2, and so forth. When a point is
    ! moved, the flag for that point is set to .false. This means that for each
    ! segment end that moves, a neighboring segment end also moves; thus, this
    ! loop is NOT PARALLEL-SAFE.
    call IO_MessageText(io, "Scanning for head-specified conditions...")
    do ibdy=1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      if (this%iBdyFlag == BDY_FREESURF) then
        call AQU_GetNeighborBdy(io, aqu, this, prev, next)
        ! If the head at the center of the element is > 0, the segment should be
        ! head-specified and both points are on the left side of the element.
        if (this%rSpecHead > max(aimag(this%cFSZ1), aimag(this%cFSZ2))) then
          this%lFSHeadSpec = .true.
          ! First end...
          this%cCPZ1 = this%cFSZ1
          this%lMoveCPZ1 = .false.
          prev%cCPZ2 = prev%cFSZ2
          prev%lMoveCPZ2 = .false.
          ! Second end...
          this%cCPZ2 = this%cFSZ2
          this%lMoveCPZ2 = .false.
          next%cCPZ1 = next%cFSZ1
          next%lMoveCPZ2 = .false.
        end if
      end if
    end do

    call IO_MessageText(io, "Moving points on free surface...")
    rMaxMove = rZERO
    do ibdy=1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      if (this%iBdyFlag == BDY_FREESURF) then
        call AQU_GetNeighborBdy(io, aqu, this, prev, next)
        ! first end
        if (this%lMoveCPZ1) then
          ph0 = rAEM_Head(io, aem, this%cCPZ1) - aimag(this%cCPZ1)
          if (ph0>rZERO .and. abs(aimag(this%cFSZ1)-aimag(this%cCPZ1)) < 1.0e-4_AE_REAL) then
            print *, "Point 1 is on boundary with positive pressure head", ibdy
          else
            rZn = rAEM_FindPressureHeadZero(io, aem, real(this%cCPZ1), aimag(this%cCPZ1), 0.01_AE_REAL, 5, aimag(this%cFSZ1))
            !print *, 'pt1', ibdy, aimag(this%cCPZ1),'zero', rzn
            rZn = 0.5_AE_Real*aimag(this%cCPZ1) + 0.5_AE_Real*rZn
            rdz = abs(rZn-aimag(this%cCPZ1))
            !print *, 'pt1', ibdy, 'moved', rdz
            if (rdz>rMaxMove) rMaxMove = rdz
            this%cCPZ1 = cmplx(real(this%cFSZ1), rZN, AE_REAL)
            this%lMoveCPZ1 = .false.
          end if
        end if
        ! Second end
        if (this%lMoveCPZ2) then
          ph0 = rAEM_Head(io, aem, this%cCPZ2) - aimag(this%cCPZ2)
          if (ph0>rZERO .and. abs(aimag(this%cFSZ2)-aimag(this%cCPZ2)) < 1.0e-4_AE_REAL) then
            print *, "Point 2 is on boundary with positive pressure head", ibdy
          else
            rZn = rAEM_FindPressureHeadZero(io, aem, real(this%cCPZ2), aimag(this%cCPZ2), 0.01_AE_REAL, 5, aimag(this%cFSZ2))
            !print *, 'pt2', ibdy, aimag(this%cCPZ2), 'zero', rzn
            rZn = 0.5_AE_Real*aimag(this%cCPZ2) + 0.5_AE_Real*rZn
            rdz = abs(rZn-aimag(this%cCPZ2))
            !print *, 'pt2', ibdy, 'moved', rdz
            if (rdz>rMaxMove) rMaxMove = rdz
            this%cCPZ2 = cmplx(real(this%cFSZ2), rZN, AE_REAL)
            this%lMoveCPZ2 = .false.
          end if
        end if
      end if
    end do
    write (IO_MessageBuffer, "('Largest vertical move on free surface: ',f10.4)") rMaxMove
    call IO_MessageText(io)

    ! Finally, check the moved elements to see if there are any seepage face elements
    do ibdy=1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      if (this%iBdyFlag == BDY_FREESURF) then
        if (.not. this%lFSHeadSpec) then
          call AQU_GetNeighborBdy(io, aqu, this, prev, next)
          if (abs(aimag(this%cFSZ1)-aimag(this%cCPZ1)) < 1.0e-4_AE_REAL .and. &
              abs(aimag(this%cFSZ2)-aimag(this%cCPZ2)) < 1.0e-4_AE_REAL) then
            write(unit=IO_MessageBuffer, fmt="('Identified seepage face at index ',i10)") ibdy
            call IO_MessageText(io)
            this%lFSHeadSpec = .true.
            this%rSpecHead = rHALF * aimag(this%cFSZ1+this%cFSZ2)
          end if
          !write(unit=*,fmt="('new', i10, l1, 1x, 2(1x,'(',f12.4,1x,f12.4,')'))") ibdy, this%lFSHeadSpec, this%cCPZ1, this%cCPZ2
          !print *,'Check', ibdy, abs(this%cCPZ1-prev%cCPZ2), aimag(this%cFSZ1-this%cCPZ1)
          !print *,'     ', ibdy, abs(this%cCPZ2-next%cCPZ1), aimag(this%cFSZ2-this%cCPZ2)
        end if
      end if
    end do

    aqu%lFSRegen = .true.
    iUpdate = 1

    return
  end function AEM_UpdateFreeSurface


  function rAEM_FindPressureHeadZero(io, aem, rX, rZ0, rTol, iMaxIter, rZLimit) result(rZn)
    !! Finds rZn, the elevation where the pressure head is zero at the x-coordinate rX
    !! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem
    real(kind=AE_REAL), intent(in) :: rX
    real(kind=AE_REAL), intent(in) :: rZ0
    real(kind=AE_REAL), intent(in) :: rTol
    real(kind=AE_REAL), intent(in) :: rZLimit
    integer(kind=AE_INT), intent(in) :: iMaxIter
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rZn
    !! [ LOCALS ]
    real(kind=AE_REAL) :: rdz, rph0, rph1, rph2, rchange, rderiv
    integer(kind=AE_INT) :: iter

    ! Set the delta for the derivative to a fraction of the tolerance
    rdz = 0.01_AE_REAL * rTol
    rZn = rz0
    do iter=1, iMaxIter
      rph0 = rAEM_Head(io, aem, cmplx(rX, rZn+rdz, AE_REAL)) - (rzn+rdz)
      rph1 = rAEM_Head(io, aem, cmplx(rX, rZn-rdz, AE_REAL)) - (rZn-rdz)
      !rph2 = rAEM_Head(io, aem, cmplx(rX, rZn+rdz, AE_REAL)) - (rZn+rdz)
      rderiv = (rph0 - rph1) / (2.0*rdz)
      rchange = rderiv * rph0
      rZn = rZn - rchange
      if (abs(rchange) < rTol) exit
    end do

    if (rZn > rZLimit) then
      rZn = rZLimit
      print *,"HIT Z LIMIT", rZLimit
    end if

    return
  end function rAEM_FindPressureHeadZero


  subroutine AEM_EnableDrawdown(io, aem)
    !! Enables the elements that are subject to the "drawdown" flag
    !! This allows the modeler to compute "drawdown" simulations based on the difference between
    !! an "unstressed" condition and a "stressed" condition. Each element module that implements
    !! the DDN operator provides an XXX_EnableDrawdown() subroutine that makes the necessary
    !! adjustments.
    !!
    !! Returns the number of changes made when enabling the drawdown simulation
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    call IO_MessageText(io, "Enabling the drawdown options")
    iChanges = 0
    iChanges = iChanges + WL0_EnableDrawdown(io, aem%wl0)
    iChanges = iChanges + CW0_EnableDrawdown(io, aem%cw0)
    write(unit=IO_MessageBuffer, fmt=*) "Number of changes: ", iChanges
    call IO_MessageText(io)

    if (iChanges > 0) then
      call AEM_Update(io, aem)
      call AEM_ComputeCheck(io, aem, .false.)
    end if

    return
  end subroutine AEM_EnableDrawdown


  subroutine AEM_ComputeCheck(io, aem, lLinearize)
    !! subroutine AEM_ComputeCheck
    !!
    !! Updates the check information for all equations
    !!
    !! Calling Sequence:
    !!    call AEM_Update(io, aem)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!
    !! If lLinearize is .true., re-linearizes the problem as necessary
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(ITERATOR_RESULT), pointer :: itr

    ! Uses the iterator mechanism to supply check information to all modules

    ! AQU module
    call AQU_ResetIterator(io, aem%aqu)
    do
      itr => AQU_NextIterator(io, aem%aqu)
      if (.not. associated(itr)) exit
      call AQU_SetIterator(io, aem%aqu, aem%fdp, itr, AEM_CheckValue(io, aem, itr), lLinearize)
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! WL0 module
    call WL0_ResetIterator(io, aem%wl0)
    do
      itr => WL0_NextIterator(io, aem%wl0)
      if (.not. associated(itr)) exit
      call WL0_SetIterator(io, aem%wl0, aem%aqu, aem%fwl, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! WL1 module
    call WL1_ResetIterator(io, aem%wl1)
    do
      itr => WL1_NextIterator(io, aem%wl1)
      if (.not. associated(itr)) exit
      call WL1_SetIterator(io, aem%wl1, aem%aqu, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! LS0 module
    call LS0_ResetIterator(io, aem%ls0)
    do
      itr => LS0_NextIterator(io, aem%ls0)
      if (.not. associated(itr)) exit
      call LS0_SetIterator(io, aem%ls0, aem%aqu, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! LS1 module
    call LS1_ResetIterator(io, aem%ls1)
    do
      itr => LS1_NextIterator(io, aem%ls1)
      if (.not. associated(itr)) exit
      call LS1_SetIterator(io, aem%ls1, aem%aqu, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! LS2 module
    call LS2_ResetIterator(io, aem%ls2)
    do
      itr => LS2_NextIterator(io, aem%ls2)
      if (.not. associated(itr)) exit
      call LS2_SetIterator(io, aem%ls2, aem%aqu, itr, AEM_CheckValue(io, aem, itr), lLinearize)
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! LS3 module
    call LS3_ResetIterator(io, aem%ls3)
    do
      itr => LS3_NextIterator(io, aem%ls3)
      if (.not. associated(itr)) exit
      call LS3_SetIterator(io, aem%ls3, aem%aqu, itr, AEM_CheckValue(io, aem, itr), lLinearize)
      deallocate(itr%cZ)
      deallocate(itr)
    end do
    ! HB0 module
    call HB0_ResetIterator(io, aem%hb0)
    do
      itr => HB0_NextIterator(io, aem%hb0)
      if (.not. associated(itr)) exit
      call HB0_SetIterator(io, aem%hb0, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
#ifndef __GPL__
    ! CW0 module
    call CW0_ResetIterator(io, aem%cw0)
    do
      itr => CW0_NextIterator(io, aem%cw0)
      if (.not. associated(itr)) exit
      call CW0_SetIterator(io, aem%cw0, itr, AEM_CheckValue(io, aem, itr))
      deallocate(itr%cZ)
      deallocate(itr)
    end do
#endif

    return
  end subroutine AEM_ComputeCheck


  function AEM_CheckValue(io, aem, itr) result(cValue)
    !! subroutine AEM_CheckValue
    !!
    !! Computes the check value requested in the ITERATOR_RESULT itr
    !!
    !! Calling Sequence:
    !!    call AEM_CheckValue(io, itr)
    !!
    !! Arguments:
    !!    (in)    type(AEM_DOMAIN), pointer :: aem
    !!              The AEM_DOMAIN object to be used
    !!   (in)    type(ITERATOR_RESULT), pointer :: itr
    !!             ITERATOR_RESULT retrieved with XXX_NextIterator()
    !!   (in)    type(IO_STATUS), pointer :: io
    !!             Tracks error conditions(if any)
    !!
    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(ITERATOR_RESULT), pointer :: itr
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cValue

    select case (itr%iValueSelector)
      case (VALUE_HEAD)
        cValue = rAEM_Head(io, aem, itr%cZ(1), .false.)
      case (VALUE_POTENTIAL)
        cValue = cAEM_Potential(io, aem, itr%cZ(1), .false.)
      case (VALUE_POTENTIALDIFF)
        cValue = cAEM_Potential(io, aem, itr%cZ(2), .false.) - cAEM_Potential(io, aem, itr%cZ(1), .false.)
      case (VALUE_DISCHARGE)
        cValue = cAEM_Discharge(io, aem, itr%cZ(1), .false.)
      case (VALUE_VELOCITY)
        cValue = cAEM_Velocity(io, aem, itr%cZ(1), .false.)
      case (VALUE_RECHARGE)
        cValue = rAEM_Recharge(io, aem, itr%cZ(1), .false.)
      case (VALUE_FLOW)
        cValue = rAEM_Flow(io, aem, itr%cZ, .false.)
      case (VALUE_SATDTHICK)
        cValue = rAQU_SatdThickness(io, aem%aqu, itr%cZ(1), real(cAEM_Potential(io, aem, itr%cZ(1), .false.)))
      case (VALUE_TRANSMISSIVITY)
        cValue = rAQU_Transmissivity(io, aem%aqu, itr%cZ(1), real(cAEM_Potential(io, aem, itr%cZ(1), .false.)))
      case (VALUE_EXTRACTION)
        cValue = rAEM_Extraction(io, aem, .false.)
    end select

    return
  end function AEM_CheckValue


  subroutine AEM_Read(io, aem)
    !! Reads an input file(including processing directives) and returns
    !! a populated AEM_DOMAIN object.

    ! [ ARGUMENTS ]
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    ! Command parsing directives
    integer(kind=AE_INT), parameter :: iEND = 1000, iDBG = 1001, &
                            iAQU = 2001, iWL0 = 2002, iPD0 = 2003, iWL1 = 2004, &
                            iLS0 = 2005, iLS1 = 2006, iLS2 = 2007, iLS3 = 2008, &
                            iHB0 = 2009, iAS0 = 20010, iOPT = 2011, &
#ifndef __GPL__
                            iCW0 = 5001, &
#endif
                            iRLX = 3001
#ifndef __GPL__
    type(DIRECTIVE), dimension(15), parameter :: dirDirectives = (/ &
                       DIRECTIVE(iCW0, 'CW0'), &
                       DIRECTIVE(iEND, 'END'), &
                       DIRECTIVE(iDBG, 'DBG'), &
                       DIRECTIVE(iAQU, 'AQU'), &
                       DIRECTIVE(iWL0, 'WL0'), &
                       DIRECTIVE(iPD0, 'PD0'), &
                       DIRECTIVE(iWL1, 'WL1'), &
                       DIRECTIVE(iLS0, 'LS0'), &
                       DIRECTIVE(iLS1, 'LS1'), &
                       DIRECTIVE(iLS2, 'LS2'), &
                       DIRECTIVE(iLS3, 'LS3'), &
                       DIRECTIVE(iHB0, 'HB0'), &
                       DIRECTIVE(iAS0, 'AS0'), &
                       DIRECTIVE(iOPT, 'OPT'), &
                       DIRECTIVE(iRLX, 'RLX')/)
#else
    type(DIRECTIVE), dimension(14), parameter :: dirDirectives = (/ &
                       DIRECTIVE(iEND, 'END'), &
                       DIRECTIVE(iDBG, 'DBG'), &
                       DIRECTIVE(iAQU, 'AQU'), &
                       DIRECTIVE(iWL0, 'WL0'), &
                       DIRECTIVE(iPD0, 'PD0'), &
                       DIRECTIVE(iWL1, 'WL1'), &
                       DIRECTIVE(iLS0, 'LS0'), &
                       DIRECTIVE(iLS1, 'LS1'), &
                       DIRECTIVE(iLS2, 'LS2'), &
                       DIRECTIVE(iLS3, 'LS3'), &
                       DIRECTIVE(iHB0, 'HB0'), &
                       DIRECTIVE(iAS0, 'AS0'), &
                       DIRECTIVE(iOPT, 'OPT'), &
                       DIRECTIVE(iRLX, 'RLX')/)
#endif
    ! Other locals
    character(len=255) :: sOptionText
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iretval
    integer(kind=AE_INT) :: iAS0Flag
    integer(kind=AE_INT) :: iNAS0

    ! Placeholders for function module test calls
    complex(kind=AE_REAL) :: cZ1, cZ2, cZC, cZE1, cZE2
    real(kind=AE_REAL) :: rR, rTol, rDPhiDX, rDPhiDY
    real(kind=AE_REAL) :: rBase, rThick, rHydCond, rPorosity, rAvgHead, rRelax
    logical :: lFlag

    call IO_MessageText(io, "Reading AEM module input")

    call IO_Assert(io, (associated(aem)), &
         "AEM_Read: the AEM_DOMAIN has not been created")

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
        case (iEND)
          ! END mark was found. Exit the file parser.
          exit
        case (iDBG)
          ! Change the io%lDebug flag
          io%lDebug = lIO_GetLogical(io, "lDebug", def=.true.)
        case (iAQU)
          ! Read the infinite aquifer properties and then enter the AQU element module
          aem%aqu => AQU_Create(io)
          call AQU_Read(io, aem%aqu)
        case (iWL0)
          ! Enter the WL0 element module
          call WL0_Alloc(io, aem%wl0)
          call WL0_Read(io, aem%wl0)
        case (iWL1)
          ! Enter the WL1 element module
          call WL1_Alloc(io, aem%wl1)
          call WL1_Read(io, aem%wl1)
        case (iPD0)
          ! Enter the PD0 element module
          call PD0_Alloc(io, aem%pd0)
          call PD0_Read(io, aem%pd0)
        case (iLS0)
          ! Enter the LS0 element module
          call LS0_Alloc(io, aem%ls0)
          call LS0_Read(io, aem%ls0)
        case (iLS1)
          ! Enter the LS1 element module
          call LS1_Alloc(io, aem%ls1)
          call LS1_Read(io, aem%ls1)
        case (iLS2)
          ! Enter the LS2 element module
          call LS2_Alloc(io, aem%ls2)
          call LS2_Read(io, aem%ls2)
        case (iLS3)
          ! Enter the LS3 element module
          call LS3_Alloc(io, aem%ls3)
          call LS3_Read(io, aem%ls3)
        case (iHB0)
          ! Enter the HB0 element module
          call HB0_Alloc(io, aem%hb0)
          call HB0_Read(io, aem%hb0)
        case (iAS0)
          ! Enter the AS0 element module(note: the iAS0Flag is AS0_TOP for the top of
          ! the aquifer and AS0_BOTTOM for the bottom of the aquifer
          !
          ! NOTE that the top/bottom options are entered backwards in ModAEM versions prior to 2.0.
          ! Switch the two arguments for backwards compatibility...
          call IO_SwitchFields(io)
          if (iIO_GetInteger(io, "top/bottom", def=0, allowed=(/0, 1/)) == 0) then
            call AS0_Alloc(io, aem%as0_top)
            call AS0_Read(io, aem%as0_top)
          else
            call AS0_Alloc(io, aem%as0_bottom)
            call AS0_Read(io, aem%as0_bottom)
          end if
        case (iOPT)
          ! Ignore OPT command
          call IO_MessageText(io, '[OPT directive was ignored]')
#ifndef __GPL__
        case (iCW0)
          ! Enter the CW0 element module
          call CW0_Alloc(io, aem%cw0)
          call CW0_Read(io, aem%cw0)
#endif
      end select
    end do

    call IO_MessageText(io, "Leaving AEM module")

    return
  end subroutine AEM_Read


  subroutine AEM_Report(io, aem)
    ! The AEMver report shows the error at the boundary conditions specified in the matrix
    ! Arguments
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    ! Locals
    integer(kind=AE_INT) :: iRow
    complex(kind=AE_REAL), dimension(10) :: cCPZ
    integer(kind=AE_INT) :: iNCP, iEqType, iElementType
    integer(kind=AE_INT) :: iElementString, iElementVertex, iElementFlag
    complex(kind=AE_REAL) :: cOrientation
    real(kind=AE_REAL) :: rGhbDistance

    call HTML_Header('Module AEM', 1)
    call HTML_Header('General Information', 2)
    call HTML_StartTable()
    call HTML_AttrInteger('Current Iteration', aem%iCurrentIteration)
    call HTML_EndTable()
    call HTML_Header('Matrix Equation Information', 3)
    if (aem%mat%iNEqn < 1) then
      call HTML_Header('No equations defined', 4)
    else
      call HTML_StartTable()
      call HTML_StartRow()
      call HTML_TableHeader((/'Row     ', 'EqType  ', 'ElemType', 'String  ', 'Vertex  ', 'Flag    '/))
      call HTML_EndRow()

      do iRow = 1, aem%mat%iNEqn
        ! Get the matrix generator information for the row
        call MAT_GetEquation(io, aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, iElementString, &
             iElementVertex, iElementFlag, cOrientation, rGhbDistance)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iRow, iEqType, iElementType, iElementString, iElementVertex, iElementFlag/))
        call HTML_EndRow()
      end do
    end if
    call HTML_EndTable()

  end subroutine AEM_Report


  subroutine AEM_Save(io, aem, sFname, mode)
    !! Saves the current solution, creating the provided file
    !! The "save" file contains only the un
    ! [ARGUMENTS]
    type(AEM_DOMAIN), pointer :: aem
    character(len=*), intent(in) :: sFname
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in) :: mode
    ! [LOCALS]
    integer(kind=AE_INT) :: istat

    if (mode == IO_MODE_BINARY) then
      open(unit=LU_SCRATCH, file=trim(sFname), status="NEW", form="UNFORMATTED", iostat=istat)
    else
      open(unit=LU_SCRATCH, file=trim(sFname), status="NEW", iostat=istat)
    end if
    call IO_Assert(io, istat==0, "Open failed on output file " // trim(sFname))

    call AQU_Save(io, aem%aqu, mode)
    call LS1_Save(io, aem%ls1, mode)
    call LS2_Save(io, aem%ls2, mode)
    call LS3_Save(io, aem%ls3, mode)
    call HB0_Save(io, aem%hb0, mode)
    call WL1_Save(io, aem%wl1, mode)
#ifndef __GPL__
    call CW0_Save(io, aem%cw0, mode)
#endif

    close(unit=LU_SCRATCH)

    return
  end subroutine AEM_Save


  subroutine AEM_Load(io, aem, sFname, mode)
    !! Loads a previously-saved preconditioning file from disk
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem
    character(len=*), intent(in) :: sFname
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat

    if (mode == IO_MODE_BINARY) then
      open(unit=LU_SCRATCH, file=trim(sFname), status="OLD", form="UNFORMATTED", iostat=istat)
    else
      open(unit=LU_SCRATCH, file=trim(sFname), status="OLD", iostat=istat)
    end if
    call IO_Assert(io, istat==0, "Open failed on intput file " // trim(sFname))

    call AQU_Load(io, aem%aqu, aem%fdp, mode)
    call LS1_Load(io, aem%ls1, aem%fwl, aem%fdp, mode)
    call LS2_Load(io, aem%ls2, aem%fwl, aem%fdp, mode)
    call LS3_Load(io, aem%ls3, aem%fwl, aem%fdp, mode)
    call HB0_Load(io, aem%hb0, aem%fdp, mode)
    call WL1_Load(io, aem%wl1, aem%fwl, mode)
#ifndef __GPL__
    call CW0_Load(io, aem%cw0, aem%fwl, aem%fdp, mode)
#endif

    close(unit=LU_SCRATCH)

    return
  end subroutine AEM_Load

end module m_aem
