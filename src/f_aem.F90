module f_aem

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

  ! AEM computational kernel for ModAEM. Holds function collections and provides
  ! potential/discharge/flow computations. All element-package orchestration has
  ! been moved to m_packages.

  use u_constants
  use u_io
  use u_matrix
  use u_domain
  use f_reference
  use f_well
  use f_dipole
  use f_pond
  use f_linesink
  use f_areasink

  implicit none

  public

  type :: AEM_DOMAIN
    !! type AEM_DOMAIN
    !!
    !! Computational kernel: holds function collections, matrix, domain, and
    !! reference-function objects. Element-package collections live in PKG_DOMAIN.
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
    integer(kind=AE_INT) :: iCW0Start
    integer(kind=AE_INT) :: iCW0NUnk

    ! True if a solution is present
    logical :: lSolutionPresent

    ! True if the drawdown elements are enabled
    logical :: lDrawdown

    ! Domain, reference function, and function-collection objects
    type(DOM_COLLECTION), pointer :: dom
    type(FRF_COLLECTION), pointer :: frf
    type(FWL_COLLECTION), pointer :: fwl
    type(FPD_COLLECTION), pointer :: fpd
    type(FDP_COLLECTION), pointer :: fdp
    type(FLS_COLLECTION), pointer :: fls
    type(FAS_COLLECTION), pointer :: fas_top
    type(FAS_COLLECTION), pointer :: fas_bottom
    type(MAT_MATRIX), pointer :: mat
  end type AEM_DOMAIN

  real(kind=AE_REAL), private, parameter :: MOVEPOINT = 1.0e-3_AE_REAL

contains


  function AEM_Create(io) result(aem)
    !! Creates a new AEM_DOMAIN object (function collections only; element
    !! collections live in PKG_DOMAIN and are created by PKG_Create).
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iStat

    allocate(aem, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AEM_Create: allocation failed")

    aem%iCurrentIteration = 0
    aem%iAQUStart = 0
    aem%iAQUNUnk = 0
    aem%iLS1Start = 0
    aem%iLS1NUnk = 0
    aem%iLS2Start = 0
    aem%iLS2NUnk = 0
    aem%iLS3Start = 0
    aem%iLS3NUnk = 0
    aem%iHB0Start = 0
    aem%iHB0NUnk = 0
    aem%iWL1Start = 0
    aem%iWL1NUnk = 0
    aem%lDrawdown = .false.
    aem%lSolutionPresent = .false.

    nullify(aem%dom)
    nullify(aem%frf)
    nullify(aem%fwl)
    nullify(aem%fpd)
    nullify(aem%fdp)
    nullify(aem%fls)
    nullify(aem%fas_top)
    nullify(aem%fas_bottom)
    aem%mat => MAT_Create(io)

    return
  end function AEM_Create


  subroutine AEM_Destroy(io, aem)
    !! Frees memory allocated for an AEM_DOMAIN object.
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iStat

    call MAT_Destroy(io, aem%mat)
    deallocate(aem, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AEM_Destroy: deallocation failed")

    return
  end subroutine AEM_Destroy


  function cAEM_Potential(io, aem, cZ, lNoCheck) result(cOmega)
    !! Returns the complex potential at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    complex(kind=AE_REAL) :: cOmega
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
        else if (lFLS_CheckPoint(io, aem%fls, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cOmega = cFRF_Potential(io, aem%frf, cZArg) + &
             cFWL_Potential(io, aem%fwl, cZArg) + &
             cFPD_Potential(io, aem%fpd, cZArg) + &
             cFDP_Potential(io, aem%fdp, cZArg) + &
             cFLS_Potential(io, aem%fls, cZArg) + &
             rFAS_InsidePotential(io, aem%fas_top, cZArg) + &
             rFAS_InsidePotential(io, aem%fas_bottom, cZArg)

    return
  end function cAEM_Potential


  function cAEM_Discharge(io, aem, cZ, lNoCheck) result(cQ)
    !! Returns the discharge vector at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    complex(kind=AE_REAL) :: cQ
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
        else if (lFLS_CheckPoint(io, aem%fls, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cQ = cFRF_Discharge(io, aem%frf, cZArg) + &
         cFWL_Discharge(io, aem%fwl, cZArg) + &
         cFPD_Discharge(io, aem%fpd, cZArg) + &
         cFDP_Discharge(io, aem%fdp, cZArg) + &
         cFLS_Discharge(io, aem%fls, cZArg) + &
         cFAS_InsideDischarge(io, aem%fas_top, cZArg) + &
         cFAS_InsideDischarge(io, aem%fas_bottom, cZArg) + &
         cZERO

    return
  end function cAEM_Discharge


  function rAEM_Recharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! Returns the net recharge rate at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rGamma
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ

    rGamma = rFWL_Recharge(io, aem%fwl, cZArg) + &
             rFPD_Recharge(io, aem%fpd, cZArg) + &
             rFDP_Recharge(io, aem%fdp, cZArg) + &
             rFAS_InsideRecharge(io, aem%fas_top, cZArg) + &
             rFAS_InsideRecharge(io, aem%fas_bottom, cZArg)

    return
  end function rAEM_Recharge


  function rAEM_TopRecharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! Returns the top recharge rate at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rGamma
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ

    rGamma = rFWL_Recharge(io, aem%fwl, cZArg) + &
             rFPD_Recharge(io, aem%fpd, cZArg) + &
             rFDP_Recharge(io, aem%fdp, cZArg) + &
             rFAS_InsideRecharge(io, aem%fas_top, cZArg)

    return
  end function rAEM_TopRecharge


  function rAEM_BottomRecharge(io, aem, cZ, lNoCheck) result(rGamma)
    !! Returns the bottom recharge rate at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rGamma
    complex(kind=AE_REAL) :: cZArg

    cZArg = cZ

    rGamma = rFAS_InsideRecharge(io, aem%fas_bottom, cZArg)

    return
  end function rAEM_BottomRecharge


  function rAEM_Extraction(io, aem, lNoCheck) result(rQ)
    !! Returns the net extraction rate of the model
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rQ

    rQ = rFWL_Extraction(io, aem%fwl) + &
         rFPD_Extraction(io, aem%fpd) + &
         rFDP_Extraction(io, aem%fdp) + &
         rFLS_Extraction(io, aem%fls) + &
         rFAS_Extraction(io, aem%fas_top) - &
         rFAS_Extraction(io, aem%fas_bottom)

    return
  end function rAEM_Extraction


  function rAEM_HeadAtWell(io, aem, cZ, pFWL, lNoCheck) result(rH)
    !! Returns the head at cZ, excluding the contribution of pFWL
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(FWL_WELL), pointer :: pFWL
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rH
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
        else if (lFLS_CheckPoint(io, aem%fls, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    rH = rDOM_PotentialToHead(io, aem%dom, &
                              real(cAEM_Potential(io, aem, cZArg) - &
                                   cFWL_Potential(io, aem%fwl, cZArg, pFWL%iIndex, 1)), &
                              cZArg)

    return
  end function rAEM_HeadAtWell


  function cAEM_DischargeAtWell(io, aem, cZ, pFWL, lNoCheck) result(cQ)
    !! Returns the discharge vector at cZ, excluding the contribution of pFWL
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(FWL_WELL), pointer :: pFWL
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    complex(kind=AE_REAL) :: cQ
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
        else if (lFLS_CheckPoint(io, aem%fls, cZArg, rVERTEXTOL, cZFix, rStrength, &
               iElementType, iElementString, iElementVertex, iElementFlag)) then
          cZArg = cZFix
          cycle
        else
          exit
        end if
      end do
    end if

    cQ = cAEM_Discharge(io, aem, cZArg) - &
         cFWL_Discharge(io, aem%fwl, cZArg, pFWL%iIndex, 1)

    return
  end function cAEM_DischargeAtWell


  function rAEM_Head(io, aem, cZ, lNoCheck) result(rHead)
    !! Returns the head at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rHead

    rHead = rDOM_PotentialToHead(io, aem%dom, real(cAEM_Potential(io, aem, cZ, lNoCheck)), cZ)

    return
  end function rAEM_Head


  function rAEM_SatdThick(io, aem, cZ, lNoCheck) result(rSatdThick)
    !! Returns the saturated thickness at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rSatdThick

    rSatdThick = rDOM_SatdThickness(io, aem%dom, cZ, real(cAEM_Potential(io, aem, cZ, lNoCheck)))

    return
  end function rAEM_SatdThick


  function rAEM_InterfaceElevation(io, aem, cZ) result(rIfcElev)
    !! Returns the aquifer base elevation at cZ (sea-water interface stub)
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    real(kind=AE_REAL) :: rIfcElev

    rIfcElev = rDOM_Base(io, aem%dom, cZ)

    return
  end function rAEM_InterfaceElevation


  function rAEM_Flow(io, aem, cZ, lNoCheck) result(rFlow)
    !! Returns the integrated flow across the path cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), dimension(:), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    real(kind=AE_REAL) :: rFlow
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
          else if (lFLS_CheckPoint(io, aem%fls, cZArg(i), rVERTEXTOL, cZFix, rStrength, &
                 iElementType, iElementString, iElementVertex, iElementFlag)) then
            cZArg(i) = cZFix
            cycle
          else
            exit
          end if
        end do
      end do
    end if

    rFlow = rFRF_Flow(io, aem%frf, cZArg) + &
            rFWL_Flow(io, aem%fwl, cZArg) + &
            rFPD_Flow(io, aem%fpd, cZArg) + &
            rFDP_Flow(io, aem%fdp, cZArg) + &
            rFLS_Flow(io, aem%fls, cZArg) + &
            rFAS_InsideFlow(io, aem%fas_top, cZArg) + &
            rFAS_InsideFlow(io, aem%fas_bottom, cZArg)
    deallocate(cZArg)

    return
  end function rAEM_Flow


  function cAEM_Velocity(io, aem, cZ, lNoCheck) result(cV)
    !! Returns the pore velocity at cZ
    type(AEM_DOMAIN), pointer :: aem
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lNoCheck
    complex(kind=AE_REAL) :: cV

    cV = cDOM_DischargeToVelocity(io, aem%dom, cAEM_Discharge(io, aem, cZ, lNoCheck), &
         cZ, real(cAEM_Potential(io, aem, cZ, lNoCheck)))

    return
  end function cAEM_Velocity


  function rAEM_FindPressureHeadZero(io, aem, rX, rZ0, rTol, iMaxIter, rZLimit) result(rZn)
    !! Finds rZn, the elevation where the pressure head is zero at the x-coordinate rX
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem
    real(kind=AE_REAL), intent(in) :: rX
    real(kind=AE_REAL), intent(in) :: rZ0
    real(kind=AE_REAL), intent(in) :: rTol
    real(kind=AE_REAL), intent(in) :: rZLimit
    integer(kind=AE_INT), intent(in) :: iMaxIter
    real(kind=AE_REAL) :: rZn
    real(kind=AE_REAL) :: rdz, rph0, rph1, rchange, rderiv
    integer(kind=AE_INT) :: iter

    rdz = 0.01_AE_REAL * rTol
    rZn = rz0
    do iter=1, iMaxIter
      rph0 = rAEM_Head(io, aem, cmplx(rX, rZn+rdz, AE_REAL)) - (rzn+rdz)
      rph1 = rAEM_Head(io, aem, cmplx(rX, rZn-rdz, AE_REAL)) - (rZn-rdz)
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


  subroutine AEM_Report(io, aem)
    !! Writes an HTML report of matrix equation information
    type(AEM_DOMAIN), pointer :: aem
    type(IO_STATUS), pointer :: io
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
        call MAT_GetEquation(io, aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, iElementString, &
             iElementVertex, iElementFlag, cOrientation, rGhbDistance)
        call HTML_StartRow()
        call HTML_ColumnInteger((/iRow, iEqType, iElementType, iElementString, iElementVertex, iElementFlag/))
        call HTML_EndRow()
      end do
    end if
    call HTML_EndTable()

  end subroutine AEM_Report


end module f_aem
