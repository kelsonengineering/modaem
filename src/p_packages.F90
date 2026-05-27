module p_packages

  ! ModAEM 2.0
  ! Copyright(c) 1995-2008 WHPA Inc. and Vic Kelson
  !
  ! Orchestration module -- holds all element-package collections and the
  ! AEM computational kernel, and provides the Read/Solve/Report lifecycle.

  use u_constants
  use u_io
  use u_matrix
  use u_domain
  use f_aem
  use p_aqu
  use p_wl0
  use p_wl1
  use p_pd0
  use p_ls0
  use p_ls1
  use p_ls2
  use p_ls3
  use p_hb0
  use p_as0
  use p_cw0

  implicit none

  public

  type :: PKG_DOMAIN
    !! type PKG_DOMAIN
    !!
    !! Orchestration container: holds the AEM computational kernel plus all
    !! element-package collections.
    !!
    ! Starting column and unknown-count for each element package in the global matrix
    integer(kind=AE_INT) :: iAQUStart, iAQUNUnk
    integer(kind=AE_INT) :: iLS1Start, iLS1NUnk
    integer(kind=AE_INT) :: iLS2Start, iLS2NUnk
    integer(kind=AE_INT) :: iLS3Start, iLS3NUnk
    integer(kind=AE_INT) :: iHB0Start, iHB0NUnk
    integer(kind=AE_INT) :: iWL1Start, iWL1NUnk
    integer(kind=AE_INT) :: iCW0Start, iCW0NUnk

    type(AEM_DOMAIN), pointer :: aem
    type(AQU_COLLECTION), pointer :: aqu
    type(WL0_COLLECTION), pointer :: wl0
    type(WL1_COLLECTION), pointer :: wl1
    type(PD0_COLLECTION), pointer :: pd0
    type(LS0_COLLECTION), pointer :: ls0
    type(LS1_COLLECTION), pointer :: ls1
    type(LS2_COLLECTION), pointer :: ls2
    type(LS3_COLLECTION), pointer :: ls3
    type(HB0_COLLECTION), pointer :: hb0
    type(AS0_COLLECTION), pointer :: as0_top
    type(AS0_COLLECTION), pointer :: as0_bottom
    type(CW0_COLLECTION), pointer :: cw0
  end type PKG_DOMAIN


contains


  function PKG_Create(io) result(pkg)
    !! Creates a new PKG_DOMAIN: allocates the AEM kernel and all element collections.
    type(IO_STATUS), pointer :: io
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT) :: iStat

    allocate(pkg, stat=iStat)
    call IO_Assert(io, (iStat == 0), "PKG_Create: allocation failed")

    pkg%iAQUStart = 0 ; pkg%iAQUNUnk = 0
    pkg%iLS1Start = 0 ; pkg%iLS1NUnk = 0
    pkg%iLS2Start = 0 ; pkg%iLS2NUnk = 0
    pkg%iLS3Start = 0 ; pkg%iLS3NUnk = 0
    pkg%iHB0Start = 0 ; pkg%iHB0NUnk = 0
    pkg%iWL1Start = 0 ; pkg%iWL1NUnk = 0
    pkg%iCW0Start = 0 ; pkg%iCW0NUnk = 0

    pkg%aem => AEM_Create(io)
    nullify(pkg%aqu)
    pkg%wl0 => WL0_Create(io)
    pkg%pd0 => PD0_Create(io)
    pkg%ls0 => LS0_Create(io)
    pkg%ls1 => LS1_Create(io)
    pkg%ls2 => LS2_Create(io)
    pkg%ls3 => LS3_Create(io)
    pkg%hb0 => HB0_Create(io)
    pkg%wl1 => WL1_Create(io)
    pkg%as0_top => AS0_Create(AS0_TOP)
    pkg%as0_bottom => AS0_Create(AS0_BOTTOM)
    pkg%cw0 => CW0_Create(io)

    return
  end function PKG_Create


  subroutine PKG_Destroy(io, pkg)
    !! Frees all memory associated with a PKG_DOMAIN.
    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iStat

    call WL0_Destroy(io, pkg%wl0)
    call PD0_Destroy(io, pkg%pd0)
    call LS0_Destroy(io, pkg%ls0)
    call LS1_Destroy(io, pkg%ls1)
    call LS2_Destroy(io, pkg%ls2)
    call LS3_Destroy(io, pkg%ls3)
    call HB0_Destroy(io, pkg%hb0)
    call WL1_Destroy(io, pkg%wl1)
    call AS0_Destroy(io, pkg%as0_top)
    call AS0_Destroy(io, pkg%as0_bottom)
    call CW0_Destroy(io, pkg%cw0)
    call AEM_Destroy(io, pkg%aem)
    deallocate(pkg, stat=iStat)
    call IO_Assert(io, (iStat == 0), "PKG_Destroy: deallocation failed")

    return
  end subroutine PKG_Destroy


  subroutine PKG_PreSolve(io, pkg)
    !! Calls all the pre-solve steps.
    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io

    call AQU_PreSolve(io, pkg%aqu)
    call HB0_PreSolve(io, pkg%hb0)
    call WL0_PreSolve(io, pkg%wl0)
    call WL1_PreSolve(io, pkg%wl1)
    call LS0_PreSolve(io, pkg%ls0)
    call LS1_PreSolve(io, pkg%ls1)
    call LS2_PreSolve(io, pkg%ls2)
    call LS3_PreSolve(io, pkg%ls3)
    call PD0_PreSolve(io, pkg%pd0)
    call AS0_PreSolve(io, pkg%as0_top, pkg%aqu)
    call AS0_PreSolve(io, pkg%as0_bottom, pkg%aqu)
    call CW0_PreSolve(io, pkg%cw0)

    return
  end subroutine PKG_PreSolve


  subroutine PKG_AllocateFunctions(io, pkg)
    !! Allocates the space in the FWL, FPD, FDP, FLS, and FAS collections and
    !! calls each element module's SetupFunctions routine.
    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iNFWL, iNFPD, iNFDP, iNFLS
    integer(kind=AE_INT) :: iNFAS_top, iNFAS_bottom
    integer(kind=AE_INT) :: i

    iNFWL = 0
    iNFPD = 0
    iNFDP = 0
    iNFLS = 0
    iNFAS_top = 0
    iNFAS_bottom = 0

    ! Given-strength elements
    iNFWL = iNFWL + iWL0_GetInfo(io, pkg%wl0, SIZE_FWL, 0)
    iNFPD = iNFPD + iWL0_GetInfo(io, pkg%wl0, SIZE_FPD, 0)
    iNFDP = iNFDP + iWL0_GetInfo(io, pkg%wl0, SIZE_FDP, 0)

    iNFWL = iNFWL + iPD0_GetInfo(io, pkg%pd0, SIZE_FWL, 0)
    iNFPD = iNFPD + iPD0_GetInfo(io, pkg%pd0, SIZE_FPD, 0)
    iNFDP = iNFDP + iPD0_GetInfo(io, pkg%pd0, SIZE_FDP, 0)

    iNFLS = iNFLS + iLS0_GetInfo(io, pkg%ls0, SIZE_FLS, 0)

    iNFWL = iNFWL + iAS0_GetInfo(io, pkg%as0_top, SIZE_FWL, 0)
    iNFPD = iNFPD + iAS0_GetInfo(io, pkg%as0_top, SIZE_FPD, 0)
    iNFDP = iNFDP + iAS0_GetInfo(io, pkg%as0_top, SIZE_FDP, 0)
    iNFAS_top = iAS0_GetInfo(io, pkg%as0_top, SIZE_FAS, 0)

    iNFWL = iNFWL + iAS0_GetInfo(io, pkg%as0_bottom, SIZE_FWL, 0)
    iNFPD = iNFPD + iAS0_GetInfo(io, pkg%as0_bottom, SIZE_FPD, 0)
    iNFDP = iNFDP + iAS0_GetInfo(io, pkg%as0_bottom, SIZE_FDP, 0)
    iNFAS_bottom = iAS0_GetInfo(io, pkg%as0_bottom, SIZE_FAS, 0)

    ! Unknown-strength elements
    iNFWL = iNFWL + iAQU_GetInfo(io, pkg%aqu, SIZE_FWL, 0)
    iNFPD = iNFPD + iAQU_GetInfo(io, pkg%aqu, SIZE_FPD, 0)
    iNFDP = iNFDP + iAQU_GetInfo(io, pkg%aqu, SIZE_FDP, 0)
    iNFLS = iNFLS + iAQU_GetInfo(io, pkg%aqu, SIZE_FLS, 0)

    iNFLS = iNFLS + iLS1_GetInfo(io, pkg%ls1, SIZE_FLS, 0)

    iNFLS = iNFLS + iLS2_GetInfo(io, pkg%ls2, SIZE_FLS, 0)

    iNFLS = iNFLS + iLS3_GetInfo(io, pkg%ls3, SIZE_FLS, 0)

    iNFWL = iNFWL + iHB0_GetInfo(io, pkg%hb0, SIZE_FWL, 0)
    iNFPD = iNFPD + iHB0_GetInfo(io, pkg%hb0, SIZE_FPD, 0)
    iNFDP = iNFDP + iHB0_GetInfo(io, pkg%hb0, SIZE_FDP, 0)

    iNFWL = iNFWL + iWL1_GetInfo(io, pkg%wl1, SIZE_FWL, 0)
    iNFPD = iNFPD + iWL1_GetInfo(io, pkg%wl1, SIZE_FPD, 0)
    iNFDP = iNFDP + iWL1_GetInfo(io, pkg%wl1, SIZE_FDP, 0)
    iNFLS = iNFLS + iCW0_GetInfo(io, pkg%cw0, SIZE_FLS, 0)

    ! Allocate function collections in the AEM kernel
    pkg%aem%fwl => FWL_Create(io, iNFWL)
    pkg%aem%fpd => FPD_Create(io, iNFPD)
    pkg%aem%fdp => FDP_Create(io, iNFDP)
    pkg%aem%fls => FLS_Create(io, max(iNFLS, 1))
    pkg%aem%fas_top => FAS_Create(io, max(iNFAS_top, 1))
    pkg%aem%fas_bottom => FAS_Create(io, max(iNFAS_bottom, 1))

    ! Call each element module's setup routine
    call AQU_SetupFunctions(io, pkg%aqu, pkg%aem%fdp, pkg%aem%fls)
    call WL0_SetupFunctions(io, pkg%wl0, pkg%aem%fwl, pkg%aqu)
    call PD0_SetupFunctions(io, pkg%pd0, pkg%aem%fpd)
    call LS0_SetupFunctions(io, pkg%ls0, pkg%aem%fls)
    call LS1_SetupFunctions(io, pkg%ls1, pkg%aem%fls)
    call LS2_SetupFunctions(io, pkg%ls2, pkg%aem%fls)
    call LS3_SetupFunctions(io, pkg%ls3, pkg%aem%fls)
    call HB0_SetupFunctions(io, pkg%hb0, pkg%aem%fdp)
    call WL1_SetupFunctions(io, pkg%wl1, pkg%aem%fwl, pkg%aqu)
    call AS0_SetupFunctions(io, pkg%as0_top, pkg%aem%fwl, pkg%aem%fdp, pkg%aem%fas_top)
    call AS0_SetupFunctions(io, pkg%as0_bottom, pkg%aem%fwl, pkg%aem%fdp, pkg%aem%fas_bottom)
    call CW0_SetupFunctions(io, pkg%cw0, pkg%aem%fls)

    return
  end subroutine PKG_AllocateFunctions


  subroutine PKG_AllocateMatrix(io, pkg, iIteration)
    !! Allocates the space in the MAT module and calls each element's SetupMatrix.
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iNEQ, iNUN
    integer(kind=AE_INT) :: i

    iNEQ = 0
    iNUN = 0

    iNEQ = iNEQ + iAQU_GetInfo(io, pkg%aqu, SIZE_EQUATIONS, iIteration)
    pkg%iAQUStart = 1
    pkg%iAQUNUnk = iAQU_GetInfo(io, pkg%aqu, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iAQUNUnk

    iNEQ = iNEQ + iLS1_GetInfo(io, pkg%ls1, SIZE_EQUATIONS, iIteration)
    pkg%iLS1Start = pkg%iAQUNUnk + 1
    pkg%iLS1NUnk = iLS1_GetInfo(io, pkg%ls1, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iLS1NUnk

    iNEQ = iNEQ + iLS2_GetInfo(io, pkg%ls2, SIZE_EQUATIONS, iIteration)
    pkg%iLS2Start = pkg%iAQUNUnk + pkg%iLS1NUnk + 1
    pkg%iLS2NUnk = iLS2_GetInfo(io, pkg%ls2, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iLS2NUnk

    iNEQ = iNEQ + iLS3_GetInfo(io, pkg%ls3, SIZE_EQUATIONS, iIteration)
    pkg%iLS3Start = pkg%iAQUNUnk + pkg%iLS1NUnk + pkg%iLS2NUnk + 1
    pkg%iLS3NUnk = iLS3_GetInfo(io, pkg%ls3, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iLS3NUnk

    iNEQ = iNEQ + iHB0_GetInfo(io, pkg%hb0, SIZE_EQUATIONS, iIteration)
    pkg%iHB0Start = pkg%iAQUNUnk + pkg%iLS1NUnk + pkg%iLS2NUnk + pkg%iLS3NUnk + 1
    pkg%iHB0NUnk = iHB0_GetInfo(io, pkg%hb0, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iHB0NUnk

    iNEQ = iNEQ + iWL1_GetInfo(io, pkg%wl1, SIZE_EQUATIONS, iIteration)
    pkg%iWL1Start = pkg%iAQUNUnk + pkg%iLS1NUnk + pkg%iLS2NUnk + &
                        pkg%iLS3NUnk + pkg%iHB0NUnk + 1
    pkg%iWL1NUnk = iWL1_GetInfo(io, pkg%wl1, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iWL1NUnk
    iNEQ = iNEQ + iCW0_GetInfo(io, pkg%cw0, SIZE_EQUATIONS, iIteration)
    pkg%iCW0Start = pkg%iAQUNUnk + pkg%iLS1NUnk + pkg%iLS2NUnk + &
                        pkg%iLS3NUnk + pkg%iHB0NUnk + pkg%iWL1NUnk + 1
    pkg%iCW0NUnk = iCW0_GetInfo(io, pkg%cw0, SIZE_UNKNOWNS, iIteration)
    iNUN = iNUN + pkg%iCW0NUnk

    call MAT_Alloc(io, pkg%aem%mat, iNEQ, iNUN)

    if (pkg%iAQUNUnk /= 0) call AQU_SetupMatrix(io, pkg%aqu, pkg%aem%mat)
    if (pkg%iLS1NUnk /= 0) call LS1_SetupMatrix(io, pkg%ls1, pkg%aem%mat)
    if (pkg%iLS2NUnk /= 0) call LS2_SetupMatrix(io, pkg%ls2, pkg%aem%mat)
    if (pkg%iLS3NUnk /= 0) call LS3_SetupMatrix(io, pkg%ls3, pkg%aem%mat)
    if (pkg%iHB0NUnk /= 0) call HB0_SetupMatrix(io, pkg%hb0, pkg%aem%mat)
    if (pkg%iWL1NUnk /= 0) call WL1_SetupMatrix(io, pkg%wl1, pkg%aqu, pkg%aem%mat)
    if (pkg%iCW0NUnk /= 0) call CW0_SetupMatrix(io, pkg%cw0, pkg%aem%mat)

    return
  end subroutine PKG_AllocateMatrix


  subroutine PKG_Update(io, pkg)
    !! Updates all element function strengths after a solution step.
    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io

    call AQU_Update(io, pkg%aqu, pkg%aem%fdp)
    call LS1_Update(io, pkg%ls1, pkg%aem%fls)
    call LS2_Update(io, pkg%ls2, pkg%aem%fls)
    call LS3_Update(io, pkg%ls3, pkg%aem%fls)
    call HB0_Update(io, pkg%hb0, pkg%aem%fdp)
    call WL0_Update(io, pkg%wl0, pkg%aem%fwl)
    call WL1_Update(io, pkg%wl1, pkg%aem%fwl)
    call CW0_Update(io, pkg%cw0, pkg%aem%fls)

    return
  end subroutine PKG_Update


  function rPKG_GetCoefficientMultiplier(io, pkg, iElementType, iElementString, &
                                         iElementVertex, iElementFlag) result(rMultiplier)
    !! Looks up the proper coefficient multiplier for the given element.
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    real(kind=AE_REAL) :: rMultiplier

    select case (iElementType)
      case (ELEM_AQU, ELEM_IN0)
        rMultiplier = rAQU_GetCoefficientMultiplier(io, pkg%aqu, iElementType, &
                                                    iElementString, iElementVertex, iElementFlag)
      case (ELEM_LS1)
        rMultiplier = rLS1_GetCoefficientMultiplier(io, pkg%ls1, iElementString, &
                                                    iElementVertex, iElementFlag)
      case (ELEM_LS2)
        rMultiplier = rLS2_GetCoefficientMultiplier(io, pkg%ls2, iElementString, &
                                                    iElementVertex, iElementFlag)
      case (ELEM_LS3)
        rMultiplier = rLS3_GetCoefficientMultiplier(io, pkg%ls3, iElementString, &
                                                    iElementVertex, iElementFlag)
      case (ELEM_HB0)
        rMultiplier = rHB0_GetCoefficientMultiplier(io, pkg%hb0, iElementString, &
                                                    iElementVertex, iElementFlag)
      case (ELEM_WL1)
        rMultiplier = rWL1_GetCoefficientMultiplier(io, pkg%wl1, iElementString, &
                                                    iElementVertex, iElementFlag)
      case (ELEM_CW0)
        rMultiplier = rCW0_GetCoefficientMultiplier(io, pkg%cw0, iElementString, &
                                                    iElementVertex, iElementFlag)
    end select

    return
  end function rPKG_GetCoefficientMultiplier


  subroutine PKG_GenerateMatrix(io, pkg, iIteration)
    !! Generates the matrix coefficients row-by-row.
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iRow
    real(kind=AE_REAL) :: rMultiplier, rGhbDistance
    real(kind=AE_REAL), dimension(:), allocatable :: rARow
    integer(kind=AE_INT) :: iStat
    complex(kind=AE_REAL), dimension(10) :: cCPZ
    integer(kind=AE_INT) :: iNCP, iEqType, ic
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag
    complex(kind=AE_REAL) :: cOrientation

    allocate(rARow(pkg%aem%mat%iNVar), stat=iStat)
    call IO_Assert(io, (iStat == 0), "PKG_GenerateMatrix: Allocation failed")

    print *, 'Solving a system of ', pkg%aem%mat%iNEqn, ' equations'
    do iRow = 1, pkg%aem%mat%iNEqn
      call MAT_GetEquation(io, pkg%aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, &
           iElementString, iElementVertex, iElementFlag, cOrientation, rGhbDistance)
      rMultiplier = rPKG_GetCoefficientMultiplier(io, pkg, iElementType, iElementString, &
                                                  iElementVertex, iElementFlag)

      if (pkg%iAQUNUnk > 0) then
        call AQU_ComputeCoefficients(io, pkg%aqu, pkg%aem%fdp, (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iAQUStart:pkg%iAQUStart+pkg%iAQUNUnk-1))
      end if

      if (pkg%iLS1NUnk > 0) then
        call LS1_ComputeCoefficients(io, pkg%ls1, pkg%aem%fls, &
             (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iLS1Start:pkg%iLS1Start+pkg%iLS1NUnk-1))
      end if

      if (pkg%iLS2NUnk > 0) then
        call LS2_ComputeCoefficients(io, pkg%ls2, pkg%aem, pkg%aem%fls, &
             (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iLS2Start:pkg%iLS2Start+pkg%iLS2NUnk-1))
      end if

      if (pkg%iLS3NUnk > 0) then
        call LS3_ComputeCoefficients(io, pkg%ls3, pkg%aem, pkg%aem%fls, &
             (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iLS3Start:pkg%iLS3Start+pkg%iLS3NUnk-1))
      end if

      if (pkg%iHB0NUnk > 0) then
        call HB0_ComputeCoefficients(io, pkg%hb0, pkg%aem%fdp, (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iHB0Start:pkg%iHB0Start+pkg%iHB0NUnk-1))
      end if

      if (pkg%iWL1NUnk > 0) then
        call WL1_ComputeCoefficients(io, pkg%wl1, pkg%aem%fwl, (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iWL1Start:pkg%iWL1Start+pkg%iWL1NUnk-1))
      end if

      if (pkg%iCW0NUnk > 0) then
        call CW0_ComputeCoefficients(io, pkg%cw0, pkg%aem, pkg%aem%fls, &
             (/(cCPZ(ic), ic=1, iNCP)/), &
             iEqType, iElementType, iElementString, iElementVertex, iElementFlag, &
             cOrientation, rGhbDistance, iIteration, rMultiplier, &
             rARow(pkg%iCW0Start:pkg%iCW0Start+pkg%iCW0NUnk-1))
      end if

      call MAT_SetRow(io, pkg%aem%mat, iRow, rARow(1:pkg%aem%mat%iNVar), rZERO)

    end do

    deallocate(rARow)

    return
  end subroutine PKG_GenerateMatrix


  subroutine PKG_GenerateRHS(io, pkg, iIteration, lDirect)
    !! Generates the right-hand side vector.
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iRow
    real(kind=AE_REAL) :: rRHS
    complex(kind=AE_REAL), dimension(10) :: cCPZ
    integer(kind=AE_INT) :: iNCP, iEqType, iElementType, ic
    integer(kind=AE_INT) :: iElementString, iElementVertex, iElementFlag
    complex(kind=AE_REAL) :: cOrientation
    real(kind=AE_REAL) :: rGhbDistance

    do iRow = 1, pkg%aem%mat%iNEqn
      call MAT_GetEquation(io, pkg%aem%mat, iRow, cCPZ, iNCP, iEqType, iElementType, &
           iElementString, iElementVertex, iElementFlag, cOrientation, rGhbDistance)

      select case (iElementType)
        case (ELEM_AQU)
          rRHS = rAQU_ComputeRHS(io, pkg%aqu, pkg%aem%fdp, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_IN0)
          rRHS = rAQU_ComputeRHS(io, pkg%aqu, pkg%aem%fdp, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_LS1)
          rRHS = rLS1_ComputeRHS(io, pkg%ls1, pkg%aem, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_LS2)
          rRHS = rLS2_ComputeRHS(io, pkg%ls2, pkg%aem, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_LS3)
          rRHS = rLS3_ComputeRHS(io, pkg%ls3, pkg%aem, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_HB0)
          rRHS = rHB0_ComputeRHS(io, pkg%hb0, iEqType, iElementType, iElementString, &
                 iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_WL1)
          rRHS = rWL1_ComputeRHS(io, pkg%wl1, pkg%aqu, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
        case (ELEM_CW0)
          rRHS = rCW0_ComputeRHS(io, pkg%cw0, pkg%aem, iEqType, iElementType, &
                 iElementString, iElementVertex, iElementFlag, iIteration, lDirect)
      end select

      call MAT_SetRHS(io, pkg%aem%mat, iRow, rRHS)
    end do

    return
  end subroutine PKG_GenerateRHS


  subroutine PKG_SolveMatrix(io, pkg, rRelaxation, lDirect)
    !! Solves the matrix and distributes results to all element modules.
    type(PKG_DOMAIN), pointer :: pkg
    real(kind=AE_REAL), intent(in) :: rRelaxation
    logical :: lDirect
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: i
    real(kind=AE_REAL) :: rValue
    integer(kind=AE_INT) :: iElementType, iElementString, iElementVertex, iElementFlag

    call MAT_Solve(io, pkg%aem%mat)

    do i = 1, pkg%aem%mat%iNEqn
      call MAT_GetVariable(io, pkg%aem%mat, i, rValue, iElementType, &
           iElementString, iElementVertex, iElementFlag)

      select case (iElementType)
        case (ELEM_AQU)
          call AQU_StoreResult(io, pkg%aqu, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_IN0)
          call AQU_StoreResult(io, pkg%aqu, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS1)
          call LS1_StoreResult(io, pkg%ls1, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS2)
          call LS2_StoreResult(io, pkg%ls2, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_LS3)
          call LS3_StoreResult(io, pkg%ls3, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_HB0)
          call HB0_StoreResult(io, pkg%hb0, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_WL1)
          call WL1_StoreResult(io, pkg%wl1, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
        case (ELEM_CW0)
          call CW0_StoreResult(io, pkg%cw0, rRelaxation*rValue, &
               iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
      end select

    end do

    return
  end subroutine PKG_SolveMatrix


  function PKG_UpdateFreeSurface(io, pkg) result(iUpdate)
    !! Adjusts boundary condition elevations for free surface BDY elements.
    type(IO_STATUS), pointer :: io
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT) :: iUpdate
    integer(kind=AE_INT) :: ibdy, iNCP, iEqnType, iElementID, iElementString, iElementVertex, iElementFlag
    real(kind=AE_REAL) :: ph0, ph1, ph2, rzn, rGhbFactor, rMaxMove, rdz
    complex(kind=AE_REAL), dimension(2) :: cCPZ
    complex(kind=AE_REAL) :: cOrientation
    type(AQU_COLLECTION), pointer :: aqu
    type(AQU_BDYELEMENT), pointer :: this, next, prev

    iUpdate = 1

    if (.not. io%lProfile) return

    aqu => pkg%aqu
    aqu%lFSRegen = .false.
    aqu%BdyElements%rSpecHead = aqu%BdyElements%rSaveSpecHead
    aqu%BdyElements%lFSHeadSpec = .false.
    aqu%BdyElements%lMoveCPZ1 = .true.
    aqu%BdyElements%lMoveCPZ2 = .true.

    do ibdy = 1, aqu%iNBdy
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

    call IO_MessageText(io, "Scanning for head-specified conditions...")
    do ibdy = 1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      if (this%iBdyFlag == BDY_FREESURF) then
        call AQU_GetNeighborBdy(io, aqu, this, prev, next)
        if (this%rSpecHead > max(aimag(this%cFSZ1), aimag(this%cFSZ2))) then
          this%lFSHeadSpec = .true.
          this%cCPZ1 = this%cFSZ1
          this%lMoveCPZ1 = .false.
          prev%cCPZ2 = prev%cFSZ2
          prev%lMoveCPZ2 = .false.
          this%cCPZ2 = this%cFSZ2
          this%lMoveCPZ2 = .false.
          next%cCPZ1 = next%cFSZ1
          next%lMoveCPZ2 = .false.
        end if
      end if
    end do

    call IO_MessageText(io, "Moving points on free surface...")
    rMaxMove = rZERO
    do ibdy = 1, aqu%iNBdy
      this => aqu%BdyElements(ibdy)
      if (this%iBdyFlag == BDY_FREESURF) then
        call AQU_GetNeighborBdy(io, aqu, this, prev, next)
        if (this%lMoveCPZ1) then
          ph0 = rAEM_Head(io, pkg%aem, this%cCPZ1) - aimag(this%cCPZ1)
          if (ph0 > rZERO .and. abs(aimag(this%cFSZ1)-aimag(this%cCPZ1)) < 1.0e-4_AE_REAL) then
            print *, "Point 1 is on boundary with positive pressure head", ibdy
          else
            rZn = rAEM_FindPressureHeadZero(io, pkg%aem, real(this%cCPZ1), aimag(this%cCPZ1), &
                                             0.01_AE_REAL, 5, aimag(this%cFSZ1))
            rZn = 0.5_AE_Real*aimag(this%cCPZ1) + 0.5_AE_Real*rZn
            rdz = abs(rZn-aimag(this%cCPZ1))
            if (rdz > rMaxMove) rMaxMove = rdz
            this%cCPZ1 = cmplx(real(this%cFSZ1), rZN, AE_REAL)
            this%lMoveCPZ1 = .false.
          end if
        end if
        if (this%lMoveCPZ2) then
          ph0 = rAEM_Head(io, pkg%aem, this%cCPZ2) - aimag(this%cCPZ2)
          if (ph0 > rZERO .and. abs(aimag(this%cFSZ2)-aimag(this%cCPZ2)) < 1.0e-4_AE_REAL) then
            print *, "Point 2 is on boundary with positive pressure head", ibdy
          else
            rZn = rAEM_FindPressureHeadZero(io, pkg%aem, real(this%cCPZ2), aimag(this%cCPZ2), &
                                             0.01_AE_REAL, 5, aimag(this%cFSZ2))
            rZn = 0.5_AE_Real*aimag(this%cCPZ2) + 0.5_AE_Real*rZn
            rdz = abs(rZn-aimag(this%cCPZ2))
            if (rdz > rMaxMove) rMaxMove = rdz
            this%cCPZ2 = cmplx(real(this%cFSZ2), rZN, AE_REAL)
            this%lMoveCPZ2 = .false.
          end if
        end if
      end if
    end do
    write (IO_MessageBuffer, "('Largest vertical move on free surface: ',f10.4)") rMaxMove
    call IO_MessageText(io)

    do ibdy = 1, aqu%iNBdy
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
        end if
      end if
    end do

    aqu%lFSRegen = .true.
    iUpdate = 1

    return
  end function PKG_UpdateFreeSurface


  subroutine PKG_EnableDrawdown(io, pkg)
    !! Enables the drawdown elements and re-checks boundary conditions.
    type(IO_STATUS), pointer :: io
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT) :: iChanges

    call IO_MessageText(io, "Enabling the drawdown options")
    iChanges = 0
    iChanges = iChanges + WL0_EnableDrawdown(io, pkg%wl0)
    iChanges = iChanges + CW0_EnableDrawdown(io, pkg%cw0)
    write(unit=IO_MessageBuffer, fmt=*) "Number of changes: ", iChanges
    call IO_MessageText(io)

    if (iChanges > 0) then
      call PKG_Update(io, pkg)
      call PKG_ComputeCheck(io, pkg, .false.)
    end if

    return
  end subroutine PKG_EnableDrawdown


  subroutine PKG_ComputeCheck(io, pkg, lLinearize)
    !! Updates check information for all element packages.
    type(PKG_DOMAIN), pointer :: pkg
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io

    call AQU_ComputeCheck(io, pkg%aqu, pkg%aem, lLinearize)
    call WL0_ComputeCheck(io, pkg%wl0, pkg%aem, pkg%aqu)
    call WL1_ComputeCheck(io, pkg%wl1, pkg%aem, pkg%aqu)
    call LS0_ComputeCheck(io, pkg%ls0, pkg%aem)
    call LS1_ComputeCheck(io, pkg%ls1, pkg%aem)
    call LS2_ComputeCheck(io, pkg%ls2, pkg%aem, lLinearize)
    call LS3_ComputeCheck(io, pkg%ls3, pkg%aem, lLinearize)
    call HB0_ComputeCheck(io, pkg%hb0, pkg%aem)
    call CW0_ComputeCheck(io, pkg%cw0, pkg%aem)

    return
  end subroutine PKG_ComputeCheck


  subroutine PKG_Solve(io, pkg, iNIter, iNPolishIter, rRelaxation)
    !! Performs the full iterative solution process.
    type(PKG_DOMAIN), pointer :: pkg
    integer(kind=AE_INT), intent(inout) :: iNIter
    integer(kind=AE_INT), intent(in) :: iNPolishIter
    real(kind=AE_REAL), intent(in) :: rRelaxation
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iIter, i, iChanges, iUpdate

    if (pkg%aem%lSolutionPresent) then
      call IO_MessageText(io, "Solution is present -- continuing solve")
    else
      call IO_MessageText(io, "Allocating space for functions and matrix")
      call PKG_PreSolve(io, pkg)
      call PKG_AllocateFunctions(io, pkg)
      pkg%aem%lSolutionPresent = .true.
    end if

    iChanges = iWL1_Prepare(io, pkg%wl1, pkg%aqu, pkg%aem%iCurrentIteration)

    do iIter = 1, iNIter
      call IO_MessageText(io, '\n')
      pkg%aem%iCurrentIteration = pkg%aem%iCurrentIteration + 1
      write (unit=IO_MessageBuffer, fmt="(""Iteration: "", i5)") pkg%aem%iCurrentIteration
      call IO_MessageText(io)

      if (iIter > 1) then
        iUpdate = iAQU_Prepare(io, pkg%aqu, pkg%aem%iCurrentIteration)
        iUpdate = iUpdate + WL0_AdjustDischarges(io, pkg%wl0, pkg%aqu)
        iUpdate = iUpdate + PKG_UpdateFreeSurface(io, pkg)

        if (iIter > 1 .and. iIter < iNIter-iNPolishIter+1) then
          iUpdate = iUpdate + &
                    iWL1_Prepare(io, pkg%wl1, pkg%aqu, pkg%aem%iCurrentIteration) + &
                    iLS1_Prepare(io, pkg%ls1, pkg%aem%iCurrentIteration) + &
                    iLS2_DoRouting(io, pkg%ls2, pkg%aem%iCurrentIteration, .true., .false.) + &
                    iLS2_Prepare(io, pkg%ls2, pkg%aem, pkg%aem%iCurrentIteration) + &
                    iLS3_DoRouting(io, pkg%ls3, pkg%aem%iCurrentIteration, .true., .false.) + &
                    iLS3_Prepare(io, pkg%ls3, pkg%aem, pkg%aem%iCurrentIteration) + &
                    iHB0_Prepare(io, pkg%hb0, pkg%aem%iCurrentIteration)
          iUpdate = iUpdate + &
                    iCW0_Prepare(io, pkg%cw0, pkg%aem%iCurrentIteration)
        end if
        if (iUpdate > 0) then
          call PKG_Update(io, pkg)
          call PKG_ComputeCheck(io, pkg, .false.)
        end if
      end if

      if (iAQU_GetInfo(io, pkg%aqu, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iWL1_GetInfo(io, pkg%wl1, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iLS1_GetInfo(io, pkg%ls1, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iLS2_GetInfo(io, pkg%ls2, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iLS3_GetInfo(io, pkg%ls3, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iHB0_GetInfo(io, pkg%hb0, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          iCW0_GetInfo(io, pkg%cw0, INFO_REGENERATE, pkg%aem%iCurrentIteration) /= 0 .or. &
          .false.) then
        call IO_MessageText(io, "Allocating matrix")
        call MAT_Clear(io, pkg%aem%mat)
        call PKG_AllocateMatrix(io, pkg, pkg%aem%iCurrentIteration)
        call IO_MessageText(io, "Generating matrix")
        call PKG_GenerateMatrix(io, pkg, pkg%aem%iCurrentIteration)
        call IO_MessageText(io, "Decomposing matrix")
        call MAT_Decompose(io, pkg%aem%mat)
        print *, 'Condition Number: ', pkg%aem%mat%rCond
      end if

      if (pkg%aem%iCurrentIteration == 1) then
        call PKG_ComputeCheck(io, pkg, .true.)
      end if

      if (io%lDebug) then
        write (IO_MessageBuffer, '(i3)') pkg%aem%iCurrentIteration
      end if

      call IO_MessageText(io, "Generating solution...")
      call PKG_GenerateRHS(io, pkg, pkg%aem%iCurrentIteration, .false.)
      if (pkg%aem%iCurrentIteration == 1) then
        call PKG_SolveMatrix(io, pkg, rONE, .false.)
      else
        call PKG_SolveMatrix(io, pkg, rRelaxation, .false.)
      end if
      call PKG_Update(io, pkg)
      call PKG_ComputeCheck(io, pkg, .true.)
      call AQU_SetPrecondition(io, pkg%aqu, .false.)

    end do

    print *, 'Performing stream routing calculations'
    iChanges = iLS2_DoRouting(io, pkg%ls2, pkg%aem%iCurrentIteration, .false., .true.)
    iChanges = iLS3_DoRouting(io, pkg%ls3, pkg%aem%iCurrentIteration, .false., .true.)
    call WL0_SolvePartialPenetration(io, pkg%wl0, pkg%aqu)

    call IO_MessageText(io)
    call IO_MessageText(io, "Solution complete")

    return
  end subroutine PKG_Solve


  subroutine PKG_Read(io, pkg)
    !! Reads the input file and populates the PKG_DOMAIN object.

    type(PKG_DOMAIN), pointer :: pkg
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), parameter :: iEND = 1000, iDBG = 1001, &
                            iAQU = 2001, iWL0 = 2002, iPD0 = 2003, iWL1 = 2004, &
                            iLS0 = 2005, iLS1 = 2006, iLS2 = 2007, iLS3 = 2008, &
                            iHB0 = 2009, iAS0 = 20010, iOPT = 2011, &
                            iCW0 = 5001, &
                            iRLX = 3001
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
    character(len=255) :: sOptionText
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iretval
    integer(kind=AE_INT) :: iAS0Flag
    integer(kind=AE_INT) :: iNAS0
    integer(kind=AE_INT) :: iNInho, iNStr
    complex(kind=AE_REAL) :: cZ1, cZ2, cZC, cZE1, cZE2
    real(kind=AE_REAL) :: rR, rTol, rDPhiDX, rDPhiDY
    real(kind=AE_REAL) :: rBase, rThick, rHydCond, rPorosity, rAvgHead, rRelax
    logical :: lFlag

    call IO_MessageText(io, "Reading AEM module input")

    call IO_Assert(io, (associated(pkg)), &
         "PKG_Read: the PKG_DOMAIN has not been created")

    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          call IO_Assert(io, .false., "PKG_Read: I/O Error")
        case (kOpFileEOF)
          call IO_Assert(io, .false., "PKG_Read: Unexpected EOF")
        case (kOpData)
          call IO_Assert(io, .false., "PKG_Read: Unexpected data record")
        case (iEND)
          exit
        case (iDBG)
          io%lDebug = lIO_GetLogical(io, "lDebug", def=.true.)
        case (iAQU)
          iNInho = iIO_GetInteger(io, 'iNInho', minimum=1)
          iNStr = iIO_GetInteger(io, 'iNStr', minimum=0)
          rBase = rIO_GetReal(io, 'rBase')
          rThick = rIO_GetReal(io, 'rThick', minimum=rTINY)
          rHydCond = rIO_GetReal(io, 'rHydCond', minimum=rZERO)
          rPorosity = rIO_GetReal(io, 'rPorosity', minimum=rTINY)
          rAvgHead = rIO_GetReal(io, 'rAvgHead', minimum=rBase)
          pkg%aem%dom => DOM_Create(io, iNInho, rBase, rThick, rHydCond, rPorosity, rAvgHead)
          pkg%aem%frf => FRF_Create(io)
          pkg%aqu => AQU_Create(io, pkg%aem%dom, pkg%aem%frf, iNStr)
          call AQU_Read(io, pkg%aqu)
        case (iWL0)
          call WL0_Alloc(io, pkg%wl0)
          call WL0_Read(io, pkg%wl0)
        case (iWL1)
          call WL1_Alloc(io, pkg%wl1)
          call WL1_Read(io, pkg%wl1)
        case (iPD0)
          call PD0_Alloc(io, pkg%pd0)
          call PD0_Read(io, pkg%pd0)
        case (iLS0)
          call LS0_Alloc(io, pkg%ls0)
          call LS0_Read(io, pkg%ls0)
        case (iLS1)
          call LS1_Alloc(io, pkg%ls1)
          call LS1_Read(io, pkg%ls1)
        case (iLS2)
          call LS2_Alloc(io, pkg%ls2)
          call LS2_Read(io, pkg%ls2)
        case (iLS3)
          call LS3_Alloc(io, pkg%ls3)
          call LS3_Read(io, pkg%ls3)
        case (iHB0)
          call HB0_Alloc(io, pkg%hb0)
          call HB0_Read(io, pkg%hb0)
        case (iAS0)
          call IO_SwitchFields(io)
          if (iIO_GetInteger(io, "top/bottom", def=0, allowed=(/0, 1/)) == 0) then
            call AS0_Alloc(io, pkg%as0_top)
            call AS0_Read(io, pkg%as0_top)
          else
            call AS0_Alloc(io, pkg%as0_bottom)
            call AS0_Read(io, pkg%as0_bottom)
          end if
        case (iOPT)
          call IO_MessageText(io, '[OPT directive was ignored]')
        case (iCW0)
          call CW0_Alloc(io, pkg%cw0)
          call CW0_Read(io, pkg%cw0)
      end select
    end do

    call IO_MessageText(io, "Leaving AEM module")

    return
  end subroutine PKG_Read


  subroutine PKG_Save(io, pkg, sFname, mode)
    !! Saves the current solution to a file.
    type(PKG_DOMAIN), pointer :: pkg
    character(len=*), intent(in) :: sFname
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in) :: mode
    integer(kind=AE_INT) :: istat

    if (mode == IO_MODE_BINARY) then
      open(unit=LU_SCRATCH, file=trim(sFname), status="NEW", form="UNFORMATTED", iostat=istat)
    else
      open(unit=LU_SCRATCH, file=trim(sFname), status="NEW", iostat=istat)
    end if
    call IO_Assert(io, istat==0, "Open failed on output file " // trim(sFname))

    call AQU_Save(io, pkg%aqu, mode)
    call LS1_Save(io, pkg%ls1, mode)
    call LS2_Save(io, pkg%ls2, mode)
    call LS3_Save(io, pkg%ls3, mode)
    call HB0_Save(io, pkg%hb0, mode)
    call WL1_Save(io, pkg%wl1, mode)
    call CW0_Save(io, pkg%cw0, mode)

    close(unit=LU_SCRATCH)

    return
  end subroutine PKG_Save


  subroutine PKG_Load(io, pkg, sFname, mode)
    !! Loads a previously-saved preconditioning file from disk.
    type(IO_STATUS), pointer :: io
    type(PKG_DOMAIN), pointer :: pkg
    character(len=*), intent(in) :: sFname
    integer(kind=AE_INT), intent(in) :: mode
    integer(kind=AE_INT) :: istat

    if (mode == IO_MODE_BINARY) then
      open(unit=LU_SCRATCH, file=trim(sFname), status="OLD", form="UNFORMATTED", iostat=istat)
    else
      open(unit=LU_SCRATCH, file=trim(sFname), status="OLD", iostat=istat)
    end if
    call IO_Assert(io, istat==0, "Open failed on input file " // trim(sFname))

    call AQU_Load(io, pkg%aqu, pkg%aem%fdp, mode)
    call LS1_Load(io, pkg%ls1, pkg%aem%fls, mode)
    call LS2_Load(io, pkg%ls2, pkg%aem%fls, mode)
    call LS3_Load(io, pkg%ls3, pkg%aem%fls, mode)
    call HB0_Load(io, pkg%hb0, pkg%aem%fdp, mode)
    call WL1_Load(io, pkg%wl1, pkg%aem%fwl, mode)
    call CW0_Load(io, pkg%cw0, pkg%aem%fls, mode)

    close(unit=LU_SCRATCH)

    return
  end subroutine PKG_Load


end module p_packages
