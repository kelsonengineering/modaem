module u_matrix

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

  ! This is a generic module which supports the generation and solution of
  ! a full MxM matrix.

  use u_constants
  use u_io

  implicit none

  public


  type, public :: MAT_EQUATION
    !! type MAT_EQUATION
    !!
    !! PUBLIC type that holds information for one equation in the matrix
    !!
    !! Members:
    !!   complex :: cCPZ(:)
    !!     The center of the well
    !!
    complex(kind=AE_REAL), dimension(:), pointer :: cCPZ
    integer(kind=AE_INT) :: iEqnType
    integer(kind=AE_INT) :: iElementID
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
    complex(kind=AE_REAL) :: cOrientation
    real(kind=AE_REAL) :: rGhbFactor
  end type MAT_EQUATION


  type, public :: MAT_VARIABLE
    !! type MAT_VARIABLE
    !!
    !! PUBLIC type that holds information for one unknown variable
    !!
    !! Members:
    !!   integer :: iElementID
    !!     The element ID(e.g. WL0, LS0, LS1, HB0)
    !!
    integer(kind=AE_INT) :: iElementID
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
    real(kind=AE_REAL) :: rValue
  end type MAT_VARIABLE


  type, public :: MAT_MATRIX
    !! type MAT_MATRIX
    !!
    !! PUBLIC type that holds information for a matrix
    !!
    !! Members:
    !!   type(MAT_EQUATION), pointer :: Equations(:)
    !!     The equation definitions
    !!
    type(MAT_EQUATION), dimension(:), pointer :: Equations
    type(MAT_VARIABLE), dimension(:), pointer :: Variables
    real(kind=AE_REAL) :: rCond
    real(kind=AE_REAL), dimension(:, :), pointer :: rMatrix
#ifdef __MATRIX_TEST__
    real(kind=AE_REAL), dimension(:, :), pointer :: rOriginalMatrix
#endif
    real(kind=AE_REAL), dimension(:, :), pointer :: rRHS
    integer(kind=AE_INT), dimension(:), pointer :: iPivot
    integer(kind=AE_INT) :: iNEqn
    integer(kind=AE_INT) :: iNVar
  end type MAT_MATRIX

  ! Iterative refinement settings for gesolve()
  integer(kind=AE_INT), private :: iMATMaxIterations = 10
  real(kind=AE_REAL), private :: rMATIterTolX = 1.0e-6_AE_REAL
  real(kind=AE_REAL), private :: rMATIterTolR = 1.0e-6_AE_REAL
  real(kind=AE_REAL), private :: rMOVEPOINT = 1.0e-3_AE_REAL


contains


  function MAT_Create(io) result(mat)
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(MAT_MATRIX), pointer :: mat
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iL, iStat

    allocate(mat, stat = iStat)
    call IO_Assert(io, (iStat == 0), "MAT_Create: Allocation failed")

    mat%iNEqn = 0
    mat%iNVar = 0
    nullify(mat%Equations, &
                             mat%Variables, &
                             mat%rMatrix, &
#ifdef __MATRIX_TEST__
                             mat%rOriginalMatrix, &
#endif
                             mat%rRHS, &
                             mat%iPivot &
                             )
    return
  end function MAT_Create


  subroutine MAT_Alloc(io, mat, iNEQ, iNUN)
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iNEQ
    integer(kind=AE_INT), intent(in) :: iNUN
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_Alloc: MAT_Create has not been called")
      call IO_Assert(io, (iNEQ > 0 .and. iNUN > 0), &
           "MAT_Alloc: Illegal dimensions")
      call IO_Assert(io, (iNEQ == iNUN), &
           "MAT_Alloc: Overspecification not supported in this version")
    end if

    allocate(mat%rMatrix(iNEQ, iNUN), &
#ifdef __MATRIX_TEST__
                                 mat%rOriginalMatrix(iNEQ, iNUN), &
#endif
                                 mat%rRHS(iNEQ, 1), &
                                 mat%Equations(iNEQ), &
                                 mat%Variables(iNUN), &
                                 mat%iPivot(iNEQ), stat = iStat &
                                 )
    call IO_Assert(io, (iStat == 0), "MAT_Alloc: Allocation failed")

    mat%rMatrix = rZERO
    mat%rRHS = rZERO
    mat%iNEqn = 0
    mat%iNVar = 0

    return
  end subroutine MAT_Alloc


  subroutine MAT_Destroy(io, mat)
    !! This subroutine de-allocates the space prevously allocated for the matrix
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    if (associated(mat%Equations)) deallocate(mat%Equations, stat = iStat)
    if (associated(mat%Variables)) deallocate(mat%Variables, stat = iStat)
    if (associated(mat%rMatrix)) deallocate(mat%rMatrix, stat = iStat)
    if (associated(mat%rRHS)) deallocate(mat%rRHS, stat = iStat)
    if (associated(mat%iPivot)) deallocate(mat%iPivot, stat = iStat)
#ifdef __MATRIX_TEST__
    if (associated(mat%rOriginalMatrix)) deallocate(mat%rOriginalMatrix, stat = iStat)
#endif
    deallocate(mat, stat = iStat)
    call IO_Assert(io, (iStat == 0), "MAT_Destroy: Deallocation failed")

    return
  end subroutine MAT_Destroy


  subroutine MAT_Clear(io, mat)
    !! This subroutine deallocates all the matrix data, leaving the matrix object intact
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iEqn

    iStat = 0
    if (.not. associated(mat)) return

    if (associated(mat%Equations)) then
      do iEqn = 1, size(mat%Equations)
        if (associated(mat%Equations(iEqn)%cCPZ)) deallocate(mat%Equations(iEqn)%cCPZ)
      end do
      deallocate(mat%Equations, stat = iStat)
    end if
    if (associated(mat%Variables)) deallocate(mat%Variables, stat = iStat)
    if (associated(mat%rMatrix)) deallocate(mat%rMatrix, stat = iStat)
    if (associated(mat%rRHS)) deallocate(mat%rRHS, stat = iStat)
    if (associated(mat%iPivot)) deallocate(mat%iPivot, stat = iStat)
#ifdef __MATRIX_TEST__
    if (associated(mat%rOriginalMatrix)) deallocate(mat%rOriginalMatrix, stat = iStat)
#endif
    call IO_Assert(io, (iStat == 0), 'MAT_Clear: Deallocation failed')

    return
  end subroutine MAT_Clear


  subroutine MAT_CreateVariable(io, mat, iElementID, iElementString, iElementVertex, iElementFlag)
    ! This function is called by the xxxSetup routines to allocate space for an unknown variable
    ! (e.g. an unknown strength coefficient).  This information is stored in the 'variables' buffer.
    ! The function returns the index into the 'variables' buffer for the unknown value
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iElementID
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(MAT_VARIABLE), pointer :: vbl

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_Alloc: MAT_Create has not been called")
      call IO_Assert(io, (mat%iNVar <= size(mat%Variables)), &
           "MAT_CreateVariable: Space exhausted")
    end if

    mat%iNVar = mat%iNVar+1
    vbl => mat%Variables(mat%iNVar)
    vbl%iElementID = iElementID
    vbl%iElementString = iElementString
    vbl%iElementVertex = iElementVertex
    vbl%iElementFlag = iElementFlag
    vbl%rValue = rZERO

    return
  end subroutine MAT_CreateVariable


  subroutine MAT_GetVariable(io, mat, iVar, rValue, iElementID, iElementString, &
               iElementVertex, iElementFlag)
    ! This function extracts the information about the specified variable
    ! It returns kOK if no error was detected, else an error code
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iVar
    real(kind=AE_REAL), intent(out) :: rValue
    integer(kind=AE_INT), intent(out) :: iElementID
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sBuf
    type(MAT_VARIABLE), pointer :: var

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_Alloc: MAT_Create has not been called")
      call IO_Assert(io, (iVar >= 1 .and. iVar <= mat%iNVar), &
           "MAT_GetVariable: Bad variable number")
    end if

    var => mat%Variables(iVar)
    rValue = var%rValue
    iElementID = var%iElementID
    iElementString = var%iElementString
    iElementVertex = var%iElementVertex
    iElementFlag = var%iElementFlag

    return
  end subroutine MAT_GetVariable


  function MAT_CreateEquation(io, mat, cCPZ, iEqnType, iElementID, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor) result(iRV)
    ! This function is called by the xxxSetup routines to create entries in the
    ! matrix generator. It returns the index of the new equation.
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    complex(kind=AE_REAL), dimension(:), intent(in) :: cCPZ
    integer(kind=AE_INT), intent(in) :: iEqnType
    integer(kind=AE_INT), intent(in) :: iElementID
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    complex(kind=AE_REAL), intent(in) :: cOrientation
    real(kind=AE_REAL), intent(in) :: rGhbFactor
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iRV
    ! [ LOCALS ]
    integer(kind=AE_INT) :: ieq, ipt
    real(kind=AE_REAL), parameter :: rCPTolerance = 1.0e-4_AE_REAL
    logical :: lFail
    type(MAT_EQUATION), pointer :: eqn

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_CreateEquation: MAT_Create has not been called")
      call IO_Assert(io, (mat%iNEqn < size(mat%Equations)), &
           "MAT_CreateEquation: Space exhausted")
    end if

    select case (iEqnType)
      case (EQN_FLOW)
        call IO_Assert(io, (size(cCPZ) > 1), &
             "MAT_CreateEquation: Path too short for EQN_FLOW")
      case (EQN_BDYGHB)
        call IO_Assert(io, (size(cCPZ) > 1), &
             "MAT_CreateEquation: Path too short for EQN_BDYGHB")
      case (EQN_POTENTIALDIFF)
        call IO_Assert(io, (size(cCPZ) == 2), &
             "MAT_CreateEquation: Exactly two control points required")
      case default
        call IO_Assert(io, (size(cCPZ) == 1), &
             "MAT_CreateEquation: Exactly one control point required")
    end select

    ! Check for overlapping control points
    do ieq = 1, mat%iNEqn
      eqn => mat%Equations(ieq)
      if ((eqn%iEqnType == iEqnType) .and. &
          eqn%iEqnType /= EQN_TOTALFLOW .and. &
          (size(cCPZ, 1) == size(eqn%cCPZ, 1))) then
        call IO_Assert(io, (maxval(abs(cCPZ-eqn%cCPZ)) > rCPTolerance), &
             "MAT_CreateEquation: Coincident control points")
      end if
    end do

    mat%iNEqn = mat%iNEqn+1
    eqn => mat%Equations(mat%iNEqn)
    allocate(eqn%cCPZ(ubound(cCPZ, 1)))
    eqn%cCPZ = cCPZ
    eqn%iEqnType = iEqnType
    eqn%iElementID = iElementID
    eqn%iElementString = iElementString
    eqn%iElementVertex = iElementVertex
    eqn%iElementFlag = iElementFlag
    eqn%cOrientation = cOrientation
    eqn%rGhbFactor = rGhbFactor
    iRV = mat%iNEqn

    return
  end function MAT_CreateEquation


  subroutine MAT_GetEquation(io, mat, iEQ, cCPZ, iNCP, iEqnType, iElementID, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor)
    ! This function extracts the information about the specified equation
    ! It returns kOK if no error was detected, else an error code
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iEQ
    complex(kind=AE_REAL), dimension(:), intent(out) :: cCPZ
    integer(kind=AE_INT), intent(out) :: iNCP
    integer(kind=AE_INT), intent(out) :: iEqnType
    integer(kind=AE_INT), intent(out) :: iElementID
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    complex(kind=AE_REAL), intent(out) :: cOrientation
    real(kind=AE_REAL), intent(out) :: rGhbFactor
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sBuf
    type(MAT_EQUATION), pointer :: eqn

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_GetEquation: MAT_Create has not been called")
      call IO_Assert(io, (iEQ <= mat%iNEqn), &
           "MAT_GetEquation: Bad equation number")
    end if

    eqn => mat%Equations(iEQ)
    iNCP = ubound(eqn%cCPZ, 1)
    cCPZ(1:iNCP) = eqn%cCPZ
    iEqnType = eqn%iEqnType
    iElementID = eqn%iElementID
    iElementString = eqn%iElementString
    iElementVertex = eqn%iElementVertex
    iElementFlag = eqn%iElementFlag
    cOrientation = eqn%cOrientation
    rGhbFactor = eqn%rGhbFactor

    return
  end subroutine MAT_GetEquation


  subroutine MAT_SetEquation(io, mat, iEQ, cCPZ, iEqnType, iElementID, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor)
    ! This function extracts the information about the specified equation
    ! It returns kOK if no error was detected, else an error code
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iEQ
    complex(kind=AE_REAL), dimension(:), intent(in) :: cCPZ
    integer(kind=AE_INT), intent(in) :: iEqnType
    integer(kind=AE_INT), intent(in) :: iElementID
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    complex(kind=AE_REAL), intent(in) :: cOrientation
    real(kind=AE_REAL), intent(in) :: rGhbFactor
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sBuf
    type(MAT_EQUATION), pointer :: eqn

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_GetEquation: MAT_Create has not been called")
      call IO_Assert(io, (iEQ <= mat%iNEqn), &
           "MAT_GetEquation: Bad equation number")
    end if

    eqn => mat%Equations(iEQ)
    ! Check to see if there's a change in the size of the CP array...
    if (ubound(eqn%cCPZ, 1) /= ubound(cCPZ, 1)) then
      deallocate(eqn%cCPZ)
      allocate(eqn%cCPZ(ubound(cCPZ,1)))
    end if
    eqn%cCPZ = cCPZ
    eqn%iEqnType = iEqnType
    eqn%iElementID = iElementID
    eqn%iElementString = iElementString
    eqn%iElementVertex = iElementVertex
    eqn%iElementFlag = iElementFlag
    eqn%cOrientation = cOrientation
    eqn%rGhbFactor = rGhbFactor

    return
  end subroutine MAT_SetEquation


  subroutine MAT_UpdateEquation(io, mat, iEQ)
    ! Updates the Check field in the specified equation
    ! It returns kOK if no error was detected, else an error code
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iEQ
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    character(len=255) :: sBuf
    type(MAT_EQUATION), pointer :: eqn

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_UpdateEquation: MAT_Create has not been called")
      call IO_Assert(io, (iEQ <= mat%iNEqn), &
           "MAT_UpdateEquation: Bad equation number")
    end if

    eqn => mat%Equations(iEQ)

    return
  end subroutine MAT_UpdateEquation


  subroutine MAT_SetRow(io, mat, iEQ, rRow, rRHS)
    ! Stores the given row data and right-hand side value into the matrix
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iEQ
    real(kind=AE_REAL), dimension(:), intent(in) :: rRow
    real(kind=AE_REAL), intent(in) :: rRHS
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_SetRow: MAT_Create has not been called")
      call IO_Assert(io, (iEQ <= mat%iNEqn), &
           "MAT_SetRow: Bad equation number")
      call IO_Assert(io, (size(rRow) == mat%iNVar), &
           "MAT_SetRow: Wrong number of terms in matrix row")
    end if

    ! All's well! Store the row
    mat%rMatrix(iEQ, :) = rRow(1:size(mat%rMatrix, 2))
#ifdef __MATRIX_TEST__
    mat%rOriginalMatrix(iEQ, :) = rRow(1:size(mat%rMatrix, 2))
#endif
    mat%rRHS(iEQ, 1) = rRHS

    return
  end subroutine MAT_SetRow


  subroutine MAT_SetRHS(io, mat, iEQ, rRHS)
    ! Stores the given row data and right-hand side value into the matrix
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT), intent(in) :: iEQ
    real(kind=AE_REAL), intent(in) :: rRHS
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_SetRHS: MAT_Create has not been called")
      call IO_Assert(io, (iEQ <= mat%iNEqn), &
           "MAT_SetRHS: Bad equation number")
    end if

    mat%rRHS(iEQ, 1) = rRHS

    return
  end subroutine MAT_SetRHS


  subroutine MAT_Decompose(io, mat)
    ! Solves the matrix
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rInfNorm, rCond
    real(kind=AE_REAL), dimension(:), allocatable :: rSums, rWork
    integer(kind=AE_INT) :: iInfo, i
    integer(kind=AE_INT), dimension(:), allocatable :: iWork
    integer(kind=AE_INT) :: iRV

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_Decompose: MAT_Create has not been called")
      call IO_Assert(io, (mat%iNEqn > 0 .and. mat%iNVar > 0), &
           "MAT_Decompose: No rows/columns exist")
      call IO_Assert(io, (mat%iNVar == mat%iNEqn), &
           "MAT_Decompose: Matrix is not square")
    end if

    ! Compute the 1-Norm of the matrix
    allocate(rSums(size(mat%rMatrix, 1)), &
                                       rWork(4*size(mat%rMatrix, 1)), &
                                       iWork(size(mat%rMatrix, 1)))
    rInfNorm = maxval((/(mat%rMatrix(i, 1), i = 1, size(mat%rMatrix, 1))/))
    ! Decompose the matrix
    call DGETRF(size(mat%rMatrix, 1), &
         size(mat%rMatrix, 2), &
         mat%rMatrix, &
         size(mat%rMatrix, 1), &
         mat%iPivot, &
         iRV)
    if (iRV /= 0) then
      write (unit=IO_MessageBuffer, fmt=*) "LAPACK solver DGETRF failed code = ", iRV
      call IO_Assert(io, .false., "")
    end if

    ! Compute the matrix condition number for testing
    call DGECON('I', &
         size(mat%rMatrix, 1), &
         mat%rMatrix, &
         size(mat%rMatrix, 2), &
         rInfNorm, &
         rCond, &
         rWork, &
         iWork, &
         iInfo)
    mat%rCond = rCond
    deallocate(rSums, rWork, iWork)

    return
  end subroutine MAT_Decompose


  subroutine MAT_Solve(io, mat)
    ! Solves the matrix
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    real(kind=AE_REAL), dimension(:), allocatable :: rX
    real(kind=AE_REAL), dimension(:, :), allocatable :: rSaveRHS
    integer(kind=AE_INT) :: iRV
    type(MAT_VARIABLE), pointer :: vbl

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_SetRHS: MAT_Create has not been called")
      call IO_Assert(io, (mat%iNEqn > 0 .and. mat%iNVar > 0), &
           "MAT_SetRHS: No rows/columns exist")
      call IO_Assert(io, (mat%iNVar == mat%iNEqn), &
           "MAT_SetRHS: Matrix is not square")
    end if

    allocate(rX(ubound(mat%rRHS, 1)), rSaveRHS(ubound(mat%rRHS, 1), 1))
    rSaveRHS = mat%rRHS
    call DGETRS('N', &
         size(mat%rMatrix, 1), &
         1, &
         mat%rMatrix, &
         size(mat%rMatrix, 1), &
         mat%iPivot, &
         mat%rRHS, &
         size(mat%rRHS, 1), &
         iRV)
    call IO_Assert(io, (iRV == 0), &
         "LAPACK solver DGETRS failed")
    print *, "Maximum change in solver: ", maxval(abs(mat%rRHS))
#ifdef __MATRIX_TEST__
    print *, 'Maximum solver error: ', maxval((matmul(mat%rMatrix, mat%rRHS)-rSaveRHS)/(matmul(mat%rMatrix, mat%rRHS)+rSaveRHS))-rOne
#endif

    ! Store the x-vector into the variables
    do i = 1, mat%iNVar
      vbl => mat%Variables(i)
      vbl%rValue = mat%rRHS(i, 1)
    end do

#ifdef __MATRIX_TEST__
    ! Print out the residuals...
    print *, 'MATRIX RESIDUALS'
    print *, '  RHS Max', maxval(rSaveRHS)
    print *, '  MIN:   ', minval(abs(matmul(mat%rOriginalMatrix, mat%rRHS)-rSaveRHS)), minloc(abs(matmul(mat%rOriginalMatrix, mat%rRHS)-rSaveRHS))
    print *, '  MAX:   ', maxval(abs(matmul(mat%rOriginalMatrix, mat%rRHS)-rSaveRHS)), maxloc(abs(matmul(mat%rOriginalMatrix, mat%rRHS)-rSaveRHS))
#endif

    deallocate(rX, rSaveRHS)
    return
  end subroutine MAT_Solve


  subroutine MAT_ComputeControlPoints(io, cZ1, cZ2, iNCP, cCPResult, rNormalOffset, rLocations)
    ! This subroutine computes a set of iNCP control points for the linear
    ! element which extends from cZ1 to cZ2.  Points are returned as follows:
    ! cCPResult(1) = cZ1, cCPResult(2:iNCP+1) are computed, cCPResult(iNCP+2) = cCZ2
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    integer(kind=AE_INT), intent(in) :: iNCP
    complex(kind=AE_REAL), dimension(:), intent(out) :: cCPResult
    real(kind=AE_REAL), intent(in) :: rNormalOffset
    type(IO_STATUS), pointer :: io
    real(kind=AE_REAL), dimension(:), intent(in), optional :: rLocations
    ! Locals
    complex(kind=AE_REAL) :: cZC, cZL, cOffset
    real(kind=AE_REAL) :: rDTheta
    integer(kind=AE_INT) :: i

    ! Check for coincident control points
    call IO_Assert(io, (abs(cZ2-cZ1) > rMOVEPOINT), &
         "MAT_ComputeControlPoints: Endpoints are coincident")
    call IO_Assert(io, (size(cCPResult) >= iNCP+2), &
         "MAT_ComputeControlPoints: Insufficient space for results")

    ! The control points are computed by turning an angle about the center
    ! of the segment.
    cZC = rHALF * (cZ1+cZ2)
    cZL = rHALF * (cZ2-cZ1)
    rDTheta = rPI / (iNCP+1)
    cOffset = cI*(cZ2-cZ1) * rNormalOffset

    if (.not. present(rLocations)) then
      ! Here we go
      cCPResult(1) = cZ1 + rTWO*rMOVEPOINT*cZL + cOffset
      do i = 1, iNCP
        cCPResult(i+1) = cZC - cos(i*rDTheta) * cZL + cOffset
      end do
      cCPResult(iNCP+2) = cZ2 - rTWO*rMOVEPOINT*cZL + cOffset
    else
      do i = 1, size(rLocations)
        cCPResult(i) = cZ1 + rLocations(i)*(cZ2-cZ1) + cOffset
      end do
    end if

    return
  end subroutine MAT_ComputeControlPoints


  subroutine MAT_Report(io, mat, label, lgridonly)
    ! This routine writes a report of all matrix information
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    character(len=*), intent(in) :: label
    logical :: lgridonly
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j, iStat, icp
    type(MAT_VARIABLE), pointer :: vbl
    type(MAT_EQUATION), pointer :: eqn

    if (io%lDebug) then
      call IO_Assert(io, (associated(mat)), &
           "MAT_SetRHS: MAT_Create has not been called")
    end if

    if (.not. lgridonly) then
      call HTML_Header('Utility module MAT', 1)
      call HTML_Header('Information about the matrix generator', 2)
      if (associated(mat%rMatrix)) then
        call HTML_StartTable()
        call HTML_AttrInteger('Number of equations', mat%iNEqn)
        call HTML_AttrInteger('Number of unknowns', mat%iNVar)
        call HTML_AttrReal('Condition number', mat%rCond)
        call HTML_EndTable()

        if (mat%iNEqn > 0) then
          call HTML_Header('Matrix equation generator data', 3)
          do i = 1, mat%iNEqn
            eqn => mat%Equations(i)
            call HTML_Header('Equation information', 4)
            call HTML_StartTable()
            call HTML_AttrInteger('Equation', i)
            call HTML_AttrInteger('BC type', eqn%iEqnType)
            call HTML_AttrInteger('Element type', eqn%iElementID)
            call HTML_AttrInteger('String', eqn%iElementString)
            call HTML_AttrInteger('Vertex', eqn%iElementVertex)
            call HTML_AttrInteger('Flag', eqn%iElementFlag)
            call HTML_EndTable()
            if (size(eqn%cCPZ) > 0) then
              call HTML_Header('Control points', 4)
              call HTML_StartTable()
              call HTML_TableHeader((/' ', 'X', 'Y'/))
              do icp = 1, size(eqn%cCPZ)
                call HTML_StartRow()
                call HTML_ColumnInteger((/i/))
                call HTML_ColumnComplex((/eqn%cCPZ(icp)/))
                call HTML_EndRow()
              end do
              call HTML_EndTable()
            else
              call HTML_Header('No control points', 4)
            end if
          end do
        else
          call HTML_Header('No equations defined', 3)
        end if

        if (mat%iNVar > 0) then
          call HTML_Header('Matrix unknown variable data', 3)
          call HTML_StartTable()
          call HTML_TableHeader((/'Unk   ', 'Elem  ', 'String', 'Vertex', 'Flag  '/))
          do i = 1, mat%iNVar
            vbl => mat%Variables(i)
            call HTML_StartRow()
            call HTML_ColumnInteger((/i, vbl%iElementID, vbl%iElementString, vbl%iElementVertex, vbl%iElementFlag/))
            call HTML_EndRow()
          end do
          call HTML_EndTable()
        end if
        call HTML_Header('Matrix values were written to mat_report.'//label//'.m', 3)

        open(unit=LU_GRID, file = "mat_report."//label//".m", iostat=iStat)
        call IO_Assert(io, (iStat == 0), "MAT_Report: Could not open grid file mat_report."//label//".m")
        write (unit=LU_GRID, fmt="('# Written by ModAEM 1.4')")
        write (unit=LU_GRID, fmt="('# name: matrix')")
        write (unit=LU_GRID, fmt="('# type: matrix')")
        write (unit=LU_GRID, fmt="('# rows: ', i10)") size(mat%rMatrix, 1)
        write (unit=LU_GRID, fmt="('# columns: ', i10)") size(mat%rMatrix, 2)
        do i = 1, mat%iNEqn
          write (unit=LU_GRID, fmt=*) (/(mat%rMatrix(i, j), j = 1, size(mat%rMatrix, 2))/)
        end do
        close(unit=LU_GRID)
      end if

    else
      call HTML_Header('Matrix is not allocated', 2)
    end if

    return
  end subroutine MAT_Report


  subroutine MAT_ControlPointReport(io, mat, iLU)
    ! This routine writes a report of control-point information
    ! [ ARGUMENTS ]
    type(MAT_MATRIX), pointer :: mat
    integer(kind=AE_INT) :: iLU
    type(IO_STATUS), pointer :: io



    return
  end subroutine MAT_ControlPointReport

end module u_matrix
