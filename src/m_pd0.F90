module m_pd0

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

  !! module m_pd0
  !!
  !! Element module for 2-D discharge specified ponds
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_pond     --  Function module for collections of ponds
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! pond elements. Each element requires a location, radius and discharge.

  use u_constants
  use u_io
  use f_pond

  implicit none

  public

  type :: PD0_POND
    !! type PD0_POND
    !!
    !! Type that holds information for one pond
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the pond
    !!   real :: rRadius
    !!     The radius of the pond
    !!   complex :: rGamma
    !!     Discharge of the pond; positive indicates extraction
    !!     and negative indicates injection.
    !!   integer :: iID
    !!     Identification label for the pond(for interaction with e.g. GUIs)
    !!   integer :: iFPDIndex
    !!     Index for the pond entry in the FPD module
    !!
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rGamma
    real(kind=AE_REAL) :: rRadius
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iFPDIndex
  end type PD0_POND

  type :: PD0_COLLECTION
    !! type PD0_COLLECTION
    !!
    !! Type that holds all ponds in a layer
    !!
    !! Members:  type(IO_STATUS), pointer :: io


    !!   type(PD0_POND), dimension(:), pointer :: Ponds
    !!     Array of PD0_POND objects for the layer; dimensioned for the maximum
    !!     number of ponds according to the input file(see PD0_Read)
    !!   integer :: iCount
    !!     The actual number of ponds in use in the layer
    !!
    type(PD0_POND), dimension(:), pointer :: Ponds
    integer(kind=AE_INT) :: iCount
  end type PD0_COLLECTION


contains

  !**pd Modified to allocate only the collection object

  function PD0_Create(io) result(pd0)
    !! function PD0_Create
    !!
    !! Creates a new PD0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    pd0 => PD0_Create()
    !!
    !! Arguments:
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(PD0_COLLECTION), pointer :: pd0
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(pd0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "PD0_Create: allocation failed")
    nullify(pd0%Ponds)
    pd0%iCount = 0

    return
  end function PD0_Create


  subroutine PD0_Alloc(io, pd0)
    !! Subroutine PD0_Alloc
    !!
    !! Allocates wells for the PD0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call PD0_Alloc(io, pd0, iNPD)
    !!
    !! Arguments:
    !!    (in)    type(PD0_COLLECTION), pointer :: pd0
    !!              The PD0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iNPD
    integer(kind=AE_INT) :: iStat

    iNPD = iIO_GetInteger(io, 'iNPD', minimum = 0)
    allocate(pd0%Ponds(iNPD), stat = iStat)
    call IO_Assert(io, (iStat == 0), "PD0_Alloc: allocation failed")

    return
  end subroutine PD0_Alloc


  subroutine PD0_Destroy(io, pd0)
    !! subroutine PD0_Destroy
    !!
    !! Frees memory allocated for an PD0 Ponds and PD0 Collection object
    !!
    !! Calling Sequence:
    !!     call PD0_Destroy(pd0)
    !!
    !! Arguments:
    !!  type(PD0_COLLECTION), pointer :: pd0
    !!              Pointer to the pd0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat


    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_Destroy: PD0_Create has not been called")
    end if

    if (associated(pd0%Ponds)) then
      deallocate(pd0%Ponds, stat = iStat)
      call IO_Assert(io, (iStat == 0), &
           "PD0_Destroy: deallocation of Ponds failed")
    end if
    deallocate(pd0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "PD0_Destroy: deallocation failed")

    return
  end subroutine PD0_Destroy


  subroutine PD0_New(io, pd0, Pond)
    !! function PD0_New
    !!
    !! Adds a new PD0_POND object to the PD0_COLLECTION 'pd0'
    !!
    !! Calling Sequence:
    !!    call PD0_New(io, pd0, Pond)
    !!
    !! Arguments:
    !!    (in)    type(PD0_COLLECTION), pointer :: pd0
    !!              The PD0_COLLECTION object to be used
    !!    (in)    type(PD0_POND), pointer :: Pond
    !!              Vector that defines the points along the barrier
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(PD0_POND) :: Pond
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_New: PD0_Create has not been called")
    end if

    call IO_Assert(io, (pd0%iCount < size(pd0%Ponds)), &
         "PD0_New: Space exhausted")

    pd0%iCount = pd0%iCount + 1
    pd0%Ponds(pd0%iCount) = Pond

    return
  end subroutine PD0_New


  subroutine PD0_PreSolve(io, pd0)
    !! subroutine PD0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call PD0_GetInfo(pd0)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(IO_STATUS), pointer :: io

    return
  end subroutine PD0_PreSolve


  function iPD0_GetInfo(io, pd0, iOption, iIteration) result(iValue)
    !! function PD0_GetInfo
    !!
    !! Returns the following sizing requirements for the PD0module
    !!
    !! Calling Sequence:
    !!    iValue = iPD0_GetInfo(io, pd0, iOption)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_GetInfo: PD0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FPD)
        iValue = pd0%iCount
      case default
        iValue = 0
    end select

    return
  end function iPD0_GetInfo


  subroutine PD0_SetupFunctions(io, pd0, fpd)
    !! subroutine PD0_Setup
    !!
    !! This routine sets up the functions in f_pond for the pond elements
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_pond by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call PD0_Setup(pd0)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             FPD_COLLECTION to be used
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(FPD_COLLECTION), pointer :: fpd
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(PD0_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_Setup: PD0_Create has not been called")
      call IO_Assert(io, (associated(fpd)), &
           "PD0_Setup: Illegal FPD_COLLECTION object")
    end if

    do i = 1, pd0%iCount
      ! Create a pond function in FPD for each pond
      pnd => pd0%Ponds(i)
      call FPD_New(io, fpd, pnd%cZ, pnd%rGamma, pnd%rRadius, ELEM_PD0, i, -1, -1, pnd%iFPDIndex)
    end do

    return
  end subroutine PD0_SetupFunctions


  subroutine PD0_ComputeRHS(io, pd0, fpd, cCPZ, iEqType, cOrientation, rRHS)
    !! subroutine PD0_ComputeRHS
    !!
    !! Computes the contribution of the PD0 elements to the right-hand-side
    !! value for the specified equation parameters.
    !!
    !! Calling Sequence:
    !!    PD0_ComputeRHS(io, pd0, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             FPD_COLLECTION to be used
    !!   (in)    complex, dimension(:) :: cCPZ
    !!             Control point(s) to be used in coefficient calculations
    !!   (in)    integer :: iEqType
    !!             The type of equation to be used
    !!   (in)    complex :: cOrientation
    !!             Orientation unit vector(for discharge-based equations)
    !!   (inout) real :: rRHS
    !!             On return, rRHS is updated by adding the contribution of
    !!             this module to the previous value.
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), dimension(:), intent(in) :: cCPZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    integer(kind=AE_INT), intent(in) :: iEqType
    real(kind=AE_REAL), intent(inout) :: rRHS
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(PD0_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_ComputeRHS: PD0_Create has not been called")
      call IO_Assert(io, (associated(fpd)), &
           "PD0_ComputeRHS: Illegal FPD_COLLECTION object")
    end if

    !  do i = 1, pd0%iCount
    !    pnd => pd0%Ponds(i)
    !    select case (iEqType)
    !      case (EQN_HEAD)
    !        rRHS = rRHS - real(cFPD_Potential(fpd, cCPZ(1), pnd%iFPDIndex, 1))
    !      case (EQN_FLOW)
    !        rRHS = rRHS - rFPD_Flow(fpd, cCPZ, pnd%iFPDIndex, 1)
    !      case (EQN_INHO)
    !        rRHS = rRHS - real(cFPD_Potential(fpd, cCPZ(1), pnd%iFPDIndex, 1))
    !      case (EQN_DISCHARGE)
    !        rRHS = rRHS + aimag(cmplx(rZERO, rONE, AE_REAL) * &
        !                      cFPD_Discharge(fpd, cCPZ(1), pnd%iFPDIndex, 1) * &
        !                      cOrientation/abs(cOrientation))
    !      case (EQN_RECHARGE)
    !        rRHS = rRHS
    !      case (EQN_CONTINUITY)
    !        rRHS = rRHS - rFPD_Extraction(fpd, pnd%iFPDIndex, 1)
    !    end select
    !  end do

    return
  end subroutine PD0_ComputeRHS


  subroutine PD0_Read(io, pd0)
    !! subroutine PD0_Read
    !!
    !! Reads the ponds for the specified PD0_COLLECTION from LU_INPUT
    !!
    !! Calling Sequence:
    !!    call PD0_Read(pd0)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be populated
    !!
    !! The format of the PD0 section of the input file appears as follows:
    !! PD0
    !!     x y rsigma id
    !!     ... Up to NPonds
    !! END
    !!
    !! NOTE: It is assumed that the PD0 line was found by the caller
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(IO_STATUS), pointer :: io
    ! [ LOCAL DIRECTIVES ]
    type(DIRECTIVE), dimension(1), parameter :: dirDirectives = (/dirEND/)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rGamma, rRad
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMaxWel
    logical :: lFlag
    type(PD0_POND), pointer :: pnd

    call IO_MessageText(io, "  Reading PD0 module input")

    call IO_Assert(io, (associated(pd0)), "PD0_Read: PD0_Create has not been called")

    ! Process input
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)

      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "PD0_Read: I/O Error")
          exit
        case (kOpFileEOF)
          ! EOF is unexpected for all ModPD0 "ifXXXRead" routines.
          ! Report the condition, but proceed as if EOD was found.
          call IO_Assert(io, .false., "PD0_Read: Unexpected EOF")
        case (kOpData)
          !****************************************************************************
          ! Here for data records
          !****************************************************************************
          call IO_Assert(io, (associated(pd0%Ponds)), "PD0_Read: No space allocated")
          call IO_Assert(io, (pd0%iCount < size(pd0%Ponds)), "PD0_Read: Space exhausted")
          cZ = cIO_GetCoordinate(io, 'cZ', extents=.true.)
          rGamma = rIO_GetReal(io, 'rGamma')
          rRad = rIO_GetReal(io, 'rRad', minimum = rTINY)
          iID = iIO_GetInteger(io, 'iID')

          pd0%iCount = pd0%iCount+1
          pnd => pd0%Ponds(pd0%iCount)
          pnd%cZ = cZ
          pnd%rRadius = rRad
          pnd%rGamma = rGamma
          pnd%iID = iID
          ! No FPD is declared here; see PD0_Setup
          pnd%iFPDIndex = -1
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
      end select
    end do

    call IO_MessageText(io, "  Leaving PD0 module")

    return
  end subroutine PD0_Read


  subroutine PD0_Inquiry(io, pd0, iLU)
    !! subroutine PD0_Inquiry
    !!
    !! Writes an inquiry report for all ponds to iLU
    !!
    !! Calling Sequence:
    !!    call PD0_Inquiry(io, pd0, iLU)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(PD0_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_Inquiry: PD0_Create has not been called")
    end if

    do i = 1, pd0%iCount
      pnd => pd0%Ponds(i)
      write (unit=iLU, &
             fmt="(""PD0"", 2("", "", i9), 4("", "", e16.8))" &
             ) pnd%iID, cIO_WorldCoords(io, pnd%cZ), pnd%rGamma, pnd%rRadius
    end do

    return
  end subroutine PD0_Inquiry


  subroutine PD0_Report(io, pd0)
    !! subroutine PD0_Report
    !!
    !! Writes a debugging report for all ponds to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call PD0_Report(pd0)
    !!
    !! Arguments:
    !!   (in)    type(PD0_COLLECTION), pointer :: pd0
    !!             PD0_COLLECTION to be used
    !!
    ! [ ARGUMENTS ]
    type(PD0_COLLECTION), pointer :: pd0
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    type(PD0_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(pd0)), &
           "PD0_Inquiry: PD0_Create has not been called")
    end if

    call HTML_Header('Module PD0', 1)
    call HTML_Header('Discharge-specified well information', 2)

    if (associated(pd0%Ponds)) then
      call HTML_StartTable()
      call HTML_AttrInteger('Number of wells', pd0%iCount)
      call HTML_AttrInteger('Number of FWL functions', iPD0_GetInfo(io, pd0, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iPD0_GetInfo(io, pd0, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iPD0_GetInfo(io, pd0, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iPD0_GetInfo(io, pd0, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iPD0_GetInfo(io, pd0, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      call HTML_Header('Wells', 4)

      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'ID    ', 'FPD # ', 'X     ', 'Y     ', 'Str   ', 'R     '/))
      do i = 1, pd0%iCount
        pnd => pd0%Ponds(i)
        call HTML_StartRow()
        call HTML_ColumnInteger((/i, pnd%iID, pnd%iFPDIndex/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, pnd%cZ)/))
        call HTML_ColumnReal((/pnd%rGamma, pnd%rRadius/))
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No ponds defined', 3)
    end if

    return
  end subroutine PD0_Report

end module m_pd0

