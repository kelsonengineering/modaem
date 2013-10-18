module m_aqu

  ! ModAEM 2.0
  ! Copyright(c) 1995-2013 Vic Kelson and Layne Christensen Company
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
  ! Contact the author by e-mail at: 
  !     victor.kelson@layne.com or vic.kelson@gmail.com
  !
  ! Or by regular mail at:
  ! 	Vic Kelson
  ! 	Chief Modeler
  ! 	Layne Christensen
  ! 	320 W 8th St
  ! 	Bloomington, IN 47401

  !! module m_layer
  !!
  !! A functioning analytic element layer, containing all the analytic
  !! functions in a model. This module simplifies the implementation of 
  !! groundwater model features, as compared to ModAEM-1.8 and earlier.
  !!
  !! Module m_layer provides a basic container object, AEM_LAYER. It contains
  !!  
  !!   1) a mapping feature (based on piecewise-constant polygons in ModAEM-2.0)
  !!	  that provides lookup capabilities for aquifer properties within then
  !!      layer.
  !!   2) Functions that convert heads and potentials, and discharges and 
  !!      velocities.
  !!   3) Collections of the analytic elements in the model, via function modules
  !!      f_well, f_pond, f_linesink, f_dipole, and f_areasink
  !!   4) Convenience functions for the creation of features in the analytic
  !!      function collections.
  !!   5) It is noted that this module does _not_ manage the creation of boundary
  !!      conditions, nor any matrix-service routines. Those are managed in the
  !!      various feature package modules.
  !!
  !! This design change is a major departure from previous ModAEM versions. It 
  !! serves several purposes:
  !!
  !!   1) Simplifies the implementation of new feature packages (they replace the
  !!      "element modules" in previous versions of ModAEM.
  !!   2) Eliminates the need for the "iterator" scheme for updating boundary
  !!      conditions (while this was elegant, it was a burden for new feature
  !!      development).
  !!   3) Facilitates "ModAEM as a module". In principle, the m_layer module and
  !!      the various function modules are a fully-functional analytic element
  !!      model library component. By the development of a fairly simple wrapper,
  !!      the library could be called from many other languages, with the calling
  !!      language (e.g. Python, C, or Ruby) implementing the matrix solver.
  !!
  !! Module use:
  !!   u_constants --  Universal ModAEM constant declarations
  !!   u_io        --  General-purpose I/O functions
  !!   u_polygon   --  Tools for polygon searching
  !!   f_well      --  Function module for collections of wells
  !!   f_pond      --  Function module for collections of circular area sinks
  !!   f_linesink  --  Function module for collections of line sinks
  !!   f_dipole    --  Function module for collections of line-dipoles
  !!

  use u_constants
  use u_io
  use u_polygon
  use f_well
  use f_pond
  use f_linesink
  use f_dipole

  implicit none

  public

  type, public :: AQU_BDYELEMENT
    !! type AQU_BDYELEMENT
    !!
    !! Type that holds a vertex for the aquifer perimeter
    !!
    !! Members:
    !!   complex :: cZ1
    !!     The complex coordinate of the first vertex(x, y)
    !!   complex :: cZ2
    !!     The complex coordinate of the second vertex(x, y)
    !!   real :: rSpecHead
    !!     The specified head for the segment
    !!   real :: rSpecFlux
    !!     The specified flux for the segment
    !!   logical :: iBdyFlag
    !!     Flag : if .false., base computations on the specified flux
    !!            if .true., base computations on the specified head
    !!   integer :: iFDPIndex
    !!     Pointer in FDP to the dipole
    !!   real :: rLength
    !!     Length of the line segment
    !!   real :: rStrength
    !!     Sink density of the line segment
    !!   real :: rDPStrength
    !!     Dipole strength of the line segment
    !!   integer :: iFWLIndex
    !!     Index into FWL for the segment
    !!
    complex(kind=AE_REAL) :: cZ1
    complex(kind=AE_REAL) :: cZ2
    complex(kind=AE_REAL) :: cCPZ1
    complex(kind=AE_REAL) :: cCPZ2
    complex(kind=AE_REAL) :: cZC
    complex(kind=AE_REAL) :: cZL
    real(kind=AE_REAL) :: rSpecHead
    real(kind=AE_REAL) :: rSpecFlux
    real(kind=AE_REAL) :: rGhbDistance
    real(kind=AE_REAL) :: rLength
    real(kind=AE_REAL) :: rStrength
    integer(kind=AE_INT) :: iBdyFlag
    real(kind=AE_REAL) :: rCheckPot
    real(kind=AE_REAL) :: rCheckFlux
    logical :: lMoveCPZ1                  ! Flag -- should CPZ1 be moved?
    logical :: lMoveCPZ2                  ! Flag -- should CPZ2 be moved?
    logical :: lFSHeadSpec                ! Flag -- is the free surface segment head-specified?
    complex(kind=AE_REAL) :: cFSZ1        ! Location of a "control point" at the first end
    complex(kind=AE_REAL) :: cFSZ2        ! Location of a "control point" at the second end
    integer(kind=AE_INT) :: iEqIndex
    real(kind=AE_REAL) :: rSaveSpecHead   ! Saves the user-specified head
  end type AQU_BDYELEMENT


  type, public :: AQU_COLLECTION
    !! type AQU_COLLECTION
    !!
    !! Type that holds information for a layer
    !!
    !! Members:
    !!   complex :: cRefPoint
    !!     Location of the reference point(if specified)
    !!   real :: rRefHead
    !!     Specified head at the reference point(if specified)
    !!   complex :: cRefUniformFlow
    !!     Infinite aquifer uniform flow discharge vector
    !!   logical :: lReference
    !!     .true. if a reference point has been specified
    !!   real :: rSolConst
    !!     The constant of integration
    !!   type(AQU_BDYELEMENT) :: Boundary(:)
    !!     Vector of vertices which make up the perimeter of the aquifer
    !!   integer :: iNBdy
    !!     Number of perimeter points currently in use
    !!
    complex(kind=AE_REAL) :: cRefPoint
    real(kind=AE_REAL) :: rRefHead
    complex(kind=AE_REAL) :: cRefUniformFlow
    logical :: lReference
    real(kind=AE_REAL) :: rSolConst
    real(kind=AE_REAL) :: rCheck
    type(AQU_BDYELEMENT), dimension(:), pointer :: BdyElements
    integer(kind=AE_INT) :: iNBdy
    logical :: lPrecondition
    type(IN0_COLLECTION), pointer :: in0
    ! Iterator information(aquifers with reference conditions)
    ! Position = -1          - > reset
    ! Position = 0           - > reference potential
    ! 0 < Position < iNBdy   - > boundary element
    ! Position > iNBdy       - > inhomogeneity element
    ! Iterator information(aquifers without reference conditions)
    ! Position = 0           - > reset
    ! 0 < Position < iNBdy   - > boundary element
    ! Position > iNBdy       - > inhomogeneity element
    integer(kind=AE_INT) :: iIterPosition
    integer(kind=AE_INT) :: iIterFlag
    ! Perimeter of the aquifer
    complex(kind=AE_REAL), dimension(:), pointer :: cPerimeter
    integer(kind=AE_INT) :: iNPerimeter
    ! Regen on free surface
    logical :: lFSRegen
    ! Boundary condition locations have been previously set (and possibly moved)
    logical :: lBoundaryInPlace
    ! Turns on the new boundary control-point option
    logical :: lNewBdy
    ! Debug flag
    logical :: lDebug
  end type AQU_COLLECTION

  ! Locally-required constants
  ! Matrix generator flag for reference point equation
  integer(kind=AE_INT), private, parameter :: kAQUReference = 1
  integer(kind=AE_INT), private, parameter :: kAQUBoundary = 2
  ! Fake a small well radius to ensure that the proper fluxes are computed
  ! in FWL_Flow
  real(kind=AE_REAL), private, parameter :: PERIMETER_TOLERANCE = 1.0e-2_AE_REAL
  real(kind=AE_REAL), private, parameter :: BDY_MOVE_FACTOR = 1.0e-4_AE_REAL


contains

  !! ELEMENT MODULE ROUTINES
  !! These routines allow AQU to behave as a ModAEM Element Module


  function AQU_Create(io) result(aqu)
    !! function AQU_Create
    !!
    !! Creates a new AQU_COLLECTION object for an infinite aquifer
    !!
    !! Calling Sequence:
    !!    aqu => AQU_Create(io)
    !!
    !! Arguments:
    !!   type(IO_STATUS), pointer :: io
    !!     All parameters are expected in the buffers in the io object
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT) :: iNInho
    integer(kind=AE_INT) :: iNStr
    real(kind=AE_REAL) :: rBase
    real(kind=AE_REAL) :: rThick
    real(kind=AE_REAL) :: rHydCond
    real(kind=AE_REAL) :: rPorosity
    real(kind=AE_REAL) :: rAvgHead
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    type(AQU_COLLECTION), pointer :: aqu
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(aqu, stat = iStat)
    call IO_Assert(io, (iStat == 0), "AQU_Create: allocation failed")

    iNInho = iIO_GetInteger(io, 'iNInho', minimum = 1)
    iNStr = iIO_GetInteger(io, 'iNStr', minimum = 0)
    rBase = rIO_GetReal(io, 'rBase')
    rThick = rIO_GetReal(io, 'rThick', minimum = rTINY)
    rHydCond = rIO_getReal(io, 'rHydCond', minimum = rZERO)
    rPorosity = rIO_getReal(io, 'rPorosity', minimum = rTINY)
    rAvgHead = rIO_getReal(io, 'rAvgHead', minimum = rBase)
    aqu%in0 => IN0_Create(io, iNInho, iNStr, rBase, rThick, rHydCond, rPorosity, rAvgHead)

    aqu%cRefPoint = cZERO
    aqu%rRefHead = rZERO
    aqu%cRefUniformFlow = cZERO
    aqu%lReference = .false.
    nullify(aqu%BdyElements)
    aqu%iNBdy = 0
    nullify(aqu%cPerimeter)
    aqu%iNPerimeter = 0
    ! The following are determined in AQU_Setup
    aqu%rSolConst = rZERO
    aqu%rCheck = rZERO
    aqu%lPrecondition = .true.
    ! First time up, no regen
    aqu%lFSRegen = .false.
    ! Boundary isn't set yet
    aqu%lBoundaryInPlace = .false.
    ! Default to the old control-point strategy
    aqu%lNewBdy = .false.

    return
  end function AQU_Create


  subroutine AQU_SetPrecondition(io, aqu, lPre)
    ! Sets the precondition flag. If it is set, it is presumed that this is the first iteration, and the
    ! saturated thickness / transmissivity are based on the "avg-head" setting.
    type(AQU_COLLECTION), pointer :: aqu
    logical, intent(in) :: lPre
    type(IO_STATUS), pointer :: io

    call IN0_SetPrecondition(io, aqu%in0, lPre)

    return
  end subroutine AQU_SetPrecondition


  subroutine AQU_SetReference(io, aqu, cRefPoint, rRefHead, cRefUniformFlow)
    !! subroutine AQU_SetReference
    !!
    !! Sets up the reference point and uniform flow
    !!
    !! Calling Sequence:
    !!    call AQU_SetReference(io, aqu, cRefPoint, rRefHead, cRefUniformFlow)
    !!
    !! Arguments:
    !!    (in)    type(AQU_COLLECTION), pointer :: aqu
    !!              The AQU_COLLECTION object to be populated
    !!    (in)    complex :: cRefPoint
    !!              Location of the reference point
    !!    (in)    real :: rRefHead
    !!              The head at the reference point
    !!    (in)    complex :: cRefUniformFlow
    !!              Far-field uniform flow
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cRefPoint
    real(kind=AE_REAL), intent(in) :: rRefHead
    complex(kind=AE_REAL), intent(in) :: cRefUniformFlow
    type(IO_STATUS), pointer :: io

    if (aqu%lDebug) then
      call IO_Assert(io, associated(aqu), "AQU_SetReference: No AQU_COLLECTION object")
    end if

    aqu%cRefPoint = cRefPoint
    aqu%rRefHead = rRefHead
    aqu%cRefUniformFlow = cRefUniformFlow
    aqu%lReference = .true.

    return
  end subroutine AQU_SetReference


  subroutine AQU_PreSolve(io, aqu)
    !! subroutine AQU_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call AQU_PreSolve(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, istat

    ! Fix up the perimeter stuff as necessary
    if (aqu%iNBdy > 0) then
      allocate(aqu%cPerimeter(aqu%iNBdy), stat = istat)
      call IO_Assert(io, (istat == 0), "Could not allocate aquifer perimeter")
      call PGN_Build(aqu%BdyElements(1:aqu%iNBdy)%cZ1, &
           aqu%BdyElements(1:aqu%iNBdy)%cZ2, &
           aqu%BdyElements(1:aqu%iNBdy)%iBdyFlag==BDY_HEAD, &
           PERIMETER_TOLERANCE, &
           aqu%cPerimeter, &
           aqu%iNPerimeter)
      call IO_Assert(io, (aqu%iNPerimeter > 0), "Could not build the aquifer perimeter")
    else
      allocate(aqu%cPerimeter(0), stat = istat)
      aqu%iNPerimeter = 0
    end if

    ! Set up the IN0 module stuff
    call IN0_PreSolve(io, aqu%in0)

    return
  end subroutine AQU_PreSolve


  function iAQU_GetInfo(io, aqu, iOption, iIteration) result(iValue)
    !! function AQU_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iAQU_GetInfo(io, aqu, iOption)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iStr

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_GetInfo: AQU_Create has not been called")
    end if

    iValue = iIN0_GetInfo(io, aqu%in0, iOption, iIteration)
    !  print *, 'IN0', iValue
    select case (iOption)
      case (SIZE_FWL)
        continue
      case (SIZE_FDP)
        continue
      case (SIZE_EQUATIONS)
        iValue = iValue + aqu%iNBdy
        if (aqu%lReference .or. aqu%iNBdy == 0) iValue = iValue + 1
      case (SIZE_UNKNOWNS)
        iValue = iValue + aqu%iNBdy
        if (aqu%lReference .or. aqu%iNBdy == 0) iValue = iValue + 1
      case (INFO_REGENERATE)
        if (iIteration < 2 .or. aqu%lFSRegen) then
          iValue = 1
        end if
      case default
        iValue = 0
    end select

    return
  end function iAQU_GetInfo


  subroutine AQU_SetupFunctions(io, aqu, fdp)
    !! subroutine AQU_SetupFunction
    !!
    !! This routine sets up the functions in f_well and f_dipole for the perimeter
    !! and sets up the equation for the reference point.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call AQU_SetupFunction(io, iL)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fwl
    !!             FDP_COLLECTION object where functions may be added
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cCPZ1, cCPZ2
    integer(kind=AE_INT) :: iBdy
    integer(kind=AE_INT) :: iNWL, iNPD, iNDP, iNEQ, iNUN
    type(AQU_BDYELEMENT), pointer :: this

    call IO_Assert(io, (associated(aqu)), "AQU_Setup: AQU_Create has not been called")

    ! Set the internal values for the BDY elements
    if (aqu%iNBdy > 0) then
      aqu%BdyElements(1:aqu%iNBdy)%cZC = rHALF * (aqu%BdyElements(1:aqu%iNBdy)%cZ2 + aqu%BdyElements(1:aqu%iNBdy)%cZ1)
      aqu%BdyElements(1:aqu%iNBdy)%cZL = rHALF * (aqu%BdyElements(1:aqu%iNBdy)%cZ2 - aqu%BdyElements(1:aqu%iNBdy)%cZ1)
    end if

    ! Set the confined potential
    call IN0_SetupFunctions(io, aqu%in0, fdp)

    return
  end subroutine AQU_SetupFunctions


  subroutine AQU_SetupMatrix(io, aqu, mat)
    !! subroutine AQU_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the module
    !! and sets up the equation for the reference point.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call AQU_SetupMatrix(io, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer :: mat
    !!             MAT_MATRIX object for the solver
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cCPZ0, cCPZ1, cCPZ2, cNormal, cN1, cN2, cUnitOffset
    real(kind=AE_REAL) :: rElev, rOffsetDistance
    integer(kind=AE_INT) :: iBdy
    integer(kind=AE_INT) :: iNWL, iNPD, iNDP, iNEQ, iNUN, iEQ
    type(AQU_BDYELEMENT), pointer :: this, prev, next

    call IO_Assert(io, (associated(aqu)), "AQU_Setup: AQU_Create has not been called")

    ! Set the confined potential
    call IN0_SetupMatrix(io, aqu%in0, mat)

    ! Build the matrix entries for the perimeter(if any)
    if (aqu%iNBdy > 0) then
      prev => aqu%BdyElements(aqu%iNBdy)
      rOffsetDistance = BDY_MOVE_FACTOR * sum(aqu%BdyElements(1:aqu%iNBdy)%rLength) / aqu%iNBdy
      do iBdy = 1, aqu%iNBdy
        this => aqu%BdyElements(iBdy)
        call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUBoundary, iBdy, 0)

        ! The locations for flux inquiries are moved "inside" the element
        cNormal = cI * (this%cZ2 - this%cZ1)
        if (.not. aqu%lBoundaryInPlace) then
          if (.not. aqu%lNewBdy) then
            this%cCPZ1 = this%cZ1 + (this%cZ2-this%cZ1)*BDY_MOVE_FACTOR + &
                             cNormal*BDY_MOVE_FACTOR*abs(this%cZ2-this%cZ1)
            this%cCPZ2 = this%cZ2 - (this%cZ2-this%cZ1)*BDY_MOVE_FACTOR +  &
                             cNormal*BDY_MOVE_FACTOR*abs(this%cZ2-this%cZ1)
          else
            ! Find the unit-normals to this element and the next...
            call AQU_GetNeighborBdy(io, aqu, this, prev, next)
            cN1 = cI * (this%cZ2 - this%cZ1) / abs(this%cZ2 - this%cZ1)
            cN2 = cI * (next%cZ2 - next%cZ1) / abs(next%cZ2 - next%cZ1)
            ! The unit vector pointing from this element's end 2 to the control point
            ! bisects the direction of cN1 and cN2
            cUnitOffset = (cN1 + cN2) / abs(cN1 + cN2)
            this%cCPZ2 = this%cZ2 + rOffsetDistance * cUnitOffset
            next%cCPZ1 = next%cZ1 + rOffsetDistance * cUnitOffset
!            write(unit=*, fmt="('bdy', 4(1x,'(',2f10.3,')'),i3)") this%cZ1, this%cZ2, &
!                & next%cZ2, this%cCPZ2, PGN_CheckPoint(aqu%cPerimeter,this%cCPZ2)
          end if
          ! Save the starting locations for the vertices
          this%cFSZ1 = this%cCPZ1
          this%cFSZ2 = this%cCPZ2
        end if

        select case (this%iBdyFlag)
          case (BDY_HEAD)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/rHALF*(this%cCPZ1+this%cCPZ2)/), &
                                                   EQN_HEAD, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
          case (BDY_FLUX)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                   EQN_FLOW, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
          case (BDY_GHB)
            this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                   EQN_BDYGHB, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, &
                 this%rGhbDistance/this%rLength)
          case (BDY_FREESURF)
            ! For free-surface elements, choose the proper initial boundary condition
            !write(unit=*,fmt="('MAT',i10,1x,i1,1x,l1,4(e16.8,1x))") ibdy, this%iBdyFlag, this%lFSHeadSpec, this%cCPZ1, this%cCPZ2
            if (.not. aqu%lBoundaryInPlace) then
              rElev = max(aimag(this%cFSZ1), aimag(this%cFSZ2))
              if (this%rSpecHead >= rElev) then
                this%lFSHeadSpec = .true.
              else
                this%lFSHeadSpec = .false.
              end if
            end if

            if (this%lFSHeadSpec) then
              ! It's a head
              this%iEqIndex = MAT_CreateEquation(io, mat, (/rHALF*(this%cCPZ1+this%cCPZ2)/), &
                                                     EQN_HEAD, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
            else
              ! It's a moveable boundary
              this%iEqIndex = MAT_CreateEquation(io, mat, (/this%cCPZ1, this%cCPZ2/), &
                                                     EQN_FLOW, ELEM_AQU, kAQUBoundary, iBdy, 0, cZERO, rZERO)
            end if
        end select
        prev => aqu%BdyElements(iBdy)
      end do
      aqu%lBoundaryInPlace = .true.

    end if
    ! Build a matrix entry for the reference point(if any)
    !print *, aqu%lReference, aqu%iNBdy
    if (aqu%lReference) then
      call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUReference, 0, 0)
      iEQ = MAT_CreateEquation(io, mat, (/aqu%cRefPoint/), EQN_HEAD, ELEM_AQU, kAQUReference, 0, 0, cZERO, rZERO)
    else if (aqu%iNBdy == 0) then
      call MAT_CreateVariable(io, mat, ELEM_AQU, kAQUReference, 0, 0)
      iEQ = MAT_CreateEquation(io, mat, (/cZERO/), EQN_CONTINUITY, ELEM_AQU, kAQUReference, 0, 0, cZERO, rZERO)
    end if

    return
  end subroutine AQU_SetupMatrix


  function iAQU_Prepare(io, aqu, iIteration) result(iChanges)
    !! subroutine AQU_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_aqu
    !!
    !! Calling Sequence:
    !!    call AQU_Prepare(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iAQU_Prepare


  subroutine AQU_GetNeighborBdy(io, aqu, this, prev, next)
    !! For profile models, returns 'prev' and 'next' as the previous and
    !! next (corresponding to cZ1 and cZ2, respectively) BDY elements adjoining
    !! the element 'this'.
    !! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    type(AQU_BDYELEMENT), pointer :: this
    type(AQU_BDYELEMENT), pointer :: prev
    type(AQU_BDYELEMENT), pointer :: next
    !! [ LOCALS ]
    integer(kind=AE_INT) :: ibdy
    type(AQU_BDYELEMENT), pointer :: bdy

    ! Find the neighbors
    prev => NULL()
    next => NULL()
    do ibdy=1,aqu%iNBdy
      bdy => aqu%BdyElements(ibdy)
      if (abs(this%cZ1 - bdy%cZ2) < PERIMETER_TOLERANCE) prev => aqu%BdyElements(ibdy)
      if (abs(this%cZ2 - bdy%cZ1) < PERIMETER_TOLERANCE) next => aqu%BdyElements(ibdy)
    end do

    call IO_Assert(io, (associated(prev) .and. associated(next)), "Failed to find neighboring BDY elements")

    return
  end subroutine AQU_GetNeighborBdy


  function rAQU_GetCoefficientMultiplier(io, aqu, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    rMultiplier = rONE

    return
  end function rAQU_GetCoefficientMultiplier


  subroutine AQU_ComputeCoefficients(io, aqu, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor, &
               iIteration, rMultiplier, rARow)
    !! subroutine HB0_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the HB0
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call AQU_ComputeCoefficients(io, aqu, fwl, fdp, mat, cPathZ, iEqType, cOrientation, rMultiplier, rRow)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer :: fwl
    !!             FWL_COLLECTION object where functions may be added
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object where functions may be added
    !!   (in)    complex :: cPathZ(:)
    !!             The control point(or control path) to be used
    !!   (in)    integer :: iEqType
    !!             The equation type
    !!   (in)    integer :: iElementType
    !!             The element type that created this equation
    !!   (in)    integer :: iElementString
    !!             The element string corresponding to this equation
    !!   (in)    integer :: iElementVertex
    !!             The element vertex corresponding to this equation
    !!   (in)    integer :: iElementFlag
    !!             The element flag(if any) for this equation
    !!   (in)    complex :: cOrientation
    !!             Orientation unit vector(for discharge-based equations)
    !!   (out)   real :: rARow(:)
    !!             The output row of coefficients(to be concatenated with
    !!             row portions for other element modules.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    complex(kind=AE_REAL), intent(in) :: cOrientation
    real(kind=AE_REAL), intent(in) :: rGhbFactor
    integer(kind=AE_INT), intent(in) :: iIteration
    real(kind=AE_REAL), intent(in) :: rMultiplier
    real(kind=AE_REAL), dimension(:), intent(out) :: rARow
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStat, iCol, iBdy, iIN0NUN
    complex(kind=AE_REAL) :: cZC, cZL, cBigZ1, cBigZ2, cUnit
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    real(kind=AE_REAL), dimension(1, 1) :: rTemp
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    type(AQU_BDYELEMENT), pointer :: bdy

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_ComputeCoefficients: No AQU_COLLECTION object")
    end if
    ! Set it up
    rARow = rZERO

    ! Build entries for the inhomogeneities
    iIN0NUN = iIN0_GetInfo(io, aqu%in0, SIZE_UNKNOWNS, 0)
    call IN0_ComputeCoefficients(io, aqu%in0, fdp, cPathZ, iEqType, iElementType, iElementString, &
         iElementVertex, iElementFlag, cOrientation, rGhbFactor, &
         iIteration, rMultiplier, rARow(1:iIN0NUN))
    iCol = iIN0NUN + 1

    ! Build entries for the boundary elements
    if (associated(aqu%BdyElements)) then
      ! Get the appropriate influence functions for the boundary condition type
      allocate(cBigZ(size(cPathZ)), stat = iStat)
      call IO_Assert(io, (iStat == 0), 'AQU_ComputeCoefficients: Allocation failed for big-Z')
      do iBdy = 1, aqu%iNBdy
        bdy => aqu%BdyElements(iBdy)
        cZC = rHALF*(bdy%cZ2+bdy%cZ1)
        cZL = rHALF*(bdy%cZ2-bdy%cZ1)
        cBigZ = (cPathZ - cZC) / cZL
        select case (iEqType)
          case (EQN_HEAD)
            cTemp = real(cILS_InfluenceP(cBigZ(1), cZL))
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_BDYGHB)
            cTemp = real(cILS_InfluenceP(rHALF*(cBigZ(1)+cBigZ(2)), cZL))
            cTemp = cTemp + rGhbFactor * real(cILS_InfluenceF(cBigZ(1), cBigZ(2), cZL))
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_FLOW)
            do i = 1, size(cBigZ)-1
              cTemp = cILS_InfluenceF(cBigZ(i), cBigZ(i+1), cZL)
              rARow(iCol) = rMultiplier * real(cTemp(1, 1))
            end do
          case (EQN_INHO)
            cTemp = cILS_InfluenceP(cBigZ(1), cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_DISCHARGE)
            cUnit = cOrientation / abs(cOrientation)
            cTemp = cILS_InfluenceW(cBigZ(1), cZL)
            rARow(iCol) = rMultiplier * real(cUnit * cTemp(1, 1))
          case (EQN_POTENTIALDIFF)
            cTemp = cILS_InfluenceP(cBigZ(1), cZL)-cILS_InfluenceP(cBigZ(2), cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_CONTINUITY)
            cBigZ1 = (bdy%cCPZ1 - cZC) / cZL
            cBigZ2 = (bdy%cCPZ2 - cZC) / cZL
            cTemp = cILS_InfluenceF(cBigZ1, cBigZ2, cZL)
            rARow(iCol) = rMultiplier * real(cTemp(1, 1))
          case (EQN_TOTALFLOW)
            rARow(iCol) = rZERO
        end select

        iCol = iCol+1
      end do
      deallocate(cBigZ)
    end if

    ! Set up the closure condition
    if (aqu%lReference .or. aqu%iNBdy == 0) then
      select case (iEqType)
        case (EQN_HEAD)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_BDYGHB)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_FLOW)
          rARow(iCol) = rZERO
        case (EQN_INHO)
          rARow(iCol) = rMultiplier * rONE
        case (EQN_DISCHARGE)
          rARow(iCol) = rZERO
        case (EQN_RECHARGE)
          rARow(iCol) = rZERO
        case (EQN_CONTINUITY)
          rARow(iCol) = rZERO
      end select
      iCol = iCol+1
    end if

    return
  end subroutine AQU_ComputeCoefficients


  function rAQU_ComputeRHS(io, aqu, fdp, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rAQU_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rAQU_ComputeRHS(io, aqu, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag, rSpecValue, rCheck)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    integer :: iElementType
    !!             Element type(either ELAM_AQU or ELEM_IN0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    !! Return Value:
    !!   real :: rRHS
    !!     The RHS value for the module
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rRHS
    type(AQU_BDYELEMENT), pointer :: bdy

    select case (iElementType)
      case (ELEM_AQU)
        select case (iElementString)
          case (kAQUBoundary)
            bdy => aqu%BdyElements(iElementVertex)
            select case (bdy%iBdyFlag)
              case (BDY_HEAD)
                if (lDirect) then
                  rRHS = rIN0_HeadToPotential(io, aqu%in0, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2))
                else
                  rRHS = rIN0_HeadToPotential(io, aqu%in0, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - bdy%rCheckPot
                end if
              case (BDY_FLUX)
                if (lDirect) then
                  rRHS = bdy%rSpecFlux
                else
                  rRHS = bdy%rSpecFlux - bdy%rCheckFlux
                end if
              case (BDY_GHB)
                rRHS = rIN0_HeadToPotential(io, aqu%in0, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - &
                       bdy%rCheckPot - bdy%rCheckFlux*bdy%rGhbDistance/bdy%rLength
              case (BDY_FREESURF)
                if (bdy%lFSHeadSpec) then
                  if (lDirect) then
                    rRHS = rIN0_HeadToPotential(io, aqu%in0, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2))
                  else
                    rRHS = rIN0_HeadToPotential(io, aqu%in0, bdy%rSpecHead, rHALF*(bdy%cZ1+bdy%cZ2)) - bdy%rCheckPot
                  end if
                else
                  if (lDirect) then
                    rRHS = bdy%rSpecFlux
                  else
                    rRHS = bdy%rSpecFlux - bdy%rCheckFlux
                  end if
                end if
            end select
          case (kAQUReference)
            if (aqu%lReference) then
              if (lDirect) then
                rRHS = rIN0_HeadToPotential(io, aqu%in0, aqu%rRefHead, aqu%cRefPoint)
              else
                rRHS = rIN0_HeadToPotential(io, aqu%in0, aqu%rRefHead, aqu%cRefPoint) - aqu%rCheck
              end if
            else
              rRHS = rZERO
            end if
        end select
      case (ELEM_IN0)
        rRHS = rIN0_ComputeRHS(io, aqu%in0, fdp, iEqType, iElementType, iElementString, iElementVertex, &
               iElementFlag, iIteration, lDirect)
    end select

    return
  end function rAQU_ComputeRHS


  subroutine AQU_StoreResult(io, aqu, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine AQU_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the LS1 module.
    !!
    !! Calling Sequence:
    !!    LS1_StoreResult(io, aqu, rValue, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(either ELAM_AQU or ELEM_IN0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(AQU_BDYELEMENT), pointer :: this

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), "AQU_Update: AQU_Create has not been called")
    end if

    select case (iElementType)
      case (ELEM_AQU)
        select case (iElementString)
          case (kAQUReference)
            if (lDirect) then
              aqu%rSolConst = rValue
            else
              aqu%rSolConst = aqu%rSolConst + rValue
            end if
          case (kAQUBoundary)
            this => aqu%BdyElements(iElementVertex)
            if (lDirect) then
              this%rStrength = rValue
            else
              this%rStrength = this%rStrength + rValue
            end if
        end select
      case (ELEM_IN0)
        call IN0_StoreResult(io, aqu%in0, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    end select

    return
  end subroutine AQU_StoreResult


  subroutine AQU_Update(io, aqu, fdp)
    !! subroutine AQU_Update
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    AQU_Update(io, aqu, fwl, fdp)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(fdp_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iBdy, iRV
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(AQU_BDYELEMENT), pointer :: this

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), "AQU_Update: AQU_Create has not been called")
    end if

    ! Update the inhomogeneities
    call IN0_Update(io, aqu%in0, fdp)

    return
  end subroutine AQU_Update


  subroutine AQU_ResetIterator(io, aqu)
    !! subroutine AQU_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call AQU_ResetIterator(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: aqu
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_ResetIterator: AQU_Create has not been called")
    end if

    aqu%iIterPosition = 0
    aqu%iIterFlag = VALUE_POTENTIAL
    call IN0_ResetIterator(io, aqu%in0)

    return
  end subroutine AQU_ResetIterator


  function AQU_NextIterator(io, aqu) result(itr)
    !! function AQU_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call AQU_NextIterator(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: aqu
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr
    type(AQU_BDYELEMENT), pointer :: bdy
    integer(kind=AE_INT) :: iStat

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_NextIterator: AQU_Create has not been called")
    end if

    aqu%iIterPosition = aqu%iIterPosition+1
    if (aqu%iIterPosition == aqu%iNBdy+1) then
      ! Handle the reference information
      if (aqu%lReference) then
        allocate(itr, stat = iStat)
        call IO_Assert(io, (iStat == 0), "AQU_NextIterator: Space Exhausted")
        itr%iElementType = ELEM_AQU
        itr%iElementString = 0
        itr%iElementVertex = 0
        itr%iElementFlag = 0
        itr%iValueSelector = VALUE_POTENTIAL
        allocate(itr%cZ(1))
        itr%cZ(1) = aqu%cRefPoint
        else!if (aqu%iNBdy == 0) then
        allocate(itr, stat = iStat)
        call IO_Assert(io, (iStat == 0), "AQU_NextIterator: Space Exhausted")
        itr%iElementType = ELEM_AQU
        itr%iElementString = 0
        itr%iElementVertex = 0
        itr%iElementFlag = 0
        itr%iValueSelector = VALUE_EXTRACTION
        allocate(itr%cZ(0))
      end if
      ! Prepare for the boundary elements(see below)
    else if (aqu%iNBdy > 0 .and. aqu%iIterPosition <= aqu%iNBdy) then
      ! Handle the perimeter elements
      bdy => aqu%BdyElements(aqu%iIterPosition)
      select case (aqu%iIterFlag)
        case (VALUE_POTENTIAL)
          ! Head condition
          allocate(itr, stat = iStat)
          call IO_Assert(io, (iStat == 0), "AQU_NextIterator: Space Exhausted")
          itr%iElementType = ELEM_AQU
          itr%iElementString = 1
          itr%iElementVertex = aqu%iIterPosition
          itr%iValueSelector = VALUE_POTENTIAL
          allocate(itr%cZ(1))
          itr%cZ(1) = rHALF * (bdy%cCPZ1+bdy%cCPZ2)
          aqu%iIterFlag = VALUE_FLOW
          ! Do not advance until the flow is extracted(see below)
          aqu%iIterPosition = aqu%iIterPosition-1
        case (VALUE_FLOW)
          ! Flux condition
          allocate(itr, stat = iStat)
          call IO_Assert(io, (iStat == 0), "AQU_NextIterator: Space Exhausted")
          itr%iElementType = ELEM_AQU
          itr%iElementString = 1
          itr%iElementVertex = aqu%iIterPosition
          itr%iValueSelector = VALUE_FLOW
          allocate(itr%cZ(2))
          itr%cZ = (/bdy%cCPZ1, bdy%cCPZ2/)
          aqu%iIterFlag = VALUE_POTENTIAL
      end select
    else
      itr => IN0_NextIterator(io, aqu%in0)
    end if

    return
  end function AQU_NextIterator


  subroutine AQU_SetIterator(io, aqu, fdp, itr, cValue, lLinearize)
    !! function AQU_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call AQU_SetIterator(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION to be used
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: aqu
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(AQU_BDYELEMENT), pointer :: bdy
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cUnit

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_NextIterator: AQU_Create has not been called")
    end if

    if (itr%iElementType == ELEM_AQU) then
      if (itr%iElementString == 0) then
        ! Handle the reference point
        aqu%rCheck = real(cValue, AE_REAL)
      else if (itr%iElementString == 1) then
        ! Handle the perimeter elements
        select case (itr%iValueSelector)
          case (VALUE_POTENTIAL)
            bdy => aqu%BdyElements(itr%iElementVertex)
            bdy%rCheckPot = real(cValue, AE_REAL)
          case (VALUE_FLOW)
            bdy => aqu%BdyElements(itr%iElementVertex)
            bdy%rCheckFlux = real(cValue, AE_REAL)
            !if (bdy%iBdyFlag == BDY_FREESURF .and. .not. bdy%lFSHeadSpec) then
            !  write(unit=*,fmt="('set', i10, 2(1x,'(',f12.4,1x,f12.4,')'),1x,e14.6)") itr%iElementVertex, itr%cZ, real(cValue,AE_REAL)
            !end if
        end select
      end if
    else if (itr%iElementType == ELEM_IN0) then
      call IN0_SetIterator(io, aqu%in0, fdp, itr, cValue, lLinearize)
    end if

    return
  end subroutine AQU_SetIterator


  subroutine AQU_Read(io, aqu)
    !! subroutine AQU_Read
    !!
    !! Reads the aquifer information for layer iL from the input LU
    !!
    !! Calling Sequence:
    !!    call AQU_Read(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             Layer number to be read
    !!   (out)   integer :: iRV
    !!             Error code or errOK for success
    !!
    !! The format of the aqu section of the input file appears as follows:
    !!
    !! AQU ninho base thickness hyd-cond porosity
    !! REF(z0) refhead(Qx0, Qy0)    [ optional ]
    !! BDY nvertices                 [ optional ]
    !!     (x1, y1) (x2, y2) head flux ghb-distance flag
    !!     (x1, y1) (x2, y2) head flux ghb-distance flag
    !!     ... up to nvertices
    !! END
    !! IN0 [optional]
    !! DOM nvertices base thickness hyd-cond porosity
    !!     (x, y)
    !!     (x, y)
    !!     ... up to nvertices
    !! DOM ...
    !! ... up to ninho
    !! END
    !!
    !! NOTE: It is assumed that the AQU line was found by the caller
    !!
    !! The BDY section(perimeter boundary) requests the specified flux or head
    !! values, based on the 'flag' value. If flag = 0, use the flux boundary; if
    !! flag = 1, uses the head boundary. Assumes that for vertex 'i', the flux
    !! or head specified is that of segment 'i', if flag = 2, a GHB is used
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! Locals -- for Directive parsing
    type(DIRECTIVE), dimension(6), parameter :: dirDirectives = &
                       (/dirEND, dirREF, dirSWI, dirBDY, dirIN0, dirOPT/)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMax
    real(kind=AE_REAL) :: rBase
    real(kind=AE_REAL) :: rThickness
    real(kind=AE_REAL) :: rHydCond
    real(kind=AE_REAL) :: rPorosity
    real(kind=AE_REAL) :: rRefHead
    complex(kind=AE_REAL) :: cRefPoint
    complex(kind=AE_REAL) :: cRefUniformFlow
    complex(kind=AE_REAL) :: cZ, cZ1, cZ2
    real(kind=AE_REAL) :: rSpecHead
    real(kind=AE_REAL) :: rSpecFlux
    real(kind=AE_REAL) :: rGhbDistance
    real(kind=AE_REAL) :: rInactiveValue
    integer(kind=AE_INT) :: iBdyFlag
    integer(kind=AE_INT) :: iNBdy
    integer(kind=AE_INT) :: iNIsl
    type(AQU_BDYELEMENT), pointer :: vtx
    ! State variables for the parser
    integer(kind=AE_INT) :: iParseMode
    integer(kind=AE_INT), parameter :: PARSE_NONE = 0
    integer(kind=AE_INT), parameter :: PARSE_BDY = 1
    type(IN0_DOMAIN), pointer :: dom
    character(len=132) :: sOption

    call IO_MessageText(io, "  Reading AQU module input")
    call IO_Assert(io, (associated(aqu)), "AQU_Read: AQU_Create has not been called")

    iParseMode = PARSE_NONE
    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "AQU_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all Mod "ifXXXRead" routines.
          call IO_Assert(io, .false., "AQU_Read: Unexpected EOF")
        case (kOpData)
          ! A data line was found. If we have a specified perimeter, add the point
          ! to the perimeter.
          select case (iParseMode)
            case (PARSE_BDY)
              call IO_Assert(io, (aqu%iNBdy < size(aqu%BdyElements)), &
                   "AQU_Read: Space exhausted")
              cZ1 = cIO_GetCoordinate(io, 'cZ1', extents=.true.)
              cZ2 = cIO_GetCoordinate(io, 'cZ2', extents=.true.)
              rSpecHead = rIO_GetReal(io, 'rSpecHead')
              rSpecFlux = rIO_GetReal(io, 'rSpecFlux')
              rGhbDistance = rIO_GetReal(io, 'rGhbDistance', minimum = rZERO)
              iBdyFlag = iIO_GetInteger(io, 'iBdyFlag')
              call IO_Assert(io, (iBdyFlag == BDY_HEAD .or. iBdyFlag == BDY_FLUX .or. &
                                  iBdyFlag == BDY_GHB .or. iBdyFlag == BDY_FREESURF), &
                   'Illegal boundary type')
              aqu%iNBdy = aqu%iNBdy+1
              vtx => aqu%BdyElements(aqu%iNBdy)
              vtx%cZ1 = cZ1
              vtx%cZ2 = cZ2
              vtx%rLength = abs(cZ2-cZ1)
              vtx%rSpecHead = rSpecHead
              vtx%rSaveSpecHead = rSpecHead
              vtx%rSpecFlux = rSpecFlux
              vtx%rGhbDistance = rGhbDistance
              vtx%iBdyFlag = iBdyFlag
          end select
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          exit
        case (kOpREF)
          ! Reference record was found. Read the reference point and uniform flow rate
          cRefPoint = cIO_GetCoordinate(io, 'cRefPoint', extents=.true.)
          rRefHead = rIO_GetReal(io, 'rRefHead')
          cRefUniformFlow = cIO_GetComplex(io, 'cUniformFlow')
          dom => IN0_FindDomain(io, aqu%in0, cRefPoint)
          call IO_Assert(io, (rRefHead > dom%rBase), "AQU_Read: Base elevation below aquifer top")
          aqu%cRefPoint = cRefPoint
          aqu%rRefHead = rRefHead
          aqu%cRefUniformFlow = cRefUniformFlow
          aqu%lReference = .true.
        case (kOpBDY)
          ! Have we already specified boundary elements?
          call IO_Assert(io, (.not. associated(aqu%BdyElements)), &
               "AQU_Read: BDY statement has already been encountered")
          iMax = iIO_GetInteger(io, 'iMax')
          call IO_Assert(io, (iMax > 0), "AQU_Read: Illegal dimension")
          allocate(aqu%BdyElements(iMax), stat = iStat)
          call IO_Assert(io, (iStat == 0), "AQU_Read: Allocation failed")
          aqu%BdyElements = AQU_BDYELEMENT(cZERO, cZERO, cZERO, cZERO, cZERO, cZERO, rZERO, &
                            rZERO, rZERO, rZERO, rZERO, BDY_FLUX, rZERO, rZERO, .false., .false., .false., &
                            cZERO, cZERO, 0, rZERO)
          aqu%iNBdy = 0
          iParseMode = PARSE_BDY
        case (kOpOPT)
          sOption = sIO_GetField(io, "sOption", allowed=(/"NEWBOUNDARY"/), force_uppercase=.true.)
          if (sOption == "NEWBOUNDARY") then
            aqu%lNewBdy = .true.
            call IO_MessageText(io, "New boundary condition logic is enabled")
          end if
        case (kOpIN0)
          call IN0_Read(io, aqu%in0)
        case default
      end select
    end do
    call IO_MessageText(io, "  Leaving AQU module")

    return
  end subroutine AQU_Read


  subroutine AQU_Report(io, aqu)
    ! Reports aquifer information for layer iL
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN, iBdy
    type(AQU_BDYELEMENT), pointer :: bdy
    real(kind=AE_REAL) :: rCheckHead

    ! Aquifer data from this module
    call HTML_Header('Module AQU', 1)
    call HTML_Header('Basic aquifer information', 2)

    if (aqu%lReference) then
      call HTML_Header('Reference aquifer', 3)
      call HTML_StartTable()
      call HTML_AttrComplex('Reference point', aqu%cRefPoint)
      call HTML_AttrReal('Reference head', aqu%rRefHead)
      call HTML_AttrComplex('Uniform flow', aqu%cRefUniformFlow)
      call HTML_AttrReal('Solution constant', aqu%rSolConst)
      call HTML_EndTable()
    else
      call HTML_Header('No reference point(continuity of flow condition)', 3)
      call HTML_StartTable()
      call HTML_AttrComplex('Reference point', cIO_WorldCoords(io, aqu%cRefPoint))
      call HTML_AttrReal('Reference head', aqu%rRefHead)
      call HTML_AttrComplex('Uniform flow', aqu%cRefUniformFlow)
      call HTML_AttrReal('Solution constant', aqu%rSolConst)
      call HTML_EndTable()
    end if

    call HTML_Header('Usage', 3)
    call HTML_StartTable()
    call HTML_AttrInteger('Number of boundary elements', aqu%iNBdy)
    call HTML_AttrInteger('Number of FWL functions', iAQU_GetInfo(io, aqu, SIZE_FWL, 0))
    call HTML_AttrInteger('Number of FPD functions', iAQU_GetInfo(io, aqu, SIZE_FPD, 0))
    call HTML_AttrInteger('Number of FDP functions', iAQU_GetInfo(io, aqu, SIZE_FDP, 0))
    call HTML_AttrInteger('Number of equations', iAQU_GetInfo(io, aqu, SIZE_EQUATIONS, 0))
    call HTML_AttrInteger('Number of unknowns', iAQU_GetInfo(io, aqu, SIZE_UNKNOWNS, 0))
    call HTML_EndTable()

    ! Report the boundary elements
    if (aqu%iNBdy > 0) then
      call HTML_Header('Boundary elements', 3)
      call HTML_StartTable()
      call HTML_TableHeader((/'Index   ', 'X1      ', 'Y1      ', 'X2      ', 'Y2      ', 'Flag    ', &
           'Spec Hd ', 'Spec Fl ', 'GHB Dis ', 'Strength', 'M Head  ', 'M Flux  ', &
           'Head Err', 'Flux Err'/))
      do iBdy = 1, aqu%iNBdy
        call HTML_StartRow()
        bdy => aqu%BdyElements(iBdy)
        rCheckHead = rAQU_PotentialToHead(io, aqu, bdy%rCheckPot, rHALF*(bdy%cZ1+bdy%cZ2))
        call HTML_ColumnInteger((/iBdy/))
        call HTML_ColumnComplex((/cIO_WorldCoords(io, bdy%cZ1), cIO_WorldCoords(io, bdy%cZ2)/))
        call HTML_ColumnInteger((/bdy%iBdyFlag/))
        call HTML_ColumnReal((/ bdy%rSpecHead, &
                                rIO_WorldLength(io, bdy%rSpecFlux, bdy%cZ2-bdy%cZ1), &
                                bdy%rGhbDistance, bdy%rStrength, &
                                rCheckHead, &
                                rIO_WorldLength(io, bdy%rCheckFlux, bdy%cZ2-bdy%cZ1) /))
        if (bdy%iBdyFlag == BDY_HEAD) then
          call HTML_ColumnReal((/rCheckHead - bdy%rSpecHead/))
          call HTML_ColumnText((/'--'/))
        else if (bdy%iBdyFlag == BDY_FLUX) then
          call HTML_ColumnText((/'--'/))
          call HTML_ColumnReal((/bdy%rCheckFlux - bdy%rSpecFlux/))
        else if (bdy%iBdyFlag == BDY_GHB) then
          call HTML_ColumnText((/'--'/))
          call HTML_ColumnText((/'--'/))
        else if (bdy%iBdyFlag == BDY_FREESURF) then
          if (bdy%lFSHeadSpec) then
            call HTML_ColumnReal((/rCheckHead - bdy%rSpecHead/))
            call HTML_ColumnText((/'--'/))
          else
            call HTML_ColumnText((/'--'/))
            call HTML_ColumnReal((/bdy%rCheckFlux - bdy%rSpecFlux/))
          end if
        end if
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No boundary elements', 3)
    end if

    ! Write out the inhomogeneities for property definitions
    call IN0_Report(io, aqu%in0)

    return
  end subroutine AQU_Report


  subroutine AQU_Inquiry(io, aqu, iLU)
    !! subroutine AQU_Inquiry
    !!
    !! Writes an inquiry report for all barriers to iLU
    !!
    !! Calling Sequence:
    !!    call AQU_Inquiry(io, aqu, iLU)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!
    ! [ aqu%lDebug ]
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    if (aqu%lDebug) then
      call IO_Assert(io, (associated(aqu)), &
           "AQU_Inquiry: AQU_Create has not been called")
    end if


    if (aqu%lReference) then
      write (unit=iLU, &
             fmt="(""#[REFERENCE]AQU, 0, 0, REF_X, REF_Y, CONST, SPEC_HEAD, MODEL_HEAD, ERROR"")")
      write (unit=iLU, &
             fmt="(""AQU"", 2("", "", i9), 5("", "", e16.8))" &
             ) 0, &
             0, &
             cIO_WorldCoords(io, aqu%cRefPoint), &
             aqu%rSolConst, &
             aqu%rRefHead, &
             rAQU_PotentialToHead(io, aqu, aqu%rCheck, aqu%cRefPoint)
    end if

    if (aqu%iNBdy > 0) then
      write (unit=iLU, &
             fmt="(""#[BOUNDARY]AQU, IBDY, FLAG, X1, Y1, X2, Y2, LENGTH, STRENGTH, SPEC_HEAD, SPEC_FLUX, " // &
             "GHB_DISTANCE, MODEL_HEAD, MODEL_FLUX"")")
      do iBdy = 1, aqu%iNBdy
        bdy => aqu%BdyElements(iBdy)
        write (unit=iLU, &
               fmt="(""BDY"", 2("", "", i9), 11("", "", e16.8))" &
               ) iBdy, &
               bdy%iBdyFlag, &
               cIO_WorldCoords(io, bdy%cZ1), &
               cIO_WorldCoords(io, bdy%cZ2), &
               bdy%rLength, bdy%rStrength, &
               bdy%rSpecHead, &
               bdy%rSpecFlux, &
               bdy%rGhbDistance, &
               rAQU_PotentialToHead(io, aqu, bdy%rCheckPot, rHALF*(bdy%cZ1+bdy%cZ2)), &
               rIO_WorldLength(io, bdy%rCheckFlux, bdy%cZ2-bdy%cZ1)
      end do
    end if

    call IN0_Inquiry(io, aqu%in0, iLU)

    return
  end subroutine AQU_Inquiry


  !! FUNCTION MODULE ROUTINES
  !! These routines allow AQU to behave as a ModAEM Function Module


  function cAQU_Potential(io, aqu, cZ) result(cOmega)
    !! complex function cAQU_Potential
    !!
    !! Computes the complex potential due to the reference aquifer at cZ
    !!
    !! Calling Sequence:
    !!    cOmega = cFDP_Potential(io, aqu, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!
    !! Return Value:
    !!           complex :: cOmega
    !!             The complex potential at the point cZ
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    if (aqu%lReference) then
      cOmega = aqu%rSolConst - conjg(aqu%cRefUniformFlow) * cZ
    else
      cOmega = aqu%rSolConst
    end if

    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cZ-bdy%cZC) / bdy%cZL
      cTemp = cILS_InfluenceP(cBigZ, bdy%cZL)
      cOmega = cOmega + cTemp(1, 1) * bdy%rStrength
    end do

  end function cAQU_Potential


  function cAQU_Discharge(io, aqu, cZ) result(cQ)
    !! complex function cAQU_Discharge
    !!
    !! Computes the complex discharge due to the reference aquifer at cZ
    !!
    !! Calling Sequence:
    !!    cOmega = cFDP_Discharge(io, aqu, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the discharge
    !!
    !! Return Value:
    !!           complex :: cW
    !!             The complex discharge at the point cZ in layer iL
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    if (aqu%lReference) then
      cQ = aqu%cRefUniformFlow
    else
      cQ = cZERO
    end if

    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cZ-bdy%cZC) / bdy%cZL
      cTemp = cILS_InfluenceW(cBigZ, bdy%cZL)
      cQ = cQ + conjg((cTemp(1, 1) * bdy%rStrength)/bdy%cZL)
    end do

    return
  end function cAQU_Discharge


  function rAQU_Flow(io, aqu, cPathZ) result(rFlow)
    !! real function cAQU_Flow
    !!
    !! Computes the integrated flow due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    rFlow = fAQU_Flow(io, aqu, cPathZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!
    !! Return Value:
    !!           real :: rFlow
    !!             The integrated flux across the path cPathZ
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy, iStat, i
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    ! The reference point(if any) and the uniform flow function satisfy
    ! Laplace equation, so we need to only look at the ends of the path
    if (aqu%lReference) then
      rFlow = aimag(cAQU_Potential(io, aqu, cPathZ(size(cPathZ))) - &
              cAQU_Potential(io, aqu, cPathZ(1)))
    else
      rFlow = rZERO
    end if


    allocate(cBigZ(size(cPathZ)), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'rAQU_Flow: Allocation failed')
    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cPathZ-bdy%cZC) / bdy%cZL
      do i = 1, size(cBigZ)-1
        cTemp = cILS_InfluenceF(cBigZ(i), cBigZ(i+1), bdy%cZL)
        rFlow = rFlow + real(cTemp(1, 1) * bdy%rStrength)
      end do
    end do
    deallocate(cBigZ)

    return
  end function rAQU_Flow


  function rAQU_Recharge(io, aqu, cZ) result(rG)
    !! complex function rAQU_Recharge
    !!
    !! Computes the recharge rate due to the aquifer
    !!
    !! Calling Sequence:
    !!    rG = rAQU_Recharge(io, aqu, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             The AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG
    ! [ LOCALS ]

    rG = rZERO

    return
  end function rAQU_Recharge


  function rAQU_Extraction(io, aqu) result(rQ)
    !! complex function rAQU_Extraction
    !!
    !! Computes the extraction due to the aquifer
    !!
    !! Calling Sequence:
    !!    rQ = rAQU_Extraction(io, aqu)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             The AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iBdy
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (.not. aqu%lReference), "AQU_Extraction: Invalid for problems with uniform flow")
    rQ = sum(aqu%BdyElements(1:aqu%iNBdy)%rCheckFlux)

    return
  end function rAQU_Extraction


  function rAQU_BranchCut(io, aqu, cPathZ) result(rBC)
    !! real function cAQU_Flow
    !!
    !! Computes the net branch cut.
    !!
    !! Calling Sequence:
    !!    rFlow = rAQU_BranchCut(io, aqu, cPathZ)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!
    !! Return Value:
    !!           real :: rFlow
    !!             The integrated flux across the path cPathZ
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBC
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: cBigZ
    complex(kind=AE_REAL), dimension(1, 1) :: cTemp
    integer(kind=AE_INT) :: iBdy, iStat, i
    type(AQU_BDYELEMENT), pointer :: bdy

    call IO_Assert(io, (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object")

    ! The reference point(if any) and the uniform flow function satisfy
    ! Laplace equation, so we need to only look at the ends of the path
    rBC = rZERO

    allocate(cBigZ(size(cPathZ)), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'rAQU_Flow: Allocation failed')
    do iBdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(iBdy)
      cBigZ = (cPathZ-bdy%cZC) / bdy%cZL
      do i = 1, size(cBigZ)-1
        cTemp = cILS_InfluenceB(cBigZ(i), cBigZ(i+1), bdy%cZL)
        rBC = rBC + aimag(cTemp(1, 1) * bdy%rStrength)
      end do
    end do
    deallocate(cBigZ)

    return
  end function rAQU_BranchCut


  function lAQU_CheckPoint(io, aqu, cZ, rTol, cZFix, rStrength, iElementType, iElementString, &
             iElementVertex, iElementFlag) result(lFound)
    !! function lAQU_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance of the
    !! perimeter of a pond.
    !! If rTol > 0, the tolerance is used in the test, otherwise the global constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rTol
    complex(kind=AE_REAL), intent(out) :: cZFix
    real(kind=AE_REAL), intent(out) :: rStrength
    integer(kind=AE_INT), intent(out) :: iElementType
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lFound
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(AQU_BDYELEMENT), pointer :: bdy
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    real(kind=AE_REAL) :: rCheckTol
    complex(kind=AE_REAL) :: cBZFix
    ! [ LOCAL PARAMETERS ]
    real(kind=AE_REAL), parameter :: rEND_FIX_POSITION = rTWO

    if (aqu%iNBdy > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a perimeter line sink
      allocate(cMapZ(aqu%iNBdy))
      cMapZ(:) = (cZ-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL

      lFound = .false.
      cZFix = cZERO
      do i = 1, aqu%iNBdy
        bdy => aqu%BdyElements(i)
        ! Check to see(1) we're within the tolerance of the extension of the line segment,
        ! then(2) we're along the(-1, +1) extension
        if (abs(aimag(cMapZ(i))) < rCheckTol) then
          ! We're near a segment. Are we at an end?
          if (abs(real(cMapZ(i))-cONE) < rCheckTol) then
            ! We're at the upper tip
            lFound = .true.
            cBZFix = cmplx(rONE-rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * bdy%cZL + bdy%cZC
            rStrength = bdy%rStrength
            iElementType = ELEM_AQU
            iElementString = kAQUBoundary
            iElementVertex = i
            iElementFlag = -1
            exit
          else if (abs(real(cMapZ(i))+cONE) < rCheckTol) then
            ! We're at the lower tip
            lFound = .true.
            cBZFix = cmplx(-rONE+rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * bdy%cZL + bdy%cZC
            rStrength = bdy%rStrength
            iElementType = ELEM_AQU
            iElementString = kAQUBoundary
            iElementVertex = i
            iElementFlag = -1
            exit
          end if
        end if
      end do
      deallocate(cMapZ)
    else
      lFound = .false.
      cZFix = cZERO
    end if

    return
  end function lAQU_CheckPoint


  function lAQU_CheckIntersection(io, aqu, cZ1, cZ2, cZInt, cZBefore, cZAfter, &
             iElementType, iElementString, iElementVertex, iElementFlag, cENormal) result(lFound)
    !! function lFDP_CheckIntersection
    !!
    !! Checks the specified line segment Z1-Z2 and returns .true. if there is an intersection.
    !! In addition, if the return value is .true., iElementType, iElementString, iElementVertex,
    !! and iElementFlag point to the feature that is intersected, cZInt is the coordinate of the
    !! intersection, and cZBefore and cZAfter are points just before and just after the
    !! intersection, normal to the dipole that is intersected.
    !!
    !! If rTol > 0, the tolerance is based on the constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), intent(out) :: cZInt, cZBefore, cZAfter
    integer(kind=AE_INT), intent(out) :: iElementType
    integer(kind=AE_INT), intent(out) :: iElementString
    integer(kind=AE_INT), intent(out) :: iElementVertex
    integer(kind=AE_INT), intent(out) :: iElementFlag
    complex(kind=AE_REAL), intent(out) :: cENormal
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lFound
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iLoc
    type(AQU_BDYELEMENT), pointer :: bdy
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1, cMapZ2, cMapZInt
    real(kind=AE_REAL), dimension(:), allocatable :: rDist
    complex(kind=AE_REAL) :: cMyZ1, cMyZ2, cZFix
    integer(kind=AE_INT) :: iEType, iEStr, iEVtx, iEFlg
    real(kind=AE_REAL) :: rStr
    real(kind=AE_REAL) :: rCheckTol
    ! [ LOCAL PARAMETERS ]
    ! This is the fraction of the dipole length that cZBefore and cZAfter will be away from
    ! the intersection point.
    real(kind=AE_REAL), parameter :: rPOINT_TOLERANCE_FACTOR = rTWO
    real(kind=AE_REAL), parameter :: rINTERSECTION_TOLERANCE = 1.0e-4_AE_REAL

    ! First, make sure we're not in tolerance
    if (lAQU_CheckPoint(io, aqu, cZ1, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ1 = cZFix
    else
      cMyZ1 = cZ1
    end if
    if (lAQU_CheckPoint(io, aqu, cZ2, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ2 = cZFix
    else
      cMyZ2 = cZ2
    end if

    if (aqu%iNBdy > 0) then
      ! Set the check tolerance
      ! Find the smallest distance to a dipole
      allocate(cMapZ1(aqu%iNBdy), cMapZ2(aqu%iNBdy), cMapZInt(aqu%iNBdy), rDist(aqu%iNBdy))
      cMapZ1(:) = (cMyZ1-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL
      cMapZ2(:) = (cMyZ2-aqu%BdyElements(1:aqu%iNBdy)%cZC) / aqu%BdyElements(1:aqu%iNBdy)%cZL
      ! Find the points of intersection with the line extensions
      cMapZInt = cmplx(rHUGE, rHUGE, AE_REAL)
      where (sign(rONE, aimag(cMapZ1)) /= sign(rONE, aimag(cMapZ2)))
        cMapZInt = cMapZ1 + abs(aimag(cMapZ1))/(abs(aimag(cMapZ1))+abs(aimag(cMapZ2))) * (cMapZ2-cMapZ1)
      end where
      ! Find the distances to line(s) that intersect
      rDist = rHUGE
      where (abs(real(cMapZInt)) <= rONE)
        rDist = abs(cMapZInt-cMapZ1)
      end where
      ! Now, find the CLOSEST intersection point
      iLoc = minloc(rDist, DIM = 1)
      if (rDist(iLoc) < rHUGE) then
        bdy => aqu%BdyElements(iLoc)
        ! Reverse-map to get the intersection
        cZInt = cMapZInt(iLoc) * bdy%cZL + bdy%cZC
        ! Compute ZBefore and ZAfter in the direction of the Z1-Z2 line segment, using the
        ! rINTERSECTION_TOLERANCE*abs(Z2-Z1) as the distance before and after the intersection
        cZBefore = cZInt - rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        cZAfter = cZInt + rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        ! If the before point is farther away from the intersection than the first point,
        ! Use the first point as the before point(and similarly for the after point)
        if (abs(cZ1-cZInt) < abs(cZBefore-cZInt)) cZBefore = cZ1
        if (abs(cZ2-cZInt) < abs(cZAfter-cZInt)) cZAfter = cZ2
        lFound = .true.
        iElementType = ELEM_AQU
        iElementString = 1
        iElementVertex = iLoc
        iElementFlag = -1
        cENormal = cI * bdy%cZL / abs(bdy%cZL)
      else
        lFound = .false.
        cZBefore = cZERO
        cZAfter = cZERO
        iElementType = -1
        iElementString = -1
        iElementVertex = -1
        iElementFlag = -1
        cENormal = cZERO
      end if
      deallocate(cMapZ1, cMapZ2, cMapZInt, rDist)
    else
      lFound = .false.
      cZBefore = cZERO
      cZAfter = cZERO
      iElementType = -1
      iElementString = -1
      iElementVertex = -1
      iElementFlag = -1
      cENormal = cZERO
    end if

    return
  end function lAQU_CheckIntersection

  !! UTILITY ROUTINES
  !! These routines provide computational aids for AEM models
  !! They make use of the IN0 equivalent calls


  function rAQU_ActiveArea(io, aqu, poly) result(rArea)
    !! real function rAQU_ComputeActiveArea
    !!
    !! Computes the area of 'poly' that lies within the active area of the
    !! aquifer 'aqu'.
    !!
    !! Calling Sequence:
    !!    rA = rAQU_ComputeActiveArea(io, aqu, poly)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex, dimension(:) :: poly
    !!             The polygon to be intersected with the aquifer
    !!
    !! Return Value:
    !!           real :: rA
    !!             The intersecting area
    !! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    type(IO_STATUS), pointer :: io
    !! [ LOCALS ]
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rArea

    rArea = PGN_Area(poly)

    return
  end function rAQU_ActiveArea


  function rAQU_HeadToPotential(io, aqu, rHead, cZ) result(rPot)
    !! real function rAQU_HeadToPotential
    !!
    !! Converts head to a discharge potential based on the in0ifer properties
    !! at cZ. Returns the(real) potential.
    !!
    !! Calling Sequence:
    !!    rPot = rAQU_HeadToPotential(io, aqu, rHead, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    real :: rHead
    !!             The head value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the aquifer where the conversion is to take place
    !!
    !! Return Value:
    !!           real :: rPot
    !!             The discharge potential corresponding to the head 'rHead'
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rHead
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPot

    rPot = rIN0_HeadToPotential(io, aqu%in0, rHead, cZ)

    return
  end function rAQU_HeadToPotential


  function rAQU_PotentialToHead(io, aqu, rPot, cZ, lTest) result(rHead)
    !! real function rAQU_HeadToPotential
    !!
    !! Converts head to a discharge potential based on the in0ifer properties
    !! at cZ. Returns the(real) potential.
    !!
    !! Calling Sequence:
    !!    rPot = rAQU_HeadToPotential(io, in0, rHead, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    real :: rHead
    !!             The head value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           real :: rHead
    !!             The head corresponding to the discharge potential 'rPot'
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lTest
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHead

    if (.not. present(lTest)) then
      rHead = rIN0_PotentialToHead(io, aqu%in0, rPot, cZ)
    else
      rHead = rIN0_PotentialToHead(io, aqu%in0, rPot, cZ, lTest)
    end if

    return
  end function rAQU_PotentialToHead


  function cAQU_DischargeToVelocity(io, aqu, cDischarge, cZ, rPot) result(cVelocity)
    !! function cAQU_DischargeToVelocity
    !!
    !! Converts Discharge to velocity, using the potential and the location cZ
    !! to compute saturated thickness
    !!
    !! Calling Sequence:
    !!    cV = cAQU_DischargeToVelocity(io, in0, rDischarge, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    real :: rDischarge
    !!             The discharge value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           complex :: cV
    !!             The velocity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cDischarge
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cVelocity

    cVelocity = cIN0_DischargeToVelocity(io, aqu%in0, cDischarge, cZ, rPot)

    return
  end function cAQU_DischargeToVelocity


  function lAQU_IsConfined(io, aqu, cZ, rPot) result(lConfined)
    !! function lAQU_IsConfined`
    !!
    !! Checks to see if the flow condition at some point is confined
    !!
    !! Calling Sequence:
    !!    if (lAQU_IsConfined(io, aqu, cZ, rPot)) then ...
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the check is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           complex :: cV
    !!             The velocity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lConfined

    lConfined = lIN0_IsConfined(io, aqu%in0, cZ, rPot)

    return
  end function lAQU_IsConfined


  function rAQU_Base(io, aqu, cZ) result(rBase)
    !! function rAQU_HydCond
    !!
    !! Computes the base elevation at the point cZ
    !!
    !! Calling Sequence:
    !!    rB = rAQU_Base(io, aqu, cZ)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQu_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           complex :: rBase
    !!             The base elevation
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBase

    rBase = rIN0_Base(io, aqu%in0, cZ)

    return
  end function rAQU_Base


  function rAQU_InterfaceElevation(io, aqu, cZ, rPot) result(rIfcElev)
    !! real function rAQU_InterfaceElevation
    !!
    !! Converts head to a discharge potential based on the in0ifer properties
    !! at cZ. Returns the(real) potential.
    !!
    !! Calling Sequence:
    !!    rPot = rAQU_InterfaceElevation(io, in0, rPot, cZ)
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rIfcElev

    call IO_MessageText(io, "Interface calculations are not yet available")
    rIfcElev = rAQU_Base(io, aqu, cZ)

    return
  end function rAQU_InterfaceElevation


  function rAQU_SatdThickness(io, aqu, cZ, rPot) result(rH)
    !! function cAQU_SatdThickness
    !!
    !! Computes the saturated thickness at point cZ in layer iL where the potential is rPot
    !!
    !! Calling Sequence:
    !!    cV = cAQU_SatdThickness(io, in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           complex :: cV
    !!             The velocity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH

    rH = rIN0_SatdThickness(io, aqu%in0, cZ, rPot)

    return
  end function rAQU_SatdThickness


  function rAQU_DefaultPotential(io, aqu, cZ) result(rP)
    !! function cAQU_SatdThickness
    !!
    !! Returns the potential at the aquifer top cZ
    !!
    !! Calling Sequence:
    !!    cV = cAQU_DefaultPotential(io, in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQU_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rP

    rP = rIN0_DefaultPotential(io, aqu%in0, cZ)

    return
  end function rAQU_DefaultPotential


  function rAQU_Transmissivity(io, aqu, cZ, rPot) result(rT)
    !! function rAQU_Transmissivity
    !!
    !! Computes the transmissivity at the point cZ where the potential is rPot
    !!
    !! Calling Sequence:
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQu_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           complex :: cV
    !!             The velocity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rT

    rT = rIN0_Transmissivity(io, aqu%in0, cZ, rPot)

    return
  end function rAQU_Transmissivity


  function lAQU_CheckActive(io, aqu, cZ) result(lActive)
    !! function lAQU_CheckActive
    !!
    !! Returns .true. if cZ is in the active area of the model
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lActive

    if (aqu%iNPerimeter==0) then
      lActive = .true.
    else
      lActive = (PGN_CheckPoint(aqu%cPerimeter(1:aqu%iNPerimeter), cZ) > 0)
    end if

    return
  end function lAQU_CheckActive


  function rAQU_HydCond(io, aqu, cZ) result(rHydCond)
    !! function rAQU_HydCond
    !!
    !! Computes the hydraulic conductivity at the point cZ
    !!
    !! Calling Sequence:
    !!    cV = cAQU_HydCond(io, aqu, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             AQu_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           complex :: rHydCond
    !!             The hydraulic conductivity
    !!
    ! [ ARGUMENTS ]
    type(AQU_COLLECTION), pointer :: aqu
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHydCond

    rHydCond = rIN0_HydCond(io, aqu%in0, cZ)

    return
  end function rAQU_HydCond


  subroutine AQU_Save(io, aqu, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: ibdy
    type(AQU_BDYELEMENT), pointer :: bdy

    ! Write the solution constant
    if (mode == IO_MODE_BINARY) then
      write (unit=LU_SCRATCH) ELEM_AQU, 0, 1, 1, aqu%rSolConst
    else
      write (unit=LU_SCRATCH, fmt=*) "AQU", 0, 1, 1, aqu%rSolConst
    end if
    ! Output records will be of the form ELEM_AQU, IBDY, 1, 1, SIGMA
    do ibdy = 1, aqu%iNBdy
      bdy => aqu%BdyElements(ibdy)
      if (mode == IO_MODE_BINARY) then
        write (unit=LU_SCRATCH) ELEM_AQU, ibdy, 1, 1, bdy%rStrength
      else
        write (unit=LU_SCRATCH, fmt=*) "AQU", ibdy, 1, 1, bdy%rStrength
      end if
    end do
    ! Save the inhomogeneities!
    call IN0_Save(io, aqu%in0, mode)

    return
  end subroutine AQU_Save


  subroutine AQU_Load(io, aqu, fdp, mode)
    !! Loads the AQU records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AQU_COLLECTION), pointer :: aqu
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, ibdy, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(AQU_BDYELEMENT), pointer :: bdy

    ! Scans the entire precondition file for the AQU data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, ibdy, ivtx, iflg, rstrength
        if (imodule /= ELEM_AQU) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, ibdy, ivtx, iflg, rstrength
        if (uppercase (trim(smodule)) /= "AQU") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      if (ibdy == 0) then
        aqu%rSolConst = rstrength
        call IO_Assert(io, ivtx == 1, "AQU vertex not found")
        call IO_Assert(io, iflg == 1, "AQU bdyength index not found")
      else
        call IO_Assert(io, ibdy > 0 .and. ibdy <= aqu%iNBdy, "AQU bdyl not found")
        bdy => aqu%BdyElements(ibdy)
        call IO_Assert(io, ivtx == 1, "AQU vertex not found")
        call IO_Assert(io, iflg == 1, "AQU bdyength index not found")
        bdy%rStrength = rstrength
      end if
    end do

    ! Load the inhomogeneities
    call IN0_Load(io, aqu%in0, fdp, mode)

    ! Now, populate the internal data bdyuctures
    call AQU_Update(io, aqu, fdp)

    return
  end subroutine AQU_Load

end module m_aqu
