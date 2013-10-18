module m_in0

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

  !! module m_in0
  !!
  !! Element module for 2-D inhomogeneities
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for inhomogeneity domains
  !! (no common boundaries at this time). The Domains are condomucted using
  !! line-doublets(implemented as line-dipoles of imaginary strength).


  ! TODO!!
  ! Replace all the polygonal geometry routines in IN0 with those in u_polygon,
  ! adding routines to u_polygon where necessary.  Current IN0 uses u_polygon
  ! only to ensure that domains are positively-oriented.

  use u_constants
  use u_io
  use u_matrix
  use u_polygon
  use f_dipole

  implicit none

  public

  private :: lIN0_PointInsideDomain, &
             lIN0_LineIntersectsDomain, &
             lIN0_DomainInsideDomain

  type :: IN0_VERTEX
    !! type IN0_VERTEX
    !!
    !! Type that holds information for one vertex along the edge of an inhomogeneity    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rStrength(3)
    !!     The doublet strength at the vertex
    !!   real :: rStrength(1)
    !!     The doublet strength at the center of the next segment
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!
    complex(kind=AE_REAL) :: cZ
    integer(kind=AE_INT) :: iFDPIndex
    real(kind=AE_REAL), dimension(3) :: rStrength
    real(kind=AE_REAL), dimension(3) :: rCheckPot
    real(kind=AE_REAL), dimension(3) :: rLeftH
    real(kind=AE_REAL), dimension(3) :: rRightH
    real(kind=AE_REAL), dimension(3) :: rLeftT
    real(kind=AE_REAL), dimension(3) :: rRightT
    real(kind=AE_REAL), dimension(3) :: rError
    complex(kind=AE_REAL), dimension(3) :: cCPZ
    complex(kind=AE_REAL), dimension(3) :: cCPZOpp
  end type IN0_VERTEX

  type :: IN0_DOMAIN
    !! type IN0_DOMAIN
    !!
    !! Type that holds information for one inhomogeneity
    !!
    !! Members:
    !!   type(IN0_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of IN0_VERTEX objects
    !!   integer :: iInsideDomain
    !!     Number of the IN0_DOMAIN that this domain is inside
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iID
    !!     The ID number for the string
    !!   real :: rBase
    !!     Base elevation
    !!   real :: rThickness
    !!     Thickness of the in0ifer(set to a large value for unconfined flow)
    !!   real :: rHydCond
    !!     Hydraulic conductivity
    !!   real :: rPorosity
    !!     Porosity
    !!   real :: rConfPot
    !!     Potential at the top of the in0ifer
    !!
    complex(kind=AE_REAL), dimension(:), pointer :: cZ
    integer(kind=AE_INT) :: iInsideDomain
    integer(kind=AE_INT) :: iOutsideDomain
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iID
    real(kind=AE_REAL) :: rBase
    real(kind=AE_REAL) :: rThickness
    real(kind=AE_REAL) :: rHydCond
    real(kind=AE_REAL) :: rPorosity
    real(kind=AE_REAL) :: rTopPot
    real(kind=AE_REAL) :: rAvgHead
  end type IN0_DOMAIN

  type :: IN0_STRING
    !! type IN0_STRING
    !!
    !! Type that holds information for one inhomogeneity string.
    !! Strings are defined separately from the domains.
    !!
    !! Members:
    !!   type(IN0_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of IN0_VERTEX objects
    !!   integer :: iLeftDomain
    !!     Number of the IN0_DOMAIN that is to the left of the string
    !!   integer :: iRightDomain
    !!     Number of the IN0_DOMAIN that is to the right of the string
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   logical :: lClosed
    !!     .true. if the string closes on itself
    !!   integer :: iID
    !!     The ID number for the string
    !!
    type(IN0_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iLeftID
    integer(kind=AE_INT) :: iRightID
    real(kind=AE_REAL) :: rLeftB
    real(kind=AE_REAL) :: rRightB
    integer(kind=AE_INT) :: iNPts
    logical :: lClosed
    integer(kind=AE_INT) :: iID
  end type IN0_STRING

  type :: IN0_COLLECTION
    !! type IN0_COLLECTION
    !!
    !! Type that holds information for all IN0 elements in a layer
    !!
    !! Members:
    !!   type(IN0_DOMAIN), dimension(:), pointer :: Domains
    !!     A vector of IN0_DOMAIN objects
    !!   integer :: iNDom
    !!     The number of Domains actually in use
    !!
    type(IN0_DOMAIN), dimension(:), pointer :: Domains
    integer(kind=AE_INT) :: iNDom
    type(IN0_STRING), dimension(:), pointer :: Strings
    integer(kind=AE_INT) :: iNStr
    integer(kind=AE_INT) :: iRegenerate
    logical :: lInitialIteration
    ! Iterator information
    integer(kind=AE_INT) :: iIterStr
    integer(kind=AE_INT) :: iIterVtx
    integer(kind=AE_INT) :: iIterFlag
    logical :: lPrecondition
  end type IN0_COLLECTION

  ! Matrix generator element flags
  integer(kind=AE_INT), private, parameter :: kIN0_Vertex = 1
  integer(kind=AE_INT), private, parameter :: kIN0_Center = 2
  integer(kind=AE_INT), private, parameter :: kIN0_Vertex2 = 3
  ! How far do I move the control points off the line segments?
  real(kind=AE_REAL), private, parameter :: IN0_NORMAL_OFFSET = -1.0e-7_AE_REAL
  real(kind=AE_REAL), private, parameter :: IN0_OPPOSITE_OFFSET = 1.0e-7_AE_REAL

  real(kind=AE_REAL) :: MOVEFACTOR = 1.0001_AE_REAL


contains


  function IN0_Create(io, iNDom, iNStr, rBase, rThickness, rHydCond, rPorosity, rAvgHead) result(in0)
    !! function IN0_Create
    !!
    !! Creates a new IN0_COLLECTION object, and builds the infinite in0ifer domain
    !!
    !! Calling Sequence:
    !!    IN0 => IN0_Create(iNDom, iNStr)
    !!
    !! Arguments:
    !!    (in)    integer :: iNDom
    !!              The number of domains to make space for; this routine will
    !!              add one domain for the "infinite" outside domain
    !!    (in)    integer :: iNStr
    !!              The number of strings to make space for; this routine will
    !!              add one domain for the "infinite" outside domain
    !!    (in)    real :: rBase
    !!              The base elevation for the infinite domain.  NOTE: Module
    !!              IN0 does not support base elevation inhomogeneities; this is
    !!              a placeholder for future expansion.
    !!    (in)    real :: rThickness
    !!              The thickness of the infinite in0ifer
    !!    (in)    real :: rHydCond
    !!              The hydraulic conductivity of the infinite in0ifer
    !!    (in)    real :: rPorosity
    !!              The porosity of the infinite in0ifer
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT), intent(in) :: iNDom
    integer(kind=AE_INT), intent(in) :: iNStr
    real(kind=AE_REAL), intent(in) :: rBase
    real(kind=AE_REAL), intent(in) :: rThickness
    real(kind=AE_REAL), intent(in) :: rHydCond
    real(kind=AE_REAL), intent(in) :: rPorosity
    real(kind=AE_REAL), intent(in) :: rAvgHead
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(IN0_COLLECTION), pointer :: in0

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(IN0_DOMAIN), pointer :: dom

    allocate(in0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "IN0_Create: allocation failed")
    allocate(in0%Domains(iNDom), in0%Strings(iNStr), stat = iStat)
    call IO_Assert(io, (iStat == 0), "IN0_Create: allocation failed")
    in0%Domains%iID = -1
    in0%Strings%iID = -1

    in0%iNDom = 1
    in0%iNStr = 0
    in0%iRegenerate = 1
    in0%lInitialIteration = .true.
    in0%lPrecondition = .true.

    ! Build the infinite aquifer
    dom => in0%Domains(1)
    dom%rBase = rBase
    dom%rThickness = rThickness
    dom%rHydCond = rHydCond
    dom%rPorosity = rPorosity
    dom%rAvgHead = rAvgHead
    dom%iInsideDomain = 0
    dom%iInsideDomain = 0
    dom%iNPts = 0
    dom%iID = 0

    return
  end function IN0_Create


  subroutine IN0_New(io, in0, cZ, iNPts, rBase, rThickness, rHydCond, rPorosity)
    !! function IN0_New
    !!
    !! Adds a new IN0_DOMAIN object to the IN0_COLLECTION 'IN0'
    !!
    !! Calling Sequence:
    !!    call IN0_New(in0, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(IN0_COLLECTION), pointer :: in0
    !!              The IN0_COLLECTION object to be used
    !!    (in)    complex :: cZ(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!    (in)    real :: rBase
    !!              Base elevation
    !!    (in)    real :: rThickness
    !!              Aquifer thickness
    !!    (in)    real :: rHydCond
    !!              Hydraulic conductivity
    !!    (in)    real :: rPorosity
    !!              Effective porosity
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), dimension(:) :: cZ
    integer(kind=AE_INT), intent(in) :: iNPts
    real(kind=AE_REAL), intent(in) :: rBase
    real(kind=AE_REAL), intent(in) :: rThickness
    real(kind=AE_REAL), intent(in) :: rHydCond
    real(kind=AE_REAL), intent(in) :: rPorosity
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iDom
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_New: IN0_Create has not been called")
    end if

    call IO_Assert(io, (in0%iNDom < size(in0%Domains)), &
         "IN0_New: Space exhausted")
    call IO_Assert(io, (iNPts <= size(dom%cZ)), &
         "IN0_New: Size of provided vertices is inconsistent")

    in0%iNDom = in0%iNDom + 1
    dom => in0%Domains(in0%iNDom)
    allocate(dom%cZ(iNPts), stat = iStat)
    call IO_Assert(io, (iStat == 0), "IN0_New: Allocation failed")
    dom%cZ = cZ(1:iNPts)
    dom%iNPts = iNPts
    dom%rBase = rBase
    dom%rThickness = rThickness
    dom%rHydCond = rHydCond
    dom%rPorosity = rPorosity

    ! This version of IN0 does not support nested inhomogeneity domains
    ! Check to ensure that no nesting or overlapping of domains has occurred
    ! nesting data structures
    do iDom = 2, in0%iNDom-1
      call IO_Assert(io, (.not. lIN0_DomainOverlapsDomain(io, in0, iDom, in0%iNDom)), &
           "IN0_New: New domain overlaps another domain")
      !    call IO_Assert(io, (.not. lIN0_DomainInsideDomain(io, in0, iDom, in0%iNDom)), &
          !                   "IN0_New: New domain is inside another domain")
    end do

    return
  end subroutine IN0_New


  subroutine IN0_SetPrecondition(io, in0, lPre)
    ! Sets the precondition flag. If it is set, it is presumed that this is the
    ! first iteration, and the saturated thickness / transmissivity are based
    ! on the "avg-head" setting.
    type(IN0_COLLECTION), pointer :: in0
    logical, intent(in) :: lPre
    type(IO_STATUS), pointer :: io

    in0%lPrecondition = lPre

    return
  end subroutine IN0_SetPrecondition


  subroutine IN0_PreSolve(io, in0)
    !! subroutine IN0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call IN0_PreSolve(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iDom, iDom2, iStr, iVtx, iNChanges
    character(len=255) :: sBuf
    type(IN0_DOMAIN), pointer :: dom, dom2, outside
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx
    complex(kind=AE_REAL), dimension(:), allocatable :: cTempZ
    real(kind=AE_REAL) :: rH, rT

    ! Compute the necessary constants
    do iDom = 1, in0%iNDom
      dom => in0%Domains(iDom)
      if (iDom > 1) then
        call PGN_MakePositivelyOriented(dom%cZ(1:dom%iNPts))
      end if
    end do

    ! Nest the domains, if necessary, using an iterative strategy
    print *, 'Checking for nested inhomogeneities...'
    ! Handle the potentials at bottom and top of all aquifers(ugh)
    in0%Domains(:)%rTopPot = rHALF * in0%Domains(:)%rHydCond * in0%Domains(:)%rThickness**2
    do
      iNChanges = 0
      do iDom = 2, in0%iNDom
        dom => in0%Domains(iDom)
        do iDom2 = 2, in0%iNDom
          dom2 => in0%Domains(iDom2)
          if (iDom /= iDom2) then
            ! Set iDom's domain to be iDom2 if (1) iDom2 contains iDom's outside
            ! domain and(2) iDom is inside iDom2
            if (lIN0_DomainInsideDomain(io, in0, iIN0_DomainIDIndex(io, in0, dom%iOutsideDomain), iDom2) .and. &
                lIN0_DomainInsideDomain(io, in0, iDom2, iDom)) then
              dom%iOutsideDomain = dom2%iID
              iNChanges = iNChanges+1
            end if
          end if
        end do

        outside => IN0_FindDomainID(io, in0, dom%iOutsideDomain)
      end do
      if (iNChanges == 0) exit
    end do

    ! If no strings were explicitly defined, create strings for the domains
    if (in0%iNStr == 0) then
      print *, 'Copying IN0 domains into strings...'
      call IO_Assert(io, (size(in0%Strings) >= in0%iNDom-1), &
           "IN0_PreSolve: Insufficient number of strings is available")
      do iDom = 2, in0%iNDom
        in0%iNStr = in0%iNStr+1
        dom => in0%Domains(iDom)
        str => in0%Strings(in0%iNStr)
        allocate(str%Vertices(dom%iNPts), stat = iStat)
        call IO_Assert(io, (iStat == 0), &
             "IN0_PreSolve: Space exhausted")
        str%iNPts = dom%iNPts
        do iVtx = 1, dom%iNPts
          vtx => str%Vertices(iVtx)
          vtx%cZ = dom%cZ(iVtx)
          vtx%rStrength = rZERO
          vtx%rCheckPot = rZERO
          vtx%iFDPIndex = -1
        end do
        str%iLeftID = dom%iInsideDomain
        str%iRightID = dom%iOutsideDomain
        str%lClosed = .true.
        str%iID = dom%iID
      end do
      !  else
      !    ! Make sure closed strings are positively oriented
      !    do iStr = 1, in0%iNStr
      !      str => in0%Strings(iStr)
      !      if (str%lClosed) then
      !        allocate(cTempZ(str%iNPts))
      !        cTempZ = str%Vertices(1:str%iNPts)%cZ
      !        call PGN_MakePositivelyOriented(io, cTempZ)
      !        str%Vertices(1:str%iNPts)%cZ = cTempZ
      !        deallocate(cTempZ)
      !      end if
      !    end do
    end if

    ! Set the RightB and LeftB base elevations for all strings
    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      dom => IN0_FindDomainID(io, in0, str%iLeftID)
      str%rLeftB = dom%rBase
      if (dom%rAvgHead > (dom%rBase+dom%rThickness)) then
        rH = dom%rThickness
      else
        rH = dom%rAvgHead - dom%rBase
      end if
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        vtx%rLeftH = rH
        vtx%rLeftT = rH * dom%rHydCond
      end do
      dom => IN0_FindDomainID(io, in0, str%iRightID)
      str%rRightB = dom%rBase
      if (dom%rAvgHead > (dom%rBase+dom%rThickness)) then
        rH = dom%rThickness
      else
        rH = dom%rAvgHead - dom%rBase
      end if
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        vtx%rRightH = rH
        vtx%rRightT = rH * dom%rHydCond
      end do
    end do


    return
  end subroutine IN0_PreSolve


  function iIN0_GetInfo(io, in0, iOption, iIteration) result(iValue)
    !! function IN0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iIN0_GetInfo(aqu, iOption)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             AQU_COLLECTION to be used
    !!   (out)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr
    type(IN0_STRING), pointer :: str

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_GetInfo: IN0_Create has not been called")
    end if

    ! Either count the domains(if in0%iNStr == 0) or the Strings
    iValue = 0
    select case (iOption)
      case (SIZE_FDP)
        do iStr = 1, in0%iNStr
          str => in0%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + str%iNPts
          else
            iValue = iValue + str%iNPts-1
          end if
        end do
      case (SIZE_EQUATIONS)
        do iStr = 1, in0%iNStr
          str => in0%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + 2*str%iNPts
          else
            iValue = iValue + 2*str%iNPts-1
          end if
        end do
      case (SIZE_UNKNOWNS)
        do iStr = 1, in0%iNStr
          str => in0%Strings(iStr)
          if (str%lClosed) then
            iValue = iValue + 2*str%iNPts
          else
            iValue = iValue + 2*str%iNPts-1
          end if
        end do
      case (INFO_REGENERATE)
        ! Forece a regen on the first two iterations
        if (iIteration < 2) then
          iValue = 1
        else
          iValue = in0%iRegenerate
        end if
      case default
        iValue = 0
    end select

    return
  end function iIN0_GetInfo


  subroutine IN0_SetupFunctions(io, in0, fdp)
    !! subroutine IN0_SetupFunctions
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call IN0_SetupFunctions(in0, fdp)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iDom, iStr, iVtx, i, ii
    real(kind=AE_REAL) :: rSum
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3
    complex(kind=AE_REAL), dimension(6) :: cCPResult1, cCPResult2
    complex(kind=AE_REAL), dimension(3) :: cCPVtx
    character(len=255) :: sBuf
    integer(kind=AE_INT) :: irv
    type(IN0_DOMAIN), pointer :: dom
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx
    type(IN0_VERTEX) :: temp_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_Setup: IN0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "IN0_Setup: Illegal FDP_COLLECTION object")
    end if

    ! Build dipoles for all segments
    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        cZ1 = vtx%cZ
        if (iVtx < str%iNPts) then
          cZ2 = str%Vertices(iVtx+1)%cZ
          call FDP_New(io, fdp, cZ1, cZ2, (/cZERO, cZERO, cZERO/), ELEM_IN0, iStr, iVtx, -1, vtx%iFDPIndex)
        else
          if (str%lClosed) then
            cZ2 = str%Vertices(1)%cZ
            call FDP_New(io, fdp, cZ1, cZ2, (/cZERO, cZERO, cZERO/), ELEM_IN0, iStr, iVtx, -1, vtx%iFDPIndex)
          end if
        end if
      end do
    end do

    return
  end subroutine IN0_SetupFunctions


  subroutine IN0_SetupMatrix(io, in0, mat)
    !! subroutine IN0_SetupMatrix
    !!
    !! This routine sets up the matrix entries for the module
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call IN0_SetupiMatrix(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!   (in)    type(MAT_MATRIX), pointer :: mat
    !!             MAT_MATRIX object to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx
    complex(kind=AE_REAL) :: cZ1, cZ2, cZ3
    complex(kind=AE_REAL), dimension(6) :: cCPResult1, cCPResult2
    complex(kind=AE_REAL), dimension(3) :: cCPVtx
    character(len=255) :: sBuf
    integer(kind=AE_INT) :: irv, iEQ
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_Setup: IN0_Create has not been called")
      call IO_Assert(io, (associated(mat)), &
           "IN0_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix entries for all segments
    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      ! Set up the unknown variables
      ! Vertex entry -- all vertices except the first and last
      ! No dipole strength at either end of the string
      do iVtx = 1, str%iNPts
        vtx => str%Vertices(iVtx)
        ! Make a doublet vertex variable
        call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kIN0_Vertex)
        call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kIN0_Center)
        if (.not. str%lClosed .and. iVtx == str%iNPts-1) then
          ! Create the 'tip' variable, only if not closed, then exit the loop!
          call MAT_CreateVariable(io, mat, ELEM_IN0, iStr, iVtx, kIN0_Vertex2)
          exit
        end if
      end do

      ! Set up control point sets and equations
      if (str%lClosed) then
        ! Closed boundary
        do iVtx = 1, str%iNPts
          vtx => str%Vertices(iVtx)
          ! Now, create the equation entries...
          ! Two equations for each segment of the perimeter path
          cZ1 = vtx%cZ
          if (iVtx < str%iNPts) then
            cZ2 = str%Vertices(iVtx+1)%cZ
          else
            cZ2 = str%Vertices(1)%cZ
          end if
          ! Compute control point locations
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, IN0_NORMAL_OFFSET, &
               (/0.001_AE_REAL, 0.501_AE_REAL/))
          ! Vertex entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(1)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kIN0_Vertex, cZ2-cZ1, rZERO)
          ! Center entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(2)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kIN0_Center, cZ2-cZ1, rZERO)
          vtx%cCPZ(1:2) = cCPResult1(1:2)
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, IN0_NORMAL_OFFSET, &
               (/0.001_AE_REAL, 0.501_AE_REAL/))
          vtx%cCPZOpp(1:2) = cCPResult1(1:2)
        end do
      else
        ! Open string
        do iVtx = 1, str%iNPts-1
          vtx => str%Vertices(iVtx)
          ! Now, create the equation entries...
          cZ1 = vtx%cZ
          cZ2 = str%Vertices(iVtx+1)%cZ
          ! Compute control point locations
          call MAT_ComputeControlPoints(io, cZ1, cZ2, 1, cCPResult1, IN0_NORMAL_OFFSET, &
               (/0.001_AE_REAL, rHALF, rONE-0.001_AE_REAL/))
          ! Vertex entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(1)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kIN0_Vertex, cZ2-cZ1, rZERO)
          vtx%cCPZ(kIN0_Vertex) = cCPResult1(1)
          ! Center entry
          iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(2)/), EQN_INHO, ELEM_IN0, &
                                   iStr, iVtx, kIN0_Center, cZ2-cZ1, rZERO)
          vtx%cCPZ(kIN0_Center) = cCPResult1(2)
          if (iVtx == str%iNPts-1) then
            ! End of an open string
            iEQ = MAT_CreateEquation(io, mat, (/cCPResult1(3)/), EQN_INHO, ELEM_IN0, &
                                     iStr, iVtx, kIN0_Vertex2, cZ2-cZ1, rZERO)
            vtx%cCPZ(kIN0_Vertex2) = cCPResult1(3)
          end if
        end do
      end if

    end do
  end subroutine IN0_SetupMatrix


  function iIN0_Prepare(io, in0, iIteration) result(iChanges)
    !! subroutine IN0_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_in0
    !!
    !! Calling Sequence:
    !!    call IN0_Setup(in0, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iIN0_Prepare


  function rIN0_GetCoefficientMultiplier(io, in0, iElementString, iElementVertex, iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier
    ! [ LOCALS ]
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx

    str => in0%Strings(iElementString)
    vtx => str%Vertices(iElementVertex)

    rMultiplier = vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)

    return
  end function rIN0_GetCoefficientMultiplier


  subroutine IN0_ComputeCoefficients(io, in0, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbFactor, &
               iIteration, rMultiplier, rARow)
    !! subroutine IN0_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the IN0
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call IN0_ComputeCoefficients(in0, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
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
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    real(kind=AE_REAL), intent(in) :: rGhbFactor
    integer(kind=AE_INT), intent(in) :: iIteration
    real(kind=AE_REAL), intent(in) :: rMultiplier
    real(kind=AE_REAL), dimension(:), intent(out) :: rARow
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iSkipCol, iVtx, iDP1, iNDP, iWhich, irv, iThisDP, i, iBaseCol, iStr
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDPW
    complex(kind=AE_REAL), dimension(1, 3, 1) :: cDPJ
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_ComputeCoefficients: IN0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "IN0_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iSkipCol = 0
    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      ! Assume: the IN0_Setup routine creates consecutive dipole entries
      iDP1 = str%Vertices(1)%iFDPIndex
      if (str%lClosed) then
        iNDP = str%iNPts
      else
        iNDP = str%iNPts-1
      end if
      allocate(cDPF(0:iNDP+1, 3, 1), cDPW(0:iNDP+1, 3, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "IN0_ComputeCoefficients: Allocation failed")

      ! Now, compute the matrix coefficients
      ! Note that xxxGetInfluence returns zero for elements 0 and iNDP+1 of the cDPFx vectors
      cDPF = cZERO
      ! Get the appropriate influence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_BDYGHB)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, (/rHALF*sum(cPathZ)/), cOrientation, cDPF(1:iNDP, :, :))
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPW(1:iNDP, :, :))
          cDPF = cDPF + rGhbFactor*cDPW
        case (EQN_FLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_INHO)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_DISCHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_W, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_RECHARGE)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_G, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_CONTINUITY)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Q, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_POTENTIALDIFF)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_D, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_TOTALFLOW)
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_Z, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
      end select

      do iVtx = 1, iNDP
        this_vtx => str%Vertices(iVtx)
        if (iVtx < iNDP) then
          next_vtx => str%Vertices(iVtx+1)
        else
          next_vtx => str%Vertices(1)
        end if

        iBaseCol = 2*iVtx-1 + iSkipCol

        ! Compute the contributions for all line-doublets
        ! Vertex 1 contribution
        rARow(iBaseCol) = rARow(iBaseCol)   - rMultiplier * aimag(cDPF(iVtx, 1, 1))
        ! Center contribution
        rARow(iBaseCol+1) = rARow(iBaseCol+1) - rMultiplier * aimag(cDPF(iVtx, 2, 1))
        ! Vertex 2 contribution
        if (str%lClosed) then
          ! For closed strings...
          if (iVtx < iNDP) then
            rARow(iBaseCol+2) = rARow(iBaseCol+2) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
          else
            rARow(1+iSkipCol) = rARow(1+iSkipCol) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
          end if
        else
          ! For open strings, store the last term
          rARow(iBaseCol+2) = rARow(iBaseCol+2) - rMultiplier * aimag(cDPF(iVtx, 3, 1))
        end if

        ! Do I compute the additional term for the inhomogeneity?
        if ((iEqType == EQN_INHO) .and. &
            (iElementType == ELEM_IN0) .and. &
            (iElementString == iStr) .and. &
            (iElementVertex == iVtx)) then
          iThisDP = str%Vertices(iVtx)%iFDPIndex
          call FDP_GetInfluence_IDP(io, fdp, INFLUENCE_J, iThisDP, 1, cPathZ(1:1), cOrientation, cDPJ)
          ! Vertex 1 contribution
          rARow(iBaseCol) = rARow(iBaseCol)   - this_vtx%rRightT(1)*aimag(cDPJ(1, 1, 1))
          ! Center contribution
          rARow(iBaseCol+1) = rARow(iBaseCol+1) - this_vtx%rRightT(2)*aimag(cDPJ(1, 2, 1))
          ! Vertex 2 contribution
          if (str%lClosed) then
            ! For closed strings...
            if (iVtx < iNDP) then
              rARow(iBaseCol+2) = rARow(iBaseCol+2) - next_vtx%rRightT(1)*aimag(cDPJ(1, 3, 1))
            else
              rARow(1+iSkipCol) = rARow(1+iSkipCol) - next_vtx%rRightT(1)*aimag(cDPJ(1, 3, 1))
            end if
          else
            ! For open strings, update the last term
            rARow(iBaseCol+2) = rARow(iBaseCol+2) - this_vtx%rRightT(3)*aimag(cDPJ(1, 3, 1))
          end if
        end if
      end do

      ! Now, advance through the buffer
      if (str%lClosed) then
        iSkipCol = iSkipCol + 2*str%iNPts
      else
        iSkipCol = iSkipCol + 2*str%iNPts - 1
      end if

      ! No memory leaks, please!
      deallocate(cDPF, cDPW)
    end do

    return
  end subroutine IN0_ComputeCoefficients


  function rIN0_ComputeRHS(io, in0, fdp, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rIN0_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rIN0_ComputeRHS(in0, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
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
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRHS
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rCorr
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx, next_vtx

    ! Assertions...
    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_ComputeRHS: IN0_Create has not been called")
      call IO_Assert(io, (iElementString <= in0%iNStr), &
           "IN0_ComputeRHS: Illegal string number")
    end if
    ! Compute the correction needed in the jump...
    str => in0%Strings(iElementString)
    if (io%lDebug) then
      call IO_Assert(io, (iElementVertex <= str%iNPts), &
           "IN0_ComputeRHS: Illegal vertex number")
      if (iElementVertex == str%iNPts .and. .not. str%lClosed) then
        call IO_Assert(io, (iElementFlag == kIN0_Vertex2), &
             "IN0_ComputeRHS: No center strength at last vertex of open string")
      end if
    end if
    vtx => str%Vertices(iElementVertex)

    ! The sign on the "Jump" is negative due to the i*i factor in the computations
    if (lDirect) then
      rRHS = - (vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)) * vtx%rCheckPot(iElementFlag) + &
             vtx%rLeftT(iElementFlag) * vtx%rRightT(iElementFlag) * &
             (str%rLeftB-str%rRightB + rHALF*(vtx%rLeftH(iElementFlag)-vtx%rRightH(iElementFlag)))
    else
      rRHS = - (vtx%rLeftT(iElementFlag) - vtx%rRightT(iElementFlag)) * vtx%rCheckPot(iElementFlag) - &
             vtx%rRightT(iElementFlag) * rFDP_PotentialJump(io, fdp, vtx%cCPZ(iElementFlag), vtx%iFDPIndex) + &
             vtx%rLeftT(iElementFlag) * vtx%rRightT(iElementFlag) * &
             (str%rLeftB-str%rRightB + rHALF*(vtx%rLeftH(iElementFlag)-vtx%rRightH(iElementFlag)))
    end if

    return
  end function rIN0_ComputeRHS


  subroutine IN0_StoreResult(io, in0, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
    !! subroutine IN0_StoreResult
    !!
    !! Stores the results of a solution for a single equation associated with
    !! the IN0 module.
    !!
    !! Calling Sequence:
    !!    IN0_StoreResult(in0, cCPZ, iEqType, cOrientation, rRHS)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    real :: rValue
    !!             The new result value from the solution vector
    !!   (in)    integer :: iElementType
    !!             Element type(always ELEM_IN0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!             For IN0, the constants kIN0_Vertex and kHBO_Center are used.
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    real(kind=AE_REAL), intent(in) :: rValue
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_StoreResult: IN0_Create has not been called")
      call IO_Assert(io, (iElementString >= 1 .and. iElementString <= in0%iNStr), &
           "IN0_StoreResult: Bad element string ID")
    end if

    str => in0%Strings(iElementString)

    if (io%lDebug) then
      call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= str%iNPts), &
           "IN0_StoreResult: Bad element vertex ID")
    end if

    ! All is well.  Store the result...
    vtx => str%Vertices(iElementVertex)
    if (lDirect) then
      vtx%rStrength = rZERO
    end if

    select case (iElementFlag)
      case (kIN0_Vertex)
        vtx%rStrength(1) = vtx%rStrength(1) + rValue
      case (kIN0_Center)
        vtx%rStrength(2) = vtx%rStrength(2) + rValue
      case (kIN0_Vertex2)
        vtx%rStrength(3) = vtx%rStrength(3) + rValue
    end select

    return
  end subroutine IN0_StoreResult


  subroutine IN0_Update(io, in0, fdp)
    !! subroutine IN0_StoreResult
    !!
    !! Updates the underlying function objects for the specified layer.
    !!
    !! Calling Sequence:
    !!    IN0_Update(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             FDP_COLLECTION object to be used
    !!   (in)    type(IO_status), pointer :: io
    !!             pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, irv
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: this_vtx, next_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_Update: IN0_Create has not been called")
      call IO_Assert(io, (associated(fdp)), &
           "IN0_Update: Illegal FDP_COLLECTION object")
    end if

    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      if (str%lClosed) then
        do iVtx = 1, str%iNPts
          this_vtx => str%Vertices(iVtx)
          cRho1 = cmplx(rZERO, this_vtx%rStrength(1), AE_REAL)
          cRho2 = cmplx(rZERO, this_vtx%rStrength(2), AE_REAL)
          if (iVtx < str%iNPts) then
            next_vtx => str%Vertices(iVtx+1)
          else
            next_vtx => str%Vertices(1)
          end if
          cRho3 = cmplx(rZERO, next_vtx%rStrength(1), AE_REAL)
          call FDP_Update(io, fdp, this_vtx%iFDPIndex, (/cRho1, cRho2, cRho3/))
        end do
      else
        do iVtx = 1, str%iNPts-1
          this_vtx => str%Vertices(iVtx)
          cRho1 = cmplx(rZERO, this_vtx%rStrength(1), AE_REAL)
          cRho2 = cmplx(rZERO, this_vtx%rStrength(2), AE_REAL)
          if (iVtx < str%iNPts-1) then
            next_vtx => str%Vertices(iVtx+1)
            cRho3 = cmplx(rZERO, next_vtx%rStrength(1), AE_REAL)
          else
            cRho3 = cmplx(rZERO, this_vtx%rStrength(3), AE_REAL)
          end if
          call FDP_Update(io, fdp, this_vtx%iFDPIndex, (/cRho1, cRho2, cRho3/))
        end do
      end if
    end do

    ! If we've done this once, it's not the initial iteration anymore
    in0%lInitialIteration = .false.

    return
  end subroutine IN0_Update


  subroutine IN0_ResetIterator(io, in0)
    !! subroutine IN0_ResetIterator
    !!
    !! Resets the module's iterator prior to traversing for check data
    !!
    !! Calling Sequence:
    !!    call IN0_ResetIterator(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: in0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(IO_STATUS), pointer :: io

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_ResetIterator: IN0_Create has not been called")
    end if

    in0%iIterStr = 1
    in0%iIterVtx = 0
    in0%iIterFlag = kIN0_Center

    return
  end subroutine IN0_ResetIterator


  function IN0_NextIterator(io, in0) result(itr)
    !! function IN0_NextIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call IN0_NextIterator(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION to be used
    !!   (in)    type(IO_STATUS), pointer :: in0
    !!             Tracks error conditions
    !!
    !! Return Value:
    !!   type(ITERATOR_RESULT), pointer :: itr
    !!     Pointer to the information for data retrieval
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(ITERATOR_RESULT), pointer :: itr
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_NextIterator: IN0_Create has not been called")
    end if

    if (in0%iIterStr > in0%iNStr) then
      nullify(itr)
      return
    end if

    ! Update the iterator(this is a mess!)
    if (in0%iIterFlag == kIN0_Vertex2) then
      ! Our last iterator step hit the end of an open string...
      in0%iIterStr = in0%iIterStr+1
      in0%iIterVtx = 1
      in0%iIterFlag = kIN0_Vertex
      if (in0%iIterStr > in0%iNStr) then
        nullify(itr)
        return
      end if
    else if (in0%iIterFlag == kIN0_Center) then
      ! We're at a center interval -- what next?
      if (in0%Strings(in0%iIterStr)%lClosed .and. &
          in0%iIterVtx+1 > in0%Strings(in0%iIterStr)%iNPts) then
        ! Aha -- end of the buffer. Handle that.
        in0%iIterFlag = kIN0_Vertex
        in0%iIterStr = in0%iIterStr+1
        in0%iIterVtx = 1
        if (in0%iIterStr > in0%iNStr) then
          nullify(itr)
          return
        end if
      else if (.not. in0%Strings(in0%iIterStr)%lClosed .and. &
             in0%iIterVtx+1 > in0%Strings(in0%iIterStr)%iNPts-1) then
        ! Open string -- Move on to the last vertex
        in0%iIterFlag = kIN0_Vertex2
      else
        ! Not at the end of the string -- move to the next vertex
        in0%iIterVtx = in0%iIterVtx + 1
        in0%iIterFlag = kIN0_Vertex
      end if
    else
      ! We were at a "first" vertex end -- move to the center of the same segment
      in0%iIterFlag = kIN0_Center
    end if

    allocate(itr)
    itr%iElementType = ELEM_IN0
    itr%iElementString = in0%iIterStr
    itr%iElementVertex = in0%iIterVtx
    itr%iElementFlag = in0%iIterFlag
    itr%iValueSelector = VALUE_POTENTIAL
    allocate(itr%cZ(1))
    itr%cZ(1) = in0%Strings(in0%iIterStr)%Vertices(in0%iIterVtx)%cCPZ(in0%iIterFlag)
    dom => IN0_FindDomain(io, in0, itr%cZ(1))

    return
  end function IN0_NextIterator


  subroutine IN0_SetIterator(io, in0, fdp, itr, cValue, lLinearize)
    !! function IN0_SetIterator
    !!
    !! Advances the module's iterator one step
    !!
    !! Calling Sequence:
    !!    call IN0_SetIterator(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION to be used
    !!   (in)    type(ITERATOR_RESULT), pointer :: itr
    !!             Pointer to the information for data retrieval
    !!   (in)    complex :: cValue
    !!             The value retrieved from the color
    !!   (in)    type(IO_STATUS), pointer :: in0
    !!             Tracks error conditions
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(FDP_COLLECTION), pointer :: fdp
    type(ITERATOR_RESULT), pointer :: itr
    complex(kind=AE_REAL), intent(in) :: cValue
    logical, intent(in) :: lLinearize
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: left, right
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx
    real(kind=AE_REAL) :: H1, H2, rHead

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_NextIterator: IN0_Create has not been called")
      call IO_Assert(io, (in0%iIterStr <= in0%iNStr), &
           "IN0_SetIterator: Iterator out of range")
    end if

    ! Clear the "regen" flag unless something interesting happens...
    in0%iRegenerate = 0

    str => in0%Strings(itr%iElementString)
    vtx => str%Vertices(itr%iElementVertex)
    if (itr%iValueSelector == VALUE_POTENTIAL) then
      !print *,"SET",itr%iElementString, itr%iElementVertex, itr%iElementFlag, real(cValue)
      vtx%rCheckPot(itr%iElementFlag) = real(cValue, AE_REAL)
      ! Do I update coefficients?
      if (lLinearize) then
        left => IN0_FindDomainID(io, in0, str%iLeftID)
        right => IN0_FindDomainID(io, in0, str%iRightID)
        rHead = rHALF * (rIN0_PotentialToHead(io, in0, vtx%rCheckPot(itr%iElementFlag), vtx%cCPZ(itr%iElementFlag)) + &
                rIN0_PotentialToHead(io, in0, vtx%rCheckPot(itr%iElementFlag), vtx%cCPZ(itr%iElementFlag)))
        H1 = rIN0_DomainThickness(io, in0, left, rHead)
        H2 = rIN0_DomainThickness(io, in0, right, rHead)
        if (H1 /= vtx%rLeftH(itr%iElementFlag) .or. H2 /= vtx%rRightH(itr%iElementFlag)) then
          in0%iRegenerate = 1
          vtx%rLeftH(itr%iElementFlag) = H1
          vtx%rRightH(itr%iElementFlag) = H2
          vtx%rLeftT(itr%iElementFlag) = H1 * left%rHydCond
          vtx%rRightT(itr%iElementFlag) = H2 * right%rHydCond
        end if
        !print *,'in0 update',rHead,str%iLeftID,str%iRightID,vtx%rLeftH(itr%iElementFlag),vtx%rRightH(itr%iElementFlag),vtx%rLeftT(itr%iElementFlag),vtx%rRightT(itr%iElementFlag)
      end if
      vtx%rError(itr%iElementFlag) = - (vtx%rLeftT(itr%iElementFlag) - &
          vtx%rRightT(itr%iElementFlag)) * vtx%rCheckPot(itr%iElementFlag) - &
          vtx%rRightT(itr%iElementFlag) * rFDP_PotentialJump(io, fdp, vtx%cCPZ(itr%iElementFlag), vtx%iFDPIndex) + &
          vtx%rLeftT(itr%iElementFlag) * vtx%rRightT(itr%iElementFlag) * &
          (str%rLeftB-str%rRightB + rHALF*(vtx%rLeftH(itr%iElementFlag)-vtx%rRightH(itr%iElementFlag)))
    end if

    return
  end subroutine IN0_SetIterator


  subroutine IN0_Inquiry(io, in0, iLU)
    !! subroutine IN0_Inquiry
    !!
    !! Writes an inquiry report for all barriers to iLU
    !!
    !! Calling Sequence:
    !!    call IN0_Inquiry(in0, iLU)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    integer :: iLU
    !!             The output LU to receive output
    !!   (in)    type(IO_status), pointer :: io
    !!             pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iLU
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStr, iVtx, i
    real(kind=AE_REAL) :: rLength, rARecip
    real(kind=AE_REAL), dimensiON(2) :: rError
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx, next

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_Inquiry: IN0_Create has not been called")
    end if

    write (unit=iLU, &
           fmt="(""#IN0, ID, VTX, X1, Y1, X2, Y2, LENGTH, STRENGTH1, STRENGTH2, POT1, POT1, ERROR1, ERROR2"")")
    do iStr = 1, in0%iNStr
      str => in0%Strings(iStr)
      do iVtx = 1, str%iNPts-1
        vtx => str%Vertices(iVtx)
        next => str%Vertices(iVtx+1)
        if (iVtx < str%iNPts) then
          rLength = abs(str%Vertices(iVtx+1)%cZ - vtx%cZ)
        else
          rLength = rZERO
        end if
        rARecip = rIN0_ARecip(io, in0, iStr)
        do i = 1, 2
          if (vtx%rStrength(i) /= rZERO) then
            rError(i) = vtx%rError(i) / vtx%rStrength(i)
          else
            rError(i) = rZERO
          end if
        end do
        write (unit=iLU, &
               fmt="('IN0', 2(', ', i9), 11(', ', e16.8))" &
               ) str%iID, iVtx, cIO_WorldCoords(io, vtx%cZ), cIO_WorldCoords(io, next%cZ), &
                 rLength, vtx%rStrength(1), vtx%rStrength(2), vtx%rCheckPot(1), vtx%rCheckPot(2), &
                 rError(1), rError(2)
      end do
    end do

    return
  end subroutine IN0_Inquiry


  subroutine IN0_Read(io, in0)
    !! subroutine IN0_Read
    !!
    !! Reads the aquifer information for layer iL from the input LU
    !!
    !! Calling Sequence:
    !!    call IN0_Read(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: aqu
    !!             Layer number to be read
    !!
    !! The format of the aqu section of the input file appears as follows:
    !!
    !! IN0
    !! DOM nvertices base thickness hyd-cond porosity
    !!     (x, y)
    !!     (x, y)
    !!     ... up to nvertices
    !! DOM ...
    !! ... up to ninho
    !! END
    !! END
    !!
    !! NOTE: It is assumed that the IN0 line was found by the caller
    !!
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(IO_STATUS), pointer :: io
    ! Locals -- for Directive parsing
    type(DIRECTIVE), dimension(3), parameter :: dirDirectives = &
                       (/dirEND, dirDOM, dirSTR/)
    ! Locals -- Input values
    integer(kind=AE_INT) :: iParseMode
    integer(kind=AE_INT) :: iVtx
    integer(kind=AE_INT) :: iOpCode
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iMax
    integer(kind=AE_INT) :: iID
    integer(kind=AE_INT) :: iLeftID
    integer(kind=AE_INT) :: iRightID
    logical :: lClosed
    real(kind=AE_REAL) :: rBase
    real(kind=AE_REAL) :: rHydCond
    real(kind=AE_REAL) :: rPorosity
    real(kind=AE_REAL) :: rThickness
    real(kind=AE_REAL) :: rRefHead
    real(kind=AE_REAL) :: rAvgHead
    complex(kind=AE_REAL) :: cZ
    logical :: lFlag
    type(IN0_DOMAIN), pointer :: dom
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx
    character(len=32) :: sTag
    ! [ CONSTANTS ]
    integer(kind=AE_INT), parameter :: PARSE_NONE = 0
    integer(kind=AE_INT), parameter :: PARSE_DOMAIN = 1
    integer(kind=AE_INT), parameter :: PARSE_STRING = 2

    iParseMode = PARSE_NONE
    call IO_MessageText(io, "  Reading IN0 module input")

    call IO_Assert(io, (associated(in0)), "IN0_Read: IN0_Create has not been called")

    do
      call IO_InputRecord(io, dirDirectives, iOpCode)
      select case (iOpCode)
        case (kOpError)
          ! A RunTime error was found during a file read operation. This
          ! condition is fatal; warn the user, and exit.
          call IO_Assert(io, .false., "IN0_Read: I/O Error")
        case (kOpFileEOF)
          ! EOF is unexpected for all Mod "ifXXXRead" routines.
          call IO_Assert(io, .false., "IN0_Read: Unexpected EOF")
        case (kOpData)
          ! A data line was found. If we have a specified perimeter, add the point
          ! to the perimeter.
          select case (iParseMode)
            case (PARSE_NONE)
              call IO_Assert(io, .false., "IN0_Read: Unexpected data record")
            case (PARSE_DOMAIN)
              call IO_Assert(io, (associated(dom%cZ)), "IN0_Read: No DOM directive")
              call IO_Assert(io, (dom%iNPts < size(dom%cZ)), &
                   "IN0_Read: Space exhausted")
              write (unit=sTag, fmt=*) 'cZ dom', in0%iNDom, dom%iNPts+1
              cZ = cIO_GetCoordinate(io, sTag, extents=.true., check_points=dom%cZ)
              dom%iNPts = dom%iNPts+1
              dom%cZ(dom%iNPts) = cZ
            case (PARSE_STRING)
              call IO_Assert(io, (associated(str%Vertices)), "IN0_Read: No STR directive")
              call IO_Assert(io, (str%iNPts < size(str%Vertices)), &
                   "IN0_Read: Space exhausted")
              write (unit=sTag, fmt=*) 'cZ str', in0%iNStr, str%iNPts+1
              cZ = cIO_GetCoordinate(io, sTag, extents=.true., check_points=str%Vertices(:)%cZ)
              str%iNPts = str%iNPts+1
              str%Vertices(str%iNPts)%cZ = cZ
          end select
        case (kOpEND)
          ! EOD mark was found. Exit the file parser.
          return
        case (kOpDOM)
          ! Start a new domain
          call IO_Assert(io, (associated(in0%Domains)), &
               "IN0_Read: No domains have been allocated")
          call IO_Assert(io, (in0%iNDom < size(in0%Domains)), &
               "IN0_Read: Space exhausted")
          iMax = iIO_GetInteger(io, 'iMax', minimum = 3)
          rBase = rIO_GetReal(io, 'rBase')
          rThickness = rIO_GetReal(io, 'rThickness', minimum = rTINY)
          rHydCond = rIO_GetReal(io, 'rHydCond', minimum = rTINY)
          rPorosity = rIO_GetReal(io, 'rPorosity', minimum = rTINY)
          rAvgHead = rIO_GetReal(io, 'rAvgHead', minimum = rBase)
          iID = iIO_GetInteger(io, 'iID', forbidden=in0%Domains%iID)

          in0%iNDom = in0%iNDom+1
          dom => in0%Domains(in0%iNDom)
          dom%iNPts = 0
          dom%rBase = rBase
          dom%rThickness = rThickness
          dom%rHydCond = rHydCond
          dom%rPorosity = rPorosity
          dom%rAvgHead = rAvgHead
          dom%iID = iID
          dom%iInsideDomain = iID
          dom%iOutsideDomain = 0
          allocate(dom%cZ(iMax), stat = iStat)
          call IO_Assert(io, (iStat == 0), "IN0_Read: Allocation failed")
          dom%cZ = cZERO
          iParseMode = PARSE_DOMAIN
        case (kOpSTR)
          ! Start a new string
          call IO_Assert(io, (associated(in0%Strings)), &
               "IN0_Read: No strings have been allocated")
          call IO_Assert(io, (in0%iNStr < size(in0%Strings)), &
               "IN0_Read: Space exhausted")
          iMax = iIO_GetInteger(io, 'iMax')
          iLeftID = iIO_GetInteger(io, 'iLeftID')
          iRightID = iIO_GetInteger(io, 'iRightID')
          lClosed = lIO_GetLogical(io, 'lClosed')
          iID = iIO_GetInteger(io, 'iID', forbidden=in0%Strings%iID)

          in0%iNStr = in0%iNStr+1
          str => in0%Strings(in0%iNStr)
          str%iNPts = 0
          str%iLeftID = iLeftID
          str%iRightID = iRightID
          str%lClosed = lClosed
          str%iID = iID
          allocate(str%Vertices(iMax), stat = iStat)
          call IO_Assert(io, (iStat == 0), "IN0_Read: Allocation failed")
          do iVtx = 1, size(str%Vertices)
            vtx => str%Vertices(iVtx)
            vtx%cZ = cZERO
            vtx%rStrength = rZERO
            vtx%rCheckPot(:) = rZERO
            vtx%iFDPIndex = -1
          end do
          iParseMode = PARSE_STRING
        case default
      end select
    end do

    call IO_MessageText(io, "Leaving IN0 module")

    return
  end subroutine IN0_Read


  subroutine IN0_Report(io, in0)
    !! subroutine IN0_Report
    !!
    !! Writes a debugging report for all line-sinks to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call IN0_Report(in0)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    type(IO_status), pointer :: io
    !!             pointer toIO_STATUS structure
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom, iStr, iVtx, i
    integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
    real(kind=AE_REAL), dimension(2) :: rError
    type(IN0_DOMAIN), pointer :: dom
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "IN0_Report: IN0_Create has not been called")
    end if

    call HTML_Header('Module IN0', 1)
    call HTML_Header('Inhomogeneity information', 2)

    if (.not. associated(in0%Domains)) then
      call HTML_Header('No domains allocated', 3)
    else
      call HTML_StartTable()
      call HTML_AttrInteger('Number of domains', in0%iNDom)
      call HTML_AttrInteger('Number of strings', in0%iNStr)
      call HTML_AttrInteger('Number of FWL functions', iIN0_GetInfo(io, in0, SIZE_FWL, 0))
      call HTML_AttrInteger('Number of FPD functions', iIN0_GetInfo(io, in0, SIZE_FPD, 0))
      call HTML_AttrInteger('Number of FDP functions', iIN0_GetInfo(io, in0, SIZE_FDP, 0))
      call HTML_AttrInteger('Number of equations', iIN0_GetInfo(io, in0, SIZE_EQUATIONS, 0))
      call HTML_AttrInteger('Number of unknowns', iIN0_GetInfo(io, in0, SIZE_UNKNOWNS, 0))
      call HTML_EndTable()

      do iDom = 1, in0%iNDom
        call HTML_Header('Domain Information', 3)
        call HTML_StartTable()
        dom => in0%Domains(iDom)
        call HTML_AttrInteger('Domain number ', iDom)
        call HTML_AttrInteger('ID', dom%iID)
        call HTML_AttrReal('Base elevation', dom%rBase)
        call HTML_AttrReal('Thickness', dom%rThickness)
        call HTML_AttrReal('Hydraulic conductivity', dom%rHydCond)
        call HTML_AttrReal('Porosity', dom%rPorosity)
        call HTML_EndTable()

        if (.not. associated(dom%cZ)) then
          call HTML_Header('No vertices -- infinite aquifer', 4)
        else
          call HTML_Header('Vertices', 4)
          call HTML_StartTable()
          call HTML_TableHeader((/'X', 'Y'/))
          do iVtx = 1, dom%iNPts
            call HTML_StartRow()
            call HTML_ColumnComplex((/cIO_WorldCoords(io, dom%cZ(iVtx))/))
            call HTML_EndRow()
          end do
          call HTML_EndTable()
        end if
      end do

      do iStr = 1, in0%iNStr
        call HTML_Header('String Information', 3)
        call HTML_StartTable()
        str => in0%Strings(iStr)
        call HTML_AttrInteger('String number', iStr)
        call HTML_AttrInteger('ID', str%iID)
        call HTML_AttrInteger('Left domain', str%iLeftID)
        call HTML_AttrInteger('Right domain', str%iRightID)
        call HTML_AttrLogical('Closed?', str%lClosed)
        call HTML_EndTable()

        call HTML_Header('Vertices', 4)
        call HTML_StartTable()
        call HTML_TableHeader((/'Vertex', 'X     ', 'Y     ', 'FDP # ', 'V1 Str', 'C Str ', 'V1 Err', 'C Err '/))
        do iVtx = 1, str%iNPts
          vtx => str%Vertices(iVtx)
          do i = 1, 2
            if (vtx%rStrength(i) /= rZERO) then
              rError(i) = vtx%rError(i) / vtx%rStrength(i)
            else
              rError(i) = rZERO
            end if
          end do
          call HTML_StartRow()
          call HTML_ColumnInteger((/iVtx/))
          call HTML_ColumnComplex((/cIO_WorldCoords(io, vtx%cZ)/))
          call HTML_ColumnInteger((/vtx%iFDPIndex/))
          call HTML_ColumnReal(vtx%rStrength(1:2), 'e13.6')
          call HTML_ColumnReal(rError, 'e13.6')
          call HTML_EndRow()
        end do
        call HTML_EndTable()
      end do

    end if

    return
  end subroutine IN0_Report


  function rIN0_ARecip(io, in0, iString) result(rARecip)
    !! function IN0_ARecip(in0, iString)
    !!
    !! Returns the factor A*, defined as the reciprocal of(k+ - k-) / k-
    !! where k+ is the inside conductivity and k- is the outside conductivity
    !!
    !! Calling Sequence:
    !!    rARecip = lIN0_ARecip(in0, iString)
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iString
    !!                String index in in0
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!    (out)     real :: rARecip
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iString
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rARecip
    ! [ LOCALS ]
    type(IN0_STRING), pointer :: str
    type(IN0_DOMAIN), pointer :: left, right

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "rIN0_Arecip: IN0_Create has not been called")
      call IO_Assert(io, (iString >= 1 .and. iString <= in0%iNDom), &
           "rIN0_ARecip: Bad domain index")
    end if

    str => in0%Strings(iString)
    left => IN0_FindDomainID(io, in0, str%iLeftID)
    right => IN0_FindDomainID(io, in0, str%iRightID)

    rARecip = right%rHydCond / (left%rHydCond - right%rHydCond)

    return
  end function rIN0_ARecip


  function lIN0_PointInsideDomain(io, in0, iDomain, cZ) result(lInside)
    !! function lIN0_PointInsideDomain
    !!
    !! Tests to see if the specified point is inside the Domain
    !!
    !! Calling Sequence:
    !!    lresult = lIN0_PointInsideDomain(in0, iDomain, cZ)
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iDomain
    !!                Domain index in in0
    !!    (in)      complex :: cZ
    !!                The point to be checked
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           logical :: lInside
    !!             .true. if the point is inside the Domain,
    !!             .false. if the point is outside the Domain
    !!
    !! NOTE:  Uses an algorithm from the comp.graphics FAQ -- THANKS!
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iDomain
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lInside
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, i, ii
    real(kind=AE_REAL) :: rSum
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_PointInsideDomain: IN0_Create has not been called")
      call IO_Assert(io, (iDomain >= 1 .and. iDomain <= in0%iNDom), &
           "lIN0_PointInsideDomain: Bad domain index")
    end if

    ! Use the "arctangent method" for now; it's slow but reliable
    dom => in0%Domains(iDomain)
    if (size(dom%cZ) == 0) then
      lInside = .true.
    else
      if (any((cZ == dom%cZ(1:dom%iNPts)), 1)) then
        lInside = .true.
      else
        rSum = rZERO
        ii = dom%iNPts
        do i = 1, dom%iNPts
          rSum = rSum + aimag(log((cZ-dom%cZ(i)) / (cZ-dom%cZ(ii))))
          ii = i
        end do
        lInside = (rSum > rONE)
      end if
    end if

    return
  end function lIN0_PointInsideDomain


  function lIN0_LineIntersectsDomain(io, in0, iDomain, cZ1, cZ2) result(lIntersects)
    !! function lIN0_LineIntersectsDomain
    !!
    !! Tests to see if the line segment cZ1-cZ2 intersects the Domain
    !!
    !! Calling Sequence:
    !!    lresult = lIN0_LineIntersectsDomain(in0, iDomain, cZ1, cZ2)
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iDomain
    !!                Domain index in in0
    !!    (in)      complex :: cZ1, cZ2
    !!                End-points of the line to be tested
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           logical :: lIntersects
    !!             .true. if the line segment intersects with any Domain edge
    !!             .false. otherwise
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iDomain
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lIntersects
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    complex(kind=AE_REAL) :: cMapZ1, cMapZ2
    real(kind=AE_REAL) :: rXInt
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_LineIntersectsDomain: IN0_Create has not been called")
      call IO_Assert(io, (iDomain >= 1 .and. iDomain <= in0%iNDom), &
           "lIN0_LineIntersectsDomain: Bad domain index")
    end if

    dom => in0%Domains(iDomain)
    lIntersects = .false.
    j = dom%iNPts
    do i = 1, dom%iNPts
      cMapZ1 = (dom%cZ(i) - rHALF * (cZ2+cZ1)) / (rHALF * (cZ2-cZ1))
      cMapZ2 = (dom%cZ(j) - rHALF * (cZ2+cZ1)) / (rHALF * (cZ2-cZ1))
      if ((aimag(cMapZ1) >= rZERO .and. aimag(cMapZ2) < rZERO) .or. &
          (aimag(cMapZ2) >= rZERO .and. aimag(cMapZ1) < rZERO)) then
        rXInt = real(cMapZ1) - aimag(cMapZ1) * ((real(cMapZ2)-real(cMapZ1)) / &
                (aimag(cMapZ2)-aimag(cMapZ1)))
        if (abs(rXInt) <= rONE) then
          lIntersects = .true.
          return
        end if
      end if
      j = i
    end do

    return
  end function lIN0_LineIntersectsDomain


  function lIN0_DomainInsideDomain(io, in0, iDomain, iDomain2) result(lInside)
    !! function lIN0_DomainInsideDomain
    !!
    !! Tests to see if the Domain 'iDomain2' is inside the Domain 'iDomain'
    !!
    !! Calling Sequence:
    !!    lresult = lIN0_DomainInsideDomain(in0, iDomain, iDomain2)
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iDomain
    !!                Domain index in in0
    !!    (in)      integer :: iDomain2
    !!                Domain index in in0
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           logical :: lInside
    !!             .true. if the Domain is completely inside the Domain,
    !!             .false. otherwise
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iDomain
    integer(kind=AE_INT), intent(in) :: iDomain2
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lInside
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    type(IN0_DOMAIN), pointer :: dom, dom2
    type(IN0_VERTEX), pointer :: vtx_i, vtx_j

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_DomainInsideDomain: IN0_Create has not been called")
      call IO_Assert(io, ((iDomain >= 1 .and. iDomain <= in0%iNDom) .and. &
           (iDomain2 >= 1 .and. iDomain2 <= in0%iNDom)), &
           "lIN0_DomainInsideDomain: Bad domain index")
    end if

    dom => in0%Domains(iDomain)
    dom2 => in0%Domains(iDomain2)
    if (size(dom%cZ, 1) == 0) then
      lInside = .true.
    else
      if (size(dom2%cZ) == 0) then
        lInside = .false.
      else
        lInside = .true.
        j = dom2%iNPts
        do i = 1, dom2%iNPts
          if (.not. lIN0_PointInsideDomain(io, in0, iDomain, dom2%cZ(i)) .or. &
              lIN0_LineIntersectsDomain(io, in0, iDomain, dom2%cZ(i), dom2%cZ(j))) then
            lInside = .false.
            return
          end if
          j = i
        end do
      end if
    end if

    return
  end function lIN0_DomainInsideDomain


  function lIN0_DomainOverlapsDomain(io, in0, iDomain, iDomain2) result(lOverlap)
    !! function lIN0_DomainOverlapsDomain
    !!
    !! Tests to see if the Domain 'iDomain2' overlaps the Domain 'iDomain'
    !!
    !! Calling Sequence:
    !!    lresult = lIN0_DomainOverlapsDomain(in0, iDomain, iDomain2)
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iDomain
    !!                Domain index in in0
    !!    (in)      integer :: iDomain2
    !!                Domain index in in0
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           logical :: lInside
    !!             .true. if Domain2 overlaps Domain,
    !!             .false. otherwise
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iDomain
    integer(kind=AE_INT), intent(in) :: iDomain2
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lOverlap
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    type(IN0_DOMAIN), pointer :: dom, dom2
    type(IN0_VERTEX), pointer :: vtx_i, vtx_j

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_DomainOverlapsDomain: IN0_Create has not been called")
      call IO_Assert(io, ((iDomain >= 1 .and. iDomain <= in0%iNDom) .and. &
           (iDomain2 >= 1 .and. iDomain2 <= in0%iNDom)), &
           "lIN0_DomainOverlapsDomain: Bad domain index")
    end if

    dom => in0%Domains(iDomain)
    dom2 => in0%Domains(iDomain2)
    if (size(dom%cZ, 1) == 0) then
      lOverlap = .true.
    else
      if (size(dom2%cZ) == 0) then
        lOverlap = .false.
      else
        lOverlap = .false.
        j = size(dom2%cZ, 1)
        do i = 1, size(dom2%cZ, 1)
          if (lIN0_LineIntersectsDomain(io, in0, iDomain, dom2%cZ(i), dom2%cZ(j))) then
            lOverlap = .true.
            return
          end if
          j = i
        end do
      end if
    end if

    return
  end function lIN0_DomainOverlapsDomain


  function IN0_FindDomain(io, in0, cZ) result(dom)
    !! function IN0_FindDomain
    !!
    !! Returns a pointer to the domain containing the point 'cZ'
    !!
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      complex :: cZ
    !!                The point in question
    !!
    !! Return Value:
    !!           type(IN0_DOMAIN), pointer :: dom
    !!             Pointer to the domain containing cZ
    !!
    !! Note:
    !!   This version does not support nested inhomogeneities
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(IN0_DOMAIN), pointer :: dom
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! Default to the infinite domain
    dom => in0%Domains(1)
    do iDom = 2, in0%iNDom
      if (lIN0_PointInsideDomain(io, in0, iDom, cZ)) then
        dom => in0%Domains(iDom)
      end if
    end do

    return
  end function IN0_FindDomain


  function iIN0_DomainIDIndex(io, in0, iID) result(iDom)
    !! function IN0_FindDomain
    !!
    !! Returns a pointer to the domain with ID 'iID'
    !!
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iID
    !!                The ID in question
    !!
    !! Return Value:
    !!           type(IN0_DOMAIN), pointer :: dom
    !!             Pointer to the domain with ID iID
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iID
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iRes
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! Default to the infinite domain
    iRes = 1
    do iDom = 1, in0%iNDom
      if (in0%Domains(iDom)%iID == iID) then
        iRes = iDom
        return
      end if
    end do
    call IO_Assert(io, .false., 'iIN0_DomainIDIndex: Could not find ID')

    return
  end function iIN0_DomainIDIndex


  function IN0_FindDomainID(io, in0, iID) result(dom)
    !! function IN0_FindDomain
    !!
    !! Returns a pointer to the domain with ID 'iID'
    !!
    !!
    !! Arguments:
    !!    (in)      type(IN0_COLLECTION), pointer :: in0
    !!                The IN0_COLLECTION object of interest
    !!    (in)      integer :: iID
    !!                The ID in question
    !!
    !! Return Value:
    !!           type(IN0_DOMAIN), pointer :: dom
    !!             Pointer to the domain with ID iID
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: iID
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(IN0_DOMAIN), pointer :: dom
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! Default to the infinite domain
    dom => in0%Domains(1)
    do iDom = 2, in0%iNDom
      if (in0%Domains(iDom)%iID == iID) then
        dom => in0%Domains(iDom)
      end if
    end do

    return
  end function IN0_FindDomainID

  !! UTILITY ROUTINES
  !! These routines provide computational aids for AEM models


  function rIN0_HeadToPotential(io, in0, rHead, cZ) result(rPot)
    !! real function rIN0_HeadToPotential
    !!
    !! Converts head to a discharge potential based on the in0ifer properties
    !! at cZ. Returns the(real) potential.
    !!
    !! Calling Sequence:
    !!    rPot = rIN0_HeadToPotential(iL, rHead, cZ)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    real :: rHead
    !!             The head value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    type(IO_status), pointer :: io
    !!             pointer toIO_STATUS structure

    !!
    !! Return Value:
    !!           real :: rPot
    !!             The discharge potential corresponding to the head 'rHead'
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    real(kind=AE_REAL), intent(in) :: rHead
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPot
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHd
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! For the general case, rK, rBase, rPerm and rThick are functions of position
    dom => IN0_FindDomain(io, in0, cZ)

    rHd = rHead - dom%rBase
    if (rHd > dom%rThickness) then
      rPot = dom%rHydCond * dom%rThickness * rHd - rHALF * dom%rHydCond * dom%rThickness**2
    else
      rPot = rHALF * dom%rHydCond * rHd**2
    end if

    return
  end function rIN0_HeadToPotential


  function rIN0_PotentialToHead(io, in0, rPot, cZ, lTest) result(rHead)
    !! real function rIN0_HeadToPotential
    !!
    !! Converts head to a discharge potential based on the in0ifer properties
    !! at cZ. Returns the(real) potential.
    !!
    !! Calling Sequence:
    !!    rPot = rIN0_HeadToPotential(in0, rHead, cZ)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    real :: rHead
    !!             The head value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!    (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           real :: rHead
    !!             The head corresponding to the discharge potential 'rPot'
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    logical, intent(in), optional :: lTest
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHead
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! For the general case, rK, rBase, rPerm and rThick are functions of position
    dom => IN0_FindDomain(io, in0, cZ)

    if (present(lTest)) print *,"P->H", rPot, cZ, dom%iID, dom%rBase, dom%rHydCond, dom%rThickness, dom%rTopPot

    if (rPot < rZERO) then
      rHead = dom%rBase
      if (present(lTest)) print *,"  P<B", rHead
    else if (rPot > dom%rTopPot) then
      rHead = (rPot + rHALF*dom%rHydCond*dom%rThickness**2) / (dom%rHydCond * dom%rThickness) + dom%rBase
      if (present(lTest)) print *,"  CNF", rHead
    else
      rHead = sqrt(rTWO * rPot / dom%rHydCond) + dom%rBase
      if (present(lTest)) print *,"  UNC", rHead
    end if

    return
  end function rIN0_PotentialToHead


  function cIN0_DischargeToVelocity(io, in0, cDischarge, cZ, rPot) result(cVelocity)
    !! function cIN0_DischargeToVelocity
    !!
    !! Converts Discharge to velocity, using the potential and the location cZ
    !! to compute saturated thickness
    !!
    !! Calling Sequence:
    !!    cV = cIN0_DischargeToVelocity(in0, rDischarge, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    real :: rDischarge
    !!             The discharge value to be converted
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer toIO_STATUS structure
    !!
    !! Return Value:
    !!           complex :: cV
    !!             The velocity
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cDischarge
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cVelocity
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rSatdThick
    type(IN0_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(in0)), &
           "lIN0_FindDomain: IN0_Create has not been called")
    end if

    ! For the general case, rK, rBase, rPerm and rThick are functions of position
    dom => IN0_FindDomain(io, in0, cZ)

    ! For the general case, rK, rBase, rPerm and rThick are functions of position, for now
    ! the position parameters are ignored
    if (rPot < rZERO) then
      cVelocity = cZERO
      return
    else
      ! Find the saturated thickness
      if (rPot > dom%rTopPot) then
        rSatdThick = dom%rThickness
      else
        rSatdThick = rIN0_PotentialToHead(io, in0, rPot, cZ) - dom%rBase
      end if

      ! Compute the velocity
      !print *, abs(cDischarge), rSatdThick, dom%rPorosity
      cVelocity = cDischarge / (rSatdThick * dom%rPorosity)
    end if

    return
  end function cIN0_DischargeToVelocity


  function lIN0_IsConfined(io, in0, cZ, rPot) result(lConfined)
    !! function lIN0_IsConfined
    !!
    !! Computes the saturated thickness at point cZ in layer iL where the potential is rPot
    !!
    !! Calling Sequence:
    !!    if (lIN0_IsConfined(in0, cZ, rPot)) then...
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           logical :: lConfined
    !!             .true. if flow is confined
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    logical :: lConfined
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(IN0_DOMAIN), pointer :: dom

    dom => IN0_FindDomain(io, in0, cZ)

    if (in0%lPrecondition) then
      rHead = dom%rAvgHead
    else
      rHead = rIN0_PotentialToHead(io, in0, rPot, cZ)
    end if

    if ((rHead - dom%rBase) < dom%rThickness) then
      lConfined = .false.
    else
      lConfined = .true.
    end if

    return
  end function lIN0_IsConfined


  function rIN0_SatdThickness(io, in0, cZ, rPot) result(rH)
    !! function cIN0_SatdThickness
    !!
    !! Computes the saturated thickness at point cZ in layer iL where the potential is rPot
    !!
    !! Calling Sequence:
    !!    cV = cIN0_SatdThickness(in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
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
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(IN0_DOMAIN), pointer :: dom

    dom => IN0_FindDomain(io, in0, cZ)

    if (in0%lPrecondition) then
      if (dom%rAvgHead > dom%rBase+dom%rThickness) then
        rHead = dom%rAvgHead + dom%rBase
      else
        rHead = dom%rAvgHead
      end if
    else
      rHead = rIN0_PotentialToHead(io, in0, rPot, cZ)
    end if

    if ((rHead - dom%rBase) < dom%rThickness) then
      rH = rHead - dom%rBase
    else
      rH = dom%rThickness
    end if

    return
  end function rIN0_SatdThickness


  function rIN0_Transmissivity(io, in0, cZ, rPot) result(rT)
    !! function rIN0_Transmissivity
    !!
    !! Computes the transmissivity at the point cZ where the potential is rPot
    !!
    !! Calling Sequence:
    !!    cV = cIN0_Transmissivity(in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
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
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rT
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(IN0_DOMAIN), pointer :: dom

    dom => IN0_FindDomain(io, in0, cZ)

    if (in0%lPrecondition) then
      rHead = dom%rAvgHead
    else
      rHead = rIN0_PotentialToHead(io, in0, rPot, cZ)
    end if

    if ((rHead - dom%rBase) < rONE_TENTH*dom%rThickness) then
      rT = rONE_TENTH * dom%rThickness
    else if ((rHead - dom%rBase) < dom%rThickness) then
      rT = dom%rHydCond * (rHead - dom%rBase)
    else
      rT = dom%rHydCond * dom%rThickness
    end if

    return
  end function rIN0_Transmissivity


  function rIN0_HydCond(io, in0, cZ) result(rHydCond)
    !! function rIN0_HydCond
    !!
    !! Computes the hydraulic conductivity at the point cZ
    !!
    !! Calling Sequence:
    !!    cV = cIN0_HydCond(in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           real :: rHydCond
    !!             The hydraulic conductivity
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHydCond
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: dom

    dom => IN0_FindDomain(io, in0, cZ)
    rHydCond = dom%rHydCond

    return
  end function rIN0_HydCond


  function rIN0_Base(io, in0, cZ) result(rBase)
    !! function rIN0_HydCond
    !!
    !! Computes the base elevation at the point cZ
    !!
    !! Calling Sequence:
    !!    cV = cIN0_Base(in0, cZ, rPot)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           real :: rBase
    !!             The base elevation
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBase
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: dom

    dom => IN0_FindDomain(io, in0, cZ)
    rBase = dom%rBase

    return
  end function rIN0_Base


  function rIN0_DefaultPotential(io, in0, cZ) result(rP)
    !! function rIN0_HydCond
    !!
    !! Computes the top potential at the point cZ
    !!
    !! Calling Sequence:
    !!    cV = cIN0_DefaultPotential(in0, cZ)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point in the in0ifer where the conversion is to take place
    !!
    !! Return Value:
    !!           real :: rP
    !!             The top potential
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rP
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: dom

    !$$ assert : allocated(Layers) : "No layers are allocated"
    !$$ assert : (iL >= 1 .and. iL <= size(Layers)) : "Wrong layer was specified"

    dom => IN0_FindDomain(io, in0, cZ)
    rP = dom%rHydCond * dom%rAvgHead

    return
  end function rIN0_DefaultPotential


  function rIN0_DomainThickness(io, in0, dom, rHead) result(rH)
    !! function rIN0DomainTransmissivity
    !!
    !! Internal helper function
    !! Computes the transmissivity in domain with ID iDom when the head is 'rHead'
    !!
    !! Calling Sequence:
    !!    cV = cIN0_DomainThickness(in0, iDom, rHead)
    !!
    !! Arguments:
    !!   (in)    type(IN0_COLLECTION), pointer :: in0
    !!             IN0_COLLECTION object to be used
    !!   (in)    integer :: iID
    !!             The ID of the domain to be queried
    !!   (in)    real :: rPot
    !!             The(real) potential at cZ
    !!
    !! Return Value:
    !!           real :: rH
    !!             The thickness
    !!
    ! [ ARGUMENTS ]
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT) :: iDom
    real(kind=AE_REAL), intent(in) :: rHead
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH
    ! [ LOCALS ]
    type(IN0_DOMAIN), pointer :: dom
    real(kind=AE_REAL) :: rMyHead

    !dom => IN0_FindDomainID(io, in0, iDom)
    if (in0%lInitialIteration) then
      rMyHead = dom%rAvgHead
    else
      rMyHead = rHead
    end if

    if (rMyHead < dom%rBase) then
      rH = 1.0e-4_AE_REAL * in0%Domains(1)%rThickness
    else if ((rMyHead - dom%rBase) < dom%rThickness) then
      rH = (rMyHead - dom%rBase)
    else
      rH = dom%rThickness
    end if
    return
  end function rIN0_DomainThickness


  subroutine IN0_Save(io, in0, mode)
    !! Saves the current solution information onto the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(IN0_COLLECTION), pointer :: in0
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istr, ivtx, iflg
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx

    ! Output records will be of the form ELEM_IN0, IWEL, IRAD, IVTX, 0, SIGMA
    do istr = 1, in0%iNStr
      str => in0%Strings(istr)
      do ivtx = 1, str%iNPts
        vtx => str%Vertices(ivtx)
        do iflg = 1, ubound(vtx%rStrength,1)
          if (mode == IO_MODE_BINARY) then
            write (unit=LU_SCRATCH) ELEM_IN0, istr, ivtx, iflg, vtx%rStrength(iflg)
          else
            write (unit=LU_SCRATCH, fmt=*) "IN0", istr, ivtx, iflg, vtx%rStrength(iflg)
          end if
        end do
      end do
    end do

    return
  end subroutine IN0_Save


  subroutine IN0_Load(io, in0, fdp, mode)
    !! Loads the IN0 records from the file on the SCRATCH LU
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(IN0_COLLECTION), pointer :: in0
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: mode
    ! [ LOCALS ]
    integer(kind=AE_INT) :: imodule, istr, ivtx, iflg, istat
    real(kind=AE_REAL) :: rstrength
    character(len=3) :: smodule
    type(IN0_STRING), pointer :: str
    type(IN0_VERTEX), pointer :: vtx
    type(IN0_VERTEX), pointer :: flg

    ! Scans the entire precondition file for the IN0 data
    rewind(unit=LU_SCRATCH)
    do
      if (mode == IO_MODE_BINARY) then
        read (unit=LU_SCRATCH, iostat=istat) imodule, istr, ivtx, iflg, rstrength
        if (imodule /= ELEM_IN0) cycle
      else
        read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, istr, ivtx, iflg, rstrength
        if (uppercase (trim(smodule)) /= "IN0") cycle
      end if
      if (istat < 0) exit
      call IO_Assert(io, istat == 0, "I/O error on precondition file")
      call IO_Assert(io, istr > 0 .and. istr <= in0%iNStr, "IN0 string not found")
      str => in0%Strings(istr)
      call IO_Assert(io, ivtx > 0 .and. ivtx <= str%iNPts, "IN0 vertex not found")
      vtx => str%Vertices(ivtx)
      call IO_Assert(io, iflg > 0 .and. iflg <= ubound(vtx%rStrength,1), "IN0 strength index not found")
      vtx%rStrength(iflg) = rstrength
    end do

    ! Now, populate the internal data structures
    call IN0_Update(io, in0, fdp)

    return
  end subroutine IN0_Load

end module m_in0
