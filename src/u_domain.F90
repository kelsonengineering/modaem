module u_domain

  ! ModAEM 1.9
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

  !! module u_domain
  !!
  !! Utility module for polygonal inhomogeneity domains with spatially varying
  !! aquifer properties (base, thickness, hydraulic conductivity, porosity).
  !!
  !! Provides DOM_DOMAIN and DOM_COLLECTION types, domain creation/reading,
  !! geometric containment tests, and aquifer-property utility functions.

  use u_constants
  use u_io
  use u_polygon

  implicit none

  public

  private :: lDOM_PointInsideDomain, &
             lDOM_LineIntersectsDomain, &
             lDOM_DomainInsideDomain

  type, public :: DOM_DOMAIN
    !! type DOM_DOMAIN
    !!
    !! Holds property and geometry information for one inhomogeneity domain.
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
  end type DOM_DOMAIN

  type, public :: DOM_COLLECTION
    !! type DOM_COLLECTION
    !!
    !! Holds a collection of DOM_DOMAIN objects plus metadata.
    !!
    type(DOM_DOMAIN), dimension(:), pointer :: Domains
    integer(kind=AE_INT) :: iNDom
    integer(kind=AE_INT) :: iRegenerate
    logical :: lInitialIteration
    logical :: lPrecondition
  end type DOM_COLLECTION

contains


  function DOM_Create(io, iNDom, rBase, rThick, rHydCond, rPorosity, rAvgHead) result(dom_coll)
    !! Creates a DOM_COLLECTION with space for iNDom domains.
    !! Domain 1 is initialized as the infinite (background) aquifer.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in) :: iNDom
    real(kind=AE_REAL), intent(in) :: rBase
    real(kind=AE_REAL), intent(in) :: rThick
    real(kind=AE_REAL), intent(in) :: rHydCond
    real(kind=AE_REAL), intent(in) :: rPorosity
    real(kind=AE_REAL), intent(in) :: rAvgHead
    ! [ RETURN VALUE ]
    type(DOM_COLLECTION), pointer :: dom_coll
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    type(DOM_DOMAIN), pointer :: dom

    allocate(dom_coll, stat=iStat)
    call IO_Assert(io, (iStat == 0), "DOM_Create: allocation failed")

    allocate(dom_coll%Domains(iNDom), stat=iStat)
    call IO_Assert(io, (iStat == 0), "DOM_Create: allocation failed")
    dom_coll%Domains%iID = -1

    dom_coll%iNDom = 1
    dom_coll%iRegenerate = 1
    dom_coll%lInitialIteration = .true.
    dom_coll%lPrecondition = .true.

    ! Initialize the infinite background domain (domain 1)
    dom => dom_coll%Domains(1)
    dom%rBase = rBase
    dom%rThickness = rThick
    dom%rHydCond = rHydCond
    dom%rPorosity = rPorosity
    dom%rAvgHead = rAvgHead
    dom%rTopPot = rHALF * rHydCond * rThick**2
    dom%iInsideDomain = 0
    dom%iOutsideDomain = 0
    dom%iNPts = 0
    dom%iID = 0
    nullify(dom%cZ)

    return
  end function DOM_Create


  subroutine DOM_Read(io, dom_coll)
    !! Reads one DOM record from the current input buffer and adds a new domain
    !! to dom_coll.  The caller is responsible for subsequently reading vertex
    !! coordinate records and appending them to dom_coll%Domains(dom_coll%iNDom)%cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iMax, iID, iStat
    real(kind=AE_REAL) :: rBase, rThickness, rHydCond, rPorosity, rAvgHead
    type(DOM_DOMAIN), pointer :: dom

    call IO_Assert(io, (associated(dom_coll)), "DOM_Read: DOM_Create has not been called")
    call IO_Assert(io, (associated(dom_coll%Domains)), &
         "DOM_Read: No domains have been allocated")
    call IO_Assert(io, (dom_coll%iNDom < size(dom_coll%Domains)), &
         "DOM_Read: Space exhausted")

    iMax = iIO_GetInteger(io, 'iMax', minimum=3)
    rBase = rIO_GetReal(io, 'rBase')
    rThickness = rIO_GetReal(io, 'rThickness', minimum=rTINY)
    rHydCond = rIO_GetReal(io, 'rHydCond', minimum=rTINY)
    rPorosity = rIO_GetReal(io, 'rPorosity', minimum=rTINY)
    rAvgHead = rIO_GetReal(io, 'rAvgHead', minimum=rBase)
    iID = iIO_GetInteger(io, 'iID', forbidden=dom_coll%Domains%iID)

    dom_coll%iNDom = dom_coll%iNDom + 1
    dom => dom_coll%Domains(dom_coll%iNDom)
    dom%iNPts = 0
    dom%rBase = rBase
    dom%rThickness = rThickness
    dom%rHydCond = rHydCond
    dom%rPorosity = rPorosity
    dom%rAvgHead = rAvgHead
    dom%iID = iID
    dom%iInsideDomain = iID
    dom%iOutsideDomain = 0
    allocate(dom%cZ(iMax), stat=iStat)
    call IO_Assert(io, (iStat == 0), "DOM_Read: Allocation failed")
    dom%cZ = cZERO

    return
  end subroutine DOM_Read


  subroutine DOM_NewDomain(io, dom_coll, cZ, iNPts, rBase, rThickness, rHydCond, rPorosity)
    !! Adds a new domain to dom_coll from caller-supplied data (no io reading).
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), dimension(:), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iNPts
    real(kind=AE_REAL), intent(in) :: rBase, rThickness, rHydCond, rPorosity
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iDom
    type(DOM_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(dom_coll)), &
           "DOM_NewDomain: DOM_Create has not been called")
    end if

    call IO_Assert(io, (dom_coll%iNDom < size(dom_coll%Domains)), &
         "DOM_NewDomain: Space exhausted")

    dom_coll%iNDom = dom_coll%iNDom + 1
    dom => dom_coll%Domains(dom_coll%iNDom)
    allocate(dom%cZ(iNPts), stat=iStat)
    call IO_Assert(io, (iStat == 0), "DOM_NewDomain: Allocation failed")
    dom%cZ = cZ(1:iNPts)
    dom%iNPts = iNPts
    dom%rBase = rBase
    dom%rThickness = rThickness
    dom%rHydCond = rHydCond
    dom%rPorosity = rPorosity

    do iDom = 2, dom_coll%iNDom - 1
      call IO_Assert(io, (.not. lDOM_DomainOverlapsDomain(io, dom_coll, iDom, dom_coll%iNDom)), &
           "DOM_NewDomain: New domain overlaps another domain")
    end do

    return
  end subroutine DOM_NewDomain


  subroutine DOM_PreSolve(io, dom_coll)
    !! Orients all domain polygons positively and resolves nesting topology.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom, iDom2, iNChanges
    type(DOM_DOMAIN), pointer :: dom, dom2, outside

    ! Ensure positive orientation for all non-infinite domains
    do iDom = 1, dom_coll%iNDom
      dom => dom_coll%Domains(iDom)
      if (iDom > 1) then
        call PGN_MakePositivelyOriented(dom%cZ(1:dom%iNPts))
      end if
    end do

    ! Compute top-of-aquifer potentials
    dom_coll%Domains(:)%rTopPot = rHALF * dom_coll%Domains(:)%rHydCond * &
                                   dom_coll%Domains(:)%rThickness**2

    ! Resolve nested containment iteratively
    print *, 'Checking for nested inhomogeneities...'
    do
      iNChanges = 0
      do iDom = 2, dom_coll%iNDom
        dom => dom_coll%Domains(iDom)
        do iDom2 = 2, dom_coll%iNDom
          dom2 => dom_coll%Domains(iDom2)
          if (iDom /= iDom2) then
            if (lDOM_DomainInsideDomain(io, dom_coll, iDOM_DomainIDIndex(io, dom_coll, dom%iOutsideDomain), iDom2) .and. &
                lDOM_DomainInsideDomain(io, dom_coll, iDom2, iDom)) then
              dom%iOutsideDomain = dom2%iID
              iNChanges = iNChanges + 1
            end if
          end if
        end do
        outside => DOM_FindDomainID(io, dom_coll, dom%iOutsideDomain)
      end do
      if (iNChanges == 0) exit
    end do

    return
  end subroutine DOM_PreSolve


  function DOM_FindDomain(io, dom_coll, cZ) result(dom)
    !! Returns a pointer to the domain containing the point cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    type(DOM_DOMAIN), pointer :: dom
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    dom => dom_coll%Domains(1)
    do iDom = 2, dom_coll%iNDom
      if (lDOM_PointInsideDomain(io, dom_coll, iDom, cZ)) then
        dom => dom_coll%Domains(iDom)
      end if
    end do

    return
  end function DOM_FindDomain


  function DOM_FindDomainID(io, dom_coll, iID) result(dom)
    !! Returns a pointer to the domain with the given ID.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iID
    ! [ RETURN VALUE ]
    type(DOM_DOMAIN), pointer :: dom
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    dom => dom_coll%Domains(1)
    do iDom = 2, dom_coll%iNDom
      if (dom_coll%Domains(iDom)%iID == iID) then
        dom => dom_coll%Domains(iDom)
      end if
    end do

    return
  end function DOM_FindDomainID


  function iDOM_DomainIDIndex(io, dom_coll, iID) result(iRes)
    !! Returns the array index of the domain with the given ID.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iID
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iRes
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDom

    iRes = 1
    do iDom = 1, dom_coll%iNDom
      if (dom_coll%Domains(iDom)%iID == iID) then
        iRes = iDom
        return
      end if
    end do
    call IO_Assert(io, .false., 'iDOM_DomainIDIndex: Could not find ID')

    return
  end function iDOM_DomainIDIndex


  function rDOM_HeadToPotential(io, dom_coll, rHead, cZ) result(rPot)
    !! Converts head to discharge potential at the domain containing cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    real(kind=AE_REAL), intent(in) :: rHead
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPot
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHd
    type(DOM_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(dom_coll)), &
           "rDOM_HeadToPotential: DOM_Create has not been called")
    end if

    dom => DOM_FindDomain(io, dom_coll, cZ)

    rHd = rHead - dom%rBase
    if (rHd > dom%rThickness) then
      rPot = dom%rHydCond * dom%rThickness * rHd - rHALF * dom%rHydCond * dom%rThickness**2
    else
      rPot = rHALF * dom%rHydCond * rHd**2
    end if

    return
  end function rDOM_HeadToPotential


  function rDOM_PotentialToHead(io, dom_coll, rPot, cZ, lTest) result(rHead)
    !! Converts discharge potential to head at the domain containing cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cZ
    logical, intent(in), optional :: lTest
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHead
    ! [ LOCALS ]
    type(DOM_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(dom_coll)), &
           "rDOM_PotentialToHead: DOM_Create has not been called")
    end if

    dom => DOM_FindDomain(io, dom_coll, cZ)

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
  end function rDOM_PotentialToHead


  function cDOM_DischargeToVelocity(io, dom_coll, cDischarge, cZ, rPot) result(cVelocity)
    !! Converts specific discharge to seepage velocity at the domain containing cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    real(kind=AE_REAL), intent(in) :: rPot
    complex(kind=AE_REAL), intent(in) :: cDischarge
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cVelocity
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rSatdThick
    type(DOM_DOMAIN), pointer :: dom

    if (io%lDebug) then
      call IO_Assert(io, (associated(dom_coll)), &
           "cDOM_DischargeToVelocity: DOM_Create has not been called")
    end if

    dom => DOM_FindDomain(io, dom_coll, cZ)

    if (rPot < rZERO) then
      cVelocity = cZERO
      return
    else
      if (rPot > dom%rTopPot) then
        rSatdThick = dom%rThickness
      else
        rSatdThick = rDOM_PotentialToHead(io, dom_coll, rPot, cZ) - dom%rBase
      end if
      cVelocity = cDischarge / (rSatdThick * dom%rPorosity)
    end if

    return
  end function cDOM_DischargeToVelocity


  function lDOM_IsConfined(io, dom_coll, cZ, rPot) result(lConfined)
    !! Returns .true. if flow at cZ is confined (head above aquifer top).
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    ! [ RETURN VALUE ]
    logical :: lConfined
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)

    if (dom_coll%lPrecondition) then
      rHead = dom%rAvgHead
    else
      rHead = rDOM_PotentialToHead(io, dom_coll, rPot, cZ)
    end if

    lConfined = ((rHead - dom%rBase) >= dom%rThickness)

    return
  end function lDOM_IsConfined


  function rDOM_Base(io, dom_coll, cZ) result(rBase)
    !! Returns the aquifer base elevation at cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBase
    ! [ LOCALS ]
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)
    rBase = dom%rBase

    return
  end function rDOM_Base


  function rDOM_SatdThickness(io, dom_coll, cZ, rPot) result(rH)
    !! Returns the saturated thickness at cZ where the potential is rPot.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)

    if (dom_coll%lPrecondition) then
      if (dom%rAvgHead > dom%rBase + dom%rThickness) then
        rHead = dom%rAvgHead + dom%rBase
      else
        rHead = dom%rAvgHead
      end if
    else
      rHead = rDOM_PotentialToHead(io, dom_coll, rPot, cZ)
    end if

    if ((rHead - dom%rBase) < dom%rThickness) then
      rH = rHead - dom%rBase
    else
      rH = dom%rThickness
    end if

    return
  end function rDOM_SatdThickness


  function rDOM_DefaultPotential(io, dom_coll, cZ) result(rP)
    !! Returns the default (initial-head) potential at cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rP
    ! [ LOCALS ]
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)
    rP = dom%rHydCond * dom%rAvgHead

    return
  end function rDOM_DefaultPotential


  function rDOM_Transmissivity(io, dom_coll, cZ, rPot) result(rT)
    !! Returns the transmissivity at cZ where the potential is rPot.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rPot
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rT
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rHead
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)

    if (dom_coll%lPrecondition) then
      rHead = dom%rAvgHead
    else
      rHead = rDOM_PotentialToHead(io, dom_coll, rPot, cZ)
    end if

    if ((rHead - dom%rBase) < rONE_TENTH * dom%rThickness) then
      rT = rONE_TENTH * dom%rThickness
    else if ((rHead - dom%rBase) < dom%rThickness) then
      rT = dom%rHydCond * (rHead - dom%rBase)
    else
      rT = dom%rHydCond * dom%rThickness
    end if

    return
  end function rDOM_Transmissivity


  function rDOM_HydCond(io, dom_coll, cZ) result(rHydCond)
    !! Returns the hydraulic conductivity at cZ.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rHydCond
    ! [ LOCALS ]
    type(DOM_DOMAIN), pointer :: dom

    dom => DOM_FindDomain(io, dom_coll, cZ)
    rHydCond = dom%rHydCond

    return
  end function rDOM_HydCond


  function rDOM_DomainThickness(io, dom_coll, dom, rHead) result(rH)
    !! Computes saturated thickness in dom given rHead (uses avg head on first iteration).
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    type(DOM_DOMAIN), pointer :: dom
    real(kind=AE_REAL), intent(in) :: rHead
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rH
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rMyHead

    if (dom_coll%lInitialIteration) then
      rMyHead = dom%rAvgHead
    else
      rMyHead = rHead
    end if

    if (rMyHead < dom%rBase) then
      rH = 1.0e-4_AE_REAL * dom_coll%Domains(1)%rThickness
    else if ((rMyHead - dom%rBase) < dom%rThickness) then
      rH = (rMyHead - dom%rBase)
    else
      rH = dom%rThickness
    end if

    return
  end function rDOM_DomainThickness


  ! ---------------------------------------------------------------------------
  ! Private geometry functions
  ! ---------------------------------------------------------------------------

  function lDOM_PointInsideDomain(io, dom_coll, iDomain, cZ) result(lInside)
    !! Tests whether cZ is inside domain iDomain.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iDomain
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    logical :: lInside
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, ii
    real(kind=AE_REAL) :: rSum
    type(DOM_DOMAIN), pointer :: dom

    dom => dom_coll%Domains(iDomain)
    if (size(dom%cZ) == 0) then
      lInside = .true.
    else
      if (any((cZ == dom%cZ(1:dom%iNPts)), 1)) then
        lInside = .true.
      else
        rSum = rZERO
        ii = dom%iNPts
        do i = 1, dom%iNPts
          rSum = rSum + aimag(log((cZ - dom%cZ(i)) / (cZ - dom%cZ(ii))))
          ii = i
        end do
        lInside = (rSum > rONE)
      end if
    end if

    return
  end function lDOM_PointInsideDomain


  function lDOM_LineIntersectsDomain(io, dom_coll, iDomain, cZ1, cZ2) result(lIntersects)
    !! Tests whether the segment cZ1-cZ2 intersects the boundary of domain iDomain.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iDomain
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    ! [ RETURN VALUE ]
    logical :: lIntersects
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    complex(kind=AE_REAL) :: cMapZ1, cMapZ2
    real(kind=AE_REAL) :: rXInt
    type(DOM_DOMAIN), pointer :: dom

    dom => dom_coll%Domains(iDomain)
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
  end function lDOM_LineIntersectsDomain


  function lDOM_DomainInsideDomain(io, dom_coll, iDomain, iDomain2) result(lInside)
    !! Tests whether domain iDomain2 is entirely inside domain iDomain.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iDomain
    integer(kind=AE_INT), intent(in) :: iDomain2
    ! [ RETURN VALUE ]
    logical :: lInside
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    type(DOM_DOMAIN), pointer :: dom, dom2

    dom => dom_coll%Domains(iDomain)
    dom2 => dom_coll%Domains(iDomain2)
    if (size(dom%cZ, 1) == 0) then
      lInside = .true.
    else
      if (size(dom2%cZ) == 0) then
        lInside = .false.
      else
        lInside = .true.
        j = dom2%iNPts
        do i = 1, dom2%iNPts
          if (.not. lDOM_PointInsideDomain(io, dom_coll, iDomain, dom2%cZ(i)) .or. &
              lDOM_LineIntersectsDomain(io, dom_coll, iDomain, dom2%cZ(i), dom2%cZ(j))) then
            lInside = .false.
            return
          end if
          j = i
        end do
      end if
    end if

    return
  end function lDOM_DomainInsideDomain


  function lDOM_DomainOverlapsDomain(io, dom_coll, iDomain, iDomain2) result(lOverlap)
    !! Tests whether domain iDomain2 overlaps domain iDomain.
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(DOM_COLLECTION), pointer :: dom_coll
    integer(kind=AE_INT), intent(in) :: iDomain
    integer(kind=AE_INT), intent(in) :: iDomain2
    ! [ RETURN VALUE ]
    logical :: lOverlap
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j
    type(DOM_DOMAIN), pointer :: dom, dom2

    dom => dom_coll%Domains(iDomain)
    dom2 => dom_coll%Domains(iDomain2)
    if (size(dom%cZ, 1) == 0) then
      lOverlap = .true.
    else
      if (size(dom2%cZ) == 0) then
        lOverlap = .false.
      else
        lOverlap = .false.
        j = size(dom2%cZ, 1)
        do i = 1, size(dom2%cZ, 1)
          if (lDOM_LineIntersectsDomain(io, dom_coll, iDomain, dom2%cZ(i), dom2%cZ(j))) then
            lOverlap = .true.
            return
          end if
          j = i
        end do
      end if
    end if

    return
  end function lDOM_DomainOverlapsDomain


end module u_domain
