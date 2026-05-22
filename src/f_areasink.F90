module f_areasink

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

  !! module f_areasink (FAS)
  !!
  !! Function module for collections of polygonal area sinks.
  !!
  !! Each FAS_AREASINK holds a closed polygon and a uniform sink density rN.
  !! The "Inside*" functions compute the contributions that arise from the
  !! unbounded area-sink solution only inside the polygon; contributions
  !! outside the polygon are handled by the dipole/well approximation
  !! registered in FDP and FWL during setup (see m_as0.F90).
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   u_io         --  Universal ModAEM I/O functions and constants
  !!   u_polygon    --  Polygon utility routines (PGN_Contains, PGN_ChopSegment, etc.)
  !!   i_areasink   --  Unbounded area-sink influence kernels

  use u_constants
  use u_io
  use u_polygon
  use i_areasink

  implicit none

  public


  type, public :: FAS_AREASINK
    !! type FAS_AREASINK
    !!
    !! Holds geometry and sink density for one polygonal area sink.
    !!
    !! Members:
    !!   complex :: poly(:)
    !!     Closed polygon vertices
    !!   real :: rN
    !!     Uniform sink density [L/T]
    !!   integer :: iNPts
    !!     Number of vertices in use
    !!   integer :: iElementType, iElementString
    !!     Bookkeeping tags identifying the element that registered this sink
    !!
    complex(kind=AE_REAL), dimension(:), pointer :: poly
    real(kind=AE_REAL) :: rN
    integer(kind=AE_INT) :: iNPts
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
  end type FAS_AREASINK


  type, public :: FAS_COLLECTION
    !! type FAS_COLLECTION
    !!
    !! Holds a collection of FAS_AREASINK objects.
    !!
    !! Members:
    !!   type(FAS_AREASINK) :: AreaSinks(:)
    !!     Pre-allocated vector of area sinks
    !!   integer :: iNAS
    !!     Number of area sinks currently in use
    !!
    type(FAS_AREASINK), dimension(:), pointer :: AreaSinks
    integer(kind=AE_INT) :: iNAS
  end type FAS_COLLECTION


contains


  function FAS_Create(io, iMax) result(fas)
    !! function FAS_Create
    !!
    !! Creates a new FAS_COLLECTION with space for iMax area sinks.
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in) :: iMax
    ! [ RETURN VALUE ]
    type(FAS_COLLECTION), pointer :: fas
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, i

    allocate(fas, stat=iStat)
    call IO_Assert(io, (iStat == 0), "FAS_Create: allocation failed")

    allocate(fas%AreaSinks(iMax), stat=iStat)
    call IO_Assert(io, (iStat == 0), "FAS_Create: allocation failed")

    do i = 1, iMax
      nullify(fas%AreaSinks(i)%poly)
      fas%AreaSinks(i)%iNPts = 0
      fas%AreaSinks(i)%rN = rZERO
      fas%AreaSinks(i)%iElementType = -1
      fas%AreaSinks(i)%iElementString = -1
    end do
    fas%iNAS = 0

    return
  end function FAS_Create


  function FAS_New(io, fas, cZ, iNPts, rN, iElementType, iElementString) result(pRV)
    !! function FAS_New
    !!
    !! Registers a new polygon area sink in the collection.
    !!
    !! Calling Sequence:
    !!    pRV => FAS_New(io, fas, cZ, iNPts, rN, iElementType, iElementString)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ(:)
    !!             Polygon vertices (at least iNPts entries)
    !!   (in)    integer :: iNPts
    !!             Number of polygon vertices
    !!   (in)    real :: rN
    !!             Uniform sink density
    !!   (in)    integer :: iElementType, iElementString
    !!             Bookkeeping tags
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    complex(kind=AE_REAL), dimension(:), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iNPts
    real(kind=AE_REAL), intent(in) :: rN
    integer(kind=AE_INT), intent(in) :: iElementType, iElementString
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FAS_AREASINK), pointer :: pRV
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    nullify(pRV)

    if (io%lDebug) then
      call IO_Assert(io, (associated(fas)), "FAS_New: FAS_Create has not been called")
      call IO_Assert(io, (associated(fas%AreaSinks)), "FAS_New: FAS_Create has not been called")
    end if

    if (fas%iNAS >= size(fas%AreaSinks)) return

    fas%iNAS = fas%iNAS + 1
    pRV => fas%AreaSinks(fas%iNAS)
    pRV%iNPts = iNPts
    pRV%rN = rN
    pRV%iElementType = iElementType
    pRV%iElementString = iElementString
    allocate(pRV%poly(iNPts), stat=iStat)
    call IO_Assert(io, (iStat == 0), "FAS_New: Allocation failed")
    pRV%poly = cZ(1:iNPts)

    return
  end function FAS_New


  function rFAS_InsidePotential(io, fas, cZ) result(rPot)
    !! function rFAS_InsidePotential
    !!
    !! Returns the sum of inside-polygon potential contributions at cZ.
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPot
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iAS
    type(FAS_AREASINK), pointer :: fas_sink

    rPot = rZERO
    do iAS = 1, fas%iNAS
      fas_sink => fas%AreaSinks(iAS)
      if (PGN_Contains(fas_sink%poly, cZ)) then
        rP = rIAS_InfluenceP(cZ, fas_sink%poly(1))
        rPot = rPot + fas_sink%rN * rP(1, 1)
      end if
    end do

    return
  end function rFAS_InsidePotential


  function cFAS_InsideDischarge(io, fas, cZ) result(cW)
    !! function cFAS_InsideDischarge
    !!
    !! Returns the sum of inside-polygon discharge contributions at cZ.
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cW
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cP(1, 1)
    integer(kind=AE_INT) :: iAS
    type(FAS_AREASINK), pointer :: fas_sink

    cW = cZERO
    do iAS = 1, fas%iNAS
      fas_sink => fas%AreaSinks(iAS)
      if (PGN_Contains(fas_sink%poly, cZ)) then
        cP = conjg(cIAS_InfluenceW(cZ, fas_sink%poly(1)))
        cW = cW + fas_sink%rN * cP(1, 1)
      end if
    end do

    return
  end function cFAS_InsideDischarge


  function rFAS_InsideRecharge(io, fas, cZ) result(rG)
    !! function rFAS_InsideRecharge
    !!
    !! Returns the sum of inside-polygon recharge contributions at cZ.
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iAS
    type(FAS_AREASINK), pointer :: fas_sink

    rG = rZERO
    do iAS = 1, fas%iNAS
      fas_sink => fas%AreaSinks(iAS)
      if (PGN_Contains(fas_sink%poly, cZ)) then
        rP = rIAS_InfluenceG(cZ, fas_sink%poly(1))
        rG = rG + fas_sink%rN * rP(1, 1)
      end if
    end do

    return
  end function rFAS_InsideRecharge


  function rFAS_InsideFlow(io, fas, cPathZ) result(rF)
    !! function rFAS_InsideFlow
    !!
    !! Returns the total inside-polygon flow across the path cPathZ(:).
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rF
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rP(1, 1)
    integer(kind=AE_INT) :: iZ, iAS, iSeg, iChop
    complex(kind=AE_REAL), dimension(:), allocatable :: cSegChop
    complex(kind=AE_REAL), dimension(2) :: seg
    type(FAS_AREASINK), pointer :: fas_sink

    rF = rZERO
    do iZ = 1, size(cPathZ)-1
      do iAS = 1, fas%iNAS
        fas_sink => fas%AreaSinks(iAS)
        allocate(cSegChop(fas_sink%iNPts+1))
        call PGN_ChopSegment(fas_sink%poly, (/cPathZ(iZ), cPathZ(iZ+1)/), cSegChop, iChop)
        do iSeg = 1, iChop-1
          seg = PGN_Segment(cSegChop, iSeg)
          if (PGN_Contains(fas_sink%poly, rHALF*(seg(1)+seg(2)))) then
            rP = -rIAS_InfluenceF(seg(1), seg(2), fas_sink%poly(1))
            rF = rF + fas_sink%rN * rP(1, 1)
          end if
        end do
        deallocate(cSegChop)
      end do
    end do

    return
  end function rFAS_InsideFlow


  function rFAS_Extraction(io, fas) result(rQ)
    !! function rFAS_Extraction
    !!
    !! Returns the total extraction rate of all area sinks.
    !! (Placeholder: returns zero until active-area tracking is implemented.)
    !!
    ! [ ARGUMENTS ]
    type(FAS_COLLECTION), pointer :: fas
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ

    rQ = rZERO

    return
  end function rFAS_Extraction


end module f_areasink
