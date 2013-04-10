module f_dipole

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

  !! module f_dipole(FDP)
  !!
  !! Written by Victor A.Kelson
  !!
  !! Revision History:
  !!   1.0.0   17 March 1999
  !!           First "source code release" version, adapted from
  !!           previous testing version.
  !!
  !! Module of data structures and functions for second order
  !! dipoles of complex strength.  The real part of the strength
  !! is the 'dipole' strength(jump in the streamfunction) and the
  !! imaginary part is the 'doublet' strength(jump in the
  !! potential)
  !!
  !! This module encapsulates containers of dipole functions
  !! in layers.  The various data structures are organized for
  !! efficient computation on vector or parallel machinery.
  !!
  !! Limitations:
  !!   At present, this module has no hooks for series expansion
  !!   block implementation.  No known computational bugs exist.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   io         --  Universal ModAEM I/O functions and constants
  !!   i_dipole   --  Influence function module for dipoles

  use u_constants
  use u_io
  use i_dipole
  use i_linesink

  implicit none

  public


  type, public :: FDP_DIPOLE
    !! type FDP_DIPOLE
    !!
    !! Type that holds information for one dipole
    !!
    !! Members:
    !!   complex :: cZC, cZL
    !!     The center and directed length vector for the dipole
    !!   complex :: cRho(3)
    !!     The complex strength coefficients for the dipole
    !!     cRho(1) is the strength at the first end, cRho(2) is the
    !!     strength at the center and cRho(3) is the strength at
    !!     the second end.
    !!
    complex(kind=AE_REAL) :: cZC
    complex(kind=AE_REAL) :: cZL
    complex(kind=AE_REAL), dimension(3) :: cRho
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
  end type FDP_DIPOLE


  type, public :: FDP_COLLECTION
    !! type FDP_COLLECTION
    !!
    !! Type that holds the dipole entries for a layer
    !!
    !! Members:
    !!   type(FDP_DIPOLE), pointer :: Dipoles(:)
    !!     FDP_DIPOLE structures that hold dipole information
    !!   integer :: iCount
    !!     Number of FDP_DIPOLE structures currently in use
    !!
    type(FDP_DIPOLE), dimension(:), pointer :: Dipoles
    integer(kind=AE_INT) :: iCount
  end type FDP_COLLECTION


contains


  function FDP_Create(io, iMax) result(fdp)
    !! function FDP_Create
    !!
    !! Creates a new FDP_COLLECTION object
    !!
    !! Calling Sequence:
    !!    fdp => FDP_Create(iMax)
    !!
    !! Arguments:
    !!    (in)      integer :: iMax
    !!                The maximum number of dipoles to be stored
    !!
    !! Return Value:
    !!   On success, fdp points to a new FDP_COLLECTION object
    !!   On failure(allocation error), fatal error
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FDP_COLLECTION), pointer :: fdp
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(fdp, stat = iStat)
    call IO_Assert(io, (iStat == 0), "FDP_Create: allocation failed")

    allocate(fdp%Dipoles(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FDP_Create: allocation failed")

    fdp%Dipoles = FDP_DIPOLE(cZERO, cmplx(rONE, rONE, AE_REAL), (/cZERO, cZERO, cZERO/), -1, -1, -1, -1)
    fdp%iCount = 0

    return
  end function FDP_Create


  subroutine FDP_New(io, fdp, cZ1, cZ2, cRho, iElementType, iElementString, iElementVertex, iElementFlag, iRV)
    !! subroutine FDP_New
    !!
    !! Makes a new dipole entry. On call, the geometry and strengths of
    !! the dipole are provided. The internal dipole structures are then
    !! set up. Returns iRV = the index into the dipole table on success or
    !! iRV < 0 on failure.
    !!
    !! Calling Sequence:
    !!    call FDP_New(io, fdp, cZ1, cZ2, cRho, iElementType, iElementString, iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    complex :: cZ1, cZ2
    !!             Complex coordinates of the ends of the dipole
    !!   (in)    complex :: cRho(3)
    !!             The complex strengths at cZ1, at the center, and
    !!             at cZ2, respectively.
    !!   (out)   integer :: iRV
    !!             On success, the index in fdp%Dipoles used
    !!
    !! Note: On failure, forces a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), dimension(:), intent(in) :: cRho
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(out) :: iRV
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_New: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_New: FDP_Alloc has not been called")
      call IO_Assert(io, (fdp%iCount < size(fdp%Dipoles)), "FDP_New: Space exhausted")
      call IO_Assert(io, (size(cRho) == 3), "FDP_New: Illegal strength vector")
      !! Note: the next assertion ensures that the caller has checked the length appropriately
      call IO_Assert(io, (abs(cZ2-cZ1) > rTINY), "FDP_New: Dipole length is too short")
    end if

    ! Compute the center and directed length parameters
    fdp%iCount = fdp%iCount + 1
    dip => fdp%Dipoles(fdp%iCount)
    dip = FDP_DIPOLE(rHALF*(cZ2+cZ1), rHALF*(cZ2-cZ1), cRho, iElementType, iElementString, iElementVertex, iElementFlag)
    iRV = fdp%iCount

    return
  end subroutine FDP_New


  subroutine FDP_Update(io, fdp, iDP, cRho)
    !! subroutine FDP_Update
    !!
    !! Updates the strength coefficients for a dipole entry.
    !!
    !! Calling Sequence:
    !!    call FDP_Update(io, fdp, iDP, cRho, iRV)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to be used
    !!   (in)    integer :: iDP
    !!             The index for the dipole
    !!   (in)    complex :: cRho(3)
    !!             The complex strengths at cZ1, at the center, and
    !!             at cZ2, respectively.
    !!
    !! Note: On failure, forces a fatal error
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: iDP
    complex(kind=AE_REAL), dimension(3), intent(in) :: cRho(3)
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
      call IO_Assert(io, (iDP <= size(fdp%Dipoles)), "FDP_Update: Space exhausted")
      call IO_Assert(io, (size(cRho) == 3), "FDP_Update: Illegal strength vector")
    end if

    dip => fdp%Dipoles(iDP)
    dip%cRho = cRho

    return
  end subroutine FDP_Update


  subroutine FDP_GetInfluence_IDP(io, fdp, iWhich, iDP1, iNDP, cPathZ, cOrientation, cF)
    !! subroutine FDP_GetInfluence
    !!
    !! Retrieves arrays of influence functions for use in matrix generation
    !!
    !! Calling Sequence:
    !!    call FDP_GetInfluence(io, fdp, iWhich, iDP1, iNDP, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION) :: fdp
    !!             The FDP_COLLECTION to be used
    !!   (in)    integer :: iWhich
    !!             The influence function to be computed;  iWhich values are
    !!                INFLUENCE_P   - Complex potential
    !!                INFLUENCE_Q   - Complex discharge
    !!                INFLUENCE_F   - Integrated flux
    !!                INFLUENCE_F_NOBC   - Integrated flux(do not compute branch cut)
    !!                INFLUENCE_G   - Areal infiltration
    !!                INFLUENCE_Q   - Extraction rate
    !!                INFLUENCE_J   - Jump magnitude
    !!                INFLUENCE_D   - Difference in potential
    !!                INFLUENCE_Z   - All zeroes
    !!   (in)    integer :: iDP1
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP
    !!             The number of consecutive dipoles to be computed
    !!   (in)    complex :: cPathZ(:)
    !!             Complex coordinates of the control point(s) to be used. For
    !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
    !!             cPathZ(1) is used; for iWhich = INFLUENCE_F, the influence
    !!             function is computed along the path cPathZ(:)
    !!   (in)    complex :: cOrientation
    !!             Orientation normal vector for iWhich = INFLUENCE_W
    !!   (out)   complex :: cF(1:iNDP, 3)
    !!             The returned influence functions.  Indexes 1:iNDP relate
    !!             to dipole indices iDP1:iDP1+iNDP-1, respectively.  Elements
    !!             cF(:, 1), cF(:, 2), cF(:, 3) are the coefficients for the first
    !!             end, center and second end of the dipole, respectively.
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: iWhich, iDP1, iNDP
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDP2, i, j, iStat
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ2
    real(kind=AE_REAL), dimension(:), allocatable :: rX
    complex(kind=AE_REAL) :: cUnit
    complex(kind=AE_REAL), dimension(3, 1) :: cP, cR, cS, cG, cQ, cJ
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_GetInfluence_IDP: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_GetInfluence_IDP: FDP_Alloc has not been called")
      call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
           "FDP_GetInfluence_IDP: Bad index range")
      call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
           "FDP_GetInfluence_IDP: Bad index range")
      call IO_Assert(io, (size(cF, 1) >= iNDP .and. size(cF, 2) >= 3 .and. size(cF, 3) >= 1), &
           "FDP_GetInfluence_IDP: Invalid result vector")
    end if

    ! It is assumed that the caller has eliminated the potential for a singularity
    ! by selecting appropriate control points
    iDP2 = iDP1+iNDP-1
    allocate(cMapZ1(iDP1:iDP2), cMapZ2(iDP1:iDP2), rX(iDP1:iDP2), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FDP_GetInfluence_IDP: Allocation failed")
    !
    select case (iWhich)
      case (INFLUENCE_P)      ! Potential
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          cP = cIDP_InfluenceP(cMapZ1(iDP1+i-1))
          cF(i, :, 1) = cP(:, 1)
        end do
      case (INFLUENCE_D)      ! Potential difference
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        cMapZ2(:) = (cPathZ(2)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          cP = cIDP_InfluenceP(cMapZ1(iDP1+i-1)) - cIDP_InfluenceP(cMapZ2(iDP1+i-1))
          cF(i, :, 1) = cP(:, 1)
        end do
      case (INFLUENCE_W)      ! Discharge
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          cUnit = cOrientation/abs(cOrientation)
          cR = cIDP_InfluenceW(cMapZ1(iDP1+i-1))
          ! Rotate the discharge vector by the unit vector in the direction of caOrientation
          ! and reverse the real and imaginary parts for use in influence functions...
          cF(i, :, 1) = (cUnit * (cR(:, 1)/fdp%Dipoles(iDP1+i-1)%cZL))
        end do
      case (INFLUENCE_F)      ! Integrated Flow
        cF = cZERO
        do j = 1, size(cPathZ)-1
          cMapZ1(:) = (cPathZ(j)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
          cMapZ2(:) = (cPathZ(j+1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
          do i = 1, iNDP
            cS = cIDP_InfluenceF(cMapZ1(iDP1+i-1), cMapZ2(iDP1+i-1))
            cF(i, :, 1) = cF(i, :, 1) + cS(:, 1)
          end do
        end do
      case (INFLUENCE_G)      ! Recharge rate
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          cG = cIDP_InfluenceG(cMapZ1(iDP1+i-1))
          cF(i, :, 1) = cG(:, 1)
        end do
      case (INFLUENCE_Q)      ! Extraction rate
        do i = 1, iNDP
          cQ = cIDP_InfluenceQ()
          cF(i, :, 1) = cQ(:, 1)
        end do
      case (INFLUENCE_J)      ! Jump
        rX(:) = real((cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL)
        do i = 1, iNDP
          cJ = cIDP_InfluenceJ(rX(iDP1+i-1))
          cF(i, :, 1) = cJ(:, 1)
        end do
      case (INFLUENCE_Z)
        cF = cZERO
    end select
    deallocate(cMapZ1, cMapZ2, rX)

    return
  end subroutine FDP_GetInfluence_IDP


  subroutine FDP_GetInfluence_ILS(io, fdp, iWhich, iDP1, iNDP, cPathZ, cOrientation, cF)
    !! subroutine FDP_GetInfluence_ILS
    !!
    !! Retrieves arrays of influence functions for use in matrix generation,
    !! using the first-order linesink function, instead of the dipole functions.
    !!
    !! Calling Sequence:
    !!    call FDP_GetInfluence(io, fdp, iWhich, iDP1, iNDP, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    integer :: iWhich
    !!             The influence function to be computed;  iWhich values are
    !!                INFLUENCE_P   - Complex potential
    !!                INFLUENCE_Q   - Complex discharge
    !!                INFLUENCE_F   - Integrated flux
    !!                INFLUENCE_G   - Areal infiltration
    !!                INFLUENCE_Q   - Extraction rate
    !!                INFLUENCE_D   - Difference in potential
    !!                INFLUENCE_Z   - All zeroes
    !!   (in)    integer :: iDP1
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP
    !!             The number of consecutive dipoles to be computed
    !!   (in)    complex :: cPathZ(:)
    !!             Complex coordinates of the control point(s) to be used. For
    !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
    !!             cPathZ(1) is used; for iWhich = INFLUENCE_F, the influence
    !!             function is computed along the path cPathZ(:)
    !!   (in)    complex :: cOrientation
    !!             Orientation normal vector for iWhich = INFLUENCE_W
    !!   (out)   complex :: cF(1:iNDP, 1)
    !!             The returned influence functions.  Indexes 1:iNDP relate
    !!             to dipole indices iDP1:iDP1+iNDP-1, respectively.  The
    !!             influence function is in terms of the sink density of the
    !!             linesink.
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: iWhich, iDP1, iNDP
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iDP2, i, j, iStat
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ2
    complex(kind=AE_REAL) :: cUnit
    complex(kind=AE_REAL), dimension(1, 1) :: cP, cR, cS, cG, cQ
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_GetInfluence_ILS: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_GetInfluence_ILS: FDP_Alloc has not been called")
      call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
           "FDP_GetInfluence_ILS: Bad index range")
      call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and.(iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
           "FDP_GetInfluence_ILS: Bad index range")
      call IO_Assert(io, (size(cF, 1) >= iNDP .and. size(cF, 2) >= 1 .and. size(cF, 3) >= 1), &
           "FDP_GetInfluence_ILS: Invalid result vector")
    end if

    ! Allocate the Map-Z arrays
    iDP2 = iDP1+iNDP-1
    allocate(cMapZ1(iDP1:iDP2), cMapZ2(iDP1:iDP2), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FDP_GetInfluence_ILS: Allocation failed")

    ! It is assumed that the caller has eliminated the potential for a singularity
    ! by selecting appropriate control points
    select case (iWhich)
      case (INFLUENCE_P)
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          dip => fdp%Dipoles(iDP1+i-1)
          cP = cILS_InfluenceP(cMapZ1(iDP1+i-1), dip%cZL)
          cF(i, :, 1) = cP(:, 1)
        end do
      case (INFLUENCE_D)
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        cMapZ2(:) = (cPathZ(2)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          dip => fdp%Dipoles(iDP1+i-1)
          cP = cILS_InfluenceP(cMapZ1(iDP1+i-1), dip%cZL) - cILS_InfluenceP(cMapZ2(iDP1+i-1), dip%cZL)
          cF(i, :, 1) = cP(:, 1)
        end do
      case (INFLUENCE_W)
        cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i = 1, iNDP
          dip => fdp%Dipoles(iDP1+i-1)
          ! Rotate the discharge vector by the unit vector in the direction of cOrientation
          cUnit = cOrientation/abs(cOrientation)
          cR = cILS_InfluenceW(cMapZ1(iDP1+i-1), dip%cZL)
          cF(i, :, 1) = cUnit * (cR(:, 1)/fdp%Dipoles(iDP1+i-1)%cZL)
        end do
      case (INFLUENCE_F)
        cF = cZERO
        do j = 1, size(cPathZ)-1
          cMapZ1(:) = (cPathZ(j)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
          cMapZ2(:) = (cPathZ(j+1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
          do i = 1, iNDP
            dip => fdp%Dipoles(iDP1+i-1)
            cS = cILS_InfluenceF(cMapZ1(iDP1+i-1), cMapZ2(iDP1+i-1), dip%cZL)
            cF(i, :, 1) = cF(i, :, 1) + cS(:, 1)
          end do
        end do
      case (INFLUENCE_G)
        do i = 1, iNDP
          cG = cILS_InfluenceG()
          cF(i, :, 1) = cG(:, 1)
        end do
      case (INFLUENCE_Q)
        do i = 1, iNDP
          dip => fdp%Dipoles(iDP1+i-1)
          cQ = cILS_InfluenceQ(dip%cZL)
          cF(i, :, 1) = cQ(:, 1)
        end do
      case (INFLUENCE_Z)
        cF = cZERO
    end select

    deallocate(cMapZ1, cMapZ2)

    return
  end subroutine FDP_GetInfluence_ILS


  function cFDP_Potential(io, fdp, cZ, iDP1, iNDP) result(cOmega)
    !! complex function cFDP_Potential
    !!
    !! Computes the complex potential due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    cOmega = cFDP_Potential(fdp, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iDP1 [OPTIONAL]
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP [OPTIONAL]
    !!             The number of consecutive dipoles to be computed
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iDP1, iNDP
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(3, 1) :: cP
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iDP1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
             "FDP_GetInfluence_ILS: Bad index range")
      end if
      iStart = iDP1
      if (present(iNDP)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
               "FDP_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNDP-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fdp%iCount
    end if

    if (fdp%iCount > 0) then
      ! Sum up the contribution of all dipoles
      allocate(cMapZ(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "cFDP_Potential: Allocation failed")
      cOmega = cZERO
      cMapZ = (cZ - fdp%Dipoles%cZC) / fdp%Dipoles%cZL
      do i = iStart, iEnd
        dip => fdp%Dipoles(i)
        cP = cIDP_InfluenceP(cMapZ(i))
        cOmega = cOmega + sum(dip%cRho*cP(:, 1))
      end do
      deallocate(cMapZ)
    else
      cOmega = cZERO
    end if

    return
  end function cFDP_Potential


  function cFDP_Discharge(io, fdp, cZ, iDP1, iNDP) result(cQ)
    !! complex function cFDP_Discharge
    !!
    !! Computes the complex discharge due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    cQ = cFDP_Discharge(fdp, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iDP1 [OPTIONAL]
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP [OPTIONAL]
    !!             The number of consecutive dipoles to be computed
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iDP1, iNDP
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(3, 1) :: cP
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iDP1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
             "FDP_GetInfluence_ILS: Bad index range")
      end if
      iStart = iDP1
      if (present(iNDP)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
               "FDP_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNDP-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fdp%iCount
    end if

    ! Sum up the contribution of all dipoles
    if (fdp%iCount > 0) then
      allocate(cMapZ(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "cFDP_Potential: Allocation failed")
      cQ = cZERO
      cMapZ = (cZ - fdp%Dipoles%cZC) / fdp%Dipoles%cZL
      do i = iStart, iEnd
        dip => fdp%Dipoles(i)
        cP = cIDP_InfluenceW(cMapZ(i))
        cQ = cQ + conjg(sum(dip%cRho*cP(:, 1) / dip%cZL))
      end do
      deallocate(cMapZ)
    else
      cQ = cZERO
    end if

    return
  end function cFDP_Discharge


  function rFDP_Flow(io, fdp, cPathZ, iDP1, iNDP) result(rFlow)
    !! real function cFDP_Flow
    !!
    !! Computes the integrated flow due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    rFlow = fFDP_Flow(fdp, cPathZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!   (in)    integer :: iDP1 [OPTIONAL]
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP [OPTIONAL]
    !!             The number of consecutive dipoles to be computed
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in), optional :: iDP1, iNDP
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(3, 1) :: cS
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1, cMapZ2
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip


    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iDP1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
             "FDP_GetInfluence_ILS: Bad index range")
      end if
      iStart = iDP1
      if (present(iNDP)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
               "FDP_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNDP-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fdp%iCount
    end if

    if (fdp%iCount > 0) then
      ! Sum up the contribution of all dipoles
      allocate(cMapZ1(iStart:iEnd), cMapZ2(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "cFDP_Potential: Allocation failed")
      rFlow = rZERO
      do j = 1, ubound(cPathZ, 1)-1
        cMapZ1 = (cPathZ(j)-fdp%Dipoles%cZC) / fdp%Dipoles%cZL
        cMapZ2 = (cPathZ(j+1)-fdp%Dipoles%cZC) / fdp%Dipoles%cZL
        do i = iStart, iEnd
          dip => fdp%Dipoles(i)
          cS = cIDP_InfluenceF(cMapZ1(i), cMapZ2(i))
          rFlow = rFlow + sum(real(dip%cRho * cS(:, 1)))
        end do
      end do
      deallocate(cMapZ1, cMapZ2)
    else
      rFlow = rZERO
    end if

    return
  end function rFDP_Flow


  function rFDP_Recharge(io, fdp, cZ, iDP1, iNDP) result(rGamma)
    !! complex function cFDP_Recharge
    !!
    !! Computes the recharge rate due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    rGamma = cFDP_Recharge(fdp, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    complex :: cZ
    !!             The point at which to determine the recharge
    !!   (in)    integer :: iDP1 [OPTIONAL]
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP [OPTIONAL]
    !!             The number of consecutive dipoles to be computed
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iDP1, iNDP
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: rGamma
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(3, 1) :: cG
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iDP1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
             "FDP_Recharge: Bad index range")
      end if
      iStart = iDP1
      if (present(iNDP)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
               "FDP_Recharge: Bad index range")
        end if
        iEnd = iStart+iNDP-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fdp%iCount
    end if

    if (fdp%iCount > 0) then
      ! Sum up the contribution of all dipoles
      allocate(cMapZ(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "FDP_Recharge: Allocation failed")
      rGamma = rZERO
      do i = iStart, iEnd
        dip => fdp%Dipoles(i)
        cG = cIDP_InfluenceG(cMapZ(i))
        rGamma = rGamma + sum(dip%cRho*cG(:, 1))
      end do
      deallocate(cMapZ)
    else
      rGamma = rZERO
    end if

    return
  end function rFDP_Recharge


  function rFDP_Extraction(io, fdp, iDP1, iNDP) result(rQ)
    !! real function rFDP_Extraction
    !!
    !! Computes the recharge rate due to the current set of dipoles.
    !!
    !! Calling Sequence:
    !!    rGamma = rFDP_Extraction(fdp, cZ, iDP1, iNDP)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    integer :: iDP1 [OPTIONAL]
    !!             The index for the first dipole to be used
    !!   (in)    integer :: iNDP [OPTIONAL]
    !!             The number of consecutive dipoles to be computed
    !! Note:
    !!   If iDP1 is not provided, all dipoles will be used
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in), optional :: iDP1, iNDP
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(3, 1) :: cQ
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iDP1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iDP1 >= lbound(fdp%Dipoles, 1) .and. iDP1 <= ubound(fdp%Dipoles, 1)), &
             "FDP_Recharge: Bad index range")
      end if
      iStart = iDP1
      if (present(iNDP)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iDP1+iNDP-1) >= lbound(fdp%Dipoles, 1) .and. (iDP1+iNDP-1) <= ubound(fdp%Dipoles, 1)), &
               "FDP_Recharge: Bad index range")
        end if
        iEnd = iStart+iNDP-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fdp%iCount
    end if

    ! Sum up the contribution of all dipoles
    rQ = rZERO
    do i = iStart, iEnd
      dip => fdp%Dipoles(i)
      cQ = cIDP_InfluenceQ()
      rQ = rQ + sum(dip%cRho*cQ(:, 1))
    end do

    return
  end function rFDP_Extraction


  function rFDP_PotentialJump(io, fdp, cZ, iDP) result(rJump)
    !! real function rFDP_PotentialJump
    !!
    !! Computes the potential jump at the position cZ in dipole iDP of FDP_COLLECTION fdp
    !!
    !! Calling Sequence:
    !!    rJump = rFDP_PotentialJump(fdp, iDP, cZ)
    !!
    !! Arguments:
    !!   (in)    type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!   (in)    integer :: iDP [OPTIONAL]
    !!             The index for the dipole to be used
    !!   (in)    complex :: cZ
    !!             The location to be investigated
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    integer(kind=AE_INT), intent(in) :: iDP
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: rJump
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    complex(kind=AE_REAL), dimension(3, 1) :: cJ
    complex(kind=AE_REAL) :: cMapZ
    complex(kind=AE_REAL) :: cA
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Jump: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Jump: FDP_Alloc has not been called")
      call IO_Assert(io, (iDP > 0 .and. iDP <= fdp%iCount), "FDP_Jump: Illegal index")
    end if

    ! Sum up the contribution of all dipoles
    dip => fdp%Dipoles(iDP)
    cMapZ = (cZ - dip%cZC) / dip%cZL
    cJ = cIDP_InfluenceJ(real(cMapZ))
    rJump = real(sum(dip%cRho*cJ(:, 1)))

    return
  end function rFDP_PotentialJump


  function lFDP_CheckPoint(io, fdp, cZ, rTol, cZFix, rStrength, iElementType, iElementString, iElementVertex, &
             iElementFlag) result(lFound)
    !! function lFDP_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance of the
    !! perimeter of a pond.
    !! If rTol > 0, the tolerance is used in the test, otherwise the global constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
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
    type(FDP_DIPOLE), pointer :: dip
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    real(kind=AE_REAL) :: rCheckTol
    complex(kind=AE_REAL) :: cBZFix
    ! [ LOCAL PARAMETERS ]
    real(kind=AE_REAL), parameter :: rEND_FIX_POSITION = rTWO

    if (fdp%iCount > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a dipole
      allocate(cMapZ(fdp%iCount))
      cMapZ(:) = (cZ-fdp%Dipoles(1:fdp%iCount)%cZC) / fdp%Dipoles(1:fdp%iCount)%cZL

      ! Assume there isn't anything there...
      lFound = .false.
      cZFix = cZERO
      do i = 1, fdp%iCount
        dip => fdp%Dipoles(i)
        ! Check to see(1) we're within the tolerance of the extension of the line segment,
        ! then(2) we're along the(-1, +1) extension
        if (abs(aimag(cMapZ(i))) < rCheckTol) then
          ! We're near a segment. Are we at an end?
          if (abs(real(cMapZ(i))-cONE) < rCheckTol) then
            ! We're at the upper tip
            lFound = .true.
            cBZFix = cmplx(rONE-rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * dip%cZL + dip%cZC
            rStrength = rFDP_PotentialJump(io, fdp, cZFix, i)
            iElementType = dip%iElementType
            iElementString = dip%iElementString
            iElementVertex = dip%iElementVertex
            iElementFlag = dip%iElementFlag
            exit
          else if (abs(real(cMapZ(i))+cONE) < rCheckTol) then
            ! We're at the lower tip
            lFound = .true.
            cBZFix = cmplx(-rONE+rEND_FIX_POSITION*rCheckTol, rZERO, AE_REAL)
            cZFix = cBZFix * dip%cZL + dip%cZC
            rStrength = rFDP_PotentialJump(io, fdp, cZFix, i)
            iElementType = dip%iElementType
            iElementString = dip%iElementString
            iElementVertex = dip%iElementVertex
            iElementFlag = dip%iElementFlag
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
  end function lFDP_CheckPoint


  function lFDP_CheckIntersection(io, fdp, cZ1, cZ2, cZInt, cZBefore, cZAfter, &
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
    type(FDP_COLLECTION), pointer :: fdp
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
    type(FDP_DIPOLE), pointer :: dip
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
    real(kind=AE_REAL), parameter :: rINTERSECTION_TOLERANCE = 1.0e-3_AE_REAL

    ! First, make sure we're not in tolerance
    if (lFDP_CheckPoint(io, fdp, cZ1, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ1 = cZFix
    else
      cMyZ1 = cZ1
    end if
    if (lFDP_CheckPoint(io, fdp, cZ2, rVERTEXTOL*rPOINT_TOLERANCE_FACTOR, cZFix, rStr, iEType, iEStr, iEVtx, iEFlg)) then
      cMyZ2 = cZFix
    else
      cMyZ2 = cZ2
    end if

    if (fdp%iCount > 0) then
      ! Set the check tolerance
      ! Find the smallest distance to a dipole
      allocate(cMapZ1(fdp%iCount), cMapZ2(fdp%iCount), cMapZInt(fdp%iCount), rDist(fdp%iCount))
      cMapZ1(:) = (cMyZ1-fdp%Dipoles(1:fdp%iCount)%cZC) / fdp%Dipoles(1:fdp%iCount)%cZL
      cMapZ2(:) = (cMyZ2-fdp%Dipoles(1:fdp%iCount)%cZC) / fdp%Dipoles(1:fdp%iCount)%cZL
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
        dip => fdp%Dipoles(iLoc)
        ! Reverse-map to get the intersection
        cZInt = cMapZInt(iLoc) * dip%cZL + dip%cZC
        ! Compute ZBefore and ZAfter in the direction of the Z1-Z2 line segment, using the
        ! rINTERSECTION_TOLERANCE*abs(Z2-Z1) as the distance before and after the intersection
        cZBefore = cZInt - rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        cZAfter = cZInt + rINTERSECTION_TOLERANCE * (cZ2-cZ1)
        ! If the before point is farther away from the intersection than the first point,
        ! Use the first point as the before point(and similarly for the after point)
        if (abs(cZ1-cZInt) < abs(cZBefore-cZInt)) cZBefore = cZ1
        if (abs(cZ2-cZInt) < abs(cZAfter-cZInt)) cZAfter = cZ2
        lFound = .true.
        iElementType = dip%iElementType
        iElementString = dip%iElementString
        iElementVertex = dip%iElementVertex
        iElementFlag = dip%iElementFlag
        cENormal = -cI * dip%cZL / abs(dip%cZL)
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
  end function lFDP_CheckIntersection


  subroutine FDP_Report(io, fdp)
    !! subroutine FDP_Report
    !!
    !! Writes a report of all dipole information to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call FDP_Report(io, fdp)
    !!
    !! Arguments:
    !!    (in)   type(FDP_COLLECTION), pointer :: fdp
    !!             The FDP_COLLECTION to use
    !!
    ! [ ARGUMENTS ]
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(FDP_DIPOLE), pointer :: dip

    if (io%lDebug) then
      call IO_Assert(io, (associated(fdp)), "FDP_Update: FDP_Create has not been called")
      call IO_Assert(io, (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called")
    end if

    call HTML_Header('Function module FDP', 1)
    call HTML_Header('Information about line-dipole functions', 2)

    if (fdp%iCount > 0) then
      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'Re(ZC)', 'Im(ZC)', 'Re(ZL)', 'Im(ZL)', &
           'Re(S1)', 'Im(S1)', 'Re(S2)', 'Im(S2)', 'Re(S3)', 'Im(S3)'/))
      do i = 1, fdp%iCount
        dip => fdp%Dipoles(i)
        call HTML_StartRow()
        call HTML_ColumnInteger((/i/))
        call HTML_ColumnComplex((/dip%cZC, dip%cZL, dip%cRho(1), dip%cRho(2), dip%cRho(3)/), 'e13.6')
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No dipole functions defined', 3)
    end if

    return
  end subroutine FDP_Report

end module f_dipole
