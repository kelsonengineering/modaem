module f_pond

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


  !! module f_pond(FPD)
  !!
  !! Written by Victor A.Kelson
  !!
  !! Revision History:
  !!   1.0.0   22 March 1999
  !!           First "source code release" version, adapted from
  !!           previous testing version.
  !!
  !! Module of data structures and functions for 2-D, steady-state
  !! ponds.
  !!
  !! This module encapsulates containers of pond functions
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
  !!   i_pond     --  Influence function module for ponds
  !!
  use u_constants
  use u_io
  use i_pond

  implicit none

  public


  type, public :: FPD_POND
    !! type FPD_POND
    !!
    !! PUBLIC type that holds information for one pond
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the pond
    !!   real :: rRadius
    !!     The radius of the pond
    !!   complex :: rGamma
    !!     Discharge of the pond; positive indicates a pond withdrawal
    !!     and negative and injection rate.
    !!
    complex(kind=AE_REAL) :: cZC
    real(kind=AE_REAL) :: rRadius
    real(kind=AE_REAL) :: rGamma
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
  end type FPD_POND


  type, public :: FPD_COLLECTION
    !! type FPD_COLLECTION
    !!
    !! PUBLIC type that holds the pond entries for a layer
    !!
    !! Members:
    !!   type(FPD_POND), pointer :: Ponds(:)
    !!     FPD_POND structures that hold pond information
    !!   integer :: iCount
    !!     Number of FPD_POND structures currently in use
    !!
    type(FPD_POND), dimension(:), pointer :: Ponds
    integer(kind=AE_INT) :: iCount
  end type


contains


  function FPD_Create(io, iMax) result(fpd)
    !! subroutine FPD_Create
    !!
    !! Creates a new FPD_COLLECTION object
    !!
    !! Calling Sequence:
    !!    fpd => FPD_Create(iMax)
    !!
    !! Arguments:
    !!    (in)      integer :: iMax
    !!                The maximum number of ponds to be stored
    !!
    !! Return Value:
    !!           type(FPD_COLLECTION), pointer :: fpd
    !!             A new FPD_COLLECTION object.  Before use, it will need to
    !!             be allocated using FPD_Alloc.
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    type(FPD_COLLECTION), pointer :: fpd
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(fpd, stat = iStat)
    call IO_Assert(io, (iStat == 0), "FPD_Create: allocation failed")

    allocate(fpd%Ponds(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FPD_Create: allocation failed")
    fpd%iCount = 0
    fpd%Ponds = FPD_POND(cZERO, rZERO, rZERO, -1, -1, -1, -1)

    return
  end function FPD_Create


  subroutine FPD_New(io, fpd, cZc, rGamma, rRadius, iElementType, iElementString, iElementVertex, iElementFlag, iRV)
    !! subroutine FPD_New
    !!
    !! Makes a new pond entry. On call, the geometry and discharge of
    !! the pond are provided. The internal pond structures are then
    !! set up. Returns iRV = the index into the pond table on success or
    !! iRV < 0 on failure.
    !!
    !! Calling Sequence:
    !!    call FPD_New(io, FPD, cZC, rGamma, rRadius, iRV)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    complex :: cZC
    !!             Complex coordinates of the center of the pond
    !!   (in)    real :: rGamma
    !!             The areal exfiltration rate of the pond
    !!   (in)    real :: rRadius
    !!             The radius of the pond
    !!   (out)   integer :: iRV
    !!             The allocated pond index on success
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), intent(in) :: cZc
    real(kind=AE_REAL), intent(in) :: rRadius
    real(kind=AE_REAL), intent(in) :: rGamma
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(out) :: iRV
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_New: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_New: FPD_Alloc has not been called")
      call IO_Assert(io, (fpd%iCount < size(fpd%Ponds)), "FPD_New: Space exhausted")
    end if

    fpd%iCount = fpd%iCount + 1
    pnd => fpd%Ponds(fpd%iCount)
    pnd = FPD_POND(cZc, rRadius, rGamma, iElementType, iElementString, iElementVertex, iElementFlag)
    iRV = fpd%iCount

    return
  end subroutine FPD_New


  subroutine FPD_Update(io, fpd, iPD, rGamma)
    !! subroutine FPD_Update
    !!
    !! Updates the discharge for a pond entry.
    !!
    !! Calling Sequence:
    !!    call FPD_Update(io, fpd, iPD, rGamma)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    integer :: iPD
    !!             The index for the pond in fpd
    !!   (in)    real :: rGamma
    !!             The new areal exfiltration rate of the pond
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    integer(kind=AE_INT), intent(in) :: iPD
    real(kind=AE_REAL), intent(in) :: rGamma
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
      call IO_Assert(io, (iPD < size(fpd%Ponds)), "FPD_Update: Space exhausted")
    end if

    pnd => fpd%Ponds(iPD)
    pnd%rGamma = rGamma

    return
  end subroutine FPD_Update


  subroutine FPD_GetInfluence(io, fpd, iWhich, iPD1, iNPD, cPathZ, cOrientation, cF)
    !! subroutine FPD_GetInfluence
    !!
    !! Retrieves arrays of influence functions for use in matrix generation
    !!
    !! Calling Sequence:
    !!    call FPD_GetInfluence(io, fpd, iWhich, iPD1, iNPD, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    integer :: iWhich
    !!             The influence function to be computed;  iWhich values are
    !!                INFLUENCE_P   - Complex potential
    !!                INFLUENCE_Q   - Complex discharge
    !!                INFLUENCE_F   - Integrated flux
    !!                INFLUENCE_G   - Areal infiltration
    !!                INFLUENCE_Q   - Extraction rate
    !!                INFLUENCE_D   - Difference in potential
    !!                INFLUENCE_Z   - All zeroes
    !!   (in)    integer :: iPD1
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD
    !!             The number of consecutive ponds to be computed
    !!   (in)    complex :: cPathZ(:)
    !!             Complex coordinates of the control point(s) to be used. For
    !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
    !!             cPathZ(1) is used; for iWhich = INFLUENCE_F, the influence
    !!             function is computed along the path cPathZ(:)
    !!   (in)    complex :: cOrientation
    !!             Orientation normal vector for iWhich = INFLUENCE_W
    !!   (out)   complex :: cF(1:iNPD, 3)
    !!             The returned influence functions.  Indexes 1:iNPD relate
    !!             to pond indices iPD1:iPD1+iNPD-1, respectively.
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    integer(kind=AE_INT), intent(in) :: iWhich, iPD1, iNPD
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iPD2, i, j
    real(kind=AE_REAL) :: rF(1, 1), rG(1, 1), rQ(1, 1)
    complex(kind=AE_REAL) :: cW(1, 1), cUnit
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_GetInfluence: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_GetInfluence: FPD_Alloc has not been called")
      call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
           "FPD_GetInfluence: Bad index range")
      call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
           "FPD_GetInfluence: Bad index range")
    end if

    select case (iWhich)
      case (INFLUENCE_P)
        do i = 1, iNPD
          pnd => fpd%Ponds(iPD1+i-1)
          cF(i, :, :) = cIPD_InfluenceP(cPathZ(1), pnd%cZC, pnd%rRadius)
        end do
      case (INFLUENCE_D)
        do i = 1, iNPD
          pnd => fpd%Ponds(iPD1+i-1)
          cF(i, :, :) = cIPD_InfluenceP(cPathZ(1), pnd%cZC, pnd%rRadius) - &
                  cIPD_InfluenceP(cPathZ(2), pnd%cZC, pnd%rRadius)
        end do
      case (INFLUENCE_W)
        do i = 1, iNPD
          pnd => fpd%Ponds(iPD1+i-1)
          cUnit = cOrientation/abs(cOrientation)
          cW = cIPD_InfluenceW(cPathZ(1), pnd%cZC, pnd%rRadius)
          cF(i, :, :) = cUnit * conjg(cW)
        end do
      case (INFLUENCE_F)
        do i = 1, iNPD
          rF = rZERO
          do j = 1, ubound(cPathZ, 1)-1
            pnd => fpd%Ponds(iPD1+i-1)
            rF = rF + cIPD_InfluenceF(cPathZ(j), cPathZ(j+1), pnd%cZC, pnd%rRadius)
          end do
          cF(i, :, :) = -cmplx(rZERO, rF, AE_REAL)
        end do
      case (INFLUENCE_G)
        do i = 1, iNPD
          pnd => fpd%Ponds(iPD1+i-1)
          rG = cIPD_InfluenceG(cPathZ(1), pnd%cZC, pnd%rRadius)
          cF(i, :, :) = rG
        end do
      case (INFLUENCE_Q)
        do i = 1, iNPD
          pnd => fpd%Ponds(iPD1+i-1)
          rQ = cIPD_InfluenceQ(pnd%rRadius)
          cF(i, :, :) = rQ
        end do
      case (INFLUENCE_Z)
        cF = cZERO
    end select

    return
  end subroutine FPD_GetInfluence


  function cFPD_Potential(io, fpd, cZ, iPD1, iNPD) result(cOmega)
    !! complex function cFPD_Potential
    !!
    !! Computes the complex potential due to the current set of ponds.
    !!
    !! Calling Sequence:
    !!    cOmega = cFPD_Potential(iL, cZ, iPD1, iNPD)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iPD1 [OPTIONAL]
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD [OPTIONAL]
    !!             The number of consecutive ponds to be computed
    !! Note:
    !!   If iPD1 is not provided, all ponds will be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iPD1, iNPD
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL) :: cP(1, 1)
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iPD1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
             "FPD_GetInfluence_ILS: Bad index range")
      end if
      iStart = iPD1
      if (present(iNPD)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
               "FPD_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNPD-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fpd%iCount
    end if

    ! Sum up the contributions of all ponds
    cOmega = cZERO
    do i = iStart, iEnd
      pnd => fpd%Ponds(i)
      cP = cIPD_InfluenceP(cZ, pnd%cZC, pnd%rRadius)
      cOmega = cOmega + pnd%rGamma*cP(1, 1)
    end do

    return
  end function cFPD_Potential


  function cFPD_Discharge(io, fpd, cZ, iPD1, iNPD) result(cQ)
    !! complex function cFPD_Discharge
    !!
    !! Computes the complex discharge due to the current set of ponds.
    !!
    !! Calling Sequence:
    !!    cQ = cFPD_Discharge(fpd, cZ, iPD1, iNPD)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iPD1 [OPTIONAL]
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD [OPTIONAL]
    !!             The number of consecutive ponds to be computed
    !! Note:
    !!   If iPD1 is not provided, all ponds will be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iPD1, iNPD

    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL) :: cW(1, 1)
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iPD1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
             "FPD_GetInfluence_ILS: Bad index range")
      end if
      iStart = iPD1
      if (present(iNPD)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
               "FPD_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNPD-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fpd%iCount
    end if

    ! Sum up the contributions of all ponds
    cQ = cZERO
    do i = iStart, iEnd
      pnd => fpd%Ponds(i)
      cW = cIPD_InfluenceW(cZ, pnd%cZC, pnd%rRadius)
      cQ = cQ - conjg(pnd%rGamma * cW(1, 1))
    end do

    return
  end function cFPD_Discharge


  function rFPD_Flow(io, fpd, cPathZ, iPD1, iNPD) result(rFlow)
    !! complex function cFPD_Discharge
    !!
    !! Computes the complex discharge due to the current set of ponds.
    !!
    !! Calling Sequence:
    !!    rFlow = fFPD_Flow(fpd, cPathZ, iPD1, iNPD)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION to be used
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!   (in)    integer :: iPD1 [OPTIONAL]
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD [OPTIONAL]
    !!             The number of consecutive ponds to be computed
    !! Note:
    !!   If iPD1 is not provided, all ponds will be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iPD1, iNPD
    real(kind=AE_REAL) :: rFlow

    ! Locals
    integer(kind=AE_INT) :: i, j, iStart, iEnd
    real(kind=AE_REAL) :: rS(1, 1)
    complex(kind=AE_REAL) :: cZC, cZL, cZMap
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iPD1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
             "FPD_GetInfluence_ILS: Bad index range")
      end if
      iStart = iPD1
      if (present(iNPD)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
               "FPD_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNPD-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fpd%iCount
    end if

    ! Sum up the contributions of all ponds
    rFlow = rZERO
    do i = iStart, iEnd
      pnd => fpd%Ponds(i)
      do j = 1, ubound(cPathZ, 1)-1
        rS = cIPD_InfluenceF(cPathZ(j), cPathZ(j+1), pnd%cZC, pnd%rRadius)
        rFlow = rFlow + pnd%rGamma*rS(1, 1)
      end do
    end do

    return
  end function rFPD_Flow


  function rFPD_Recharge(io, fpd, cZ, iPD1, iNPD) result(rG)
    !! complex function cFPD_Discharge
    !!
    !! Computes the recharge rate due to the current set of ponds.
    !!
    !! Calling Sequence:
    !!    rG = rFPD_Recharge(fpd, cZ, iPD1, iNPD)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    complex :: cZ
    !!             The point at which to determine the potential
    !!   (in)    integer :: iPD1 [OPTIONAL]
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD [OPTIONAL]
    !!             The number of consecutive ponds to be computed
    !! Note:
    !!   If iPD1 is not provided, all ponds will be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iPD1, iNPD
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL) :: cG(1, 1)
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iPD1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
             "FPD_GetInfluence_ILS: Bad index range")
      end if
      iStart = iPD1
      if (present(iNPD)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
               "FPD_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNPD-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fpd%iCount
    end if

    ! Sum up the contributions of all ponds
    rG = rZERO
    do i = iStart, iEnd
      pnd => fpd%Ponds(i)
      cG = cIPD_InfluenceG(cZ, pnd%cZC, pnd%rRadius)
      rG = rG + pnd%rGamma * cG(1, 1)
    end do

    return
  end function rFPD_Recharge


  function rFPD_Extraction(io, fpd, iPD1, iNPD) result(rQ)
    !! complex function cFPD_Discharge
    !!
    !! Computes the recharge rate due to the current set of ponds.
    !!
    !! Calling Sequence:
    !!    rG = rFPD_Extraction(fpd, iPD1, iNPD)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!   (in)    integer :: iPD1 [OPTIONAL]
    !!             The index for the first pond to be used
    !!   (in)    integer :: iNPD [OPTIONAL]
    !!             The number of consecutive ponds to be computed
    !! Note:
    !!   If iPD1 is not provided, all ponds will be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT), intent(in), optional :: iPD1, iNPD

    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL) :: cQ(1, 1)
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    ! Set the range, using the optional arguments
    if (present(iPD1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iPD1 >= lbound(fpd%Ponds, 1) .and. iPD1 <= ubound(fpd%Ponds, 1)), &
             "FPD_GetInfluence_ILS: Bad index range")
      end if
      iStart = iPD1
      if (present(iNPD)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iPD1+iNPD-1) >= lbound(fpd%Ponds, 1) .and. (iPD1+iNPD-1) <= ubound(fpd%Ponds, 1)), &
               "FPD_GetInfluence_ILS: Bad index range")
        end if
        iEnd = iStart+iNPD-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fpd%iCount
    end if

    ! Sum up the contributions of all ponds
    rQ = rZERO
    do i = iStart, iEnd
      pnd => fpd%Ponds(i)
      cQ = cIPD_InfluenceQ(pnd%rRadius)
      rQ = rQ + pnd%rGamma * cQ(1, 1)
    end do

    return
  end function rFPD_Extraction


  function lFPD_CheckPoint(io, fpd, cZ, rTol, cZFix, rStrength, iElementType, iElementString, &
             iElementVertex, iElementFlag) result(lFound)
    !! function lfpd_CheckPoint
    !!
    !! Checks the specified point and returns .true. if the point is within a tolerance of the
    !! perimeter of a pond.
    !! If rTol > 0, the tolerance is used in the test, otherwise the global constant rVERTEXTOL
    !! is used. If the point is within the tolerance, cZFix is set to a point on the well bore.
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
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
    type(FPD_POND), pointer :: pnd
    real(kind=AE_REAL), dimension(:), allocatable :: rDist
    real(kind=AE_REAL) :: rCheckTol

    lFound = .false.
    cZFix = cZERO
    if (fpd%iCount > 0) then
      ! Set the check tolerance
      if (rTol <= rZERO) then
        rCheckTol = rVERTEXTOL
      else
        rCheckTol = rTol
      end if
      ! Find the smallest distance to a well
      allocate(rDist(fpd%iCount))
      rDist = abs(fpd%Ponds(1:fpd%iCount)%rRadius - abs(cZ-fpd%Ponds(1:fpd%iCount)%cZC))
      i = minloc(rDist, DIM = 1)
      if (rDist(i) < rCheckTol**2) then
        lFound = .true.
        pnd => fpd%Ponds(i)
        cZFix = pnd%cZC + pnd%rRadius - rCheckTol * rMOVE_FACTOR
        rStrength = pnd%rGamma
        iElementType = pnd%iElementType
        iElementString = pnd%iElementString
        iElementVertex = pnd%iElementVertex
        iElementFlag = pnd%iElementFlag
      else
        rStrength = rZERO
        iElementType = -1
        iElementString = -1
        iElementVertex = -1
        iElementFlag = -1
        cZFix = cZERO
        lFound = .false.
      end if
      deallocate(rDist)
    else
      rStrength = rZERO
      iElementType = -1
      iElementString = -1
      iElementVertex = -1
      iElementFlag = -1
      cZFix = cZERO
      lFound = .false.
    end if

    return
  end function lFPD_CheckPoint


  subroutine FPD_Report(io, fpd)
    !! subroutine FPD_Report
    !!
    !! Writes a report of all pond information to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call FPD_Report(io, fpd)
    !!
    !! Arguments:
    !!   (in)    type(FPD_COLLECTION), pointer :: fpd
    !!             The FPD_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(FPD_COLLECTION), pointer :: fpd
    type(IO_STATUS), pointer :: io

    ! Locals
    integer(kind=AE_INT) :: i
    type(FPD_POND), pointer :: pnd

    if (io%lDebug) then
      call IO_Assert(io, (associated(fpd)), "FPD_Update: FPD_Create has not been called")
      call IO_Assert(io, (associated(fpd%Ponds)), "FPD_Update: FPD_Alloc has not been called")
    end if

    call HTML_Header('Function module FPD', 1)
    call HTML_Header('Information about pond functions', 2)

    if (fpd%iCount > 0) then
      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'XC    ', 'YC    ', 'Radius', 'Gamma '/))
      do i = 1, fpd%iCount
        pnd => fpd%Ponds(i)
        call HTML_StartRow()
        call HTML_ColumnInteger((/i/))
        call HTML_ColumnComplex((/pnd%cZC/))
        call HTML_ColumnReal((/pnd%rRadius, pnd%rGamma/), 'e13.6')
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No pond functions defined', 3)
    end if

    return
  end subroutine FPD_Report

end module f_pond

