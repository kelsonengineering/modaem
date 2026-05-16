module f_linesink

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

  !! module f_linesink (FLS)
  !!
  !! Written by Victor A. Kelson
  !!
  !! Module of data structures and functions for first-order uniform
  !! linesinks of complex strength.  The real part of the strength
  !! is the 'linesink' strength (jump in the streamfunction, i.e. the
  !! sink density sigma) and the imaginary part is the 'doublet'
  !! strength (jump in the potential).
  !!
  !! Unlike the second-order dipole (f_dipole), the linesink strength
  !! is spatially uniform along each element, so only a single complex
  !! coefficient per element is stored.
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   u_io         --  Universal ModAEM I/O functions and constants
  !!   i_linesink   --  First-order linesink influence function kernels

  use u_constants
  use u_io
  use i_linesink

  implicit none

  public


  type, public :: FLS_LINESINK
    !! type FLS_LINESINK
    !!
    !! Type that holds information for one first-order linesink
    !!
    !! Members:
    !!   complex :: cZC
    !!     Center of the linesink: 0.5*(z2+z1)
    !!   complex :: cZL
    !!     Directed half-length: 0.5*(z2-z1)
    !!   complex :: cSigma
    !!     Complex sink density.  Real part is the linesink (flux) strength;
    !!     imaginary part is the doublet (potential jump) strength.
    !!
    complex(kind=AE_REAL) :: cZC
    complex(kind=AE_REAL) :: cZL
    complex(kind=AE_REAL) :: cSigma
    integer(kind=AE_INT) :: iElementType
    integer(kind=AE_INT) :: iElementString
    integer(kind=AE_INT) :: iElementVertex
    integer(kind=AE_INT) :: iElementFlag
    integer(kind=AE_INT) :: iIndex
  end type FLS_LINESINK


  type, public :: FLS_COLLECTION
    !! type FLS_COLLECTION
    !!
    !! Type that holds the linesink entries for a layer
    !!
    !! Members:
    !!   type(FLS_LINESINK), pointer :: Linesinks(:)
    !!     FLS_LINESINK structures that hold linesink information
    !!   integer :: iCount
    !!     Number of FLS_LINESINK structures currently in use
    !!
    type(FLS_LINESINK), dimension(:), pointer :: Linesinks
    integer(kind=AE_INT) :: iCount
  end type FLS_COLLECTION


contains


  function FLS_Create(io, iMax) result(fls)
    !! function FLS_Create
    !!
    !! Creates a new FLS_COLLECTION object
    !!
    !! Calling Sequence:
    !!    fls => FLS_Create(io, iMax)
    !!
    !! Arguments:
    !!    (in)    integer :: iMax
    !!              The maximum number of linesinks to be stored
    !!
    !! Return Value:
    !!   On success, fls points to a new FLS_COLLECTION object.
    !!   On failure (allocation error), fatal error.
    !!
    ! [ ARGUMENTS ]
    integer(kind=AE_INT), intent(in) :: iMax
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FLS_COLLECTION), pointer :: fls
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(fls, stat = iStat)
    call IO_Assert(io, (iStat == 0), "FLS_Create: allocation failed")

    allocate(fls%Linesinks(iMax), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FLS_Create: allocation failed")

    fls%Linesinks = FLS_LINESINK(cZERO, cmplx(rONE, rONE, AE_REAL), cZERO, -1, -1, -1, -1, -1)
    fls%iCount = 0

    return
  end function FLS_Create


  function FLS_New(io, fls, cZ1, cZ2, cSigma, iElementType, iElementString, iElementVertex, iElementFlag) result(pRV)
    !! function FLS_New
    !!
    !! Makes a new linesink entry.  On call, the geometry and sink density of
    !! the linesink are provided.  The internal linesink structures are then
    !! set up.  Returns pRV pointing to the new FLS_LINESINK on success,
    !! or null if space is exhausted or the element is too short.
    !!
    !! Calling Sequence:
    !!    pRV => FLS_New(io, fls, cZ1, cZ2, cSigma, iElementType, iElementString, iElementVertex, iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    complex :: cZ1, cZ2
    !!             Complex coordinates of the ends of the linesink
    !!   (in)    complex :: cSigma
    !!             The complex sink density (real part: flux strength;
    !!             imaginary part: doublet strength)
    !!   (in)    integer :: iElementType, iElementString, iElementVertex, iElementFlag
    !!             Element bookkeeping tags stored verbatim in the FLS_LINESINK record
    !!
    !! Return value:
    !!   Pointer to the new FLS_LINESINK entry, or null if space is exhausted
    !!   or the element length is less than rTINY
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), intent(in) :: cSigma
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FLS_LINESINK), pointer :: pRV

    nullify(pRV)

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "FLS_New: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "FLS_New: FLS_Create has not been called")
    end if

    if (fls%iCount >= size(fls%Linesinks)) return
    if (abs(cZ2-cZ1) <= rTINY) return

    fls%iCount = fls%iCount + 1
    pRV => fls%Linesinks(fls%iCount)
    pRV = FLS_LINESINK(rHALF*(cZ2+cZ1), rHALF*(cZ2-cZ1), cSigma, &
                       iElementType, iElementString, iElementVertex, iElementFlag, fls%iCount)

    return
  end function FLS_New


  subroutine FLS_GetInfluence(io, fls, iWhich, pLS1, iNLS, cPathZ, cOrientation, cF)
    !! subroutine FLS_GetInfluence
    !!
    !! Retrieves arrays of influence functions for use in matrix generation,
    !! using the first-order linesink kernels from i_linesink.
    !!
    !! Calling Sequence:
    !!    call FLS_GetInfluence(io, fls, iWhich, pLS1, iNLS, cPathZ, cOrientation, cF)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    integer :: iWhich
    !!             The influence function to be computed; iWhich values are
    !!                INFLUENCE_P   - Complex potential
    !!                INFLUENCE_W   - Complex discharge
    !!                INFLUENCE_F   - Integrated flux
    !!                INFLUENCE_G   - Areal infiltration (always zero for linesinks)
    !!                INFLUENCE_Q   - Extraction rate
    !!                INFLUENCE_D   - Difference in potential
    !!                INFLUENCE_Z   - All zeroes
    !!   (in)    type(FLS_LINESINK), pointer :: pLS1
    !!             Pointer to the first linesink to be used
    !!   (in)    integer :: iNLS
    !!             The number of consecutive linesinks to be computed
    !!   (in)    complex :: cPathZ(:)
    !!             Complex coordinates of the control point(s) to be used.  For
    !!             iWhich = (INFLUENCE_P, INFLUENCE_W, INFLUENCE_G) only
    !!             cPathZ(1) is used; for iWhich = INFLUENCE_F, the influence
    !!             function is computed along the path cPathZ(:)
    !!   (in)    complex :: cOrientation
    !!             Orientation normal vector for iWhich = INFLUENCE_W
    !!   (out)   complex :: cF(1:iNLS, 1, 1)
    !!             The returned influence functions.  Index 1:iNLS relates
    !!             to the iNLS consecutive linesinks starting at pLS1%iIndex.
    !!             Only cF(i, 1, 1) is populated; the second and third dimensions
    !!             are present for calling-convention compatibility with other
    !!             GetInfluence routines.
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    integer(kind=AE_INT), intent(in) :: iWhich, iNLS
    type(FLS_LINESINK), pointer :: pLS1
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    complex(kind=AE_REAL), dimension(:, :, :), intent(out) :: cF
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iLS1, iLS2, i, j, iStat
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ2
    complex(kind=AE_REAL) :: cUnit
    complex(kind=AE_REAL), dimension(1, 1) :: cP, cR, cS, cG, cQ
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "FLS_GetInfluence: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "FLS_GetInfluence: FLS_Create has not been called")
      call IO_Assert(io, (associated(pLS1)), "FLS_GetInfluence: null pointer pLS1")
      call IO_Assert(io, (pLS1%iIndex >= lbound(fls%Linesinks, 1) .and. pLS1%iIndex <= ubound(fls%Linesinks, 1)), &
           "FLS_GetInfluence: Bad index range")
      call IO_Assert(io, ((pLS1%iIndex+iNLS-1) >= lbound(fls%Linesinks, 1) .and. &
                          (pLS1%iIndex+iNLS-1) <= ubound(fls%Linesinks, 1)), &
           "FLS_GetInfluence: Bad index range")
      call IO_Assert(io, (size(cF, 1) >= iNLS .and. size(cF, 2) >= 1 .and. size(cF, 3) >= 1), &
           "FLS_GetInfluence: Invalid result array")
    end if

    iLS1 = pLS1%iIndex
    iLS2 = iLS1+iNLS-1
    allocate(cMapZ1(iLS1:iLS2), cMapZ2(iLS1:iLS2), stat = iStat)
    call IO_Assert(io, (iStat == 0), "FLS_GetInfluence: Allocation failed")

    ! It is assumed that the caller has eliminated the potential for a singularity
    ! by selecting appropriate control points
    select case (iWhich)
      case (INFLUENCE_P)
        cMapZ1(:) = (cPathZ(1)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
        do i = 1, iNLS
          ls => fls%Linesinks(iLS1+i-1)
          cP = cILS_InfluenceP(cMapZ1(iLS1+i-1), ls%cZL)
          cF(i, 1, 1) = cP(1, 1)
        end do
      case (INFLUENCE_D)
        cMapZ1(:) = (cPathZ(1)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
        cMapZ2(:) = (cPathZ(2)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
        do i = 1, iNLS
          ls => fls%Linesinks(iLS1+i-1)
          cP = cILS_InfluenceP(cMapZ1(iLS1+i-1), ls%cZL) - cILS_InfluenceP(cMapZ2(iLS1+i-1), ls%cZL)
          cF(i, 1, 1) = cP(1, 1)
        end do
      case (INFLUENCE_W)
        cMapZ1(:) = (cPathZ(1)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
        do i = 1, iNLS
          ls => fls%Linesinks(iLS1+i-1)
          cUnit = cOrientation/abs(cOrientation)
          cR = cILS_InfluenceW(cMapZ1(iLS1+i-1), ls%cZL)
          cF(i, 1, 1) = cUnit * (cR(1, 1) / ls%cZL)
        end do
      case (INFLUENCE_F)
        cF = cZERO
        do j = 1, size(cPathZ)-1
          cMapZ1(:) = (cPathZ(j)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
          cMapZ2(:) = (cPathZ(j+1)-fls%Linesinks(iLS1:iLS2)%cZC) / fls%Linesinks(iLS1:iLS2)%cZL
          do i = 1, iNLS
            ls => fls%Linesinks(iLS1+i-1)
            cS = cILS_InfluenceF(cMapZ1(iLS1+i-1), cMapZ2(iLS1+i-1), ls%cZL)
            cF(i, 1, 1) = cF(i, 1, 1) + cS(1, 1)
          end do
        end do
      case (INFLUENCE_G)
        ! Linesinks satisfy Laplace's equation; the recharge influence is zero
        do i = 1, iNLS
          cG = cILS_InfluenceG()
          cF(i, 1, 1) = cG(1, 1)
        end do
      case (INFLUENCE_Q)
        do i = 1, iNLS
          ls => fls%Linesinks(iLS1+i-1)
          cQ = cILS_InfluenceQ(ls%cZL)
          cF(i, 1, 1) = cQ(1, 1)
        end do
      case (INFLUENCE_Z)
        cF = cZERO
    end select

    deallocate(cMapZ1, cMapZ2)

    return
  end subroutine FLS_GetInfluence


  function cFLS_Potential(io, fls, cZ, iLS1, iNLS) result(cOmega)
    !! complex function cFLS_Potential
    !!
    !! Computes the complex potential due to the specified linesinks.
    !!
    !! Calling Sequence:
    !!    cOmega = cFLS_Potential(io, fls, cZ, iLS1, iNLS)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    complex :: cZ
    !!             The point at which to evaluate the potential
    !!   (in)    integer :: iLS1 [OPTIONAL]
    !!             Index of the first linesink to include
    !!   (in)    integer :: iNLS [OPTIONAL]
    !!             Number of consecutive linesinks to include
    !!
    !! Note:
    !!   If iLS1 is not provided, all linesinks in the collection are used.
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iLS1, iNLS
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(1, 1) :: cP
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "cFLS_Potential: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "cFLS_Potential: FLS_Create has not been called")
    end if

    if (present(iLS1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iLS1 >= lbound(fls%Linesinks, 1) .and. iLS1 <= ubound(fls%Linesinks, 1)), &
             "cFLS_Potential: Bad index range")
      end if
      iStart = iLS1
      if (present(iNLS)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iLS1+iNLS-1) >= lbound(fls%Linesinks, 1) .and. &
                              (iLS1+iNLS-1) <= ubound(fls%Linesinks, 1)), &
               "cFLS_Potential: Bad index range")
        end if
        iEnd = iStart+iNLS-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fls%iCount
    end if

    if (fls%iCount > 0) then
      allocate(cMapZ(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "cFLS_Potential: Allocation failed")
      cOmega = cZERO
      cMapZ(iStart:iEnd) = (cZ - fls%Linesinks(iStart:iEnd)%cZC) / fls%Linesinks(iStart:iEnd)%cZL
      do i = iStart, iEnd
        ls => fls%Linesinks(i)
        cP = cILS_InfluenceP(cMapZ(i), ls%cZL)
        cOmega = cOmega + ls%cSigma * cP(1, 1)
      end do
      deallocate(cMapZ)
    else
      cOmega = cZERO
    end if

    return
  end function cFLS_Potential


  function cFLS_Discharge(io, fls, cZ, iLS1, iNLS) result(cQ)
    !! complex function cFLS_Discharge
    !!
    !! Computes the complex discharge due to the specified linesinks.
    !!
    !! Calling Sequence:
    !!    cQ = cFLS_Discharge(io, fls, cZ, iLS1, iNLS)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    complex :: cZ
    !!             The point at which to evaluate the discharge
    !!   (in)    integer :: iLS1 [OPTIONAL]
    !!             Index of the first linesink to include
    !!   (in)    integer :: iNLS [OPTIONAL]
    !!             Number of consecutive linesinks to include
    !!
    !! Note:
    !!   If iLS1 is not provided, all linesinks in the collection are used.
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in), optional :: iLS1, iNLS
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(1, 1) :: cW
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "cFLS_Discharge: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "cFLS_Discharge: FLS_Create has not been called")
    end if

    if (present(iLS1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iLS1 >= lbound(fls%Linesinks, 1) .and. iLS1 <= ubound(fls%Linesinks, 1)), &
             "cFLS_Discharge: Bad index range")
      end if
      iStart = iLS1
      if (present(iNLS)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iLS1+iNLS-1) >= lbound(fls%Linesinks, 1) .and. &
                              (iLS1+iNLS-1) <= ubound(fls%Linesinks, 1)), &
               "cFLS_Discharge: Bad index range")
        end if
        iEnd = iStart+iNLS-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fls%iCount
    end if

    if (fls%iCount > 0) then
      allocate(cMapZ(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "cFLS_Discharge: Allocation failed")
      cQ = cZERO
      cMapZ(iStart:iEnd) = (cZ - fls%Linesinks(iStart:iEnd)%cZC) / fls%Linesinks(iStart:iEnd)%cZL
      do i = iStart, iEnd
        ls => fls%Linesinks(i)
        cW = cILS_InfluenceW(cMapZ(i), ls%cZL)
        cQ = cQ + conjg(ls%cSigma * cW(1, 1) / ls%cZL)
      end do
      deallocate(cMapZ)
    else
      cQ = cZERO
    end if

    return
  end function cFLS_Discharge


  function rFLS_Flow(io, fls, cPathZ, iLS1, iNLS) result(rFlow)
    !! real function rFLS_Flow
    !!
    !! Computes the integrated flow across cPathZ due to the specified linesinks.
    !!
    !! Calling Sequence:
    !!    rFlow = rFLS_Flow(io, fls, cPathZ, iLS1, iNLS)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    complex :: cPathZ(:)
    !!             The path across which the flow is desired
    !!   (in)    integer :: iLS1 [OPTIONAL]
    !!             Index of the first linesink to include
    !!   (in)    integer :: iNLS [OPTIONAL]
    !!             Number of consecutive linesinks to include
    !!
    !! Note:
    !!   If iLS1 is not provided, all linesinks in the collection are used.
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    integer(kind=AE_INT), intent(in), optional :: iLS1, iNLS
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, j, iStart, iEnd, iStat
    complex(kind=AE_REAL), dimension(1, 1) :: cS
    complex(kind=AE_REAL), dimension(:), allocatable :: cMapZ1, cMapZ2
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "rFLS_Flow: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "rFLS_Flow: FLS_Create has not been called")
    end if

    if (present(iLS1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iLS1 >= lbound(fls%Linesinks, 1) .and. iLS1 <= ubound(fls%Linesinks, 1)), &
             "rFLS_Flow: Bad index range")
      end if
      iStart = iLS1
      if (present(iNLS)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iLS1+iNLS-1) >= lbound(fls%Linesinks, 1) .and. &
                              (iLS1+iNLS-1) <= ubound(fls%Linesinks, 1)), &
               "rFLS_Flow: Bad index range")
        end if
        iEnd = iStart+iNLS-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fls%iCount
    end if

    if (fls%iCount > 0) then
      allocate(cMapZ1(iStart:iEnd), cMapZ2(iStart:iEnd), stat = iStat)
      call IO_Assert(io, (iStat == 0), "rFLS_Flow: Allocation failed")
      rFlow = rZERO
      do j = 1, ubound(cPathZ, 1)-1
        cMapZ1(iStart:iEnd) = (cPathZ(j)-fls%Linesinks(iStart:iEnd)%cZC) / fls%Linesinks(iStart:iEnd)%cZL
        cMapZ2(iStart:iEnd) = (cPathZ(j+1)-fls%Linesinks(iStart:iEnd)%cZC) / fls%Linesinks(iStart:iEnd)%cZL
        do i = iStart, iEnd
          ls => fls%Linesinks(i)
          cS = cILS_InfluenceF(cMapZ1(i), cMapZ2(i), ls%cZL)
          rFlow = rFlow + real(ls%cSigma * cS(1, 1))
        end do
      end do
      deallocate(cMapZ1, cMapZ2)
    else
      rFlow = rZERO
    end if

    return
  end function rFLS_Flow


  function rFLS_Extraction(io, fls, iLS1, iNLS) result(rQ)
    !! real function rFLS_Extraction
    !!
    !! Computes the total volumetric extraction rate due to the specified linesinks.
    !!
    !! Calling Sequence:
    !!    rQ = rFLS_Extraction(io, fls, iLS1, iNLS)
    !!
    !! Arguments:
    !!   (in)    type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to use
    !!   (in)    integer :: iLS1 [OPTIONAL]
    !!             Index of the first linesink to include
    !!   (in)    integer :: iNLS [OPTIONAL]
    !!             Number of consecutive linesinks to include
    !!
    !! Note:
    !!   If iLS1 is not provided, all linesinks in the collection are used.
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    integer(kind=AE_INT), intent(in), optional :: iLS1, iNLS
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rQ
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, iStart, iEnd
    complex(kind=AE_REAL), dimension(1, 1) :: cQ
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "rFLS_Extraction: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "rFLS_Extraction: FLS_Create has not been called")
    end if

    if (present(iLS1)) then
      if (io%lDebug) then
        call IO_Assert(io, (iLS1 >= lbound(fls%Linesinks, 1) .and. iLS1 <= ubound(fls%Linesinks, 1)), &
             "rFLS_Extraction: Bad index range")
      end if
      iStart = iLS1
      if (present(iNLS)) then
        if (io%lDebug) then
          call IO_Assert(io, ((iLS1+iNLS-1) >= lbound(fls%Linesinks, 1) .and. &
                              (iLS1+iNLS-1) <= ubound(fls%Linesinks, 1)), &
               "rFLS_Extraction: Bad index range")
        end if
        iEnd = iStart+iNLS-1
      else
        iEnd = iStart
      end if
    else
      iStart = 1
      iEnd = fls%iCount
    end if

    rQ = rZERO
    do i = iStart, iEnd
      ls => fls%Linesinks(i)
      cQ = cILS_InfluenceQ(ls%cZL)
      rQ = rQ + real(ls%cSigma * cQ(1, 1))
    end do

    return
  end function rFLS_Extraction


  subroutine FLS_Report(io, fls)
    !! subroutine FLS_Report
    !!
    !! Writes a report of all linesink information to LU_OUTPUT
    !!
    !! Calling Sequence:
    !!    call FLS_Report(io, fls)
    !!
    !! Arguments:
    !!    (in)   type(FLS_COLLECTION), pointer :: fls
    !!             The FLS_COLLECTION to report on
    !!
    ! [ ARGUMENTS ]
    type(FLS_COLLECTION), pointer :: fls
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i
    type(FLS_LINESINK), pointer :: ls

    if (io%lDebug) then
      call IO_Assert(io, (associated(fls)), "FLS_Report: FLS_Create has not been called")
      call IO_Assert(io, (associated(fls%Linesinks)), "FLS_Report: FLS_Create has not been called")
    end if

    call HTML_Header('Function module FLS', 1)
    call HTML_Header('Information about first-order linesink functions', 2)

    if (fls%iCount > 0) then
      call HTML_StartTable()
      call HTML_TableHeader((/'      ', 'Re(ZC)', 'Im(ZC)', 'Re(ZL)', 'Im(ZL)', 'Re(S) ', 'Im(S) '/))
      do i = 1, fls%iCount
        ls => fls%Linesinks(i)
        call HTML_StartRow()
        call HTML_ColumnInteger((/i/))
        call HTML_ColumnComplex((/ls%cZC, ls%cZL, ls%cSigma/), 'e13.6')
        call HTML_EndRow()
      end do
      call HTML_EndTable()
    else
      call HTML_Header('No linesink functions defined', 3)
    end if

    return
  end subroutine FLS_Report


end module f_linesink
