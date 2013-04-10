module m_cw0

  ! ModAEM 1.8
  ! Copyright(c) 2003-2008 WHPA Inc. and Vic Kelson
  !
  ! This code is the property of WHPA, Inc.  Do not distribute.
  !

  !! module m_cw0
  !!
  !! Element module for 2-D discharge specified collector wells with resistance
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! line-sink elements. Elements are defined as strings of points, with the
  !! discharge rate(per length) specified for each segment.
  !!
  !! Note: ModAEM uses the "traditional" line-sink function (io, Strack, 1989)
  !! for matrix generation, but uses strings of dipoles terminated by a well
  !! for computational performance, once a solution is achieved.

  use u_constants
  use u_io
  use f_well
  use f_dipole
  use u_matrix
  use m_aqu

  implicit none

  public

  type :: CW0_VERTEX
    !! type CW0_VERTEX    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   complex :: cZCPZ
    !!     The complex coordinate of the segment's control point
    !!   real :: rHeadDifference
    !!     The difference in head between the vertex and next vertex in a string
    !!   real :: rDPStrength
    !!     The dipole strength at the vertex
    !!   real :: rLength
    !!     The segment length
    !!   real :: rSpecValue
    !!     The 'specified value' (potential difference between this vertex and the
    !!     next one) for this iteration
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !!   real :: rCheck
    !!     Heads at the segment's control-point
    !!
    complex(kind=AE_REAL) :: cZ
    complex(kind=AE_REAL) :: cCPZ
    real(kind=AE_REAL) :: rHeadDifference
    real(kind=AE_REAL) :: rDPStrength
    real(kind=AE_REAL) :: rStrength
    real(kind=AE_REAL) :: rLength
    real(kind=AE_REAL) :: rSpecValue
    integer(kind=AE_INT) :: iFDPIndex
    real(kind=AE_REAL) :: rCheck
    real(kind=AE_REAL) :: rTransmissivity
    real(kind=AE_REAL) :: rResistanceTerm
  end type CW0_VERTEX

  type :: CW0_RADIAL
    !! type CW0_RADIAL
    !!
    !! Type that holds information for one radial arm of a collector well
    !!
    !! Members:
    !!   complex :: cZEnd
    !!     Complex coordinate of the end of the radial
    !!   type(CW0_VERTEX), dimension(:), pointer :: Vertices
    !!     A vector of CW0_VERTEX objects
    !!   integer :: iFWLIndex
    !!     Points to the FWL function associated with the(inner) end of the radial
    !!
    real(kind=AE_REAL) :: rResistance
    real(kind=AE_REAL) :: rWidth
    real(kind=AE_REAL) :: rBlankLength
    complex(kind=AE_REAL) :: cZEnd
    type(CW0_VERTEX), dimension(:), pointer :: Vertices
    integer(kind=AE_INT) :: iFWLIndex
  end type CW0_RADIAL

  type :: CW0_WELL
    !! type CW0_WELL
    !!
    !! Type that holds information for one collector well
    !!
    !! Members:
    !!   complex :: cZC
    !!     Coordinate of the center of the well
    !!   real :: rC
    !!     Radius of the well caisson
    !!   real :: rQ
    !!     Pumping rate
    !!   real :: rResistance
    !!     Resistance of the radial arms, defined as
    !!                    (resistance of the arms)/(perimeter of the arms)
    !!   integer :: iNRad
    !!     Number of radial arms
    !!   integer :: iResolution
    !!     Number of linesinks to be manufactured per radial arm
    !!   type(CW0_RADIAL), dimension(:), pointer :: Radials
    !!     A vector of CW0_RADIAL objects for the radial arms
    !!   real :: rCheck
    !!     The total pumping rate of the well at this iteration
    !!
    integer(kind=AE_INT) :: iNRad
    complex(kind=AE_REAL) :: cZC
    real(kind=AE_REAL) :: rRc
    real(kind=AE_REAL) :: rQ
    integer(kind=AE_INT) :: iResolution
    integer(kind=AE_INT) :: iID
    type(CW0_RADIAL), dimension(:), pointer :: Radials
    real(kind=AE_REAL) :: rCheck
    logical :: lHeadSpec
    ! For Drawdown simulations
    logical :: lDdn
    logical :: lDDnHeadSpec
    real(kind=AE_REAL) :: rDdnValue
  end type CW0_WELL

  type :: CW0_COLLECTION
    !! type CW0_COLLECTION
    !!
    !! Type that holds information for all CW0 elements in a layer
    !!
    !! Members:
    !!   type(CW0_WELL), dimension(:), pointer :: Wells
    !!     A vector of CW0_WELL objects
    !!   integer :: iCount
    !!     The number of wells actually in use
    !!   integer :: iFluxSpecified
    !!     Nonzero if the element is to be flux-specified(first iteration)
    !!
    type(CW0_WELL), dimension(:), pointer :: Wells
    integer(kind=AE_INT) :: iCount
    integer(kind=AE_INT) :: iRegenerate
    ! Iterator Information
    integer(kind=AE_INT) :: iIterWell
    integer(kind=AE_INT) :: iIterRad
    integer(kind=AE_INT) :: iIterVtx
  end type CW0_COLLECTION


contains


  function CW0_Create(io) result(cw0)
    !! function CW0_Create
    !!
    !! Creates a new CW0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    cw0 => CW0_Create()
    !!
    !! Arguments:
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    type(CW0_COLLECTION), pointer :: cw0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(cw0, stat = iStat)
    call IO_Assert(io, (iStat == 0), "CW0_Create: allocation failed")
    nullify(cw0%Wells)
    cw0%iCount = 0
    cw0%iRegenerate = 1

    return
  end function CW0_Create


  subroutine CW0_Alloc(io, cw0)
    !! Subroutine CW0_Alloc
    !!
    !! Allocates Strings for the CW0_COLLECTION object
    !!
    !! Calling Sequence:
    !!    call CW0_Alloc(io, cw0, iCount)
    !!
    !! Arguments:
    !!    (in)    type(CW0_COLLECTION), pointer :: ls0
    !!              The CW0_COLLECTION object to be used
    !!    (in)    integer :: iCount
    !!              The number of strings to make space for
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iCount
    integer(kind=AE_INT) :: iStat

    iCount = iIO_GetInteger(io, "iCount", minimum=0)
    allocate(cw0%Wells(iCount), stat = iStat)
    call IO_Assert(io, (iStat == 0), "CW0_Alloc: allocation failed")
    cw0%Wells(:)%iID = -1

    return
  end subroutine CW0_Alloc


  subroutine CW0_Destroy(io, cw0)
    !! subroutine CW0_Destroy
    !!
    !! Frees memory allocated for cw0 Linesinks and strings of vertices
    !! and the cw0 Collection object
    !!
    !! Calling Sequence:
    !!     call CW0_Destroy(cw0)
    !!
    !! Arguments:
    !!  type(cw0_COLLECTION), pointer :: cw0
    !!              Pointer to the cw0_COLLECTION object to be used
    !!
    !! Return Value:
    !!
    ! [ ARGUMENTS ]
    type(cw0_COLLECTION), pointer :: cw0
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    integer(kind=AE_INT) :: iWel, iRad
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad

    if (io%lDebug) then
      call IO_Assert(io, (associated(cw0)), &
           "CW0_Destroy: CW0_Create has not been called")
    end if

    if (cw0%iCount > 0) then
      do iWel = 1, cw0%iCount
        wel => cw0%Wells(iWel)
        do iRad = 1, wel%iNRad
          rad => wel%Radials(iRad)
          deallocate(rad%Vertices, stat = iStat)
          call IO_Assert(io, (iStat == 0), "CW0_Destroy: deallocation of vertices failed")
        end do
        deallocate(wel%Radials, stat = iStat)
        call IO_Assert(io, (iStat == 0), "CW0_Destroy: deallocation of radials failed")
      end do
      deallocate(cw0%Wells, stat = iStat)
      call IO_Assert(io, (iStat == 0), "CW0_Destroy: deallocation of wells failed")
      ! Now deallocate the collection
      deallocate(cw0, stat = iStat)
      call IO_Assert(io, (iStat == 0), "CW0_Destroy: deallocation failed")
    end if

    return
  end subroutine CW0_Destroy


  subroutine CW0_New(cw0)
    !! function CW0_New
    !!
    !! Adds a new CW0_WELL object to the CW0_COLLECTION 'cw0'
    !!
    !! Calling Sequence:
    !!    call CW0_New(io, cw0, Vertices, iNPt)
    !!
    !! Arguments:
    !!    (in)    type(CW0_COLLECTION), pointer :: cw0
    !!              The CW0_COLLECTION object to be used
    !!    (in)    type(CW0_VERTEX) :: Vertices(:)
    !!              Vector that defines the points along the barrier
    !!    (in)    integer :: iNPt
    !!              The number of vertices in the string
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]

    call IO_Assert(io, (.false.), "CW0_New is not yet implemented")

    return
  end subroutine CW0_New


  function iCW0_GetID(io, cw0, iIndex) result(iID)
    !! Returns the ID number for the well at index 'iIndex'
    type(CW0_COLLECTION), pointer :: cw0
    integer(kind=AE_INT), intent(in) :: iIndex
    type(IO_STATUS), pointer :: io
    integer(kind=AE_INT) :: iID

    call IO_Assert(io, (iIndex > 0 .and. iIndex <= cw0%iCount), "Internal error -- no such index")
    iID = cw0%Wells(iIndex)%iID

    return
  end function iCW0_GetID


  subroutine CW0_PreSolve(io, cw0)
    !! subroutine CW0_PreSolve
    !!
    !! Steps to be executed prior to beginning the solution process
    !! This routine adjusts elements as necessary, and allocates internal buffers
    !!
    !! Calling Sequence:
    !!    call CW0_PreSolve(cw0)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer :: cw0
    !!             CW0_COLLECTION to be used
    !!   (in)    type(IO_status), pointer :: io
    !!              pointer to IO_Status structure
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(IO_STATUS), pointer :: io

    return
  end subroutine CW0_PreSolve


  function iCW0_GetInfo(io, cw0, iOption, iIteration) result(iValue)
    !! function CW0_GetInfo
    !!
    !! Returns the following sizing requirements for the WL0module
    !!
    !! Calling Sequence:
    !!    iValue = iCW0_GetInfo(io, cw0, iOption)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer :: cw0
    !!             CW0_COLLECTION to be used
    !!   (in)   integer :: iOption
    !!             The(see u_constants.f90) to be retrieved
    !!   (in)   integer :: iIteration
    !!             The current iteration
    !!
    !! Return Value:
    !!   integer :: iOption
    !!     The requested information for the object. Note: Unrecognized options
    !!     should always return zero; (via 'case default' in 'select' structure)
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    integer(kind=AE_INT), intent(in) :: iOption
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iValue
    integer(kind=AE_INT) :: iWel, iRad
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad

    if (io%lDebug) then
      call IO_Assert(io, (associated(cw0)), &
           "CW0_GetInfo: CW0_Create has not been called")
    end if

    iValue = 0
    select case (iOption)
      case (SIZE_FWL)
        iValue = 0
        do iWel = 1, cw0%iCount
          wel => cw0%Wells(iWel)
          iValue = iValue + wel%iNRad
        end do
      case (SIZE_FDP)
        iValue = 0
        do iWel = 1, cw0%iCount
          wel => cw0%Wells(iWel)
          do iRad = 1, wel%iNRad
            rad => wel%Radials(iRad)
            iValue = iValue + wel%iResolution
          end do
        end do
      case (SIZE_EQUATIONS)
        iValue = 0
        do iWel = 1, cw0%iCount
          wel => cw0%Wells(iWel)
          do iRad = 1, wel%iNRad
            rad => wel%Radials(iRad)
            iValue = iValue + wel%iResolution
          end do
        end do
      case (SIZE_UNKNOWNS)
        iValue = 0
        do iWel = 1, cw0%iCount
          wel => cw0%Wells(iWel)
          do iRad = 1, wel%iNRad
            rad => wel%Radials(iRad)
            iValue = iValue + wel%iResolution
          end do
        end do
      case (INFO_REGENERATE)
        iValue = 0
        if (iIteration == 1) then
          iValue = 1
        end if
      case default
        iValue = 0
    end select

    return
  end function iCW0_GetInfo


  subroutine CW0_SetupFunctions(io, cw0, fwl, fdp)
    !! subroutine CW0_Setup
    !!
    !! This routine sets up the functions in f_well and f_dipole for the line-sinks
    !! Since this module creates given-strength elements, the strengths of
    !! all functions are computed at set-up time.
    !!
    !! Note: This routine assumes that sufficient space has been allocated
    !! in f_well and in f_dipole by SOL_Alloc.
    !!
    !! Calling Sequence:
    !!    call CW0_Setup(cw0)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer
    !!             CW0_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
    !!             FDP_COLLECTION object to be used
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iSeg, iWel, iRad, iVtx, i, iDP, iWL, iStat
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead, rLen, rSumRes, &
                         rTotalLength, rAvgStrength
    complex(kind=AE_REAL) :: cZ1, cZ, cDZi, cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad
    type(CW0_VERTEX), pointer :: this, next, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(cw0)), &
           "CW0_Setup: CW0_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "CW0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "CW0_Setup: Illegal FDP_COLLECTION object")
    end if

    do iWel = 1, cw0%iCount
      wel => cw0%Wells(iWel)
      ! Create the vertices for the sub-segments on the radials
      rTotalLength = rZERO
      do iRad = 1, wel%iNRad
        ! Build dipoles for all segments of each radial arm
        rad => wel%Radials(iRad)
        ! The length of the radial arm is the distance from the tip of the arm
        ! to the center of the radial collector, less the radius of the caisson
        !rLen = abs(rad%cZEnd - wel%cZC) - wel%rRC
        rLen = abs(rad%cZEnd - wel%cZC) - rad%rBlankLength
        ! Unit vector along the arm(from tip towards center)
        cDZi = (wel%cZC-rad%cZEnd)/abs(wel%cZC-rad%cZEnd)
        allocate(rad%Vertices(wel%iResolution+1), stat = iStat)
        call IO_Assert(io, (iStat == 0), "CW0_SetupFunctions: Allocation failed")

        ! OK.  Now we'll put the vertices along the arm, using constant spacing,
        do iVtx = 1, wel%iResolution+1
          this => rad%Vertices(iVtx)
          this%cZ = rad%cZEnd + cDZI * rLen * float(iVtx-1) / float(wel%iResolution)
          this%rHeadDifference = rZERO
          this%rDPStrength = rZERO
          this%rStrength = rZERO
          this%rLength = rZERO
          this%rSpecValue = rZERO
          this%rCheck = rZERO
        end do

        ! Woo hoo! Now, make the line dipoles and a well for each arm
        do iVtx = 1, wel%iResolution
          this => rad%Vertices(iVtx)
          next => rad%Vertices(iVtx+1)
          this%rLength = abs(next%cZ - this%cZ)
          rTotalLength = rTotalLength + this%rLength
          this%rDPStrength = rZERO
          call FDP_New(io, fdp, this%cZ, next%cZ, (/cZERO, cZERO, cZERO/), ELEM_CW0, &
               iWel, iRad, iWel, this%iFDPIndex)
          ! For this version, no overspecification -- control points at centers of segments
          this%cCPZ = rHALF * (this%cZ + next%cZ)
        end do
        ! Put a well at the end of each arm
        last_vtx => rad%Vertices(wel%iResolution+1)
        call FWL_New(io, fwl, last_vtx%cZ, rZERO, rZERO, ELEM_CW0, iWel, iRad, iVtx, iWL)
        rad%iFWLIndex = iWL
      end do

      ! Now, set all the arms' strengths to the average pumping rate per unit length
      !    if (.not. wel%lHeadSpec) then
      !      rAvgStrength = wel%rQ / rTotalLength
      !      do iRad = 1, wel%iNRad
      !        do iVtx = 1, wel%iResolution
      !          call CW0_StoreResult(io, cw0, rAvgStrength, ELEM_CW0, iWel, iRad, iVtx, .true.)
      !        end do
      !      end do
      !    end if
    end do

    ! Make sure that the functions are updated
    call CW0_Update(io, cw0, fwl, fdp)

    return
  end subroutine CW0_SetupFunctions


  subroutine CW0_SetupMatrix(io, cw0, aqu, mat)
    !! subroutine CW0_Setup
    !!
    !!
    !! Calling Sequence:
    !!    call CW0_Setup(cw0)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer
    !!             CW0_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer
    !!             AQU_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(AQU_COLLECTION), pointer :: aqu
    type(MAT_MATRIX), pointer :: mat
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iWel, iRad, iVtx, i, iDP, iWL, iEQ
    real(kind=AE_REAL) :: rStrength, rDisch, rHead1, rHead2, rHead, rHydCond
    complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
    complex(kind=AE_REAL), dimension(3) :: cCPResult
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad
    type(CW0_VERTEX), pointer :: this, next, last_vtx

    if (io%lDebug) then
      call IO_Assert(io, (associated(cw0)), &
           "CW0_Setup: CW0_Create has not been called")
      call IO_Assert(io, (associated(mat)), &
           "CW0_Setup: Illegal MAT_MATRIX object")
    end if

    ! Build matrix generator entries for all segments
    do iWel = 1, cw0%iCount
      ! Set up the unknown variables -- one per segment
      wel => cw0%Wells(iWel)
      do iRad = 1, wel%iNRad
        rad => wel%Radials(iRad)
        do iVtx = 1, wel%iResolution
          call MAT_CreateVariable(io, mat, ELEM_CW0, iWel, iRad, iVtx)
        end do
      end do

      ! Set up equations -- One equation per segment
      ! last_vtx is the last segment of the last radial
      ! I know this is grotesque looking, but...
      last_vtx => wel%Radials(wel%iNRad)%Vertices(wel%iResolution)
      do iRad = 1, wel%iNRad
        rad => wel%Radials(iRad)
        do iVtx = 1, wel%iResolution
          this => rad%Vertices(iVtx)
          next => rad%Vertices(iVtx+1)
          if (wel%lHeadSpec) then
            iEQ = MAT_CreateEquation(io, mat, (/this%cCPZ/), EQN_HEAD, &
                                     ELEM_CW0, iWel, iRad, iVtx, cZERO, rZERO)
          else
            if (iVtx < wel%iResolution) then
              iEQ = MAT_CreateEquation(io, mat, (/this%cCPZ, next%cCPZ/), EQN_POTENTIALDIFF, &
                                       ELEM_CW0, iWel, iRad, iVtx, cZERO, rZERO)
            else
              if (iRad < wel%iNRad) then
                ! For all but the last radial, make a delta-Phi equation...
                iEQ = MAT_CreateEquation(io, mat, (/this%cCPZ, last_vtx%cCPZ/), EQN_POTENTIALDIFF, &
                                         ELEM_CW0, iWel, iRad, iVtx, cZERO, rZERO)
              else
                ! For the last radial, make a total discharge equation
                iEQ =  MAT_CreateEquation(io, mat, (/cZERO/), EQN_TOTALFLOW, &
                                          ELEM_CW0, iWel, iRad, iVtx, cZERO, rZERO)
              end if
            end if
          end if
        end do
      end do

    end do

    return
  end subroutine CW0_SetupMatrix


  function iCW0_Prepare(io, cw0, aqu, iIteration) result(iChanges)
    !! subroutine CW0_Prepare
    !!
    !! Prepares the module for a new iteration
    !!
    !! Do-nothing for m_cw0
    !!
    !! Calling Sequence:
    !!    call CW0_Setup(io, cw0, aqu, mat)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer
    !!             CW0_COLLECTION object to be used
    !!   (in)    type(MAT_MATRIX), pointer
    !!             MAT_MATRIX object to be used
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iIteration
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iChanges

    iChanges = 0

    return
  end function iCW0_Prepare


  function rCW0_GetCoefficientMultiplier(io, cw0, iElementString, iElementVertex, &
             iElementFlag) result(rMultiplier)
    !! Returns the coefficient multiplier
    !! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    type(IO_STATUS), pointer :: io
    !! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier

    rMultiplier = rONE

    return
  end function rCW0_GetCoefficientMultiplier


  subroutine CW0_SetRegenerate(io, cw0, iRegenerate)
    !! Sets the regenerate flag to the provided value
    type(CW0_COLLECTION), pointer :: cw0
    integer(kind=AE_INT), intent(in) :: iRegenerate
    type(IO_STATUS), pointer :: io

    cw0%iRegenerate = iRegenerate

    return
  end subroutine CW0_SetRegenerate


  subroutine CW0_ComputeCoefficients(io, cw0, aqu, fwl, fdp, cPathZ, iEqType, iElementType, iElementString, &
               iElementVertex, iElementFlag, cOrientation, rGhbResistance, &
               iIteration, rMultiplier, rARow)
    !! subroutine CW0_ComputeCoefficients
    !!
    !! Computes a row of matrix coefficients(with no corrections) for the CW0
    !! elements in layer iL.
    !!
    !! Calling Sequence:
    !!    call CW0_ComputeCoefficients(io, cw0, cPathZ, iEqType, cOrientation, rRow)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer
    !!             CW0_COLLECTION object to be used
    !!   (in)    type(FWL_COLLECTION), pointer
    !!             FWL_COLLECTION object to be used
    !!   (in)    type(FDP_COLLECTION), pointer
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
    !!   (in)    integer :: iIteration
    !!             The element flag(if any) for this equation
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(AQU_COLLECTION), pointer :: aqu
    type(FWL_COLLECTION), pointer :: fwl
    type(FDP_COLLECTION), pointer :: fdp
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    complex(kind=AE_REAL), intent(in) :: cOrientation
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    real(kind=AE_REAL), intent(in) :: rGhbResistance
    integer(kind=AE_INT), intent(in) :: iIteration
    real(kind=AE_REAL), intent(in) :: rMultiplier
    real(kind=AE_REAL), dimension(:), intent(out) :: rARow
    type(IO_STATUS), pointer :: io
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat, iCol, iWel, iVtx, iDP1, iNDP, iWhich, iVtxCol, &
                           iLastCol, iNextCol, iBaseCol
    complex(kind=AE_REAL), dimension(:, :, :), allocatable :: cDPF, cDP, cDPW
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad
    type(CW0_VERTEX), pointer :: first_vtx, this, next, last_vtx
    real(kind=AE_REAL) :: rT1, rT2

    if (io%lDebug) then
      call IO_Assert(io, (associated(cw0)), &
           "CW0_ComputeCoefficients: CW0_Create has not been called")
      call IO_Assert(io, (associated(fwl)), &
           "CW0_Setup: Illegal FWL_COLLECTION object")
      call IO_Assert(io, (associated(fdp)), &
           "CW0_ComputeCoefficients: Illegal FDP_COLLECTION object")
    end if

    iCol = 0
    rARow = rZERO
    do iWel = 1, cw0%iCount
      wel => cw0%Wells(iWel)
      rad => wel%Radials(wel%iNRad)
      last_vtx => rad%Vertices(wel%iResolution)
      first_vtx => wel%Radials(1)%Vertices(1)
      iDP1 = first_vtx%iFDPIndex
      iNDP = wel%iNRad * wel%iResolution
      iBaseCol = iCol

      allocate(cDPF(0:iNDP+1, 1, 1), stat = iStat)
      call IO_Assert(io, (iStat == 0), "CW0_ComputeCoefficients: Allocation failed")

      ! Get the appropriate incluence functions for the boundary condition type
      select case (iEqType)
        case (EQN_HEAD)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_BDYGHB)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, (/rHALF*sum(cPathZ)/), cOrientation, cDPF(1:iNDP, :, :))
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPW(1:iNDP, :, :))
          cDPF = cDPF + rGhbResistance*cDPW
        case (EQN_FLOW)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_F, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_INHO)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_P, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_DISCHARGE)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_W, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_RECHARGE)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_G, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_CONTINUITY)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_Q, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
        case (EQN_TOTALFLOW)
          if (iElementType == ELEM_CW0 .and. iElementString == iWel) then
            call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_Q, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
          else
            call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_Z, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
          end if
        case (EQN_POTENTIALDIFF)
          call FDP_GetInfluence_ILS(io, fdp, INFLUENCE_D, iDP1, iNDP, cPathZ, cOrientation, cDPF(1:iNDP, :, :))
      end select

      do iVtx = 1, iNDP
        iCol = iCol + 1
        rARow(iCol) = real(cDPF(iVtx, 1, 1))
      end do

      ! Look up the matrix column for the specified vertex and adjust the coefficients
      if (iElementType == ELEM_CW0 .and. iEqType == EQN_POTENTIALDIFF .and. iElementString == iWel) then
        rad => wel%Radials(iElementVertex)
        this => rad%Vertices(iElementFlag)
        iVtxCol = iBaseCol + this%iFDPIndex - iDP1 + 1
        if (iElementFlag /= wel%iResolution) then
          next => rad%Vertices(iElementFlag+1)
          iNextCol = iVtxCol+1
          rT1 = rAQU_Transmissivity(io, aqu, this%cCPZ, this%rCheck)
          rT2 = rAQU_Transmissivity(io, aqu, next%cCPZ, next%rCheck)
          this%rResistanceTerm = (rHALF * (rT1+rT2) * rad%rResistance) / rad%rWidth
          rARow(iVtxCol) = rARow(iVtxCol) - this%rResistanceTerm
          rARow(iNextCol) = rARow(iNextCol) + this%rResistanceTerm
        else
          next => last_vtx
          iNextCol = iCol
          rT1 = rAQU_Transmissivity(io, aqu, this%cCPZ, this%rCheck)
          rT2 = rAQU_Transmissivity(io, aqu, next%cCPZ, next%rCheck)
          this%rResistanceTerm = (rHALF * (rT1+rT2) * rad%rResistance) / rad%rWidth
          rARow(iVtxCol) = rARow(iVtxCol) - this%rResistanceTerm
          rARow(iNextCol) = rARow(iNextCol) + this%rResistanceTerm
        end if
      else if (iElementType == ELEM_CW0 .and. iEqType == EQN_TOTALFLOW .and. iElementString == iWel) then
        rT1 = rAQU_Transmissivity(io, aqu, last_vtx%cCPZ, last_vtx%rCheck)
        last_vtx%rResistanceTerm = (rT1 * rad%rResistance) / rad%rWidth
      else if (iElementType == ELEM_CW0 .and. iEqType == EQN_HEAD .and. iElementString == iWel) then
        rad => wel%Radials(iElementVertex)
        this => rad%Vertices(iElementFlag)
        rT1 = rAQU_Transmissivity(io, aqu, this%cCPZ, this%rCheck)
        this%rResistanceTerm = (rT1 * rad%rResistance) / rad%rWidth
        iVtxCol = iBaseCol + this%iFDPIndex - iDP1 + 1
        rARow(iVtxCol) = rARow(iVtxCol) - this%rResistanceTerm
      end if

      deallocate(cDPF)
    end do

    ! Use the multiplier
    rARow = rMultiplier * rARow

    return
  end subroutine CW0_ComputeCoefficients


  function rCW0_ComputeRHS(io, cw0, aqu, iEqType, iElementType, iElementString, iElementVertex, &
             iElementFlag, iIteration, lDirect) result(rRHS)
    !! function rCW0_ComputeRHS
    !!
    !! Computes the right-hand side value for the solution
    !!
    !! Calling Sequence:
    !!   rRHS = rCW0_ComputeRHS(io, cw0, rValue, iElementType, iElementString, iElementVertex, &
         !!                          iElementFlag)
    !!
    !! Arguments:
    !!   (in)    type(CW0_COLLECTION), pointer :: cw0
    !!             CW0_COLLECTION object to be used
    !!   (in)    type(AQU_COLLECTION), pointer :: aqu
    !!             The AQU_COLLECTION used for head-to-potential conversions
    !!   (in)    integer :: iElementType
    !!             Element type(either ELAM_AQU or ELEM_IN0)
    !!   (in)    integer :: iElementString
    !!             Element string number
    !!   (in)    integer :: iElementVertex
    !!             Element vertex number
    !!   (in)    integer :: iElementFlag
    !!             Element flag(e.g. for vertices which yield more than one equation)
    !!   (in)    type(IO_STATUS), pointer :: io
    !!             Tracking for IO module
    !!
    !! Return Value:
    !!   real :: rRHS
    !!     The RHS value for the module
    !!
    ! [ ARGUMENTS ]
    type(CW0_COLLECTION), pointer :: cw0
    type(AQU_COLLECTION), pointer :: aqu
    integer(kind=AE_INT), intent(in) :: iEqType
    integer(kind=AE_INT), intent(in) :: iElementType
    integer(kind=AE_INT), intent(in) :: iElementString
    integer(kind=AE_INT), intent(in) :: iElementVertex
    integer(kind=AE_INT), intent(in) :: iElementFlag
    integer(kind=AE_INT), intent(in) :: iIteration
    logical, intent(in) :: lDirect
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRHS
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rThisPot, rNextPot
    type(CW0_WELL), pointer :: wel
    type(CW0_RADIAL), pointer :: rad
    type(CW0_VERTEX), pointer :: this, next, last_vtx

    call IO_Assert(io, (associated(cw0)), "CW0_ComputeRHS: CW0_COLLECTION not available")
    call IO_Assert(io, (iElementString > 0 .and. iElementString <= cw0%iCount), &
         "CW0_ComputeRHS: Bad well index")

    wel => cw0%Wells(iElementString)
    rad => wel%Radials(wel%iNRad)
    last_vtx => rad%Vertices(wel%iResolution)

    call IO_Assert(io, (iElementVertex > 0 .and. iElementVertex <= wel%iNRad), &
         "CW0_ComputeRHS: Bad radial index")
    rad => wel%Radials(iElementVertex)

    call IO_Assert(io, (iElementFlag > 0 .and. iElementFlag <= wel%iResolution), &
         "CW0_ComputeRHS: Bad vertex index")
    this => rad%Vertices(iElementFlag)

    if (iEqType == EQN_POTENTIALDIFF) then
      if (iElementFlag < wel%iResolution-1) then
        next => rad%Vertices(iElementFlag+1)
      else
        next => wel%Radials(wel%iNRad)%Vertices(wel%iResolution)
      end if
      ! Now, add the head correction terms, the pipe loss(if implemented), and convert to
      ! potentials; use these potentials to compute the delta-potential value for the RHS
      if (iElementFlag < wel%iResolution) then
        if (lDirect) then
          rRHS = rZERO
        else
          rRHS = this%rResistanceTerm * (this%rStrength-next%rStrength) - &
                 (this%rCheck-next%rCheck)
          endif
        else
          if (lDirect) then
            rRHS = rZERO
          else
            rRHS = this%rResistanceTerm * (this%rStrength-last_vtx%rStrength) - &
                   (this%rCheck-last_vtx%rCheck)
          end if
        end if
        !    print *, 'CW0 DH RHS ', iElementVertex, iElementString, rRHS
      else if (iEqType == EQN_HEAD) then
        ! Closure condition -- specified head
        if (lDirect) then
          rRHS = rAQU_HeadToPotential(io, aqu, wel%rQ, this%cCPZ)
        else
          rRHS = rAQU_HeadToPotential(io, aqu, wel%rQ, this%cCPZ) - this%rCheck
          rRHS = rRHS + this%rResistanceTerm * this%rStrength
        end if
      else if (iEqType == EQN_TOTALFLOW) then
        ! Closure condition -- specified pumping rate
        if (lDirect) then
          rRHS = wel%rQ
        else
          rRHS = wel%rQ - rCW0_ComputeTotalFlow(io, cw0, iElementString)
        end if
      end if

      return
    end function rCW0_ComputeRHS


    subroutine CW0_StoreResult(io, cw0, rValue, iElementType, iElementString, iElementVertex, iElementFlag, lDirect)
      !! subroutine CW0_StoreResult
      !!
      !! Stores the results of a solution for a single equation associated with
      !! the CW0 module.
      !!
      !! Calling Sequence:
      !!    CW0_StoreResult(io, cw0, cCPZ, iEqType, cOrientation, rRHS)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!   (in)    real :: rValue
      !!             The new result value from the solution vector
      !!   (in)    integer :: iElementType
      !!             Element type(always ELEM_CW0)
      !!   (in)    integer :: iElementString
      !!             Element string number
      !!   (in)    integer :: iElementVertex
      !!             Element vertex number
      !!   (in)    integer :: iElementFlag
      !!             Element flag(e.g. for vertices which yield more than one equation)
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      real(kind=AE_REAL), intent(in) :: rValue
      integer(kind=AE_INT), intent(in) :: iElementType
      integer(kind=AE_INT), intent(in) :: iElementString
      integer(kind=AE_INT), intent(in) :: iElementVertex
      integer(kind=AE_INT), intent(in) :: iElementFlag
      logical, intent(in) :: lDirect
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: this, next

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_StoreResult: CW0_Create has not been called")
        call IO_Assert(io, (iElementString >= 1 .and. iElementString <= cw0%iCount), &
             "CW0_StoreResult: Bad element string ID")
      end if

      wel => cw0%Wells(iElementString)

      if (io%lDebug) then
        call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= wel%iNRad), &
             "CW0_StoreResult: Bad radial ID")
      end if
      rad => wel%Radials(iElementVertex)

      if (io%lDebug) then
        call IO_Assert(io, (iElementFlag >= 1 .and. iElementFlag <= wel%iResolution), &
             "CW0_StoreResult: Bad segment ID")
      end if

      this => rad%Vertices(iElementFlag)
      next => rad%Vertices(iElementFlag+1)
      if (lDirect) then
        this%rStrength = rValue
      else
        this%rStrength = this%rStrength + rValue
      end if
      next%rDPStrength = this%rDPStrength + this%rStrength * this%rLength

      return
    end subroutine CW0_StoreResult


    subroutine CW0_ComputeCheck(io, cw0, aqu, iEqType, iElementString, iElementVertex, iElementFlag, rHead)
      !! function rCW0_ComputeCheck
      !!
      !! Returns the check value for the specified domain and vertex at the point cZ
      !!
      !! Calling Sequence:
      !!    call CW0_ComputeCheck(in0, iElementString, iElementVertex, iElementFlag, rHead)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION object to be used
      !!   (in)    integer :: iElementString
      !!             The radial collector index
      !!   (in)    integer :: iElementVertex
      !!             The vertex index
      !!   (in)    integer :: iElementFlag
      !!             The flag value
      !!   (in)    integer :: iEqType
      !!             The equation type(EQ_POTENTIALDIFF or EQ_TOTALFLOW)
      !!   (in)    real :: rHead1
      !!             Head at this segment's control point
      !!   (in)    real :: rHead2
      !!             Head at the next segment's control point
      !!   (in)    type(IO_STATUS) :: io
      !!             Tracking for IO module
      !!
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(AQU_COLLECTION), pointer :: aqu
      integer(kind=AE_INT), intent(in) :: iEqType
      integer(kind=AE_INT), intent(in) :: iElementString
      integer(kind=AE_INT), intent(in) :: iElementVertex
      integer(kind=AE_INT), intent(in) :: iElementFlag
      real(kind=AE_REAL), intent(in) :: rHead
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: this, next

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_ComputeCheck: CW0_Create has not been called")
        call IO_Assert(io, (iElementString >= 1 .and. iElementString <= cw0%iCount), &
             "CW0_ComputeCheck: Bad element string ID")
      end if
      wel => cw0%Wells(iElementString)

      if (io%lDebug) then
        call IO_Assert(io, (iElementVertex >= 1 .and. iElementVertex <= wel%iNRad), &
             "CW0_ComputeCheck: Bad radial ID")
      end if
      rad => wel%Radials(iElementVertex)

      if (io%lDebug) then
        call IO_Assert(io, (iElementFlag >= 1 .and. iElementFlag <= wel%iResolution), &
             "CW0_ComputeCheck: Bad segment ID")
      end if

      ! Note: this code only stores the modeled head at the segment; CW0_ComputeRHS computes the
      ! actual delta-potential information at the beginning of the next iteration
      this => rad%Vertices(iElementFlag)
      this%rCheck = rHead

      if (iEqType == EQN_TOTALFLOW) then
        wel%rCheck = rCW0_ComputeTotalFlow(io, cw0, iElementString)
      end if

      return
    end subroutine CW0_ComputeCheck


    function rCW0_ComputeTotalFlow(io, cw0, iElementString) result(rCheck)
      !! function rCW0_ComputeTotalFlow
      !!
      !! Returns the total flow for the well with index iElementString
      !!
      !! Calling Sequence:
      !!    rCheck = rCW0_ComputeTotalFlow(in0, iElementString)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION object to be used
      !!   (in)    integer :: iElementString
      !!             The index of the well
      !!
      !! Return Value:
      !!   (out)   real :: rCheck
      !!             The check value
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      integer(kind=AE_INT), intent(in) :: iElementString
      type(IO_STATUS), pointer :: io
      ! [ RETURN VALUE ]
      real(kind=AE_REAL) :: rCheck
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iRad, iVtx
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_ComputeCheck: CW0_Create has not been called")
      end if

      rCheck = rZERO
      wel => cw0%Wells(iElementString)
      do iRad = 1, wel%iNRad
        rad => wel%Radials(iRad)
        do iVtx = 1, wel%iResolution
          vtx => rad%Vertices(iVtx)
          rCheck = rCheck + vtx%rLength * vtx%rStrength
        end do
      end do

      return
    end function rCW0_ComputeTotalFlow


    function rCW0_ComputeCaissonHead(io, cw0, aqu, iWel) result(rCaissonHead)
      !! function rCW0_ComputeCaissonHead
      !!
      !! Returns the head in the caisson
      !!
      !! Calling Sequence:
      !!    h = rCW0_ComputeCaissonHead(io, cw0, aqu, iWel)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!   (in)    type(AQU_COLLECTION), pointer
      !!             AQU_COLLECTION object for head/potential conversion
      !!   (in)    integer :: iLU
      !!             The output LU to receive output
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(AQU_COLLECTION), pointer :: aqu
      integer(kind=AE_INT), intent(in) :: iWel
      type(IO_STATUS), pointer :: io
      ! [ RETURN VALUE ]
      real(kind=AE_REAL) :: rCaissonHead
      ! [ LOCALS ]
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "rCW0_ComputeCaissonHead: CW0_Create has not been called")
      end if

      wel => cw0%Wells(iWel)
      rad => wel%Radials(wel%iNRad)
      vtx => rad%Vertices(wel%iResolution)
      rCaissonHead = rAQU_PotentialToHead(io, aqu, vtx%rCheck, vtx%cCPZ) - vtx%rStrength * rad%rResistance / rad%rWidth

      return
    end function rCW0_ComputeCaissonHead


    function CW0_EnableDrawdown(io, cw0) result(iChanges)
      !! Enables the elements that are subject to the "drawdown" flag
      !! This allows the modeler to compute "drawdown" simulations based on the difference between
      !! an "unstressed" condition and a "stressed" condition. Each element module that implements
      !! the DDN operator provides an XXX_EnableDrawdown() subroutine that makes the necessary
      !! adjustments.
      !!
      !! Returns the number of changes made when enabling the drawdown simulation
      !!
      ! [ ARGUMENTS ]
      type(IO_STATUS), pointer :: io
      type(CW0_COLLECTION), pointer :: cw0
      ! [ RETURN VALUE ]
      integer(kind=AE_INT) :: iChanges
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iwel
      type(CW0_WELL), pointer :: wel

      iChanges = 0
      do iwel = 1, cw0%iCount
        wel => cw0%Wells(iwel)
        if (wel%lDdn) then
          iChanges = iChanges + 1
          if (wel%lDdnHeadSpec) then
            wel%rQ = wel%rDdnValue
            wel%lHeadSpec = .true.
          else
            wel%rQ = wel%rDdnValue
            wel%lHeadSpec = .false.
          end if
        end if
      end do

      return
    end function CW0_EnableDrawdown


    subroutine CW0_Update(io, cw0, fwl, fdp)
      !! subroutine CW0_Update
      !!
      !! Updates the underlying function objects for the specified layer.
      !!
      !! Calling Sequence:
      !!    CW0_Update(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!   (in)    type(FWL_COLLECTION), pointer
      !!             FWL_COLLECTION object to be used
      !!   (in)    type(FDP_COLLECTION), pointer
      !!             FDP_COLLECTION object to be used
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(FWL_COLLECTION), pointer :: fwl
      type(FDP_COLLECTION), pointer :: fdp
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iWel, iRad, iVtx
      complex(kind=AE_REAL) :: cRho1, cRho2, cRho3
      real(kind=AE_REAL) :: rTotalExtraction
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: this, next


      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_Update: CW0_Create has not been called")
        call IO_Assert(io, (associated(fdp)), &
             "CW0_Update: Illegal FDP_COLLECTION object")
      end if

      do iWel = 1, cw0%iCount
        rTotalExtraction = rZERO
        wel => cw0%Wells(iWel)
        do iRad = 1, wel%iNRad
          rad => wel%Radials(iRad)
          do iVtx = 1, wel%iResolution
            this => rad%Vertices(iVtx)
            next => rad%Vertices(iVtx+1)
            cRho1 = cmplx(this%rDPStrength, rZERO, AE_REAL)
            cRho3 = cmplx(next%rDPStrength, rZERO, AE_REAL)
            cRho2 = rHALF * (cRho1 + cRho3)
            call FDP_Update(io, fdp, this%iFDPIndex, (/cRho1, cRho2, cRho3/))
            rTotalExtraction = rTotalExtraction + this%rLength*this%rStrength
          end do
          ! Put a well at the end of the string
          call FWL_Update(io, fwl, rad%iFWLIndex, real(cRho3))
        end do
        wel%rCheck = rTotalExtraction
        !    print *, 'ext ', iWel, rTotalExtraction
      end do

      return
    end subroutine CW0_Update


    function lCW0_CheckPoint(io, cw0, cZ, rTol) result(lSing)
      !! logical function lCW0_CheckPoint
      !!
      !! Checks for a singularity near cZ
      !!
      !! Calling Sequence:
      !!    lSing = lCW0_CheckPoint(io, cw0, cZ, rTol)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             The CW0_COLLECTION object to be used
      !!   (in)    complex :: cZ
      !!             The point to be checked
      !!   (in)    real :: rTol
      !!             Tolerance about cZ to be checked
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      complex(kind=AE_REAL), intent(in) :: cZ
      real(kind=AE_REAL), intent(in) :: rTol
      type(IO_STATUS), pointer :: io
      ! [ RETURN VALUE ]
      logical :: lSing
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iWel, iRad, iVtx
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      ! NOTE: Singularities are found at the end points of the linesinks making up the radials
      lSing = .false.
      do iWel = 1, cw0%iCount
        wel => cw0%Wells(iWel)
        do iRad = 1, wel%iNRad
          rad => wel%Radials(iRad)
          do iVtx = 1, wel%iResolution+1
            vtx => rad%Vertices(iVtx)
            if (abs(real(vtx%cZ-cZ, AE_REAL)) < rTol .and. abs(aimag(vtx%cZ-cZ)) < rTol) then
              lSing = .true.
              exit
            end if
          end do
        end do
      end do

      return
    end function lCW0_CheckPoint


    subroutine CW0_FindStringPointer(io, cw0, iCWID, CWWell, lfound)
      !! subroutine CW0_FindStringPointer
      !!
      !! Finds the linesink string specified by the ID and returns a pointer to it
      !!
      !! Calling Sequence:
      !!    call CW0_FindStringPointer(ls1, iLSID, LSString, lfound)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be used
      !!   (in)    integer :: iCWID
      !!             The linesink string ID number
      !!   (out)   type(CW0_WELL) :: CWWell
      !!             Pointer to the well
      !!   (out)   logical :: lFound
      !!             .true. if the well was found
      !!             .false. if the well was not found
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      integer(kind=AE_INT), intent(in) :: iCWID
      type(CW0_WELL), pointer :: CWWell
      logical :: lFound
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      integer(kind=AE_INT) :: i

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_FindStringPointer: CW0_Create has not been called")
      end if

      lFound = .false.
      do i = 1, cw0%iCount
        CWWell => cw0%Wells(i)
        if (CWWell%iID == iCWID) then
          lFound = .true.
          return
        end if
      end do

      return
    end subroutine CW0_FindStringPointer


    function lCW0_CheckProximity(io, cw0, cZ, rProximity, iDir, iElementID, iElementVtx) result(lFound)
      !! function lCW0_CheckProximity
      !!
      !! Checks to see if the point specified is near a line-sink, then
      !! moves the point a small distance if it is.  Sets lChange = .true. if a
      !! change is made to cZ.  NOTE: lChange is preset by the caller, so this
      !! routine only changes it to .true. if necessary.
      !!
      !! An example application is found in cAEM_Potential()
      !!
      !! Calling Sequence:
      !!    lFound = CW0_CheckProximity(io, cw0, cZ, rProximity, iDir, iElementID) result(lFound)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!   (in)    complex :: cZ
      !!             The point to be checked
      !!   (in)    real :: rProximity
      !!             The tolerance about cZ to be examined
      !!   (in)    integer :: iDir
      !!             Flag for checking line-sinks:
      !!               iDir = 0   -  Check without regard to the direction of travel
      !!               iDir = < 0   -  Check for back tracing(ignore gaining line-sinks)
      !!               iDir => 0   -  Check for formard tracing(ignore losing line-sinks)
      !!   (out)   integer :: iElementID
      !!             The element(string) ID for the line-sink found(if lFound is .true.)
      !!   (out)   integer :: iElementVtx
      !!             The vertex associated with the string segment which was found
      !!   (out)   logical :: lFound
      !!             .true. if a line-sink is found
      !!             .false. if a line-sink is not found
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      complex(kind=AE_REAL), intent(in) :: cZ
      integer(kind=AE_INT), intent(in) :: iDir
      real(kind=AE_REAL), intent(in) :: rProximity
      integer(kind=AE_INT), intent(out) :: iElementID, iElementVtx
      type(IO_STATUS), pointer :: io
      ! [ RETURN VALUE ]
      logical :: lFound
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iWel, iRad, iVtx
      real(kind=AE_REAL) :: rDist, rMinDist
      complex(kind=AE_REAL) :: cZL, cZC, cBigZ
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: this, next

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_CheckProximity: CW0_Create has not been called")
      end if

      iElementID = 0
      lFound = .false.
      rMinDist = HUGE(AE_REAL)
      do iWel = 1, cw0%iCount
        wel => cw0%Wells(iWel)
        do iRad = 1, wel%iNRad
          rad => wel%Radials(iRad)
          do iVtx = 1, wel%iResolution
            this => rad%Vertices(iVtx)
            next => rad%Vertices(iVtx+1)
            ! iDir < 0 -- skip gaining line sinks
            if (iDir < 0 .and. this%rStrength > rZERO) cycle
            ! iDir > 0 -- skip losing line sinks
            if (iDir > 0 .and. this%rStrength < rZERO) cycle
            ! Compute the mapped Z-value and then the distance from the linesink
            cZC = rHALF * (this%cZ + &
                  next%cZ)
            cZL = rHALF * (this%cZ - &
                  next%cZ)
            cBigZ = (cZ-cZC) / cZL
            rDist = abs(aimag(cBigZ) * abs(cZL))
            if (rDist < rProximity .and. abs(real(cBigZ)) <= rONE) then
              if (rDist < rMinDist) then
                rMinDist = rDist
                iElementID = wel%iID
                iElementVtx = iVtx
                lFound = .true.
              end if
            end if
          end do
        end do
      end do

      return
    end function lCW0_CheckProximity


    subroutine CW0_ResetIterator(io, cw0)
      !! subroutine CW0_ResetIterator
      !!
      !! Resets the module's iterator prior to traversing for check data
      !!
      !! Calling Sequence:
      !!    call CW0_ResetIterator(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be used
      !!   (in)    type(IO_STATUS), pointer :: cw0
      !!             Tracks error conditions
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(IO_STATUS), pointer :: io

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_ResetIterator: CW0_Create has not been called")
      end if

      cw0%iIterWell = 1
      cw0%iIterRad = 1
      cw0%iIterVtx = 0

      return
    end subroutine CW0_ResetIterator


    function CW0_NextIterator(io, cw0) result(itr)
      !! function CW0_NextIterator
      !!
      !! Advances the module's iterator one step
      !!
      !! Calling Sequence:
      !!    call CW0_NextIterator(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be used
      !!   (in)    type(IO_STATUS), pointer :: cw0
      !!             Tracks error conditions
      !!
      !! Return Value:
      !!   type(ITERATOR_RESULT), pointer :: itr
      !!     Pointer to the information for data retrieval
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(IO_STATUS), pointer :: io
      ! [ RETURN VALUE ]
      type(ITERATOR_RESULT), pointer :: itr
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_NextIterator: CW0_Create has not been called")
      end if

      if (cw0%iIterWell > cw0%iCount) then
        nullify(itr)
        return
      end if

      wel => cw0%Wells(cw0%iIterWell)
      rad => wel%Radials(cw0%iIterRad)
      cw0%iIterVtx = cw0%iIterVtx + 1
      if (cw0%iIterVtx > wel%iResolution) then
        cw0%iIterRad = cw0%iIterRad + 1
        cw0%iIterVtx = 1
        if (cw0%iIterRad > wel%iNRad) then
          cw0%iIterWell = cw0%iIterWell + 1
          cw0%iIterRad = 1
          cw0%iIterVtx = 1
          if (cw0%iIterWell > cw0%iCount) then
            nullify(itr)
            return
          end if
          wel => cw0%Wells(cw0%iIterWell)
        end if
        rad => wel%Radials(cw0%iIterRad)
      end if
      vtx => rad%Vertices(cw0%iIterVtx)

      allocate(itr)

      itr%iElementType = ELEM_CW0
      itr%iElementString = cw0%iIterWell
      itr%iElementVertex = cw0%iIterRad
      itr%iElementFlag = cw0%iIterVtx
      itr%iValueSelector = VALUE_POTENTIAL
      allocate(itr%cZ(1))
      itr%cZ(1) = vtx%cCPZ

      return
    end function CW0_NextIterator


    subroutine CW0_SetIterator(io, cw0, itr, cValue)
      !! function CW0_SetIterator
      !!
      !! Advances the module's iterator one step
      !!
      !! Calling Sequence:
      !!    call CW0_SetIterator(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be used
      !!   (in)    type(ITERATOR_RESULT), pointer :: itr
      !!             Pointer to the information for data retrieval
      !!   (in)    complex :: cValue
      !!             The value retrieved from the color
      !!   (in)    type(IO_STATUS), pointer :: cw0
      !!             Tracks error conditions
      !!
      !! Return Value:
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(ITERATOR_RESULT), pointer :: itr
      complex(kind=AE_REAL), intent(in) :: cValue
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_NextIterator: CW0_Create has not been called")
      end if

      cw0%Wells(itr%iElementString)%Radials(itr%iElementVertex)%Vertices(itr%iElementFlag)%rCheck = real(cValue, AE_REAL)

      return
    end subroutine CW0_SetIterator


    subroutine CW0_FindWell(io, cw0, iWellID, cZWell, rDischarge, lFound)
      !! subroutine CW0_FindWell
      !!
      !! Finds the well specified by the Well ID and returns its parameters
      !!
      !! Calling Sequence:
      !!    call CW0_FindWell(wl0, iWellID, cZWell, rDischarge)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be used
      !!   (in)    integer :: iWellID
      !!             The well ID number
      !!   (out)   complex :: cZWell
      !!             Location of the well caisson center
      !!   (out)   real :: rDischarge
      !!             The discharge of the well
      !!   (out)   logical :: lFound
      !!             .true. if the well was found
      !!             .false. if the well was not found
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      integer(kind=AE_INT), intent(in) :: iWellID
      complex(kind=AE_REAL), intent(out) :: cZWell
      real(kind=AE_REAL), intent(out) :: rDischarge
      logical :: lFound
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      integer(kind=AE_INT) :: i
      type(CW0_WELL), pointer :: wel

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_FindWell: CW0_Create has not been called")
      end if

      lFound = .false.
      do i = 1, cw0%iCount
        wel => cw0%Wells(i)
        if (wel%iID == iWellID) then
          lFound = .true.
          cZWell = wel%cZC
          rDischarge = wel%rQ
        end if
      end do

      return
    end subroutine CW0_FindWell


    subroutine CW0_Read(io, cw0)
      !! subroutine CW0_Read
      !!
      !! Populates an CW0_COLLECTION using data from LU_INPUT
      !!
      !! Calling Sequence:
      !!    call CW0_Read(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer :: cw0
      !!             CW0_COLLECTION to be populated
      !!
      !! The format of the CW0 section of the input file appears as follows:
      !! cw0
      !!   wel max_radials zc rc q c
      !!     (x, y)     # tip of radial #1
      !!     (x, y)     # tip of radial #2
      !!     ...  up to max_radials
      !!   wel...
      !!
      !! NOTE: It is assumed that the CW0 line was found by the caller

      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(IO_STATUS), pointer :: io
      ! [ LOCAL DIRECTIVES ]
      type(DIRECTIVE), dimension(3), parameter :: dirDirectives = (/ dirEND, dirWEL, dirDDN /)
      ! [ LOCALS ]
      real(kind=AE_REAL) :: rRc
      real(kind=AE_REAL) :: rQ
      real(kind=AE_REAL) :: rResistance
      real(kind=AE_REAL) :: rWidth
      real(kind=AE_REAL) :: rBlankLength
      real(kind=AE_REAL) :: rDdnValue
      complex(kind=AE_REAL) :: cZ
      complex(kind=AE_REAL) :: cZC
      integer(kind=AE_INT) :: iOpCode
      integer(kind=AE_INT) :: iStat
      integer(kind=AE_INT) :: iMaxRad
      integer(kind=AE_INT) :: iID
      integer(kind=AE_INT) :: iResolution
      logical :: lFlag, lHeadSpec, lDdnHeadSpec
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad

      call IO_MessageText(io, "  Reading CW0 module input")

      call IO_Assert(io, (associated(cw0)), "CW0_Read: CW0_Create has not been called")

      ! Use IO_InputRecord to process the model input file.
      nullify(wel, rad)
      do
        call IO_InputRecord(io, dirDirectives, iOpCode)
        select case (iOpCode)
          case (kOpError)
            ! A RunTime error was found during a file read operation. This
            ! condition is fatal; warn the user, and exit.
            call IO_Assert(io, .false., "CW0_Read: I/O Error")
          case (kOpFileEOF)
            ! EOF is unexpected for all ModCW0 "ifXXXRead" routines.
            call IO_Assert(io, .false., "CW0_Read: Unexpected EOF")
          case (kOpData)
            !****************************************************************************
            ! Here for data records
            !****************************************************************************
            call IO_Assert(io, (associated(wel)), "CW0_Read: No current well")
            call IO_Assert(io, (wel%iNRad < size(wel%Radials)), "CW0_Read: Space exhausted")
            cZ = cIO_GetCoordinate(io, "cZ")
            rWidth = rIO_GetReal(io, "rWidth", minimum=rTINY, def=rONE)
            rResistance = rIO_GetReal(io, "rResistance", minimum=rZERO)
            rBlankLength = rIO_GetReal(io, "rBlankLength", minimum=rZERO)
            wel%iNRad = wel%iNRad+1
            rad => wel%Radials(wel%iNRad)
            rad%cZEnd = cZ
            rad%rWidth = rWidth
            rad%rResistance = rResistance
            rad%rBlankLength = rBlankLength
          case (kOpEND)
            ! EOD mark was found. Exit the file parser.
            exit
          case (kOpWEL)
            !****************************************************************************
            ! Here for the WEL command -- create a new collector well
            ! the maximum number of radials is in the input record
            !****************************************************************************
            call IO_Assert(io, (cw0%iCount < size(cw0%Wells)), "CW0_Read: Space exhausted")
            ! Retrive the number of vertices desired...
            iMaxRad = iIO_GetInteger(io, "iMaxRad", minimum=1)
            cZC = cIO_GetCoordinate(io, "cZC")
            rRc = rIO_GetReal(io, "rRc", minimum=rTINY)
            rQ = rIO_GetReal(io, "rQ")
            iResolution = iIO_GetInteger(io, "iResolution", minimum=1)
            lHeadSpec = lIO_GetLogical(io, "lHeadSpec")
            iID = iIO_GetInteger(io, "iID", forbidden=cw0%Wells(:)%iID)
            ! OKAY! Allocate the vertices...
            cw0%iCount = cw0%iCount+1
            wel => cw0%Wells(cw0%iCount)
            allocate(wel%Radials(iMaxRad), stat = iStat)
            call IO_Assert(io, (iStat == 0), "CW0_Read: Allocation failed")
            ! Made it!
            wel%iNRad = 0
            wel%cZC = cZC
            wel%rRc = rRc
            wel%rQ = rQ
            wel%iResolution = iResolution
            wel%lHeadSpec = lHeadSpec
            wel%iID = iID
          case (kOpDDN)
            call IO_Assert(io, (cw0%iCount>0), "No current well to modify")
            wel => cw0%Wells(cw0%iCount)
            rDdnValue = rIO_GetReal(io, "rDdnValue")
            lDdnHeadSpec = lIO_GetLogical(io, "lDdnHeadSpec", def=.false.)
            wel%lDdn = .true.
            wel%rDdnValue = rDdnValue
            wel%lDdnHeadSpec = lDdnHeadSpec
        end select
      end do

      call IO_MessageText(io, "  Leaving CW0 module")

    end subroutine CW0_Read


    subroutine CW0_Inquiry(io, cw0, aqu, iLU)
      !! subroutine CW0_Inquiry
      !!
      !! Writes an inquiry report for all line-sinks to iLU
      !!
      !! Calling Sequence:
      !!    call CW0_Inquiry(io, cw0, iLU)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!   (in)    integer :: iLU
      !!             The output LU to receive output
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(AQU_COLLECTION), pointer :: aqu
      integer(kind=AE_INT), intent(in) :: iLU
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iWel, iRad, iVtx
      real(kind=AE_REAL) :: rLength
      real(kind=AE_REAL) :: rCaissonHead, rModeledQ
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx, this, next

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_Inquiry: CW0_Create has not been called")
      end if

    write (unit=iLU, &
           fmt="(""#CW0, ID, HEAD_SPEC, XC, YC, Q, MODEL_Q, CAISSON_HEAD"")")
      do iWel = 1, cw0%iCount
        wel => cw0%Wells(iWel)
        ! Write the overall collector-well record CW0
        rad => wel%Radials(wel%iNRad)
        vtx => rad%Vertices(wel%iResolution)
        rModeledQ = rCW0_ComputeTotalFlow(io, cw0, iWel)
        rCaissonHead = rCW0_ComputeCaissonHead(io, cw0, aqu, iWel)
        write (unit=iLU, &
               fmt="(""CW0"", 1("", "", i9), 1("", "", l9), 5("", "", e16.8))" &
               ) wel%iID, wel%lHeadSpec, cIO_WorldCoords(io, wel%cZC), wel%rQ, rModeledQ, rCaissonHead
        ! Write records for all radials(CWR)
        write (unit=iLU, &
               fmt="(""#CWR, WELL_ID, RADIAL, VERTEX, X0, Y0, X1, Y1, LENGTH, STRENGTH, MODEL_HEAD, ERROR"")")
        do iRad = 1, wel%iNRad
          rad => wel%Radials(iRad)
          do iVtx = 1, wel%iResolution-1
            this => rad%Vertices(iVtx)
            next => rad%Vertices(iVtx+1)
            rLength = abs(next%cZ - this%cZ)
            write (unit=iLU, &
                   fmt="(""CWR"", 3("", "", i9), 8("", "", e17.10))") &
                   wel%iID, &
                   iRad, &
                   iVtx, &
                   cIO_WorldCoords(io, this%cZ), &
                   cIO_WorldCoords(io, next%cZ), &
                   rLength, &
                   this%rStrength, &
                   rAQU_PotentialToHead(io, aqu, this%rCheck, this%cCPZ), &
                   rAQU_PotentialToHead(io, aqu, this%rCheck - &
                   this%rResistanceTerm*this%rStrength, this%cCPZ)
          end do
        end do
      end do

      return
    end subroutine CW0_Inquiry


    subroutine CW0_Report(io, cw0, aqu)
      !! subroutine CW0_Report
      !!
      !! Writes a debugging report for all line-sinks to LU_OUTPUT
      !!
      !! Calling Sequence:
      !!    call CW0_Report(cw0)
      !!
      !! Arguments:
      !!   (in)    type(CW0_COLLECTION), pointer
      !!             CW0_COLLECTION object to be used
      !!
      ! [ ARGUMENTS ]
      type(CW0_COLLECTION), pointer :: cw0
      type(AQU_COLLECTION), pointer :: aqu
      type(IO_STATUS), pointer :: io
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iWel, iRad, iVtx
      integer(kind=AE_INT) :: nWL, nPD, nDP, nEQ, nUN
      real(kind=AE_REAL) :: rCaissonHead, rModeledQ
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      if (io%lDebug) then
        call IO_Assert(io, (associated(cw0)), &
             "CW0_Report: CW0_Create has not been called")
      end if

      call HTML_Header('Module CW0', 1)
      call HTML_Header('Discharge-specified area sink information', 2)

      if (associated(cw0%Wells)) then
        call HTML_StartTable()
        call HTML_AttrInteger('Number of wells', cw0%iCount)
        call HTML_AttrInteger('Number of FWL functions', iCW0_GetInfo(io, cw0, SIZE_FWL, 0))
        call HTML_AttrInteger('Number of FPD functions', iCW0_GetInfo(io, cw0, SIZE_FPD, 0))
        call HTML_AttrInteger('Number of FDP functions', iCW0_GetInfo(io, cw0, SIZE_FDP, 0))
        call HTML_AttrInteger('Number of equations', iCW0_GetInfo(io, cw0, SIZE_EQUATIONS, 0))
        call HTML_AttrInteger('Number of unknowns', iCW0_GetInfo(io, cw0, SIZE_UNKNOWNS, 0))
        call HTML_EndTable()

        ! Write the string information for all wells...
        do iWel = 1, cw0%iCount
          wel => cw0%Wells(iWel)
          rad => wel%Radials(wel%iNRad)
          vtx => rad%Vertices(wel%iResolution)
          rModeledQ = rCW0_ComputeTotalFlow(io, cw0, iWel)
          rCaissonHead = rCW0_ComputeCaissonHead(io, cw0, aqu, iWel)
          call HTML_Header('Well information', 3)
          call HTML_StartTable()
          call HTML_AttrInteger('Well #', iWel)
          call HTML_AttrInteger('ID', wel%iID)
          call HTML_AttrInteger('Number of radials', wel%iNRad)
          call HTML_AttrInteger('Resolution', wel%iResolution)
          call HTML_AttrReal('Center X', real(wel%cZC), '(f12.1)')
          call HTML_AttrReal('Center Y', aimag(wel%cZC), '(f12.1)')
          call HTML_AttrReal('Caisson radius', wel%rRc)
          call HTML_AttrReal('Pumping rate', wel%rQ)
          call HTML_AttrReal('Modeled pumping rate', rModeledQ)
          call HTML_AttrReal('Modeled caisson head', rCaissonHead)
          call HTML_EndTable()

          do iRad = 1, wel%iNRad
            rad => wel%Radials(iRad)
            call HTML_Header('Radial information', 4)
            call HTML_StartTable()
            call HTML_AttrInteger('Radial #', iRad)
            call HTML_AttrInteger('FWL index', rad%iFWLIndex)
            call HTML_AttrReal('Tip X', real(rad%cZEnd), '(f12.1)')
            call HTML_AttrReal('Tip Y', aimag(rad%cZEnd), '(f12.1)')
            call HTML_AttrReal('Resistance', rad%rResistance)
            call HTML_AttrReal('Radial width', rad%rWidth)
            call HTML_AttrReal('Blank length', rad%rBlankLength)
            call HTML_EndTable()

            call HTML_Header('Radial vertices', 4)
            call HTML_StartTable()
            call HTML_TableHeader((/'      ', 'FDP # ', 'X     ', 'Y     ', 'Sigma ', 'Length', 'Head  ','Inside'/))
            do iVtx = 1, wel%iResolution+1
              vtx => rad%Vertices(iVtx)
              call HTML_StartRow()
              call HTML_ColumnInteger((/ iVtx, vtx%iFDPIndex /))
              call HTML_ColumnComplex((/ vtx%cZ /))
              call HTML_ColumnReal((/ vtx%rStrength /))
              call HTML_ColumnReal((/ vtx%rLength /))
              call HTML_ColumnReal((/ rAQU_PotentialToHead(io, aqu, vtx%rCheck, vtx%cCPZ) /))
              call HTML_ColumnReal((/ rAQU_PotentialToHead(io, aqu, vtx%rCheck, vtx%cCPZ) - &
                                      vtx%rStrength*rad%rResistance/rad%rWidth /))
              call HTML_EndRow()
            end do
            call HTML_EndTable()
          end do
        end do
      else
        call HTML_Header('No collector wells defined', 3)
      end if

      return
    end subroutine CW0_Report


    subroutine CW0_Save(io, cw0, mode)
      !! Saves the current solution information onto the SCRATCH LU
      ! [ ARGUMENTS ]
      type(IO_STATUS), pointer :: io
      type(CW0_COLLECTION), pointer :: cw0
      integer(kind=AE_INT), intent(in) :: mode
      ! [ LOCALS ]
      integer(kind=AE_INT) :: iwel, irad, ivtx
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      ! Output records will be of the form ELEM_CW0, IWEL, IRAD, IVTX, 0, SIGMA
      do iwel = 1, cw0%iCount
        wel => cw0%Wells(iwel)
        do irad = 1, wel%iNRad
          rad => wel%Radials(irad)
          do ivtx = 1, wel%iResolution
            vtx => rad%Vertices(ivtx)
            if (mode == IO_MODE_BINARY) then
              write (unit=LU_SCRATCH) ELEM_CW0, iwel, irad, ivtx, vtx%rStrength
            else
              write (unit=LU_SCRATCH, fmt=*) "CW0", iwel, irad, ivtx, vtx%rStrength
            end if
          end do
        end do
      end do

      return
    end subroutine CW0_Save


    subroutine CW0_Load(io, cw0, fwl, fdp, mode)
      !! Loads the CW0 records from the file on the SCRATCH LU
      ! [ ARGUMENTS ]
      type(IO_STATUS), pointer :: io
      type(CW0_COLLECTION), pointer :: cw0
      type(FDP_COLLECTION), pointer :: fdp
      type(FWL_COLLECTION), pointer :: fwl
      integer(kind=AE_INT), intent(in) :: mode
      ! [ LOCALS ]
      integer(kind=AE_INT) :: imodule, iwel, irad, ivtx, istat
      real(kind=AE_REAL) :: rstrength
      character(len=3) :: smodule
      type(CW0_WELL), pointer :: wel
      type(CW0_RADIAL), pointer :: rad
      type(CW0_VERTEX), pointer :: vtx

      ! Scans the entire precondition file for the CW0 data
      rewind(unit=LU_SCRATCH)
      do
        if (mode == IO_MODE_BINARY) then
          read (unit=LU_SCRATCH, iostat=istat) imodule, iwel, irad, ivtx, rstrength
          if (imodule /= ELEM_CW0) cycle
        else
          read (unit=LU_SCRATCH, fmt=*, iostat=istat) smodule, iwel, irad, ivtx, rstrength
          if (trim(uppercase (smodule)) /= "CW0") cycle
        end if
        if (istat < 0) exit
        call IO_Assert(io, istat == 0, "I/O error on precondition file")
        call IO_Assert(io, iwel > 0 .and. iwel <= cw0%iCount, "CW0 well not found")
        wel => cw0%Wells(iwel)
        call IO_Assert(io, irad > 0 .and. irad <= wel%iNRad, "CW0 radial not found")
        rad => wel%Radials(irad)
        call IO_Assert(io, ivtx > 0 .and. irad <= wel%iResolution, "CW0 vertex not found")
        vtx => rad%Vertices(ivtx)
        vtx%rStrength = rstrength
      end do

      ! Now, populate the internal data structures
      call CW0_Update(io, cw0, fwl, fdp)

      return
    end subroutine CW0_Load

  end module m_cw0
