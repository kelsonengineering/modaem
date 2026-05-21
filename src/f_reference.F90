module f_reference

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

  !! module f_reference
  !!
  !! Function module for the reference flow field (singleton).
  !!
  !! Provides the complex potential, discharge, and integrated flow
  !! influence functions for:
  !!   - A constant of integration rSolConst (the matrix unknown)
  !!   - A prescribed uniform flow cRefUniformFlow = Qx + i*Qy
  !!   - An optional reference point cRefPoint at which head is specified
  !!
  !! Physics:
  !!   Omega_ref(z) = rSolConst - conjg(W_0) * z   (reference point specified)
  !!   Omega_ref(z) = rSolConst                     (no reference point)
  !!   W_ref(z)     = W_0                           (constant everywhere)
  !!
  !! Because only one reference flow field exists per model, FRF_COLLECTION
  !! holds a single record rather than an array.
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   u_io         --  Universal ModAEM I/O functions and constants

  use u_constants
  use u_io

  implicit none

  public


  type, public :: FRF_COLLECTION
    !! type FRF_COLLECTION
    !!
    !! Holds the singleton reference flow field for the model.
    !!
    !! Members:
    !!   complex :: cRefPoint
    !!     Location of the reference point (if specified)
    !!   real :: rRefHead
    !!     Specified head at the reference point (if specified)
    !!   complex :: cRefUniformFlow
    !!     Prescribed uniform flow vector W_0 = Qx + i*Qy
    !!   logical :: lReference
    !!     .true. if a reference point has been specified
    !!   real :: rSolConst
    !!     The constant of integration (the matrix unknown)
    !!   real :: rCheck
    !!     Check value of the potential at the reference point
    !!
    complex(kind=AE_REAL) :: cRefPoint
    real(kind=AE_REAL) :: rRefHead
    complex(kind=AE_REAL) :: cRefUniformFlow
    logical :: lReference
    real(kind=AE_REAL) :: rSolConst
    real(kind=AE_REAL) :: rCheck
  end type FRF_COLLECTION


contains


  function FRF_Create(io) result(frf)
    !! function FRF_Create
    !!
    !! Creates and initializes a new FRF_COLLECTION object.
    !!
    !! Calling Sequence:
    !!    frf => FRF_Create(io)
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    type(FRF_COLLECTION), pointer :: frf
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    allocate(frf, stat=iStat)
    call IO_Assert(io, (iStat == 0), "FRF_Create: allocation failed")

    frf%cRefPoint = cZERO
    frf%rRefHead = rZERO
    frf%cRefUniformFlow = cZERO
    frf%lReference = .false.
    frf%rSolConst = rZERO
    frf%rCheck = rZERO

    return
  end function FRF_Create


  subroutine FRF_SetReference(io, frf, cRefPoint, rRefHead, cRefUniformFlow)
    !! subroutine FRF_SetReference
    !!
    !! Stores the reference point location, head, and uniform flow vector.
    !!
    !! Calling Sequence:
    !!    call FRF_SetReference(io, frf, cRefPoint, rRefHead, cRefUniformFlow)
    !!
    !! Arguments:
    !!   (in)    complex :: cRefPoint
    !!             Complex coordinate of the reference point
    !!   (in)    real :: rRefHead
    !!             Specified head at the reference point
    !!   (in)    complex :: cRefUniformFlow
    !!             Prescribed uniform flow W_0 = Qx + i*Qy
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cRefPoint
    real(kind=AE_REAL), intent(in) :: rRefHead
    complex(kind=AE_REAL), intent(in) :: cRefUniformFlow
    type(IO_STATUS), pointer :: io

    call IO_Assert(io, associated(frf), "FRF_SetReference: No FRF_COLLECTION object")

    frf%cRefPoint = cRefPoint
    frf%rRefHead = rRefHead
    frf%cRefUniformFlow = cRefUniformFlow
    frf%lReference = .true.

    return
  end subroutine FRF_SetReference


  function cFRF_Potential(io, frf, cZ) result(cOmega)
    !! complex function cFRF_Potential
    !!
    !! Returns the reference flow field contribution to the complex potential at cZ.
    !!
    !!   Omega_ref(z) = rSolConst - conjg(W_0) * z   (reference point specified)
    !!   Omega_ref(z) = rSolConst                     (no reference point)
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cOmega

    call IO_Assert(io, associated(frf), "cFRF_Potential: No FRF_COLLECTION object")

    if (frf%lReference) then
      cOmega = frf%rSolConst - conjg(frf%cRefUniformFlow) * cZ
    else
      cOmega = frf%rSolConst
    end if

    return
  end function cFRF_Potential


  function cFRF_Discharge(io, frf, cZ) result(cQ)
    !! complex function cFRF_Discharge
    !!
    !! Returns the reference flow field contribution to the complex discharge at cZ.
    !!
    !!   W_ref(z) = W_0   (constant everywhere when reference point specified)
    !!   W_ref(z) = 0     (no reference point)
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), intent(in) :: cZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ

    call IO_Assert(io, associated(frf), "cFRF_Discharge: No FRF_COLLECTION object")

    if (frf%lReference) then
      cQ = frf%cRefUniformFlow
    else
      cQ = cZERO
    end if

    return
  end function cFRF_Discharge


  function rFRF_Flow(io, frf, cPathZ) result(rFlow)
    !! real function rFRF_Flow
    !!
    !! Returns the reference flow field contribution to the integrated flow
    !! along the path cPathZ. Since the reference field is analytic, only the
    !! endpoint potentials are needed.
    !!
    ! [ ARGUMENTS ]
    type(FRF_COLLECTION), pointer :: frf
    complex(kind=AE_REAL), dimension(:), intent(in) :: cPathZ
    type(IO_STATUS), pointer :: io
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rFlow

    call IO_Assert(io, associated(frf), "rFRF_Flow: No FRF_COLLECTION object")

    if (frf%lReference) then
      rFlow = aimag(cFRF_Potential(io, frf, cPathZ(size(cPathZ))) - &
                    cFRF_Potential(io, frf, cPathZ(1)))
    else
      rFlow = rZERO
    end if

    return
  end function rFRF_Flow


end module f_reference
