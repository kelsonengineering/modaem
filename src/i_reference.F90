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

!> Module of element customization for the 2-D reference flow field.
!>
!> This module contains only the influence functions, to complete the specification of the
!> ELEMENT class.
!>
!> Module use:
!>   u_constants  --  Universal ModAEM constant declarations
!>   u_element    --  Base element class specifications
!>
module i_reference

  use u_constants
  use u_element
  implicit none

  !> Contains a reference flow field. The reference flow field has three components:
  !>    1) The constant of integration
  !>    2) The far-field total discharge in the x direction (Q_x0)
  !>    3) The far-field total discharge in the y direction (Q_y0)
  type, public, extends(ELEMENT) :: IRF_REFERENCE
    complex(kind=AE_REAL) :: cZ0 !< Complex coordinate of the reference point.

  contains
    !> Derived-class initialization routine.
    procedure, pass :: Configure => IRF_Configure
    !> Derived-class generic creation routine.
    procedure, nopass :: New => IRF_New

    ! Implement the base-class methods
    !> Type binding for the unit potential
    procedure, pass :: UnitPotential => cIRF_UnitPotential
    !> Type binding for the unit discharge
    procedure, pass :: UnitDischarge => cIRF_UnitDischarge
    !> Type binding for the unit integrated flow
    procedure, pass :: UnitFlow => rIRF_UnitFlow
    !> Type binding for the unit recharge
    procedure, pass :: UnitRecharge => rIRF_UnitRecharge
    !> Type binding for the unit extraction
    procedure, pass :: UnitExtraction => rIRF_UnitExtraction

  end type IRF_REFERENCE

  public

contains

  !> Configures a new IRF_Well object. Typically, this is used to configure an entry in
  !> a collection in another package, e.g. an array of IRF_REFERENCE objects in module f_well.
  !> @param cZ0 Complex coordinate of the reference point.
  !> @param rConstant The initial reference field constant of integration.
  !> @param cQ0 The complex far-field total discharge vector.
  !> @returns Allocation status. Nonzero indicates an internal allocation failure.
  function IRF_Configure(el, cZ0, rConstant, cQ0) result(iStat)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(inout) :: el
    complex(kind=AE_REAL), intent(in) :: cZ0
    real(kind=AE_REAL), intent(in) :: rConstant
    complex(kind=AE_REAL), intent(in) :: cQ0
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iStat

    el%iOrder = 3
    allocate(el%rStrength(el%iOrder), stat=iStat)
    if (iStat == 0) then
      el%cZ0 = cZ0
      el%rStrength(1) = rConstant
      el%rStrength(2) = real(cQ0, AE_REAL)
      el%rStrength(3) = aimag(cQ0)
    end if

    return
  end function IRF_Configure


  !> Creates a new IRF_REFERENCE object and returns a standalone ELEMENT_REFERENCE object.
  !> @param cZ0 Complex coordinate of the reference point.
  !> @param rConstant The initial reference field constant of integration.
  !> @param cQ0 The complex far-field total discharge vector.
  !> @returns An ELEMENT_REFERENCE containing a pointer to the IRF_REFERENCE.
  function IRF_New(cZ0, rConstant, cQ0) result(eref)
    complex(kind=AE_REAL), intent(in) :: cZ0
    real(kind=AE_REAL), intent(in) :: rConstant
    complex(kind=AE_REAL), intent(in) :: cQ0
    ! [ RETURN VALUE ]
    type(ELEMENT_REFERENCE) :: eref
    ! [ LOCALS ]
    type(IRF_REFERENCE), pointer :: irf
    integer(kind=AE_INT) :: iStat

    allocate(irf, stat=iStat)
    if (iStat == 0) then
      if (irf%Configure(cZ0, rConstant, cQ0) == 0) then
        eref%el => irf
      else
        eref%el => NULL()
      end if
    else
      eref%el => NULL()
    endif

    return
  end function IRF_New


  !> Computes the unit complex potential for a well at the provided
  !> complex coordinate.
  !> @param el An IRF_REFERENCE object.
  !> @param cZ The complex coordinate in question.
  !> @returns A single-element array containing the unit potential for the well.
  pure function cIRF_UnitPotential(el, cZ, iOrder) result(cRV)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iOrder
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cRV(iOrder)
    ! [ LOCALS ]

    cRV(1) = rONE
    cRV(2) = -rONE
    cRV(3) = -rONE

    return
  end function cIRF_UnitPotential


  !> Computes the unit complex discharge for a well at the provided
  !> complex coordinate. The returned value is the complex conjugate of the
  !> discharge vector per unit strength of the well.
  !> @param el An IRF_REFERENCE object.
  !> @param cZ The complex coordinate in question.
  !> @returns A single-element array containing the unit discharge for the well.
  pure function cIRF_UnitDischarge(el, cZ, iOrder) result(cRV)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iOrder
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cRV(iOrder)
    ! [ LOCALS ]

    ! Note: the ...InfluenceW functions return the CONJUGATE of the discharge per unit strength
    cRV(1) = cZERO
    cRV(2) = rONE
    cRV(3) = -rONE

    return
  end function cIRF_UnitDischarge


  !> Computes the unit flow for a well over the provided line segment. The returned value
  !> is the complex conjugate of the integrated flow per unit strength of the well.
  !> @param el An IRF_REFERENCE object.
  !> @param cZ The complex coordinate in question.
  !> @returns A single-element array containing the unit integrated flow for the well over the
  !>    line segment CZ1-CZ2.
  pure function rIRF_UnitFlow(el, cZ1, cZ2, iOrder) result(rRV)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ1
    complex(kind=AE_REAL), intent(in) :: cZ2
    integer(kind=AE_INT), intent(in) :: iOrder
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRV(iOrder)
    ! [ LOCALS ]

    rRV(1) = rZERO
    rRV(2) = -aimag(cZ2 - cZ1)
    rRV(3) = real(cZ2 - cZ1)

    return
  end function rIRF_UnitFlow


  !> Computes the unit recharge for a well at the provided complex coordinate. For the
  !> reference field, this value is zero everywhere.
  !> @param el An IRF_REFERENCE object.
  !> @param cZ The complex coordinate in question.
  !> @returns A single-element array containing the unit recharge for the well.
  pure function rIRF_UnitRecharge(el, cZ, iOrder) result(rRV)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(in) :: iOrder
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRV(iOrder)
    ! [ LOCALS ]

    rRV = rZERO

    return
  end function rIRF_UnitRecharge


  !> Computes the unit extraction rate for a well. For the reference field, the value
  !> is infinite for all three coefficients; the reference point does not enter into
  !> the model water budget.
  !> @param el An IRF_REFERENCE object.
  !> @returns A single-element array containing the unit recharge for the well.
  pure function rIRF_UnitExtraction(el, iOrder) result(rRV)
    ! [ ARGUMENTS ]
    class(IRF_REFERENCE), intent(in) :: el
    integer(kind=AE_INT), intent(in) :: iOrder
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rRV(iOrder)
    ! [ LOCALS ]

    rRV = HUGE(AE_REAL)

    return
  end function rIRF_UnitExtraction


end module i_reference
