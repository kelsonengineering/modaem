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

!> Module of declarations for generic element classes. Makes it possible for generic
!> use of element types in ModAEM-2 packages and features.
module u_element

  use u_constants
  implicit none

  public

  !> Abstract base type for 2-D analytic elements. All element types have the same
  !> configuration: one or more strength coefficients. There are 7 bound "unit" functions
  !> for the complex potential, discharge, and other features. Each "unit" function
  !> returns a vector of coefficients, either real or complex; each entry corresponds
  !> to an entry in the strength-coefficient vector for the element. For example, a
  !> well has a single strength and a single entry in each unit function, while a
  !> line-dipole has three of each.
  !>
  !> In addition, each element has a set of corresponding bound "contribution" functions
  !> that return the dot product of the strength coefficient vector and unit functions.
  type, abstract :: ELEMENT
    !> A vector of strength coefficients
    integer(kind=AE_INT) :: iOrder
    real(kind=AE_REAL), dimension(:), allocatable :: rStrength

  contains

    !> Unit complex potential.
    procedure(Unit1CC), pass, deferred :: UnitPotential
    !> Unit complex discharge.
    procedure(Unit1CC), pass, deferred :: UnitDischarge
    !> Unit integrated discharge.
    procedure(Unit2CR), pass, deferred :: UnitFlow
    !> Unit recharge at a point.
    procedure(Unit1CR), pass, deferred :: UnitRecharge
    !> Unit extraction rate.
    procedure(Unit0R), pass, deferred :: UnitExtraction

    !> Complex potential contribution for the element.
    procedure, pass :: Potential => ElementPotential
    !> Complex discharge contribution for the element.
    procedure, pass :: Discharge => ElementDischarge
    !> Integrated discharge contribution for the element.
    procedure, pass :: Flow => ElementFlow
    !> Recharge at a point contribution for the element.
    procedure, pass :: Recharge => ElementRecharge
    !> Total extraction rate contribution for the element.
    procedure, pass :: Extraction => ElementExtraction

  end type ELEMENT


  !> Base type for line elements, such as line-sink and line-dipole elements. Adds
  !> functionality for looking up values along the line segment.
  type, abstract, extends(ELEMENT) :: LINE_ELEMENT

  contains

    !> Unit potential jump at a position along a line element. Override this function
    !> in types that require it.
    procedure(LE_Unit1RR), pass, deferred :: UnitPotentialJump
    !> Unit streamfunction jump at a position along a line element.  Override this function
    !> in types that require it.
    procedure(LE_Unit1RR), pass, deferred :: UnitStreamfunctionJump

    !> Discharge potential jump at a position along a line element.
    procedure, pass :: PotentialJump => LineElementPotentialJump
    !> Unit streamfunction jump at a position along a line element.
    procedure, pass :: StreamfunctionJump => LineElementStreamfunctionJump

  end type LINE_ELEMENT


  !> An object that contains a pointer to an ELEMENT object. This makes it possible
  !> to make arrays of generic ELEMENT objects, facilitating parallelism and
  !> simplifying application code.
  type, public :: ELEMENT_REFERENCE
    class(ELEMENT), pointer :: el
  end type ELEMENT_REFERENCE


  ! Abstract interfaces for base element unit functions.
  abstract interface
    !> An element unit function that takes no arguments, save for the element, and
    !> returns a vector of real unit functions.
    !> @param el The element instance.
    !> @returns A vector of real results.
    pure function Unit0R(el, iOrder) result(rRV)
      import :: ELEMENT
      import :: AE_REAL, AE_INT
      class(ELEMENT), intent(in) :: el
      integer(kind=AE_INT), intent(in) :: iOrder
      real(kind=AE_REAL), dimension(iOrder) :: rRV
    end function Unit0R

    !> An element unit bound function that takes a single complex coordinate and
    !> returns a vector of real unit functions.
    !> @param el The element instance.
    !> @param cZ A complex coordinate.
    !> @returns A vector of real results.
    pure function Unit1CR(el, cZ, iOrder) result(rRV)
      import :: ELEMENT
      import :: AE_REAL, AE_INT
      class(ELEMENT), intent(in) :: el
      complex(kind=AE_REAL), intent(in) :: cZ
      integer(kind=AE_INT), intent(in) :: iOrder
      real(kind=AE_REAL), dimension(iOrder) :: rRV
    end function Unit1CR

    !> An element unit bound function that takes a single complex coordinate and
    !> returns a vector of complex unit functions.
    !> @param el The element instance.
    !> @param cZ A complex coordinate.
    !> @returns A vector of real results.
    pure function Unit1CC(el, cZ, iOrder) result(cRV)
      import :: ELEMENT
      import :: AE_REAL, AE_INT
      class(ELEMENT), intent(in) :: el
      complex(kind=AE_REAL), intent(in) :: cZ
      integer(kind=AE_INT), intent(in) :: iOrder
      complex(kind=AE_REAL), dimension(iOrder) :: cRV
    end function Unit1CC

    !> An element unit bound function that takes a pair of complex coordinates and
    !> returns a vector of real unit functions.
    !> @param el The element instance.
    !> @param cZ1 The beginning complex coordinate.
    !> @param cZ2 The ending complex coordinate.
    !> @returns A vector of real results.
    pure function Unit2CR(el, cZ1, cZ2, iOrder) result(rRV)
      import :: ELEMENT
      import :: AE_REAL, AE_INT
      class(ELEMENT), intent(in) :: el
      complex(kind=AE_REAL), intent(in) :: cZ1
      complex(kind=AE_REAL), intent(in) :: cZ2
      integer(kind=AE_INT), intent(in) :: iOrder
      real(kind=AE_REAL), dimension(iOrder) :: rRV
    end function Unit2CR
  end interface

  !> Abstract interfaces for line element unit functions.
  abstract interface
    !> An element unit bound function that takes a single real position (along an element)
    !> and returns a vector of real unit functions.
    !> @param el The element instance.
    !> @param cZ A complex coordinate.
    !> @returns A vector of  results.
    pure function LE_Unit1RR(el, rX, iOrder) result(rRV)
      import :: LINE_ELEMENT
      import :: AE_REAL, AE_INT
      class(LINE_ELEMENT), intent(in) :: el
      real(kind=AE_REAL), intent(in) :: rX
      integer(kind=AE_INT), intent(in) :: iOrder
      real(kind=AE_REAL), dimension(:), allocatable :: rRV
    end function LE_Unit1RR

  end interface

contains

  !> Implements the Potential bound method.
  !> @param el An ELEMENT object
  !> @param cZ A complex coordinate point.
  !> @returns The complex potential contribution of the element at cZ.
  pure function ElementPotential(el, cZ) result(cRV)
    class(ELEMENT), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    complex(kind=AE_REAL) :: cRV

    cRV = dot_product(el%UnitPotential(cZ, el%iOrder), el%rStrength)

    return
  end function ElementPotential

  !> Implements the Discharge bound method.
  !> @param el An ELEMENT object
  !> @param cZ A complex coordinate point.
  !> @returns The complex discharge contribution of the element at cZ.
  pure function ElementDischarge(el, cZ) result(cRV)
    class(ELEMENT), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    complex(kind=AE_REAL) :: cRV

    cRV = dot_product(el%UnitDischarge(cZ, el%iOrder), el%rStrength)

    return
  end function ElementDischarge

  !> Implements the Flow bound method.
  !> @param el An ELEMENT object
  !> @param cZ1 Complex coordinate of the first end of a line segment.
  !> @param cZ2 Complex coordinate of the second end of a line segment.
  !> @returns The integrated discharge contribution of the element at cZ.
  pure function ElementFlow(el, cZ1, cZ2) result(rRV)
    class(ELEMENT), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ1
    complex(kind=AE_REAL), intent(in) :: cZ2
    real(kind=AE_REAL) :: rRV

    rRV = dot_product(el%UnitFlow(cZ1, cZ2, el%iOrder), el%rStrength)

    return
  end function ElementFlow

  !> Implements the Recharge bound method.
  !> @param el An ELEMENT object
  !> @param cZ A complex coordinate point.
  !> @returns The recharge contribution of the element at cZ.
  pure function ElementRecharge(el, cZ) result(rRV)
    class(ELEMENT), intent(in) :: el
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL) :: rRV

    rRV = dot_product(el%UnitRecharge(cZ, el%iOrder), el%rStrength)

    return
  end function ElementRecharge

  !> Implements the Extraction bound method.
  !> @param el An ELEMENT object
  !> @returns The extraction rate of the element.
  pure function ElementExtraction(el) result(rRV)
    class(ELEMENT), intent(in) :: el
    real(kind=AE_REAL) :: rRV

    rRV = dot_product(el%UnitExtraction(el%iOrder), el%rStrength)

    return
  end function ElementExtraction

  !> Implements the potential jump bound method for line elements.
  !> @param el An ELEMENT object
  !> @param rX A fractional position (range 0-1) along the element.
  !> @returns The discharge potential jump along the element at position rX.
  pure function LineElementPotentialJump(el, rX) result(rRV)
    class(LINE_ELEMENT), intent(in) :: el
    real(kind=AE_REAL), intent(in) :: rX
    real(kind=AE_REAL) :: rRV

    rRV = real(dot_product(el%UnitPotentialJump(rX, el%iOrder), el%rStrength), AE_REAL)

    return
  end function LineElementPotentialJump

  !> Implements the streamfunction jump bound method for line elements.
  !> @param el An ELEMENT object
  !> @param rX A fractional position (range 0-1) along the element.
  !> @returns The streamfunction jump along the element at position rX.
  pure function LineElementStreamfunctionJump(el, rX) result(rRV)
    class(LINE_ELEMENT), intent(in) :: el
    real(kind=AE_REAL), intent(in) :: rX
    real(kind=AE_REAL) :: rRV

    rRV = real(dot_product(el%UnitPotentialJump(rX, el%iOrder), el%rStrength), AE_REAL)

    return
  end function LineElementStreamfunctionJump

end module u_element
