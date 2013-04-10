!> module i_areasink
!!
!! Module of influence function for 1-D unbounded unbounded area-sinks
!!
!! This module contains only the influence functions used in
!! modules m_as0 and others.
!!
module i_areasink

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


  !! module i_areasink
  !!
  !! Module of influence function for 1-D unbounded unbounded area-sinks
  !!
  !! This module contains only the influence functions used in
  !! modules m_as0 and others.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!
  use u_constants

  implicit none

  ! Public functions -- called by other ModAEM routines
  private

  public :: rIAS_InfluenceP, &
           cIAS_InfluenceW, &
           rIAS_InfluenceF, &
           rIAS_InfluenceG

  ! NOTE: This module does not contain the "InfluenceQ" function -- it
  ! is only relevant for closed polygons, and thus must be handled by the
  ! m_as* modules explicitly.


contains

  pure function rIAS_InfluenceP(cZ, cZC) result(rP)
    !! complex function rIAS_InfluenceP
    !!
    !! Computes the complex potential influence function for a single
    !! unbounded area-sink with origin cZC at the complex coordinate cZ.
    !!
    !! For the unbounded area-sink, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cP = rIAS_InfluenceP(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             Returns NaN for cZ = cZC
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the origin of the unbounded area-sink
    !!
    !! Returns:
    !!           complex :: cP
    !!             Influence function for complex potential at the point
    !!             cZ. Both the real and imaginary parts exist for all
    !!             cZ /= cZC
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rP(1, 1)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rDX

    rDX = real(cZ-cZC)
    rP(1, 1) = rHALF * rDX*rDX

    return
  end function rIAS_InfluenceP

  pure function cIAS_InfluenceW(cZ, cZC) result(cW)
    !! complex function rIAS_InfluenceQ
    !!
    !! Computes the complex discharge influence function for a single
    !! unbounded area-sink with origin cZC at the complex coordinate cMapZ. If the unbounded area-sink
    !! function is used to represent a unbounded area-sink element, it is up to the
    !! caller to move the point to the perimeter of the unbounded area-sink,
    !!
    !! For the unbounded area-sink, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cW = rIAS_InfluenceW(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             Returns -Inf for cZ = cZC
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the origin of the unbounded area-sink
    !!
    !! Returns:
    !!           complex :: cW
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cW(1, 1)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rDX

    ! Note: the ...InfluenceW routines return the CONJUGATE of the discharge per unit strength
    rDX = real(cZ-cZC)
    cW(1, 1) = -conjg(cmplx(rDX, rZERO, AE_REAL))

    return
  end function cIAS_InfluenceW

  pure function rIAS_InfluenceF(cZ1, cZ2, cZC) result(rF)
    !! complex function rIAS_InfluenceF
    !!
    !! Computes the integrated flow influence function for a single
    !! unbounded area-sink with origin cZC between the points cZ1 and cZ2.
    !!
    !! For the unbounded area-sink, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cF = rIAS_InfluenceF(cZ1, cZ2, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinates of the ends of the line
    !!             segment in question. Returns -NaN for cZ1 = cZC or
    !!             cZ2 = cZC.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the origin of the unbounded area-sink
    !!
    !! Returns:
    !!           complex :: cF
    !!             Influence function for total flow across the segment cZ1-cZ2
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), intent(in) :: cZC
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rF(1, 1)
    ! [ LOCALS ]

    rF(1, 1) = - aimag(cZ2-cZ1) * rHALF * real(cZ2+cZ1-2*cZC)

    return
  end function rIAS_InfluenceF

  pure function rIAS_InfluenceG(cZ, cZC) result(rG)
    !! complex function rIAS_InfluenceG
    !!
    !! Computes the recharge influence function for a single
    !! unbounded area-sink with origin cZC at the complex coordinate cMapZ.
    !!
    !! For the unbounded area-sink, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state unbounded area-sink, rG = 1 for all Z
    !!
    !! Calling Sequence:
    !!    rG = rIAS_InfluenceG(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the origin of the unbounded area-sink
    !!
    !! Returns:
    !!           complex :: cQ
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rG(1, 1)
    ! [ LOCALS ]

    rG = rONE

    return
  end function rIAS_InfluenceG

end module i_areasink

