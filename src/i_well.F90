module i_well

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


  !! module i_well
  !!
  !! Module of influence function for 2-D wells
  !!
  !! This module contains only the influence functions used in
  !! modules f_well and others.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!
  use u_constants

  implicit none

  ! Public functions -- called by other ModAEM routines
  private

  public :: cIWL_InfluenceP, &
           cIWL_InfluenceW, &
           cIWL_InfluenceF, &
           cIWL_InfluenceG, &
           cIWL_InfluenceQ


contains

  pure function cIWL_InfluenceP(cZ, cZC) result(cP)
    !! complex function cIWL_InfluenceP
    !!
    !! Computes the complex potential influence function for a single
    !! well with center cZC at the complex coordinate cMapZ. If the well
    !! function is used to represent a well element, it is up to the
    !! caller to move the point to the perimeter of the well,
    !!
    !! For the well, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cP = cIWL_InfluenceP(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             Returns NaN for cZ = cZC
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the well
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
    complex(kind=AE_REAL) :: cP(1, 1)
    ! [ LOCALS ]

    cP(1, 1) = log(cZ-cZC)/(rTWO * rPI)

    return
  end function cIWL_InfluenceP

  pure function cIWL_InfluenceW(cZ, cZC) result(cW)
    !! complex function cIWL_InfluenceQ
    !!
    !! Computes the complex discharge influence function for a single
    !! well with center cZC at the complex coordinate cMapZ. If the well
    !! function is used to represent a well element, it is up to the
    !! caller to move the point to the perimeter of the well,
    !!
    !! For the well, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cW = cIWL_InfluenceW(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             Returns -Inf for cZ = cZC
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the well
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

    ! Note: the ...InfluenceW functions return the CONJUGATE of the discharge per unit strength
    cW(1, 1) = -rONE / ((cZ-cZC)*(rTWO * rPI))

    return
  end function cIWL_InfluenceW

  pure function cIWL_InfluenceF(cZ1, cZ2, cZC) result(cF)
    !! complex function cIWL_InfluenceF
    !!
    !! Computes the integrated flow influence function for a single
    !! well with center cZC between the points cZ1 and cZ2. If the well
    !! function is used to represent a well element, it is up to the
    !! caller to move point(s) to the perimeter of the well,
    !!
    !! For the well, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cF = cIWL_InfluenceF(cZ1, cZ2, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinates of the ends of the line
    !!             segment in question. Returns -NaN for cZ1 = cZC or
    !!             cZ2 = cZC.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the well
    !!
    !! Returns:
    !!           complex :: cQ
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2
    complex(kind=AE_REAL), intent(in) :: cZC
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cF(1, 1)
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cP1(1, 1), cP2(1, 1)
    real(kind=AE_REAL) :: rX1, rX2, rY1, rY2, rXi, rBC(1, 1)

    ! Get the complex potentials first -- they generate the Psi values.
    cP1 = cIWL_InfluenceP(cZ1, cZC)
    cP2 = cIWL_InfluenceP(cZ2, cZC)
    ! Now, look for a branch cut!
    rBC = rZERO
    rX1 = real(cZ1-cZC)
    rY1 = aimag(cZ1-cZC)
    rX2 = real(cZ2-cZC)
    rY2 = aimag(cZ2-cZC)
    ! There MIGHT be a branch cut...
    if ((rY1 >= rZERO .and. rY2 < rZERO) .or. &
        (rY1 < rZERO .and. rY2 >= rZERO)) then
      rXi = rX1 - rY1*((rX2-rX1)/(rY2-rY1))
      if (rXi < rZERO) then
        ! BRANCH CUT!  Check for sign on branch cut...
        if (rY1 < 0) then
          rBC = -rONE
        else
          rBC = rONE
        end if
      end if
    end if

    ! Compute the influence function as the difference in Psi plus the branch cut
    cF = aimag(cP2-cP1) + rBC

    return
  end function cIWL_InfluenceF

  pure function cIWL_InfluenceG(cZ, cZC) result(cG)
    !! complex function rIWL_InfluenceG
    !!
    !! Computes the recharge influence function for a single
    !! well with center cZC at the complex coordinate cMapZ. If the well
    !! function is used to represent a well element, it is up to the
    !! caller to move the point to the perimeter of the well,
    !!
    !! For the well, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rG = 0 for all Z
    !!
    !! Calling Sequence:
    !!    rG = rIWL_InfluenceG(cZ, cZC)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the well
    !!
    !! Returns:
    !!           complex :: cQ
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cG(1, 1)
    ! [ LOCALS ]

    cG = rZERO

    return
  end function cIWL_InfluenceG

  pure function cIWL_InfluenceQ() result(cQ)
    !! complex function cIWL_InfluenceQ
    !!
    !! Computes the extraction rate influence function for a single
    !! well.
    !!
    !! For the well, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rQ = 1
    !!
    !! Calling Sequence:
    !!    cQ = cIWL_InfluenceQ()
    !!
    !! Arguments:
    !!           (none)
    !!
    !! Returns:
    !!           complex :: cQ
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ(1, 1)
    ! [ LOCALS ]

    cQ = rONE

    return
  end function cIWL_InfluenceQ

end module i_well

