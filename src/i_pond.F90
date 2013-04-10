module i_pond

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

  !! module i_pond
  !!
  !! Module of influence functions for first-order
  !! circular recharge areas("ponds").
  !!
  !! This module contains only the influence functions used in
  !! modules f_pond and others.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!
  use u_constants

  implicit none

  public
  
  ! Private functions
  private :: rIPD_InfFlowInside, &
             rIPD_InfFlowOutside


contains

  pure function cIPD_InfluenceP(cZ, cZC, rRad) result(cP)
    !! complex function cIPD_InfluenceP
    !!
    !! Computes the complex potential influence function for a single
    !! pond with center cZC and radius rRad at the complex coordinate cMapZ.
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cP = cIPD_InfluenceP(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           complex :: cP
    !!             Influence function for complex potential at the point
    !!             cZ.  Note that for cZ inside the pond, no streamfunction
    !!             exists; aimag(cP) is zero. The streamfunction is
    !!             calculated for cZ outside the pond.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cP
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rR

    rR = abs(cZ-cZC)
    if (rR <= rRad) then
      ! Inside function
      cP(1, 1) = cmplx(rONE_FOURTH * (rR*rR - rRad*rRad), rZERO, AE_REAL)
    else
      ! Outside function -- treat as a well
      cP(1, 1) = rHALF * rRad*rRad * log((cZ-cZC)/rRad)
    end if

    return
  end function cIPD_InfluenceP

  pure function cIPD_InfluenceW(cZ, cZC, rRad) result(cW)
    !! complex function cIPD_InfluenceW
    !!
    !! Computes the complex discharge influence functions for a single
    !! pond with center cZC and radius rRad at the coordinate cZ.
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    cW = cIPD_InfluenceW(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ
    !!             The complex coordinate of the point in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           complex :: cW
    !!             Influence function for complex discharge at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cW

    ! Note: the ...InfluenceW functions return the CONJUGATE of the discharge per unit strength
    if (abs(cZ-cZC) <= rRad) then
      ! Inside function
      cW(1, 1) = -rHALF * conjg(cZ-cZC)
    else
      ! Outside function -- treat as a well
      cW(1, 1) = -rHALF * rRad*rRad / (cZ-cZC)
    end if

    return
  end function cIPD_InfluenceW

  pure function cIPD_InfluenceF(cZ1, cZ2, cZC, rRad) result(cF)
    !! real function cIPD_InfluenceF
    !!
    !! Computes the integrated flow influence function for a single
    !! pond with center cZC and radius rRad between the complex coordinates
    !! cZ1 and cZ2.
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    rF = cIPD_InfluenceP(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinate of the points in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           real :: rS
    !!             Influence function for integrated flow between cZ1 and
    !!             cZ2.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2, cZC
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cF
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cZ, cZ0, cZi, cR1, cR2, cV12, cV21, cZInt1, &
                            cZInt2, cZCZ0, cZChord
    real(kind=AE_REAL) :: rR1, rR2, rA, rB, rC, rD, rZ0ZChord, rDZ

    ! The segment cZ1-cZ2 has five possiblities:
    !   1) Both points outside, no intersection with pond
    !   2) Both points outside, two intersections with pond
    !   3) First point inside, second outside, one intersection
    !   4) First point outside, second inside, one intersection
    !   5) Both points inside, no intersections with pond
    ! First, one must figure out which case it is.

    ! Disclaimer: this section of the code is very inelegant, but
    ! it works.  It was extracted from the research code SLWL without
    ! modification -- this part of the code can use a good working
    ! over.

    ! Get the radii for the two points
    rR1 = abs(cZ1-cZC)
    rR2 = abs(cZ2-cZC)
    ! Both points inside?
    if (rR1 <= rRad .and. rR2 <= rRad) then
      ! From in- > in with no intersections -- use inside function
      cF = rIPD_InfFlowInside(cZ1, cZ2, cZC)
    else if (rR1 < rRad .and. rR2 > rRad) then
      ! From in- > out with one intersection
      cR1 = cZ1-cZC
      cV12 = cZ2-cZ1
      rA = abs((real(cR1)*real(cV12)+aimag(cR1)*aimag(cV12)) / abs(cV12))
      rB = abs(cR1)
      rC = sqrt(rB**2-rA**2)
      rD = sqrt(rRad**2-rC**2)
      cZInt1 = (rD-rA)*cV12/abs(cV12)+cZ1
      ! Use inside function from Z1- > ZInt1 and Outside Function from ZInt1- > Z2
      cF = rIPD_InfFlowInside(cZ1, cZInt1, cZC) + &
           rIPD_InfFlowOutside(cZInt1, cZ2, cZC, rRad)
    else if (rR2 < rRad .and. rR1 > rRad) then
      ! From out- > in with one intersection
      cR2 = cZ2-cZC
      cV21 = cZ1-cZ2
      rA = abs((real(cR2)*real(cV21)+aimag(cR2)*aimag(cV21)) / abs(cV21))
      rB = abs(cR2)
      rC = sqrt(rB**2-rA**2)
      rD = sqrt(rRad**2-rC**2)
      cZInt1 = (rD-rA)*cV21 / abs(cV21)+cZ2
      ! Use outside function from Z1- > ZInt1 and inside Function from ZInt1- > Z2
      cF = rIPD_InfFlowOutside(cZ1, cZInt1, cZC, rRad) + &
           rIPD_InfFlowInside(cZInt1, cZ2, cZC)
    else
      ! Both points outside
      ! Check to see if there is a possibility that there COULD be an intersection
      ! Compute the directed distance between the center of the segment Z1-Z2 and
      ! the center of the circle. The imaginary part of this quantity is the
      ! distance from the center to the segment Z1-Z2. If this distance is greater
      ! than the circle radius, no intersections. Otherwise, compute the
      ! location of the intersections.
      cZi = (cZ2-cZ1) / abs(cZ2-cZ1)        ! Unit vector from Z1-Z2
      cZ0 = rHALF*(cZ2+cZ1)                 ! Center of segment Z1-Z2
      cZ = (cZC - rHALF*(cZ2+cZ1))/cZi
      if (abs(aimag(cZ)) >= rRad .or. abs(real(cZ)) > rONE) then
        ! From out- > out with no intersections -- Use Outside function
        cF = rIPD_InfFlowOutside(cZ1, cZ2, cZC, rRad)
      else
        ! From out- > out with two intersections
        cZCZ0 = cZ0-cZC                        ! Vector from Zc to Z0
        ! The distance from the center of the segment Z1-Z2 to the center of the
        ! chord of the circle subtended by Z1-Z2 is determined by the pythagorean
        ! theorem
        rZ0ZChord = sqrt(abs(cZCZ0)**2 - aimag(cZ)**2)
        ! Compute the location of the chord center by moving along the unit vector
        ! Z1-Z2 from the center of Z1-Z2. The direction is determined by the
        ! relative distances Zc-Z1 and Zc-Z2.
        if (abs(cZC-cZ1) <= abs(cZC-cZ1)) then
          cZChord = cZ0 - rZ0ZChord*cZi
        else
          cZChord = cZ0 + rZ0ZChord*cZi
        end if
        ! We are almost there! The distance from the center of the chord to the
        ! intersections is now determined by the pythagorean theorem
        rDZ = sqrt(rRad*rRad - aimag(cZ)**2)
        ! Move both directions from the chord center to the intersections
        cZInt1 = cZChord - rDZ*cZi
        cZInt2 = cZChord + rDZ*cZi
        ! Use Outside function from Z1- > ZInt1, Inside function from ZInt1- > ZInt2,
        ! Outside Function from ZInt2- > Z2.  PHEW!
        cF = rIPD_InfFlowOutside(cZ1, cZInt1, cZC, rRad) + &
             rIPD_InfFlowInside(cZInt1, cZInt2, cZC) + &
             rIPD_InfFlowOutside(cZInt2, cZ2, cZC, rRad)
      end if
    end if

    return
  end function cIPD_InfluenceF

  pure function rIPD_InfFlowInside(cZ1, cZ2, cZC) result(rF)
    !! real function cIPD_InfFlowInside
    !!
    !! Computes the integrated flow influence function for a single
    !! pond with center cZC and radius rRad between the complex coordinates
    !! cZ1 and cZ2, where both cZ1 and cZ2 are INSIDE the pond.
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    rF = cIPD_InfFlowInside(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinates of the points in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           real :: rS
    !!             Influence function for integrated flow between cZ1 and
    !!             cZ2.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2, cZC
    ! [ RETURN VALUE ]
    real(kind=AE_REAL), dimension(1, 1) :: rF
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cR, cS
    real(kind=AE_REAL) :: rDR1, rDR2, rDS1, rDS2, rDW3

    ! Vector cR is the vector from cZ1 to the center of the recharge circle
    cR = cZC-cZ1                        ! To center of circle
    rDR1 = real(cR)
    rDR2 = aimag(cR)
    ! Vector cS is the vector from cZ1 to cZ2
    cS = cZ2-cZ1
    rDS1 = real(cS)
    rDS2 = aimag(cS)
    ! rDW3 is the z-component of the curl of R and S
    rDW3 = rDR1*rDS2 - rDR2*rDS1
    ! The signed area bounded by R and S is 1/2 of W3
    rF(1, 1) = -rHALF*rDW3

    return
  end function rIPD_InfFlowInside

  pure function rIPD_InfFlowOutside(cZ1, cZ2, cZC, rRad) result(rF)
    !! real function cIPD_InfFlowInside
    !!
    !! Computes the integrated flow influence function for a single
    !! pond with center cZC and radius rRad between the complex coordinates
    !! cZ1 and cZ2, where both cZ1 and cZ2 are OUTSIDE the pond.
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    rS = cIPD_InfFlowOutside(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinates of the points in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           real :: rS
    !!             Influence function for integrated flow between cZ1 and
    !!             cZ2.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ1, cZ2, cZC
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    real(kind=AE_REAL), dimension(1, 1) :: rF
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cP1, cP2
    real(kind=AE_REAL) :: rX1, rX2, rY1, rY2, rXi, rBC

    ! Get the complex potentials first -- they generate the Psi values.
    cP1 = rHALF * rRad*rRad * log((cZ1-cZC)/rRad)
    cP2 = rHALF * rRad*rRad * log((cZ2-cZC)/rRad)
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
          rBC = -rPI * rRad*rRad
        else
          rBC = rPI * rRad*rRad
        end if
      end if
    end if

    ! The influence function is the difference in Psi plus the branch cut
    rF(1, 1) = aimag(cP2-cP1) + rBC

    return
  end function rIPD_InfFlowOutside

  pure function cIPD_InfluenceG(cZ, cZC, rRad) result(cG)
    !! real function cIPD_InfFlowInside
    !!
    !! Computes the recharge influence function for a single
    !! pond with center cZC and radius rRad at the coordinate cZ
    !!
    !! For the pond, one influence function is computed.
    !!
    !! Calling Sequence:
    !!    rS = cIPD_InfluenceG(cZ, cZC, rRad)
    !!
    !! Arguments:
    !!   (in)    complex :: cZ1, cZ2
    !!             The complex coordinates of the points in question.
    !!             No singularities exist got the pond function.
    !!   (in)    complex :: cZC
    !!             The complex coordinate of the center of the pond
    !!   (in)    real :: rRad
    !!             The radius of the pond
    !!
    !! Returns:
    !!           real :: rG
    !!             Influence function for recharge at cZ
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZ, cZC
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cG
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rR

    rR = abs(cZ-cZC)
    if (rR <= rRad) then
      cG = rONE
    else
      cG = rZERO
    end if

    return
  end function cIPD_InfluenceG

  pure function cIPD_InfluenceQ(rRad) result(cQ)
    !! complex function rIPD_InfluenceQ
    !!
    !! Computes the extraction rate influence function for a single
    !! well.
    !!
    !! For the well, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rQ = 1
    !!
    !! Calling Sequence:
    !!    rQ = rIPD_InfluenceQ()
    !!
    !! Arguments:
    !!    (in)   real :: rRad
    !!             Radius of the well
    !!
    !! Returns:
    !!           complex :: rQ
    !!             Influence function for extraction rate at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    real(kind=AE_REAL), intent(in) :: rRad
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ(1, 1)
    ! [ LOCALS ]

    cQ = rPI * rRad * rRad

    return
  end function cIPD_InfluenceQ

end module i_pond

