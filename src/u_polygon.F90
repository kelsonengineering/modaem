module u_polygon

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

  !! module u_poly
  !!
  !! Utility module for polygons
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!
  !! This module contains routines computations based on polygonal
  !! areas in 2-D space.

  use u_constants

  implicit none

  private

  public :: PGN_Perimeter, &
           PGN_SignedArea, &
           PGN_Area, &
           PGN_IsPositivelyOriented, &
           PGN_MakePositivelyOriented, &
           PGN_CheckPoint, &
           PGN_Contains, &
           PGN_Segment, &
           PGN_ChopSegment, &
           PGN_LeftmostIndex, &
           PGN_OrientLeftmost, &
           PGN_ClipArea, &
           PGN_Build

  real(kind=AE_REAL) :: rPGNTempArea


contains


  function PGN_Perimeter(poly) result(rPerim)
    ! Returns the perimeter of the polygon
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rPerim
    ! [ LOCALS ]
    integer(kind=AE_INT) :: current, previous

    if (size(poly) < 3) then
      rPerim = HUGE(AE_REAL)
    else
      rPerim = rZERO
      previous = size(poly)
      do current = 1, size(poly)
        rPerim = rPerim + abs(poly(current)-poly(previous))
        previous = current
      end do
    end if

    return
  end function PGN_Perimeter


  function PGN_SignedArea(poly) result(rArea)
    ! Returns the(signed) area of the poly
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rArea
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: localPoly
    real(kind=AE_REAL) :: rSum
    integer(kind=AE_INT) :: i, n

    n = size(poly)
    if (n < 3) then
      rArea = HUGE(AE_REAL)
    else
      allocate(localPoly(0:n+1))
      localPoly(0:n-1) = poly
      localPoly(n:n+1) = localPoly(0:1)
      rSum = rZERO
      do i = 1, n
        rSum = rSum + real(localPoly(i)) * aimag(localPoly(i+1)-localPoly(i-1))
      end do
      rArea = rHALF * rSum
      deallocate(localPoly)
    end if

    return
  end function PGN_SignedArea


  function PGN_Area(poly) result(rArea)
    ! Returns the(signed) area of the poly
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rArea

    rArea = abs(PGN_SignedArea(poly))

    return
  end function PGN_Area


  function PGN_IsPositivelyOriented(poly) result(lPositive)
    ! Returns .true. if the polygon is positively-oriented
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    ! [ RETURN VALUE ]
    logical :: lPositive

    lPositive = (PGN_SignedArea(poly) > rZERO)

    return
  end function PGN_IsPositivelyOriented


  subroutine PGN_MakePositivelyOriented(poly)
    ! Changes the polygon in-place(if necessary) to make it positively-oriented
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(inout) :: poly
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cTemp
    integer(kind=AE_INT) :: i, n

    if (.not. PGN_IsPositivelyOriented(poly)) then
      n = size(poly)
      do i = 1, n/2
        cTemp = poly(i)
        poly(i) = poly(n-i+1)
        poly(n-i+1) = cTemp
      end do
    end if

    return
  end subroutine PGN_MakePositivelyOriented


  function PGN_CheckPoint(poly, cZ) result(iResult)
    ! Checks the point cZ against the polygon.
    !
    ! Returns: -1 if the point is outside the polygon
    !           0 if the point is on an edge on at a vertex
    !          +1 if the point is inside the polygon
    !
    ! Adapted from the 1970 FORTRAN code of W.R. Franklin, found at
    ! http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iResult
    ! [ LOCALS ]
    real(kind=AE_REAL), dimension(:), allocatable :: x, y
    real(kind=AE_REAL) :: rchk
    integer(kind=AE_INT) :: iStat, i, j
    logical :: mx, nx, my, ny

    allocate(x(size(poly)), y(size(poly)), stat = iStat)
    if (iStat /= 0) then
      print *, "ALLOCATION FAILED IN PGN_CheckPoint"
      stop 1
    end if

    iResult = -1
    x = real(poly-cZ, AE_REAL)
    y = aimag(poly-cZ)
    do i = 1, size(poly)
      ! This check guarantees that the vertices work
      ! (it failed in some cases in the original)
      if (x(i) == rZERO.and.y(i) == rZERO) then
        iResult = 0
        exit
      end if
      j = 1+mod(i, size(poly))
      mx = x(i) >= rZERO
      nx = x(j) >= rZERO
      my = y(i) >= rZERO
      ny = y(j) >= rZERO
      if (.not.((my.or.ny).and.(mx.or.nx)).or.(mx.and.nx)) then
        continue
      else if (.not.(my.and.ny.and.(mx.or.nx).and..not.(mx.and.nx))) then
        rchk = (y(i)*x(j) - x(i)*y(j)) / (x(j) - x(i))
        if (rchk == rZERO) then
          iResult = 0
          exit
        else if (rchk > rZERO) then
          iResult = -iResult
        end if
      else
        iResult = -iResult
      end if
    end do

    deallocate(x, y)
    return
  end function PGN_CheckPoint


#ifndef __PNPLOGS__
  ! Checks for point-in-polygon using the fast algorithm

  function PGN_Contains(poly, cZ) result(lInside)
    ! Returns .true. if the point cZ lies inside the polygon
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    logical :: lInside

    lInside = PGN_CheckPoint(poly, cZ) >= 0

    return
  end function PGN_Contains

#else
  !! Checks for point-in-Polygon using logs(slow)

  function PGN_Contains(poly, cZ) result(lInside)
    ! Returns .true. if the point cZ lies inside the polygon
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    complex(kind=AE_REAL), intent(in) :: cZ
    ! [ RETURN VALUE ]
    logical :: lInside
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, ii
    real(kind=AE_REAL) :: rSum

    if (size(poly) < 3) then
      lInside = .true.
    else if (any((cZ == poly), 1)) then
      lInside = .true.
    else
      rSum = rZERO
      ii = size(poly)
      do i = 1, size(poly)
        rSum = rSum + aimag(log((cZ-poly(i)) / (cZ-poly(ii))))
        ii = i
      end do
      lInside = (abs(rSum) > rONE)
    end if

    return
  end function PGN_Contains
#endif


  function PGN_Segment(poly, iSeg) result(polySeg)
    ! Returns an array of two points that define segment "iSeg" of the polygon
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    integer(kind=AE_INT), intent(in) :: iSeg
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(2) :: polySeg
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istart, iend

    if (iSeg < 1 .or. iSeg > size(poly)) then
      !    print *, "PGN_Segment: Illegal segment"
      stop
    end if

    istart = iSeg
    iend=iSeg+1
    if (iend > size(poly)) iend=1
    polySeg = (/poly(istart), poly(iend)/)

    return
  end function PGN_Segment


  subroutine PGN_ChopSegment(poly, polySeg, polyChop, iout)
    ! Returns an array of points that define the sub-segments of "polySeg" when
    ! intersected with "poly"
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    complex(kind=AE_REAL), dimension(2), intent(in) :: polySeg
    complex(kind=AE_REAL), dimension(:), intent(inout) :: polyChop
    integer(kind=AE_INT), intent(out) :: iout
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(2) :: seg
    complex(kind=AE_REAL) :: A, B, C, D, P
    real(kind=AE_REAL) :: r, s, denom, r_num, s_num, distance
    integer(kind=AE_INT) :: i, ii, n

    polyChop(1:2) = polySeg
    iout = 2
    n = size(poly)
    A = polySeg(1)
    B = polySeg(2)
    do i = 1, n
      seg = PGN_Segment(poly, i)
      C = seg(1)
      D = seg(2)
      denom = real(B-A)*aimag(D-C) - aimag(B-A)*real(D-C)
      r_num = aimag(A-C)*real(D-C) - real(A-C)*aimag(D-C)
      s_num = aimag(A-C)*real(B-A) - real(A-C)*aimag(B-A)
      if (denom == rZERO) then
        ! Inserts NO points -- this could be trouble in extreme cases...
        cycle
      end if
      r = r_num / denom
      s = s_num / denom

      if ((rZERO < r .and. r < rONE) .and. (rZERO < s .and. s < rONE)) then
        iout = iout+1
        P = A + r*(B-A)
        distance = abs(P-polyChop(1))
        ! Avoid the silly zero-pass DO loop
        ii = 2
        do
          if (distance < abs(polyChop(ii)-polyChop(1))) then
            polyChop(ii+1:iout) = polyChop(ii:iout-1)
            polyChop(ii) = P
            exit
          end if
          ii = ii+1
          if (ii > iout-1) exit
        end do
      end if
    end do

    return
  end subroutine PGN_ChopSegment


  function PGN_LeftmostIndex(poly) result(iLeft)
    ! Returns the index of the left-most point in the polygon
    ! If more than one point share the smallest x, the smaller y wins
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly
    ! [ RETURN VALUE ]
    integer(kind=AE_INT) :: iLeft
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i

    iLeft = 1
    do i = 2, size(poly)
      if (real(poly(i)) < real(poly(iLeft))) then
        iLeft = i
      else if (real(poly(i)) == real(poly(iLeft)) .and. &
             aimag(poly(i)) < aimag(poly(iLeft))) then
        iLeft = i
      end if
    end do

    return
  end function PGN_LeftmostIndex


  subroutine PGN_OrientLeftmost(poly)
    ! Performs all the steps needed to make sure that the polygon is oriented for
    ! convenience in AEM codes. That is, (1) make it positively-oriented and(2)
    ! put the points in order such that the leftmost point in the polygon is the
    ! first in the buffer.
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), dimension(:), intent(inout) :: poly
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(:), allocatable :: localPoly
    integer(kind=AE_INT) :: iLeft, iTo, iFrom, n

    call PGN_MakePositivelyOriented(poly)
    iLeft = PGN_LeftmostIndex(poly)
    if (iLeft /= 1) then
      ! Copy into a new buffer to fix the bugger up
      n = size(poly)
      allocate(localPoly(n))
      localPoly = poly
      !    print *, "local ", localPoly
      iFrom = iLeft
      do iTo = 1, n
        poly(iTo) = localPoly(iFrom)
        iFrom = iFrom+1
        if (iFrom > n) iFrom = 1
      end do
      deallocate(localPoly)
    end if

    return
  end subroutine PGN_OrientLeftmost


  function PGN_ClipArea(poly, clip) result(rArea)
    ! Clips the source polygon 'poly' against the clipping polygon 'clip'
    ! and computes the area of the 'poly' that is within 'clip'.
    !
    ! Uses the clipping routines in polypack.f
    !
    complex(kind=AE_REAL), dimension(:), intent(in) :: poly, clip
    real(kind=AE_REAL) :: rArea
    ! Arrays for calls to PPINTR
    real(kind=AE_REAL), dimension(:), allocatable :: xccp, yccp, xcsp, ycsp, rwrk
    integer(kind=AE_INT), dimension(:), allocatable :: iwrk
    integer(kind=AE_INT) :: nccp, ncsp, nwrk, ierr, nbts
    external :: PPPPAP, PPINTR
    integer(kind=AE_INT) :: i

    if (size(clip) == 0) then
      rArea = PGN_Area(poly)
    else
      nccp = size(clip, 1)
      ncsp = size(poly, 1)
      nbts = 24
      nwrk = 20*(nccp+ncsp)
      allocate(xccp(nccp), yccp(nccp), xcsp(ncsp), ycsp(ncsp), rwrk(nwrk), iwrk(nwrk))
      ! Build the clip polygon
      xccp = real(clip, AE_REAL)
      yccp = aimag(clip)
      call PPPPAP(xccp, yccp, nccp, nbts)
      ! Build the subject polygon
      xcsp = real(poly, AE_REAL)
      ycsp = aimag(poly)
      call PPPPAP(xcsp, ycsp, ncsp, nbts)
      ! Now, sum up the area!
      rPGNTempArea = rZERO
      call PPINTR(xccp, yccp, nccp, xcsp, ycsp, ncsp, rwrk, iwrk, nwrk, PGN_AccumulateArea, ierr)
      if (ierr /= 0) then
        print *, 'INTERNAL ERROR -- PPINTR RETURNED IERR = ', ierr
        stop 1
      end if
      rArea = rPGNTempArea
      deallocate(xccp, yccp, xcsp, ycsp, rwrk, iwrk)
    end if

    return
  end function PGN_ClipArea


  subroutine PGN_AccumulateArea(xcbl, xcbr, ycob, dxle, dxre, ycot)
    ! Accumulates the area from ppintr, one trapezoid at a time
    real(kind=AE_REAL), intent(in) :: xcbl, xcbr, ycob, dxle, dxre, ycot
    real(kind=AE_REAL) :: xul, xur, a

    xul = xcbl+dxle*(ycot-ycob)
    xur = xcbr+dxre*(ycot-ycob)
    a = 0.5 * ((xcbr-xcbl)+(xur-xul)) * (ycot-ycob)
    rPGNTempArea = rPGNTempArea + a

    return
  end subroutine PGN_AccumulateArea


  subroutine PGN_Build(z1, z2, allow_swap, chktol, poly, n)
    ! Builds a polygon from a collection of line segments. The line segments
    ! are provided in the arrays 'z1' and 'z2', which must be of equal length
    ! and are expected to map out a closed polygon. The internal area must be
    ! on the left-hand side of each segment z1(i)-- > z2(i). The segments need
    ! not be linked end-to-end; this code finds the proper endpoints, checking
    ! within the tolerance 'chktol'. Segments of length shorter than chktol
    ! will be ignored in the construction of the polygon.
    !
    ! If allow_swap(i) is true, the segment may be swapped during the check.
    ! If not, the orientation is required to be as specified.
    !
    ! The resulting polygon is provided in the array 'poly', with the number
    ! of points in the polygon in 'n'. The value 'n' can also specify that an
    ! error occurred:
    !  n = -1 : some of the segments are not on the perimeter
    !  n = -2 : The polygon did not close!
    !  n = -3 : Sizes of z1 and z2 are different
    !  n = -4 : Output buffer poly is too small
    !
    ! This code is intended for constructing the active area of a model in
    ! module AQU of ModAEM.
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in), dimension(:) :: z1, z2
    logical, intent(in), dimension(:) :: allow_swap
    real(kind=AE_REAL), intent(in) :: chktol
    complex(kind=AE_REAL), intent(inout), dimension(:) :: poly
    integer(kind=AE_INT), intent(out) :: n
    ! [ LOCALS ]
    logical, dimension(:), allocatable :: lAvail
    integer(kind=AE_INT) :: i, iStatus
    ! iStatus: 0 = no change this time; 1 = change this time; 2 = closed polygon;

    ! Checks...
    if (size(z1) /= size(z2)) then
      n = -3
      return
    end if
    if (size(poly) < size(z1)) then
      n = -4
      return
    end if
    allocate(lAvail(size(z1)))
    lAvail = abs(z2-z1) > chktol

    ! Find the first active
    poly = cZERO
    n = 0
    ! We'll loop through the data, building segment-by-segment
    do
      iStatus = 0
      do i = 1, size(z1)
        if (.not. lAvail(i)) then
          continue
        else if (n == 0) then
          poly(1) = z1(i)
          poly(2) = z2(i)
          lAvail(i) = .false.
          n = 2
          iStatus = 1
        else if (abs(z1(i)-poly(n)) < chktol) then
          if (abs(z2(i)-poly(1)) < chktol) then
            ! All done!
            iStatus = 2
            lAvail(i) = .false.
            exit
          else
            n = n+1
            poly(n) = z2(i)
            lAvail(i) = .false.
            iStatus = 1
          end if
        else if ( allow_swap(i) ) then
          if ( abs(z2(i)-poly(n)) < chktol ) then
            ! Matched with swapped ends!
            if (abs(z1(i)-poly(1)) < chktol) then
              ! All done!
              iStatus = 2
              lAvail(i) = .false.
              exit
            else
              n = n+1
              poly(n) = z1(i)
              lAvail(i) = .false.
              iStatus = 1
            end if
          end if
        end if
      end do
      ! End of this pass through the data
      if (iStatus == 2) then
        if (any(lAvail)) then
          n = -1      ! Some segments unused!
        end if
        exit
      else if (iStatus == 0) then
        n = -2        ! No closure!
        exit
      end if
    end do

    deallocate(lAvail)

    return
  end subroutine PGN_Build

end module u_polygon
#ifdef __UNITTEST__

program unittest

  use u_constants
  use s_unittest
  use u_polygon

  ! Main program -- checks the functions in module u_polygon
  complex(kind=AE_REAL), dimension(4) :: poly1, poly2, polyChop, testChop
  complex(kind=AE_REAL), dimension(2) :: seg, testseg
  complex(kind=AE_REAL), dimension(5) :: z1, z2, pbuild, pcheck
  logical :: lFail
  real(kind=AE_REAL) :: rV
  integer(kind=AE_INT) :: iV
  logical :: lV

  lFail = .false.
  poly1 = (/cmplx(rZERO, rZERO), cmplx(rONE, rZERO), cmplx(rONE, rONE), cmplx(rZERO, rONE)/)
  poly2 = (/cmplx(-rONE, rZERO), cmplx(rZERO, rONE), cmplx(rONE, rZERO), cmplx(rZERO, -rONE)/)

  print *, "** PERFORMING UNIT TESTS FOR MODULE u_polygon **"
  print *

  print *, "Testing PGN_Perimeter..."
  rV = PGN_Perimeter(poly1)
  call realValueCheck("Return value(poly1): ", rV, 4.0_AE_REAL, lFail)

  print *, "Testing PGN_SignedArea..."
  rV = PGN_SignedArea(poly1)
  call realValueCheck("Return value(poly1): ", rV, 1.0_AE_REAL, lFail)
  rV = PGN_SignedArea(poly2)
  call realValueCheck("Return value(poly2): ", rV, -2.0_AE_REAL, lFail)

  print *, "Testing PGN_Area..."
  rV = PGN_Area(poly2)
  call realValueCheck("Return value(poly2): ", rV, 2.0_AE_REAL, lFail)

  print *, "Testing PGN_IsPositivelyOriented..."
  lV = PGN_IsPositivelyOriented(poly1)
  call logicalValueCheck("Return value(poly1) = ", lV, .true., lFail)
  lV = PGN_IsPositivelyOriented(poly2)
  call logicalValueCheck("Return value(poly2) = ", lV, .false., lFail)

  print *, "Testing PGN_MakePositivelyOriented..."
  call PGN_MakePositivelyOriented(poly2)
  print *, "  Called PGN_MakePositivelyOriented(poly2)"
  lV = PGN_IsPositivelyOriented(poly2)
  call logicalValueCheck("PGN_IsPositivelyOriented return value(poly2) = ", lV, .true., lFail)

  print *, "Testing PGN_CheckPoint..."
  iV = PGN_CheckPoint(poly2, cmplx(rZERO, rZERO, AE_REAL))
  call integerValueCheck("Return value(center of poly2) = ", iV, 1, lFail)
  iV = PGN_CheckPoint(poly2, cmplx(rZERO, -rONE, AE_REAL))
  call integerValueCheck("Return value(vertex of poly2) = ", iV, 0, lFail)
  ! Note that the following test _should_ be guaranteed to work, since
  ! one-half is excactly representable in binary floating-point.
  iV = PGN_CheckPoint(poly2, cmplx(-rHALF, rHALF, AE_REAL))
  call integerValueCheck("Return value(edge of poly2) = ", iV, 0, lFail)
  iV = PGN_CheckPoint(poly2, cmplx(-rTWO, rZERO, AE_REAL))
  call integerValueCheck("Return value(left of poly2) = ", iV, -1, lFail)
  iV = PGN_CheckPoint(poly2, cmplx(rTWO, rZERO, AE_REAL))
  call integerValueCheck("Return value(right of poly2) = ", iV, -1, lFail)
  iV = PGN_CheckPoint(poly2, cmplx(rZERO, rTWO, AE_REAL))
  call integerValueCheck("Return value(above poly2) = ", iV, -1, lFail)
  iV = PGN_CheckPoint(poly2, cmplx(rZERO, -rTWO, AE_REAL))
  call integerValueCheck("Return value(below poly2) = ", iV, -1, lFail)

  print *, "Testing PGN_Segment..."
  seg = PGN_Segment(poly1, 2)
  testseg = (/cmplx(rONE, rZERO, AE_REAL), cmplx(rONE, rONE, AE_REAL)/)
  call complexArrayCheck('Segment 2(poly1)', seg, testseg, 2, lFail)

  print *, "Testing PGN_ChopSegment..."
  testseg = (/cmplx(-rTWO, -rTWO, AE_REAL), cmplx(rTWO, rTWO, AE_REAL)/)
  call PGN_ChopSegment(poly2, testseg, polyChop, iV)
  call integerValueCheck("Number of points in chopped polygon = ", iV, 4, lFail)
  testChop = (/cmplx(-rTWO, -rTWO, AE_REAL), cmplx(-rHALF, -rHALF, AE_REAL), &
             cmplx(rHALF, rHALF, AE_REAL), cmplx(rTWO, rTWO, AE_REAL)/)
  call complexArrayCheck('Segment 2(poly1)', polyChop, testChop, iV, lFail)

  print *, "Testing PGN_ClipArea..."
  rV = PGN_ClipArea(poly1, poly2)
  call realValueCheck("Clipped area(poly1 from poly2): ", rV, 0.5_AE_REAL, lFail)

  print *, "Testing PGN_Build..."
  ! This should work...
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0), (0.0, 0.0)/)
  call PGN_Build(z1, z2, 1.0E-8_AE_REAL, pbuild, iV)
  call integerValueCheck("Build a 5-point polygon n = ", iv, 5, lFail)
  pcheck = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  call complexArrayCheck("New polygon:", pbuild, pcheck, iv, lFail)
  ! This should fail(missing internal connection)
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 0.0), (0.0, 0.0)/)
  call PGN_Build(z1, z2, 1.0E-8_AE_REAL, pbuild, iV)
  call integerValueCheck("Missing internal connection n = ", iv, -1, lFail)
  ! This should fail(unused segment)
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 0.0), (-1.0, 0.0)/)
  call PGN_Build(z1, z2, 1.0E-8_AE_REAL, pbuild, iV)
  call integerValueCheck("Missing closure test n = ", iv, -1, lFail)
  ! This should fail(no closure)
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.5, 4.0), (0.0, 0.0), (0.0, 0.0)/)
  call PGN_Build(z1, z2, 1.0E-8_AE_REAL, pbuild, iV)
  call integerValueCheck("Unused segment test n = ", iv, -2, lFail)
  ! This should fail(z1 not same length as z2)
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.5, 4.0), (0.0, 0.0), (0.0, 0.0)/)
  call PGN_Build(z1(1:4), z2, 1.0E-8_AE_REAL, pbuild, iV)
  call integerValueCheck("Wrong size test n = ", iv, -3, lFail)
  ! This should fail(output buffer too small)
  z1 = (/(0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0)/)
  z2 = (/(1.0, 0.0), (2.0, 3.0), (1.0, 4.0), (0.0, 3.0), (0.0, 0.0)/)
  call PGN_Build(z1, z2, 1.0E-8_AE_REAL, pbuild(1:4), iV)
  call integerValueCheck("Output too small test n = ", iv, -4, lFail)

  if (lFail) then
    print *, "One or more unit tests failed"
  else
    print *, "*** ALL TESTS SUCCESSFUL ***"
  end if

  stop
end program unittest
#endif
