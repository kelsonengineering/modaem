module u_math

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

  ! This module contains useful computational routines and special functions required
  ! by certain modules

  use u_constants

  implicit none

  public :: rExpIntE1, &
           rPolynomial

  ! Constants
  real(kind=AE_DBL), public, parameter :: rGAMMA = .5772156649_AE_DBL
  ! The cutoff for exp(-rU) ~ 0.0 was adopted from Haitjema - GFLOW
  real(kind=AE_DBL), public, parameter :: rU_CUTOFF = 150.0_AE_REAL
  
  ! Bessel function constants
  integer (kind=AE_INT),public,parameter :: NBesAprox = 8
  real (kind=AE_REAL),parameter,dimension(0:NBesAprox) :: rRange = &
      (/0,1,2,3,4,5,6,7,8/)
  real (kind=AE_REAL),public,parameter :: rRBesConv = 7
  real (kind=AE_REAL),public,parameter :: rRBesConvSq = rRBesConv**2
  
  real (kind=AE_REAL),public,parameter,dimension(0:NBesAprox) :: raBes = &
    (/ -0.500004211065677_AE_REAL, &
       -0.124989431448700_AE_REAL, &
       -0.781685695640975e-2_AE_REAL, &
       -0.216324415010288e-3_AE_REAL, &
       -0.344525393452639e-5_AE_REAL, &
       -0.315133836828774e-7_AE_REAL, &
       -0.296636186265427e-9_AE_REAL, &
       -0.313689942474032e-12_AE_REAL, &
       -0.112031912249579e-13_AE_REAL /)
  real (kind=AE_REAL),public,parameter,dimension(0:NBesAprox) :: rbBes = &
    (/  0.115956920789028_AE_REAL, &
        0.278919134951974_AE_REAL, &
        0.252752008621167e-1_AE_REAL, &
        0.841879407543506e-3_AE_REAL, &
        0.152425102734818e-4_AE_REAL, &
        0.148292488399579e-6_AE_REAL, &
        0.157622547107156e-8_AE_REAL, &
        0.117975437124933e-11_AE_REAL, &
        0.655345107753534e-13_AE_REAL /)

contains


  function rExpIntE1(rU) result(rE1)
    ! Returns the exponential integral of the first kind, E1(U)
    ! From Abramowicz and Stegun, 5.1.53 and 5.1.56
    real(kind=AE_REAL), intent(in) :: rU
    real(kind=AE_REAL) :: rE1
    ! Locals
    real(kind=AE_REAL) :: rU0

    ! Local computations are performed in double precision, return value is single
    ! Coefficient arrays for polynomial computations
    ! Polynomial for 0 <= x <= 1 [A&S 5.1.53]
    real(kind=AE_DBL), dimension(0:5) :: &
                         rA_SmallX = &
                         (/0.00107857_AE_DBL, &
                         -0.00976004_AE_DBL, &
                         0.05519968_AE_DBL, &
                         -0.24991055_AE_DBL, &
                         0.99999193_AE_DBL, &
                         -0.57721566_AE_DBL/)
    ! Polynomial for the numerator of the rational expression for x > 1 [A&S 5.1.56]
    real(kind=AE_DBL), dimension(0:4) :: &
                         rA_LargeX_Numerator = &
                         (/1.0000000000_AE_DBL, &
                         8.5733287401_AE_DBL, &
                         18.0590169730_AE_DBL, &
                         8.6347608925_AE_DBL, &
                         0.2677737343_AE_DBL/)
    ! Polynomial for the denominator of the rational expression for x > 1 [A&S 5.1.56]
    real(kind=AE_DBL), dimension(0:4) :: &
                         rA_LargeX_Denominator = &
                         (/1.0000000000_AE_DBL, &
                         9.5733223454_AE_DBL, &
                         25.6329561486_AE_DBL, &
                         21.0996530827_AE_DBL, &
                         3.9584969228_AE_DBL/)

    ! Check the argument
    if (rU <= rZERO) then
      write (unit=*, fmt="(""rfExpIntE1: Error - argument is less than zero"")")
      rE1 = -HUGE(AE_REAL)
      return
    end if

    ! Compute the E1(U) function
    if (rU == rZERO) then
      rE1 = -HUGE(AE_REAL)
    else if (rU < rONE) then
      rE1 = rPolynomial(rU, rA_SmallX, 5) - log(rU)
    else
      ! The cutoff for exp(-rU) ~ 0.0 was adopted from Haitjema - GFLOW
      if (rU < rU_CUTOFF) then
        rU0 = rU
      else
        rU0 = rZERO
      end if
      rE1 = rPolynomial(rU, rA_LargeX_Numerator, 4) /   &
            (rPolynomial(rU, rA_LargeX_Denominator, 4) * rU * exp(rU0))
    end if

    return
  end function rExpIntE1


  function rPolynomial(rX, rA, iOrder) result(rY)
    ! Computes the iOrder order polynomial rA(0)x^iOrder + rA(1)x^(iOrder-1) + ...
    ! rA = vector of coefficients, dimension(0:iOrder)
    ! Allcomputations are in double precision
    ! Arguments
    real(kind=AE_DBL), intent(in) :: rX
    real(kind=AE_DBL), intent(in), dimension(0:) :: rA
    integer(kind=AE_INT), intent(in) :: iOrder
    real(kind=AE_DBL) :: rY
    ! Locals
    integer(kind=AE_INT) :: i

    ! Range checking
    if (ubound(rA, 1) < iOrder) then
      rY = -HUGE(AE_DBL)
      return
    end if

    rY = rA(0)
    do i = 1, iOrder
      rY = rY*rX + rA(i)
    end do

    return
  end function rPolynomial
  
  
  function rBesselK0(rolam) result(rK0)
    !! Returns the modified Bessel function of order zero and the second kind
    !! for r/lambda where r and lambda are given
    ! [ ARGUMENTS ]
    real(kind=AE_REAL), intent(in) :: rolam
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rK0
    ! [ LOCALS ]
    real(kind=AE_REAL),dimension(0:NBesAprox) :: rsqolamvec
    real(kind=AE_REAL) :: rsqolam
    integer(kind=AE_INT) :: n
    
    rsqolam = rolam * rolam
    if (rsqolam > rRBesConvSq) then
      rK0 = rZERO
    else
      rsqolamvec(0)=rONE
      do n=1,NBesAprox
        rsqolamvec(n) = rsqolamvec(n-1) * rsqolam
      end do
      rK0 = sum ( ( raBes * log(rsqolam) + rbBes ) * rsqolamvec )
    end if
    
    return
  end function rBesselK0


end module u_math

