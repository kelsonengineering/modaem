module i_dipole

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

  !! module i_dipole
  !!
  !! Module of data structures and functions for second order
  !! dipoles of complex strength.  The real part of the strength
  !! is the 'dipole' strength(jump in the streamfunction) and the
  !! imaginary part is the 'doublet' strength(jump in the
  !! potential)
  !!
  !! This module encapsulates containers of dipole elements,
  !! organized for efficient computation on vector or parallel
  !! machinery.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!
  use u_constants

  implicit none

  private

  public :: cIDP_InfluenceP, &
           cIDP_InfluenceW, &
           cIDP_InfluenceF, &
           cIDP_InfluenceG, &
           cIDP_InfluenceQ, &
           cIDP_InfluenceJ


contains

  pure function cIDP_InfluenceP(cMapZ) result(cP)
    !! complex function cIDP_InfluenceP
    !!
    !! Computes the complex potential influence functions for a single
    !! second-order dipole at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the dipole element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the second-order dipole, three influence functions are computed,
    !! one for the strengths at each end of the dipole and one for the
    !! strength at the center. They are stored in the complex array cP.
    !!
    !! Calling Sequence:
    !!    cP = cIDP_InfluenceP(cMapZ)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ
    !!             The complex coordinates of the point in question
    !!             NOTE: This version does not handle the case where the
    !!             potential is requested at a vertex along a string of
    !!             dipoles(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!
    !! Returns:
    !!           complex :: cP(3, 1)
    !!             Row-vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(3, 1) :: cP
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cF, cFStar, cG, cCauchy

    ! Compute the functions from [Strack, 1989]
    cCauchy = log((cMapZ-cONE)/(cMapZ+cONE))
    cF = -rHALF * (cMapZ-cONE) * cCauchy - cONE
    cFStar = (cMapZ*cMapZ - cONE) * cCauchy + rTWO*cMapZ
    cG = rHALF * (cMapZ+cONE) * cCauchy + cONE
    ! Compute the influence functions Pi
    cP(1, 1) = -(cF + rHALF*cFStar) / (rTWO*rPI)
    cP(2, 1) = cFStar / (rTWO*rPI)
    cP(3, 1) = -(cG + rHALF*cFStar) / (rTWO*rPI)

    return
  end function cIDP_InfluenceP

  pure function cIDP_InfluenceW(cMapZ) result(cW)
    !! complex function cIDP_InfluenceW
    !!
    !! Computes the complex discharge influence functions for a single
    !! second-order dipole at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the dipole element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the second-order dipole, three influence functions are computed,
    !! one for the strengths at each end of the dipole and one for the
    !! strength at the center. They are stored in the complex array cP.
    !!
    !! Since this function performs its computations in the mapped-Z
    !! domain, it is up to the caller to convert the discharge back to the
    !! model domain:  Q_i = conjg(strength * mapW / (.5 * z2-z1))
    !!
    !! Calling Sequence:
    !!    cW = cIDP_InfluenceW(cMapZ)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ
    !!             The complex coordinates of the point in question
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             dipoles(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!
    !! Returns:
    !!           complex :: cW(3, 1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(3, 1) :: cW
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cDF, cDFStar, cDG, cCauchy

    ! Compute the derivatives of the functions from [Strack, 1989]
    cCauchy = log((cMapZ-cONE)/(cMapZ+cONE))
    cDF = -rHALF * (cONE - (cMapZ-cONE)/(cMapZ+cONE) + cCauchy)
    cDFStar = (rTWO*cMapZ*cCauchy + rFOUR)
    cDG = rHALF * ((cMapZ+cONE)/(cMapZ-cONE) - cONE + cCauchy)
    ! Compute the influence functions Ri
    cW(1, 1) = (cDF + rHALF*cDFStar) / (rTWO*rPI)
    cW(2, 1) = -cDFStar / (rTWO*rPI)
    cW(3, 1) = (cDG + rHALF*cDFStar) / (rTWO*rPI)

    return
  end function cIDP_InfluenceW

  pure function cIDP_InfluenceF(cMapZ1, cMapZ2) result(cF)
    !! complex function cIDP_InfluenceF
    !!
    !! Computes the integrated flux influence functions for a single
    !! second-order dipole at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the dipole element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the second-order dipole, three influence functions are computed,
    !! one for the strengths at each end of the dipole and one for the
    !! strength at the center. They are stored in the complex array cP.
    !!
    !! Calling Sequence:
    !!    cP = cIDP_InfluenceQ(cMapZ1, cMapZ2)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ1, cMapZ2
    !!             The complex coordinates of the ends of the line segment
    !!             for the integrated flux computation.
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             dipoles(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!
    !! Returns:
    !!           complex :: cP(3, 1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ1, cMapZ2
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(3, 1) :: cF
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(3, 1) :: cP1, cP2
    real(kind=AE_REAL) :: rXInt
    complex(kind=AE_REAL), dimension(3, 1) :: cBC

    ! Get the discharge potential values first
    cP1 = cIDP_InfluenceP(cMapZ1)
    cP2 = cIDP_InfluenceP(cMapZ2)
    ! Now, compute the branch cut, if any.  The branch cut will only occur
    ! in the imaginary part of the potential.
    cBC = cZERO
    if ((aimag(cMapZ1) >= rZERO .and. aimag(cMapZ2) < rZERO) .or. &
        (aimag(cMapZ2) >= rZERO .and. aimag(cMapZ1) < rZERO)) then
      rXInt = real(cMapZ1) - aimag(cMapZ1) * ((real(cMapZ2)-real(cMapZ1)) / &
              (aimag(cMapZ2)-aimag(cMapZ1)))
      if (rXInt > -rONE .and. rXInt < rONE) then
        if (aimag(cMapZ1) < rZERO) then
          cBC(1, 1) = cmplx(rZERO, (rHALF * (-(rXInt-rONE) - (rXInt*rXInt-rONE))), AE_REAL)
          cBC(2, 1) = cmplx(rZERO, (rXInt*rXInt - rONE), AE_REAL)
          cBC(3, 1) = cmplx(rZERO, (rHALF * ((rXInt+rONE) - (rXInt*rXInt-rONE))), AE_REAL)
        else
          cBC(1, 1) = -cmplx(rZERO, (rHALF * (-(rXInt-rONE) - (rXInt*rXInt-rONE))), AE_REAL)
          cBC(2, 1) = -cmplx(rZERO, (rXInt*rXInt - rONE), AE_REAL)
          cBC(3, 1) = -cmplx(rZERO, (rHALF * ((rXInt+rONE) - (rXInt*rXInt-rONE))), AE_REAL)
        end if
      end if
    end if
    cF = (cP2-cP1) + cBC
    ! The caller expects that the real part of cF will contain the contribution of
    ! the real part of the dipole strength to the integrated flux and the imaginary
    ! part, the contribution of the imaginary strength to the integrated flux.  This
    ! requires that the real and imaginary parts of Delta-Psi be switched...
    cF = cmplx(aimag(cF), -real(cF), AE_REAL)

    return
  end function cIDP_InfluenceF

  pure function cIDP_InfluenceG(cMapZ) result(cR)
    !! complex function cIDP_InfluenceG
    !!
    !! Computes the recharge influence functions for a single
    !! second-order dipole at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the dipole element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! Since this function satisfies Laplace equation, this influence
    !! function is zero for all values of mapZ
    !!
    !! Calling Sequence:
    !!    cP = cIDP_InfluenceG(cMapZ)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ
    !!             The complex coordinates of the point in question
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             dipoles(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!
    !! Returns:
    !!           complex :: cP(3, 1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(3, 1) :: cR
    ! [ LOCALS ]

    cR = cZERO

    return
  end function cIDP_InfluenceG

  pure function cIDP_InfluenceQ() result(cQ)
    !! complex function rIDP_InfluenceQ
    !!
    !! Computes the extraction rate influence function for a single
    !! well.
    !!
    !! For the well, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rQ = 1
    !!
    !! Calling Sequence:
    !!    rQ = rIDP_InfluenceQ()
    !!
    !! Arguments:
    !!           (none)
    !! Returns:
    !!           real :: cQ
    !!             Influence function for extraction rate at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ(3, 1)
    ! [ LOCALS ]

    cQ = rZERO

    return
  end function cIDP_InfluenceQ

  pure function cIDP_InfluenceJ(rX) result(cJ)
    !! complex function rIDP_InfluenceJ
    !!
    !! Computes the potential/streamfunction jump influence function for a single
    !! dipole.
    !!
    !! For the dipole, three influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rQ = 1
    !!
    !! Calling Sequence:
    !!    rJ = rIDP_InfluenceJ(rX)
    !!
    !! Arguments:
    !!   (in)    real :: rX
    !!             Real part of the mapped 'Z' coordinate along the element
    !!             (assumes that the imaginary part is zero)
    !! Returns:
    !!           real :: cJ
    !!             Influence functions for the jumps at the point rX.
    !!
    ! [ ARGUMENTS ]
    real(kind=AE_REAL), intent(in) :: rX
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cJ(3, 1)
    ! [ LOCALS ]
    real(kind=AE_REAL) :: rJF, rJFStar, rJG

    rJF = -rHALF * (rX - rONE)
    rJFStar = rX*rX - rONE
    rJG = rHALF * (rX + rONE)

    cJ(1, 1) = cmplx(rJF+rHALF*rJFStar, rJF+rHALF*rJFStar, AE_REAL)
    cJ(2, 1) = cmplx(-rJFStar, -rJFStar, AE_REAL)
    cJ(3, 1) = cmplx(rJG+rHALF*rJFStar, rJG+rHALF*rJFStar, AE_REAL)

    return
  end function cIDP_InfluenceJ

end module i_dipole
