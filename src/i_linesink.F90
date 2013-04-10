module i_linesink

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


  !! module i_linesink
  !!
  !! Module of data structures and functions for second order
  !! linesinks of complex strength.  The real part of the strength
  !! is the 'linesink' strength(jump in the streamfunction) and the
  !! imaginary part is the 'doublet' strength(jump in the
  !! potential)
  !!
  !! This module encapsulates containers of linesink elements,
  !! organized for efficient computation on vector or parallel
  !! machinery.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!
  use u_constants

  implicit none

  private

  public :: cILS_InfluenceP, &
           cILS_InfluenceW, &
           cILS_InfluenceF, &
           cILS_InfluenceG, &
           cILS_InfluenceQ, &
           cILS_InfluenceB

  real(kind=AE_REAL) :: END_TOLERANCE = 1.0e-6_AE_REAL


contains

  pure function cILS_InfluenceP(cMapZ, cZL) result(cP)
    !! complex function cILS_InfluenceP
    !!
    !! Computes the complex potential influence functions for a single
    !! first-order linesink at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the linesink element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the first-order linesink, one influence function is computed,
    !! in terms of the sink density of the element. It is stored in the
    !! complex array cP.
    !!
    !! Calling Sequence:
    !!    cP = cILS_InfluenceP(cMapZ)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ
    !!             The complex coordinates of the point in question
    !!             NOTE: This version does not handle the case where the
    !!             potential is requested at a vertex along a string of
    !!             linesinks(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!   (in)    complex :: cZL
    !!             The complex quantity(z_2-z_1)/2, where z1 and z2 are
    !!             the model-coordinates of the ends of the linesink
    !!
    !! Returns:
    !!           complex :: cP(1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ, cZL
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cP
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cOm1, cOm2

    ! Compute the function from [Strack, 1989]
    ! NOTE: this version does not use far-field expansions at this time
    cOm1 = cZERO
    cOm2 = cZERO
    if (abs(cMapZ+rONE) > rZERO) then
      cOm1 = (cMapZ+rONE)*log(cMapZ+rONE)
    end if
    if (abs(cMapZ-rONE) > rZERO) then
      cOm2 = (cMapZ-rONE)*log(cMapZ-rONE)
    end if
    cP(1, 1) = (rHALF*abs(cZL) / rPI) * &
            (cOm1 - cOm2 + rTWO*log(cZL) - rTWO)
    return
  end function cILS_InfluenceP

  pure function cILS_InfluenceW(cMapZ, cZL) result(cW)
    !! complex function cILS_InfluenceW
    !!
    !! Computes the complex discharge influence functions for a single
    !! first-order linesink at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the linesink element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the first-order linesink, one influence function is computed,
    !! in terms of the sink density of the element. It is stored in the
    !! complex array cP.
    !!
    !! Since this function performs its computations in the mapped-Z
    !! domain, it is up to the caller to convert the discharge back to the
    !! model domain:  Q_i = conjg(strength * mapW / (.5 * z2-z1))
    !!
    !! Calling Sequence:
    !!    cW = cILS_InfluenceW(cMapZ)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ
    !!             The complex coordinates of the point in question
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             linesinks(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!   (in)    complex :: cZL
    !!             The complex quantity(z_2-z_1)/2, where z1 and z2 are
    !!             the model-coordinates of the ends of the linesink
    !!
    !! Returns:
    !!           complex :: cW(1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ, cZL
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cW
    ! [ LOCALS ]
    complex(kind=AE_REAL) :: cDF, cDFStar, cDG, cCauchy

    cW(1, 1) = -(rHALF*abs(cZL) / rPI) * log((cMapZ+rONE)/(cMapZ-rONE))

    ! Compute the derivatives of the function from [Strack, 1989]
    ! NOTE: this version does not use far-field expansions at this time

    return
  end function cILS_InfluenceW

  pure function cILS_InfluenceF(cMapZ1, cMapZ2, cZL) result(cF)
    !! complex function cILS_InfluenceF
    !!
    !! Computes the integrated flux influence functions for a single
    !! first-order linesink at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the linesink element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the first-order linesink, one influence function is computed,
    !! in terms of the sink density of the element. It is stored in the
    !! complex array cP.
    !!
    !! Calling Sequence:
    !!    cP = cILS_InfluenceQ(cMapZ1, cMapZ2, cZL)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ1, cMapZ2
    !!             The complex coordinates of the ends of the line segment
    !!             for the integrated flux computation.
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             linesinks(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!   (in)    complex :: cZL
    !!             The complex quantity(z_2-z_1)/2, where z1 and z2 are
    !!             the model-coordinates of the ends of the linesink
    !!
    !! Returns:
    !!           complex :: cP(1)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ1, cMapZ2, cZL
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cF
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(1, 1) :: cP1, cP2

    ! Get the discharge potential values first
    cP1 = cILS_InfluenceP(cMapZ1, cZL)
    cP2 = cILS_InfluenceP(cMapZ2, cZL)
    ! Now, compute the branch cut, if any.  The branch cut will only occur
    ! in the imaginary part of the potential.
    cF = (cP2-cP1) + cILS_InfluenceB(cMapZ1, cMapZ2, cZL)
    ! The caller expects that the real part of cF will contain the contribution of
    ! the real part of the linesink strength to the integrated flux and the imaginary
    ! part, the contribution of the imaginary strength to the integrated flux.  This
    ! requires that the real and imaginary parts of Delta-Psi be switched...
    cF = cmplx(aimag(cF), real(cF), AE_REAL)

    return
  end function cILS_InfluenceF

  pure function cILS_InfluenceG() result(cR)
    !! complex function cILS_InfluenceG
    !!
    !! Computes the recharge influence functions for a single
    !! first-order linesink at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the linesink element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! Since this function satisfies Laplace equation, this influence
    !! function is zero for all values of mapZ
    !!
    !! Calling Sequence:
    !!    cP = cILS_InfluenceG()
    !!
    !! Arguments:
    !!
    !! Returns:
    !!           complex :: cP(3)
    !!             Vector of three influence functions; element 1 has the
    !!             influence of the strength at mapZ = -1, element 2 has the
    !!             influence of the center strength and element 3 has the
    !!             influence of the strength at mapZ = +1
    !!
    ! [ ARGUMENTS ]
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cR
    ! [ LOCALS ]

    cR = cZERO

    return
  end function cILS_InfluenceG

  pure function cILS_InfluenceQ(cZL) result(cQ)
    !! complex function rILS_InfluenceQ
    !!
    !! Computes the extraction rate influence function for a single
    !! well.
    !!
    !! For the well, one influence function is computed.
    !!
    !! NOTE: for the 2-D steady-state well, rQ = 1
    !!
    !! Calling Sequence:
    !!    rQ = rILS_InfluenceQ()
    !!
    !! Arguments:
    !!    (in)   complex :: cZL
    !!             The quantity 0.5*(cZ2-cZ1) where cZ2 and cZ1 are the
    !!             complex coordinates of the ends of the line-sink(as
    !!             stored in e.g. LS0 and LS1.
    !!
    !! Returns:
    !!           real :: rQ
    !!             Influence function for extraction rate at the point
    !!             cZ.
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cZL
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cQ(1, 1)
    ! [ LOCALS ]

    cQ = rTWO * abs(cZL)

    return
  end function cILS_InfluenceQ

  pure function cILS_InfluenceB(cMapZ1, cMapZ2, cZL) result(cB)
    !! complex function cILS_InfluenceB
    !!
    !! Computes the branch cut (in streamfunction) influence function for a single
    !! first-order linesink at the complex coordinate cMapZ. The position
    !! is in the mapped domain, with the linesink element on the real axis
    !! extending from mapZ = (-1, 0) to mapZ = (+1, 0).
    !!
    !! For the first-order linesink, one influence function is computed,
    !! in terms of the sink density of the element. It is stored in the
    !! complex array cP.
    !!
    !! Calling Sequence:
    !!    cP = cILS_InfluenceQ(cMapZ1, cMapZ2, cZL)
    !!
    !! Arguments:
    !!   (in)    complex :: cMapZ1, cMapZ2
    !!             The complex coordinates of the ends of the line segment
    !!             for the integrated flux computation.
    !!             NOTE: This version does not handle the case where the
    !!             discharge is requested at a vertex along a string of
    !!             linesinks(see Strack). It is expected that the calling
    !!             routine has eliminated this possibility.
    !!   (in)    complex :: cZL
    !!             The complex quantity(z_2-z_1)/2, where z1 and z2 are
    !!             the model-coordinates of the ends of the linesink
    !!
    !! Returns:
    !!           complex :: cP(1)
    !!             Vector of (one) influence function
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cMapZ1, cMapZ2, cZL
    ! [ RETURN VALUE ]
    complex(kind=AE_REAL), dimension(1, 1) :: cB
    ! [ LOCALS ]
    complex(kind=AE_REAL), dimension(1, 1) :: cP1, cP2
    real(kind=AE_REAL) :: rXInt
    complex(kind=AE_REAL), dimension(1, 1) :: cBC

    ! Compute the branch cut, if any.  The branch cut will only occur
    ! in the imaginary part of the potential.
    cBC = cZERO
    if ((aimag(cMapZ1) >= rZERO .and. aimag(cMapZ2) < rZERO) .or. &
        (aimag(cMapZ2) >= rZERO .and. aimag(cMapZ1) < rZERO)) then
      rXInt = real(cMapZ1) - aimag(cMapZ1) * ((real(cMapZ2)-real(cMapZ1)) / &
              (aimag(cMapZ2)-aimag(cMapZ1)))
      if (rXInt < -rONE) then
        if (aimag(cMapZ1) < rZERO) then
          cBC = -cmplx(rZERO, rTWO*abs(cZL), AE_REAL)
        else
          cBC = cmplx(rZERO, rTWO*abs(cZL), AE_REAL)
        end if
      else if (rXInt < rONE) then
        if (aimag(cMapZ1) <= rZERO) then
          cBC = cmplx(rZERO, (rXInt-rONE)*abs(cZL), AE_REAL)
        else
          cBC = -cmplx(rZERO, (rXInt-rONE)*abs(cZL), AE_REAL)
        end if
      else
        cBC = cZERO
      end if
    end if
    cB = cBC

    return
  end function cILS_InfluenceB

end module i_linesink

