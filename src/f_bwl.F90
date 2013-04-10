module f_bwl


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
  !
  ! Contains code objects for the computation of influence functions for Bessel wells

  use u_constants
  use u_io
  use u_math

  implicit none

  public

  type :: BWL_WELL
    !! type BWL_WELL
    !!
    !! Type that holds information for one partially-penetrating well
    complex(kind=AE_REAL) :: cZ
    real(kind=AE_REAL) :: rRadius
    real(kind=AE_REAL) :: rAqBot
    real(kind=AE_REAL) :: rAqTop
    real(kind=AE_REAL) :: rPrevAqTop
    real(kind=AE_REAL) :: rScrBot
    real(kind=AE_REAL) :: rScrTop
    real(kind=AE_REAL) :: rK
    real(kind=AE_REAL) :: rKhKv
    integer(kind=AE_INT) :: iNLay
    integer(kind=AE_INT) :: iLayer
    logical :: lPpWell
    real(kind=AE_REAL), dimension(:), allocatable :: rT
    real(kind=AE_REAL), dimension(:), allocatable :: rC
    real(kind=AE_REAL), dimension(:), allocatable :: rLambda
    real(kind=AE_REAL), dimension(:,:), allocatable :: rEigVec
    real(kind=AE_REAL), dimension(:), allocatable :: rBCoef
    real(kind=AE_REAL), dimension(:), allocatable :: rInfl
  end type BWL_WELL

contains

  function BWL_New(io, cZ, rRadius, rScrBot, rScrTop, rAqBot, rAqTop, rK, rKhKv) result(bwl)
    !! Builds the internal data structures for a partially-penetrating well
    !! Side effects -- allocates space as necessary for the ppwell functionality.
    !! Returns a pointer to a BWL_WELL object
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    complex(kind=AE_REAL), intent(in) :: cZ
    real(kind=AE_REAL), intent(in) :: rRadius
    real(kind=AE_REAL), intent(in) :: rAqBot
    real(kind=AE_REAL), intent(in) :: rAqTop
    real(kind=AE_REAL), intent(in) :: rScrBot
    real(kind=AE_REAL), intent(in) :: rScrTop
    real(kind=AE_REAL), intent(in) :: rK
    real(kind=AE_REAL), intent(in) :: rKhKv
    ! [ LOCALS ]
    integer(kind=AE_INT) :: istat
    ! [ RETURN VALUE ]
    type(BWL_WELL), pointer :: bwl
    
    ! Build the new object
    allocate(bwl, stat=istat)
    call IO_Assert(io, (istat==0), "Allocation failed for BWL_WELL object")    
    bwl%cZ = cZ
    bwl%rRadius = rRadius    
    bwl%rAqBot = rAqBot
    bwl%rAqTop = rAqTop
    bwl%rPrevAqTop = -rHUGE  ! Forces recalc
    bwl%rScrBot = rScrBot
    bwl%rScrTop = rScrTop
    bwl%rK = rK
    bwl%rKhKv = rKhKv

    ! Set up the object
    if (bwl%rScrBot <= bwl%rAqBot .and. bwl%rScrTop >= bwl%rAqTop) then
      bwl%lPpWell = .false.
    else if (bwl%rScrBot > bwl%rAqBot .and. bwl%rScrTop >= bwl%rAqTop) then
      bwl%lPpWell = .true.
      bwl%iNLay = 2
      bwl%iLayer = 1
    else if (bwl%rScrBot <= bwl%rAqBot .and. bwl%rScrTop < bwl%rAqTop) then
      bwl%lPpWell = .true.
      bwl%iNLay = 2
      bwl%iLayer = 2
    else
      bwl%lPpWell = .true.
      bwl%iNLay = 3
      bwl%iLayer = 2
    end if
    
    ! Make space in the arrays
    allocate (bwl%rT(bwl%iNLay), bwl%rC(bwl%iNLay+1), bwl%rInfl(bwl%iNLay), bwl%rLambda(bwl%iNLay-1), &
              bwl%rEigVec(bwl%iNlay, bwl%iNLay-1), bwl%rBCoef(bwl%iNLay-1), stat=istat)
    call IO_Assert(io, (istat==0), "Allocation failed creating partially-penetrating well")

    return
  end function BWL_New


  subroutine BWL_Solve(io, bwl, rAqTop, lChange)
    !! Solves a single partially-penetrating well, updating the elevation of the aquifer top 
    !! As necessary. Note: if a previous solution is present for the well with the same aquifer
    !! top elevation, no action is needed.
    !! 
    !! Returns an array of influence functions B_n; n=1,3. If the well results in a two-layer
    !! solution, only two influence functions are returned.
    !!
    !! If the new aquifer top value forces a re-solve, lChange is .true.
    !!
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(BWL_WELL), pointer :: bwl
    real(kind=AE_REAL), intent(in) :: rAqTop
    logical, intent(out) :: lChange
    ! [ LOCALS ]
    integer(kind=AE_INT) :: i, istat, lwork, info, icount, n, m, irow
    integer(kind=AE_INT), dimension(1) :: imin
    real(kind=AE_REAL), dimension(:,:), allocatable :: rA, vl, vr, rAA
    real(kind=AE_REAL), dimension(:), allocatable :: wi, wr, work, rBB
    integer(kind=AE_INT), dimension(:), allocatable :: iPiv
    character (len=1) :: jobvl, jobvr
    
    ! Is the well configured for PP but fully-penetrating?
    if (.not. bwl%lPpWell) then
      bwl%rInfl = rZERO
      lChange = .false.
      return
    end if
    
    ! Do I have a valid solution?
    if (abs(rAqTop-bwl%rPrevAqTop) <= 0.001_AE_REAL) then
      lChange = .false.
      return
    end if
    
    ! Get the current solution values
    lChange = .true.
    bwl%rAqTop = rAqTop
    bwl%rPrevAqTop = bwl%rAqTop
    ! Set up the layer properties
    if (bwl%iNLay == 2 .and. bwl%iLayer == 1) then
      ! Screened at top
      bwl%rT(1) = bwl%rK * (bwl%rAqTop- bwl%rScrBot)
      bwl%rC(1) = rHUGE
      bwl%rT(2) = bwl%rK * (bwl%rScrBot - bwl%rAqBot)
      bwl%rC(2) = rHALF * (bwl%rAqTop-bwl%rAqBot) / (bwl%rK / bwl%rKhKV)
      bwl%rC(3) = rHUGE
    else if (bwl%iNLay == 2 .and. bwl%iLayer == 2) then
      ! Screened at bottom
      bwl%rT(1) = bwl%rK * (bwl%rAqTop- bwl%rScrTop)
      bwl%rC(1) = rHUGE
      bwl%rT(2) = bwl%rK * (bwl%rScrTop - bwl%rAqBot)
      bwl%rC(2) = rHALF * (bwl%rAqTop-bwl%rAqBot) / (bwl%rK / bwl%rKhKV)
      bwl%rC(3) = rHUGE
    else 
      ! Screened in the middle
      bwl%rC(1) = rHUGE
      bwl%rT(1) = bwl%rK * (bwl%rAqTop- bwl%rScrTop)
      bwl%rC(2) = rHALF * (bwl%rAqTop-bwl%rScrBot) / (bwl%rK / bwl%rKhKV)
      bwl%rT(2) = bwl%rK * (bwl%rScrTop - bwl%rScrBot)
      bwl%rC(3) = rHALF * (bwl%rScrTop - bwl%rAqBot) / (bwl%rK / bwl%rKhKV)
      bwl%rT(3) = bwl%rK * (bwl%rScrBot - bwl%rAqBot)
      bwl%rC(4) = rHUGE
    end if
    ! Construct the system matrix for the well
    allocate(rA(bwl%iNLay, bwl%iNLay), vr(bwl%iNLay, bwl%iNLay), vl(bwl%iNLay, bwl%iNLay), &
             wi(bwl%iNLay), wr(bwl%iNLay), work(10*bwl%iNLay), &
             rAA(bwl%iNLay-1, bwl%iNLay-1), rBB(bwl%iNLay-1), &
             iPiv(bwl%iNLay-1), stat=istat)
    call IO_Assert(io, (istat==0), "Could not allocate system matrix")
    rA = rZERO
    do i = 1, bwl%iNLay
      if (i>1) rA(i,i-1) = -rONE/(bwl%rC(I)*bwl%rT(i-1))
      rA(i,i) = rONE/(bwl%rC(i)*bwl%rT(i)) + rONE/(bwl%rC(i+1)*bwl%rT(i))
      if (i<bwl%iNLay) rA(i,i+1) = -rONE/(bwl%rC(i+1)*bwl%rT(i+1))
    end do

    ! Compute eigenvalues and eigenvectors
    jobvl = 'N'
    jobvr = 'V'
    lwork = 10*bwl%iNLay
    call DGEEV(jobvl, jobvr, bwl%iNLay, rA, bwl%iNLay, wr, wi, vl, bwl%iNLay, vr, bwl%iNLay, work, lwork, info)
    
    ! Now, store the lambdas and the eigenvector terms
    imin = minloc(wr)
    icount = 0
    do i = 1, bwl%iNLay
      if (i /= imin(1)) then
        icount = icount+1
        bwl%rLambda(icount) = rONE/sqrt(wr(i))
        bwl%rEigVec(:, icount) = vr(:, i)
      end if    
    end do

    ! Now, build the matrix and solve for the Bessel strength coefficients
    rAA = rZERO
    irow = 0
    do n = 1, bwl%iNLay
      if (n /= bwl%iLayer) then
        irow = irow+1
        do m = 1, bwl%iNLay-1
          rAA(irow,m) = rAA(irow,m) + bwl%rEigVec(n, m)
          rBB(irow) = bwl%rT(n) / (rTWO * rPI * sum(bwl%rT))
        end do
      end if
    end do
    call DGESV(bwl%iNLay-1, 1, rAA, bwl%iNLay-1, iPiv, rBB, bwl%iNLay-1, info)
    call IO_Assert(io, (INFO == 0), "Matrix solution failed for Bessel strengths")
    bwl%rBCoef = rBB

    ! Calculate the influence functions in the layers
    bwl%rInfl = rZERO
    do n=1, bwl%iNLay
      do m=1, bwl%iNLay-1
        bwl%rInfl(n) = bwl%rInfl(n) - bwl%rBCoef(m) * &
                                      rBesselK0(bwl%rRadius/bwl%rLambda(m)) * bwl%rEigVec(n,m) / bwl%rT(n)
      end do
    end do
    
    deallocate(rA, vr, vl, wi, wr, work, rAA, rBB, iPiv)
    return
  end subroutine BWL_Solve
  
  
  function BWL_GetCoefficientMultiplier(io, bwl) result(rMultiplier)
    !! Returns the matrix equation multiplier for the well bwl. This is computed based on the 
    !! ratio of the transmissivity of the layer with the well in it to the total transmissivity
    !! at the well. 
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(BWL_WELL), pointer :: bwl
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rMultiplier
    
    rMultiplier = bwl%rT(bwl%iLayer) / sum(bwl%rT)
    
    return    
  end function BWL_GetCoefficientMultiplier


  function BWL_ComputeCoefficient(io, bwl) result(rAji)
    !! Returns the contribution of bwl to the matrix coefficient, based on the total transmissivity 
    !! at the well. 
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(BWL_WELL), pointer :: bwl
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rAji
    
    rAji = -bwl%rInfl(bwl%iLayer) * sum(bwl%rT)
    
    return    
  end function BWL_ComputeCoefficient


  function BWL_ComputeRHS(io, bwl, rQ) result(rBj)
    !! Returns the contribution of bwl to the RHS, based on the total transmissivity 
    !! at the well. 
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(BWL_WELL), pointer :: bwl
    real(kind=AE_REAL), intent(in) :: rQ
    ! [ RETURN VALUE ]
    real(kind=AE_REAL) :: rBj
    
    rBj = -rQ * bwl%rInfl(bwl%iLayer) * sum(bwl%rT)
    
    return    
  end function BWL_ComputeRHS

end module f_bwl
