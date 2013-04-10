module u_grid

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

  ! module u_grdd(GRD)
  ! Data structures and convenient methods for regularly grddded data

  use u_constants
  use u_io

  implicit none

  public

  type :: GRD_GRID
    !! type GRD_GRID
    !!
    !! PUBLIC type that holds the information for regularly-grddded data
    !!
    !! Members:
    !!   real :: rMinX
    !!     Minimum X value
    !!   real :: rMaxX
    !!     Maximum X value
    !!   real :: rMinY
    !!     Minimum Y value
    !!   real :: rMaxY
    !!     Maximum Y value
    !!   real :: rMinZ
    !!     Minimum Z value
    !!   real :: rMaxZ
    !!     Maximum Z value
    !!   integer :: iResX
    !!     Number of points along the X axis
    !!   integer :: iResY
    !!     Number of points along the Y axis
    !!   real :: rDelta
    !!     The spacing between grdd points
    !!   real :: rValues(:, :)
    !!     The values of points in the grd
    !!
    real(kind=AE_REAL) :: rMinX
    real(kind=AE_REAL) :: rMaxX
    real(kind=AE_REAL) :: rMinY
    real(kind=AE_REAL) :: rMaxY
    real(kind=AE_REAL) :: rMinZ
    real(kind=AE_REAL) :: rMaxZ
    integer(kind=AE_INT) :: iResX
    integer(kind=AE_INT) :: iResY
    real(kind=AE_REAL) :: rDelta
    real(kind=AE_REAL), dimension(:, :), pointer :: rValues
  end type GRD_GRID

  ! Define the smalles legal grid window here
  real(kind=AE_REAL), private, parameter :: rSMALL_GRID = 1.0e-3_AE_REAL


contains


  function GRD_Create(io, cLL, cUR, iRes) result(grd)
    !! function GRD_Create
    !!
    !! Creates a new GRD_GRID object
    !!
    !! Calling Sequence:
    !!    grd => GRD_Create(cLL, cUR, iRes)
    !!
    !! Arguments:
    !!    (in)    complex :: cLL
    !!              Lower-left corner of the grdd
    !!    (in)    complex :: cUR
    !!              Upper-rigth corner of the grdd
    !!    (in)    integer :: iRes
    !!              Number of points along the longer axis
    !!
    ! [ ARGUMENTS ]
    complex(kind=AE_REAL), intent(in) :: cLL
    complex(kind=AE_REAL), intent(in) :: cUR
    integer(kind=AE_INT), intent(in) :: iRes
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    type(GRD_GRID), pointer :: grd

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    real(kind=AE_REAL) :: rSizeX, rSizeY

    call IO_Assert(io, (abs(cUR-cLL) > rSMALL_GRID), 'GRD_Create: Bad window')
    call IO_Assert(io, (iRes > 2), 'GRD_Create: Bad resolution')

    allocate(grd, stat = iStat)
    call IO_Assert(io, (iStat == 0), 'GRD_Create: Allocation failed')

    ! Compute the actual data sizes for the grdd and allocate space
    rSizeX = real(cUR-cLL)
    rSizeY = aimag(cUR-cLL)
    if (rSizeX >= rSizeY) then
      grd%iResX = iRes
      grd%rDelta = rSizeX/real(iRes-1, AE_REAL)
      grd%iResY = int(rSizeY/grd%rDelta)+1
    else
      grd%iResY = iRes
      grd%rDelta = rSizeY/real(iRes-1, AE_REAL)
      grd%iResX = int(rSizeX/grd%rDelta)+1
    end if
    allocate(grd%rValues(grd%iResY, grd%iResX), stat = iStat)
    call IO_Assert(io, (iStat == 0), 'MakeGrid: Allocation failed')

    ! Make the real window span the center of the specified window
    grd%rMinX = rHALF * (real(cLL+cUR) - (grd%iResX-1)*grd%rDelta)
    grd%rMaxX = grd%rMinX + (grd%iResX-1)*grd%rDelta
    grd%rMinY = rHALF * (aimag(cLL+cUR) - (grd%iResY-1)*grd%rDelta)
    grd%rMaxY = grd%rMinY + (grd%iResY-1)*grd%rDelta
    grd%rValues = rZERO

    return
  end function GRD_Create


  function GRD_CreateFromGRD(io, sFileName) result(grd)
    !! function GRD_CreateFromGRD
    !!
    !! Creates a new GRD_GRID object from a SURFER(tm)-compatible ASCII .GRD file
    !!
    !! Calling Sequence:
    !!    grd => GRD_CreateFromGRD(sFileName)
    !!
    !! Arguments:
    !!    (in)    character(len=*) :: sFileName
    !!              Name of the SURFER(tm) GRD file
    !!
    ! [ ARGUMENTS ]
    character(len=*), intent(in) :: sFileName
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    type(GRD_GRID), pointer :: grd
    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat
    character(len=4) :: sDSAA
    integer(kind=AE_INT) :: iNX, iNY
    real(kind=AE_REAL) :: rDX, rDY, rMinX, rMaxX, rMinY, rMaxY, rMinZ, rMaxZ

    ! Open up the GRD file
    open(unit=LU_GRID, file = trim(sFileName), status="OLD", iostat=iStat)
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: Could not open " // sFileName // " for input")

    ! Read the grid header
    read (unit=LU_GRID, fmt="(a4)", iostat=iStat) sDSAA
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: I/O Error")
    call IO_Assert(io, (sDSAA == "DSAA"), "GRD_CreateFromGRD: Illegal header record")
    ! Dimensions and limits
    read (unit=LU_GRID, fmt=*, iostat=iStat) iNX, iNY
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: I/O Error")
    read (unit=LU_GRID, fmt=*, iostat=iStat) rMinX, rMaxX
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: I/O Error")
    read (unit=LU_GRID, fmt=*, iostat=iStat) rMinY, rMaxY
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: I/O Error")
    read (unit=LU_GRID, fmt=*, iostat=iStat) rMinZ, rMaxZ
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: I/O Error")
    ! Made it this far... build the object and populate it
    allocate(grd, stat = iStat)
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: Allocation failed")
    grd%iResX = iNX
    grd%iResY = iNY
    grd%rMinX = rMinX
    grd%rMaxX = rMaxX
    grd%rMinY = rMinY
    grd%rMaxY = rMaxY
    grd%rMinZ = rMinZ
    grd%rMaxZ = rMaxZ
    allocate(grd%rValues(iNY, iNX), stat = iStat)
    call IO_Assert(io, (iStat == 0), "GRD_CreateFromGRD: Allocation failed")
    read (unit=LU_GRID, fmt=*, iostat=iStat) grd%rValues

    return
  end function GRD_CreateFromGRD


  subroutine GRD_WriteToGRD(io, grd, sFileName)
    !! function GRD_CreateFromGRD
    !!
    !! Creates a new GRD_GRID object from a SURFER(tm)-compatible ASCII .GRD file
    !!
    !! Calling Sequence:
    !!    grd => GRD_CreateFromGRD(grd, sFileName)
    !!
    !! Arguments:
    !!    (in)    type(GRD_GRID), pointer :: grd
    !!              GRD_GRID object to be written
    !!    (in)    character(len=*) :: sFileName
    !!              Name of the SURFER(tm) GRD file
    !!
    ! [ ARGUMENTS ]
    type(GRD_GRID), pointer :: grd
    character(len=*), intent(in) :: sFileName
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: j, iStat

    open(unit=LU_GRID, file = trim(sFileName), status="NEW", iostat=iStat)
    call IO_Assert(io, (iStat == 0), "MakeGrid: Could not open " // sFileName // " for output")

    call IO_MessageText(io, "  Writing grid to " // trim(sFileName))
    write (unit=LU_GRID, &
           fmt="('DSAA')" &
           )
    write (unit=LU_GRID, &
           fmt="(2(1x, i5))" &
           ) grd%iResX, grd%iResY
    write (unit=LU_GRID, &
           fmt="(2(1x, g12.5))" &
           ) grd%rMinX, grd%rMaxX
    write (unit=LU_GRID, &
           fmt="(2(1x, g12.5))" &
           ) grd%rMinY, grd%rMaxY
    write (unit=LU_GRID, &
           fmt="(2(1x, g12.5))" &
           ) minval(grd%rValues), maxval(grd%rValues)
    do j = 1, grd%iResY
      write (unit=LU_GRID, &
             fmt=* &
             ) grd%rValues(j, :)
    end do
    ! That's it
    close(unit=LU_GRID)

    return
  end subroutine GRD_WriteToGRD


  subroutine GRD_WriteToMatlab(io, grd, sFileName)
    !! function GRD_WriteToMatlab
    !!
    !! Creates a matlab '.m' file for contouring
    !!
    !! Calling Sequence:
    !!    grd => GRD_CreateFromMatlab(grd, sFileName)
    !!
    !! Arguments:
    !!    (in)    type(GRD_GRID), pointer :: grd
    !!              GRD_GRID object to be written
    !!    (in)    character(len=*) :: sFileName
    !!              Name of the Matlab(tm) '.m' file
    !!
    ! [ ARGUMENTS ]
    type(GRD_GRID), pointer :: grd
    character(len=*), intent(in) :: sFileName
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: j, iStat

    open(unit=LU_GRID, file = trim(sFileName), status="NEW", iostat=iStat)
    call IO_Assert(io, (iStat == 0), "MakeGrid: Could not open " // sFileName // " for output")

    call IO_MessageText(io, "  Writing grid to " // trim(sFileName))
    write (unit=LU_GRID, &
           fmt="('[x, y] = meshgrid(', g12.5, ':', g12.5, ':', g12.5, ', ', g12.5, ':', g12.5, ':', g12.5, ')')" &
           ) grd%rMinX, (grd%rMaxX-grd%rMinX)/float(grd%iResX-1), grd%rMaxX, &
           grd%rMinY, (grd%rMaxY-grd%rMinY)/float(grd%iResY-1), grd%rMaxY
    write (unit=LU_GRID, &
           fmt="('f = [')" &
           )
    do j = 1, grd%iResY
      write (unit=LU_GRID, &
             fmt=* &
             ) grd%rValues(j, :)
    end do
    write (unit=LU_OUTPUT, &
           fmt="(']')" &
           )
    ! That's it
    close(unit=LU_GRID)

    return
  end subroutine GRD_WriteToMatlab


  subroutine GRD_Destroy(io, grd)
    !! subroutine GRD_Destroy
    !!
    !! Destroys a GRD_GRID object
    !!
    !! Usage:
    !!        call GRD_Destroy(io, grd)
    !!
    !! Arguments:
    !!   (in)     type(GRD_GRID), pointer :: grd
    !!              The grid to be destroyed
    !!
    ! [ ARGUMENTS ]
    type(GRD_GRID), pointer :: grd
    type(IO_STATUS), pointer :: io

    ! [ LOCALS ]
    integer(kind=AE_INT) :: iStat

    call IO_Assert(io, (associated(grd)), "GRD_Destroy: No GRD_GRID object")
    call IO_Assert(io, (associated(grd%rValues)), "GRD_Destroy: GRD_GRID object has no data")

    deallocate(grd%rValues, stat = iStat)
    call IO_Assert(io, (iStat == 0), "GRD_Destroy: Deallocation error")
    deallocate(grd, stat = iStat)
    call IO_Assert(io, (iStat == 0), "GRD_Destroy: Deallocation error")

    return
  end subroutine GRD_Destroy


  function cGRD_GetZ(io, grd, iRow, iCol) result(cZ)
    !! function cGRD_GetZ
    !!
    !! Gets the complex coordinate of the center of the cell(iRow, iCol)
    !!
    !! Calling Sequence:
    !!    cZ = cGRD_GetZ(grd, iRow, iCol)
    !!
    !! Arguments:
    !!     (in)     type(GRD_GRID), pointer :: grd
    !!                GRD_GRID object to be queried
    !!     (in)     integer :: iRow
    !!                Row number
    !!     (in)     integer :: iCol
    !!                Column number
    !!
    !! Return Value:
    !!              complex :: cZ
    !!                Complex coordinate of the center of cell(iRow, iCol)
    !!                returns HUGE if outside the grid.
    !!
    ! [ ARGUMENTS ]
    type(GRD_GRID), pointer :: grd
    integer(kind=AE_INT), intent(in) :: iRow
    integer(kind=AE_INT), intent(in) :: iCol
    type(IO_STATUS), pointer :: io

    ! [ RETURN VALUE ]
    complex(kind=AE_REAL) :: cZ

    call IO_Assert(io, (associated(grd)), "GRD_GetZ: No GRD_GRID object")

    if (iRow < 1 .or. &
        iRow > grd%iResY .or. &
        iCol < 1 .or. &
        iCol > grd%iResX) then
      cZ = cmplx(HUGE(AE_REAL), HUGE(AE_REAL), AE_REAL)
    else
      cZ = cmplx(grd%rMinY+(grd%iResY-iRow-1)*grd%rDelta, &
           grd%rMinX+(iCol-1)*grd%rDelta, &
           AE_REAL)
    end if

    return
  end function cGRD_GetZ


  subroutine GRD_Lookup(io, grd, cZ, iRow, iCol)
    !! subroutine GRD_Lookup
    !!
    !! Looks up the row, col coordinate of the point cZ
    !!
    !! Calling Sequence:
    !!    call GRD_Lookup(io, grd, cZ, iRow, iCol)
    !!
    !! Arguments:
    !!     (in)     type(GRD_GRID), pointer :: grd
    !!                GRD_GRID object to be queried
    !!     (in)     complex :: cZ
    !!                Complex coordinate
    !!     (out)    integer :: iRow
    !!                Row number or -1 if outside the grid
    !!     (out)    integer :: iCol
    !!                Column number or -1 if outside the grid
    !!
    !!
    ! [ ARGUMENTS ]
    type(GRD_GRID), pointer :: grd
    complex(kind=AE_REAL), intent(in) :: cZ
    integer(kind=AE_INT), intent(out) :: iRow
    integer(kind=AE_INT), intent(out) :: iCol
    type(IO_STATUS), pointer :: io

    call IO_Assert(io, (associated(grd)), "GRD_Lookup: No GRD_GRID object")

    if (real(cZ) < grd%rMinX .or. &
        real(cZ) > grd%rMaxX .or. &
        aimag(cZ) < grd%rMinY .or. &
        aimag(cZ) > grd%rMaxY) then
      iRow = -1
      iCol = -1
    else
      ! Compute iRow and iCol(adjust for off-by one exactly on boundaries)
      iRow = int(grd%rMaxY-aimag(cZ))/grd%rDelta + 1
      if (iRow > grd%iResY) iRow = grd%iResY
      iCol = int(real(cZ)-grd%rMinX)/grd%rDelta + 1
      if (iCol > grd%iResX) iRow = grd%iResX
    end if

    return
  end subroutine GRD_Lookup

end module u_grid

