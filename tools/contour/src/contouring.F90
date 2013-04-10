module contouring

  ! ModAEM 1.4pre1
  ! Copyright (c) 2001 WHPA Inc.
  !
  ! This program is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License
  ! as published by the Free Software Foundation; either version 2
  ! of the License, or (at your option) any later version.
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


! This is an F90 wrapper for the CONTUR routine from SLWL

use u_constants

implicit none 

private

  integer (4),parameter :: MGSZ = 400
  real (8) :: rgrx1,rgry1                ! world coordinates of lower left of graphics window
  real (8) :: rgrx2,rgry2                ! world coordinates of upper right of graphics window
  real (8) :: rgrmax,rgrmin              ! max and min of any scalar field plotted on the grid
  integer (4) :: ngrx,ngry               ! number of grid points in x and y directions
  real (8),dimension(MGSZ) :: rgrxg      ! grid x-coordinates
  real (8),dimension(MGSZ) :: rgryg      ! grid y-coordinates
  real (8),dimension(MGSZ,MGSZ) :: rgrid ! grid z-values
  real (8) :: rgrdx,rgrdy                ! grid point spacing in x and y directions
  integer (4),dimension(1+(MGSZ*MGSZ-1)/31) :: iSon

public :: CTR_ReadGRD,CTR_DrawContours

contains

subroutine CTR_ReadGRD(iInLU)
  implicit none
  ! Reads the specified SURFER-compatible grid file from iInLU
  !
  ! VAKelson, October, 1998
  integer (kind=AE_INT),intent(in) :: iInLU
  character (len=4) s4Temp
  integer (kind=AE_INT) :: i,j

  read(unit=iInLU,fmt="(a4)") s4Temp
  if (s4Temp /= "DSAA" ) then
    print *,"Illegal grid file"
    stop
  endif

  ! Read the grid dimensions and then the extents in x,y and z
  read (unit=iInLU,fmt=*) ngrx,ngry
  read (unit=iInLU,fmt=*) rgrx1,rgrx2
  read (unit=iInLU,fmt=*) rgry1,rgry2
  read (unit=iInLU,fmt=*) rgrmin,rgrmax
  ! Compute the x- and y-xoordinates of the grid rows and columns
  rgrdx = (rgrx2-rgrx1)/(ngrx-1)
  do i=1,ngrx
    rgrxg(i) = rgrx1+(i-1)*rgrdx
  end do
  rgrdy = (rgry2-rgry1)/(ngry-1)
  do j=1,ngry
    rgryg(i) = rgry1+(i-1)*rgrdy
  end do

  ! Now, read the z data
  do j=1,ngry
    read (unit=iInLU,fmt=*) (rgrid(i,j),i=1,ngrx)
  end do

  return
end subroutine CTR_ReadGRD

subroutine CTR_DrawContours(iOutLU,rZMin,rZIncr,rZMax)
  integer (kind=AE_INT),intent(in) :: iOutLU
  real (kind=AE_REAL),intent(in) :: rZMin,rZIncr,rZMax
  ! Writes a contour file on iOutLU, using subroutine DrawContours
  real (4) :: rLevel

  rLevel = rZMin
  print *,"Computing contours... "
  do
    call Contur(rGrid,MGSZ,ngrx,ngry,rgrx1,rgry1,rgrdx,rgrdy,rLevel,iSon,iOutLU)
    if ( rZIncr <= rZERO ) exit
    rLevel = rLevel+rZIncr
    if ( rLevel > rZMax ) exit
  end do
  print *,"Contouring complete"

  return
end subroutine CTR_DrawContours

end module contouring


