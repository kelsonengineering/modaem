module a_stdio

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

  !! Provides STDIO support for external programs

  use u_constants
  use u_io
  use m_aem

  implicit none

  public

  integer(kind=AE_INT), parameter :: STDIN=5
  integer(kind=AE_INT), parameter :: STDOUT=6
  integer(kind=AE_INT), parameter :: STDERR=7

contains


  subroutine STD_IO(io, aem)
    !! Top-level STDIO accessor routine
    ! [ ARGUMENTS ]
    type(IO_STATUS), pointer :: io
    type(AEM_DOMAIN), pointer :: aem

    call IO_MessageText(io, "[ENTERING STDIO]")
    ! Someday there'll be code in here <grin>
    call IO_MessageText(io, "[LEAVING STDIO]")

    return
  end subroutine STD_IO

end module a_stdio
