module coordsmod

use parametersmod, only : dp,i4

implicit none

public :: parsecoords

type index
  real(dp)    :: minlon
  real(dp)    :: maxlon
  real(dp)    :: minlat
  real(dp)    :: maxlat
  integer(i4) :: startx
  integer(i4) :: starty
  integer(i4) :: endx
  integer(i4) :: endy
  integer(i4) :: countx
  integer(i4) :: county
end type index

contains

!-----------------------------------------------------------------------------------------------

subroutine parsecoords(coordstring,id)

!subroutine to parse a coordinate string

use parametersmod, only : dp

implicit none

!arguments

character(*), intent(in)  :: coordstring
type(index),  intent(out) :: id

!local variables

real(dp), dimension(4) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

!----

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

id%minlon = val(1)
id%maxlon = val(2)
id%minlat = val(3)
id%maxlat = val(4)

end subroutine parsecoords

!-----------------------------------------------------------------------------------------------

end module coordsmod
