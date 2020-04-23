module utilitiesmod

!recoded from Numerical Recipes

use parametersmod, only : i4,sp,dp

implicit none

!----------------------------------------------------------------------------------------
!interfaces for overloaded functions

interface iminloc

  module procedure iminloc_i4,iminloc_sp,iminloc_dp

end interface iminloc

!--

interface imaxloc

  module procedure imaxloc_i4,imaxloc_sp,imaxloc_dp

end interface imaxloc

!--

interface assert_eq

  module procedure assert_eq_i4,assert_eq_sp
    
end interface assert_eq

!--

interface outerprod

  module procedure outerprod_sp,outerprod_dp

end interface outerprod

!--

interface swap

  module procedure swap_i4,swap_sp,swap_i4_vect,swap_sp_vect

end interface swap

!----------------------------------------------------------------------------------------

contains

!-----------

integer(i4) function iminloc_i4(arr)

  implicit none

  integer(i4), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imin

  !--

  imin = minloc(arr)

  iminloc_i4 = imin(1)

end function iminloc_i4

!-----------

integer(i4) function iminloc_sp(arr)

  implicit none

  real(sp), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imin

  !--

  imin = minloc(arr)

  iminloc_sp = imin(1)

end function iminloc_sp

!-----------

integer(i4) function iminloc_dp(arr)

  implicit none

  real(dp), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imin

  !--

  imin = minloc(arr)

  iminloc_dp = imin(1)

end function iminloc_dp

!-----------

integer(i4) function imaxloc_i4(arr)

  implicit none

  integer(i4), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imax

  !--

  imax = maxloc(arr)

  imaxloc_i4 = imax(1)

end function imaxloc_i4

!-----------

integer(i4) function imaxloc_sp(arr)

  implicit none

  real(sp), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imax

  !--

  imax = maxloc(arr)

  imaxloc_sp = imax(1)

end function imaxloc_sp

!-----------

integer(i4) function imaxloc_dp(arr)

  implicit none

  real(dp), dimension(:), intent(in) :: arr

  integer(i4), dimension(1) :: imax

  !--

  imax = maxloc(arr)

  imaxloc_dp = imax(1)

end function imaxloc_dp

!-----------

integer(i4) function ifirstloc(mask)
  
  implicit none
  
  logical, dimension(:), intent(in) :: mask
  
  integer(i4) :: i
  
  !--
  
  do i = 1,size(mask)
    if (mask(i)) then
      ifirstloc = i
      return
    end if
  end do

  ifirstloc = i

end function ifirstloc

!-----------

integer(i4) function assert_eq_i4(val1,val2,errormsg)
  
  implicit none
  
  integer(i4),  intent(in) :: val1
  integer(i4),  intent(in) :: val2
  character(*), intent(in) :: errormsg

  !--

  if (val1 == val2) then
    assert_eq_i4 = val1
  else
    write(0,*)errormsg
    stop
  end if

end function assert_eq_i4

!-----------

real(sp) function assert_eq_sp(val1,val2,errormsg)
  
  implicit none
  
  real(sp),     intent(in) :: val1
  real(sp),     intent(in) :: val2
  character(*), intent(in) :: errormsg
  
  !--
  
  if (val1 == val2) then
    assert_eq_sp = val1
  else
    write(0,*)errormsg
    stop
  end if

end function assert_eq_sp

!-----------

function outerprod_sp(a,b)
  
  implicit none
  
  real(sp), dimension(:), intent(in) :: a
  real(sp), dimension(:), intent(in) :: b
  
  real(sp), dimension(size(a),size(b)) :: outerprod_sp
  
  !--
  
  outerprod_sp = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
  
end function outerprod_sp

!-----------

function outerprod_dp(a,b)
  
  implicit none
  
  real(dp), dimension(:), intent(in) :: a
  real(dp), dimension(:), intent(in) :: b
  
  real(dp), dimension(size(a),size(b)) :: outerprod_dp
  
  !--
  
  outerprod_dp = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
  
end function outerprod_dp
  
!-----------

subroutine swap_i4(a,b)

implicit none

integer(i4), intent(inout) :: a
integer(i4), intent(inout) :: b

integer(i4) :: dum

!--

dum = a

a = b

b = dum

end subroutine swap_i4

!-----------

subroutine swap_sp(a,b)

implicit none

real(sp), intent(inout) :: a
real(sp), intent(inout) :: b

real(sp) :: dum

!--

dum = a

a = b

b = dum

end subroutine swap_sp

!-----------

subroutine swap_i4_vect(a,b,mask)

implicit none

integer(i4), dimension(:), intent(inout) :: a
integer(i4), dimension(:), intent(inout) :: b
logical,     dimension(:), intent(in), optional    :: mask

integer(i4), dimension(size(a)) :: dum

!--

where(mask)

  dum = a

  a = b

  b = dum

end where

end subroutine swap_i4_vect

!-----------

subroutine swap_sp_vect(a,b,mask)

implicit none

real(sp), dimension(:), intent(inout) :: a
real(sp), dimension(:), intent(inout) :: b
logical,     dimension(:), intent(in), optional    :: mask

real(sp), dimension(size(a)) :: dum

!--

where(mask)

  dum = a

  a = b

  b = dum

end where

end subroutine swap_sp_vect

!-----------

end module utilitiesmod
