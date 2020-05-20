module statsmod

use parametersmod, only : i2,i4,sp,dp,d2r,r2d

implicit none

!module to calculate median and standard deviation of an input vector

!----------------------------------------------------------------------------------------
!interfaces for overloaded functions

interface median

  module procedure median_i2,median_sp

end interface median

!--

interface stdev

  module procedure stdev_i2,stdev_sp

end interface stdev

!----------------------------------------------------------------------------------------

contains

!-----------

integer(i2) function median_i2(vals)

  implicit none

  integer(i2), dimension(:) :: vals

  !--

  median_i2 = nint(median_sp(real(vals)))

end function median_i2

!-----------

real(sp) function median_sp(vals)

  use sortmod, only : sort

  implicit none

  real(sp), dimension(:) :: vals

  integer(i4) :: n
  integer(i4) :: i

  !--

  n = size(vals)
  i = 1 + n / 2

  call sort(vals)

  if (mod(n,2) == 1) then !there are a odd number of values

    median_sp = vals(i)

  else

    median_sp = 0.5 * (vals(i-1) + vals(i))

  end if

end function median_sp

!-----------

integer(i2) function stdev_i2(vals)

  implicit none

  integer(i2), dimension(:), intent(in) :: vals

  !--

  stdev_i2 = nint(stdev_sp(real(vals)))

end function stdev_i2

!-----------

real(sp) function stdev_sp(vals)

  implicit none

  real(sp), dimension(:), intent(in) :: vals

  integer(i4) :: n

  real(sp) :: rn
  real(sp) :: mean

  !--

  n = size(vals)

  rn = 1. / real(n)

  mean = rn * sum(vals)

  stdev_sp = sqrt(rn * sum((vals - mean)**2))

end function stdev_sp

!-----------

real(sp) function vector_mean(vals)

  implicit none

  real(sp), dimension(:), intent(in) :: vals

  real(sp) :: sin_sum
  real(sp) :: cos_sum

  sin_sum = sum(sin(vals * d2r))

  cos_sum = sum(cos(vals * d2r))

  vector_mean = atan2(sin_sum, cos_sum)

  vector_mean = vector_mean * r2d + 180.

  vector_mean = max(vector_mean, 0.) !Correct for rounding error

end function vector_mean

!-----------

end module statsmod
