module statsmod

use parametersmod, only : i2,i4,sp,dp,d2r,r2d,missing_sp

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

real(sp) function circle_mean(vals)

  !Reference: (Wikipedia) Mean of circular quantities

  implicit none

  real(sp), dimension(:), intent(in) :: vals

  integer(i4) :: n

  real(sp) :: sin_sum
  real(sp) :: cos_sum

  n = size(vals)

  sin_sum = sum(sin(vals * d2r))

  cos_sum = sum(cos(vals * d2r))

  !---

  if (sin_sum == 0. .and. cos_sum == 0.) then ! all vectors cancalled out, atan2(0,0) is undefined

    circle_mean = missing_sp

  else

    circle_mean = atan2((sin_sum/n), (cos_sum/n))

    circle_mean = circle_mean * r2d

    !---

    if (circle_mean < 0.) then

      circle_mean = circle_mean + 360. !Convert from atan2 to degrees north

    end if

    !---

  end if

  !circle_mean = max(circle_mean, 0.) !Correct for rounding error

end function circle_mean

!-----------

real(sp) function circle_stdev(vals)

  ! Reference: (Wikipedia) Directional statistics & wrapped normal distribution; https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf

  implicit none

  real(sp), dimension(:), intent(in) :: vals

  integer(i4) :: n
  real(sp) :: r

  real(sp) :: sin_sum
  real(sp) :: cos_sum

  n = size(vals)

  sin_sum = sum(sin(vals * d2r))

  cos_sum = sum(cos(vals * d2r))

  r = sqrt((sin_sum ** 2) + (cos_sum ** 2))

  !---

  if (r == 0.) then

    circle_stdev = missing_sp

  else

    r = -2. * log(r / real(n))

    r = max(r, 0.) !Correct for rounding error

    r = sqrt(r) * r2d

    circle_stdev = r

  end if

end function circle_stdev

!-----------

end module statsmod
