module statsmod

use parametersmod, only : i2,i4,sp,dp

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

end module statsmod
