module calcfracmod

use parametersmod, only : i2,i4,sp,dp
use calcstatsmod,  only : classbin,nclasses_slope,nclasses_aspect,nclasses_cti

contains

!-----------------

subroutine calcfrac(llim)

use parametersmod, only : i2,i4,sp,dp

type(classbin), intent(inout) :: llim
integer :: i

  !Slope bin

  llim%slope_bin(1:2) = 0.

  do i = 3, nclasses_slope
    llim%slope_bin(i) = 10.**(0.2 * real(i-1) - 2.4)
  end do

  !Aspect bin

  do i = 1, nclasses_aspect
    llim%aspect_bin(i) = real((360/nclasses_aspect) * i)
  end do

  !CTI bin

  do i = 1, 10
    llim%cti_bin(i) = real(i)
  end do

  llim%cti_bin(11) = 15.
  llim%cti_bin(12:13) = 20.

end subroutine calcfrac

end module calcfracmod
