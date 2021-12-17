module calcfracmod

use parametersmod, only : i2,i4,sp,dp
use calcstatsmod,  only : classbin,nclasses_slope,nclasses_aspect,nclasses_cti,nclasses_hand

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

  !Aspect bin (by 8 compass directions)

  llim%aspect_bin(1) = 22.5

  do i = 2, nclasses_aspect
    llim%aspect_bin(i) = real(22.5 + 45. * (i-1))
  end do

  !CTI bin (Fig.2 in Marthews et al. 2015)

  do i = 1, 10
    llim%cti_bin(i) = real(i)
  end do

  llim%cti_bin(11) = 15.
  llim%cti_bin(12:13) = 20.

  ! HAND bin (using exponential values from exp(-2) to exp(6))

  llim%hand_bin(1) = 0.

  do i = 2, nclasses_hand-1
    llim%hand_bin(i) = exp(real(i) - 4.)
  end do

  llim%hand_bin(nclasses_hand) = exp(6.)

end subroutine calcfrac

end module calcfracmod
