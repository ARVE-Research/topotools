module parametersmod

use iso_fortran_env, only : int16,int32,int64,real32,real64

implicit none

integer, parameter :: i2 = int16
integer, parameter :: i4 = int32
integer, parameter :: i8 = int64
integer, parameter :: sp = real32
integer, parameter :: dp = real64

real(dp), parameter :: pi  =    3.14159265358979323846_dp
real(dp), parameter :: d2r = pi / 180._dp

end module parametersmod
