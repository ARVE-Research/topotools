module parametersmod

integer, parameter :: i2 = selected_int_kind(4)
integer, parameter :: i4 = selected_int_kind(8)
integer, parameter :: i8 = selected_int_kind(13)
integer, parameter :: sp = selected_real_kind(4)   !4 for single precision
integer, parameter :: dp = selected_real_kind(13)

real(dp), parameter :: pi     =    3.14159265359d0
real(dp), parameter :: d2r    = pi / 180.d0

end module parametersmod
