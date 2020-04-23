module netcdfmod

implicit none

public :: handle_err

integer :: ncstat

contains

!--------------------------------------------------------------------------------

subroutine handle_err(ncstat)

use netcdf

implicit none

!Internal subroutine - checks error status after each netcdf call,
!prints out text message each time an error code is returned. 

integer, intent (in) :: ncstat
    
if(ncstat /= nf90_noerr) then 
  write(0,*)trim(nf90_strerror(ncstat))
  stop
end if

end subroutine handle_err

!--------------------------------------------------------------------------------

end module netcdfmod
