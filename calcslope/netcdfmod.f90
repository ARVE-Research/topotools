module netcdfmod

implicit none

public :: handle_err
public :: makeoutfile

integer :: ncstat

contains

!--------------------------------------------------------------------------------

subroutine handle_err(ncstat)

use netcdf

implicit none

!Internal subroutine - checks error status after each netcdf call,
!prints out text message each time an error code is returned. 

integer, intent (in) :: ncstat

! --

write(0,*)trim(nf90_strerror(ncstat))
stop

end subroutine handle_err

!--------------------------------------------------------------------------------

subroutine makeoutfile(filename,xlen,ylen,lonrange,latrange,ofid)

use parametersmod, only : i4,sp,dp
use netcdf
use typesizes

implicit none

character(*), intent(in) :: filename

integer(i4), intent(in)  :: xlen
integer(i4), intent(in)  :: ylen

real(dp), dimension(2), intent(in) :: lonrange
real(dp), dimension(2), intent(in) :: latrange

integer(i4), intent(out) :: ofid

character(8)  :: today
character(10) :: now

integer(i4), allocatable, dimension(:) :: dimids
integer(i4), allocatable, dimension(:) :: chunks

real(sp), parameter :: missing_sp = -9999.

integer(i4) :: dimid
integer(i4) :: varid

!---

ncstat = nf90_create(filename,nf90_hdf5,ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','DEM slope file')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'created',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----

allocate(dimids(2))

ncstat = nf90_def_dim(ofid,'lon',xlen,dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',lonrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'lat',ylen,dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',latrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----

allocate(chunks(2))

chunks = [600,600]

!----

ncstat = nf90_def_var(ofid,'slope',nf90_float,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','slope')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m m-1')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)
  
end subroutine makeoutfile

!--------------------------------------------------------------------------------

end module netcdfmod
