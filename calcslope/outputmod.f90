module outputmod

use parametersmod, only : i2,i4,sp

implicit none

public :: genoutfile
public :: putlonlat

real(sp),    parameter :: missing = -9999.

contains

!-------------------------------------------------------------------------------------------------

subroutine genoutfile(outfile,id,chunks,ofid)

use parametersmod, only : i4,sp,dp
use netcdf
use typesizes
use netcdfmod,     only : ncstat,handle_err
use coordsmod,     only : index

implicit none

character(*), intent(in)  :: outfile
type(index),  intent(in)  :: id
integer, dimension(:), intent(in) :: chunks
integer,      intent(out) :: ofid

!local variables

integer(i4) :: dimid
integer(i4) :: varid

!real(dp), allocatable, dimension(:) :: lon
!real(dp), allocatable, dimension(:) :: lat

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

!integer :: i

character(8)  :: today
character(10) :: now

integer(i4), allocatable, dimension(:) :: dimids

integer, parameter :: ndims = 2

!---

xrange = [id%minlon,id%maxlon]
yrange = [id%minlat,id%maxlat]

!---

!write(0,'(a,a)')'creating output file: ',trim(outfile)

!write(0,*)'create',id%countx,id%county

ncstat = nf90_create(outfile,nf90_hdf5,ofid)
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

!write(0,*)'added global atts'

!-----------
!dimensions

allocate(dimids(ndims))

!----
! lon

ncstat = nf90_def_dim(ofid,'lon',id%countx,dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

dimids(1) = dimid

ncstat = nf90_def_var(ofid,'lon',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
! lat

ncstat = nf90_def_dim(ofid,'lat',id%county,dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

dimids(2) = dimid

ncstat = nf90_def_var(ofid,'lat',nf90_double,dimid,varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
! dem (lon,lat)

ncstat = nf90_def_var(ofid,'elev',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','elevation above mean sea level')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
! slope (lon,lat)

ncstat = nf90_def_var(ofid,'slope',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','terrain slope')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m m-1')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
! aspect (lon,lat)

ncstat = nf90_def_var(ofid,'aspect',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','terrain aspect; flat ground set to missing')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
! topographic wetness index CTI (lon,lat)

ncstat = nf90_def_var(ofid,'cti',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','topographic wetness index (CTI)')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','index')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'method','using GA2 algorithm from Marthews et al., HESS, 2015')
if (ncstat/=nf90_noerr) call handle_err(ncstat)


ncstat = nf90_put_att(ofid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

end subroutine genoutfile

!-------------------------------------------------------------------------------------------------

subroutine putlonlat(ofid,id_olon,id_olat,lon,lat,srtx,srty)

use parametersmod, only : i4,sp,dp
use netcdf
use typesizes
use netcdfmod,     only : ncstat,handle_err

implicit none

integer, intent(in) :: ofid
integer, intent(in) :: id_olon
integer, intent(in) :: id_olat
real(dp), dimension(:), intent(in) :: lon
real(dp), dimension(:), intent(in) :: lat
integer, intent(in) :: srtx
integer, intent(in) :: srty

!---

ncstat = nf90_put_var(ofid,id_olon,lon,start=[srtx])
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,id_olat,lat,start=[srty])
if (ncstat/=nf90_noerr) call handle_err(ncstat)    

end subroutine putlonlat 

!-------------------------------------------------------------------------------------------------

end module outputmod

