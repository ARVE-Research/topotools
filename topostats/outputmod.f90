module outputmod

use parametersmod, only : i2,i4,sp

implicit none

public :: genoutfile
public :: putlonlat

real(sp),    parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

contains

!-------------------------------------------------------------------------------------------------

subroutine genoutfile(outfile,id,llim,chunks,ofid)

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
real(sp), dimension(:), intent(in) :: llim

!local variables

integer(i4) :: dimid
integer(i4) :: varid

!real(dp), allocatable, dimension(:) :: lon
!real(dp), allocatable, dimension(:) :: lat

real(dp), dimension(2) :: xrange
real(dp), dimension(2) :: yrange

integer :: nclass

character(8)  :: today
character(10) :: now

integer(i4), allocatable, dimension(:) :: dimids

integer, parameter :: ndims = 3

!---

nclass = size(llim)

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
!lon

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
!lat

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
!slope classes

ncstat = nf90_def_dim(ofid,'slopeclass',nclass,dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

dimids(3) = dimid

ncstat = nf90_def_var(ofid,'slopeclass',nf90_float,dimid,varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','lower limit of slope class threshold')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m m-1')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
!elevation (lon,lat)

ncstat = nf90_def_var(ofid,'elev',nf90_short,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','median elevation above mean sea level')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
!slope (lon,lat)

ncstat = nf90_def_var(ofid,'slope',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','median terrain slope')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m m-1')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
!standard deviation of elevation (lon,lat)

ncstat = nf90_def_var(ofid,'elev_stdev',nf90_short,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','standard deviation of elevation')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'scale_factor',0.1)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'add_offset',0.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
!standard deviation of slope (lon,lat)

ncstat = nf90_def_var(ofid,'slope_stdev',nf90_float,dimids(1:2),varid,chunksizes=chunks(1:2),deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','standard deviation of terrain slope')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','m m-1')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!----
!slope class fraction (lon,lat,slopeclass)

ncstat = nf90_def_var(ofid,'classfrac',nf90_float,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.true.)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','fraction of the gridcell in each slope class')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','fraction')
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'_FillValue',missing_sp)
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

