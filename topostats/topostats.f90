program topostats

use parametersmod, only : i1,i2,i4,sp,dp,missing_i2,missing_sp
use netcdfmod,     only : ncstat,handle_err
use netcdf
use typesizes
use coordsmod,     only : index,parsecoords
use statsmod,      only : median,stdev
use utilitiesmod,  only : iminloc,cellarea
use outputmod,     only : genoutfile
use calcstatsmod,  only : nclasses,llim,statvals,calcstats

implicit none

character(100) :: topofile
character(100) :: landmaskfile
character(100) :: joboptions
character(100) :: outfile

namelist / infiles / topofile, landmaskfile

integer(i4) :: ifid
integer(i4) :: lfid
integer(i4) :: ofid

integer(i4) :: dimid
integer(i4) :: xvarid
integer(i4) :: yvarid

integer(i4) :: id_landmask
integer(i4) :: id_elev_in
integer(i4) :: id_slope_in
integer(i4) :: id_elev_out
integer(i4) :: id_elev_std
integer(i4) :: id_slope_out
integer(i4) :: id_slope_std
integer(i4) :: id_classfrac

integer(i4) :: xpos_out
integer(i4) :: ypos_out

integer(i4) :: ompchunk = 8

integer(i4), parameter :: maxchunklen = 2896

!integer(i4) :: nv
!real(sp)    :: rv

real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat

real(dp),    allocatable, dimension(:)   :: lon_out
real(dp),    allocatable, dimension(:)   :: lat_out

real(sp),    allocatable, dimension(:,:) :: area

integer(i1), allocatable, dimension(:,:) :: landmask
integer(i2), allocatable, dimension(:,:) :: elev_in
real(sp),    allocatable, dimension(:,:) :: slope_in

real(sp),    allocatable, dimension(:,:,:) :: classfrac

real(dp), dimension(2) :: hpx

integer(i4), parameter :: pixchunkmaxf = 10 !max number of input chunks to read at once (along one side)
integer(i4) :: pixchunkmax

integer, dimension(2) :: pixpix
!integer, dimension(2) :: pixchunk

integer(i4) :: a,b
integer(i4) :: i,j
integer(i4) :: x,y
! integer(i4) :: g,h
integer(i4) :: m,n
integer(i4) :: xlen
integer(i4) :: ylen
integer(i4) :: p,q,r,s
integer(i4) :: k

integer(i4) :: varid

type(statvals), allocatable, dimension(:,:) :: stats

integer, dimension(2) :: sblock_in
integer, dimension(2) :: nblock_out  !number of output block
integer, dimension(2) :: sblock_out  !size of the output block

integer, dimension(3) :: schunk_out  !size of the output chunk for netcdf
integer, dimension(3) :: nchunk_out  !number of chunks in the output file

real(sp), dimension(2) :: resratio

integer(i4) :: xoffset
integer(i4) :: yoffset

integer(i4) :: startx
integer(i4) :: starty

character(80) :: coordstring

type(index) :: in
type(index) :: out

!resolutions in seconds

real(dp) :: xres_in
real(dp) :: yres_in
real(dp) :: xres_out
real(dp) :: yres_out

real(dp), dimension(2) :: lonrange
real(dp), dimension(2) :: latrange

character(80) :: status_line

character(10) :: outres

integer(i4), dimension(2) :: nblkchk

!-------------------------------------------------------------------------------------
!set up the bounds for the slope classes

llim(1:2) = 0.

do i = 3,nclasses
  llim(i) = 10.**(0.2 * real(i-1) - 2.4)
end do

!--------

call getarg(1,coordstring)

call parsecoords(coordstring,in)

call getarg(2,outres)

read(outres,*)xres_out

yres_out = xres_out

!---

call getarg(3,joboptions)

open(10,file=joboptions,status='old')

read(10,nml=infiles)

close(10)

ncstat = nf90_open(topofile,nf90_nowrite,ifid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_open(landmaskfile,nf90_nowrite,lfid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---

ncstat = nf90_inq_dimid(ifid,'lon',dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)
  
ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_dimid(ifid,'lat',dimid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)
  
ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---

ncstat = nf90_inq_varid(ifid,'lon',xvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)
  
ncstat = nf90_get_att(ifid,xvarid,'actual_range',lonrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'lat',yvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,yvarid,'actual_range',latrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

xres_in = 3600. * (lonrange(2) - lonrange(1)) / real(xlen)  !units arc seconds
yres_in = 3600. * (latrange(2) - latrange(1)) / real(ylen)

!---
!get indices then repurpose lon and lat

allocate(lon(xlen))
allocate(lat(ylen))

ncstat = nf90_get_var(ifid,xvarid,lon)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(ifid,yvarid,lat)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

in%startx = iminloc(abs(in%minlon - lon))  !start-x
in%starty = iminloc(abs(in%minlat - lat))  !start-y

in%endx = iminloc(abs(in%maxlon - lon))    !end-x
in%endy = iminloc(abs(in%maxlat - lat))    !end-y

!adjust the nodes to select if the nearest node automatically selected is outside the area of interest

if (lon(in%startx) < in%minlon) in%startx = in%startx + 1
if (lon(in%endx)   > in%maxlon) in%endx   = in%endx   - 1

if (lat(in%starty) < in%minlat) in%starty = in%starty + 1
if (lat(in%endy)   > in%maxlat) in%endy   = in%endy   - 1

in%countx = 1 + in%endx - in%startx  !add one to include the start and end points
in%county = 1 + in%endy - in%starty

deallocate(lon)
deallocate(lat)

!---

ncstat = nf90_inq_varid(lfid,'land',id_landmask)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'elev',id_elev_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'slope',id_slope_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_elev_in,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_slope_in,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_variable(ifid,id_elev_in,chunksizes=sblock_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!maximum number of input pixels in one output chunk (along one dimension)
!should be even multiple of the input file chunk size

pixchunkmax = sblock_in(1) * pixchunkmaxf

resratio(1) = xres_in / xres_out
resratio(2) = yres_in / yres_out

pixpix = 1. / resratio

out%minlon = in%minlon
out%maxlon = in%maxlon
out%minlat = in%minlat
out%maxlat = in%maxlat

out%countx = in%countx * resratio(1)
out%county = in%county * resratio(2)

!limit the block size for the input data to a maximum of 12,000 x 12,000 input pixels (about 1Gb)

sblock_in(1) = min(in%countx,pixchunkmax)
sblock_in(2) = min(in%county,pixchunkmax)

!the output chunk size will be set to the input block size

sblock_out(1:2) = sblock_in * resratio
!sblock_out(3)   = 1

nblock_out(1) = out%countx / sblock_out(1)
nblock_out(2) = out%county / sblock_out(2)

!---

! write(0,*)missing_i2,missing_sp
write(0,*)'total input pixels:     ',in%countx,in%county
write(0,*)'total output pixels:    ',out%countx,out%county
write(0,*)'input read block size:  ',sblock_in
write(0,*)'output block size:      ',sblock_out
write(0,*)'input pixels per opixel:',pixpix
write(0,*)'pixchunkmax:            ',pixchunkmax
write(0,*)'number of output blocks:',nblock_out

if (out%countx > maxchunklen) then
  schunk_out(1) = out%countx / (1 + out%countx / maxchunklen)
else
  schunk_out(1) = out%countx
end if

if (out%county > maxchunklen) then
  schunk_out(2) = out%county / (1 + out%county / maxchunklen)
else
  schunk_out(2) = out%county
end if

schunk_out(3) = 1

nchunk_out(1) = out%countx / schunk_out(1)
nchunk_out(2) = out%county / schunk_out(2)
nchunk_out(3) = nclasses

write(0,*)'output chunk size:      ',schunk_out
write(0,*)'num output chunks:      ',nchunk_out

!number of output blocks per output chunk

nblkchk = schunk_out(1:2) / sblock_out

write(0,*)'num outblocks per chunk:',nblkchk

write(0,*)'press return to continue'
read(*,*)

!---
!generate output file

call getarg(4,outfile)

call genoutfile(outfile,out,llim,schunk_out,ofid)

ncstat = nf90_inq_varid(ofid,'elev',id_elev_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slope',id_slope_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'elev_stdev',id_elev_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slope_stdev',id_slope_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'classfrac',id_classfrac)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---
!calculate median elevation and slope and standard deviations

allocate(lat(sblock_in(2)))
allocate(area(sblock_in(1),sblock_in(2)))
allocate(elev_in(sblock_in(1),sblock_in(2)))
allocate(slope_in(sblock_in(1),sblock_in(2)))
allocate(landmask(sblock_in(1),sblock_in(2)))

allocate(stats(sblock_out(1),sblock_out(2)))
allocate(classfrac(sblock_out(1),sblock_out(2),nclasses))

do y = 1,nblock_out(2)
  
  yoffset = (y-1) * sblock_out(2) * pixpix(2)
  
  do x = 1,nblock_out(1)

    write(status_line,'(a,i5,a,i5)')' working on row: ',y,' col: ',x
    call overprint(trim(status_line))
    
    xoffset = (x-1) * sblock_out(1) * pixpix(1)
    
    !output grid offset
    
    b = 1 + (y-1) * sblock_out(2)
    a = 1 + (x-1) * sblock_out(1)
    
    !set output structure to default values
    
    stats%elev_med  = missing_i2
    stats%elev_std  = missing_i2
    stats%slope_med = missing_sp
    stats%slope_std = missing_sp

    do k = 1,nclasses
      stats%classfrac(k) = missing_sp
    end do
    
    !---
    
    startx = in%startx + xoffset
    starty = in%starty + yoffset
    
    !---
    !read block of input data
    
!     write(0,*)'read'

    ncstat = nf90_get_var(ifid,yvarid,lat,start=[starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(ifid,id_elev_in,elev_in,start=[startx,starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(ifid,id_slope_in,slope_in,start=[startx,starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(lfid,id_landmask,landmask,start=[startx,starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)
      
    !---

!     write(0,*)'calc',size(slope_in),size(landmask)
    
    do j = 1,sblock_in(1)
      do i = 1,sblock_in(2)
        if (landmask(i,j) /= 1) then  !slope_in(i,j) == 0. .and. 

          elev_in(i,j)  = missing_i2
          slope_in(i,j) = missing_sp

	end if
      end do
    end do

    if (all(elev_in == missing_i2)) cycle

    !calculate input cell area

    do k = 1,sblock_in(2)
      area(:,k) = cellarea(lat(k),[xres_in,yres_in])
    end do
    
    !---
    !calculate stats for each output pixel in this input block
    
    ompchunk = min(8,sblock_out(2))

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,n,p,q,r,s,xpos_out,ypos_out)
    !$OMP DO SCHEDULE(DYNAMIC,ompchunk)
    
    do n = 1,sblock_out(2)
      do m = 1,sblock_out(1)
        
        p = 1 + pixpix(1) * (m-1)
        q = p + pixpix(1) - 1
        r = 1 + pixpix(2) * (n-1)
        s = r + pixpix(2) - 1
        
        if (all(elev_in(p:q,r:s) == missing_i2)) cycle

        call calcstats(area(p:q,r:s),elev_in(p:q,r:s),slope_in(p:q,r:s),stats(m,n))
        
      end do
    end do

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
  
    !---
    !write data if output buffer is full

    !write the output chunk
    !should try to optimize for fewer writes, but this is as far as I have got

    !write(0,*)'write'

    ncstat = nf90_put_var(ofid,id_elev_out,stats%elev_med,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_elev_std,stats%elev_std,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_slope_out,stats%slope_med,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_slope_std,stats%slope_std,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)
      
    do k = 1,nclasses
      classfrac(:,:,k) = stats%classfrac(k)
    end do

    ncstat = nf90_put_var(ofid,id_classfrac,classfrac,start=[a,b,1])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

  end do    !output chunks
end do

write(0,*)

!write the actual ranges

!write the output longitude and latitude and slope classes

allocate(lon_out(out%countx))
allocate(lat_out(out%county))

!write(0,*)out%minlon,out%minlat

xres_out = xres_out / 3600.  !convert to degrees from seconds
yres_out = yres_out / 3600.

hpx = [0.5 * xres_out, 0.5 * yres_out]

lon_out = [(out%minlon + real(i-1) * xres_out + hpx(1),i=1,out%countx)]
lat_out = [(out%minlat + real(i-1) * yres_out + hpx(2),i=1,out%county)]

ncstat = nf90_inq_varid(ofid,'lon',varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lon_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,lat_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slopeclass',varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,llim)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---
!close

ncstat = nf90_close(ifid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(lfid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

end program topostats
