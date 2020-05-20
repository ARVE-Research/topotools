program topostats

use parametersmod, only : i1,i2,i4,sp,dp,missing_sp
use netcdfmod,     only : ncstat,handle_err
use netcdf
use typesizes
use coordsmod,     only : index,parsecoords
use statsmod,      only : median,stdev
use utilitiesmod,  only : iminloc,cellarea
use outputmod,     only : genoutfile
use calcstatsmod,  only : nclasses_slope,nclasses_aspect,nclasses_cti,classbin,statvals,calcstats,calcstats_aspect
use calcfracmod,   only : calcfrac

implicit none

character(100) :: infile
character(100) :: outfile

integer(i4) :: ifid
integer(i4) :: ofid

integer(i4) :: dimid
integer(i4) :: xvarid
integer(i4) :: yvarid

integer(i4) :: id_elev_in
integer(i4) :: id_slope_in
integer(i4) :: id_aspect_in
integer(i4) :: id_cti_in
integer(i4) :: id_elev_out
integer(i4) :: id_elev_std
integer(i4) :: id_slope_out
integer(i4) :: id_slope_std
integer(i4) :: id_aspect_out
integer(i4) :: id_aspect_std
integer(i4) :: id_cti_out
integer(i4) :: id_cti_std
integer(i4) :: id_classfrac_slope
integer(i4) :: id_classfrac_aspect
integer(i4) :: id_classfrac_cti

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

real(sp),    allocatable, dimension(:,:) :: elev_in
real(sp),    allocatable, dimension(:,:) :: slope_in
real(sp),    allocatable, dimension(:,:) :: aspect_in
real(sp),    allocatable, dimension(:,:) :: cti_in

real(sp),    allocatable, dimension(:,:,:) :: classfrac_slope
real(sp),    allocatable, dimension(:,:,:) :: classfrac_aspect
real(sp),    allocatable, dimension(:,:,:) :: classfrac_cti

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

type(classbin) :: llim
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
!set up the bounds for the variable bin classes

call calcfrac(llim)

!--------

call getarg(1,coordstring)

call parsecoords(coordstring,in)

call getarg(2,outres)

read(outres,*)xres_out

yres_out = xres_out

!---

call getarg(3,infile)

ncstat = nf90_open(infile,nf90_nowrite,ifid)
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

ncstat = nf90_inq_varid(ifid,'elev',id_elev_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'slope',id_slope_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'aspect',id_aspect_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'cti',id_cti_in)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_elev_in,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_slope_in,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_aspect_in,'missing_value',missing_sp)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,id_cti_in,'missing_value',missing_sp)
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

! write(0,*)missing_sp,missing_sp
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
nchunk_out(3) = nclasses_slope

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

ncstat = nf90_inq_varid(ofid,'aspect',id_aspect_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'cti',id_cti_out)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'elev_stdev',id_elev_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slope_stdev',id_slope_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'aspect_stdev',id_aspect_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'cti_stdev',id_cti_std)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slope_classfrac',id_classfrac_slope)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'aspect_classfrac',id_classfrac_aspect)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'cti_classfrac',id_classfrac_cti)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---
!nclasses_slopeulate median elevation and slope and standard deviations

allocate(lat(sblock_in(2)))
allocate(area(sblock_in(1),sblock_in(2)))
allocate(elev_in(sblock_in(1),sblock_in(2)))
allocate(slope_in(sblock_in(1),sblock_in(2)))
allocate(aspect_in(sblock_in(1),sblock_in(2)))
allocate(cti_in(sblock_in(1),sblock_in(2)))

allocate(stats(sblock_out(1),sblock_out(2)))
allocate(classfrac_slope(sblock_out(1),sblock_out(2),nclasses_slope))
allocate(classfrac_aspect(sblock_out(1),sblock_out(2),nclasses_aspect))
allocate(classfrac_cti(sblock_out(1),sblock_out(2),nclasses_cti))

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

    stats%elev_med   = missing_sp
    stats%elev_std   = missing_sp
    stats%slope_med  = missing_sp
    stats%slope_std  = missing_sp
    stats%aspect_med = missing_sp
    stats%aspect_std = missing_sp
    stats%cti_med    = missing_sp
    stats%cti_std    = missing_sp

    do k = 1,nclasses_slope
      stats%classfrac_slope(k) = missing_sp
    end do

    do k = 1,nclasses_aspect
      stats%classfrac_aspect(k) = missing_sp
    end do

    do k = 1,nclasses_cti
      stats%classfrac_cti(k) = missing_sp
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

    ncstat = nf90_get_var(ifid,id_aspect_in,aspect_in,start=[startx,starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(ifid,id_cti_in,cti_in,start=[startx,starty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    !---

!     write(0,*)'nclasses_slope',size(slope_in)

    if (all(elev_in == missing_sp)) cycle

    !nclasses_slopeulate input cell area

    do k = 1,sblock_in(2)
      area(:,k) = cellarea(lat(k),[xres_in,yres_in])
    end do

    !---
    !nclasses_slopeulate stats for each output pixel in this input block

    ompchunk = min(8,sblock_out(2))

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,n,p,q,r,s,xpos_out,ypos_out)
    !$OMP DO SCHEDULE(DYNAMIC,ompchunk)

    do n = 1,sblock_out(2)
      do m = 1,sblock_out(1)

        p = 1 + pixpix(1) * (m-1)
        q = p + pixpix(1) - 1
        r = 1 + pixpix(2) * (n-1)
        s = r + pixpix(2) - 1

        if (all(elev_in(p:q,r:s) == missing_sp)) cycle

          if(all(aspect_in(p:q,r:s) == missing_sp)) then

            call calcstats(area(p:q,r:s),elev_in(p:q,r:s),slope_in(p:q,r:s),cti_in(p:q,r:s),stats(m,n),llim) !Avoid floating invalid in std_dev calculation
            !Aspect remain as "missing" in this subroutine
          else

            call calcstats_aspect(area(p:q,r:s),elev_in(p:q,r:s),slope_in(p:q,r:s),aspect_in(p:q,r:s),cti_in(p:q,r:s),stats(m,n),llim)

          end if

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

    ncstat = nf90_put_var(ofid,id_aspect_out,stats%aspect_med,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_aspect_std,stats%aspect_std,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_cti_out,stats%cti_med,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_cti_std,stats%cti_std,start=[a,b])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    do k = 1,nclasses_slope
      classfrac_slope(:,:,k) = stats%classfrac_slope(k)
    end do

    do k = 1,nclasses_aspect
      classfrac_aspect(:,:,k) = stats%classfrac_aspect(k)
    end do

    do k = 1,nclasses_cti
      classfrac_cti(:,:,k) = stats%classfrac_cti(k)
    end do

    ncstat = nf90_put_var(ofid,id_classfrac_slope,classfrac_slope,start=[a,b,1])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_classfrac_aspect,classfrac_aspect,start=[a,b,1])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_classfrac_cti,classfrac_cti,start=[a,b,1])
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

ncstat = nf90_put_var(ofid,varid,llim%slope_bin)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'aspectclass',varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,llim%aspect_bin)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'cticlass',varid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,llim%cti_bin)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---
!close

ncstat = nf90_close(ifid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

print *, llim

end program topostats
