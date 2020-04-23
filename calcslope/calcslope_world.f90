program calcslope

!calculate slopes from a 3" DEM using the method of Song & Shan (Photogrammetric Eng. & Remote Sens., 75(3), 281-290)

use parametersmod, only : i2,i4,sp,dp
use netcdfmod,     only : ncstat,handle_err,makeoutfile
use distmod,       only : sphdist
use coordsmod,     only : index,parsecoords
use utilitiesmod,  only : iminloc
use outputmod,     only : genoutfile,putlonlat
use netcdf

implicit none

character(100) :: infile
character(100) :: outfile

integer(i4) :: ifid
integer(i4) :: ofid

integer :: dimid
integer :: xvarid
integer :: yvarid
integer :: zvarid

integer :: id_olon
integer :: id_olat
integer :: id_dem
integer :: id_slope

real(dp),    allocatable, dimension(:)   :: all_lon
real(dp),    allocatable, dimension(:)   :: lon
real(dp),    allocatable, dimension(:)   :: lat
integer(i2), allocatable, dimension(:,:) :: dem
real(sp),    allocatable, dimension(:,:) :: slope

real(dp),    dimension(8) :: snbr
real(sp),    dimension(8) :: dist
integer(i2), dimension(8) :: elev

integer(i4), dimension(8,2) :: idx

real(dp), dimension(2)   :: ll0
real(dp), dimension(8,2) :: llnbr

integer :: chunk
integer :: i,j
integer :: x,y
integer :: x0,y0
integer :: n
integer :: xlen
integer :: ylen
integer :: l

real(sp), dimension(8) :: dz

integer(i2) :: missing_i2

real(sp),    parameter :: missing_sp = -9999.

logical, dimension(8) :: nbr

real(dp) :: dzdxc
real(dp) :: dzdyc
real(dp) :: dzdxd
real(dp) :: dzdyd

real(dp) :: Sd8
real(dp) :: Sfd

real(dp) :: w

real(dp), dimension(4) :: ddz

integer, dimension(2) :: cs

integer :: nblkx
integer :: nblky

integer :: blklenx
integer :: blkleny

integer :: srtx
integer :: srty
integer :: srtx_o
integer :: srty_o

character(80) :: coordstring

type(index) :: id

real(dp) :: xres
real(dp) :: yres

real(dp), dimension(2) :: lonrange
real(dp), dimension(2) :: latrange

character(80) :: status_line

logical, dimension(8) :: nsf

!--------

idx(1,:) = [-1, 1]
idx(2,:) = [ 0, 1]
idx(3,:) = [ 1, 1]
idx(4,:) = [-1, 0]
idx(5,:) = [ 1, 0]
idx(6,:) = [-1,-1]
idx(7,:) = [ 0,-1]
idx(8,:) = [ 1,-1]

!---

open(99,file='highslope.log',status='unknown')

!define an input neighborhood that is n9 around the 0.5 degree block of interest

call getarg(1,coordstring)

call parsecoords(coordstring,id)

!---

call getarg(2,infile)

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

!write(0,*)'infile size: ',xlen,ylen

!---

ncstat = nf90_inq_varid(ifid,'lon',xvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)
  
ncstat = nf90_get_att(ifid,xvarid,'valid_range',lonrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ifid,'lat',yvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,yvarid,'valid_range',latrange)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

xres = (lonrange(2) - lonrange(1)) / real(xlen)
yres = (latrange(2) - latrange(1)) / real(ylen)

!write(0,*)xres,yres

!---
!get indices then repurpose lon and lat

allocate(lon(xlen))
allocate(lat(ylen))

ncstat = nf90_get_var(ifid,xvarid,lon)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(ifid,yvarid,lat)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

id%startx = iminloc(abs(id%minlon - lon))  !start-x
id%starty = iminloc(abs(id%minlat - lat))  !start-y

id%endx = iminloc(abs(id%maxlon - lon))    !end-x
id%endy = iminloc(abs(id%maxlat - lat))    !end-y

!adjust the nodes to select if the nearest node automatically selected is outside the area of interest

if (lon(id%startx) < id%minlon) id%startx = id%startx + 1
if (lon(id%endx)   > id%maxlon) id%endx   = id%endx   - 1

if (lat(id%starty) < id%minlat) id%starty = id%starty + 1
if (lat(id%endy)   > id%maxlat) id%endy   = id%endy   - 1

id%countx = 1 + id%endx - id%startx  !add one to include the start and end points
id%county = 1 + id%endy - id%starty

write(0,*)minval(lon),maxval(lon)
write(0,*)minval(lat),maxval(lat)

write(0,*)id%startx,id%endx,id%countx
write(0,*)id%starty,id%endy,id%county

!read(*,*)

allocate(all_lon(xlen))

all_lon = lon

deallocate(lon)
deallocate(lat)

!---

ncstat = nf90_inq_varid(ifid,'elv',zvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,zvarid,'missing_value',missing_i2)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_variable(ifid,zvarid,chunksizes=cs)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

write(0,*)'chunksizes:',cs

!---

write(0,*)id%startx,id%countx
write(0,*)id%starty,id%county

!---
!generate output file

call getarg(3,outfile)

call genoutfile(outfile,id,cs,ofid)

ncstat = nf90_inq_varid(ofid,'lon',id_olon)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',id_olat)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'elev',id_dem)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'slope',id_slope)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

!---
!divide desired job up into blocks depending on the input chunksize

if (mod(id%countx,cs(1)) == 0) then
  nblkx = id%countx / cs(1)
else
  nblkx = id%countx / cs(1) + 1
end if
  
if (mod(id%county,cs(2)) == 0) then
  nblky = id%county / cs(2)
else
  nblky = id%county / cs(2) + 1
end if

!write(0,'(i10,a,2i6)')nblkx*nblky,' blocks to calculate ',nblkx,nblky

!define the computational block size

blklenx = min(id%countx,cs(1))  ! + 2  !buffer one column right and left
blkleny = min(id%county,cs(2))  ! + 2  !buffer one row above and below

if (blkleny > 8) then
  chunk = blkleny / 8
else
  chunk = blkleny
end if

allocate(lon(0:blklenx+1))
allocate(lat(0:blkleny+1))
allocate(dem(0:blklenx+1,0:blkleny+1))
allocate(slope(blklenx,blkleny))

!l = 0

do j = 1,nblky
  
  srty = id%starty + (j-1) * blkleny
 
  srty_o = 1 + (j-1) * blkleny

  write(status_line,'(a,i5,a,i5)')' working on row: ',j,' of ',nblky
  call overprint(trim(status_line))
  
  !--------------------------------------
  ! get latitude
  
  if (srty == 1) then  !at bottom edge of global grid
    
    !ncstat = nf90_get_var(ifid,yvarid,lat(0),start=[srty+1])
    !if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(ifid,yvarid,lat(1:),start=[srty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    call putlonlat(ofid,id_olon,id_olat,all_lon(id%startx:id%endx),lat(1:blkleny),1,srty_o)

    cycle !skip it for now

  else if (srty + blkleny + 1 > ylen) then  !at top edge of global grid

    !ncstat = nf90_get_var(ifid,yvarid,lat(blkleny+1),start=[srty-1])
    !if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_get_var(ifid,yvarid,lat(:blkleny),start=[srty-1])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    call putlonlat(ofid,id_olon,id_olat,all_lon(id%startx:id%endx),lat(1:blkleny),1,srty_o)

    cycle !skip it for now (no land anyway)

  else  !normal situation
    
    ncstat = nf90_get_var(ifid,yvarid,lat,start=[srty-1])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

  end if
  
  !--------------------------------------
  ! get longitude
  
  do i = 1,nblkx

    srtx_o = 1 + (i-1) * blklenx

    !l = l + 1

    srtx = id%startx + (i-1) * blklenx

    !write(status_line,'(a,3i5,4i8)')' working on: ',l,i,j,srtx,xlen,srty,ylen
    !call overprint(trim(status_line))
    
    if (srtx == 1) then  !at left edge of global grid
      
      !write(0,*)'case left edge'
      
      ncstat = nf90_get_var(ifid,xvarid,lon(0),start=[xlen])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)
        
      ncstat = nf90_get_var(ifid,xvarid,lon(1:),start=[srtx],count=[blklenx+1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(0,:),start=[xlen,srty-1],count=[1,blkleny+2])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(1:,:),start=[srtx,srty-1],count=[blklenx+1,blkleny+2])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

    else if (srtx + blklenx + 1 > xlen) then  !at right edge of global grid

      !write(0,*)'case right edge'
      
      ncstat = nf90_get_var(ifid,xvarid,lon(blklenx+1),start=[1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,xvarid,lon(:blklenx),start=[srtx-1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(blklenx+1,:),start=[1,srty-1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(:blklenx,:),start=[srtx-1,srty-1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

    else  !normal situation
    
      ncstat = nf90_get_var(ifid,xvarid,lon,start=[srtx-1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem,start=[srtx-1,srty-1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

    end if
        
    where (dem < -1000) dem = missing_i2

    if (all(dem == missing_i2)) then

      call putlonlat(ofid,id_olon,id_olat,lon(1:blklenx),lat(1:blkleny),srtx_o,srty_o)

      cycle

    end if
    
    slope = missing_sp

    !calculate slopes in the block that is in the center of the superblock
    
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(blklenx,blkleny,missing_i2,dem,idx,slope,lon,lat,chunk)
    !$OMP DO SCHEDULE(DYNAMIC,chunk)
    
    do y0 = 1,blkleny
      do x0 = 1,blklenx
        
        if (dem(x0,y0) == missing_i2) cycle
        
        ll0 = [lon(x0),lat(y0)]
        
        nbr = .true.
        nsf = .false.
        
        do n = 1,8
          
          x = x0 + idx(n,1)
          y = y0 + idx(n,2)
            
          if (dem(x,y) == dem(x0,y0)) nsf(n) = .true.

        end do
        
        if (all(nsf)) then  !all cells have the same elevation as center pixel - slope zero
          
          slope(x0,y0) = 0.
          
          cycle
       
        else
       
          do n = 1,8
          
            x = x0 + idx(n,1)
            y = y0 + idx(n,2)
            
            if (dem(x,y) == missing_i2) then                    !or otherwise missing neighbor

              nbr(n) = .false.
           
            else
          
              elev(n) = dem(x,y)

              dz(n) = real(elev(n) - dem(x0,y0))

              llnbr(n,:) = [lon(x),lat(y)]

              !calculate distance to neighbors

              dist(n) = sphdist(ll0,llnbr(n,:))

              !calculate slope to neighbors

              snbr(n) = dz(n) / dist(n)

            end if
          end do
        end if

        Sd8 = max(maxval(snbr,mask=nbr),0.d0)

        if (all(nbr)) then

          dzdxc = real(elev(2) - elev(6)) / (dist(2) + dist(6))
          dzdyc = real(elev(8) - elev(4)) / (dist(8) + dist(4))
          dzdxd = real(elev(1) - elev(5)) / (dist(1) + dist(5))
          dzdyd = real(elev(7) - elev(3)) / (dist(7) + dist(3))

          ddz = [dzdxc,dzdyc,dzdxd,dzdyd]

          if (sum(ddz) == 0.d0) then

            slope(x0,y0) = Sd8        

          else

            Sfd = sqrt(dzdxc**2 + dzdyc**2)

            w = 0.5d0 * (sqrt(dzdxd**2 + dzdyd**2) - sqrt(dzdxc**2 + dzdyc**2))**2 / (dzdxd**2 + dzdyd**2 + dzdxc**2 + dzdyc**2)

            slope(x0,y0) = w * Sd8 + (1.d0 - w) * Sfd

          end if

        else

          slope(x0,y0) = Sd8

        end if
        
        if (slope(x0,y0) > 4.) then

          write(99,'(2f15.10,2i8,f9.4,9i5)')ll0,x0+srtx-1,y0+srty-1,slope(x0,y0),dem(x0,y0),elev
          
        end if
        
      end do    !block x
    end do      !block y
    
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !write out block
    
    !write(0,*)'writing'
    
    call putlonlat(ofid,id_olon,id_olat,lon(1:blklenx),lat(1:blkleny),srtx_o,srty_o)

    ncstat = nf90_put_var(ofid,id_dem,dem(1:blklenx,1:blkleny),start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)    

    ncstat = nf90_put_var(ofid,id_slope,slope,start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)    

    !write(0,*)'done'

  end do
end do

write(0,*)

close(99)

ncstat = nf90_close(ifid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

end program calcslope
