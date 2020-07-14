program calcslope

!calculate slopes from a 3" DEM using the method of Song & Shan (Photogrammetric Eng. & Remote Sens., 75(3), 281-290)

use parametersmod, only : i4,sp,dp,d2r
use netcdfmod,     only : ncstat,handle_err
use distmod,       only : sphdist
use coordsmod,     only : index,parsecoords
use utilitiesmod,  only : iminloc,imaxloc
use outputmod,     only : genoutfile,putlonlat
use netcdf

implicit none

real(sp), parameter :: minslope = 0.001  ! minimum slope for calculating aspect and CTI: 0.1% or 1/1000 m m-1 after Marthews et al (2015)

character(100) :: infile
character(100) :: outfile

integer(i4) :: ifid
integer(i4) :: ofid

integer :: dimid
integer :: xvarid
integer :: yvarid
integer :: zvarid
integer :: favarid

integer :: id_olon
integer :: id_olat
integer :: id_dem
integer :: id_slope
integer :: id_aspect
integer :: id_cti
integer :: id_curvature_profile
integer :: id_curvature_plan
integer :: id_curvature_tangent

real(dp), allocatable, dimension(:)   :: all_lon
real(dp), allocatable, dimension(:)   :: lon
real(dp), allocatable, dimension(:)   :: lat
real(sp), allocatable, dimension(:,:) :: dem
real(sp), allocatable, dimension(:,:) :: flowacc
real(sp), allocatable, dimension(:,:) :: slope
real(sp), allocatable, dimension(:,:) :: aspect
real(sp), allocatable, dimension(:,:) :: cti
real(sp), allocatable, dimension(:,:) :: curvature_profile
real(sp), allocatable, dimension(:,:) :: curvature_plan
real(sp), allocatable, dimension(:,:) :: curvature_tangent

real(dp), dimension(8) :: snbr
real(sp), dimension(8) :: dist
real(sp), dimension(8) :: elev

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
integer :: counting = 0

real(sp), dimension(8) :: dz

real(sp), dimension(2) :: range_elv
real(sp), dimension(2) :: range_flowacc

real(sp) :: missing

logical, dimension(8) :: nbr
logical, dimension(8) :: ucell

real(dp) :: dzdxc
real(dp) :: dzdyc
real(dp) :: dzdxd
real(dp) :: dzdyd
real(dp) :: d2zdx2
real(dp) :: d2zdy2
real(dp) :: d2zdxdy
real(dp) :: p
real(dp) :: q

real(dp) :: rc, qc, pc, tc, sc, w_dis ! Curvature calculation parameters

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

integer, dimension(1) :: dir

real(sp) :: Ad8
real(sp) :: Afd

real(sp) :: cellsize_ns
real(sp) :: cellsize_ew
real(sp) :: pixlen
real(sp) :: c_all
real(sp) :: c_in
real(sp) :: cout
real(sp) :: cdiag
real(sp) :: msu
real(sp) :: dfltsink
real(sp) :: aspec
real(sp) :: lnpixlen

integer :: flowdir

!--------
! neighbor node numbering convention used in Wilson and Gallant (Terrain Analysis, 2000) and in other papers

idx(1,:) = [ 1, 1]    ! northeast
idx(2,:) = [ 1, 0]    ! east
idx(3,:) = [ 1,-1]    ! southeast
idx(4,:) = [ 0,-1]    ! south
idx(5,:) = [-1,-1]    ! southwest
idx(6,:) = [-1, 0]    ! west
idx(7,:) = [-1, 1]    ! northwest
idx(8,:) = [ 0, 1]    ! north

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

write(0,*)'infile size: ',xlen,ylen

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

if (id%startx /= id%endx) then

  if (lon(id%startx) < id%minlon) id%startx = id%startx + 1
  if (lon(id%endx)   > id%maxlon) id%endx   = id%endx   - 1

end if

if (id%starty /= id%endy) then

  if (lat(id%starty) < id%minlat) id%starty = id%starty + 1
  if (lat(id%endy)   > id%maxlat) id%endy   = id%endy   - 1

end if

id%countx = 1 + id%endx - id%startx  !add one to include the start and end points
id%county = 1 + id%endy - id%starty

write(0,*)'lon range: ',minval(lon),maxval(lon)
write(0,*)'lat_range: ',minval(lat),maxval(lat)

write(0,*)'x range: ',id%startx,id%endx,id%countx
write(0,*)'y range: ',id%starty,id%endy,id%county

!read(*,*)

allocate(all_lon(xlen))

all_lon = lon

deallocate(lon)
deallocate(lat)

!---

ncstat = nf90_inq_varid(ifid,'elv',zvarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,zvarid,'missing_value',missing)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,zvarid,'valid_range',range_elv)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_variable(ifid,zvarid,chunksizes=cs)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

cs(1) = min(cs(1),id%countx)
cs(2) = min(cs(2),id%county)

write(0,*)'chunksizes:',cs
write(0,*)'range elv: ',range_elv

ncstat = nf90_inq_varid(ifid,'flowacc',favarid)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ifid,favarid,'valid_range',range_flowacc)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

write(0,*)'range flowacc: ',range_flowacc

!---

! write(0,*)id%startx,id%countx
! write(0,*)id%starty,id%county

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

ncstat = nf90_inq_varid(ofid,'aspect',id_aspect)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'cti',id_cti)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'curvature_profile',id_curvature_profile)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'curvature_plan',id_curvature_plan)
if (ncstat/=nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'curvature_tangent',id_curvature_tangent)
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
allocate(flowacc(blklenx,blkleny))
allocate(slope(blklenx,blkleny))
allocate(aspect(blklenx,blkleny))
allocate(cti(blklenx,blkleny))
allocate(curvature_profile(blklenx, blkleny))
allocate(curvature_plan(blklenx, blkleny))
allocate(curvature_tangent(blklenx, blkleny))

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

    cycle !skip it for now (no Antarctica data anyway)

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

    if (srtx == 1) then  !at left edge of global grid, first column copied from last column in input

      !write(0,*)'case left edge'

      ncstat = nf90_get_var(ifid,xvarid,lon(0),start=[xlen])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,xvarid,lon(1:),start=[srtx],count=[blklenx+1])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(0,:),start=[xlen,srty-1],count=[1,blkleny+2])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

      ncstat = nf90_get_var(ifid,zvarid,dem(1:,:),start=[srtx,srty-1],count=[blklenx+1,blkleny+2])
      if (ncstat/=nf90_noerr) call handle_err(ncstat)

    else if (srtx + blklenx + 1 > xlen) then  !at right edge of global grid, last column copied from first column in input

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

    where (dem < range_elv(1)) dem = missing

    if (all(dem == missing)) then

      call putlonlat(ofid,id_olon,id_olat,lon(1:blklenx),lat(1:blkleny),srtx_o,srty_o)

      cycle

    end if

    ! get the flow accumulation (for calculating CTI)

    ncstat = nf90_get_var(ifid,favarid,flowacc,start=[srtx,srty])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    slope  = missing
    aspect = missing
    cti = missing
    curvature_profile = missing
    curvature_plan = missing
    curvature_tangent = missing

    !calculate slopes in the block that is in the center of the superblock

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(blklenx,blkleny,missing,dem,idx,slope,aspect,flowacc,cti,curvature_profile,curvature_plan,curvature_tangent,lon,lat,chunk)
    !$OMP DO SCHEDULE(DYNAMIC,chunk)

    do y0 = 1,blkleny
      do x0 = 1,blklenx

        if (dem(x0,y0) == missing) cycle

        ll0 = [lon(x0),lat(y0)]

        nbr = .true.
        nsf = .false.

        snbr = 0.

        do n = 1,8

          x = x0 + idx(n,1)
          y = y0 + idx(n,2)

          if (dem(x,y) == dem(x0,y0)) nsf(n) = .true.

        end do

!         if (all(nsf)) then  !all cells have the same elevation as center pixel - slope zero
!
!           slope(x0,y0) = 0.
!
!           cycle
!
!         else

          do n = 1,8

            x = x0 + idx(n,1)
            y = y0 + idx(n,2)

            if (dem(x,y) == missing) then    ! or otherwise missing neighbor

              nbr(n) = .false.

              llnbr(n,:) = [lon(x),lat(y)]

              !calculate distance to neighbors

              dist(n) = sphdist(ll0,llnbr(n,:)) ! still need to assign dist(n) to avoid pixlen == 0. in CTI calculation

            else

              elev(n) = dem(x,y)

              dz(n) = dem(x0,y0) - elev(n)     ! uphill slopes will be negative; convention as in Wilson and Gallant (2000)

              llnbr(n,:) = [lon(x),lat(y)]

              !calculate distance to neighbors

              dist(n) = sphdist(ll0,llnbr(n,:))

              !calculate slope to neighbors

              snbr(n) = dz(n) / dist(n)

            end if

          end do

          if (dist(2) == 0.) then
            print *, 'x0,y0: ', x0, y0
            print *, 'lon0, lat0: ', lon(x0), lat(y0)
            print *, 'x, y : ', x0 + idx(2,1), y0 + idx(2,2)
            print *, 'lon1, lat1: ', lon(x0 + idx(2,1)), lat(y0 + idx(2,2))
            print *, 'dist: ', dist
            !stop "dist is 0"
          end if

!         end if

        ! D8 slope calculation

        Sd8 = max(maxval(snbr,mask=nbr),0._dp)

        ! D8 aspect

        dir = maxloc(snbr,mask=nbr)

        Ad8 = 45. * real(dir(1))

        if (all(nbr)) then

          ! all neighbor gridcells are present so we can calculate finite differences

          dzdxc = real(elev(2) - elev(6)) / (dist(2) + dist(6))  ! east-west   (dz/dx cardinal)
          dzdyc = real(elev(8) - elev(4)) / (dist(8) + dist(4))  ! north-south (dz/dy cardinal)
          dzdxd = real(elev(1) - elev(5)) / (dist(1) + dist(5))  ! northeast-southwest (dz/dx diagonal)
          dzdyd = real(elev(7) - elev(3)) / (dist(7) + dist(3))  ! northwest-southeast (dz/dy diagonal)

          ddz = [dzdxc,dzdyc,dzdxd,dzdyd]

          ! 1st method: Curvature parameters from Evan's method (Florinsky, 1988) http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.86.9066&rep=rep1&type=pdf

          w_dis = sum(dist) / 8. ! Average cell length to calculate parameters: rc, tc, sc, pc, qc

          rc = (elev(7) + elev(1) + elev(6) + elev(2) + elev(5) + elev(3) - 2. * (elev(7) + elev(7) + dem(x0,y0))) / (3. * (w_dis ** 2))
          tc = (elev(7) + elev(8) + elev(1) + elev(5) + elev(4) + elev(3) - 2. * (elev(6) + elev(2) + dem(x0,y0))) / (3. * (w_dis ** 2))
          sc = (elev(1) + elev(5) - elev(7) - elev(3)) / (4. * (w_dis ** 2))
          pc = (elev(1) + elev(2) + elev(3) - elev(7) - elev(6) - elev(5)) / (6. * w_dis)
          qc = (elev(7) + elev(8) + elev(1) - elev(5) - elev(4) - elev(3)) / (6. * w_dis)

          ! Formulea for for profile and plan of curvature
          ! Hengi & Reuter (2009) Geomorphmetry: Concepts, Software and Applications (P.151)

          curvature_profile(x0,y0) = -((pc ** 2) * rc + 2 * pc * qc * rc * sc + (qc ** 2) * tc) / (((pc ** 2) + (qc ** 2)) * ((1 + (pc ** 2) + (qc ** 2)) ** 1.5))
          curvature_plan(x0,y0) = -((qc ** 2) * rc - 2 * pc * qc * sc + (pc ** 2) * tc) / (((1 + (pc ** 2) + (qc ** 2)) ** 1.5))
          curvature_tangent(x0,y0) = -((qc ** 2) * rc - 2 * pc * qc * sc + (pc ** 2) * tc) / (((pc ** 2) + (qc ** 2)) * ((1 + (pc ** 2) + (qc ** 2)) ** 0.5))

          ! Convert to radian per 100m
          curvature_profile(x0,y0) = curvature_profile(x0,y0) * 100.
          curvature_plan(x0,y0) = curvature_plan(x0,y0) * 100.
          curvature_tangent(x0,y0) = curvature_tangent(x0,y0) * 100.

          ! ---

          ! 2nd method: curvature from second derivative calculations from Wilson & Gallant (not working)

          ! second derivatives for curvature calculations
          !d2zdx2 = real((elev(2) - 2. * elev(9) + elev(6)) / (((dist(2) + dist(6)) / 2.) ** 2))
          !d2zdy2 = real((elev(8) - 2. * elev(9) + elev(4)) / (((dist(8) + dist(4)) / 2.) ** 2))
          !d2zdxdy = real((elev(1) + elev(5) - elev(7) - elev(3)) / ((((dist(2) + dist(4) + dist(6) + dist(8)) / 4.) ** 2) * 4.))

          !p = (dzdxc ** 2) + (dzdyc ** 2)       !TAPES-G formulea from Wilson & Gallant p.57
          !q = p + 1._dp

          ! Convert to radian per 100m
          !curvature_profile(x0,y0) = curvature_profile(x0,y0) * 100.
          !curvature_plan(x0,y0) = curvature_tangent(x0,y0) * 100.


          !curvature_profile(x0,y0) = ((d2zdx2 * (dzdxc ** 2)) + (2. * d2zdxdy * dzdxc * dzdyc) + (d2zdy2 * (dzdyc ** 2))) / (p * (q ** 1.5)) !Profile curvature
          !curvature_plan(x0,y0) = ((d2zdx2 * (dzdyc ** 2)) - (2. * d2zdxdy * dzdxc * dzdyc) + (d2zdy2 * (dzdxc ** 2))) / (p ** 1.5) !Plan curvature
          !curvature_tangent(x0,y0) = ((d2zdx2 * (dzdyc ** 2)) - (2. * d2zdxdy * dzdxc * dzdyc) + (d2zdy2 * (dzdxc ** 2))) / (p * (q ** 0.5)) !Tangent curvature
          !curvature_total(x0,y0) = ((d2zdx2 ** 2) + (2. * (d2zdxdy ** 2)) + (d2zdy2 ** 2)) ! Formula for total curvature Wilson & Gallant p.57

          ! --

          if (sum(ddz) == 0._dp) then

            slope(x0,y0) = Sd8 ! Sd8 = 0 if land is flat, which is the case for inland waterbodies

            if(slope(x0,y0) == 0.) then

              curvature_profile(x0,y0) = missing
              curvature_plan(x0,y0) = missing
              curvature_tangent(x0,y0) = missing

            end if

            if (slope(x0,y0) > minslope) then

              aspect(x0,y0) = Ad8  ! Aspect will remain "missing" if slope(x0,y0) < minslope,

            end if

          else

            ! finite difference slope calculation

            Sfd = sqrt(dzdxc**2 + dzdyc**2)

            ! Song and Shan (2009) weight factor

            w = 0.5_dp * (sqrt(dzdxd**2 + dzdyd**2) - sqrt(dzdxc**2 + dzdyc**2))**2 / (dzdxd**2 + dzdyd**2 + dzdxc**2 + dzdyc**2)

            ! final calculation of slope (Song and Shan, 2009)

            slope(x0,y0) = w * Sd8 + (1._dp - w) * Sfd

            !finite difference aspect

            if (slope(x0,y0) > minslope) then

              if (dzdxc == 0.) then !Avoid returning NaN

                aspect(x0,y0) = Ad8

              else

                Afd = 180. - (atan(dzdyc / dzdxc)) / d2r + 90. * (dzdxc / abs(dzdxc))

                aspect(x0,y0) = w * Ad8 + (1. - w) * Afd

              end if

            end if

          end if


        else  ! missing some neighbor gridcells

          slope(x0,y0) = Sd8

          aspect(x0,y0) = Ad8

          curvature_profile(x0,y0) = missing
          curvature_plan(x0,y0) = missing
          curvature_tangent(x0,y0) = missing

        end if

        ! Avoid unphysically large values for curvature when slope is extremely small
        if (slope(x0,y0) < minslope) then

          curvature_profile(x0,y0) = missing
          curvature_plan(x0,y0) = missing
          curvature_tangent(x0,y0) = missing

        end if

        ! topographic index calculations after Marthews et al., 2015 (Figure A1)

        ! part 1, contour length

        cellsize_ew = dist(2)
        cellsize_ns = dist(8)

        pixlen = 0.5 * (cellsize_ew + cellsize_ns)  ! mean pixel side length

                  if (pixlen == 0.) then
                    print *, x0, y0
                    print *, 'dist:', dist
                    print *, "infinity"
                    !stop
                  end if

        lnpixlen = log(pixlen)

        cdiag = 0.354 * pixlen  ! this is an estimate, would not work for very large pixels (several degrees)

        c_all = cellsize_ew + cellsize_ns + 4. * cdiag

        if (all(snbr <= 0.)) then  ! the cell is a sink or the whole neighborhood is flat (or if its water surface, then all(snbr == 0))

          cout = 0.

          msu = max(sum(abs(snbr)) / 8.,minslope) ! mean slope across the non-outflow contour

          if (msu > 0. .and. flowacc(x0,y0) > 0.) then

            dfltsink = log(flowacc(x0,y0) * 1.e6 / (2. * pixlen * msu))

            cti(x0,y0) = max(dfltsink - lnpixlen,0.)

            !cti(x0,y0) = min(cti(x0,y0), 25.004) !this is working...quick fix but need to check pixlen (Infinity = -log 0)

          else

           cti(x0,y0) = missing

          end if

          !print*, 'cti:',cti(x0,y0)

        else  ! normal case

          flowdir = imaxloc(snbr)

          select case(flowdir)
          case(1,3,5,7)  ! diagonal
            cout = cdiag
          case(2,6)      ! east or west
            cout = 0.5 * cellsize_ew
          case(4,8)      ! north or south
            cout = 0.5 * cellsize_ns
          end select

          c_in = c_all - cout


          if (flowacc(x0,y0) > 0.) then

            ! part 2, specific catchment area

            aspec = flowacc(x0,y0) * 1.e6 / cout

            ! part 3, slopes

            ucell = .false.

            where (nbr) ucell = .true.

            ucell(flowdir) = .false.

            msu = max(sum(abs(snbr),mask=ucell), minslope) / count(ucell)  ! mean slope across the non-outflow contour

            ! part 4, topographic index

            cti(x0,y0) = max(log(aspec / Sd8) - lnpixlen,0.) ! max(Sd8, minslope) needs checking, in case Sd8 == 0  but its NOT a dfltsink???

            !cti(x0,y0) = min(cti(x0,y0), 25.)

          else

            cti(x0,y0) = missing

          end if

          !print *, cti(x0,y0)

        end if

        if (slope(x0,y0) > 4.) then

          write(99,'(2f15.10,2i8,f9.4,9i5)')ll0,x0+srtx-1,y0+srty-1,slope(x0,y0),dem(x0,y0),elev

        end if

        ! write(0,*)'D8 slope',Sd8
        ! write(0,*)'FD slope',Sfd
        ! write(0,*)'D8 aspect',Ad8
        ! write(0,*)'FD aspect',Afd
        ! write(0,*)'final slope ',slope(x0,y0)
        ! write(0,*)'final aspect',aspect(x0,y0)

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

    ncstat = nf90_put_var(ofid,id_aspect,aspect,start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_cti,cti,start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_curvature_profile,curvature_profile,start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_curvature_plan,curvature_plan,start=[srtx_o,srty_o],count=[blklenx,blkleny])
    if (ncstat/=nf90_noerr) call handle_err(ncstat)

    ncstat = nf90_put_var(ofid,id_curvature_tangent,curvature_tangent,start=[srtx_o,srty_o],count=[blklenx,blkleny])
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

write(0,*) ' '
write(0,*) '   Done!   '
write(0,*) ' '

end program calcslope
