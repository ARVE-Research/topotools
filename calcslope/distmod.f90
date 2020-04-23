module distmod

!calculate the distance between two points in space, vector and scalar versions

use parametersmod, only : sp,dp,pi,d2r

implicit none

public :: sphdist
public :: cardist
public :: area

interface sphdist
  module procedure sphdist_scalar,sphdist_vector
end interface sphdist

real(dp), parameter :: earthrad = 6378.137d0 !km, WGS-84 spherical approximation

contains

!-------------------------------------------------------------------------------------------------

function sphdist_scalar(p0,p1)
  
!spherical distance calculation; returns distance between lat-lon points in meters
!spherical approximation form of the Vincenty formula
!coded from formula in Wikipedia article "Great-circle distance"

!scalar version

use parametersmod

implicit none

!arguments

real(dp), dimension(:), intent(in) :: p0  !lon-lat of pos 0
real(dp), dimension(:), intent(in) :: p1  !lon-lat pair of points to calculate distance to

real(sp) :: sphdist_scalar

!parameters

real(dp), parameter :: rm = earthrad * 1000.d0  !convert to meters

!local variables

real(dp) :: lon0
real(dp) :: lat0

real(dp) :: lon1
real(dp) :: lat1
real(dp) :: dlon
real(dp) :: num
real(dp) :: den
real(dp) :: dang

!---

lon0 = d2r * p0(1)
lon1 = d2r * p1(1)
lat0 = d2r * p0(2)
lat1 = d2r * p1(2)

dlon = lon1 - lon0

num = sqrt((cos(lat1) * sin(dlon))**2 + (cos(lat0) * sin(lat1) - sin(lat0) * cos(lat1) * cos(dlon))**2)
den = sin(lat0) * sin(lat1) + cos(lat0) * cos(lat1) * cos(dlon)

dang = atan2(num,den)

sphdist_scalar = rm * dang

end function sphdist_scalar

!-------------------------------------------------------------------------------------------------

function sphdist_vector(p0,p1)

!spherical distance calculation; returns distance between lat-lon points in meters
!spherical approximation form of the Vincenty formula
!coded from formula in Wikipedia article "Great-circle distance"

!this version operates on a vector

use parametersmod

implicit none

!arguments

real(dp), dimension(:),   intent(in) :: p0  !lon-lat of pos 0
real(dp), dimension(:,:), intent(in) :: p1  !lon-lat vector of points to calculate distance to

real(sp), dimension(size(p1,dim=1)) :: sphdist_vector

!parameters

real(dp), parameter :: rm     = earthrad * 1000.d0  !convert to meters

!local variables


real(dp) :: lon0
real(dp) :: lat0

real(dp), dimension(size(p1,dim=1)) :: lon1
real(dp), dimension(size(p1,dim=1)) :: lat1
real(dp), dimension(size(p1,dim=1)) :: dlon
real(dp), dimension(size(p1,dim=1)) :: num
real(dp), dimension(size(p1,dim=1)) :: den
real(dp), dimension(size(p1,dim=1)) :: dang

!---

lon0 = d2r * p0(1)
lon1 = d2r * p1(:,1)
lat0 = d2r * p0(2)
lat1 = d2r * p1(:,2)

dlon = lon1 - lon0

num = sqrt((cos(lat1) * sin(dlon))**2 + (cos(lat0) * sin(lat1) - sin(lat0) * cos(lat1) * cos(dlon))**2)
den = sin(lat0) * sin(lat1) + cos(lat0) * cos(lat1) * cos(dlon)

dang = atan2(num,den)

sphdist_vector = rm * dang

end function sphdist_vector

!-------------------------------------------------------------------------------------------------

real(sp) function cardist(p0,p1)

!cartesian distance calculation

use parametersmod

implicit none

real(dp), dimension(2) :: p0  !xycoords of pos 0
real(dp), dimension(2) :: p1  !xycoords of pos 1

!---

cardist = real(sqrt((p0(1) - p1(1))**2 + (p0(2) - p1(2))**2))

end function cardist

!-------------------------------------------------------------------------------------------------

real(sp) function area(lat,gridres)

!this function returns the size of a regular grid cell in square meters.

use parametersmod

implicit none

!arguments

real(dp), intent(in)               :: lat      !(decimal degrees)
real(dp), intent(in), dimension(2) :: gridres  !(decimal degrees)

!local variables

real(dp) :: lonres
real(dp) :: latres
real(dp) :: elevation
real(dp) :: deltalat
real(dp) :: deltalon
real(dp) :: cellarea

!-----------------

lonres = gridres(1)
latres = gridres(2)

elevation = d2r * (lat + 0.5d0 * latres)

deltalon = d2r * lonres
deltalat = d2r * latres

cellarea = 2.d0 * earthrad**2 * deltalon * cos(elevation) * sin(0.5d0 * deltalat)

area = cellarea * 1.e6

end function area

!-------------------------------------------------------------------------------------------------

end module distmod
