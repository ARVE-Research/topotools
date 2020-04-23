module calcstatsmod

use parametersmod, only : i2,i4,sp,dp

implicit none

integer(i4), parameter :: nclasses = 14

real(sp), dimension(nclasses) :: llim

type statvals
  integer(i2) :: elev_med
  integer(i2) :: elev_std
  real(sp)    :: slope_med
  real(sp)    :: slope_std
  real(sp), dimension(nclasses) :: classfrac
end type statvals

contains

!-----------------

subroutine calcstats(area,elev,slope,stats)

use parametersmod, only : missing_i2,missing_sp
use statsmod,      only : median,stdev

implicit none

real(sp),       dimension(:,:), intent(in)  :: area
integer(i2),    dimension(:,:), intent(in)  :: elev
real(sp),       dimension(:,:), intent(in)  :: slope
type(statvals),                 intent(out) :: stats

!--

logical,  allocatable, dimension(:,:) :: valid
real(sp), allocatable, dimension(:)   :: blockvect
real(sp), allocatable, dimension(:)   :: areavect

integer :: xlen
integer :: ylen
integer :: nv

real(sp) :: validarea

integer :: i

!--

xlen = size(elev,dim=1)
ylen = size(elev,dim=2)

nv = count(elev /= missing_i2)

allocate(blockvect(nv))
allocate(areavect(nv))

allocate(valid(xlen,ylen))

where (elev /= missing_i2)
  valid = .true.
elsewhere
  valid = .false.
end where

!elevation

blockvect = pack(real(elev),mask=valid)

stats%elev_med = nint(median(blockvect))

stats%elev_std = nint(10.*stdev(blockvect))

!slope

blockvect = pack(slope,mask=valid)

stats%slope_med = median(blockvect)

stats%slope_std = stdev(blockvect)

areavect = pack(area,mask=valid)

validarea = sum(areavect)

!calculate slope classes

!zero slope

stats%classfrac(1) = sum(areavect,mask=blockvect == 0.) / validarea

!all other classes

do i = 2,nclasses-1
  stats%classfrac(i) = sum(areavect,mask=blockvect > llim(i) .and. blockvect <= llim(i+1)) / validarea
end do

!top class

i=nclasses
stats%classfrac(i) = sum(areavect,mask=blockvect > llim(i)) / validarea

!correct for rounding error

stats%classfrac = min(1.,stats%classfrac)

end subroutine calcstats

end module calcstatsmod
