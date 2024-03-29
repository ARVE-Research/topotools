module calcstatsmod

use parametersmod, only : i2,i4,sp,dp

implicit none

integer(i4), parameter :: nclasses_slope  = 14
integer(i4), parameter :: nclasses_aspect = 8
integer(i4), parameter :: nclasses_cti    = 13
integer(i4), parameter :: nclasses_hand   = 11

!------------------------------------------------

type classbin
  real(sp), dimension(nclasses_slope)  :: slope_bin
  real(sp), dimension(nclasses_aspect) :: aspect_bin
  real(sp), dimension(nclasses_cti)    :: cti_bin
  real(sp), dimension(nclasses_hand)   :: hand_bin
end type classbin

!------------------------------------------------

type statvals
  real(sp) :: elev_med
  real(sp) :: elev_std
  real(sp) :: slope_med
  real(sp) :: slope_std
  real(sp) :: aspect_med
  real(sp) :: aspect_std
  real(sp) :: cti_med
  real(sp) :: cti_std
  real(sp) :: hand_med
  real(sp) :: hand_std
  real(sp) :: areafrac
  real(sp), dimension(nclasses_slope)  :: classfrac_slope
  real(sp), dimension(nclasses_aspect) :: classfrac_aspect
  real(sp), dimension(nclasses_cti)    :: classfrac_cti
  real(sp), dimension(nclasses_cti)    :: classfrac_hand
end type statvals

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine calcstats(area,elev,slope,cti,hand,stats,llim)

use parametersmod, only : missing_sp
use statsmod,      only : median,stdev

implicit none

real(sp), dimension(:,:), intent(in)  :: area
real(sp), dimension(:,:), intent(in)  :: elev
real(sp), dimension(:,:), intent(in)  :: slope
real(sp), dimension(:,:), intent(in)  :: cti
real(sp), dimension(:,:), intent(in)  :: hand
type(classbin),           intent(in)  :: llim
type(statvals),           intent(out) :: stats

logical,  allocatable, dimension(:,:) :: valid
real(sp), allocatable, dimension(:)   :: blockvect
real(sp), allocatable, dimension(:)   :: areavect

integer :: xlen
integer :: ylen
integer :: nv

real(sp) :: validarea

integer :: i

!---

xlen = size(elev,dim=1)
ylen = size(elev,dim=2)

nv = count(elev /= missing_sp)

allocate(blockvect(nv))
allocate(areavect(nv))

allocate(valid(xlen,ylen))

where (elev /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where


areavect = pack(area,mask=valid)

validarea = sum(areavect)

!areafrac-----------

stats%areafrac = validarea / sum(area)

if (all(valid)) stats%areafrac = 1. !Adjust for rounding error


!elevation-----------

blockvect = pack(elev,mask=valid)

stats%elev_med = median(blockvect)

stats%elev_std = 10.*stdev(blockvect)

!slope-----------

blockvect = pack(slope,mask=valid)

stats%slope_med = median(blockvect)

stats%slope_std = stdev(blockvect)

  !calculate slope classes

  !zero slope

  stats%classfrac_slope(1) = sum(areavect,mask=blockvect == 0.) / validarea

  !all other classes

  do i = 2,nclasses_slope-1
    stats%classfrac_slope(i) = sum(areavect,mask=blockvect > llim%slope_bin(i) .and. blockvect <= llim%slope_bin(i+1)) / validarea
  end do

  !top class

  i=nclasses_slope
  stats%classfrac_slope(i) = sum(areavect,mask=blockvect > llim%slope_bin(i)) / validarea

  !correct for rounding error

  stats%classfrac_slope = min(1.,stats%classfrac_slope)


deallocate(blockvect)


!cti-----------

nv = count(cti /= missing_sp)

where (cti /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

allocate(blockvect(nv))

blockvect = pack(cti,mask=valid)

stats%cti_med = median(blockvect)

stats%cti_std = stdev(blockvect)

  !first cti class

  stats%classfrac_cti(1) = sum(areavect,mask=blockvect >= 0. .and. blockvect <= llim%hand_bin(1)) / validarea

  !all other classes

  do i = 2,nclasses_cti-1
    stats%classfrac_cti(i) = sum(areavect,mask=blockvect > llim%cti_bin(i-1) .and. blockvect <= llim%cti_bin(i)) / validarea
  end do

  stats%classfrac_cti(13) = sum(areavect,mask=blockvect > llim%cti_bin(12)) / validarea

deallocate(blockvect)

!hand-----------

nv = count(hand /= missing_sp)

where (hand /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

allocate(blockvect(nv))

blockvect = pack(hand,mask=valid)

stats%hand_med = median(blockvect)

stats%hand_std = stdev(blockvect)

  !first hand class

  stats%classfrac_hand(1) = sum(areavect,mask=blockvect >= 0. .and. blockvect <= llim%hand_bin(1)) / validarea

  !all other classes

  do i = 2, nclasses_hand-1
    stats%classfrac_hand(i) = sum(areavect,mask=blockvect > llim%hand_bin(i-1) .and. blockvect <= llim%hand_bin(i)) / validarea
  end do

  stats%classfrac_hand(11) = sum(areavect,mask=blockvect > llim%hand_bin(10)) / validarea


end subroutine calcstats

!-------------------------------------------------------------------------------

subroutine calcstats_aspect(area,elev,slope,aspect,cti,hand,stats,llim)

use parametersmod, only : missing_sp
use statsmod,      only : median,stdev,circle_mean,circle_stdev

implicit none

real(sp), dimension(:,:), intent(in)  :: area
real(sp), dimension(:,:), intent(in)  :: elev
real(sp), dimension(:,:), intent(in)  :: slope
real(sp), dimension(:,:), intent(in)  :: aspect
real(sp), dimension(:,:), intent(in)  :: cti
real(sp), dimension(:,:), intent(in)  :: hand
type(classbin),           intent(in)  :: llim
type(statvals),           intent(out) :: stats

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

nv = count(elev /= missing_sp)

allocate(blockvect(nv))
allocate(areavect(nv))

allocate(valid(xlen,ylen))

where (elev /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

areavect = pack(area,mask=valid)

validarea = sum(areavect)

!areafrac-----------

stats%areafrac = validarea / sum(area)

if (all(valid)) stats%areafrac = 1.


!elevation-----------

blockvect = pack(elev,mask=valid)

stats%elev_med = median(blockvect)

stats%elev_std = 10.*stdev(blockvect)

!slope-----------

blockvect = pack(slope,mask=valid)

stats%slope_med = median(blockvect)

stats%slope_std = stdev(blockvect)

  !calculate slope classes

  !zero slope

  stats%classfrac_slope(1) = sum(areavect,mask=blockvect == 0.) / validarea

  !all other classes

  do i = 2,nclasses_slope-1
    stats%classfrac_slope(i) = sum(areavect,mask=blockvect > llim%slope_bin(i) .and. blockvect <= llim%slope_bin(i+1)) / validarea
  end do

  !top class

  i=nclasses_slope
  stats%classfrac_slope(i) = sum(areavect,mask=blockvect > llim%slope_bin(i)) / validarea

  !correct for rounding error

  stats%classfrac_slope = min(1.,stats%classfrac_slope)


deallocate(blockvect)



!aspect-----------

nv = count(aspect /= missing_sp)

allocate(blockvect(nv))

where (aspect /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

blockvect = pack(aspect,mask=valid)

areavect = pack(area,mask=valid)

validarea = sum(areavect)

!write(0,*) 'valid aspect:',  nv
!write(0,*) 'aspect values:', aspect
!write(0,*) 'valid area: ', validarea

stats%aspect_med = circle_mean(blockvect)

stats%aspect_std = circle_stdev(blockvect)


  !first aspect class

  stats%classfrac_aspect(1) = (sum(areavect,mask=blockvect >= 0. .and. blockvect <= llim%aspect_bin(1)) &
                              + sum(areavect,mask=blockvect > llim%aspect_bin(8) .and. blockvect <= 360.)) / validarea

  !all other classes

  do i = 2,nclasses_aspect
    stats%classfrac_aspect(i) = sum(areavect,mask=blockvect > llim%aspect_bin(i-1) &
                              .and. blockvect <= llim%aspect_bin(i)) / validarea
  end do

  !correct for rounding error

  stats%classfrac_slope = min(360.,stats%classfrac_slope)

deallocate(blockvect)

!cti-----------

nv = count(cti /= missing_sp)

where (cti /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

allocate(blockvect(nv))

blockvect = pack(cti,mask=valid)

stats%cti_med = median(blockvect)

!print*, '                 '
!print*, nv, blockvect, slope

stats%cti_std = stdev(blockvect)

  !first cti class

  stats%classfrac_cti(1) = sum(areavect,mask=blockvect >= 0. .and. blockvect <= llim%cti_bin(1)) / validarea

  !all other classes

  do i = 2,nclasses_cti-1
    stats%classfrac_cti(i) = sum(areavect,mask=blockvect > llim%cti_bin(i-1) .and. blockvect <= llim%cti_bin(i)) / validarea
  end do

  stats%classfrac_cti(13) = sum(areavect,mask=blockvect > llim%cti_bin(12)) / validarea

deallocate(blockvect)

!hand-----------

nv = count(hand /= missing_sp)

where (hand /= missing_sp)
  valid = .true.
elsewhere
  valid = .false.
end where

allocate(blockvect(nv))

blockvect = pack(hand,mask=valid)

stats%hand_med = median(blockvect)

stats%hand_std = stdev(blockvect)

  !first hand class

  stats%classfrac_hand(1) = sum(areavect,mask=blockvect >= 0. .and. blockvect <= llim%hand_bin(1)) / validarea

  !all other classes

  do i = 2, nclasses_hand-1
    stats%classfrac_hand(i) = sum(areavect,mask=blockvect > llim%hand_bin(i-1) .and. blockvect <= llim%hand_bin(i)) / validarea
  end do

  stats%classfrac_hand(11) = sum(areavect,mask=blockvect > llim%hand_bin(10)) / validarea


end subroutine calcstats_aspect

!-------------------------------------------------------------------------------

end module calcstatsmod
