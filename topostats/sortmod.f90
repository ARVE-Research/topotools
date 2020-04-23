module sortmod

!Modified from: Brainerd, W.F., (year) Guide to Fortran 2003 Programming, 

use parametersmod, only : i4,sp

implicit none

public  :: sortstruct
public  :: sort
private :: quicksort
private :: insertsort

type sortvalues
  real(sp)    :: vals
  integer(i4) :: indx
end type sortvalues

contains

!--------------------------------------------

subroutine sortstruct(struct,key,dir)

!sort a data structure of any type using a real vector

implicit none

class(*), dimension(:), intent(inout) :: struct   !data structure to be sorted
real(sp), dimension(:), intent(inout) :: key      !vector of values that are the sort key
character(10),          intent(in)    :: dir      !direction of the sort

integer(i4), allocatable, dimension(:) :: index

integer(i4) :: i

integer(i4) :: n

class(*), allocatable, dimension(:) :: tmp

!---

n = size(key)

allocate(index(n))

allocate(tmp(n),source=struct)

call sort(key,index)

do i = 1,n
!  struct(i) = tmp(index(i))    !this doesn't quite work
end do

end subroutine sortstruct

!--------------------------------------------

subroutine sort(values,index)

!if the optional argument index is provided, leave the values vector in place,
!and return a vector of indices in the sort order
!this subroutine sorts ascending by default

implicit none

real(sp),    dimension(:), intent(inout)           :: values
integer(i4), dimension(:), intent(inout), optional :: index

type(sortvalues), allocatable, dimension(:) :: vs

!-----

allocate(vs(size(values)))

vs%vals = values

if (present(index)) vs%indx = index

call quicksort(vs)

if (present(index)) then
  index  = vs%indx
else
  values = vs%vals
end if

end subroutine sort

!--------------------------------------------

recursive subroutine quicksort(vs)

implicit none

type(sortvalues), dimension(:), intent(inout) :: vs

integer(i4), parameter :: insert = 6

integer(i4) :: i
integer(i4) :: j
integer(i4) :: n

real(sp) :: chosen

type(sortvalues) :: temp

!-----

n = size(vs)

if (n <= insert) then

  call insertsort(vs)

else

  chosen = vs(n/2)%vals

  i = 0
  j = n + 1

  do

    do
      i = i + 1
      !write(0,*)i,vs(i)%vals,chosen
      if (vs(i)%vals >= chosen) exit
    end do

    do
      j = j - 1
      if (vs(j)%vals <= chosen) exit
    end do

    if (i < j) then

      temp  = vs(i)
      vs(i) = vs(j)
      vs(j) = temp

    else if (i == j) then

      i = i + 1
      exit

    else

      exit

    end if
  end do
  
  if (1 < j) call quicksort(vs(:j))
  if (i < n) call quicksort(vs(i:))

end if

end subroutine quicksort

!--------------------------------------------

subroutine insertsort(vs)

implicit none

type(sortvalues), dimension(:), intent(inout) :: vs
  
integer(i4) :: i
integer(i4) :: j
integer(i4) :: n

type(sortvalues) :: temp

!-----

n = size(vs)

do i = 1,n-1
  do j = i+1,n
    if (vs(i)%vals > vs(j)%vals) then
      temp  = vs(i)
      vs(i) = vs(j)
      vs(j) = temp
    end if
  end do
end do

end subroutine insertsort

!--------------------------------------------

end module sortmod
