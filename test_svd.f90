program test_eigb_f
 use tt_lib
 use ttaux_lib
 use ttop_lib
 use ttio_lib
!  use ttm_lib
 use time_lib
 implicit none


 real(8), allocatable :: arr(:)
 type(dtt), target :: tt
 integer(8), allocatable :: n(:)
 integer(8) :: d, nel, rmax
 integer(8) :: i
 integer(8), pointer :: r(:), nn(:)
 real(8), pointer :: current_core(:, :, :) !They are all 3d-dimensional, including the border ones
 rmax = 20

 
 !Create an array
 d = 5 !5-dim array
 allocate(n(d))
 n(:) = 2
 allocate(arr(product(n)))
 nel = product(n)
 do i = 1, nel
     arr(i) = i
 end do
 !Call the function
 call svd(n, arr, tt, 1d-8)
 !call say(tt)

 !The TT-structure is:
 ! tt % r(0:d+1) -> the ranks
 ! tt % n(1:d) -> mode sizes
 ! tt % u(1:d) -> 3d arrays o
 r => tt % r(0:d)
 nn => tt % n(1:d)
 print *, 'ranks:', r(1:d+1)
 print *, 'mode sizes:', nn
 print *, 'printing cores'
 do i = 1, d
     current_core => tt % u(i) % p(1:r(i), 1:nn(i), 1:r(i+1))
     print *, 'The core', i, ' is', current_core
 end do
 deallocate(n)
 deallocate(arr)
 !call dealloc(tt)
end program
