program test_eigb_f
 use tt_lib
 use ttaux_lib
 use ttop_lib
 use ttio_lib
!  use ttm_lib
 use time_lib
 implicit none


 real(8), allocatable :: arr(:)
 type(dtt) :: tt
 integer(8), allocatable :: n(:)
 integer(8) :: d, nel, rmax
 integer(8) :: i
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
 call say(tt)
 deallocate(n)
 deallocate(arr)

end program
