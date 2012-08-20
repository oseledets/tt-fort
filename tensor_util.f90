module tensor_util_lib
use time_lib
contains
  subroutine next_ind(ind,n,d)
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(inout) :: ind(d)
    integer k
    k=1
    ind(1) = ind(1) + 1
    do while (ind(k) > n(k) .and.  k < d)
       ind(k)=ind(k)-n(k)
       ind(k+1) = ind(k+1) + 1
       k = k + 1
    end do
   !print *,'get_ind called, i=',i
   !print *,'get_ind called, d=',d
   !print *,'get_ind called, ind=',ind
  end subroutine next_ind

  subroutine permute(a,M,n,d,b,prm)
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: M
    integer, dimension(d), intent(in) :: prm
    double precision, intent(in) :: a(M)
    double precision, intent(out) :: b(M)
    double precision t1,t2
    integer i,j
    integer :: ind(d) 
    integer :: tmp(d)
    !A really stupid code:  array (b) is just a permutation of the array a.
    !A one-to-one index correspondence is given as
    !a(i(1) + (i2-1)*n1 + (i3-1)*n1*n2 + ... + (id-1)*n1*...*n(d-1)) = 
    !b(i(prm(1))+ (i(prm(2)))*n(prm(1))....
    !for all multidimensional indices, do:
    !b(dot(ind(prm), sz1(prm)))=a(dot(ind,sz1))
    tmp(2:d)=n(prm(1:d-1))
    tmp(1)=1
    do i=2,d
       tmp(i)=tmp(i-1)*tmp(i)
    end do
    ind(1:d)=1
    do i=1,M
       !call get_ind(i,ind,n,d)
       j=dot_product(tmp,ind(prm)-1)+1 !tmp = 1,2 (4x2) ind=(2,2) ind-1=(1,1)
       !1+2+1 = 4
       !ind = 1,4 (?) ind(prm) = 4,1 ind(prm-1)=3,0 3+2*0 + 1 =4 (?????)
       !2,2 -> 2 + 1*2 = 4
       !print *,'ind:',ind(1:d),'tmp:',tmp
       call next_ind(ind,n,d)
       !print *,'tmp',tmp
       !print *,'ind:',ind
       !print *,'prm:',prm
       !print *,'ind(prm):',ind(prm)
       !print *,'i=',i,'j=',j,'prm:',prm,'sz:',n
       b(j)=a(i)
    end do
  end subroutine permute
  subroutine get_ind(i,ind,n,d)
    integer, intent(in) :: d

    integer, intent(in) :: n(d)
    integer, intent(out) :: ind(d)
    integer, intent(in) :: i
    integer j
    j = i
    do k=1,d
       !print *,'inside get_ind loop,j=',j
       !print *,'inside get_ind loop,n=',n
       ind(k) = mod(j-1,n(k))+1
       j=(j-ind(k))/n(k)+1
    end do
   !print *,'get_ind called, i=',i
   !print *,'get_ind called, d=',d
   !print *,'get_ind called, ind=',ind

  end subroutine get_ind

end module tensor_util_lib
