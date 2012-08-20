!The Fortran module (to be wrapped with f2py)
!for doing certain TT-matrix stuff 
module tt_matrix
  use tensor_util_lib
contains
  subroutine tt_mv_full(n,d,m,r,ps,crm,crm_size,x,x_size,rb,y,y_size)
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: ps(d+1)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: m(d)
    integer, intent(in) :: crm_size
    integer, intent(in) :: x_size,y_size
    integer, intent(in) :: rb !In the case we multiply by a matrix
    double precision, intent(in) :: crm(crm_size)
    double precision, intent(in) :: x(x_size,rb)
    double precision, intent(out) :: y(y_size)  !We have to provide the size of the output array here
    !Local variables
    integer i,max_mem,c_size,c0
    double precision, allocatable :: b(:),c(:)
    double precision, allocatable :: cr(:),cr1(:)
    integer :: sz1(4),sz2(2)
    integer :: prm1(4),prm2(2)
    max_mem=0
    !print *,'d=',d,'n=',n,'m=',m
    do i=1,d
       max_mem = max(max_mem,r(i)*n(i)*m(i)*r(i+1))
    end do
    !print *,'max_mem1=',max_mem
    allocate(cr(max_mem))
    allocate(cr1(max_mem))
    max_mem = 0
    c_size=product(m)*rb
    do i=1,d
       max_mem = max(max_mem,c_size)
       c_size = n(i)*r(i+1)*c_size/(r(i)*m(i)) 
    end do 
    !print *,'max_mem2=',max_mem

    allocate(c(max_mem))
    allocate(b(max_mem))
    prm1(1) = 2
    prm1(2) = 4
    prm1(3) = 1
    prm1(4) = 3
    prm2(1) = 2
    prm2(2) = 1
    c_size=product(m)*rb

    call dcopy(c_size,x,1,c,1)
    do i = 1, d  
       cr(1:r(i)*n(i)*m(i)*r(i+1))=crm(ps(i):ps(i+1)-1)
       sz1(1) = r(i)
       sz1(2) = n(i)
       sz1(3) = m(i)
       sz1(4) = r(i+1)
       call permute(cr,product(sz1),sz1,4,cr1,prm1)
       c0 = c_size/(r(i)*m(i))
       !print *,'c before:',c(1:8)
       !print *,'cr:',cr1(:)
       call dgemm('n','n',n(i)*r(i+1),c0,r(i)*m(i),1.d0,cr1,n(i)*r(i+1),c,r(i)*m(i),0.d0,b,n(i)*r(i+1))
       !print *,'b after:',b(1:8)
       !The size of b is 
       !n(i)*r(i+1)*c_size/(r(i)*m(i))
       c_size = n(i)*r(i+1)*c0
       sz2(1) = n(i)
       sz2(2) =  r(i+1)*c0
       call permute(b,product(sz2),sz2,2,c,prm2)
       !print *,'c after:',c(1:8)
    end do
    call dcopy(c_size,c,1,y,1)
    !cr=cra(ps(k):ps(k+1)-1);
    !  cr=reshape(cr,[r(k),n(k),m(k),r(k+1)]);
    !  cr=permute(cr,[2,4,1,3]); cr=reshape(cr,[n(k)*r(k+1),r(k)*m(k)]);
    !  M=numel(c);
    !  c=reshape(c,[r(k)*m(k),M/(r(k)*m(k))]);
    !  c=cr*c; c=reshape(c,[n(k),numel(c)/n(k)]);
    !  c=permute(c,[2,1]);
    !end
    !c=c(:); c=reshape(c,[rb,numel(c)/rb]);
    !c=c.';
    deallocate(c)
    deallocate(b)
    deallocate(cr)
    deallocate(cr1)
    
  end subroutine tt_mv_full
!    subroutine tt_mv(n,d,m,r,ps,crm,crm_size,rx,psx,crx,x_size,ry,psy,cry,y_size)
!    integer, intent(in) :: d
!    integer, intent(in) :: n(d)
!    integer, intent(in) :: ps(d+1)
!    integer, intent(in) :: r(d+1)
!    integer, intent(in) :: m(d)
!    integer, intent(in) :: crm_size
!    integer, intent(in) :: x_size,y_size
!    integer, intent(in) :: rb !In the case we multiply by a matrix
!    double precision, intent(in) :: crm(crm_size)
!    double precision, intent(in) :: crx(x_size)
!    double precision, intent(out) :: cry(y_size)  !We have to provide the size of the output array here
    !Local variables
!    integer i,max_mem,c_size,c0
!    double precision, allocatable :: b(:),c(:)
!    double precision, allocatable :: cr(:),cr1(:)
!    integer :: sz1(4),sz2(2)
!    integer :: prm1(4),prm2(2)
!    end subroutine tt_mv
end module tt_matrix
