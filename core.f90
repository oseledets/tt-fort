!This file is a very simple attempt to write down 
!simple TT stuff (i.e., matrix-by-matrix or matrix-by-vector product) in the TT-format
!Many of the procedures are already in the SDV stuff (adding, dot product, etc, but
!we do not yet have the tt-matrix structure; we can (in principle) use it like
!mat_mat(n,m,k,tt1,tt2) (and put it inside tt.f90); then the call would be quite simple;  
module core
  real(8),  allocatable :: result_core(:)
  complex(8), allocatable :: zresult_core(:)
 contains
 
 subroutine dealloc()
   if ( allocated(result_core) ) then
      deallocate(result_core)
   end if
   if ( allocated(zresult_core) ) then
      deallocate(zresult_core)
   end if
 end subroutine dealloc
 

 subroutine dmat_mat_loc(n,m,k,r1,r2,p1,p2,core1,core2,res_core) 
     use matrix_util, only : dtransp, dperm321
     integer, intent(in) :: n,m,k,r1,r2,p1,p2
     real(8), intent(in) :: core1(*), core2(*)
     real(8), intent(inout) :: res_core(*)
     real(8) :: tmp_core1(r1*n*m*r2), tmp_core2(p1*m*k*p2)
     
     !core1(a1,i1,j1,a2) -> core1(a2,a1,i1,j1) * core2(j1,k1,b2,b1) -> res_core(a2,a1,i1,k1,b2,b1) -> 
     call dtransp(r1*n*m,r2,core1,tmp_core1)
     call dtransp(p1,m*k*p2,core2,tmp_core2)
     call dgemm('n','n',r2*r1*n,k*p2*p1,m,1d0,tmp_core1,r2*r1*n,tmp_core2,m,0d0,res_core,r2*r1*n)
     !We have a2,a1,i1,k1,b2,b1 -> but! we need (?) a1, b1, 
     call dperm321(r2,r1*n*k*p2,p1,res_core)
 end subroutine dmat_mat_loc
 
 subroutine zmat_mat_loc(n,m,k,r1,r2,p1,p2,core1,core2,res_core) 
     use matrix_util, only : ztransp, zperm321
     integer, intent(in) :: n,m,k,r1,r2,p1,p2
     complex(8), intent(in) :: core1(*), core2(*)
     complex(8), intent(inout) :: res_core(*)
     complex(8) :: tmp_core1(r1*n*m*r2), tmp_core2(p1*m*k*p2)
     
     !core1(a1,i1,j1,a2) -> core1(a2,a1,i1,j1) * core2(j1,k1,b2,b1) -> res_core(a2,a1,i1,k1,b2,b1) -> 
     call ztransp(r1*n*m,r2,core1,tmp_core1)
     call ztransp(p1,m*k*p2,core2,tmp_core2)
     call zgemm('n','n',r2*r1*n,k*p2*p1,m,(1d0, 0d0),tmp_core1,r2*r1*n,tmp_core2,m,(0d0, 0d0),res_core,r2*r1*n)
     call zperm321(r2,r1*n*k*p2,p1,res_core)
 end subroutine zmat_mat_loc


 subroutine dmat_mat(d,n,m,k,cr1,cr1size,cr2,cr2size,r1,r2,rres)
    integer, intent(in) :: d
    integer, intent(in) :: cr1size, cr2size
    integer, intent(in) :: n(d)
    integer, intent(in) :: m(d)
    integer, intent(in) :: k(d)
    integer, intent(in) :: r1(d+1)
    integer, intent(in) :: r2(d+1)
    integer, intent(out) :: rres(d+1)
    real(8), intent(in) :: cr1(cr1size)
    real(8), intent(in) :: cr2(cr2size)
    integer :: i,pos,pos1,pos2,mem
    
    rres(1:d+1) = r1(1:d+1) * r2(1:d+1)
    mem = 0
    do i = 1,d
       mem = mem + n(i) * k(i) * rres(i) * rres(i+1) 
    end do 
    if ( allocated(result_core) ) then
       if ( size(result_core) < mem ) then
          deallocate(result_core)
          allocate(result_core(mem))
       end if
    else
        allocate(result_core(mem))
    end if
    pos1 = 1
    pos2 = 1
    pos = 1
    do i = 1,d
       call dmat_mat_loc(n(i),m(i),k(i),r1(i),r1(i+1),r2(i),r2(i+1),cr1(pos1),cr2(pos2),result_core(pos))
       pos1 = pos1 + r1(i) * n(i) * m(i) * r1(i+1)
       pos2 = pos2 + r2(i) * m(i) * k(i) * r2(i+1)
       pos = pos + rres(i) * n(i) * k(i) * rres(i+1)
    end do 
 end subroutine dmat_mat

 subroutine zmat_mat(d,n,m,k,cr1,cr1size,cr2,cr2size,r1,r2,rres)
    integer, intent(in) :: d
    integer, intent(in) :: cr1size, cr2size
    integer, intent(in) :: n(d)
    integer, intent(in) :: m(d)
    integer, intent(in) :: k(d)
    integer, intent(in) :: r1(d+1)
    integer, intent(in) :: r2(d+1)
    integer, intent(out) :: rres(d+1)
    complex(8), intent(in) :: cr1(cr1size)
    complex(8), intent(in) :: cr2(cr2size)
    integer :: i,pos,pos1,pos2,mem
    
    rres(1:d+1) = r1(1:d+1) * r2(1:d+1)
    mem = 0
    do i = 1,d
       mem = mem + n(i) * k(i) * rres(i) * rres(i+1) 
    end do 
    if ( allocated(zresult_core) ) then
       if ( size(zresult_core) < mem ) then
          deallocate(zresult_core)
          allocate(zresult_core(mem))
       end if
    else
        allocate(zresult_core(mem))
    end if
    pos1 = 1
    pos2 = 1
    pos = 1
    do i = 1,d
       call zmat_mat_loc(n(i),m(i),k(i),r1(i),r1(i+1),r2(i),r2(i+1),cr1(pos1),cr2(pos2),zresult_core(pos))
       pos1 = pos1 + r1(i) * n(i) * m(i) * r1(i+1)
       pos2 = pos2 + r2(i) * m(i) * k(i) * r2(i+1)
       pos = pos + rres(i) * n(i) * k(i) * rres(i+1)
    end do 
 end subroutine zmat_mat
end module core
