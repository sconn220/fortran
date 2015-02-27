program NL_Routines

  implicit none

  real(8), Dimension(2,2) :: A
  real(8), Dimension(2,2) :: Q,R

  A = 0
  Q = 0
  R = 0

  call qr(A,Q,R)

end program 


  subroutine qr(A,Q,R)
         
      integer :: n,m
      integer :: arrsize
      real(8),Dimension(2,2) :: A,Q,R
      real(8),Dimension(2,1) :: v
      m = size(A(:,1))
      n = size(A(1,:))
     
      do ii = 1,n
          v = A(:,ii)
 !         do jj = 1,ii-1
 !             R(jj,ii) = Q(:,jj)*A(:,ii)
 !             v = v - R(jj,ii)*Q(:,jj)
 !         end do
 !         R(ii,ii) = norm2(v)
 !         Q(:,ii) = v/R(ii,ii)
      end do
  
  !print*, R
  !:print*, Q

  end subroutine
