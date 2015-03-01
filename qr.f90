program NL_Routines

  implicit none

  real(8), Dimension(3,3) :: A,Q,R


  A(1,:) = (/2,4,2/)
  A(2,:) = (/-1,0,-4/)
  A(3,:) = (/2,2,-1/)
  Q = 0
  R = 0

  call qr(A,Q,R,3,3)

end program 


  subroutine qr(A,Q,R,m,n)!{{{
  !---------------------------------------------------
  !subroutine-Computes qr of a matrix                !
  ! -See Numerical Linear Algebra for algoritm       !
  !                                                  !
  !Input:  A - Matrix to be operated on              !
  !        m - Number of rows                        !
  !        n - Number of Columns                     ! 
  !                                                  !
  !Output: Q - Q matrix of orthgonal Vectors         !
  !        R- Upper triangular matrix of coeffients  !
  !---------------------------------------------------
      !Declaration of Variables
      integer :: n,m
      real(8),Dimension(m,n),intent(IN) :: A
      real(8),Dimension(m,n),intent(OUT):: Q,R
      real(8),Dimension(m) :: v
    
      !initalize Variables
      v = 0
      
      !QR Algorithm
      do ii = 1,n
          v = A(:,ii)
          do jj = 1,ii-1
              R(jj,ii) = Dot_product(Q(:,jj),A(:,ii))
              v = v - R(jj,ii)*Q(:,jj)
          end do
          R(ii,ii) = norm2(v)
          Q(:,ii) = v/R(ii,ii)
      end do
  
  end subroutine !}}}
!----------------------------------------------------


