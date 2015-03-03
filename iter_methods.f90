program NL_Routines

  implicit none

  real(8), Dimension(3,3) :: A,Q,R
  real(8), Dimension(3) :: x
  integer :: iters
  !------QR testing--------!
  A(1,:) = (/2,4,2/)
  A(2,:) = (/-1,0,-4/)
  A(3,:) = (/2,2,-1/)
  !Q = 0
  !R = 0
  !call qr(A,Q,R,3,3)

  !------Arnoldi testing----!
  x = (/1,2,3/)
  iters = 2
  call arnoldi(A,x,3,iters)
end program 

!----------------------------------------------------!
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
              R(jj,ii) = Dot_Product(Q(:,jj),A(:,ii))
              v = v - R(jj,ii)*Q(:,jj)
          end do
          R(ii,ii) = norm2(v)
          Q(:,ii) = v/R(ii,ii)
      end do
  
  end subroutine !}}}
!----------------------------------------------------!
 subroutine arnoldi(A,x,n,iters)!{{{

 !--------------------------------------------------!
 !Subroutine: arnoldi - Computes the orthonormal    !
 !            basis for the Krylov subspace of      !
 !            matrix A                              !
 !Input: A- an nxn matrix                           !
 !                                                  !
 !Output: H-Hessenberg Matrix resulting from method !
 !                                                  !
 !--------------------------------------------------!
    
    !Declaration of Variables
    !Input Variables
    integer,intent(IN) :: n      !Size of matrix A
    integer,intent(IN) :: iters  !Number of iterations in method
    real(8), Dimension(n,n), intent(IN) :: A !Matrix A
    real(8), Dimension(n), intent(IN) :: x   !Initial Guess
    
    !Internal Variables
    real(8), Dimension(n) :: q   !Current iteration's orthonormal vectors
    real(8), Dimension(n) :: r   !
    
    !Output Variables
    real(8), Dimension(n+1,n) :: H          !Hessenberg Matrix
    real(8), Dimension(n,iters) :: Q_full !Set of all Orthonormal Vecotrs
    
    if (iters >= n) then 
        print*,'Error: Number of iteration exceeds matrix size'
    else
        q = x/norm2(x)
        Q_full(:,1) = q
        !Number of iteration to solve
        do jj = 1,iters
            r = matmul(A,q)
            !Stablized Gram-Schmidt Orthogonalization
            do ii = 1,jj
                H(ii,jj) = dot_product(q,r)
                r = r - q*H(ii,jj)
            end do
            H(jj+1,jj) = norm2(r)
            Q_full(:,jj) = q
            q = r/H(jj+1,jj)
        end do
    end if
 end subroutine!}}}
!----------------------------------------------------!
subroutine newton_gmres(F,X)!{{{
 !----------------------------------------------------------
 !Subroutine: Generalized Minimal Residual Method (GMRES)
 !            Approximates solution to a nonsymmetrix
 !            system of equations system of linear equations
 !            by looking for the minimium residue in a krylov
 !            subspace. 
 !            This method is combined with newton's method to 
 !            Solve non-linear systems
 !
 !Input: F - The Non-linear Function from R^n to R^n
 !       J - Jacobian
 !       X - inital condition 
 !       delta - Inital Guess
 !       R - Residue
 !
 !Ouput:
 !----------------------------------------------------------
  
  !Declear Variables
  
  !Initiate Variables

  !Arnoldi Process

  !Least Squares Using QR Process


end subroutine!}}}
!----------------------------------------------------!
