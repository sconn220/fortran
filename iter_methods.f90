program NL_Routines

  implicit none

  real(8), Dimension(3,3) :: J
  real(8), Dimension(3) :: F,u,delta
  integer :: iters

  F = 0.0     !Function Result
  J = 0.0     !Jacobian
  u = (/0.0,0.0,0.0/)     !Input Variables
  delta = 0.0 !inital guess for results
  iters = 3

  call nl_j(J,u)
  call nl_f(F,u)
  call newton_gmres(F,J,u,delta,3,iters)

!--------------Test Arnoldi---------------!
!{{{ 
!real(8),dimension(5,5) :: A,H,Q_full
!  real(8),dimension(5) :: x
!  integer :: ii,jj
!  H = 0
!  Q_full = 0
  
!  A(1,:) = (/0.8147,0.0975,0.1576,0.1419,0.6557/)
!  A(2,:) = (/0.9058,0.2785,0.9706,0.4218,0.0357/)
!  A(3,:) = (/0.1270,0.5469,0.9572,0.9157,0.8491/)
!  A(4,:) = (/0.9134,0.9575,0.4854,0.7922,0.9340/)
!  A(5,:) = (/0.6324,0.9649,0.8003,0.9595,0.6787/)

!  x = (/0.7577,0.7431,0.3922,0.6555,0.1712/)
!  call arnoldi(A,x,5,3,H,Q_full)
  
  !write h to file
!  open (unit = 1, file = 'h.out', status = 'unknown')
!  do ii = 1,5
!      do jj = 1,5
!          write (1,*),ii,jj,H(ii,jj)
!      end do 
!  end do}}}
!------------------------------------------!  
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
 subroutine arnoldi(A,x,n,iters,H,Q_full)!{{{

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
    integer,intent(in) :: n      !Size of matrix A
    integer,intent(in) :: iters  !Number of iterations in method
    real(8), Dimension(n,n), intent(in) :: A !Matrix A
    real(8), Dimension(n), intent(in) :: x   !Initial Guess
    
    !Internal Variables
    real(8), Dimension(n) :: r   !
    
    !Output Variables
    real(8), Dimension(n,n) :: H          !Hessenberg Matrix
    real(8), Dimension(n,iters) :: Q_full !Set of all Orthonormal Vecotrs
    print*,'Performing Arnoldi'
    if (iters .gt. n) then 
        print*,'Error: Number of iteration exceeds matrix size'
    else
        Q_full(:,1) = x/norm2(x)
        !Number of iteration to solve
        do jj = 1,iters
            r = matmul(A,Q_full(:,jj))
            !Stablized Gram-Schmidt Orthogonalization
            do ii = 1,jj
                H(ii,jj) = dot_product(Q_full(:,ii),r)
                r = r - H(ii,jj)*Q_Full(:,ii)
            end do
            if (jj .lt. n) then 
                H(jj+1,jj) = norm2(r)
                Q_full(:,jj+1) = r/H(jj+1,jj)
            end if 
        end do
    end if
 end subroutine!}}}
!----------------------------------------------------!
subroutine newton_gmres(F,J,x,delta,n,iters)!{{{
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
 !       n - size of F
 !
 !Ouput: Y - Aproximate solution to current time step.
 !----------------------------------------------------------

  !Declear Variables
  !Input Variables
  real(8),dimension(n,n),intent(in) :: J   !Function and Jacobian
  real(8),dimension(n),intent(in) :: F,x   !Function, Input Variables
  real(8),dimension(n),intent(in) :: delta !Guess to solution
  integer,intent(in) :: iters, n           !Max Iterations

  !Internal Variables
  real(8),dimension(n) :: r                !Residue 
  real(8),dimension(n,iters) :: Q_Full     !Set of Orthonormal Basis Functions
  real(8),dimension(n+1,n) :: H            !Hessenberg Matrix 
  integer :: jj, ii
  !Initiate Variables
  r = 0
  Q_Full = 0
  H = 0 
  
  !Arnoldi Process:
  r = -F-matmul(J,delta)  !compute residue
 
  print*,'Performing Arnoldi'
  if (iters .gt. n) then 
      print*,'Error: Number of iteration exceeds matrix size'
  else
      Q_full(:,1) = x/norm2(x)
      !Number of iteration to solve
      do jj = 1,iters
          r = matmul(A,Q_full(:,jj))
          !Stablized Gram-Schmidt Orthogonalization
          do ii = 1,jj
              H(ii,jj) = dot_product(Q_full(:,ii),r)
              r = r - H(ii,jj)*Q_Full(:,ii)
          end do
          if (jj .lt. n) then 
              H(jj+1,jj) = norm2(r)
              Q_full(:,jj+1) = r/H(jj+1,jj)
          end if 
      end do
  end if

end subroutine!}}}
!----------------------------------------------------!
subroutine nl_j(J,u) !{{{
  real(8),Dimension(3,3) :: J
  real(8),Dimension(3) :: u
  
  J(1,1) = 0.001
  J(1,2) = -0.001
  J(1,3) = 1.0
  J(2,1) = -0.001
  J(2,2) = 0.001+0.001*3.0*u(2)**2
  J(2,3) = 0.0
  J(3,1) = 1.0
  J(3,2) = 0.0
  J(3,3) = 0.0

end subroutine !}}}
!----------------------------------------------------!
subroutine nl_f(F,u) !{{{
  real(8),Dimension(3):: F,u
  F(1) = 0.001*(u(1)-u(2))+u(3)
  F(2) = 0.001*(u(2)-u(1))+0.001*(u(2)**3)
  F(3) = u(3)-1
end subroutine !}}}
