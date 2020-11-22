!Example sum prod
 
 
        use matforf90
        use linear_operator
 
        implicit none
 
        integer, parameter :: n=4
        integer :: i
                                     ! DP rank-1 vector
        real (dp) , dimension(n)     ::   x
 
                                     ! DP rank-2 matrix (n,n)
        real (dp) , dimension(n,1)   ::   y
        real (dp) , dimension(n,n)   ::   z
 
 
         z = reshape( (/ 11,12,13,14,                   &
                         21,22,23,24,                   &
                         31,32,33,34,                   &
                         41,42,43,44 /)                 &
                        ,(/4,4/) )
 
 
!SUM            PRODUCT works the same way as SUM
!           returns in a rank-1 vector or scalar
 
        print *, "matrix z"; call prnt(z)
 
        x = sumc(z)   ! sum along rows of each column
        print *,"rank-1 vector of sum on columns "
        call prnt(x)
        print *,"total", sumc(x)
        print *,"shape(x)", shape(x)
 
        x = sumc(.t. z)  ! sum along columns of each row
        print *,"rank-1 vector sum on rows "
        call prnt(x)
        print *,"total", sumc(x)
        print *,"shape(x)", shape(x)
 
                                  ! How to get the sums into a rank-2 vector.
        y(:,1) = sumc(z)   ! sum along dimension 1
        print *,"rank-2 vector sum over columns"
        call prnt(y)
        print *,"total", sumc(y)
        print *,"shape(y)", shape(y)
 
 
!product                product works the same way as sum
!               returns in a rank-1 vector or scalar
 
        print *, "matrix z"; call prnt(z)
 
        x = prodc(z)   ! sum along rows of each column
        print *,"rank-1 vector of sum on columns "
        call prnt(x)
        print *,"total", prodc(x)
        print *,"shape(x)", shape(x)
 
        x = prodc(.t. z)  ! sum along columns of each row
        print *,"rank-1 vector sum on rows "
        call prnt(x)
        print *,"total", prodc(x)
        print *,"shape(x)", shape(x)
 
                                  ! How to get the products into a rank-2 vector.
        y(:,1) = prodc(z)   ! sum along dimension 1
        print *,"rank-2 vector sum over columns"
        call prnt(y)
        print *,"total", prodc(y)
        print *,"shape(y)", shape(y)
 
 
end
 