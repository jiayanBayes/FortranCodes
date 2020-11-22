!Driver matproducts

        use matforf90
        use linear_operator

        implicit none

        integer, parameter :: n=4
        integer :: i
 
        real (dp) , dimension(n)     ::   x
                                             ! DP rank-2 matrix (n,n)
        real (dp) , dimension(n,1)   ::   y
        real (dp) , dimension(:,:), allocatable   ::   y1
        real (dp) , dimension(n,n)   ::   z
        real (dp) , dimension(:,:), allocatable   ::   z1

        x = reshape( (/ 1,2,3,4 /) , (/4 /) )
        y = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

        z = reshape( (/ 11,12,13,14,                   &
                        21,22,23,24,                   &
                        31,32,33,34,                   &
                        41,42,43,44 /)                 &
                      ,(/4,4/) )

!MULTIPLICATION

        allocate(y1(n,1))

                                    ! mprd
        y1 = mprd(z,y)
        print *,"z*y";  call prnt(y1)

                                    ! tprd
        print*, "Matmul"
        print*, matmul(transpose(z),y)
        y1 = tprd(z,y)
        print *,"z'y";  call prnt(y1)

                                    ! prdt
        print*, "Matmul"
        print*, matmul(z,transpose(z))

        allocate(z1(4,4))
        z1 = prdt(z,z)
        print *,"z * z'"; call prnt(z1)

                        ! linear operator

        y1 = z .x. y
        print *,"operator z .x. y";  call prnt(y1)

        y1 = z .tx. y
        print *,"operator z .tx. y ";  call prnt(y1)

        z1 = z .xt. z
        print *,"operator z .xt. z";  call prnt(z1)

		print *, "dot product", dot_product(y(:,1),x)
end