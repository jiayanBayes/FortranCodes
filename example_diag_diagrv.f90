!Driver diag diagrv

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


			!diag diagrv

	allocate(z1(n,n), y1(n,1))

	x = diag(z)
	print *, "x = diag(z)";	call prnt(x)

	z1 = diagrv(z,y)
	print *, "diagrv(z,y)";	call prnt(z1)

	z1 = diagrv(z,x)
	print *, "diagrv(z,x)";	call prnt(z1)


end