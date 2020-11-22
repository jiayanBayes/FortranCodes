!Driver for cdfninv

	use matforf90
	use linear_operator
	implicit none

	integer, parameter ::  n=4
	integer :: iseed

							!  DP objects
	real (dp) , dimension(n)  ::   x
    real (dp) , dimension(:), allocatable  ::   x1

	real (dp) , dimension(n,1)  ::   y

	real (dp) , dimension(n,n)::   z
	real (dp) , dimension(:,:), allocatable  ::   z1


    x = reshape( (/ 1,2,3,4 /) , (/4 /) )
    y = reshape( (/ 10,20,30,40 /) , (/4,1 /) )

    z = reshape( (/ 11,12,13,14,                   &
                    21,22,23,24,                   &
                    31,32,33,34,                   &
                    41,42,43,44 /)                 &
                  ,(/4,4/) )

!Start of program

	iseed = 123457
	y = rndus(n,1,iseed)
	x = y(:,1)
	z = rndus(n,n,iseed)

	allocate(x1(n), z1(n,n))

				! cdfninv
	x1 = cdfninv(x)
	print *,"x vector of unif numbers"; call prnt(x)
	print *,"cdfninv of rank-1 vector"; call prnt(x1)

	z1 = cdfninv(z)
	print *,"z matrix of unif numbers"; call prnt(z)
	print *,"cdfninv of rank-2 matrix"; call prnt(z1)
	

end

