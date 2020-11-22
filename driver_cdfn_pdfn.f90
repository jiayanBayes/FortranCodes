!Driver for cdfn pdfn

	use matforf90
	use linear_operator
	implicit none

	integer, parameter ::  n=4

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


	z = z .div. y

	print *, "z: ";call prnt(z)

	allocate(x1(n), z1(n,n))

				! cdfn
	x1 = cdfn(x)
	print *,"cdfn of rank-1 vector"; call prnt(x1)

	z1 = cdfn(z)
	print *,"cdfn of rank-2 matrix"; call prnt(z1)
	
				! pdfn
	x1 = pdfn(x)
	print *,"pdfn of rank-1 vector"; call prnt(x1)

	z1 = pdfn(z)
	print *,"pdfn of rank-2 matrix"; call prnt(z1)
	

end

