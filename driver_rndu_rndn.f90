!Driver rndn rndu


	use matforf90
	use linear_operator
	implicit none

	integer, parameter :: n=4

							!  DP objects
	real (dp) , dimension(:,:),   allocatable ::   x

	integer :: nlines, ncols, iseed

	nlines = 10
	ncols  = 4
	iseed  = 123457

	allocate(x(nlines,ncols))

		!rndns
	print *,"rndns"
	x = rndns(nlines,ncols,iseed)
	call prnt(x)

    print *,"rndnls"
	x = rndnls(nlines,ncols,iseed)
	call prnt(x)

		!rndns
	print *,"rndus"
	x = rndus(nlines,ncols,iseed)
	call prnt(x)

    print *,"rnduls"
	x = rnduls(nlines,ncols,iseed)
	call prnt(x)


end

