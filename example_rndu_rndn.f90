!Example rndn rndu


	use matforf90
	use linear_operator
	implicit none

	integer, parameter :: n=4

							!  DP objects
	real (dp) , dimension(:,:),   allocatable ::   x
	real (dp) , dimension(n,n)                ::   z

	integer :: nlines, ncols, iseed, iseedkeep

       z = reshape( (/ 11,12,13,14,                    &
                       21,22,23,24,                    &
                       31,32,33,34,                    &
                       41,42,43,44 /)                  &
                      ,(/4,4/) )

	nlines = 10
	ncols  = 4
	iseed  = 123457

	allocate(x(nlines,ncols))

		!rndns
	print *,"rndns"
	x = rndns(nlines,ncols,iseed);  call prnt(x)

    print *,"rndnls"
	x = rndnls(nlines,ncols,iseed)
	call prnt(x)

		!rndns
	print *,"rndus"
	x = rndus(nlines,ncols,iseed);	call prnt(x)

	iseedkeep = iseed
	x = z .xt. rndus(nlines,ncols,iseed)
	print *, " z * rndus' "; call prnt(x)

	iseed = iseedkeep
	x = mprd(z , .t. rndus(nlines,ncols,iseed) )
	print *, " z * rndus' "; call prnt(x)



end

