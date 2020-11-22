!Driver random numbers

	use matforf90
	use linear_operator
	implicit none

	integer, parameter :: n = 4

							!  DP objects
	real (dp) , dimension(:,:),   allocatable ::   x
	real (dp) , dimension(:),   allocatable ::   x1

	integer :: nlines, ncols, iseed
	integer :: df 
	character(2) :: achar

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

	deallocate(x)
	allocate(x(n,n))

		! wishart
	df = 7
	x = wishart(df, n, iseed)
	print *, "  n x n cov matrix from Inverted Wishart"
	call prnt( invpd ( x .xt. x) )

		! chi2
	achar = cvtis(n)
	x = chi2dev(n, n, df, iseed)
	print *, "  ",achar , " x ", achar, " matrix of chi2 random deviates "; call prnt(x)

	allocate( x1(n*n))
	x1 = vec(x)
	print *,"vec form";call prnt(x1)

	print *,"mean of chi2 random deviates drawn", meanc(x1)

end

