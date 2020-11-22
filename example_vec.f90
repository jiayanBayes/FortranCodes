!Example vec


	use matforf90
	use linear_operator
	implicit none

	integer, parameter :: n=4

							!  DP objects
	real (dp) , dimension(:),     allocatable ::   x1
	real (dp) , dimension(n,n)                ::   z

	integer :: nlines, ncols, iseed, iseedkeep
	integer :: df 
	character(2) :: achar


       z = reshape( (/ 11,12,13,14,                    &
                       21,22,23,24,                    &
                       31,32,33,34,                    &
                       41,42,43,44 /)                  &
                      ,(/4,4/) )

! Start of program



	allocate( x1(n*n))
	x1 = vec(z)
	print *,"vec form of z matrix";call prnt(x1)


end

