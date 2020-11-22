!Driver for eye

	use matforf90
	use linear_operator
	implicit none

	integer, parameter ::  n=4


	real (dp) , dimension(n,n)::   z

	
	z = eye(n)
	print *,"Identity matrix"; call prnt(z)

	z = 0
	print *,"Matrix of zeros"; call prnt(z)

	z = 1
	print *,"Matrix of ones"; call prnt(z)


end

