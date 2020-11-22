!Driver invpd chol
 
 
	use matforf90
	use linear_operator
	implicit none
 
	integer, parameter :: n=4
	integer :: i
 
					! DP matrix (n,n)
	real (dp) , dimension(n,n)   ::   z
	real (dp) , dimension(:,:), allocatable ::   z1
 
 
	z = reshape( (/ 11,12,13,14,        &				
                        12,22,23,24,        &
                        13,23,33,34,        &
                        14,24,34,44 /)      &
                        ,(/4,4/) )
 
!invpd and chol
 
	allocate(z1(n,n))
	z1 = invpd(z)
	print *,"z inverse"; 	call prnt(z1)
	print *,"check: z*zinv = I "; call prnt( z .x. z1)
 
 
	z1 = chol(z)
	print *,"chol(z)";	call prnt(z1)
 
	print *,"check P'P = Z"; call prnt( chol(z) .tx. chol(z))
 
end
 
 