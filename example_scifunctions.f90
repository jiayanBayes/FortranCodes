        use matforf90
        use linear_operator
        implicit none
 
        integer, parameter :: n=4
        integer :: i
 
                                         ! DP matrix (n,n)
        real (dp) , dimension(n,n)   ::   z
        real (dp) , dimension(n,n)   ::   z1
 
 
        z = reshape( (/ 11,12,13,14,                    &
                        12,22,23,24,                    &
                        13,23,33,34,                    &
                        14,24,34,44 /)                  &
                       ,(/4,4/) )
 
 ! Start of program

		z1 = log(z)
		print *, "z"; call prnt(z)
		print *, "log(z)"; call prnt(z1)
		print *, "exp(log(z))"; call prnt(exp(z1))
		print *, "sqrt(z)"; call prnt(sqrt(z))
 
 
end