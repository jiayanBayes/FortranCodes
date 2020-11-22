        use matforf90
        use linear_operator
        implicit none
 
        integer, parameter :: n=4
        integer :: i, j
 
        real (dp) , dimension(n,1)   ::   x

                                         ! DP matrix (n,n)
        real (dp) , dimension(n,n)   ::   z
        logical   , dimension(n,n)   ::   z1
 
 
    	x = reshape( (/ 11,12,13,14 /) , (/4,1/) )

        z = reshape( (/ 11,12,13,14,                    &
                        12,22,23,24,                    &
                        13,23,33,34,                    &
                        14,24,34,44 /)                  &
                       ,(/4,4/) )
 
 ! Start of program

		z1 = (z <= 13)
 
		print *,""
		do i=1,n
			print *, (z1(i,j),j=1,n)
		enddo
					! element by element

		z1 = (z == z)
 
		print *,""
		do i=1,n
			print *, (z1(i,j),j=1,n)
		enddo

					! useless but interesting code
		z1 = prdt(x,x)
		z1 = (z == z1)
 
		print *,""
		do i=1,n
			print *, (z1(i,j),j=1,n)
		enddo

end