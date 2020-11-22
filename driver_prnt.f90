! driver for prnt


    use matforf90
    implicit none

	integer, parameter :: n=4

							! Integer vector (n)
	integer , dimension(n)   ::  ix

							! Integer vector (n,1)
	integer , dimension(n,1) ::   iy

							! Integer matrix (n,n)
	integer , dimension(n,n) ::   iz


							! DP rank-1 vector
	real (dp) , dimension(n)     ::   x

							! DP rank-2 matrix (n,n)
	real (dp) , dimension(n,1)   ::   y
	real (dp) , dimension(n,n)   ::   z

							! Character vectors
	character(nmaxchar), dimension(n) :: namevars


! Initialisation

  	 ix = reshape( (/ 1,2,3,4 /) , (/4/) )
     iy = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

 	 iz = reshape( (/ 11,12,13,14,			&
				 	  21,22,23,24,			&
				 	  31,32,33,34,			&
					  41,42,43,44 /)		&
				  ,(/4,4/) )

   	 x = reshape( (/ 1,2,3,4 /) , (/4/) )
     y = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

 	 z = reshape( (/ 11,12,13,14,			&
					21,22,23,24,			&
					31,32,33,34,			&
					41,42,43,44 /)			&
				  ,(/4,4/) )

									! for some reason, the "" must have the same
									! number of characters when inputed like this.
									! Individual elements may be changed later on.

	namevars  = (/"Hello", "How  ", "are  ", "you. "/)


! Start of program


    print *, "rank-1 ix vector"; call prnt(ix)
    print *, "rank-2 iy vector"; call prnt(iy)
    print *, "izmatrix"; call prnt(iz)

    print *, "rank-1 x vector"; call prnt(x)
    print *, "rank-2 y vector"; call prnt(y)
    print *, "zmatrix"; call prnt(z)

	print *, "character vector"; call prnt(namevars)

end

