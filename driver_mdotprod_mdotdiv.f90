!Driver mdotprod mdotdiv

	use matforf90
	use linear_operator

 	implicit none

	integer, parameter :: n=4

	integer ,  dimension(n)  :: index
	integer ,  dimension(:), allocatable  :: index1

							! Integer vector (n)
	integer , dimension(n)   ::  ix
	integer , dimension(:), allocatable   ::   ix1

							! Integer vector (n,1)
	integer , dimension(n,1) ::   iy
	integer , dimension(:,:), allocatable   ::   iy1

							! Integer matrix (n,n)
	integer , dimension(n,n) ::   iz
	integer , dimension(:,:), allocatable   ::   iz1


							! DP rank-1 vector
	real (dp) , dimension(n)     ::   x
	real (dp) , dimension(:), allocatable     ::   x1

							! DP rank-2 matrix (n,n)
	real (dp) , dimension(n,1)   ::   y
	real (dp) , dimension(:,:), allocatable   ::   y1
	real (dp) , dimension(n,n)   ::   z
	real (dp) , dimension(:,:), allocatable   ::   z1


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


! Start of program
	
	print *, "Double precision version"
	print *, "x"; call prnt(x)
	print *, "z"; call prnt(z)
				
	allocate(z1(n,n), x1(n), y1(n,1))

				!MDOTPROD  dp

	print *,"dp version"
	print *,"MDOTPROD"

	z1 = mdotprod(z,z)
	print *,"z.*z"; call prnt(z1)

	z1 = mdotprod(z,y)
	print *,"z.*y"; call prnt(z1)

	z1 = mdotprod(z, .t. y)
	print *,"z.*y'"; call prnt(z1)

	y1 = mdotprod(y,y)
	print *,"y.*y"; call prnt(y1)

				! with rank 1
	print *,"rank-1 x"
	z1 = mdotprod(z,x)
	print *,"z.*x"; call prnt(z1)

		!linear operator
	print *, "using linear operator .prod."
	z1 = z .prod. x
	print *,"z.prod.x"; call prnt(z1)


	print *,"MDODIV"

			!MDOTDIV dp

	z1 = mdotdiv(z,z)
	print *,"z./z"; call prnt(z1)

	z1 = mdotdiv(z,y)
	print *,"z./y"; call prnt(z1)

	z1 = mdotdiv(z, .t. y)
	print *,"z./y'"; call prnt(z1)

	y1 = mdotdiv(y,y)
	print *,"y./y"; call prnt(y1)

				! with rank 1
	print *,"rank-1 x"
	z1 = mdotdiv(z,x)
	print *,"z./x"; call prnt(z1)

		!linear operator
	print *, "using linear operator .div."
	z1 = z .div. x
	print *,"z.div.x"; call prnt(z1)

				! In place calculation

	z = mdotdiv(z,x)
	print *,"z./x in place"; call prnt(z)

    z = mdotdiv(z,z)
	print *,"z./z in place"; call prnt(z)


	print *,"integer version"

	allocate(iz1(n,n), ix1(n), iy1(n,1))

		 !MDOTPROD  int

	print *,"MDOTPROD "

	iz1 = mdotprod(iz,iz)
	print *,"iz.*iz"; call prnt(iz1)

	iz1 = mdotprod(iz,iy)
	print *,"iz.*iy"; call prnt(iz1)

	iz1 = mdotprod(iz, .t. iy)
	print *,"iz.*iy'"; call prnt(iz1)

	iy1 = mdotprod(iy,iy)
	print *,"iy.*iy"; call prnt(iy1)

				! with rank 1
	print *,"rank-1 ix"
	iz1 = mdotprod(iz,ix)
	print *,"iz.*ix"; call prnt(iz1)

		!linear operator
	print *, "using linear operator .prod."
	iz1 = iz .prod. ix
	print *,"iz.prod.ix"; call prnt(iz1)


	print *,"MDOTDIV"

			!MDOTDIV in

	iz1 = mdotdiv(iz,iz)
	print *,"iz./iz"; call prnt(iz1)

	iz1 = mdotdiv(iz,iy)
	print *,"iz./iy"; call prnt(iz1)

	iz1 = mdotdiv(iz, .t. iy)
	print *,"iz./iy'"; call prnt(iz1)

	iy1 = mdotdiv(iy,iy)
	print *,"iy./iy"; call prnt(iy1)

				! with rank 1
	print *,"rank-1 ix"
	iz1 = mdotdiv(iz,ix)
	print *,"iz./ix"; call prnt(iz1)

		!linear operator
	print *, "using linear operator .div."
	iz1 = iz .div. ix
	print *,"iz.div.ix"; call prnt(iz1)


				! In place calculation

	iz = mdotdiv(iz,ix)
	print *,"iz./ix in place"; call prnt(iz)

    iz = mdotdiv(iz,iz)
	print *,"iz./iz in place"; call prnt(iz)


end