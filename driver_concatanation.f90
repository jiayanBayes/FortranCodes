!Driver concatanation

	use matforf90
	use linear_operator

 	implicit none

	integer, parameter :: n=4

							! Integer vector (n,1)
	integer , dimension(n,1) ::   iy
	integer , dimension(:,:), allocatable   ::   iy1

							! Integer matrix (n,n)
	integer , dimension(n,n) ::   iz
	integer , dimension(:,:), allocatable   ::   iz1
	integer , dimension(:,:), allocatable   ::   iz2

							! DP rank-2 matrix (n,n)
	real (dp) , dimension(n,1)   ::   y
	real (dp) , dimension(:,:), allocatable   ::   y1
	real (dp) , dimension(n,n)   ::   z
	real (dp) , dimension(:,:), allocatable   ::   z1
	real (dp) , dimension(:,:), allocatable   ::   z2

							! Character vectors
	character(nmaxchar), dimension(n) :: namevars
	character(nmaxchar), dimension(1) :: namevars1
	character(nmaxchar), dimension(:), allocatable :: namevars2


! Initialisation

     iy = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

 	 iz = reshape( (/ 11,12,13,14,			&
				 	  21,22,23,24,			&
				 	  31,32,33,34,			&
					  41,42,43,44 /)		&
				  ,(/4,4/) )

     y = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

 	 z = reshape( (/ 11,12,13,14,			&
					21,22,23,24,			&
					31,32,33,34,			&
					41,42,43,44 /)			&
				  ,(/4,4/) )

	namevars  = (/"Hello", "How  ", "are  ", "you. "/)
    namevars1 = (/"Very well."/)

! Start of code

				! Double precision

	print *, "z"; call prnt(z)
	print *, "y"; call prnt(y)

				! Concatanation ~ = horizontal  | = vertical
				! result of z~y
	allocate(z1(n,n+1))
	z1 = hcon(z,y)
	print *,"z~y "; call prnt(z1)
	z1 = 0
	z1 = z .hc. y
	print *,"z .hc. y "; call prnt(z1)
	deallocate(z1)

					! result of z|x

	allocate(y1(1,n))
	y1 = .t. y			! transposed of y

	allocate(z1(n+1,n))
	z1 = vcon(z,y1)
	print *,"z|y' "; call prnt(z1) 
	z1 = 0
	z1 = z .vc. y1
	print *,"z .vc. y' "; call prnt(z1)
	deallocate(z1)

					! four blocks

	allocate(z1(n+1,n+1))
	allocate(z2(1,1))		! rank-2 scalar

	z2 = 5

	z1  =      ( z    .hc. y    )     &
		  .vc. ( y1   .hc. z2)	

	print *,"z1 "; call prnt(z1)
	deallocate(z1)


				! Integer

	print *, "iz"; call prnt(iz)
	print *, "iy"; call prnt(iy)

				! Concatanation ~ = horiizontal  | = vertical
				! result of iz~iy
	allocate(iz1(n,n+1))
	iz1 = hcon(iz,iy)
	print *,"iz~iy "; call prnt(iz1)
	iz1 = 0
	iz1 = iz .hc. iy
	print *,"iz .hc. iy "; call prnt(iz1)
	deallocate(iz1)

					! result of iz|iy

	allocate(iy1(1,n))
	iy1 = .t. iy			! transposed of iy

	allocate(iz1(n+1,n))
	iz1 = vcon(iz,iy1)
	print *,"iz|iy' "; call prnt(iz1) 
	iz1 = 0
	iz1 = iz .vc. iy1
	print *,"iz .vc. iy' "; call prnt(iz1)
	deallocate(iz1)

					! four blocks

	allocate(iz1(n+1,n+1))
	allocate(iz2(1,1))		! rank-2 scalar

	iz2 = 5

	iz1  =      ( iz    .hc. iy    )     &
		   .vc. ( iy1   .hc. iz2)	

	print *,"iz1 "; call prnt(iz1)
	deallocate(iz1)

					! character vectors

	allocate(namevars2(n+1))
	namevars2 = namevars .vc. namevars1
	call prnt(namevars2)
	


end