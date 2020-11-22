!Driver for selif and delif

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

							! Character vectors
	character(nmaxchar), dimension(n) :: namevars
	character(nmaxchar), dimension(:), allocatable :: namevars1

	logical , dimension(n) ::  b
	integer i, nbr

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


	b=(ix == 1)   .or. (ix == 2)
	print *, "b vector", b
	
	print *, "rank-1 integer" 		
 				! Delif
    nbr    = count( .not. b)            ! counts the number of false elements in b
    if (allocated(ix1) ) deallocate(ix1)
	allocate( ix1(nbr) )                ! allocate ix1 to be subvector from ix.
	ix1  = delif(ix,  b)                ! delete row if b is true
	print *, "Delif:"
	print *, "ix  vector", ix
	print *, "ix1 subvector"; call prnt(ix1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(ix1) ) deallocate(ix1)
	allocate( ix1(nbr) )                ! allocate ix1 to be subvector from ix.
	ix1  = selif(ix,  b)                ! select row if b is true
	print *, "Selif:"
	print *, "ix1 subvector"; call prnt(ix1)


	print *," "
	print *, "rank-2 integer" 		
 				! Delif
    nbr    = count( .not. b)            ! counts the number of false elements in b
    if (allocated(iz1) ) deallocate(iz1)
	allocate( iz1(nbr,cols(iz)) )       ! allocate iz1 to be submatrix from iz.
	iz1  = delif(iz,  b)                ! delete row if b is true
	print *, "Delif:"
	print *, "iz  matrix"; call prnt(iz)
	print *, "iz1 submatrix"; call prnt(iz1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(iz1) ) deallocate(iz1)
	allocate( iz1(nbr,cols(iz)) )       ! allocate iz1 to be submatrix from iz.
	iz1  = selif(iz,  b)                ! select row if b is true
	print *, "Selif:"
	print *, "iz1 submatrix"; call prnt(iz1)

				! DP

	print *," "
	print *, "rank-1 dp" 		
 				! Delif
    nbr    = count( .not. b)            ! counts the number of false elements in b
    if (allocated(x1) ) deallocate(x1)
	allocate( x1(nbr) )                 ! allocate ix1 to be subvector from ix.
	x1  = delif(x,  b)                  ! delete row if b is true
	print *, "Delif:"
	print '(a10,(10f8.2))',  "xvector ", x
	print *, "x1 subvector"; call prnt(x1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(x1) ) deallocate(x1)
	allocate( x1(nbr) )                 ! allocate ix1 to be subvector from ix.
	x1  = selif(x,  b)                  ! select row if b is true
	print *, "Selif:"
	print *, "x1 subvector"; call prnt(x1)


	print *," "
	print *, "rank-2 dp" 		
 				! Delif
    nbr    = count( .not. b)            ! counts the number of false elements in b
    if (allocated(z1) ) deallocate(z1)
	allocate( z1(nbr,cols(iz)) )        ! allocate iz1 to be submatrix from iz.
	z1  = delif(z,  b)                  ! delete row if b is true
	print *, "Delif:"
	print *, "z  matrix"; call prnt(z)
	print *, "z1 submatrix"; call prnt(z1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(z1) ) deallocate(z1)
	allocate( z1(nbr,cols(iz)) )        ! allocate iz1 to be submatrix from iz.
	z1  = selif(z,  b)                  ! select row if b is true
	print *, "Selif:"
	print *, "z1 submatrix"; call prnt(z1)

				! Character

	print *," "
	print *, "rank-1 character vector" 		
 				! Delif
    nbr    = count( .not. b)            ! counts the number of false elements in b
    if (allocated(namevars1) ) deallocate(namevars1)
	allocate( namevars1(nbr) )          ! allocate namevars1 to be subvector from namevars.
	namevars1  = delif(namevars,  b)    ! delete row if b is true
	print *, "Delif:"
	print '(a10,(10A8))',  "namevars : ", namevars
	print *, "namevars1 subvector"; call prnt(namevars1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(namevars1) ) deallocate(namevars1)
	allocate( namevars1(nbr) )          ! allocate namevars1 to be subvector from namevars.
	namevars1  = selif(namevars,  b)    ! select row if b is true
	print *, "Selif:"
	print *, "namevars1 subvector"; call prnt(namevars1)

end