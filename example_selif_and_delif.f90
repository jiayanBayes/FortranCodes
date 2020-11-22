!Example for selif and delif

	use matforf90
	use linear_operator

 	implicit none

	integer, parameter :: n=4

	integer ,  dimension(n)  :: index
	integer ,  dimension(:), allocatable  :: index1

							! DP rank-1 vector
	real (dp) , dimension(n)     ::   x
	real (dp) , dimension(:), allocatable     ::   x1

							! DP rank-2 matrix (n,n)
	real (dp) , dimension(n,n)   ::   z
	real (dp) , dimension(:,:), allocatable   ::   z1

							! Character vectors
	character(nmaxchar), dimension(n) :: namevars
	character(nmaxchar), dimension(:), allocatable :: namevars1

	logical , dimension(n) ::  b
	integer i, nbr

! Initialisation

 
   	 x = reshape( (/ 1,2,3,4 /) , (/4/) )
  
 	 z = reshape( (/ 11,12,13,14,			&
					21,22,23,24,			&
					31,32,33,34,			&
					41,42,43,44 /)			&
				  ,(/4,4/) )

	namevars  = (/"Hello", "How  ", "are  ", "you. "/)


! Start of program


	b=(x == 1)   .or. (x == 2)
	print *, "b vector", b
	
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
	allocate( z1(nbr,cols(z)) )         ! allocate iz1 to be submatrix from iz.
	z1  = delif(z,  b)                  ! delete row if b is true
	print *, "Delif:"
	print *, "z  matrix"; call prnt(z)
	print *, "z1 submatrix"; call prnt(z1)


 				! Selif
    nbr    = count( b)                  ! counts the number of true elements in b
    if (allocated(z1) ) deallocate(z1)
	allocate( z1(nbr,cols(z)) )         ! allocate iz1 to be submatrix from iz.
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