! Example of subseting


    use matforf90
    implicit none

	integer, parameter :: n=4

	integer ,  dimension(n)  :: index
	integer ,  dimension(:), allocatable  :: index1

							! DP rank-1 vector
	real (dp) , dimension(n)     ::   x
	real (dp) , dimension(:), allocatable     ::   x1

							! DP rank-2 matrix 
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

	index  = (/(i,i=1,n)/)

					! We want to keep all rows where x = 1 or 2. 
	b=(x == 1.0)   .or. (x == 2.0)     ! logical vector with T's and F's

    nbr    = count(b)                  ! counts the number of true elements in b
	allocate( index1(nbr) )
	index1 = pack((/(i,i=1,n)/), b)            ! nbr-vector with positions of T's in b vector  				

	print *, "x vector"; call prnt(x)
	print *, "b vector", b
	print *, "index ",  index
	print *, "index1", index1

	allocate( x1(nbr) )                ! allocate x1 to be subvector from x.
	x1  = x(index1)
	print *, "x1 subvector"; call prnt(x1)

				! alternative way, using delif
	print *,"alternative way"
	b=(x == 1.0)   .or. (x == 2.0)     ! logical vector with T's and F's
    nbr    = count( .not. b)           ! counts the number of false elements in b
 
    if (allocated(x1) ) deallocate(x1)
	allocate( x1(nbr) )                ! allocate x1 to be subvector from x.
	x1  = delif(x, .not. b)
	print *, "x1 subvector"; call prnt(x1)



end

