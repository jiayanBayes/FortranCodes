!Driver random numbers

	use matforf90
	use linear_operator
	implicit none

	integer, parameter :: n = 4
    integer istart, incr, iend, nbrterms

							!  integer objects
	integer , dimension(n)    ::   ix
	integer , dimension(:),   allocatable   ::   ix1

							!  DP objects
	real (dp) , dimension(n)  ::   x
	real (dp) , dimension(:),   allocatable  ::   x1
	real (dp) , dimension(:,:), allocatable  ::   y


						! seqa, seqam, seqamax return in rank-1 vectors
	
	istart = 1
	incr   = 2
	iend   = 10
	nbrterms = 10

						! seqamax
	allocate(ix1(1+(iend-istart)/incr))	
	ix1 = seqamax(istart,incr,iend)
	print *, "sequence from istart to iend increasing by incr"; call prnt(ix1) 

						! seqa
	allocate(y(nbrterms,1))
	y(:,1) = seqa(istart,incr,nbrterms)
	print *,"sequence starting at istart, incrasing by incr and containing nbrterms terms&
	& stored in rank-2 vector"
	call prnt(y)


	y(:,1) = seqam(istart,incr,nbrterms)
	print *,"multiplicative sequence starting at istart, incrasing by incr & 
			&and containing nbrterms terms stored in rank-2 vector"
	call prnt(y)

end

