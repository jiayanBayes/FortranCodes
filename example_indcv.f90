!Driver selected objects


	use matforF90
	use linear_operator
	implicit none

	integer, parameter :: n=4
	integer :: i
							! Character vectors
	character(nmaxchar), dimension(n) :: namevars
	character(nmaxchar), dimension(3) ::  namevars1
	character(nmaxchar), dimension(:), allocatable ::	 namevars2

	integer, dimension(:), allocatable :: index

	namevars  = (/"Hello", "How  ", "are  ", "you. "/)
	namevars1 = (/"How", "are", "nos"/)


!indcv

	allocate(index(rows(namevars1)-1))
	index = indcv(namevars1(1:2), namevars)
	print *,"indcv"; call prnt(index)

	deallocate(index)
	allocate(index(rows(namevars1)))
	index = indcv(namevars1, namevars)
	print *,"indcv"; call prnt(index)


end

