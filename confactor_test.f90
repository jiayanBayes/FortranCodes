module mysubs
	implicit none
	contains

		subroutine truncated_normal(u, v, a, ind, dim, x)
			! Random generator from truncated normal
			! The ouput is the a vector with dimension of dim
			! 05/03/2003 Jia Yan, ChoiceStream Inc.
			use matforF90
			use linear_operator
			use numerical_libraries
			implicit none

			! Decalre inputs and output
			integer, intent(in) :: dim
			real(dp), dimension(dim), intent(in) :: u		!mean of the normal distribution
			real(dp), intent(in) :: v						!variance of the normal distribution
			real(dp), intent(in) :: a						!truncation point
			logical,  intent(in)  :: ind					!indicator, truncated below if true and truncated above if false
			real(dp), dimension(dim), intent(out) :: x		!the draw from truncated normal


			! Declare local variables
			real(dp), dimension(dim) :: point						!the point which standard normal c.d.f is evaluated at
			real(dp), dimension(dim) :: cdf							 
			real(dp), dimension(dim) :: uniform						!uniform random number
			real(dp), dimension(dim) :: temp
			real(dp) :: i 
			
			! Body of the procedure
			point = (a-u)/sqrt(v)
			cdf = 0.
			do i = 1, dim
				cdf(i) = dnordf(point(i))
			end do 
			
			call drnun(dim,uniform)
	
			if (ind) then
				temp = cdf + uniform * (1 - cdf)
			else
				temp = uniform * cdf	
			end if
			
			x = 0.
			do i = 1, dim
				x(i) = u(i) + sqrt(v) * dnorin(temp(i))
			end do
	end subroutine truncated_normal

end module

program lcfa
	! This program implements Bayesian factor analysis for binary data. 
	! For details, see LCFA.doc by Jia Yan
	! 04/02/2003 Jia Yan, ChoiceStream Inc.
	use matforF90
	use mysubs
	use numerical_libraries
	use linear_operator
	implicit none
	
	! Decalre local variables
	integer, parameter :: Nitems = 1331					 !number of items 
	integer, parameter :: Nattrs = 760					 !dimenison of attributes of each item
	integer, parameter :: Nfactors = 10 				 !desired number of factors
	integer, parameter :: ifinput = 15
	
	integer :: i
	real(dp), parameter :: tpoint = 0.
	real(dp), parameter :: v = 1.
	real(dp), dimension(Nfactors, 1) :: mu 
	real(dp), dimension(Nfactors, 1) :: tn
	real(dp), dimension(Nfactors, 1) :: tnp   
	character(20), dimension(Nitems, 1) :: id 	
	integer, dimension(Nitems, Nattrs) :: attributes 

	character(80) :: inputfile = "biodata.asc"				!inputed data set
	
	! Body part begins from here
	open( ifinput,      FILE=inputfile,			   STATUS='UNKNOWN')

	! Read data. The data format should be: the 1st column is the item id, and from the second column to the 
	! end are the attributes 
	do i = 1, Nitems  
		read(ifinput,*) id(i,1), attributes(i, :)
		!if (i==1) then
		!	call prnt(attributes(i, :))
		!end if
	end do
	
	mu = 1. 
	call truncated_normal(mu, v, tpoint, .TRUE., Nfactors, tnp)
	call truncated_normal(mu, v, tpoint, .FALSE., Nfactors, tn)
	call prnt(tn)
	!call prnt(tn)
		
	
			!Calculating the Schwarz or Bayesian Information Criterion (BIC)
		!likelihood = 0. 
		!do i = 1, Nitems
		!	part1 = 1
		!	do j = 1, Nattrs
		!		xb = -1 * alphha_final(j,1) - dot_product(betta_final(j, :), f_final(i, :))
		!		if (attributes(i,j) == 1) then
		!			ncf = 1 - dnordf(xb)
		!		else
		!			ncf = dnordf(xb)
		!		end if
		!		part1 = part1 * ncf  
		!	end do
		!	part2 = product( (1 / sqrt(2*pi)) * exp( -1 * (f_final(i, :)**2) / 2 ) )
		!	likelihood = likelihood + log(part1 * part2)
		!end do
		!n = Nitems
		!pk = Nattrs * Nfactors - Nfactors * (Nfactors - 1) / 2
		!bic = 2 * likelihood + log(n) * pk  
	!end if
	!print *, bic
	!call prnt(attributes)
end program

