module mysubs
	implicit none
	contains

		subroutine truncated_normal_vector(u, v, a, ind, dim, x)
			! Random generator from truncated normal
			! The ouput is the a vector with dimension of dim
			! 05/03/2003 Jia Yan, ChoiceStream Inc.
			use matforF90
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
		end subroutine truncated_normal_vector

		subroutine truncated_normal_univariate(u, v, a, ind, x)
			! Random generator from truncated normal
			! 03/06/2003 Jia Yan, ChoiceStream Inc.
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and output
			real(dp), intent(in) :: u				!mean of the normal distribution
			real(dp), intent(in) :: v				!variance of the normal distribution
			real(dp), intent(in) :: a				!truncation point
			logical,  intent(in)  :: ind			!indicator, truncated below if true and truncated above if false
			real(dp), intent(out) :: x				!the draw from truncated normal


			! Declare local variables
			real(dp) :: point						!the point which standard normal c.d.f is evaluated at
			real(dp) :: cdf							 
			real(dp) :: uniform						!uniform random number
			real(dp) :: temp 
			
			! Body of the procedure
			point = (a-u)/sqrt(v)
			cdf = dnordf(point)
			call drnun(1,uniform)
	
			if (ind) then
				temp = cdf + uniform * (1 - cdf)
			else
				temp = uniform * cdf	
			end if
		
			x = u + sqrt(v) * dnorin(temp)
	
	end subroutine truncated_normal_univariate

end module

program lcfa
	! This program implements Bayesian factor analysis for binary data. 
	! For details, see LCFA.doc by Jia Yan
	! 04/02/2003 Jia Yan, ChoiceStream Inc.
	use matforF90
	use mysubs
	use numerical_libraries
	implicit none
	
	! Decalre local variables
	integer, parameter :: Nitems = 1331					 !number of items 
	integer, parameter :: Nattrs = 760					 !dimenison of attributes of each item
	integer, parameter :: Nfactors = 5 				 !desired number of factors
	integer, parameter :: ifinput = 15
	integer, parameter :: ifloading = 16
	integer, parameter :: iffactor = 18 
	integer, parameter :: ifalpha = 19
	integer, parameter :: ifbic = 20
	integer, parameter :: Nloops = 8000					 !number of Gibbs draws
	integer, parameter :: truncation = 1000
	
	logical, parameter :: check = .TRUE.

	real(dp), parameter :: tpoint = 0.
	real(dp), parameter :: v = 1.
	real(dp), parameter :: pi = 3.1415926
	integer :: i, iloops, j
	real(dp) :: b0, invvb, tn, xb, s, likelihood, ncf, part1, part2, n, pk, bic	 

	character(20), dimension(Nitems, 1) :: id 	
	integer, dimension(Nitems, Nattrs) :: attributes 
	real(dp), dimension(Nattrs,1) :: alphha_ini, alphha, alphha_final
	real(dp), dimension(Nfactors + 1, 1) :: delta0
	real(dp), dimension(Nfactors + 1, Nfactors + 1) :: invvdelta 
	real(dp), dimension(Nitems, Nfactors) :: f_ini, f, f_final
	real(dp), dimension(Nattrs, Nfactors) :: betta_ini, betta, betta_final
	real(dp), dimension(Nfactors, Nfactors) :: BB
	real(dp), allocatable, dimension(:, :) :: I_star, I_temp, I_i, capD, litd, randn, f_temp, &
	 										  phi0, gphi, x, phik, y_temp, cons, xx, deltak, &
											  tmat1, tmat2
			  
	
	character(80) :: inputfile = "biodata.asc"							!inputed data set
	character(80) :: outfile_loading = "loading_bio_5.asc"				!result for loading matrix
	character(80) :: outfile_factor = "factor_bio_5.asc"					!result for factor matrix
	character(80) :: outfile_alpha = "alpha_bio_5.asc"					!result for the intercepts
	character(80) :: outfile_bic = "bic_5.asc"

	! Body part begins from here
	open( ifinput,    FILE=inputfile,       STATUS='UNKNOWN')
	open( ifloading,  FILE=outfile_loading, STATUS='UNKNOWN', RECL=8000)
	open( iffactor,   FILE=outfile_factor,  STATUS='UNKNOWN', RECL=8000)
	open( ifalpha,    FILE=outfile_alpha,   STATUS='UNKNOWN', RECL=8000)
	open( ifbic,      FILE=outfile_bic,     STATUS='UNKNOWN', RECL=8000)
	! Read data. The data format should be: the 1st column is the item id, and from the second column to the 
	! end are the attributes 
	do i = 1, Nitems  
		read(ifinput,*) id(i,1), attributes(i, :)
		!if (i==1) then
		!	call prnt(attributes(i, :))
		!end if
	end do
		
	if (check) then
		! Initialize the priors for hyperparameters of the two blocks: theta2 and theta3
		b0	  = 0.									!prior mean of the diaganol elements 				
		invvb = 0.25 								!prior precision of the diaganol elements	 
		delta0 = 0.									!prior mean of the last P-q rows 
		invvdelta = 0.25 * eye(Nfactors + 1)		!prior precision of the last P-q rows
	
		f_ini = 0.
		! Set initial values for loading matrix and mean vector
		betta_ini = 0.
		alphha_ini = 0.

		! Gibbs Sampling
		f_final = 0.				!f_final stores the final results of factors 	
		alphha_final = 0.			!alphha_final stores the final results of population mean
		betta_final = 0.			!betta_final stores the final results of loading matrix 			
		do iloops = 1, Nloops
			!Update parameters
			f = f_ini
			betta = betta_ini
			alphha = alphha_ini	
			!calculating values will be used later
			BB = betta .tx. betta
		
			allocate( I_star(Nitems, Nattrs) ) 
			do i = 1, Nitems
				 !ith row of latent index
				do j = 1, Nattrs
					xb = alphha(j,1) + dot_product( betta(j, :), f(i, :) ) 
					if (attributes(i,j) == 1) then
						call truncated_normal_univariate(xb, v, tpoint, .TRUE. ,tn)
						I_star(i,j) = tn
					else
						call truncated_normal_univariate(xb, v, tpoint, .FALSE., tn)
						I_star(i,j) = tn
					end if
				end do
			end do
		
			!I_star = attributes
			!Draw latent indicator and factors
			do i = 1, Nitems
				!ith row of latent factors
				allocate( capD(Nfactors, Nfactors) )
				allocate( I_temp(1, Nattrs) )
				allocate( I_i(Nattrs, 1) )
				allocate( litd(Nfactors, 1) )
				allocate( randn(Nfactors, 1) )
				allocate( tmat1(Nfactors, 1) )
				allocate( tmat2(Nfactors, 1) )
				allocate( f_temp(1, Nfactors) )
				capD = invpd( BB + eye(Nfactors) )
				I_temp(1, :) = I_star(i, :)
				I_i = transpose( I_temp )
				I_i = I_i - alphha
				litd = betta .tx. I_i
				call drnnor(Nfactors, randn)
				tmat1 = capD .x. litd
				tmat2 = chol(capD) .tx. randn
				f_temp = transpose(tmat1 + tmat2) 
				f_ini(i,:) = f_temp(1, :) 
				deallocate(randn)
				deallocate(f_temp)  
				deallocate(capD)
				deallocate(I_temp)
				deallocate(I_i)
				deallocate(litd)
				deallocate(tmat1)
				deallocate(tmat2)
			end do  	  
		
			!Draw the first and second blocks of loading matrix
			allocate( I_temp(Nitems, 1) )
			allocate( y_temp(Nitems, 1) ) 
			do i = 1, Nfactors
				allocate( phi0(i, 1) )
				allocate( gphi(i, i) ) 
				allocate( x(Nitems, i) )
				allocate( capD(i, i) )
				allocate( litd(i, 1) )
				allocate( randn(i, 1) )
				allocate( phik(i, 1) )
				phi0 = 0.					!prior mean of the first block
				gphi = 0.25 * eye(i)		!prior precision of the first block
				if (i == 1) then
					x = 1.
				else
					x(:, 1) = 1.
					do j = 1, i - 1
						x(:, j + 1) = f_ini(:, j)						
					end do	
				end if
				I_temp(:, 1) = I_star(:, i) - betta(i, i) * f_ini(:, i)
				capD = x .tx. x
				capD = capD + gphi
				capD = invpd(capD)
				litd = (x .tx. I_temp) + phi0
				call drnnor(i, randn)
				phik = chol(capD) .tx. randn
				phik = (capD .x. litd) + phik 
				alphha_ini(i, 1) = phik(1,1) 
				if (i >= 2 ) then
					do j = 1, i - 1
						betta_ini(i, j) = phik(j+1, 1)
					end do
				end if
			
				!Draw the diagonal elements: the second block of loading matrix 
				s = 1 / (dot_product( f_ini(:,i), f_ini(:,i)) + invvb)
				I_temp = x .x. phik 
				y_temp(:, 1) = I_star(:, i) - I_temp(:, 1)
				xb = s * ( dot_product(f_ini(:, i), y_temp(:, 1) ) + b0 )
				call truncated_normal_univariate(xb, s, tpoint, .TRUE., tn) 
				betta_ini(i,i) = tn
			
				deallocate(phi0)
				deallocate(gphi)
				deallocate(x)
				deallocate(capD)
				deallocate(litd)
				deallocate(randn) 
				deallocate(phik)
			end do
			deallocate(I_temp) 
			deallocate(y_temp)
	
			!Draw the third block of loading matrix
			allocate( I_temp(Nitems, 1) )
			allocate( x(Nitems, Nfactors + 1) ) 
			allocate( xx(Nfactors + 1, Nfactors + 1) )
			allocate( cons(Nitems, 1) )
			allocate( capD(Nfactors + 1, Nfactors + 1) )
			allocate( litd(Nfactors + 1, 1) ) 
			allocate( randn(Nfactors + 1, 1) )
			allocate( tmat1(Nfactors + 1, 1) )
			allocate( tmat2(Nfactors + 1, 1) ) 
			allocate( deltak(1, Nfactors + 1) )
			cons = 1.
			x = cons .hc. f_ini
			xx = x .tx. x
			capD = invpd( xx + invvdelta)
			deallocate(xx)
			do i = 1, Nattrs-Nfactors
				I_temp(:,1) = I_star(:, Nfactors + i)
				litd = (x .tx. I_temp) + delta0
				call drnnor(Nfactors + 1, randn)
				tmat1 = capD .x. litd
				tmat2 = chol(capD) .tx. randn
				deltak = transpose( tmat1 + tmat2 )
				alphha_ini(Nfactors + i, 1) = deltak(1,1)
				betta_ini(Nfactors + i, :) = deltak(1, 2:Nfactors + 1 )
			end do
			!Now we have obtained all elements in loading matrix
			deallocate(I_temp)
			deallocate(x)
			deallocate(cons)
			deallocate(capD)
			deallocate(litd)
			deallocate(randn)
			deallocate(tmat1)
			deallocate(tmat2)
			deallocate(deltak)
			deallocate(I_star)
			!Storing results 
			if (iloops > truncation) then
				f_final = (f_final + f_ini) 
				betta_final = (betta_final + betta_ini) 
				alphha_final = (alphha_final + alphha_ini) 
			end if
			write(*,*) iloops
		end do
		alphha_final = alphha_final / (Nloops - truncation)
		betta_final = betta_final / (Nloops - truncation)
		f_final = f_final / (Nloops - truncation)
		call prnt(alphha_final) 
		call prnt(betta_final)

		!Write results onto hard disk
		do i = 1, Nitems
			write(iffactor,*) id(i,1), achar(9), (f_final(i,j), achar(9), j = 1, Nfactors -1), f_final(i, Nfactors)
			!100 format(1x, 20F20.4)
		end do

		do i = 1, Nattrs
			write(ifloading,*) (betta_final(i,j), achar(9), j = 1, Nfactors -1), betta_final(i, Nfactors)
			write(ifalpha,*)   alphha_final(i,1)
			!200 format(1x, 20F20.4)
		end do

		!Calculating the Schwarz or Bayesian Information Criterion (BIC)
		likelihood = 0. 
		do i = 1, Nitems
			part1 = 1
			do j = 1, Nattrs
				xb = -1 * alphha_final(j,1) - dot_product(betta_final(j, :), f_final(i, :))
				if (attributes(i,j) == 1) then
					ncf = 1 - dnordf(xb)
				else
					ncf = dnordf(xb)
				end if
				part1 = part1 * ncf  
			end do
			part2 = product( (1 / sqrt(2*pi)) * exp( -1 * (f_final(i, :)**2) / 2 ) )
			likelihood = likelihood + log(part1 * part2)
		end do
		n = Nitems
		pk = Nattrs * Nfactors - Nfactors * (Nfactors - 1) / 2
		bic = 2 * likelihood + log(n) * pk  
	end if
	write(ifbic,*) bic
	!call prnt(attributes)
end program

