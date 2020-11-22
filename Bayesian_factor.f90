module mysubs
	implicit none
	contains

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
	integer, parameter :: Nitems = 10000					 !number of items 
	integer, parameter :: Nattrs = 7					 !dimenison of attributes of each item
	integer, parameter :: Nfactors = 3	 				 !desired number of factors
	integer, parameter :: ifinput = 15
	integer, parameter :: ifloading = 16
	integer, parameter :: iffactor = 18 
	integer, parameter :: ifalpha = 19
	integer, parameter :: ifbic = 20
	integer, parameter :: Nloops = 1000					 !number of Gibbs draws
	integer, parameter :: truncation = 500
	
	logical, parameter :: check = .TRUE.

	real(dp), parameter :: tpoint = 0.
	real(dp), parameter :: v = 1.
	integer :: i, iloops, j
	real(dp) :: b0, invvb, tn, xb, s, rho, R, deg, chi2, scale 

	character(20), dimension(Nattrs, 1) :: id 	
	real(dp), dimension(Nitems, Nattrs) :: attributes 
	real(dp), dimension(Nattrs,1) :: alphha_ini, alphha, alphha_final, ommega, ommega_ini, ommega_final
	real(dp), dimension(Nfactors + 1, 1) :: delta0
	real(dp), dimension(Nfactors + 1, Nfactors + 1) :: invvdelta 
	real(dp), dimension(Nitems, Nfactors) :: f_ini, f, f_final
	real(dp), dimension(Nattrs, Nfactors) :: betta_ini, betta, betta_final, betta_temp 
	real(dp), dimension(Nfactors, Nfactors) :: BB
	real(dp), allocatable, dimension(:, :) :: I_temp, I_i, capD, litd, randn, f_temp, &
	 										  phi0, gphi, x, phik, y_temp, cons, xx, deltak, &
											  tmat1, tmat2, epsi, datta
			  
	
	character(80) :: inputfile = "gene.asc"							!inputed data set
	character(80) :: outfile_loading = "loading.asc"				!result for loading matrix
	character(80) :: outfile_factor = "factor.asc"					!result for factor matrix
	character(80) :: outfile_alpha = "alpha.asc"					!result for the intercepts
	character(80) :: outfile_bic = "bic.asc"

	! Body part begins from here
	open( ifinput,      FILE=inputfile,			   STATUS='UNKNOWN')
	open( ifloading,    FILE=outfile_loading,      STATUS='UNKNOWN', RECL=8000)
	open( iffactor,     FILE=outfile_factor,       STATUS='UNKNOWN', RECL=8000)
	open( ifalpha,      FILE=outfile_alpha,        STATUS='UNKNOWN', RECL=8000)
	open( ifbic,        FILE=outfile_bic,          STATUS='UNKNOWN', RECL=8000)
	! Read data. The data format should be: the 1st column is the item id, and from the second column to the 
	! end are the attributes 
	allocate( datta(Nattrs,Nitems) )
	do i = 1, Nattrs  
		read(ifinput,*) id(i,1), datta(i, :)
	end do
	attributes = transpose(datta)
	deallocate(datta)
	
		
	if (check) then
		! Initialize the priors for hyperparameters of the two blocks: theta2 and theta3
		b0	  = 0.									!prior mean of the diaganol elements 				
		invvb = 0.25 								!prior precision of the diaganol elements	 
		delta0 = 0.									!prior mean of the last P-q rows 
		invvdelta = 0.25 * eye(Nfactors + 1)		!prior precision of the last P-q rows
		rho = 2
		R = 1

		f_ini = 0.
		! Set initial values for loading matrix and mean vector
		betta_ini = 0.
		alphha_ini = 0.

		! Gibbs Sampling
		f_final = 0.				!f_final stores the final results of factors 	
		alphha_final = 0.			!alphha_final stores the final results of population mean
		betta_final = 0.			!betta_final stores the final results of loading matrix 
		ommega_final = 0.			
		ommega_ini = 1.
		do iloops = 1, Nloops
			!Update parameters
			f = f_ini
			betta = betta_ini
			alphha = alphha_ini	
			ommega = ommega_ini
			!calculating values will be used later
			do i = 1, Nattrs
				do j = 1, Nfactors
					betta_temp(i,j) = betta(i,j) * ( 1 / ommega(i,1) )
				end do
			end do
			BB = betta_temp .tx. betta
			!Draw latent factors
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
				I_temp(1, :) = attributes(i, :)
				I_i = transpose( I_temp )
				I_i = I_i - alphha
				litd = betta_temp .tx. I_i
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
				I_temp(:, 1) = attributes(:, i) - betta(i, i) * f_ini(:, i)
				capD = x .tx. x
				capD = capD * (1/ommega(i,1))
				capD = capD + gphi
				capD = invpd(capD)
				litd = x .tx. I_temp
				litd = litd * (1/ommega(i,1))
				litd = litd + phi0
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
				s = 1 / (dot_product( f_ini(:,i), f_ini(:,i)) * ( 1/ommega(i,1) ) + invvb)
				I_temp = x .x. phik 
				y_temp(:, 1) = attributes(:, i) - I_temp(:, 1)
				xb = s * ( dot_product(f_ini(:, i), y_temp(:, 1) ) * ( 1/ommega(i,1) )+ b0 )
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
			do i = 1, Nattrs-Nfactors
				xx = x .tx. x
				xx = xx * (1/ommega(Nfactors+i, 1))
				capD = invpd(xx + invvdelta)	
				I_temp(:,1) = attributes(:, Nfactors + i)
				litd = x .tx. I_temp
				litd = litd * (1/ommega(Nfactors+i, 1))
				litd = litd + delta0
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
			deallocate(xx)

			! Draw ommega matrix
			allocate(epsi(Nitems,Nattrs))
			epsi = 0.
			do i = 1, Nitems
				 !ith row of latent index
				do j = 1, Nattrs
					xb = alphha_ini(j,1) + dot_product( betta_ini(j, :), f_ini(i, :) ) 
					epsi(i,j) = attributes(i,j) - xb
				end do
			end do
			
			do i = 1, Nattrs
				scale = dot_product(epsi(:,i), epsi(:,i))
				scale = rho * R + scale
				deg = rho + Nitems
				call drnchi(1, deg, chi2)
				ommega_ini(i,1) = scale / chi2
			end do
			deallocate(epsi)

			!Storing results 
			if (iloops > truncation) then
				f_final = (f_final + f_ini) 
				betta_final = (betta_final + betta_ini) 
				alphha_final = (alphha_final + alphha_ini) 
				ommega_final = ommega_final + ommega_ini
			end if
			write(*,*) iloops
		end do
		alphha_final = alphha_final / (Nloops - truncation)
		betta_final = betta_final / (Nloops - truncation)
		f_final = f_final / (Nloops - truncation)
		ommega_final = ommega_final / (Nloops - truncation)
		!call prnt(alphha_final)
		call prnt(betta_final)
		!call prnt(ommega_final)
		!Write results onto hard disk
		do i = 1, Nitems
			write(iffactor,*) (f_final(i,j), achar(9), j = 1, Nfactors -1), f_final(i, Nfactors)
			!100 format(1x, 20F20.4)
		end do

		do i = 1, Nattrs
			write(ifloading,*) (betta_final(i,j), achar(9), j = 1, Nfactors -1), betta_final(i, Nfactors)
			write(ifalpha,*)   alphha_final(i,1)
			!200 format(1x, 20F20.4)
		end do
	end if
end program

