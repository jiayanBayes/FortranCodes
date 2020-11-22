
module tools
	
	implicit none
	contains

		subroutine truncated_normal_2side_univariate (mu, v, a, b, x, iseed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure draws a random number from truncated normal,
			! which is truncated below a and above b.
			! The inputs are: 
			! mu -- mean of the normal distribution being truncated;
			! v -- variance of the normal distribution;
			! a, b -- vector of truncation points;
		
			! The ouput is x.
			! Jia Yan, Nov. 05, 2003 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and outputs
			real(dp), intent(in) :: mu, v, a, b
			integer, intent(in), optional :: iseed
			real(dp), intent(out) :: x

			! Decalre local variables
			real(dp) :: points1, points2, uniform, cdf1, cdf2, temp
		
			! Body begins from here
			if (present(iseed)) then
				call rnset(iseed)	
			end if
			call drnun(1, uniform)
		
			points1 = a - mu
			points1 = points1 / sqrt(v)
			points2 = b - mu
			points2 = points2 / sqrt(v) 
			x = 0.
		
			cdf1 = dnordf(points1)
			cdf2 = dnordf(points2)
			temp = cdf1 + uniform * (cdf2 - cdf1)
			x = mu + sqrt(v) * dnorin(temp)
		
		end subroutine truncated_normal_2side_univariate


		subroutine truncated_normal_2side_vec (dim, mu, v, a, b, x, iseed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure draws a vector of random number from truncated normal,
			! which is truncated below a and above b.
			! The 5 required inputs are: dim -- dimesnion of draws;
			! mu -- mean of the normal distribution being truncated;
			! v -- diagonal vector of covariance matrix of the normal distribution;
			! a, b -- vector of truncation points;
		
			! The ouput is x.
			! Jia Yan, Nov. 05, 2003 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and outputs
			integer, intent(in) :: dim
			real(dp), intent(in), dimension(dim,1) :: mu
			real(dp), intent(in), dimension(dim,1) :: v
			real(dp), intent(in), dimension(dim,1) :: a, b
			integer, intent(in), optional :: iseed
			real(dp), intent(out), dimension(dim,1) :: x

			! Decalre local variables
			real(dp) :: cdf1, cdf2, temp
			real(dp), dimension(dim,1) :: points1, points2, uniform
			integer :: i

			! Body begins from here
			if (present(iseed)) then
				call rnset(iseed)	
			end if
			call drnun(dim, uniform)
		
			points1 = a - mu
			points1 = points1 / sqrt(v)
			points2 = b - mu
			points2 = points2 / sqrt(v) 
			x = 0.
			do i = 1, dim
				cdf1 = dnordf(points1(i,1))
				cdf2 = dnordf(points2(i,1))
				temp = cdf1 + uniform(i,1) * (cdf2 - cdf1)
				x(i,1) = mu(i,1) + sqrt(v(i,1)) * dnorin(temp)
			end do
	
		end subroutine truncated_normal_2side_vec

	
		subroutine truncated_normal_1side_univariate(mu, v, a, ind, x, iseed)
			! Random generator from truncated normal, which 
			! is trncated below or above point a.
			! The inputs are:
			! 1. mu -- mean of the normal being truncated;
			! 2. v -- variance of the normal being truncated;
			! 3. a -- truncation point;
			! 4. ind -- 1 if truncated below; 0 if truncated above;
		
			! The output is x. 
			! 03/06/2003 Jia Yan, ChoiceStream Inc.
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and output
			real(dp), intent(in) :: mu				
			real(dp), intent(in) :: v				
			real(dp), intent(in) :: a				
			integer,  intent(in)  :: ind
			integer, intent(in), optional :: iseed 			
			real(dp), intent(out) :: x				


			! Declare local variables
			real(dp) :: point						
			real(dp) :: cdf							 
			real(dp) :: uniform						
			real(dp) :: temp 
			
			! Body of the procedure
			point = (a-mu)/sqrt(v)
			cdf = dnordf(point)

			if (present(iseed)) then
				call rnset(iseed)
			end if
		
			call drnun(1, uniform)
	
			if (ind == 1) then
				temp = cdf + uniform * (1 - cdf)
			else
				temp = uniform * cdf	
			end if
		
			x = mu + sqrt(v) * dnorin(temp)
	
		end subroutine truncated_normal_1side_univariate


		subroutine truncated_normal_1side_vec (dim, mu, v, a, ind, x, iseed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure draws a vector of random number from truncated normal.
			! The 5 required inputs are: dim -- dimesnion of draws;
			! mu -- mean of the normal distribution being truncated;
			! v -- diagonal vector of covariance matrix of the normal distribution;
			! a -- vector of truncation points;
			! ind -- 1 or 0, 1 if truncated below and 0 if truncated above;
			! The ouput is x.
			! Jia Yan, Sep. 27, 2003 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and outputs
			integer, intent(in) :: dim
			real(dp), intent(in), dimension(dim,1) :: mu
			real(dp), intent(in), dimension(dim,1) :: v
			real(dp), intent(in), dimension(dim,1) :: a
			integer, intent(in), dimension(dim,1) :: ind
			integer, intent(in), optional :: iseed
			real(dp), intent(out), dimension(dim,1) :: x

			! Decalre local variables
			real(dp) :: cdf, temp
			real(dp), dimension(dim,1) :: points, uniform
			integer :: i

			! Body begins from here
			if (present(iseed)) then
				call rnset(iseed)	
			end if
			call drnun(dim, uniform)
		
			points = a - mu
			points = points / sqrt(v)
			x = 0.
			do i = 1, dim
				cdf = dnordf(points(i,1))
				if (ind(i,1) == 1) then
					temp = cdf + uniform(i,1) * (1 - cdf)
				else
					temp = uniform(i,1) * cdf				
				end if	 
				x(i,1) = mu(i,1) + sqrt(v(i,1)) * dnorin(temp)
			end do
	
		end subroutine truncated_normal_1side_vec

		subroutine wishart_jia (dim, mat, deg, w, iseed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!This procudure generates wishart random draws w
			!with scale matrix mat, degree deg, dimension dim
			!based on seed seed
			!Jia Yan, Sep. 02, 2003
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre arguments
			integer, intent(in) :: dim
			integer, intent(in) :: deg
			integer, intent(in), optional :: iseed
			real(dp), intent(in), dimension(dim, dim) :: mat
			real(dp), intent(out), dimension(dim, dim) :: w

			! Decalre local variables
			real(dp), dimension(dim,1) :: randn, temp1
			real(sp), dimension(dim, dim) :: temp2
			integer :: i
	
			! Body begins from here
			if (present(iseed)) then
				call rnset(iseed)
			end if

			w = 0.
			do i = 1, deg
				call drnnor(dim, randn)
				temp1 = chol(mat) .tx. randn
				temp2 = temp1 .xt. temp1
				w = w + temp2
			end do

		end subroutine wishart_jia 

		subroutine dirichlet (dim, v, draw, seed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure generates random draws with dimension of d
			! from Dirichlet distribution
			! Jia Yan, Sep. 02, 2003
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
			use matforF90
			use numerical_libraries
			implicit none

			! Decalre inputs and outputs
			integer, intent(in) :: dim
			integer, intent(in), optional :: seed
			real(dp), intent(in), dimension(dim,1) :: v
			real(dp), intent(out), dimension(dim,1) :: draw
	
			! Declare local variables 
			integer :: i
			real(dp) :: s, g
			real(dp), dimension(dim ,1) :: gdraws
	
			! Body from here
			if (present(seed)) then
				call rnset(seed)
			end if
			gdraws = 0.
			do i = 1, dim
				s = v(i, 1)
				call drngam(1, s, g)
				gdraws(i ,1) = g
			end do 

			draw = 0.
			do i = 1, dim
				s = gdraws(i ,1)
				draw(i, 1) = s / sum(gdraws) 
			end do

		end subroutine dirichlet

		subroutine mnormal_density (dim, x, mu, var, fval)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure evaluates the density of multivariate
			! normal with dimension dim, mean vector mu, as well as
			! covariance matrix var, at x
			! Jia Yan, Sep. 04, 2003
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			use matforF90
			use numerical_libraries
			implicit none
		
			! Decalre inputs and outputs
			integer, intent(in) :: dim
			real(dp), intent(in), dimension(dim,1) :: mu, x
			real(dp), intent(in), dimension(dim,dim) :: var
			real(dp), intent(out) :: fval

			! Declqre local variables
			real(dp), dimension(dim,1) :: tempmat1
			real(dp), dimension(dim, dim) :: tempmat2, luvar
			real(dp), dimension(1,dim) :: tempmat3
			real(dp), dimension(1,1) :: tempmat4
			integer, dimension(dim) :: ipvt
			real(dp) :: tempscalar, det1, det2, det3

			! Body begins from here
			tempmat1 = x - mu
			tempmat2 = invpd(var)
			tempmat3 = tempmat1 .tx. tempmat2
			tempmat4 = tempmat3 .x. tempmat1

			tempscalar = tempmat4(1,1)

			call dlftrg(dim, var, dim, luvar, dim, ipvt)
			call dlfdrg(dim, luvar, dim, ipvt, det1, det2)
			det3 = det1 * (10**det2)
			fval = (1/sqrt(det3)) * exp(-0.5 * tempscalar)

		end subroutine mnormal_density
	
		subroutine multinomial (dim, n, p, draw, seed)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! This procedure generate random draws, draw, from
			! multinomial distribution with dimension, dim, parameters
			! n (independent trials) and p (probability vector)
			! Jia Yan, Sep. 09, 2003
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use matforF90
			use numerical_libraries
			implicit none
		
			! Declare inputs and outputs
			integer, intent(in) :: dim
			integer, intent(in) :: n
			real(dp), intent(in), dimension(dim,1) :: p
			integer, intent(in), optional :: seed
			integer, intent(out), dimension(1,dim) :: draw
		
			! Decalre local varaibales
			real(dp), allocatable, dimension(:,:) :: bp
			logical, dimension(dim,1) :: mask
			integer :: i, j, theta, thetasum, nelement
			real(dp) :: pelement, psum
				
			! Body begines from here
			mask = .TRUE.
			psum = 1
			thetasum = 0
			draw = 0
		
			if (present(seed)) then
				call rnset(seed)
			end if
		
			do i = 1, dim - 1
				nelement = n - thetasum
				pelement = p(i,1) / psum
				allocate(bp(nelement,1))
				call drnun(nelement, bp)
				theta = 0
				do j = 1, nelement
					if (bp(j,1) <= pelement) then
							theta = theta + 1
					end if
				end do
				deallocate(bp)
				draw(1,i) = theta
				mask(i,1) = .FALSE.
				psum = sum(p,mask)
				thetasum = thetasum + theta
				if (thetasum == n) exit
			end do
			draw(1,dim) = n - thetasum
		end subroutine multinomial
		
		subroutine relabel (n, arr, order)
			! This procudure implements relabeling solving label-switching problem
		      	! assoicated with Bayesian estimation for latent class model
			! The inputs of this procedure are:
			! 1. n -- size of the vector being sorted;
			! 2. arr -- rank 1 vector being sorted;
			! The output of the procedure is rank 2 (n by 1) integer mat -- order
			
			! Jia Yan, Nov 04, 2003   
			use matforF90
			implicit none

			! Declare inputs and outputs
			integer,  intent(in) :: n 
			real(dp), intent(inout),  dimension(n) :: arr
			integer,  intent(out), dimension(n,1) :: order

			! Declare local variables
			integer :: i, j, iptr
			real(dp) :: temp
						
			! Body begins from here
			order = 1
			do i = 1, n-1
				iptr = i
				do j = i+1, n
					if (arr(j) < arr(iptr)) then
						iptr = j
					end if
				end do
				order(i,1) = iptr
	
				if (i /= iptr) then
					temp = arr(i)
					arr(i) = arr(iptr)
					arr(iptr) = temp
				end if
			end do
		end subroutine relabel
end module tools

