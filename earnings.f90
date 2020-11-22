
module mysub
	implicit none
	contains

	subroutine mixture_reg (nitems, nf, ng, cons, nloop, burnin, dataset, mmat, vmat, pmat, bic)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! This procedure estimates regression model with random parameters
		! in which the distribution of parametesr has a discrete type. 
		! The inputs are:
		! ntimes -- number of obs;
		! nf -- dimension of the parameter space;
		! ng -- number of groups;
		! cons -- optional input inputed as integer number 1 (constarined variance) or 0 (unconstarined). 
		! nloop -- number of Gibbs draws; 
		! burnin -- number of draws considered as burn in period;
		! dataset -- name of data set (ASCII file);

		! The four outputs are: 
		! mmat -- a nf by ng matrix, and each column represents the mean of parameters of each group;
		! vmat -- a ng by 1 vector recording the variance of error of each group;
		! pmat -- probability of each item belonging to each group
		! bic -- Bayesian Information Criterion for model selection 
			
		! Jia Yan, Hong Kong Polytechnic,  June 03, 2006    	 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		use tools
		use matforF90
		use numerical_libraries
		implicit none
	
		! Declare inputs and outputs
		integer, intent(in)							:: nitems, nf, ng, nloop, burnin, cons
		character(40), intent(in)					:: dataset
		real(dp), intent(out), dimension(nf, ng)    :: mmat
		real(dp), intent(out), dimension(ng, 1)     :: vmat
		real(dp), intent(out), dimension(ng,1)      :: pmat
		real(dp), intent(out)						:: bic
			
		! Declare local variables
		integer,  parameter					:: ifdata = 15
		integer,  parameter					:: ifout  = 16	
		integer,  parameter					:: pi = 3.1415926535897932
		integer,  dimension(nitems, ng)		:: member
		integer,  dimension(ng, 1)			:: label, counts	
		integer,  dimension(1, ng)			:: ctemp

		character(40)						:: outfile = "out.asc"
		real(dp), dimension(nitems, 1)		:: y, yg, likelihood
		real(dp), dimension(nitems, nf)		:: x, xg 
		real(dp), dimension(1, nf)			:: xtemp
		real(dp), dimension(nf, 1)			:: mbprior, btemp, xy, randn
		real(dp), dimension(nf, nf)			:: vbprior, xx, vcb
		real(dp), dimension(nf, ng)			:: mmat_ini
		real(dp), dimension(ng)				:: array
		real(dp), dimension(ng, 1)			:: pmat_ini, vmat_ini, amat, mvector, epsi2
		real(dp), dimension(1, 1)			:: error

		logical, dimension(nitems,ng)		:: mask

		integer  :: iloop, i, j, igroup, istat, num, bidon
		real(dp) ::	vtemp, scalar, fv, e, degree, rho, tao

		! Body begins from here
		open(ifdata,FILE=dataset,STATUS='OLD',FORM='FORMATTED',ACTION='READ',IOSTAT=istat)
		open(ifout, FILE=outfile, STATUS='UNKNOWN', FORM='FORMATTED',ACTION='WRITE', RECL = 80000)
	
		DataCheck: if (istat /= 0) then
			write (*,*) 'File Open Failed'
		else	
			! Read in data
			do j = 1, nitems
				read(ifdata,*) bidon, y(j,1), x(j,:)		
			end do

			close(ifdata, STATUS = 'KEEP')

			! Set priors for parameters
			mbprior = 0.
			vbprior = 0.01 * eye(nf)
			rho = 4+2
			tao = 4.	
			amat = 0.2

			! Set initial values for parameters
			member = 0
			array = 0.
			label = 1
			pmat_ini = 1./ng
			mmat_ini = 1.
			vmat_ini = 1.
			mvector = 0.
		 
			mmat = 0.
			vmat = 0.
			pmat = 0.

			!! Calculate moment matrix for ng = 1
			xx = x .tx. x
			xy = x .tx. y 
			Gibbs: do iloop = 1, nloop
				write(*,*)"LOOP:", iloop

				if (ng > 1) then

					! Updating parameters here
					mask = .FALSE.

					! Relabeling based on the first element of each column of mmat
					! with ascending order
								
					if (iloop == 1) then
						do igroup = 1, ng
							label(igroup,1) = igroup
						end do
					else
						array(:) = mmat_ini(1,:)
						call relabel(ng,array,label)
					end if
						
					! Draw member
					do j = 1, nitems
						xtemp(1,:) = x(j,:)
						do igroup = 1, ng
							num = label(igroup,1)
							btemp(:,1) = mmat_ini(:, num)
							vtemp = vmat_ini(num, 1)
							error = xtemp .x. btemp
							e = y(j,1) - error(1,1) 
							e = e ** 2
							scalar = -0.5 * (e/vtemp)
							fv = 1 / sqrt(2 * pi * vtemp)
							fv = fv * exp(scalar)    						 
							mvector(igroup,1) = fv * pmat_ini(num,1)
						end do
						fv = sum(mvector)
						mvector = mvector / fv
						call multinomial(ng, 1, mvector, ctemp)
						member(j,:) = ctemp(1,:)
						do igroup = 1,ng
							if (ctemp(1,igroup) == 1) then
								mask(j,igroup) = .TRUE.
							end if
						end do
					end do

					! Draw pmat and mmat
					counts = 0
					do igroup = 1, ng
						! for pmat
						num = count(mask(:, igroup))
						counts(igroup,1) = num
						mvector(igroup,1) = amat(igroup,1) + real(num)

						! for mmat
						num = label(igroup,1)
						scalar = vmat_ini(num,1)
						call drnnor(nf, randn)

						yg(:,1) = y(:,1) * member(:,igroup) 
						do j = 1, nf
							xg(:,j) = x(:,j) * member(:,igroup) 
						end do

						xx = xg .tx. xg
						xy = xg .tx. yg
	
						vcb = invpd( (xx/scalar) + vbprior)
						btemp = vcb .x. ( (xy/scalar) + (vbprior .x. mbprior) )
						btemp = btemp + ( chol(vcb) .tx. randn )

						mmat_ini(:,igroup) = btemp(:,1)

						yg = yg - (xg .x. btemp)	
						scalar = dot_product( yg(:,1), yg(:,1) )
						epsi2(igroup,1) = scalar

					end do
					call dirichlet(ng, mvector, pmat_ini)

					! Draw vmat
					if (cons == 1) then
						scalar = sum(epsi2) + rho * tao
						degree = rho + nitems
						call drnchi(1, degree, fv)
						vmat_ini = scalar / fv
					else
						do igroup = 1, ng
							scalar = epsi2(igroup,1) + rho * tao
							degree = rho + real(count(mask(:,igroup)))
							call drnchi(1, degree, fv)
							vmat_ini(igroup,1) = scalar / fv
						end do
					end if


				else
					pmat_ini = 1.

					call drnnor(nf, randn)
					scalar = vmat_ini(1,1)
					vcb = invpd( (xx/scalar) + vbprior)
					mmat_ini = vcb .x. ( (xy/scalar) + (vbprior .x. mbprior) )
					mmat_ini = mmat_ini + ( chol(vcb) .tx. randn )

					yg = y - (x .x. mmat_ini)	
					scalar = dot_product( yg(:,1), yg(:,1) )
					scalar = scalar + rho * tao
					degree = rho + nitems
					call drnchi(1, degree, fv)
					vmat_ini(1,1) = scalar / fv

				end if

				! Storing Results	
				if (iloop > burnin) then
					mmat = mmat + mmat_ini
					vmat = vmat + vmat_ini
					pmat = pmat + pmat_ini
				
					! write results to the hard disk for future use
					write(ifout,*) iloop, achar(9), ( ( (mmat_ini(i, j), achar(9)), i = 1, nf), j = 1, ng ), &
					( (vmat_ini(i, 1), achar(9)), i = 1, ng), ( (pmat_ini(i, 1), achar(9)), i = 1, ng-1), &
					pmat_ini(ng,1)

				end if
					
			end do Gibbs

			mmat = mmat / (nloop - burnin)
			vmat = vmat / (nloop - burnin)
			pmat = pmat / (nloop - burnin)

			!!! Need one more part for calculating BIC !!!
			likelihood = 0.
			if (ng > 1) then
				if (cons == 1) then
					num = ng * nf + ng
				else
					num = ng * nf + 1
				end if

				do igroup = 1, ng
					btemp(:,1) = mmat(:, igroup)
					vtemp = vmat(igroup, 1)
					yg = y - (x .x. btemp)
					yg = yg ** 2
					yg = yg / vtemp
					yg = -0.5 * yg 
					fv = 1 / sqrt(2 * pi * vtemp)
					yg = fv * exp(yg)
					yg = yg * pmat(igroup,1)
					likelihood = likelihood + yg  
				end do
			else
				num = nf + 1
				vtemp = vmat(1,1)
				yg = y - (x .x. mmat)
				yg = yg ** 2
				yg = yg / vtemp
				yg = -0.5 * yg 
				fv = 1 / sqrt(2 * pi * vtemp)
				yg = fv * exp(yg)
				likelihood = yg 
			end if
			likelihood = log(likelihood)
			fv = 2 * sum(likelihood)
			bic = fv - num * log(real(nitems))

			write(ifout,*) bic

		end if DataCheck

	end subroutine mixture_reg

end module mysub

program main
	use mysub
	use tools
	use matforF90
	implicit none

	! Decalre local variables
	integer:: nobs, dim, group, const, ndraws, throw 
	character(40) :: infile, buffer
	
	real(dp), allocatable, dimension(:,:) :: mean, var, p
	real(dp) :: BIC

	! read in arguments
	call getarg(1, buffer)
	read (buffer,*) nobs
	call getarg(2, buffer)
	read (buffer,*) dim
	call getarg(3, buffer)
	read (buffer,*) group
	call getarg(4, buffer)
	read (buffer,*) const
	call getarg(5, buffer)
	read (buffer,*) ndraws
	call getarg(6, buffer)
	read (buffer,*) throw
	call getarg(7, buffer)
	infile = buffer

	allocate(mean(dim,group), var(group,1), p(group,1))	

	! Do Bayesian model-based clustering
	call mixture_reg(nobs, dim, group, const, ndraws, throw,infile, mean,var,p,BIC)

	write(*,*) "Probability Vector:"
	call prnt(p)
	write(*,*) "Beta:" 
	call prnt(mean)
	write(*,*) "Variance:"
	call prnt(var)
	write(*,*) "BIC:"
	print *, BIC	 

end program main
