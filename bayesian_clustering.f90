
module mysub
	implicit none
	contains
	subroutine bayesian_clustering (nitems, nf, ng, nloop, burnin, dataset, mmat, vmat, pmat, prob, cons)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! This procedure implements Bayesian clustering based on 
		! the assumption that the data is mixture normal. The inputs
		! are: nf -- dimension of the data; ng -- number of groups
		! of clustering; nloop -- number of Gibbs draws; 
		! burnin -- number of draws considered as burn in period;
		! dataset -- name of data set (ASCII file);
		! cons -- optional input inputed as integer number 1 or 2. If inputed as 1,
		! the procedure constrians the covariances of different groups are
		! the same; if inputed as 2, the procedure constrain that all covariance
		! matrices are diagonal. 
		! The three outputs are: mmat -- a nf by ng matrix, and each column
		! represents the mean of each group; vmat -- a nf*ng by ng matrix, and
		! each ng by ng submatrix represents the covariance matrix of each group;
		! member -- membership of each item in the data set, and it is a
		! # of items by nf matrix
	
		! Jia Yan, Sep. 12, 2003    	 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		use tools
		use matforF90
		use numerical_libraries
		implicit none
	
		! Declare inputs and outputs
		integer, intent(in) :: nitems, nf, ng, nloop, burnin
		integer, intent(in), optional :: cons
		character(40), intent(in) :: dataset
		real(dp), intent(out), dimension(nf, ng) :: mmat
		real(dp), intent(out), dimension(nf*ng, nf) :: vmat
		real(dp), intent(out), dimension(ng,1) :: pmat
		real(dp), intent(out), dimension(nitems,ng) :: prob
			
		! Declare local variables
		integer, parameter :: ifdata = 15
		real(dp), dimension(ng,1) :: pmat_ini, pmat_final, amat, mvector
		real(dp), dimension(nf,1) :: x, mutemp, xsum
		real(dp), dimension(1,nf) :: ftemp
		integer, dimension(1,ng) :: ctemp
		real(dp), dimension(1,ng) :: probtemp
		real(dp), dimension(nf,ng) :: munod, mmat_ini, mmat_final
		real(dp), dimension(nf*ng,nf) :: vmat_ini, vmat_final
		real(dp), dimension(nf,nf) :: gnod, smat, vtemp, wmat, xtemp
		real(dp), dimension(nitems,nf) :: factors
		real(dp), dimension(nitems,ng) :: prob_inter
		integer, dimension(nitems,ng) :: member_ini
		logical, dimension(nitems,ng) :: mask
		real(dp) :: fv, cd, ld, chirho, chiscale, chideg, chi2
		integer :: iloop, igroup, ifactor, j, indic1, indic2, rho, istat, degree
		character(40) :: bidon
			
		! Body begins from here
		open(ifdata,FILE=dataset,STATUS='OLD',FORM='FORMATTED',ACTION='READ',IOSTAT=istat)
	
		if (istat /= 0) then
			write (*,*) 'File Open Failed'
		else	
			! Read in data
			do j = 1, nitems
				read(ifdata,*) bidon, factors(j,:)		
			end do

			close(ifdata, STATUS = 'KEEP')

			! Set priors for parameters
			amat = 2
			munod = 0.
			gnod = eye(nf)

			if (present(cons)) then
				if (cons == 1) then
					rho = nf + 10
					smat = eye(nf)
				else 
					chirho = 20.
					smat = 4*eye(nf)
				end if
			else
				rho = nf + 10
				smat = eye(nf)	
			end if

			! Set initial values for parameters
			member_ini = 0
			pmat_ini = 1./ng
			mmat_ini = 1.
			vmat_ini = 0.
			do igroup = 1,ng
				indic1 = (igroup - 1) * nf + 1
				indic2 = (igroup - 1) * nf + nf
				vmat_ini(indic1:indic2,:) = eye(nf) 
			end do
	 
			mmat_final = 0.
			vmat_final = 0.
			pmat_final = 0.
			prob = 0.
			do iloop = 1, nloop
				print *, iloop

				! Updating parameters here
				mask = .FALSE.
				pmat = pmat_ini
				mmat = mmat_ini
				vmat = vmat_ini
				
				! Relabeling based on the first element of each column of mmat
				! with ascending order
				call matsort(nf,ng,1,mmat)
				
				! Draw member
				do j = 1, nitems
					ftemp(1,:) = factors(j,:)
					x = transpose(ftemp)
					mvector = 0.
					do igroup = 1, ng
						indic1 = (igroup - 1) * nf + 1
						indic2 = (igroup - 1) * nf + nf
						mutemp(:,1) = mmat(:,igroup)
						vtemp(:,:) = vmat(indic1:indic2,:)
						call mnormal_density(nf, x, mutemp, vtemp, fv)
						mvector(igroup,1) = fv * pmat(igroup,1)
					end do
					fv = sum(mvector)
					mvector = mvector / fv
					call multinomial(ng, 1, mvector, ctemp)
					member_ini(j,:) = ctemp(1,:)
					do igroup = 1,ng
						if (member_ini(j,igroup) == 1) then
							mask(j,igroup) = .TRUE.
						end if
					end do
					
					probtemp = transpose(mvector)
					prob_inter(j,:) = probtemp(1,:)
				end do

				! Draw pmat
				mvector = 0.
				do igroup = 1, ng
					mvector(igroup,1) = amat(igroup,1) + count(mask(:,igroup))
				end do
				call dirichlet(ng, mvector, pmat_ini)
				
				! Draw mmat
				do igroup = 1, ng
					do j = 1, nf
						fv = sum(factors(:,j), mask(:,igroup))
						cd = 1 / ( (count(mask(:,igroup)) / vmat(j,j)) + &
							(1/gnod(j,j)) )
					
						ld = (fv / vmat(j, j)) + &
								(munod(j,igroup) / gnod(j, j))

						call drnnor(1, fv)	
						mmat_ini(j,igroup) = cd*ld + sqrt(cd)*fv
					end do
				end do
				
				! Draw vmat
				if (present(cons)) then
					if (cons == 1) then
						! This part implements the constraint that
						! all groups have the same covariance matrix.
						vtemp = 0.
						do igroup = 1,ng
							do j = 1,nitems
								if (member_ini(j,igroup) == 1) then
									ftemp(1,:) = factors(j,:)
									x = transpose(ftemp)
									x(:,1) = x(:,1) - mmat_ini(:,igroup)
									xtemp = x .xt. x
									vtemp = vtemp + xtemp
								end if
							end do
						end do
						wmat = vtemp + rho * smat
						wmat = invpd(wmat)
						degree = rho + nitems
						call wishart_jia(nf, wmat, degree, vtemp)
						vtemp = invpd(vtemp)
						do igroup = 1,ng
							indic1 = (igroup - 1) * nf + 1
							indic2 = (igroup - 1) * nf + nf
							vmat_ini(indic1:indic2,:) = vtemp(:,:)
						end do
					else
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! This part (cons == 2) implements the constraint that 
						! covariance matrices are all diagonal. ( To be added) 
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						do igroup = 1,ng
							indic1 = (igroup - 1) * nf + 1
							indic2 = (igroup - 1) * nf + nf
							vtemp = 0.
							xsum = 0.	
							do j = 1,nitems
								if (member_ini(j,igroup) == 1) then
									ftemp(1,:) = factors(j,:)
									x = transpose(ftemp)
									x(:,1) = x(:,1) - mmat_ini(:,igroup)
									do ifactor = 1, nf
										x(ifactor, 1) = x(ifactor,1) ** 2
									end do 
									xsum = xsum + x
								end if
							end do

							do ifactor = 1, nf
								chiscale = chirho * smat(ifactor,ifactor) + xsum(ifactor,1)
								chideg = chirho + count(mask(:,igroup)) 
								call drnchi(1, chideg, chi2)
								vtemp(ifactor,ifactor) = chiscale / chi2
							end do
							vmat_ini(indic1:indic2,:) = vtemp(:,:)
						end do
					end if
				else
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! This part estimates covariance matrcies 
					! of each group without any constraints.
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					do igroup = 1,ng
						indic1 = (igroup - 1) * nf + 1
						indic2 = (igroup - 1) * nf + nf
						vtemp = 0.
						do j = 1,nitems
							if (member_ini(j,igroup) == 1) then
								ftemp(1,:) = factors(j,:)
								x = transpose(ftemp)
								x(:,1) = x(:,1) - mmat_ini(:,igroup)
								xtemp = x .xt. x
						        vtemp = vtemp + xtemp
							end if
						end do
						wmat = vtemp + (rho * smat)
						wmat = invpd(wmat)
						degree = rho + count(mask(:,igroup))
						call wishart_jia(nf, wmat, degree, vtemp)
						vtemp = invpd(vtemp)
						vmat_ini(indic1:indic2,:) = vtemp(:,:)
					end do
				end if

				! Storing Results	
				if (iloop > burnin) then
					mmat_final = mmat_final + mmat_ini
					vmat_final = vmat_final + vmat_ini
					pmat_final = pmat_final + pmat_ini
					prob = prob + prob_inter
				end if
					
			end do

			mmat = mmat_final / (nloop - burnin)
			vmat = vmat_final / (nloop - burnin)
			pmat = pmat_final / (nloop - burnin)
			prob = prob / (nloop - burnin)

		end if
	end subroutine bayesian_clustering
end module mysub

program testuse
	use mysub
	use tools
	use matforF90
	implicit none

	! Decalre local variables
	integer, parameter :: iuin = 16
	integer, parameter :: iuout = 17
	integer, parameter :: nobs = 1331
	integer, parameter :: dim = 5
	integer, parameter :: group = 10
	integer, parameter :: ndraws = 6000
	integer, parameter :: throw = 1000
	character(40) :: infile = 'factor_bio_5.asc'
	character(40) :: outfile = 'output10.asc'
	
	real(dp), dimension(dim,group) :: mean
	real(dp), dimension(dim*group,dim) :: var
	real(dp), dimension(group,1) :: p
	integer, dimension(nobs,group) :: cluster
	real(dp), dimension(nobs,group) :: probability
	real(dp), dimension(nobs,dim) :: mydata
	character(40), dimension(nobs,1) :: id
	real(dp), dimension(dim,1) :: x, mutemp
	real(dp), dimension(1,dim) :: ftemp
	real(dp), dimension(dim,dim) :: vtemp
	real(dp) :: fv, loglikelihood, lli, BIC
	integer :: j,i,igroup,istat,indic1,indic2

	! Body begins from here
	open(iuin,  FILE=infile,  STATUS='OLD',     FORM='FORMATTED',ACTION='READ',IOSTAT=istat)	
	open(iuout, FILE=outfile, STATUS='UNKNOWN', FORM='FORMATTED',ACTION='WRITE', RECL = 8000)
	
	! Read in item id for the purpose of output results
	do j =1, nobs
		read(iuin,*) id(j,1), mydata(j,:)
	end do 
	! Close input data file
	close(iuin, STATUS = 'KEEP')

	! Do Bayesian model-based clustering
	call bayesian_clustering(nobs,dim,group,ndraws,throw,infile, mean,var,p,probability,1)

	! Calculate BIC 
	! Draw member
	loglikelihood = 0.
	do j = 1, nobs
		ftemp(1,:) = mydata(j,:)
		x = transpose(ftemp)
		lli = 0.
		do igroup = 1, group
			indic1 = (igroup - 1) * dim + 1
			indic2 = (igroup - 1) * dim + dim
			mutemp(:,1) = mean(:,igroup)
			vtemp(:,:) = var(indic1:indic2,:)
			call mnormal_density(dim, x, mutemp, vtemp, fv)
			lli = lli + fv * p(igroup,1)
		end do
		loglikelihood = loglikelihood + log(lli)
	end do

	print *, loglikelihood
	do j = 1, nobs
		write(iuout,*) id(j,1), achar(9), probability(j,1), achar(9), ((probability(j,i)), achar(9), i=2,group)
	end do 
	call prnt(var)
	 
end program
