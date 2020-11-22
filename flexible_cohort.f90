
module mysub
	implicit none
	contains
	subroutine bayesian_clustering (ntotal, nusers, nitems, capK, capL, ng, nloop, burnin, dataset, &
									bmat, vmat, amat, pmat, prob, cons)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! This procedure implements Bayesian clustering based on 
		! the assumption that the data is mixture normal. The inputs
		! are: nf -- dimension of the data; ng -- number of groups
		! of clustering; nloop -- number of Gibbs draws; 
		! burnin -- number of draws considered as burn in period;
		! dataset -- name of data set (ASCII file);
		! cons -- optional input inputed as integer number 1 or 2. If inputed as 1,
		! the procedure constrains the covariances of different groups are
		! the same; if inputed as 2, the procedure constrain that all covariance
		! matrices are diagonal. 
		! The four outputs are: mmat -- a nf by ng matrix, and each column
		! represents the mean of each group; vmat -- a nf*ng by ng matrix, and
		! each ng by ng submatrix represents the covariance matrix of each group;
		! pmat -- ng by 1 matrix recording the share of items in each group.
		! prob -- nitems by ng matrix, recording the probability of each item in
		! each group.
	
		! Jia Yan, Sep. 29, 2003    	 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		use tools
		use matforF90
		use numerical_libraries
		implicit none
	
		! Declare inputs and outputs
		integer,	   intent(in)							:: ntotal, nusers, nitems, capK, &
															   capL, ng, nloop, burnin
		integer,	   intent(in), optional					:: cons
		character(40), intent(in)							:: dataset
		real(dp),	   intent(out), dimension(capK, ng)		:: bmat
		real(dp),	   intent(out), dimension(ng, 1)		:: vmat
		real(dp),      intent(out), dimension(nitems, ng)	:: amat
		real(dp),      intent(out), dimension(ng, 1)		:: pmat
		real(dp),      intent(out), dimension(nusers,ng)	:: prob
			
		! Declare local variables
		integer,	   parameter						:: ifdata = 15
		character(40), dimension(nusers,1)				:: userid
		integer,	   dimension(nitems,1)				:: itemid
		real(dp),      dimension(ntotal,capK+1)			:: mydata
		real(dp),	   dimension(nusers,capL)			:: z
		integer,       dimension(nusers,1)				:: times_user
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!Revised up to here 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
		real(dp), dimension(ng,1) :: pmat_ini, pmat_final, dmat, mvector
		real(dp), dimension(capK,ng) :: bmat_ini, bmat_final
		real(dp), dimension(capK,1) :: bnod, btemp, randn
		real(dp), dimension(capK,capK) :: invbnod, xx
		real(dp), dimension(nitems,ng) :: alphha, alphha_ini, alphha_final 
		real(dp), dimension(1,nf) :: ftemp
		integer, dimension(1,ng) :: ctemp
		real(dp), dimension(1,ng) :: probtemp
		real(dp), dimension(nf,ng) :: munod, mmat_ini, mmat_final
		real(dp), dimension(nf*ng,nf) :: vmat_ini, vmat_final
		real(dp), dimension(nf,nf) :: gnod, smat, vtemp, wmat, xtemp
		real(dp), dimension(nitems,nf) :: factors
		real(dp), dimension(nusers,ng) :: prob_inter
		integer, dimension(nusers,ng) :: member_ini
		logical, dimension(nusers,ng) :: mask
		logical, dimension(ntotal)    :: ind
		logical						  :: ltemp
		real(dp) :: fv,cd,ld,rhoepsi,rhoalphha,sepsi,salphha,chiscale,chideg,chi2, stemp
		integer :: num, iloop, igroup, iuser, ifactor, j, l, indic1, indic2, rho, istat, degree, kstart, &
				   kend, nobsuser	
		character(40) :: bidon

		real(dp), allocatable, dimension(:,:) :: yuser, xuser, itemuser, yalp, xy, group
			
		! Body begins from here
		open(ifdata,FILE=dataset,STATUS='OLD',FORM='FORMATTED',ACTION='READ',IOSTAT=istat)
	
		if (istat /= 0) then
			write (*,*) 'File Open Failed'
		else	
			! Read in data
			num = 0
			do iuser = 1, nusers
				read(ifdata,*) userid(iuser,1), itemid(num+1,1), mydata(num+1,:), (z(iuser,l), l = 1, capL), &
							   times_user(iuser,1)

				do j = 1, (times_user(iuser,1) - 1)
					read(ifdata,*) bidon, itemid(num+j+1,1), mydata(num+j+1,:)
				end do
				num = num + times_user(iuser,1)		
			end do

			close(ifdata, STATUS = 'KEEP')
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Set priors for parameters
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Drichlet prior 
			dmat = 0.2
			! Prior  for beta
			bnod = 0.
			invbnod = 0.01 * eye(capK)
			! Priro for sigma of epsilon
			rhoepsi = 1.
			sepsi = 1.
			! Prior for sigma of alphha
			rhoalphha = 1.
			salphha = 1.

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Initialize parameters
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			member_ini = 0
			pmat_ini = 1./ng
			bmat_ini = 1.
			vmat_ini = 1.
			alphha_ini = 0.
			 
			bmat_final = 0.
			vmat_final = 0.
			pmat_final = 0.
			alphha_final = 0.
			prob = 0.
			do iloop = 1, nloop
				print *, iloop

				! Updating parameters here
				mask = .FALSE.
				pmat = pmat_ini
				bmat = bmat_ini
				vmat = vmat_ini
				alphha = alphha_ini
								
				! Relabeling based on the first element of each column of mmat
				! with ascending order
				call matsort(capK,ng,1,bmat)
				
				! Draw member
				num = 0
				do iuser = 1, nusers
					nobsuser = times_user(iuser,1)
					kstart = num + 1
					kend = num + nobsuser
					allocate(yuser(nobsuser,1))
					allocate(xuser(nobsuser,capK))
					allocate(itemuser(nobsuser,1))
					yuser(:,1) = mydata(kstart:kend, 1)
					xuser(:,:) = mydata(kstart:kend, 2:(capK+1) )
					itemuser(:,1) = itemid(kstart:kend,1) 
					mvector = 0.
					do igroup = 1, ng
						btemp(:,1) = bmat(:,igroup)
						yuser = yuser - (xuser .x. btemp)
						fv = 1.
						do j = 1, nobsuser
							stemp = yuser(j,1) - alphha(itemuser(j,1),igroup)
							stemp = stemp / sqrt(vmat(igroup,1))
							fv = fv * dnordf(stemp)
						end do
						mvector(igroup,1) = fv * pmat(igroup,1)
					end do
					fv = sum(mvector)
					mvector = mvector / fv
					call multinomial(ng, 1, mvector, ctemp)
					member_ini(iuser,:) = ctemp(1,:)
					do igroup = 1,ng
						if (member_ini(iuser,igroup) == 1) then
							mask(iuser,igroup) = .TRUE.
						end if
					end do
					
					probtemp = transpose(mvector)
					prob_inter(iuser,:) = probtemp(1,:)
					deallocate(yuser)
					deallocate(xuser)
					deallocate(itemuser)
					num = num + times_user(iuser,1)
				end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This part will be raplaced by a MNL model
				! Draw pmat
				mvector = 0.
				do igroup = 1, ng
					mvector(igroup,1) = dmat(igroup,1) + count(mask(:,igroup))
				end do
				call dirichlet(ng, mvector, pmat_ini)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! REVISED UP TO HERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
				do igroup = 1, ng
					num = 0
					do iuser = 1, nusers
						nobsuser = times_user(iuser,1)
						kstart = num + 1
						kend = num + nobsuser	
						ltemp = mask(iuser,igroup)
						ind(kstart:kend) = ltemp
						num = num + nobsuser
					end do
					
					! Draw bmat
					num = count(ind)
					allocate(yuser(num,1))
					allocate(yalp(num,1))
					allocate(xuser(num,capK))
					allocate(itemuser(num,1))
					allocate(xy(num,1))
					yuser(:,1) = selif(mydata(:,1), ind)
					xuser(:,:) = selif(mydata(:,2:(capK+1)), ind)
					itemuser(:,1) = selif(itemid(:,1), ind)
					
					yalp = 0.
					do j = 1, num
						yalp(j,1) = yuser(j,1) - alphha(itemuser(j,1), igroup)
					end do
					
					call drnnor(capK, randn)
					xx = xuser .tx. xuser
					xx = xx / vmat(igroup,1)
					xx = xx + invbnod
					xx = invpd(xx)

					xy = xuser .tx. yalp 
					xy = xy / vmat(igroup,1)
					xy = xy + (invbnod .x. bnod)
					
					btemp = xx .x. xy
					btemp = btemp + (chol(xx) .tx. randn) 
					bmat_ini(:,igroup) = btemp(:,1)
					
					! Draw alphha
					num = count(mask(:,igroup))
					allocate(group(num,1))
					yalp = yuser - (xuser .x.  btemp)
					! STOP HERE, Sep. 29
										 	
					deallocate(yuser)
					deallocate(yalp)
					deallocate(xuser)
					deallocate(itemuser) 
					deallocate(xy)
					deallocate(group)
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
	integer, parameter :: nobs = 2000
	integer, parameter :: dim = 2
	integer, parameter :: group = 2
	integer, parameter :: ndraws = 2000
	integer, parameter :: throw = 500
	character(40) :: infile = 'simulatedata.asc'
	character(40) :: outfile = 'output.asc'
	
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
	call bayesian_clustering(nobs,dim,group,ndraws,throw,infile, mean,var,p,probability,2)

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
