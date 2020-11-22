
module mysub1
	implicit none
	contains

	subroutine getmat(dim, t, x, matout)
		! This procedure is to generate the covariance matrix
		! of the errors which are AR(1) correlated
		
		! Jia, Jan. 13, 2007		
		use matforF90
		implicit none

		! define inputs and outputs
		integer, intent(in) :: dim
		integer, intent(in), dimension(dim,1) :: t
		real(dp), intent(in) :: x
		real(dp), intent(out), dimension(dim, dim) :: matout

		! define locals
		integer :: j, i, k, size
		real(dp), allocatable, dimension(:) :: junk  
		matout = 0.

	   j = 0	
	   Outer: do i = 1, dim
	        size = dim - j
				
			allocate(junk(size))
			Inner: do k = 1, size
				junk(k) = x ** (t(k+j,1) - t(i,1))
			end do Inner

			matout( (j+1):dim,i) = junk(:)
			matout( i,(j+1):dim) = junk(:)

			deallocate(junk)
			j = j + 1
	   end do Outer
	
	end subroutine getmat

	subroutine timelag(dim, yin, yout)
		use matforF90
		implicit none

		! declare inputs and outputs
		integer, intent(in) :: dim
		real(dp), intent(in), dimension(dim,1) :: yin
		real(dp), intent(out), dimension(dim-1,2) :: yout

		! declare locals
		integer :: i
		
		! code starts from here
		do i = 1, dim - 1
			yout(i,1) = yin(i+1,1)
			yout(i,2) = yin(i,1)
		end do	
	
	end subroutine timelag	 

end module mysub1

module mysub2
	
	implicit none
	contains

	subroutine pwreg (ntotal, nfirm, nmax, capK, nloop, dataset)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! This procedure is to implement the Bayesian Prais-Winsten type regression on the unbalanced 
		! and interupted panel data. 
		 	
		! Jia Yan, Hong Kong Polytechnic,  Jan. 12, 2007    	 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		use mysub1
		use tools
		use matforF90
		use numerical_libraries
		implicit none
	
		! Declare inputs and outputs
		integer, intent(in)					:: ntotal, nfirm, nmax, capK, nloop
		character(40), intent(in)			:: dataset
			
		! Declare local variables
		integer,  parameter					:: ifdata = 15
		integer,  parameter					:: ifout  = 16	

		character(40)						:: outfile = "out.asc"
		character(40), dimension(nfirm)		:: id				
		integer,  dimension(nfirm)			:: times, times_max
		integer,  dimension(ntotal,1)		:: tseries, gap
		real(dp), dimension(ntotal, 1)		:: yvec, error
		real(dp), dimension(nmax-nfirm,2)  :: elag
		real(dp), dimension(ntotal, capK)	:: xmat 
		real(dp), dimension(capK, 1)		:: b0, b_ini, b, xy, randn
		real(dp), dimension(capK, capK)		:: vb, xx, vcb

		integer,  allocatable, dimension(:,:) :: tjunk
		real(dp), allocatable, dimension(:,:) :: yjunk, xjunk, covmat

		integer  :: istat, ifirm, iloop, i, j, k, m, bidon, num, kstart, kend, counter, nlag, maxn
		real(dp) ::	degree, rho_mu, tao_mu, corr0, vcorr, corr_ini, corr, sigmu_ini, sigmu, ts1, ts2, rscalar, &
					ms, vs

		! Body begins from here
		open(ifdata,FILE=dataset,STATUS='OLD',FORM='FORMATTED',ACTION='READ',IOSTAT=istat)
		open(ifout, FILE=outfile, STATUS='UNKNOWN', FORM='FORMATTED',ACTION='WRITE', RECL = 80000)
	
		DataCheck: if (istat /= 0) then
			write (*,*) 'File Open Failed'
		else	
			! Read data into memory once for all
			num = 0
			dataLoad: do ifirm=1,nfirm
				read(ifdata,*) id(ifirm), yvec(num+1, 1), xmat(num+1,:), tseries(num+1,1), gap(num+1,1), &
							   times(ifirm), times_max(ifirm)
				do j=1,(times(ifirm)-1)
					read(ifdata,*) bidon, yvec(num+j+1,1), xmat(num+j+1,:), tseries(num+j+1,1), gap(num+j+1,1)
				end do
				num = num + times(ifirm)
			end do dataLoad

			close(ifdata, STATUS = 'KEEP')

			! Set priors for parameters
		    ! priors for betta
			b0 = 0
			vb = .01 * eye(capK)

			! priors on corr
			corr0 = 0.5
			vcorr = 0.01

			! priors for sig_mu		
			rho_mu  = 1
			tao_mu  = 1

			! Set initial values for parameters
			b_ini = 0.
			corr_ini = 0.5
			sigmu_ini = 1.

			Gibbs: do iloop = 1, nloop
				write(*,*)"LOOP:", iloop

				! Update
				b = b_ini
				corr = corr_ini
				sigmu = sigmu_ini

				!! matrix accummulation !!
				counter = 0
				xx = 0.
				xy = 0.
				accum: do ifirm = 1, nfirm
					num = times(ifirm)
					kstart = counter + 1
					kend = counter + num

					allocate(yjunk(num,1))
					allocate(xjunk(num, capK))
					allocate(tjunk(num,1))
					allocate(covmat(num, num))

					yjunk(:,1) = yvec(kstart:kend, 1)
					xjunk(:,:) = xmat(kstart:kend, :) 
					tjunk(:,1) = tseries(kstart:kend,1)
					call getmat(num, tjunk, corr, covmat)
					
					covmat = (sigmu / (1 - corr ** 2)) * covmat 
					covmat = invpd(covmat)
					covmat = chol(covmat)
					yjunk = covmat .tx. yjunk
					xjunk = covmat .tx. xjunk

					xx = xx + (xjunk .tx. xjunk)
					xy = xy + (xjunk .tx. yjunk)

					deallocate(yjunk)
					deallocate(xjunk)
					deallocate(tjunk) 
					deallocate(covmat)

					counter = counter + num
				end do accum
				
				! update betta
				call drnnor(capK, randn)
				vcb = invpd(xx + vb)
				b_ini = vcb .x. ( xy + (vb .x. b0) )
				b_ini = b_ini + ( chol(vcb) .tx. randn )

				! update corr
				error = yvec - (xmat .x. b_ini)
				elag = 0.
				counter = 0
				nlag = 0
				ts1 = 0.
				ts2 = 0.
				lag: do ifirm = 1, nfirm
					num = times(ifirm)
					maxn = times_max(ifirm)

					allocate(xjunk(maxn,1))
					allocate(yjunk(maxn-1,2))
					
					!!!! filling in the missing points in an interrupted series
					bidon = 0
					fillin: do j = 1, num
						bidon = bidon + 1

						! to remeber that gap of the first observation is always 1
						if (gap(counter+j,1) == 1 ) then	
							xjunk(bidon,1) = error(counter + j,1)
						else
							do k = 1, (gap(counter+j,1) - 1)
								call drnnor(1, rscalar)
								xjunk(bidon, 1) = corr * xjunk(bidon - 1, 1) + sqrt(sigmu) * rscalar
								bidon = bidon + 1
							end do
							xjunk(bidon,1) = error(counter + j,1)

						end if
					end do fillin 

					call timelag(maxn, xjunk, yjunk)

					ts1 = ts1 + dot_product(yjunk(:,2), yjunk(:,2))
					ts2 = ts2 + dot_product(yjunk(:,2), yjunk(:,1))

					kstart = nlag + 1
					kend = nlag + maxn - 1
					elag(kstart:kend, :) = yjunk(:,:)
					
					deallocate(xjunk)
					deallocate(yjunk)
					counter = counter + num
					nlag = nlag + maxn - 1
				end do lag
				ts1 = ts1 / sigmu
				ts2 = ts2 / sigmu
				call drnnor(1, rscalar)
				corr_ini = (vcorr * corr0 + ts2) / (ts1 + vcorr)
				corr_ini = corr_ini + sqrt(1 / (ts1 + vcorr)) * rscalar

				elag(:,1) = elag(:,1) - corr_ini * elag(:,2)				
				ts1  = dot_product(elag(:,1), elag(:,1))
	
				ts2 = ts1 + rho_mu * tao_mu
				degree = rho_mu + nmax - nfirm
				call drnchi(1, degree, rscalar)
				sigmu_ini = ts2 / rscalar
				
				! write results to the hard disk for future use
				write(ifout,*) iloop, achar(9), ((b_ini(j,1), achar(9)),j=1,capK), corr_ini, achar(9), &
							   sigmu_ini 	 

			
			end do Gibbs


		end if DataCheck

	end subroutine pwreg

end module mysub2

program main
	use mysub2
	use tools
	use matforF90
	implicit none

	! Decalre local variables
	integer:: nobs, nf, mx, dim, ndraws 
	character(40) :: infile, buffer

	!nobs = 21516
	!nf = 938
	!mx = 26861
	!dim = 4
	!ndraws = 100
	!infile = "datain.asc"
	
	! read in arguments
	call getarg(1, buffer)
	read (buffer,*) nobs
	call getarg(2, buffer)
	read (buffer,*) nf
	call getarg(3, buffer)
	read (buffer,*) mx
	call getarg(4, buffer)
	read (buffer,*) dim
	call getarg(5, buffer)
	read (buffer,*) ndraws
	call getarg(6, buffer)
	infile = buffer


	! Do Bayesian model-based clustering
	call pwreg(nobs, nf, mx, dim, ndraws, infile)

end program main
