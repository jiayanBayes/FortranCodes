
! simple version.  Next version will be more flexible since dynamic allocation
! will be used.  This is for testing.
!
!	Writes two tab separated ascii output files. 
!		output.asc contains loop, gamma1,...,gammaL, sigDelta1,...,sigDeltaK, sigmaepsi
!       betan.asc  conains  loop, userid, beta1,...,betak
! Denis Bolduc, December 18, 2002. 

program MCMC

	use numerical_libraries
	use matforF90
	use linear_operator
    implicit none

	integer, parameter :: nusers      = 2000          ! max number of users
	integer, parameter :: capK        = 3             ! number of expl. variables
	integer, parameter :: capL        = 2             ! number of z variables in generation
	                                                  ! of beta_n's.
	integer, parameter :: nbrLoops    = 200           ! nbr loops for Gibbs sampling   

	integer, parameter :: ifinput         = 10            ! unit input file
	integer, parameter :: ifoutput        = 11            ! unit output file
	integer, parameter :: ifoutputbeta    = 12            ! unit output file betan's

	character(80)  :: infile      = "dataset.asc"
	character(80)  :: outfile     = "output.asc"
	character(80)  :: outfilebeta = "betan.asc"
	character(200) :: string

	integer,  allocatable, dimension(:)         :: userid
	! times records the number of observations for each user. This setting can handle unbalanced panel data 
	integer,  dimension(nusers)		 		    :: times 	
	real(dp), allocatable, dimension(:,:)       :: y, epsi
	real(dp), allocatable, dimension(:,:)    :: x
    real(dp), dimension(nusers,capL) :: z
	   
   	real(dp), dimension(nusers, capK)      :: b, bres

	real(dp), dimension(capL,capL) :: vcg, zz
	real(dp), dimension(capL * capK, capL * capK) :: vg  
	real(dp), dimensaion
	real(dp), dimension(capK * capL,1)    :: g0, g, accg, gmean, gnew, g_ini, randl, accZsbs

	real(dp), dimension(capK,capK) :: sigmabres, sigmabres_ini, precbres, precbresHalf, XX, vcbn
    real(dp), dimension(capK,1)    :: zng, bn, Xy, bmean, randn, sigbres, rho, R
	
 
    real(dp) :: k1, k2, sigmaepsi, sigmaepsi_ini, epsi2, scale, chi2, deg, degree, Smat, tempmat
	integer  :: i,j,k,l, iloop, iuser, bidon, nobs, ik, kstart, kend
	
	character(8) :: aa

			! order on dataset is:
			! userid	y	x1	x2	x3	z_n11	z_n12	z_n21	z_n22	z_n31	z_n32
			! always keep z type data at right of x, and times at the end of dataset

! Code starts here

	open( ifinput,      FILE=infile,       STATUS='UNKNOWN')
	open( ifoutput,     FILE=outfile,      STATUS='UNKNOWN', RECL=8000)
	open( ifoutputbeta, FILE=outfilebeta,  STATUS='UNKNOWN', RECL=8000)


! Fill zn matrix once and for all. Only first observation of userspecific data is stored.
! z matrix is stored as KxL for first two dimension.  Third dimension rolls over individuals.  
	allocate(userid(1))
	allocate(y(1,1))
	allocate(x(1,capK))
	do iuser=1,nusers
		read(ifinput,*) userid(1), y(1,1), x(1,:), (z(l,iuser),l=1,capL), times(iuser)	
		do j=1,(times(iuser)-1)
			read(ifinput,*) bidon
		end do
	end do
	deallocate(userid)
	deallocate(y)
	deallocate(x)
				
! initialize priors for hyperparameters: the 1st layer
	g0 = 0
	vg = .001*eye(capL)	
	k1  = 0
	k2  = 1
! initialize priors for hyperparameters: the 2nd layer
	rho = 0
	R = 1

! write namefields on output datasets
						! Note: could have used achar(9) which is tab character
	string = "Loop"  // achar(9) 
	do j=1,capL
		aa = cvtis(j)
		string = string(1:len_trim(string)) // " gamma" // aa(1:len_trim(aa))  // achar(9) 
	end do
	do j=1,capK
		aa = cvtis(j)
		string = string(1:len_trim(string)) // " delta" // aa(1:len_trim(aa)) // aa(1:len_trim(aa)) // achar(9)  
	end do
	string = string(1:len_trim(string)) // " sigmaepsi"
	write(ifoutput,*) string(1:len_trim(string))

! write namefields on betan output datasets
	string = "Loop" // achar(9) // "userid" // achar(9) 
	do j=1,capK
		aa = cvtis(j)
		if (j .eq. capK) then
			string = string(1:len_trim(string)) // "beta" // aa(1:len_trim(aa))
		else
			string = string(1:len_trim(string)) // "beta" // aa(1:len_trim(aa)) // achar(9)  
		endif
	end do
	write(ifoutputbeta,*) string(1:len_trim(string))
	
! set initial values for parameters
! gamma, sig_eps and Sigma_delta
	g_ini = 0
	sigmaepsi_ini = 1
	sigmabres_ini = eye(capK)
	accg = 0

	iseed = 123457     ! reset seed ot keep same sequence of random numbers

	do iloop = 1,nbrLoops
		print *,"Loop : ", iloop
				! rewind dataset
		rewind ifinput

	! update
		g = g_ini
		sigmabres = sigmabres_ini
		precbres = invpd(sigmabres)
							! keep as lower triangular
		precbresHalf = .t. chol(precbres)
		sigmaepsi = sigmaepsi_ini

	! draw beta_n
					! reset accumulators for gamma later on.
		epsi2 = 0
		accZsbs  = 0
		accZsZs  = 0
		do iuser=1,nusers
			allocate(userid(times(iuser)))
			allocate(y(times(iuser),1))
			allocate(x(times(iuser),capK))
			allocate(epsi(times(iuser),1))
			
		
			
			call drnnor(capK, randn)
			! read in data specific to user: iuser
			
			do j=1,times(iuser)
				read(ifinput,*) userid(j), y(j,1), x(j,:)
			end do
			XX = x .tx. x
			Xy = x .tx. y

			vcbn   = invpd(precbres + XX/sigmaepsi)
			bmean  = vcbn .x. ( (precbres .x. zng) + Xy/sigmaepsi)			 
			bn     = bmean + ( chol(vcbn) .tx. randn )
			b(iuser,:)     = .t. bn(:,1)
			
			! accumulate zstar'zstar and zstar*betanstar
			zz   = zz + (z(iuser,:) .tx. z(iuser,:) )
			do ik = 1, capK
				kstart = (ik - 1) * capL + 1
				kend = (ik - 1) * capL + capL 
				zng(ik,1) = z(iuser,:) .x. g(kstart:kend,1)
			end do
		
	
			! write on betan output file
			write(ifoutputbeta,*) iloop, achar(9), iuser, achar(9), ((bn(j,1), achar(9)),j=1,capK-1), bn(capK,1) 

			! accumulate epsilon2
			epsi  = y - (x .x. bn)
			epsi2 = epsi2 + dot_product(epsi(:,1),epsi(:,1))
			
			deallocate(userid)
			deallocate(y)
			deallocate(x)
			deallocate(epsi)

		end do
		
	! draw sigmaepsi
		scale  = k1*k2 + epsi2
		degree = k1 + sum(times)
		call drnchi(1, degree, chi2)
		sigmaepsi_ini = scale / chi2

	! draw gamma
		call drnnor(capL,randl)
		do ik = 1, capK
			kstart = (ik - 1) * capL + 1
			kend = (ik -1) * capL + capL
			vcg = invpd(vg(kstart:kend, kstart:kend) + zz)
			zb = z .tx. b(:,ik)
			gmean = vcg .x. ( (vg(kstart:kend, kstart:kend) .x. g0(kstart:kend,1)) + zb ) 
			gnew(kstart:kend,1)  = gmean + ( chol(vcg) .tx. randl)
			bres(:,ik) = b(:,ik) - z * gnew(kstart:kend,1)
		end do
		accg  = accg + gnew
		print *, "mean gamma";call prnt(.t. (accg/iloop))
		

		! draw new sigma matrix			
				! see Johnson, p. 204
				! 3 key results:
				! wishart(int(deg), capK, iseed) makes draws from standard Wishart: Z ~ W(m,I), where Z = T*T'.
				! If Z ~ W(m,A), GZG'~W(m,GAG').  A=I corresonds to standard Wishart.
				! If Z ~ W(m,A), then inv(Z) ~ IW(m,inv(A)).

		!Smat    = .t. chol( invpd( eye(capK)*rho + (bres .xt. bres) ))
		!deg     = rho + nusers
		!Tmat    = wishart(int(deg), capK, iseed)
		!tempmat = Smat .x. Tmat
		!resMat  = tempmat .xt. tempmat
		
		do ik = 1, capK
			tempmat = dot_product(bres(:,ik), bres(:,ik))
			Smat = rho(ik,1) * R(ik,1) + tempmat
			deg = rho(ik,1) + nusers
			call drnchi(1, deg, chi2)
			sigmabres_ini(ik,ik) = Smat / chi2
		end do
			 
	! revise Sigma and update
		!sigmabres_ini = invpd(resMat)
		g_ini = gnew
		sigbres(:,1) = diag(sigmabres_ini)

	! write on output file
		write(ifoutput,*) iloop, achar(9),((gnew(j,1), achar(9)),j=1,capL), & 
		   ((sigbres(j,1), achar(9)),j=1,capK-1), sigbres(capK,1), achar(9), sigmaepsi_ini  

	end do

end
