module sharedata

	! To define the globals to be used in later procedures

	! Jia Yan, 11/27/2004
	use matforF90
	implicit none
	save

	! parameters of the price random process
	integer,  parameter :: DIMP = 10
	real(dp), parameter :: PUP = 30
	real(dp), parameter :: PLOW = 16
	real(dp), parameter :: LMD0 = 0.06
	real(dp), parameter :: LMD1 = 0.98
	real(dp), parameter :: SIGP = 0.000394

	! parameters of the inventory
	integer,  parameter :: DIMQ = 30
	real(dp), parameter :: QUP = 100
	real(dp), parameter :: QLOW = 0
	
	! parameters of the demand
	real(dp), parameter :: DUP = 500
	real(dp), parameter :: DLOW = 0.001 
	real(dp), parameter :: APHA0 = 5.70
	real(dp), parameter :: APHA1 = 0.7
	real(dp), parameter :: SIGD = 0.75

	! parameters of the cost
	real(dp), parameter :: CUP = 7
	real(dp), parameter :: CLOW = 1
	real(dp), parameter :: SIGC = 0.1
	real(dp), parameter :: MUC = 4

	! OTHERS
	integer,  parameter :: NDRAWS = 2000
	real(dp), parameter :: DRATE = 0.995
	real(dp), parameter :: PBASE = 200
	real(dp), parameter :: PENALTY = 0
	real(dp), parameter :: HCOST = 2
	real(dp), parameter :: PCAP  = 100
	real(dp), parameter :: PI = 3.1415926
	 
end module sharedata

module mysub
	
	! This module contains five procedures:
	! dstate: to generate the discretized random space;
	! pmat:   to generate the tranistion probability matrix over the discretized random space;
	! ufucntion: to evaluate the one-period profit function
	! vfunction: to evaluate the value function
	! optfunc:   to define the objective function for maixmization in policy improvement step
	    
	!! Jia Yan, Dec. 06, 2004, Hong Kong Polytechnic University
	implicit none
	contains

		subroutine dstate (mout)
			! The output of this procedure is a DIMP*DIMQ by 2 vector, in which
			! the first column is the price, and the 2nd column is the inventory.
			! The random space is discretized randomly, with the following format
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! p1	z1  
			! p2	z1
			! p1	z2
			! p2	z2

			! This one is fixed (Jia, 13/12/2004)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use sharedata 
			use matforF90
			use numerical_libraries 
			implicit none

			! Define inputs and outputs
			real(dp), intent(out), dimension(DIMQ * DIMP, 2) :: mout 
	
			! Define locals
			integer :: i, num, kstart, kend
			real(dp), dimension(DIMP-2,1) :: randx, junkx
			real(dp), dimension(DIMQ-2,1) :: randy, junky
			real(dp), dimension(DIMP, 1)  :: disx
			real(dp), dimension(DIMQ, 1)  :: disy 
			! Code starts here
	
			! initializing
			mout = 0.
			disx = 0.
			disy = 0.

			! draw 
			call drnun(DIMP-2, randx)
			call drnun(DIMQ-2, randy)
			randx = PLOW + (PUP - PLOW) * randx
			randy = QLOW + (QUP - QLOW) * randy

			call dsvrgn(DIMP-2, randx, junkx)
			call dsvrgn(DIMQ-2, randy, junky)

			! to get final discretization
			disx(1, 1)			= PLOW
			disy(1, 1)          = QLOW
			disx(2:(DIMP-1), 1) = junkx(:,1)
			disy(2:(DIMQ-1), 1) = junky(:,1)
			disx(DIMP,1) = PUP
			disy(DIMQ,1) = QUP  

			num = 0
			do i = 1, DIMQ
				kstart = num + 1
				kend  =  num + DIMP
				mout(kstart:kend, 1) = disx(:,1)
				mout(kstart:kend, 2) = disy(i,1)
				num = num + DIMP
			end do

		end subroutine dstate

		subroutine pmat(index, order, smat, sdraw, mout)
			! To genrate the transition probability from one state to all others.
			! The inputs are:
			! index -- integer indicating the posiiton of initial state in discrete state space;
			! order -- deicsion on the state;
			! smat  -- DIMP*DIMQ by 2 matrix representing the discretized random space;
			! sdraw -- NDRAWS by 2 matrix; the 1st column is the random draw for price, and
			!          the 2nd column is the random draw for inventory, all conditional on the
			!          indexed initial state; 

			! The output is:
			! mout -- DIMP*DIMQ by 1 vector.  

			! Jia, 13/12/2004
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! Declare inputs and outputs
			integer,  intent(in)						  :: index	
			real(dp), intent(in)						  :: order
			real(dp), intent(in),  dimension(DIMP*DIMQ,2) :: smat
			real(dp), intent(in),  dimension(NDRAWS, 2)   :: sdraw
			real(dp), intent(out), dimension(DIMP*DIMQ,1) :: mout
	
			! Declare locals
			integer :: i, j, num, kstart, kend
			real(dp) :: cut0, cut1, tscalar 
			real(dp), dimension(DIMP,1) :: marginp
			real(dp), dimension(DIMQ,1) :: marginq
			real(dp), dimension(DIMP*DIMQ, 2) :: junk
			real(dp), dimension(NDRAWS,1) :: pdraw, ddraw 
			logical,  dimension(NDRAWS,1) :: ind
		
			! initialization
			mout = 0.	
			ind = .FALSE.
			pdraw(:,1) = sdraw(:,1)
			ddraw(:,1) = sdraw(:,2)
	
			! We first calculate the marginal transition probability of price
			do i = 1, DIMP
				if (i == 1) then
					cut0 = 0
					cut1 = smat(i,1)
				else
					cut0 = smat(i-1,1)
					cut1 = smat(i,1)
				end if
		
				forall (j = 1:NDRAWS, pdraw(j,1) <= cut1 .AND. pdraw(j,1) > cut0 )
					ind(j,1) = .TRUE.
				end forall
	
				marginp(i,1) =  count(ind)
				ind = .FALSE.
			end do
			marginp(:,1) = marginp(:,1) / real(NDRAWS)
			  
			! Then we calculate the conditional transition probability of inventory, and
			! thus the joint transition probability
			ddraw = smat(index, 2) + order - ddraw

			forall (i = 1:NDRAWS, ddraw(i,1) < 0)
					ddraw(i,1) = 0
			end forall

			num = 0
			junk = 0
			do i = 1, DIMQ
				kstart = num + 1
				kend   = num + DIMP
				
				! Get the conditional transition probability of inventory
				if (i == 1) then
					cut0 = -1
					cut1 = smat(1,2)
				else
					cut0 = smat(num,2)
					cut1 = smat(num+1,2)
				end if
		
				forall (j = 1:NDRAWS, ddraw(j,1) <= cut1 .AND. ddraw(j,1) > cut0 )
					ind(j,1) = .TRUE.
				end forall
	
				tscalar =  count(ind)
				tscalar = tscalar / real(NDRAWS)

				! The product of conditional and marginal is the joint transition probability
				mout(kstart:kend,1) = tscalar * marginp(:,1)
				ind = .FALSE.
				num = num + DIMP
			end do

		end subroutine pmat

		subroutine ufunction(order, zstate, sdraw, srent, uval)
			! To evaluate one-period profit function at a given state. The inputs are:
			! order  -- the deicsion on the state;
			! zstate -- inventory level of current state
			! sdraw  -- NDRAWS by 3 matrix including random draws for price, inventory, and cost
			! srent  -- the payment to the supplier
			! purpose -- 1 to calculate the optimal order at a given state and 0 to calculate (S,s)
			!			 bands

			! The output is:
			! uval -- the current profit at the state   

			!! Jia Yan, Jan. 31, 2005
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! declare inputs and outputs
			real(dp), intent(in)					  :: order, zstate, srent
			real(dp), intent(in), dimension(NDRAWS,3) :: sdraw
			real(dp), intent(out)					  :: uval

			! declare locals
			integer :: i, j, tpoint
			real(dp) :: meanp, cut, tscalar, pstockout, meandup, meandlow, wpart
			real(dp), dimension(NDRAWS,1)    :: pdraw, ddraw, udraw, stock
			logical,  dimension(NDRAWS,1)    :: ind

			! Body starts here -- assigning random draws 
			pdraw(:,1) = sdraw(:,1)
			ddraw(:,1) = sdraw(:,2)
			udraw(:,1) = sdraw(:,3)
			
			! get the mean price
			meanp = sum(pdraw)
			meanp = meanp / real(NDRAWS)
			
			! get the mean quantity demanded above/below a certain stock level   		
			cut = zstate + order
			ind = .FALSE.
			stock = 0.
			forall (i = 1:NDRAWS, ddraw(i,1) > cut)
				ind(i,1) = .TRUE.
!				stock(i,1) = ddraw(i,1)
				stock(i,1) = ddraw(i,1) - cut
			end forall
			tpoint = count(ind)
			pstockout = real(tpoint) / real(NDRAWS)

			if (pstockout > 0) then
				meandup = sum(stock)
				meandup = meandup / real(tpoint)
			else
				meandup = 0.	
			end if

			if (cut == 0) then
				meandlow = 0
			else
				ind = .FALSE.
				stock = 0.
				forall(i=1:NDRAWS, ddraw(i,1) <= cut)
					ind(i,1) = .TRUE.
					stock(i,1) = ddraw(i,1)
				end forall	
				tpoint = count(ind)
				if (tpoint == 0) then
					meandlow = 0
				else
					meandlow = sum(stock)
					meandlow = meandlow / real(tpoint)
				end if
			end if 

			if (pstockout > 1) then
				write (*,*) "wrong"
			end if

			! Now we can get the wpart
			wpart = meanp * meandlow + meanp * pstockout * cut - pstockout * PENALTY - &
			        HCOST * cut

			if (order > 0) then
				uval = (wpart - srent * order) * 10 - PBASE
			else
				uval = wpart * 10
			end if

		end subroutine ufunction

		subroutine vfunction(ntimes, pmat, u, value)
			! To approximate the value function using a geometric series. 
			! The inputs are:
			! ntimes -- integer indicating the # of terms of the series;
			! pmat      -- DIMP*DIMQ by DIMP*DIMQ matrix representing the transition matrix
			! u      -- DIMP*DIMQ by 1 vector representing the one-period profit

			! The output is:
			! value  -- DIMP*DIMQ by 1 vector 
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! Define input and output
			integer,  intent(in)                                   :: ntimes
			real(dp), intent(in),  dimension(DIMP*DIMQ,1)          :: u
			real(dp), intent(in),  dimension(DIMP*DIMQ, DIMP*DIMQ) :: pmat 
			real(dp), intent(out), dimension(DIMP*DIMQ,1)          :: value

			! Define locals
			integer :: i  
			real(dp), dimension(DIMP*DIMQ,1) :: tv
			real(dp), dimension(DIMP*DIMQ, DIMP*DIMQ) :: junk

			! Body starts here
			junk = pmat
			value = u
			do i = 1, ntimes
				if (i == 1) then
					tv = pmat .tx. u
					value = value + DRATE * tv
				else
					junk = junk .x. pmat	
					junk = ( DRATE**(i) ) * junk
					tv = junk .tx. u
					value = value + tv
				end if 
			end do
		end subroutine vfunction

end module mysub


program main
	! This program implements the policy iteration algorithm to solve a 
	! dynamic programming model
	use sharedata
	use mysub
	use matforF90
	use numerical_libraries
	use tools
	implicit none

	integer,  parameter :: iftemp = 21
	integer,  parameter :: ifout1 = 22
	integer,  parameter :: ifout2 = 23
	integer,  parameter :: ifout3 = 24
	integer,  parameter :: nsums  = 50
	integer,  parameter :: maxfcn = 2010
	real(dp), parameter :: tolerance = 1 * 10 ** (-3)
	
	integer  :: i, j, k, m, n, loop, group, ncounter, num, inter
	integer,  dimension(NDRAWS,1) :: sindex, intp, intz
	integer,  dimension(1)	  :: location		
	logical,  dimension(DIMP*DIMQ,1) :: ind
	real(dp), dimension(NDRAWS,1) :: mu, vp, vd, pt1, pt2, dt1, dt2, pdraw, ddraw, udraw, zlevel, &
									 mcost, vcost, ct1, ct2, hrate 	
	real(dp), dimension(NDRAWS,2) :: statedraw
	real(dp), dimension(NDRAWS,3) :: rdraw
	real(dp), dimension(DIMP*DIMQ, 1) :: order, ordernew, u, v, vjunk, xdiff, svalue
	real(dp), dimension(DIMP*DIMQ, 2) :: dsmat
	real(dp), dimension(DIMP*DIMQ,DIMP*DIMQ) :: tmat
	real(dp), dimension(maxfcn, DIMP) :: junk_ss
	real(dp), dimension(DIMP, 3) :: ssband
	real(dp), dimension(DIMP, 1) :: obj
	real(dp), dimension(maxfcn) :: junk, jval
	real(dp) :: cutpl, cutpu, cutzl, cutzu, action, z, tscalar, f, aub, atry, cutc, tvar1, &
				tvar2, tvar3, payment
	logical :: check
	
	! Body starts here
	vp = SIGP
	vd = SIGD
	pt1 = log(PLOW)
	pt2 = log(PUP)
	dt1 = log(DLOW)
	dt2 = log(DUP)
	ct1 = CLOW
	ct2 = CUP
	mcost = MUC
	vcost = SIGC

	open(iftemp,FILE='rjunk.dat',STATUS='UNKNOWN',FORM='UNFORMATTED', RECL=80000) 
	open(ifout1, FILE='tdata_HCOST_2_SIGC_01_PBASE_200_PCAP_100.dat',STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE', RECL=80000) 
	open(ifout2, FILE='rdraws.dat',STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE', RECL=80000) 
	open(ifout3, FILE='ssband.dat',STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE', RECL=80000) 


	! Discretize the continues random space, and write to disk for future use
	call dstate(dsmat)

	! generate random draws for production cost from a truncated normal distribution
	call tnormal_2side_vec(NDRAWS, mcost, vcost, ct1, ct2, udraw)

	! calculate the hazard rate of truncated normal distribution
	tvar1 = (CLOW - MUC) / sqrt(SIGC)
	tvar2 = (CUP -  MUC) / sqrt(SIGC)
	tvar3 = dnordf(tvar2) - dnordf(tvar1)
	hrate = 0.
	do i = 1, NDRAWS
		cutc = udraw(i,1)
		tvar1 = (cutc - MUC) / sqrt(SIGC)
		ncounter = 0
		do k = 1, NDRAWS
			if (udraw(k,1) <= cutc) then
				ncounter = ncounter + 1
			end if
		end do
		tscalar = ( 1 / sqrt(2*PI*SIGC) ) * exp(-0.5 * (tvar1 ** 2 ) )
		tscalar = tscalar / tvar3
		hrate(i,1) = real(ncounter)
		hrate(i,1) = hrate(i,1) / real(NDRAWS)
		hrate(i,1) = hrate(i,1) / tscalar
	end do
	payment = sum(hrate + udraw)
	payment = payment / real(NDRAWS)
	print *, payment

	! Set initial values for optimal order
	order(:,1) = 0.
	ssband = 0.

	! Start the policy iteration
	loop = 0
	do 
		loop = loop + 1
		print *, loop

		if (loop > 30 )  then
			write(*,*) 'maximum number of func evaluation exceeded'
			exit
		end if

		! Loop over state space to get the transition probability, and u-value at each state		
		num = 0
		do i = 1, DIMP
			num = num + 1

			if (loop == 1) then
				mu = LMD0 + LMD1 * log(dsmat(i,1))   
				call tnormal_2side_vec(NDRAWS, mu, vp, pt1, pt2, pdraw)
				pdraw = exp(pdraw)

			
				! mapping price draws into the discrete price categories, or states of price
				intp = 0
				do j = 1, DIMP
					if (j == 1) then
						cutpl = 0
						cutpu = dsmat(j, 1)
					else
						cutpl = dsmat(j-1, 1)
						cutpu = dsmat(j, 1) 
					end if

					forall (k = 1:NDRAWS, pdraw(k,1) > cutpl .AND. pdraw(k,1) <= cutpu)
						intp(k,1) = j
					end forall
				end do

				! random draws of demand conditional on drawed price
				mu = APHA0 - APHA1 * log(pdraw)
				call tnormal_2side_vec(NDRAWS, mu, vd, dt1, dt2, ddraw)
				ddraw = exp(ddraw)

				! write random draws to the disk for future use
				write(iftemp) pdraw, ddraw, intp
			else
				read(iftemp) pdraw, ddraw, intp
			end if
			statedraw(:,1) = pdraw(:,1)
			statedraw(:,2) = ddraw(:,1)

			ncounter = 0
			do j = 1, DIMQ
				inter = num + ncounter
				action = order(inter, 1)
				z = dsmat(inter, 2)

				! Get the transition probability matrix
				call pmat(inter, action, dsmat, statedraw, vjunk)
				tmat(:, inter) = vjunk(:,1)


				rdraw(:,1) = pdraw(:,1)
				rdraw(:,2) = ddraw(:,1)
				rdraw(:,3) = udraw(:,1)

				! Get the values of one-period profit for each state
				call ufunction(action, z, rdraw, payment, tscalar)
				u(inter,1) = tscalar

				ncounter = ncounter + DIMP
			end do

		end do
		rewind iftemp
			
		! Now we can approximate the value function
		call vfunction(nsums, tmat, u, v)

		!!! Policy improvement step for calculating optimal order for each state !!! 
		num = 0
		do i = 1, DIMP
			num = num + 1

			! Read in the data of random draws
			read(iftemp) pdraw, ddraw, intp
			rdraw(:,1) = pdraw(:,1)
			rdraw(:,2) = ddraw(:,1)
			rdraw(:,3) = udraw(:,1) 

			ncounter = 0
			do j = 1, DIMQ
				inter = num + ncounter
	
				z = dsmat(inter, 2)
				if (PCAP > QUP - z) then
					aub = QUP - z
				else
					aub = PCAP
				end if


				atry = 0
				junk = -1000000
				jval = 0

				do k = 1,maxfcn
					if (k == 1) then
						atry = 0.
					else
						atry = atry + 0.1
					end if

					if (atry > aub + 0.1) exit

					if (atry >= aub) then
						atry = aub
					end if

					! Get the one-period profit at the new order
					call ufunction(atry, z, rdraw, payment, tscalar)

					!!!! Update the transition probability under the new order !!!!
					! Get inventory level given the realizations of random demand
					zlevel = z + atry - ddraw
					! The inventory can not be negative
					forall (m = 1:NDRAWS, zlevel(m,1) < 0)
						zlevel(m,1) = 0
					end forall

					! mapping the inventory levels into discrete categories, or inventory states
					intz = 0
					group = 0
					do m = 1, DIMQ
						if (m == 1) then
							cutzl = -1
							cutzu = dsmat(group+1, 2)
						else
							cutzl = dsmat(group, 2)
							cutzu = dsmat(group+1, 2) 
						end if


						forall (n = 1:NDRAWS, zlevel(n,1) > cutzl .AND. zlevel(n,1) <= cutzu)
							intz(n,1) = m
						end forall

						group = group + DIMP
					end do

					sindex = intp + (intz - 1) * DIMP
				
					!!! Calculate the average value function !!!
					f = 0
					do m = 1, NDRAWS
						f = f + v(sindex(m,1),1)
					end do	
					f = f / NDRAWS
					junk(k) = tscalar + DRATE * f
					jval(k) = atry

				end do
			
				!!! Calculate (S,s) bands !!!
				if (j == 1) then
					junk_ss(:,i) = junk
				end if

				!!! Define optimal qunatity if the minimum of optimals !!!
				f = maxval(junk)
				do k = 1, maxfcn
					if (junk(k) >= f) then
						ordernew(inter, 1) = jval(k)
						exit
					end if
				end do

				!!! Calculate (S,s) bands !!!
				if (j == 1) then
					obj(i,1) = f
				end if 
				ncounter = ncounter + DIMP
			end do

		end do
		rewind iftemp

		! Check convergency
		xdiff = abs(ordernew - order)
		if (maxval(xdiff) <= tolerance) then
			write(*,*) 'policy iteration converges successfully'

			!!!!Caluclate the profit of the supplier
			tscalar = sum(hrate)
			tscalar = tscalar / NDRAWS 
			svalue = (PBASE + ordernew * tscalar) / (1-DRATE)
			do i = 1, DIMP*DIMQ
				write(ifout1,*) dsmat(i,1),achar(9),dsmat(i,2),achar(9),u(i,1),achar(9),v(i,1), &
							   achar(9), ordernew(i,1), achar(9), svalue(i,1)
			end do

			do i = 1, NDRAWS
				write(ifout2,*) udraw(i,1), achar(9), hrate(i,1)
			end do

			!! Calculate (S,s) bands !!
			obj = obj - PBASE
			ssband(:,1) = dsmat(1:DIMP,1)
			ssband(:,2) = ordernew(1:DIMP,1)
			do i = 1, DIMP
				atry = 0.
				do j = 1, maxfcn
					if (j == 1) then
						atry = 0.
					else
						atry = atry + 0.1
					end if

					if (junk_ss(j, i) >= obj(i,1)) then
						ssband(i, 3) = atry

						exit
					end if
				end do

				!! Write ssband data onto the hard disk !!!
				write(ifout3,*) ssband(i,1), achar(9), ssband(i,2), achar(9), ssband(i,3)

			end do
		
			call prnt(ssband)
			exit
		end if
		
		print *, maxval(xdiff), tolerance			
		order = ordernew
	end do

end program main
