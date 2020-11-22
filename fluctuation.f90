module sharedata

	! To define the globals to be used in later procedures

	! Jia Yan, 11/23/2005
	use matforF90
	implicit none
	save

	real(dp), parameter :: PRICE = 20.
	real(dp), parameter :: COST = 10.

	! parameters of the inventory
	integer,  parameter :: DIMQ = 20
	real(dp), parameter :: QUP = 100.
	real(dp), parameter :: QLOW = -100.
	
	! parameters of the demand
	real(dp), parameter :: DUP = 100.
	real(dp), parameter :: DLOW = 0.

	! OTHERS
	integer,  parameter :: NDRAWS = 2000
	real(dp), parameter :: DRATE = 0.995
	real(dp), parameter :: PENALTY = 5 - DRATE * (PRICE - COST) 
	real(dp), parameter :: HCOST = 1.
	real(dp), parameter :: PI = 3.1415926535

end module sharedata

module mysub
	
	! This module contains five procedures:
	! dstate: to generate the discretized random space;
	! pmat:   to generate the tranistion probability matrix over the discretized random space;
	! ufucntion: to evaluate the one-period profit function
	! vfunction: to evaluate the value function
	! optfunc:   to define the objective function for maixmization in policy improvement step
	    
	!! Jia Yan, Dec. 06, 2005, Hong Kong Polytechnic University
	implicit none
	contains

		subroutine dstate (mout)
			! The output of this procedure is a DIMQ-1 by 1 vector, representing inventory.
			! The random space is discretized randomly, with the following format

			! This one is fixed (Jia, 11/23/2005)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use sharedata 
			use matforF90
			use numerical_libraries 
			implicit none

			! Define inputs and outputs
			real(dp), intent(out), dimension(DIMQ, 1) :: mout 
	
			! Define locals
			integer :: i, num, kstart, kend
			real(dp), dimension(DIMQ-1,1) :: randy, junky

			! Code starts here

			! initializing
			mout = 0.

			! draw 
			call drnun(DIMQ-1, randy)
			randy = QLOW + (QUP - QLOW) * randy

			! sort in descending order	
			call dsvrgn(DIMQ-1, randy, junky)

			mout(1:DIMQ-1,1) = junky(:,1)
			mout(DIMQ,1) = QUP

		end subroutine dstate

		subroutine pmat(index, order, smat, sdraw, mout)
			! To genrate the transition probability from one state to all others.
			! The inputs are:
			! index -- integer indicating the posiiton of initial state in discrete state space;
			! order -- deicsion on the state;
			! smat  -- DIMQ by 1 matrix representing the discretized random space;
			! sdraw -- NDRAWS by 1 matrix; the random draw for inventory, conditional on the
			!          indexed initial state; 

			! The output is:
			! mout -- DIMQ by 1 vector.  

			! Jia, 23/11/2005
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! Declare inputs and outputs
			integer,  intent(in)						  :: index	
			real(dp), intent(in)						  :: order
			real(dp), intent(in),  dimension(DIMQ,1)	  :: smat
			real(dp), intent(in),  dimension(NDRAWS, 1)   :: sdraw
			real(dp), intent(out), dimension(DIMQ,1)      :: mout
	
			! Declare locals
			integer :: i, j
			real(dp) :: cut, tscalar 
			real(dp), dimension(DIMQ,1) :: junk
			real(dp), dimension(NDRAWS,1) :: demand 
			logical,  dimension(NDRAWS,1) :: ind
		
			! initialization
			mout = 0.	
			ind = .FALSE.
	
			! calculate the conditional transition probability of inventory
			demand = smat(index, 1) + order - sdraw

			junk = 0
			do i = 1, DIMQ

				cut = smat(i,1)

				forall (j = 1:NDRAWS, demand(j,1) <= cut )
					ind(j,1) = .TRUE.
				end forall
	
				tscalar =  count(ind)
				tscalar = tscalar / real(NDRAWS)
				junk(i,1) = tscalar

				if (i == 1) then
					mout(i,1) = tscalar
				else
					mout(i,1) = tscalar - junk(i-1,1)
				end if

				ind = .FALSE.

			end do

		end subroutine pmat

		subroutine ufunction(order, zstate, d, uval)
			! To evaluate one-period profit function at a given state. The inputs are:
			! order  -- the deicsion on the state;
			! zstate -- inventory level of current state
			! ddraw  -- NDRAWS by 1 random draws for inventory

			! The output is:
			! uval -- the current profit at the state   

			!! Jia Yan, Nov. 23, 2005
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! declare inputs and outputs
			real(dp), intent(in)					  :: order, zstate
			real(dp), intent(in), dimension(NDRAWS,1) :: d
			real(dp), intent(out)					  :: uval

			! declare locals
			integer :: i, tpoint
			real(dp) :: cut, pstockout, meandlow, meanstock, meanout
			real(dp), dimension(NDRAWS,1)    :: stockout, stock
			logical,  dimension(NDRAWS,1)    :: ind

			! Body starts here -- assigning random draws 
			
			! get the mean quantity demanded above/below a certain stock level   		
			cut = zstate + order
			ind = .FALSE.
			stockout = 0.
			stock = 0.

			meandlow = 0.
			do i = 1, NDRAWS
			
				if (d(i,1) > cut) then
					ind(i,1) = .TRUE.
					stockout(i,1) = d(i,1) - cut 
				else
					meandlow = meandlow + d(i,1)  
					stock(i,1) = cut - d(i,1)
				end if
				
			end do

			tpoint = count(ind) ! # of stockout
			pstockout = real(tpoint) / real(NDRAWS)

			if (pstockout > 1) then
				write (*,*) "wrong: probability can't be greater than 1"
			end if

			if (pstockout == 1) then
				meanstock = 0
			else
				meandlow = meandlow / real(NDRAWS - tpoint)
				meanstock = sum(stock)
				meanstock = meanstock / real(NDRAWS - tpoint)
			end if

			if (pstockout == 0) then
				meanout = 0
			else
				meanout = sum(stockout)
				meanout = meanout / real(tpoint)
			end if

			uval = PRICE * (meandlow + pstockout * cut) - PENALTY * meanout - &
			        HCOST * meanstock - COST * order

		end subroutine ufunction

		subroutine vfunction(ntimes, pmat, u, value)
			! To approximate the value function using a geometric series. 
			! The inputs are:
			! ntimes -- integer indicating the # of terms of the series;
			! pmat      -- DIMQ by DIMQ matrix representing the transition matrix
			! u      -- DIMQ by 1 vector representing the one-period profit

			! The output is:
			! value  -- DIMQ by 1 vector 
			use sharedata
			use matforF90
			use numerical_libraries
			implicit none

			! Define input and output
			integer,  intent(in)                                   :: ntimes
			real(dp), intent(in),  dimension(DIMQ,1)			   :: u
			real(dp), intent(in),  dimension(DIMQ, DIMQ)		   :: pmat 
			real(dp), intent(out), dimension(DIMQ,1)               :: value

			! Define locals
			integer :: i  
			real(dp), dimension(DIMQ,1) :: tv
			real(dp), dimension(DIMQ, DIMQ) :: junk

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

	integer,  parameter :: ifout1 = 22
	integer,  parameter :: ifout2 = 23
	integer,  parameter :: nsums  = 50
	real(dp), parameter :: tolerance = 1
	
	integer  :: i, j, k, m, n, loop, num, ibtype, maxfcn
	logical,  dimension(DIMQ,1) :: ind
	real(dp), dimension(DIMQ, 1) :: order, ordernew, u, vjunk, xdiff, DSMAT, VVAL
	real(dp), dimension(DIMQ,DIMQ) :: tmat
	real(dp), dimension(NDRAWS,1)  :: DDRAW, sindex, zlevel
	real(dp) :: aub, alow, atry, action, z, tscalar, f, ftemp, quant, ZPOS, intercut

	logical :: check

	!external :: obj
	
	! Body starts here
	num = 1
	alow = 0.
	ibtype = 0
	maxfcn = 30000

	open(ifout1, FILE='tdata_p00.dat', STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE', RECL=80000) 
	open(ifout2, FILE='rdraws.dat',STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE', RECL=80000) 

	check = .TRUE.

	! Discretize the continues random space, and write to disk for future use
	call dstate(DSMAT)

	! to get the random demand draws
	call drnun(NDRAWS, DDRAW)
	DDRAW = DLOW + (DUP - DLOW) * DDRAW

	do i = 1, NDRAWS
		write(ifout2, *) DDRAW(i,1)
	end do

	if (check) then

	! Set initial values for optimal order
	order(:,1) = 0.

	! Start the policy iteration
	loop = 0
	do 
		loop = loop + 1
		print *, loop

		if (loop > 10 )  then
			write(*,*) 'maximum number of func evaluation exceeded'
			exit
		end if

		! Loop over state space to get the transition probability, and u-value at each state		
		do i = 1, DIMQ
			quant = order(i, 1)
			ZPOS = DSMAT(i, 1)

			! Get the transition probability matrix
			call pmat(i, quant, DSMAT, DDRAW, vjunk)
			tmat(:, i) = vjunk(:,1)

			! Get the values of one-period profit for each state
			call ufunction(quant, ZPOS, DDRAW, tscalar)
			u(i,1) = tscalar
		end do
		
		! Now we can approximate the value function
		call vfunction(nsums, tmat, u, VVAL)

		!!! Policy improvement step for calculating optimal order for each state !!! 
		do i = 1, DIMQ

			! Read in the data of random draws
			ZPOS = DSMAT(i, 1)
			aub = QUP - ZPOS  ! Define the contraint correspondence

			action = 0.
			ftemp = -99999999
			atry = 0.
			do j = 1, 300000000

				if(atry > aub) exit

				call ufunction(atry, ZPOS, DDRAW, tscalar)
	
				zlevel = ZPOS + atry - DDRAW

				! The inventory can not be negative
				!do m = 1, NDRAWS
				!	if (zlevel(m,1) < 0) then
				!		zlevel(m,1) = 0
				!	end if
				!end do

				! mapping the inventory levels into discrete categories, or inventory states
				sindex = 0
				do m = 1, DIMQ

					intercut = DSMAT(m,1)

					do n = 1, NDRAWS
						if (zlevel(n,1) <= intercut .AND. sindex(n,1) == 0) then
							sindex(n,1) = m
						end if
					end do
				end do

				!!! Calculate the average value function !!!
				f = 0
				do m = 1, NDRAWS
					f = f + VVAL(sindex(m,1),1)
				end do	
				f = f / NDRAWS
				f = tscalar + DRATE * f

				if (f > ftemp) then
					action = atry
					ftemp = f
				end if

				atry = atry + 0.0001
			end do

			!call dbcpol(obj, num, atry, ibtype, alow, aub, tolerance, maxfcn, action, ftemp)
			ordernew(i,1) = action
			print *, ZPOS, action
		end do

		! Check convergency and output results
		xdiff = abs(ordernew - order)
		if (maxval(xdiff) <= 0.001) then
			write(*,*) 'policy iteration converges successfully'

			do i = 1, DIMQ
				write(ifout1,*) DSMAT(i,1),achar(9), u(i,1),achar(9),VVAL(i,1), &
							   achar(9), ordernew(i,1)
			end do

			exit
		end if
		
		order = ordernew

		print *,"tolerance", maxval(xdiff), tolerance			

	end do

	end if

end program main


