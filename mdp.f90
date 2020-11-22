module sharedata

	! Define the globals to be used in later procedures

	! Jia Yan
	use matforF90
	implicit none
	save

	! State parameters
	real(dp) :: kstate, bstate, estate, brate

	! parameters of the operating profit function
	real(dp), parameter :: apha0 = 0.
	real(dp), parameter :: apha1 = 0.75
	real(dp), parameter :: rho = 0.66
	real(dp), parameter :: sigmu = 0.1221

	! parameters of the capital adjustment function
	real(dp), parameter :: lamd0 = 0.039
	real(dp), parameter :: lamd1 = 0.049
	
	! parameters of the equity cost function
	real(dp), parameter :: theta0 = 0.389
	real(dp), parameter :: theta1 = 0.053
	real(dp), parameter :: theta2 = 0.0002

	! OTHERS
	real(dp), parameter :: taxc = 0.30
	real(dp), parameter :: taxd = 0.30
	real(dp), parameter :: drate = 0.995
	real(dp), parameter :: rr = 0.04
	real(dp), parameter :: rc = 0.10 
	real(dp), parameter :: dpre = 0.15
	integer,  parameter :: ndraws = 500
	integer,  parameter :: kdim = 3
	integer,  parameter :: bdim  = 3
	integer,  parameter :: edim = 3
	real(dp), parameter :: kup = 1000
	real(dp), parameter :: klow = 0
	real(dp), parameter :: bup = 500
	real(dp), parameter :: blow = -500
	real(dp), parameter :: eup = 0.65
	real(dp), parameter :: elow = -0.65
	integer,  parameter :: tperiods = 1000

	! Random draws for Monte-Carlo Integration
	real(dp), dimension(ndraws,1) :: edraw

	! States
	real(dp), dimension(kdim*bdim*edim, 3) :: state

	! value-to-go vector
	real(dp), dimension(kdim*bdim*edim, 1) :: vvec
		 
end module sharedata

module mysub1
	
	! This module contains five procedures:
	! dstate: to generate the discretized random space;
    
	!! Jia Yan, Feb. 06, 2008, WSU
	implicit none
	contains

		subroutine dstate (mout)
			! The output of this procedure is a kdim*bdim*edim by 3 vector, in which
			! the first column is capital; the 2nd column is debt, and the third column is shock.
			! The random space is discretized randomly, with the following format
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! k1	b1  e1
			! k2	b1	e1
			! k1	b2	e1
			! k2	b2	e1
			! k1	b1  e2
			! k2	b1	e2
			! k1	b2	e2
			! k2	b2	e2
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			use sharedata 
			use matforF90
			use numerical_libraries 
			implicit none

			! Define inputs and outputs
			real(dp), intent(out), dimension(kdim * bdim * edim, 3) :: mout 
	
			! Define locals
			integer :: i, j, num1, num2, kstart, kend
			real(dp), dimension(kdim-2,1) :: randk, disk
			real(dp), dimension(bdim-2,1) :: randb, disb
			real(dp), dimension(edim-2,1) :: rande, dise

			real(dp), dimension(kdim,1) :: kpoints
			real(dp), dimension(bdim,1) :: bpoints
			real(dp), dimension(edim,1) :: epoints

			! initializing
			mout = 0.
			disk = 0.
			disb = 0.
			dise = 0.

			! draw 
			call drnun(kdim-2, randk)
			call drnun(bdim-2, randb)
			call drnun(edim-2, rande)
			randk = klow + (kup - klow) * randk
			randb = blow + (bup - blow) * randb
			rande = elow + (eup - elow) * rande

			call dsvrgn(kdim-2, randk, disk)
			call dsvrgn(bdim-2, randb, disb)
			call dsvrgn(edim-2, rande, dise)

			kpoints(1,1) = klow
			kpoints(kdim,1) = kup
			kpoints(2:(kdim-1),1) = disk(:,1)

			bpoints(1,1) = blow
			bpoints(bdim,1) = bup
			bpoints(2:(bdim-1),1) = disb(:,1)

			epoints(1,1) = elow
			epoints(edim,1) = eup
			epoints(2:(edim-1),1) = dise(:,1)

			! order
			num1 = 0
			do i = 1, edim
				do j = 1, bdim
					kstart = num1 + 1
					kend  =  num1 + kdim
					mout(kstart:kend, 1) = kpoints(:,1)
					mout(kstart:kend, 2) = bpoints(j,1)
					mout(kstart:kend, 3) = epoints(i,1)
				
					num1 = num1 + kdim
				end do
			end do

		end subroutine dstate

end module mysub1

subroutine dividend(dim,x ,fval)
	
	use matforf90
	use mysub1
	use sharedata
	implicit none

	integer,  intent(in) :: dim
	real(dp), intent(in), dimension(dim) :: x
	real(dp), intent(out) :: fval

	integer :: i, j
	real(dp) :: y, c, scalar, knew, bnew, enew, r

	! current period
	y = apha0 + apha1 * log(kstate) + estate
	y = exp(y)

	knew = x(1)
	bnew = x(2)
	if (bnew > 0) then
		r = rr
	else
		r = rc
	end if

	scalar = knew - (1-dpre) * kstate

	if (abs(scalar) > 0.0001) then
		c = (lamd0 + 0.5 * lamd1 * (scalar/kstate)**2) * kstate
	else
		c = 0.
	end if

	fval = (1-taxc) * y - c - scalar + bstate - (bnew/(1+r*(1-taxc)) )

end subroutine dividend

subroutine objlast(dim, x, fval)

	use matforf90
	use mysub1
	use sharedata
	implicit none

	integer,  intent(in) :: dim
	real(dp), intent(in), dimension(dim) :: x
	real(dp), intent(out) :: fval

	integer :: i, j  
	real(dp) :: y, c,  div, knew, bnew, enew, v, vsum, scalar, r

	! current period	
	y = apha0 + apha1 * log(kstate) + estate
	y = exp(y)
	knew = x(1)
	bnew = x(2)
	if (bnew > 0) then
		r = rr
	else
		r = rc
	end if

	scalar = knew - (1-dpre)* kstate
	if (abs(scalar) > 0.00001) then
		c = (lamd0  + 0.5 * lamd1 * (scalar/kstate)**2) * kstate
	else
		c = 0
	end if


	div = (1-taxc) * y - c - scalar + bstate - (bnew/(1+r*(1-taxc)) )
	if (div < 0) then
		div = div - (theta0 - theta1 * div + 0.5 * theta2*(div**2))
	end if

	! next period
	! Monte-Carlo interation
	vsum = 0.	
	do i = 1, ndraws
		enew = rho * estate + sigmu * edraw(i,1)

		y = apha0 + apha1 * log(knew) + enew
		y = exp(y)
		v = y + knew + bnew
		vsum = vsum + v
	end do
	vsum = vsum / ndraws

	fval = div + brate * drate * vsum
end subroutine objlast

subroutine objt(dim, x, xindex, fval)

	use matforf90
	use mysub1
	use sharedata
	implicit none

	integer,  intent(in) :: dim
	real(dp), intent(in), dimension(dim) :: x
	integer,  intent(in), dimension(dim) :: xindex
	real(dp), intent(out) :: fval 

	integer :: i, j, intk, intb, inte, sindex, group 
	real(dp) :: y, c, div, knew, bnew, enew, vsum, cutl, cutu, diff1, diff2, scalar, r

	! current period	
	y = apha0 + apha1 * log(kstate) + estate
	y = exp(y)
	
	if (x(2) > 0) then
		r = rr
	else
		r = rc
	end if

	scalar = x(1) - (1-dpre)*kstate
	if (abs(scalar) > 0.0001) then
			c = (lamd0 + 0.5 * lamd1 * (scalar/kstate)**2) * kstate
	else
			c = 0.
	end if
	div = (1-taxc) * y - c - scalar + bstate - (x(2)/(1+r*(1-taxc)) )
	if (div < 0) then
		div = div - (theta0 - theta1 * div + 0.5 * theta2 * (div**2)) 
	end if

	! next period
	intk = xindex(1)
	intb = xindex(2)
	! Monte-Carlo interation
	vsum = 0.	
	do i = 1, ndraws
		enew = rho * estate + sigmu * edraw(i,1)

		group = 0
		inte = 0
		do j = 1, edim

			if (j==1) then
				cutl = -100000000000
				cutu = state(1,3)
			else
				cutl = state(group-kdim*bdim+1, 3)
				cutu = state(group+1, 3) 
			end if 	

			if (enew > cutl .AND. enew <= cutu) then
				diff1 = enew - cutl
				diff2 = cutu - enew
				if (diff1 > diff2) then
					inte = j
				else 
					inte = j-1
				end if
				exit
			end if

			group = group + kdim*bdim
		end do

		if (inte == 0) then
			inte = edim			
		end if	
		
		sindex = intk + (intb-1)*kdim + (inte-1)*kdim*bdim
		vsum = vsum + vvec(sindex,1)
	end do
	vsum = vsum / ndraws
			
	fval = div + drate * brate * vsum

end subroutine objt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main

	use matforf90
	use numerical_libraries
	use mysub1
	use sharedata
	use tools
	implicit none

	integer,  parameter :: ifstate = 23
	integer,  parameter :: ifout = 24 
	integer,  parameter :: nvar = 2

	logical :: check			
	integer :: iloop, i, num, maxfcn, ibtype, j, group	
	real(dp) :: fout, tol, mu, vp, rscalar, lbd, ubd, x1, c, y
	integer,  dimension(2) :: location, locationjunk
 	real(dp), dimension(nvar) :: xguess, xout, xlb, xub 
	real(dp), dimension(1) :: bini, bout, blb, bub
	real(dp), dimension(kdim*bdim*edim,8) :: results
	real(dp), dimension(kdim*bdim*edim,1) :: value
	real(dp), dimension(kdim*bdim,1) :: sortD
	real(dp), dimension(kdim*bdim,2) :: indjunk 

	! Body starts here	
	check = .FALSE.

	mu = 0.
	vp =0.1221 
	open(ifstate,FILE='sjunk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
	open(ifout,FILE='results.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
	
	! genearet random draws used in subroutines to do Monte-Carlo integration; edraw is a macro used by other procedures	
	edraw = 0.
	do i = 1, ndraws
		call tnormal_2side_univariate(mu, vp, elow, eup, rscalar)
		edraw(i,1) = rscalar
	end do

	! discretize the state space randomly; state is the macro used by other procedures
	! also, get the data points for interpolation	
	call dstate(state)
	call prnt(state)

	! initializing value-to-go vector
	vvec = 0.

	debug: if (check) then
	! starting the bakward induction
	results = 0.
	Time: do iloop = 1, tperiods					
		print *, "TIME", iloop 

		! starting from the last period	
		if (iloop == 1) then
			innerlast: do i = 1, kdim*bdim*edim

				! Write states onto disk for checking 	
				if (iloop == 1) then
						write(ifstate,*) state(i,1), achar(9), state(i,2), achar(9), state(i,3)
				end if

				! Assigning values of state macros
				brate = 0.75
				kstate = state(i,1)
				bstate = state(i,2)
				estate = state(i,3)

				do j = 1, kdim*bdim
					xguess(1) = state(j,1)
					xguess(2) = state(j,2)
					call objlast(nvar, xguess, fout)
					sortD(j,1) = fout
				end do

				! finding maximum
				location = maxloc(sortD)
				xguess(1) = state(location(1), 1)
				xguess(2) = state(location(1), 2)
				
				! now calculate value function recursively
				brate = 1.
				call objlast(nvar, xguess, fout)
				vvec(i,1) = fout
				
				! calculate dividend at the solutions
				call dividend(nvar, xguess, fout)
				 
				results(i,1) = kstate
				results(i,2) = bstate
				results(i,3) = estate
				results(i,4) = vvec(i,1)
				results(i,5) = xguess(1)	
				results(i,6) = xguess(2)
				results(i,7) = fout
				results(i,8) = iloop
				
			end do innerlast
		else
			innerperiod: do i = 1, kdim*bdim*edim

				! Assigning values of state macros
				brate = 0.75
				kstate = state(i,1)
				bstate = state(i,2)
				estate = state(i,3)

				group = 0
				do j = 1, bdim
					do num = 1, kdim
						group = group + 1
						xguess(1) = state(group,1)
						xguess(2) = state(group,2)
						location(1) = num
						location(2) = j
						call objt(nvar, xguess, location, fout)

						sortD(group,1) = fout
						indjunk(group,1) = num
						indjunk(group,2) = j 
					end do

				end do

				! finding maximum
				location = maxloc(sortD)
				xguess(1) = state(location(1),1)
				xguess(2) = state(location(1),2)
				locationjunk(1) = indjunk(location(1),1)
				locationjunk(2) = indjunk(location(1),2)

				! now calculate value function recursively
				brate = 1.
				call objt(nvar, xguess, locationjunk, fout)
				value(i,1) = fout

				! calculate dividend at the solutions
				call dividend(nvar, xguess, fout)
				 
				results(i,1) = kstate
				results(i,2) = bstate
				results(i,3) = estate
				results(i,4) = vvec(i,1)
				results(i,5) = xguess(1)
				results(i,6) = xguess(2)	
				results(i,7) = fout
				results(i,8) = iloop

			end do innerperiod

			! update vvec
			vvec = value			
		end if
		
			
	end do Time		

	! write out results
	do iloop = 1, kdim*bdim*edim
		write(ifout,*) ((results(iloop, i), achar(9)), i = 1, 8) 
	end do

	end if debug
end program main
