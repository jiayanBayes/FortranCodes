module sharedata

	! To define the globals to be used in later procedures

	! Jia Yan
	use matforF90
	implicit none
	save

	! parameters of the operating profit function
	real(dp), parameter :: apha0 = -2.2898
	real(dp), parameter :: apha1 = 0.6242
	real(dp), parameter :: apha2 = 0.0735
	real(dp), parameter :: apha3 = -0.0045
	real(dp), parameter :: rho = 0.6242
	real(dp), parameter :: sigmu = 0.4339

	! parameters of the capital adjustment function
	real(dp), parameter :: lamd0 = 0.039
	real(dp), parameter :: lamd1 = 0.049
	
	! parameters of the interest function
	real(dp), parameter :: eta0 = 500
	real(dp), parameter :: eta1 = 0.001 
	real(dp), parameter :: eta2 = 0.01

	! parameters of the equity cost function
	real(dp), parameter :: theta0 = 0.389
	real(dp), parameter :: theta1 = 0.053
	real(dp), parameter :: theta2 = 0.0002

	! OTHERS
	real(dp), parameter :: brate = 0.9
	real(dp), parameter :: drate = 0.9
	real(dp), parameter :: dpre = 0.15
	integer,  parameter :: ndraws = 1000
	integer,  parameter :: kdim = 2
	integer,  parameter :: bdim  = 2
	integer,  parameter :: edim = 3
	real(dp), parameter :: kup = 2000000
	real(dp), parameter :: klow = 0
	real(dp), parameter :: bup = 500
	real(dp), parameter :: blow = 0
	real(dp), parameter :: eup = 100
	real(dp), parameter :: elow = 0
	integer,  parameter :: tperiods = 1
		 
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
			real(dp), dimension(kdim,1) :: randk, disk
			real(dp), dimension(bdim,1) :: randb, disb
			real(dp), dimension(edim,1) :: rande, dise

			! initializing
			mout = 0.
			disk = 0.
			disb = 0.
			dise = 0.

			! draw 
			call drnun(kdim, randk)
			call drnun(bdim, randb)
			call drnun(edim, rande)
			randk = klow + (kup - klow) * randk
			randb = blow + (bup - blow) * randb
			rande = elow + (eup - elow) * rande

			call dsvrgn(kdim, randk, disk)
			call dsvrgn(bdim, randb, disb)
			call dsvrgn(edim, rande, dise)

			! order
			num1 = 0
			do i = 1, edim
				do j = 1, bdim
					kstart = num1 + 1
					kend  =  num1 + kdim
					mout(kstart:kend, 1) = disk(:,1)
					mout(kstart:kend, 2) = disb(j,1)
					mout(kstart:kend, 3) = dise(i,1)
				
					num1 = num1 + kdim
				end do
			end do

		end subroutine dstate

end module mysub1

subroutine objlast(nv, x, f)
	
	use matforf90
	use mysub1
	use sharedata
	implicit none

	integer, intent(in) :: nv
	real(dp), intent(in), dimension(nv) :: x
	real(dp), intent(out) :: f

	integer, parameter:: iftemp = 11
	integer :: i, istat 
	real(dp) :: kstate, bstate, estate, y, c, equity, div, r, knew, bnew, enew, v, vsum
	real(dp), dimension(ndraws,1) :: edraw

	open(iftemp, FILE='junk.asc',STATUS='OLD', FORM='FORMATTED', ACTION='READ', IOSTAT=istat) 
	if (istat /= 0) then
		write(*,*) "file open failed"
	else
		read(iftemp, *) kstate, bstate, estate
	end if
	print *, kstate, bstate, estate 
	close(iftemp)
	kstate = 100000000
	bstate = 10000000
	estate = 0.1

	! current period	
	y = apha0 + apha1 * log(kstate) + apha2 * ((log(kstate)) ** 2) + apha3 * ((log(kstate)) ** 3) + estate
	y = exp(y)
	c = (lamd0 * x(1) + 0.5 * lamd1 * (x(1)/kstate)**2) * kstate
	r = 0.15
	div = y + x(2) - x(1) - c - r * bstate

	if (div < 0) then
		equity = theta0 + theta1 * abs(div) + 0.5 * theta2 * (div ** 2)
		div = div - equity 
	end if

	! next period
	knew = (1 - dpre) * kstate + x(1)
	bnew = (1 - r) * bstate + x(2)
	
	call drnnor(ndraws, edraw)
	! Monte-Carlo interation
	vsum = 0.	
	do i = 1, ndraws
		enew = rho * estate + sqrt(sigmu) * edraw(i,1)
		y = apha0 + apha1 * log(knew) + apha2 * ((log(knew)) ** 2) + apha3 * ((log(knew)) ** 3) + enew
		y = exp(y)
		v = y - (1+r) * bnew
		vsum = vsum + v
	end do
	vsum = vsum / ndraws

	f = div + drate * brate * vsum
	f = -1 * f
end subroutine objlast


subroutine objt(nv, x, f)
	
	use matforf90
	use mysub1
	use sharedata
	implicit none

	integer, intent(in) :: nv
	real(dp), intent(in), dimension(nv) :: x
	real(dp), intent(out) :: f

	integer, parameter:: iftemp1 = 25
	integer, parameter:: iftemp2 = 20
	integer, parameter:: iftemp3 = 21
	integer :: i, j, intk, intb, inte, sindex, group 
	real(dp) :: kstate, bstate, estate, y, c, equity, div, r, knew, bnew, enew, vsum, cutl, cutu, diff1, diff2
	real(dp), dimension(ndraws,1) :: edraw
	real(dp), dimension(kdim*bdim*edim, 1):: vvec
	real(dp), dimension(kdim*bdim*edim, 3):: state


	open(iftemp1,FILE='junk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
	open(iftemp2,FILE='vjunk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
	open(iftemp3,FILE='sjunk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 

	read(iftemp1,*) kstate, bstate, estate
	read(iftemp2,*) vvec
	read(iftemp3,*) state
	close(iftemp1)
	close(iftemp2)
	close(iftemp3)

	! current period	
	y = apha0 + apha1 * log(kstate) + apha2 * ((log(kstate)) ** 2) + apha3 * ((log(kstate)) ** 3) + estate
	y = exp(y)
	c = (lamd0 * x(1) + 0.5 * lamd1 * (x(1)/kstate)**2) * kstate
	r = 0.15
	div = y + x(2) - x(1) - c - r * bstate

	if (div < 0) then
		equity = theta0 + theta1 * abs(div) + 0.5 * theta2 * (div ** 2)
		div = div - equity 
	end if

	! next period
	knew = (1 - dpre) * kstate + x(1)
	bnew = (1 - r) * bstate + x(2)

	do i = 1, kdim
		if (i == 1) then
			if (knew <= state(i,1)) then
				intk = i
			end if
		else if (i == kdim) then
			if (knew > state(i,1)) then
				intk = i
			end if
		else
			cutl = state(i-1, 1)
			cutu = state(i, 1) 
			
			if (knew > cutl .AND. knew <= cutu) then
				diff1 = knew - cutl
				diff2 = cutu - knew
				if (diff1 > diff2) then
					intk = i
				else 
					intk = i-1
				end if

			end if
				
		end if

	end do

	group = 0
	do i = 1, bdim
		if (i == 1) then
			if (bnew <= state(i,2)) then
				intb = i
			end if
		else if (i == bdim) then
			if (bnew > state(kdim*bdim,2)) then
				intb = i
			end if
		else
			cutl = state(group-kdim+1, 2)
			cutu = state(group+1, 2) 
			
			if (bnew > cutl .AND. bnew <= cutu) then
				diff1 = bnew - cutl
				diff2 = cutu - bnew
				if (diff1 > diff2) then
					intb = i
				else 
					intb = i-1
				end if

			end if
				
		end if
		group = group + kdim
	end do

	call drnnor(ndraws, edraw)
	! Monte-Carlo interation
	vsum = 0.	
	do i = 1, ndraws
		enew = rho * estate + sqrt(sigmu) * edraw(i,1)

		group = 0
		do j = 1, edim
			if (j == 1) then
				if (enew <= state(j,3)) then
					inte = j
				end if
			else if (j == edim) then
				if (enew > state(kdim*bdim*edim,3)) then
					inte = j
				end if
			else
				cutl = state(group-kdim*bdim+1, 3)
				cutu = state(group+1, 3) 
			
				if (enew > cutl .AND. enew <= cutu) then
					diff1 = enew - cutl
					diff2 = cutu - enew
					if (diff1 > diff2) then
						inte = j
					else 
						inte = j-1
					end if
				end if
			end if
			group = group + kdim*bdim
		end do

		sindex = intk + (intb-1)*kdim + (inte-1)*kdim*bdim
		vsum = vsum + vvec(sindex,1)
	end do
	vsum = vsum / ndraws

	f = div + drate * brate * vsum

end subroutine objt


program main

	use matforf90
	use numerical_libraries
	use mysub1
	use sharedata
	implicit none

	integer, parameter :: ifstate = 23
	integer, parameter :: ifvalue = 24
	integer, parameter :: ifcurrent = 25
	integer, parameter :: dim = 2
	integer, parameter :: maxfcn = 30000
	real(dp), parameter :: tol = 0.001
	real(dp), parameter :: simp = 0.
				
	integer :: iloop, i	
	real(dp) :: fval
	real(dp), dimension(kdim*bdim*edim, 3):: rstate
	real(dp), dimension(kdim*bdim*edim, 1):: value
	real(dp), dimension(dim) :: inv, invguess
	real(dp), dimension(kdim*bdim*edim,6) :: results
	logical :: check
	external :: objlast, objt

	! Body starts here	
	check = .TRUE.
	open(ifstate,FILE='sjunk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 

	! discretize the state space randomly	
	call dstate(rstate)
	do i = 1, kdim*bdim*edim	
		write(ifstate,*) rstate(i,1), achar(9), rstate(i,2), achar(9), rstate(i,3)
	end do
	close(ifstate)

	
	if (check) then	
	! starting the backward induction
	invguess = 0.
	value = 0.
	results = 0.
	do iloop = 1, tperiods
		print *, iloop
		
		if (iloop == 1) then
			! the last period
			
			! loop over states
			do i = 1, kdim*bdim*edim
				! write current state onto disk
				open(ifcurrent,FILE='junk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
				write(ifcurrent,*) rstate(i,1), achar(9), rstate(i,2), achar(9), rstate(i,3)
				close(ifcurrent)	

				! solve optimal investment plans at current state
				call dumpol(objlast, dim, invguess, simp, tol, maxfcn, inv, fval)
				!call objlast(dim, invguess, fval)
				value(i,1) = fval
			end do	
			! Now we have got the invesment palns for all states.
		else

			! loop over states			
			do i = 1, kdim*bdim*edim
				! write current state onto disk
				open(ifcurrent,FILE='junk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
				write(ifcurrent,*) rstate(i,1), achar(9), rstate(i,2), achar(9), rstate(i,3)
				close(ifcurrent)					

				! solve optimal investment plans at current state
				call dumpol(objt, dim, invguess, simp, tol, maxfcn, inv, fval)
				value(i,1) = fval
			
				if (iloop == tperiods) then
					results(i, 1) = inv(1)
					results(i, 2) = inv(2)
					results(i, 3) = fval
					results(i, 4) = rstate(i,1)
					results(i, 5) = rstate(i,2)
					results(i, 6) = rstate(i,3)
				end if
			end do	

		end if

		! write "value" onto disk 
		open(ifvalue,FILE='vjunk.asc',STATUS='UNKNOWN',FORM='FORMATTED', RECL=80000) 
		do i = 1, kdim*bdim*edim
			write(ifvalue,*) value(i,1)	
		end do
		close(ifvalue)

	end do

	! interpolation here

	end if
end program main
