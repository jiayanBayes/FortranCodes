
Module shareddata
	
	use matforF90
	implicit none
	save

	real(dp), parameter :: r = 0.0
	real(dp), parameter :: b = 0.15
	real(dp), parameter :: phi = 0.995
	real(dp), parameter :: deta = 3
	real(dp), parameter :: a0 = 50
	real(dp), parameter :: apha = 0.3
	real(dp), parameter :: g0 = 100 
	real(dp), parameter :: gma = -0.01
	real(dp), parameter :: kini = 200

end module shareddata

subroutine functions (x, f, num)

	use shareddata
	use matforF90
	use numerical_libraries
	implicit none

	! Decalre inputs and outputs
	integer, intent(in) :: num
	real(dp), intent(in), dimension(num) :: x
	real(dp), intent(out), dimension(num) :: f

	! declare locals 
	integer :: period, i
	real(dp) :: z, rho1, rho2, ins, ratio, pro, dpro, dins
	real(dp), allocatable, dimension(:)  :: h, m

	period = num / 2
	allocate(h(period))
	allocate(m(period))
	do i = 1, period
		 if (i==1) then
			pro = a0 * kini ** apha
			z = pro - x(i) + (1 - b) * kini -x(period+i)
		 else
		    pro = a0 * x(i-1) ** apha
			z = pro - x(i) + (1-b) * x(i-1) - x(period+i)
		 end if
		 ins = g0 * x(period+i) ** gma
		 rho1 = (ins / z) ** deta
		 rho2 = (ins / z) ** (deta-1)

		 if (z <= ins) then
			rho1 = 1.
			rho2 = 1.
		 end if
			 	
         h(i) = 1 + r * rho1
	     m(i) = 1 + r * rho2
	end do

	do i = 1, period
		dpro = a0 * apha * x(i) ** (apha-1)
		dins = g0 * gma * x(period+i) ** (gma-1)
		if (i < period) then
			f(i) = dpro + (1-b) - h(i) / (phi * h(i+1))
		else
			f(i) = dpro - h(i)/phi
		end if   
		f(period+i) = (deta / (1 - deta)) * dins - (h(i) / m(i))
	end do

	deallocate(h)
	deallocate(m)	

end subroutine functions

program main
	use shareddata
	use matforF90
	use tools
	use numerical_libraries
	implicit none
	
	integer, parameter :: ifout = 30
	character(80), parameter :: outfile = "xout.asc"
	integer, parameter :: dim = 2
	integer, parameter :: imax = 1000000
	real(dp), parameter :: error = 0.00001
	real(dp) :: fnorm, z, ins, output, prob, profit, cinv
	real(dp), dimension(dim) :: xguess, xout, fout, xstar
	real(dp), dimension(dim/2) :: urv
	external :: functions
	integer :: i

	open(ifout,  FILE=outfile,  STATUS='UNKNOWN', RECL=80000)

	xguess(1:(dim/2)) = 30
	!xguess(1:1)= 10
	!xguess((dim/2-1):(dim/2)) =1
	xguess((1 + dim/2):dim) = 0.1
	 
	call dneqnf(functions, error, dim, imax, xguess, xout, fnorm)
	call prnt(xout)
	xstar = xout
	do i = 1, dim/2
		if (i == 1) then
			output = a0 * kini ** apha
			cinv = xout(i) - (1-b) * kini
			cinv = max(0., cinv)
			xout(i) = (1-b) * kini + cinv
			z = output - cinv - xout((dim/2)+i)
		else
			output = a0 * xout(i-1) ** apha
			cinv = xout(i) - (1-b) * xout(i-1)
			cinv = max(0., cinv)
			xout(i) = (1-b) * xout(i-1) + cinv
			z = output - cinv - xout((dim/2)+i)
		end if
		ins = g0 * xout((dim/2)+i) ** gma
		prob = (ins/z) ** deta
		prob = min(1., prob)
		if (z <= ins) then
			prob = 1.
		end if
		profit = z - (deta / (deta - 1)) * ins - (r/(deta-1)) * prob * z
		write(ifout,*) i, achar(9), profit, achar(9), cinv, achar(9), xout(i), achar(9),xstar(i), achar(9), xout((dim/2)+i), achar(9), prob, achar(9), &
						output, achar(9), ins, achar(9), z
	end do
end program main
