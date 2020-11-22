
program ols

	use matforF90
	use linear_operator
    implicit none

	integer, parameter :: ifinput = 10
	integer, parameter :: n=2000, k=4
	character(80)  :: infile   = "dataset.asc"
	integer :: i
	integer userid

	real(dp), dimension(n)   :: y, resid
	real(dp), dimension(n,k) :: x
    real(dp), dimension(k,k) :: xx, xxtmp, xxi, varcov
	real(dp), dimension(k)   :: beta, xy, xytmp, std, tstud
	real(dp) :: sig2
	character(8), dimension(k) :: namevec 

! Code starts here

	namevec = (/ "x1  ","x2  ","x3  ","cons" /)


	open( ifinput,  FILE=infile,   STATUS='UNKNOWN')   
						! read all data into vectors and matrix
	do i=1,n
		read (ifinput,*), userid, y(i), x(i,:)
	end do

	xx = x .tx. x
	xy = x .tx. y
	xxi  = invpd(xx)
	beta = xxi .x. xy

	resid  = y - mprd(x,beta)
	sig2   = dot_product(resid,resid)/(n-k)
	varcov = sig2*xxi
	std    = sqrt(diag(varcov))
	tstud  = beta .div. std

		! print final results
	print *,""
    print *,"coef           beta      st.dev.   t-student"
	print *,"---------------------------------------------"
	do i=1,k
		print 100, namevec(i), beta(i), std(i), tstud(i)
	end do
	print *,"---------------------------------------------"

100 format(1x, a8,2f12.6, f12.2)

end