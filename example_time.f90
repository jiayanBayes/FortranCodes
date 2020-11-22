
    use matforf90

    implicit none

        integer, parameter :: n=4
		integer    :: stime
		integer    :: istart, iend , iesttime
		integer    :: i,j
        real(dp)   :: x 
		                  
        character(40)   :: timestring
		character(24)   :: greeting

 
        call fdate(greeting)             
        print 100, greeting
100     format(34x,'Welcome : ' , a24 /)
 

		istart = stime()
		x=0
		do i=1,10000
			do j=1,100000
			x = x*x
			end do
		end do
		iend = stime()
		iesttime = iend-istart
		call etstr(iesttime,timestring)
		print *, "elapsed time : ", timestring

end

