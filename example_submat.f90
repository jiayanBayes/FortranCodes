        use matforf90
        use linear_operator
        implicit none
 
        integer, parameter :: n=4
        integer :: i
 

                                         ! Integer vector (n)
        integer , dimension(:), allocatable ::   retain
 
                                         ! DP matrix (n,n)
        real (dp) , dimension(n,n)   ::   z
        real (dp) , dimension(:,:), allocatable ::   ztmp
 
 
        z = reshape( (/ 11,12,13,14,                    &
                        12,22,23,24,                    &
                        13,23,33,34,                    &
                        14,24,34,44 /)                  &
                       ,(/4,4/) )
 
 ! Start of program
                        ! take submatrix of z
 
        allocate(ztmp(n-1,n-1), retain(n-1))
        retain = seqa(1,1,n-1)
        retain(3) = 4
        print *, "retain vector"; call prnt(retain)
 
        ztmp = z(retain,retain)
        print *,"z matrix"; call prnt(z)
        print *,"ztmp";call prnt(ztmp)
 
                        ! in place substitution
 
        print *,"in place substitution  z(seqa(1,1,3),seqa(1,1,3)) = z(retain,retain) "
        z(seqa(1,1,3),seqa(1,1,3)) = z(retain,retain)
        print *,"modified z matrix"; call prnt(z)
 
 
 
end