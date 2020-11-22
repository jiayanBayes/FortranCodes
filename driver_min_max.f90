!Driver selected objects

        use matforf90
        use linear_operator

        implicit none
        integer, parameter :: n=4

                                          !  INT vector (n)
        integer , dimension(n)   ::   index

                                          !  DP vector (n)
        real (dp) , dimension(n)   ::   x

                                          ! DP matrix (n,n)
        real (dp) , dimension(n,n)   ::   z

        z = reshape( (/ 11,21,31,41,                    &
                        12,22,32,42,                    &
                        13,23,33,43,                    &
                        14,24,34,44 /)                  &
                        ,(/4,4/) )

         print *,"z"; call prnt(z)

                        !MIN
                                ! output is automatically rank-1 vector
         x  = minc(z)
         print *, "minc" ; call prnt(x)

         index = minindc(z)
         print *,"index";  call prnt(index)

                                ! output is automatically rank-1 vector
         x  = maxc(z)
         print *, "maxc" ; call prnt(x)

         index = maxindc(z)
         print *,"index";  call prnt(index)



end