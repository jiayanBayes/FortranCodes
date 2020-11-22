!Example mdotprod mdotdiv

        use matforf90
        use linear_operator

        implicit none
        integer, parameter :: n=4

                                       ! DP rank-1 vector
        real (dp) , dimension(n)     ::   x
        real (dp) , dimension(:), allocatable     ::   x1

                                       ! DP rank-2 matrix (n,n
        real (dp) , dimension(n,1)   ::   y
        real (dp) , dimension(:,:), allocatable   ::   y1
        real (dp) , dimension(n,n)   ::   z
        real (dp) , dimension(:,:), allocatable   ::   z1


! Initialisation

         x = reshape( (/ 1,2,3,4 /) , (/4/) )
         y = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

         z = reshape( (/ 11,12,13,14,                   &
                         21,22,23,24,                    &
                         31,32,33,34,                    &
                         41,42,43,44 /)                  &
                        ,(/4,4/) )

! Start of program

        print *, "Double precision version"
        print *, "x"; call prnt(x)
        print *, "z"; call prnt(z)

        allocate(z1(n,n), x1(n), y1(n,1))

                                !MDOTPROD  dp

        print *,"MDOTPROD"

        z1 = mdotprod(z,z)
        print *,"z.*z"; call prnt(z1)

        z1 = mdotprod(z,y)
        print *,"z.*y"; call prnt(z1)

        z1 = mdotprod(z, .t. y)
        print *,"z.*y'"; call prnt(z1)


                                ! with rank 1
        print *,"rank-1 x"
        z1 = mdotprod(z,x)
        print *,"z.*x"; call prnt(z1)

                !linear operator
        print *, "using linear operator .prod."
        z1 = z .prod. x
        print *,"z.prod.x"; call prnt(z1)


        print *,"MDOTDIV"        !MDOTDIV dp

        z1 = mdotdiv(z,z)
        print *,"z./z"; call prnt(z1)

                !linear operator
        print *, "using linear operator .div."
        z1 = z .div. x
        print *,"z.div.x"; call prnt(z1)

                                ! In place calculation

        z = mdotdiv(z,x)
        print *,"z./x in place"; call prnt(z)

        z = mdotdiv(z,z)
        print *,"z./z in place"; call prnt(z)



end