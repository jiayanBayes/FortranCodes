
    use matforf90
    implicit none

        integer, parameter :: n=4
                               ! DP rank-1 vector
        real (dp) , dimension(n)     ::   x

                               ! string
        character(80)   :: astring
        character(4)    :: a
        integer :: i
! Initialisation

        x = reshape( (/ 1,2,3,4 /) , (/4/) )

        astring  = "23  .  33.23 "
        a = "1432"
        i = 1234

! Start of program

        print *, "rows of x", rows(x), "   cols of x", cols(x)
        x = 0
        call parse_numbers(astring, x)
        print *,"parsed string"; call prnt(x)
        print *, "string to integer    ", cvtsi(a)
        print *, "integer to string    ", cvtis(i)


end

