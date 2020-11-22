!Example concatanation

        use matforf90
        use linear_operator

        implicit none

        integer, parameter :: n=4

                                          ! Integer vector (n,1)
        integer , dimension(n,1) ::   iy
        integer , dimension(:,:), allocatable   ::   iy1

                                          ! Integer matrix (n,n)
        integer , dimension(n,n) ::   iz
        integer , dimension(:,:), allocatable   ::   iz1
        integer , dimension(:,:), allocatable   ::   iz2


                                                        ! Character vectors
        character(nmaxchar), dimension(n) :: namevars
        character(nmaxchar), dimension(1) :: namevars1
        character(nmaxchar), dimension(:), allocatable :: namevars2


! Initialisation

     iy = reshape( (/ 1,2,3,4 /) , (/4,1 /) )

         iz = reshape( (/ 11,12,13,14,                  &
                          21,22,23,24,                  &
                          31,32,33,34,                  &
                          41,42,43,44 /)                &
                        ,(/4,4/) )



        namevars  = (/"Hello", "How  ", "are  ", "you. "/)
        namevars1 = (/"Very well."/)

! Start of code

                                ! Integer

        print *, "iz"; call prnt(iz)
        print *, "iy"; call prnt(iy)

                            ! Concatanation ~ = horizontal  | = vertical
                            ! result of iz~iy
        allocate(iz1(n,n+1))
        iz1 = hcon(iz,iy)
        print *,"iz~iy "; call prnt(iz1)
        iz1 = 0
        iz1 = iz .hc. iy
        print *,"iz .hc. iy "; call prnt(iz1)
        deallocate(iz1)

                                        ! result of iz|iy

        allocate(iy1(1,n))
        iy1 = .t. iy                    ! transposed of iy

        allocate(iz1(n+1,n))
        iz1 = vcon(iz,iy1)
        print *,"iz|iy' "; call prnt(iz1)
        iz1 = 0
        iz1 = iz .vc. iy1
        print *,"iz .vc. iy' "; call prnt(iz1)
        deallocate(iz1)

                                        ! four blocks

        allocate(iz1(n+1,n+1))
        allocate(iz2(1,1))              ! rank-2 scalar

        iz2 = 5

        iz1  =      ( iz    .hc. iy    )     &
               .vc. ( iy1   .hc. iz2)

        print *,"iz1 "; call prnt(iz1)
        deallocate(iz1)

                                        ! character vectors

        allocate(namevars2(n+1))
        namevars2 = namevars .vc. namevars1
        call prnt(namevars2)



end