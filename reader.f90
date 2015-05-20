program reader
implicit none
real(8) :: dCV, dCA, signif


    open(unit=10,file="./grids/CMS_LO040_0217.dat",form='formatted',access='sequential')
    do while ( .not. eof(10) )
    

! -12.00000   -5.000       5.186608     2.589E+02

        read(unit=10,fmt="(1X,F9.5,3X,F6.3,5X,F8.6)") dCV, dCA, signif 
        print *, "xxx"

    enddo



end program