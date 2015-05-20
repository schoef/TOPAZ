module ModIntDipoles_eeTTB
use ModIntDipoles
implicit none

integer, parameter, private :: dp = selected_real_kind(15)

double precision, private :: MomDK(0:3,1:6)

 contains

 
  SUBROUTINE IntDip_eettb(p,IDip)
    use ModParameters
    use ModMisc
    implicit none
    real(8) :: p(4,4)
    real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sHat
    real(8) :: dipsoft,dipfini,dipplus,epcorr
    real(8) :: CF 
    integer :: n,emi


    CF = 4d0/3d0


    z = 0.2d0 ! dummy

    IDip(1:3) = 0d0
    do n=1,2
       if(n.eq.1) then
          dipsoft =ff_qq(m_Top,m_Top,p,3,4,z,1) * CF
          dipfini =ff_qq(m_Top,m_Top,p,3,4,z,2) * CF!  is zero
          dipplus =ff_qq(m_Top,m_Top,p,3,4,z,3) * CF!  is zero
          emi = 3
       elseif(n.eq.2) then
          dipsoft =ff_qq(m_Top,m_Top,p,4,3,z,1) * CF
          dipfini =ff_qq(m_Top,m_Top,p,4,3,z,2) * CF!  is zero
          dipplus =ff_qq(m_Top,m_Top,p,4,3,z,3) * CF!  is zero
          emi = 3
       endif

       ! this is for check against virtual amplitude
!       dipplus = 0d0; print *, 'dipoles check'
!       dipfini = 0d0; print *, 'dipoles check'

       if(emi.eq.1) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)
       elseif(emi.eq.2) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(3) = IDip(3) + (dipfini+dipplus)
       elseif(emi.eq.3) then
          IDip(1) = IDip(1) + (dipsoft-dipplus)
          IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0!  is zero
          IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0!  is zero
       endif
    enddo

    
     !print *, epinv
     !print *, "IntDip",IDip(2:3)

    ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
!     epcorr=epinv


!     APsoft= 3d0/2d0*CF     * epcorr
!     APfini= (-1d0-z)*CF    * epcorr
!     APplus= 2d0*CF/(1d0-z) * epcorr

! this is for check against virtual amplitude
!    APplus = 0d0; print *, 'dipoles check'
!    APfini = 0d0; print *, 'dipoles check'

!     IDip(1) = IDip(1) + (APsoft - APplus)*2d0
!     IDip(2) = IDip(2) + (APfini + APplus)
!     IDip(3) = IDip(3) + (APfini + APplus)

     !print *, "AP",(APsoft - APplus),(APfini + APplus)
     !print *, "sum",IDip(2:3)
     !pause

  END SUBROUTINE 
  

END MODULE
