

! make Opt=Yes; ./TOPAZ Collider=6 TopDK=1 ObsSet=71 Correction=0 NLOParam=1 Process=91 VegasNc0=100000 VegasNc1=100000

MODULE ModAmplitudes_eeTTB
  use ModMisc
  implicit none
  
  
  

  CONTAINS

  
  
  
  SUBROUTINE Tree_ee_tbt(LO_Res_Pol)
  use ModProcess
  use ModParameters
  implicit none
    complex(8), intent(out) :: LO_Res_Pol
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4)
    complex(8) :: p1(4), p2(4), p3(4), p4(4)
    complex(8) :: lqcurrL(4), lqcurrR(4)   ! Light quark currents
    complex(8) :: hqcurrZL(4),hqcurrZR(4)  ! BSM Heavy quark current for Z
    complex(8) :: hqcurrGaL(4),hqcurrGaR(4)! BSM Heavy quark current for photon
    complex(8) :: propfact_Z,propfact_Pho
    logical,parameter :: removePhoton=.false.
    
    
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4)  = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4)  = ExtParticle(1)%Pol(1:4)
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)
    
   
    propfact_Z   = -cI/(2d0*sc_(p1,p2)-m_Z**2 + cI * Ga_ZExp * M_Z)
    if( removePhoton ) then
        propfact_Pho = (0d0,0d0)
        print *, "WARNING: no intermediate photon in Tree_ee_tbt"
    else
        propfact_Pho = -cI/(2d0*sc_(p1,p2))
    endif

    ! light quark current
    lqcurrL(:) = (-cI) * vbqq2(ub1,chir(.false.,v2)) 
    lqcurrR(:) = (-cI) * vbqq2(ub1,chir(.true.,v2))
    

    ! heavy quark current for SM and BSM Z coupling;  (-cI) is already included inside vbqV
    hqcurrZL(:)   = vbqV(ub3,lqcurrL,   &                                   
                         dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                         dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                         p1(1:4)+p2(1:4),2)

    hqcurrZR(:)   = vbqV(ub3,lqcurrR,   &
                         dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                         dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                         p1(1:4)+p2(1:4),2)

    

    ! heavy quark current for SM and BSM photon coupling;  (-cI) is already included inside vbqV
    hqcurrGaL(:)  = vbqV(ub3,lqcurrL,   &                                   
                         dcmplx(Q_top),dcmplx(Q_top),   &
                         dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                         p1(1:4)+p2(1:4),2)

    hqcurrGaR(:)  = vbqV(ub3,lqcurrR,   &
                         dcmplx(Q_top),dcmplx(Q_top),   &
                         dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                         p1(1:4)+p2(1:4),2)


    LO_Res_Pol = psp1_(hqcurrZL,v4) * couplZEE_left * propfact_Z   &
               + psp1_(hqcurrZR,v4) * couplZEE_right* propfact_Z   &
               + psp1_(hqcurrGaL,v4)* Q_el          * propfact_Pho &
               + psp1_(hqcurrGaR,v4)* Q_el          * propfact_Pho 



!    print *, "tree ",LO_Res_Pol,cdabs(LO_Res_Pol)
!    print *, "checker",(lqcurrR(:))*couplZEE_right
!    print *, "checker",(lqcurrL(:))*couplZEE_left
!    print *, "checker",propfact_Z
!    print *, "checker",couplZTT_left_dyn,couplZTT_right_dyn,couplZTT_left2_dyn,couplZTT_right2_dyn
!   print *, "check e-",-p2(:)
!   print *, "check e+",-p1(:)
!   print *, "check t ",p3(:)
!   print *, "check tb",p4(:)
!   print *, "check e-",v2(:)
!   print *, "check e+",ub1(:)
!   print *, "check t ",ub3(:)
!   print *, "check tb",v4(:)
!   print *, ""
!   pause
 

  RETURN
  END SUBROUTINE 


  
  
  


  


! my version of the virtual subroutine
SUBROUTINE VirtAmp_ee_Z_ttb(VirtAmpPol_ee_Z_ttb)
use ModParameters
! use ModKinematics
use ModProcess
use ModMisc
implicit none
real(8) :: MomExt(1:4,1:4),MomZ(1:4)
integer :: Heli(1:4),xe
integer, parameter :: Hplus=+1, Hminus=-1, Hdummy=0
complex(8) ::VirtAmpPol_ee_Z_ttb,VirtAmp2(Hminus:Hplus),TreeAmp,FF1,FF2,FF3,FF4
complex(8) :: Spi(1:4),BarSpi(1:4,Hminus:Hplus),ZPol1(1:4,Hminus:Hplus),ZPol2(1:4),ZProp,PProp,PPol1(1:4),PPol2(1:4)
real(8) :: sHat,IZV2,IZA2,IPV2,IPA2
logical,parameter :: removePhoton=.false.


    VirtAmpPol_ee_Z_ttb = (0d0,0d0)
    xe=0!   -2, -1, 0 

    MomExt(1:4,1) =-dble( ExtParticle(3)%Mom(1:4) )
    MomExt(1:4,2) =-dble( ExtParticle(4)%Mom(1:4) )
    MomExt(1:4,3) = dble( ExtParticle(2)%Mom(1:4) )
    MomExt(1:4,4) = dble( ExtParticle(1)%Mom(1:4) )
    MomZ(1:4) = MomExt(1:4,1)+MomExt(1:4,2)
    sHat = MomZ(1:4).dot.MomZ(1:4)


!   tree for initial state
    ZProp = -cI/dcmplx((MomZ(1:4).dot.MomZ(1:4))-M_Z**2, Ga_ZExp*M_Z)
    if( removePhoton ) then
        PProp = (0d0,0d0)
        print *, "WARNING: no intermediate photon in VirtAmp_ee_Z_ttb"
    else
        PProp = -cI/dcmplx((MomZ(1:4).dot.MomZ(1:4)), 0d0)
    endif
    
    
    Spi(1:4)    = ExtParticle(3)%Pol(1:4)
    BarSpi(1:4,Hdummy) = ExtParticle(4)%Pol(1:4)
    BarSpi(1:4,Hplus)  = Chir(.true.,BarSpi(1:4,Hdummy))
    BarSpi(1:4,Hminus) = Chir(.false.,BarSpi(1:4,Hdummy))
    ZPol1(1:4,Hplus)   = vbqq(4,BarSpi(1:4,Hplus),Spi(1:4))  * ZProp * couplZEE_left  *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq
    ZPol1(1:4,Hminus)  = vbqq(4,BarSpi(1:4,Hminus),Spi(1:4)) * ZProp * couplZEE_right *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq
    PPol1(1:4) = (ZPol1(1:4,Hplus)/couplZEE_left*(-Q_el)  + ZPol1(1:4,Hminus)/couplZEE_right*(-Q_el) )/ZProp*PProp




! --------------------------------------------


!   1-loop SM Z form factors for final state
    call CalcFormFactor2(xe,sHat,FF2,FF3,FF4)
    if(xe.eq.-2) FF2 = FF2 + 0d0! WFR  counter terms
    if(xe.eq.-1) FF2 = FF2 - 3d0! WFR  counter terms
    if(xe.eq. 0) FF2 = FF2 - 4d0 - 3d0*dlog(MuRen**2/m_top**2)! WFR  counter terms
! print *, "WARNING: switching off WFR CT"


!   FF2
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = (Chir(.true.,BarSpi(1:4,Hdummy)) *couplZTT_left    &
                       +  Chir(.false.,BarSpi(1:4,Hdummy))*couplZTT_right) * FF2
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    ZPol2(1:4) = vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq

!   FF3
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = (Chir(.true.,BarSpi(1:4,Hdummy)) *couplZTT_right    &
                       +  Chir(.false.,BarSpi(1:4,Hdummy))*couplZTT_left ) * FF3
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    ZPol2(1:4) = ZPol2(1:4) + vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq

!   FF4
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * (couplZTT_left+couplZTT_right) * FF4     
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    ZPol2(1:4) = ZPol2(1:4) + psp1_(BarSpi(1:4,Hdummy),Spi(1:4))*MomExt(1:4,4)   *(0d0,-1d0)!  multipy by -I because all the vbqq have this factor, too

    VirtAmp2(Hplus)  = ( (ZPol1(1:4,Hplus).dot.ZPol2(1:4))  ) 
    VirtAmp2(Hminus) = ( (ZPol1(1:4,Hminus).dot.ZPol2(1:4)) ) 
    VirtAmpPol_ee_Z_ttb = VirtAmp2(Hplus) + VirtAmp2(Hminus) 




!    print *, "virt ",VirtAmpPol_ee_Z_ttb, cdabs(VirtAmpPol_ee_Z_ttb)
! print *, "checker",(ZPol1(1:4,Hminus)) * (couplZEE_right)
! print *, "checker",(ZPol1(1:4,Hplus))  * (couplZEE_left)
! print *, "checker",ZProp
!    print *, "checker",couplZTT_left_dyn,couplZTT_right_dyn,couplZTT_left2_dyn,couplZTT_right2_dyn
!   print *, "check e-",momext(1:4,1)
!   print *, "check e+",momext(1:4,2)
!   print *, "check t ",momext(1:4,3)
!   print *, "check tb",momext(1:4,4)
!   print *, "check e-",ExtParticle(3)%Pol(1:4)
!   print *, "check e+",ExtParticle(4)%Pol(1:4)
!   print *, "check t ",ExtParticle(2)%Pol(1:4)
!   print *, "check tb",ExtParticle(1)%Pol(1:4)
! pause



! --------------------------------------------


!   1-loop SM Photon form factors for final state
    call CalcFormFactor2(xe,sHat,FF2,FF3,FF4)
    if(xe.eq.-2) FF2 = FF2 + 0d0! WFR  counter terms
    if(xe.eq.-1) FF2 = FF2 - 3d0! WFR  counter terms
    if(xe.eq. 0) FF2 = FF2 - 4d0 - 3d0*dlog(MuRen**2/m_top**2)! WFR  counter terms
! print *, "WARNING: switching off WFR CT"


!   FF2
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * (-Q_top) * FF2
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    PPol2(1:4) = vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq

!   FF3
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * (-Q_top) * FF3
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    PPol2(1:4) = PPol2(1:4) + vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq

!   FF4
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * (-Q_top) * FF4     
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    PPol2(1:4) = PPol2(1:4) + psp1_(BarSpi(1:4,Hdummy),Spi(1:4))*MomExt(1:4,4)   *(0d0,-1d0)! multipy by -I because all the vbqq have this factor, too

    VirtAmp2(Hplus)  = ( (PPol1(1:4).dot.PPol2(1:4))  ) 
    VirtAmpPol_ee_Z_ttb = VirtAmpPol_ee_Z_ttb + VirtAmp2(Hplus)




! --------------------------------------------


! make Opt=No; ./TOPAZ Collider=6 TopDK=0 ObsSet=70 Correction=1 NLOParam=2 Process=91 RelDelF1V=-1.0 RelDelF1A=-1.0 RelDelF2V=0.0 RelDelF2A=1.0
! todo: where is the log(mu) from the WFR ?? 
!       check CT against ttb+z 
! question: do SM and BSM piece interfere ???? 

!   1-loop BSM Z form factors for final state
    IZV2 = 1d0/2d0*( couplZTT_left2 + couplZTT_right2 )
    IZA2 = 1d0/2d0*( couplZTT_left2 - couplZTT_right2 )
    call CalcFormFactor2XX(xe,sHat,FF1,FF4)
    if(xe.eq.-1) then 
          FF1 = FF1 +(-3d0+1d0) *m_top!   WFR counter terms
          FF4 = FF4 +(-3d0+1d0)       !   WFR counter terms
    endif
    if(xe.eq. 0) then 
          FF1 = FF1 + (-4d0 - 3d0*dlog(MuRen**2/m_top**2) + 1d0*dlog(MuRen**2/TTBZ_MassScale**2)) *m_top!   WFR counter terms
          FF4 = FF4 + (-4d0 - 3d0*dlog(MuRen**2/m_top**2) + 1d0*dlog(MuRen**2/TTBZ_MassScale**2))       !   WFR counter terms
    endif
! print *, "WARNING: switching off WFR CT"


    FF1 = FF1 * 2d0! normalize such that FF1=1*mtop, FF4=1 gives the tree
    FF4 = FF4 * 2d0
!   FF1
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * FF1 * IZV2
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    ZPol2(1:4) = vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq


!   FF4
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = ( +Chir(.true.,BarSpi(1:4,Hdummy)) + Chir(.false.,BarSpi(1:4,Hdummy))  ) * IZV2 * FF4   &
                       + ( +Chir(.true.,BarSpi(1:4,Hdummy)) - Chir(.false.,BarSpi(1:4,Hdummy))  ) * (0d0,-1d0)*IZA2 * FF4 
    ZPol2(1:4) = ZPol2(1:4) + psp1_(BarSpi(1:4,Hdummy),Spi(1:4))*MomExt(1:4,4)   *(0d0,-1d0)

    VirtAmp2(Hplus)  = ( (ZPol1(1:4,Hplus).dot.ZPol2(1:4))  )     *(-1d0)! phase convention to match SM part
    VirtAmp2(Hminus) = ( (ZPol1(1:4,Hminus).dot.ZPol2(1:4)) )     *(-1d0)! phase convention to match SM part 
    VirtAmpPol_ee_Z_ttb  = VirtAmpPol_ee_Z_ttb + VirtAmp2(Hplus) + VirtAmp2(Hminus) 
 

!    print *, "virt ",VirtAmpPol_ee_Z_ttb, cdabs(VirtAmpPol_ee_Z_ttb)
! print *, "checker",(ZPol1(1:4,Hminus)) * (couplZEE_right)
! print *, "checker",(ZPol1(1:4,Hplus))  * (couplZEE_left)
! print *, "checker",ZProp
!    print *, "checker",couplZTT_left_dyn,couplZTT_right_dyn,couplZTT_left2_dyn,couplZTT_right2_dyn
!   print *, "check e-",momext(1:4,1)
!   print *, "check e+",momext(1:4,2)
!   print *, "check t ",momext(1:4,3)
!   print *, "check tb",momext(1:4,4)
!   print *, "check e-",ExtParticle(3)%Pol(1:4)
!   print *, "check e+",ExtParticle(4)%Pol(1:4)
!   print *, "check t ",ExtParticle(2)%Pol(1:4)
!   print *, "check tb",ExtParticle(1)%Pol(1:4)
! pause



! --------------------------------------------

!   1-loop BSM photon form factors for final state
    IPV2 = 1d0/2d0*( couplGaTT_left2 + couplGaTT_right2 )
    IPA2 = 1d0/2d0*( couplGaTT_left2 - couplGaTT_right2 )
    call CalcFormFactor2XX(xe,sHat,FF1,FF4)
    if(xe.eq.-1) then 
          FF1 = FF1 +(-3d0+1d0) *m_top!   WFR counter terms
          FF4 = FF4 +(-3d0+1d0)       !   WFR counter terms
    endif
    if(xe.eq. 0) then 
          FF1 = FF1 + (-4d0 - 3d0*dlog(MuRen**2/m_top**2) + 1d0*dlog(MuRen**2/TTBZ_MassScale**2)) *m_top!   WFR counter terms
          FF4 = FF4 + (-4d0 - 3d0*dlog(MuRen**2/m_top**2) + 1d0*dlog(MuRen**2/TTBZ_MassScale**2))       !   WFR counter terms
    endif
! print *, "WARNING: switching off WFR CT"


    FF1 = FF1 * 2d0
    FF4 = FF4 * 2d0
!   FF1
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = BarSpi(1:4,Hdummy) * FF1 * IPV2
    Spi(1:4) = ExtParticle(1)%Pol(1:4)
    PPol2(1:4) = vbqq(4,BarSpi(1:4,Hdummy),Spi(1:4)) *dsqrt(2d0)!  compensates 1/sqrt(2) from vbqq


!   FF4
    BarSpi(1:4,Hdummy) = ExtParticle(2)%Pol(1:4)
    BarSpi(1:4,Hdummy) = ( +Chir(.true.,BarSpi(1:4,Hdummy)) + Chir(.false.,BarSpi(1:4,Hdummy))  ) * IPV2 * FF4   &
                       + ( +Chir(.true.,BarSpi(1:4,Hdummy)) - Chir(.false.,BarSpi(1:4,Hdummy))  ) * (0d0,-1d0)*IPA2 * FF4 
    PPol2(1:4) = PPol2(1:4) + psp1_(BarSpi(1:4,Hdummy),Spi(1:4))*MomExt(1:4,4)   *(0d0,-1d0)

    VirtAmp2(Hplus)  = ( (PPol1(1:4).dot.PPol2(1:4))  )      *(-1d0)! phase convention to match SM part
    VirtAmpPol_ee_Z_ttb = VirtAmpPol_ee_Z_ttb + VirtAmp2(Hplus)


RETURN
END SUBROUTINE


  




SUBROUTINE CalcFormFactor2(xe,sHat,FF2,FF3,FF4)
use ModParameters
use ModMisc
implicit none
complex(8) :: FF2,FF3,FF4
real(8) :: xs,sHat,beta,qp,qm,DDILOG
integer :: xe
complex(8) :: SI1finite,SI2finite,SI3finite,SI4finite
complex(8) :: qlI1,qlI2,qlI3,SI(1:4),C(1:1)
logical,parameter :: onlyLO=.false.

    if( onlyLO ) then! disable WFR CT for this
        print *, "WARNING: returning only tree in CalcFormFactor2"
        FF2 = 1d0
        FF3 = 0d0
        FF4 = 0d0
    RETURN
    endif


    C(1) = (2d0,0d0)
    if( xe.ne.0 ) C(:)=(0d0,0d0)

    SI(1) = qlI1(m_Top**2,MuRen**2,xe)
    SI(2) = qlI2(m_Top**2,0d0,m_Top**2,MuRen**2,xe)
    SI(3) = qlI2(sHat,m_Top**2,m_Top**2,MuRen**2,xe)
    SI(4) = qlI3(shat,m_top**2,m_Top**2,m_Top**2,m_Top**2,0d0,MuRen**2,xe)
    
    FF2 = ( 16d0*m_Top**4*SI(4) + sHat*(C(1) - 4d0*SI(2) + 3d0*SI(3) + 2d0*sHat*SI(4)) -4d0*m_Top**2*(C(1) - 3d0*SI(2) + 2d0*SI(3) + 3d0*sHat*SI(4)) )/(4d0*m_Top**2 - sHat)
    FF3 = 4d0*m_Top**2*(SI(2) - SI(3))/(4d0*m_Top**2 - sHat)
    FF4 = -2d0*(SI(1) + m_Top**2*(C(1)/2d0 - 2d0*SI(2) + SI(3)))/(4d0*m_Top**3 - m_Top*sHat)  !*Dot(ep3,p5)

END SUBROUTINE


  
  

SUBROUTINE CalcFormFactor2XX(xe,sHat,FF1,FF4)
use ModParameters
use ModMisc
implicit none
complex(8) :: FF1,FF4
real(8) :: xs,sHat
integer :: xe
complex(8) :: qlI1,qlI2,qlI3,SI(1:3)
logical,parameter :: onlyLO=.false.

    if( onlyLO ) then! disable WFR CT for this
        print *, "WARNING: returning only tree in CalcFormFactor2XX"
        FF1 = m_Top
        FF4 = 1d0         
        RETURN
    endif

    SI(1) = qlI2(m_Top**2,0d0,m_Top**2,MuRen**2,xe)
    SI(2) = qlI2(sHat,m_Top**2,m_Top**2,MuRen**2,xe)
    SI(3) = qlI3(shat,m_top**2,m_Top**2,m_Top**2,m_Top**2,0d0,MuRen**2,xe)

    FF1 = 2d0*m_Top*(SI(1) - SI(2) + (2d0*m_Top**2 - sHat)*SI(3))
    FF4 = 2d0*(2d0*m_Top**2 - sHat)*(2d0*SI(1) - 2d0*SI(2) + (4d0*m_Top**2 - sHat)*SI(3))/(4d0*m_Top**2 - sHat)   ! * Dot(ep3,p5)


END SUBROUTINE


  
  


  




  SUBROUTINE Tree_ee_tbtg_f(LO_Res_Pol) ! Emission from the final state
  use ModProcess
  use ModParameters
  use ModMisc
  implicit none
    complex(8), intent(out) :: LO_Res_Pol
    complex(8) :: ub1(4), v2(4), ub3(4), v4(4), e5(4), sptmp(4)
    complex(8) :: p1(4), p2(4), p3(4), p4(4), p5(4), p35(4), p45(4)
    complex(8) :: lqcurrL(4), lqcurrR(4) ! Light quark currents
    complex(8) :: hqcurrZL(1:2,1:4),hqcurrZR(1:2,1:4)  ! BSM Heavy quark current for Z
    complex(8) :: hqcurrGaL(1:2,1:4),hqcurrGaR(1:2,1:4)! BSM Heavy quark current for photon
    complex(8) :: propfact_Z,propfact_Pho
    logical,parameter :: removePhoton=.false.

   
    ub1(1:4) = ExtParticle(4)%Pol(1:4)
    v2(1:4) = ExtParticle(3)%Pol(1:4)
    ub3(1:4) = ExtParticle(2)%Pol(1:4)
    v4(1:4) = ExtParticle(1)%Pol(1:4)
    e5(1:4) = ExtParticle(5)%Pol(1:4)
!   e5(1:4) = ExtParticle(5)%Mom(1:4); print *, "checking gauge invariance"
    
    p1(1:4) = ExtParticle(4)%Mom(1:4)
    p2(1:4) = ExtParticle(3)%Mom(1:4)
    p3(1:4) = ExtParticle(2)%Mom(1:4)
    p4(1:4) = ExtParticle(1)%Mom(1:4)
    p5(1:4) = ExtParticle(5)%Mom(1:4)

    p35 = p3 + p5
    p45 = p4 + p5

    propfact_Z   = -cI/(2d0*sc_(p1,p2)-m_Z**2 + cI * Ga_ZExp * M_Z)
    if( removePhoton ) then
        propfact_Pho = (0d0,0d0)
        print *, "WARNING: no intermediate photon in Tree_ee_tbt"
    else
        propfact_Pho = -cI/(2d0*sc_(p1,p2))
    endif    



    ! light quark current
    lqcurrL(:) = (-cI) * vbqq2(ub1,chir(.false.,v2)) 
    lqcurrR(:) = (-cI) * vbqq2(ub1,chir(.true.,v2))
    


    ! heavy quark current for SM and BSM Z coupling;  (-cI) is already included inside vbqV
    ! First diagram: emission from quark line
    ! Phase factor: -cI (Vertex) * cI (Propagator) = 1
    sptmp(:) = spb2_(ub3,e5)
    sptmp(:) = (spb2_(sptmp,p35)+m_Top * sptmp)/(sc_(p35,p35)-m_Top**2)
    hqcurrZL(1,:) = vbqV(sptmp,lqcurrL,   &                                   
                         dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                         dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                         p1(1:4)+p2(1:4),2)
    hqcurrZR(1,:) = vbqV(sptmp,lqcurrR,   &
                         dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                         dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                         p1(1:4)+p2(1:4),2)
    ! Second diagram: emission from anti-quark line
    sptmp(:) = spi2_(e5,v4)
    sptmp(:) = (-spi2_(p45,sptmp)+m_Top * sptmp)/(sc_(p45,p45)-m_Top**2)    
    hqcurrZL(2,:) = vVq(lqcurrL,sptmp,   &                                   
                        dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                        dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                        p1(1:4)+p2(1:4),2)
    hqcurrZR(2,:) = vVq(lqcurrR,sptmp,   &
                        dcmplx(couplZTT_left),dcmplx(couplZTT_right),   &
                        dcmplx(couplZTT_left2),dcmplx(couplZTT_right2), &
                        p1(1:4)+p2(1:4),2)






    ! heavy quark current for SM and BSM photon coupling;  (-cI) is already included inside vbqV
    ! First diagram: emission from quark line
    ! Phase factor: -cI (Vertex) * cI (Propagator) = 1
    sptmp(:) = spb2_(ub3,e5)
    sptmp(:) = (spb2_(sptmp,p35)+m_Top * sptmp)/(sc_(p35,p35)-m_Top**2)
    hqcurrGaL(1,:) = vbqV(sptmp,lqcurrL,   &                                   
                          dcmplx(Q_top),dcmplx(Q_top),   &
                          dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                          p1(1:4)+p2(1:4),2)
    hqcurrGaR(1,:) = vbqV(sptmp,lqcurrR,   &
                          dcmplx(Q_top),dcmplx(Q_top),   &
                          dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                          p1(1:4)+p2(1:4),2)
    ! Second diagram: emission from anti-quark line
    sptmp(:) = spi2_(e5,v4)
    sptmp(:) = (-spi2_(p45,sptmp)+m_Top * sptmp)/(sc_(p45,p45)-m_Top**2)    
    hqcurrGaL(2,:) = vVq(lqcurrL,sptmp,   &                                   
                         dcmplx(Q_top),dcmplx(Q_top),   &
                         dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                         p1(1:4)+p2(1:4),2)
    hqcurrGaR(2,:) = vVq(lqcurrR,sptmp,   &
                         dcmplx(Q_top),dcmplx(Q_top),   &
                         dcmplx(couplGaTT_left2),dcmplx(couplGaTT_right2), &
                         p1(1:4)+p2(1:4),2)





    LO_Res_Pol =   &
!                diagram 1
               + psp1_(hqcurrZL(1,:),v4) * couplZEE_left * propfact_Z   &
               + psp1_(hqcurrZR(1,:),v4) * couplZEE_right* propfact_Z   &
               + psp1_(hqcurrGaL(1,:),v4)* Q_el          * propfact_Pho &
               + psp1_(hqcurrGaR(1,:),v4)* Q_el          * propfact_Pho &
!                diagram 2
               + psp1_(ub3,hqcurrZL(2,:)) * couplZEE_left * propfact_Z   &
               + psp1_(ub3,hqcurrZR(2,:)) * couplZEE_right* propfact_Z   &
               + psp1_(ub3,hqcurrGaL(2,:))* Q_el          * propfact_Pho &
               + psp1_(ub3,hqcurrGaR(2,:))* Q_el          * propfact_Pho 


  RETURN
  END SUBROUTINE 



  

END MODULE 
