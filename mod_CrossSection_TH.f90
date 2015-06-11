MODULE ModCrossSection_TH
  use ModTopDecay
  implicit none
  
  integer,private,parameter :: NumMaxHisto=45


contains




  FUNCTION EvalCS_LO_tdubH_MPI(yRnd,VgsWgt,res)
    implicit none
    integer :: EvalCS_LO_tdubH_MPI
    real(8) ::  yRnd(*),res(*),VgsWgt

    res(1) = EvalCS_LO_tdubH(yRnd,VgsWgt)
    EvalCS_LO_tdubH_MPI=0
    RETURN
  END FUNCTION EvalCS_LO_tdubH_MPI





  
  FUNCTION EvalCS_LO_tdubH(yRnd,VgsWgt)
! RR May 29 2015
! Routine for production of H(p3)+t(p4)+jet(p5)
! Looks very different from usual TOPAZ routines because LO amplitudes computed analytically
! Many formulae taken from implementation in MCFM, see hep-ph:/1302.3856 and hep-ph:/1204.1513
    use ModProcess
    use ModParameters
    use ModKinematics
    use ModMisc
    use ModSingleTopHAmps
    implicit none
    real(8) :: EvalCS_LO_tdubH,yRnd(1:VegasMxDim),VgsWgt
    real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles,PSWgt,PSWgt2,PSWgt3,pdf(-6:6,1:2)
    real(8) ::PObs(1:NumMaxHisto) 
    real(8) :: MomExt(1:4,1:11),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol,s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
    complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),CoupFac,decay_amp(1:2)
    real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
    integer :: NBin(1:NumMaxHisto),NHisto,j
    logical :: applyPSCut
    
    EvalCS_LO_tdubH = 0d0

    call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    call SmearExternal(yRnd(18),M_H,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)                                                                                                      
    MH_Inv = M_H
   if( EHat.le.m_Top+MH_Inv ) then
      EvalCS_LO_tdubH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ISFac = MomCrossing(MomExt)
   call setPDFs(eta1,eta2,MuFac,pdf)


!!!!! MCFM phase space point to check
 !! undecayed tops
!   MomExt(1:4,1)=  (/12.804565518142741d0,   0.0000000000000000d0,  0.0000000000000000d0, 12.804565518142741d0 /)     
!   MomExt(1:4,2)=  (/6853.0396776786965d0,   0.0000000000000000d0,  0.0000000000000000d0, -6853.0396776786965d0 /)     
!   MomExt(1:4,3)=  (/ 2151.2050020915885d0,   59.848489328696687d0, -121.56893937712810d0, -2143.2326314069742d0 /)     
!   MomExt(1:4,4)=  (/ 4272.0664685769034d0,   10.381035193439132d0,  36.524127576457872d0, -4268.3932731359500d0 /)     
!   MomExt(1:4,5)=  (/ 442.57277252834774d0,  -70.229524522135819d0,  85.044811800670232d0, -428.60920761762918d0 /)
!   MomExt=MomExt*GeV
!!!!!end MCFM

   CoupFac=2d0*g_weak**2/Vev*ci
   ColFac=9d0
   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
      PSWgt = PSWgt * PSWgt2
! usual decay top(p4) --> b(p6) + e+(p7) + nu(p8)
!      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
   ENDIF


   IF( HDECAYS.NE.0 .AND. HDECAYS.NE.-2 ) THEN
      call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,9:10),PSWgt3)
      PSWgt = PSWgt * PSWgt3   !* MHJacobi                                                                                                 
   ENDIF

! MCFM point -- decayed tops
!   MomExt(1:4,1)= (/  221.80711450036728d0,  0.0000000000000000d0,   0.0000000000000000d0,  221.80711450036728d0/)
!   MomExt(1:4,2)= (/  1566.2816147059500d0,  0.0000000000000000d0,   0.0000000000000000d0,   -1566.2816147059500d0/)
!   MomExt(1:4,3)= (/   718.44658318060931d0, -384.77702987430052d0,  -164.14263304640136d0,  -570.34491860767741d0/)
!   MomExt(1:4,4)= (/   359.92426242009503d0,  47.101295335937778d0,  -181.18432631335395d0,  -254.10663618784588d0/)
!   MomExt(1:4,5)= (/   709.71788360561300d0,  337.67573453836275d0,   345.32695935975534d0,  -520.02294541005938d0/)
!   MomExt(1:4,6)= (/   72.712116618267146d0, -38.834935190585639d0,  -3.7286676369820384d0,  -61.359569339300407d0/)
!   MomExt(1:4,7)= (/   250.71614048335576d0,  99.022624769419508d0,  -147.13355105314275d0,  -177.21405428784917d0/)
!   MomExt(1:4,8)= (/   36.496005318472129d0, -13.086394242896088d0,  -30.322107623229172d0,  -15.533012560696307d0/)
!   MomExt=MomExt*GeV
! end MCFM

   call Kinematics_TH(0,MomExt(1:4,1:11),(/4,5,3,1,2,0,6,7,8,9,10/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      EvalCS_LO_tdubH = 0d0
      return
   endif

  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt



! setup momenta for spinor helicity products -- undecayed tops
      p4Dp5=MomExt(1,5)*MomExt(1,4)-MomExt(2,5)*MomExt(2,4)-MomExt(3,5)*MomExt(3,4)-MomExt(4,5)*MomExt(4,4)
      p4Dp7=MomExt(1,7)*MomExt(1,4)-MomExt(2,7)*MomExt(2,4)-MomExt(3,7)*MomExt(3,4)-MomExt(4,7)*MomExt(4,4)
      p2Dp3=MomExt(1,2)*MomExt(1,3)-MomExt(2,2)*MomExt(2,3)-MomExt(3,2)*MomExt(3,3)-MomExt(4,2)*MomExt(4,3)
      MomExtFlat(1,1:4)=MomExt(1:4,1)
      MomExtFlat(2,1:4)=MomExt(1:4,2)
      MomExtFlat(3,1:4)=m_H**2/2d0/p2Dp3*MomExt(1:4,2)
      MomExtFlat(4,1:4)=MomExt(1:4,3)-MomExtFlat(3,1:4)
      MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,5)/2d0/p4Dp5
      MomExtFlat(6,1:4)=MomExt(1:4,4)-MomExtFlat(5,1:4)
      MomExtFlat(7,1:4)=MomExt(1:4,5)

    ! use different flattened momenta for top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,7)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,4)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,6)
         MomExtFlatDK(9,1:4)=MomExt(1:4,7)
         MomExtFlatDK(10,1:4)=MomExt(1:4,8)
      ENDIF
      

! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

   IF (CORRECTION .EQ. 0) THEN
      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call tdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
         call ubhtdamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Up_,Bot_,1:2))
         call ubhtdamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Bot_,Up_,1:2))        
         call ubhtdamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(ADn_,Bot_,1:2))
         call ubhtdamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(Bot_,ADn_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp=LOAmp*CoupFac


!      IF (TOPDECAYS .EQ. 0) THEN
!         print *, "u_b",(abs(LOAmp(Up_,Bot_,1))**2+abs(LOAmp(Up_,Bot_,2))**2) * ColFac * ISFac * GeV**2
!         print *, "b_u",(abs(LOAmp(Bot_,Up_,1))**2+abs(LOAmp(Bot_,Up_,2))**2) * ColFac * ISFac * GeV**2
!         print *, "ad_b",(abs(LOAmp(ADn_,Bot_,1))**2+abs(LOAmp(ADn_,Bot_,2))**2) * ColFac * ISFac * GeV**2
!         print *, "b_adn",(abs(LOAmp(Bot_,ADn_,1))**2+abs(LOAmp(Bot_,ADn_,2))**2) * ColFac * ISFac * GeV**2
!         pause
!      ELSE
!         print *, "u_b",(abs(LOAmp(Up_,Bot_,1))**2+abs(LOAmp(Up_,Bot_,2))**2) * ColFac * ISFac * GeV**6
!         print *, "b_u",(abs(LOAmp(Bot_,Up_,1))**2+abs(LOAmp(Bot_,Up_,2))**2) * ColFac * ISFac * GeV**6
!         print *, "ad_b",(abs(LOAmp(ADn_,Bot_,1))**2+abs(LOAmp(ADn_,Bot_,2))**2) * ColFac * ISFac * GeV**6
!         print *, "b_adn",(abs(LOAmp(Bot_,ADn_,1))**2+abs(LOAmp(Bot_,ADn_,2))**2) * ColFac * ISFac * GeV**6
!         pause
!      ENDIF
      

      
      LO_Res_UnPol= &
           + (abs(LOAmp(Up_,Bot_,1))**2+abs(LOAmp(Up_,Bot_,2))**2)   * (pdf(Up_,1)*pdf(Bot_,2)  + pdf(Chm_,1)*pdf(Bot_,2)) &
           + (abs(LOAmp(Bot_,Up_,1))**2+abs(LOAmp(Bot_,Up_,2))**2)   * (pdf(Bot_,1)*pdf(Up_,2)  + pdf(Bot_,1)*pdf(Chm_,2)) &
           + (abs(LOAmp(ADn_,Bot_,1))**2+abs(LOAmp(ADn_,Bot_,2))**2) * (pdf(ADn_,1)*pdf(Bot_,2) + pdf(AStr_,1)*pdf(Bot_,2)) &
           + (abs(LOAmp(Bot_,ADn_,1))**2+abs(LOAmp(Bot_,ADn_,2))**2) * (pdf(Bot_,1)*pdf(ADn_,2) + pdf(Bot_,1)*pdf(AStr_,2))

      LO_Res_UnPol = LO_Res_UnPol * ColFac * ISFac * WidthExpansion
      EvalCS_LO_tdubH = LO_Res_Unpol * PreFac

   ENDIF

   if( IsNan(EvalCS_LO_tdubH) ) then
        print *, "NAN:",EvalCS_LO_tdubH
        print *, yRnd(:)
        print *, LO_Res_UnPol

        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_LO_tdubH = 0d0
!         pause                                                                                                                                 
        return
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_LO_tdubH,BinValue=PObs(NHisto))                                                                            
!      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggH)
   enddo
   EvalCounter = EvalCounter + 1


   EvalCS_LO_tdubH = EvalCS_LO_tdubH/VgsWgt


 end FUNCTION EvalCS_LO_tdubH


  FUNCTION EvalCS_LO_tbardubbarH_MPI(yRnd,VgsWgt,res)
    implicit none
    integer :: EvalCS_LO_tbardubbarH_MPI
    real(8) ::  yRnd(*),res(*),VgsWgt

    res(1) = EvalCS_LO_tbardubbarH(yRnd,VgsWgt)
    EvalCS_LO_tbardubbarH_MPI=0
    RETURN
  END FUNCTION EvalCS_LO_tbardubbarH_MPI





 FUNCTION EvalCS_LO_tbardubbarH(yRnd,VgsWgt)
    use ModProcess
    use ModParameters
    use ModKinematics
    use ModMisc
    use ModSingleTopHAmps
    implicit none
    real(8) :: EvalCS_LO_tbardubbarH,yRnd(1:VegasMxDim),VgsWgt
    real(8) :: Ehat,MH_Inv,eta1,eta2,ISFac,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles,PSWgt,PSWgt2,PSWgt3,pdf(-6:6,1:2)
    real(8) ::PObs(1:NumMaxHisto) 
    real(8) :: MomExt(1:4,1:11),MomExtFlat(1:7,1:4),p4Dp5,p4Dp7,LO_Res_UnPol,s(10,10),p2Dp3,MomExtFlatDK(1:10,1:4)
    complex(8) :: za(10,10),zb(10,10),LOAmp(-6:6,-6:6,1:2),CoupFac,decay_amp(1:2)
    real(8) :: s12,s13,s1e4,s1k4,s15,s23,s2e4,s2k4,s25,s3e4,s3k4,s35,se45,sk45,se4k4,ColFac
    integer :: NBin(1:NumMaxHisto),NHisto,j,k
    logical :: applyPSCut



    EvalCS_LO_tbardubbarH = 0d0

    call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    call SmearExternal(yRnd(18),M_H,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)
    MH_Inv = M_H

   if( EHat.le.m_Top+MH_Inv ) then
      EvalCS_LO_tbardubbarH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3ArbMass(EHat,(/MH_Inv,M_Top,0d0/),yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ISFac = MomCrossing(MomExt)
   call setPDFs(eta1,eta2,MuFac,pdf)



!!!!!! MCFM phase space point to check                                                                                                                                         
!   MomExt(1:4,1)=  (/12.804565518142741d0,   0.0000000000000000d0,  0.0000000000000000d0, 12.804565518142741d0 /)
!   MomExt(1:4,2)=  (/6853.0396776786965d0,   0.0000000000000000d0,  0.0000000000000000d0, -6853.0396776786965d0 /)
!   MomExt(1:4,3)=  (/ 2151.2050020915885d0,   59.848489328696687d0, -121.56893937712810d0, -2143.2326314069742d0 /)
!   MomExt(1:4,4)=  (/ 4272.0664685769034d0,   10.381035193439132d0,  36.524127576457872d0, -4268.3932731359500d0 /)
!   MomExt(1:4,5)=  (/ 442.57277252834774d0,  -70.229524522135819d0,  85.044811800670232d0, -428.60920761762918d0 /)
!   MomExt=MomExt*GeV
!!!!!end MCFM                                                                                                                                                       
!   print *, g_weak**2
   CoupFac=2d0*g_weak**2/Vev*ci
   ColFac=9d0

   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
      PSWgt = PSWgt * PSWgt2
! usual decay atop(p4) --> ab(p6) + e-(p7) + nubar(p8)
!      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
   ENDIF
   IF( HDECAYS.NE.0 .AND. HDECAYS.NE.-2 ) THEN
      call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,9:10),PSWgt3)
      PSWgt = PSWgt * PSWgt3   !* MHJacobi                                                                            
   ENDIF


! MCFM point -- decayed tops
!   MomExt(1:4,1)= (/ 3650.0144958499450d0,   0.0000000000000000d0,   0.0000000000000000d0, 3650.0144958499450d0/)       
!   MomExt(1:4,2)= (/ 4416.8360375648035d0,   0.0000000000000000d0,   0.0000000000000000d0,  -4416.8360375648035d0/)       
!   MomExt(1:4,3)= (/  2524.5820457377790d0,   1740.5566521568635d0,   1463.5730674468346d0, -1089.0614870328973d0/)       
!   MomExt(1:4,4)= (/  2820.6965893355869d0,  -441.49038581856871d0,  -1901.3596045977481d0, -2028.8711497586719d0/)       
!   MomExt(1:4,5)= (/  2721.5718983413826d0,  -1299.0662663382948d0,   437.78653715091355d0,  2351.1110950767106d0/)       
!   MomExt(1:4,6)= (/  584.28028358321626d0,  -116.34103819613028d0,  -429.22474142090778d0, -378.96481890498103d0/)       
!   MomExt(1:4,7)= (/  2142.4030726221968d0,  -306.25535180144135d0,  -1400.3724496098203d0, -1592.1857892368578d0/)   
!   MomExt(1:4,8)= (/  94.013233130173830d0,  -18.893995820997077d0,  -71.762413567020076d0, -57.720541616833088d0/)       
!
!   MomExt=MomExt*GeV
! end MCFM

   call Kinematics_TBH(0,MomExt(1:4,1:11),(/4,5,3,1,2,0,6,7,8,9,10/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      EvalCS_LO_tbardubbarH = 0d0
      return
   endif

  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt

! setup momenta for spinor helicity products                                                                                                                                 
   p4Dp5=MomExt(1,5)*MomExt(1,4)-MomExt(2,5)*MomExt(2,4)-MomExt(3,5)*MomExt(3,4)-MomExt(4,5)*MomExt(4,4)
   p4Dp7=MomExt(1,7)*MomExt(1,4)-MomExt(2,7)*MomExt(2,4)-MomExt(3,7)*MomExt(3,4)-MomExt(4,7)*MomExt(4,4)
   p2Dp3=MomExt(1,2)*MomExt(1,3)-MomExt(2,2)*MomExt(2,3)-MomExt(3,2)*MomExt(3,3)-MomExt(4,2)*MomExt(4,3)
   MomExtFlat(1,1:4)=MomExt(1:4,1)
   MomExtFlat(2,1:4)=MomExt(1:4,2)
   MomExtFlat(3,1:4)=m_H**2/2d0/p2Dp3*MomExt(1:4,2)
   MomExtFlat(4,1:4)=MomExt(1:4,3)-MomExtFlat(3,1:4)
   MomExtFlat(5,1:4)=m_Top**2*MomExt(1:4,5)/2d0/p4Dp5
   MomExtFlat(6,1:4)=MomExt(1:4,4)-MomExtFlat(5,1:4)
   MomExtFlat(7,1:4)=MomExt(1:4,5)

    ! use different flattened momenta for anti-top decays
      IF (TOPDECAYS .NE. 0) THEN 
         MomExtFlatDK(1:7,1:4)=MomExtFlat(1:7,1:4)
         ! overwrite flattened top momenta         
         MomExtFlatDK(5,1:4)=m_Top**2*MomExt(1:4,7)/2d0/p4Dp7
         MomExtFlatDK(6,1:4)=MomExt(1:4,4)-MomExtFlatDK(5,1:4)
         ! top decay products 
         MomExtFlatDK(8,1:4)=MomExt(1:4,6)
         MomExtFlatDK(9,1:4)=MomExt(1:4,7)
         MomExtFlatDK(10,1:4)=MomExt(1:4,8)

      ENDIF



! Get spinor products
      za=0d0
      zb=0d0
      s=0d0
      IF (TOPDECAYS .EQ. 0) THEN
         do j=1,7
            call convert_to_MCFM(MomExtFlat(j,1:4))
         enddo
         MomExtFlat(1,1:4)=-MomExtFlat(1,1:4)
         MomExtFlat(2,1:4)=-MomExtFlat(2,1:4)

         call spinoru(7,MomExtFlat,za,zb,s)
      ELSE
         do j=1,10
            call convert_to_MCFM(MomExtFlatDK(j,1:4))
         enddo
         MomExtFlatDK(1,1:4)=-MomExtFlatDK(1,1:4)
         MomExtFlatDK(2,1:4)=-MomExtFlatDK(2,1:4)

         call spinoru(10,MomExtFlatDK,za,zb,s)
      ENDIF

   IF (CORRECTION .EQ. 0) THEN
      LOAmp=(0d0,0d0)
      IF (TOPDECAYS .EQ. 0) THEN
         decay_amp(1)=dcmplx(1d0,0d0)
         decay_amp(2)=dcmplx(1d0,0d0)
      ELSE
         call atdecay(5,6,8,9,10,za,zb,decay_amp)
         
      ENDIF
      call dbbarhtbaruamp(1,2,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(Dn_,ABot_,1:2))
      call dbbarhtbaruamp(2,1,3,4,5,6,7,za,zb,s,decay_amp,LOAmp(ABot_,Dn_,1:2))
      call dbbarhtbaruamp(7,2,3,4,5,6,1,za,zb,s,decay_amp,LOAmp(AUp_,ABot_,1:2))
      call dbbarhtbaruamp(7,1,3,4,5,6,2,za,zb,s,decay_amp,LOAmp(ABot_,AUp_,1:2))

      ! coupling factors in decay incl in tdecay function
      LOAmp=LOAmp*CoupFac

!      IF (TOPDECAYS .EQ. 0) THEN
!      print *, "d_ab",(abs(LOAmp(Dn_,ABot_,1))**2+abs(LOAmp(Dn_,ABot_,2))**2) * ColFac * ISFac * GeV**2
!      print *, "ab_d",(abs(LOAmp(ABot_,Dn_,1))**2+abs(LOAmp(ABot_,Dn_,2))**2) * ColFac * ISFac * GeV**2
!      print *, "aup_ab",(abs(LOAmp(AUp_,ABot_,1))**2+abs(LOAmp(AUp_,ABot_,2))**2) * ColFac * ISFac * GeV**2
!      print *, "ab_au",(abs(LOAmp(ABot_,AUp_,1))**2+abs(LOAmp(ABot_,AUp_,2))**2) * ColFac * ISFac * GeV**2
!      pause
!      ELSE
!      print *, "d_ab",(abs(LOAmp(Dn_,ABot_,1))**2+abs(LOAmp(Dn_,ABot_,2))**2) * ColFac * ISFac * GeV**6
!      print *, "ab_d",(abs(LOAmp(ABot_,Dn_,1))**2+abs(LOAmp(ABot_,Dn_,2))**2) * ColFac * ISFac * GeV**6
!      print *, "aup_ab",(abs(LOAmp(AUp_,ABot_,1))**2+abs(LOAmp(AUp_,ABot_,2))**2) * ColFac * ISFac * GeV**6
!      print *, "ab_au",(abs(LOAmp(ABot_,AUp_,1))**2+abs(LOAmp(ABot_,AUp_,2))**2) * ColFac * ISFac * GeV**6
!      pause
!      ENDIF

      LO_Res_UnPol= &
           + (abs(LOAmp(Dn_,ABot_,1))**2+abs(LOAmp(Dn_,ABot_,2))**2)   * (pdf(Dn_,1)*pdf(ABot_,2)  + pdf(Str_,1)*pdf(ABot_,2)) &
           + (abs(LOAmp(ABot_,Dn_,1))**2+abs(LOAmp(ABot_,Dn_,2))**2)   * (pdf(ABot_,1)*pdf(Dn_,2)  + pdf(ABot_,1)*pdf(Str_,2)) &
           + (abs(LOAmp(AUp_,ABot_,1))**2+abs(LOAmp(AUp_,ABot_,2))**2) * (pdf(AUp_,1)*pdf(ABot_,2) + pdf(AChm_,1)*pdf(ABot_,2)) &
           + (abs(LOAmp(ABot_,AUp_,1))**2+abs(LOAmp(ABot_,AUp_,2))**2) * (pdf(ABot_,1)*pdf(AUp_,2) + pdf(ABot_,1)*pdf(AChm_,2))

      LO_Res_UnPol = LO_Res_UnPol * ColFac * ISFac * WidthExpansion
      EvalCS_LO_tbardubbarH = LO_Res_Unpol * PreFac

   ENDIF

   if( IsNan(EvalCS_LO_tbardubbarH) ) then
        print *, "NAN:",EvalCS_LO_tbardubbarH
        print *, yRnd(:)
        print *, LO_Res_UnPol

        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_LO_tbardubbarH = 0d0
!         pause                                                                                                                                                                 
        return
   endif

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),EvalCS_LO_tbardubbarH,BinValue=PObs(NHisto))
!      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggH)                                                                                                                     
   enddo
   EvalCounter = EvalCounter + 1

   EvalCS_LO_tbardubbarH = EvalCS_LO_tbardubbarH/VgsWgt




 end FUNCTION EvalCS_LO_tbardubbarH

!! SUBROUTINE TDECAY(k4,e4,b,ep,nu,za,zb,dkamp)
!!    use ModParameters
!!    implicit none
!!    integer :: k4,e4,b,ep,nu
!!    complex(8) :: za(:,:),zb(:,:),dkamp(1:2)
!!    real(8) :: NWAFactor_Top,NWAFactor_W
!!    complex(8) :: WProp
!!    
!!! if one flattens the top wrt to e, then amp(2) = 0.
!!    dkamp(1) = za(b,nu)*zb(ep,e4)
!!    dkamp(2) = m_top * za(b,nu)*zb(ep,k4)/za(e4,k4)
!!
!!!    NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top(0)*m_Top)
!!!    NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W(0)*m_W)
!!! RR remove  -- compare with MCFM
!!    NWAFactor_Top = 1d0/(dsqrt(2d0*Ga_TopExp*m_Top))
!!    NWAFactor_W   = 1d0/(dsqrt(2d0*Ga_WExp*m_W))
!!
!!    WProp = (0d0,-1d0)*NWAFactor_W
!!
!!!    dkamp = dkamp * WProp * NWAFactor_Top * g_weak**2
!!    dkamp = dkamp * WProp * NWAFactor_Top * g_weak**2
!!
!!
!!  end SUBROUTINE TDECAY

end MODULE ModCrossSection_TH

      


