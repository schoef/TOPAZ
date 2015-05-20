MODULE ModCrossSection_eeTTB
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45



 CONTAINS






FUNCTION EvalCS_1L_ee_ttb_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ee_ttb_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ee_ttb(yRnd,VgsWgt)
EvalCS_1L_ee_ttb_MPI=0
RETURN
END FUNCTION





FUNCTION EvalCS_1L_ee_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_eeTTB
use ModParameters
use ModIntDipoles_eeTTB
implicit none
real(8) ::  EvalCS_1L_ee_ttb,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,Virt_Res_Pol
real(8) :: LO_Res_Unpol,Virt_Res_Unpol
integer :: iHel
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,ISFac
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:6)
logical :: applyPSCut
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2), z
real(8) :: IDip(3), IDipAmp,MadGraphRes
integer :: NHisto,NBin(1:NumMaxHisto),npdf
! real(8) :: MadGraphRes
! include 'misc/global_import'
include "vegas_common.f"

  EvalCS_1L_ee_ttb = 0d0
  yrnd(1:2) = 1d0! fix Ehat to ECollider


  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! Flat mapping
  if( EHat.le.2d0*m_Top) then
     EvalCS_1L_ee_ttb = 0d0
     return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt)
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,3),yRnd(5:8),.false.,MomDK(1:4,1:3),PSWgt2)
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(9:12),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
  ENDIF


  call Kinematics_eeTTBAR(.false.,MomExt,MomDK,applyPSCut,NBin)
  if( applyPSCut ) then
     EvalCS_1L_ee_ttb = 0d0
     return
  endif


  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  RunFactor = RunAlphaS(NLOParam,MuRen)
  
  LO_Res_Unpol   = 0d0
  Virt_Res_Unpol = 0d0
  !------------ LO --------------
  IF( CORRECTION.EQ.0 ) THEN

     ISFac = MomCrossing(MomExt)
     IF( TOPDECAYS.GE.1 ) THEN
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF

     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_ee_tbt(LO_Res_Pol)
        LO_Res_UnPol = LO_Res_UnPol + dreal(LO_Res_Pol*dconjg(LO_Res_Pol))
     enddo
     EvalCS_1L_ee_ttb = ISFac*PreFac  * WidthExpansion  * LO_Res_Unpol * alpha4Pi**2 * 3d0!=Nc = delta_ii

     
!  call coupsm(0)
!  call SEMEP_TTB((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,4),MomExt(1:4,3)/)*100d0,MadGraphRes)
!  print *, "TOPAZ", ISFac*LO_Res_Unpol *alpha4Pi**2* 3d0
!  print *, "Madgr",MadGraphRes   *(alpha**2)/(137d0**-2)
!  print *, "ratio",ISFac*LO_Res_Unpol*alpha4Pi**2* 3d0 /(MadGraphRes *(alpha**2)/(137d0**-2))
!  pause
     

  ELSEIF(CORRECTION.EQ.1) THEN
  
     ISFac = MomCrossing(MomExt)
     IF( TOPDECAYS.GE.1 ) THEN
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF


     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_ee_tbt(LO_Res_Pol)
        call VirtAmp_ee_Z_ttb(Virt_Res_Pol)
        LO_Res_UnPol   = LO_Res_UnPol   + dreal(LO_Res_Pol*dconjg(LO_Res_Pol))
        Virt_Res_UnPol = Virt_Res_UnPol + dreal(LO_Res_Pol*dconjg(Virt_Res_Pol))
     enddo
     EvalCS_1L_ee_ttb = ISFac*PreFac * Virt_Res_Unpol  * 4d0/3d0*3d0!=  Tr[T^aT^a]
     EvalCS_1L_ee_ttb = EvalCS_1L_ee_ttb * (alpha_sOver2Pi *RunFactor) * alpha4Pi**2

     
     
  ELSEIF(CORRECTION.EQ.3) THEN

! print *, "Virt.correction",EvalCS_1L_ee_ttb

     ISFac = MomCrossing(MomExt)
     IF( TOPDECAYS.GE.1 ) THEN
        call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
        call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
     ENDIF

     LO_Res_Unpol   = 0d0
     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_ee_tbt(LO_Res_Pol)
        LO_Res_UnPol = LO_Res_UnPol + dreal(LO_Res_Pol*dconjg(LO_Res_Pol))
     enddo
     
     call IntDip_eettb(MomExt(1:4,1:4), IDip)
     IDip = IDip * (alpha_s*RunFactor)/(2d0*Pi) 

     IDipAmp = IDip(1) * LO_Res_UnPol* alpha4Pi**2  * ISFac * 3d0!=Nc = delta_ii
     EvalCS_1L_ee_ttb = IDipAmp * PreFac

! print *, "Int.dipole", EvalCS_1L_ee_ttb
! pause
  ENDIF

  
  do NHisto=1,NumHistograms
     call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ee_ttb)
  enddo
  EvalCS_1L_ee_ttb = EvalCS_1L_ee_ttb/VgsWgt

  
RETURN
END FUNCTION EvalCS_1L_ee_ttb






FUNCTION EvalCS_Real_ee_ttbg(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_eeTTB
use ModParameters
use ModDipoles_eeTTB
implicit none
real(8) ::  EvalCS_Real_ee_ttbg,yRnd(1:VegasMxDim),VgsWgt,DipoleResult,EvalCS_Dips_ee_ttb
complex(8) :: Res_Po, prim2
real(8) :: Res_UnPol
integer :: iHel
real(8) :: EHat,eta1,eta2,RunFactor,PSWgt,PS,PSWgt2,PSWgt3,xFrag, sHatJacobi
real(8) :: FluxFac, PreFac, ISFac,MadGraphRes
real(8) :: MomExt(1:4,1:5),MomDK(1:4,1:6), MomExtTd(1:4,1:4), MomDKTd(1:4,1:6)
integer :: NBin(1:NumMaxHisto),NHisto
logical :: applyPSCut,applySingCut
real(8) :: colf
real(8),parameter :: CF = 4d0/3d0
integer :: nDip
real(8) :: resdip, dipoles
real(8) :: resLO
complex(8) :: ampLO
include "vegas_common.f"
integer :: i1,i2,i3,i4,i5


  EvalCS_Real_ee_ttbg = 0d0
  EvalCS_Dips_ee_ttb = 0d0
  NumExtParticles = 5

  
  yrnd(1:2) = 1d0! fix Ehat to ECollider
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi) ! no mapping 
  
  if( EHat.le.2d0*m_Top ) then
      EvalCS_Real_ee_ttbg = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to3(EHat,yRnd(3:7),MomExt(1:4,1:5),PSWgt)

  call CheckSing(MomExt,applySingCut)
  if( applySingCut ) then
     EvalCS_Real_ee_ttbg = 0d0
     SkipCounter = SkipCounter + 1
     return
  endif
  
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomDK(1:4,1:3),PSWgt2) ! anti-top
     call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomDK(1:4,4:6),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
  ENDIF

  call Kinematics_eeTTBAR(.true.,MomExt,MomDK,applyPSCut,NBin)
  PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
  colf =  CF * 3d0! = TR[T^a.T^a]

  RunFactor = RunAlphaS(NLOParam,MuRen)
  ISFac = MomCrossing(MomExt)

  IF( TOPDECAYS.GE.1 ) THEN
     call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
     call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
  ENDIF
  if( applyPSCut ) then
     EvalCS_Real_ee_ttbg = 0d0
  else

     ISFac = MomCrossing(MomExt)
     Res_UnPol = 0d0
     do iHel=1,NumHelicities
        call HelCrossing(Helicities(iHel,1:NumExtParticles))
        call SetPolarizations()
        call Tree_ee_tbtg_f(prim2)
        Res_UnPol = Res_UnPol +  dreal(prim2*dconjg(prim2))
     enddo!helicity loop

     EvalCS_Real_ee_ttbg = Res_UnPol * (4d0*Pi*alpha_s*RunFactor)*alpha4Pi**2 * PreFac * ISFac * colf
     
! call coupsm(0)
! call SEMEP_TTBG((/MomExt(1:4,1),MomExt(1:4,2),MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,3)/)*100d0,MadGraphRes)
! print *, "TOPAZ", ISFac*Res_Unpol * (4d0*Pi*alpha_s*RunFactor) *alpha4Pi**2 * colf
! print *, "Madgr",MadGraphRes *1d4  *(alpha**2)/(137d0**-2)    *(alpha_s)/(0.13d0)
! print *, "ratio", (ISFac*Res_Unpol * (4d0*Pi*alpha_s*RunFactor) *alpha4Pi**2 * colf)     /(MadGraphRes*1d4 *(alpha**2)/(137d0**-2)*(alpha_s)/(0.13d0))
! pause

     do NHisto=1,NumHistograms
        call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ee_ttbg)
     enddo
     EvalCounter = EvalCounter + 1

  endif!applyPSCut

  
  

     dipoles = 0d0
     colf =  CF * 3d0  * 2d0
     NumExtParticles = 4
     do nDip = 1,2
        call Dipoles_eettb(nDip, MomExt(1:4,1:5), MomExtTd(1:4,1:4), resdip)
        if (resdip .EQ. 0d0) cycle
        Crossing(:) = (/3,4,-1,-2,0/) ! For dipoles
        ISFac = MomCrossing(MomExtTd)
        Crossing(:) = (/4,5,-1,-2,3/)        
        IF( TOPDECAYS.NE.0 ) THEN
           call EvalPhasespace_TopDecay(MomExtTd(1:4,3),yRnd(8:11),.false.,MomDKTd(1:4,1:3),PSWgt2)
           call EvalPhasespace_TopDecay(MomExtTd(1:4,4),yRnd(12:15),.false.,MomDKTd(1:4,4:6),PSWgt3)
           PSWgt = PSWgt * PSWgt2*PSWgt3! this is better never be used anywhere below
           
           call Kinematics_eeTTBAR(.false.,MomExtTd,MomDKTd,applyPSCut,NBin)

           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
              call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
              
              resLO = 0d0              
              i1=0;  i2=0;
              do i3 = -1,1,2
                 ExtParticle(1)%Helicity = i1
                 ExtParticle(2)%Helicity = i2
                 ExtParticle(3)%Helicity = i3
                 ExtParticle(4)%Helicity =-i3
                 call SetPolarizations()
                 call Tree_ee_tbt(ampLO)
                 resLO = resLO +  dreal(ampLO * dconjg(ampLO))
              enddo
              
           endif
           
        ELSE
           
           call Kinematics_eeTTBAR(.false.,MomExtTd,MomDK,applyPSCut,NBin)
           
           if ( applyPSCut ) then
              resdip = 0d0
           else
              
              resLO = 0d0
              do i1 = -1,1,2
                 do i2 = -1,1,2
                    do i3 = -1,1,2
                       ExtParticle(1)%Helicity = i1
                       ExtParticle(2)%Helicity = i2
                       ExtParticle(3)%Helicity = i3
                       ExtParticle(4)%Helicity = -i3
                       call SetPolarizations()
                       call Tree_ee_tbt(ampLO)
                       resLO = resLO +  dreal(ampLO * dconjg(ampLO))
                    enddo
                 enddo
              enddo
              
           endif
           
        ENDIF !TopDecays
        
        resdip = resdip * resLO
        resdip = resdip * (4d0*Pi*alpha_s*RunFactor)*alpha4Pi**2 * PreFac * ISFac * colf
        
        do NHisto=1,NumHistograms
           call intoHisto(NHisto,NBin(NHisto),resdip)
        enddo
     
        dipoles = dipoles + resdip
     
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !-- dipoles check (only one soft singularity)
!  print *, 'inv    ', (MomExt(1:4,1).dot.MomExt(1:4,3))/Ehat**2,(MomExt(1:4,2).dot.MomExt(1:4,3))/Ehat**2
!  print *, 'real   ', EvalCS_Real_ee_ttbg
!  print *, 'dipoles', dipoles
!  print *, 'dipoles/real+1', dipoles/EvalCS_Real_ee_ttbg + 1d0
!  pause 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  EvalCS_Real_ee_ttbg = (EvalCS_Real_ee_ttbg+dipoles)/VgsWgt

  RETURN

END FUNCTION EvalCS_Real_ee_ttbg







FUNCTION EvalCS_NLODK_ee_ttb(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes_eeTTB
use ModParameters
use ModHadrWDecay
implicit none
#DEFINE TheReal_Top 2233
#DEFINE TheDipole_Top 2244
#DEFINE TheDipole_ATop 2255
#DEFINE TheWm_decay 2260
#DEFINE TheReal_Top_Wdecay 2266
#DEFINE TheDipole_Top_Wdecay 2277
#DEFINE TheWp_decay 2280
#DEFINE TheDipole_ATop_Wdecay 2288
#DEFINE TheEnd 2299
real(8) ::  EvalCS_NLODK_ee_ttb, yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,NLO_Res_Pol, NLO_Res_Unpol,Dip_Res_Unpol
integer :: iHel,jHel,kHel,GluHel,iPrimAmp,jPrimAmp,ndip
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac,dip_res_w
real(8) :: MomExt(1:4,1:NumExtParticles),MomDK(1:4,1:7),MomDKTd(1:4,1:6),MomDKx(1:4,1:7)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,RunFactor
integer :: NBin(1:NumMaxHisto),NHisto
real(8) :: pbDpg,ptDpg,ptDpb,z,omz,Dipole,rsq,y
real(8), parameter :: CF=4d0/3d0,PhotonCouplCorr=2d0
real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4)
integer,parameter :: up=1,dn=2,glu=1
include "vegas_common.f"




  EvalCS_NLODK_ee_ttb = 0d0

  yrnd(1:2) = 1d0! fix Ehat to ECollider
  call PDFMapping(1,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
  if( EHat.le.2d0*m_Top) then
      EvalCS_NLODK_ee_ttb = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)

  call EvalPhaseSpace_2to2(EHat,yRnd(3:4),MomExt(1:4,1:4),PSWgt1)


IF( CORRECTION.EQ.4 ) THEN
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8) ,MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(9:12),MomDK(1:4,4:6),PSWgt3)

   call Kinematics_eeTTBAR(.false.,MomExt,MomDK,applyPSCut,NBin)!,xJPsiFrag=xFrag)
   if( applyPSCut ) then
      EvalCS_NLODK_ee_ttb = 0d0
      return
   endif

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)

!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_ee_tbt(LO_Res_Pol)

      call TopDecay(ExtParticle(1),DK_1L_T,MomDK(1:4,1:3))
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + dreal(LO_Res_Pol*dconjg(NLO_Res_Pol))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii

   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

   

!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_ee_tbt(LO_Res_Pol)

      call TopDecay(ExtParticle(2),DK_1L_T,MomDK(1:4,4:6))
      call Tree_ee_tbt(NLO_Res_Pol)

      !--F By looking at Markus' implementation I assume that the factor 2 in 2*dreal(Alo*conjg(Avirt)) is already taken care of
      NLO_Res_UnPol= NLO_Res_UnPol + dreal(LO_Res_Pol*dconjg(NLO_Res_Pol))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo




IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.4 ) THEN
!----------------------------------------
! one loop correction to W- decay |
!----------------------------------------
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_ee_tbt(LO_Res_Pol)

      call TopDecay(ExtParticle(1),DK_1L_Q,MomDK(1:4,1:3))
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + dreal(LO_Res_Pol*dconjg(NLO_Res_Pol))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
ENDIF


IF( TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.3 ) THEN
!-------------------------------------
! one loop correction to W+ decay   |
!-------------------------------------
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,4:6))
      call Tree_ee_tbt(LO_Res_Pol)

      call TopDecay(ExtParticle(2),DK_1L_Q,MomDK(1:4,4:6))
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + dreal(LO_Res_Pol*dconjg(NLO_Res_Pol))
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
ENDIF



ELSEIF( CORRECTION.EQ.5 ) THEN
!------------------------------------
! real correction to Anti-top decay |
!------------------------------------
   call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,3),yRnd(5:11),MomDK(1:4,1:4),PSWgt2) 
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(12:15),MomDK(1:4,5:7),PSWgt3) 
!   write(6,*) dsqrt(dabs(MomExt(1,4)**2-MomExt(2,4)**2-MomExt(3,4)**2-MomExt(4,4)**2))
!   pause

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)


   MomDKx(1:4,1:3) = MomDK(1:4,1:3) 
   MomDKx(1:4,4:6) = MomDK(1:4,5:7) 
   MomDKx(1:4,7) = MomDK(1:4,4) 
   call Kinematics_eeTTBAR(.true.,MomExt,MomDKx,applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_ATop

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2 
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_T,MomDK(1:4,1:4),GluonHel=GluHel) 
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7)) 
      call Tree_ee_tbt(NLO_Res_Pol)
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop 
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
! print *, "real  ", dble(NLO_Res_Unpol)


TheDipole_ATop continue

   call WTransform(MomDK(1:4,1:4),MomDKTd(1:4,1:3),pbDpg,ptDpg,ptDpb)
   MomDKTd(1:4,4:6) = MomDK(1:4,5:7)
   omz=ptDpg/(ptDpb+ptDpg-pbDpg)
   rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
   z=1d0-omz
   y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
   Dipole = - alpha_s*4d0*Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
   Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   call Kinematics_eeTTBAR(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) goto TheReal_Top
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol * Dipole * alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
! print *, "dipole", dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,1).dot.MomDK(1:4,4))/m_Top**2
! pause

TheReal_Top continue

!------------------------------------
! real correction to Top decay |
!------------------------------------
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_BG_W,MomExt(1:4,4),yRnd(9:15),MomDK(1:4,4:7),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   call Kinematics_eeTTBAR(.true.,MomExt,MomDK(1:4,1:7),applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_Top

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_RE_T,MomDK(1:4,4:7),GluonHel=GluHel)
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop
   enddo!helicity loop


   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
!   print *, 'real  ', dble(NLO_Res_Unpol)

TheDipole_Top continue

   call WTransform(MomDK(1:4,4:7),MomDKTd(1:4,4:6),pbDpg,ptDpg,ptDpb)
   MomDKTd(1:4,1:3) = MomDK(1:4,1:3)
   omz=ptDpg/(ptDpb+ptDpg-pbDpg)
   rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
   z=1d0-omz
   y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
   Dipole = - alpha_s*4d0*Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
   Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   call Kinematics_eeTTBAR(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) goto TheWm_decay
   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6))
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol *Dipole * alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo
! print *, "dipole", dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,4).dot.MomDK(1:4,7))/m_Top**2
! pause

TheWm_decay continue

IF (TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.4) THEN
!------------------------------
! real correction to W- decay |
!------------------------------
   call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,3),yRnd(5:11),MomDK(1:4,1:4),PSWgt2) 
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,4),yRnd(12:15),MomDK(1:4,5:7),PSWgt3) 
!   write(6,*) dsqrt(dabs(MomExt(1,4)**2-MomExt(2,4)**2-MomExt(3,4)**2-MomExt(4,4)**2))
!   pause

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)

   MomDKx(1:4,1:3) = MomDK(1:4,1:3) 
   MomDKx(1:4,4:6) = MomDK(1:4,5:7) 
   MomDKx(1:4,7) = MomDK(1:4,4) 
   call Kinematics_eeTTBAR(.true.,MomExt,MomDKx,applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_ATop_Wdecay

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2 
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_Q,MomDK(1:4,1:4),GluonHel=GluHel) 
      call TopDecay(ExtParticle(2),DK_LO,MomDK(1:4,5:7)) 
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop 
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!print *, "real              ", dble(NLO_Res_Unpol)

TheDipole_ATop_Wdecay continue

do ndip=1,2
   call wdec_trans(ndip,MomDK(1:4,1:4),MomDKTd(1:4,1:3),alpha_DKWff,dip_res_w)
   if( dip_res_w.eq.0d0 ) goto 20
   Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w

   MomDKTd(1:4,4:6) = MomDK(1:4,5:7)

   call Kinematics_eeTTBAR(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) then
      goto 20
      return
   endif

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3)) 
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6)) 
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   NLO_Res_Unpol = NLO_Res_Unpol * Dipole
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!   print *, "dipole", ndip, dble(NLO_Res_Unpol)
!   print *, "sing", (MomDK(1:4,2).dot.MomDK(1:4,4))/m_Top**2

20 continue
enddo ! dipole loop

ENDIF ! Wm decay real correction


TheWp_decay continue

IF (TOPDECAYS.EQ.2 .OR. TOPDECAYS.EQ.3) THEN
!------------------------------
! real correction to W+ decay |
!------------------------------
   call EvalPhasespace_TopDK(T_B_W,MomExt(1:4,3),yRnd(5:8),MomDK(1:4,1:3),PSWgt2)
   call EvalPhasespace_TopDK(T_B_WG,MomExt(1:4,4),yRnd(9:15),MomDK(1:4,4:7),PSWgt3)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)
   ISFac = MomCrossing(MomExt)
   call Kinematics_eeTTBAR(.true.,MomExt,MomDK(1:4,1:7),applyPSCut,NBin)
   if( applyPSCut ) goto TheDipole_Top_Wdecay

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
   do GluHel=-1,1,2
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDK(1:4,1:3))
      call TopDecay(ExtParticle(2),DK_RE_Q,MomDK(1:4,4:7),GluonHel=GluHel)
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop
   enddo!helicity loop


   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol) 
   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

!print *, "real              ", dble(NLO_Res_Unpol)

TheDipole_Top_Wdecay continue

do ndip=1,2
   call wdec_trans(ndip,MomDK(1:4,4:7),MomDKTd(1:4,4:6),alpha_DKWff,dip_res_w)
   if( dip_res_w.eq.0d0 ) goto 22
   Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w

   MomDKTd(1:4,1:3) = MomDK(1:4,1:3)

   call Kinematics_eeTTBAR(.false.,MomExt,MomDKTd,applyPSCut,NBin)
   if( applyPSCut ) then
      goto 22
      return
   endif

   NLO_Res_UnPol = (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomDKTd(1:4,1:3)) 
      call TopDecay(ExtParticle(2),DK_LO,MomDKTd(1:4,4:6)) 
      call Tree_ee_tbt(NLO_Res_Pol)

      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol*dconjg(NLO_Res_Pol)
   enddo!helicity loop

   !  normalization
   NLO_Res_Unpol = ISFac * PreFac * NLO_Res_Unpol* alpha4Pi**2 * 3d0!=Nc = delta_ii
   NLO_Res_Unpol = NLO_Res_Unpol * Dipole
   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol))
   enddo

! print *, "dipole", ndip, dble(NLO_Res_Unpol)
! print *, "sing", (MomDK(1:4,4).dot.MomDK(1:4,7))/m_Top**2

22 continue
enddo ! dipole loop

! pause

ENDIF ! Wp decay real correction

TheEnd continue


ENDIF

   EvalCS_NLODK_ee_ttb = EvalCS_NLODK_ee_ttb/VgsWgt
RETURN
END FUNCTION






END MODULE


