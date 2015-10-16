MODULE ModCrossSection_TTBH
use ModTopDecay
implicit none

integer,private,parameter :: NumMaxHisto=45


contains




FUNCTION EvalCS_1L_ttbggH_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ttbggH_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ttbggH(yRnd,VgsWgt)
EvalCS_1L_ttbggH_MPI=0
RETURN
END FUNCTION






FUNCTION EvalCS_1L_ttbggH(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_new
use ModUCuts_128
use ModUCuts_128_new
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHDecay
use ModIntDipoles_GGTTBGH
use ModWeighting
implicit none
real(8) ::  EvalCS_1L_ttbggH,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),FermionLoopPartAmp(1:4,-2:1)
complex(8) :: propH,BosonicPartAmp(1:3,-2:1),mydummy,ZPolVec(1:4),BarSpi(1:4),Spi(1:4)
integer :: iHel,jHel,kHel,iPrimAmp,jPrimAmp,BPrimAmp,APrimAmp,ListPrimAmps(14)
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,ZDKMatel,Msq_T_BWENU
real(8) :: MomExt(1:4,1:14),pHsq
logical :: applyPSCut,nonrenormcoupl
logical, save :: first=.true.
real(8) :: Col1Lf_ttbggH(2,4), Col1L_ttbggH(2,3)
real(8) :: MG_MOM(0:3,1:5),tmpmom(1:4)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0, Cf=4d0/3d0
real(8) :: eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac,AccPoles
real(8) :: pdf(-6:6,1:2),pdf_z(-6:6,1:2),xE,HOp(1:3),MH_Inv,PObs(1:NumMaxHisto)
real(8) :: QPtol,DPtol,couplZLL,couplGLL
real(8) :: tmpVcoupl,tmpAcoupl,tmpV2coupl,tmpA2coupl,LightLoopCoupl
integer :: NBin(1:NumMaxHisto),NHisto,PhotonCouplCorr=2d0,nHel(1:2)
integer :: ZQcoupl,jj,lastSister
integer :: QPredo(1:8,1:5)
complex(8) :: couplZQQ_left_dyn_old,couplZQQ_right_dyn_old,couplZTT_left_dyn_old,couplZTT_right_dyn_old,couplZTT_left_old,couplZTT_right_old,couplZQQ_left2_dyn_old,couplZQQ_right2_dyn_old,couplZTT_left2_dyn_old,couplZTT_right2_dyn_old,couplZTT_left2_old,couplZTT_right2_old
complex(8) :: LOPartAmp(1:NumBornAmps),RenormAmp(1:NumBornAmps)
real(8) :: Ren_Res_Pol,Ren_Res_UnPol,R_V,R_A
include 'misc/global_import'
include 'vegas_common.f'
complex(8) :: tmpBornResults(14),FullBornAmps(1:NumPrimAmps),RenormAmps(1:NumPrimAmps),HOO_RenormAmps(1:NumPrimAmps),HOO_RenormPartAmp(1:3)
real(8)   :: HOO_Ren_Res_Pol,HOO_Ren_Res_UnPol
integer :: MaxBins=2500,NumReWeightHist
integer, allocatable :: ReWeightHisto(:)
character :: reweightfilename*(80)
complex(8) :: LO_Res_UnPol_RW

EvalCS_1L_ttbggH = 0d0


  DPtol=1d-4
  QPtol=1d-3


   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    call SmearExternal(yRnd(18),M_H,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)
   MH_Inv = M_H
   if( EHat.le.2d0*m_Top+MH_Inv ) then
      EvalCS_1L_ttbggH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to3M(EHat,MH_Inv,yRnd(3:7),MomExt(1:4,1:5),PSWgt)
   call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

   ISFac = MomCrossing(MomExt)

   IF( TOPDECAYS.NE.0 ) THEN
      call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
      call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
      PSWgt = PSWgt * PSWgt2*PSWgt3
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
   ENDIF
   IF( HDECAYS.NE.0 .AND. HDECAYS.NE.-2 ) THEN
      call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
      PSWgt = PSWgt * PSWgt4   !* MHJacobi
   ENDIF





   call Kinematics_TTBARH(0,MomExt(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      EvalCS_1L_ttbggH = 0d0
      return
   endif


! read up the weights from the appropriate histograms
  if (Correction .eq. 0 .and. LO_ReWeighting) then
      NumReWeightHist=2
      allocate(ReWeightHisto(1:NumReWeightHist))
      ReWeightHisto=0
      ReWeightHisto(1)=8
      ReWeightHisto(2)=16
      reweightfilename="LHC14_82_2.20/101/LHC.101.kfactor_2DHist.dat"
      if (first) then
         allocate(Weights(1:NumReWeightHist,1:2500))
         ! this uses the max histo bins (2500) from mod_Kinematics -- can change this if needed, but need to make sure they match up
         call ReadWeights(ReWeightHisto,reweightfilename)
      print *, "weights read in"
      endif
     first = .false.
  endif




! ! we still need this for the massless qqZ couplings
!    pHsq=MomExt(1,3)*MomExt(1,3)-MomExt(2,3)*MomExt(2,3)-MomExt(3,3)*MomExt(3,3)-MomExt(4,3)*MomExt(4,3)
!    if ( HDecays.lt.10 .and. HDecays.gt.0 ) then
!       propH = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)
!    elseif (HDecays .gt. 10) then 
!       propH=(1d0,0d0)/(pHsq-m_Z**2+ci*Ga_ZExp*m_Z)
!    endif
!    propH  = pHsq * propH! the pZ^2 will cancel against a 1/kZ^2 term in the Z decay matrix element

   call SetPropagators()
   ExtParticle(5)%Pol(:) = 0d0
   ExtParticle(5)%Pol(1) = 1d0
   
   call SetPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   nHel(1:2) = getHelicity(0d0)

   
   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   Ren_Res_Pol=0d0
   Ren_Res_UnPol=0d0   
   HOO_Ren_Res_Pol=0d0
   HOO_Ren_Res_UnPol=0d0   
!------------ LO --------------
IF( Correction.EQ.0 ) THEN

   do iHel=nHel(1),nHel(2)
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      if( HDecays.gt.0 ) then
          call HDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
      endif
      do iPrimAmp=1,2!  set this explicitely to avoid a problem with warm-up run for Correction=1 
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      LO_Res_Pol = (0d0,0d0)
      do jPrimAmp=1,2
      do iPrimAmp=1,2
          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
      enddo
      enddo
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!------------ 1 LOOP --------------
ELSEIF( Correction.EQ.1 ) THEN

   do iHel=nHel(1),nHel(2)
       QPredo=-1
       call HelCrossing(Helicities(iHel,1:NumExtParticles))
       call SetPolarizations()
       if( HDecays.gt.0 ) then
           call HDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
       endif
! 
       do iPrimAmp=1,NumBornAmps
           call EvalTree(BornAmps(iPrimAmp))
! !! now can overwrite BornAmps with various renorm amps...
! !          print *, 'born before',iPrimAmp,BornAmps(iPrimAmp)%Result
!           FullBornAmps(iPrimAmp)=BornAmps(iPrimAmp)%Result
       enddo
! 
        LO_Res_Pol = (0d0,0d0)   ! MARKUS: has to be disabled for now because NumBornAmps counts more than just 1,2.
        do jPrimAmp=1,2
           do iPrimAmp=1,2
            LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)
        enddo
        enddo
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
! 
! 
! !!!! RR Higher Order Operator -- i.e. sigma_{mu,nu} q^{nu} --  renorm term ??????
! if (nonrenormcoupl) then
!    call SigmaRenorm_gg(MomExt(1:4,12:13),BornAmps,HOO_RenormAmps,HOO_Ren_Res_Pol)
!    HOO_Ren_Res_UnPol=HOO_Ren_Res_UnPol + HOO_Ren_Res_Pol
! endif
! 
        ListPrimAmps=(/1,2,3,5,7,10,13,14,15,18,21,22,23,25/)
! 
! ! ----------------- bosonic loops ---------------------------
! !------------------------------------------------------------
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.1 ) THEN
       do iPrimAmp=1,12
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=1,12
           call SetKirill(PrimAmps(iPrimAmp))
           call QuadCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=1,12
           call SetKirill(PrimAmps(iPrimAmp))
           call TripCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=1,12
           call SetKirill(PrimAmps(iPrimAmp))
           call DoubCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=1,12
           call SetKirill(PrimAmps(iPrimAmp))
           call SingCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=1,12
          call SetKirill(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
        enddo
 
 
 ! now combine into real prims:
        PrimAmps(3)%Result = PrimAmps(3)%Result+PrimAmps(4)%Result
        PrimAmps(5)%Result = PrimAmps(5)%Result+PrimAmps(6)%Result
        PrimAmps(7)%Result = PrimAmps(7)%Result+PrimAmps(8)%Result   + PrimAmps(9)%Result
        PrimAmps(10)%Result= PrimAmps(10)%Result+PrimAmps(11)%Result + PrimAmps(12)%Result
! 
! 
 ! check on poles -- bosonic loops only
        do iPrimAmp=1,6
           APrimAmp=ListPrimAmps(iPrimAmp)
           call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
           PrimAmps(APrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
           call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
! 
! ! RR  -- now we modify rdiv for this process --- 
! ! this one has been around for ages          
           rdiv(1)=rdiv(1)+3d0
! ! now a modification for the non-renorm couplings
!           if (nonrenormcoupl) then
!              rdiv(1)=rdiv(1)-1d0/2d0*HOO_RenormAmps(APrimAmp)/BornAmps(APrimAmp)%Result
!           endif
! !! --- end
!           
! 
! 
! ! RR REMOVE -- overwrite poles with analytic value to check SP
! !          print *, 'overwriting single poles with analytic result'
! !          PrimAmps(APrimAmp)%Result(-1)=rdiv(1)*BornAmps(APrimAmp)%Result
! !          print *, 'LO',BornAmps(APrimAmp)%Result
! !          print *, 'diff', APrimAmp,PrimAmps(APrimAmp)%Result(-1)/BornAmps(APrimAmp)%Result-rdiv(1)
! !          pause
! 
! 
           AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
! 
! ! QP 
           if ( AccPoles .gt. DPtol ) then
!           if (.true.) then
!              print *, "pole accuracy=",AccPoles
!              print *, "using QP"
              useQP=useQP+1
! !             QPredo(iPrimAmp,1)=APrimAmp
! !             QPredo(iPrimAmp,2:PrimAmps(APrimAmp)%NumSisters+1)=PrimAmps(APrimAmp)%Sisters(1:PrimAmps(APrimAmp)%NumSisters)
              if (PrimAmps(APrimAmp)%NumSisters .gt. 0) then
                 lastSister=PrimAmps(APrimAmp)%Sisters(PrimAmps(APrimAmp)%NumSisters)
              else
                 lastSister=APrimAmp
              endif
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call PentCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call QuadCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call TripCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call DoubCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call SingCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
              enddo
 
 !  sum up sisters to one prim
              do jPrimAmp=APrimAmp+1,lastSister
                 PrimAmps(APrimAmp)%Result=PrimAmps(APrimAmp)%Result+PrimAmps(jPrimAmp)%Result
              enddo
 
              call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
              PrimAmps(APrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
! 
              AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
!              print *, 'accpoles after QP', accpoles 
              if ( AccPoles .gt. QPtol) then
                  print *, 'QP fails: ', AccPoles,APrimAmp
                 PrimAmps(APrimAmp)%Result=0d0
                 pole_skipped=pole_skipped+1
                 SkipCounter = SkipCounter + 1
                 RETURN ! reject the whole event instead of just this primamp
              endif
           endif! QP
 
        enddo! iPrimAmp
 
 
        BosonicPartAmp(1,-2:1) = PrimAmps(PrimAmp1_15234)%Result(-2:1) &
             &                 - 1d0/Nc**2 *(  PrimAmps(PrimAmp1_15432)%Result(-2:1) ) 
        BosonicPartAmp(2,-2:1) = PrimAmps(PrimAmp1_15243)%Result(-2:1) &
             &                 - 1d0/Nc**2 *(  PrimAmps(PrimAmp1_15342)%Result(-2:1) )
 
        BosonicPartAmp(3,-2:1) = PrimAmps(PrimAmp1_15234)%Result(-2:1) &
                               + PrimAmps(PrimAmp1_15243)%Result(-2:1) &
                               + PrimAmps(PrimAmp1_13524)%Result(-2:1) &
                               + PrimAmps(PrimAmp1_14523)%Result(-2:1) &
                               + PrimAmps(PrimAmp1_15342)%Result(-2:1) &
                               + PrimAmps(PrimAmp1_15432)%Result(-2:1) 
 
       Col1L_ttbggH = 0d0
       Col1L_ttbggH(1,1)= 4d0 * Cf**2 * Nc**2  ! = 64
       Col1L_ttbggH(1,2)= - 2d0 * Cf * Nc      ! =-8
       Col1L_ttbggH(1,3)= 2d0 * Cf * Nc        ! = 8
       Col1L_ttbggH(2,1)= Col1L_ttbggH(1,2)
       Col1L_ttbggH(2,2)= Col1L_ttbggH(1,1)
       Col1L_ttbggH(2,3)= Col1L_ttbggH(1,3)
       NLO_Res_Pol(-2:1) = (0d0,0d0)
       do jPrimAmp=1,3
          do iPrimAmp=1,2
           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1L_ttbggH(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result * dconjg(BosonicPartAmp(jPrimAmp,-2:1)) )
          enddo
       enddo
       NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)
! 
! 
! !! print *, ""
! !! print *, "1-loop amplitudes:"
! !! print *, "bosonic loops"
! !! print *, "|Ta*Tb|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(1,1)*dreal( BornAmps(1)%Result * dconjg( BosonicPartAmp(1,-2:1) ))
! !! print *, "|Tb*Ta|^2*Re[M0*cong(M1)]:", PhotonCouplCorr*ColLO_ttbgg(2,2)*dreal( BornAmps(2)%Result * dconjg( BosonicPartAmp(2,-2:1) ))
! !! print *, "check",cdabs( BosonicPartAmp(3,-2:1) )
! 
    ENDIF
! 
! 
! 
! 
! 
! !------------ fermionic loops -------------------------------
! !------------------------------------------------------------
! 
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.2 ) THEN
       do iPrimAmp=13,26
           call SetKirill(PrimAmps(iPrimAmp))
           call PentCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=13,26
           call SetKirill(PrimAmps(iPrimAmp))
           call QuadCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=13,26
           call SetKirill(PrimAmps(iPrimAmp))
           call TripCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=13,26
           call SetKirill(PrimAmps(iPrimAmp))
           call DoubCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=13,26
           call SetKirill(PrimAmps(iPrimAmp))
           call SingCut_new(PrimAmps(:),iPrimAmp)
        enddo
 
       do iPrimAmp=13,26
          call SetKirill(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
       enddo
 
 
! the fermion loops are combined into gauge invariant prims
        PrimAmps(PrimAmp2m_12534)%Result = PrimAmps(PrimAmp2m_12534)%Result +PrimAmps(PrimAmp2m_12354)%Result +PrimAmps(PrimAmp2m_12345)%Result 
        PrimAmps(PrimAmp2m_12543)%Result = PrimAmps(PrimAmp2m_12543)%Result +PrimAmps(PrimAmp2m_12453)%Result +PrimAmps(PrimAmp2m_12435)%Result 
        PrimAmps(PrimAmp2m_13254)%Result = PrimAmps(PrimAmp2m_13254)%Result +PrimAmps(PrimAmp2m_13245)%Result           
        PrimAmps(PrimAmp2m_14253)%Result = PrimAmps(PrimAmp2m_14253)%Result +PrimAmps(PrimAmp2m_14235)%Result   
! 
!        
! ! check on poles
        do iPrimAmp=7,14
           APrimAmp=ListPrimAmps(iPrimAmp)
           call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
           PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
           call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
 
           rdiv(1)=rdiv(1)+1.0d0/3.0d0
           AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
 
 
 ! QP 
           if ( AccPoles .gt. DPtol ) then
!          if (.true.) then
              useQP=useQP+1
              if (PrimAmps(APrimAmp)%NumSisters .gt. 0) then
                 lastSister=PrimAmps(APrimAmp)%Sisters(PrimAmps(APrimAmp)%NumSisters)
              else
                 lastSister=APrimAmp
              endif
 
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call PentCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call QuadCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call TripCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call DoubCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call SingCut_128_new(PrimAmps(:),jPrimAmp)
              enddo
              do jPrimAmp=APrimAmp,lastSister
                 call SetKirill(PrimAmps(jPrimAmp))
                 call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
              enddo
 
 !  sum up sisters to one prim
              do jPrimAmp=APrimAmp+1,lastSister
                 PrimAmps(APrimAmp)%Result=PrimAmps(APrimAmp)%Result+PrimAmps(jPrimAmp)%Result
              enddo
 
              call RenormalizeUV(PrimAmps(APrimAmp),BornAmps(APrimAmp),MuRen**2)
              PrimAmps(APrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(APrimAmp)%Result(-2:1)
 
              AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
              if ( AccPoles .gt. QPtol) then
                  print *, 'QP fails: ', AccPoles,APrimAmp
                 PrimAmps(APrimAmp)%Result=0d0
                 pole_skipped=pole_skipped+1
                 SkipCounter = SkipCounter + 1
                 RETURN ! reject the whole event instead of just this primamp                
              endif
           endif! QP
         enddo! iPrimAmp
 
 
!       if (HDecays .gt. 0) then
!          PrimAmps(PrimAmp2_12534)%Result(:) = PrimAmps(PrimAmp2_12534)%Result(:) * propH * couplZLL
!          PrimAmps(PrimAmp2_12543)%Result(:) = PrimAmps(PrimAmp2_12543)%Result(:) * propH * couplZLL
!          PrimAmps(PrimAmp2_13254)%Result(:) = PrimAmps(PrimAmp2_13254)%Result(:) * propH * couplZLL
!          PrimAmps(PrimAmp2_14253)%Result(:) = PrimAmps(PrimAmp2_14253)%Result(:) * propH * couplZLL
!       endif
 
 ! now combine into fermloops partials, with appropriate Z-coupls
       FermionLoopPartAmp=0d0
! 

           FermionLoopPartAmp(1,-2:1)=0d0 &
                & + nf_light*PrimAmps(PrimAmp2_15234)%Result &
                & + PrimAmps(PrimAmp2m_15234)%Result &
                & + PrimAmps(PrimAmp2m_12534)%Result
           
           FermionLoopPartAmp(2,-2:1)=0d0 &
                & + nf_light*PrimAmps(PrimAmp2_15243)%Result &
                & + PrimAmps(PrimAmp2m_15243)%Result &
                & + PrimAmps(PrimAmp2m_12543)%Result
 
           FermionLoopPartAmp(3,-2:1)=&               
                &  +PrimAmps(PrimAmp2m_13254)%Result(-2:1) 
 
           FermionLoopPartAmp(4,-2:1)=&
                & +PrimAmps(PrimAmp2m_14253)%Result(-2:1) 

 
       Col1Lf_ttbggH = 0d0
       Col1Lf_ttbggH(1,1) = 4d0 * Cf**2 * Nc - 2d0*Cf  !  56/3 = 64/3 -8/3
       Col1Lf_ttbggH(1,2) = -4d0*Cf                    ! -16/3 = -8/3 -8/3
       Col1Lf_ttbggH(1,3) = -2d0*Cf
       Col1Lf_ttbggH(1,4) = Col1Lf_ttbggH(1,3) 

       Col1Lf_ttbggH(2,2) = Col1Lf_ttbggH(1,1)
       Col1Lf_ttbggH(2,1) = Col1Lf_ttbggH(1,2)
       Col1Lf_ttbggH(2,3) = Col1Lf_ttbggH(1,3)
       Col1Lf_ttbggH(2,4) = Col1Lf_ttbggH(1,3)
 
       NLO_Res_Pol(-2:1) = (0d0,0d0)
       do jPrimAmp=1,4
          do iPrimAmp=1,2
           NLO_Res_Pol(-2:1) = NLO_Res_Pol(-2:1) + Col1Lf_ttbggH(iPrimAmp,jPrimAmp) * dreal( BornAmps(iPrimAmp)%Result*dconjg(FermionLoopPartAmp(jPrimAmp,-2:1)) )
          enddo
       enddo
       NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)
 
 
ENDIF
 
 
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.3 ) THEN
    call Gamma5Renorm_gg(MomExt(1:4,12:13),BornAmps,RenormAmps,Ren_Res_Pol)
    Ren_Res_UnPol=Ren_Res_UnPol + Ren_Res_Pol
 ENDIF
 
     enddo! helicity loop
 ENDIF


IF( Correction.EQ.0 ) THEN
   
! reweighting
   if (LO_Reweighting) then
      call ReWeight(LO_Res_UnPol,NBin,ReWeightHisto,LO_Res_UnPol_RW)
      LO_Res_UnPol = LO_Res_UnPol_RW
   endif


!  normalization
   
   
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * WidthExpansion
   EvalCS_1L_ttbggH = LO_Res_Unpol * PreFac


ELSEIF( Correction.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
!  CT contributions

IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.1 ) THEN

!!   print *, "before renorm"
!!   print *, "bosonic corrections"
!!   print *, "DP",NLO_Res_UnPol(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!!   print *, "fermionic corrections"
!!   print *, "DP",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol_Ferm(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol_Ferm(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol

                        ! beta        !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0 )*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu) contrib. from  top WFRC
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) - (-2d0/3d0*Nf_light)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol  ! shift alpha_s^DR --> alpha_s^MSbar
! Yukawa coupling renormalization
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) - 3d0*4d0/3d0*LO_Res_UnPol                                             ! -3/eps*CF*
   NLO_Res_UnPol(0)  = NLO_Res_UnPol(0)  + (-5d0-3d0*2d0*dlog(MuRen/m_top))*4d0/3d0*LO_Res_UnPol                ! -5*CF-3*CF*log(mu^2/mt^2)

!   if (nonrenormcoupl) then
!      NLO_Res_UnPol(-1)=NLO_Res_UnPol(-1)+HOO_Ren_Res_UnPol
!      NLO_Res_UnPol(0) =NLO_Res_UnPol(0) +dlog(MuRen**2/TTBZ_MassScale**2)*HOO_Ren_Res_UnPol
!   endif

!!   print *, "after renorm"
!!   print *, "bosonic corrections"
!!   print *, "DP",NLO_Res_UnPol(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!!   print *, "log term from Yukawa renorm", (-3d0*2d0*dlog(MuRen/m_top))*4d0/3d0*LO_Res_UnPol 
!!   print *, "fermionic corrections"
!!   print *, "DP",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol_Ferm(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol_Ferm(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!!   print *, "both corrections"
!!   print *, "DP",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol
!!   print *, "SP",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol
!!   print *, "cc",(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0))/LO_Res_UnPol
!!   print *, "rat",(NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol(0)+NLO_Res_UNpol(1)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!!
!!   stop

ENDIF
IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.3 ) THEN
! Raoul's gamma5 ren -- CF now included in subroutine
! NB -- additional factor of two relative to ttb+Z case -- see Larin, hep-ph/9302240; compare Eq. 10 for non-singlet axial-vector current, 
! Eq. 33 for singlet axial-vector current, and Eq. 15 for pseudoscalar current
   NLO_Res_UnPol( 0) = NLO_Res_UnPol(0)  - Ren_Res_UnPol*4d0   
ENDIF

!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2                            
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor 
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2 * alpha_sOver2Pi*RunFactor 
   EvalCS_1L_ttbggH = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN


!! comment this out
! NLO_Res_UnPol      = NLO_Res_UnPol      * PreFac
! NLO_Res_UnPol_Ferm = NLO_Res_UnPol_Ferm * PreFac
! LO_Res_Unpol       = LO_Res_Unpol       * PreFac

   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8)
   ENDIF

!! xe=0.3d0

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)

   call EvalIntDipoles_GGTTBGH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,HOp(1:3))
   HOp(1:3) = HOp(1:3)*RunFactor**3 * PreFac
   EvalCS_1L_ttbggH = HOp(1)    * pdf(0,1)  * pdf(0,2)   &
                    + HOp(2)/xE * pdf_z(0,1)* pdf(0,2)   &
                    + HOp(3)/xE * pdf(0,1)  * pdf_z(0,2)

! print *, "1L check",NLO_Res_UnPol(-2)                           /(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "1L check",( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "ID check", EvalCS_1L_ttbggH/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
!! pause
! stop



ENDIF





! !      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
! !      ./TOPAZ Collider=1 TopDK=0 HDK=0 Process=101 Correction=0 NLOParam=0 PDFSet=2 MTop=1.73 ObsSet=81
! !      MADGRAPH CHECK: gg->ttbH, mt=173, alpha_s=0.13 mH=125.00, vev=250.618249228543
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBH(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_H,vev
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!       pause
       



   if( IsNan(EvalCS_1L_ttbggH) ) then
        print *, "NAN:",EvalCS_1L_ttbggH
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac, sHatJacobi !, MHJacobi
        print *, eta1,eta2,MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbggH = 0d0
!         pause
        return
   endif

   do NHisto=1,NumHistograms
!      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggH,BinValue=PObs(NHisto))
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbggH)
   enddo
   EvalCounter = EvalCounter + 1 


   EvalCS_1L_ttbggH = EvalCS_1L_ttbggH/VgsWgt


RETURN
END FUNCTION










FUNCTION EvalCS_1L_ttbqqbH_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_1L_ttbqqbH_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_1L_ttbqqbH(yRnd,VgsWgt)
EvalCS_1L_ttbqqbH_MPI=0
RETURN
END FUNCTION






FUNCTION EvalCS_1L_ttbqqbH(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModUCuts
use ModUCuts_128
use ModIntegrals
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHDecay
use ModIntDipoles_QQBTTBGH
use ModIntDipoles_QGTTBQH
use ModIntDipoles_QBGTTBQBH
implicit none
real(8) ::  EvalCS_1L_ttbqqbH,yRnd(1:VegasMxDim),VgsWgt,xE
complex(8) :: rdiv(1:2),LO_Res_Pol,LO_Res_Unpol,NLO_Res_Pol(-2:1),NLO_Res_UnPol(-2:1),NLO_Res_Unpol_Ferm(-2:1),rdiv5(1:2)
complex(8) :: BosonicPartAmp(1:2,-2:1),FermionPartAmp(1:2,-2:1),mydummy(1:2),LOPartAmp(1:2)
integer :: iHel,iPrimAmp,jPrimAmp,tmphel,APrimAmp,ListPrimAmps(1:7),LocalSisters(1:7),npdfmax
real(8) :: EHat,RunFactor,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,AccPoles,HOp(1:2,1:3),pdf_z(-6:6,1:2)
real(8) :: MomExt(1:4,1:14),MomZl(1:4), MomZa(1:4),pHsq,MH_Inv,ZDKMatel,Msq_T_BWENU,p12Z(1:4),p34Z(1:4)
complex(8) :: s12Z,s34Z,fermloop_fin(1:2)
complex(8) :: propH,ZPolVec(1:4),BarSpi(1:4),Spi(1:4),light_quark_coupl(1:2)
logical :: applyPSCut,nonrenormcoupl
logical, save :: first=.true.
real(8) :: DPtol, QPtol,PObs(1:NumMaxHisto)
! real(8) :: couplZUU,couplZDD,couplZLL,couplGUU,couplGDD,couplGLL,couplZTT,couplGTT
real(8) :: MG_MOM(0:3,1:NumExtParticles)
real(8) :: MadGraph_tree
real(8),parameter :: Nc=3d0
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a(1:2),PDFFac_b(1:2),PDFFac(1:2),pdf(-6:6,1:2)
integer :: NHisto,NBin(1:NumMaxHisto),npdf,qf
integer,parameter :: up=1,dn=2
complex(8) :: couplZQQ_left_dyn_old,couplZQQ_right_dyn_old,couplZTT_left_dyn_old,couplZTT_right_dyn_old,RenormAmp(1:2),couplZQQ_left2_dyn_old,couplZQQ_right2_dyn_old,couplZTT_left2_dyn_old,couplZTT_right2_dyn_old,couplZTT_left2_old,couplZTT_right2_old,couplZTT_left_old,couplZTT_right_old
real(8) :: Ren_Res_Pol,Ren_Res_UnPol,R_V,R_A,prim_opp_err(1:10)
include 'misc/global_import'
include "vegas_common.f"
real(8) :: tmpVcoupl,tmpAcoupl,tmpV2coupl,tmpA2coupl
complex(8) :: tmpBornResults(14),RenormAmps(14),HOORenormPartAmp(1:2),HOO_RenormAmp(1:NumPrimAmps),HOO_RenormPartAmp(1:2),HOO_Ren_Res_Pol,HOO_Ren_Res_UnPol


  EvalCS_1L_ttbqqbH = 0d0
  DPtol=1d-4
  QPtol=1d-3
  opp_err=0d0
  prim_opp_err=0d0

  
  call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!   call SmearExternal(yRnd(18),M_H,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)
  MH_Inv = M_H
  if( EHat.le.2d0*m_Top+MH_Inv ) then
     EvalCS_1L_ttbqqbH = 0d0
     return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
  call EvalPhaseSpace_2to3M(EHat,MH_Inv,yRnd(3:7),MomExt(1:4,1:5),PSWgt)! q q G tb t
  call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))


  ISFac = MomCrossing(MomExt)   
  IF( TOPDECAYS.NE.0 ) THEN
     call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
     call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
     PSWgt = PSWgt * PSWgt2*PSWgt3
     call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
     call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))
  ENDIF

!   IF( HDECAYS.NE.0 .AND. HDECAYS.NE.-2 ) THEN
!      call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)
!      NRndHel=NRndHel+2
!      PSWgt = PSWgt * PSWgt4 * MHJacobi
!   ENDIF


   call Kinematics_TTBARH(0,MomExt(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      EvalCS_1L_ttbqqbH = 0d0
      return
   endif

! ! we still need this for the massless qqZ couplings
!    pHsq=MomExt(1,3)*MomExt(1,3)-MomExt(2,3)*MomExt(2,3)-MomExt(3,3)*MomExt(3,3)-MomExt(4,3)*MomExt(4,3)
!    if ( HDecays.lt.10 .and. HDecays.gt.0 ) then
!       propH = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)
!    elseif (HDecays .gt. 10) then 
!       propH=(1d0,0d0)/(pHsq-m_Z**2+ci*Ga_ZExp*m_Z)
!    endif
!    propH  = pHsq*propH! the pZ^2 will cancel against a 1/kZ^2 term in the Z decay matrix element
   
   call SetPropagators()
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
   PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
   PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
   PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
 
! !    ! this is for the Madgraph check
!    PDFFac_a(up) = 1d0
!    PDFFac_a(dn) = 0d0
!    PDFFac_b(up) = 0d0
!    PDFFac_b(dn) = 0d0
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)

   LO_Res_Unpol             = (0d0,0d0)
   NLO_Res_Unpol(-2:1)      = (0d0,0d0)
   NLO_Res_Unpol_Ferm(-2:1) = (0d0,0d0)
   Ren_Res_Pol=0d0
   Ren_Res_UnPol=0d0
   HOO_Ren_Res_UnPol=0d0
   HOO_Ren_Res_Pol=0d0

!------------ LO --------------
IF( CORRECTION.EQ.0 ) THEN
   do npdf=1,2
      if (npdf .eq. 1) then
         PDFFac(1:2)=PDFFac_a(1:2) 
      elseif (npdf .eq. 2) then
         PDFFac(1:2)=PDFFac_b(1:2)
         call swapMom(MomExt(1:4,1),MomExt(1:4,2))
         ISFac = MomCrossing(MomExt)
     endif
         
    ISFac = MomCrossing(MomExt)
    call SetPropagators()

   do iHel=1,NumHelicities
       call HelCrossing(Helicities(iHel,1:NumExtParticles))
       call SetPolarizations()
       if( ExtParticle(3)%Helicity*ExtParticle(4)%Helicity.eq.+1 ) cycle 
! RR added this line 04/04/2015 -- should never be true
       if( ExtParticle(3)%Helicity*ExtParticle(4)%Helicity.ne.-1 ) cycle 
!          if( HDecays.gt.0 ) then
!             if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be more elegantly done in mod_process
!             call HDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
! !             call ZGamLCoupl(1,Helicities(iHel,5),couplZLL,couplGLL)  ! charged lept
!          else
! !             couplZTT_left_dyn  = couplZTT_left
! !             couplZTT_right_dyn = couplZTT_right
!          endif

!        if (npdf .eq. 2) then! change helicities of the massless quarks for the couplings to Z          
!           Helicities(iHel,3)=-Helicities(iHel,3)
!           Helicities(iHel,4)=-Helicities(iHel,4)
!        endif
!        call ZGamQcoupl(Up_,Helicities(iHel,3),couplZUU,couplGUU)
!        call ZGamQcoupl(Dn_,Helicities(iHel,3),couplZDD,couplGDD)
!        couplZQQ_left_dyn=one
!        couplZQQ_right_dyn=one                  
       do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
       enddo

       LO_Res_Pol = (0d0,0d0)
!        if (Zdecays .eq. 0) then
         LOPartAmp(up)=BornAmps(1)%Result !+ BornAmps(2)%Result*couplZUU
         LOPartAmp(dn)=BornAmps(1)%Result !+ BornAmps(2)%Result*couplZDD
!        elseif( HDecays .ge. 1 .and. Zdecays.lt.10 ) then
!          LOPartAmp(up)=BornAmps(1)%Result + BornAmps(2)%Result*couplZUU*propH*couplZLL 
!          LOPartAmp(dn)=BornAmps(1)%Result + BornAmps(2)%Result*couplZDD*propH*couplZLL  
!        elseif( Zdecays.gt.10 ) then
!          LOPartAmp(up)s=BornAmps(1)%Result + BornAmps(2)%Result*( couplZUU*propH*couplZLL + couplGUU*propPh*couplGLL )
!          LOPartAmp(dn)=BornAmps(1)%Result + BornAmps(2)%Result*( couplZDD*propH*couplZLL + couplGDD*propPh*couplGLL )
!        endif
        LO_Res_Pol =  ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
        LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
     enddo!helicity loop
  enddo! npdf loop
  !   call swapMom(MomExt(1:4,1),MomExt(1:4,2))   ! swap back to original order, for ID below
  !print * , "remove swap"


!------------ 1 LOOP --------------
ELSEIF( CORRECTION.EQ.1 ) THEN
    npdfmax=2 
    do npdf=1,npdfmax
       if (npdf .eq. 1) then
          PDFFac(1:2)=PDFFac_a(1:2)
!          if( npdfmax.eq.1 ) PDFFac(1:2) = PDFFac(1:2) + PDFFac_b(1:2)
! !                                          PDFFac(1:2)= (/1d0,0d0 /)
       elseif (npdf .eq. 2) then
          PDFFac(1:2)=PDFFac_b(1:2)
          call swapMom(MomExt(1:4,1),MomExt(1:4,2))
          ISFac = MomCrossing(MomExt)
       endif
       ISFac = MomCrossing(MomExt)
       call SetPropagators()
! 
!       !   do iHel=nHel(1),nHel(2)
       do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          if ( Helicities(iHel,3)*Helicities(iHel,4) .eq. 1) cycle
! RR added this line 04/04/2015 -- should never be true
          if ( Helicities(iHel,3)*Helicities(iHel,4) .ne. -1) cycle
          call SetPolarizations()

!          if( HDecays.gt.0 ) then
!             if( ExtParticle(5)%Helicity.eq.0 ) cycle!   this can be done more elegantly in mod_process
!             call HDecay(ExtParticle(5),DK_LO,MomExt(1:4,12:13))
! !             call ZGamLCoupl(1,Helicities(iHel,5),couplZLL,couplGLL)  ! charged lept
!          else
!             couplZTT_left_dyn=couplZTT_left
!             couplZTT_right_dyn=couplZTT_right
!          endif
! 
! 
!        if (npdf .eq. 2) then
!           ! change helicities of the massless quarks for the couplings to Z          
!           Helicities(iHel,3)=-Helicities(iHel,3)
!           Helicities(iHel,4)=-Helicities(iHel,4)
!        endif
! !        call ZGamQcoupl(Up_,Helicities(iHel,3),couplZUU,couplGUU)
! !        call ZGamQcoupl(Dn_,Helicities(iHel,3),couplZDD,couplGDD)
!        couplZQQ_left_dyn=one
!        couplZQQ_right_dyn=one
! 
        do iPrimAmp=1,NumBornAmps
           call EvalTree(BornAmps(iPrimAmp))
        enddo
! 
! 
!         if (Zdecays .eq. 0) then
            LOPartAmp(up)=BornAmps(1)%Result!+BornAmps(2)%Result*( couplZUU )
            LOPartAmp(dn)=BornAmps(1)%Result!+BornAmps(2)%Result*( couplZDD )
! !         elseif( Zdecays.lt.10 .and. Zdecays .gt. 0 ) then
! !            LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*couplZUU*propH*couplZLL 
! !            LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*couplZDD*propH*couplZLL 
! !         elseif( Zdecays.gt.10 ) then
! !            LOPartAmp(up)=BornAmps(1)%Result + BornAmps(2)%Result*( couplZUU*propH*couplZLL + couplGUU*propPh*couplGLL )
! !            LOPartAmp(dn)=BornAmps(1)%Result + BornAmps(2)%Result*( couplZDD*propH*couplZLL + couplGDD*propPh*couplGLL )
! !         elseif (Zdecays .eq. -2)  then    ! this is the photon comparison
! !            LOPartAmp(up)=BornAmps(1)%Result+BornAmps(2)%Result*Q_up
! !            LOPartAmp(dn)=BornAmps(1)%Result+BornAmps(2)%Result*Q_dn
!         endif
! 
       LO_Res_Pol = (0d0,0d0)
       LO_Res_Pol =  ColLO_ttbqqb(1,1) * ( LOPartAmp(up)*dconjg(LOPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(LOPartAmp(dn))*PDFFac(dn))
       LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
! 
! !!!! RR Higher Order Operator -- i.e. sigma_{mu,nu} q^{nu} --  renorm term ??????
!       if (nonrenormcoupl) then
!          call SigmaRenorm_qqb(MomExt(1:4,12:13),BornAmps,HOO_RenormAmp,HOO_RenormPartAmp)
!          !  incl CF factor from vertex correction.
!          HOO_Ren_Res_Pol = 4d0/3d0*ColLO_ttbqqb(1,1) * dreal( LOPartAmp(up)*dconjg(HOO_RenormPartAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(HOO_RenormPartAmp(dn))*PDFFac(dn) )
!          HOO_Ren_Res_UnPol=HOO_Ren_Res_UnPol + HOO_Ren_Res_Pol
!       endif
! 
!       if (Zdecays .eq. 0) then
!          light_quark_coupl(up)= couplZUU 
!          light_quark_coupl(dn)= couplZDD 
! !       elseif( Zdecays.lt.10 .and. Zdecays .gt. 0 ) then
! !          light_quark_coupl(up)= couplZUU*propH*couplZLL 
! !          light_quark_coupl(dn)= couplZDD*propH*couplZLL 
! !       elseif( Zdecays.gt.10 ) then
! !          light_quark_coupl(up)= couplZUU*propH*couplZLL + couplGUU*propPh*couplGLL 
! !          light_quark_coupl(dn)= couplZDD*propH*couplZLL + couplGDD*propPh*couplGLL 
! !       elseif (Zdecays .eq. -2)  then    ! this is the photon comparison
! !          light_quark_coupl(up)= Q_up
! !          light_quark_coupl(dn)= Q_dn
!       endif
! 
       ListPrimAmps=(/1,2,3,5,6,8,9/)
       LocalSisters=(/0,0,1,0,1,0,0/)
! 
! 
! ! ------------ bosonic loops --------------------------------
! !------------------------------------------------------------
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.1 ) THEN
 
      do iPrimAmp=1,5
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
 
 ! only UV renorm those with virt gluons on the top line
          if (iPrimAmp .eq. 1 .or. iPrimAmp .eq. 2 .or. iPrimAmp .eq. 3  .or. iPrimAmp .eq. 5) then
             call RenormalizeUV(PrimAmps(iPrimAmp),BornAmps(iPrimAmp),MuRen**2)
          endif
          PrimAmps(iPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(iPrimAmp)%Result(-2:1)
          call OneLoopDiv(PrimAmps(iPrimAmp),MuRen**2,3,rdiv(2),rdiv(1)) 
          prim_opp_err(iPrimAmp)=opp_err
       enddo
 
 
 ! combine into real primitives      
       PrimAmps(PrimAmp3_15432)%Result=PrimAmps(PrimAmp3_15432)%Result+PrimAmps(PrimAmp3_14352)%Result
       BornAmps(PrimAmp3_15432)%Result=BornAmps(PrimAmp3_15432)%Result+BornAmps(PrimAmp3_14352)%Result
!       HOO_RenormAmp(PrimAmp3_15432)=HOO_RenormAmp(PrimAmp3_15432)+HOO_RenormAmp(PrimAmp3_14352)
       prim_opp_err(PrimAmp3_15432)=prim_opp_err(PrimAmp3_15432) + prim_opp_err(PrimAmp3_14352)
       
 
       do iPrimAmp=1,4
          APrimAmp=ListPrimAmps(iPrimAmp)
          call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
! 
! ! RR  -- now we modify rdiv for this process --- 
          if (APrimAmp .le. 2) then
             rdiv(1)=rdiv(1)+3d0 
          elseif (APrimAmp .eq. 3) then
             rdiv(1)=rdiv(1)+1.5d0 
          endif
!          if (nonrenormcoupl) then
!             if (APrimAmp .eq. 1 .or. APrimAmp .eq. 3 .or. APrimAmp .eq. 5) then
!                rdiv(1)=rdiv(1)-1d0/2d0*HOO_RenormAmp(APrimAmp)/BornAmps(APrimAmp)%Result
!             endif
!          endif
! 
          AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
          
          if ( AccPoles .gt. DPtol .or. prim_opp_err(APrimAmp) .gt. 1d-2) then
!          if (.true.) then
!             print *, "pole accuracy=",AccPoles
!             print *,"using QP"
             useQP=useQP+1
             PrimAmps(APrimAmp)%Result=(0d0,0d0)
 
             
             do jPrimAmp=APrimAmp,APrimAmp+LocalSisters(iPrimAmp)
                call SetKirill(PrimAmps(jPrimAmp))
                call PentCut_128(PrimAmps(jPrimAmp))
                call QuadCut_128(PrimAmps(jPrimAmp))
                call TripCut_128(PrimAmps(jPrimAmp))
                call DoubCut_128(PrimAmps(jPrimAmp))
                call SingCut_128(PrimAmps(jPrimAmp))
                call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
                if (jPrimAmp .eq. 1 .or. jPrimAmp .eq. 2 .or. jPrimAmp .eq. 3  .or. jPrimAmp .eq. 5) then
                   call RenormalizeUV(PrimAmps(jPrimAmp),BornAmps(jPrimAmp),MuRen**2)
                endif
                PrimAmps(jPrimAmp)%Result(-2:1) = (0d0,1d0) * PrimAmps(jPrimAmp)%Result(-2:1)               
             enddo
 
 
             do jPrimAmp=1,LocalSisters(iPrimAmp)
                PrimAmps(APrimAmp)%Result(-2:1)=PrimAmps(APrimAmp)%Result(-2:1)+PrimAmps(APrimAmp+jPrimAmp)%Result(-2:1)
             enddo
 
             AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
             prim_opp_err(APrimAmp)=opp_err
             if ( AccPoles .gt. QPtol .or. prim_opp_err(APrimAmp) .gt. 1d-2) then
                print *, 'QP fails: ', AccPoles, prim_opp_err(APrimAmp), APrimAmp
                PrimAmps(APrimAmp)%Result=0d0
                pole_skipped=pole_skipped+1               
                SkipCounter = SkipCounter + 1
                RETURN ! reject the whole event instead of just this primamp
              endif
           endif
        enddo
        
 
       BosonicPartAmp(up,-2:1)=  &
                              +   Nc *  PrimAmps(PrimAmp1_15234)%Result(-2:1) &
                              - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1)+ PrimAmps(PrimAmp1_15243)%Result(-2:1) ) &
                              - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1) &
                              & + PrimAmps(PrimAmp4_15234)%Result(-2:1))

       BosonicPartAmp(dn,-2:1)=  &
                              +   Nc *  PrimAmps(PrimAmp1_15234)%Result(-2:1) &
                              - 2d0/Nc*( PrimAmps(PrimAmp1_15234)%Result(-2:1) + PrimAmps(PrimAmp1_15243)%Result(-2:1) ) &
                              - 1d0/Nc*( PrimAmps(PrimAmp3_15432)%Result(-2:1)  + PrimAmps(PrimAmp4_15234)%Result(-2:1) )

       NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(BosonicPartAmp(up,-2:1))) *PDFFac(up) &
            + dreal(LOPartAmp(dn)*dconjg(BosonicPartAmp(dn,-2:1)))*PDFFac(dn) )
 
       NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1) + NLO_Res_Pol(-2:1)
! 
! 
    ENDIF
! 
! ! ------------ fermionic loops ------------------------------
! !------------------------------------------------------------
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.2 ) THEN
! 

! ! this provides an analytic form for Prim 8
       p12Z(1:4)=ExtParticle(1)%Mom(1:4)+ExtParticle(2)%Mom(1:4)+ExtParticle(5)%Mom(1:4)
       s12Z=p12Z(1)**2-p12Z(2)**2-p12Z(3)**2-p12Z(4)**2
       fermloop_fin(1) = -2d0/3d0*log(-MuRen**2/s12Z)-10d0/9d0

 
       do iPrimAmp=6,9
          if (iPrimAmp .eq. 8) then
             PrimAmps(iPrimAmp)%Result(-1)=-2d0/3d0 * BornAmps(1)%Result
             PrimAmps(iPrimAmp)%Result(0)=fermloop_fin(1) * BornAmps(1)%Result
             PrimAmps(iPrimAmp)%Result(1)= (0d0,0d0)
             cycle
          endif
          call SetKirill(PrimAmps(iPrimAmp))
          call PentCut(PrimAmps(iPrimAmp))
          call QuadCut(PrimAmps(iPrimAmp))
          call TripCut(PrimAmps(iPrimAmp))
          call DoubCut(PrimAmps(iPrimAmp))
          call SingCut(PrimAmps(iPrimAmp))
          call EvalMasterIntegrals(PrimAmps(iPrimAmp),MuRen**2)
          PrimAmps(iPrimAmp)%Result(-2:1) = -(0d0,1d0)*PrimAmps(iPrimAmp)%Result(-2:1) !minus is from closed fermion loop
       enddo
 
! 
       PrimAmps(PrimAmp2m_12534)%Result = PrimAmps(PrimAmp2m_12534)%Result + PrimAmps(PrimAmp2m_12345)%Result
! hack to get correct Born Amp to compare with heavy loop with H on loop
       BornAmps(PrimAmp2m_12534)%Result = BornAmps(1)%Result
       do iPrimAmp=5,7
          APrimAmp=ListPrimAmps(iPrimAmp)
          call OneLoopDiv(PrimAmps(APrimAmp),MuRen**2,3,rdiv(2),rdiv(1))
          rdiv(1)=rdiv(1)+1.0d0/3.0d0         
          ! this is a hack to give me the "correct pole for the Z on the top loop
          ! ideally, of course, this will be correctly calculated by oneloopdiv, and we can remove these lines...
          if (iPrimAmp .eq. 5) then
             rdiv=0d0
          endif
 
 
          AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
 
          if (AccPoles .gt. DPtol) then
!          if (.true.) then
!             print *, "pole accuracy=",AccPoles
!             print *,"using QP"
             useQP=useQP+1
             PrimAmps(APrimAmp)%Result=(0d0,0d0)
 
             do jPrimAmp=APrimAmp,APrimAmp+LocalSisters(iPrimAmp)
                call SetKirill(PrimAmps(jPrimAmp))
                call PentCut_128(PrimAmps(jPrimAmp))
                call QuadCut_128(PrimAmps(jPrimAmp))
                call TripCut_128(PrimAmps(jPrimAmp))
                call DoubCut_128(PrimAmps(jPrimAmp))
                call SingCut_128(PrimAmps(jPrimAmp))
                call EvalMasterIntegrals(PrimAmps(jPrimAmp),MuRen**2)
!                if (jPrimAmp .eq. 1 .or. jPrimAmp .eq. 3 .or. jPrimAmp .eq. 5  .or. jPrimAmp .eq. 10) then
!                   call RenormalizeUV(PrimAmps(jPrimAmp),BornAmps(jPrimAmp),MuRen**2)
!                endif
                PrimAmps(jPrimAmp)%Result(-2:1) = -(0d0,1d0) * PrimAmps(jPrimAmp)%Result(-2:1)               
             enddo
 
 
             do jPrimAmp=1,LocalSisters(iPrimAmp)
                PrimAmps(APrimAmp)%Result(-2:1)=PrimAmps(APrimAmp)%Result(-2:1)+PrimAmps(APrimAmp+jPrimAmp)%Result(-2:1)
             enddo
 
             AccPoles = CheckPoles(PrimAmps(APrimAmp),BornAmps(APrimAmp),rdiv(1:2))
             if ( AccPoles .gt. QPtol ) then
                print *, 'QP fails: ', AccPoles,APrimAmp
                PrimAmps(APrimAmp)%Result=0d0
                pole_skipped=pole_skipped+1
                SkipCounter = SkipCounter + 1
                RETURN ! reject the whole event instead of just this primamp
             endif
          endif! QP         
       enddo
 

       FermionPartAmp(up,-2:1) = Nf_light*PrimAmps(PrimAmp2_15234)%Result(-2:1) &
            + PrimAmps(PrimAmp2m_15234)%Result(-2:1)  &
            + PrimAmps(PrimAmp2m_12534)%Result(-2:1) 
 
       FermionPartAmp(dn,-2:1) = Nf_light*PrimAmps(PrimAmp2_15234)%Result(-2:1) &
            + PrimAmps(PrimAmp2m_15234)%Result(-2:1) &
            + PrimAmps(PrimAmp2m_12534)%Result(-2:1) 

!       
!       if (HDecays .gt. 0) then
!          PrimAmps(PrimAmp2_12534)%Result(:) = PrimAmps(PrimAmp2_12534)%Result(:) * propH * couplZLL
! !          if (ZQcoupl .ne. 3) then
! !             call Error("WARNING : Fermion loop with decaying Z on light quark loop -- propagator only calculated using AXIAL VECTOR coupling")
! !          endif
!       endif
! 

      
      NLO_Res_Pol(-2:1) = Col1L_ttbqqb(1,1) *( dreal(LOPartAmp(up)*dconjg(FermionPartAmp(up,-2:1)))*PDFFac(up) &
                                             + dreal(LOPartAmp(dn)*dconjg(FermionPartAmp(dn,-2:1)))*PDFFac(dn) )
      NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) + NLO_Res_Pol(-2:1)
 
 ENDIF
 
 
 
 
! !     GAMMA-5 RENORMALIZATION -- SEE RR NOTES, 23 AUG 2013                        !
 IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.3 ) THEN
    call Gamma5Renorm_qqb(MomExt(1:4,12:13),BornAmps,RenormAmp)
! 4/3=CF is vertex color factor
    Ren_Res_Pol = 4d0/3d0*ColLO_ttbqqb(1,1) * dreal( LOPartAmp(up)*dconjg(RenormAmp(up))*PDFFac(up) + LOPartAmp(dn)*dconjg(RenormAmp(dn))*PDFFac(dn) )     
    Ren_Res_UnPol=Ren_Res_UnPol + Ren_Res_Pol
ENDIF

! 
! 
enddo!helicity loop
  enddo ! npdf
ENDIF


IF( CORRECTION.EQ.0 ) THEN
!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * WidthExpansion
   EvalCS_1L_ttbqqbH = LO_Res_Unpol * PreFac

ELSEIF( CORRECTION.EQ.1 ) THEN
!  overall normalization: (4*Pi)^eps/Gamma(1-eps)
IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.1 ) THEN
!  CT contributions                           ! beta           !top WFRC
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) + (-11d0/3d0*3d0 - 3d0*4d0/3d0)*LO_Res_Unpol
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-3d0*4d0/3d0)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from  top WFRC
   NLO_Res_UnPol_Ferm(-1) = NLO_Res_UnPol_Ferm(-1) + (+2d0/3d0*Nf_light+2d0/3d0*Nf_heavy)*LO_Res_Unpol
   NLO_Res_UnPol_Ferm( 0) = NLO_Res_UnPol_Ferm( 0) + (2d0/3d0*Nf_heavy)*2d0*dlog(MuRen/m_top)*LO_Res_Unpol  ! finite log(mu2) contrib. from heavy flavor in alpha_s ren.
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + (-5d0/2d0*8d0/3d0 )*LO_Res_Unpol   ! finite contribution from top WFRC's
   NLO_Res_UnPol( 0) = NLO_Res_UnPol( 0) + LO_Res_Unpol ! shift alpha_s^DR --> alpha_s^MSbar
! Yukawa coupling renormalization                                                                                                                                             
   NLO_Res_UnPol(-1) = NLO_Res_UnPol(-1) - 3d0*4d0/3d0*LO_Res_UnPol                                             ! -3/eps*CF*                                                   
   NLO_Res_UnPol(0)  = NLO_Res_UnPol(0)  + (-5d0-3d0*2d0*dlog(MuRen/m_top))*4d0/3d0*LO_Res_UnPol                ! -5*CF-3*CF*log(mu^2/mt^2)     

!!   print *, "after renorm"
!!   print *, "LO", LO_Res_UnPol
!!   print *, "DP",NLO_Res_UnPol(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol(0)+NLO_Res_UnPol(1))/LO_Res_UnPol
!!   print *, "log term from Yukawa renorm", (-3d0*2d0*dlog(MuRen/m_top))*4d0/3d0*LO_Res_UnPol
!!   print *, "fermionic corrections"
!!   print *, "DP",NLO_Res_UnPol_Ferm(-2)/LO_Res_UnPol
!!   print *, "SP",NLO_Res_UnPol_Ferm(-1)/LO_Res_UnPol
!!   print *, "cc",NLO_Res_UnPol_Ferm(0)/LO_Res_UnPol
!!   print *, "rat",NLO_Res_UnPol_Ferm(1)/LO_Res_UnPol
!!   print *, "fin",(NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol
!!   print *, "both corrections"
!!   print *, "DP",(NLO_Res_UnPol(-2)+NLO_Res_UnPol_Ferm(-2))/LO_Res_UnPol
!!   print *, "SP",(NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1))/LO_Res_UnPol                                                                                                  
!!   print *, "cc",(NLO_Res_UnPol(0)+NLO_Res_UnPol_Ferm(0))/LO_Res_UnPol
!!   print *, "rat",(NLO_Res_UnPol(1)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol 
!!   print *, "fin",(NLO_Res_UnPol(0)+NLO_Res_UNpol(1)+NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1))/LO_Res_UnPol                                                               
!!   stop         

ENDIF

IF( TTBZ_DebugSwitch.EQ.0 .OR. TTBZ_DebugSwitch.EQ.3 ) THEN
! Raoul's gamma5 ren
! NB -- additional factor of two relative to ttb+Z case -- see Larin, hep-ph/9302240; compare Eq. 10 for non-singlet axial-vector current, 
! Eq. 33 for singlet axial-vector current, and Eq. 15 for pseudoscalar current
   NLO_Res_UnPol( 0)=NLO_Res_UnPol(0) - Ren_Res_UnPol*4d0  
ENDIF

!   if (nonrenormcoupl) then
!      print *, "incl counterterms for non-renorm sigma couplings"
!      print *, " NOTE: double check finite counterterms !!!"
!      NLO_Res_UnPol(-1)=NLO_Res_UnPol(-1)+HOO_Ren_Res_UnPol
!      NLO_Res_UnPol(0)=NLO_Res_UnPol(0)+dlog(MuRen**2/TTBZ_MassScale)*HOO_Ren_Res_UnPol
!endif



!  normalization
   LO_Res_Unpol = LO_Res_Unpol                         * ISFac * (alpha_s4Pi*RunFactor)**2 
   NLO_Res_UnPol(-2:1) = NLO_Res_UnPol(-2:1)           * ISFac * (alpha_s4Pi*RunFactor)**2  * alpha_sOver2Pi*RunFactor
   NLO_Res_UnPol_Ferm(-2:1) = NLO_Res_UnPol_Ferm(-2:1) * ISFac * (alpha_s4Pi*RunFactor)**2  * alpha_sOver2Pi*RunFactor
   EvalCS_1L_ttbqqbH = ( NLO_Res_UnPol(0)+NLO_Res_UnPol(1) + NLO_Res_UnPol_Ferm(0)+NLO_Res_UnPol_Ferm(1) ) * PreFac

ELSEIF( CORRECTION.EQ.3 ) THEN

! NLO_Res_UnPol      = NLO_Res_UnPol      * PreFac
! NLO_Res_UnPol_Ferm = NLO_Res_UnPol_Ferm * PreFac
! LO_Res_Unpol       = LO_Res_Unpol       * PreFac
! call swapMom(MomExt(1:4,1),MomExt(1:4,2))


   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   IF( TOPDECAYS.GE.1 ) THEN
       xE = yRnd(16)
   ELSEIF( TOPDECAYS.EQ.0 ) THEN
       xE = yRnd(8)
   ENDIF
   
!! xe=0.3d0!; print *, "fixed xE"

   call setPDFs(eta1/xE,eta2/xE,MuFac,pdf_z)
   IF( PROCESS.EQ.106 ) THEN

      npdf=1
      call EvalIntDipoles_QQBTTBGH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH= HOp(1,1)    * (pdf(Up_,1)*pdf(AUp_,2)+pdf(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(AUp_,2)+pdf_z(Chm_,1)*pdf(AChm_,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(AUp_,2)+pdf(Chm_,1)*pdf_z(AChm_,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(ADn_,2)+pdf(Str_,1)*pdf(AStr_,2)+pdf(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(ADn_,2)+pdf_z(Str_,1)*pdf(AStr_,2)+pdf_z(Bot_,1)*pdf(ABot_,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(ADn_,2)+pdf(Str_,1)*pdf_z(AStr_,2)+pdf(Bot_,1)*pdf_z(ABot_,2) )

      npdf=2
      call EvalIntDipoles_QQBTTBGH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH= EvalCS_1L_ttbqqbH          &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(AUp_,1)+pdf(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(AUp_,1)+pdf_z(Chm_,2)*pdf(AChm_,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(AUp_,1)+pdf(Chm_,2)*pdf_z(AChm_,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(ADn_,1)+pdf(Str_,2)*pdf(AStr_,1)+pdf(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(ADn_,1)+pdf_z(Str_,2)*pdf(AStr_,1)+pdf_z(Bot_,2)*pdf(ABot_,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(ADn_,1)+pdf(Str_,2)*pdf_z(AStr_,1)+pdf(Bot_,2)*pdf_z(ABot_,1) )

! print *, "1L check",NLO_Res_UnPol(-2)                           /(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "1L check",( NLO_Res_UnPol(-1)+NLO_Res_UnPol_Ferm(-1) )/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
! print *, "ID check", xe,EvalCS_1L_ttbqqbH/(alpha_sOver2Pi*RunFactor)/(LO_Res_Unpol)
!
! stop


   ELSEIF( PROCESS.EQ.103 ) THEN


! for int dipole check -- REMOVE!
!     PreFac=1d0

      npdf=1
      call EvalIntDipoles_QGTTBQH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH=HOp(1,1)    *  (pdf(Up_,1)*pdf(0,2)+pdf(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,1)*pdf(0,2)+pdf_z(Chm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(Up_,1)*pdf_z(0,2)+pdf(Chm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(Dn_,1)*pdf(0,2)+pdf(Str_,1)*pdf(0,2)+pdf(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,1)*pdf(0,2)+pdf_z(Str_,1)*pdf(0,2)+pdf_z(Bot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,1)*pdf_z(0,2)+pdf(Str_,1)*pdf_z(0,2)+pdf(Bot_,1)*pdf_z(0,2) )

      npdf=2
      call EvalIntDipoles_QGTTBQH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH= EvalCS_1L_ttbqqbH  &
                        +HOp(1,1)    * (pdf(Up_,2)*pdf(0,1)+pdf(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(Up_,2)*pdf(0,1)+pdf_z(Chm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(Up_,2)*pdf_z(0,1)+pdf(Chm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(Dn_,2)*pdf(0,1)+pdf(Str_,2)*pdf(0,1)+pdf(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(Dn_,2)*pdf(0,1)+pdf_z(Str_,2)*pdf(0,1)+pdf_z(Bot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(Dn_,2)*pdf_z(0,1)+pdf(Str_,2)*pdf_z(0,1)+pdf(Bot_,2)*pdf_z(0,1) )




   ELSEIF( PROCESS.EQ.104 ) THEN

! for int dipole check -- REMOVE!
!      PreFac=1d0

      npdf=1
      call EvalIntDipoles_QBGTTBQBH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH= HOp(1,1)    * (pdf(AUp_,1)*pdf(0,2)+pdf(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,1)*pdf(0,2)+pdf_z(AChm_,1)*pdf(0,2) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,1)*pdf_z(0,2)+pdf(AChm_,1)*pdf_z(0,2) ) &
                        +HOp(2,1)    * (pdf(ADn_,1)*pdf(0,2)+pdf(AStr_,1)*pdf(0,2)+pdf(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,1)*pdf(0,2)+pdf_z(AStr_,1)*pdf(0,2)+pdf_z(ABot_,1)*pdf(0,2) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,1)*pdf_z(0,2)+pdf(AStr_,1)*pdf_z(0,2)+pdf(ABot_,1)*pdf_z(0,2) )

      npdf=2
      call EvalIntDipoles_QBGTTBQBH((/MomExt(1:4,4),MomExt(1:4,3),MomExt(1:4,5),-MomExt(1:4,1),-MomExt(1:4,2)/),MomExt(1:4,6:13),xE,npdf,HOp(1:2,1:3))
      HOp(1:2,1:3) = HOp(1:2,1:3)*RunFactor**3 * PreFac
      EvalCS_1L_ttbqqbH= EvalCS_1L_ttbqqbH &
                        +HOp(1,1)    * (pdf(AUp_,2)*pdf(0,1)+pdf(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,2)/xE * (pdf_z(AUp_,2)*pdf(0,1)+pdf_z(AChm_,2)*pdf(0,1) ) &
                        +HOp(1,3)/xE * (pdf(AUp_,2)*pdf_z(0,1)+pdf(AChm_,2)*pdf_z(0,1) ) &
                        +HOp(2,1)    * (pdf(ADn_,2)*pdf(0,1)+pdf(AStr_,2)*pdf(0,1)+pdf(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,2)/xE * (pdf_z(ADn_,2)*pdf(0,1)+pdf_z(AStr_,2)*pdf(0,1)+pdf_z(ABot_,2)*pdf(0,1) ) &
                        +HOp(2,3)/xE * (pdf(ADn_,2)*pdf_z(0,1)+pdf(AStr_,2)*pdf_z(0,1)+pdf(ABot_,2)*pdf_z(0,1) )


! print *, "ID check", xe,EvalCS_1L_ttbqqbH/(alpha_sOver2Pi*RunFactor)
!STOP

   ENDIF

ENDIF



! !      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
! !      ./TOPAZ Collider=1 TopDK=0 HDK=0 Process=102 Correction=0 NLOParam=0 PDFSet=2 MTop=1.73 ObsSet=81
! !      MADGRAPH CHECK: gg->ttbH, mt=173, alpha_s=0.13 mH=125.00, vev=250.618249228543
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SUUB_TTBH(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_H,vev
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**2
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**2)
!       pause
       




   if( IsNan(EvalCS_1L_ttbqqbH) ) then
        print *, "NAN:",EvalCS_1L_ttbqqbH
        print *, yRnd(:)
        print *, NLO_Res_UnPol(0),NLO_Res_UnPol(1),NLO_Res_UnPol_Ferm(0),NLO_Res_UnPol_Ferm(1)
        print *, PSWgt , VgsWgt , PDFFac_a,PDFFac_b, sHatJacobi
        print *, eta1,eta2,((1d0-eta1)*xE+eta1),((1d0-eta2)*xE+eta2),MuFac,EHat
        print *, "Mom"
        print *, MomExt(1:4,:)
        print *, "SKIP EVENT!!!!!"
        EvalCS_1L_ttbqqbH = 0d0
        return
   endif
    EvalCounter = EvalCounter + 1 


   do NHisto=1,NumHistograms
!      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqbH,BinValue=PObs(NHisto))
      call intoHisto(NHisto,NBin(NHisto),EvalCS_1L_ttbqqbH)
   enddo



   EvalCS_1L_ttbqqbH = EvalCS_1L_ttbqqbH/VgsWgt

RETURN
END FUNCTION











FUNCTION EvalCS_Real_ttbgggH_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_Real_ttbgggH_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_Real_ttbgggH(yRnd,VgsWgt)
EvalCS_Real_ttbgggH_MPI=0
RETURN
END FUNCTION





FUNCTION EvalCS_Real_ttbgggH(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModMisc
use ModProcess
use ModDipoles_GGTTBGH
use ModHDecay
implicit none
real(8) ::  EvalCS_Real_ttbgggH,yRnd(1:VegasMxDim),VgsWgt,DipoleResult
complex(8) :: LO_Res_Pol,LO_Res_Unpol,PartAmp(1:4),PolExt3(1:4),The2ndResult(1:6)
integer :: iHel,jPrimAmp,iPrimAmp,NHisto,NBin(1:NumMaxHisto)
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,RunFactor,PreFac,sij
real(8) :: eta1,eta2,MHJacobi,sHatJacobi,FluxFac,PDFFac,MH_Inv
real(8) :: MomExt(1:4,1:14),pdf(-6:6,1:2),PObs(1:NumMaxHisto)
real(8) :: MG_MOM(0:3,1:6)
real(8) :: MadGraph_tree
logical :: applyPSCut,applySingCut
include "vegas_common.f"

yRnd( 1)=  0.3385585941088194d0
yRnd( 2)=  0.2799513116385563d0
yRnd( 3)=  0.012473622342792d0
yRnd( 4)=  0.2879364093709448d0
yRnd( 5)=  0.1334328211068331d0
yRnd( 6)=  0.7829718273519412d0
yRnd( 7)=  0.3479862101366653d0
yRnd( 8)=  0.1332233664734401d0
yRnd( 9)=  0.2332185946559626d0
yRnd(10)=  0.7471774192280964d0
print *, "fixing yrnds"
print *, yRnd


   EvalCS_Real_ttbgggH= 0d0
   DipoleResult = 0d0
   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    call SmearExternal(yRnd(21),M_H,Ga_H,Zero,EHat,MH_Inv,MHJacobi)
   MH_Inv = M_H
   if( EHat.le.2d0*m_Top+MH_Inv) then
      EvalCS_Real_ttbgggH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)

   call EvalPhaseSpace_2to4M(EHat,(/0d0,MH_Inv,m_Top,m_Top/),yRnd(3:10),MomExt(1:4,1:6),PSWgt)! gluon gluon gluon H tb t
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))
   ISFac = MomCrossing(MomExt)


   PSWgt2 = 1d0
   PSWgt3 = 1d0
   PSWgt4 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomExt(1:4,7:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
   call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,7:9))
   call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
ENDIF
! IF( HDECAYS.NE.0 ) THEN
!    call EvalPhasespace_HDecay(MomExt(1:4,4),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 * MHJacobi
! ENDIF

   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbgggH = 0d0
       SkipCounter = SkipCounter + 1
       return
   endif

   call Kinematics_TTBARH(1,MomExt(1:4,1:14),(/5,6,4,1,2,3,7,8,9,10,11,12,13,14/),applyPSCut,NBin,PObs)
   call setPDFs(eta1,eta2,MuFac,pdf)
   PDFFac = pdf(0,1) * pdf(0,2)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt * PDFFac
   RunFactor = RunAlphaS(NLOParam,MuRen)
   if( applyPSCut ) then
       EvalCS_Real_ttbgggH = 0d0
   else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
          call SetPolarizations(3,PolExt3(1:4))
!           if( HDecays.ne.0 ) then
!               call HDecay(ExtParticle(6),DK_LO,MomExt(1:4,13:14))
!           endif
          do iPrimAmp=1,NumBornAmps
                 call EvalTree(BornAmps(iPrimAmp),PolExt3(1:4),The2ndResult(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,6
          do iPrimAmp=1,6
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result * dconjg(BornAmps(jPrimAmp)%Result)
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbggg(iPrimAmp,jPrimAmp) * The2ndResult(iPrimAmp)    * dconjg(The2ndResult(jPrimAmp))
          enddo
          enddo
          LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 
        EvalCS_Real_ttbgggH = LO_Res_Unpol * PreFac
        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),EvalCS_Real_ttbgggH,BinValue=PObs(NHisto))
        enddo

        EvalCounter = EvalCounter + 1
endif!applyPSCut

     PreFac = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3  /PSWgt2/PSWgt3/PSWgt4
     call EvalDipoles_GGTTBGH((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:20),PreFac,DipoleResult)

     
     
     sij = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))                            ! collinear
!       sij = MomExt(1,3)**2                                                ! soft
     print *,  sij/EHat**2,EvalCS_Real_ttbgggH,DipoleResult,(1d0+EvalCS_Real_ttbgggH/DipoleResult)
     pause




!      careful, alphas=0.13 only for NLOParam=0 PDFSet=2
!      ./TOPAZ Collider=1 TopDK=0 HDK=0 Process=105 Correction=2 NLOParam=0 PDFSet=2 MTop=1.73 ObsSet=81
!      MADGRAPH CHECK: gg->ttbHg, mt=173, alpha_s=0.13 mH=125.00, vev=250.618249228543
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!       call coupsm(0)
!       call SGG_TTBHG(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_h
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**4
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**4)
!       pause

       

    EvalCS_Real_ttbgggH = (EvalCS_Real_ttbgggH + DipoleResult)/VgsWgt
RETURN
END FUNCTION








FUNCTION EvalCS_Real_ttbqqbgH_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_Real_ttbqqbgH_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_Real_ttbqqbgH(yRnd,VgsWgt)
EvalCS_Real_ttbqqbgH_MPI=0
RETURN
END FUNCTION



FUNCTION EvalCS_Real_ttbqqbgH(yRnd,VgsWgt)
use ModParameters
use ModKinematics
use ModAmplitudes
use ModProcess
use ModMisc
use ModHDecay
use ModDipoles_QQBTTBGH
use ModDipoles_QGTTBQH
use ModDipoles_QBGTTBQBH
implicit none
real(8) ::  EvalCS_Real_ttbqqbgH,EvalCS_Dips_ttbqqbgH,yRnd(1:VegasMxDim),VgsWgt,DipoleResult(1:2)
complex(8) :: LO_Res_Pol,LO_Res_Unpol,PropH,LOPartAmp(1:4,1:2),PolExt3(1:4),The2ndResult(1:6)
integer :: iHel,jPrimAmp,iPrimAmp,NHisto,NBin(1:NumMaxHisto),NPDF
real(8) :: EHat,PSWgt,PSWgt2,PSWgt3,PSWgt4,ISFac,RunFactor,PreFac,PreFacDip
real(8) :: eta1,eta2,MHJacobi,sHatJacobi,FluxFac,PDFFac(1:2),PDFFac_a(1:2),PDFFac_b(1:2),PObs(1:NumMaxHisto)
real(8) :: MomExt(1:4,1:14),sij,pdf(-6:6,1:2),MH_Inv,couplZLL,couplGLL,couplZUU,couplGUU,couplZDD,couplGDD
real(8) :: MG_MOM(0:3,1:6),MadGraph_tree
logical :: applyPSCut,applySingCut
include "vegas_common.f"


! yrnd(1)=0.1d0
! yrnd(2)=0.03d0
! yrnd(3:10)=0.7d0
! print *, "fixing yrnd"

   EvalCS_Real_ttbqqbgH= 0d0
   EvalCS_Dips_ttbqqbgH= 0d0

   call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!    call SmearExternal(yRnd(21),M_H,Ga_H,Zero,EHat,MH_Inv,MHJacobi)
   MH_inv = M_H
   if( EHat.le.2d0*m_Top+MH_Inv ) then
      EvalCS_Real_ttbqqbgH = 0d0
      return
   endif
   FluxFac = 1d0/(2d0*EHat**2)
   call EvalPhaseSpace_2to4M(EHat,(/0d0,MH_Inv,m_Top,m_Top/),yRnd(3:10),MomExt(1:4,1:6),PSWgt)! q qb gluon H tb t / q g q H tb t
   call boost2Lab(eta1,eta2,6,MomExt(1:4,1:6))

   PSWgt2 = 1d0
   PSWgt3 = 1d0
   PSWgt4 = 1d0
IF( TopDecays.GE.1 ) THEN
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(11:14),.false.,MomExt(1:4,7:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,6),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
   PSWgt = PSWgt * PSWgt2*PSWgt3
ENDIF
! IF( HDECAYS.GT.0 ) THEN
!    call EvalPhasespace_HDecay(MomExt(1:4,4),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
!    PSWgt = PSWgt * PSWgt4 * MHJacobi
! ENDIF


   call CheckSing(MomExt,applySingCut)
   if( applySingCut ) then
       EvalCS_Real_ttbqqbgH = 0d0
       SkipCounter = SkipCounter + 1
       return
   endif

   call Kinematics_TTBARH(1,MomExt(1:4,1:14),(/5,6,4,1,2,3,7,8,9,10,11,12,13,14/),applyPSCut,NBin,PObs)
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt * VgsWgt
   RunFactor = RunAlphaS(NLOParam,MuRen)


   call setPDFs(eta1,eta2,MuFac,pdf)
   IF( PROCESS.EQ.106 ) THEN
      PDFFac_a(Up_) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(Dn_) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(Up_) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(Dn_) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
   ELSEIF( PROCESS.EQ.103 ) THEN
      PDFFac_a(Up_) = pdf(Up_,1)*pdf(0,2) + pdf(Chm_,1)*pdf(0,2)
      PDFFac_a(Dn_) = pdf(Dn_,1)*pdf(0,2) + pdf(Str_,1)*pdf(0,2) + pdf(Bot_,1)*pdf(0,2)
      PDFFac_b(Up_) = pdf(Up_,2)*pdf(0,1) + pdf(Chm_,2)*pdf(0,1)
      PDFFac_b(Dn_) = pdf(Dn_,2)*pdf(0,1) + pdf(Str_,2)*pdf(0,1) + pdf(Bot_,2)*pdf(0,1)
   ELSEIF( PROCESS.EQ.104 ) THEN
      PDFFac_a(Up_) = pdf(AUp_,1)*pdf(0,2) + pdf(AChm_,1)*pdf(0,2)
      PDFFac_a(Dn_) = pdf(ADn_,1)*pdf(0,2) + pdf(AStr_,1)*pdf(0,2) + pdf(ABot_,1)*pdf(0,2)
      PDFFac_b(Up_) = pdf(AUp_,2)*pdf(0,1) + pdf(AChm_,2)*pdf(0,1)
      PDFFac_b(Dn_) = pdf(ADn_,2)*pdf(0,1) + pdf(AStr_,2)*pdf(0,1) + pdf(ABot_,2)*pdf(0,1)
   ENDIF

   ! for Madgraph check
!    PDFFac_a(Up_) = 0d0; print *, "disabled pdfs"
!    PDFFac_a(Dn_) = 0d0        
!    PDFFac_b(Up_) = 1d0
!    PDFFac_b(Dn_) = 0d0


   DO NPDF=1,2
        if(npdf.eq.1) then
              PDFFac(1:2) = PDFFac_a(1:2)
        elseif(npdf.eq.2) then
              PDFFac(1:2) = PDFFac_b(1:2)
              call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        endif

        ISFac = MomCrossing(MomExt)
        IF( TopDecays.GE.1 ) THEN
          call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,7:9))
          call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
        ENDIF

if( applyPSCut ) then
        EvalCS_Real_ttbqqbgH = 0d0
else
        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities
          call HelCrossing(Helicities(iHel,1:NumExtParticles))
!           call SetPolarizations()
          call SetPolarizations(3,PolExt3(1:4))
!           if( Process.eq.106 .and. ExtParticle(3)%Helicity.eq.ExtParticle(4)%Helicity ) cycle
!           if( HDecays.gt.0 ) then
!               if( ExtParticle(6)%Helicity.eq.0 ) cycle!   this can be more elegantly done in mod_process
!               call HDecay(ExtParticle(6),DK_LO,MomExt(1:4,13:14))
!           endif

          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp),PolExt3(1:4),The2ndResult(iPrimAmp))
!               call EvalTree(BornAmps(iPrimAmp))
          enddo
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,4
          do iPrimAmp=1,4
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbqqbg(iPrimAmp,jPrimAmp) * &
                           (  BornAmps(iPrimAmp)%Result *dconjg(BornAmps(jPrimAmp)%Result) &
                            + The2ndResult(iPrimAmp)    *dconjg(The2ndResult(jPrimAmp))    &
                           )
          enddo
          enddo
          LO_Res_Unpol = LO_Res_Unpol + LO_Res_Pol * ( PDFFac(Up_)+PDFFac(Dn_) )
        enddo!helicity loop
        
        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**3 
        EvalCS_Real_ttbqqbgH = EvalCS_Real_ttbqqbgH + dble(LO_Res_Unpol*PreFac)
        do NHisto=1,NumHistograms
               call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol*PreFac),BinValue=PObs(NHisto))
        enddo
        EvalCounter = EvalCounter + 1
endif!applyPSCut


    PreFacDip = PreFac * ISFac * (alpha_s4Pi*RunFactor)**3  /PSWgt2/PSWgt3/PSWgt4
    IF( PROCESS.EQ.106 ) THEN
        call EvalDipoles_QQBTTBGH((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),-MomExt(1:4,2),MomExt(1:4,3)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),npdf,DipoleResult)
    ELSEIF( PROCESS.EQ.103 ) THEN
        call EvalDipoles_QGTTBQH((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),-MomExt(1:4,1),MomExt(1:4,3),-MomExt(1:4,2)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),npdf,DipoleResult)
    ELSEIF( PROCESS.EQ.104 ) THEN
        call EvalDipoles_QBGTTBQBH((/MomExt(1:4,5),MomExt(1:4,4),MomExt(1:4,6),MomExt(1:4,3),-MomExt(1:4,1),-MomExt(1:4,2)/),yRnd(11:20),(/PreFacDip*PDFFac(1),PreFacDip*PDFFac(2)/),npdf,DipoleResult)
    ENDIF
    EvalCS_Dips_ttbqqbgH = EvalCS_Dips_ttbqqbgH + DipoleResult(1) + DipoleResult(2)

  ENDDO! loop over a<-->b pdfs


!      sij = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
! !      sij = MomExt(1,3)**2
!      print *,  sij/EHat**2,EvalCS_Real_ttbqqbgH,EvalCS_Dips_ttbqqbgH,(1d0+EvalCS_Real_ttbqqbgH/EvalCS_Dips_ttbqqbgH)
!      pause


    EvalCS_Real_ttbqqbgH = (EvalCS_Real_ttbqqbgH + EvalCS_Dips_ttbqqbgH) /VgsWgt



!       !careful, alphas=0.13 only for NLOParam=0 PDFSet=2
!       !./TOPAZ Collider=1 TopDK=0 HDK=0 Process=106 Correction=2 NLOParam=0 PDFSet=2 MTop=1.73 ObsSet=81
!       !MADGRAPH CHECK: qqb->ttbHg, mt=173, alpha_s=0.13 mH=125.00, vev=250.618249228543
!       MG_MOM(0:3,1) = MomExt(1:4,1)*100d0
!       MG_MOM(0:3,2) = MomExt(1:4,2)*100d0
!       MG_MOM(0:3,3) = MomExt(1:4,6)*100d0
!       MG_MOM(0:3,4) = MomExt(1:4,5)*100d0
!       MG_MOM(0:3,5) = MomExt(1:4,4)*100d0
!       MG_MOM(0:3,6) = MomExt(1:4,3)*100d0
!       call coupsm(0)
! !       call SUUB_TTBHG(MG_MOM,MadGraph_tree)
! !       call SUG_TTBHU(MG_MOM,MadGraph_tree)
!       call SUBG_TTBHUB(MG_MOM,MadGraph_tree)
!       print *, ""
!       print *, alpha_s*RunFactor,m_top,m_h
!       print *, "My tree:         ", LO_Res_Unpol/(100d0)**4
!       print *, "MadGraph hel.amp:", MadGraph_tree
!       print *, "MG/ME ratio: ", MadGraph_tree/(dble(LO_Res_Unpol)/(100d0)**4)
!       pause
    
    
END FUNCTION










FUNCTION EvalCS_NLODK_ttbH_MPI(yRnd,VgsWgt,res)
implicit none
integer :: EvalCS_NLODK_ttbH_MPI
real(8) ::  yRnd(*),res(*),VgsWgt

res(1) = EvalCS_NLODK_ttbH(yRnd,VgsWgt)
EvalCS_NLODK_ttbH_MPI=0
RETURN
END FUNCTION






FUNCTION EvalCS_NLODK_ttbH(yRnd,VgsWgt)
use ModProcess
use ModKinematics
use ModAmplitudes
use ModMyRecurrence
use ModParameters
use ModHadrWDecay
use ModHDecay
implicit none
real(8) ::  EvalCS_NLODK_ttbH,yRnd(1:VegasMxDim),VgsWgt
complex(8) :: LO_Res_Pol,LO_Res_Unpol,Dip_Res_Unpol,NLO_Res_Pol,NLO_Res_UnPol
complex(8) :: TreeResult(1:NumBornAmps),DKResult(1:NumBornAmps)
integer :: iHel,jHel,kHel,GluHel,iPrimAmp,jPrimAmp,ndip
real(8) :: EHat,PSWgt1,PSWgt2,PSWgt3,ISFac,dip_res_w,MH_INV
real(8) :: MomExt(1:4,1:14),MomExtTd(1:4,1:14)
logical :: applyPSCut,applySingCut
real(8) :: tau,eta1,eta2,sHatJacobi,PreFac,FluxFac,PDFFac_a(1:2),PDFFac_b(1:2),PDFFac(1:2),RunFactor
real(8) :: pdf(-6:6,1:2)
! real(8) :: couplZUU,couplZDD,couplZLL,couplGUU,couplGDD,couplGLL,couplZTT,couplGTT
integer :: NBin(1:NumMaxHisto),NHisto,npdf
real(8) :: pbDpg,ptDpg,ptDpb,z,omz,Dipole,rsq,y,PObs(1:NumMaxHisto)
real(8), parameter :: CF=4d0/3d0,PhotonCouplCorr=2d0
real(8) :: MomBoost(1:4),MomLep1(1:4),MomLep2(1:4)
integer,parameter :: up=1,dn=2,glu=1
real(8) :: PSWgt4
include "vegas_common.f"



! yrnd(1:20) = 0.3d0; print *, "fixed yrnd"


  EvalCS_NLODK_ttbH = 0d0
  call PDFMapping(2,yRnd(1:2),eta1,eta2,Ehat,sHatJacobi)
!   IF( CORRECTION.EQ.4 ) call SmearExternal(yRnd(18),M_Z,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)
!   IF( CORRECTION.EQ.5 ) call SmearExternal(yRnd(21),M_Z,Ga_ZExp,Zero,EHat,MH_Inv,MHJacobi)
  MH_Inv = m_H
  if( EHat.le.2d0*m_Top+MH_Inv ) then
      EvalCS_NLODK_ttbH = 0d0
      return
  endif
  FluxFac = 1d0/(2d0*EHat**2)
!   if ( HDecays .lt. 10) then
!       propH = (1d0,0d0)/dsqrt(2d0*Ga_Zexp*m_Z)
!   elseif (HDecays .gt. 10) then
!       propH=(1d0,0d0)/(MH_Inv**2-m_Z**2+ci*Ga_ZExp*m_Z)
!   endif
!   propH=MH_Inv**2*propH


  call EvalPhaseSpace_2to3M(EHat,MH_Inv,yRnd(3:7),MomExt(1:4,1:5),PSWgt1)
  call boost2Lab(eta1,eta2,5,MomExt(1:4,1:5))

  call setPDFs(eta1,eta2,MuFac,pdf)
  IF( PROCESS.EQ.101 ) THEN
      PDFFac(glu)   = pdf(0,1) * pdf(0,2)
      PDFFac_a(glu) = PDFFac(glu)
      PDFFac_b(glu) = 0d0
  ELSEIF( PROCESS.EQ.102 ) THEN
      PDFFac_a(up) = pdf(Up_,1)*pdf(AUp_,2) + pdf(Chm_,1)*pdf(AChm_,2)
      PDFFac_a(dn) = pdf(Dn_,1)*pdf(ADn_,2) + pdf(Str_,1)*pdf(AStr_,2) + pdf(Bot_,1)*pdf(ABot_,2)
      PDFFac_b(up) = pdf(Up_,2)*pdf(AUp_,1) + pdf(Chm_,2)*pdf(AChm_,1)
      PDFFac_b(dn) = pdf(Dn_,2)*pdf(ADn_,1) + pdf(Str_,2)*pdf(AStr_,1) + pdf(Bot_,2)*pdf(ABot_,1)
  ENDIF
  PSWgt4=1d0


IF( CORRECTION.EQ.4 ) THEN
!----------------------------------------
! one loop correction to Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:15),.false.,MomExt(1:4,9:11),PSWgt3)
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(16:17),MomExt(1:4,12:13),PSWgt4)

   call Kinematics_TTBARH(0,MomExt(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      EvalCS_NLODK_ttbH = 0d0
      return
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3*PSWgt4 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
   if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
   elseif(npdf.eq.2) then
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        if( Process.eq.101 ) cycle
   endif
   ISFac = MomCrossing(MomExt)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(1),DK_1L_T,MomExt(1:4,6:8))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      if(PROCESS.EQ.101) then
          NLO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
         NLO_Res_Pol = ColLO_ttbqqb(1,1) * dreal( TreeResult(1)*dconjg(DKResult(1)) ) * ( PDFFac(up)+PDFFac(dn) )
      endif
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol


      if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.4 ) then  !  virt.corr. to hadronic W decay
          call TopDecay(ExtParticle(1),DK_1L_Q,MomExt(1:4,6:8))
          do iPrimAmp=1,NumBornAmps
              call EvalTree(BornAmps(iPrimAmp))
              DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
          enddo
          if(PROCESS.EQ.101) then
              NLO_Res_Pol = (0d0,0d0)
              do jPrimAmp=1,NumBornAmps
              do iPrimAmp=1,NumBornAmps
                  NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
              enddo
              enddo
          elseif(PROCESS.EQ.102) then
              NLO_Res_Pol = ColLO_ttbqqb(1,1) * dreal( TreeResult(1)*dconjg(DKResult(1)) ) * ( PDFFac(up)+PDFFac(dn) )
          endif
          NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
      endif!  virt.corr. to hadronic W decay

   enddo!helicity loop


!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol),BinValue=PObs(NHisto))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back






!-------------------------------------
! one loop correction to top-decay   |
!-------------------------------------
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3*PSWgt4 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        PDFFac(1:2) = PDFFac_b(1:2)
        if( Process.eq.101 ) cycle
    endif
    ISFac = MomCrossing(MomExt)

   LO_Res_Unpol = (0d0,0d0)
   NLO_Res_UnPol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,9:11))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          TreeResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo

      call TopDecay(ExtParticle(2),DK_1L_T,MomExt(1:4,9:11))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
          DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
      enddo
      if(PROCESS.EQ.101) then
          NLO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
          NLO_Res_Pol = ColLO_ttbqqb(1,1) * dreal( TreeResult(1)*dconjg(DKResult(1)) ) * ( PDFFac(up)+PDFFac(dn) )
      endif
      NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol


      if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.3 ) then  !  virt.corr. to hadronic W decay
        call TopDecay(ExtParticle(2),DK_1L_Q,MomExt(1:4,9:11))
        do iPrimAmp=1,NumBornAmps
            call EvalTree(BornAmps(iPrimAmp))
            DKResult(iPrimAmp) = BornAmps(iPrimAmp)%Result
        enddo
        if(PROCESS.EQ.101) then
            NLO_Res_Pol = (0d0,0d0)
            do jPrimAmp=1,NumBornAmps
            do iPrimAmp=1,NumBornAmps
                NLO_Res_Pol = NLO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal( TreeResult(iPrimAmp)*dconjg(DKResult(jPrimAmp)) )*PDFFac(glu)
            enddo
            enddo
        elseif(PROCESS.EQ.102) then
            NLO_Res_Pol = ColLO_ttbqqb(1,1) * dreal( TreeResult(1)*dconjg(DKResult(1)) ) * ( PDFFac(up)+PDFFac(dn) )
        endif
        NLO_Res_UnPol = NLO_Res_UnPol + NLO_Res_Pol
      endif!  virt.corr. to hadronic W decay

   enddo!helicity loop

!  normalization
   NLO_Res_Unpol = NLO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(NLO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(NLO_Res_Unpol),BinValue=PObs(NHisto))
   enddo

enddo! npdf loop
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back





ELSEIF( CORRECTION.EQ.5 ) THEN
if( DKRE_switch.eq.0 .or. DKRE_switch.eq.1 ) then
!----------------------------------------
! real gluon emission for Anti-top decay |
!----------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
   call CheckSing(MomExt(1:4,6:9),applySingCut)
   if( applySingCut) then
      EvalCS_NLODK_ttbH = 0d0
      SkipCounter = SkipCounter + 1
      goto 13
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3*PSWgt4 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        PDFFac(1:2) = PDFFac_b(1:2)
        if( Process.eq.101 ) cycle
    endif
    ISFac = MomCrossing(MomExt)
    call Kinematics_TTBARH(1,MomExt(1:4,1:14),(/4,5,3,1,2,9,6,7,8,10,11,12,13,14/),applyPSCut,NBin,PObs)
    if( applyPSCut ) then
      goto 14
    endif
   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_RE_T,MomExt(1:4,6:9),GluonHel=GluHel)
      call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))
      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo

      if(PROCESS.EQ.101) then
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result) * PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
          LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
      endif
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop

!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(LO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol),BinValue=PObs(NHisto))
   enddo
   EvalCounter = EvalCounter + 1


14 continue

!-------------------------------------
! dipole subtraction for Atop-decay |
!-------------------------------------
   call WTransform(MomExt(1:4,6:9),MomExtTd(1:4,6:8),pbDpg,ptDpg,ptDpb)
   omz=ptDpg/(ptDpb+ptDpg-pbDpg)  !  for some reason this is not (1-z) as defined in the paper...
   rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
   z=1d0-omz
   y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
   Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
   Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   MomExtTd(1:4,1:5)  = MomExt(1:4,1:5)
   MomExtTd(1:4,9:13) = MomExt(1:4,10:14)
   call Kinematics_TTBARH(0,MomExtTd(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
   if( applyPSCut ) cycle ! = goto next npdf

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      if(PROCESS.EQ.101) then
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
          LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
      endif
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + Dip_Res_Unpol

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol),BinValue=PObs(NHisto))
   enddo

!   print *, npdf,(MomExt(1:4,6).dot.MomExt(1:4,9))/m_top**2
!   print *, MomExt(1,9)**2/m_top**2
!!   print *, dble(LO_Res_Unpol),dble(Dip_Res_Unpol),dble(LO_Res_Unpol)/dble(Dip_Res_Unpol) + 1d0
!   pause
enddo! npdf loop


13 continue
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back






if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.4 ) then
          call EvalPhasespace_TopDecay2(MomExt(1:4,4),yRnd(8:14),.true.,MomExt(1:4,6:9),PSWgt2)
          call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(15:18),.false.,MomExt(1:4,10:12),PSWgt3)
!           call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
          if( PSWgt2 .eq. 0d0 ) then  ! this rejects "too singular" events, similar to CheckSing
              EvalCS_NLODK_ttbH = 0d0
              SkipCounter = SkipCounter + 1
              goto 19
          endif

          PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3*PSWgt4 * VgsWgt
          RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
          if(npdf.eq.1) then
              PDFFac(1:2) = PDFFac_a(1:2)
          elseif(npdf.eq.2) then
              call swapMom(MomExt(1:4,1),MomExt(1:4,2))
              PDFFac(1:2) = PDFFac_b(1:2)
              if( Process.eq.101 ) cycle
          endif
          ISFac = MomCrossing(MomExt)
          call Kinematics_TTBARH(1,MomExt(1:4,1:14),(/4,5,3,1,2,9,6,7,8,10,11,12,13,14/),applyPSCut,NBin,PObs)
          if( applyPSCut ) then
            goto 16
          endif


        LO_Res_Unpol = (0d0,0d0)
        do iHel=1,NumHelicities ! loop over initial state chiralities
        do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
            call HelCrossing(Helicities(iHel,1:NumExtParticles))
            call SetPolarizations()
            call TopDecay(ExtParticle(1),DK_RE_Q,MomExt(1:4,6:9),GluonHel=GluHel)
            call TopDecay(ExtParticle(2),DK_LO,MomExt(1:4,10:12))

            do iPrimAmp=1,NumBornAmps
                call EvalTree(BornAmps(iPrimAmp))
            enddo
            if(PROCESS.EQ.101) then
                LO_Res_Pol = (0d0,0d0)
                do jPrimAmp=1,NumBornAmps
                do iPrimAmp=1,NumBornAmps
                    LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
                enddo
                enddo
            elseif(PROCESS.EQ.102) then
                LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
            endif
            LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
        enddo!helicity loop
        enddo!helicity loop

      !  normalization
        LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * PreFac
        EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(LO_Res_Unpol)
        do NHisto=1,NumHistograms
            call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol),BinValue=PObs(NHisto))
        enddo


16 continue


          do ndip=1,2!   there are two dipoles
              call wdec_trans(ndip,MomExt(1:4,6:9),MomExtTd(1:4,6:8),alpha_DKWff,dip_res_w)
              if( dip_res_w.eq.0d0 ) cycle
              Dipole = - alpha_s4Pi*RunFactor * CF * dip_res_w
              MomExtTd(1:4,1:5)  = MomExt(1:4,1:5)
              MomExtTd(1:4,9:13) = MomExt(1:4,10:14)
              call Kinematics_TTBARH(0,MomExtTd(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
              if( applyPSCut ) cycle! = goto next dipole

              Dip_Res_Unpol= (0d0,0d0)
              do iHel=1,NumHelicities
                  call HelCrossing(Helicities(iHel,1:NumExtParticles))
                  call SetPolarizations()
                  call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
                  call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

                  do iPrimAmp=1,NumBornAmps
                      call EvalTree(BornAmps(iPrimAmp))
                  enddo
                  if(PROCESS.EQ.101) then
                      LO_Res_Pol = (0d0,0d0)
                      do jPrimAmp=1,NumBornAmps
                      do iPrimAmp=1,NumBornAmps
                          LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
                      enddo
                      enddo
                  elseif(PROCESS.EQ.102) then
                      LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
                  endif
                  Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
              enddo!helicity loop

!             normalization
              Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * Dipole * PreFac
              EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + Dip_Res_Unpol
              do NHisto=1,NumHistograms
                  call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol),BinValue=PObs(NHisto))
              enddo
!              print *, 'dipoles in W decay'
!              print *, MomExt(1,9)**2/m_W**2
!!                         print *, npdf,(MomExt(1:4,7).dot.MomExt(1:4,9))/m_W**2,(MomExt(1:4,8).dot.MomExt(1:4,9))/m_W**2
!                         print *, ndip,dble(LO_Res_Unpol),dble(Dip_Res_Unpol),dble(LO_Res_Unpol)/dble(Dip_Res_Unpol) + 1d0
!                         pause
          enddo!dipole loop


enddo! npdf loop
19 continue
call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back
endif! correction on W



endif! DKRE_switch





if( DKRE_switch.eq.0 .or. DKRE_switch.eq.2 ) then
!-------------------------------------
! real gluon emission for top-decay  |
!-------------------------------------
   call EvalPhasespace_TopDecay(MomExt(1:4,4),yRnd(8:11),.false.,MomExt(1:4,6:8),PSWgt2)
   call EvalPhasespace_TopDecay(MomExt(1:4,5),yRnd(12:18),.true.,MomExt(1:4,9:12),PSWgt3)
!    call EvalPhasespace_HDecay(MomExt(1:4,3),yRnd(19:20),MomExt(1:4,13:14),PSWgt4)
   call CheckSing(MomExt(1:4,9:12),applySingCut)
   if( applySingCut ) then
      SkipCounter = SkipCounter + 1
      goto 17
   endif
   PreFac = fbGeV2 * FluxFac * sHatJacobi * PSWgt1*PSWgt2*PSWgt3*PSWgt4 * VgsWgt
   RunFactor = RunAlphaS(2,MuRen)

do npdf=1,2
    if(npdf.eq.1) then
        PDFFac(1:2) = PDFFac_a(1:2)
    elseif(npdf.eq.2) then
        PDFFac(1:2) = PDFFac_b(1:2)
        call swapMom(MomExt(1:4,1),MomExt(1:4,2))
        if( Process.eq.101 ) cycle
    endif
    ISFac = MomCrossing(MomExt)
    call Kinematics_TTBARH(1,MomExt(1:4,1:14),(/4,5,3,1,2,12,6,7,8,9,10,11,13,14/),applyPSCut,NBin,PObs)
    if( applyPSCut ) then
      goto 15
    endif

   LO_Res_Unpol = (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
   do GluHel=1,-1,-2 ! loop over additional gluon chiralities from decay
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExt(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_RE_T,MomExt(1:4,9:12),GluonHel=GluHel)

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      if(PROCESS.EQ.101) then
           LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
          LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
      endif
      LO_Res_UnPol = LO_Res_UnPol + LO_Res_Pol
   enddo!helicity loop
   enddo!helicity loop



!  normalization
   LO_Res_Unpol = LO_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(LO_Res_Unpol)

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(LO_Res_Unpol),BinValue=PObs(NHisto))
   enddo
   EvalCounter = EvalCounter + 1


15 continue


!-------------------------------------
! dipole subtraction for top-decay   |
!-------------------------------------
  call WTransform(MomExt(1:4,9:12),MomExtTd(1:4,9:11),pbDpg,ptDpg,ptDpb)
  omz=ptDpg/(ptDpb+ptDpg-pbDpg)
  rsq = 1d0 - 2d0/m_top**2*(ptDpb+ptDpg-pbDpg)
  z=1d0-omz
  y=pbDpg*2d0/m_top**2/(1d0-dsqrt(rsq))**2
  Dipole = - alpha_s4Pi*RunFactor * CF * ( 1d0/pbDpg*(2d0/omz-1d0-z) - (m_Top/ptDpg)**2 )
  Dipole = Dipole * (1d0 - StepFunc(1d0-alpha_DKTfi-z) * StepFunc(y-alpha_DKTfi*(1d0+dsqrt(rsq))**2*z*omz/(z+rsq*omz)) )

   MomExtTd(1:4,1:8) = MomExt(1:4,1:8)
   MomExtTd(1:4,12:13) = MomExt(1:4,13:14)
   call Kinematics_TTBARH(0,MomExtTd(1:4,1:14),(/4,5,3,1,2,0,6,7,8,9,10,11,12,13/),applyPSCut,NBin,PObs)
   if( applyPSCut ) then
      goto 17
   endif

   Dip_Res_Unpol= (0d0,0d0)
   do iHel=1,NumHelicities ! loop over initial state chiralities
      call HelCrossing(Helicities(iHel,1:NumExtParticles))
      call SetPolarizations()
      call TopDecay(ExtParticle(1),DK_LO,MomExtTd(1:4,6:8))
      call TopDecay(ExtParticle(2),DK_LO,MomExtTd(1:4,9:11))

      do iPrimAmp=1,NumBornAmps
          call EvalTree(BornAmps(iPrimAmp))
      enddo
      if(PROCESS.EQ.101) then
          LO_Res_Pol = (0d0,0d0)
          do jPrimAmp=1,NumBornAmps
          do iPrimAmp=1,NumBornAmps
              LO_Res_Pol = LO_Res_Pol + ColLO_ttbgg(iPrimAmp,jPrimAmp) * BornAmps(iPrimAmp)%Result*dconjg(BornAmps(jPrimAmp)%Result)*PDFFac(glu)
          enddo
          enddo
      elseif(PROCESS.EQ.102) then
          LO_Res_Pol = ColLO_ttbqqb(1,1) * BornAmps(1)%Result*dconjg(BornAmps(1)%Result)  * ( PDFFac(up)+PDFFac(dn) )
      endif
      Dip_Res_UnPol = Dip_Res_UnPol + LO_Res_Pol
   enddo!helicity loop

!  normalization
   Dip_Res_Unpol = Dip_Res_Unpol * ISFac * (alpha_s4Pi*RunFactor)**2 * Dipole  * PreFac
   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH + dble(Dip_Res_Unpol)

!             print *, npdf,(MomExt(1:4,9).dot.MomExt(1:4,12))/m_top**2
!             print *, MomExt(1,9)**2/m_top**2
!             print *, dble(LO_Res_Unpol),dble(Dip_Res_Unpol),dble(LO_Res_Unpol)/dble(Dip_Res_Unpol) + 1d0
!             pause

   do NHisto=1,NumHistograms
      call intoHisto(NHisto,NBin(NHisto),dble(Dip_Res_Unpol),BinValue=PObs(NHisto))
   enddo

enddo! npdf loop
! call swapMom(MomExt(1:4,1),MomExt(1:4,2))! swap back

if( TOPDECAYS.eq.2 .or. TOPDECAYS.eq.3 ) then
    call Error("Real correction on W+ for TopDecay=2,3 is not yet implemented.")
endif

endif! DKRE_switch
17 continue
ENDIF


   EvalCS_NLODK_ttbH = EvalCS_NLODK_ttbH/VgsWgt

RETURN
END FUNCTION








!!   subroutine SigmaRenorm_gg(plept,TheBornAmps,RenormAmps,Renorm_Res)
!! ! UV counterterm to renormalize sigma-q couplings in gg channel
!! ! Input  : array of BornAmps
!! !        : momenta of leptons, needed only for calling of ZDecay (which redefines the dyn coupl)
!! ! Output : array of counterterm amps (= LO amps with C1V=C1A=0)
!! !          interference between CTamps and BornAmps, incl color factor CF
!! !NB: this last NOT summed over helicity! So this routine should be called INSIDE a helicity sum.
!!     use ModHDecay
!!     use ModProcess
!!     use ModParameters
!!     use ModAmplitudes
!!     implicit none
!!     type(BornAmplitude),target :: TheBornAmps(1:NumBornAmps)
!!     complex(8)  :: RenormAmps(1:NumBornAmps)
!!     real(8)     :: Renorm_Res,plept(1:4,1:2)
!!     integer     :: iPrimAmp,jPrimAmp
!!        
!! 
!!     call StoreTopZCouplings
!! 
!! ! now set the SM-like couplings to zero, and recalculate the LO amps
!!     couplZTT_left  = 0d0
!!     couplZTT_right = 0d0
!!     Q_Top=0d0
!!     
!!     if( HDecays.gt.0 ) then
!!        call HDecay(ExtParticle(5),DK_LO,plept(1:4,1:2))
!!     else
!!        couplZTT_left_dyn=couplZTT_left
!!        couplZTT_right_dyn=couplZTT_right
!!     endif
!!     
!!     do iPrimAmp=1,NumBornAmps
!!        call EvalTree2(BornAmps(iPrimAmp)%TreeProc,RenormAmps(iPrimAmp))
!!     enddo
!!     
!!     Renorm_Res=0d0
!!     do jPrimAmp=1,2
!!        do iPrimAmp=1,2
!!           !  incl CF factor from vertex correction.
!!           Renorm_Res =Renorm_Res + 4d0/3d0*ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal(BornAmps(iPrimAmp)%Result*dconjg(RenormAmps(jPrimAmp)))
!!        enddo
!!     enddo
!!     
!!     call RetrieveTopZCouplings 
!! 
!!   end subroutine SigmaRenorm_gg
!!     
!! 
!!     
!!     
!!     
!!   subroutine SigmaRenorm_qqb(plept,TheBornAmps,RenormAmps,RenormPartAmps)
!! ! UV counterterm to renormalize sigma-q couplings in qqb channel
!! ! Input  : array of BornAmps
!! !         : array of BornAmps incl qqbZ factors etc
!! !        : momenta of leptons, needed only for calling of ZDecay (which redefines the dyn coupl)
!! ! Output : array of counterterm amps (= LO amps with C1V=C1A=0)
!! ! differences with gg channel: combination into Z on qqb and ttb line; color factors; interference with Born not computed
!!     use ModHDecay
!!     use ModProcess
!!     use ModParameters
!!     use ModAmplitudes
!!     implicit none
!!     type(BornAmplitude),target :: TheBornAmps(1:NumBornAmps)
!!     integer,parameter :: up=1,dn=2
!!     complex(8)  :: RenormAmps(1:NumBornAmps),RenormPartAmps(1:2)
!!     real(8)     :: plept(1:4,1:2)
!!     integer     :: iPrimAmp,jPrimAmp
!!        
!!     call StoreTopZCouplings
!! 
!! ! now set the SM-like couplings to zero, and recalculate the LO amps
!!     couplZTT_left  = 0d0
!!     couplZTT_right = 0d0
!!     Q_Top=0d0
!!     
!!     if( HDecays.gt.0 ) then
!!        call HDecay(ExtParticle(5),DK_LO,plept(1:4,1:2))
!!     else
!!        couplZTT_left_dyn=couplZTT_left
!!        couplZTT_right_dyn=couplZTT_right
!!     endif
!!     
!!     do iPrimAmp=1,NumBornAmps
!!        call EvalTree2(BornAmps(iPrimAmp)%TreeProc,RenormAmps(iPrimAmp))
!!     enddo
!! 
!! ! don't include amplitudes with Z attached to qqb line
!!     RenormPartAmps(up)=RenormAmps(1)
!!     RenormPartAmps(dn)=RenormAmps(1)
!!     
!!     call RetrieveTopZCouplings 
!!       
!!   end subroutine SigmaRenorm_qqb
!!   
!!   
!!   
  

  subroutine Gamma5Renorm_gg(plept,TheBornAmps,RenormAmps,Renorm_Res)
! Counterterm to renormalize gamma-5 in D-dim
! Input  : array of BornAmps
!        : momenta of leptons, needed only for calling of ZDecay (which redefines the dyn coupl)
! Output : array of counterterm amps (= LO amps with C1V=C1A=0)
!          interference between CTamps and BornAmps, incl color factor CF
!NB: this last NOT summed over helicity! So this routine should be called INSIDE a helicity sum.
    use ModHDecay
    use ModProcess
    use ModParameters
    use ModAmplitudes
    implicit none
    type(BornAmplitude),target :: TheBornAmps(1:NumBornAmps)
    complex(8)  :: RenormAmps(1:NumBornAmps)
    real(8)     :: Renorm_Res,plept(1:4,1:2)
    integer     :: iPrimAmp,jPrimAmp
       

    call StoreTopHCouplings

! now set all vector couplings to zero (i.e. L=-R=A) and recompute Born ampls

    if( HDecays.gt.0 ) then
       call HDecay(ExtParticle(5),DK_LO,plept(1:4,1:2))
    else
       couplZTT_left_dyn = m_top/vev * (0d0,1d0)*kappaTTBH_tilde
       couplZTT_right_dyn=-m_top/vev * (0d0,1d0)*kappaTTBH_tilde
    endif

    do iPrimAmp=1,2
       call EvalTree2(BornAmps(iPrimAmp)%TreeProc,RenormAmps(iPrimAmp))
    enddo
    Renorm_Res=0d0
    do jPrimAmp=1,2
       do iPrimAmp=1,2
          !  incl CF factor from vertex correction.
          Renorm_Res =Renorm_Res + 4d0/3d0*ColLO_ttbgg(iPrimAmp,jPrimAmp) * dreal(BornAmps(iPrimAmp)%Result*dconjg(RenormAmps(jPrimAmp)))
       enddo
    enddo

      
    call RetrieveTopHCouplings 

  end subroutine Gamma5Renorm_gg


  
  
  subroutine Gamma5Renorm_qqb(plept,TheBornAmps,RenormPartAmps)
! Counterterm to renormalize gamma-5 in D-dim
! Input  : array of BornAmps
!        : momenta of leptons, needed only for calling of ZDecay (which redefines the dyn coupl)
! Output : array of counterterm amps (= LO amps with C1V=C1A=0)
!          interference between CTamps and BornAmps, incl color factor CF
!NB: this last NOT summed over helicity! So this routine should be called INSIDE a helicity sum.
    use ModHDecay
    use ModProcess
    use ModParameters
    use ModAmplitudes
    implicit none
    type(BornAmplitude),target :: TheBornAmps(1:NumBornAmps)
    integer,parameter :: up=1,dn=2
    complex(8)  :: RenormAmps(1:NumBornAmps),RenormPartAmps(1:2)
    real(8)     :: plept(1:4,1:2)
    integer     :: iPrimAmp,jPrimAmp
    
    call StoreTopHCouplings
      

! now set all vector couplings to zero (i.e. L=-R=A) and recompute Born ampls
    if( HDecays.gt.0 ) then
       call HDecay(ExtParticle(5),DK_LO,plept(1:4,1:2))
    else
       couplZTT_left_dyn = m_top/vev * (0d0,1d0)*kappaTTBH_tilde
       couplZTT_right_dyn=-m_top/vev *(0d0,1d0)*kappaTTBH_tilde
    endif
    
    iPrimAmp=1
    call EvalTree2(BornAmps(iPrimAmp)%TreeProc,RenormAmps(iPrimAmp))
        
! don't include amplitudes with Z attached to qqb line
    RenormPartAmps(up)=RenormAmps(1)
    RenormPartAmps(dn)=RenormAmps(1)
      
    call RetrieveTopHCouplings 

  end subroutine Gamma5Renorm_qqb



  subroutine StoreTopHCouplings
    use ModParameters
    implicit none
    
    couplZTT_left_dyn_store  = couplZTT_left_dyn
    couplZTT_right_dyn_store = couplZTT_right_dyn

  end subroutine STORETOPHCOUPLINGS
  


  subroutine RetrieveTopHCouplings
    use ModParameters
    implicit none
    
    couplZTT_left_dyn   = couplZTT_left_dyn_store 
    couplZTT_right_dyn  = couplZTT_right_dyn_store 
    
  end subroutine RetrieveTopHCouplings
  










END MODULE


