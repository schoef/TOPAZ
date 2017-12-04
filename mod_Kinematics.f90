MODULE ModKinematics
implicit none
save


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
!     
    logical :: BinSmearing=.false.
    real(8) :: SmearSigma=0.1d0   
end type



type,extends(Histogram) :: HistogramMultiDim
    integer :: HistoDim

    integer :: NBins2
    real(8) :: BinSize2
    real(8) :: LowVal2
    real(8) :: SetScale2
    
    integer :: NBins3
    real(8) :: BinSize3
    real(8) :: LowVal3
    real(8) :: SetScale3
end type





integer,public :: it_sav

type(Histogram),allocatable   :: Histo(:)
type(HistogramMultiDim),allocatable       :: Histo3D(:),Histo2D(:)


real(8) :: pT_jet_cut, pT_bjet_cut, pT_lep_cut, Rsep_jet, Rsep_LepJet, pT_miss_cut, eta_sepa_cut, MInv_jets_cut, eta_lep_cut, eta_jet_cut, eta_bjet_cut, HT_cut, pT_hardestjet_cut
real(8) :: pT_pho_cut,Rsep_Pj,Rsep_Pbj,Rsep_Plep,eta_pho_cut,MTW_cut, Mttbar_cut,Rsep_jetlep,Rsep_leplep
real(8) :: pT_lepZ_cut,pt_lept_cut,pT_ll_cut,HT_jet_cut,Frac_sep_jetlep,MZ_window

real(8),public :: MInv_LB, MInv_T2







!DEC$ IF(_UseMPIVegas.EQ.1)
integer,public,parameter :: NUMHISTO=40       ! this has to match the constants in pvegas_mpi.c
integer,public,parameter :: MXHISTOBINS=500
type, BIND(C) :: ReducedHistogram
    real(8) :: Value(1:MXHISTOBINS)
    real(8) :: Value2(1:MXHISTOBINS)
    integer :: Hits(1:MXHISTOBINS)
end type
type(ReducedHistogram)  :: RedHisto(1:NUMHISTO)
public :: RedHisto,getRedHisto,transferHisto,clearRedHisto
!DEC$ ENDIF



contains





!DEC$ IF(_UseMPIVegas.EQ.1)
INTEGER FUNCTION getRedHisto(TheHisto,NHisto)
implicit none
type(ReducedHistogram) :: TheHisto
integer NHisto,NBin


  do NBin=1,MXHISTOBINS
    TheHisto%Value(NBin)  = RedHisto(NHisto)%Value(NBin)
    TheHisto%Value2(NBin) = RedHisto(NHisto)%Value2(NBin)
    TheHisto%Hits(NBin)   = RedHisto(NHisto)%Hits(NBin)
  enddo

getRedHisto=0
RETURN
END FUNCTION



INTEGER FUNCTION transferHisto(TheHisto,NHisto)
use ModParameters
implicit none
type(ReducedHistogram) :: TheHisto
integer NHisto,NBin
transferHisto=0

  if( NHisto.gt.NumHistograms ) return! this is required because in pvegas we loop until max.number of histograms (NUMHISTO)
  do NBin=1,Histo(NHisto)%NBins
    Histo(NHisto)%Value(NBin)  = TheHisto%Value(NBin)
    Histo(NHisto)%Value2(NBin) = TheHisto%Value2(NBin)
    Histo(NHisto)%Hits(NBin)   = TheHisto%Hits(NBin)
  enddo

RETURN
END FUNCTION



SUBROUTINE clearRedHisto()
implicit none
integer NHisto,NBin

  do NHisto=1,NUMHISTO
  do NBin=1,MXHISTOBINS
    RedHisto(NHisto)%Value(NBin)  = 0d0
    RedHisto(NHisto)%Value2(NBin) = 0d0
    RedHisto(NHisto)%Hits(NBin)   = 0
  enddo
  enddo

RETURN
END SUBROUTINE
!DEC$ ENDIF






SUBROUTINE InitPSCuts()
use ModMisc
use ModParameters
implicit none


pT_jet_cut        = 1d100
pT_bjet_cut       = 1d100
pT_lep_cut        = 1d100
Rsep_jet          = 1d100
Rsep_LepJet       = 1d100
Rsep_LepLep       = 1d100
pT_miss_cut       = 1d100
eta_sepa_cut      = 1d100
MInv_jets_cut     = 1d100
eta_lep_cut       = 1d100
eta_jet_cut       = 1d100
eta_bjet_cut      = 1d100
HT_cut            = 1d100
pT_hardestjet_cut = 1d100
pT_pho_cut        = 1d100
Rsep_Pj           = 1d100
Rsep_Pbj          = 1d100
Rsep_Plep         = 1d100
eta_pho_cut       = 1d100
MZ_window         = 1d100


IF( ObsSet.EQ.0 ) THEN! set of observables for ttb production without decays at Tevatron

ELSEIF( ObsSet.EQ.1 ) THEN! set of observables for ttb production without decays at LHC

ELSEIF( ObsSet.EQ.2 ) THEN! set of observables for ttb production as signal process at Tevatron (di-lept. decay)
    Rsep_jet    = 0.4d0         !*0d0
    pT_bjet_cut = 20d0*GeV      !*0d0
    eta_bjet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 20d0*GeV      !*0d0
    pT_miss_cut = 25d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*1d2

ELSEIF( ObsSet.EQ.3 ) THEN! set of observables for ttb production as signal process at LHC (di-lept. decay)
    Rsep_jet    = 0.4d0         !*0d0
    pT_bjet_cut = 25d0*GeV      !*0d0
    eta_bjet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 25d0*GeV      !*0d0
    pT_miss_cut = 50d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*1d2

ELSEIF( ObsSet.EQ.4 ) THEN! ! set of observables for ttb production with hadr. Atop, lept. top decay
    Rsep_jet    = 0.5d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.0d0
    pT_jet_cut  = 20d0*GeV
    eta_jet_cut = 2.0d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.0d0
    pT_miss_cut = 20d0*GeV
    HT_cut      = 220d0*GeV

ELSEIF( ObsSet.EQ.5 ) THEN! set of observables for ttb production with hadr. top, lept. Atop decay at LHC

!   these are the cuts for muons
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0

    pT_bjet_cut = 20d0*GeV
    pT_jet_cut  = 20d0*GeV
    eta_bjet_cut= 2.0d0
    eta_jet_cut = 2.5d0

    pT_miss_cut = 20d0*GeV
!   Mwt cut is hard coded below

    Rsep_LepJet = 0.0d0
    Rsep_jet    = 0.4d0

!     HT_cut      = 220d0*GeV
    HT_cut      = 380d0*GeV  ! this is for ratios to ttb+gamma and ttb+Z
!     HT_cut      = 420d0*GeV  ! this is mtt_cut for ratios to ttb+gamma


ELSEIF( ObsSet.EQ.6 ) THEN! set of observables for ttb production with lept. top, hadr. Atop decay at LHC
    Rsep_jet    = 0.4d0        ! *0d0
    pT_bjet_cut = 20d0*GeV     ! *0d0
    eta_bjet_cut= 2.5d0        ! *1d2
    pT_jet_cut  = 20d0*GeV     ! *0d0
    eta_jet_cut = 2.5d0        ! *1d2
    pT_lep_cut  = 20d0*GeV     ! *0d0
    eta_lep_cut = 2.5d0        ! *1d2
    pT_miss_cut = 20d0*GeV     ! *0d0


ELSEIF( ObsSet.EQ.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC
    Rsep_jet    = 0.5d0
    HT_cut      = 100d0*GeV
    pT_jet_cut  = 20d0*GeV
    pT_lep_cut  = 20d0*GeV

ELSEIF( ObsSet.EQ.8 ) THEN! set of observables for ttb spin correlations at LHC (di-lept. decay)
    pT_bjet_cut = 25d0*GeV
    pT_lep_cut  = 20d0*GeV
    pT_miss_cut = 40d0*GeV
    eta_lep_cut = 2.5d0
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.9 ) THEN! this is for the factorization check


ELSEIF( ObsSet.EQ.10 ) THEN! set of observables for ttbjet production without decays at Tevatron
    pT_jet_cut  = 20d0*GeV
    Rsep_jet    = 1d0

ELSEIF( ObsSet.EQ.11 ) THEN! set of observables for ttbjet production without decays at LHC
    pT_jet_cut  = 50d0*GeV
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.12 ) THEN! set of observables for ttbjet production as signal process at Tevatron (hadr.Atop, lept.top decay)
    pT_jet_cut  = 20d0*GeV
    pT_bjet_cut = pT_jet_cut
    eta_jet_cut = 2d0            *100d0
    eta_bjet_cut= eta_jet_cut    *100d0
    pT_lep_cut  = 20d0*GeV       *0d0
    eta_lep_cut = 1d0            *100d0
    pT_miss_cut = 20d0*GeV       *0d0
    HT_cut      = 220d0*GeV      *0d0
    Rsep_jet    = 0.5d0

ELSEIF( ObsSet.EQ.13 ) THEN! set of observables for ttbjet production as signal process at LHC (di-lept. decay)
    pT_jet_cut  = 25d0*GeV
    pT_bjet_cut = 25d0*GeV
    eta_jet_cut = 2.5d0
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 50d0*GeV
    Rsep_jet    = 0.4d0

ELSEIF( ObsSet.EQ.14 ) THEN! set of observables for ttbjet production as background process to VBF at LHC (di-lept. decay)
    pT_hardestjet_cut = 40d0*GeV
    pT_jet_cut   = 20d0*GeV
    pT_bjet_cut = pT_jet_cut
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2d0
    eta_sepa_cut = 3.0d0  ! 2.5d0 !  3.8d0
    MInv_jets_cut= 550d0*GeV
!     eta_jet_cut  = 2.5d0  ! REMOVED!!  this value is used in the opposite way, i.e. if(eta_jet .lt. eta_jet_cut) then cut
    eta_bjet_cut = eta_jet_cut
    Rsep_jet     = 0.5d0

ELSEIF( ObsSet.EQ.15 ) THEN! set of observables for ttbjet production as signal process at LHC (hadr.Atop, lept.top decay)
    pT_jet_cut  = 25d0*GeV
    pT_bjet_cut = pT_jet_cut
    eta_jet_cut = 2.5d0
    eta_bjet_cut= eta_jet_cut
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 30d0*GeV
    Rsep_jet    = 0.4d0
!   added an additional cut on MT(W)(+ET_miss) for ATLAS analysis, see kinematics_ttbjet subroutine


ELSEIF( ObsSet.EQ.19 ) THEN! for checks of ttbjet
    pT_jet_cut  = 20d0*GeV
    Rsep_jet    = 0.4d0


ELSEIF( ObsSet.EQ.20 ) THEN! set of observables for ttbgamma production without decays at Tevatron
!     Rsep_jet    = 1d0
    pT_pho_cut  = 20d0*GeV
    Rsep_Pj      = 0.4d0
!     Rsep_Plep    = 0.4d0


ELSEIF( ObsSet.EQ.21 ) THEN! set of observables for ttbgamma production without decays at LHC
!     Rsep_jet    = 1d0
    pT_pho_cut  = 20d0*GeV
    Rsep_Pj      = 0.4d0
!     Rsep_Plep    = 0.4d0


ELSEIF( ObsSet.EQ.22 ) THEN! set of observables for ttbgamma production with di-lept.decays at Tevatron
    Rsep_jet    = 0.4d0
    pT_pho_cut  = 10d0*GeV
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0
    Rsep_Plep   = 0.4d0

    pT_bjet_cut = 25d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 40d0*GeV

ELSEIF( ObsSet.EQ.23 ) THEN! set of observables for ttbgamma production with di-lept.decays at LHC
    pT_pho_cut  = 10d0*GeV
    eta_pho_cut = 5.0d0
    Rsep_Pj     = 0.2d0
    Rsep_Pbj    = 0.2d0
    Rsep_Plep   = 0.5d0

    pT_bjet_cut = 10d0*GeV  ! added this to define a jet 
    eta_bjet_cut= 5.0d0     ! added this to define a jet
    Rsep_jet    = 0.4d0

    pT_lep_cut  = 15d0*GeV
    eta_lep_cut = 5d0
    pT_miss_cut = 0d0*GeV

ELSEIF( ObsSet.EQ.24 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at Tevatron

    pT_pho_cut  = 10d0*GeV
    eta_pho_cut = 1.1d0
    Rsep_Plep   = 0.4d0
    Rsep_Pj     = 0.4d0
    Rsep_Pbj    = 0.4d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 15d0*GeV
    pT_jet_cut  = 15d0*GeV
    eta_bjet_cut= 2d0
    eta_jet_cut = 2d0

    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 1.1d0

    pT_miss_cut = 20d0*GeV
    HT_cut      = 200d0*GeV



! ObsSet=25 is found now together with ObsSet=26,27,28
! ObsSet=25 is the same as ObsSet=26 but with TopDK=4 instead of TopDK=3. 
! For the LHC asymmetry A_C=|y_t|-|y_tbar| one has to run both, ObsSet=25 and 26, in order to average 
! over the asymmetry induced by asymmetric cuts on hadr. and lept. decaying tops.
ELSEIF( ObsSet.EQ.25 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC
! 
! 
! ATLAS cuts 8TeV
!   these are the cuts for muons

    pT_lep_cut  = 15d0*GeV
    eta_lep_cut = 5d0

    pT_bjet_cut = 10d0*GeV ! added this to define a jet 
    pT_jet_cut  = 10d0*GeV ! added this to define a jet 
    eta_bjet_cut= 5.0d0    ! added this to define a jet 
    eta_jet_cut = 5.0d0    ! added this to define a jet 

    pT_miss_cut = 00d0*GeV
!   Mwt cut is hard coded below: removed

    pT_pho_cut  = 10d0*GeV
    eta_pho_cut = 5.0d0
!   cracks for photon are hard coded below: removed

    Rsep_LepJet = 0.0d0
    Rsep_jet    = 0.4d0
    Rsep_Pj     = 0.2d0
    Rsep_Pbj    = 0.2d0
    Rsep_Plep   = 0.5d0
    HT_cut      = 0d0*GeV

! 
! ! ATLAS cuts 7TeV
!  
! !   these are the cuts for muons
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
! 
!     pT_bjet_cut = 25d0*GeV
!     pT_jet_cut  = 25d0*GeV
!     eta_bjet_cut= 2.5d0
!     eta_jet_cut = 2.5d0
! 
!     pT_pho_cut  = 15d0*GeV
!     eta_pho_cut = 2.37d0
! !   cracks for photon are hard coded below
! 
!     pT_miss_cut = 25d0*GeV
! !   Mwt cut is hard coded below
! 
!     pT_pho_cut  = 15d0*GeV
!     eta_pho_cut = 2.37d0
! !   cracks for photon are hard coded below
! 
!     Rsep_LepJet = 0.4d0
!     Rsep_jet    = 0.4d0
!     Rsep_Pj     = 0.5d0
!     Rsep_Pbj    = 0.5d0
!     Rsep_Plep   = 0.4d0!  not specified




!   below is a copy of the ObsSet=28 case
!     pT_pho_cut  = 20d0*GeV
!     eta_pho_cut = 2.5d0
!     Rsep_Plep   = 0.4d0
!     Rsep_Pj     = 0.4d0
!     Rsep_Pbj    = 0.4d0
! 
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 20d0*GeV
!     pT_jet_cut  = 20d0*GeV
!     eta_bjet_cut= 2d0
!     eta_jet_cut = 2.5d0
! 
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
! 
!     pT_miss_cut = 20d0*GeV
!     HT_cut      = 200d0*GeV





! ELSEIF( ObsSet.EQ.26 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at Tevatron with Qtop cuts
! 
!     pT_pho_cut  = 20d0*GeV
!     eta_pho_cut = 2.5d0
!     Rsep_Plep   = 0.4d0
!     Rsep_Pj     = 0.4d0
!     Rsep_Pbj    = 0.4d0
! 
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 20d0*GeV
!     pT_jet_cut  = 20d0*GeV
!     eta_bjet_cut= 2d0
!     eta_jet_cut = 2.5d0
! 
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
! 
!     pT_miss_cut = 20d0*GeV
!     HT_cut      = 200d0*GeV


ELSEIF( ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28) THEN
! ELSEIF( ObsSet.EQ.25 .OR. ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC
                                                                                ! 26=no suppression cuts, 27=suppress PR, 28=suppress DK
    pT_pho_cut  = 20d0*GeV
    eta_pho_cut = 2.5d0         !*100d0
    Rsep_Plep   = 0.4d0         !*0d0
    Rsep_Pj     = 0.4d0         !*0d0
    Rsep_Pbj    = 0.4d0         !*0d0

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV      !*0d0
    pT_jet_cut  = 20d0*GeV      !*0d0
    eta_bjet_cut= 2.0d0         !*100d0
    eta_jet_cut = 2.5d0         !*100d0

    pT_lep_cut  = 20d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*100d0

    pT_miss_cut = 20d0*GeV      !*0d0
    HT_cut      = 200d0*GeV     !*0d0

! ELSEIF( ObsSet.EQ.28 ) THEN! set of observables for ttbgamma production with semi-lept.decays(hadr.Atop, lept.top decay) at LHC
! 
!     pT_pho_cut  = 20d0*GeV
!     eta_pho_cut = 2.5d0
!     Rsep_Plep   = 0.4d0
!     Rsep_Pj     = 0.4d0
!     Rsep_Pbj    = 0.4d0
! 
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 20d0*GeV
!     pT_jet_cut  = 20d0*GeV
!     eta_bjet_cut= 2d0
!     eta_jet_cut = 2.5d0
! 
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
! 
!     pT_miss_cut = 20d0*GeV
!     HT_cut      = 200d0*GeV

ELSEIF( ObsSet.EQ.29 ) THEN! this is for the factorization check

    pT_pho_cut  = 20d0*GeV
    Rsep_Pj     = 0.4d0




ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for HTHTbar + A0/BH production (stable)


ELSEIF( ObsSet.EQ.32 ) THEN! set of observables for HTHTbar + A0/BH production (di-lept. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV 
    eta_bjet_cut= 2.5d0 

    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0 
    pT_miss_cut = 25d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.33 ) THEN! set of observables for HTHTbar + A0/BH production (semi-hadr. tops)
    Rsep_jet    = 0.4d0

    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 30d0*GeV 
    eta_jet_cut= 2.5d0     
 
    pT_lep_cut  = 20d0*GeV  
    eta_lep_cut = 2.5d0      
    pT_miss_cut = 150d0*GeV! note that this is ET and not pT

    MTW_cut = 120d0*GeV



ELSEIF( ObsSet.EQ.34 ) THEN! set of observables for HTHTbar + A0/BH production (di-lept. tops) without acceptance cuts
    Rsep_jet    = 0d0

    pT_bjet_cut = 0d0*GeV 
    eta_bjet_cut= 100d0 

    pT_lep_cut  = 0d0*GeV  
    eta_lep_cut = 100d0 
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.35 ) THEN! set of observables for HTHTbar + A0/BH production (semi-hadr. tops) without acceptance cuts
    Rsep_jet    = 0d0

    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_jet_cut =  0d0*GeV 
    eta_jet_cut=  100d0     
 
    pT_lep_cut  = 0d0*GeV  
    eta_lep_cut = 100d0      
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT

    MTW_cut = 150d0*GeV  * 0d0



ELSEIF( ObsSet.EQ.36 ) THEN! set of observables for HTHTbar + A0/BH production  (semi-hadr. tops) for mtop determination


    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 20d0*GeV
    eta_jet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 20d0*GeV! note that this is ET and not pT
    MTW_cut = 0d0*GeV




ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar + Chi production (stable)



ELSEIF( ObsSet.EQ.42 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)
    Rsep_jet    = 0.4d0 
    pT_bjet_cut = 25d0*GeV 
    eta_bjet_cut= 2.5d0     
    pT_lep_cut  = 20d0*GeV       
    eta_lep_cut = 2.5d0            
    pT_miss_cut = 80d0*GeV      ! note that this is ET and not pT


ELSEIF( ObsSet.EQ.43 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops)

! !   these are the cuts for mstop/chi = 500/100 GeV analysis at 8TeV
!  if( Collider_Energy.eq.8000d0*GeV .or. Collider_Energy.eq.14000d0*GeV ) then

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 30d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 30d0*GeV
    eta_jet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 150d0*GeV! note that this is ET and not pT
    MTW_cut = 120d0*GeV

!  elseif( Collider_Energy.eq.7000d0*GeV ) then
! !   these are the cuts for mstop/chi = 300/100 GeV analysis at 7TeV
!     Rsep_jet    = 0.4d0
!     pT_bjet_cut = 30d0*GeV
!     eta_bjet_cut= 2.5d0
!     pT_jet_cut = 25d0*GeV
!     eta_jet_cut= 2.5d0
!     pT_lep_cut  = 20d0*GeV
!     eta_lep_cut = 2.5d0
!     pT_miss_cut = 125d0*GeV! note that this is ET and not pT
!     MTW_cut = 0d0*GeV
!  else
!     call Error("This ObsSet only supports Collider=1,11,12")
!  endif


ELSEIF( ObsSet.EQ.44 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops) without acceptance cuts

    Rsep_jet    = 0d0
    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_lep_cut  = 0d0*GeV
    eta_lep_cut = 100d0
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT


ELSEIF( ObsSet.EQ.45 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops) without acceptance cuts

    Rsep_jet    = 0d0
    pT_bjet_cut = 0d0*GeV
    eta_bjet_cut= 100d0
    pT_jet_cut = 0d0*GeV
    eta_jet_cut= 100d0
    pT_lep_cut  = 0d0*GeV
    eta_lep_cut = 100d0
    pT_miss_cut = 0d0*GeV! note that this is ET and not pT

    MTW_cut = 150d0*GeV  * 0d0



ELSEIF( ObsSet.EQ.46 ) THEN! set of observables for STSTbar + Chi production (semi-hadr. tops) for mtop determination


    Rsep_jet    = 0.4d0
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut = 20d0*GeV
    eta_jet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    eta_lep_cut = 2.5d0
    pT_miss_cut = 20d0*GeV! note that this is ET and not pT
    MTW_cut = 0d0*GeV


ELSEIF( ObsSet.EQ.48 ) THEN! set of observables for STOP width



ELSEIF( ObsSet.EQ.51 ) THEN! set of observables for ttb+Z (stable tops)

ELSEIF( ObsSet.EQ.52 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )


    Rsep_jet    = 0.4d0 
    pT_bjet_cut = 20d0*GeV
    eta_bjet_cut= 2.5d0
    pT_lep_cut  = 20d0*GeV
    pT_miss_cut = 40d0*GeV
    eta_lep_cut = 2.5d0



ELSEIF( ObsSet.EQ.53 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay )

    Rsep_jet    = 0.4d0
    pT_bjet_cut = 40d0*GeV
    eta_bjet_cut= 2.4d0
    pT_jet_cut  = 40d0*GeV
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 20d0*GeV
    pT_miss_cut = 20d0*GeV
    eta_lep_cut = 2.4d0
    Rsep_jetlep = 0.4d0
    MZ_window   = 15d0*GeV

ELSEIF( ObsSet.EQ.54 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay from 7 TeV CMS data)

   Rsep_jet    = 0.5d0
   pT_bjet_cut = 20d0*GeV
   eta_bjet_cut= 2.4d0
   pT_jet_cut  = 20d0*GeV
   eta_jet_cut = 2.4d0

   eta_lep_cut = 2.5d0
   pT_lepZ_cut  = 20d0*GeV
   pT_lept_cut  = 10d0*GeV
   pT_ll_cut  = 35d0*GeV
   HT_jet_cut = 120*GeV  

   Rsep_jetlep = 0.3d0
   Frac_sep_jetlep=0.15d0


ELSEIF( ObsSet.EQ.55 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay ) same as 52 but no cuts


    Rsep_jet    = 0.4d0         *0d0
    pT_bjet_cut = 25d0*GeV      *0d0
    eta_bjet_cut= 2.5d0         *1d2
    pT_lep_cut  = 25d0*GeV      *0d0
    pT_miss_cut = 50d0*GeV      *0d0
    eta_lep_cut = 2.5d0         *1d2


ELSEIF( ObsSet.EQ.56 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay ) same as 53 but no cuts

    Rsep_jet    = 0.4d0         *0d0
    pT_bjet_cut = 20d0*GeV      *0d0
    eta_bjet_cut= 2.5d0         *1d2
    pT_jet_cut  = 20d0*GeV      *0d0
    eta_jet_cut = 2.5d0         *1d2

    pT_lep_cut  = 15d0*GeV      *0d0 
    pT_miss_cut = 20d0*GeV      *0d0 
    eta_lep_cut = 2.5d0         *1d2
    Rsep_jetlep = 0.4d0         *0d0
    MZ_window   = 50d0*GeV      ! finite number needed to remove photon singularity at MZ^2 = 0       



ELSEIF( ObsSet.EQ.57 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay ) at Tevatron


    Rsep_jet    = 0.4d0
    pT_bjet_cut = 15d0*GeV
    eta_bjet_cut= 2.5d0
    pT_jet_cut  = 15d0*GeV
    eta_jet_cut = 2.5d0

    pT_lep_cut  = 15d0*GeV
    pT_miss_cut = 20d0*GeV
    eta_lep_cut = 2.5d0
    Rsep_jetlep = 0.4d0


ELSEIF (ObsSet.EQ.58) then    
! cuts correspoding to CMS 20 invfb analysis, cf. hep-ex/1406.7830

   Rsep_jet     = 0.5d0
   pT_jet_cut   = 30d0*GeV
   pT_bjet_cut  = pT_jet_cut
   eta_jet_cut  = 2.4d0
   eta_bjet_cut = eta_jet_cut
   
   pT_lep_cut   = 20d0*GeV
   eta_lep_cut  = 2.4d0
   pT_miss_cut  = 0d0*GeV
   MZ_window    = 20d0*GeV
   
   Rsep_jetlep  = 0.3d0
   Rsep_LepLep  = 0.3d0
   


ELSEIF ( ObsSet.EQ.60 ) THEN ! Zprime, stable top

   Mttbar_cut = 500d0*GeV

ELSEIF ( ObsSet.EQ.61 ) THEN ! Zprime, top decay to dileptons

   Rsep_jet = 0.5d0             !*0d0 !this removes all cuts for fact.checks

   pT_lep_cut  = 20d0*GeV       !*0d0
   eta_lep_cut = 2.5d0          !*1d6

   pT_bjet_cut = 30d0*GeV       !*0d0
   eta_bjet_cut = 2.5d0         !*1d6
   
ELSEIF ( ObsSet.EQ.62 ) THEN ! Zprime, fully hadronic top decay



ELSEIF ( ObsSet.EQ.64 ) THEN ! Zprime, semi-hadronic top decay (for factorization checks)
   
ELSEIF ( ObsSet.EQ.65 ) THEN ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)
   
!   this is for electrons and muons
    pT_lep_cut  = 25d0*GeV
    eta_lep_cut = 2.47d0

    Rsep_LepJet = 0.4d0


    pT_miss_cut = 30d0*GeV! note that this is ET and not pT
    MTW_cut = 30d0*GeV

    pT_jet_cut = 300d0*GeV
    eta_jet_cut= 2.0d0
!   more jet cuts are defined inside KinematicsZprimeTTB subroutine   




ELSEIF ( ObsSet.EQ.66 ) THEN ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)

!   this is for electrons
    pT_lep_cut  = 35d0*GeV
    eta_lep_cut = 2.5d0

!   this is for muons
!    pT_lep_cut  = 45d0*GeV
!    eta_lep_cut = 2.1d0

    Rsep_LepJet = 0.5d0
!   pTrel is defined inside KinematicsZprimeTTB subroutine

    pT_miss_cut = 50d0*GeV! note that this is ET and not pT
    HT_cut = 150d0*GeV

    Rsep_jet = 0.5d0
    pT_jet_cut = 150d0*GeV
    eta_jet_cut= 2.4d0
!   more jet cuts are defined inside KinematicsZprimeTTB subroutine   

   

ELSEIF ( ObsSet.EQ.67 ) THEN ! SM Z boson, stable top



ELSEIF ( ObsSet.EQ.70 ) THEN ! ee-->ttb with  stable tops


ELSEIF ( ObsSet.EQ.71 ) THEN ! ee-->ttb with di-leptonic tops
    Rsep_jet    = 0.4d0         !*0
    pT_bjet_cut = 20d0*GeV      !*0
    eta_bjet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 20d0*GeV      !*0
    pT_miss_cut = 40d0*GeV      !*0
    eta_lep_cut = 2.5d0         !*1d2


    
    

ELSEIF ( ObsSet.EQ.72 ) THEN ! ee-->ttb with semi-leptonic tops
    Rsep_jet    = 0.4d0         !*0
    pT_bjet_cut = 20d0*GeV      !*0
    pT_jet_cut = 20d0*GeV      !*0
    eta_bjet_cut= 2.5d0         !*1d2
    eta_jet_cut= 2.5d0         !*1d2
    pT_lep_cut  = 20d0*GeV      !*0
    pT_miss_cut = 20d0*GeV      !*0
    eta_lep_cut = 2.5d0         !*1d2




ELSEIF( ObsSet.EQ.81 ) THEN! set of observables for ttb+H (stable tops)


ELSEIF( ObsSet.EQ.82 ) THEN! set of observables for ttb+H (di-leptonic tops)


    Rsep_jet    = 0.4d0


ELSEIF( ObsSet.EQ.83 ) THEN! set of observables for ttb+H (semi-leptonic tops)


    Rsep_jet    = 0.4d0         !*0d0     
    pT_bjet_cut = 20d0*GeV      !*0d0
    pT_jet_cut  = 20d0*GeV      !*0d0
    eta_bjet_cut= 2.5d0         !*1d2
    eta_jet_cut = 2.5d0         !*1d2
    pT_lep_cut  = 20d0*GeV      !*0d0
    pT_miss_cut = 20d0*GeV      !*0d0
    eta_lep_cut = 2.5d0         !*1d2


ELSEIF( ObsSet.EQ.91 ) THEN! set of observables for t+H (stable tops)

ELSEIF( ObsSet.EQ.92 ) THEN! set of observables for tb+H (stable tops)

ELSEIF( ObsSet.EQ.93 ) THEN! set of observables for t+H leptonic top decay)
    Rsep_jet    = 0.5d0
ELSEIF( ObsSet.EQ.94 ) THEN! set of observables for tb+H (leptonic top decay)
    Rsep_jet    = 0.5d0
ENDIF



END SUBROUTINE





SUBROUTINE InitHisto()
use ModMisc
use ModParameters
implicit none
integer :: AllocStatus,NHisto,i,j

it_sav = 1
! RR 2015-04-15 -- modify in each ObsSet if needed
Num2DHistograms=0


IF( ObsSet.EQ.0 ) THEN! set of observables for ttb production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

ELSEIF( ObsSet.EQ.1 ) THEN! set of observables for ttb production without decays at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

ELSEIF( ObsSet.EQ.2 ) THEN! set of observables for ttb production as signal process at Tevatron (di-lept. decay)

          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be 1!")
          NumHistograms = 20
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 50d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_lep_bquark"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "Minv_lep_bquark"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

!          Histo(7)%Info   = "Psi"
!          Histo(7)%NBins  = 60
!          Histo(7)%BinSize= 0.1d0
!          Histo(7)%LowVal = -3d0
!          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "DeltaPhi"
          Histo(8)%NBins  = 65
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -3d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "E_bj1+E_bj2"
          Histo(9)%NBins  = 200
          Histo(9)%BinSize= 5d0*GeV
          Histo(9)%LowVal = 50d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "ET_miss"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 5d0*GeV
          Histo(10)%LowVal = 20d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "E_lep1+E_lep2"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 5d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pt_bj1"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 5d0*GeV
          Histo(12)%LowVal = 20d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pt_bj2"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 5d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "E_bj1"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "E_bj2"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 5d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "M_lj"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 8d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "ET_leptons"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 5d0*GeV
          Histo(17)%LowVal = 20d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "m_T"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 20d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "etaFB_lepMinus"
          Histo(19)%NBins  = 2
          Histo(19)%BinSize= 5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "etaFB_lepPlus"
          Histo(20)%NBins  = 2
          Histo(20)%BinSize= 5d0
          Histo(20)%LowVal =-5.0d0
          Histo(20)%SetScale= 1d0

ELSEIF( ObsSet.EQ.3 ) THEN! set of observables for ttb production as signal process at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be |1|!")
          NumHistograms = 20
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 25d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_lep_bquark"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "Psi"
          Histo(7)%NBins  = 60
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -3d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "DeltaPhi"
          Histo(8)%NBins  = 65
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -3d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "E_bj1+E_bj2"
          Histo(9)%NBins  = 200
          Histo(9)%BinSize= 5d0*GeV
          Histo(9)%LowVal = 50d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "ET_miss"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 5d0*GeV
          Histo(10)%LowVal = 20d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "E_lep1+E_lep2"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 5d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pt_bj1"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 0d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pt_bj2"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 5d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "E_bj1"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "E_bj2"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 5d0*GeV
          Histo(15)%LowVal = 20d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "Minv_leptons"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 10d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "ET_leptons"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 5d0*GeV
          Histo(17)%LowVal = 20d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "m_T"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 20d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "etaFB_lepMinus"
          Histo(19)%NBins  = 2
          Histo(19)%BinSize= 5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "etaFB_lepPlus"
          Histo(20)%NBins  = 2
          Histo(20)%BinSize= 5d0
          Histo(20)%LowVal =-5.0d0
          Histo(20)%SetScale= 1d0

ELSEIF( ObsSet.EQ.4 ) THEN! set of observables for ttb production with semi hadronic decay
          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Lep"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept+miss)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 25d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(lep+bjet)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pT_ttbar"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 15d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0



ELSEIF( ObsSet.EQ.5 ) THEN! set of observables for ttb production with hadr. top, lept. Atop decay
          if(Collider.ne.1 .and. Collider.ne.13 )  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "(pT_Top+pT_ATop)/2"
          Histo(1)%NBins  = 4
          Histo(1)%BinSize= 60d0*GeV
          Histo(1)%LowVal = 40d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "(y_ATop+y_Top)/2"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "(pT_Top+pT_ATop)/2 for MTT>660"
          Histo(3)%NBins  = 8
          Histo(3)%BinSize= 60d0*GeV
          Histo(3)%LowVal = 40d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "(y_ATop+y_Top)/2 for MTT>660"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "y_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 50d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(ttbar)"
          Histo(11)%NBins  = 90
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 340d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "<Ehat>"
          Histo(12)%NBins  = 1
          Histo(12)%BinSize= 10000d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale=0.01d0





ELSEIF( ObsSet.EQ.6 ) THEN! set of observables for ttb production with lept. top, hadr. Atop decay
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_LepP"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_LepP"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "HT(jets+lept)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 50d0*GeV
          Histo(10)%LowVal = 150d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "m(lep+bjet)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 10d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "pT_ttbar"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 15d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0



ELSEIF( ObsSet.EQ.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.6) call Error("TopDecays needs to be 6!")
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          Histo(1)%Info   = "M_LJPsi"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 5d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_JPsi"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_LepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "<M_LJPsi>"
          Histo(4)%NBins  = 1
          Histo(4)%BinSize= 10000d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale=1d0

          Histo(5)%Info   = "<M_LJPsi^2>"
          Histo(5)%NBins  = 1
          Histo(5)%BinSize= 10000d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale=1d0

          Histo(6)%Info   = "<x>"
          Histo(6)%NBins  = 1
          Histo(6)%BinSize= 10000d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale=1d0


ELSEIF( ObsSet.EQ.8 ) THEN! set of observables for ttb spin correlations at LHC (di-lept. decay)

!           if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(abs(TopDecays).ne.1) call Error("TopDecays needs to be 1!")
          NumHistograms = 4
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          Histo(1)%Info   = "r"
          Histo(1)%NBins  = 20
          Histo(1)%BinSize= 0.05d0
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 1d0

          Histo(2)%Info    = "pT_lepMinus"
          Histo(2)%NBins   = 40
          Histo(2)%BinSize = 50d0*GeV
          Histo(2)%LowVal  = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "DeltaPhi"
          Histo(4)%NBins  = 65
          Histo(4)%BinSize= 0.1d0
          Histo(4)%LowVal = -3d0
          Histo(4)%SetScale= 1d0

ELSEIF( ObsSet.EQ.9 ) THEN! this is for the factorization check
          NumHistograms = 6
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0


ELSEIF( ObsSet.EQ.10 ) THEN! set of observables for ttbjet production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.11 ) THEN! set of observables for ttbjet production without decays at LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.12 ) THEN! set of observables for ttbjet production as signal process at Tevatron (hadr.Atop, lept.top decay)
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          NumHistograms = 44
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
!          if( .not.allocated(Histo2D) ) then
!                allocate( Histo2D(1:3), stat=AllocStatus  )
!                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo2D")
!          endif


          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 20
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 20
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pt_5th_jet"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_5th_jet"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 0.4d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_miss"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "HT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 25d0*GeV
          Histo(8)%LowVal = 100d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m_lb"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "phi(l,b)"
          Histo(10)%NBins  = 15
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "R(l,b)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta_FBlepPlus"
          Histo(12)%NBins  = 2
          Histo(12)%BinSize= 5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "pT(ttbar)"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 10d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0

! new stuff for A_FB analysis, this is for ideally reconstructed tops
          Histo(14)%Info   = "pT(ttbar) FWD"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 10d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "pT(ttbar) BWD"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 10d0*GeV
          Histo(15)%LowVal =  0d0*GeV
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "m(ttbar) FWD"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 50d0*GeV
          Histo(16)%LowVal =  300d0*GeV
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "m(ttbar) BWD"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 50d0*GeV
          Histo(17)%LowVal =  300d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "y(ttbar) FWD"
          Histo(18)%NBins  = 40
          Histo(18)%BinSize= 0.5d0
          Histo(18)%LowVal =-5.0d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "y(ttbar) BWD"
          Histo(19)%NBins  = 40
          Histo(19)%BinSize= 0.5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "dy(tops) FWD"
          Histo(20)%NBins  = 40
          Histo(20)%BinSize= 0.25d0
          Histo(20)%LowVal = 0.0d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "dy(tops) BWD"
          Histo(21)%NBins  = 40
          Histo(21)%BinSize= 0.25d0
          Histo(21)%LowVal = 0.0d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "pT(ttbar) FWD lept"
          Histo(22)%NBins  = 50
          Histo(22)%BinSize= 10d0*GeV
          Histo(22)%LowVal =  0d0*GeV
          Histo(22)%SetScale= 100d0

          Histo(23)%Info   = "pT(ttbar) BWD lept"
          Histo(23)%NBins  = 50
          Histo(23)%BinSize= 10d0*GeV
          Histo(23)%LowVal =  0d0*GeV
          Histo(23)%SetScale= 100d0

          Histo(24)%Info   = "y_FB(top)"
          Histo(24)%NBins  = 2
          Histo(24)%BinSize= 5d0
          Histo(24)%LowVal =-5.0d0
          Histo(24)%SetScale= 1d0

          Histo(25)%Info   = "pT(top)"
          Histo(25)%NBins  = 40
          Histo(25)%BinSize= 20d0*GeV
          Histo(25)%LowVal = 0d0
          Histo(25)%SetScale= 100d0

          Histo(26)%Info   = "phi(ttbar) FWD" 
          Histo(26)%NBins  = 15
          Histo(26)%BinSize= 0.25d0
          Histo(26)%LowVal = 0d0
          Histo(26)%SetScale= 1d0

          Histo(27)%Info   = "phi(ttbar) BWD" 
          Histo(27)%NBins  = 15
          Histo(27)%BinSize= 0.25d0
          Histo(27)%LowVal = 0d0
          Histo(27)%SetScale= 1d0


! new stuff for A_FB analysis, this is for realistically reconstructed tops
          Histo(28)%Info   = "pT(ttbar)"
          Histo(28)%NBins  = 50
          Histo(28)%BinSize= 10d0*GeV
          Histo(28)%LowVal = 0d0
          Histo(28)%SetScale= 100d0

          Histo(29)%Info   = "pT(ttbar) FWD"
          Histo(29)%NBins  = 50
          Histo(29)%BinSize= 10d0*GeV
          Histo(29)%LowVal =  0d0*GeV
          Histo(29)%SetScale= 100d0

          Histo(30)%Info   = "pT(ttbar) BWD"
          Histo(30)%NBins  = 50
          Histo(30)%BinSize= 10d0*GeV
          Histo(30)%LowVal =  0d0*GeV
          Histo(30)%SetScale= 100d0

          Histo(31)%Info   = "m(ttbar) FWD"
          Histo(31)%NBins  = 50
          Histo(31)%BinSize= 50d0*GeV
          Histo(31)%LowVal =  300d0*GeV
          Histo(31)%SetScale= 100d0

          Histo(32)%Info   = "m(ttbar) BWD"
          Histo(32)%NBins  = 50
          Histo(32)%BinSize= 50d0*GeV
          Histo(32)%LowVal =  300d0*GeV
          Histo(32)%SetScale= 100d0

          Histo(33)%Info   = "y(ttbar) FWD"
          Histo(33)%NBins  = 40
          Histo(33)%BinSize= 0.5d0
          Histo(33)%LowVal =-5.0d0
          Histo(33)%SetScale= 1d0

          Histo(34)%Info   = "y(ttbar) BWD"
          Histo(34)%NBins  = 40
          Histo(34)%BinSize= 0.5d0
          Histo(34)%LowVal =-5.0d0
          Histo(34)%SetScale= 1d0

          Histo(35)%Info   = "dy(tops) FWD"
          Histo(35)%NBins  = 40
          Histo(35)%BinSize= 0.25d0
          Histo(35)%LowVal = 0.0d0
          Histo(35)%SetScale= 1d0

          Histo(36)%Info   = "dy(tops) BWD"
          Histo(36)%NBins  = 40
          Histo(36)%BinSize= 0.25d0
          Histo(36)%LowVal = 0.0d0
          Histo(36)%SetScale= 1d0

          Histo(37)%Info   = "pT(ttbar) FWD lept"
          Histo(37)%NBins  = 50
          Histo(37)%BinSize= 10d0*GeV
          Histo(37)%LowVal =  0d0*GeV
          Histo(37)%SetScale= 100d0

          Histo(38)%Info   = "pT(ttbar) BWD lept"
          Histo(38)%NBins  = 50
          Histo(38)%BinSize= 10d0*GeV
          Histo(38)%LowVal =  0d0*GeV
          Histo(38)%SetScale= 100d0

          Histo(39)%Info   = "y_FB(top)"
          Histo(39)%NBins  = 2
          Histo(39)%BinSize= 5d0
          Histo(39)%LowVal =-5.0d0
          Histo(39)%SetScale= 1d0

          Histo(40)%Info   = "pT(top)"
          Histo(40)%NBins  = 40
          Histo(40)%BinSize= 20d0*GeV
          Histo(40)%LowVal = 0d0
          Histo(40)%SetScale= 100d0

          Histo(41)%Info   = "phi(ttbar) FWD" 
          Histo(41)%NBins  = 15
          Histo(41)%BinSize= 0.25d0
          Histo(41)%LowVal = 0d0
          Histo(41)%SetScale= 1d0

          Histo(42)%Info   = "phi(ttbar) BWD" 
          Histo(42)%NBins  = 15
          Histo(42)%BinSize= 0.25d0
          Histo(42)%LowVal = 0d0
          Histo(42)%SetScale= 1d0

          Histo(43)%Info   = "cos(theta_top*)" 
          Histo(43)%NBins  = 20
          Histo(43)%BinSize= 0.1d0
          Histo(43)%LowVal = -1d0
          Histo(43)%SetScale= 1d0

          Histo(44)%Info   = "beta_top*cos(theta_top*)" 
          Histo(44)%NBins  = 20
          Histo(44)%BinSize= 0.1d0
          Histo(44)%LowVal = -1d0
          Histo(44)%SetScale= 1d0

! 2D histograms for realistically reconstructed tops      
!
!          Histo2D(1)%Info   = "m_ttbar over pT_jet"  
!          Histo2D(1)%NBins(1)  = 50
!          Histo2D(1)%BinSize(1)= 50d0
!          Histo2D(1)%LowVal(1) = 300d0
!          Histo2D(1)%SetScale(1)= 100d0
!          Histo2D(1)%NBins(2)  = 50
!          Histo2D(1)%BinSize(2)= 10d0
!          Histo2D(1)%LowVal(2) = 0d0
!          Histo2D(1)%SetScale(2)= 100d0
!
!          Histo2D(2)%Info   = "m_ttbar over pT_ttbar FWD"  
!          Histo2D(2)%NBins(1)  = 50
!          Histo2D(2)%BinSize(1)= 50d0
!          Histo2D(2)%LowVal(1) = 300d0
!          Histo2D(2)%SetScale(1)= 100d0
!          Histo2D(2)%NBins(2)  = 50
!          Histo2D(2)%BinSize(2)= 20d0
!          Histo2D(2)%LowVal(2) = 0d0
!          Histo2D(2)%SetScale(2)= 100d0
!
!          Histo2D(3)%Info   = "m_ttbar over pT_ttbar BWD"  
!          Histo2D(3)%NBins(1)  = 50
!          Histo2D(3)%BinSize(1)= 50d0
!          Histo2D(3)%LowVal(1) = 300d0
!          Histo2D(3)%SetScale(1)= 100d0
!          Histo2D(3)%NBins(2)  = 50
!          Histo2D(3)%BinSize(2)= 20d0
!          Histo2D(3)%LowVal(2) = 0d0
!          Histo2D(3)%SetScale(2)= 100d0





ELSEIF( ObsSet.EQ.13 ) THEN! set of observables for ttbjet production as signal process at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1) call Error("TopDecays needs to be 1!")
          if(Collider_Energy.ne.7000d0*GeV) call Error("wrong collider energy for ObsSet=13!")
          NumHistograms = 16
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 20d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "HT"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 100d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "Minv_leptons"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT_jet"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_jet"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_miss"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "pT_ttbar"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 15d0*GeV
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "phi(l+,l-)"
          Histo(11)%NBins  = 15
          Histo(11)%BinSize= 0.25d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "m_l+l-"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "pT_bjet"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "eta_bjet"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 0.4d0
          Histo(14)%LowVal =-5.0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "m_bb"
          Histo(15)%NBins  = 50
          Histo(15)%BinSize= 20d0*GeV
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 100d0

          Histo(16)%Info   = "m_bj"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 20d0*GeV
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 100d0




ELSEIF( ObsSet.EQ.14 ) THEN! set of observables for ttbjet production as background process to VBF at LHC (di-lept. decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1) call Error("TopDecays needs to be 1!")
          if(Collider_Energy.ne.14000d0*GeV) call Error("wrong collider energy for ObsSet=14!")
          NumHistograms = 11
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info    = "pT_lepMinus"
          Histo(1)%NBins   = 40
          Histo(1)%BinSize = 50d0*GeV
          Histo(1)%LowVal  = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepMinus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_lepPlus"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_lepPlus"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "Minv_leptons"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT_jet"
          Histo(6)%NBins  = 40
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "phi_leptons"
          Histo(7)%NBins  = 20
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "mT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta(hardest jet)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "eta(2nd hardest jet)"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "eta_Zeppi"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0



ELSEIF( ObsSet.EQ.15 ) THEN! set of observables for ttbjet production as signal process at LHC (hadr.Atop, lept.top decay)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4) call Error("TopDecays needs to be 4!")
          if(Collider_Energy.ne.7000d0*GeV) call Error("wrong collider energy for ObsSet=15!")
          NumHistograms = 13
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pt_5th_jet"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_5th_jet"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.4d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_miss"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "HT"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 100d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m_lb"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "phi(l,b)"
          Histo(10)%NBins  = 15
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal = 0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "R(l,b)"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta_FBlepPlus"
          Histo(12)%NBins  = 2
          Histo(12)%BinSize= 5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "pT_ttbar"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 15d0*GeV
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 100d0





ELSEIF( ObsSet.EQ.19 ) THEN! for checks of ttbjet
          NumHistograms = 10
          if(TopDecays.eq.0) call Error("TopDecays needs to be >0!")
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lepPlus"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_lepPlus"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_jet"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_jet"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.4d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_miss"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "HT"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "m_lb"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "phi(l,b)"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "R(l,b)"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal = 0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "eta_FBlepPlus"
          Histo(10)%NBins  = 2
          Histo(10)%BinSize= 5d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0


ELSEIF( ObsSet.EQ.20 ) THEN! set of observables for ttbgamma production without decays at Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "etaFB_CP"
          Histo(9)%NBins  = 2
          Histo(9)%BinSize= 10d0
          Histo(9)%LowVal =-10.0d0
          Histo(9)%SetScale= 0.1d0


ELSEIF( ObsSet.EQ.21 ) THEN! set of observables for ttbgamma production without decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "etaFB_CP"
          Histo(9)%NBins  = 2
          Histo(9)%BinSize= 10d0
          Histo(9)%LowVal =-10.0d0
          Histo(9)%SetScale= 0.1d0


ELSEIF( ObsSet.EQ.22 ) THEN! set of observables for ttbgamma production di-lept. decays at the Tevatron
          if(Collider.ne.2)  call Error("Collider needs to be TEV!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal = 20d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "ET_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 20d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho+miss)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "m(LepP+bjet)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "phi(LepP,LepM)"
          Histo(14)%NBins  = 20
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal = 0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "etaFB_CP"
          Histo(15)%NBins  = 2
          Histo(15)%BinSize= 10d0
          Histo(15)%LowVal =-10.0d0
          Histo(15)%SetScale= 0.1d0


ELSEIF( ObsSet.EQ.23 ) THEN! set of observables for ttbgamma production di-lept. decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1 ) call Error("TopDecays needs to be 1!")
          NumHistograms = 14
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 200
          Histo(7)%BinSize= 5d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 2d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "ET_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal = 0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho+miss)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 100d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "m(LepP+bjet)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 20d0*GeV
          Histo(13)%LowVal = 20d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "phi(LepP,LepM)"
          Histo(14)%NBins  = 20
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal = 0d0
          Histo(14)%SetScale= 1d0




ELSEIF( ObsSet.EQ.24 ) THEN! set of observables for ttbgamma production semi-lept. decays at the TEV
          if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          NumHistograms = 16
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_ATop"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 50d0*GeV
          Histo(3)%LowVal = 0d0
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 25d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 25d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 50d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 90
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "etaFB_CP"
          Histo(16)%NBins  = 2
          Histo(16)%BinSize= 10d0
          Histo(16)%LowVal =-10.0d0
          Histo(16)%SetScale= 0.1d0




ELSEIF( ObsSet.EQ.25 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4 .and. TopDecays.ne.3) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
          NumHistograms = 15
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0
                    
          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "etaFB_ATop"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 5d0
          Histo(5)%LowVal =-5.0d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 200
          Histo(7)%BinSize= 5d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 20d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 20d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal = 100d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 90
          Histo(14)%BinSize= 5d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0



! ELSEIF( ObsSet.EQ.26 ) THEN! set of observables for ttbgamma production semi-lept. decays at the TEV
!           if(Collider.ne.2)  call Error("Collider needs to be Tevatron!")
!           if(TopDecays.ne.4  .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
!           NumHistograms = 15
!           if( .not.allocated(Histo) ) then
!                 allocate( Histo(1:NumHistograms), stat=AllocStatus  )
!                 if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!           endif
! 
!           Histo(1)%Info   = "pT_ATop"
!           Histo(1)%NBins  = 40
!           Histo(1)%BinSize= 50d0*GeV
!           Histo(1)%LowVal = 0d0
!           Histo(1)%SetScale= 100d0
! 
!           Histo(2)%Info   = "eta_ATop"
!           Histo(2)%NBins  = 40
!           Histo(2)%BinSize= 0.25d0
!           Histo(2)%LowVal =-5.0d0
!           Histo(2)%SetScale= 1d0
! 
!           Histo(3)%Info   = "pT_Top"
!           Histo(3)%NBins  = 40
!           Histo(3)%BinSize= 50d0*GeV
!           Histo(3)%LowVal = 0d0
!           Histo(3)%SetScale= 100d0
! 
!           Histo(4)%Info   = "eta_Top"
!           Histo(4)%NBins  = 40
!           Histo(4)%BinSize= 0.25d0
!           Histo(4)%LowVal =-5.0d0
!           Histo(4)%SetScale= 1d0
! 
!           Histo(5)%Info   = "etaFB_ATop"
!           Histo(5)%NBins  = 2
!           Histo(5)%BinSize= 5d0
!           Histo(5)%LowVal =-5.0d0
!           Histo(5)%SetScale= 1d0
! 
!           Histo(6)%Info   = "etaFB_Top"
!           Histo(6)%NBins  = 2
!           Histo(6)%BinSize= 5d0
!           Histo(6)%LowVal =-5.0d0
!           Histo(6)%SetScale= 1d0
! 
!           Histo(7)%Info   = "pT_Photon"
!           Histo(7)%NBins  = 40
!           Histo(7)%BinSize= 20d0*GeV
!           Histo(7)%LowVal = 20d0*GeV
!           Histo(7)%SetScale= 100d0
! 
!           Histo(8)%Info   = "eta_Photon"
!           Histo(8)%NBins  = 40
!           Histo(8)%BinSize= 0.25d0
!           Histo(8)%LowVal =-5.0d0
!           Histo(8)%SetScale= 1d0
! 
!           Histo(9)%Info   = "pT_LepP"
!           Histo(9)%NBins  = 40
!           Histo(9)%BinSize= 20d0*GeV
!           Histo(9)%LowVal = 20d0*GeV
!           Histo(9)%SetScale= 100d0
! 
!           Histo(10)%Info   = "eta_LepP"
!           Histo(10)%NBins  = 40
!           Histo(10)%BinSize= 0.25d0
!           Histo(10)%LowVal =-5.0d0
!           Histo(10)%SetScale= 1d0
! 
!           Histo(11)%Info   = "pT_miss"
!           Histo(11)%NBins  = 40
!           Histo(11)%BinSize= 20d0*GeV
!           Histo(11)%LowVal = 20d0*GeV
!           Histo(11)%SetScale= 100d0
! 
!           Histo(12)%Info   = "HT(jets+lept+pho)"
!           Histo(12)%NBins  = 40
!           Histo(12)%BinSize= 20d0*GeV
!           Histo(12)%LowVal = 150d0*GeV
!           Histo(12)%SetScale= 100d0
! 
!           Histo(13)%Info   = "R(pho,bjet)"
!           Histo(13)%NBins  = 25
!           Histo(13)%BinSize= 0.2d0
!           Histo(13)%LowVal = 0d0
!           Histo(13)%SetScale= 1d0
! 
!           Histo(14)%Info   = "m(lep+bjet)"
!           Histo(14)%NBins  = 40
!           Histo(14)%BinSize= 20d0*GeV
!           Histo(14)%LowVal = 20d0*GeV
!           Histo(14)%SetScale= 100d0
! 
!           Histo(15)%Info   = "phi(photon,lept)"
!           Histo(15)%NBins  = 20
!           Histo(15)%BinSize= 0.2d0
!           Histo(15)%LowVal = 0d0
!           Histo(15)%SetScale= 1d0


! ELSEIF( ObsSet.EQ.25 .OR. ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28 ) THEN
ELSEIF( ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC for Q_top measurement
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if( (ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28) .and. TopDecays.ne.4) call Error("TopDecays needs to be 4 for ObsSet=26,27,28")
          if( (ObsSet.EQ.25) .and. TopDecays.ne.3) call Error("TopDecays needs to be 3 for ObsSet=25")
          if( (ObsSet.EQ.26 .OR. ObsSet.EQ.27 .OR. ObsSet.EQ.28) .and. (Q_top.ne.Q_up .and. TopDecays.ne.3) ) call Error("TopDecays needs to be 3 for Qt=-4/3")
          NumHistograms = 20
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "(pT_Top+pT_ATop)/2"
          Histo(1)%NBins  = 4
          Histo(1)%BinSize= 60d0*GeV
          Histo(1)%LowVal = 40d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "(y_ATop+y_Top)/2"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT_Top"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 40
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0
          
          Histo(5)%Info   = "etaFB_CP"
          Histo(5)%NBins  = 2
          Histo(5)%BinSize= 10d0
          Histo(5)%LowVal =-10.0d0
          Histo(5)%SetScale= 0.1d0

          Histo(6)%Info   = "etaFB_Top"
          Histo(6)%NBins  = 2
          Histo(6)%BinSize= 5d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "pT_Photon"
          Histo(7)%NBins  = 40
          Histo(7)%BinSize= 25d0*GeV
          Histo(7)%LowVal = 0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "eta_Photon"
          Histo(8)%NBins  = 40
          Histo(8)%BinSize= 0.25d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "pT_LepP"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "eta_LepP"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 0.25d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "pT_miss"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 25d0*GeV
          Histo(11)%LowVal =  0d0*GeV
          Histo(11)%SetScale= 100d0

          Histo(12)%Info   = "HT(jets+lept+pho)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 50d0*GeV
          Histo(12)%LowVal = 150d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "R(pho,bjet)"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.25d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "m(lep+bjet)"
          Histo(14)%NBins  = 40
          Histo(14)%BinSize= 10d0*GeV
          Histo(14)%LowVal = 20d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "phi(photon,lept)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal = 0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "mjjb"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 4d0*GeV
          Histo(16)%LowVal = 70d0*GeV
          Histo(16)%SetScale= 100d0

          Histo(17)%Info   = "mjj"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 4d0*GeV
          Histo(17)%LowVal = 50d0*GeV
          Histo(17)%SetScale= 100d0

          Histo(18)%Info   = "mTbP"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 4d0*GeV
          Histo(18)%LowVal = 70d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "mTlp"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 4d0*GeV
          Histo(19)%LowVal = 70d0*GeV
          Histo(19)%SetScale= 100d0

          Histo(20)%Info   = "<Ehat>"
          Histo(20)%NBins  = 1
          Histo(20)%BinSize= 10000d0*GeV
          Histo(20)%LowVal = 0d0
          Histo(20)%SetScale=0.01d0
          
! ELSEIF( ObsSet.EQ.28 ) THEN! set of observables for ttbgamma production semi-lept. decays at the LHC for Q_top measurement
!           if(Collider.ne.1)  call Error("Collider needs to be LHC!")
!           if(TopDecays.ne.4 .and. TopDecays.ne.3) call Error("TopDecays needs to be 3 (for Qt=-4/3) or 4 (for Qt=Qup)")
!           NumHistograms = 15
!           if( .not.allocated(Histo) ) then
!                 allocate( Histo(1:NumHistograms), stat=AllocStatus  )
!                 if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!           endif
! 
!           Histo(1)%Info   = "pT_ATop"
!           Histo(1)%NBins  = 40
!           Histo(1)%BinSize= 50d0*GeV
!           Histo(1)%LowVal = 0d0
!           Histo(1)%SetScale= 100d0
! 
!           Histo(2)%Info   = "eta_ATop"
!           Histo(2)%NBins  = 40
!           Histo(2)%BinSize= 0.25d0
!           Histo(2)%LowVal =-5.0d0
!           Histo(2)%SetScale= 1d0
! 
!           Histo(3)%Info   = "pT_Top"
!           Histo(3)%NBins  = 40
!           Histo(3)%BinSize= 50d0*GeV
!           Histo(3)%LowVal = 0d0
!           Histo(3)%SetScale= 100d0
! 
!           Histo(4)%Info   = "eta_Top"
!           Histo(4)%NBins  = 40
!           Histo(4)%BinSize= 0.25d0
!           Histo(4)%LowVal =-5.0d0
!           Histo(4)%SetScale= 1d0
! 
!           Histo(5)%Info   = "etaFB_ATop"
!           Histo(5)%NBins  = 2
!           Histo(5)%BinSize= 5d0
!           Histo(5)%LowVal =-5.0d0
!           Histo(5)%SetScale= 1d0
! 
!           Histo(6)%Info   = "etaFB_Top"
!           Histo(6)%NBins  = 2
!           Histo(6)%BinSize= 5d0
!           Histo(6)%LowVal =-5.0d0
!           Histo(6)%SetScale= 1d0
! 
!           Histo(7)%Info   = "pT_Photon"
!           Histo(7)%NBins  = 40
!           Histo(7)%BinSize= 25d0*GeV
!           Histo(7)%LowVal = 0d0*GeV
!           Histo(7)%SetScale= 100d0
! 
!           Histo(8)%Info   = "eta_Photon"
!           Histo(8)%NBins  = 40
!           Histo(8)%BinSize= 0.25d0
!           Histo(8)%LowVal =-5.0d0
!           Histo(8)%SetScale= 1d0
! 
!           Histo(9)%Info   = "pT_LepP"
!           Histo(9)%NBins  = 40
!           Histo(9)%BinSize= 25d0*GeV
!           Histo(9)%LowVal =  0d0*GeV
!           Histo(9)%SetScale= 100d0
! 
!           Histo(10)%Info   = "eta_LepP"
!           Histo(10)%NBins  = 40
!           Histo(10)%BinSize= 0.25d0
!           Histo(10)%LowVal =-5.0d0
!           Histo(10)%SetScale= 1d0
! 
!           Histo(11)%Info   = "pT_miss"
!           Histo(11)%NBins  = 40
!           Histo(11)%BinSize= 25d0*GeV
!           Histo(11)%LowVal =  0d0*GeV
!           Histo(11)%SetScale= 100d0
! 
!           Histo(12)%Info   = "HT(jets+lept+pho)"
!           Histo(12)%NBins  = 40
!           Histo(12)%BinSize= 50d0*GeV
!           Histo(12)%LowVal = 150d0*GeV
!           Histo(12)%SetScale= 100d0
! 
!           Histo(13)%Info   = "R(pho,bjet)"
!           Histo(13)%NBins  = 25
!           Histo(13)%BinSize= 0.25d0
!           Histo(13)%LowVal = 0d0
!           Histo(13)%SetScale= 1d0
! 
!           Histo(14)%Info   = "m(lep+bjet)"
!           Histo(14)%NBins  = 90
!           Histo(14)%BinSize= 5d0*GeV
!           Histo(14)%LowVal = 20d0*GeV
!           Histo(14)%SetScale= 100d0
! 
!           Histo(15)%Info   = "phi(photon,lept)"
!           Histo(15)%NBins  = 20
!           Histo(15)%BinSize= 0.2d0
!           Histo(15)%LowVal = 0d0
!           Histo(15)%SetScale= 1d0


ELSEIF( ObsSet.EQ.29 ) THEN! set of observables for ttbgamma production without decays at Tevatron
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 25d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "eta_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "eta_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "pT_Photon"
          Histo(5)%NBins  = 40
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "eta_Photon"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.2d0
          Histo(6)%LowVal =-5.0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "etaFB_ATop"
          Histo(7)%NBins  = 2
          Histo(7)%BinSize= 5d0
          Histo(7)%LowVal =-5.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "etaFB_Top"
          Histo(8)%NBins  = 2
          Histo(8)%BinSize= 5d0
          Histo(8)%LowVal =-5.0d0
          Histo(8)%SetScale= 1d0



ELSEIF( ObsSet.EQ.31 ) THEN! set of observables for HTHTbar + A0 production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(XTopDecays.ne.0)  call Error("XTopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 25d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0



ELSEIF( ObsSet.EQ.32 .OR. ObsSet.EQ.34 ) THEN! set of observables for HTHTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          if(XTopDecays.ne.1 .and. XTopDecays.ne.2 ) call Error("XTopDecays needs to be 1(BH) or 2(A0)!")
          NumHistograms = 7
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "phi(l+,l-)"
          Histo(4)%NBins  = 15
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "m_l+l-"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0


ELSEIF( ObsSet.EQ.33 .OR. ObsSet.EQ.35 .OR. ObsSet.EQ.36 ) THEN! set of observables for HTHTbar + A0/BH production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          if(XTopDecays.ne.1 ) call Error("XTopDecays needs to be 1!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "HT"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 50d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(softest jet)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 40d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT_Top"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta_Top"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.25d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "Log(S)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 1.0d0
          Histo(10)%LowVal =-10.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "r_pT"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-2.0d0
          Histo(11)%SetScale= 1d0


          Histo(12)%Info   = "m(lep+bjet)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 10d0*GeV
          Histo(12)%LowVal = 20d0*GeV
          Histo(12)%SetScale= 100d0

   

ELSEIF( ObsSet.EQ.41 ) THEN! set of observables for STSTbar (stable stops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(XTopDecays.ne.0)  call Error("XTopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0



ELSEIF( ObsSet.EQ.42 .OR. ObsSet.EQ.44 ) THEN! set of observables for STSTbar + Chi production (di-lept. tops)

          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          if(XTopDecays.ne.3  ) call Error("XTopDecays needs to be 3!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 40d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "phi(l+,l-)"
          Histo(4)%NBins  = 15
          Histo(4)%BinSize= 0.25d0
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "m_l+l-"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "HT"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT(W)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 40d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "MT^eff"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 100d0*GeV
          Histo(8)%LowVal = 0d0
          Histo(8)%SetScale= 100d0


ELSEIF( ObsSet.EQ.43 .OR. ObsSet.EQ.45  .OR. ObsSet.EQ.46 ) THEN! set of observables for STSTbar + Chi production
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          if(XTopDecays.ne.3 ) call Error("XTopDecays needs to be 3!")
          NumHistograms = 12
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_LepP"
          Histo(1)%NBins  = 40
          Histo(1)%BinSize= 40d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta_LepP"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "ET_miss"
          Histo(3)%NBins  = 40
          Histo(3)%BinSize= 40d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "HT"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 50d0*GeV
          Histo(4)%LowVal = 0d0
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(softest jet)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "MT(W)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 40d0*GeV
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "MT^eff"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 50d0*GeV
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT_Top"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 50d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "eta_Top"
          Histo(9)%NBins  = 40
          Histo(9)%BinSize= 0.25d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "Log(S)"
          Histo(10)%NBins  = 40
          Histo(10)%BinSize= 1.0d0
          Histo(10)%LowVal =-10.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "r_pT"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal =-2.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "m(lep+bjet)"
          Histo(12)%NBins  = 40
          Histo(12)%BinSize= 10d0*GeV
          Histo(12)%LowVal = 20d0*GeV
          Histo(12)%SetScale= 100d0

   

ELSEIF( ObsSet.EQ.48 ) THEN! set of observables for STOP width
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          if(XTopDecays.ne.3 .and. XTopDecays.ne.1) call Error("XTopDecays needs to be 1 or 3!")
!           NumHistograms = 0
!           if( .not.allocated(Histo) ) then
!                 allocate( Histo(1:NumHistograms), stat=AllocStatus  )
!                 if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
!           endif



ELSEIF( ObsSet.EQ.51 ) THEN! set of observables for ttb+Z (stable tops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0)  call Error("TopDecays needs to be 0")
          NumHistograms = 1
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(top)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0





ELSEIF( ObsSet.EQ.52 .or. ObsSet.EQ.55 ) THEN! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )
          if(abs(TopDecays).ne.1)  call Error("TopDecays needs to be 1")
          if(abs(ZDecays).ne.1)    call Error("ZDecays needs to be 1")
          NumHistograms = 30
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(lep+)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT(lep-))"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal =  0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT(Z)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(top)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(Atop)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT(j1)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT(j2)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT(mu+)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal =  0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "pT(e-)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal =  0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "pT(miss)"
          Histo(10)%NBins  = 50
          Histo(10)%BinSize= 10d0*GeV
          Histo(10)%LowVal =  0d0*GeV
          Histo(10)%SetScale= 100d0

          Histo(11)%Info   = "eta(lep+)"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.5d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "eta(lep-)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.5d0
          Histo(12)%LowVal =-5.0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "eta(Z)"
          Histo(13)%NBins  = 30
          Histo(13)%BinSize= 0.2d0
          Histo(13)%LowVal =-3.0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "eta(top)"
          Histo(14)%NBins  = 30
          Histo(14)%BinSize= 0.2d0
          Histo(14)%LowVal =-3.0d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "eta(atop)"
          Histo(15)%NBins  = 30
          Histo(15)%BinSize= 0.2d0
          Histo(15)%LowVal =-3.0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "eta(j1)"
          Histo(16)%NBins  = 50
          Histo(16)%BinSize= 0.5d0
          Histo(16)%LowVal =-5.0d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "eta(j2)"
          Histo(17)%NBins  = 50
          Histo(17)%BinSize= 0.5d0
          Histo(17)%LowVal =-5.0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "eta(mu+)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 0.5d0
          Histo(18)%LowVal =-5.0d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "eta(e-)"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 0.5d0
          Histo(19)%LowVal =-5.0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "phi(l+,l-)"
          Histo(20)%NBins  = 15    *4d0
          Histo(20)%BinSize= 0.25d0/4d0
          Histo(20)%LowVal = 0d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "pseudo(Z)"
          Histo(21)%NBins  = 30
          Histo(21)%BinSize= 0.2d0
          Histo(21)%LowVal =-3.0d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "pseudo(top)"
          Histo(22)%NBins  = 30
          Histo(22)%BinSize= 0.2d0
          Histo(22)%LowVal =-3.0d0
          Histo(22)%SetScale= 1d0

          Histo(23)%Info   = "pseudo(atop)"
          Histo(23)%NBins  = 30
          Histo(23)%BinSize= 0.2d0
          Histo(23)%LowVal =-3.0d0
          Histo(23)%SetScale= 1d0

          Histo(24)%Info   = "Minv_lep_bquark"
          Histo(24)%NBins  = 40
          Histo(24)%BinSize= 10d0*GeV
          Histo(24)%LowVal = 0d0
          Histo(24)%SetScale= 100d0

          Histo(25)%Info   = "m_T2"
          Histo(25)%NBins  = 50
          Histo(25)%BinSize= 10d0*GeV
          Histo(25)%LowVal = 0d0*GeV
          Histo(25)%SetScale= 100d0

          Histo(26)%Info   = "<m(lep+bjet)>"
          Histo(26)%NBins  = 2
          Histo(26)%BinSize= 10000d0*GeV
          Histo(26)%LowVal = 0d0*GeV
          Histo(26)%SetScale= 1d0/Histo(26)%BinSize

          Histo(27)%Info   = "<mT2>"
          Histo(27)%NBins  = 2
          Histo(27)%BinSize= 10000d0*GeV
          Histo(27)%LowVal = 0d0*GeV
          Histo(27)%SetScale= 1d0/Histo(27)%BinSize

          Histo(28)%Info   = "sqrt(shat)"
          Histo(28)%NBins  = 10
          Histo(28)%BinSize= 100d0*GeV
          Histo(28)%LowVal =  350d0*GeV
          Histo(28)%SetScale= 100d0

          Histo(29)%Info   = "pT(top)"
          Histo(29)%NBins  = 10
          Histo(29)%BinSize= 60d0*GeV
          Histo(29)%LowVal =  0d0*GeV
          Histo(29)%SetScale= 100d0

          Histo(30)%Info   = "alpha(t,tbar)"
          Histo(30)%NBins  = 10
          Histo(30)%BinSize= 0.31d0
          Histo(30)%LowVal =  0d0*GeV
          Histo(30)%SetScale= 1d0


ELSEIF( ObsSet.EQ.53 .or. ObsSet.EQ.56 .or. ObsSet.EQ.58  ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay )
          if(abs(TopDecays).ne.4)  call Error("TopDecays needs to be 4")
          if(abs(ZDecays).ne.1 .and. abs(ZDecays).ne.11)    call Error("ZDecays needs to be 1 or 11")
          NumHistograms = 28
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

!           Histo(1)%Info   = "pT(lep+)"
!           Histo(1)%NBins  = 50
!           Histo(1)%BinSize= 20d0*GeV
!           Histo(1)%LowVal =  0d0*GeV
!           Histo(1)%SetScale= 100d0
! 
!           Histo(2)%Info   = "pT(mu-)"
!           Histo(2)%NBins  = 50
!           Histo(2)%BinSize= 20d0*GeV
!           Histo(2)%LowVal =  0d0*GeV
!           Histo(2)%SetScale= 100d0

          Histo(1)%Info   = "(pT_Top+pT_ATop)/2"
          Histo(1)%NBins  = 4
          Histo(1)%BinSize= 60d0*GeV
          Histo(1)%LowVal = 40d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "(y_ATop+y_Top)/2"
          Histo(2)%NBins  = 40
          Histo(2)%BinSize= 0.25d0
          Histo(2)%LowVal =-5.0d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT(mu+)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(b1)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(b2)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT(j1)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT(j2)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT(miss)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal =  0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "y(lep+)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 0.5d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y(mu-)"
          Histo(10)%NBins  = 20
          Histo(10)%BinSize= 0.5d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "y(mu+)"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.5d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "pT(Z)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal =  0d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "y(Z)"
          Histo(13)%NBins  = 20
          Histo(13)%BinSize= 0.5d0
          Histo(13)%LowVal =-5.0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "pT(top)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 20d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "y(top)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.5d0
          Histo(15)%LowVal =-5.0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "yFB(Top)"
          Histo(16)%NBins  = 2
          Histo(16)%BinSize= 5d0
          Histo(16)%LowVal =-5.0d0
          Histo(16)%SetScale= 1d0

!           Histo(16)%Info   = "y(antitop)"
!           Histo(16)%NBins  = 20
!           Histo(16)%BinSize= 0.5d0
!           Histo(16)%LowVal =-5.0d0
!           Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "phi(mu-,mu+)"
          Histo(17)%NBins  = 26
          Histo(17)%BinSize= 0.25d0/2d0
          Histo(17)%LowVal = 0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "CosAlpha*(Z,mu-)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 2d0/50d0 
          Histo(18)%LowVal = -1d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "phi(Z,antit)"
          Histo(19)%NBins  = 26
          Histo(19)%BinSize= 0.25d0/2d0
          Histo(19)%LowVal = 0d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "phi(t,antit)"
          Histo(20)%NBins  = 26
          Histo(20)%BinSize= 0.25d0/2d0
          Histo(20)%LowVal = 0d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "N(jets)"
          Histo(21)%NBins  = 5
          Histo(21)%BinSize= 1d0
          Histo(21)%LowVal = 1d0
          Histo(21)%SetScale= 1d0

          ! distributions 17..20 but with bin smearing and double the bin size
          Histo(22)%Info   = "phi(mu-,mu+) smeared"
          Histo(22)%NBins  = 26
          Histo(22)%BinSize= 0.25d0/2d0
          Histo(22)%LowVal = 0d0
          Histo(22)%SetScale= 1d0
          Histo(22)%BinSmearing=.true.
          Histo(22)%SmearSigma=1d0

          Histo(23)%Info   = "CosAlpha*(Z,mu-) smeared"
          Histo(23)%NBins  = 50
          Histo(23)%BinSize= 2d0/50d0
          Histo(23)%LowVal = -1d0
          Histo(23)%SetScale= 1d0
          Histo(22)%BinSmearing=.true.
          Histo(22)%SmearSigma=1d0

          Histo(24)%Info   = "phi(Z,antit) smeared"
          Histo(24)%NBins  = 26
          Histo(24)%BinSize= 0.25d0/2d0
          Histo(24)%LowVal = 0d0
          Histo(24)%SetScale= 1d0
          Histo(24)%BinSmearing=.true.
          Histo(24)%SmearSigma=1d0

          Histo(25)%Info   = "phi(t,antit) smeared"
          Histo(25)%NBins  = 26
          Histo(25)%BinSize= 0.25d0/2d0
          Histo(25)%LowVal = 0d0
          Histo(25)%SetScale= 1d0
          Histo(25)%BinSmearing=.true.
          Histo(25)%SmearSigma=1d0
          
          Histo(26)%Info   = "M(Z)"
          Histo(26)%NBins  = 40
          Histo(26)%BinSize= 2.5d0*GeV
          Histo(26)%LowVal = 50d0*GeV
          Histo(26)%SetScale= 100d0
          Histo(26)%BinSmearing=.false.

          Histo(27)%Info   = "<Ehat>"
          Histo(27)%NBins  = 1
          Histo(27)%BinSize= 10000d0*GeV
          Histo(27)%LowVal = 0d0
          Histo(27)%SetScale=0.01d0

          Histo(28)%Info   = "CosThetaStar"
          Histo(28)%NBins  = 50
          Histo(28)%BinSize= 2d0/50d0 
          Histo(28)%LowVal = -1d0
          Histo(28)%SetScale= 1d0


!       ELSEIF( ObsSet.EQ.54 .or. ObsSet.EQ.58 ) THEN
       ELSEIF( ObsSet.EQ.54 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay ) observed at CMS at sqrt(s)=7 TeV
          if(abs(TopDecays).ne.4)  call Error("TopDecays needs to be 4")
          if(abs(ZDecays).ne.1)    call Error("ZDecays needs to be 1")
!          print *, Collider
!          if (Collider .ne. 11) call Error("Collider needs to be 11 for sqrt{s}=7 TeV")
          NumHistograms = 25
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif
          Histo(1)%Info   = "pT(lep+)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT(mu-)"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal =  0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT(mu+)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(b1)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "pT(b2)"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal =  0d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "pT(j1)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 20d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "pT(j2)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 20d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "pT(miss)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 20d0*GeV
          Histo(8)%LowVal =  0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "y(lep+)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 0.5d0
          Histo(9)%LowVal =-5.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y(mu-)"
          Histo(10)%NBins  = 20
          Histo(10)%BinSize= 0.5d0
          Histo(10)%LowVal =-5.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "y(mu+)"
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.5d0
          Histo(11)%LowVal =-5.0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "pT(Z)"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 20d0*GeV
          Histo(12)%LowVal =  0d0*GeV
          Histo(12)%SetScale= 100d0

          Histo(13)%Info   = "y(Z)"
          Histo(13)%NBins  = 20
          Histo(13)%BinSize= 0.5d0
          Histo(13)%LowVal =-5.0d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "pT(top)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 20d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "y(top)"
          Histo(15)%NBins  = 20
          Histo(15)%BinSize= 0.5d0
          Histo(15)%LowVal =-5.0d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "y(antitop)"
          Histo(16)%NBins  = 20
          Histo(16)%BinSize= 0.5d0
          Histo(16)%LowVal =-5.0d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "phi(mu-,mu+)"
          Histo(17)%NBins  = 26
          Histo(17)%BinSize= 0.25d0/2d0
          Histo(17)%LowVal = 0d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "HTTOT"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 20d0*GeV
          Histo(18)%LowVal = 0d0*GeV
          Histo(18)%SetScale= 100d0



ELSEIF( ObsSet.EQ.57 ) THEN! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay )
          if(abs(TopDecays).ne.4)  call Error("TopDecays needs to be 4")
          if(abs(ZDecays).ne.1)    call Error("ZDecays needs to be 1")
          NumHistograms = 4
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "y(top)"
          Histo(1)%NBins  = 30
          Histo(1)%BinSize= 0.2d0
          Histo(1)%LowVal =  -3d0
          Histo(1)%SetScale= 1d0

          Histo(2)%Info   = "y(tbar)"
          Histo(2)%NBins  = 30
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal =  -3d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "y(top) FWD"
          Histo(3)%NBins  = 1
          Histo(3)%BinSize= 10d0
          Histo(3)%LowVal = -5d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y(top) BWD"
          Histo(4)%NBins  = 1
          Histo(4)%BinSize= 10d0
          Histo(4)%LowVal = -5d0
          Histo(4)%SetScale= 1d0




ELSEIF( ObsSet.EQ.60  ) THEN! set of observables for Zprime, stable tops
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.5d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.5d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.1d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "M_TTbar+jet"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 40d0*GeV
          Histo(9)%LowVal = 350d0*GeV
          Histo(9)%SetScale= 100d0

ELSEIF( ObsSet.EQ.61 ) THEN! set of observables for Zprime, top decaying to dileptons
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 13
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lep1"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal = 20d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_bjet1"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal = 30d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lep2"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal = 20d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT_bjet2"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 20d0*GeV
          Histo(4)%LowVal = 30d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 20
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 1025d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "M_eff"
          Histo(6)%NBins  = 20
          Histo(6)%BinSize= 50d0*GeV
          Histo(6)%LowVal = 1025d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "y_lep1"
          Histo(7)%NBins  = 30
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal =-3.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "y_bjet1"
          Histo(8)%NBins  = 30
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-3.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "y_lep2"
          Histo(9)%NBins  = 30
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-3.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y_bjet2"
          Histo(10)%NBins  = 30
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-3.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "deltaPhi_LL"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.08d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "dPhiMinus"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.126d0
          Histo(12)%LowVal = -DblPi
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "dPhiPlus"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 0.063d0
          Histo(13)%LowVal = 0
          Histo(13)%SetScale= 1d0


ELSEIF( ObsSet.EQ.62 ) THEN! set of observables for Zprime, fully hadronic top decay
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.2  ) call Error("TopDecays needs to be 2!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0


ELSEIF( ObsSet.EQ.64 .OR. ObsSet.EQ.65 ) THEN ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4 .and. TopDecays.ne.3 ) call Error("TopDecays needs to be 3 or 4!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0



ELSEIF( ObsSet.EQ.66 ) THEN ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)
          if(Collider.ne.12)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 4!")
          NumHistograms = 8
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.2d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.2d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 40d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.08d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.04d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.04d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0



ELSEIF( ObsSet.EQ.67  ) THEN! set of observables for SM Z, stable tops
!           if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.5d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top (FW BWD)"
          Histo(4)%NBins  = 2
          Histo(4)%BinSize= 10d0
          Histo(4)%LowVal =-10.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 20d0*GeV
          Histo(5)%LowVal = 50d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.1d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "M_TTbar+jet"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 40d0*GeV
          Histo(9)%LowVal = 350d0*GeV
          Histo(9)%SetScale= 100d0



ELSEIF( ObsSet.EQ.70  ) THEN! set of observables for ee->ttb, stable tops
          if( .not. (Collider.ge.5 .and. Collider.le.7) )  call Error("Collider needs to be ee!")
          if(TopDecays.ne.0  ) call Error("TopDecays needs to be 0!")
          NumHistograms = 9
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif


          Histo(1)%Info   = "pT_ATop"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 50d0*GeV
          Histo(1)%LowVal = 0d0
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_Top"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 50d0*GeV
          Histo(2)%LowVal = 0d0
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "y_ATop"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 0.5d0
          Histo(3)%LowVal =-5.0d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "y_Top"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 0.5d0
          Histo(4)%LowVal =-5.0d0
          Histo(4)%SetScale= 1d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 50d0*GeV
          Histo(5)%LowVal = 350d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "deltaPhi_TTbar"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 0.1d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0

          Histo(7)%Info   = "CosTheta_scatter"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 0.1d0
          Histo(7)%LowVal = -1d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "CosTheta_star"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "M_TTbar+jet"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 40d0*GeV
          Histo(9)%LowVal = 350d0*GeV
          Histo(9)%SetScale= 100d0

          
ELSEIF( ObsSet.EQ.71 ) THEN! set of observables for ee->ttb, top decaying to dileptons
          if( .not. (Collider.ge.5 .and. Collider.le.7) )  call Error("Collider needs to be ee!")
          if(TopDecays.ne.1  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 19
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lep1"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_bjet1"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lep2"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT_bjet2"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal = 00d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 300d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "M_eff"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 300d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "y_lep1"
          Histo(7)%NBins  = 30
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal =-3.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "y_bjet1"
          Histo(8)%NBins  = 30
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-3.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "y_lep2"
          Histo(9)%NBins  = 30
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-3.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y_bjet2"
          Histo(10)%NBins  = 30
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-3.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "deltaPhi_LL"
          Histo(11)%NBins  = 50
          Histo(11)%BinSize= 0.08d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "dPhiMinus"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.126d0
          Histo(12)%LowVal = -DblPi
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "dPhiPlus"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 0.126d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0          

          Histo(14)%Info   = "CosTheta(e-,t)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 0.1d0
          Histo(14)%LowVal = -1d0
          Histo(14)%SetScale= 1d0
          
          

          
          
          Histo(15)%Info   = "dPhiMinus"
          Histo(15)%NBins  = 7
          Histo(15)%BinSize= 0.9d0
          Histo(15)%LowVal = -DblPi
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "dPhiPlus"
          Histo(16)%NBins  = 7
          Histo(16)%BinSize= 0.45d0
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 1d0          
          
          Histo(17)%Info   = "dPhiMinus vs. dPhiPlus"
          Histo(17)%NBins  = 49
          Histo(17)%BinSize= 0.9d0
          Histo(17)%LowVal = -DblPi
          Histo(17)%SetScale= 1d0
          

          
ELSEIF( ObsSet.EQ.72 ) THEN! set of observables for ee->ttb, tops decaying semi-leptonically
          if( .not. (Collider.ge.5 .and. Collider.le.7) )  call Error("Collider needs to be ee!")
          if(TopDecays.ne.4  ) call Error("TopDecays needs to be 1!")
          NumHistograms = 19
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT_lep1"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "pT_bjet1"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 10d0*GeV
          Histo(2)%LowVal = 0d0*GeV
          Histo(2)%SetScale= 100d0

          Histo(3)%Info   = "pT_lep2"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT_bjet2"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal = 00d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "M_TTbar"
          Histo(5)%NBins  = 50
          Histo(5)%BinSize= 10d0*GeV
          Histo(5)%LowVal = 300d0*GeV
          Histo(5)%SetScale= 100d0

          Histo(6)%Info   = "M_eff"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 300d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "y_lep1"
          Histo(7)%NBins  = 30
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal =-3.0d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "y_bjet1"
          Histo(8)%NBins  = 30
          Histo(8)%BinSize= 0.2d0
          Histo(8)%LowVal =-3.0d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "y_lep2"
          Histo(9)%NBins  = 30
          Histo(9)%BinSize= 0.2d0
          Histo(9)%LowVal =-3.0d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y_bjet2"
          Histo(10)%NBins  = 30
          Histo(10)%BinSize= 0.2d0
          Histo(10)%LowVal =-3.0d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "deltaPhi_LL"
          Histo(11)%NBins  = 40
          Histo(11)%BinSize= 0.08d0
          Histo(11)%LowVal = 0d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "dPhiMinus"
          Histo(12)%NBins  = 50
          Histo(12)%BinSize= 0.126d0
          Histo(12)%LowVal = -DblPi
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "dPhiPlus"
          Histo(13)%NBins  = 25
          Histo(13)%BinSize= 0.126d0
          Histo(13)%LowVal = 0d0
          Histo(13)%SetScale= 1d0          

          Histo(14)%Info   = "CosTheta(e-,t)"
          Histo(14)%NBins  = 25
          Histo(14)%BinSize= 0.1d0
          Histo(14)%LowVal = -1d0
          Histo(14)%SetScale= 1d0
          
         
          Histo(15)%Info   = "dPhiMinus"
          Histo(15)%NBins  = 7
          Histo(15)%BinSize= 0.9d0
          Histo(15)%LowVal = -DblPi
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "dPhiPlus"
          Histo(16)%NBins  = 7
          Histo(16)%BinSize= 0.45d0
          Histo(16)%LowVal = 0d0
          Histo(16)%SetScale= 1d0          
          
          Histo(17)%Info   = "dPhiMinus vs. dPhiPlus"
          Histo(17)%NBins  = 49
          Histo(17)%BinSize= 0.9d0
          Histo(17)%LowVal = -DblPi
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "y_top"
          Histo(18)%NBins  = 40
          Histo(18)%BinSize= 0.25d0
          Histo(18)%LowVal = -5d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "y_top (FB)"
          Histo(19)%NBins  = 2
          Histo(19)%BinSize= 10d0
          Histo(19)%LowVal = -10d0
          Histo(19)%SetScale= 1d0
                                       




ELSEIF( ObsSet.EQ.81 ) THEN! set of observables for ttb+H (stable tops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.0)  call Error("TopDecays needs to be 0")
          NumHistograms = 2
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(top)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0
          
          Histo(2)%Info   = "pT(Higgs)"
          Histo(2)%NBins  = 50
          Histo(2)%BinSize= 20d0*GeV
          Histo(2)%LowVal =  0d0*GeV
          Histo(2)%SetScale= 100d0


          
          
ELSEIF( ObsSet.EQ.82 ) THEN! set of observables for ttb+H (di-leptonic tops)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
          if(TopDecays.ne.1)  call Error("TopDecays needs to be 1")
          Num2DHistograms=11
          NumHistograms = 7
          NumHistograms = NumHistograms + Num2DHistograms

          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          if( .not.allocated(Histo2D) ) then
                allocate( Histo2D(1:Num2DHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

         Histo(1)%Info   = "pT(top)"
         Histo(1)%NBins  = 20
         Histo(1)%BinSize= 50d0*GeV
         Histo(1)%LowVal =  0d0*GeV
         Histo(1)%SetScale= 100d0

         Histo(2)%Info   = "pT(Atop)"
         Histo(2)%NBins  = 20
         Histo(2)%BinSize= 50d0*GeV
         Histo(2)%LowVal =  0d0*GeV
         Histo(2)%SetScale= 100d0

         Histo(3)%Info   = "y(top)"
         Histo(3)%NBins  = 20
         Histo(3)%BinSize= 0.5d0
         Histo(3)%LowVal = -5d0
         Histo(3)%SetScale= 1d0

         Histo(4)%Info   = "y(ATop)"
         Histo(4)%NBins  = 20
         Histo(4)%BinSize= 0.5d0
         Histo(4)%LowVal =-5.0d0
         Histo(4)%SetScale= 1d0

         Histo(5)%Info   = "M_TTbar"
         Histo(5)%NBins  = 20
         Histo(5)%BinSize= 40d0*GeV
         Histo(5)%LowVal = 350d0*GeV
         Histo(5)%SetScale= 100d0

         Histo(6)%Info   = "pT(Higgs)"
         Histo(6)%NBins  = 20
         Histo(6)%BinSize= 50d0*GeV
         Histo(6)%LowVal =  0d0*GeV
         Histo(6)%SetScale= 100d0

         Histo(7)%Info   = "y(Higgs)"
         Histo(7)%NBins  = 20
         Histo(7)%BinSize= 0.5d0
         Histo(7)%LowVal =-5.0d0
         Histo(7)%SetScale= 1d0

!! 2-d histogram
!! NB -- 2-dim histograms MUST have "2D_" as the first three characters of the info
!! -- this is to ensure that they are properly written out
         i=1
         j=2
         Histo2D(1)%Info  = "2D_pT(top),pT(ATop)"
         Histo2D(1)%HistoDim=2
         Histo2D(1)%NBins2=Histo(i)%NBins
         Histo2D(1)%NBins3=Histo(j)%NBins
         Histo2D(1)%NBins=(Histo2D(1)%NBins2)*(Histo2D(1)%NBins3)
         Histo2D(1)%BinSize= 1d0
         Histo2D(1)%LowVal =  0d0*GeV
         Histo2D(1)%SetScale =  1d0
         Histo2D(1)%BinSize2= Histo(i)%BinSize
         Histo2D(1)%LowVal2 =  Histo(i)%LowVal
         Histo2D(1)%SetScale2= Histo(i)%SetScale
         Histo2D(1)%BinSize3= Histo(j)%BinSize
         Histo2D(1)%LowVal3 =  Histo(j)%LowVal
         Histo2D(1)%SetScale3= Histo(j)%SetScale

         i=1
         j=3
         Histo2D(2)%Info  = "2D_pT(top),y(top)"
         Histo2D(2)%HistoDim=2
         Histo2D(2)%NBins2=Histo(i)%NBins
         Histo2D(2)%NBins3=Histo(j)%NBins
         Histo2D(2)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(2)%BinSize= 1d0
         Histo2D(2)%LowVal =  0d0*GeV
         Histo2D(2)%SetScale =  1d0
         Histo2D(2)%BinSize2= Histo(i)%BinSize
         Histo2D(2)%LowVal2 =  Histo(i)%LowVal
         Histo2D(2)%SetScale2= Histo(i)%SetScale
         Histo2D(2)%BinSize3= Histo(j)%BinSize
         Histo2D(2)%LowVal3 =  Histo(j)%LowVal
         Histo2D(2)%SetScale3= Histo(j)%SetScale

         i=1
         j=4
         Histo2D(3)%Info  = "2D_pT(top),y(ATop)"
         Histo2D(3)%HistoDim=2
         Histo2D(3)%NBins2=Histo(i)%NBins
         Histo2D(3)%NBins3=Histo(j)%NBins
         Histo2D(3)%NBins=(Histo2D(1)%NBins2)*(Histo2D(1)%NBins3)
         Histo2D(3)%BinSize= 1d0
         Histo2D(3)%LowVal =  0d0*GeV
         Histo2D(3)%SetScale =  1d0
         Histo2D(3)%BinSize2= Histo(i)%BinSize
         Histo2D(3)%LowVal2 =  Histo(i)%LowVal
         Histo2D(3)%SetScale2= Histo(i)%SetScale
         Histo2D(3)%BinSize3= Histo(j)%BinSize
         Histo2D(3)%LowVal3 =  Histo(j)%LowVal
         Histo2D(3)%SetScale3= Histo(j)%SetScale
         
         i=1
         j=6
         Histo2D(4)%Info  = "2D_pT(top),pT(H)"
         Histo2D(4)%HistoDim=2
         Histo2D(4)%NBins2=Histo(i)%NBins
         Histo2D(4)%NBins3=Histo(j)%NBins
         Histo2D(4)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(4)%BinSize= 1d0
         Histo2D(4)%LowVal =  0d0*GeV
         Histo2D(4)%SetScale =  1d0
         Histo2D(4)%BinSize2= Histo(i)%BinSize
         Histo2D(4)%LowVal2 =  Histo(i)%LowVal
         Histo2D(4)%SetScale2= Histo(i)%SetScale
         Histo2D(4)%BinSize3= Histo(j)%BinSize
         Histo2D(4)%LowVal3 =  Histo(j)%LowVal
         Histo2D(4)%SetScale3= Histo(j)%SetScale

         i=1
         j=7
         Histo2D(5)%Info  = "2D_pT(top),y(H)"
         Histo2D(5)%HistoDim=2
         Histo2D(5)%NBins2=Histo(i)%NBins
         Histo2D(5)%NBins3=Histo(j)%NBins
         Histo2D(5)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(5)%BinSize= 1d0
         Histo2D(5)%LowVal =  0d0*GeV
         Histo2D(5)%SetScale =  1d0
         Histo2D(5)%BinSize2= Histo(i)%BinSize
         Histo2D(5)%LowVal2 =  Histo(i)%LowVal
         Histo2D(5)%SetScale2= Histo(i)%SetScale
         Histo2D(5)%BinSize3= Histo(j)%BinSize
         Histo2D(5)%LowVal3 =  Histo(j)%LowVal
         Histo2D(5)%SetScale3= Histo(j)%SetScale

         i=3
         j=5
         Histo2D(6)%Info  = "2D_y(top),M(tt)"
         Histo2D(6)%HistoDim=2
         Histo2D(6)%NBins2=Histo(i)%NBins
         Histo2D(6)%NBins3=Histo(j)%NBins
         Histo2D(6)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(6)%BinSize= 1d0
         Histo2D(6)%LowVal =  0d0*GeV
         Histo2D(6)%SetScale =  1d0
         Histo2D(6)%BinSize2= Histo(i)%BinSize
         Histo2D(6)%LowVal2 =  Histo(i)%LowVal
         Histo2D(6)%SetScale2= Histo(i)%SetScale
         Histo2D(6)%BinSize3= Histo(j)%BinSize
         Histo2D(6)%LowVal3 =  Histo(j)%LowVal
         Histo2D(6)%SetScale3= Histo(j)%SetScale

         i=3
         j=6
         Histo2D(7)%Info  = "2D_y(top),pT(H)"
         Histo2D(7)%HistoDim=2
         Histo2D(7)%NBins2=Histo(i)%NBins
         Histo2D(7)%NBins3=Histo(j)%NBins
         Histo2D(7)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(7)%BinSize= 1d0
         Histo2D(7)%LowVal =  0d0*GeV
         Histo2D(7)%SetScale =  1d0
         Histo2D(7)%BinSize2= Histo(i)%BinSize
         Histo2D(7)%LowVal2 =  Histo(i)%LowVal
         Histo2D(7)%SetScale2= Histo(i)%SetScale
         Histo2D(7)%BinSize3= Histo(j)%BinSize
         Histo2D(7)%LowVal3 =  Histo(j)%LowVal
         Histo2D(7)%SetScale3= Histo(j)%SetScale

         i=3
         j=7
         Histo2D(8)%Info  = "2D_y(top),y(H)"
         Histo2D(8)%HistoDim=2
         Histo2D(8)%NBins2=Histo(i)%NBins
         Histo2D(8)%NBins3=Histo(j)%NBins
         Histo2D(8)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(8)%BinSize= 1d0
         Histo2D(8)%LowVal =  0d0*GeV
         Histo2D(8)%SetScale =  1d0
         Histo2D(8)%BinSize2= Histo(i)%BinSize
         Histo2D(8)%LowVal2 =  Histo(i)%LowVal
         Histo2D(8)%SetScale2= Histo(i)%SetScale
         Histo2D(8)%BinSize3= Histo(j)%BinSize
         Histo2D(8)%LowVal3 =  Histo(j)%LowVal
         Histo2D(8)%SetScale3= Histo(j)%SetScale

         i=5
         j=6
         Histo2D(9)%Info  = "2D_M(tt),pt(H)"
         Histo2D(9)%HistoDim=2
         Histo2D(9)%NBins2=Histo(i)%NBins
         Histo2D(9)%NBins3=Histo(j)%NBins
         Histo2D(9)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(9)%BinSize= 1d0
         Histo2D(9)%LowVal =  0d0*GeV
         Histo2D(9)%SetScale =  1d0
         Histo2D(9)%BinSize2= Histo(i)%BinSize
         Histo2D(9)%LowVal2 =  Histo(i)%LowVal
         Histo2D(9)%SetScale2= Histo(i)%SetScale
         Histo2D(9)%BinSize3= Histo(j)%BinSize
         Histo2D(9)%LowVal3 =  Histo(j)%LowVal
         Histo2D(9)%SetScale3= Histo(j)%SetScale

         i=5
         j=7
         Histo2D(10)%Info  = "2D_M(tt),y(H)"
         Histo2D(10)%HistoDim=2
         Histo2D(10)%NBins2=Histo(i)%NBins
         Histo2D(10)%NBins3=Histo(j)%NBins
         Histo2D(10)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(10)%BinSize= 1d0
         Histo2D(10)%LowVal =  0d0*GeV
         Histo2D(10)%SetScale =  1d0
         Histo2D(10)%BinSize2= Histo(i)%BinSize
         Histo2D(10)%LowVal2 =  Histo(i)%LowVal
         Histo2D(10)%SetScale2= Histo(i)%SetScale
         Histo2D(10)%BinSize3= Histo(j)%BinSize
         Histo2D(10)%LowVal3 =  Histo(j)%LowVal
         Histo2D(10)%SetScale3= Histo(j)%SetScale

         i=6
         j=7
         Histo2D(11)%Info  = "2D_pt(H),y(H)"
         Histo2D(11)%HistoDim=2
         Histo2D(11)%NBins2=Histo(i)%NBins
         Histo2D(11)%NBins3=Histo(j)%NBins
         Histo2D(11)%NBins=(Histo2D(2)%NBins2)*(Histo2D(2)%NBins3)
         Histo2D(11)%BinSize= 1d0
         Histo2D(11)%LowVal =  0d0*GeV
         Histo2D(11)%SetScale =  1d0
         Histo2D(11)%BinSize2= Histo(i)%BinSize
         Histo2D(11)%LowVal2 =  Histo(i)%LowVal
         Histo2D(11)%SetScale2= Histo(i)%SetScale
         Histo2D(11)%BinSize3= Histo(j)%BinSize
         Histo2D(11)%LowVal3 =  Histo(j)%LowVal
         Histo2D(11)%SetScale3= Histo(j)%SetScale




         do NHisto=NumHistograms-Num2DHistograms+1,NumHistograms
            Histo(NHisto)%Info    = Histo2D(NHisto-NumHistograms+Num2DHistograms)%Info
            Histo(NHisto)%NBins   = Histo2D(NHisto-NumHistograms+Num2DHistograms)%NBins
            Histo(NHisto)%BinSize = Histo2D(NHisto-NumHistograms+Num2DHistograms)%BinSize
            Histo(NHisto)%Lowval  = Histo2D(NHisto-NumHistograms+Num2DHistograms)%Lowval
            Histo(NHisto)%Setscale= Histo2D(NHisto-NumHistograms+Num2DHistograms)%Setscale
         enddo

            
            

                                       




ELSEIF( ObsSet.EQ.83 ) THEN! set of observables for ttb+H (semi-leptonic tops, stable Higgs)
          if(Collider.ne.1)  call Error("Collider needs to be LHC!")
!          if(TopDecays.ne.1 .and. TopDecays.ne.3 .and. TopDecays.ne.4)  call Error("TopDecays needs to be 1 or 3 or 4")
          if(TopDecays.ne.3 .and. TopDecays.ne.4 )  call Error("TopDecays needs to be 3 or 4")
          NumHistograms = 14
          if( .not.allocated(Histo) ) then
                allocate( Histo(1:NumHistograms), stat=AllocStatus  )
                if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
          endif

          Histo(1)%Info   = "pT(top)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 20d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta(top)"
          Histo(2)%NBins  = 20
          Histo(2)%BinSize= 0.5d0
          Histo(2)%LowVal = -5d0
          Histo(2)%SetScale= 1d0
          
          Histo(3)%Info   = "pT(Higgs)"
          Histo(3)%NBins  = 50
          Histo(3)%BinSize= 20d0*GeV
          Histo(3)%LowVal =  0d0*GeV
          Histo(3)%SetScale= 100d0
          
          Histo(4)%Info   = "eta(Higgs)"
          Histo(4)%NBins  = 20
          Histo(4)%BinSize= 0.5d0
          Histo(4)%LowVal = -5d0
          Histo(4)%SetScale= 1d0
          
          Histo(5)%Info   = "delta_eta(t,tbar)"
          Histo(5)%NBins  = 20
          Histo(5)%BinSize= 0.25d0
          Histo(5)%LowVal = 0d0
          Histo(5)%SetScale= 1d0
                   
          Histo(6)%Info   = "delta_eta(l+,l-)"
          Histo(6)%NBins  = 20
          Histo(6)%BinSize= 0.25d0
          Histo(6)%LowVal = 0d0
          Histo(6)%SetScale= 1d0
                   
          Histo(7)%Info   = "delta_eta(b,bbar)"
          Histo(7)%NBins  = 20
          Histo(7)%BinSize= 0.25d0
          Histo(7)%LowVal = 0d0
          Histo(7)%SetScale= 1d0                   
         
          Histo(8)%Info   = "cos(theta_ll)" 
          Histo(8)%NBins  = 20
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -1d0
          Histo(8)%SetScale= 1d0          
         
          Histo(9)%Info   = "cos(theta_bb)" 
          Histo(9)%NBins  = 20
          Histo(9)%BinSize= 0.1d0
          Histo(9)%LowVal = -1d0
          Histo(9)%SetScale= 1d0
          
          Histo(10)%Info   = "cos(theta_H)" 
          Histo(10)%NBins  = 20
          Histo(10)%BinSize= 0.1d0
          Histo(10)%LowVal = -1d0
          Histo(10)%SetScale= 1d0
         
          Histo(11)%Info   = "cos(theta_t)" 
          Histo(11)%NBins  = 20
          Histo(11)%BinSize= 0.1d0
          Histo(11)%LowVal = -1d0
          Histo(11)%SetScale= 1d0
          
          Histo(12)%Info   = "D_cp" 
          Histo(12)%NBins  = 20
          Histo(12)%BinSize= 0.05d0
          Histo(12)%LowVal = 0d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "pT(jet)"
          Histo(13)%NBins  = 50
          Histo(13)%BinSize= 10d0*GeV
          Histo(13)%LowVal =  0d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "N(jets)"
          Histo(14)%NBins  = 6
          Histo(14)%BinSize= 1d0
          Histo(14)%LowVal =  0d0
          Histo(14)%SetScale= 1d0     



 ELSEIF( ObsSet.EQ.91 ) THEN! set of observables for tH (stable tops)                                                                                                         
         if(Collider.ne.1)  call Error("Collider needs to be LHC!")
         if(TopDecays.ne.0)  call Error("TopDecays needs to be 0")
         NumHistograms = 22
         if( .not.allocated(Histo) ) then
            allocate( Histo(1:NumHistograms), stat=AllocStatus  )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
         endif

          Histo(1)%Info   = "pT(top)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta(top)"
          Histo(2)%NBins  = 60
          Histo(2)%BinSize= 0.1d0
          Histo(2)%LowVal = -3d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "y(top)"
          Histo(3)%NBins  = 60
          Histo(3)%BinSize= 0.1d0
          Histo(3)%LowVal = -3d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "pT(jet)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "y(jet)"
          Histo(5)%NBins  = 60
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -6d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "pT(Higgs)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "m(tj)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "y(Higgs)"
          Histo(8)%NBins  = 100
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -5d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "eta(Higgs)"
          Histo(9)%NBins  = 100
          Histo(9)%BinSize= 0.1d0
          Histo(9)%LowVal = -5d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y(tj)"
          Histo(10)%NBins  = 60
          Histo(10)%BinSize= 0.1d0
          Histo(10)%LowVal = -3d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "eta(tj)"
          Histo(11)%NBins  = 60
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = -6d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "Deltay(t,j)"
          Histo(12)%NBins  = 100
          Histo(12)%BinSize= 0.1d0
          Histo(12)%LowVal = -5d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "m(Htop)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 50d0*GeV
          Histo(13)%LowVal = 0d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "m(Hj)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 10d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "y(Ht)"
          Histo(15)%NBins  = 60
          Histo(15)%BinSize= 0.1d0
          Histo(15)%LowVal = -3d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "eta(Ht)"
          Histo(16)%NBins  = 60
          Histo(16)%BinSize= 0.2d0
          Histo(16)%LowVal = -6d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "y(Hj)"
          Histo(17)%NBins  = 80
          Histo(17)%BinSize= 0.1d0
          Histo(17)%LowVal = -4d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "eta(Hj)"
          Histo(18)%NBins  = 100
          Histo(18)%BinSize= 0.1d0
          Histo(18)%LowVal = -5d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "Deltay(H,t)"
          Histo(19)%NBins  = 100
          Histo(19)%BinSize= 0.1d0
          Histo(19)%LowVal = -5d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "Delta eta(H,T)"
          Histo(20)%NBins  = 100
          Histo(20)%BinSize= 0.1d0
          Histo(20)%LowVal = -5d0
          Histo(20)%SetScale= 1d0

          Histo(21)%Info   = "Deltay(H,j)"
          Histo(21)%NBins  = 100
          Histo(21)%BinSize= 0.1d0
          Histo(21)%LowVal = -5d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "Delta eta(H,j)"
          Histo(22)%NBins  = 100
          Histo(22)%BinSize= 0.1d0
          Histo(22)%LowVal = -5d0
          Histo(22)%SetScale= 1d0


ELSEIF( ObsSet.EQ.92 ) THEN! set of observables for tb+H (stable tops)                                                                                             
                                                                                                                                                                         
         if(Collider.ne.1)  call Error("Collider needs to be LHC!")
         if(TopDecays.ne.0)  call Error("TopDecays needs to be 0")
         NumHistograms = 22
         if( .not.allocated(Histo) ) then
            allocate( Histo(1:NumHistograms), stat=AllocStatus  )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
         endif

          Histo(1)%Info   = "pT(Atop)"
          Histo(1)%NBins  = 50
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal =  0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta(Atop)"
          Histo(2)%NBins  = 60
          Histo(2)%BinSize= 0.1d0
          Histo(2)%LowVal = -3d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "y(Atop)"
          Histo(3)%NBins  = 60
          Histo(3)%BinSize= 0.1d0
          Histo(3)%LowVal = -3d0
          Histo(3)%SetScale= 1d0

          Histo(4)%Info   = "pT(jet)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal =  0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "y(jet)"
          Histo(5)%NBins  = 60
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -6d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "pT(Higgs)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal =  0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "m(tbarj)"
          Histo(7)%NBins  = 50
          Histo(7)%BinSize= 10d0*GeV
          Histo(7)%LowVal =  0d0*GeV
          Histo(7)%SetScale= 100d0

          Histo(8)%Info   = "y(Higgs)"
          Histo(8)%NBins  = 100
          Histo(8)%BinSize= 0.1d0
          Histo(8)%LowVal = -5d0
          Histo(8)%SetScale= 1d0

          Histo(9)%Info   = "eta(Higgs)"
          Histo(9)%NBins  = 100
          Histo(9)%BinSize= 0.1d0
          Histo(9)%LowVal = -5d0
          Histo(9)%SetScale= 1d0

          Histo(10)%Info   = "y(tbarj)"
          Histo(10)%NBins  = 60
          Histo(10)%BinSize= 0.1d0
          Histo(10)%LowVal = -3d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "eta(tbarj)"
          Histo(11)%NBins  = 60
          Histo(11)%BinSize= 0.2d0
          Histo(11)%LowVal = -6d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "Deltay(tbar,j)"
          Histo(12)%NBins  = 100
          Histo(12)%BinSize= 0.1d0
          Histo(12)%LowVal = -5d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "m(Htop)"
          Histo(13)%NBins  = 40
          Histo(13)%BinSize= 50d0*GeV
          Histo(13)%LowVal = 0d0*GeV
          Histo(13)%SetScale= 100d0

          Histo(14)%Info   = "m(Hj)"
          Histo(14)%NBins  = 50
          Histo(14)%BinSize= 10d0*GeV
          Histo(14)%LowVal =  0d0*GeV
          Histo(14)%SetScale= 100d0

          Histo(15)%Info   = "y(Htbar)"
          Histo(15)%NBins  = 60
          Histo(15)%BinSize= 0.1d0
          Histo(15)%LowVal = -3d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "eta(Htbar)"
          Histo(16)%NBins  = 60
          Histo(16)%BinSize= 0.2d0
          Histo(16)%LowVal = -6d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "y(Hj)"
          Histo(17)%NBins  = 80
          Histo(17)%BinSize= 0.1d0
          Histo(17)%LowVal = -4d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "eta(Hj)"
          Histo(18)%NBins  = 100
          Histo(18)%BinSize= 0.1d0
          Histo(18)%LowVal = -5d0
          Histo(18)%SetScale= 1d0

          Histo(19)%Info   = "Deltay(H,tbar)"
          Histo(19)%NBins  = 100
          Histo(19)%BinSize= 0.1d0
          Histo(19)%LowVal = -5d0
          Histo(19)%SetScale= 1d0

          Histo(20)%Info   = "Delta eta(H,Tbar)"
          Histo(20)%NBins  = 100
          Histo(20)%BinSize= 0.1d0
          Histo(20)%LowVal = -5d0
          Histo(20)%SetScale= 1d0


          Histo(21)%Info   = "Deltay(H,j)"
          Histo(21)%NBins  = 100
          Histo(21)%BinSize= 0.1d0
          Histo(21)%LowVal = -5d0
          Histo(21)%SetScale= 1d0

          Histo(22)%Info   = "Delta eta(H,j)"
          Histo(22)%NBins  = 100
          Histo(22)%BinSize= 0.1d0
          Histo(22)%LowVal = -5d0
          Histo(22)%SetScale= 1d0
       
ELSEIF( ObsSet.EQ.93 .or. ObsSet .eq. 95 ) THEN! set of observables for tH (leptonic top decay)                                                            
         if(Collider.ne.1)  call Error("Collider needs to be LHC!")
         if(TopDecays.ne.1)  call Error("TopDecays needs to be 1")
         NumHistograms = 37
         if( .not.allocated(Histo) ) then
            allocate( Histo(1:NumHistograms), stat=AllocStatus  )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
         endif


          Histo(1)%Info   = "pT(e+)"
          Histo(1)%NBins  = 20
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta(e+)"
          Histo(2)%NBins  = 22
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal = -2.2d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT(miss)"
          Histo(3)%NBins  = 20
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(j)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "eta(j)"
          Histo(5)%NBins  = 60
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -6d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "pT(b)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "eta(b)"
          Histo(7)%NBins  = 60
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal = -6d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "pT(H)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m(W,trans)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal = 0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "y(H)"
          Histo(10)%NBins  = 100
          Histo(10)%BinSize= 0.1d0
          Histo(10)%LowVal = -5d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "y(b,ep)"
          Histo(11)%NBins  = 100
          Histo(11)%BinSize= 0.1d0
          Histo(11)%LowVal = -5d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "y(j,ep)"
          Histo(12)%NBins  = 100
          Histo(12)%BinSize= 0.1d0
          Histo(12)%LowVal = -5d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "eta(H)"
          Histo(13)%NBins  = 100
          Histo(13)%BinSize= 0.1d0
          Histo(13)%LowVal = -5d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "eta(W)"
          Histo(14)%NBins  = 100
          Histo(14)%BinSize= 0.1d0
          Histo(14)%LowVal = -5d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "Delta eta(j,b)"
          Histo(15)%NBins  = 100
          Histo(15)%BinSize= 0.1d0
          Histo(15)%LowVal = -5d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "Delta R(b,e+)"
          Histo(16)%NBins  = 100
          Histo(16)%BinSize= 0.1d0
          Histo(16)%LowVal = -5d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "Delta R(j,e+)"
          Histo(17)%NBins  = 100
          Histo(17)%BinSize= 0.1d0
          Histo(17)%LowVal = -5d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "m(H,b)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 0d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "m(W,b)"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 10d0*GeV
          Histo(19)%LowVal = 0d0*GeV
          Histo(19)%SetScale= 100d0

          Histo(20)%Info   = "m(W,j)"
          Histo(20)%NBins  = 50
          Histo(20)%BinSize= 10d0*GeV
          Histo(20)%LowVal = 0d0*GeV
          Histo(20)%SetScale= 100d0

          Histo(21)%Info   = "m(H,j)"
          Histo(21)%NBins  = 30
          Histo(21)%BinSize= 5d0*GeV
          Histo(21)%LowVal = 125d0*GeV
          Histo(21)%SetScale= 100d0

          Histo(22)%Info   = "y(H,j)"
          Histo(22)%NBins  = 80
          Histo(22)%BinSize= 0.1d0
          Histo(22)%LowVal = -4d0
          Histo(22)%SetScale= 1d0

          Histo(23)%Info   = "eta(H,j)"
          Histo(23)%NBins  = 80
          Histo(23)%BinSize= 0.1d0
          Histo(23)%LowVal = -4d0
          Histo(23)%SetScale= 1d0

          Histo(24)%Info   = "eta(H,b)"
          Histo(24)%NBins  = 80
          Histo(24)%BinSize= 0.1d0
          Histo(24)%LowVal = -4d0
          Histo(24)%SetScale= 1d0

          Histo(25)%Info   = "eta(W,b)"
          Histo(25)%NBins  = 80
          Histo(25)%BinSize= 0.1d0
          Histo(25)%LowVal = -4d0
          Histo(25)%SetScale= 1d0

          Histo(26)%Info   = "eta(W,j)"
          Histo(26)%NBins  = 80
          Histo(26)%BinSize= 0.1d0
          Histo(26)%LowVal = -4d0
          Histo(26)%SetScale= 1d0

          Histo(27)%Info   = "Delta eta(H,j)"
          Histo(27)%NBins  = 100
          Histo(27)%BinSize= 0.1d0
          Histo(27)%LowVal = -5d0
          Histo(27)%SetScale= 1d0

          Histo(28)%Info   = "Delta y(H,j)"
          Histo(28)%NBins  = 100
          Histo(28)%BinSize= 0.1d0
          Histo(28)%LowVal = -5d0
          Histo(28)%SetScale= 1d0

          Histo(29)%Info   = "Delta eta(H,b)"
          Histo(29)%NBins  = 100
          Histo(29)%BinSize= 0.1d0
          Histo(29)%LowVal = -5d0
          Histo(29)%SetScale= 1d0

          Histo(30)%Info   = "Delta eta(W,j)"
          Histo(30)%NBins  = 100
          Histo(30)%BinSize= 0.1d0
          Histo(30)%LowVal = -5d0
          Histo(30)%SetScale= 1d0

          Histo(31)%Info   = "Delta eta(W,b)"
          Histo(31)%NBins  = 100
          Histo(31)%BinSize= 0.1d0
          Histo(31)%LowVal = -5d0
          Histo(31)%SetScale= 1d0

          Histo(32)%Info   = "Delta eta(H,ep)"
          Histo(32)%NBins  = 100
          Histo(32)%BinSize= 0.1d0
          Histo(32)%LowVal = -5d0
          Histo(32)%SetScale= 1d0

          Histo(33)%Info   = "m(H,b,ep)"
          Histo(33)%NBins  = 50
          Histo(33)%BinSize= 10d0*GeV
          Histo(33)%LowVal = 0d0*GeV
          Histo(33)%SetScale= 100d0

          Histo(34)%Info   = "m(H,j,ep)"
          Histo(34)%NBins  = 50
          Histo(34)%BinSize= 10d0*GeV
          Histo(34)%LowVal = 0d0*GeV
          Histo(34)%SetScale= 100d0

          Histo(35)%Info   = "Delta eta(t,j)"
          Histo(35)%NBins  = 100
          Histo(35)%BinSize= 0.1d0
          Histo(35)%LowVal = -5d0
          Histo(35)%SetScale= 1d0

          Histo(36)%Info   = "Delta eta(t,H)"
          Histo(36)%NBins  = 100
          Histo(36)%BinSize= 0.1d0
          Histo(36)%LowVal = -5d0
          Histo(36)%SetScale= 1d0

          Histo(37)%Info   = "cos theta(l,j1)"
          Histo(37)%NBins  = 20
          Histo(37)%BinSize= 0.1d0
          Histo(37)%LowVal = -1d0
          Histo(37)%SetScale= 1d0


ELSEIF( ObsSet.EQ.94 .or. ObsSet .eq. 96) THEN! set of observables for tbH (leptonic top decay)                                                         

                                          
         if(Collider.ne.1)  call Error("Collider needs to be LHC!")
         if(TopDecays.ne.1)  call Error("TopDecays needs to be 1")
         NumHistograms = 37
         if( .not.allocated(Histo) ) then
            allocate( Histo(1:NumHistograms), stat=AllocStatus  )
            if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
         endif


          Histo(1)%Info   = "pT(-)"
          Histo(1)%NBins  = 20
          Histo(1)%BinSize= 10d0*GeV
          Histo(1)%LowVal = 0d0*GeV
          Histo(1)%SetScale= 100d0

          Histo(2)%Info   = "eta(e-)"
          Histo(2)%NBins  = 22
          Histo(2)%BinSize= 0.2d0
          Histo(2)%LowVal = -2.2d0
          Histo(2)%SetScale= 1d0

          Histo(3)%Info   = "pT(miss)"
          Histo(3)%NBins  = 20
          Histo(3)%BinSize= 10d0*GeV
          Histo(3)%LowVal = 0d0*GeV
          Histo(3)%SetScale= 100d0

          Histo(4)%Info   = "pT(j)"
          Histo(4)%NBins  = 50
          Histo(4)%BinSize= 10d0*GeV
          Histo(4)%LowVal = 0d0*GeV
          Histo(4)%SetScale= 100d0

          Histo(5)%Info   = "eta(j)"
          Histo(5)%NBins  = 60
          Histo(5)%BinSize= 0.2d0
          Histo(5)%LowVal = -6d0
          Histo(5)%SetScale= 1d0

          Histo(6)%Info   = "pT(bbar)"
          Histo(6)%NBins  = 50
          Histo(6)%BinSize= 10d0*GeV
          Histo(6)%LowVal = 0d0*GeV
          Histo(6)%SetScale= 100d0

          Histo(7)%Info   = "eta(bbar)"
          Histo(7)%NBins  = 60
          Histo(7)%BinSize= 0.2d0
          Histo(7)%LowVal = -6d0
          Histo(7)%SetScale= 1d0

          Histo(8)%Info   = "pT(H)"
          Histo(8)%NBins  = 50
          Histo(8)%BinSize= 10d0*GeV
          Histo(8)%LowVal = 0d0*GeV
          Histo(8)%SetScale= 100d0

          Histo(9)%Info   = "m(W,trans)"
          Histo(9)%NBins  = 50
          Histo(9)%BinSize= 10d0*GeV
          Histo(9)%LowVal = 0d0*GeV
          Histo(9)%SetScale= 100d0

          Histo(10)%Info   = "y(H)"
          Histo(10)%NBins  = 100
          Histo(10)%BinSize= 0.1d0
          Histo(10)%LowVal = -5d0
          Histo(10)%SetScale= 1d0

          Histo(11)%Info   = "y(b,em)"
          Histo(11)%NBins  = 100
          Histo(11)%BinSize= 0.1d0
          Histo(11)%LowVal = -5d0
          Histo(11)%SetScale= 1d0

          Histo(12)%Info   = "y(j,em)"
          Histo(12)%NBins  = 100
          Histo(12)%BinSize= 0.1d0
          Histo(12)%LowVal = -5d0
          Histo(12)%SetScale= 1d0

          Histo(13)%Info   = "eta(H)"
          Histo(13)%NBins  = 100
          Histo(13)%BinSize= 0.1d0
          Histo(13)%LowVal = -5d0
          Histo(13)%SetScale= 1d0

          Histo(14)%Info   = "eta(W)"
          Histo(14)%NBins  = 100
          Histo(14)%BinSize= 0.1d0
          Histo(14)%LowVal = -5d0
          Histo(14)%SetScale= 1d0

          Histo(15)%Info   = "Delta eta(j,bbar)"
          Histo(15)%NBins  = 100
          Histo(15)%BinSize= 0.1d0
          Histo(15)%LowVal = -5d0
          Histo(15)%SetScale= 1d0

          Histo(16)%Info   = "Delta R(bbar,e-)"
          Histo(16)%NBins  = 100
          Histo(16)%BinSize= 0.1d0
          Histo(16)%LowVal = -5d0
          Histo(16)%SetScale= 1d0

          Histo(17)%Info   = "Delta R(j,e-)"
          Histo(17)%NBins  = 100
          Histo(17)%BinSize= 0.1d0
          Histo(17)%LowVal = -5d0
          Histo(17)%SetScale= 1d0

          Histo(18)%Info   = "m(H,bbar)"
          Histo(18)%NBins  = 50
          Histo(18)%BinSize= 10d0*GeV
          Histo(18)%LowVal = 0d0*GeV
          Histo(18)%SetScale= 100d0

          Histo(19)%Info   = "m(W,bbar)"
          Histo(19)%NBins  = 50
          Histo(19)%BinSize= 10d0*GeV
          Histo(19)%LowVal = 0d0*GeV
          Histo(19)%SetScale= 100d0

          Histo(20)%Info   = "m(W,j)"
          Histo(20)%NBins  = 50
          Histo(20)%BinSize= 10d0*GeV
          Histo(20)%LowVal = 0d0*GeV
          Histo(20)%SetScale= 100d0

          Histo(21)%Info   = "m(H,j)"
          Histo(21)%NBins  = 30
          Histo(21)%BinSize= 5d0*GeV
          Histo(21)%LowVal = 125d0*GeV
          Histo(21)%SetScale= 100d0

          Histo(22)%Info   = "y(H,j)"
          Histo(22)%NBins  = 80
          Histo(22)%BinSize= 0.1d0
          Histo(22)%LowVal = -4d0
          Histo(22)%SetScale= 1d0

          Histo(23)%Info   = "eta(H,j)"
          Histo(23)%NBins  = 80
          Histo(23)%BinSize= 0.1d0
          Histo(23)%LowVal = -4d0
          Histo(23)%SetScale= 1d0

          Histo(24)%Info   = "eta(H,bbar)"
          Histo(24)%NBins  = 80
          Histo(24)%BinSize= 0.1d0
          Histo(24)%LowVal = -4d0
          Histo(24)%SetScale= 1d0

          Histo(25)%Info   = "eta(W,bbar)"
          Histo(25)%NBins  = 80
          Histo(25)%BinSize= 0.1d0
          Histo(25)%LowVal = -4d0
          Histo(25)%SetScale= 1d0

          Histo(26)%Info   = "eta(W,j)"
          Histo(26)%NBins  = 80
          Histo(26)%BinSize= 0.1d0
          Histo(26)%LowVal = -4d0
          Histo(26)%SetScale= 1d0

          Histo(27)%Info   = "Delta eta(H,j)"
          Histo(27)%NBins  = 100
          Histo(27)%BinSize= 0.1d0
          Histo(27)%LowVal = -5d0
          Histo(27)%SetScale= 1d0

          Histo(28)%Info   = "Delta y(H,j)"
          Histo(28)%NBins  = 100
          Histo(28)%BinSize= 0.1d0
          Histo(28)%LowVal = -5d0
          Histo(28)%SetScale= 1d0

          Histo(29)%Info   = "Delta eta(H,bbar)"
          Histo(29)%NBins  = 100
          Histo(29)%BinSize= 0.1d0
          Histo(29)%LowVal = -5d0
          Histo(29)%SetScale= 1d0

          Histo(30)%Info   = "Delta eta(W,j)"
          Histo(30)%NBins  = 100
          Histo(30)%BinSize= 0.1d0
          Histo(30)%LowVal = -5d0
          Histo(30)%SetScale= 1d0

          Histo(31)%Info   = "Delta eta(W,bbar)"
          Histo(31)%NBins  = 100
          Histo(31)%BinSize= 0.1d0
          Histo(31)%LowVal = -5d0
          Histo(31)%SetScale= 1d0

          Histo(32)%Info   = "Delta eta(H,em)"
          Histo(32)%NBins  = 100
          Histo(32)%BinSize= 0.1d0
          Histo(32)%LowVal = -5d0
          Histo(32)%SetScale= 1d0

          Histo(33)%Info   = "m(H,bbar,em)"
          Histo(33)%NBins  = 50
          Histo(33)%BinSize= 10d0*GeV
          Histo(33)%LowVal = 0d0*GeV
          Histo(33)%SetScale= 100d0

          Histo(34)%Info   = "m(H,j,em)"
          Histo(34)%NBins  = 50
          Histo(34)%BinSize= 10d0*GeV
          Histo(34)%LowVal = 0d0*GeV
          Histo(34)%SetScale= 100d0

          Histo(35)%Info   = "Delta eta(tbar,j)"
          Histo(35)%NBins  = 100
          Histo(35)%BinSize= 0.1d0
          Histo(35)%LowVal = -5d0
          Histo(35)%SetScale= 1d0

          Histo(36)%Info   = "Delta eta(tbar,H)"
          Histo(36)%NBins  = 100
          Histo(36)%BinSize= 0.1d0
          Histo(36)%LowVal = -5d0
          Histo(36)%SetScale= 1d0

          Histo(37)%Info   = "cos theta(l,j1)"
          Histo(37)%NBins  = 20
          Histo(37)%BinSize= 0.1d0
          Histo(37)%LowVal = -1d0
          Histo(37)%SetScale= 1d0
     
          
          
ELSE
    call Error("This ObsSet is not available ",ObsSet)
ENDIF


!DEC$ IF(_UseMPIVegas.EQ.1)
    if( NumHistograms.gt.NUMHISTO ) call Error("Number of histograms exceeds limit for MPI implemenation. Increase NUMHISTO in .f90 and .c",NumHistograms)
!DEC$ ENDIF


do NHisto=1,NumHistograms
      if( .not.allocated(Histo(NHisto)%Value) ) then
        allocate( Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Value2) ) then
        allocate( Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo(NHisto)%Hits) ) then
        allocate( Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      Histo(NHisto)%Value(0:Histo(NHisto)%NBins+1) = 0d0
      Histo(NHisto)%Value2(0:Histo(NHisto)%NBins+1)= 0d0
      Histo(NHisto)%Hits(0:Histo(NHisto)%NBins+1)  = 0
!DEC$ IF(_UseMPIVegas.EQ.1)
    if( Histo(NHisto)%NBins.gt.MXHISTOBINS ) call Error("Number of histo bins exceeds limit for MPI implemenation. Increase MXHISTOBINS in .f90 and .c",NHisto)
!DEC$ ENDIF
enddo


! RR 2015-04-15
do NHisto=1,Num2DHistograms
      if( .not.allocated(Histo2D(NHisto)%Value) ) then
        allocate( Histo2D(NHisto)%Value(0:Histo2D(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo2D(NHisto)%Value2) ) then
         allocate( Histo2D(NHisto)%Value2(0:Histo2D(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      if( .not.allocated(Histo2D(NHisto)%Hits) ) then
        allocate( Histo2D(NHisto)%Hits(0:Histo2D(NHisto)%NBins+1), stat=AllocStatus  )
        if( AllocStatus .ne. 0 ) call Error("Memory allocation in Histo")
      endif
      Histo2D(NHisto)%Value(0:Histo2D(NHisto)%NBins+1) = 0d0
      Histo2D(NHisto)%Value2(0:Histo2D(NHisto)%NBins+1)= 0d0
      Histo2D(NHisto)%Hits(0:Histo2D(NHisto)%NBins+1)  = 0
!DEC$ IF(_UseMPIVegas.EQ.1)                                                                                                                                                     
    if( Histo2D(NHisto)%NBins.gt.MXHISTOBINS ) call Error("Number of histo bins exceeds limit for MPI implemenation. Increase MXHISTOBINS in .f90 and .c",NHisto)
!DEC$ ENDIF                                                                                                                                                                     
enddo



RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_HTopDK(Topol,HTopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: HTopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)


   if( Topol.eq.HT_A0_T ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_A0,.false.,MomDK,PSWgt)
   elseif( Topol.eq.HT_BH_T ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_BH,.false.,MomDK,PSWgt)
   elseif( Topol.eq.HT_BH_T_G ) then
      call EvalPhasespace_HTopDecay(HTopMom,xRndPS,m_BH,.true.,MomDK,PSWgt)
   else
      call Error("EvalPS not yet implemented")
   endif


RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_HTopDecay(HTopMom,xRndPS,Mass,GluonRad,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,Mass
real(8) :: HTopMom(1:4),A0Mom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
      call genps(2,m_HTop,xRndPS(1:2),(/Mass,m_SMTop/),MomDK(1:4,1:2),PSWgt2)! Htop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,2),HTopMom(1:4),m_HTop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission

     call genps(3,m_HTop,xRndPS(1:5),(/Mass,m_SMTop,0d0/),MomDK(1:4,1:3),PSWgt2)! Htop decay with gluon

!          Pcol1= 5 -1
!          Pcol2= 5 -1
!          SingDepth = 1e-10
!          Steps = 20
!          call gensing(3,m_HTop,(/Mass,m_SMTop,0d0/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "running gensing"
!          PSWgt2=1d0

!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,2),HTopMom(1:4),m_HTop)
      call boost(MomDK(1:4,3),HTopMom(1:4),m_HTop)
      PSWgt = PSWgt2*PiWgt3

    endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_STopDK(Topol,STopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: STopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)


   if( Topol.eq.ST_Chi0_T ) then
      call EvalPhasespace_STopDecay(STopMom,xRndPS,.false.,MomDK,PSWgt)

   elseif( Topol.eq.ST_Chi0_T_G ) then
      call EvalPhasespace_STopDecay(STopMom,xRndPS,.true.,MomDK,PSWgt)
   else
      call Error("EvalPS not yet implemented")
   endif


RETURN
END SUBROUTINE



SUBROUTINE EvalPhasespace_STopDecay(STopMom,xRndPS,GluonRad,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2
real(8) :: STopMom(1:4),ChiMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
!     MomDK(1:4,i): i= 1:Chi, 2:top
      call genps(2,m_STop,xRndPS(1:2),(/m_Chi,m_Top/),MomDK(1:4,1:2),PSWgt2)! Stop decay
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,2),STopMom(1:4),m_STop)
      PSWgt = PSWgt2*PiWgt2

    else! extra gluon emission
!     MomDK(1:4,i): i= 1:Chi, 2:top, 3:gluon
      call genps(3,m_STop,xRndPS(1:5),(/m_Chi,m_Top,0d0/),MomDK(1:4,1:3),PSWgt2)! Stop decay with gluon

!          Pcol1= 5 -1
!          Pcol2= 5 -1
!          SingDepth = 1e-10
!          Steps = 20
!          call gensing(3,m_STop,(/m_Chi,m_Top,0d0/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "running gensing"
!          PSWgt2=1d0

!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,2),STopMom(1:4),m_STop)
      call boost(MomDK(1:4,3),STopMom(1:4),m_STop)
      PSWgt = PSWgt2*PiWgt3

      soft = dabs(MomDK(1,3))/m_STop
      if( soft.lt.1d-6  ) PSWgt = 0d0

    endif


RETURN
END SUBROUTINE







SUBROUTINE EvalPhasespace_TopDK(Topol,TopMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: Topol
real(8) :: PSWgt
real(8) :: TopMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)



   if( Topol.eq.T_B_W ) then
      call EvalPhasespace_TopDecay(TopMom,xRndPS,.false.,MomDK,PSWgt)
   elseif( Topol.eq.T_BG_W ) then
      call EvalPhasespace_TopDecay(TopMom,xRndPS,.true.,MomDK,PSWgt)
   elseif( Topol.eq.T_B_WG ) then
      call EvalPhasespace_TopDecay2(TopMom,xRndPS,.true.,MomDK,PSWgt)
   elseif( Topol.eq.T_BGG_W ) then
      call EvalPhasespace_TopDecay3(TopMom,xRndPS,MomDK,PSWgt)
   elseif( Topol.eq.T_BG_WG) then
      call EvalPhasespace_TopDecay4(TopMom,xRndPS,MomDK,PSWgt)
   elseif( Topol.eq.T_B_WGG) then
      call EvalPhasespace_TopDecay5(TopMom,xRndPS,MomDK,PSWgt)
   endif

RETURN
END SUBROUTINE






SUBROUTINE EvalPhasespace_TopDecay(TopMom,xRndPS,GluonRad,MomDK,PSWgt)!  top quark decay phase space with/without additional massless particle
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll
logical,save :: flip=.true.

    if( .not.GluonRad ) then!  no extra gluon radiation
!     MomDK(1:4,i): i= 1:bottom, 2:lepton, 3:neutrino
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      call genps(2,m_W,xRndPS(3:4),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      PSWgt = PSWgt2*PiWgt2 * PSWgt3*PiWgt2

    else! extra gluon emission
!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon
!      call genps(3,m_Top,xRndPS(1:5),(/0d0,0d0,m_W/),MomDK(1:4,1:3),PSWgt2)! top decay with additional gluon
!      WMom(1:4) = MomDK(1:4,3)
!      MomDK(1:4,4) = MomDK(1:4,2)
!      PSWgt = PSWgt2*PiWgt3

      call yeti3(m_Top,xRndPS(1:5),(/0d0,m_W,0d0/),MomDK(1:4,1:3),PSWgt2)
      WMom(1:4) = MomDK(1:4,2)
      MomDK(1:4,4) = MomDK(1:4,3)
      PSWgt = PSWgt2

! flip=.not.flip
! if( flip ) then!  every second call the singular event is generated
!         Pcol1= 3 -1
!         Pcol2= 4 -1
!         SingDepth = 1d-5
!         Steps = 10
!         call gensing(3,m_Top,(/0d0,0d0,m_W/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"
!         PSWgt2=1d0
!         WMom(1:4) = MomDK(1:4,3)
!         MomDK(1:4,4) = MomDK(1:4,2)
! endif

      soft = dabs(MomDK(1,4))/m_Top
      coll = dabs(MomDK(1:4,1).dot.MomDK(1:4,4))/m_Top**2
      if( soft.lt.1d-6 .or. coll.lt.1d-10 ) PSWgt = 0d0

      call genps(2,m_W,xRndPS(6:7),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
      PSWgt = PSWgt * PSWgt3*PiWgt2
    endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay2(TopMom,xRndPS,GluonRad,MomDK,PSWgt)!  top quark decay phase space with additional massless particle from W decay
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:3)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,s_glu_qb,s_glu_q,E_glu
logical,save :: flip=.true.

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      PSWgt = PSWgt2*PiWgt2

!       call genps(3,m_W,xRndPS(3:7),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! W decay
!       PSWgt = PSWgt * PSWgt3*PiWgt3
      call yeti3(m_W,xRndPS(3:7),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! W decay
      PSWgt = PSWgt * PSWgt3


! flip=.not.flip
! if( flip ) then!  every second call the singular event is generated
!         Pcol1= 5 -1
!         Pcol2= 5 -1
!         SingDepth = 1e-5
!         Steps = 5
!         call gensing(3,m_W,(/0d0,0d0,0d0/),MomDK(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"
!         PSWgt3=1d0
! endif

!     boost leptons to the W frame:
      call boost(MomDK(1:4,2),WMom(1:4),m_W)
      call boost(MomDK(1:4,3),WMom(1:4),m_W)
      call boost(MomDK(1:4,4),WMom(1:4),m_W)
!     boost all guys to the top frame:
      call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
      call boost(MomDK(1:4,4),TopMom(1:4),m_Top)

!     cut out very collinear/soft configurations
      s_glu_qb = dabs(MomDK(1:4,4).dot.MomDK(1:4,2))/m_W**2
      s_glu_q  = dabs(MomDK(1:4,4).dot.MomDK(1:4,3))/m_W**2
      E_glu    = dabs(MomDK(1,4))/m_W
      if( s_glu_qb.lt.1d-10 .or. s_glu_q.lt.1d-10 .or. E_glu.lt.1d-5) then
          PSWgt = 0d0
          return
      endif

RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_TopDecay3(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth, soft, coll

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon, 5: gluon/photon
    call genps(4,m_Top,xRndPS(1:8),(/m_W,0d0,0d0,0d0/),MomDK(1:4,1:4),PSWgt2)

!        Pcol1= 4 -1
!        Pcol2= 5 -1
!        SingDepth = 1e-10
!        Steps = 20
!        call gensing(4,m_Top,(/m_W,0d0,0d0,0d0/),MomDK(1:4,1:4),Pcol1,Pcol2,SingDepth,Steps)
!        PSWgt2=1d0

       WMom(1:4) = MomDK(1:4,1)
       MomDK(1:4,1) = MomDK(1:4,2)
       MomDK(1:4,5) = MomDK(1:4,4)
       MomDK(1:4,4) = MomDK(1:4,3)
       call genps(2,m_W,xRndPS(9:10),(/0d0,0d0/),MomDK(1:4,2:3),PSWgt3)! W decay
       !     boost leptons to the W frame:
       call boost(MomDK(1:4,2),WMom(1:4),m_W)
       call boost(MomDK(1:4,3),WMom(1:4),m_W)
       !     boost all guys to the top frame:
       call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
       call boost(MomDK(1:4,5),TopMom(1:4),m_Top)

! cut to singular events
   soft = dmin1(dabs(MomDK(1,4)),dabs(MomDK(1,5)))/TopMom(1)
   coll = dmin1(dabs(MomDk(1:4,1).dot.MomDk(1:4,4)),dabs(MomDk(1:4,1).dot.MomDk(1:4,5)),dabs(MomDk(1:4,4).dot.MomDk(1:4,5)))/TopMom(1)**2
! print *, "soft/coll",soft,coll
   if(soft.gt.1.d-6 .and. coll.gt.1.d-10 ) then
       PSWgt = PSWgt2*PiWgt4 * PSWgt3*PiWgt2
   else
!        MomDk(1:4,1:5) = 0.d0!   why is that necessary??  it leads to NaN in sq.mat.elements
       PSWgt = 0.d0
   endif


! if( MomDK(1,4)/M_top.lt.1d-3 ) PSWgt = 0d0
! if( MomDK(1,5)/M_top.lt.1d-3 ) PSWgt = 0d0
! if( dsqrt((MomDk(1:4,1).dot.MomDk(1:4,4))/m_top**2).lt.1d-3 ) PSWgt = 0d0!  has zero effect!!
! if( dsqrt((MomDk(1:4,1).dot.MomDk(1:4,5))/m_top**2).lt.1d-3 ) PSWgt = 0d0!  has zero effect!!
! if( dsqrt((MomDk(1:4,4).dot.MomDk(1:4,5))/m_top**2).lt.1d-3 ) PSWgt = 0d0


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay4(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles from top and W line, resp.
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4),MomChk(1:4,1:4)
real(8) :: MomDK(1:4,1:5)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth, soft, coll1, coll2


!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon/photon 5: gluon
    call genps(3,m_Top,xRndPS(1:5),(/0d0,0d0,m_W/),MomDK(1:4,1:3),PSWgt2)! top decay with additional gluon
!        Pcol1= 4 -1
!        Pcol2= 4 -1
!        SingDepth = 1e-15
!        Steps = 20
!        call gensing(3,m_Top,(/0d0,0d0,m_W/),MomDK(1:4,1:3),Pcol1,Pcol2,SingDepth,Steps)
!        PSWgt2=1d0

    soft = dabs(MomDK(1,2))/M_top
    coll1= dabs(MomDk(1:4,2).dot.MomDk(1:4,1))/m_top**2
    if(soft.gt.1d-6 .and. coll1.gt.1d-10) then!  reject too soft/collinear gluon configurations
        WMom(1:4) = MomDK(1:4,3)
        MomDK(1:4,5) = MomDK(1:4,2)
!         call genps(3,m_W,xRndPS(6:10),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)
!         PSWgt = PSWgt2*PiWgt3 * PSWgt3*PiWgt3
        call yeti3(m_W,xRndPS(6:10),(/0d0,0d0,0d0/),MomDK(1:4,2:4),PSWgt3)! highest sensitivity: 3-1-2
        PSWgt = PSWgt2*PiWgt3 * PSWgt3
!               Pcol1= 5 -1
!               Pcol2= 5 -1
!               SingDepth = 1e-15
!               Steps = 20
!               call gensing(3,m_W,(/0d0,0d0,0d0/),MomDK(1:4,2:4),Pcol1,Pcol2,SingDepth,Steps)
!               PSWgt3=1d0

        soft = dabs(MomDK(1,4))/m_W
        coll1= dabs(MomDk(1:4,2).dot.MomDk(1:4,4))/m_W**2
        coll2= dabs(MomDk(1:4,3).dot.MomDk(1:4,4))/m_W**2
        if(soft.gt.1d-5 .and. coll1.gt.1d-10 .and. coll2.gt.1d-10) then!  reject too soft/collinear gluon configurations
!            boost leptons to the W frame:
             call boost(MomDK(1:4,2),WMom(1:4),m_W)
             call boost(MomDK(1:4,3),WMom(1:4),m_W)
             call boost(MomDK(1:4,4),WMom(1:4),m_W)
!            boost all guys to the top frame:
             call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
             call boost(MomDK(1:4,5),TopMom(1:4),m_Top)
        else
             PSWgt = 0d0
        endif
   else
         PSWgt = 0d0
   endif


RETURN
END SUBROUTINE





SUBROUTINE EvalPhasespace_TopDecay5(TopMom,xRndPS,MomDK,PSWgt)!  top quark decay phase space with two additional massless particles from W decay
use ModProcess
use ModMisc
use ModParameters
implicit none
logical GluonRad
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: TopMom(1:4),WMom(1:4)
real(8) :: MomDK(:,:)
real(8) :: xRndPS(:)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
real(8),parameter :: N4=4, PiWgt4 = (2d0*Pi)**(4-N4*3) * (4d0*Pi)**(N4-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,soft,coll

!     MomDK(1:4,i): i=1:bottom, 2:lepton, 3:neutrino, 4: gluon, 5: gluon
      call genps(2,m_Top,xRndPS(1:2),(/0d0,m_W/),MomDK(1:4,1:2),PSWgt2)! top decay
      WMom(1:4) = MomDK(1:4,2)
      call genps(4,m_W,xRndPS(3:10),(/0d0,0d0,0d0,0d0/),MomDK(1:4,2:5),PSWgt3)! W decay

      soft = dmin1(dabs(MomDK(1,4)),dabs(MomDK(1,5)))/M_W
      coll = dmin1(dabs(MomDk(1:4,4).dot.MomDk(1:4,5)),dabs(MomDk(1:4,2).dot.MomDk(1:4,4)),dabs(MomDk(1:4,2).dot.MomDk(1:4,5)), &
                   dabs(MomDk(1:4,3).dot.MomDk(1:4,4)),dabs(MomDk(1:4,3).dot.MomDk(1:4,5)))/m_W**2
      if(soft.gt.1.d-6 .and. coll.gt.1.d-11) then!  reject too soft/collinear gluon configurations
!               Pcol1= 4 -1
!               Pcol2= 6 -1
!               SingDepth = 1e-15
!               Steps = 15
!               call gensing(4,m_W,(/0d0,0d0,0d0,0d0/),MomDK(1:4,2:5),Pcol1,Pcol2,SingDepth,Steps)
!               PSWgt3=1d0
!         boost leptons to the W frame:
          call boost(MomDK(1:4,2),WMom(1:4),m_W)
          call boost(MomDK(1:4,3),WMom(1:4),m_W)
          call boost(MomDK(1:4,4),WMom(1:4),m_W)
          call boost(MomDK(1:4,5),WMom(1:4),m_W)
!         boost all guys to the top frame:
          call boost(MomDK(1:4,1),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,2),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,3),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,4),TopMom(1:4),m_Top)
          call boost(MomDK(1:4,5),TopMom(1:4),m_Top)
          PSWgt = PSWgt2*PiWgt2 * PSWgt3*PiWgt4
       else
          PSWgt = 0d0
       endif



RETURN
END SUBROUTINE







SUBROUTINE EvalPhasespace_ZDecay(ZMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt
real(8) :: ZMom(1:4)
real(8) :: MomDK(1:4,1:2)
real(8) :: xRndPS(1:2),mZ_inv
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


        mZ_inv = dsqrt( dabs((ZMom(1:4).dot.ZMom)) )
        call genps(2,mZ_inv,xRndPS(1:2),(/0d0,0d0/),MomDK(1:4,1:2),PSWgt)! top decay

!       boost leptons to the Z frame:
        call boost(MomDK(1:4,1),ZMom(1:4),mZ_inv)
        call boost(MomDK(1:4,2),ZMom(1:4),mZ_inv)
        PSWgt = PSWgt*PiWgt2


RETURN
END SUBROUTINE




SUBROUTINE EvalPhasespace_HDecay(ZMom,xRndPS,MomDK,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: PSWgt
real(8) :: ZMom(1:4)
real(8) :: MomDK(1:4,1:2)
real(8) :: xRndPS(1:2),mZ_inv
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


        mZ_inv = dsqrt( dabs((ZMom(1:4).dot.ZMom)) )
        call genps(2,mZ_inv,xRndPS(1:2),(/0d0,0d0/),MomDK(1:4,1:2),PSWgt)! top decay

!       boost leptons to the Z frame:
        call boost(MomDK(1:4,1),ZMom(1:4),mZ_inv)
        call boost(MomDK(1:4,2),ZMom(1:4),mZ_inv)
        PSWgt = PSWgt*PiWgt2


RETURN
END SUBROUTINE










SUBROUTINE EvalPhaseSpace_2tobbbWW(EHat,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:12)
real(8) :: Mom(1:4,1:10),MomAux(1:4,3:8)
real(8),parameter :: PiWgtPr4 = (2d0*Pi)**(4-4*3) * (4d0*Pi)**(4-1)
real(8),parameter :: PiWgtPr2 = (2d0*Pi)**(4-2*3) * (4d0*Pi)**(2-1)


    call genps(4,Ehat,xRndPS(1:8),(/0d0,m_W,0d0,m_W/),MomAux(1:4,3:6),PSWgt)
    PSWgt = PSWgt*PiWgtPr4
    Mom(1:4,5)=MomAux(1:4,3)
    Mom(1:4,8)=MomAux(1:4,5)


    call genps(2,m_W,xRndPS(9:10),(/0d0,0d0/),Mom(1:4,6:7),PSWgt2)
    PSWgt = PSWgt*PSWgt2*PiWgtPr2
    call boost(Mom(1:4,6),MomAux(1:4,4),m_W)
    call boost(Mom(1:4,7),MomAux(1:4,4),m_W)

    call genps(2,m_W,xRndPS(11:12),(/0d0,0d0/),Mom(1:4,9:10),PSWgt3)
    PSWgt = PSWgt*PSWgt3*PiWgtPr2
    call boost(Mom(1:4,9),MomAux(1:4,6),m_W)
    call boost(Mom(1:4,10),MomAux(1:4,6),m_W)

    Mom(1:4,3)=Mom(1:4,5)+Mom(1:4,6)+Mom(1:4,7)
    Mom(1:4,4)=Mom(1:4,8)+Mom(1:4,9)+Mom(1:4,10)


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0





RETURN
END SUBROUTINE












SUBROUTINE EvalPhasespace_2to2(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
! integer :: NPart,i
! real(8) :: vel,theta
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_Top,m_Top/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!   call SetKinVars(4,Mom(1:4,1:4),(/0d0,0d0,m_Top,m_Top/))

!     include "kinpointDK.MCFM.f90"

return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to2HT(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_HTop,m_HTop/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2to3HT(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:5),MomW(1:4),xRndPS(1:5)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8) :: s13,s23


   call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Htop,m_Htop/),Mom(1:4,3:5),PSWgt)
   PSWgt = PSWgt*PiWgt3


!     Pcol1= 2 -1
!     Pcol2= 3 -1
!     SingDepth = 1e-10
!     Steps = 15
!     PSWgt = 1d0
!     call gensing(3,EHat,(/0d0,m_HTop,m_HTop/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "generating singular point"


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   s13 = Mom(1:4,1).dot.Mom(1:4,3)
   s23 = Mom(1:4,2).dot.Mom(1:4,3)
   if( abs(s13)/EHat**2.lt.1d-9 .or. abs(s23)/EHat**2.lt.1d-9 ) PSWgt=0d0
   if( abs(Mom(1,3)/EHat).lt.1d-5  ) PSWgt=0d0


return
END SUBROUTINE








SUBROUTINE EvalPhasespace_2to2Stops(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
integer,save :: it=1
real(8) :: beta,t,u,cos13,phi


   call genps(2,Ehat,xRndPS(1:2),(/m_Stop,m_Stop/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! print *, "INPUT MOMENTA FOR COMPARISON WITH RADJA"
! phi = 0.231231d0
! if(it.eq.1) then
!  EHat=dsqrt(1000000.0000000d0) *GeV
!  t = -836410.161513775d0 *GeV**2
!  u = -143589.838486225d0*GeV**2
! elseif(it.eq.2) then
!  EHat=dsqrt(1440000.00000000d0) *GeV
!  t = -1324817.04595758d0*GeV**2
!  u = -95182.9540424241d0*GeV**2
! elseif(it.eq.3) then
!  EHat=dsqrt(1960000.00000000d0) *GeV
!  t = -1454974.22611929d0*GeV**2
!  u = -485025.773880714d0*GeV**2
! endif
! 
!  PSWgt=1d0
!  m_Stop = dsqrt(0.5d0*(EHat**2+t+u))
! 
!  beta=dsqrt(1d0-m_stop**2/(EHat*0.5d0)**2)
!  cos13=1d0/beta*(1d0+(t-m_stop**2)/(EHat**2*0.5d0))
! 
! 
!  Mom(1,4) = EHat*0.5d0
!  Mom(2,4) = EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dsin(phi)
!  Mom(3,4) = EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dcos(phi)
!  Mom(4,4) = EHat*0.5d0*beta*cos13
! 
!  Mom(1,3) = EHat*0.5d0
!  Mom(2,3) =-EHat*beta*0.5d0*dsqrt(1d0-cos13**2)*dsin(phi)
!  Mom(3,3) =-EHat*0.5d0*beta*dsqrt(1d0-cos13**2)*dcos(phi)
!  Mom(4,3) =-EHat*beta*0.5d0*cos13
!  it=it+1
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 




!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


! print *, "check",((Mom(1:4,3)).dot.(Mom(1:4,3)))-m_stop**2
! print *, "check",((Mom(1:4,4)).dot.(Mom(1:4,4)))-m_stop**2
! print *, "check",Mom(1,1)+Mom(1,2)-Mom(1,3)-Mom(1,4)
! print *, "check",Mom(2,1)+Mom(2,2)-Mom(2,3)-Mom(2,4)
! print *, "check",Mom(3,1)+Mom(3,2)-Mom(3,3)-Mom(3,4)
! print *, "check",Mom(4,1)+Mom(4,2)-Mom(4,3)-Mom(4,4)
! print *, "check",ehat**2+t+u-2d0*m_stop**2
! print *, "t",(Mom(1:4,1)-Mom(1:4,3)).dot.(Mom(1:4,1)-Mom(1:4,3))
! print *, "u",(Mom(1:4,1)-Mom(1:4,4)).dot.(Mom(1:4,1)-Mom(1:4,4))
! pause


return
END SUBROUTINE






SUBROUTINE EvalPhasespace_2to3Stops(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:5),MomW(1:4),xRndPS(1:5)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8) :: s13,s23


   call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Stop,m_Stop/),Mom(1:4,3:5),PSWgt)
   PSWgt = PSWgt*PiWgt3

!     Pcol1= 3 -1
!     Pcol2= 3 -1
!     SingDepth = 1e-10
!     Steps = 15
!     PSWgt = 1d0
!     call gensing(3,EHat,(/0d0,m_sTop,m_sTop/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "generating singular point"

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   s13 = Mom(1:4,1).dot.Mom(1:4,3)
   s23 = Mom(1:4,2).dot.Mom(1:4,3)
   if( abs(s13)/EHat**2.lt.1d-9 .or. abs(s23)/EHat**2.lt.1d-9 ) PSWgt=0d0
   if( abs(Mom(1,3)/EHat).lt.1d-5  ) PSWgt=0d0


return
END SUBROUTINE






SUBROUTINE EvalPhasespace_2to3(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(3,Ehat,xRndPS(1:5),(/0d0,m_Top,m_Top/),Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,0d0/),Mom(1:4,3:5),PSWgt)
!   TmpMom(1:4) = Mom(1:4,3)
!   Mom(1:4,3)  = Mom(1:4,5)
!   Mom(1:4,5)  = TmpMom(1:4)

!      Pcol1= 3 -1
!      Pcol2= 3 -1
!      SingDepth = 1e-10
!      Steps = 10
!      PSWgt = 1d0
!      call gensing(3,EHat,(/0d0,m_Top,m_Top/),Mom(1:4,3:5),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing activated"

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!   call SetKinVars(5,Mom(1:4,1:5),(/0d0,0d0,0d0,m_Top,m_Top/))


!      print *, "using DUW kinpoint with mtop=174"
!      include './misc/DUWkinpoint'        ! Dittm.,Uwer,Weinzierl p1+p2 --> p3+p4+p5

!     print *, "using kinpoint3."
!     include 'misc/kinpoint3.'

!    include "kinpoint."        ! 0 --> p1+..+p4
!    Mom(1:4,1) = -Mom(1:4,1)     ! p1+p2 --> p3+p4
!    Mom(1:4,2) = -Mom(1:4,2)


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to3M(EHat,Mass,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Mass
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massive(M) + massive(anti-top) + massive(top)
  call genps(3,Ehat,xRndPS(1:5),(/Mass,m_Top,m_Top/),Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,Mass/),Mom(1:4,3:5),PSWgt)
!   TmpMom(1:4) = Mom(1:4,3)
!   Mom(1:4,3)  = Mom(1:4,5)
!   Mom(1:4,5)  = TmpMom(1:4)

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to3ArbMass(EHat,Mass,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Mass(1:3)
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)                                                                                                                                                     
! integer :: NPart,i                                                                                                                                                            
! real(8) :: vel,parx,theta ! for checks                                                                                                                                        
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=3, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


  call genps(3,Ehat,xRndPS(1:5),Mass,Mom(1:4,3:5),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   call yeti3(Ehat,xRndPS(1:5),(/m_Top,m_Top,Mass/),Mom(1:4,3:5),PSWgt)                                                                                                        
!   TmpMom(1:4) = Mom(1:4,3)                                                                                                                                                    
!   Mom(1:4,3)  = Mom(1:4,5)                                                                                                                                                    
!   Mom(1:4,5)  = TmpMom(1:4)                                                                                                                                                   

!  particles on the beam axis:                                                                                                                                                  
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE


SUBROUTINE EvalPhasespace_2to4M(EHat,Masses,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,Masses(1:4)
real(8) :: xRndPS(1:5)
real(8) :: Mom(1:4,1:5),TmpMom(1:4)
! real(8) :: MomDK(1:4,1:6)
! integer :: NPart,i
! real(8) :: vel,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth,velo,parx
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(M) + massive(anti-top) + massive(top)
  call genps(4,Ehat,xRndPS(1:8),Masses(1:4),Mom(1:4,3:6),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!   Pcol1= 1 -1
!   Pcol2= 3 -1
!   SingDepth = 1d-16
!   Steps = 20
!   PSWgt = 1d0
!   call gensing(4,EHat,Masses(1:4),Mom(1:4,3:6),Pcol1,Pcol2,SingDepth,Steps); print *, "gensing"


!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE



SUBROUTINE SmearExternal(xRnd,Mass,Width,MinEnergy,MaxEnergy,invMass,Jacobi)
use ModParameters
real(8) :: xRnd,Mass,Width,invMass,Jacobi,MinEnergy,MaxEnergy
real(8) :: r,rmin,rmax,BW


IF( Width.lt.(1d-6)*GeV ) THEN
   invMass = Mass
   Jacobi  = 1d0
ELSE
   rmin=1d0/(Width*Mass) * datan( (MinEnergy**2-Mass**2)/(Width*Mass)  )    
   rmax=1d0/(Width*Mass) * datan( (MaxEnergy**2-Mass**2)/(Width*Mass)  )
   r = xRnd*(rmax-rmin) + rmin
   
   invMass = dsqrt(dabs( Mass*Width * dtan( Mass*Width*r )  +  Mass**2 ))
!    BW = (invMass**2 - Mass**2)**2 + Mass**2 * Width**2            ! Breit-Wiegner propagator 
   BW = (Mass**2 * Width**2) * ( dtan( Mass*Width*r )**2 + 1d0 )    ! equivalent to above

   Jacobi  = (rmax-rmin) * BW   /(2d0*DblPi)                        ! this Jacobian has [GeV^2]; it becomes [GeV^-2] when hitting the sq. propagator 
ENDIF                                                               ! factor 2*Pi comes from integr. measure


return
END SUBROUTINE




!!! Zprime section !!!

SUBROUTINE EvalPhasespaceBWMapp(EHat,Masses,xRndPS,Mom,PSWgt)                                                                                                
use ModParameters                                                                                                                                            
use ModMisc                                                                                                                                                  
implicit none                                                                                                                                                
real(8) :: EHat,Masses(1:3),xRndPS(1:3*3-4)                                                                                                                  
real(8) :: Mom(1:4,1:5),PSWgt                                                                                                                                
real(8) :: PiWgt2,SingDepth                                                                                                                                  
integer :: N,NVar,Pcol1,Pcol2,Steps                                                                                                                          
                                                                                                                                                             
   N=3                                                                                                                                                       
   NVar=3*N-4                                                                                                                                                
   PiWgt2 = (2d0*Pi)**(-NVar) * (4d0*Pi)**(N-1)                                                                                                              
   call genpszttbg(N,Ehat,xRndPS(1:NVar),M_Zpr,Ga_Zpr,m_Top,(/0d0,m_top,m_top/),Mom(1:4,3:5),PSWgt)                                                          
   PSWgt = PSWgt*PiWgt2                                                                                                                                      
                                                                                                                                                             
   Mom(1,1) =  EHat*0.5d0                                                                                                                                    
   Mom(2,1) =  0d0                                                                                                                                           
   Mom(3,1) =  0d0                                                                                                                                           
   Mom(4,1) = +EHat*0.5d0                                                                                                                                    
                                                                                                                                                             
   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

   call swapmom(Mom(1:4,3),Mom(1:4,5))

   if( N.eq.3 ) then
      if( dmin1( (Mom(1:4,1).dot.Mom(1:4,5))/EHat**2,(Mom(1:4,2).dot.Mom(1:4,5))/EHat**2,(Mom(1,5)/EHat)**2 ).lt.1d-10 ) PSWgt = 0d0
   endif


return
END SUBROUTINE


!!! End Zprime section !!!



SUBROUTINE EvalPhasespace_2to4(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:6)!,MomTmp(1:4,3:6)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massless + massive(anti-top) + massive(top)
    call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),PSWgt)
! call genps(4,Ehat,xRndPS(1:8),(/m_Top,m_Top,0d0,0d0/),MomTmp(1:4,3:6),PSWgt)
    PSWgt = PSWgt*PiWgtPr

! Mom(1:4,3)=MomTmp(1:4,5)
! Mom(1:4,4)=MomTmp(1:4,6)
! Mom(1:4,5)=MomTmp(1:4,3)
! Mom(1:4,6)=MomTmp(1:4,4)

!    Pcol1= 1 -1
!    Pcol2= 3 -1
!    SingDepth = 1d-10
!    Steps = 10
!    PSWgt = 1d0
!    call gensing(4,EHat,(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),Pcol1,Pcol2,SingDepth,Steps)

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0

!     write(*,"(4F25.16)") Mom(1:4,1)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,2)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,3)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,4)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,5)  *1d0
!     write(*,"(4F25.16)") Mom(1:4,6)  *1d0
!
! !     print *, (Mom(1:4,5).dot.Mom(1:4,6))/(Mom(1:4,1).dot.Mom(1:4,2))
!     print *, (Mom(1:4,3).dot.Mom(1:4,4))/(Mom(1:4,1).dot.Mom(1:4,2))
!     pause
!


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to5(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:11)
real(8) :: Mom(1:4,1:7)
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=4, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massless + massive(anti-top) + massive(top)
 call genps(5,Ehat,xRndPS(1:11),(/0d0,0d0,0d0,m_Top,m_Top/),Mom(1:4,3:7),PSWgt)
 PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to6(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
use ifport
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:20)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: NPart,i
real(8) :: velo,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8),parameter :: NPr=6, PiWgtPr = (2d0*Pi)**(4-NPr*3) * (4d0*Pi)**(NPr-1)


!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(6,Ehat,xRndPS(1:14),(/0d0,0d0,0d0,0d0,0d0,0d0/),Mom(1:4,3:8),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE




SUBROUTINE EvalPhasespace_2toN(N,EHat,xRndPS,Mom,Mass,PSWgt)
use ModProcess
use ModMisc
use ModParameters
use ifport
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:(3*N+4))
real(8) :: Mom(1:4,1:N),Mass(1:N)
integer :: NPart,i,N
real(8) :: velo,parx,theta ! for checks
integer :: Pcol1,Pcol2,Steps
real(8) :: SingDepth
real(8) :: PiWgtPr


   PiWgtPr= (2d0*Pi)**(4-N*3) * (4d0*Pi)**(N-1)
!  generate PS: massless + massless --> massless + massive(anti-top) + massive(top)
  call genps(N,Ehat,xRndPS(1:(3*N+4)),Mass,Mom(1:4,3:N+2),PSWgt)
  PSWgt = PSWgt*PiWgtPr

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE







SUBROUTINE EvalPhasespace_2to4ML(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:5),MomDK(1:4,1:6)
integer :: NPart,i

!  generate PS: massless + massless --> massless +  massless +  massless +  massless
   call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,0d0,0d0/),Mom(1:4,3:6),PSWgt)
!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE



SUBROUTINE EvalPhasespace_2to4ZDK(EHat,xRndPS,Mom,PSWgt)
use ModProcess
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3
real(8) :: xRndPS(1:8)
real(8) :: Mom(1:4,1:6),MomDK(1:4,1:6)
integer :: NPart,i

!  generate PS: massless + massless --> massless +  massless +  massive(anti-top) +  massive(Top)
   call genps(4,Ehat,xRndPS(1:8),(/0d0,0d0,m_Top,m_Top/),Mom(1:4,3:6),PSWgt)
!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE EvalPhasespace_2to4ZDK











SUBROUTINE Kinematics_TTBARJET(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS
real(8) :: MomExt(1:4,1:5+NPlus1PS),MomDK(1:4,1:6),MomJet(1:4,1:8),zeros(1:10),MomJet_ordered(1:4,1:8)
real(8) :: MomHadr(1:4,1:8),MomLept(1:4,1:4)
real(8) :: MomBoost(1:4),MomW(1:4),MomTops(1:4,1:2)
logical :: applyPSCut
integer :: NBin(:),PartList(1:8),JetList(1:8),NJet,NObsJet,n,NObsJet_Tree,nWJets
real(8) :: pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top,HT,m_lb,R_lb,m_bb,m_bj
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,eta_miss,beta,costheta
real(8) :: pT_jet(1:8),eta_jet(1:8),eta_sepa,eta_Zeppi,s34,s35,s36,s45,s46,s56,mTopHadr,mTopLept
real(8) :: R_bb,MinvJets,MinvLept,phi_Lept,pT_lept,ET_lept,ET_miss,mT,pT_x,pT_y,MTW
real(8) :: MomTTbar(1:4),pT_ttbar,m_ttbar,y_top,y_Atop,y_ttbar,dy_tops,dphi_ttbar


!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if( NPlus1PS.eq.0 ) then
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   else
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5) - MomExt(1:4,6)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,4) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   endif
   if( any(abs(zeros(1:4)/MomExt(1,1)).gt.1d-6) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:4)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( NPlus1PS.eq.0 ) then
        zeros(1) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(2) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
   else
        zeros(1) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
        zeros(2) = (MomExt(1:4,6).dot.MomExt(1:4,6)) - m_Top**2
        zeros(10)=  MomExt(1:4,4).dot.MomExt(1:4,4)
   endif
   zeros(3) =  MomExt(1:4,3).dot.MomExt(1:4,3)
   zeros(4) =  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(5) =  MomDK(1:4,2).dot.MomDK(1:4,2)
   zeros(6) =  MomDK(1:4,3).dot.MomDK(1:4,3)
   zeros(7) =  MomDK(1:4,4).dot.MomDK(1:4,4)
   zeros(8) =  MomDK(1:4,5).dot.MomDK(1:4,5)
   zeros(9) =  MomDK(1:4,6).dot.MomDK(1:4,6)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:10)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:10)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:10)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARJET(",NPlus1PS,"): ",zeros(1:10)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
!DEC$ ENDIF







! NumPart = particles in the final state
! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light partons, N-1=ATop, N=Top
!                          MomDK(:,1:6) : 1=ABot, 2=lep-/q, 3=ANeu/qbar, 4=Bot, 5=lep+/qbar, 6=Neu/q

applyPSCut = .false.
NBin(1:NumHistograms) = 0

MomHadr(1:4,1:8) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:8)=(/1,2,3,4,5,6,7,8/)


! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   MomHadr(1:4,1) = MomExt(1:4,3)   ! q/qbar/glu
   if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,2) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
   else
      NumHadr = 1
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   endif
!-------------------------------------------------------
elseif( TopDecays.eq.1 ) then  ! full leptonic decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu

  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,4) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomHadr(1:4,4) = MomDK(1:4,2)  ! q
  MomHadr(1:4,5) = MomDK(1:4,3)  ! qbar
  MomHadr(1:4,6) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,7) = MomDK(1:4,6)  ! q
  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,8) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 8
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,4) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,5) = MomDK(1:4,6)  ! q
  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,6) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
  MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
  MomHadr(1:4,3) = MomExt(1:4,3) ! q/qbar/glu
  MomHadr(1:4,4) = MomDK(1:4,2)  ! q
  MomHadr(1:4,5) = MomDK(1:4,3)  ! qbar
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu

  if( NPlus1PS.eq.1 ) then
      MomHadr(1:4,6) = MomExt(1:4,4)   ! q/qbar/glu
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,5)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,6)   ! Top
  else
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  endif

else
  call Error("this TopDecay is not implemented in Kinematics_TTBARJET")
endif


!---------------------- (anti) kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero

    NJet=0
    MomJet(1:4,1:8) = MomHadr(1:4,1:8)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets

!-------------------------------------------------------------------------
pT_jet(1:8)  = 0d0
eta_jet(1:8) = 0d0

if( ObsSet.eq.10 .or. ObsSet.eq.11 ) then! set of observables for ttbjet production without decays at Tevatron & LHC

    call pT_order(NJet,MomJet(1:4,1:NJet))! pT ordering of jet momenta

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

    pT_jet(1) = get_PT(MomJet(1:4,1))
    eta_jet(1)= get_ETA(MomJet(1:4,1))

    pT_jet(2) = get_PT(MomJet(1:4,2))
    eta_jet(2)= get_ETA(MomJet(1:4,2))

! check cuts
    NObsJet_Tree = 1
    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!     if( pT_Top.lt.800d0*GeV ) then!   this is for the boosted observable
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_jet(1))
    NBin(8) = WhichBin(8,eta_jet(1))


!-------------------------------------------------------
elseif( ObsSet.eq.12 ) then! set of observables for ttbjet production as signal process at Tevatron (semi-lept decay)

    NObsJet_Tree = 5

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_miss  = get_PT(MomLept(1:4,4))

    phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))

    HT = pT_lepP + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo




!    check cuts
     if( pT_lepP.lt.pT_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( abs(eta_lepP).gt.eta_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( pT_miss.lt.pT_miss_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( HT.lt.HT_cut ) then
         applyPSCut = .true.
         RETURN
     endif


! binning

    NBin(:) = 0

    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,pT_jet(3))
    NBin(4) = WhichBin(4,eta_jet(3))

    MomJet_ordered(1:4,1:NObsJet) = MomJet(1:4,1:NObsJet)
    call pT_order(NObsJet,MomJet_ordered(1:4,1:NObsJet))! pT ordering of jet momenta for b AND non-b jets
    NBin(5) = WhichBin(5,get_pT(MomJet_ordered(1:4,5)))
    NBin(6) = WhichBin(6,get_eta(MomJet_ordered(1:4,5)))

    NBin(7) = WhichBin(7,pT_miss)
    NBin(8) = WhichBin(8,HT)
    NBin(9) = WhichBin(9,m_lb)
    NBin(10)= WhichBin(10,phi_Lept)
    NBin(11)= WhichBin(11,R_lb)
    NBin(12)= WhichBin(12,eta_lepP)




! additional histograms for A_FB analysis
! ** ideally reconstructed tops
    MomTTbar(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    pT_ttbar = get_PT(MomTTbar(1:4))
    m_ttbar = get_MInv(MomTTbar(1:4))
    y_ttbar = get_ETA(MomTTbar(1:4))
    y_top = get_ETA(MomTops(1:4,2))
    y_ATop = get_ETA(MomTops(1:4,1))
    dy_tops  = y_top-y_ATop
    pT_top = get_PT(MomTops(1:4,2))
    dphi_ttbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2)) )
    if( dphi_ttbar.gt.Pi ) dphi_ttbar=2d0*Pi-dphi_ttbar
    beta = dsqrt( abs(1d0-m_top**2/MomTops(1,2)**2) )
    costheta = Get_CosTheta(MomTops(1:4,2))


    NBin(13)= WhichBin(13,pT_ttbar)
    if(dy_tops.ge.0d0) NBin(14)= WhichBin(14,pT_ttbar)
    if(dy_tops.lt.0d0) NBin(15)= WhichBin(15,pT_ttbar)
    if(dy_tops.ge.0d0) NBin(16)= WhichBin(16,m_ttbar)
    if(dy_tops.lt.0d0) NBin(17)= WhichBin(17,m_ttbar)
    if(dy_tops.ge.0d0) NBin(18)= WhichBin(18,y_ttbar)
    if(dy_tops.lt.0d0) NBin(19)= WhichBin(19,y_ttbar)
    if(dy_tops.ge.0d0) NBin(20)= WhichBin(20,abs(dy_tops))
    if(dy_tops.lt.0d0) NBin(21)= WhichBin(21,abs(dy_tops))
    if(eta_lepP.ge.0d0) NBin(22)= WhichBin(22,pT_ttbar)
    if(eta_lepP.lt.0d0) NBin(23)= WhichBin(23,pT_ttbar)
    NBin(24)= WhichBin(24,y_top)
    NBin(25)= WhichBin(25,pT_top)
    if(dy_tops.ge.0d0) NBin(26)= WhichBin(26,dphi_ttbar)
    if(dy_tops.lt.0d0) NBin(27)= WhichBin(27,dphi_ttbar)

    NBin(43) = WhichBin(43,costheta)
    NBin(44) = WhichBin(44,beta*costheta)



! additional histograms for A_FB analysis
! ** realistically reconstructed tops


! reconstruct top (decays leptonically)
! reconstruct anti-top (decays hadronically)
!   find two non-b jets that are closest to MW mass
    if( NObsJet.eq.6 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s36= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,6))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s46= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,6))-M_W )
        s56= dabs( get_MInv(MomJet(1:4,5)+MomJet(1:4,6))-M_W )
    elseif( NObsJet.eq.5 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s36=1d10; s46=1d10; s56=1d10
    else
        call Error("this should not happen")
    endif
    nWJets=minloc((/s34,s35,s45,s36,s46,s56/),1)

!   construct hadr. W momentum
    if(nWJets.eq.1) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,4)
    elseif(nWJets.eq.2) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,5)
    elseif(nWJets.eq.3) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,5)
    elseif(nWJets.eq.4) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,6)
    elseif(nWJets.eq.5) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,6)
    elseif(nWJets.eq.6) then
        MomW(1:4) = MomJet(1:4,5)+MomJet(1:4,6)
    else
        MomW(1:4) = 0d0
    endif

    if( get_R(MomJet(1:4,1),MomLept(1:4,3)) .lt. get_R(MomJet(1:4,2),MomLept(1:4,3))  ) then ! find smaller R-distance between lepton and bjet
        MomTops(1:4,2) = MomJet(1:4,1) + MomLept(1:4,3)+MomLept(1:4,4)
        MomTops(1:4,1) = MomJet(1:4,2) + MomW(1:4)
    else
        MomTops(1:4,2) = MomJet(1:4,2) + MomLept(1:4,3)+MomLept(1:4,4)
        MomTops(1:4,1) = MomJet(1:4,1) + MomW(1:4)
    endif

! print *, "check j",NObsJet
! print *, "check W",get_MInv(MomW(1:4))*100d0
! print *, "check T",get_MInv(MomTops(1:4,1))*100d0,get_MInv(MomTops(1:4,2))*100d0
! print *, "xxx",get_MInv(MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4))*100d0
! pause

!   require a 30 GeV window around M_W
!   require a 50 GeV window around M_Top
    if( dabs(get_MInv(MomW(1:4))-M_W).lt.30d0*GeV .and. & 
        dabs(get_MInv(MomTops(1:4,1))-M_Top).lt.50d0*GeV .and. dabs(get_MInv(MomTops(1:4,2))-M_Top).lt.50d0*GeV) then

        MomTTbar(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
        pT_ttbar = get_PT(MomTTbar(1:4))
        m_ttbar = get_MInv(MomTTbar(1:4))
        y_ttbar = get_ETA(MomTTbar(1:4))
        y_top = get_ETA(MomTops(1:4,2))
        y_ATop = get_ETA(MomTops(1:4,1))
        dy_tops  = y_top-y_ATop
        pT_top = get_PT(MomTops(1:4,2))
        dphi_ttbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2)) )
        if( dphi_ttbar.gt.Pi ) dphi_ttbar=2d0*Pi-dphi_ttbar

        NBin(28)= WhichBin(28,pT_ttbar)
        if(dy_tops.ge.0d0) NBin(29)= WhichBin(29,pT_ttbar)
        if(dy_tops.lt.0d0) NBin(30)= WhichBin(30,pT_ttbar)
        if(dy_tops.ge.0d0) NBin(31)= WhichBin(31,m_ttbar)
        if(dy_tops.lt.0d0) NBin(32)= WhichBin(32,m_ttbar)
        if(dy_tops.ge.0d0) NBin(33)= WhichBin(33,y_ttbar)
        if(dy_tops.lt.0d0) NBin(34)= WhichBin(34,y_ttbar)
        if(dy_tops.ge.0d0) NBin(35)= WhichBin(35,abs(dy_tops))
        if(dy_tops.lt.0d0) NBin(36)= WhichBin(36,abs(dy_tops))
        if(eta_lepP.ge.0d0) NBin(37)= WhichBin(37,pT_ttbar)
        if(eta_lepP.lt.0d0) NBin(38)= WhichBin(38,pT_ttbar)
        NBin(39)= WhichBin(39,y_top)
        NBin(40)= WhichBin(40,pT_top)
        if(dy_tops.ge.0d0) NBin(41)= WhichBin(41,dphi_ttbar)
        if(dy_tops.lt.0d0) NBin(42)= WhichBin(42,dphi_ttbar)


    endif






!-------------------------------------------------------
elseif( ObsSet.eq.13 ) then! set of observables for ttbjet production as signal process at LHC (di-lept decay)

    NObsJet_Tree = 3

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_lepM  = get_PT(MomLept(1:4,1))
    eta_lepM = get_eta(MomLept(1:4,1))

    pT_miss  = get_PT(MomLept(1:4,2)+MomLept(1:4,4))


    HT = pT_lepP + pT_lepM + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo

    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    m_bb = get_MInv(MomJet(1:4,1)+MomJet(1:4,2))
    m_bj = get_MInv(MomJet(1:4,1)+MomJet(1:4,3))

    phi_Lept = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    pT_Top = get_PT( MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,1)+MomLept(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) )



! check cuts

    if( pT_lepP.lt.pT_lep_cut .or. pT_lepM.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut .or. abs(eta_lepM).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepM)
    NBin(2) = WhichBin(2,eta_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,HT)
    NBin(6) = WhichBin(6,MinvLept)
    NBin(7) = WhichBin(7,pT_jet(3))
    NBin(8) = WhichBin(8,eta_jet(3))
    NBin(9) = WhichBin(9,pT_miss)
    NBin(10) = WhichBin(10,pT_Top)
    NBin(11) = WhichBin(11,phi_Lept)
    NBin(12) = WhichBin(12,MInvLept)
    NBin(13) = WhichBin(13,pt_jet(1))
    NBin(14) = WhichBin(14,eta_jet(1))
    NBin(15) = WhichBin(15,m_bb)
    NBin(16) = WhichBin(16,m_bj)



!-------------------------------------------------------
elseif( ObsSet.eq.14 ) then! set of observables for ttbjet production as background process to VBF at LHC (di-lept decay)

print *, "STOP! proper jet selectin needs to be implemented first! see ObsSet=12,13"; stop

    NObsJet_Tree = 3
!   request at least two b-jets and one non-b-jet
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    call pT_order(NJet,MomJet(1:4,1:NJet))! pT ordering of jet momenta for b- and non-b-jets, Note: b-jets are no longer in position 1 & 2

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_lepM  = get_PT(MomLept(1:4,1))
    eta_lepM = get_eta(MomLept(1:4,1))

    do n=1,NJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
    enddo

    MinvJets = Get_MInv(MomJet(1:4,1)+MomJet(1:4,2))
    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!     phi_Lept = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
!     if( phi_Lept.gt.Pi ) phi_Lept=2d0*DblPi-phi_Lept
    phi_Lept = dacos( (MomLept(2,1)*MomLept(2,3)+MomLept(3,1)*MomLept(3,3))/pT_lepP/pT_lepM )  ! dacos el [0,Pi]
!PRINT *, "CHECK EQUIVALENCE!"
!STOP

    pT_lept = get_PT(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    ET_lept = dsqrt(pT_lept**2 + MinvLept**2)
    ET_miss = dsqrt(pT_miss**2 + MinvLept**2)
    pT_x = MomLept(2,1)+MomLept(2,2)+MomLept(2,3)+MomLept(2,4)
    pT_y = MomLept(3,1)+MomLept(3,2)+MomLept(3,3)+MomLept(3,4)
    mT = dsqrt( (ET_lept+ET_miss)**2 - (pT_x**2+pT_y**2) )
    eta_sepa = abs(eta_jet(1)-eta_jet(2))


! check cuts
    if( pT_jet(1).lt. pT_hardestjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo

    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

    if( eta_jet(1)*eta_jet(2).gt.0d0 .or. eta_sepa.lt.eta_sepa_cut ) then
        applyPSCut = .true.
        RETURN
    endif


!   a "veto jet" must lie in between the two tagged jets
    if( abs(eta_jet(3)-eta_jet(1)).lt.eta_sepa .and. abs(eta_jet(3)-eta_jet(2)).lt.eta_sepa ) then  ! jet(3) is "veto jet"
        eta_Zeppi = eta_jet(3) - 0.5d0*(eta_jet(1)+eta_jet(2))
    else  ! there is no "veto jet"
        eta_Zeppi = -100d0
    endif

!     do n=1,NJet
!         if( pT_jet(n).gt.pT_jet_cut .and. abs(eta_jet(n)).lt.eta_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     enddo

    if( pT_lepP.lt.pT_lep_cut .or. pT_lepM.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut .or. abs(eta_lepM).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( MinvJets.lt.Minv_jets_cut ) then
        applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepM)
    NBin(2) = WhichBin(2,eta_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,MinvLept)
    NBin(6) = WhichBin(6,pT_jet(1))
    NBin(7) = WhichBin(7,phi_Lept)
    NBin(8) = WhichBin(8,mT)
    NBin(9) = WhichBin(9,eta_jet(1))
    NBin(10) = WhichBin(10,eta_jet(2))
    NBin(11) = WhichBin(11,eta_Zeppi)




elseif( ObsSet.eq.15 ) then! set of observables for ttbjet production as signal process at LHC (semi-lept decay)

    NObsJet_Tree = 5

!   check that there are at least NObsJet_Tree resolved jets
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

!   check that b jets pass pT and y cut
    if( get_pT(MomJet(1:4,1)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,1))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( get_pT(MomJet(1:4,2)).lt.pT_bjet_cut .or. abs(get_eta(MomJet(1:4,2))).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   calculate number of observable jets
    NObsJet = 2 != two b-jets that already passed all criteria
    do n=3,NJet
        if( get_pT(MomJet(1:4,n)).gt.pT_jet_cut .and. abs(get_eta(MomJet(1:4,n))).lt.eta_jet_cut ) then
             NObsJet = NObsJet +1
            if( n.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,n)
        endif
    enddo
    MomJet(1:4,NObsJet+1:NJet) = 0d0!  purging the remaining array

!   check that there are at least NObsJet_Tree observable jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif





!   evaluate kinematic variables

    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

    pT_miss  = get_PT(MomLept(1:4,4))

    MTW = dsqrt( (MomLept(2,3)+MomLept(2,4))**2 + (MomLept(3,3)+MomLept(3,4))**2 )

    phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
    if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept

    m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
    R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))

    HT = pT_lepP + pT_miss
    do n=1,NObsJet
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
      HT = HT + pT_jet(n)
    enddo



!    check cuts
     if( pT_lepP.lt.pT_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( abs(eta_lepP).gt.eta_lep_cut ) then
         applyPSCut = .true.
         RETURN
     endif

     if( pT_miss.lt.pT_miss_cut ) then
         applyPSCut = .true.
         RETURN
     endif

!     if( pT_miss+MTW.lt.60d0*GeV ) then!   additional cut for ATLAS muon analysis
     if( MTW.lt.30d0*GeV ) then!   additional cut for ATLAS electron analysis
         applyPSCut = .true.
         RETURN
     endif


!   find two non-b jets that are closest to MW mass
    if( NObsJet.eq.6 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s36= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,6))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s46= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,6))-M_W )
        s56= dabs( get_MInv(MomJet(1:4,5)+MomJet(1:4,6))-M_W )
    elseif( NObsJet.eq.5 ) then
        s34= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W )
        s35= dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,5))-M_W )
        s45= dabs( get_MInv(MomJet(1:4,4)+MomJet(1:4,5))-M_W )
        s36=1d10; s46=1d10; s56=1d10
    endif
    nWJets=minloc((/s34,s35,s45,s36,s46,s56/),1)

!   construct hadr. W momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4)
    if(nWJets.eq.1) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,4)
    elseif(nWJets.eq.2) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,5)
    elseif(nWJets.eq.3) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,5)
    elseif(nWJets.eq.4) then
        MomW(1:4) = MomJet(1:4,3)+MomJet(1:4,6)
    elseif(nWJets.eq.5) then
        MomW(1:4) = MomJet(1:4,4)+MomJet(1:4,6)
    elseif(nWJets.eq.6) then
        MomW(1:4) = MomJet(1:4,5)+MomJet(1:4,6)
    else
        MomW(1:4) = 0d0
    endif
    MomTops(1:4,1) = MomTops(1:4,1) + MomW(1:4)!   construct the t+bar system
    if( dmin1(s34,s35,s45,s36,s46,s56).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,pT_jet(3))
    NBin(4) = WhichBin(4,eta_jet(3))

    call pT_order(NObsJet,MomJet(1:4,1:NObsJet))! pT ordering of jet momenta for b AND non-b jets
    NBin(5) = WhichBin(5,get_pT(MomJet(1:4,5)))
    NBin(6) = WhichBin(6,get_eta(MomJet(1:4,5)))

    NBin(7) = WhichBin(7,pT_miss)
    NBin(8) = WhichBin(8,HT)
    NBin(9) = WhichBin(9,m_lb)
    NBin(10)= WhichBin(10,phi_Lept)
    NBin(11)= WhichBin(11,R_lb)
    NBin(12)= WhichBin(12,eta_lepP)
    NBin(13)= WhichBin(13,pT_Top)




!-------------------------------------------------------
elseif( ObsSet.eq.19 ) then! for checks of ttbjet



!DEC$ IF(_FactCheck .EQ.1)


!    this is for "checkCD"
! if( NPlus1PS.eq.0) then
!     pT_jet(3)  = get_pT(MomExt(1:4,3))  ! leading non-b jet
!     if( pT_jet(3).lt.pT_jet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! else
!     if( get_R(MomExt(1:4,3),MomExt(1:4,4)).lt.0.4d0 ) then
!         pT_jet(3)  = get_pT(MomExt(1:4,3)+MomExt(1:4,4))
!         if( pT_jet(3).lt.pT_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     else
!         pT_jet(3)  = get_pT(MomExt(1:4,3))
!         pT_jet(4)  = get_pT(MomExt(1:4,4))
!         if( pT_jet(3).lt.pT_jet_cut .and. pT_jet(4).lt.pT_jet_cut ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     endif
! endif

RETURN!     this is needed to avoid the cuts below
!DEC$ ENDIF






    NObsJet_Tree = 3

    if( TopDecays.eq.2 ) NObsJet_Tree = NObsJet_Tree + 4
    if( TopDecays.eq.3 .or. TopDecays.eq.4 ) NObsJet_Tree = NObsJet_Tree + 2
    if( NJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif
    pT_lepP  = get_PT(MomLept(1:4,3))
    eta_lepP = get_eta(MomLept(1:4,3))

!     pT_miss  = get_PT(MomLept(1:4,4))
!     phi_Lept = dabs( Get_PHI(MomLept(1:4,3)) - Get_PHI(MomJet(1:4,1)) )
!     if( phi_Lept.gt.Pi ) phi_Lept=2d0*Pi-phi_Lept
!     m_lb = get_MInv(MomLept(1:4,3)+MomJet(1:4,1))
!     R_lb = get_R(MomLept(1:4,3),MomJet(1:4,1))


    do n=1,NJet! first two jets are always b-jets
      pT_jet(n)  = get_pT( MomJet(1:4,n))
      eta_jet(n) = get_eta(MomJet(1:4,n))
    enddo

!     HT = pT_lepP + pT_miss
!     do n=1,NJet
!         HT = HT + pT_jet(n)
!     enddo


! check cuts
    NObsJet = 0
    do n=1,NJet
        if( pT_jet(n).gt.pT_jet_cut ) NObsJet = NObsJet +1
    enddo
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif
    if( TOPDECAYS.eq.1 ) then
        pT_lepP  = get_PT(MomLept(1:4,3))
        eta_lepP = get_eta(MomLept(1:4,3))
        pT_lepM  = get_PT(MomLept(1:4,1))
        eta_lepM = get_eta(MomLept(1:4,1))
        pT_miss  = get_PT (MomLept(1:4,2)+MomLept(1:4,4))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepP.lt.0.2d0 .or. pT_lepM.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepP).gt.2.1d0 .or. abs(eta_lepM).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif



    if( TOPDECAYS.eq.3 ) then
        pT_lepM  = get_PT(MomLept(1:4,1))
        eta_lepM = get_eta(MomLept(1:4,1))
        pT_miss  = get_PT (MomLept(1:4,2))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepM.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepM).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif



    if( TOPDECAYS.eq.4 ) then
        pT_lepP  = get_PT(MomLept(1:4,3))
        eta_lepP = get_eta(MomLept(1:4,3))
        pT_miss  = get_PT (MomLept(1:4,4))
        if( pT_miss.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( pT_lepP.lt.0.2d0 ) then
            applyPSCut = .true.
            RETURN
        endif
        if( abs(eta_lepP).gt.2.1d0 ) then
            applyPSCut = .true.
            RETURN
        endif
    endif


! binning
    NBin(1) = WhichBin(1,pT_lepP)
    NBin(2) = WhichBin(2,eta_lepP)
    NBin(3) = WhichBin(3,get_pT(MomJet(1:4,1)))
    NBin(4) = WhichBin(4,get_eta(MomJet(1:4,1)))
    NBin(5) = WhichBin(5,pT_miss)
    NBin(6) = WhichBin(6,HT)
    NBin(7) = WhichBin(7,m_lb)
    NBin(8)= WhichBin(8,phi_Lept)
    NBin(9)= WhichBin(9,R_lb)
    NBin(10)= WhichBin(10,eta_lepP)





!-------------------------------------------------------
else
  print *, "ObsSet not implemented TTBJET",ObsSet
  stop
endif


return
END SUBROUTINE







SUBROUTINE Kinematics_TTBARPHOTON(NPlus1PS,Mom,MomOrder,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS,MomOrder(1:12)
real(8) :: Mom(1:4,1:12),zeros(1:12)
real(8) :: MomJet(1:4,1:7),MomJet_CHECK(1:4,1:7)
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4)
logical :: applyPSCut,isolated
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree,NJet_CHECK
real(8) :: pT_lepM,pT_lepP,ET_miss,pT_ATop,pT_Top,HT,ET_bjet,eta_CP
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,m_lb,m_jj,mTblP,m_jjb,m_jjbP,mT_lp
real(8) :: pT_jet(1:7),eta_jet(1:7),eta_sepa,pt_Pho,eta_Pho,Rphobjet,mT_bln(1:2),mT_blnp(1:2)
real(8) :: R_Pj(1:5),R_lj(1:5),R_PlepP,R_PlepM,pT_lept,ET_lept,mT,MInvPb1jj,mTb2lP,MInvPb2jj,mTb1lP,Phi_LP,Phi_LL
integer :: tbar,t,pho,inLeft,inRight,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,L,N




! momentum ordering
  tbar    = MomOrder(1)
  t       = MomOrder(2)
  pho     = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu


!DEC$ IF(_CheckMomenta .EQ.1)
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,tbar) - Mom(1:4,t) - Mom(1:4,pho)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,pho) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:4)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   zeros(3) =  Mom(1:4,pho).dot.Mom(1:4,pho)
   zeros(4) =  Mom(1:4,bbar).dot.Mom(1:4,bbar)
   zeros(5) =  Mom(1:4,lepM).dot.Mom(1:4,lepM)
   zeros(6) =  Mom(1:4,nubar).dot.Mom(1:4,nubar)
   zeros(7) =  Mom(1:4,b).dot.Mom(1:4,b)
   zeros(8) =  Mom(1:4,lepP).dot.Mom(1:4,lepP)
   zeros(9) =  Mom(1:4,nu).dot.Mom(1:4,nu)
   if( NPlus1PS.eq.1 ) zeros(10)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:10)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:10)
      print *, Mom(1:4,1:2)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:10)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARPHOTON(): ",zeros(1:10)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
!DEC$ ENDIF


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)

MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
    if(NPlus1PS.eq.0) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   L = LepM
   N = nubar
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   L = LepP
   N = nu
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF



!------------------ Frixione photon isolation ----------------------------

    isolated = FrixioneIsolated(Mom(1:4,pho),Rsep_Pj,NumHadr,MomHadr(1:4,1:NumHadr))
    if( .not. isolated ) then
        applyPSCut = .true.
        RETURN
    endif


!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))


!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.20 .or. ObsSet.eq.21) then! ttb+photon production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    eta_CP = dabs(eta_Top)-dabs(eta_ATop)

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

!     if(NumHadr.eq.1) then
!         if( dabs(dacos((MomHadr(2,1)*Mom(2,pho)+MomHadr(3,1)*Mom(3,pho)+MomHadr(4,1)*Mom(4,pho))/MomHadr(1,1)/Mom(1,pho))).lt.3d0/360d0*2d0*DblPi ) then !   for Chinese check
!             applyPSCut = .true.
!             RETURN
!         endif
!         R_Pj(1)  = get_R(Mom(1:4,pho),MomHadr(1:4,1))
!         if( R_Pj(1).lt.Rsep_Pj ) then
!             applyPSCut = .true.
!             RETURN
!         endif
!     endif

! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,pT_Pho)
    NBin(6) = WhichBin(6,eta_Pho)
    NBin(7) = WhichBin(7,eta_ATop)
    NBin(8) = WhichBin(8,eta_Top)
    NBin(9) = WhichBin(9,eta_CP)

!-------------------------------------------------------
elseif( ObsSet.eq.22 .or. ObsSet.eq.23 ) then! set of observables for ttb+gamma production with di-lept. decays at the Tevatron & LHC
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))

!   determine observable jets
    NObsJet = 0
    do k=1,NJet
        pT_jet(k)  = get_PT(MomJet(1:4,k))
        eta_jet(k) = get_ETA(MomJet(1:4,k))
        
        if( pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
        endif
    enddo

!     NObsJet_Tree = 4! request two b-jets and at least two light jets
!     if( NObsJet.lt.NObsJet_Tree ) then
!         applyPSCut = .true.
!         RETURN
!     endif




    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepM = get_R(Mom(1:4,pho),Mom(1:4,lepM))
    R_PlepP = get_R(Mom(1:4,pho),Mom(1:4,lepP))

    pT_lepM = get_PT(Mom(1:4,lepM))
    pT_lepP = get_PT(Mom(1:4,lepP))
    eta_lepM = get_ETA(Mom(1:4,lepM))
    eta_lepP = get_ETA(Mom(1:4,lepP))

ET_miss=0.1d0
HT= 0.1d0
m_lb=0.1d0
phi_ll=0.5d0

! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.eta_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    do k=1,NObsJet
        if( get_R(MomJet(1:4,k),Mom(1:4,pho)).lt.Rsep_Pj  )then
            applyPSCut = .true.
            RETURN
        endif
    enddo


!    if( pT_lepM.gt.pT_lepP .and. pT_lepM.lt.pT_lep_cut ) then
!        applyPSCut = .true.
!        RETURN
!    endif
!    if( pT_lepP.gt.pT_lepM .and. pT_lepP.lt.pT_lep_cut ) then
!        applyPSCut = .true.
!        RETURN
!    endif
    if( pT_lepM.lt.pT_lep_cut .and. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif


    if( abs(eta_lepM).gt.eta_lep_cut .or. abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( R_PlepM.lt.Rsep_Plep .or. R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif


!     if( ET_miss.lt.pT_miss_cut ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!     if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!     if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
!     NBin(5) = WhichBin(5,eta_ATop)! Tevatron
    NBin(5) = WhichBin(5,dabs(eta_Top)-dabs(eta_ATop)  )! LHC
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11)= WhichBin(11,ET_miss)
    NBin(12)= WhichBin(12,HT)
    NBin(13)= WhichBin(13,m_lb)
    NBin(14)= WhichBin(14,Phi_LL)



elseif( ObsSet.eq.24 .or. ObsSet.eq.25 ) then! set of observables for ttb+gamma production with semi-lept. at the Tevatron/LHC
! elseif( ObsSet.eq.24  ) then! set of observables for ttb+gamma production with semi-lept. at the Tevatron/LHC


    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))

!   determine observable jets
    NObsJet = 0
    do k=1,NJet
        pT_jet(k)  = get_PT(MomJet(1:4,k))
        eta_jet(k) = get_ETA(MomJet(1:4,k))
        
        if( pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
        endif
    enddo

!     NObsJet_Tree = 4! request two b-jets and at least two light jets
!     if( NObsJet.lt.NObsJet_Tree ) then
!         applyPSCut = .true.
!         RETURN
!     endif




    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepP = get_R(Mom(1:4,pho),Mom(1:4,L))

    pT_lepP = get_PT(Mom(1:4,L))
    eta_lepP = get_ETA(Mom(1:4,L))

ET_miss=0.1d0
HT= 0.1d0
m_lb=0.1d0
phi_ll=0.5d0

! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.eta_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    do k=1,NObsJet
        if( get_R(MomJet(1:4,k),Mom(1:4,pho)).lt.Rsep_Pj  )then
            applyPSCut = .true.
            RETURN
        endif
    enddo
   
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif


!     if( ET_miss.lt.pT_miss_cut ) then
!        applyPSCut = .true.
!         RETURN
!     endif

!     if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!     if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
!     NBin(5) = WhichBin(5,eta_ATop)! Tevatron
    NBin(5) = WhichBin(5,dabs(eta_Top)-dabs(eta_ATop)  )! LHC
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11)= WhichBin(11,ET_miss)
    NBin(12)= WhichBin(12,HT)
    NBin(13)= WhichBin(13,m_lb)
    NBin(14)= WhichBin(14,Phi_LL)



elseif( ObsSet.eq.26 .or. ObsSet.eq.27 .or. ObsSet.eq.28 ) then
! elseif( ObsSet.eq.25 .or. ObsSet.eq.26 .or. ObsSet.eq.27 .or. ObsSet.eq.28 ) then! set of observables for ttb+gamma production with semi-lept. decays(hadr.Atop, lept.top decay) at the Tevatron/LHC
                                                                                 ! ObsSet 26,27 include suppression cuts for photon radiation from top decay


                                                               
                                                                                                                            
!   request two separated b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    
    do k=1,NJet
      pt_jet(k)  = get_PT(MomJet(1:4,k))
      eta_jet(k) = get_eta(MomJet(1:4,k))
      R_Pj(k)    = get_R(MomJet(1:4,k),Mom(1:4,pho))
    enddo
    
    

!   request b-jets to be outside the Frixione cone
    if( R_Pj(1).lt.Rsep_Pj .or. R_Pj(2).lt.Rsep_Pj ) then
        applyPSCut = .true.
        RETURN
    endif

!   check if b-jets pass cuts
    if(  pT_jet(1).lt.pT_bjet_cut .or. abs(eta_jet(1)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_jet(2).lt.pT_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_jet(1) + pT_jet(2)


!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( R_Pj(k).gt.Rsep_Pj .and. pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside Frixione cone and outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + pT_jet(k)
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))
    eta_CP = dabs(eta_Top)-dabs(eta_ATop)


    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

    R_PlepP  = get_R(Mom(1:4,pho),Mom(1:4,L))

    Phi_LP = dabs( Get_PHI(Mom(1:4,pho)) - Get_PHI(Mom(1:4,L)) )
    if( Phi_LP.gt.Pi ) Phi_LP=2d0*Pi-Phi_LP

    pT_lepP  = get_PT(Mom(1:4,L))
    eta_lepP = get_ETA(Mom(1:4,L))

!     MomObs(1:4) = MomObs(1:4) + Mom(1:4,pho) + Mom(1:4,L)
!     MomMiss(1:4) = Mom(1:4,inLeft) + Mom(1:4,inRight) - MomObs(1:4)
    MomMiss(1:4) = Mom(1:4,N)
    ET_miss  = get_ET(MomMiss(1:4))
    HT = HT + pT_lepP + pT_Pho + ET_miss


    m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))    ! these are the pT-hardest b-jets
    mT = get_MT(Mom(1:4,L),MomMiss(1:4))
    Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))


!     if( dabs(get_Minv(Mom(1:4,pho)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,pho)+MomJet(1:4,2))) ) then
!         m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
!     else
!         m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,2))
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))
!     endif

!     Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))!  min.pT b-jet
!     m_lb = get_Minv(Mom(1:4,L)+MomJet(1:4,1))
!     if( get_R(Mom(1:4,pho),MomJet(1:4,1)).lt.get_R(Mom(1:4,pho),MomJet(1:4,2)) ) then
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,1))
!     else
!         Rphobjet = get_R(Mom(1:4,pho),MomJet(1:4,2))
!     endif


! check cuts
    if(pt_Pho.lt.pT_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if(abs(eta_Pho).gt.eta_pho_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( R_PlepP.lt.Rsep_Plep ) then
        applyPSCut = .true.
        RETURN
    endif

    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    m_jj = get_Minv(MomJet(1:4,3)+MomJet(1:4,4))! 3 4 is the pair
    
    ! using m_lb to associate lepton and b-jets
    if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
        mTblP = get_MT(MomJet(1:4,1)+Mom(1:4,pho)+Mom(1:4,L),MomMiss(1:4))
        m_jjbP= get_mInv(MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,4)+Mom(1:4,pho))
        m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,4))
        
        if( get_Minv(MomJet(1:4,1)+Mom(1:4,L)+MomMiss(1:4)+Mom(1:4,pho)) .lt. get_Minv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,4)+Mom(1:4,pho))  ) then
            R_Pj(1) = get_R(MomJet(1:4,1),Mom(1:4,pho))
        else
            R_Pj(1) = get_R(MomJet(1:4,2),Mom(1:4,pho))
        endif
    else
        mTblP = get_MT(MomJet(1:4,2)+Mom(1:4,pho)+Mom(1:4,L),MomMiss(1:4))
        m_jjbP= get_mInv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,4)+Mom(1:4,pho))
        m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,4))
        if( get_Minv(MomJet(1:4,2)+Mom(1:4,L)+MomMiss(1:4)+Mom(1:4,pho)) .lt. get_Minv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,4)+Mom(1:4,pho))  ) then
            R_Pj(1) = get_R(MomJet(1:4,2),Mom(1:4,pho))
        else
            R_Pj(1) = get_R(MomJet(1:4,1),Mom(1:4,pho))
        endif
    endif

    if(NObsJet.eq.5) then
         if( dabs(get_Minv(MomJet(1:4,3)+MomJet(1:4,5))-m_W) .lt. dabs(m_jj-m_W) ) then
              m_jj = get_Minv(MomJet(1:4,3)+MomJet(1:4,5))! 3 5 is the pair
              if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
                  m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,5))
              else
                  m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,3)+MomJet(1:4,5))
              endif
         endif
         if( dabs(get_Minv(MomJet(1:4,4)+MomJet(1:4,5))-m_W) .lt. dabs(m_jj-m_W) ) then
              m_jj = get_Minv(MomJet(1:4,4)+MomJet(1:4,5))! 4 5 is the pair
              if( dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,1))).lt.dabs(get_Minv(Mom(1:4,L)+MomJet(1:4,2))) ) then
                  m_jjb = get_mInv(MomJet(1:4,2)+MomJet(1:4,4)+MomJet(1:4,5))
              else
                  m_jjb = get_mInv(MomJet(1:4,1)+MomJet(1:4,4)+MomJet(1:4,5))
              endif
         endif
    endif

    
   ! this is Ehat for <EHat> in NHisto=12
   MInv_LB= get_MInv(Mom(1:4,inLeft)+Mom(1:4,inRight))


IF( OBSSET.EQ.27 ) THEN!   these are the cuts to suppress photon radiation from top quarks and W bosons
    
!     if( m_jjb.lt.160d0*GeV .or. m_jjb.gt.180d0*GeV) then
!          applyPSCut = .true.
!         RETURN
!     endif
!     if( (mTblP.lt.180d0*GeV) ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! 
! 
     if( m_jj.lt.75d0*GeV .or. m_jj.gt.86d0*GeV) then
             applyPSCut = .true.
             RETURN
     endif

     m_jj = get_Minv(Mom(1:4,L)+MomMiss(1:4))
     if( m_jj.lt.75d0*GeV .or. m_jj.gt.86d0*GeV) then
             applyPSCut = .true.
             RETURN
     endif

     if( mTblP.gt.175d0*GeV ) then!   this is to enhance the decay contribution
         applyPSCut = .true.
         RETURN
     endif



ELSEIF( OBSSET.EQ.28 ) THEN!   these are the cuts to suppress photon radiation from top quarks and W bosons
    
!     if( m_jjb.lt.160d0*GeV .or. m_jjb.gt.180d0*GeV) then
!          applyPSCut = .true.
!         RETURN
!     endif
!     if( (mTblP.lt.180d0*GeV) ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! 
! 
     if( m_jj.lt.75d0*GeV .or. m_jj.gt.86d0*GeV) then
             applyPSCut = .true.
             RETURN
     endif

     m_jj = get_Minv(Mom(1:4,L)+MomMiss(1:4))
     if( m_jj.lt.75d0*GeV .or. m_jj.gt.86d0*GeV) then
             applyPSCut = .true.
             RETURN
     endif

     if( mTblP.lt.175d0*GeV ) then!   this is to suppress the decay contribution
         applyPSCut = .true.
         RETURN
     endif

     
          
     
     
ENDIF





! binning
    NBin(1:4) = 0
    if( 0.5d0*(pT_ATop+pT_Top).lt.800d0*GeV ) then
        NBin(1) = WhichBin(1,0.5d0*(pT_ATop+pT_Top)  )
        NBin(2) = WhichBin(2,0.5d0*(eta_ATop+eta_Top))
    endif
!     NBin(1) = WhichBin(1,pT_ATop)
!     NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_CP)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_Pho)
    NBin(8) = WhichBin(8,eta_Pho)
    NBin(9) = WhichBin(9,pT_lepP)
    NBin(10)= WhichBin(10,eta_lepP)
    NBin(11) = WhichBin(11,ET_miss)
    NBin(12) = WhichBin(12,HT)
    NBin(13) = WhichBin(13,R_Pj(1)) !Rphobjet)
    NBin(14) = WhichBin(14,m_lb)
    NBin(15) = WhichBin(15,Phi_LP)
    NBin(16) = WhichBin(16,m_jjbP)
    NBin(17) = WhichBin(17,m_jj)
    NBin(18) = WhichBin(18,mTblP)
    NBin(19) = WhichBin(19,mT)


elseif( ObsSet.eq.29) then! ttb+photon production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    eta_ATop = get_ETA(Mom(1:4,tbar))
    eta_Top  = get_ETA(Mom(1:4,t))

    pT_Pho  = get_PT(Mom(1:4,pho))
    eta_Pho = get_ETA(Mom(1:4,pho))

! check cuts
!     if(pt_Pho.lt.pT_pho_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif

!     R_Pj(1)  = get_R(Mom(1:4,pho),Mom(1:4,realp))
!     if( (Process.eq.24 .or. Process.eq.26) .and. R_Pj(1).lt.Rsep_Pj ) then
!        applyPSCut = .true.
!        RETURN
!     endif


!     if( dabs(Mom(1:4,bbar).dot.Mom(1:4,pho))/m_Top**2.lt.1d-3 .or. dabs(Mom(1:4,b).dot.Mom(1:4,pho))/m_Top**2.lt.1d-3) then
!         applyPSCut = .true.
!         RETURN
!     endif


! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,pT_Pho)
    NBin(6) = WhichBin(6,eta_Pho)
    NBin(7) = WhichBin(7,eta_ATop)
    NBin(8) = WhichBin(8,eta_Top)


!-------------------------------------------------------
else
  print *, "ObsSet not implemented TTBPHOTON",ObsSet
  stop
endif


return
END SUBROUTINE





SUBROUTINE Kinematics_TTBARZ(NPlus1PS,Mom,MomOrder,applyPSCut,NBin,PObs)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,NPlus1PS,MomOrder(1:14)
real(8) :: Mom(1:4,1:14),zeros(1:14)
real(8) :: MomJet(1:4,1:7) !,MomJet_CHECK(1:4,1:7)
real(8) :: MomHadr(1:4,0:8),MomZ(1:4),MomFermZframe(1:4)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4)
logical :: applyPSCut
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree,leptj(1:3),i,j,ii,jj
real(8),optional :: PObs(:)
real(8) :: nZLept,s12,s13,s14,s23,s24,s34
real(8) :: pT_lep(4),ET_miss,PT_miss,pT_ATop,pT_Top,HT,ET_bjet
real(8) :: eta_ATop,eta_Top,eta_lep(1:4),pseudo_Z, pseudo_top, pseudo_tbar
real(8) :: pT_jet(1:7),eta_jet(1:7),eta_sepa,mT_bln(1:2),pT_Z,eta_Z
real(8) :: R_lj(1:5),R_PlepM,pT_lept,ET_lept,mT,dPhiLL,CosTheta1,DphiZt,Dphittbar
integer :: tbar,t,Zbos,inLeft,inRight,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,L,N,Zl,Za,ferm_Z,Aferm_Z,jlabel
real(8) :: pT_ll,HT_jet,WithinCone(1:3),RLept,Minv_Z,sqrtshat
integer :: iLept,jLept,jJet,JetIndex(1:4),LepIndex(1:3)
real(8) :: mT2,pA(2:4),pB(2:4),pTInvis(2:4),mA,mB,mInvis! this is for MT2 calculation
real(8) :: pseudoEta_Z,cosTheta,gamma_boost,beta_boost,cosThetaStar 

applyPSCut = .false.
if( Process.eq.81 ) return!  return for Z => photon


! momentum ordering
  tbar    = MomOrder(1)
  t       = MomOrder(2)
  Zbos    = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)
  ferm_Z  = MomOrder(13)!      fermion from Z decay (lep-, q, nu)
  Aferm_Z = MomOrder(14)! anti-fermion from Z decay (lep+, qbar, nubar)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu


!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,tbar) - Mom(1:4,t) - Mom(1:4,ZBos)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ferm_Z) - Mom(1:4,Aferm_Z) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZ(): ",NPlus1PS,zeros(1:4)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   if (ZDecays .le. 10) then   ! Z onshell
      zeros(3) = (Mom(1:4,ZBos).dot.Mom(1:4,ZBos)) - M_Z**2
   else ! Z off shell, so this check is meaningless
      zeros(3)=0d0
   endif
   zeros(4) =  Mom(1:4,bbar).dot.Mom(1:4,bbar)
   zeros(5) =  Mom(1:4,lepM).dot.Mom(1:4,lepM)
   zeros(6) =  Mom(1:4,nubar).dot.Mom(1:4,nubar)
   zeros(7) =  Mom(1:4,b).dot.Mom(1:4,b)
   zeros(8) =  Mom(1:4,lepP).dot.Mom(1:4,lepP)
   zeros(9) =  Mom(1:4,nu).dot.Mom(1:4,nu)
   zeros(10) = (Mom(1:4,ferm_Z).dot.Mom(1:4,ferm_Z))
   zeros(11) = (Mom(1:4,Aferm_Z).dot.Mom(1:4,Aferm_Z))


   if( NPlus1PS.eq.1 ) zeros(12)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:3)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:3)
      print *, Mom(1:4,1:3)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:12)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:12)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
!DEC$ ENDIF


!!! RR 2013/03/16 : havent made any changes for Z decay from here on - it looks like it shoul just be top decay though (other than some cuts on Z leptons -- much later)

NBin(1:NumHistograms) = 0
MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)
MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
    if(NPlus1PS.eq.0) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 .OR. TOPDECAYS.eq.-1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   L = LepM
   N = nubar
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   L = LepP
   N = nu
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF





!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))





!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.51) then! ttb+Z production without top decays at Tevatron & LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))


! binning
    NBin(1) = WhichBin(1,pT_Top)




elseif( ObsSet.eq.52 .or. ObsSet.eq.55 ) then! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay )


!     s12 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))-M_Z )
!     s13 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,ferm_Z))-M_Z )
!     s14 = dabs( get_MInv(Mom(1:4,lepM)+Mom(1:4,Aferm_Z))-M_Z )
!     s23 = dabs( get_MInv(Mom(1:4,lepP)+Mom(1:4,ferm_Z))-M_Z )
!     s24 = dabs( get_MInv(Mom(1:4,lepP)+Mom(1:4,Aferm_Z))-M_Z )
!     s34 = dabs( get_MInv(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))-M_Z )
!     nZLept=minloc((/s12,s13,s14,s23,s24,s34/),1)


    pT_jet(1) = get_PT(Mom(1:4,b))
    pT_jet(2) = get_PT(Mom(1:4,bbar))
    eta_jet(1) = get_ETA(Mom(1:4,b))
    eta_jet(2) = get_ETA(Mom(1:4,bbar))

    eta_Z = get_ETA(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))
    pT_top = get_PT(Mom(1:4,t))
    pT_atop = get_PT(Mom(1:4,tbar))
    eta_top = get_ETA(Mom(1:4,t))
    eta_atop = get_ETA(Mom(1:4,tbar))

    pseudo_Z=get_pseudoETA(Mom(1:4,Zbos))
    pseudo_top=get_pseudoETA(Mom(1:4,t))
    pseudo_tbar=get_pseudoETA(Mom(1:4,tbar))

    pT_Lep(1)  = get_PT(Mom(1:4,ferm_Z))
    pT_Lep(2)  = get_PT(Mom(1:4,Aferm_Z))
    pT_Lep(3)  = get_PT(Mom(1:4,LepP))
    pT_Lep(4)  = get_PT(Mom(1:4,LepM))
    eta_Lep(1)  = get_ETA(Mom(1:4,ferm_Z))
    eta_Lep(2)  = get_ETA(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,LepP))
    eta_Lep(4) = get_ETA(Mom(1:4,LepM))

    pT_miss = get_PT(Mom(1:4,nu)+Mom(1:4,nubar))
!!!!
!!!!    pT_miss = get_PT(Mom(1:4,nu)+Mom(1:4,nubar) + Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))! assuming Z-->inv
!!!!

    pT_Z  = get_PT(Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z))

    DphiLL = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,Aferm_Z))  )
    if( DphiLL.gt.Pi ) DphiLL=2d0*Pi-DphiLL



 CosTheta1 = Get_CosAlpha( Mom(1:4,t),Mom(1:4,Zbos) )
 sqrtshat = dsqrt( 2d0*(Mom(1:4,1).dot.Mom(1:4,2)) )



! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_lep(1).lt.pT_lep_cut .OR. pT_lep(2).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif




!     mA = get_MInv(Mom(1:4,LepP)+MomJet(1:4,1))
!     mB = get_MInv(Mom(1:4,LepP)+MomJet(1:4,2))
!     if( mA.lt.mB ) then 
!         MInv_LB = mA
!         mA = get_MInv(Mom(1:4,LepM)+MomJet(1:4,2))
!         mB = get_MInv(Mom(1:4,LepP)+MomJet(1:4,1))
!         pA(2:4) = (/Mom(2,LepM)+MomJet(2,2),Mom(3,LepM)+MomJet(3,2), mA /)*100d0
!         pB(2:4) = (/Mom(2,LepP)+MomJet(2,1),Mom(3,LepP)+MomJet(3,1), mB /)*100d0
!     else
!         MInv_LB = mB
!         mA = get_MInv(Mom(1:4,LepM)+MomJet(1:4,1))
!         mB = get_MInv(Mom(1:4,LepP)+MomJet(1:4,2))
!         pA(2:4) = (/Mom(2,LepM)+MomJet(2,1),Mom(3,LepM)+MomJet(3,1), mA /)*100d0
!         pB(2:4) = (/Mom(2,LepP)+MomJet(2,2),Mom(3,LepP)+MomJet(3,2), mB /)*100d0
!     endif
!!   ! calculate MT2 with pA,pB calculated above
!     mInvis = 0d0
!     pTInvis(2:4) = (/Mom(2,nu)+Mom(2,nubar)+Mom(2,ferm_Z)+Mom(2,Aferm_Z),Mom(3,nu)+Mom(3,nubar)+Mom(2,ferm_Z)+Mom(2,Aferm_Z), mInvis/) *100d0
!!     call calcMT2( pA(2:3), pB(2:3), pTInvis(2:4),MInv_T2 ) !   this is the call to the full MT2 library
!     MInv_T2=MInv_T2/100d0




! binning
    NBin(1) = WhichBin(1,pT_Lep(1))
    NBin(2) = WhichBin(2,pT_Lep(2))
    NBin(3) = WhichBin(3,pT_Z)
    NBin(4) = WhichBin(4,pT_top)
    NBin(5) = WhichBin(5,pT_atop)
    NBin(6) = WhichBin(6,pT_jet(1))
    NBin(7) = WhichBin(7,pT_jet(2))
    NBin(8) = WhichBin(8,pT_Lep(3))
    NBin(9) = WhichBin(9,pT_Lep(4))
    NBin(10) = WhichBin(10,pT_miss)
    NBin(11) = WhichBin(11,eta_Lep(1))
    NBin(12) = WhichBin(12,eta_Lep(2))
    NBin(13) = WhichBin(13,eta_Z)
    NBin(14) = WhichBin(14,eta_top)
    NBin(15) = WhichBin(15,eta_atop)
    NBin(16) = WhichBin(16,eta_jet(1))
    NBin(17) = WhichBin(17,eta_jet(2))
    NBin(18) = WhichBin(18,eta_Lep(3))
    NBin(19) = WhichBin(19,eta_Lep(4))
    NBin(20) = WhichBin(20,dphill)
    NBin(21) = WhichBin(21,pseudo_Z)
    NBin(22) = WhichBin(22,pseudo_top)
    NBin(23) = WhichBin(23,pseudo_tbar)
    NBin(24) = WhichBin(24,MInv_LB)
    NBin(25) = WhichBin(25,MInv_T2)


    NBin(28) = WhichBin(28,sqrtshat)
    NBin(29) = WhichBin(29,pT_top)
    NBin(30) = WhichBin(30, dacos(CosTheta1) )



    if( present(PObs) ) then
      PObs(1) = pT_Lep(1)
      PObs(2) = pT_Lep(2)
      PObs(3) = pT_Z
      PObs(4) = pT_top
      PObs(5) = pT_atop
      PObs(6) = pT_jet(1)
      PObs(7) = pT_jet(2)
      PObs(8) = pT_jet(3)
      PObs(9) = pT_jet(4)
      PObs(10) = pT_miss
      PObs(11) = eta_Lep(1)
      PObs(12) = eta_Lep(2)
      PObs(13) = eta_Z
      PObs(14) = eta_top
      PObs(15) = eta_atop
      PObs(16) = eta_jet(1)
      PObs(17) = eta_jet(1)
      PObs(18) = eta_Lep(3)
      PObs(19) = eta_Lep(4)
      PObs(20) = dphill
      PObs(21) = pseudo_Z
      PObs(22) = pseudo_top
      PObs(23) = pseudo_tbar
      PObs(24) = MInv_LB
      PObs(25) = MInv_T2

      PObs(28) = sqrtshat
      PObs(29) = pT_top
      PObs(30) = dacos(CosTheta1)
    endif



elseif( ObsSet.EQ.53 .or. ObsSet.EQ.56 .or. ObsSet.EQ.58 ) then! set of observables for ttb+Z ( semi-lept. ttbar decays and di-lept. Z decay )

! request at least four jets where two are b-jets
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    pT_jet(3) = get_PT(MomJet(1:4,3))
    pT_jet(4) = get_PT(MomJet(1:4,4))
    eta_jet(1) = get_ETA(MomJet(1:4,1))
    eta_jet(2) = get_ETA(MomJet(1:4,2))
    eta_jet(3) = get_ETA(MomJet(1:4,3))
    eta_jet(4) = get_ETA(MomJet(1:4,4))


    pT_Lep(1)  = get_PT(Mom(1:4,LepP))
    eta_Lep(1) = get_ETA(Mom(1:4,LepP))

    pT_Lep(2)  = get_PT(Mom(1:4,ferm_Z))
    eta_Lep(2) = get_ETA(Mom(1:4,ferm_Z))

    pT_Lep(3)  = get_PT(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,Aferm_Z))

    pT_miss = get_PT(Mom(1:4,nu))

    MomZ(1:4) = Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z)
    pT_Z  = get_PT(MomZ(1:4))
    eta_Z = get_eta(MomZ(1:4))
    Minv_Z = get_MInv(MomZ(1:4))

    pT_Top   = get_PT(Mom(1:4,t))
    pT_ATop  = get_PT(Mom(1:4,tbar))
    eta_top  = get_ETA(Mom(1:4,t))
    eta_atop = get_ETA(Mom(1:4,tbar))

!   cosThetaStar 
    pseudoEta_Z = get_PseudoETA(MomZ(1:4))
    cosTheta    = Get_CosAlpha(Mom(1:4,ferm_Z), MomZ(1:4)) 
    gamma_boost = dsqrt(1. + (pT_Z/Minv_Z)**2*cosh(pseudoEta_Z)**2)
    beta_boost  = dsqrt(1. - 1./gamma_boost**2)
    cosThetaStar= (-beta_boost + cosTheta) / (1-beta_boost*cosTheta)

    DphiLL = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,Aferm_Z))  )
    if( DphiLL.gt.Pi ) DphiLL=dabs(2d0*Pi-DphiLL)



    MomBoost(1)   = +MomZ(1)
    MomBoost(2:4) = -MomZ(2:4)
    MomFermZframe(1:4) = Mom(1:4,ferm_Z)
    call boost(MomFermZframe(1:4),MomBoost(1:4), get_MInv(MomZ(1:4)) )
    CosTheta1 = Get_CosAlpha( MomFermZframe(1:4),MomZ(1:4) ) !  = angle between: fermion from Z in the rest frame of the Z and the direction of flight of the Z
!     CosTheta1 = Get_CosAlpha( MomFermZframe(1:4),MomBoost(1:4) ) !  = angle between: fermion from Z in the rest frame of the Z and the direction of flight of the Z

    DphiZt = dabs( Get_PHI(Mom(1:4,Zbos)) - Get_PHI(Mom(1:4,t))  )
    if( DphiZt.gt.Pi ) DphiZt=2d0*Pi-DphiZt

    Dphittbar = dabs( Get_PHI(Mom(1:4,t)) - Get_PHI(Mom(1:4,tbar))  )
    if( Dphittbar.gt.Pi ) Dphittbar=2d0*Pi-Dphittbar




! check cuts
!    if (Minv_Z .lt. 5d0*GeV) then
!       applyPSCut=.true.
!       return
!    endif
    if( dabs(Minv_Z-m_Z).gt.MZ_window ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .OR. pT_jet(4).lt.pT_jet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(3)).gt.eta_jet_cut .OR. abs(eta_jet(4)).gt.eta_jet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    NObsJet=0
    do i=1,NJet
        if( get_PT(MomJet(1:4,i)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,i))).lt.eta_jet_cut ) NObsJet=NObsJet+1
    enddo

    if( pT_lep(1).lt.pT_lep_cut .OR. pT_lep(2).lt.pT_lep_cut .OR. pT_lep(3).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut  .OR. abs(eta_lep(3)).gt.eta_lep_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    leptj = (/ LepP,ferm_Z,Aferm_Z /)
    do i=1,4
    do j=1,3
      if( get_R( MomJet(1:4,i),Mom(1:4,leptj(j)) ) .lt. Rsep_jetlep ) then 
        applyPSCut = .true.
        RETURN
      endif
    enddo
    enddo


! apply lepton-lepton isolation only for ObsSet=58 (CMS 8 TeV analysis)
    if (ObsSet.EQ.58) then
       do i=1,3
          do j=i+1,3
             ii = leptj(i)
             jj = leptj(j)
             if (get_R( Mom(1:4,ii),Mom(1:4,jj)) .lt. Rsep_LepLep) then
                applyPSCut = .true.
                RETURN
             endif
          enddo
       enddo
    endif

   ! this is Ehat for <EHat> in NHisto=12
   MInv_LB= get_MInv(Mom(1:4,inLeft)+Mom(1:4,inRight))


! binning
    NBin(1:2) = 0
    if( 0.5d0*(pT_ATop+pT_Top).lt.800d0*GeV ) then
        NBin(1) = WhichBin(1,0.5d0*(pT_ATop+pT_Top)  )
        NBin(2) = WhichBin(2,0.5d0*(eta_ATop+eta_Top))
    endif
!     NBin(1) = WhichBin(1,pT_Lep(1))
!     NBin(2) = WhichBin(2,pT_Lep(2))
    NBin(3) = WhichBin(3,pT_Lep(3))
    NBin(4) = WhichBin(4,pT_jet(1))
    NBin(5) = WhichBin(5,pT_jet(2))
    NBin(6) = WhichBin(6,pT_jet(3))
    NBin(7) = WhichBin(7,pT_jet(4))
    NBin(8) = WhichBin(8,pT_miss)
    NBin(9) = WhichBin(9,eta_Lep(1))
    NBin(10) = WhichBin(10,eta_Lep(2))
    NBin(11) = WhichBin(11,eta_Lep(3))
    NBin(12) = WhichBin(12,pT_Z)
    NBin(13) = WhichBin(13,eta_Z)
    NBin(14) = WhichBin(14,pT_top)
    NBin(15) = WhichBin(15,eta_top)
    NBin(16) = WhichBin(16,dabs(eta_Top)-dabs(eta_ATop)  )
!     NBin(16) = WhichBin(16,eta_atop)
    NBin(17) = WhichBin(17,DphiLL)
    NBin(18) = WhichBin(18,CosTheta1)
    NBin(19) = WhichBin(19,DphiZt)
    NBin(20) = WhichBin(20,Dphittbar)
    NBin(21) = WhichBin(21,dble(NObsJet))
    NBin(22) = WhichBin(22,DphiLL)
    NBin(23) = WhichBin(23,CosTheta1)
    NBin(24) = WhichBin(24,DphiZt)
    NBin(25) = WhichBin(25,Dphittbar)
    NBin(26) = WhichBin(26,Minv_Z)
    NBin(28) = WhichBin(28,CosThetaStar)


    if( present(PObs) ) then
      PObs(1) = pT_Lep(1)
      PObs(2) = pT_Lep(2)
      PObs(3) = pT_Lep(3)
      PObs(4) = pT_jet(1)
      PObs(5) = pT_jet(2)
      PObs(6) = pT_jet(3)
      PObs(7) = pT_jet(4)
      PObs(8) = pT_miss
      PObs(9) = eta_Lep(1)
      PObs(10) = eta_Lep(2)
      PObs(11) = eta_Lep(3)
      PObs(12) = pT_Z
      PObs(13) = eta_Z
      PObs(14) = pT_top
      PObs(15) = eta_top
      PObs(16) = dabs(eta_Top)-dabs(eta_ATop)
!       PObs(16) = eta_atop
      PObs(17) = DphiLL
      PObs(18) = CosTheta1
      PObs(19) = DphiZt
      PObs(20) = Dphittbar
      PObs(21) = dble(NObsJet)
      PObs(22) = DphiLL
      PObs(23) = CosTheta1
      PObs(24) = DphiZt
      PObs(25) = Dphittbar
      PObs(26) = Minv_Z
      PObs(28) = CosThetaStar
    endif


elseif (ObsSet .eq. 54) then    ! this is for the observed CMS set at 7 TeV

   ! request at least three jets where two are b-jets
   NObsJet_Tree = 3
   ! Ask Markus: how does this enforce 2 bjet requirement???
   if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
      applyPSCut = .true.
      RETURN
    endif

    LepIndex=(/LepP,ferm_Z,Aferm_Z/)
    JetIndex=(/1,2,3,4/)
    
    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    pT_jet(3) = get_PT(MomJet(1:4,3))
    pT_jet(4) = get_PT(MomJet(1:4,4))
    HT_jet=pT_jet(1)+pT_jet(2)+pT_jet(3)+pT_jet(4)
    eta_jet(1) = get_ETA(MomJet(1:4,1))
    eta_jet(2) = get_ETA(MomJet(1:4,2))
    eta_jet(3) = get_ETA(MomJet(1:4,3))
    eta_jet(4) = get_ETA(MomJet(1:4,4))


    pT_Lep(1)  = get_PT(Mom(1:4,LepP))
    eta_Lep(1) = get_ETA(Mom(1:4,LepP))
    pT_Lep(2)  = get_PT(Mom(1:4,ferm_Z))
    eta_Lep(2) = get_ETA(Mom(1:4,ferm_Z))
    pT_Lep(3)  = get_PT(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,Aferm_Z))

    pT_ll=get_PT( (Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z)) )

    pT_miss = get_PT(Mom(1:4,nu))
    
    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
       applyPSCut = .true.
       RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
       RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .OR. pT_jet(4).lt.pT_jet_cut ) then
       applyPSCut = .true.
       RETURN
    endif

    if( abs(eta_jet(3)).gt.eta_jet_cut .OR. abs(eta_jet(4)).gt.eta_jet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut  .OR. abs(eta_lep(3)).gt.eta_lep_cut ) then
       applyPSCut = .true.
       RETURN
    endif

    if (pT_lep(2).lt.pT_lepZ_cut .OR. pT_lep(3).lt.pT_lepZ_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lep(1).lt.pT_lept_cut ) then 
        applyPSCut = .true.
        RETURN
    endif

    if (pT_ll .le. pT_ll_cut) then
       applyPSCut = .true.
       RETURN
    endif


    if (HT_jet .le. HT_jet_cut) then
       applyPSCut = .true.
       RETURN
    endif

    WithinCone=0d0
    do iLept=1,3
       do jLept=1,3
          if (jLept .eq. iLept) cycle
          RLept=Get_R(Mom(1:4,LepIndex(iLept)),Mom(1:4,LepIndex(jLept)))
          if (Rlept .le. Rsep_jetlep ) then
             WithinCone(iLept)=WithinCone(iLept)+get_PT(Mom(1:4,LepIndex(jLept))) &
                  &+Mom(1,LepIndex(jLept))
          endif
       enddo
       do jJet=1,NJet
          RLept=Get_R(Mom(1:4,LepIndex(iLept)),MomJet(1:4,Jetindex(jJet)))
          if (Rlept .le. Rsep_jetlep ) then
             WithinCone(iLept)=WithinCone(iLept)+get_PT(MomJet(1:4,JetIndex(jJet))) &
                  &+MomJet(1,JetIndex(jJet))
          endif
       enddo
       if ( WithinCone(iLept) .ge. Frac_sep_jetlep*get_PT(Mom(1:4,LepIndex(iLept))) ) then
          applyPSCut = .true.
          RETURN
       endif
    enddo



elseif( ObsSet.eq.57 ) then! set of observables for ttb+Z ( di-lept. ttbar decays and di-lept. Z decay ) at TEvatron

! request at least four jets where two are b-jets
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif


    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    pT_jet(3) = get_PT(MomJet(1:4,3))
    pT_jet(4) = get_PT(MomJet(1:4,4))
    eta_jet(1) = get_ETA(MomJet(1:4,1))
    eta_jet(2) = get_ETA(MomJet(1:4,2))
    eta_jet(3) = get_ETA(MomJet(1:4,3))
    eta_jet(4) = get_ETA(MomJet(1:4,4))


    pT_Lep(1)  = get_PT(Mom(1:4,LepP))
    eta_Lep(1) = get_ETA(Mom(1:4,LepP))

    pT_Lep(2)  = get_PT(Mom(1:4,ferm_Z))
    eta_Lep(2) = get_ETA(Mom(1:4,ferm_Z))

    pT_Lep(3)  = get_PT(Mom(1:4,Aferm_Z))
    eta_Lep(3) = get_ETA(Mom(1:4,Aferm_Z))

    pT_miss = get_PT(Mom(1:4,nu))

    MomZ(1:4) = Mom(1:4,ferm_Z)+Mom(1:4,Aferm_Z)
    pT_Z  = get_PT(MomZ(1:4))
    eta_Z = get_eta(MomZ(1:4))

    pT_top = get_PT(Mom(1:4,t))
    eta_top = get_ETA(Mom(1:4,t))
    eta_atop = get_ETA(Mom(1:4,tbar))



    DphiLL = dabs( Get_PHI(Mom(1:4,ferm_Z)) - Get_PHI(Mom(1:4,Aferm_Z))  )
    if( DphiLL.gt.Pi ) DphiLL=dabs(2d0*Pi-DphiLL)



    MomBoost(1)   = +MomZ(1)
    MomBoost(2:4) = -MomZ(2:4)
    MomFermZframe(1:4) = Mom(1:4,ferm_Z)
    call boost(MomFermZframe(1:4),MomBoost(1:4), get_MInv(MomZ(1:4)) )
    CosTheta1 = Get_CosAlpha( MomFermZframe(1:4),MomZ(1:4) ) !  = angle between: fermion from Z in the rest frame of the Z and the direction of flight of the Z
!     CosTheta1 = Get_CosAlpha( MomFermZframe(1:4),MomBoost(1:4) ) !  = angle between: fermion from Z in the rest frame of the Z and the direction of flight of the Z

    DphiZt = dabs( Get_PHI(Mom(1:4,Zbos)) - Get_PHI(Mom(1:4,t))  )
    if( DphiZt.gt.Pi ) DphiZt=2d0*Pi-DphiZt

    Dphittbar = dabs( Get_PHI(Mom(1:4,t)) - Get_PHI(Mom(1:4,tbar))  )
    if( Dphittbar.gt.Pi ) Dphittbar=2d0*Pi-Dphittbar




! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .OR. pT_jet(4).lt.pT_jet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_jet(3)).gt.eta_jet_cut .OR. abs(eta_jet(4)).gt.eta_jet_cut) then
       applyPSCut = .true.
        RETURN
    endif

    NObsJet=0
    do i=1,NJet
        if( get_PT(MomJet(1:4,i)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,i))).lt.eta_jet_cut ) NObsJet=NObsJet+1
    enddo

    if( pT_lep(1).lt.pT_lep_cut .OR. pT_lep(2).lt.pT_lep_cut .OR. pT_lep(3).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut .OR. abs(eta_lep(2)).gt.eta_lep_cut  .OR. abs(eta_lep(3)).gt.eta_lep_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    leptj = (/ LepP,ferm_Z,Aferm_Z /)
    do i=1,4
    do j=1,3
      if( get_R( MomJet(1:4,i),Mom(1:4,leptj(j)) ) .lt. Rsep_jetlep ) then 
        applyPSCut = .true.
        RETURN
      endif
    enddo
    enddo



! binning
    NBin(1) = WhichBin(1,eta_top)
    NBin(2) = WhichBin(2,eta_atop)
    if( eta_top.gt.0d0 ) NBin(3) = WhichBin(3,eta_top)
    if( eta_top.le.0d0 ) NBin(4) = WhichBin(4,eta_top)


    if( present(PObs) ) then
      PObs(1) = pT_Lep(1)
      PObs(2) = pT_Lep(2)
      PObs(3) = pT_Lep(3)
      PObs(4) = pT_jet(1)
      PObs(5) = pT_jet(2)
      PObs(6) = pT_jet(3)
      PObs(7) = pT_jet(4)
      PObs(8) = pT_miss
      PObs(9) = eta_Lep(1)
      PObs(10) = eta_Lep(2)
      PObs(11) = eta_Lep(3)
      PObs(12) = pT_Z
      PObs(13) = eta_Z
      PObs(14) = pT_top
      PObs(15) = eta_top
      PObs(16) = eta_atop
      PObs(17) = DphiLL
      PObs(18) = CosTheta1
      PObs(19) = DphiZt
      PObs(20) = Dphittbar
      PObs(21) = dble(NObsJet)
      PObs(22) = DphiLL
      PObs(23) = CosTheta1
      PObs(24) = DphiZt
      PObs(25) = Dphittbar
    endif

!-------------------------------------------------------
else
  print *, "ObsSet not implemented TTBZ",ObsSet
  stop
endif


return
END SUBROUTINE






SUBROUTINE Kinematics_TTBARH(NPlus1PS,Mom,MomOrder,applyPSCut,NBin,PObs)
use ModMisc
use ModParameters
#if _UseJHUGenMELA==1
use ModTTBH
#endif
implicit none
integer :: NumHadr,NPlus1PS,MomOrder(1:14)
real(8) :: Mom(1:4,1:14),zeros(1:14)
real(8) :: MomJet(1:4,1:7) 
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4)
logical :: applyPSCut
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree,i,j
real(8),optional :: PObs(:)
real(8) :: pT_lep(4),ET_miss,PT_miss,pT_ATop,pT_Top,pT_Higgs,HT,ET_bjet,eta_Higgs
real(8) :: eta_ATop,eta_Top,eta_lep(1:4),m_ttbar,delta_eta_T,delta_eta_B,cos_thetaLL
real(8) :: pT_jet(1:7),eta_jet(1:7),MELA_rescale
integer :: tbar,t,Hig,inLeft,inRight,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,L,N,HDK1,HDK2
real(8) :: pT_ll,HT_jet,RLept,CosTheta1,CosTheta2,CosThetaStar,delta_eta_L,cos_thetaBB
integer :: iLept,jLept,jJet
real(8) :: MomMELA(1:4,1:13),MatElSq_H0,MatElSq_H1,D_0minus,MomAux1(1:4),MomAux2(1:4),M_aux,MomRest(1:4)
logical,save :: FirstTime=.true.


applyPSCut = .false.


! momentum ordering
  tbar    = MomOrder(1)
  t       = MomOrder(2)
  Hig     = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)
  HDK1    = MomOrder(13)
  HDK2    = MomOrder(14)
  
  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu

!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,tbar) - Mom(1:4,t) - Mom(1:4,Hig)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,Hig) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZ(): ",NPlus1PS,zeros(1:4)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
   
   if(TopDecays.ne.0 .and. Correction.ne.5) then
      zeros(1:4) = Mom(1:4,tbar) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar)
      zeros(5:8) = Mom(1:4,t) - Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
   endif
   if( any(abs(zeros(1:8)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZ(): ",NPlus1PS,zeros(1:8)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
      
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   zeros(3) = 0d0
   zeros(4) =  Mom(1:4,bbar).dot.Mom(1:4,bbar)
   zeros(5) =  Mom(1:4,lepM).dot.Mom(1:4,lepM)
   zeros(6) =  Mom(1:4,nubar).dot.Mom(1:4,nubar)
   zeros(7) =  Mom(1:4,b).dot.Mom(1:4,b)
   zeros(8) =  Mom(1:4,lepP).dot.Mom(1:4,lepP)
   zeros(9) =  Mom(1:4,nu).dot.Mom(1:4,nu)
!    zeros(10)=  ((Mom(1:4,lepM)+Mom(1:4,nubar)).dot.(Mom(1:4,lepM)+Mom(1:4,nubar))) - m_W**2
!    zeros(11)=  ((Mom(1:4,lepP)+Mom(1:4,nu)).dot.(Mom(1:4,lepP)+Mom(1:4,nu))) - m_W**2
   if( NPlus1PS.eq.1 ) zeros(12)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   

   if( TopDecays.eq.0 .and. any(abs(zeros(1:3)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:3)
      print *, Mom(1:4,1:3)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:12)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZ(): ",zeros(1:12)
      print *, "momenta dump:"
      print *, Mom(1:4,1:12)
   endif
!DEC$ ENDIF




NBin(1:NumHistograms) = 0
MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)
MomMiss(1:4) = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
    if(NPlus1PS.eq.0) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 .OR. TOPDECAYS.eq.-1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    L = LepM  ! ambiguously defined for di-leptons
    N = nubar
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   L = LepM
   N = nubar
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   L = LepP
   N = nu
   if(NPlus1PS.eq.0) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF





!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    if( ObsSet.ne.83 ) then
        call pT_order(2,MomJet(1:4,1:2))
        call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))  ! remove this for Obsset=83 due to MELA momenta assignment (should be unordered)
    endif

! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))





!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.81) then! ttb+H production without top decays at LHC

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))

    pT_Higgs= get_PT(Mom(1:4,Hig))

! binning
    NBin(1) = WhichBin(1,pT_Top)
    NBin(2) = WhichBin(2,pT_Higgs)
    if( present(PObs) ) then
      PObs(1) = pT_Top
      PObs(2) = pT_Higgs
   endif





!------------------------ cuts and binning --------------------------------
elseif( ObsSet.eq.82) then! ttb+H production with di-leptonic tops


    NObsJet_Tree = 2
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif


!     if( njet.eq.3 .and. get_PT(MomJet(1:4,3)).gt.10d0*GeV  ) then! reject events with resolved jet
!         applyPSCut = .true.
!         RETURN        
!     endif



    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))
    eta_top = get_ETA(Mom(1:4,t))
    eta_Atop = get_ETA(Mom(1:4,tbar))

    pT_Higgs= get_PT(Mom(1:4,Hig))
    eta_Higgs = get_ETA(Mom(1:4,Hig))
    m_ttbar = get_MInv(Mom(1:4,t)+Mom(1:4,tbar))

! binning
    NBin(1) = WhichBin(1,pT_Top)
    NBin(2) = WhichBin(2,pT_ATop)
    NBin(3) = WhichBin(3,eta_top)
    NBin(4) = WhichBin(4,eta_Atop)
    NBin(5) = WhichBin(5,m_ttbar)
    NBin(6) = WhichBin(6,pT_Higgs)
    NBin(7) = WhichBin(7,eta_Higgs)

! for the multidimensional histogram, use the above bins
    NBin(8) = Convert2DHist(NBin(1),NBin(2),Histo2D(1)%NBins2,Histo2D(1)%NBins3) 
    NBin(9) = Convert2DHist(NBin(1),NBin(3),Histo2D(2)%NBins2,Histo2D(2)%NBins3) 
    NBin(10) = Convert2DHist(NBin(1),NBin(4),Histo2D(3)%NBins2,Histo2D(3)%NBins3) 
    NBin(11) = Convert2DHist(NBin(1),NBin(6),Histo2D(4)%NBins2,Histo2D(4)%NBins3) 
    NBin(12) = Convert2DHist(NBin(1),NBin(7),Histo2D(5)%NBins2,Histo2D(5)%NBins3) 
    NBin(13) = Convert2DHist(NBin(3),NBin(5),Histo2D(6)%NBins2,Histo2D(6)%NBins3) 
    NBin(14) = Convert2DHist(NBin(3),NBin(6),Histo2D(7)%NBins2,Histo2D(7)%NBins3) 
    NBin(15) = Convert2DHist(NBin(3),NBin(7),Histo2D(8)%NBins2,Histo2D(8)%NBins3) 
    NBin(16) = Convert2DHist(NBin(5),NBin(6),Histo2D(9)%NBins2,Histo2D(9)%NBins3) 
    NBin(17) = Convert2DHist(NBin(5),NBin(7),Histo2D(10)%NBins2,Histo2D(10)%NBins3) 
    NBin(18) = Convert2DHist(NBin(6),NBin(7),Histo2D(11)%NBins2,Histo2D(11)%NBins3) 

    if( present(PObs) ) then
      PObs(1) = pT_Top
      PObs(2) = pT_ATop
      PObs(3) = eta_top
      PObs(4) = eta_Atop
      PObs(5) = m_ttbar
      PObs(6) = pT_Higgs
      PObs(7) = eta_Higgs
   endif

!------------------------ cuts and binning --------------------------------
elseif( ObsSet.eq.83) then! ttb+H production with semi-leptonic tops

      if( .not.( any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
          applyPSCut = .true.
          RETURN
     endif

    pT_ATop = get_PT(Mom(1:4,tbar))
    pT_Top  = get_PT(Mom(1:4,t))
    eta_top = get_ETA(Mom(1:4,t))
    eta_Atop = get_ETA(Mom(1:4,tbar))
    pT_Higgs= get_PT(Mom(1:4,Hig))
    eta_Higgs = get_ETA(Mom(1:4,Hig))

    do j = 1,NJet
      pT_jet(j)  = get_PT(MomJet(1:4,j))
      eta_jet(j) = get_ETA(MomJet(1:4,j))
    enddo
    
    
    
    
    
!   check if b-jets pass cuts
    if(  pT_jet(1).lt.pT_bjet_cut .or. abs(eta_jet(1)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_jet(2).lt.pT_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( pT_jet(k).gt.pT_jet_cut .and. abs(eta_jet(k)).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif
    
    
    pT_lep(1) = get_PT(Mom(1:4,L))
    eta_lep(1) = get_ETA(Mom(1:4,L))
    pT_miss   = get_PT(Mom(1:4,N))
    delta_eta_T = abs( eta_top - eta_Atop )
    delta_eta_B = abs( eta_jet(1) - eta_jet(2) )
    delta_eta_L = 0d0
    
    cos_thetaLL = 0d0
    cos_thetaBB = 0d0


! construct cos(theta1)=cos(theta_tt):
!       MomRest(1:4)  = Mom(1:4,Hig)
!       MomBoost(1)   = +MomRest(1)
!       MomBoost(2:4) = -MomRest(2:4)
!       MomAux(1:4) = Mom(1:4,t)+Mom(1:4,tbar)
!       M_Aux = get_MInv(MomAux(1:4))
!       call boost(MomAux(1:4),MomBoost(1:4),M_Aux)
      CosTheta1 = 0d0




! construct cos(theta*)=Cos(theta_H):    scattering angle of Higgs in ttb+H rest frame
      MomRest(1:4)= Mom(1:4,tbar) + Mom(1:4,t) + Mom(1:4,Hig)
      MomBoost(1)   = +MomRest(1)
      MomBoost(2:4) = -MomRest(2:4)
      
      MomAux1(1:4) = Mom(1:4,Hig)
      MomAux2(1:4) = (/0d0,0d0,0d0,1d0/) ! z-axis
      
      M_Aux = dabs( get_MInv(MomRest(1:4)) + 1d-12 )
      call boost(MomAux1(1:4),MomBoost(1:4),M_Aux)
      call boost(MomAux2(1:4),MomBoost(1:4),M_Aux)
      
      CosThetaStar = Get_CosAlpha( MomAux1(1:4),MomAux2(1:4) )

      
      
! construct cos(theta2)=cos(theta_t):
      MomRest(1:4)  = Mom(1:4,t)+Mom(1:4,tbar)
      MomBoost(1)   = +MomRest(1)
      MomBoost(2:4) = -MomRest(2:4)
      
!       MomAux1(1:4) = Mom(1:4,Hig)
      MomAux1(1:4) = Mom(1:4,t)
!       MomAux2(1:4) = (/0d0,0d0,0d0,1d0/) ! z-axis
      MomAux2(1:4) = Mom(1:4,Hig)


      M_Aux = dabs( get_MInv(MomRest(1:4)) + 1d-12 )
      call boost(MomAux1(1:4),MomBoost(1:4),M_Aux)
      call boost(MomAux2(1:4),MomBoost(1:4),M_Aux)
      
      CosTheta2 = -Get_CosAlpha( MomAux1(1:4),MomAux2(1:4) )



! check cuts
!     if( pT_jet(1).lt.pT_bjet_cut .OR. pT_jet(2).lt.pT_bjet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! 
!     if( abs(eta_jet(1)).gt.eta_bjet_cut .OR. abs(eta_jet(2)).gt.eta_bjet_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif
! 
!     if( pT_jet(3).lt.pT_jet_cut .OR. pT_jet(4).lt.pT_jet_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif
! 
!     if( abs(eta_jet(3)).gt.eta_jet_cut .OR. abs(eta_jet(4)).gt.eta_jet_cut) then
!         applyPSCut = .true.
!         RETURN
!     endif

    if( pT_lep(1).lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lep(1)).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif


     
#if _UseJHUGenMELA==1
   
    if( FirstTime ) then
      call NNPDFDriver("./PDFS/NNPDF30_lo_as_0130.LHgrid")
!       call NNPDFDriver("./PDFS/NNPDF30_nlo_as_011.LHgrid")
      call NNinitPDF(0)
      call InitProcess_TTBH(m_H,m_top)
      FirstTime = .false.
    endif
    MomMELA(1:4,1) = -Collider_Energy/2d0 *  (/+1d0, 0d0, 0d0, 1d0 /)
    MomMELA(1:4,2) = -Collider_Energy/2d0 *  (/+1d0, 0d0, 0d0,-1d0 /)  
    MomMELA(1:4,3) = Mom(1:4,Hig)
    MomMELA(1:4,4) = Mom(1:4,tbar)
    MomMELA(1:4,5) = Mom(1:4,t)
    MomMELA(1:4,6) = MomJet(1:4,1)! bbar
    MomMELA(1:4,9) = MomJet(1:4,2)! b
    MomMELA(1:4,7) = MomJet(1:4,4)! this ordering co-incides with the exact LO choice
    MomMELA(1:4,8) = MomJet(1:4,3)
    MomMELA(1:4,10)= Mom(1:4,lepP)
    MomMELA(1:4,11)= Mom(1:4,nu)
    MomMELA(1:4,12:13) = 0d0
    
    call EvalXSec_PP_TTBH(MomMELA(1:4,1:13),(/(1d0,0d0),(0d0,0d0)/),TopDecays,2,MatElSq_H0)
    call EvalXSec_PP_TTBH(MomMELA(1:4,1:13),(/(0d0,0d0),(1d0,0d0)/),TopDecays,2,MatElSq_H1)
    
    MELA_rescale = 2.537649d0
    
    D_0minus = MatElSq_H0/(MatElSq_H0 + MELA_rescale*MatElSq_H1 )

#endif    
    

! binning
    NBin(1) = WhichBin(1,pT_Top)
    NBin(2) = WhichBin(2,eta_Top)
    NBin(3) = WhichBin(3,pT_Higgs)
    NBin(4) = WhichBin(4,eta_Higgs)
    NBin(5) = WhichBin(5,delta_eta_T)
    NBin(6) = WhichBin(6,delta_eta_L)
    NBin(7) = WhichBin(7,delta_eta_B)
    NBin(8) = WhichBin(8,cos_thetaLL)
    NBin(9) = WhichBin(9,cos_thetaBB)
    NBin(10)= WhichBin(10,CosThetaStar)
    NBin(11)= WhichBin(11,CosTheta2)
    NBin(12)= WhichBin(12,D_0minus)
    if( NObsJet.eq.5 ) then 
!        NBin(13)= WhichBin(13,get_PT(MomJet(1:4,5)))! if 5 jets then plot softest
        NBin(13)= WhichBin(13,get_PT(MomJet(1:4,1)))! if 5 jets then plot hardest
    else
        NBin(13)= 0
    endif
    NBin(14)= WhichBin(14,dble(NObsJet))

    if( present(PObs) ) then
      PObs(1) = pT_Top
      PObs(2) = eta_Top
      PObs(3) = pT_Higgs
      PObs(4) = eta_Higgs
      PObs(5) = delta_eta_T
      PObs(6) = delta_eta_L
      PObs(7) = delta_eta_B
      PObs(8) = cos_thetaLL
      PObs(9) = cos_thetaBB
      PObs(10)= CosThetaStar
      PObs(11)= CosTheta2
      PObs(12)= D_0minus
      if( NObsJet.eq.5 ) then 
!        PObs(13)= get_PT(MomJet(1:4,5))! if 5 jets then plot softest
        PObs(13)= get_PT(MomJet(1:4,1))! if 5 jets then plot hardest
      else
        PObs(13)= 0d0
      endif
      PObs(14)= dble(NObsJet)
   endif



!-------------------------------------------------------
else
  print *, "ObsSet not implemented TTBH",ObsSet
  stop
endif


return
END SUBROUTINE










SUBROUTINE Kinematics_TH(NPlus1PS,Mom,MomOrder,applyPSCut,NBin,PObs)
  use ModMisc
  use ModParameters
  implicit none
  integer :: NPlus1PS,MomOrder(1:11)
  real(8) :: Mom(1:4,1:11),zeros(1:9)
  integer :: t,ljet,Hig,inLeft,inRight,realp,b,lepP,nu,Hdk1,Hdk2,k
  real(8) :: MomJet(1:4,1:7) !,MomJet_CHECK(1:4,1:7)                                                                                                                    
  real(8) :: MomHadr(1:4,0:8),MomMiss(1:4)
  real(8) :: pT_top,pT_higgs,pT_j,eta_j,eta_top,eta_Higgs,y_top,y_Higgs,y_j,eta_Hj,Dy_jt,Dy_Ht,Dy_Hj,Deta_Hj,Deta_Ht,eta_Ht,y_Ht,y_tj,m_Ht,eta_tj,y_Hj,m_tj,m_Hj
  real(8) :: pT_miss,pT_ep,pT_b,eta_ep,eta_b,y_bep,y_jep,eta_W,eta_Hb,eta_Wj,eta_Wb,Deta_jb,Deta_Hb,Deta_Wj,Deta_Wb,Deta_tj,m_Hb,m_Wb,m_Wj,m_Hbep,m_Hjep,DR_bep,mtrans_W,HTtot,DR_jep,Deta_Hep, costheta_j,costheta_l,costheta_lj,DR_bj
  real(8),optional :: PObs(:)
  logical :: applyPSCut
  integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NumHadr


  applyPSCut = .false.

! momentum ordering
  t      = MomOrder(1)
  ljet   = MomOrder(2) 
  Hig     = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  b    = MomOrder(7)
  lepP    = MomOrder(8)
  nu   = MomOrder(9)
  Hdk1    = MomOrder(10)
  Hdk2    = MomOrder(11)

!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ljet) - Mom(1:4,t) - Mom(1:4,Hig)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ljet) - Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu) - Mom(1:4,Hdk1)-Mom(1:4,Hdk2)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TH(): ",NPlus1PS,zeros(1:4)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
   zeros(1) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_Top**2
   zeros(2) = (Mom(1:4,Hig).dot.Mom(1:4,Hig)) - M_H**2
   zeros(3) = (Mom(1:4,ljet).dot.Mom(1:4,ljet))
   zeros(4) = (Mom(1:4,b).dot.Mom(1:4,b))
   zeros(5) = (Mom(1:4,lepP).dot.Mom(1:4,lepP))
   zeros(6) = (Mom(1:4,nu).dot.Mom(1:4,nu))
   zeros(7) = (Mom(1:4,Hdk1).dot.Mom(1:4,Hdk1))
   zeros(8) = (Mom(1:4,Hdk2).dot.Mom(1:4,Hdk2))


   if( NPlus1PS.eq.1 ) zeros(9)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:3)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TH(): ",zeros(1:3)
      print *, Mom(1:4,1:3)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:9)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TH(): ",zeros(1:9)
      print *, "momenta dump:"
      print *, Mom(1:4,1:9)
   endif
   print *, "CHECK MOMENTA PASSED!"
!DEC$ ENDIF  


NBin(1:NumHistograms) = 0
MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)
MomMiss(1:4) = 0d0



IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays
   MomHadr(1:4,1) = Mom(1:4,ljet)
    if(NPlus1PS.eq.0) then
        NumHadr = 1
    else
        NumHadr = 2
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 .OR. TOPDECAYS.eq.-1 ) THEN  ! di-leptonic decay                                                                                                         
!    MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,1) = Mom(1:4,b)     ! Bot                                                                                                                                       
   MomHadr(1:4,2) = Mom(1:4,ljet)
 
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF




!---------------------- kT jet algorithm ---------------------------------                                                                                                     
! important: b-quarks need to be at the first two positions of MomHadr(:,:)                                                                                                    
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.                                                                                        
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order                                                                                                 

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets   
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))                                                                                                                           
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA                               
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))                                                              

    
    if( ObsSet.eq.91) then! t+H production without top decays at LHC                                                                                                      


       pT_Top  = get_PT(Mom(1:4,t))
       pT_Higgs= get_PT(Mom(1:4,Hig))
       pT_j= get_PT(Mom(1:4,ljet))

       y_top=get_eta(Mom(1:4,t))
       y_Higgs=get_eta(Mom(1:4,Hig))
       y_j=get_eta(Mom(1:4,ljet))

       eta_top=get_pseudoeta(Mom(1:4,t))
       eta_Higgs=get_pseudoeta(Mom(1:4,Hig))
       eta_j=get_pseudoeta(Mom(1:4,ljet))

       eta_tj=get_pseudoeta(Mom(1:4,t)+Mom(1:4,ljet))
       eta_Hj=get_pseudoeta(Mom(1:4,Hig)+Mom(1:4,ljet))
       eta_Ht=get_pseudoeta(Mom(1:4,Hig)+Mom(1:4,t))

       y_tj=get_eta(Mom(1:4,t)+Mom(1:4,ljet))       
       y_Hj=get_eta(Mom(1:4,Hig)+Mom(1:4,ljet))
       y_Ht=get_eta(Mom(1:4,Hig)+Mom(1:4,t))
       
       m_Ht=get_MInv(Mom(1:4,Hig)+Mom(1:4,t))
       m_Hj=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet))
       m_tj=get_MInv(Mom(1:4,t)+Mom(1:4,ljet))
       
       Dy_Hj=y_Higgs-y_j
       Dy_Ht=y_Higgs-y_top
       Dy_jt=y_j-y_top

       Deta_Hj=eta_Higgs-eta_j
       Deta_Ht=eta_Higgs-eta_top

       ! binning                                                                                                                                             
       NBin(1) = WhichBin(1,pT_Top)
       NBin(2) = WhichBin(2,eta_Top)
       NBin(3) = WhichBin(3,y_Top)
       NBin(4) = WhichBin(4,pT_j)
       NBin(5) = WhichBin(5,y_j)
       NBin(6) = WhichBin(6,pT_Higgs)
       NBin(7) = WhichBin(7,m_tj)
       NBin(8) = WhichBin(8,y_Higgs)
       NBin(9) = WhichBin(9,eta_Higgs)
       NBin(10) = WhichBin(10,y_tj)
       NBin(11) = WhichBin(11,eta_tj)
       NBin(12) = WhichBin(12,Dy_jt)
       NBin(13) = WhichBin(13,m_Ht)
       NBin(14) = WhichBin(14,m_Hj)
       NBin(15) = WhichBin(15,y_Ht)
       NBin(16) = WhichBin(16,eta_Ht)
       NBin(17) = WhichBin(17,y_Hj)
       NBin(18) = WhichBin(18,eta_Hj)
       NBin(19) = WhichBin(19,Dy_Ht)
       NBin(20) = WhichBin(20,Deta_Ht)
       NBin(21) = WhichBin(21,Dy_Hj)
       NBin(22) = WhichBin(22,Deta_Hj)

       if( present(PObs) ) then
          PObs(1) = pT_Top
          PObs(2) = eta_Top
          PObs(3) = y_Top
          PObs(4) = pt_j
          PObs(5) = y_j
          PObs(6) = pT_Higgs
          PObs(7) = m_tj
          PObs(8) = y_Higgs
          PObs(9) = eta_Higgs
          PObs(10) = y_tj
          PObs(11) = eta_tj
          POBs(12) = Dy_jt
          PObs(13) = m_Ht
          PObs(14) = m_Hj
          PObs(15) = y_Ht
          PObs(16) = eta_Ht
          PObs(17) = y_Hj
          PObs(18) = eta_Hj
          POBs(19) = Dy_Ht
          POBs(20) = Deta_Ht
          POBs(21) = Dy_Hj
          POBs(22) = Deta_Hj

       endif
       !-------------------------------------------------------             
    elseif (ObsSet .eq. 93 .or. ObsSet .eq. 95) then      ! leptonic top decay
       pT_ep = get_PT(Mom(1:4,lepP))
       pT_miss = get_PT(Mom(1:4,nu))
       pT_Higgs= get_PT(Mom(1:4,Hig))
       pT_j = get_PT(Mom(1:4,ljet))
       pT_b = get_PT(Mom(1:4,b))
       pT_Top= get_PT(Mom(1:4,t))
       
       eta_ep=get_pseudoeta(Mom(1:4,lepP))
       eta_b=get_pseudoeta(Mom(1:4,b))
       eta_top=get_pseudoeta(Mom(1:4,t))
       eta_Higgs=get_pseudoeta(Mom(1:4,Hig))
       eta_j=get_pseudoeta(Mom(1:4,ljet))
       eta_W=get_pseudoeta(Mom(1:4,lepP)+Mom(1:4,nu))
       eta_Hj=get_pseudoeta(Mom(1:4,ljet)+Mom(1:4,Hig))
       eta_Hb=get_pseudoeta(Mom(1:4,b)+Mom(1:4,Hig))
       eta_Wj=get_pseudoeta(Mom(1:4,lepP)+Mom(1:4,nu)+Mom(1:4,ljet))
       eta_Wb=get_pseudoeta(Mom(1:4,lepP)+Mom(1:4,nu)+Mom(1:4,b))

       DR_bj=get_R(Mom(1:4,b),Mom(1:4,ljet))
       DR_bep=get_R(Mom(1:4,b),Mom(1:4,lepP))
       DR_jep=get_R(Mom(1:4,ljet),Mom(1:4,lepP))


       IF (ObsSet .eq. 95) THEN
          if (pT_ep .le. pT_lep_cut .or. pT_b .le. pT_bjet_cut .or. pT_j .le. pT_jet_cut) then
             applyPSCut = .true.
             RETURN
          endif

          if (abs(eta_ep) .ge. eta_lep_cut .or. abs(eta_j) .ge. eta_jet_cut .or. abs(eta_b) .ge. eta_bjet_cut) then
             applyPSCut = .true.
             RETURN
          endif

          if (DR_bj .le. Rsep_jet .or. DR_bep .le. Rsep_lepjet .or. DR_jep .le. Rsep_lepjet) then
             applyPSCut = .true.
             RETURN
          endif
       ENDIF


 !      y_top=get_eta(Mom(1:4,t))
       y_Higgs=get_eta(Mom(1:4,Hig))
       y_bep=get_eta(Mom(1:4,b)+Mom(1:4,lepP))
       y_jep=get_eta(Mom(1:4,ljet)+Mom(1:4,lepP))
       y_Hj=get_eta(Mom(1:4,ljet)+Mom(1:4,Hig))
       y_j=get_eta(Mom(1:4,ljet))
       
       costheta_j = get_CosTheta(Mom(1:4,ljet))
       costheta_l = get_CosTheta(Mom(1:4,lepP))
!       costheta_lj=costheta_l-costheta_j
       costheta_lj=get_cosalpha(Mom(1:4,lepP),Mom(1:4,ljet))

       Deta_jb=eta_j-eta_b
       Deta_Hj=eta_Higgs-eta_j
       Deta_Hb=eta_Higgs-eta_b
       Deta_Wj=eta_W-eta_j
       Deta_Wb=eta_W-eta_b
       Deta_Ht=eta_Higgs-eta_top
       Deta_Hep=eta_Higgs-eta_ep
       Deta_tj=eta_top-eta_j
       Dy_Hj=y_Higgs-y_j

       m_Hb=get_MInv(Mom(1:4,Hig)+Mom(1:4,b))
       m_Hj=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet))
       m_Wb=get_MInv(Mom(1:4,lepP)+Mom(1:4,nu)+Mom(1:4,b))
       m_Wj=get_MInv(Mom(1:4,lepP)+Mom(1:4,nu)+Mom(1:4,ljet))
       m_Hbep=get_MInv(Mom(1:4,Hig)+Mom(1:4,b)+Mom(1:4,lepP))
       m_Hjep=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet)+Mom(1:4,lepP))
       
       mtrans_W=get_MT(Mom(1:4,lepP),Mom(1:4,nu))

      
       NBin(1) = WhichBin(1,pT_ep)
       NBin(2) = WhichBin(2,eta_ep)
       NBin(3) = WhichBin(3,pt_miss)
       NBin(4) = WhichBin(4,pt_j)
       NBin(5) = WhichBin(5,eta_j)
       NBin(6) = WhichBin(6,pt_b)
       NBin(7) = WhichBin(7,eta_b)
       NBin(8) = WhichBin(8,pt_Higgs)
       NBin(9) = WhichBin(9,mtrans_W)
       NBin(10) = WhichBin(10,y_Higgs)
       NBin(11) = WhichBin(11,y_bep)
       NBin(12) = WhichBin(12,y_jep)
       NBin(13) = WhichBin(13,eta_Higgs)
       NBin(14) = WhichBin(14,eta_W)
       NBin(15) = WhichBin(15,Deta_jb)
       NBin(16) = WhichBin(16,DR_bep)
       NBin(17) = WhichBin(17,DR_jep)
       NBin(18) = WhichBin(18,m_Hb)
       NBin(19) = WhichBin(19,m_Wb)
       NBin(20) = WhichBin(20,m_Wj)
       NBin(21) = WhichBin(21,m_Hj)
       NBin(22) = WhichBin(22,y_Hj)
       NBin(23) = WhichBin(23,eta_Hj)
       NBin(24) = WhichBin(24,eta_Hb)
       NBin(25) = WhichBin(25,eta_Wb)
       NBin(26) = WhichBin(26,eta_Wj)
       NBin(27) = WhichBin(27,Dy_Hj)
       NBin(28) = WhichBin(28,Deta_Hj)
       NBin(29) = WhichBin(29,Deta_Hb)
       NBin(30) = WhichBin(30,Deta_Wj)
       NBin(31) = WhichBin(31,Deta_wb)
       NBin(32) = WhichBin(32,Deta_Hep)
       NBin(33) = WhichBin(33,m_Hbep)
       NBin(34) = WhichBin(34,m_Hjep)
       NBin(35) = WhichBin(35,Deta_tj)
       NBin(36) = WhichBin(36,Deta_Ht)
       NBin(37) = WhichBin(37,costheta_lj)


      if( present(PObs) ) then


       Pobs(1) = pT_ep
       PObs(2) = eta_ep
       PObs(3) = pt_miss
       PObs(4) = pt_j
       PObs(5) = eta_j
       PObs(6) = pt_b
       PObs(7) = eta_b
       PObs(8) = pt_Higgs
       PObs(9) = mtrans_W
       PObs(10) = y_Higgs
       PObs(11) = y_bep
       PObs(12) = y_jep
       PObs(13) = eta_Higgs
       PObs(14) = eta_W
       PObs(15) = Deta_jb
       PObs(16) = DR_bep
       PObs(17) = DR_jep
       PObs(18) = m_Hb
       PObs(19) = m_Wb
       PObs(20) = m_Wj
       PObs(21) = m_Hj
       PObs(22) = y_Hj
       PObs(23) = eta_Hj
       PObs(24) = eta_Hb
       PObs(25) = eta_Wb
       PObs(26) = eta_Wj
       PObs(27) = Dy_Hj
       PObs(28) = Deta_Hj
       PObs(29) = Deta_Hb
       PObs(30) = Deta_Wj
       PObs(31) = Deta_wb
       PObs(32) = Deta_Hep
       PObs(33) = m_Hbep
       PObs(34) = m_Hjep
       PObs(35) = Deta_tj
       PObs(36) = Deta_Ht
       PObs(37) = costheta_lj
endif


    else
       print *, "ObsSet not implemented TH",ObsSet
       stop
    endif


return
END SUBROUTINE Kinematics_TH




SUBROUTINE Kinematics_TBH(NPlus1PS,Mom,MomOrder,applyPSCut,NBin,PObs)
  use ModMisc
  use ModParameters
  implicit none
  integer :: NPlus1PS,MomOrder(1:11)
  real(8) :: Mom(1:4,1:11),zeros(1:9)
  integer :: tbar,ljet,Hig,inLeft,inRight,realp,bbar,lepM,nubar,Hdk1,Hdk2,k
  real(8) :: MomJet(1:4,1:7) !,MomJet_CHECK(1:4,1:7)                                                                
  real(8) :: MomHadr(1:4,0:8),MomMiss(1:4)
  real(8) :: pT_ATop,pT_higgs,pT_j,eta_j,eta_atop,eta_Higgs,y_atop,y_Higgs,y_j,eta_Hj,Dy_jtbar,Dy_Htbar,Dy_Hj,Deta_Hj,Deta_Htbar,y_Htbar,y_tbarj,m_Htbar,eta_tbarj,y_Hj,m_tbarj,m_Hj,eta_Htbar
  real(8) :: pT_miss,pT_em,pT_bbar,eta_em,eta_bbar,y_bbarem,y_jem,eta_W,eta_Hbbar,eta_Wj,eta_Wbbar,Deta_jbbar,Deta_Hbbar,Deta_Wj,Deta_Wbbar,Deta_tbarj,m_Hbbar,m_Wbbar,m_Wj,m_Hbbarem,m_Hjem,DR_bbarem,mtrans_W,HTtot,DR_jem,Deta_Hem,costheta_j,costheta_l,costheta_lj,DR_bbarj
  real(8),optional :: PObs(:)
  logical :: applyPSCut
  integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NumHadr


  applyPSCut = .false.

! momentum ordering                                                                                                                                                        
  tbar      = MomOrder(1)
  ljet   = MomOrder(2)
  Hig     = MomOrder(3)
  inLeft  = MomOrder(4)
  inRight = MomOrder(5)
  realp   = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  Hdk1    = MomOrder(10)
  Hdk2    = MomOrder(11)


!DEC$ IF(_CheckMomenta .EQ.1)                                                                                                                                              
   zeros(:) = 0d0
   if(TopDecays.eq.0) then
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ljet) - Mom(1:4,tbar) - Mom(1:4,Hig)
   else
      zeros(1:4) = Mom(1:4,inLeft)+Mom(1:4,inRight) - Mom(1:4,ljet) - Mom(1:4,bbar)-Mom(1:4,lepM)-Mom(1:4,nubar) - Mom(1:4,Hdk1)-Mom(1:4,Hdk2)
   endif
   if( NPlus1PS.eq.1 ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/Mom(1,inLeft)).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TBH(): ",NPlus1PS,zeros(1:4)
      print *, "momenta dump:"
      do k=1,12
         print *,k, Mom(1:4,k)
      enddo
   endif
   zeros(1) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_Top**2
   zeros(2) = (Mom(1:4,Hig).dot.Mom(1:4,Hig)) - M_H**2
   zeros(3) = (Mom(1:4,ljet).dot.Mom(1:4,ljet))
   zeros(4) = (Mom(1:4,bbar).dot.Mom(1:4,bbar))
   zeros(5) = (Mom(1:4,lepM).dot.Mom(1:4,lepM))
   zeros(6) = (Mom(1:4,nubar).dot.Mom(1:4,nubar))
   zeros(7) = (Mom(1:4,Hdk1).dot.Mom(1:4,Hdk1))
   zeros(8) = (Mom(1:4,Hdk2).dot.Mom(1:4,Hdk2))


   if( NPlus1PS.eq.1 ) zeros(9)=  Mom(1:4,realp).dot.Mom(1:4,realp)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:3)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TH(): ",zeros(1:3)
      print *, Mom(1:4,1:3)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:9)/Mom(1,inLeft)**2).gt.1d-8) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TH(): ",zeros(1:9)
      print *, "momenta dump:"
      print *, Mom(1:4,1:9)
   endif
   print *, "CHECK MOMENTA PASSED!"
!DEC$ ENDIF                                                                                                          
NBin(1:NumHistograms) = 0
MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)
MomMiss(1:4) = 0d0

IF( TOPDECAYS.EQ.0 ) THEN  ! no top decays    
   MomHadr(1:4,1) = Mom(1:4,ljet)                                                                           
    if(NPlus1PS.eq.0) then
        NumHadr = 1
    else
        NumHadr = 2
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 .OR. TOPDECAYS.eq.-1 ) THEN  ! di-leptonic decay                                                                                                         
   MomHadr(1:4,1) = Mom(1:4,bbar)     ! anti-b                                      
   MomHadr(1:4,2) = Mom(1:4,ljet)
 
    if(NPlus1PS.eq.0) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF





!---------------------- kT jet algorithm ---------------------------------                                                                                                     
! important: b-quarks need to be at the first two positions of MomHadr(:,:)                                                                                                    
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.                                                                                        
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order                                                                                                 

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))                                                                                                                           
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA                    
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))                                                                                                            
if( ObsSet.eq.92) then! tb+H production without top decays at LHC
       pT_ATop  = get_PT(Mom(1:4,tbar))
       pT_Higgs= get_PT(Mom(1:4,Hig))
       pT_j= get_PT(Mom(1:4,ljet))

       y_atop=get_eta(Mom(1:4,tbar))
       y_Higgs=get_eta(Mom(1:4,Hig))
       y_j=get_eta(Mom(1:4,ljet))

       eta_atop=get_pseudoeta(Mom(1:4,tbar))
       eta_Higgs=get_pseudoeta(Mom(1:4,Hig))
       eta_j=get_pseudoeta(Mom(1:4,ljet))

       eta_Hj=get_pseudoeta(Mom(1:4,Hig)+Mom(1:4,ljet))

       Dy_Hj=y_Higgs-y_j
       Dy_Htbar=y_Higgs-y_atop
       Dy_jtbar=y_j-y_atop

       Deta_Hj=eta_Higgs-eta_j
       Deta_Htbar=eta_Higgs-eta_atop

       eta_tbarj=get_pseudoeta(Mom(1:4,tbar)+Mom(1:4,ljet))
       eta_Hj=get_pseudoeta(Mom(1:4,Hig)+Mom(1:4,ljet))
       eta_Htbar=get_pseudoeta(Mom(1:4,Hig)+Mom(1:4,tbar))

       y_tbarj=get_eta(Mom(1:4,tbar)+Mom(1:4,ljet))
       y_Hj=get_eta(Mom(1:4,Hig)+Mom(1:4,ljet))
       y_Htbar=get_eta(Mom(1:4,Hig)+Mom(1:4,tbar))

       m_Htbar=get_MInv(Mom(1:4,Hig)+Mom(1:4,tbar))
       m_Hj=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet))
       m_tbarj=get_MInv(Mom(1:4,tbar)+Mom(1:4,ljet))


      ! binning                                                                                                
       NBin(1) = WhichBin(1,pT_ATop)
       NBin(2) = WhichBin(2,eta_ATop)
       NBin(3) = WhichBin(3,y_ATop)
       NBin(4) = WhichBin(4,pT_j)
       NBin(5) = WhichBin(5,y_j)
       NBin(6) = WhichBin(6,pT_Higgs)
       NBin(7) = WhichBin(7,m_tbarj)
       NBin(8) = WhichBin(8,y_Higgs)
       NBin(9) = WhichBin(9,eta_Higgs)
       NBin(10) = WhichBin(10,y_tbarj)
       NBin(11) = WhichBin(11,eta_tbarj)
       NBin(12) = WhichBin(12,Dy_jtbar)
       NBin(13) = WhichBin(13,m_Htbar)
       NBin(14) = WhichBin(14,m_Hj)
       NBin(15) = WhichBin(15,y_Htbar)
       NBin(16) = WhichBin(16,eta_Htbar)
       NBin(17) = WhichBin(17,y_Hj)
       NBin(18) = WhichBin(18,eta_Hj)
       NBin(19) = WhichBin(19,Dy_Htbar)
       NBin(20) = WhichBin(20,Deta_Htbar)
       NBin(21) = WhichBin(21,Dy_Hj)
       NBin(22) = WhichBin(22,Deta_Hj)

       if( present(PObs) ) then
          PObs(1) = pT_ATop
          PObs(2) = eta_ATop
          PObs(3) = y_ATop
          PObs(4) = pT_j
          PObs(5) = y_j
          PObs(6) = pT_Higgs
          PObs(7) = m_tbarj
          PObs(8) = y_Higgs
          PObs(9) = eta_Higgs
          PObs(10) = y_tbarj
          PObs(11) = eta_tbarj
          PObs(12) = Dy_jtbar
          PObs(13) = m_Htbar
          PObs(14) = m_Hj
          PObs(15) = y_Htbar
          PObs(16) = eta_Htbar
          PObs(17) = y_Hj
          PObs(18) = eta_Hj
          PObs(19) = Dy_Htbar
          PObs(20) = Deta_Htbar
          PObs(21) = Dy_Hj
          PObs(22) = Deta_Hj
       endif
       !-------------------------------------------------------                                                                   
    elseif (ObsSet .eq. 94 .or. ObsSet .eq. 96) then      ! leptonic top decay
       pT_em = get_PT(Mom(1:4,lepM))
       pT_miss = get_PT(Mom(1:4,nubar))
       pT_Higgs= get_PT(Mom(1:4,Hig))
       pT_j = get_PT(Mom(1:4,ljet))
       pT_bbar = get_PT(Mom(1:4,bbar))
       pT_ATop= get_PT(Mom(1:4,tbar))
       
       eta_em=get_pseudoeta(Mom(1:4,lepM))
       eta_bbar=get_pseudoeta(Mom(1:4,bbar))
       eta_atop=get_pseudoeta(Mom(1:4,tbar))
       eta_Higgs=get_pseudoeta(Mom(1:4,Hig))
       eta_j=get_pseudoeta(Mom(1:4,ljet))
       eta_W=get_pseudoeta(Mom(1:4,lepM)+Mom(1:4,nubar))
       eta_Hj=get_pseudoeta(Mom(1:4,ljet)+Mom(1:4,Hig))
       eta_Hbbar=get_pseudoeta(Mom(1:4,bbar)+Mom(1:4,Hig))
       eta_Wj=get_pseudoeta(Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,ljet))
       eta_Wbbar=get_pseudoeta(Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,bbar))

       DR_bbarj=get_R(Mom(1:4,bbar),Mom(1:4,ljet))
       DR_bbarem=get_R(Mom(1:4,bbar),Mom(1:4,lepM))
       DR_jem=get_R(Mom(1:4,ljet),Mom(1:4,lepM))


       IF (ObsSet .eq. 96) THEN
          if (pT_em .le. pT_lep_cut .or. pT_bbar .le. pT_bjet_cut .or. pT_j .le. pT_jet_cut) then
             applyPSCut = .true.
             RETURN
          endif

          if (abs(eta_em) .ge. eta_lep_cut .or. abs(eta_j) .ge. eta_jet_cut .or. abs(eta_bbar) .ge. eta_bjet_cut) then
             applyPSCut = .true.
             RETURN
          endif

          if (DR_bbarj .le. Rsep_jet .or. DR_bbarem .le. Rsep_lepjet .or. DR_jem .le. Rsep_lepjet) then
             applyPSCut = .true.
             RETURN
          endif
       ENDIF


 !      y_atop=get_eta(Mom(1:4,t))
       y_Higgs=get_eta(Mom(1:4,Hig))
       y_bbarem=get_eta(Mom(1:4,bbar)+Mom(1:4,lepM))
       y_jem=get_eta(Mom(1:4,ljet)+Mom(1:4,lepM))
       y_Hj=get_eta(Mom(1:4,ljet)+Mom(1:4,Hig))
       y_j=get_eta(Mom(1:4,ljet))

       costheta_j = get_CosTheta(Mom(1:4,ljet))
       costheta_l = get_CosTheta(Mom(1:4,lepM))
!       costheta_lj=costheta_l-costheta_j
       costheta_lj=get_cosalpha(Mom(1:4,lepM),Mom(1:4,ljet))

       Deta_jbbar=eta_j-eta_bbar
       Deta_Hj=eta_Higgs-eta_j
       Deta_Hbbar=eta_Higgs-eta_bbar
       Deta_Wj=eta_W-eta_j
       Deta_Wbbar=eta_W-eta_bbar
       Deta_Htbar=eta_Higgs-eta_atop
       Deta_Hem=eta_Higgs-eta_em
       Deta_tbarj=eta_atop-eta_j
       Dy_Hj=y_Higgs-y_j

       m_Hbbar=get_MInv(Mom(1:4,Hig)+Mom(1:4,bbar))
       m_Hj=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet))
       m_Wbbar=get_MInv(Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,bbar))
       m_Wj=get_MInv(Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,ljet))
       m_Hbbarem=get_MInv(Mom(1:4,Hig)+Mom(1:4,bbar)+Mom(1:4,lepM))
       m_Hjem=get_MInv(Mom(1:4,Hig)+Mom(1:4,ljet)+Mom(1:4,lepM))
       
       mtrans_W=get_MT(Mom(1:4,lepM),Mom(1:4,nubar))
       
       NBin(1) = WhichBin(1,pT_em)
       NBin(2) = WhichBin(2,eta_em)
       NBin(3) = WhichBin(3,pt_miss)
       NBin(4) = WhichBin(4,pt_j)
       NBin(5) = WhichBin(5,eta_j)
       NBin(6) = WhichBin(6,pt_bbar)
       NBin(7) = WhichBin(7,eta_bbar)
       NBin(8) = WhichBin(8,pt_Higgs)
       NBin(9) = WhichBin(9,mtrans_W)
       NBin(10) = WhichBin(10,y_Higgs)
       NBin(11) = WhichBin(11,y_bbarem)
       NBin(12) = WhichBin(12,y_jem)
       NBin(13) = WhichBin(13,eta_Higgs)
       NBin(14) = WhichBin(14,eta_W)
       NBin(15) = WhichBin(15,Deta_jbbar)
       NBin(16) = WhichBin(16,DR_bbarem)
       NBin(17) = WhichBin(17,DR_jem)
       NBin(18) = WhichBin(18,m_Hbbar)
       NBin(19) = WhichBin(19,m_Wbbar)
       NBin(20) = WhichBin(20,m_Wj)
       NBin(21) = WhichBin(21,m_Hj)
       NBin(22) = WhichBin(22,y_Hj)
       NBin(23) = WhichBin(23,eta_Hj)
       NBin(24) = WhichBin(24,eta_Hbbar)
       NBin(25) = WhichBin(25,eta_Wbbar)
       NBin(26) = WhichBin(26,eta_Wj)
       NBin(27) = WhichBin(27,Dy_Hj)
       NBin(28) = WhichBin(28,Deta_Hj)
       NBin(29) = WhichBin(29,Deta_Hbbar)
       NBin(30) = WhichBin(30,Deta_Wj)
       NBin(31) = WhichBin(31,Deta_wbbar)
       NBin(32) = WhichBin(32,Deta_Hem)
       NBin(33) = WhichBin(33,m_Hbbarem)
       NBin(34) = WhichBin(34,m_Hjem)
       NBin(35) = WhichBin(35,Deta_tbarj)
       NBin(36) = WhichBin(36,Deta_Htbar)
       NBin(37) = WhichBin(37,costheta_lj)


      if( present(PObs) ) then


       Pobs(1) = pT_em
       PObs(2) = eta_em
       PObs(3) = pt_miss
       PObs(4) = pt_j
       PObs(5) = eta_j
       PObs(6) = pt_bbar
       PObs(7) = eta_bbar
       PObs(8) = pt_Higgs
       PObs(9) = mtrans_W
       PObs(10) = y_Higgs
       PObs(11) = y_bbarem
       PObs(12) = y_jem
       PObs(13) = eta_Higgs
       PObs(14) = eta_W
       PObs(15) = Deta_jbbar
       PObs(16) = DR_bbarem
       PObs(17) = DR_jem
       PObs(18) = m_Hbbar
       PObs(19) = m_Wbbar
       PObs(20) = m_Wj
       PObs(21) = m_Hj
       PObs(22) = y_Hj
       PObs(23) = eta_Hj
       PObs(24) = eta_Hbbar
       PObs(25) = eta_Wbbar
       PObs(26) = eta_Wj
       PObs(27) = Dy_Hj
       PObs(28) = Deta_Hj
       PObs(29) = Deta_Hbbar
       PObs(30) = Deta_Wj
       PObs(31) = Deta_wbbar
       PObs(32) = Deta_Hem
       PObs(33) = m_Hbbarem
       PObs(34) = m_Hjem
       PObs(35) = Deta_tbarj
       PObs(36) = Deta_Htbar
       PObs(37) = costheta_lj
endif


                                                                                                                                                                               
    else
       print *, "ObsSet not implemented TBH",ObsSet
       stop
    endif


return
END SUBROUTINE Kinematics_TBH










SUBROUTINE Kinematics_TTBAR(NPlus1PS,MomExt,MomDK,applyPSCut,NBin,xJPsiFrag)
use ModMisc
use ModParameters
use ModSpinCorrel
implicit none
real(8), optional :: xJPsiFrag
integer :: NumHadr
real(8) :: MomExt(:,:),MomDK(:,:),MomJet(1:4,1:8)
real(8) :: MomLepP(1:4),MomLepM(1:4),MomBoost(1:4),MomZ(1:4),MomTops(1:4,1:2),MomJPsi(1:4),MomMiss(1:4),MomLepTR(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,iJet,NObsJet_Tree,NObsJet
real(8) :: pT_jet1,pT_jet2,pT_bjet1,pT_bjet2,pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,eta_bjet1,eta_bjet2,R_bb
real(8) :: MInv_lepP_bjet,CosPhi_LepLep,CosPsiT_LepLep,CosPsi_LepLep,Psi_LepLep,DeltaPhi
real(8) :: MInv_Tops,CosPhi_LepPZ,InvM_Lep,CosPhi_LepPlanes,HT,E_bjet1,E_bjet2,ET_miss
real(8) :: MInv_lepP_lepM,pT_lept,ET_lept,pT_x,pT_y,mT,MinvLept,E_lept,m_lb
real(8),save :: EMin=1d12
real(8) :: MomHadr(1:4,1:7),MomLept(1:4,1:4),MomJet_aux(1:4,1:7)
real(8) :: second_pair, first_pair, pT_bjet_min, mtop_reconstr,ptaux,pt_cand(4),pij(4)
real(8) :: ptmax, xij, MomBBar(1:4),diff,cosLeBb,Ebjets,r_sc,mT_lp
integer :: i1,i,j,n,k
real(8) :: MomLeptOrd(1:4,1:2), MomJetOrd(1:4,1:8), M_eff, Minv_Lept, ASV_ttbar
!-- for Zprime background
real(8) :: costheta_star, MomAux1(1:4), MomAux2(1:4), y_Top, y_Atop, m_ttbar, Mttb_cut, dphi_ttbar, costheta_scatter
real(8) :: pt_lep1, pt_lep2, eta_lep1, eta_lep2, dphi_LL, MassAux, MomTopsCMS(1:4,1:2)
real(8) :: MomBeam1(1:4), MomBeam2(1:4), MomBeam(1:4), nx(2:4), ny(2:4), nz(2:4), MomLeptTRF(1:4)
real(8) :: cosPhi, sinPhi, cosPhiBar, sinPhiBar, deltaSin, dPhiMinus, dPhiPlus

! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light particles, N-1=ATop, N=Top    ! HERE WAS A BUG: tops need to be at the end
!                          MomDK(:,1:7) : 1=ABot, 2=lep-, 3=ANeu, 4=Bot, 5=lep+, 6=Neu, 7=(Glu)
! MomLept(:,1:4): 1=lep-, 2=ANeu, 3=lep+, 4=Neu


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)




! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   NumHadr = 0
   if( NPlus1PS ) then
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   else
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
   endif
!-------------------------------------------------------
elseif( abs(TopDecays).eq.1 ) then  ! full leptonic decay
  MomLept(1:4,1) = MomDK(1:4,2)   ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)   ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)   ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)   ! Neu

  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,3) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,3) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,3) = MomDK(1:4,2)   ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)   ! q
  MomHadr(1:4,5) = MomDK(1:4,5)   ! qbar
  MomHadr(1:4,6) = MomDK(1:4,6)   ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,7) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,7) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.5 ) then  ! lept. Atop, hadr. top decay with J/Psi fragmentation

  MomHadr(1:4,1) = (1d0-xJPsiFrag) * MomDK(1:4,1)  ! b-jet remnand
  MomJPsi(1:4)   = xJPsiFrag * MomDK(1:4,1)        ! J/Psi momentum

  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.6) then  ! hadr. Atop, lept. top decay with J/Psi fragmentation


! MomHad(1:4,2) is the b-quark momentum because J/Psi and b-remnand jet are always recombined

  MomJPsi(1:4)   = xJPsiFrag * MomDK(1:4,4)        ! J/Psi momentum
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q

  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif

endif



!---------------------- kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)

! print *, PartList(1:NumHadr)
! print *, MomJet(1:4,1)
! print *, MomJet(1:4,2)
! print *, MomJet(1:4,3)
! print *, MomJet(1:4,4)
! print *, MomJet(1:4,5)

    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))

! print *, "NJet",NJet
! print *, JetList(1:NJet)
! print *, MomJet(1:4,1)
! print *, MomJet(1:4,2)
! print *, MomJet(1:4,3)
! print *, MomJet(1:4,4)
! print *, MomJet(1:4,5)
! pause

!--------------------------------------------------------------------------

if( ObsSet.eq.0 .or. ObsSet.eq.1 .or. ObsSet.eq.9) then! set of observables for ttb production without decays at Tevatron & LHC
!  eval kinematic variables
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

!     if( pT_Top.lt.800d0*GeV .or. pT_ATop.lt.800*GeV ) then!   this is for the boosted observable
!         applyPSCut = .true.
!         RETURN
!     endif

! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,eta_ATop)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)


!-------------------------------------------------------
elseif( ObsSet.eq.2 .or. ObsSet.eq.3) then! set of observables for ttb production with di-lept. decays at TEV and LHC

! request at least two b-jets
    if( .not.(NJet.ge.2 .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
   call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))
   call pT_order(2,MomJet(1:4,1:2))


! evaluate kinematic variables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))
    E_bjet1 = MomJet(1,1)
    E_bjet2 = MomJet(1,2)

    MinvLept = Get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    pT_lept = get_PT(MomLept(1:4,1)+MomLept(1:4,3))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    ET_lept = dsqrt(pT_lept**2 + MinvLept**2)
    ET_miss = get_ET(MomLept(1:4,2)+MomLept(1:4,4))
    pT_x = MomLept(2,1)+MomLept(2,2)+MomLept(2,3)+MomLept(2,4)
    pT_y = MomLept(3,1)+MomLept(3,2)+MomLept(3,3)+MomLept(3,4)
    mT = dsqrt( (ET_lept+ET_miss)**2 - (pT_x**2+pT_y**2) )

    E_lept = MomLept(1,1)+MomLept(1,3)

    pT_lepM = get_PT(MomLept(1:4,1))
    pT_lepP = get_PT(MomLept(1:4,3))

    eta_lepM = get_ETA(MomLept(1:4,1))
    eta_lepP = get_ETA(MomLept(1:4,3))

    MInv_lepP_lepM = get_MInv( MomLept(1:4,1)+MomLept(1:4,3) )

    if( get_MInv(MomLept(1:4,3)+MomJet(1:4,1)).lt.get_MInv( MomLept(1:4,3)+MomJet(1:4,2)) ) then 
        MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomJet(1:4,1) )
    else
        MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomJet(1:4,2) )
    endif
!MInv_lepP_bjet = get_MInv( MomLept(1:4,3)+MomDK(1:4,4) )

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))
    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))

    HT = pT_bjet1 + pT_bjet2 + pT_miss + pT_lepM + pT_lepP
    if( NJet.gt.2 ) HT = HT + get_PT(MomJet(1:4,3))

    MInv_Tops = get_MInv( MomTops(1:4,1)+MomTops(1:4,2) )

!     MomLepM(1:4)  = MomLept(1:4,1)
!     MomBoost(1)   =+MomTops(1,1)
!     MomBoost(2:4) =-MomTops(2:4,1)
!     call boost(MomLepM(1:4),MomBoost(1:4),m_top)
!     MomLepP(1:4)  = MomLept(1:4,3)
!     MomBoost(1)   =+MomTops(1,2)
!     MomBoost(2:4) =-MomTops(2:4,2)
!     call boost(MomLepP(1:4),MomBoost(1:4),m_top)
!     CosPhi_LepLep = (MomLepP(2)*MomLepM(2)+MomLepP(3)*MomLepM(3)+MomLepP(4)*MomLepM(4))/MomLepP(1)/MomLepM(1)

!     CosPsiT_LepLep = ( MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) )/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2)
    Psi_LepLep  = dacos( ( MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) + MomLept(4,1)*MomLept(4,3))/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2+MomLept(4,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2+MomLept(4,3)**2) )

   DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi


!                     MomBoost(1)   =+MomTops(1,1)
!                     MomBoost(2:4) =-MomTops(2:4,1)
!                     call boost(MomLept(1:4,1),MomBoost(1:4),m_top)
!                     MomBoost(1)   =+MomTops(1,2)
!                     MomBoost(2:4) =-MomTops(2:4,2)
!                     call boost(MomLept(1:4,3),MomBoost(1:4),m_top)
!                     Psi_LepLep = dacos( VectorProd(MomTops(2:4,1),MomLept(2:4,1))/MomLept(1,1)/dsqrt((MomTops(1,1))**2-m_top**2) )
!                     DeltaPhi   = dacos( VectorProd(MomTops(2:4,2),MomLept(2:4,3))/MomLept(1,3)/dsqrt((MomTops(1,2))**2-m_top**2) )

!                       if( get_eta(MomTops(1:4,1)).lt.0.8d0 ) then
!                           applyPSCut = .true.
!                           RETURN
!                       endif
!                     Psi_LepLep = ( MomLept(4,1)/MomLept(1,1) )
!                     DeltaPhi   = ( MomLept(4,3)/MomLept(1,3) )


   DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
   Psi_LepLep  = dacos( ( MomLept(4,1)*MomLept(4,3) +MomLept(2,1)*MomLept(2,3) + MomLept(3,1)*MomLept(3,3) )/dsqrt(MomLept(2,1)**2+MomLept(3,1)**2+MomLept(4,1)**2)/dsqrt(MomLept(2,3)**2+MomLept(3,3)**2+MomLept(4,3)**2) )




   Ebjets = MomJet(1,1)+MomJet(1,2)

! check cuts
    if( pT_bjet1.lt.pT_bjet_cut .OR. pT_bjet2.lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_bjet1).gt.eta_bjet_cut .OR. abs(eta_bjet2).gt.eta_bjet_cut) then
       applyPSCut = .true.
        RETURN
    endif


    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepM).gt.eta_lep_cut .OR. abs(eta_lepP).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif






!  additional cuts to observe spin correlations
!     if( pT_lepM.gt.50d0*GeV .OR. pT_lepP.gt.50d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!
!     if( MInv_lepP_lepM.gt.100d0*GeV .OR. MInv_lepP_lepM.gt.100d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif

!    if( MInv_Tops.gt.400d0*GeV ) then
!      applyPSCut = .true.
!       RETURN
!    endif

!     if( pT_bjet1.gt.100d0*GeV .OR. pT_bjet2.gt.100d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif





! !  additional cuts for W+dijet search
!
!     if( pT_miss.lt.25d0*GeV ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( pT_bjet1.lt.30d0*GeV .OR. pT_bjet2.lt.30d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif
!
!     if( abs(eta_bjet1).gt.2.4d0 .OR. abs(eta_bjet2).gt.2.4d0 ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( Get_PT( MomJet(1:4,1)+MomJet(1:4,2) ).lt.40d0*GeV ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!     if( dabs(Get_ETA(MomJet(1:4,1))  - Get_ETA(MomJet(1:4,2)) ).gt.2.5d0 ) then
!        applyPSCut = .true.
!         RETURN
!     endif
!
!    if( get_pt(MomJet(1:4,1)).gt.get_pt(MomJet(1:4,2)) ) then
!         DeltaPhi = dabs( Get_PHI(MomJet(1:4,1)) - Get_PHI(MomLept(1:4,2)+MomLept(1:4,4))  )
!         if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
!         if( abs(DeltaPhi).lt.0.4d0 ) then
!               applyPSCut = .true.
!               RETURN
!         endif
!    else
!         DeltaPhi = dabs( Get_PHI(MomJet(1:4,2)) - Get_PHI(MomLept(1:4,2)+MomLept(1:4,4))  )
!         if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
!         if( abs(DeltaPhi).lt.0.4d0 ) then
!               applyPSCut = .true.
!               RETURN
!         endif
!    endif
!
!
!
! ! now, pick the lepton
!    if( get_pt(MomLept(1:4,1)).gt.20d0*GeV .AND. get_pt(MomLept(1:4,3)).lt.10d0*GeV ) then! the lepton is MomLept(1:4,1)
!               if( get_ETA(MomLept(1:4,1)).gt.1d0 ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).gt.76d0*GeV .and. get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).lt.106d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_R(MomLept(1:4,1),MomJet(1:4,1)).le.0.52d0 .or. get_R(MomLept(1:4,1),MomJet(1:4,2)).le.0.52d0) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!
!               MomMiss(1:4) = MomLept(1:4,2)+MomLept(1:4,4)
!               MomMiss(4) = 0d0
!               MomLepTR(1:4) = MomLept(1:4,1)
!               MomLepTR(4) = 0d0
!               mT = dsqrt( 2d0*get_PT(MomLept(1:4,1))*pT_miss*(1d0-Get_CosAlpha(MomMiss(1:4),MomLepTR(1:4))) )
!               if( mT.lt.30d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!    elseif( get_pt(MomLept(1:4,3)).gt.20d0*GeV .AND. get_pt(MomLept(1:4,1)).lt.10d0*GeV ) then! the lepton is MomLept(1:4,2)
!               if( get_ETA(MomLept(1:4,3)).gt.1d0 ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).gt.76d0*GeV .and. get_MInv(MomLept(1:4,1)+MomLept(1:4,3)).lt.106d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!               if( get_R(MomLept(1:4,3),MomJet(1:4,1)).le.0.52d0 .or. get_R(MomLept(1:4,3),MomJet(1:4,2)).le.0.52d0) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!
!               MomMiss(1:4) = MomLept(1:4,2)+MomLept(1:4,4)
!               MomMiss(4) = 0d0
!               MomLepTR(1:4) = MomLept(1:4,3)
!               MomLepTR(4) = 0d0
!               mT = dsqrt( 2d0*get_PT(MomLept(1:4,3))*pT_miss*(1d0-Get_CosAlpha(MomMiss(1:4),MomLepTR(1:4))) )
!               if( mT.lt.30d0*GeV ) then
!                   applyPSCut = .true.
!                   RETURN
!               endif
!    else
!        applyPSCut = .true.
!        RETURN
!    endif
!   MInv_lepP_lepM = get_MInv( MomJet(1:4,1)+MomJet(1:4,2) )



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,eta_lepP)
    NBin(5) = WhichBin(5,HT)
    NBin(6) = WhichBin(6,MInv_lepP_bjet)
    NBin(7) = WhichBin(7,get_MInv(MomLept(1:4,3)+MomDK(1:4,4))) ! careful! MomDK(:,4) is not IR save
!    NBin(7) = WhichBin(7,Psi_LepLep)
    NBin(8) = WhichBin(8,DeltaPhi)
    NBin(9) = WhichBin(9,Ebjets)

    NBin(10) = WhichBin(10,ET_miss)
    NBin(11) = WhichBin(11,E_lept)
    NBin(12) = WhichBin(12,pT_bjet1)
    NBin(13) = WhichBin(13,pT_bjet2)
    NBin(14) = WhichBin(14,E_bjet1)
    NBin(15) = WhichBin(15,E_bjet2)
    NBin(16) = WhichBin(16,MInv_lepP_lepM)
    NBin(17) = WhichBin(17,ET_lept)
    NBin(18) = WhichBin(18,mT)

    NBin(19) = WhichBin(19,eta_lepM)
    NBin(20) = WhichBin(20,eta_lepP)



!-------------------------------------------------------
elseif( ObsSet.eq.4 ) then! set of observables for ttb production with hadr. Atop, lept. top decay

!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
!     if( NObsJet.ne.NObsJet_Tree ) then!!!    CAREFUL: this cuts out the additional hard jet: only for combination with ttbjet
        applyPSCut = .true.
        RETURN
    endif


    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif


!   construct ttb momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_lepP)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,m_lb)
    NBin(12)= WhichBin(12,pT_Top)



!-------------------------------------------------------
elseif( ObsSet.eq.5 ) then! set of observables for ttb production with hadr. Atop, lept. top decay


!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut .or. get_R(MomLept(1:4,3),MomJet(1:4,1)).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut .or. get_R(MomLept(1:4,3),MomJet(1:4,2)).lt.Rsep_LepJet ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut .and. get_R(MomLept(1:4,3),MomJet(1:4,k)).gt.Rsep_LepJet ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif



    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    mT_lp = get_MT(MomLept(1:4,3),MomLept(1:4,4))! this is the transverse W mass

    MInv_Tops = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    ASV_ttbar = Get_Mom3Abs((/0d0,MomTops(2,1)+MomTops(2,2),MomTops(3,1)+MomTops(3,2),0d0/))/(pT_Top+pT_ATop + 1d-10)


    if( ASV_ttbar.gt.0.4d0 ) then
        applyPSCut = .true.
        RETURN
    endif




!     if( ET_miss+mT_lp.lt.60*GeV  ) then
!         applyPSCut = .true.
!         RETURN
!     endif


!    if( HT.lt.HT_cut ) then
!        applyPSCut = .true.
!        RETURN
!    endif

    if(  MInv_Tops.lt.HT_cut ) then!  --> 460 pushes <Ehat>=621     ! 660 pushes <Ehat>=835          ! 630 860
        applyPSCut = .true.
        RETURN
    endif


!   construct hadr. W momentum
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif

   
   ! this is Ehat for <EHat> in NHisto=12
   MInv_LB= get_MInv(MomExt(1:4,1)+MomExt(1:4,2))

! binning
    NBin(1:4) = 0
    if( 0.5d0*(pT_ATop+pT_Top).lt.800d0*GeV ) then
        NBin(1) = WhichBin(1,0.5d0*(pT_ATop+pT_Top)  )
        NBin(2) = WhichBin(2,0.5d0*(eta_ATop+eta_Top))
        if( MInv_Tops.gt.660d0*GeV ) then
            NBin(3) = WhichBin(3,0.5d0*(pT_ATop+pT_Top)  )
            NBin(4) = WhichBin(4,0.5d0*(eta_ATop+eta_Top))    
        endif
    endif
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,MInv_Tops)  !m_lb)





!-------------------------------------------------------
elseif( ObsSet.eq.6 ) then! set of observables for ttb production with hadr. Atop, lept. top decay


!   request at least two b-jets
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))

    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))
    eta_bjet1 = get_ETA(MomJet(1:4,1))
    eta_bjet2 = get_ETA(MomJet(1:4,2))

!   check if b-jets pass cuts
    if(  pT_bjet1.lt.pT_bjet_cut .or. abs(eta_bjet1).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( pT_bjet2.lt.pT_bjet_cut .or. abs(eta_bjet2).gt.eta_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    HT = pT_bjet1 + pT_bjet2

!   determine observable light jets
    NObsJet = 2! these are the b-jets
    do k=3,NJet
        if( get_PT(MomJet(1:4,k)).gt.pT_jet_cut .and. abs(get_ETA(MomJet(1:4,k))).lt.eta_jet_cut ) then! count jets outside beam pipe
            NObsJet = NObsJet +1
            if( k.ne.NObsJet ) MomJet(1:4,NObsJet) = MomJet(1:4,k)
            HT = HT + get_PT(MomJet(1:4,k))
        endif
    enddo

    NObsJet_Tree = 4! request two b-jets and at least two light jets
    if( NObsJet.lt.NObsJet_Tree ) then
!     if( NObsJet.ne.NObsJet_Tree ) then!!!    CAREFUL: this cuts out the additional hard jet: only for combination with ttbjet
        applyPSCut = .true.
        RETURN
    endif



    pT_lepP = get_PT(MomLept(1:4,3))
    eta_lepP = get_ETA(MomLept(1:4,3))

    ET_miss  = get_ET(MomLept(1:4,4))
    HT = HT + pT_lepP + ET_miss

    if( pT_bjet1.gt.pT_bjet2 ) then
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,1))
    else
        m_lb = get_Minv(MomLept(1:4,3)+MomJet(1:4,2))
    endif

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    eta_ATop = get_ETA(MomTops(1:4,1))
    eta_Top  = get_ETA(MomTops(1:4,2))


! check cuts
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepP).gt.eta_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( ET_miss.lt.pT_miss_cut ) then
        applyPSCut = .true.
        RETURN
    endif

!     if( HT.lt.HT_cut ) then
!         applyPSCut = .true.
!         RETURN
!     endif



!   construct hadr. ttb pT
    MomTops(1:4,1) = MomJet(1:4,1)+MomJet(1:4,2)+MomLept(1:4,3)+MomLept(1:4,4) + MomJet(1:4,3)+MomJet(1:4,4)
    if( dabs( get_MInv(MomJet(1:4,3)+MomJet(1:4,4))-M_W ).lt.20d0*GeV ) then!   require a 20GeV window around M_W
        pT_Top = get_pT(MomTops(1:4,1))
    else
        pT_Top   = -1d0
    endif



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,eta_ATop)
    NBin(3) = WhichBin(3,pT_Top)
    NBin(4) = WhichBin(4,eta_Top)
    NBin(5) = WhichBin(5,eta_ATop)
    NBin(6) = WhichBin(6,eta_Top)
    NBin(7) = WhichBin(7,pT_lepP)
    NBin(8) = WhichBin(8,eta_lepP)
    NBin(9) = WhichBin(9,ET_miss)
    NBin(10)= WhichBin(10,HT)
    NBin(11)= WhichBin(11,m_lb)
    NBin(12)= WhichBin(12,pT_Top)





!-------------------------------------------------------
elseif( ObsSet.eq.7 ) THEN! set of observables for ttb production with lept. top and J/Psi fragmentation, hadr. Atop decay at LHC

    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))

    NObsJet_Tree = 4
    NObsJet = 0
    HT = 0d0
    do n=1,NJet
        pT_jet1 = get_pT(MomJet(1:4,n))
        if( pT_jet1.gt.pT_jet_cut ) NObsJet = NObsJet +1
        HT = HT + pT_jet1
    enddo

    if( NObsJet.lt.NObsJet_Tree ) then
        applyPSCut = .true.
        RETURN
    endif

    if( HT.lt.HT_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_lepP = get_PT(MomLept(1:4,3))
    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
! print *,"removed cuts",  NBin(1:3)


    MInv_LB = get_MInv( xJPsiFrag*MomDK(1:4,4)+MomLept(1:4,3) )
    pT_lepP = get_PT(MomLept(1:4,3))
    pT_bjet1 = get_pT(MomJPsi(1:4)) ! pT(J/Psi)


! binning
    NBin(1) = WhichBin(1,MInv_LB)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lepP)




!-------------------------------------------------------
elseif( ObsSet.eq.8 ) then! set of observables for ttb spin correlations at LHC (di-lept. decay)

! request at least two b-jets
    if( .not.(NJet.ge.2 .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif
   call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))


! evaluate kinematic variables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_bjet2 = get_PT(MomJet(1:4,2))

    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))
    pT_miss = get_PT(MomLept(1:4,2)+MomLept(1:4,4))

    pT_lepM = get_PT(MomLept(1:4,1))
    pT_lepP = get_PT(MomLept(1:4,3))

    eta_lepM = get_ETA(MomLept(1:4,1))
    eta_lepP = get_ETA(MomLept(1:4,3))

    DeltaPhi = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi

    MInv_Tops = get_MInv( MomTops(1:4,1)+MomTops(1:4,2) )

! check cuts
!     if( pT_lepM.gt.50d0*GeV .OR. pT_lepP.gt.50d0*GeV ) then
!         applyPSCut = .true.
!         RETURN
!     endif

    if( pT_bjet1.lt.pT_bjet_cut .OR. pT_bjet2.lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then
       applyPSCut = .true.
        RETURN
    endif

    if( abs(eta_lepM).gt.eta_lep_cut .OR. abs(eta_lepP).gt.eta_lep_cut) then
       applyPSCut = .true.
        RETURN
    endif

!    if( MInv_Tops.gt.400d0*GeV ) then
!      applyPSCut = .true.
!       RETURN
!    endif


   MomHadr(1:4,1) = MomJet(1:4,1)
   MomHadr(1:4,2:3) = MomDK(1:4,2:3)

   MomHadr(1:4,4) = MomJet(1:4,2)
   MomHadr(1:4,5:6) = MomDK(1:4,5:6)
   if( Collider.eq.1 ) then
     r_sc = calc_rgg(MomExt(1:4,1:4),MomHadr(1:4,1:6))
   else
     r_sc = calc_rqq(MomExt(1:4,1:4),MomHadr(1:4,1:6))
   endif

! binning
    NBin(1) = WhichBin(1,r_sc)
    NBin(2) = WhichBin(2,pT_lepM)
    NBin(3) = WhichBin(3,pT_lepP)
    NBin(4) = WhichBin(4,DeltaPhi)

elseif (ObsSet.EQ.60) then ! set of observables for ttb production without decays, for Zprime background

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    if (M_TTbar.lt.Mttb_cut) then
        applyPSCut = .true.
        RETURN
    endif


    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))

!-------------------------------------------------------
elseif( ObsSet.eq.61 ) then! set of observables for ttb production with di-lept. decays, for Zprime background

   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

   call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:8) = MomJet(1:4,1:8)
   call pT_order(NumHadr,MomJetOrd(1:4,1:NumHadr))


!   M_LL = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Phase space cuts

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

    if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi
!pause

    else
       dPhiPlus=20d0
       dPhiMinus=20d0

    endif

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)

endif

return
END SUBROUTINE










SUBROUTINE Kinematics_TTbarETmiss(NPlus1PS,Mom,MomOrder,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr,MomOrder(1:13)
real(8) :: Mom(1:4,1:15),zeros(1:13)
real(8) :: MomJet(1:4,1:7)
real(8) :: MomHadr(1:4,0:8)
real(8) :: MomBoost(1:4),MomMiss(1:4),MomObs(1:4),MomAux(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet,k,NObsJet_Tree
real(8) :: pT_lepM,pT_lepP,pT_miss,pT_ATop,pT_Top
real(8) :: eta_ATop,eta_Top,eta_lepM,eta_lepP,m_lb
real(8) :: pT_jet(1:7),eta_jet(1:7),m_X0,phi_ll,m_ll,MTW,MTeff
integer :: Htbar,Ht,X0bar,X0,tbar,t,realp,bbar,lepM,nubar,b,lepP,nu,qdn,qbup,qbdn,qup,Lep,Neu,ib_lep,ib_had
real(8) :: Stopness,rpT,H_T



! momentum ordering
  Htbar   = MomOrder(1)
  Ht      = MomOrder(2)
  X0bar   = MomOrder(3)
  X0      = MomOrder(4)
  tbar    = MomOrder(5)
  t       = MomOrder(6)
  bbar    = MomOrder(7)
  lepM    = MomOrder(8)
  nubar   = MomOrder(9)
  b       = MomOrder(10)
  lepP    = MomOrder(11)
  nu      = MomOrder(12)
  realp   = MomOrder(13)

  qdn    = lepM
  qbup   = nubar
  qbdn   = lepP
  qup    = nu

IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 .or. ObsSet.eq.34 .or. ObsSet.eq.35 ) then
   IF(XTOPDECAYS.EQ.1) m_X0 = m_BH
   IF(XTOPDECAYS.EQ.2) m_X0 = m_A0
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 .or. ObsSet.eq.44 .or. ObsSet.eq.45 ) then
   m_X0 = m_Chi
endif

!DEC$ IF(_CheckMomenta .EQ.1)
IF(XTOPDECAYS.NE.0) THEN
   zeros(1:4) = Mom(1:4,1)+Mom(1:4,2) - (Mom(1:4,X0bar)+Mom(1:4,X0)+Mom(1:4,bbar)+Mom(1:4,lepM)+Mom(1:4,nubar)+Mom(1:4,b)+Mom(1:4,lepP)+Mom(1:4,nu))
   if( NPlus1PS ) zeros(1:4) = zeros(1:4) - Mom(1:4,realp)
   if( any(abs(zeros(1:4)/m_HTop).gt.1d-8) ) then
      print *, "ERROR: energy-momentum violation 1 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:4),NPlus1PS
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif

   if( Correction.ne.5 .or. .not.NPlus1PS ) then
          zeros(1:4) = Mom(1:4,Htbar)-Mom(1:4,X0bar)-Mom(1:4,tbar)
          zeros(5:8) = Mom(1:4,Ht)-Mom(1:4,X0)-Mom(1:4,t)
          if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 2 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:8),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif

          zeros(1:4) = Mom(1:4,tbar)-Mom(1:4,bbar)-Mom(1:4,LepM)-Mom(1:4,nubar)
          zeros(5:8) = Mom(1:4,t)-Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu)
          if( any(abs(zeros(1:8)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 3 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:8),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif
    else
          zeros(1:4) = Mom(1:4,Htbar)-Mom(1:4,X0bar) + Mom(1:4,Ht)-Mom(1:4,X0) & 
                     - Mom(1:4,bbar)-Mom(1:4,LepM)-Mom(1:4,nubar) -Mom(1:4,b)-Mom(1:4,lepP)-Mom(1:4,nu) - Mom(1:4,realp)
          if( any(abs(zeros(1:4)/m_HTop).gt.1d-8) ) then
              print *, "ERROR: energy-momentum violation 3 in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:4),NPlus1PS
              print *, "momentum dump:"
              print *, Mom(:,:)
          endif
    endif



   zeros(:) = 0d0
IF( ObsSet.eq.31 .or. ObsSet.eq.32 .or. ObsSet.eq.33 .or. ObsSet.eq.34 .or. ObsSet.eq.35 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_HTop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_HTop**2
elseif( ObsSet.eq.41 .or. ObsSet.eq.42 .or. ObsSet.eq.43 .or. ObsSet.eq.44 .or. ObsSet.eq.45 ) then
   zeros(1) = (Mom(1:4,HTbar).dot.Mom(1:4,HTbar))  - m_STop**2
   zeros(2) = (Mom(1:4,HT).dot.Mom(1:4,HT)) - m_STop**2
endif
   zeros(3) = (Mom(1:4,X0bar).dot.Mom(1:4,X0bar)) - m_X0**2
   zeros(4) = (Mom(1:4,X0).dot.Mom(1:4,X0)) - m_X0**2
   zeros(5) = (Mom(1:4,tbar).dot.Mom(1:4,tbar)) - m_SMTop**2
   zeros(6) = (Mom(1:4,t).dot.Mom(1:4,t)) - m_SMTop**2
   zeros(7) = (Mom(1:4,bbar).dot.Mom(1:4,bbar))
   zeros(8) = (Mom(1:4,b).dot.Mom(1:4,b))
   zeros(9) = (Mom(1:4,lepm).dot.Mom(1:4,lepm))
   zeros(10) = (Mom(1:4,lepp).dot.Mom(1:4,lepp))
   zeros(11) = (Mom(1:4,nubar).dot.Mom(1:4,nubar))
   zeros(12) = (Mom(1:4,nu).dot.Mom(1:4,nu))
   if(NPlus1PS) zeros(13) = (Mom(1:4,realp).dot.Mom(1:4,realp))

   if( any(abs(zeros(1:13)/1d0).gt.1d-5) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTbarETmiss: ",zeros(1:13),NPlus1PS
      print *, "momentum dump:"
      print *, Mom(:,:)
   endif
ENDIF
!DEC$ ENDIF



applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)

MomMiss(1:4) = 0d0
MomObs(1:4)  = 0d0

! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
!-------------------------------------------------------

IF( TOPDECAYS.EQ.0 ) THEN  ! no decays
    if(.not.NPlus1PS) then
        NumHadr = 0
    else
        NumHadr = 1
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.1 ) THEN  ! di-leptonic decay
    MomHadr(1:4,1) = Mom(1:4,bbar)
    MomHadr(1:4,2) = Mom(1:4,b)     ! Bot
    if(.not.NPlus1PS) then
        NumHadr = 2
    else
        NumHadr = 3
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
    endif

ELSEIF( TOPDECAYS.EQ.3 ) THEN  ! lept. Atop, hadr. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbdn)
   MomHadr(1:4,4) = Mom(1:4,qup)
   Lep = LepM
   Neu = nubar
   if(.not.NPlus1PS) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)   ! q/qbar/glu
   endif

ELSEIF( TOPDECAYS.EQ.4 ) THEN  ! hadr. Atop, lept. top decay
   MomHadr(1:4,1) = Mom(1:4,bbar)
   MomHadr(1:4,2) = Mom(1:4,b)
   MomHadr(1:4,3) = Mom(1:4,qbup)
   MomHadr(1:4,4) = Mom(1:4,qdn)
   Lep = LepP
   Neu = nu
   if(.not.NPlus1PS) then
        NumHadr = 4
   else
        NumHadr = 5
        MomHadr(1:4,NumHadr) = Mom(1:4,realp)
   endif

ELSE
  call Error("this decay is not yet implemented")
ENDIF



!---------------------- kT jet algorithm ---------------------------------
! important: b-quarks need to be at the first two positions of MomHadr(:,:)
! recombination of MomHadr(1:4,i) and MomHadr(1:4,j) results in MomJet(1:4,i) with i<j.
! but take care: some MomJet(1:4,i) can be zero, this is why we apply pt_order

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr)) ! hard protojets in beam pipe are counted as jets
    call pT_order(2,MomJet(1:4,1:2))
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))
! call SwitchEnergyComponent(MomHadr(1:4,1:NumHadr))
! call fastjetppgenkt(MomHadr(1:4,1:NumHadr),NumHadr,Rsep_jet,-1d0,MomJet(1:4,1:NumHadr),NJet)! 4th argument:  +1d0=kT,  -1d0=anti-kT,   0d0=CA
! call SwitchEnergyComponentBack(MomJet_CHECK(1:4,1:NumHadr))

!------------------------ cuts and binning --------------------------------
if( ObsSet.eq.31 .or. ObsSet.eq.41) then! no decays


    pT_jet(1) = get_PT(Mom(1:4,HTbar))


! binning
    NBin(1)= WhichBin(1,pT_jet(1))


!-------------------------------------------------------
elseif( ObsSet.eq.32 .or. ObsSet.eq.34 .or. ObsSet.eq.42 .or. ObsSet.eq.44 )  then! set of observables for TTbar -> ttbar + ETmiss  in di-lept. top decays

! request at least two b-jets
    NObsJet_Tree = 2
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_jet(1) = get_PT(MomJet(1:4,1))
    pT_jet(2) = get_PT(MomJet(1:4,2))
    eta_jet(1)= get_ETA(MomJet(1:4,1))
    eta_jet(2)= get_ETA(MomJet(1:4,2))

    pT_lepM = get_PT(Mom(1:4,lepM))
    pT_lepP = get_PT(Mom(1:4,lepP))
    eta_lepM = get_ETA(Mom(1:4,lepM))
    eta_lepP = get_ETA(Mom(1:4,lepP))

    MomMiss(1:4) = Mom(1:4,nu)+Mom(1:4,nubar)+Mom(1:4,X0)+Mom(1:4,X0bar)
    pT_miss = get_ET( MomMiss(1:4) )! note that this is ET and not pT


    phi_ll = dabs( Get_PHI(Mom(1:4,lepM)) - Get_PHI(Mom(1:4,lepP)) )
    if( phi_ll.gt.Pi ) phi_ll=2d0*Pi-phi_ll

    m_ll = get_MInv(Mom(1:4,lepM)+Mom(1:4,lepP))


    MTW = dsqrt(  2d0*pT_lepP*get_ET( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lepP),MomMiss(1:4))) )! let's define MTW with ET instead of pT in accordance with ATLAS

    MTeff = pt_jet(1) + pt_jet(2) + pT_lepP + pT_lepM + pT_miss
    if( NJet.eq.3 ) MTeff = MTeff + pt_jet(3)

    H_T = pt_jet(1) + pt_jet(2) + MTW ! definition in accordance with ATLAS


! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_lepM.lt.pT_lep_cut .OR. pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_lepM).gt.eta_lep_cut .or. abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then! note that this is ET and not pT
       applyPSCut = .true.
        RETURN
    endif


! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,phi_ll)
    NBin(5)= WhichBin(5,m_ll)
    NBin(6)= WhichBin(6,H_T)
    NBin(7)= WhichBin(7,MTW)
    NBin(8)= WhichBin(8,MTeff)





!-------------------------------------------------------
elseif( ObsSet.eq.33 .or. ObsSet.eq.35 .or. ObsSet.eq.36 .or.ObsSet.eq.43 .or. ObsSet.eq.45 .or. ObsSet.eq.46 ) then! set of observables for TTbar -> ttbar + ETmiss  in semi-hadr. top decays
! request at least two b-jets and 2 jets from hadronic W decay
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_jet(1) = get_PT( MomJet(1:4,1))
    pT_jet(2) = get_PT( MomJet(1:4,2))
    pT_jet(3) = get_PT( MomJet(1:4,3))
    pT_jet(4) = get_PT( MomJet(1:4,4))
    if( NJet.eq.5 ) pT_jet(5) = get_PT( MomJet(1:4,5))
    eta_jet(1)= get_ETA(MomJet(1:4,1))
    eta_jet(2)= get_ETA(MomJet(1:4,2))
    eta_jet(3)= get_ETA(MomJet(1:4,3))
    eta_jet(4)= get_ETA(MomJet(1:4,4))
    if( NJet.eq.5 ) eta_jet(5)= get_ETA(MomJet(1:4,5))

    pT_lepP  = get_PT(Mom(1:4,lep))
    eta_lepP = get_ETA(Mom(1:4,lep))

    MomMiss(1:4) = Mom(1:4,nu)+Mom(1:4,X0)+Mom(1:4,X0bar)
    pT_miss = get_ET( MomMiss(1:4) )! note that this is ET and not pT

    pT_Top = get_PT(Mom(1:4,t))
    eta_Top = get_ETA(Mom(1:4,t))


!    MTW = Get_MT(Mom(1:4,lep),MomMiss(1:4))!    ==  dsqrt(  2d0*pT_lepP*get_PT( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lep),MomMiss(1:4))) )
   MTW = dsqrt(  2d0*pT_lepP*get_ET( MomMiss(1:4) )*(1d0-Get_CosPhi(Mom(1:4,lep),MomMiss(1:4))) )! let's define MTW with ET instead of pT in accordance with ATLAS

   MTeff = pt_jet(1) + pt_jet(2) + pt_jet(3) + pt_jet(4) + pT_lepP + pT_miss
   if( NJet.eq.5 ) MTeff = MTeff + pt_jet(5)

   H_T = pt_jet(1) + pt_jet(2) + pt_jet(3) + pt_jet(4) + MTW ! definition in accoradance with ATLAS


! check cuts
    if( pT_jet(1).lt.pT_bjet_cut .or. pT_jet(2).lt.pT_bjet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(1)).gt.eta_bjet_cut .or. abs(eta_jet(2)).gt.eta_bjet_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_jet(3).lt.pT_jet_cut .or. pT_jet(4).lt.pT_jet_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_jet(3)).gt.eta_jet_cut .or. abs(eta_jet(4)).gt.eta_jet_cut) then
        applyPSCut = .true.
        RETURN
    endif


    if( pT_lepP.lt.pT_lep_cut ) then
        applyPSCut = .true.
        RETURN
    endif
    if( abs(eta_lepP).gt.eta_lep_cut) then
        applyPSCut = .true.
        RETURN
    endif

    if( pT_miss.lt.pT_miss_cut ) then! note that this is ET and not pT
       applyPSCut = .true.
        RETURN
    endif
  
   if( MTW.lt.MTW_cut ) then
      applyPSCut = .true.
       RETURN
   endif


!  topness: step 1
   if( get_MInv2(MomJet(1:4,1)+Mom(1:4,lep)) .lt. get_MInv2(MomJet(1:4,2)+Mom(1:4,lep)) ) then! select the b-jet that belongs to the leptonic side
        ib_lep=1
        ib_had=2
   else
        ib_lep=2
        ib_had=1
   endif
   Stopness = (M_W**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)))**2/(5*GeV)**4  & 
            + (M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,ib_lep)))**2/(15*GeV)**4

!  topness: step 2
   if( NJet.eq.4 ) then! if there are 4 jets, we assume that all come from the decay process
      Stopness =  Stopness  +  &
                + (M_top**2 - Get_MInv2(MomJet(1:4,3)+MomJet(1:4,4)+MomJet(1:4,ib_had)))**2/(15*GeV)**4 
   elseif( NJet.eq.5 ) then! if there are 5 jets, we assume that the hardest non-bjet comes from the production process
      Stopness =  Stopness  +  &
                + (M_top**2 - Get_MInv2(MomJet(1:4,4)+MomJet(1:4,5)+MomJet(1:4,ib_had)))**2/(15*GeV)**4 
   endif

!  topness: step 3
   if( NJet.eq.4 ) then! if there are 4 jets, we assume that all come from the decay process
      Stopness =  Stopness  +  &
               + (4*M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,1)+MomJet(1:4,2)+MomJet(1:4,3)+MomJet(1:4,4)))**2/(1000*GeV)**4
   elseif( NJet.eq.5 ) then! if there are 5 jets, we assume that the hardest non-bjet comes from the production process
      Stopness =  Stopness  +  &
               + (4*M_top**2 - Get_MInv2(Mom(1:4,lep)+MomMiss(1:4)+MomJet(1:4,1)+MomJet(1:4,2)+MomJet(1:4,4)+MomJet(1:4,5)))**2/(1000*GeV)**4
   endif


    
    rpT = (pT_jet(1)-pT_lepP)/(pT_jet(1)+pT_lepP)
    



! binning
    NBin(1)= WhichBin(1,pT_lepP)
    NBin(2)= WhichBin(2,eta_lepP)
    NBin(3)= WhichBin(3,pT_miss)
    NBin(4)= WhichBin(4,H_T)
    NBin(5)= WhichBin(5,pt_jet(NJet))! softest jet
    NBin(6)= WhichBin(6,MTW)
    NBin(7)= WhichBin(7,MTeff)
    NBin(8)= WhichBin(8,pT_Top)
    NBin(9)= WhichBin(9,eta_Top)
    NBin(10)= WhichBin(10,dlog(Stopness))
    NBin(11)= WhichBin(11,rpT)




!-------------------------------------------------------
else
  print *, "ObsSet not implemented TTBETmiss",ObsSet
  stop
endif


return
END SUBROUTINE






SUBROUTINE Kinematics_TTBARZprime(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr
real(8) :: MomExt(:,:),MomDK(:,:),MomJet(1:4,1:8),MomAux1(1:4),MomAux2(1:4)
real(8) :: MomTops(1:4,1:2),MomBoost(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet_Tree,NObsJet
real(8) :: y_Top,y_ATop,pT_Top,pT_ATop,M_TTbar,Dphi_TTbar,pT_Lept,y_Lept,Dphi_LL,M_LL
real(8) :: pT_LeptP, pT_LeptM, eta_LeptP, eta_LeptM, pT_b, pT_bbar, eta_b, eta_bbar
real(8) :: CosTheta_scatter,CosTheta_soper,CosTheta_star,Phi,Phibar,dPhiPlus,dPhiMinus,MassAux
real(8) :: MomHadr(1:4,1:7),MomLept(1:4,1:4),zeros(1:9),cosPhi,cosPhibar
real(8) :: MomTopsCMS(1:4,1:2),nx(2:4),ny(2:4),nz(2:4),MomLeptTRF(1:4),MomBeam(2:4)
integer :: i,j,n,k
real(8) :: sinPhi, sinPhibar, deltaSin, MomBeam1(1:4), MomBeam2(1:4), Mttb_cut,MTW
real(8) :: pT_bjet1, pT_bjet2, pT_lep1, pT_lep2, eta_bjet1, eta_bjet2, eta_lep1, eta_lep2,R_LepJet, eta_jet1, eta_jet2
real(8) :: MomLeptOrd(1:4,1:2), MomJetOrd(1:4,1:8), pT_miss, ET_miss, pT_jet1, pT_jet2, M_eff, HT, Minv_Lept,MomFatJet(1:4)

!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if( .not.NPlus1PS ) then
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   else
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5)
        elseif( TopDecays.ne.0 .and. Correction.eq.2 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        elseif( TopDecays.ne.0 .and. Correction.eq.5 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6) - MomDK(1:4,7)
        endif
   endif
   if( any(abs(zeros(1:4)/MomExt(1,1)).gt.1d-6) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:4)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( .not. NPlus1PS ) then
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3) = 0d0
   elseif( NPlus1PS .and. Correction.eq.2 ) then 
        zeros(1) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(2) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
        zeros(3)=  MomExt(1:4,3).dot.MomExt(1:4,3)
   elseif( NPlus1PS .and. Correction.eq.5 ) then 
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3)=  MomDK(1:4,7).dot.MomDK(1:4,7)
   endif
   zeros(4) =  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(5) =  MomDK(1:4,2).dot.MomDK(1:4,2)
   zeros(6) =  MomDK(1:4,3).dot.MomDK(1:4,3)
   zeros(7) =  MomDK(1:4,4).dot.MomDK(1:4,4)
   zeros(8) =  MomDK(1:4,5).dot.MomDK(1:4,5)
   zeros(9) =  MomDK(1:4,6).dot.MomDK(1:4,6)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_TTBARZprime(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:7)
   endif
!DEC$ ENDIF
 


! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light particles, N-1=ATop, N=Top
!                          MomDK(:,1:7) : 1=ABot, 2=lep-/q, 3=ANeu/qbar, 4=Bot, 5=lep+/qbar, 6=Neu/q, 7=(Glu)
! MomLept(:,1:4): 1=lep-, 2=ANeu, 3=lep+, 4=Neu


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)




! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   NumHadr = 0
   if( NPlus1PS ) then
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   else
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
   endif
!-------------------------------------------------------
elseif( TopDecays.eq.1 ) then  ! full leptonic decay
  MomLept(1:4,1) = MomDK(1:4,2)   ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)   ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)   ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)   ! Neu

  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,3) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,3) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,3) = MomDK(1:4,2)   ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)   ! q
  MomHadr(1:4,5) = MomDK(1:4,5)   ! qbar
  MomHadr(1:4,6) = MomDK(1:4,6)   ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,7) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,7) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
endif





if( ObsSet.eq.60 ) then! set of observables for ttb production without decays

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    if (M_TTbar.lt.Mttb_cut) then
        applyPSCut = .true.
        RETURN
    endif


    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))

!-------------------------------------------------------
elseif( ObsSet.eq.61 ) then! set of observables for ttb production with di-lept. decays

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------



   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

   call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:8) = MomJet(1:4,1:8)
   call pT_order(NumHadr,MomJetOrd(1:4,1:NumHadr))


!   M_LL = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--F Phase space cuts

!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

    if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi
!pause

 else
    dPhiPlus = 20d0
    dPhiMinus = 20d0
 endif

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)

!-------------------------------------------------------
elseif( ObsSet.eq.62 ) then! set of observables for ttb production with hadr. Atop, hadr. top decay

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------


    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)



elseif( ObsSet.eq.64 ) then ! Zprime, semi-hadronic top decay (for factorization checks)

   NBin(:) = 1

elseif( ObsSet.eq.65 ) then ! Zprime, semi-hadronic top decay (for ATLAS analysis: James Ferrando)


!---------------------- (anti) kT jet algorithm ---------------------------------

    Rsep_Jet = 1.0d0

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

    if( NJet.lt.1 ) then
       applypscut = .true.
       RETURN
    endif
    pT_jet1 = get_PT(MomJet(1:4,1))
    if( pT_jet1.lt.pT_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_jet1 = get_eta(MomJet(1:4,1))
    if( abs(eta_jet1).gt.eta_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 

    MomFatJet(1:4) = MomJet(1:4,1)



!---------------------- (anti) kT jet algorithm ---------------------------------

    Rsep_Jet = 0.4d0

    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

   k=0
   do i=1,NJet
       if( get_R(MomJet(1:4,i),MomFatJet(1:4)).gt.1.5d0 .and. get_PT(MomJet(1:4,i)).gt.25d0*GeV .and. abs(get_eta(MomJet(1:4,i))).lt.2.5d0 ) then
            k=k+1
       endif
   enddo
   if( k.lt.1 ) then
       applypscut = .true.
       RETURN
    endif


    pT_lep1 = get_PT(MomLept(1:4,3))
    if( pT_lep1.lt.pT_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_lep1 = get_ETA(MomLept(1:4,3))
    if( abs(eta_lep1).lt.eta_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 

    pT_miss = get_pT(MomLept(1:4,4))
    if( pT_miss.lt.pT_miss_cut ) then
       applypscut = .true.
       RETURN
    endif 

    MTW = Get_MT(MomLept(1:4,3),MomLept(1:4,4))
    if( MTW.lt.MTW_cut ) then
       applypscut = .true.
       RETURN
    endif 

    do i=1,NJet
      R_LepJet = get_R(MomJet(1:4,i),MomLept(1:4,3))
      if( R_LepJet.lt.Rsep_LepJet ) then
         applypscut = .true.
         RETURN
      endif 
    enddo






!   ------------------------
!   this is just for binning
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)





elseif( ObsSet.eq.66 ) then ! Zprime, semi-hadronic top decay (for CMS analysis: Roman Kogler)


!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(NumHadr,MomJet(1:4,1:NumHadr))! pT-order jets
!--------------------------------------------------------------------------

   if( NJet.lt.2 ) then
       applypscut = .true.
       RETURN
    endif

    pT_jet1 = get_PT(MomJet(1:4,1))
    pT_jet2 = get_PT(MomJet(1:4,2))
    if( pT_jet1.lt.pT_jet_cut ) then
       applypscut = .true.
       RETURN
    endif 
    if( pT_jet2.lt.50d0*GeV ) then
       applypscut = .true.
       RETURN
    endif
   

    eta_jet1 = get_eta(MomJet(1:4,1))
    eta_jet2 = get_eta(MomJet(1:4,2))
    if (dabs(eta_jet1).gt.eta_jet_cut .or. dabs(eta_jet2).gt.eta_jet_cut) then
       applypscut = .true.
       RETURN
    endif

    pT_lep1 = get_PT(MomLept(1:4,3))
    if( pT_lep1.lt.pT_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 
    eta_lep1 = get_ETA(MomLept(1:4,3))
    if( abs(eta_lep1).lt.eta_lep_cut ) then
       applypscut = .true.
       RETURN
    endif 

    pT_miss = get_pT(MomLept(1:4,4))
    if( pT_miss.lt.pT_miss_cut ) then
       applypscut = .true.
       RETURN
    endif 

    HT = pT_lep1 + pT_miss
    if( HT.lt.HT_cut ) then
       applypscut = .true.
       RETURN
    endif 

    do i=1,NJet
      R_LepJet = get_R(MomJet(1:4,i),MomLept(1:4,3))
      if( R_LepJet.lt.Rsep_LepJet ) then
         applypscut = .true.
         RETURN
      endif 
    enddo

!   --------------------------------
!   this is just for binning
    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)










elseif (ObsSet.EQ.67) then ! set of observables for ttb production without decays, for SM Z


    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))



    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))


endif

return
END SUBROUTINE







SUBROUTINE Kinematics_eeTTBAR(NPlus1PS,MomExt,MomDK,applyPSCut,NBin)
use ModMisc
use ModParameters
implicit none
integer :: NumHadr
real(8) :: MomExt(:,:),MomDK(:,:),MomJet(1:4,1:8),MomAux1(1:4),MomAux2(1:4)
real(8) :: MomTops(1:4,1:2),MomBoost(1:4)
logical :: applyPSCut,NPlus1PS
integer :: NBin(:),PartList(1:7),JetList(1:7),NJet,NObsJet_Tree,NObsJet
real(8) :: y_Top,y_ATop,pT_Top,pT_ATop,M_TTbar,Dphi_TTbar,pT_Lept,y_Lept,Dphi_LL,M_LL
real(8) :: pT_LeptP, pT_LeptM, eta_LeptP, eta_LeptM, pT_b, pT_bbar, eta_b, eta_bbar
real(8) :: CosTheta_scatter,CosTheta_soper,CosTheta_star,Phi,Phibar,dPhiPlus,dPhiMinus,MassAux
real(8) :: MomHadr(1:4,1:7),MomLept(1:4,1:4),zeros(1:9),cosPhi,cosPhibar
real(8) :: MomTopsCMS(1:4,1:2),nx(2:4),ny(2:4),nz(2:4),MomLeptTRF(1:4),MomBeam(2:4)
integer :: i,j,n,k
real(8) :: sinPhi, sinPhibar, deltaSin, MomBeam1(1:4), MomBeam2(1:4), Mttb_cut,MTW
real(8) :: pT_bjet1, pT_bjet2, pT_lep1, pT_lep2, eta_bjet1, eta_bjet2, eta_lep1, eta_lep2,R_LepJet, eta_jet1, eta_jet2
real(8) :: MomLeptOrd(1:4,1:2), MomJetOrd(1:4,1:8), pT_miss, ET_miss, pT_jet1, pT_jet2, M_eff, HT, Minv_Lept,MomFatJet(1:4)

!DEC$ IF(_CheckMomenta .EQ.1)
   zeros(:) = 0d0
   if( .not.NPlus1PS ) then
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4)
        else
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        endif
   else
        if(TopDecays.eq.0) then
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomExt(1:4,4) - MomExt(1:4,5)
        elseif( TopDecays.ne.0 .and. Correction.eq.2 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomExt(1:4,3) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6)
        elseif( TopDecays.ne.0 .and. Correction.eq.5 ) then 
            zeros(1:4) = MomExt(1:4,1)+MomExt(1:4,2) - MomDK(1:4,1)-MomDK(1:4,2)-MomDK(1:4,3)-MomDK(1:4,4)-MomDK(1:4,5)-MomDK(1:4,6) - MomDK(1:4,7)
        endif
   endif
   if( any(abs(zeros(1:4)/MomExt(1,1)).gt.1d-6) ) then
      print *, "ERROR: energy-momentum violation in SUBROUTINE Kinematics_eeTTBAR(",NPlus1PS,"): ",zeros(1:4)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:5+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( .not. NPlus1PS ) then
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3) = 0d0
   elseif( NPlus1PS .and. Correction.eq.2 ) then 
        zeros(1) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(2) = (MomExt(1:4,5).dot.MomExt(1:4,5)) - m_Top**2
        zeros(3)=  MomExt(1:4,3).dot.MomExt(1:4,3)
   elseif( NPlus1PS .and. Correction.eq.5 ) then 
        zeros(1) = (MomExt(1:4,3).dot.MomExt(1:4,3)) - m_Top**2
        zeros(2) = (MomExt(1:4,4).dot.MomExt(1:4,4)) - m_Top**2
        zeros(3)=  MomDK(1:4,7).dot.MomDK(1:4,7)
   endif
   zeros(4) =  MomDK(1:4,1).dot.MomDK(1:4,1)
   zeros(5) =  MomDK(1:4,2).dot.MomDK(1:4,2)
   zeros(6) =  MomDK(1:4,3).dot.MomDK(1:4,3)
   zeros(7) =  MomDK(1:4,4).dot.MomDK(1:4,4)
   zeros(8) =  MomDK(1:4,5).dot.MomDK(1:4,5)
   zeros(9) =  MomDK(1:4,6).dot.MomDK(1:4,6)
   if( TopDecays.eq.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_eeTTBAR(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:6)
   endif
   if( TopDecays.ne.0 .and. any(abs(zeros(1:9)/MomExt(1,1)**2).gt.1d-6) ) then
      print *, "ERROR: onshell-ness violation in SUBROUTINE Kinematics_eeTTBAR(",NPlus1PS,"): ",zeros(1:9)
      print *, "momenta dump:"
      print *, MomExt(1:4,1:4+NPlus1PS)
      print *, MomDK(1:4,1:7)
   endif
!DEC$ ENDIF
 


! required momentum order: MomExt(:,:): 1=In_left, 2=In_right, 3,4,...=light particles, N-1=ATop, N=Top
!                          MomDK(:,1:7) : 1=ABot, 2=lep-/q, 3=ANeu/qbar, 4=Bot, 5=lep+/qbar, 6=Neu/q, 7=(Glu)
! MomLept(:,1:4): 1=lep-, 2=ANeu, 3=lep+, 4=Neu


applyPSCut = .false.
NBin(1:NumHistograms) = 0


MomHadr(1:4,1:7) = 0d0
MomLept(1:4,1:4) = 0d0
PartList(1:7)=(/1,2,3,4,5,6,7/)




! separating momenta into hadron momenta and lepton momenta to which various procedures can be applied
MomHadr(1:4,1) = MomDK(1:4,1)  ! ABot
MomHadr(1:4,2) = MomDK(1:4,4)  ! Bot
!-------------------------------------------------------
if( TopDecays.eq.0 ) then  ! no top decays
   NumHadr = 0
   if( NPlus1PS ) then
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
   else
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
   endif
!-------------------------------------------------------
elseif( TopDecays.eq.1 ) then  ! full leptonic decay
  MomLept(1:4,1) = MomDK(1:4,2)   ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)   ! ANeu
  MomLept(1:4,3) = MomDK(1:4,5)   ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)   ! Neu

  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,3) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,3) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 3
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 2
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.2 ) then  ! full hadronic decay
  MomHadr(1:4,3) = MomDK(1:4,2)   ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)   ! q
  MomHadr(1:4,5) = MomDK(1:4,5)   ! qbar
  MomHadr(1:4,6) = MomDK(1:4,6)   ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,7) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,7) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 7
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 6
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.3 ) then  ! lept. Atop, hadr. top decay
  MomLept(1:4,1) = MomDK(1:4,2)  ! lep-
  MomLept(1:4,2) = MomDK(1:4,3)  ! ANeu
  MomHadr(1:4,3) = MomDK(1:4,5)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,6)  ! q
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
!-------------------------------------------------------
elseif( TopDecays.eq.4 ) then  ! hadr. Atop, lept. top decay
  MomHadr(1:4,3) = MomDK(1:4,2)  ! qbar
  MomHadr(1:4,4) = MomDK(1:4,3)  ! q
  MomLept(1:4,3) = MomDK(1:4,5)  ! lep+
  MomLept(1:4,4) = MomDK(1:4,6)  ! Neu
  if( NPlus1PS .and. CORRECTION.eq.2 ) then
      MomHadr(1:4,5) = MomExt(1:4,3)   ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,4)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,5)   ! Top
  elseif( NPlus1PS .and. CORRECTION.eq.5 ) then
      MomHadr(1:4,5) = MomDK(1:4,7)    ! q/qbar/glu
      NumHadr = 5
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  else
      NumHadr = 4
      MomTops(1:4,1) = MomExt(1:4,3)   ! ATop
      MomTops(1:4,2) = MomExt(1:4,4)   ! Top
  endif
endif





if( ObsSet.eq.70 ) then! set of observables for ttb production without decays

    pT_ATop = get_PT(MomTops(1:4,1))
    pT_Top  = get_PT(MomTops(1:4,2))

    y_ATop = get_ETA(MomTops(1:4,1))
    y_Top  = get_ETA(MomTops(1:4,2))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    if (M_TTbar.lt.Mttb_cut) then
        applyPSCut = .true.
        RETURN
    endif


    Dphi_TTbar = dabs( Get_PHI(MomTops(1:4,1)) - Get_PHI(MomTops(1:4,2))  )
    if( Dphi_TTbar.gt.Pi ) Dphi_TTbar=2d0*Pi-Dphi_TTbar

!   scattering angle of the top quark in lab frame
    CosTheta_scatter = get_CosTheta( MomTops(1:4,2) )

!   see chinese paper, approximates scattering angle in ttb rest frame
    MomAux1(1:4) = MomTops(1:4,2)
    MomAux2(1:4) = MomExt(1:4,1)
    call boost(MomAux1(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)! this seems wrong!!!
    call boost(MomAux2(1:4),MomTops(1:4,1)+MomTops(1:4,2),m_top)
    CosTheta_star = get_CosAlpha(MomAux1(1:4),MomAux2(1:4))



! binning
    NBin(1) = WhichBin(1,pT_ATop)
    NBin(2) = WhichBin(2,pT_Top)
    NBin(3) = WhichBin(3,y_ATop)
    NBin(4) = WhichBin(4,y_Top)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,Dphi_TTbar)
    NBin(7) = WhichBin(7,CosTheta_scatter)
    NBin(8) = WhichBin(8,CosTheta_star)

    if( NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)+MomExt(1:4,3)))
    if( .not. NPlus1PS ) NBin(9) = WhichBin(9,get_MInv(MomTops(1:4,1)+MomTops(1:4,2)))

!-------------------------------------------------------
elseif( ObsSet.eq.71 ) then! set of observables for ttb production with di-lept. decays

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------



   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

!    call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:8) = MomJet(1:4,1:8)
   call pT_order(NumHadr,MomJetOrd(1:4,1:NumHadr))



   
!   check that there are two b jets    
    if( .not.(any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

!     if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi

!  else
!     dPhiPlus = 20d0
!     dPhiMinus = 20d0
!  endif


!   scattering angles for comparisons
    CosTheta_scatter = VectorProd(MomTops(2:4,2),MomExt(2:4,1))/MomTops(1,2)/MomExt(1,1) ! checked against literature (Peskin), maximum at cos()=+1, minimum at cos()=-1
!     CosTheta_scatter = (MomDK(1:4,4).dot.MomDK(1:4,5))*4d0/(m_top**2-m_w**2)-1d0! checked against literature (Peskin): zero at cos()=+1, non-zero at cos()=-1
    
    

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)
    NBin(14) = WhichBin(14,CosTheta_scatter)

    
    NBin(15) = WhichBin(15,dPhiMinus)
    NBin(16) = WhichBin(16,dPhiPlus)
    NBin(17) = Histo(15)%NBins * ( NBin(16) -1  ) + NBin(15) 
    
    
!-------------------------------------------------------
elseif( ObsSet.eq.72 ) then! set of observables for ttb production with di-lept. decays

!---------------------- (anti) kT jet algorithm ---------------------------------
    NJet=0
    MomJet(1:4,1:7) = MomHadr(1:4,1:7)
    call JetAlgo_kt(Rsep_jet,PartList(1:NumHadr),MomHadr(1:4,1:NumHadr),NJet,JetList(1:NumHadr),MomJet(1:4,1:NumHadr))
    call pT_order(2,MomJet(1:4,1:2))! pT-order b-jets first
    call pT_order(NumHadr-2,MomJet(1:4,3:NumHadr))! pT-order non-b jets
!--------------------------------------------------------------------------


   !-- define ordering for leptons, without modifying existing lepton array used below for angles
   MomLeptOrd(1:4,1) = MomLept(1:4,1)
   MomLeptOrd(1:4,2) = MomLept(1:4,3)

!    call pT_order(2,MomLeptOrd)! pT-order for leptons

   !-- order all jets, irrespective whether they are b-jets or not
   MomJetOrd(1:4,1:2) = MomJet(1:4,1:2)
   call pT_order(2,MomJetOrd(1:4,1:2))
   MomJetOrd(1:4,3:NJet) = MomJet(1:4,3:NJet)
   call pT_order(NJet-2,MomJetOrd(1:4,3:NJet))



   
!   check that there are two b jets    
    NObsJet_Tree = 4
    if( .not.(NJet.ge.NObsJet_Tree .and. any(JetList(1:NJet).eq.1) .and. any(JetList(1:NJet).eq.2)) ) then
        applyPSCut = .true.
        RETURN
    endif

    pT_bjet2 = get_PT(MomJet(1:4,2))
    pT_lep2  = get_PT(MomLeptOrd(1:4,2))

    eta_bjet1 = get_eta(MomJet(1:4,1))
    eta_bjet2 = get_eta(MomJet(1:4,2))

    eta_lep1 = get_eta(MomLept(1:4,1))
    eta_lep2 = get_eta(MomLeptOrd(1:4,2))

    y_Top  = get_ETA(MomTops(1:4,2))

    pT_jet2 = get_PT(MomJet(1:4,4))
    eta_jet1 = get_eta(MomJet(1:4,3))
    eta_jet2 = get_eta(MomJet(1:4,4))

    if (pT_lep2.lt.pt_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_lep1).gt.eta_lep_cut .or. dabs(eta_lep2).gt.eta_lep_cut) then
       applypscut = .true.
       RETURN
    endif

    if (pT_bjet2.lt.pt_bjet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_bjet1).gt.eta_bjet_cut .or. dabs(eta_bjet2).gt.eta_bjet_cut) then
       applypscut = .true.
       RETURN
    endif


    if (pT_jet2.lt.pt_jet_cut) then
       applypscut = .true.
       RETURN
    endif

    if (dabs(eta_jet1).gt.eta_jet_cut .or. dabs(eta_jet2).gt.eta_jet_cut) then
       applypscut = .true.
       RETURN
    endif


!-- If we pass the cuts, compute other observables
    pT_bjet1 = get_PT(MomJet(1:4,1))
    pT_lep1 = get_PT(MomLeptOrd(1:4,1))

    M_TTbar = get_MInv(MomTops(1:4,1)+MomTops(1:4,2))

    Dphi_LL = dabs( Get_PHI(MomLept(1:4,1)) - Get_PHI(MomLept(1:4,3))  )
    if( Dphi_LL.gt.Pi ) Dphi_LL=2d0*Pi-Dphi_LL

    pT_miss = get_pT(MomLept(1:4,2) + MomLept(1:4,4))
    Minv_Lept = get_MInv(MomLept(1:4,1)+MomLept(1:4,3))
    ET_miss = dsqrt(pT_miss**2 + Minv_Lept**2)

    pT_jet1 = get_PT(MomJetOrd(1:4,1))
    pT_jet2 = get_PT(MomJetOrd(1:4,2))

    HT = pT_lep1 + pT_lep2 + pT_jet1 + pT_jet2

    M_eff = ET_miss + HT

!   construct Baumgart-Tweedie angles // for non-zero transverse momentum of TTBAR, use Collins-Soper construction

!     if (dabs(M_TTbar-M_Zpr).le.(100d0*GeV)) then ! compute angles only close to the resonance

!   step 1: boost top momenta into a TTBAR CMS frame ( pT+pTbar=0 )
    MomAux1(1:4) = MomTops(1:4,1)+MomTops(1:4,2)
    MomAux1(2:4) =-MomAux1(2:4)
    MassAux = dsqrt( (MomAux1(1:4).dot.MomAux1(1:4)) +1d-12 )
    MomTopsCMS(1:4,1:2) = MomTops(1:4,1:2)
    call boost(MomTopsCMS(1:4,1),MomAux1(1:4),MassAux)
    call boost(MomTopsCMS(1:4,2),MomAux1(1:4),MassAux)

!   step 2: construct coordinate system in a TTBAR CMS frame

    MomBeam1(1:4) = (/1d0,0d0,0d0,1d0/)
    call boost(MomBeam1(1:4),MomAux1(1:4),MassAux)

    MomBeam2(1:4) = (/1d0,0d0,0d0,-1d0/)
    call boost(MomBeam2(1:4),MomAux1(1:4),MassAux)

    MomBeam(2:4) = MomBeam1(2:4)-MomBeam2(2:4)
    MomBeam(2:4) = MomBeam(2:4)/dsqrt(MomBeam(2)**2+MomBeam(3)**2+MomBeam(4)**2)     ! Collins-Soper direction

    ny(2:4) = MomTopsCMS(2:4,2).cross.MomBeam(2:4)
    ny(2:4) = ny(2:4)/dsqrt(ny(2)**2+ny(3)**2+ny(4)**2)
    nz(2:4) = MomTopsCMS(2:4,2)/dsqrt(MomTopsCMS(2,2)**2+MomTopsCMS(3,2)**2+MomTopsCMS(4,2)**2)
    nx(2:4) = ny(2:4).cross.nz(2:4)

!   step 3: boost lepton momenta first into CMS frame, then into frame with ptop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,3)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)
    MomAux2(1:4) = MomTopsCMS(1:4,2)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with nx

    cosPhi = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhi = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2  & 
    + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )


!    Phi = dacos(cosPhi)

!------- repeat the same for the anti-top

!   step 3: boost lepton momenta into frame with patop=0 and call it MomLeptTRF
    MomLeptTRF(1:4) = MomLept(1:4,1)
    call boost(MomLeptTRF(1:4),MomAux1(1:4),MassAux)

    MomAux2(1:4) = MomTopsCMS(1:4,1)
    MomAux2(2:4) =-MomAux2(2:4)
    call boost(MomLeptTRF(1:4),MomAux2(1:4),m_top)

!   step 4: project MomLeptTRF onto the nx-ny plane and calculate its angle with n
    cosPhibar = VectorProd(MomLeptTRF(2:4),nx(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

    sinPhibar = VectorProd(MomLeptTRF(2:4),ny(2:4))/dsqrt( VectorProd(MomLeptTRF(2:4),nx(2:4))**2 & 
      + VectorProd(MomLeptTRF(2:4),ny(2:4))**2 )

!    Phibar = dacos(cosPhibar)

    deltaSin = sinPhi * cosPhibar - cosPhi * sinPhibar

    if (deltaSin .gt. 0d0) then
       dPhiMinus = dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    else
       dPhiMinus = -dacos( cosPhi*cosPhibar + sinPhi*sinPhibar)
    endif

    dPhiPlus = dabs(dacos(cosPhi*cosPhibar - sinPhi*sinPhibar))

!    dPhiMinus = Phi-Phibar          ! -pi..pi
!    dPhiPlus  = dabs(Phi+Phibar)    ! 0...2pi

!    print *, phi/DblPi,phibar/DblPi
!    print *, dPhiMinus/DblPi,dPhiPlus/DblPi

!  else
!     dPhiPlus = 20d0
!     dPhiMinus = 20d0
!  endif


!   scattering angles for comparisons
    CosTheta_scatter = VectorProd(MomTops(2:4,2),MomExt(2:4,1))/MomTops(1,2)/MomExt(1,1) ! checked against literature (Peskin), maximum at cos()=+1, minimum at cos()=-1
!     CosTheta_scatter = (MomDK(1:4,4).dot.MomDK(1:4,5))*4d0/(m_top**2-m_w**2)-1d0! checked against literature (Peskin): zero at cos()=+1, non-zero at cos()=-1
    
    

! binning
    NBin(1) = WhichBin(1,pT_lep1)
    NBin(2) = WhichBin(2,pT_bjet1)
    NBin(3) = WhichBin(3,pT_lep2)
    NBin(4) = WhichBin(4,pT_bjet2)
    NBin(5) = WhichBin(5,M_TTbar)
    NBin(6) = WhichBin(6,M_eff)
    NBin(7) = WhichBin(7,eta_lep1)
    NBin(8) = WhichBin(8,eta_bjet1)
    NBin(9) = WhichBin(9,eta_lep2)
    NBin(10) = WhichBin(10,eta_bjet2)
    NBin(11) = WhichBin(11,dPhi_LL)
    NBin(12) = WhichBin(12,dPhiMinus)
    NBin(13) = WhichBin(13,dPhiPlus)
    NBin(14) = WhichBin(14,CosTheta_scatter)
    NBin(18) = WhichBin(18,y_Top)
    NBin(19) = WhichBin(19,y_Top)
    
    NBin(15) = WhichBin(15,dPhiMinus)
    NBin(16) = WhichBin(16,dPhiPlus)
    NBin(17) = Histo(15)%NBins * ( NBin(16) -1  ) + NBin(15) 
    

    
endif

return
END SUBROUTINE















RECURSIVE SUBROUTINE JetAlgo_kt(Rsep_jet,PartonList,MomParton,NJet,JetList,MomJet)  ! initial call must have NJet=0 and MomJet(1:4,:) = MomPartons(1:4,:)
use ModMisc
use ModParameters
implicit none
integer :: PartonList(:), JetList(:)
real(8) :: MomParton(:,:)
integer :: NJet
real(8) :: MomJet(:,:)
real(8) :: Rsep_jet
integer :: NParton,i,j,k,dii_minp,dij_minp,ij(1:120,1:2)  ! max. partons=8, (8-1)! =5040  ! max. partons=6, (6-1)! =120
real(8) :: dii(1:6),dij(1:120),Rij,dii_min,dij_min!,eta(1:5),phi(1:5)


NParton = size(PartonList)
if( NParton.eq.0 ) then
    return
elseif( NParton.eq.1 ) then
!   print *, "HARD JET", PartonList(NParton)
   NJet = NJet +1
   JetList(NJet) = PartonList(NParton)
   return
endif

!generate dii, dij
do i=1,NParton
   dii(i) = get_PT2(MomParton(1:4,PartonList(i)))**AlgoType
enddo

k=0
do i=1,NParton-1
do j=i+1,NParton
   k = k+1
   Rij = get_R( MomParton(1:4,PartonList(i)), MomParton(1:4,PartonList(j)) )
   dij(k) = dmin1(dii(i),dii(j)) * (Rij/Rsep_jet)**2
   ij(k,1)=i
   ij(k,2)=j
enddo
enddo
!!print *, dii(1:NParton)
!!print *, dij(1:k)


! find minima
dii_min=dii(1)
dii_minp=1
do i=2,NParton
  if( dii(i).lt.dii_min ) then
    dii_min = dii(i)
    dii_minp= i
  endif
enddo

dij_min=dij(1)
dij_minp=1
do i=2,k
  if( dij(i).lt.dij_min ) then
    dij_min = dij(i)
    dij_minp= i
  endif
enddo


if( dii_min.lt.dij_min ) then ! no recombination
   NJet = NJet +1
   k=0
!   print *, "HARD JET", PartonList(dii_minp)
   JetList(NJet) = PartonList(dii_minp)
   MomJet(1:4,PartonList(dii_minp)) = MomParton(1:4,PartonList(dii_minp))
   do i=1,NParton!  remove momenta dii_minp from parton list
      if( i.eq.dii_minp ) cycle
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomParton(:,:),NJet,JetList(:),MomJet(:,:))
else ! recombination of dij(dij_min)
   if( RecombPrescr.eq.1 ) then
        MomJet(1:4,PartonList(ij(dij_minp,1))) = EllisSoperComb(MomJet(1:4,PartonList(ij(dij_minp,1))),MomJet(1:4,PartonList(ij(dij_minp,2))))   ! Ellis-Soper combination
   else
        MomJet(1:4,PartonList(ij(dij_minp,1))) = MomJet(1:4,PartonList(ij(dij_minp,1))) + MomJet(1:4,PartonList(ij(dij_minp,2)))                 ! four vector addition
   endif
   MomJet(1:4,PartonList(ij(dij_minp,2))) = 0d0
   k=0
!   print *, "RECOMB.",PartonList(ij(dij_minp,1)),PartonList(ij(dij_minp,2))
   do i=1,NParton!  remove momenta j from parton list
      if( i.eq.ij(dij_minp,2) ) then
         cycle
      endif
      k=k+1
      PartonList(k) = PartonList(i)
   enddo
   call JetAlgo_kt(Rsep_jet,PartonList(1:k),MomJet(:,:),NJet,JetList(:),MomJet(:,:))
endif

return
END SUBROUTINE




FUNCTION EllisSoperComb(pi,pj)
use ModMisc
implicit none
real(8), intent(in) :: pi(1:4),pj(1:4)
real(8) :: EllisSoperComb(1:4)
real(8) :: ET_i,ET_j,ET_k,eta_i,eta_j,eta_k,phi_i,phi_j,phi_k,theta_k,E_k

   ET_i = get_PT(pi(1:4))   ! ET=pT as defined in hep-ph/9305266; Ellis,Soper
   ET_j = get_PT(pj(1:4))
   eta_i = get_PseudoEta(pi(1:4))  ! check pseudo_eta = eta for m=0
   eta_j = get_PseudoEta(pj(1:4))


   phi_i = get_Phi(pi(1:4))
   phi_j = get_Phi(pj(1:4))

   ET_k = ET_i + ET_j
   eta_k = (ET_i*eta_i + ET_j*eta_j)/ET_k
   phi_k = (ET_i*phi_i + ET_j*phi_j)/ET_k
   theta_k = 2d0*datan( dexp(-eta_k) )   ! datan( (0,inf) ) = (0,pi/2)  --> theta_k = (0,pi)
   E_k = ET_k/dsin(theta_k)

   EllisSoperComb(1) = E_k
   EllisSoperComb(2) = E_k * ( dsin(theta_k) * dcos(phi_k) )     ! check E_k prefactor
   EllisSoperComb(3) = E_k * ( dsin(theta_k) * dsin(phi_k) )
   EllisSoperComb(4) = E_k * ( dcos(theta_k) )

return
END FUNCTION




FUNCTION FrixioneIsolated(MomPho,Riso,NumParton,MomParton)! Frixione photon isolation cut, arXiv:hep-ph/9801442
use ModMisc
implicit none
logical :: FrixioneIsolated
integer :: NumParton,i,j
real(8) :: MomParton(1:4,1:NumParton),MomPho(1:4),Riso
real(8) :: Chi,Esum,OneMinusCosRip,OneMinusCosRiso,Rip,Rjp

  FrixioneIsolated = .true.
  do i=1,NumParton
      Rip = dmin1(get_R(MomParton(1:4,i),MomPho(1:4)),Riso)
      OneMinusCosRip  = Rip**2/2d0*(1d0 - Rip**2/12d0 + Rip**4/360d0 - Rip**6/20160d0)
      OneMinusCosRiso = Riso**2/2d0*(1d0 - Riso**2/12d0 + Riso**4/360d0 - Riso**6/20160d0)
      Chi = get_ET(MomPho(1:4))*(OneMinusCosRip)/(OneMinusCosRiso)

      Esum = 0d0
      do j=1,NumParton
         Rjp = get_R(MomParton(1:4,j),MomPho(1:4))
         if( Rjp.gt.Rip ) cycle
         if( Rjp.gt.Riso) cycle
         Esum = Esum + get_ET(MomParton(1:4,j))
      enddo

      if( Chi.lt.Esum ) then
         FrixioneIsolated = .false.
         return
      endif
  enddo


!12345 continue
!if(  IsolateFrix(MomPho,Riso,NumParton,MomParton(1:4,1:NumParton)).eq.FrixioneIsolated   ) then
!  print *, "Andreas", .not.IsolateFrix(MomPho,Riso,NumParton,MomParton(1:4,1:NumParton))
!  print *, "Markus ",FrixioneIsolated
!  print *, Chi,Esum
!  pause
!endif

return
END FUNCTION




function IsolateFrix(kpho,delta0,numHad,MomHad)! Andreas' implementation
use ModMisc
  implicit none
  integer numHad, i, numH
  real(8) :: MomHad(4,numHad), kpho(4), delta, ET(numHad),RiPho(numHad)
  real(8) :: delta0, LHS, RHS
  logical :: IsolateFrix
  IsolateFrix = .false.

  do i=0,1000
    delta = dble(i)/1000.d0*delta0
    LHS =0.d0
    do numH=1,numHad
        ET(numH) = dsqrt(MomHad(2,numH)**2 + MomHad(3,numH)**2)
        RiPho(numH)  = get_R(MomHad(1:4,numH),kpho)
        LHS = LHS + ET(numH)*StepFunc(delta-RiPho(numH))
    enddo
    RHS = dsqrt(kpho(2)**2 + kpho(3)**2)*((1.d0-dcos(delta))/(1.d0-dcos(delta0)))
    if(LHS .gt. RHS) then
        IsolateFrix = .true.
        return
    endif
  enddo

  return
end function







FUNCTION getHelicity(yRnd)
use ModParameters
use ModProcess
implicit none
integer :: getHelicity(1:2)
real(8) :: yRnd

  IF(HelSampling) THEN
    getHelicity(1:2)=yRnd*dble(NumHelicities) + 1
  ELSE
    getHelicity(1)=1
    getHelicity(2)=NumHelicities
  ENDIF

return
END FUNCTION



FUNCTION MomCrossing(Mom)
use ModProcess
use ModParameters
implicit none
real(8) :: Mom(1:4,1:NumExtParticles),MomCrossing
integer :: NPart
!real(8) :: QuarkCrossing=-1d0, SpinAvg=1d0/4d0, QuarkColAvg=1d0/3d0, GluonColAvg=1d0/8d0

do NPart=1,NumExtParticles
   ExtParticle(NPart)%Mom(1:4) = sign(1,Crossing(NPart)) * Mom(1:4,abs(Crossing(NPart)))
   ExtParticle(NPart)%Mom(5:8) = (0d0,0d0)
enddo
MomCrossing = AvgFactor

return
END FUNCTION





SUBROUTINE HelCrossing(Hel)
use ModProcess
use ModParameters
implicit none
integer :: Hel(1:NumExtParticles),n


do n=1,NumExtParticles
  ExtParticle(n)%Helicity = Hel(n)
enddo


return
END SUBROUTINE








SUBROUTINE CheckSing(MomExt,applySingCut)
use ModMisc
use ModParameters
use ModProcess
implicit none
logical :: applySingCut
real(8) :: MomExt(:,:),s12,s13,s23,s14,s24,s34,E3,E4


!     applySingCut=.false.
!     if( dsqrt(MomExt(2,3)**2 + MomExt(3,3)**2).lt.1d0 ) then
!         applySingCut=.true.
!         return
!     endif

IF( NumExtParticles.EQ.5 .AND. CORRECTION.EQ.2 ) THEN
    s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
    s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
    s23 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
    E3 = MomExt(1,3)
    applySingCut=.false.
    if( dabs(s13/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s23/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E3/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif
ELSEIF( NumExtParticles.EQ.6 .AND. CORRECTION.EQ.2 ) THEN
    s12 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,2))
    s13 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,3))
    s23 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,3))
    s14 = 2d0*(MomExt(1:4,1).dot.MomExt(1:4,4))
    s24 = 2d0*(MomExt(1:4,2).dot.MomExt(1:4,4))
    s34 = 2d0*(MomExt(1:4,3).dot.MomExt(1:4,4))
    E3 = MomExt(1,3)
    E4 = MomExt(1,4)
    applySingCut=.false.
    if( dabs(s13/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s23/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s14/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s24/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(s34/s12).lt.1d-9 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E3/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif
    if( dabs(E4/(2d0*MomExt(1,1))).lt.1d-5 ) then
        applySingCut=.true.
        return
    endif

ELSEIF( CORRECTION.EQ.5 ) THEN
    s14 = MomExt(1:4,1).dot.MomExt(1:4,4)
    E4 = dabs(MomExt(1,4))
    applySingCut=.false.
    if( dabs(s14)/m_Top**2 .lt. 1d-8 ) then
!     if( dabs(s14)/m_Top**2 .lt. 1d-4 ) then
        applySingCut=.true.
        return
    endif
    if( E4/m_Top .lt. 1d-4 ) then
!     if( E4/m_Top .lt. 1d-2 ) then
        applySingCut=.true.
        return
    endif
ENDIF

return
END SUBROUTINE






FUNCTION WhichBin(NHisto,Value)
implicit none
integer :: WhichBin,NHisto
real(8) :: Value

   WhichBin = (Value-Histo(NHisto)%LowVal)/Histo(NHisto)%BinSize + 1
   if( WhichBin.lt.0 ) then
      WhichBin = 0
   elseif( WhichBin.gt.Histo(NHisto)%NBins ) then
      WhichBin = Histo(NHisto)%NBins+1
   endif

RETURN
END FUNCTION







! SUBROUTINE IntoHisto(NHisto,NBin,Value)!    older version without compatibility for bin smearing
! use ModMisc
! implicit none
! integer :: NHisto,NBin
! real(8) :: Value
! 
!      if( IsNaN(Value) ) return
! 
! !DEC$ IF(_UseMPIVegas.EQ.0)
!      Histo(NHisto)%Value(NBin) = Histo(NHisto)%Value(NBin)  + Value
!      Histo(NHisto)%Value2(NBin)= Histo(NHisto)%Value2(NBin) + Value**2
!      Histo(NHisto)%Hits(NBin)  = Histo(NHisto)%Hits(NBin)+1
! !DEC$ ELSE
!      if( nbin.le.0 ) return
!      if( nbin.gt.MXHISTOBINS ) return
!      RedHisto(NHisto)%Hits(NBin)   = RedHisto(NHisto)%Hits(NBin)+1
!      RedHisto(NHisto)%Value(NBin)  = RedHisto(NHisto)%Value(NBin)+Value
!      RedHisto(NHisto)%Value2(NBin) = RedHisto(NHisto)%Value2(NBin)+Value**2
! !DEC$ ENDIF
! 
! RETURN
! END SUBROUTINE


SUBROUTINE IntoHisto(NHisto,NBin,Value,BinValue)
use ModMisc
implicit none
integer :: NHisto,NBin,NeighbBin,i
real(8) :: Value
real(8),optional :: BinValue
real(8) :: LowerBinValue,UpperBinValue,NeighbBinValue,ErrorFunct


    if( IsNaN(Value) ) return
    if( (.not. Histo(NHisto)%BinSmearing) .or. NBin.eq.0 .or. NBin.eq.Histo(NHisto)%NBins+1 ) then
!DEC$ IF(_UseMPIVegas.EQ.0)
        Histo(NHisto)%Value(NBin)  = Histo(NHisto)%Value(NBin)  + Value
        Histo(NHisto)%Value2(NBin) = Histo(NHisto)%Value2(NBin) + Value**2
        Histo(NHisto)%Hits(NBin)   = Histo(NHisto)%Hits(NBin)+1
!DEC$ ELSE
        if( nbin.le.0 ) return
        if( nbin.gt.MXHISTOBINS ) return
        RedHisto(NHisto)%Value(NBin)  = RedHisto(NHisto)%Value(NBin)  + Value
        RedHisto(NHisto)%Value2(NBin) = RedHisto(NHisto)%Value2(NBin) + Value**2
        RedHisto(NHisto)%Hits(NBin)   = RedHisto(NHisto)%Hits(NBin)+1
!DEC$ ENDIF
    else
        if( .not. present(BinValue) ) then 
             call Error("Argument BinValue is missing in call to IntoHisto.")
        else
             NBin = WhichBin(NHisto,BinValue)
        endif
        LowerBinValue=(NBin-1)*Histo(NHisto)%BinSize + Histo(NHisto)%LowVal
        UpperBinValue=LowerBinValue + Histo(NHisto)%BinSize
        if( BinValue.gt.LowerBinValue+Histo(NHisto)%BinSize/2d0 ) then
           NeighbBinValue=UpperBinValue
           NeighbBin=NBin+1
           if( NeighbBin.gt.Histo(NHisto)%NBins  ) NeighbBin=Histo(NHisto)%NBins
        else
           NeighbBinValue=LowerBinValue
           NeighbBin=NBin-1
           if( NeighbBin.le.0 ) NeighbBin=1
        endif

! print *, NBin,NeighbBin
! pause

           ! ErrorFunct ranges between 0...+1
           ! erf(0)=0, erf(1/sqrt2)=0.68, erf(2/sqrt2)=0.95, erf(3/sqrt2)=0.99
           ! --> the smaller SmearSigma the less leakage into the other bin
           ! --> SmearSigma=0.5 means that sigma of the gauss distribution is 0.5 of half the bin size
           ErrorFunct=erf( dabs(NeighbBinValue-BinValue)/(Histo(NHisto)%SmearSigma*0.5d0*Histo(NHisto)%BinSize)/dsqrt(2d0) )
!DEC$ IF(_UseMPIVegas.EQ.0)
           Histo(NHisto)%Value(NBin)      = Histo(NHisto)%Value(NBin)       + 0.5d0*(1d0+ErrorFunct)*Value
           Histo(NHisto)%Value(NeighbBin) = Histo(NHisto)%Value(NeighbBin)  + 0.5d0*(1d0-ErrorFunct)*Value
           Histo(NHisto)%Value2(NBin)     = Histo(NHisto)%Value2(NBin)      +(0.5d0*(1d0+ErrorFunct)*Value)**2
           Histo(NHisto)%Value2(NeighbBin)= Histo(NHisto)%Value2(NeighbBin) +(0.5d0*(1d0-ErrorFunct)*Value)**2
           Histo(NHisto)%Hits(NBin)       = Histo(NHisto)%Hits(NBin)+1
!DEC$ ELSE
           RedHisto(NHisto)%Value(NBin)      = RedHisto(NHisto)%Value(NBin)       + 0.5d0*(1d0+ErrorFunct)*Value
           RedHisto(NHisto)%Value(NeighbBin) = RedHisto(NHisto)%Value(NeighbBin)  + 0.5d0*(1d0-ErrorFunct)*Value
           RedHisto(NHisto)%Value2(NBin)     = RedHisto(NHisto)%Value2(NBin)      +(0.5d0*(1d0+ErrorFunct)*Value)**2
           RedHisto(NHisto)%Value2(NeighbBin)= RedHisto(NHisto)%Value2(NeighbBin) +(0.5d0*(1d0-ErrorFunct)*Value)**2
           RedHisto(NHisto)%Hits(NBin)       = RedHisto(NHisto)%Hits(NBin)+1
!DEC$ ENDIF
    endif

RETURN
END SUBROUTINE








SUBROUTINE setPolarizations(iSel,SecondPol)
use ModProcess
use ModMisc
use ModParameters
implicit none
integer :: NPart
complex(8),optional :: SecondPol(1:4)
integer,optional :: iSel


   do NPart=1,NumExtParticles
      if( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.lt.0 .and. abs(ExtParticle(NPart)%Helicity).eq.1 ) then
         call vSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
         if( present(SecondPol) .and. iSel.eq.NPart ) call vSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,-ExtParticle(NPart)%Helicity,SecondPol(1:4))

      elseif( IsAQuark(ExtParticle(NPart)%PartType) .and. ExtParticle(NPart)%PartType.gt.0 .and. abs(ExtParticle(NPart)%Helicity).eq.1 ) then
         call ubarSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
         if( present(SecondPol) .and. iSel.eq.NPart ) call ubarSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,-ExtParticle(NPart)%Helicity,SecondPol(1:4))

      elseif( ExtParticle(NPart)%PartType.eq.Glu_ ) then
         call pol_mless(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
         if( present(SecondPol) .and. iSel.eq.NPart ) call pol_mless(ExtParticle(NPart)%Mom(1:4),-ExtParticle(NPart)%Helicity,SecondPol(1:4))
!          if( present(SecondPol) .and. present(iSel) ) call pol_mless(ExtParticle(iSel)%Mom(1:4),-ExtParticle(iSel)%Helicity,SecondPol(1:4))

      elseif ( ExtParticle(NPart)%PartType.eq.Z0_ .and. .not. ZDecays.gt.0 ) then 
         call pol_massSR(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))! on-shell Z-boson polarization
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)

      elseif( ExtParticle(NPart)%PartType.eq.Pho_  ) then
         call pol_mless(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)

      elseif( ExtParticle(NPart)%PartType.eq.ElP_ ) then
         call vSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)

      elseif( ExtParticle(NPart)%PartType.eq.ElM_ ) then
         call ubarSpi(ExtParticle(NPart)%Mom(1:4),ExtParticle(NPart)%Mass,ExtParticle(NPart)%Helicity,ExtParticle(NPart)%Pol(1:4))
         ExtParticle(NPart)%Pol(5:16) = (0d0,0d0)
      endif
   enddo

!       check gauge invariance
!        ExtParticle(4)%Pol(1:4) = ExtParticle(4)%Mom(1:4);       print *, "gauge invariance check"

END SUBROUTINE



SUBROUTINE SetPropagators()
use ModProcess
implicit none
integer :: NPrimAmp,Prop

   do NPrimAmp=1,NumPrimAmps
         PrimAmps(NPrimAmp)%IntPart(1)%Mom(1:5) = (0d0,0d0)
         do Prop=2,NumExtParticles
            PrimAmps(NPrimAmp)%IntPart(Prop)%Mom(1:4) = PrimAmps(NPrimAmp)%IntPart(Prop-1)%Mom(1:4)  &
                                                      + dcmplx( ExtParticle( PrimAmps(NPrimAmp)%ExtLine(Prop-1) )%Mom(1:4) )
            PrimAmps(NPrimAmp)%IntPart(Prop)%Mom( 5 ) = (0d0,0d0)
         enddo
   enddo
END SUBROUTINE






SUBROUTINE boost2Lab(x1,x2,NumPart,Mom)
implicit none
real(8) Mom(1:4,1:NumPart)
real(8) x1,x2
real(8) gamma,betagamma,MomTmp1,MomTmp4
integer :: i,NumPart
!   beta  = (x2-x1)/(x1+x2)
!   gamma = 1d0/dsqrt(1d0-beta**2)
!   betagamma = beta*gamma

  gamma     = (x1+x2)/2d0/dsqrt(x1*x2)
  betagamma = (x2-x1)/2d0/dsqrt(x1*x2)

!   Mom(1,1)= x1*ColliderEnergy/2d0
!   Mom(2,1)= 0d0
!   Mom(3,1)= 0d0
!   Mom(4,1)= x1*ColliderEnergy/2d0
!
!   Mom(1,2)= x2*ColliderEnergy/2d0
!   Mom(2,2)= 0d0
!   Mom(3,2)= 0d0
!   Mom(4,2)=-x2*ColliderEnergy/2d0

  do i=1,NumPart
      MomTmp1=Mom(1,i)
      MomTmp4=Mom(4,i)
      Mom(1,i)= gamma*MomTmp1 - betagamma*MomTmp4
      Mom(4,i)= gamma*MomTmp4 - betagamma*MomTmp1
  enddo

RETURN
END SUBROUTINE






SUBROUTINE PDFMapping(MapType,yRnd,eta1,eta2,Ehat,sHatJacobi)
use ModParameters
use ModMisc
implicit none
integer :: MapType
real(8) :: yRnd(1:2),eta1,eta2,EHat,sHatJacobi,tau,nPotMap,z,sbar,fmax

  if( MapType.eq.1 ) then!  no mapping
      eta1 = yRnd(1)
      eta2 = yRnd(2)
      sHatJacobi = 1d0
  elseif( MapType.eq.2 ) then!  exponential mapping
      tau = (2d0*m_Top/Collider_Energy)**2
      eta1 = tau**yRnd(1)
      eta2 = tau**( (1d0-yRnd(1))*yRnd(2) )
      sHatJacobi = dlog(tau)**2*(1d0-yRnd(1))*eta1*eta2
  elseif( MapType.eq.3 ) then!  linear mapping
      tau = (2d0*m_Top/Collider_Energy)**2
      eta1 = (1d0-tau)*yRnd(1) + tau
      eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
      sHatJacobi = (1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)
  elseif( MapType.eq.4 ) then!  MCFM mapping
      tau = dexp(dlog(((2d0*m_Top/Collider_Energy)**2))*yRnd(1))
      eta1 = dsqrt(tau)*dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      eta2 = dsqrt(tau)/dexp(0.5d0*dlog(tau)*(1d0-2d0*yRnd(2)))
      sHatJacobi = dlog(((2d0*m_Top/Collider_Energy)**2))*tau*dlog(tau)
  elseif( MapType.eq.5 ) then!  nPotMap mapping
      nPotMap = 0.5d0
      tau = (2d0*m_Top/Collider_Energy)**2
      yRnd(1) = yRnd(1)**nPotMap
      yRnd(2) = yRnd(2)**nPotMap
      eta1 = (1d0-tau) * yRnd(1) + tau
      eta2 = ((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*yRnd(2) + tau/((1d0-tau)*yRnd(1)+tau)
      sHatJacobi=(1d0-tau)*((1d0-tau)*yRnd(1))/((1d0-tau)*yRnd(1)+tau)*nPotMap**2*((yRnd(1)*yRnd(2))**(1d0/nPotMap))**(nPotMap-1d0)
  elseif( MapType.eq.10 ) then!  Breit-Wigner mapping
!       fmax = 1d0/m_Grav/Ga_Grav * ( datan((Collider_Energy**2-m_Grav**2)/m_Grav/Ga_Grav) - datan(-m_Grav/Ga_Grav) )
!       sbar = m_Grav*Ga_Grav * dtan(fmax*yRnd(1)*m_Grav*Ga_Grav - atan(m_Grav/Ga_Grav) ) + m_Grav**2
!       z = sbar/Collider_Energy**2
!       eta1 = z + (1d0-z)*yRnd(2)
!       eta2 = z/eta1
!       sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (Collider_Energy**2*eta1*eta2 - m_Grav**2)**2 + m_Grav**2*Ga_Grav**2 )

!!! Zprime section !!!
  elseif ( MapType.eq.62) then ! BW for Zprime
     fmax = 1d0/M_Zpr/Ga_Zpr * ( datan((Collider_Energy**2-M_Zpr**2)/M_Zpr/Ga_Zpr) - datan(-M_Zpr/Ga_Zpr) )
     sbar = M_Zpr*Ga_Zpr * dtan(fmax*yRnd(1)*M_Zpr*Ga_Zpr - atan(M_Zpr/Ga_Zpr) ) + M_Zpr**2
     z = sbar/Collider_Energy**2
     eta1 = z + (1d0-z)*yRnd(2)
     eta2 = z/eta1
     sHatJacobi = fmax/Collider_Energy**2 * (1d0-z)/eta1  * ( (Collider_Energy**2*eta1*eta2 - M_Zpr**2)**2 + M_Zpr**2*Ga_Zpr**2 )
!!! End Zprime section !!!

  else

      call Error("PDF mapping not available")
  endif
  EHat = Collider_Energy*dsqrt(eta1*eta2)

RETURN
END SUBROUTINE


FUNCTION Convert2DHist(nx,ny,totx,toty)
  implicit none
  integer nx,ny,totx,toty,Convert2DHist
  
  if (nx .gt. totx .or. ny .gt. toty .or. nx .eq. 0 .or. ny .eq. 0) then
! send these events to the last bin
     Convert2DHist=totx*toty+1
  else
     Convert2DHist=nx+(ny-1)*totx
  endif


end FUNCTION Convert2DHist


SUBROUTINE setPDFs(x1,x2,MuFacIn,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale,MuFacIn
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),sbar(1:2),chm(1:2),cbar(1:2),bot(1:2),bbar(1:2),glu(1:2),phot
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2),NNpdf(1:2,-6:7)

        PDFScale=MuFacIn*100d0


#if _UseLHAPDF==1

        call evolvePDF(x1,PDFScale,NNpdf(1,-6:7))
        call evolvePDF(x2,PDFScale,NNpdf(2,-6:7))
        NNpdf(1,-6:7) = NNpdf(1,-6:7)/x1
        NNpdf(2,-6:7) = NNpdf(2,-6:7)/x2
        
        pdf(Up_,1)   = NNpdf(1,+2)         * swPDF_u
        pdf(AUp_,1)  = NNpdf(1,-2)         * swPDF_u
        pdf(Dn_,1)   = NNpdf(1,+1)         * swPDF_d
        pdf(ADn_,1)  = NNpdf(1,-1)         * swPDF_d
        pdf(Chm_,1)  = NNpdf(1,+4)         * swPDF_c
        pdf(AChm_,1) = NNpdf(1,-4)         * swPDF_c
        pdf(Str_,1)  = NNpdf(1,+3)         * swPDF_s
        pdf(AStr_,1) = NNpdf(1,-3)         * swPDF_s
        pdf(Bot_,1)  = NNpdf(1,+5)         * swPDF_b
        pdf(ABot_,1) = NNpdf(1,-5)         * swPDF_b
        pdf(0,1)     = NNpdf(1,+0)         * swPDF_g            
            
        pdf(Up_,2)   = NNpdf(2,+2)         * swPDF_u
        pdf(AUp_,2)  = NNpdf(2,-2)         * swPDF_u
        pdf(Dn_,2)   = NNpdf(2,+1)         * swPDF_d
        pdf(ADn_,2)  = NNpdf(2,-1)         * swPDF_d
        pdf(Chm_,2)  = NNpdf(2,+4)         * swPDF_c
        pdf(AChm_,2) = NNpdf(2,-4)         * swPDF_c
        pdf(Str_,2)  = NNpdf(2,+3)         * swPDF_s
        pdf(AStr_,2) = NNpdf(2,-3)         * swPDF_s
        pdf(Bot_,2)  = NNpdf(2,+5)         * swPDF_b
        pdf(ABot_,2) = NNpdf(2,-5)         * swPDF_b
        pdf(0,2)     = NNpdf(2,+0)         * swPDF_g            

#else



!   MRSW PDFS
IF( PDFSET.EQ.1 .AND. NLOPARAM.LE.1) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines, where eta/z appears
!             call mrstlo(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
! RR added -- these are needed if using MRST LO as above
!             sbar(1)=str(1)
!             cbar(1)=chm(1)
!             bbar(1)=bot(1)
!             call mrs96(x1,PDFScale,2,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            call GetAllPDFs("mstw2008lo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 ) then
!             call mrstlo(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
! RR added -- these are needed if using MRST LO as above
!             sbar(2)=str(2)
!             cbar(2)=chm(2)
!             bbar(2)=bot(2)

!             call mrs96(x2,PDFScale,2,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
            call GetAllPDFs("mstw2008lo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif
ELSEIF( PDFSET.EQ.1 .AND. NLOPARAM.EQ.2) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines, where eta/z appears
!             call mrst2004(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
!             call mrst2001(x1,PDFScale,1,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
! RR added -- needed if using MRST2001 as above
!            sbar(1)=str(1)
!            cbar(1)=chm(1)
!            bbar(1)=bot(1)
            call GetAllPDFs("mstw2008nlo",0,x1,PDFScale,upv(1),dnv(1),usea(1),dsea(1),str(1),sbar(1),chm(1),cbar(1),bot(1),bbar(1),glu(1),phot)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 ) then
!             call mrst2004(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
!             call mrst2001(x2,PDFScale,1,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
! RR added -- needed if using MRST2001 as above
!            sbar(2)=str(2)
!            cbar(2)=chm(2)
!            bbar(2)=bot(2)
            call GetAllPDFs("mstw2008nlo",0,x2,PDFScale,upv(2),dnv(2),usea(2),dsea(2),str(2),sbar(2),chm(2),cbar(2),bot(2),bbar(2),glu(2),phot)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif





!   CTEQ PDFS
ELSEIF( PDFSET.EQ.2 ) THEN
        if( x1.lt.1d0 ) then ! this is needed for integrated dipole routines
            call cteq10(x1,PDFScale,99,upv(1),dnv(1),usea(1),dsea(1),str(1),chm(1),bot(1),glu(1))
            sbar(1) = str(1)
            cbar(1) = chm(1)
            bbar(1) = bot(1)
        else
            upv(1) = 0d0
            dnv(1) = 0d0
            usea(1)= 0d0
            dsea(1)= 0d0
            str(1) = 0d0
            sbar(1)= 0d0
            chm(1) = 0d0
            cbar(1)= 0d0
            bot(1) = 0d0
            bbar(1)= 0d0
            glu(1) = 0d0
        endif
        if( x2.lt.1d0 )then
            call cteq10(x2,PDFScale,99,upv(2),dnv(2),usea(2),dsea(2),str(2),chm(2),bot(2),glu(2))
            sbar(2) = str(2)
            cbar(2) = chm(2)
            bbar(2) = bot(2)
        else
            upv(2) = 0d0
            dnv(2) = 0d0
            usea(2)= 0d0
            dsea(2)= 0d0
            str(2) = 0d0
            sbar(2)= 0d0
            chm(2) = 0d0
            cbar(2)= 0d0
            bot(2) = 0d0
            bbar(2)= 0d0
            glu(2) = 0d0
        endif
ENDIF




IF( COLLIDER.EQ.1 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = cbar(1)             * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = sbar(1)             * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bbar(1)             * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       PROTON CONTENT
        pdf(Up_,2)   = (upv(2) + usea(2))  * swPDF_u / x2
        pdf(AUp_,2)  = usea(2)             * swPDF_u / x2
        pdf(Dn_,2)   = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(ADn_,2)  = dsea(2)             * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = cbar(2)             * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = sbar(2)             * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bbar(2)             * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2

ELSEIF( COLLIDER.EQ.2 ) THEN
!       PROTON CONTENT
        pdf(Up_,1)   = (upv(1) + usea(1))  * swPDF_u / x1
        pdf(AUp_,1)  = usea(1)             * swPDF_u / x1
        pdf(Dn_,1)   = (dnv(1) + dsea(1))  * swPDF_d / x1
        pdf(ADn_,1)  = dsea(1)             * swPDF_d / x1
        pdf(Chm_,1)  = chm(1)              * swPDF_c / x1
        pdf(AChm_,1) = cbar(1)             * swPDF_c / x1
        pdf(Str_,1)  = str(1)              * swPDF_s / x1
        pdf(AStr_,1) = sbar(1)             * swPDF_s / x1
        pdf(Bot_,1)  = bot(1)              * swPDF_b / x1
        pdf(ABot_,1) = bbar(1)             * swPDF_b / x1
        pdf(0,1)     = glu(1)              * swPDF_g / x1

!       ANTI-PROTON CONTENT
        pdf(Up_,2)   = usea(2)             * swPDF_u / x2
        pdf(AUp_,2)  = (upv(2) + usea(2))  * swPDF_u / x2
        pdf(Dn_,2)   = dsea(2)             * swPDF_d / x2
        pdf(ADn_,2)  = (dnv(2) + dsea(2))  * swPDF_d / x2
        pdf(Chm_,2)  = chm(2)              * swPDF_c / x2
        pdf(AChm_,2) = cbar(2)             * swPDF_c / x2
        pdf(Str_,2)  = str(2)              * swPDF_s / x2
        pdf(AStr_,2) = sbar(2)             * swPDF_s / x2
        pdf(Bot_,2)  = bot(2)              * swPDF_b / x2
        pdf(ABot_,2) = bbar(2)             * swPDF_b / x2
        pdf(0,2)     = glu(2)              * swPDF_g / x2
ENDIF

#endif


RETURN
END SUBROUTINE








SUBROUTINE WriteLHEvent_TTB(Mom,InFlav,EventWeight,nPhoRad)
use ModParameters
use ModMisc
implicit none
real(8) :: Mom(1:4,1:13),DKRnd(1:2)
real(8),optional :: EventWeight
integer,optional :: nPhoRad
integer :: InFlav(1:2)
integer :: ICOLUP(1:2,1:13),LHE_IDUP(1:13),ISTUP(1:13),MOTHUP(1:2,1:13)
integer :: NUP,IDPRUP,i
real(8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,Lifetime,Spin,TheMass
character(len=*),parameter :: Fmt1 = "(6X,I3,2X,I3,3X,I2,3X,I2,2X,I3,2X,I3,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,X,1PE18.11,1PE18.11,X,1F3.1)"
integer, parameter :: inLeft=1,inRight=2,Xbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
integer, parameter :: io_LHEOutFile=17

! For description of the LHE format see http://arxiv.org/abs/hep-ph/0109068 and http://arxiv.org/abs/hep-ph/0609017
! The LHE numbering scheme can be found here: http://pdg.lbl.gov/mc_particle_id_contents.html and http://lhapdf.hepforge.org/manual#tth_sEcA


      IDPRUP=100
      SCALUP=MuFac * 100d0
      AQEDUP=alpha
      AQCDUP=alpha_s*RunAlphaS(NLOParam,MuRen)

      MOTHUP(1:2,inLeft) = (/0,0/);             ISTUP(inLeft) = -1
      MOTHUP(1:2,inRight)= (/0,0/);             ISTUP(inRight)= -1
      MOTHUP(1:2,tbar)   = (/inLeft,inRight/);  ISTUP(tbar)   = +2
      MOTHUP(1:2,t)      = (/inLeft,inRight/);  ISTUP(t)      = +2
      MOTHUP(1:2,bbar)   = (/tbar,tbar/);       ISTUP(bbar)   = +1
      MOTHUP(1:2,Wm)     = (/tbar,tbar/);       ISTUP(Wm)     = +2
      MOTHUP(1:2,lepM)   = (/Wm,Wm/);           ISTUP(lepM)   = +1
      MOTHUP(1:2,nubar)  = (/Wm,Wm/);           ISTUP(nubar)  = +1
      MOTHUP(1:2,b)      = (/t,t/);             ISTUP(b)      = +1
      MOTHUP(1:2,Wp)     = (/t,t/);             ISTUP(Wp)     = +2
      MOTHUP(1:2,lepP)   = (/Wp,Wp/);           ISTUP(lepP)   = +1
      MOTHUP(1:2,nu)     = (/Wp,Wp/);           ISTUP(nu)     = +1
      if( Process.ge.1 .and. Process.le.6 ) MOTHUP(1:2,6:13)=MOTHUP(1:2,6:13)-1
      
      if( .not.present(nPhoRad) .or. nPhoRad.eq.0 ) then 
          MOTHUP(1:2,Xbos) = (/t,tbar/);    ISTUP(XBos) = +1
      elseif( nPhoRad.eq.1 ) then
          MOTHUP(1:2,Xbos) = (/t,t/);       ISTUP(XBos) = +1
      elseif( nPhoRad.eq.2 ) then
          MOTHUP(1:2,Xbos) = (/Wp,Wp/);     ISTUP(XBos) = +1
      elseif( nPhoRad.eq.3 ) then
          MOTHUP(1:2,Xbos) = (/tbar,tbar/); ISTUP(XBos) = +1
      elseif( nPhoRad.eq.4 ) then
          MOTHUP(1:2,Xbos) = (/Wm,Wm/);     ISTUP(XBos) = +1
      endif


      LHE_IDUP(inLeft) = InFlav(1)
      LHE_IDUP(inRight)= InFlav(2)
      LHE_IDUP(Xbos)   = LHE_Pho_
      LHE_IDUP(tbar)   = -LHE_Top_
      LHE_IDUP(t)      = LHE_Top_
      LHE_IDUP(bbar)   = -LHE_Bot_
      LHE_IDUP(Wm)     = -LHE_Wp_
      LHE_IDUP(b)      = LHE_Bot_
      LHE_IDUP(Wp)     = LHE_Wp_
      
      if( InFlav(1).eq.LHE_Glu_ ) then
          ICOLUP(1:2,inLeft)  = (/501,510/)
          ICOLUP(1:2,inRight) = (/510,502/)    
      elseif( InFlav(1).gt.0 ) then
          ICOLUP(1:2,inLeft)  = (/501,000/)
          ICOLUP(1:2,inRight) = (/000,502/)
      else
          ICOLUP(1:2,inLeft)  = (/000,502/)
          ICOLUP(1:2,inRight) = (/501,000/)
      endif
      ICOLUP(1:2,XBos) = (/000,000/)
      ICOLUP(1:2,tbar) = (/000,502/)
      ICOLUP(1:2,t)    = (/501,000/)
      ICOLUP(1:2,bbar) = (/000,502/)
      ICOLUP(1:2,b)    = (/501,000/)
      ICOLUP(1:2,Wm)   = (/000,000/)
      ICOLUP(1:2,Wp)   = (/000,000/)
      
      call random_number(DKRnd)      
      if( TopDecays.eq.1 ) then
          if( DKRnd(1).lt.1d0/3d0 ) then
            LHE_IDUP(lepM) = LHE_ElM_
            LHE_IDUP(nubar)=-LHE_NuE_
          elseif( DKRnd(1).lt.2d0/3d0 ) then
            LHE_IDUP(lepM) = LHE_MuM_
            LHE_IDUP(nubar)=-LHE_NuM_ 
          else
            LHE_IDUP(lepM) = LHE_TaM_
            LHE_IDUP(nubar)=-LHE_NuT_
          endif

          if( DKRnd(2).lt.1d0/3d0 ) then
            LHE_IDUP(lepP)=-LHE_ElM_
            LHE_IDUP(nu)  = LHE_NuE_
          elseif( DKRnd(2).lt.2d0/3d0 ) then
            LHE_IDUP(lepP)=-LHE_MuM_ 
            LHE_IDUP(nu)  = LHE_NuM_
          else
            LHE_IDUP(lepP)= -LHE_TaM_
            LHE_IDUP(nu)  = LHE_NuT_
          endif
          ICOLUP(1:2,lepM) = (/000,000/)
          ICOLUP(1:2,nubar)= (/000,000/)
          ICOLUP(1:2,lepP) = (/000,000/)
          ICOLUP(1:2,nu)   = (/000,000/)
      
      elseif( TopDecays.eq.4 ) then
          if( DKRnd(1).lt.1d0/2d0 ) then
            LHE_IDUP(lepM) = LHE_Dn_
            LHE_IDUP(nubar)=-LHE_Up_
          else
            LHE_IDUP(lepM) = LHE_Str_
            LHE_IDUP(nubar)=-LHE_Chm_
          endif
          if( DKRnd(2).lt.1d0/3d0 ) then
            LHE_IDUP(lepP)=-LHE_ElM_
            LHE_IDUP(nu)  = LHE_NuE_
          elseif( DKRnd(2).lt.2d0/3d0 ) then
            LHE_IDUP(lepP)=-LHE_TaM_ 
            LHE_IDUP(nu)  = LHE_NuM_
          else
            LHE_IDUP(lepP)= -LHE_TaM_
            LHE_IDUP(nu)  = LHE_NuT_
          endif
          ICOLUP(1:2,lepM) = (/601,000/)
          ICOLUP(1:2,nubar)= (/000,601/)
          ICOLUP(1:2,lepP) = (/000,000/)
          ICOLUP(1:2,nu)   = (/000,000/)

          
      elseif( TopDecays.eq.3 ) then
          if( DKRnd(1).lt.1d0/3d0 ) then
            LHE_IDUP(lepM) = LHE_ElM_
            LHE_IDUP(nubar)=-LHE_NuE_
          elseif( DKRnd(1).lt.2d0/3d0 ) then
            LHE_IDUP(lepM) = LHE_MuM_
            LHE_IDUP(nubar)=-LHE_NuM_ 
          else
            LHE_IDUP(lepM) = LHE_TaM_
            LHE_IDUP(nubar)=-LHE_NuT_
          endif
          if( DKRnd(2).lt.1d0/2d0 ) then
            LHE_IDUP(lepP)=-LHE_Dn_
            LHE_IDUP(nu)  = LHE_Up_
          else
            LHE_IDUP(lepP)=-LHE_Str_
            LHE_IDUP(nu)  = LHE_Chm_
          endif
          ICOLUP(1:2,lepM) = (/000,000/)
          ICOLUP(1:2,nubar)= (/000,000/)
          ICOLUP(1:2,lepP) = (/000,701/)
          ICOLUP(1:2,nu)   = (/701,000/)
      
      else! TopDecays.eq.2
          if( DKRnd(1).lt.1d0/2d0 ) then
            LHE_IDUP(lepM) = LHE_Dn_
            LHE_IDUP(nubar)=-LHE_Up_
          else
            LHE_IDUP(lepM) = LHE_Str_
            LHE_IDUP(nubar)=-LHE_Chm_
          endif
          if( DKRnd(2).lt.1d0/2d0 ) then
            LHE_IDUP(lepP)=-LHE_Dn_
            LHE_IDUP(nu)  = LHE_Up_
          else
            LHE_IDUP(lepP)=-LHE_Str_
            LHE_IDUP(nu)  = LHE_Chm_
          endif
          ICOLUP(1:2,lepM) = (/601,000/)
          ICOLUP(1:2,nubar)= (/000,601/)
          ICOLUP(1:2,lepP) = (/000,701/)
          ICOLUP(1:2,nu)   = (/701,000/)
      endif
      

      if( TopDecays.eq.0 ) then
        NUP = 5
        ISTUP(tbar) = +1
        ISTUP(t)    = +1 
      else
        NUP=13
      endif

      if( present(EventWeight) ) then
          XWGTUP=EventWeight
      else
          XWGTUP=1.0d0
      endif
      Lifetime = 0.0d0
      Spin     = 0.1d0


      write(io_LHEOutFile,"(A)") "<event>"
      write(io_LHEOutFile,"(I2,X,I3,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7,2X,1PE14.7)") NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      ! in order of appearance:
      ! (*) number of particles in the event
      ! (*) process ID (user defined)
      ! (*) weighted or unweighted events: +1=unweighted, otherwise= see manual
      ! (*) pdf factorization scale in GeV
      ! (*) alpha_QED coupling for this event 
      ! (*) alpha_s coupling for this event

      do i=1,NUP
          if( i.eq.3 .and. Process.ge.1 .and. Process.le.6 ) cycle
      !      TheMass = GetMass( MY_IDUP(i) )*100d0
          TheMass = get_Minv(Mom(:,i))
          if( TheMass/Mom(1,i).lt.1d-7 ) TheMass = 0.0d0
          write(io_LHEOutFile,fmt1) LHE_IDUP(i),ISTUP(i), MOTHUP(1,i),MOTHUP(2,i), ICOLUP(1,i),ICOLUP(2,i),Mom(2:4,i)*100d0,Mom(1,i)*100d0,TheMass*100d0,Lifetime,Spin
      enddo
      write(io_LHEOutFile,"(A)") "</event>"

RETURN
END SUBROUTINE











SUBROUTINE TTbar_OffShellProjection(MomIn,MomOut,Jacobian)
use modParameters
use modMisc
implicit none
real(8) :: MomIn(1:4,1:13),MomOut(:,:),MomTmp(1:4),Jacobian
real(8) :: xRndWidth(2:5),BW_Mass(2:5),BW_Jacobi(2:5)
integer, parameter :: inLeft=1,inRight=2,Xbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13

    call random_number(xRndWidth)
    
    call SmearExternal(xRndWidth(2),m_top,Ga_TopExp,m_top-6d0*Ga_TopExp,m_top+6d0*Ga_TopExp,BW_Mass(2),BW_Jacobi(2))
    call SmearExternal(xRndWidth(3),m_top,Ga_TopExp,m_top-6d0*Ga_TopExp,m_top+6d0*Ga_TopExp,BW_Mass(3),BW_Jacobi(3))
    call SmearExternal(xRndWidth(4),m_W,Ga_WExp,m_W-6d0*Ga_WExp,m_W+6d0*Ga_WExp,BW_Mass(4),BW_Jacobi(4))
    call SmearExternal(xRndWidth(5),m_W,Ga_WExp,m_W-6d0*Ga_WExp,m_W+6d0*Ga_WExp,BW_Mass(5),BW_Jacobi(5))
    Jacobian = BW_Jacobi(2) * BW_Jacobi(3) * BW_Jacobi(4) * BW_Jacobi(5)    

    call ShiftMass(MomIn(1:4,tbar),MomIn(1:4,t),BW_Mass(2),BW_Mass(3),MomOut(1:4,tbar),MomOut(1:4,t))
    
    MomTmp(1:4) = MomOut(1:4,t) - MomIn(1:4,Wp)
    call ShiftMass(MomTmp,MomIn(1:4,Wp),m_BotExp,BW_Mass(4),MomOut(1:4,b),MomOut(1:4,Wp))
    
    MomTmp(1:4) = MomOut(1:4,tbar) - MomIn(1:4,Wm)
    call ShiftMass(MomTmp,MomIn(1:4,Wm),m_BotExp,BW_Mass(5),MomOut(1:4,bbar),MomOut(1:4,Wm))
    
    MomTmp(1:4) = MomOut(1:4,Wp) - MomIn(1:4,lepP)
    MomOut(1:4,nu)   = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP)) * MomIn(1:4,lepP)
    MomOut(1:4,lepP) = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepP))) * MomIn(1:4,lepP)
    
    MomTmp(1:4) = MomOut(1:4,Wm) - MomIn(1:4,lepM)
    MomOut(1:4,nubar) = MomTmp(1:4) - (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM)) * MomIn(1:4,lepM)
    MomOut(1:4,lepM)  = (1d0 + (MomTmp(1:4).dot.MomTmp(1:4))/2d0/(MomTmp(1:4).dot.MomIn(1:4,lepM))) * MomIn(1:4,lepM)

return
END SUBROUTINE





SUBROUTINE ShiftMass(p1,p2,m1,m2,p1hat,p2hat)
use ModMisc
implicit none
real(8),intent(in) :: p1(1:4),p2(1:4)
real(8) :: m1,m2,p1hat(1:4),p2hat(1:4)
real(8) :: xi,eta,a,b,c,p1sq,p2sq,p1p2

  p1sq = p1(1:4).dot.p1(1:4)
  p2sq = p2(1:4).dot.p2(1:4)
  p1p2 = p1(1:4).dot.p2(1:4)

  a = ( p1sq*p2(1:4) - p2sq*p1(1:4) + p1p2*(p2(1:4)-p1(1:4)) ).dot.( p1sq*p2(1:4) - p2sq*p1(1:4) + p1p2*(p2(1:4)-p1(1:4)) )
  b = ( p1sq+p2sq+2d0*p1p2+m2**2-m1**2 ) * ( p1p2**2 - p1sq*p2sq )
  c = 0.25d0*( p1sq+p2sq+2d0*p1p2+m2**2-m1**2 )**2*p1sq - (p1sq+p1p2)**2*m2**2
  eta = 1d0/2d0/a * ( -b - dsqrt( dabs(b**2 -4d0*a*c) ) )
  xi = ( p1sq+p2sq+2d0*p1p2 + m2**2 - m1**2 - 2d0*eta*(p2sq+p1p2) )/2d0/( p1sq + p1p2 )

  p2hat(1:4) = xi*p1(1:4) + eta*p2(1:4)
  p1hat(1:4) = (1d0-xi)*p1(1:4) + (1d0-eta)*p2(1:4)


RETURN
END SUBROUTINE







END MODULE
