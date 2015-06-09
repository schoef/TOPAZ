# TOPAZ

TOPAZ is a Fortran90 program that allows precise calculations of various top quark processes at hadron and lepton colliders.
Production and decay dynamics are described at NLO QCD through the narrow-width approximation for top quarks.

The code requires the intel fortran compiler (ifort) and can optionally be run using MPI parallelization. Interfaces to FastJet and LHAPDFS are available but not required.
The QCDLoop-1.9 library is included to evaluate the 1-loop master integrals. CTEQ6L1, CTEQ10, MSTW8 and NNPDFS3.0 pdf sets are included in the package.


The following processes are available:

 - ttbar
 - ttbar+jet
 - ttbar+photon
 - ttbar+Z
 - ttbar+Higgs
 - single top + H
 - Zprime --> ttbar
 - fourth generation TT --> ttbar+ETmiss
 - scalar top partners stops --> ttbar+ETmiss 


Basic parameters in TOPAZ are set in mod_Parameters.f90. 
Run-time parameters are set via the command line.
Kinematic cuts and histograms are defined in mod_Kinematics.



