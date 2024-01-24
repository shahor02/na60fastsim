baseDir=$PWD
fileDir=$baseDir/files

seed=$RANDOM
nev=1000
#ProcType='rho'    
ProcType='omega2Body'
#ProcType='omegaDalitz'
#ProcType='etaDalitz'
#ProcType='eta2Body'
#ProcType='phi'
#ProcType='etaPrime'
#ProcType='jpsi'

Energy=40 
Resolution='10um'
#MagnetType='Toroid1R-B0.5-Rmin2cm-Z1.5m-DeadZonesOn' 
MagnetType='Toroid1R-B0.1Tm-Rmin30cm-Rmax300cm-Z3.0m' 
#MagnetType='Toroid1R-B0.2Tm-Rmin30cm-Rmax280cm-Z3.0m' 
#MagnetType='without_toroid_without_dipole'
#MagnetType='B0.3Tm-ACM' 
Chi2Cut='Chi21.5'       
Wall='Wall180'
MinITSHits='MinITSHits5'
#MagnetDeadZones='ToroidDeadZones'
        
localdir=$baseDir/runDiMu-$ProcType-$Resolution-$Chi2Cut-$nev-ev-${Energy}GeV-$MagnetType-$Wall-$MinITSHits-Tree_${seed}

mkdir $localdir
cd $localdir
cp $fileDir/*.{h,cxx} ./
# overwrite magnetic field implementation
cp $fileDir/KMCUtilsTLowEnergy_2reg_setupII.cxx $localdir/KMCUtils.cxx
cp $fileDir/macros/runDiMuGenLMR.C . 

#cp $fileDir/setup-$Resolution-toroid-$Wall.txt $localdir/setup.txt
cp $fileDir/setup-$Resolution-itssa.txt $localdir/setup.txt

#gdb aliroot

aliroot <<EOF >runDiMuGenLMR.out 2>runDiMuGenLMR.err
  gROOT->ProcessLine(".L ./KMCUtils.cxx+g");
  gROOT->ProcessLine(".L ./NaMaterial.cxx+g");
  gROOT->ProcessLine(".L ./KMCProbeFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCClusterFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCLayerFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCPolyLayer.cxx+g");
  gROOT->ProcessLine(".L ./KMCMSStation.cxx+g");
  gROOT->ProcessLine(".L ./KMCDetectorFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCFlukaParser.cxx+g");
  gROOT->ProcessLine(".L ./KMCMagnetBuilder.cxx+g");
  gROOT->ProcessLine(".L ./GenMUONLMR.cxx+g");

  gROOT->ProcessLine(".L LocLog.cxx+g");
  gROOT->ProcessLine(".L TrackPar.cxx+g");
  gROOT->ProcessLine(".L TLocTreeStream.cxx+g");
    

.L runDiMuGenLMR.C+g
runDiMuGenLMR($nev,"$ProcType",$seed,$Energy) 

EOF
