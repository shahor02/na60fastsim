void load(){

  gROOT->ProcessLine(".L ./NaMaterial.cxx+g");
  gROOT->ProcessLine(".L ./KMCUtils.cxx+g");
  gROOT->ProcessLine(".L ./KMCProbeFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCClusterFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCLayerFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCDetectorFwd.cxx+g");
  gROOT->ProcessLine(".L ./KMCFlukaParser.cxx+g");
  gROOT->ProcessLine(".L ./GenMUONLMR.cxx+g");
  
}
