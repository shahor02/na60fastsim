#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include "KMCProbeFwd.h"
#include "KMCDetectorFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "./HFUtils.C"
#endif

constexpr double kC = 0.0299792458;
constexpr double kTau = 274.;
constexpr double kCTau = kTau * kC;
constexpr double kHyperMass = (3727.380 + 1115.683 - 3.102) * 1.e-3;
constexpr double kDaughterMasses[] = {3.727380, 0.938272, 0.139570};
constexpr double kDaughterCharges[] = {2., 1., -1.};
constexpr int kNdaughters = 3;
constexpr int seed = 1;

struct MiniH
{
  float px, py, pz;
  float x, y, z;
  float m;
  float ProngPvDCA[kNdaughters];
  float pxMC, pyMC, pzMC;
  float xMC, yMC, zMC;
  unsigned char NClusters[kNdaughters];
  char reconstructed = 0;
  char fake = 0;
};

// Track cuts
double ChiTot = 1.5;
int minITShits = 4;

KMCDetectorFwd* createDetectorFwd(double Eint, const char *setup) {
  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

  det->SetMinITSHits(TMath::Min(minITShits, det->GetNumberOfActiveLayersITS())); // NA60+
  // det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); // NA60+
  // det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(0);
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); // NA60+
  // det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  //

  TVirtualMagField *fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class())
  {
    MagField *mag = (MagField *)fld;
    int BNreg = mag->GetNReg();
    const double *BzMin = mag->GetZMin();
    const double *BzMax = mag->GetZMax();
    const double *BVal;
    printf("*************************************\n");
    printf("number of magnetic field regions = %d\n", BNreg);
    for (int i = 0; i < BNreg; i++)
    {
      BVal = mag->GetBVals(i);
      printf("*** Field region %d ***\n", i);
      if (i == 0)
      {
        printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
      else if (i == 1)
      {
        printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
    }
  }
  return det;
}

void runHyperNucleiCandidates(int nevents = 14000000,
                              double Eint = 160.,
                              const char *setup = "setup.txt",
                              const char *filNamPow = "AnalysisResults.root",
                              bool writeNtuple = true,
                              bool simulateBg = true)
{
  int refreshBg = 10000;
  gRandom->SetSeed(seed);

  // MUSIC input
  printf("--> pt and y shape of the hypernuclei from %s\n", filNamPow);
  TFile *filPow = new TFile(filNamPow);
  TH3D *ptYetaShape = (TH3D *)filPow->Get("ptYeta");
  TH1D *ptShape = (TH1D *)ptYetaShape->ProjectionX("ptShape");
  TH1D *yShape = (TH1D *)ptYetaShape->ProjectionY("yShape");

  KMCDetectorFwd* det = createDetectorFwd(Eint, setup);

  // prepare decays
  MiniH hyper;
  TGenPhaseSpace decay;
  TLorentzVector parentgen, daugen[3], parent, daurec[3];
  KMCProbeFwd recProbe[3];

  TFile fnt("Hypernuclei-Signal-ntuple.root", "recreate");
  TTree tree("hyper", "hyper");
  tree.Branch("h", &hyper);

  for (int iev = 0; iev < nevents; iev++)
  {
    double vprim[3] = {0, 0, 0};
    if (iev % 100 == 0)
      printf(" ***************  ev = %d \n", iev);
    double pxyz[3];

    if (simulateBg && (iev % refreshBg) == 0)
      det->GenBgEvent(0., 0., 0.);
    double ptGenD = ptShape->GetRandom(); // get Lc distribution from file
    double yGenD = yShape->GetRandom();
    double phi = gRandom->Rndm() * 2 * TMath::Pi();
    hyper.pxMC = ptGenD * std::cos(phi);
    hyper.pyMC = ptGenD * std::sin(phi);

    double mt = std::hypot(ptGenD, kHyperMass);
    hyper.pzMC = mt * TMath::SinH(yGenD);
    double totGenMom{std::hypot(ptGenD, hyper.pzMC)};
    double en = mt * TMath::CosH(yGenD);
    hyper.reconstructed = 0;
    hyper.fake = 0;

    parentgen.SetPxPyPzE(hyper.pxMC, hyper.pyMC, hyper.pzMC, en);

    double decayLenght = gRandom->Exp(kCTau) * parentgen.Beta() * parentgen.Gamma();
    hyper.zMC = decayLenght * hyper.pzMC / totGenMom;
    hyper.xMC = decayLenght * ptGenD * std::cos(phi) / totGenMom;
    hyper.yMC = decayLenght * ptGenD * std::sin(phi) / totGenMom;
    decay.SetDecay(parentgen, 3, kDaughterMasses);
    decay.Generate();

    for (int iD{0}; iD < kNdaughters; ++iD)
    {
      daugen[iD] = *(decay.GetDecay(iD));
      if (!det->SolveSingleTrack(daugen[iD].Pt(), daugen[iD].Rapidity(), daugen[iD].Phi(), kDaughterMasses[iD], kDaughterCharges[iD], hyper.xMC, hyper.yMC, hyper.zMC, 0, 1, 99))
        continue;
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw || trw->GetNormChi2(kTRUE) > ChiTot)
        continue;
      hyper.reconstructed++;

      hyper.NClusters[iD] = trw->GetNITSHits();
      hyper.fake += trw->GetNFakeITSHits();
      recProbe[iD] = *trw;
    }

    if (hyper.reconstructed < 3) {
      tree.Fill();
      continue;
    }

    for (int iD{0}; iD < kNdaughters; ++iD) {
      hyper.ProngPvDCA[iD] = std::hypot(recProbe[iD].GetX(), recProbe[iD].GetY()) * (recProbe[iD].GetX() > 0 ? 1. : -1.);
    }

    double xV, yV, zV;
    double sigmaVert;
    ComputeVertex(recProbe[0], recProbe[1], recProbe[2], vprim[2], xV, yV, zV, sigmaVert);
    hyper.x = xV;
    hyper.y = yV;
    hyper.z = zV;

    // Get daughter track momentum at decay vertex
    for (int idau = 0; idau < kNdaughters; idau++)
    {
      recProbe[idau].PropagateToZBxByBz(zV);
      recProbe[idau].GetPXYZ(pxyz);
      daurec[idau].SetXYZM(pxyz[0], pxyz[1], pxyz[2], kDaughterMasses[idau]);
    }

    parent = daurec[0];
    parentgen = daugen[0];
    for (int iD = 1; iD < kNdaughters; ++iD) {
      parent += daurec[iD];
      parentgen += daugen[iD];
    }

    hyper.px = parent.Px();
    hyper.py = parent.Py();
    hyper.pz = parent.Pz();
    hyper.m = parent.M();

    tree.Fill();
  }

  fnt.cd();
  tree.Write();
  fnt.Close();
}
