#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

double vX = 0, vY = 0, vZ = 0; // event vertex


TDatime dt;

void runBkgVT(Int_t nevents = 100, 
	      double Eint = 160.,
	      const char *setup = "setup-10um-itssa_Eff1.txt",
	      bool simulateBg=kTRUE)
{

  int refreshBg = 100;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  //gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8210/xmldoc")); // check if pythia8 path is set correctly !!!!
  
  
  TH3F *h3DPiBkg = new TH3F("h3DPiBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DKBkg = new TH3F("h3DKBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH3F *h3DPBkg = new TH3F("h3DPBkg", "pt,y,phi pions", 50, 0., 5., 50, 0., 5., 50, 0, 2 * TMath::Pi());
  TH1F *hNevents = new TH1F("hNevents", "", 1, 0, 1);
  TH1F* hGenStat = new TH1F("hGenStat","",18,0.5,18.5);
  hGenStat->GetXaxis()->SetBinLabel(1,"#pi to gen");
  hGenStat->GetXaxis()->SetBinLabel(2,"#pi bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(3,"#pi bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(4,"#pi bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(5,"#pi in tree");
  hGenStat->GetXaxis()->SetBinLabel(6,"#pi with fake clusters");
  hGenStat->GetXaxis()->SetBinLabel(7,"K to gen");
  hGenStat->GetXaxis()->SetBinLabel(8,"K bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(9,"K bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(10,"K bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(11,"K in tree");
  hGenStat->GetXaxis()->SetBinLabel(12,"K with fake clusters");
  hGenStat->GetXaxis()->SetBinLabel(13,"p to gen");
  hGenStat->GetXaxis()->SetBinLabel(14,"p bad SolveSingleTrack");
  hGenStat->GetXaxis()->SetBinLabel(15,"p bad GetWinnerMCTrack");
  hGenStat->GetXaxis()->SetBinLabel(16,"p bad chi2");
  hGenStat->GetXaxis()->SetBinLabel(17,"p in tree");
  hGenStat->GetXaxis()->SetBinLabel(18,"p with fake clusters");
  
  
  //Magnetic field and detector parameters
  MagField *mag = new MagField(1);
  int BNreg = mag->GetNReg();
  const double *BzMin = mag->GetZMin();
  const double *BzMax = mag->GetZMax();
  const double *BVal;
  printf("*************************************\n");
  printf("number of magnetic field regions = %d\n", BNreg);
  for (int i = 0; i < BNreg; i++){
    BVal = mag->GetBVals(i);
    printf("*** Field region %d ***\n", i);
    if (i == 0){
      printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }else if (i == 1){
      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }
  }

  KMCDetectorFwd *det = new KMCDetectorFwd();
  printf("Setup file = %s\n",setup);
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT
  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
  // we don't need MS part here, even if it is in the setup
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); //NA60+
  //det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
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
  det->SetMinP2Propagate(1); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  //
  det->BookControlHistos();

  // Get Pi, K, P spectral shapes
  TF1* fdNdYPi=det->GetdNdYPi();
  TF1* fdNdYK=det->GetdNdYK();
  TF1* fdNdYP=det->GetdNdYP();
  TF1* fdNdPtPi=det->GetdNdPtPi();
  TF1* fdNdPtK=det->GetdNdPtK();
  TF1* fdNdPtP=det->GetdNdPtP();





  TFile *f = new TFile("treeBkgEvents.root", "RECREATE");
  TTree *tree = new TTree("tree", "tree Bkg");
  TClonesArray *arrtr = new TClonesArray("KMCProbeFwd");
  TClonesArray &aarrtr = *arrtr;
  tree->Branch("tracks", &arrtr);
  
  for (Int_t iev = 0; iev < nevents; iev++){
    aarrtr.Clear();
    printf(" ***************  ev = %d \n", iev);
    double pxyz[3];
    hNevents->Fill(0.5);
    if (simulateBg && (iev % refreshBg) == 0)
      det->GenBgEvent(vX, vY, vZ);
    
    double ntrPi = gRandom->Poisson(det->GetNChPi());
    //printf("fNChPi=%f ntrPi=%f\n", det->GetNChPi(), ntrPi);

    double yrap, pt, phi;
    int charge;
    double mass;
    Int_t icount = 0;
    for (int itr = 0; itr < ntrPi; itr++){	
      yrap = fdNdYPi->GetRandom();
      pt = fdNdPtPi->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.52 ? 1 : -1;
      mass = KMCDetectorFwd::kMassPi;
      h3DPiBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(1);

      TLorentzVector *ppi = new TLorentzVector(0., 0., 0., 0.);
      ppi->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      if (!det->SolveSingleTrack(ppi->Pt(), ppi->Rapidity(), ppi->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
	hGenStat->Fill(2);
	continue;
      }
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw){
	hGenStat->Fill(3);
	continue;
      }
      if (trw->GetNormChi2(kTRUE) > ChiTot){
	hGenStat->Fill(4);
	continue;
      }
      trw->GetPXYZ(pxyz);
      //printf("charge = %d, %f \n", charge, trw->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw);
      hGenStat->Fill(5);
      if(trw->GetNFakeITSHits()>0) hGenStat->Fill(6);
      icount++;
    }

    // kaons
    double ntrK = gRandom->Poisson(det->GetNChK());
    //printf("fNChK=%f ntrK=%f\n", det->GetNChK(), ntrK);
    
    for (int itr = 0; itr < ntrK; itr++){
      yrap = fdNdYK->GetRandom();
      pt = fdNdPtK->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = gRandom->Rndm() > 0.3 ? 1 : -1;
      mass = KMCDetectorFwd::kMassK;
      h3DKBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(7);
      
      TLorentzVector *pk = new TLorentzVector(0., 0., 0., 0.);
      pk->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      if (!det->SolveSingleTrack(pk->Pt(), pk->Rapidity(), pk->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
	hGenStat->Fill(8);
	continue;
      }
      KMCProbeFwd *trw2 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw2){
	hGenStat->Fill(9);
	continue;
      }
      if (trw2->GetNormChi2(kTRUE) > ChiTot){
	hGenStat->Fill(10);
	continue;
      }
      trw2->GetPXYZ(pxyz);
      //printf("charge = %d, %f \n", charge, trw2->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw2);
      hGenStat->Fill(11);
      if(trw2->GetNFakeITSHits()>0) hGenStat->Fill(12);
      icount++;
    }
    
    // protons
    double ntrP = gRandom->Poisson(det->GetNChP());
    //printf("fNChP=%f ntrP=%f\n", det->GetNChP(), ntrP);
    for (int itr = 0; itr < ntrP; itr++){
      yrap = fdNdYP->GetRandom();
      pt = fdNdPtP->GetRandom();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      charge = 1;
      mass = KMCDetectorFwd::kMassP;
      h3DPBkg->Fill(pt, yrap, phi);
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + mass * mass) * TMath::SinH(yrap)};
      hGenStat->Fill(13);
      
      TLorentzVector *pp = new TLorentzVector(0., 0., 0., 0.);
      pp->SetXYZM(pxyz[0], pxyz[1], pxyz[2], mass);
      if (!det->SolveSingleTrack(pp->Pt(), pp->Rapidity(), pp->Phi(), mass, charge, vX, vY, vZ, 0, 1, 99)){
	hGenStat->Fill(14);
	continue;
      }
      KMCProbeFwd *trw3 = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw3){
	hGenStat->Fill(15);
	continue;
      }
      if (trw3->GetNormChi2(kTRUE) > ChiTot){
	hGenStat->Fill(16);
	continue;
      }
      trw3->GetPXYZ(pxyz);
      //printf("charge = %d, %f \n", charge, trw3->GetCharge());
      new (aarrtr[icount]) KMCProbeFwd(*trw3);
      hGenStat->Fill(17);
      if(trw3->GetNFakeITSHits()>0) hGenStat->Fill(18);
      icount++;
    }
    printf("Pions+Kaons+Protons in array = %d out of %.0f \n",icount,ntrPi+ntrK+ntrP);
    tree->Fill();
  }
  f->cd();
  tree->Write();
  f->Close();
  
  TFile *outfile = new TFile("bkgdistributions.root", "recreate");
  outfile->cd();
  hNevents->Write();
  hGenStat->Write();
  h3DPiBkg->Write();
  h3DKBkg->Write();
  h3DPBkg->Write();
  outfile->Close();
}
