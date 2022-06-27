#if !defined(__CINT__) || defined(__MAKECINT__)
#define _NOALIROOT_
#include "GenMUONLMR.h"
#include "KMCDetectorFwd.h"
#include "KMCFlukaParser.h"
#include "KMCProbeFwd.h"
#include "TF1.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"
#include "TH2.h"
#include "TLocTreeStream.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TTree.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

double vX = 0, vY = 0, vZ = 0; // event vertex
KMCDetectorFwd *det = 0;

void runDiMu_Geant(
    int nev = 30000,    // n events to generate
    int refreshBg = 10, // generate new bg event for each refreshBg-th, 1 to refresh for every signal
    const char *setup = "setups/setup-EHN1_40GeV_5pixel.txt", // setup.txt", // setup to load
    const char *geantList = "g.lst",
    // arguments below are at the moment not to be used
    const char *flukaBGList = "", //"fluka.lst", // optional fluka background file, if empty, then use parametric background
    const char *interactionSource = "" // optional primary interaction volume in  fluka files, if empty, take all
) {
  TLocTreeSRedirector outStream("dimuGenLMR.root"); // the output stream trees will go here
  TString flukaBG = flukaBGList;
  TString signalGeant = geantList;
  MagField *mag = new MagField(1);
  int BNreg = mag->GetNReg();
  const double *BzMin = mag->GetZMin();
  const double *BzMax = mag->GetZMax();
  const double *BVal;
  printf("*************************************\n");
  printf("number of magnetic field regions = %d\n", BNreg);
  for (int i = 0; i < BNreg; i++) {
    BVal = mag->GetBVals(i);
    printf("*** Field region %d ***\n", i);
    if (i == 0) {
      printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    } else if (i == 1) {
      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
    }
  }

  int outN = 1;//nev / 10;
  if (outN < 1) outN = 1;
  //
  //
  det = new KMCDetectorFwd();
  // det->SetUseRPhiErrorMS(true);
  det->ReadSetup(setup, setup);
  det->SetExternalInput(kTRUE); // !!! important
  
  
  KMCFlukaParser flukaParserSig;
  KMCFlukaParser *flukaParser;

  flukaParserSig.SetInpList(geantList);

  /*
  if (!flukaBG.IsNull()) {
    printf("Will use Fluka background from %s\n", flukaBG.Data());
    flukaParser = new KMCFlukaParser();
    flukaParser->SetInpList(flukaBG.Data());
  }
  */

  // set the min N tracker hits needed to validate the track
  printf("min number of hits in MS = %d\n", det->GetNumberOfActiveLayersMS());
  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); // NA60+
  det->SetMinMSHits(det->GetNumberOfActiveLayersMS());   // NA60+
  det->SetMinTRHits(det->GetNumberOfActiveLayersTR());
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);

  // set chi2 cuts
  det->SetMaxChi2Cl(40.);  // max track to cluster chi2
  det->SetMaxChi2NDF(13.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx

  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); // NA60+
  // det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //
  det->SetApplyBransonPCorrection(); // Branson correction

  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parent, dimuGen, muGen[2], muRec[2], muRecMS[2], dimuRecMS,
      dimuRec; // eli
  int npix[2] = {0}, nfake[2] = {0};
  float chiGlo[2] = {999., 999.};
  float chiGloITS[2] = {999., 999.};

  Double_t prodM[2] = {KMCDetectorFwd::kMassMu, KMCDetectorFwd::kMassMu};
  const int crgMu[2] = {-1, 1};
  //
  det->BookControlHistos();

  for (int iev = 0; iev < nev; iev++) {
    //
    if ((iev % outN) == 0) printf("Done %d out of %d\n", iev, nev);
    auto &simEv = flukaParserSig.getSimEvent();
    simEv.clear();
    if (!flukaParserSig.readNextTrackGeant() || !flukaParserSig.readNextTrackGeant()) {
      return;
    }
    if (simEv.signal.size() < 2) {
      printf("Got %zu signal particles at event %d\n", simEv.signal.size(), iev);
      return;
    }
    /*
    if (flukaParser) {
      det->ImposeFlukaBackground(flukaParser, interactionSource, true); // allow
    rewind
    }
    */
    det->ImposeBackgroundHits(simEv.bgHits);

    int nrec = 0;
    int nfakeHits = 0;
    double pxyz[3];
    const TParticle *fMu[2] = {&simEv.signal[0], &simEv.signal[1]};

    for (int imu = 0; imu < 2; imu++) {
      muGen[imu].SetXYZM(fMu[imu]->Px(), fMu[imu]->Py(), fMu[imu]->Pz(), 0.105658369);
    }
    dimuGen = muGen[0];
    dimuGen += muGen[1];

    for (int imu = 0; imu < 2; imu++) {
      TParticlePDG *particle = fMu[imu]->GetPDG();
      int crg = particle->Charge() / 3;
      //      if
      //      (!det->SolveSingleTrack(fMu[imu]->Pt(),fMu[imu]->Y(),fMu[imu]->Phi(),particle->Mass(),particle->Charge()/3,
      //      fMu[imu]->Vx(), fMu[imu]->Vy(),fMu[imu]->Vz(), 0, 1, 99))
      //      continue;
      /* test
      det->SetExternalInput(kFALSE);
      det->SolveSingleTrack(fMu[imu]->Pt(),fMu[imu]->Y(),fMu[imu]->Phi(),fMu[imu]->GetMass(),crg,fMu[imu]->Vx(), fMu[imu]->Vy(), fMu[imu]->Vz(), 0,1,99);
      det->SetExternalInput(kTRUE);
      */
      det->ImposeSignalHits(simEv.signalHits[imu]);
      det->CreateProbe(det->GetProbe(), fMu[imu]->Pt(), fMu[imu]->Y(), fMu[imu]->Phi(), fMu[imu]->GetMass(), crg, fMu[imu]->Vx(), fMu[imu]->Vy(), fMu[imu]->Vz());
      det->SetLastActiveLayerTracked(det->GetLastActiveLayer());
      if (!det->SolveSingleTrackViaKalmanMC(999)) break;

      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw) break;
      printf("Doing mu = %d\n", imu);
      trw->Print();

      nfakeHits += trw->GetNFakeITSHits();
      trw->GetPXYZ(pxyz);
      muRec[imu].SetXYZM(pxyz[0], pxyz[1], pxyz[2], prodM[imu]);
      npix[imu] = trw->GetNITSHits();
      nfake[imu] = trw->GetNFakeITSHits();
      chiGlo[imu] = trw->GetNormChi2(kTRUE);
      chiGloITS[imu] = trw->GetNormChi2ITS(kTRUE);

      KMCProbeFwd *muMS = det->GetMuBransonCorrVtx();
      if (muMS) {
        muMS->GetPXYZ(pxyz);
        muRecMS[imu].SetXYZM(pxyz[0], pxyz[1], pxyz[2],
                             KMCDetectorFwd::kMassMu);
      }

      KMCProbeFwd *seed0 = det->GetLayerTR(1)->GetMCTrack(0);
      KMCProbeFwd *seed1 = det->GetLayerTR(1)->GetMCTrack(1);
      if (!seed0 || !seed1) {
        printf("failed to get seed at L11: %p %p\n", seed0, seed1);
        break;
      }
      outStream << "probeRec"
                << "genMu=" << &muGen[imu] << "recMu=" << &muRec[imu]
                << "recMuMS=" << &muRecMS[imu] << "npix=" << npix[imu]
                << "nfake=" << nfake[imu] << "chi=" << chiGlo[imu]
                << "chiITS=" << chiGloITS[imu] << "probe=" << trw
                << "pr11=" << seed0 << "pr11s=" << seed1 << "\n";

      if (trw->GetNormChi2(kTRUE) > ChiTot)
        continue;
      nrec++;
    }

    if (nrec < 2)
      continue;

    dimuRec = muRec[0];
    dimuRec += muRec[1];
    dimuRecMS = muRecMS[0];
    dimuRecMS += muRecMS[1];

    outStream << "genrecAcc"
              << "gen=" << &dimuGen << "rec=" << &dimuRec
              << "recMS=" << &dimuRecMS << "genMu0=" << &muGen[0]
              << "genMu1=" << &muGen[1] << "recMu0=" << &muRec[0]
              << "recMu1=" << &muRec[1] << "recMuMS0=" << &muRecMS[0]
              << "recMuMS1=" << &muRecMS[1] << "npix0=" << npix[0]
              << "npix1=" << npix[1] << "nfake0=" << nfake[0]
              << "nfake1=" << nfake[1] << "chi0=" << chiGlo[0]
              << "chi1=" << chiGlo[1] << "chiITS0=" << chiGloITS[0]
              << "chiITS1=" << chiGloITS[1] << "\n";
  }
  outStream.Close();
}
