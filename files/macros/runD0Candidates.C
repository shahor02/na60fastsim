#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "KMCProbeFwd.h"
#include "KMCDetectorFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "./HFUtils.C"
#endif




// settings for signal generation
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double ptminSG = 0.;
double ptmaxSG = 10; //Elena's change, it was 3 GeV/c

double vX = 0, vY = 0, vZ = 0; // event vertex

THnSparseF* CreateSparse();
TDatime dt;
Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk);

void GenerateD0SignalCandidates(Int_t nevents = 100000, 
				double Eint = 160., 
				const char *setup = "setup-10um-itssa_Eff1.txt", 
				const char *filNamPow="/home/prino/na60plus/POWHEG/pp20/Charm1dot5/pp0_frag-PtSpectra-Boost.root",
				const char *privateDecayTable = "/home/prino/na60plus/decaytables/USERTABD0.DEC",
				int optPartAntiPart=3,
				int minITShits=4,
				double chi2Cut = 1.5,
				double minTrackP = 1.,
				bool writeNtuple = kFALSE, 
				bool simulateBg=kTRUE,
				bool optLastLayClean=kFALSE){
  
  // Generate D0->Kpi signals and simulate detector response for decay tracks
  // Need input D0 pt and y ditribution from POWHEG


  int refreshBg = 1000;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);
  gSystem->Load("$ALICE_ROOT/lib/libEvtGen.so");
  gSystem->Load("$ALICE_ROOT/lib/libEvtGenExternal.so");
  gSystem->Load("$ALICE_ROOT/lib/libTEvtGen.so");
  

  //  PYTHIA input -> not used
  // TFile *fin = new TFile("Mergedfkine.root");
  // // TFile *fin = new TFile("fkineNew.root");
  // TH1D *hD0pt = (TH1D *)fin->Get("hD0pt");
  // TH1D *hD0y = (TH1D *)fin->Get("hD0y");
  // hD0y->Rebin(4);
  
  //  POWHEG+PYTHIA input 
  printf("--> pt and y shape of D0 from %s\n",filNamPow);
  TFile *filPow=new TFile(filNamPow);
  TH3D* h3Dpow=(TH3D*)filPow->Get("hptyeta421");
  TH1D *hD0pt = (TH1D*)h3Dpow->ProjectionX("hD0pt");
  TH1D *hD0y = (TH1D*)h3Dpow->ProjectionY("hD0y");
  TH3D* h3Dbarpow=(TH3D*)filPow->Get("hptyetam421");
  if(h3Dbarpow){
    TH1D *hD0barpt = (TH1D*)h3Dbarpow->ProjectionX("hD0barpt");
    TH1D *hD0bary = (TH1D*)h3Dbarpow->ProjectionY("hD0bary");
    if(optPartAntiPart==3){
      hD0pt->Add(hD0barpt);
      hD0y->Add(hD0bary);
    }else if(optPartAntiPart==2){
      hD0pt=hD0barpt;
      hD0y=hD0bary;
    }
  }
  TH2F *hptK = new TH2F("hptK", "kaons from D0 decays", 50,0.,10.,50, 0., 10.);
  TH2F *hptPi = new TH2F("hptPi", "pions from D0 decays", 50, 0.,10.,50,0., 10.);
  TH2F *hmomK = new TH2F("hmomK", "kaons from D0 decays ; #eta ; p (GeV/c)", 50,0.,5.,50, 0., 20.);
  TH2F *hmomPi = new TH2F("hmomPi", "pions from D0 decays ; #eta ; p (GeV/c)", 50, 0.,5.,50,0., 20.);
  TH1D *hyK = new TH1D("hyK", "y kaons from D0 decays", 50, 0., 5.);
  TH1D *hyPi = new TH1D("hyPi", "y pions from D0 decays", 50, 0., 5.);
  TH2F *hyPiK = new TH2F("hyPiK", "y pions vs y Kaons from D0 decays", 50, 0., 5., 50, 0., 5.);
  TFile *fout = new TFile("D0-Signal-histos.root", "recreate");

  //Magnetic field and detector parameters

  //int outN = nev/10;
  //if (outN<1) outN=1;

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT
  if(optLastLayClean) det->setLastITSLayerClean(true);
    
  det->SetMinITSHits(TMath::Min(minITShits,det->GetNumberOfActiveLayersITS())); //NA60+
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
  det->SetMinP2Propagate(0.01); //NA60+
  //det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  det->BookControlHistos();
  //
  
  TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* mag = (MagField*) fld;
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
  }

  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parentgen, daugen[2], parent, daurec[2], parentrefl, daurecswapmass[2]; 
  KMCProbeFwd recProbe[2];  
  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  bool privTab=kFALSE;
  if (strlen(privateDecayTable)>0){
    if(gSystem->Exec(Form("ls -l %s",privateDecayTable))==0){
      fDecayer->SetDecayTablePath((char*)privateDecayTable);
      fDecayer->ReadDecayTable();
      printf("-- Use D0 decay table from file %s\n",privateDecayTable);
      privTab=kTRUE;
    }
  }
  if(!privTab){
    printf("-- Use existing decay modes in aliroot\n");
    fDecayer->SetForceDecay(kHadronicD); 
  }
  fDecayer->ForceDecay();
  
  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  TLorentzVector *mom = new TLorentzVector();
  
  // define mother particle
  Int_t pdgParticle = 421;
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "Y-Pt corr match", 20, 1., 5., 40, ptminSG, ptmaxSG);
  TH1D* hPtGen = new TH1D("hPtGen", "Pt gen", 40, ptminSG, ptmaxSG);
  TH1D* hYGen = new TH1D("hYGen", "Y full phase space", 20, 1., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 20, 1., 5., 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Reconstructed Pt all match", 40, ptminSG, ptmaxSG);  
  TH1D* hPtGenRecoAll = new TH1D("hPtGenRecoAll", "Generated Pt all match", 40, ptminSG, ptmaxSG);
  TH2F* hPtRecoVsGenAll = new TH2F("hPtRecoVsGenAll"," ; Generated p_{T} ; Reconstructed p_{T}",40, ptminSG, ptmaxSG,40, ptminSG, ptmaxSG);
  TH2F* hDiffPtRecoGenAll = new TH2F("hDiffPtRecoGenAll"," ; Generated p_{T} ; Reco p_{T} - Gen p_{T}",40, ptminSG, ptmaxSG,100,-0.2,0.2);

  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Reconstructed Y all match", 20, 1., 5.);
  TH1D* hYGenRecoAll = new TH1D("hYGenRecoAll", "Generated Y all match", 20, 1.,5.);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "Y-Pt fake match", 20,1., 5., 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "Pt fake match", 40, ptminSG, ptmaxSG);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 200, 1., 3.5);
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match", 200, 1., 3.5);
  TH1D* hMassRefl = new TH1D("hMassRefl", "Mass reflections", 200, 1., 3.5);

  TH2D* hKaonDauClu = new TH2D("hKaonDauClu", "N hits kaon daughter ; N hits ; N fake hits ", 11, -0.5, 10.5, 11, -0.5, 10.5);
  TH2D* hPionDauClu = new TH2D("hPionDauClu", "N hits pion daughter ; N hits ; N fake hits ", 11, -0.5, 10.5, 11, -0.5, 10.5);
  TH1D* hPionDauHitPat = new TH1D("hPionDauHitPat", "Hits on layer ; Layer", 11, -0.5, 10.5);
  TH1D* hKaonDauHitPat = new TH1D("hKaonDauHitPat", "Hits on layer ; Layer", 11, -0.5, 10.5);
  TH2D* hDauChi2 = new TH2D("hDauChi2", " ;  N fake hits ; chi2", 11, -0.5, 10.5, 50,0.,5.);
  TH2F *hDistXY = new TH2F("hDistXY", " ", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 300, -1, 1, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", "", 200, 1.5, 2.5, 6, 0, 3);
  TH2F *hMassVsY = new TH2F("hMassVsY", "", 200, 1.5, 2.5, 10, 0, 5);
  TH2F *hMassReflVsPt = new TH2F("hMassReflVsPt", "", 200, 1.5, 2.5, 6, 0, 3);
  TH2F *hMassReflVsY = new TH2F("hMassReflVsY", "", 200, 1.5, 2.5, 10, 0, 5);
  
  TH2F *hResVx = new TH2F("hResVx", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResVxVsY = new TH2F("hResVxVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsYfake00 = new TH2F("hResVzVsYfake00", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsYfake10 = new TH2F("hResVzVsYfake10", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsYfake11 = new TH2F("hResVzVsYfake11", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsYfake5X = new TH2F("hResVzVsYfake5X", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsYfake55 = new TH2F("hResVzVsYfake55", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", "", 100, -1, 1, 50, 0, 5); //for Kaons
  TH2F *hResPyVsY = new TH2F("hResPyVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", "", 2, 0, 2);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntD0cand = 0x0;
  if (writeNtuple)
    {
      fnt = new TFile("fntSig.root", "recreate");
      ntD0cand = new TNtuple("ntD0cand", "ntD0cand", "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax:d0D:costhst", 32000);
    }
  Float_t arrnt[13];
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    if(iev%100==0) printf(" ***************  ev = %d \n", iev);
    int nrec = 0;
    int nfake = 0;
    double pxyz[3];
    
    if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);
    Double_t ptGenD = hD0pt->GetRandom(); // get D0 distribution from file
    Double_t yGenD = hD0y->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGenD = ptGenD * TMath::Cos(phi);
    Double_t pyGenD = ptGenD * TMath::Sin(phi);
    
    Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
    Double_t mt = TMath::Sqrt(ptGenD * ptGenD + mass * mass);
    Double_t pzGenD = mt * TMath::SinH(yGenD);
    Double_t en = mt * TMath::CosH(yGenD);

    Double_t massKaon = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    Double_t massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    mom->SetPxPyPzE(pxGenD, pyGenD, pzGenD, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);
    
    Int_t arrpdgdau[2];
    Double_t ptK = -999.;
    Double_t ptPi = -999.;
    Double_t yK=-999.;
    Double_t yPi = -999.;
    Int_t icount = 0;
    Double_t secvertgenK[3]={0.,0.,0.};
    Double_t secvertgenPi[3]={0.,0.,0.};
    
    // loop on decay products
    for (int i = 0; i < np; i++) { 
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = TMath::Abs(iparticle1->GetPdgCode());
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
      if (kf == pdgParticle){
	// D0 particle
	hYGen->Fill(iparticle1->Y());
	hPtGen->Fill(iparticle1->Pt());
	hYPtGen->Fill(iparticle1->Y(), iparticle1->Pt());
	//printf("mother part = %d, code=%d, pt=%f, y=%f \n", i, kf, iparticle1->Pt(), iparticle1->Y());
      }	
      if (kf == 321 || kf == 211){
	// daughters that can be reconstructed: K and pi
	Double_t e = iparticle1->Energy();
	Double_t px = iparticle1->Px();
	Double_t py = iparticle1->Py();
	Double_t pz = iparticle1->Pz();	    
	TLorentzVector *pDecDau = new TLorentzVector(0., 0., 0., 0.);
	pDecDau->SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	Int_t crg=1;
	if(iparticle1->GetPdgCode()<0) crg=-1;
	if (!det->SolveSingleTrack(pDecDau->Pt(), pDecDau->Rapidity(), pDecDau->Phi(), iparticle1->GetMass(), crg, vX, vY, vZ, 0, 1, 99)) continue;
	KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
	if (!trw) continue;
	if (trw->GetNormChi2(kTRUE) > chi2Cut) continue;
	if (trw->GetP() < minTrackP) continue;
	nrec++;
	    
	nfake += trw->GetNFakeITSHits();
	trw->GetPXYZ(pxyz);
	if (kf == 321){
	  // Kaon daughter
	  ptK = iparticle1->Pt();
	  yK = iparticle1->Y();
	  hptK->Fill(ptGenD,ptK);
	  hmomK->Fill(iparticle1->Eta(),iparticle1->P());
	  hyK->Fill(iparticle1->Y());
	  secvertgenK[0] = iparticle1->Vx();
	  secvertgenK[1] = iparticle1->Vy();
	  secvertgenK[2] = iparticle1->Vz();
	  daugen[0].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
	  daurecswapmass[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massPion);
	  recProbe[0] = *trw;
	}else if (kf == 211){
	  // Pion daughter
	  ptPi = iparticle1->Pt();
	  yPi = iparticle1->Y();
	  hptPi->Fill(ptGenD,ptPi);
	  hmomPi->Fill(iparticle1->Eta(),iparticle1->P());
	  hyPi->Fill(iparticle1->Y());
	  secvertgenPi[0] = iparticle1->Vx();
	  secvertgenPi[1] = iparticle1->Vy();
	  secvertgenPi[2] = iparticle1->Vz();
	  daugen[1].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],  iparticle1->GetMass());
	  daurecswapmass[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2],massKaon);
	  recProbe[1] = *trw;
	}
      }
    }
    if (ptK > 0 && ptPi > 0) hyPiK->Fill(yPi, yK);
    if (nrec < 2) continue;
    hNevents->Fill(1.5);
    
    recProbe[0].PropagateToDCA(&recProbe[1]);
    
    parent = daurec[0];
    parent += daurec[1];
    parentgen = daugen[0];
    parentgen += daugen[1];
    parentrefl = daurecswapmass[0];
    parentrefl += daurecswapmass[1];
      
    Double_t  ptRecD=parent.Pt();
    Double_t  massRecD=parent.M();
    Double_t  massRecReflD=parentrefl.M();
    Double_t yRecD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
    hYPtRecoAll->Fill(yRecD, ptRecD);
    hPtRecoAll->Fill(ptRecD);
    hPtGenRecoAll->Fill(ptGenD);
    hPtRecoVsGenAll->Fill(ptGenD,ptRecD);
    hDiffPtRecoGenAll->Fill(ptGenD,(ptRecD-ptGenD));
    hYRecoAll->Fill(yRecD);
    hYGenRecoAll->Fill(yGenD);
    hMassAll->Fill(massRecD);
    hMassRefl->Fill(massRecReflD);
    if (nfake > 0){
      hYPtRecoFake->Fill(yRecD, ptRecD);
      hPtRecoFake->Fill(ptRecD);
      hMassFake->Fill(massRecD);
    }
    hMassVsPt->Fill(massRecD,ptRecD);
    hMassVsY->Fill(massRecD,yRecD);
    hMassReflVsPt->Fill(massRecReflD,ptRecD);
    hMassReflVsY->Fill(massRecReflD,yRecD);
    hDauChi2->Fill(recProbe[0].GetNFakeITSHits(),recProbe[0].GetNormChi2(kTRUE));
    hDauChi2->Fill(recProbe[1].GetNFakeITSHits(),recProbe[1].GetNormChi2(kTRUE));
    hKaonDauClu->Fill(recProbe[0].GetNHits(),recProbe[0].GetNFakeITSHits());
    hPionDauClu->Fill(recProbe[1].GetNHits(),recProbe[1].GetNFakeITSHits());
    UInt_t hmap0=recProbe[0].GetHitsPatt();
    UInt_t hmap1=recProbe[1].GetHitsPatt();
    for(Int_t jl=0; jl<=10; jl++){
      if( hmap0 & (1<<jl) ) hKaonDauHitPat->Fill(jl);
      if( hmap1 & (1<<jl) ) hPionDauHitPat->Fill(jl);
    }
    Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
    Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
    Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
    Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
    hDCA->Fill(dca, ptRecD);
    hDCAx->Fill(d1, ptRecD);
    hDCAy->Fill(d2, ptRecD);
    hDCAz->Fill(d3, ptRecD);
    
    //      Double_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
    //      Double_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
    //      Double_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;
    
    Double_t xP, yP, zP;
    ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
    Double_t residVx=10000.*(xP - secvertgenK[0]);
    Double_t residVy=10000.*(yP - secvertgenK[1]);
    Double_t residVz=10000.*(zP - secvertgenK[2]);
    hResVx->Fill(residVx, ptRecD);
    hResVy->Fill(residVy, ptRecD);
    hResVz->Fill(residVz, ptRecD);
    hResVxVsY->Fill(residVx, yRecD);
    hResVyVsY->Fill(residVy, yRecD);
    hResVzVsY->Fill(residVz, yRecD);
    if(recProbe[0].GetNFakeITSHits()==0 && recProbe[1].GetNFakeITSHits()==0)  hResVzVsYfake00->Fill(residVz, yRecD);
    else if((recProbe[0].GetNFakeITSHits()==1 && recProbe[1].GetNFakeITSHits()==0) || (recProbe[1].GetNFakeITSHits()==1 && recProbe[0].GetNFakeITSHits()==0))  hResVzVsYfake10->Fill(residVz, yRecD);
    else if(recProbe[0].GetNFakeITSHits()==1 && recProbe[1].GetNFakeITSHits()==1)  hResVzVsYfake11->Fill(residVz, yRecD);
    else if((recProbe[0].GetNFakeITSHits()==5 && recProbe[1].GetNFakeITSHits()<5) || (recProbe[1].GetNFakeITSHits()==5 && recProbe[0].GetNFakeITSHits()<5))  hResVzVsYfake5X->Fill(residVz, yRecD);
    else if(recProbe[0].GetNFakeITSHits()==5 && recProbe[1].GetNFakeITSHits()==5)  hResVzVsYfake55->Fill(residVz, yRecD);

    hResPx->Fill(daurec[0].Px() - daugen[0].Px(), ptRecD);
    hResPy->Fill(daurec[0].Py() - daugen[0].Py(), ptRecD);
    hResPz->Fill(daurec[0].Pz() - daugen[0].Pz(), ptRecD);
    hResPx->Fill(daurec[1].Px() - daugen[1].Px(), ptRecD);
    hResPy->Fill(daurec[1].Py() - daugen[1].Py(), ptRecD);
    hResPz->Fill(daurec[1].Pz() - daugen[1].Pz(), ptRecD);

    hResPxVsY->Fill(daurec[0].Px() - daugen[0].Px(), yRecD);
    hResPyVsY->Fill(daurec[0].Py() - daugen[0].Py(), yRecD);
    hResPzVsY->Fill(daurec[0].Pz() - daugen[0].Pz(), yRecD);
    hResPxVsY->Fill(daurec[1].Px() - daugen[1].Px(), yRecD);
    hResPyVsY->Fill(daurec[1].Py() - daugen[1].Py(), yRecD);
    hResPzVsY->Fill(daurec[1].Pz() - daugen[1].Pz(), yRecD);
    
    // cout << "secvert generated Pion: " << secvertgenPi[0] << "  " << secvertgenPi[1] << "  " << secvertgenPi[2] << endl;
    // cout << "Reco Vert  Pion: " << xP << "  " << yP << "  " << zP << endl;
    
    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
    Float_t distgen = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1] + secvertgenPi[2] * secvertgenPi[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1]);
    // printf("dist = %f , distXY=%f , dx=%f, dy=%f, dz=%f z1=%f, z2=%f \n", dist, distXY, xP, yP, zP, recProbe[0].GetZ(), recProbe[1].GetZ());
    // printf("distgen = %f , distgenXY=%f \n", distgen, distgenXY);
    
    Double_t vsec[3] = {xP, yP, zP};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t cts = CosThetaStar(parent,daurec[0]);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRecD);
    // printf(" ***** ***** cos point = %f \n", cosp);
    //if (cosp < -0.98)
    //    printf("SMALL COSPOINT");
    
    hResDist->Fill(dist - distgen, ptRecD);
    hResDistXY->Fill(distXY - distgenXY, ptRecD);
    
    //recProbe[0].PropagateToDCA(&recProbe[1]);
    
    // hYPtAll->Fill(parent.Y(), ptRecD);
    // hPtAll->Fill(ptRecD);
    hDistXY->Fill(distXY, ptRecD);
    hDist->Fill(dist, ptRecD);
    hDistgenXY->Fill(distgenXY, ptRecD);
    hDistgen->Fill(distgen, ptRecD);
      
    //AliExternalTrackParam *track1 = (AliExternalTrackParam *)recProbe[0].GetTrack();
    //AliExternalTrackParam *track2 = (AliExternalTrackParam *)recProbe[1].GetTrack();
    recProbe[0].PropagateToZBxByBz(0);
    Double_t d0x1 = recProbe[0].GetX();
    Double_t d0y1 = recProbe[0].GetY();
    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
    if (d0x1 < 0)
      d0xy1 *= -1;
    
    recProbe[1].PropagateToZBxByBz(0);
    Double_t d0x2 = recProbe[1].GetX();
    Double_t d0y2 = recProbe[1].GetY();
    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
    if (d0x2 < 0)
      d0xy2 *= -1;
    
    // printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
    
    hd0XYprod->Fill(d0xy1 * d0xy2, ptRecD);
    hd0XY1->Fill(d0xy1, ptRecD);
    hd0XY2->Fill(d0xy2, ptRecD);
      
    arrsp[0] = massRecD;
    arrsp[1] = ptRecD;
    arrsp[2] = yRecD;
    arrsp[3] = dist;
    arrsp[4] = cosp;
    arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
    arrsp[6] = d0xy1 * d0xy2;
    arrsp[7] = dca;
    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
    arrsp[9] = TMath::Abs(ipD);	    
    arrsp[10] = cts;      
    hsp->Fill(arrsp);
    
    if (ntD0cand){
      arrnt[0] = massRecD;
      arrnt[1] = ptRecD;
      arrnt[2] = yRecD;
      arrnt[3] = dist;
      arrnt[4] = cosp;
      arrnt[5] = d0xy1;
      arrnt[6] = d0xy2;
      arrnt[7] = d0xy1 * d0xy2;
      arrnt[8] = dca;
      arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      arrnt[11] = TMath::Abs(ipD);
      arrnt[12] = cts;
      ntD0cand->Fill(arrnt);
    }
  } //event loop
  
  hMassAll->SetLineColor(kBlue);
  hMassAll->Draw();
  hMassAll->SetMinimum(0.1);
  hMassFake->SetLineColor(kRed);
  hMassFake->Draw("same");
  
  fout->cd();  
  hMassAll->Write();
  hMassFake->Write();
  hMassRefl->Write();
  hMassVsPt->Write();
  hMassVsY->Write();
  hMassReflVsPt->Write();
  hMassReflVsY->Write();
  hYPtGen->Write();
  hPtGen->Write();
  hYGen->Write();
  hYPtRecoAll->Write();  
  hYPtRecoFake->Write();
  hPtRecoAll->Write();
  hPtGenRecoAll->Write();
  hPtRecoVsGenAll->Write();
  hDiffPtRecoGenAll->Write();
  hYRecoAll->Write();
  hYGenRecoAll->Write();
  hPtRecoFake->Write();
  hDauChi2->Write();
  hKaonDauClu->Write();
  hPionDauClu->Write();
  hKaonDauHitPat->Write();
  hPionDauHitPat->Write();
  hDistXY->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hResVx->Write();
  hResVy->Write();
  hResVz->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResVzVsYfake00->Write();
  hResVzVsYfake10->Write();
  hResVzVsYfake11->Write();
  hResVzVsYfake5X->Write();
  hResVzVsYfake55->Write();

  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hNevents->Write();
  hsp->Write();
  if (ntD0cand){
    fnt->cd();
    ntD0cand->Write();
    fnt->Close();
  }
  // TCanvas *ccdau = new TCanvas();
  // ccdau->Divide(3, 2);
  // ccdau->cd(1)->SetLogy();
  // hPtGen->Draw();
  // ccdau->cd(2)->SetLogy();
  // hptK->Draw();
  // ccdau->cd(3)->SetLogy();
  // hptPi->Draw();
  // ccdau->cd(4)->SetLogy();
  // hYGen->Draw();
  // ccdau->cd(5)->SetLogy();
  // hyK->Draw();
  // ccdau->cd(6)->SetLogy();
  // hyPi->Draw();
  
  TFile fout2("DecayHistos.root", "RECREATE");
  // TFile fout2("DecayHistostest.root", "RECREATE");
  hPtGen->Write();
  hptK->Write();
  hptPi->Write();
  hmomK->Write();
  hmomPi->Write();
  hyK->Write();
  hyPi->Write();
  hYGen->Write();
  hyPiK->Write();
  fout2.Close();

  fout->Close();
}



void MakeD0CombinBkgCandidates(const char *setup = "setup-10um-itssa_Eff1.txt",
			       const char* trackTreeFile="treeBkgEvents.root",
			       Int_t nevents = 999999,
			       int minITShits=4,
			       double chi2Cut = 1.5,
			       double minTrackP = 1.,
			       Int_t writeNtuple = kFALSE){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create D0 combinatorial background candidates (= OS pairs of tracks)
  // Store in THnSparse and (optionally) TNtuple

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  
  TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* mag = (MagField*) fld;
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
  }

  TFile *filetree = new TFile(trackTreeFile);
  TTree *tree = (TTree *)filetree->Get("tree");
  TClonesArray *arr = 0;
  tree->SetBranchAddress("tracks", &arr);
  Int_t entries = tree->GetEntries();
  printf("Number of events in tree = %d\n",entries);
  if(nevents>entries) nevents=entries;
  else printf(" --> Analysis performed on first %d events\n",nevents);

  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  // define mother particle
  Int_t pdgParticle = 421;
  
  TFile *fout = new TFile("D0-Bkg-histos.root", "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 0., 2.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDCAx = new TH2F("hDCAx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAy = new TH2F("hDCAy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hDCAz = new TH2F("hDCAz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XYprod = new TH2F("hd0xyprod", "", 100, -0.01, 0.01, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  
  TH1D *hVx = new TH1D("hVx", "", 200, -0.05, 0.05);
  TH1D *hVy = new TH1D("hVy", "", 200, -0.05, 0.05);
  TH1D *hVz = new TH1D("hVz", "", 200, -0.05, 0.05);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  TH2F *hCospXY = new TH2F("hCospXY", "", 100, -1, 1, 30, 0, 3);
  TH2F *hCosThStVsMass = new TH2F("hCosThStVsMass", "", 50, 1.5, 2.5, 40, -1, 1);
  
  TH2F *hmomK = new TH2F("hmomK", "kaons from D0 decays ; #eta ; p (GeV/c)", 50,0.,5.,50, 0., 20.);
  TH2F *hmomPi = new TH2F("hmomPi", "pions from D0 decays ; #eta ; p (GeV/c)", 50, 0.,5.,50,0., 20.);

  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
    
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntD0cand = 0x0;
  Float_t arrnt[13];
  if (writeNtuple){
    fnt = new TFile("fntBkg.root", "recreate");
    ntD0cand = new TNtuple("ntD0cand", "ntD0cand",  "mass:pt:y:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax:d0D:costhst", 32000);
  }

  KMCProbeFwd recProbe[2];
  TLorentzVector parent, daurec[2];
  Double_t massD0 = 1.864;
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    tree->GetEvent(iev);
    Int_t arrentr = arr->GetEntriesFast();
    
    for (Int_t itr = 0; itr < arrentr; itr++){
      KMCProbeFwd *tr1 = (KMCProbeFwd *)arr->At(itr);
      // cout << "tr P=" << tr1->GetP() << endl;
      if(tr1->GetNHits()<minITShits) continue;
      if (tr1->GetNormChi2(kTRUE) > chi2Cut) continue;
      if (tr1->GetP() < minTrackP) continue;
      Float_t ch1 = tr1->GetCharge();
      for (Int_t itr2 = itr; itr2 < arrentr; itr2++){
	KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
	if(tr2->GetNHits()<minITShits) continue;
	if (tr2->GetNormChi2(kTRUE) > chi2Cut) continue;
	if (tr2->GetP() < minTrackP) continue;
	Float_t ch2 = tr2->GetCharge();
	if (ch1 * ch2 > 0) continue;
	if (ch1 < 0){ //convention: first track negative
	  recProbe[0] = *tr1;
	  recProbe[1] = *tr2;
	}else if (ch2 < 0){
	  recProbe[0] = *tr2;
	  recProbe[1] = *tr1;
	}
	Bool_t ok=recProbe[0].PropagateToDCA(&recProbe[1]);
	Double_t pxyz[3];
	recProbe[0].GetPXYZ(pxyz);
	
	Double_t pxyz2[3];
	recProbe[1].GetPXYZ(pxyz2);

	Double_t xP, yP, zP;
	ComputeVertex(recProbe[0],recProbe[1],xP,yP,zP);
	
	Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
	Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
	Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();

	recProbe[0].PropagateToZBxByBz(0);
	recProbe[1].PropagateToZBxByBz(0);
	
	for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
	  // mass hypothesis: Kpi, piK
	  Int_t iKaon=-1;
	  Int_t iPion=-1;
	  if(iMassHyp==0){
	    daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], KMCDetectorFwd::kMassK);
	    daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassPi);
	    iKaon=0;
	    iPion=1;
	  }else{
	    daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2], KMCDetectorFwd::kMassPi);
	    daurec[1].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassK);
	    iKaon=1;
	    iPion=0;
	  }
	  parent = daurec[0];
	  parent += daurec[1];
	  countCand++;
	  Float_t ptD=parent.Pt();
	  Float_t invMassD=parent.M();
	  Float_t yD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
	  hYPtRecoAll->Fill(yD, ptD);
	  hPtRecoAll->Fill(ptD);
	  hMassAll->Fill(invMassD);
	  if(invMassD>1.65  && invMassD<2.15){
	    // range to fill histos
	    if(invMassD>1.805 && invMassD<1.925) countCandInPeak++;

	    Float_t dca = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
	    
	    //printf(" DCA = %f\n", sqrt(d1 * d1 + d2 * d2 + d3 * d3));
	    hDCA->Fill(dca, ptD);
	    hDCAx->Fill(d1, ptD);
	    hDCAy->Fill(d2, ptD);
	    hDCAz->Fill(d3, ptD);
	    // Float_t xP = (recProbe[1].GetX() + recProbe[0].GetX()) / 2.;
	    // Float_t yP = (recProbe[1].GetY() + recProbe[0].GetY()) / 2.;
	    // Float_t zP = (recProbe[1].GetZ() + recProbe[0].GetZ()) / 2.;

	    hVx->Fill(xP);
	    hVy->Fill(yP);
	    hVz->Fill(zP);

	    Float_t dist = TMath::Sqrt(xP * xP + yP * yP + zP * zP);
	    Float_t distXY = TMath::Sqrt(xP * xP + yP * yP);
	    Double_t vsec[3] = {xP, yP, zP};
	    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
	    Double_t cospxy = CosPointingAngleXY(vprim, vsec, parent);
	    Double_t cts = CosThetaStar(parent,daurec[iKaon]);
	    Double_t ipD = ImpParXY(vprim, vsec, parent);
	    Double_t momK=daurec[iKaon].P();
	    Double_t momPi=daurec[iPion].P();
	    Double_t etaK=daurec[iKaon].Eta();
	    Double_t etaPi=daurec[iPion].Eta();
	    hmomK->Fill(etaK,momK);
	    hmomPi->Fill(etaPi,momPi);
	    hCosp->Fill(cosp, ptD);
	    hCospXY->Fill(cospxy, ptD);
	    hCosThStVsMass->Fill(invMassD,cts);
	    //	    printf(" ***** ***** cos point = %f vprim=%f %f %f  vsec=%f %f %f\n", cosp,vprim[0],vprim[1],vprim[2],vsec[0],vsec[1],vsec[2]);
	    hDistXY->Fill(distXY, ptD);
	    hDist->Fill(dist, ptD);
    
	    Double_t d0x1 = recProbe[0].GetX();
	    Double_t d0y1 = recProbe[0].GetY();
	    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
	    if (d0x1 < 0) d0xy1 *= -1;
	    
	    Double_t d0x2 = recProbe[1].GetX();
	    Double_t d0y2 = recProbe[1].GetY();
	    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
	    if (d0x2 < 0) d0xy2 *= -1;
	  	      
	    //printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
	    
	    hd0XYprod->Fill(d0xy1 * d0xy2, ptD);
	    hd0XY1->Fill(d0xy1, ptD);
	    hd0XY2->Fill(d0xy2, ptD);
	    if(cosp>0.97 && (d0xy1*d0xy2)<0.0001){
	      arrsp[0] = invMassD;
	      arrsp[1] = ptD;
	      arrsp[2] = yD;
	      arrsp[3] = dist;
	      arrsp[4] = cosp;
	      arrsp[5] = TMath::Min(TMath::Abs(d0xy1),TMath::Abs(d0xy2));
	      arrsp[6] = d0xy1 * d0xy2;
	      arrsp[7] = dca;
	      arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      arrsp[9] = TMath::Abs(ipD);	    
	      arrsp[10] = cts;
	      hsp->Fill(arrsp);
	    }
	    if (ntD0cand && cosp>0.97 && (d0xy1*d0xy2)<0.0001){
	      arrnt[0] = invMassD;
	      arrnt[1] = ptD;
	      arrnt[2] = yD;
	      arrnt[3] = dist;
	      arrnt[4] = cosp;
	      arrnt[5] = d0xy1;
	      arrnt[6] = d0xy2;
	      arrnt[7] = d0xy1 * d0xy2;
	      arrnt[8] = dca;
	      arrnt[9] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      arrnt[10] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      arrnt[11] = TMath::Abs(ipD);	    
	      arrnt[12] = cts;
	      ntD0cand->Fill(arrnt);
	    }
	    
	  } // check on inv mass
	} // loop on mass hypothesis
      } // loop on first track
    } // loop on second track
    
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot D0 candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hmomK->Write();
  hmomPi->Write();
  hDistXY->Write();
  hDist->Write();
  hDCA->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAz->Write();
  hVx->Write();
  hVy->Write();
  hVz->Write();
  hCosp->Write();
  hCospXY->Write();
  hCosThStVsMass->Write();
  hd0XYprod->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hsp->Write();
  fout->Close();
  if (ntD0cand){
    fnt->cd();
    ntD0cand->Write();
    fnt->Close();
  }
}


Double_t CosThetaStar(TLorentzVector &parent, TLorentzVector &dauk) {

  Double_t massMoth = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t massK = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  Double_t pStar = TMath::Sqrt((massMoth*massMoth-massK*massK-massPi*massPi)*(massMoth*massMoth-massK*massK-massPi*massPi)-4.*massK*massK*massPi*massPi)/(2.*massMoth);

  Double_t pMoth=parent.P();
  Double_t e=TMath::Sqrt(massMoth*massMoth+pMoth*pMoth);
  Double_t beta = pMoth/e;
  Double_t gamma = e/massMoth;
  TVector3 momDau(dauk.Px(),dauk.Py(),dauk.Pz());
  TVector3 momMoth(parent.Px(),parent.Py(),parent.Pz());
  Double_t qlProng=momDau.Dot(momMoth)/momMoth.Mag();
  Double_t cts = (qlProng/gamma-beta*TMath::Sqrt(pStar*pStar+massK*massK))/pStar;

  return cts;
}

THnSparseF* CreateSparse(){
  const Int_t nAxes=11;
  TString axTit[nAxes]={"Inv. mass (GeV/c^{2})","p_{T} (GeV/c)","y",
			"Dec Len (cm)","cos(#vartheta_{p})",
			"d_0^{min} (cm)",
			"d_0*d_0 (cm^{2})","DCA",
			"p_{T}^{min} (GeV/c)",
			"d_0^{D} (cm)","cos(#theta*)"};
  Int_t bins[nAxes] =   {100,   10, 20, 30,  20,   10,   10,      12,    6,  16,  10}; 
  Double_t min[nAxes] = {1.65,  0., 1., 0., 0.98, 0.,   -0.0006,  0.0,   0.,  0.,  -1.};
  Double_t max[nAxes] = {2.15,  5., 5., 0.3, 1.,   0.05, 0.,      0.03,  3.,  0.04, 1.};  
  THnSparseF *hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
  for(Int_t iax=0; iax<nAxes; iax++) hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  return hsp;
}




