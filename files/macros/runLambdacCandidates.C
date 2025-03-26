#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH1D.h>
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
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "./HFUtils.C"
#endif

// Track cuts
double ChiTot = 1.5;
int minITShits=4;
double mind0xyDau=-9999;
double cutCosPointCand=-1.01;
double cutImpParProd=99999.;

// settings for signal generation
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double ptminSG = 0.;
double ptmaxSG = 10; //Elena's change, it was 3 GeV/c

double vX = 0, vY = 0, vZ = 0; // event vertex

// THnSparse axis config
Int_t nBinsDecLen=20;
Double_t minDecLen=0.;
Double_t maxDecLen=0.1;
Int_t nBinsDecLenXY=16;
Double_t minDecLenXY=0.;
Double_t maxDecLenXY=0.04;
Int_t nBinsCosPoi=20;
Double_t minCosPoi=0.98;
Double_t maxCosPoi=1.;
Int_t nBinsImpDau=10;
Double_t minImpDau=0.;
Double_t maxImpDau=0.05;
Int_t nBinsImpPro=20;
Double_t minImpPro=-0.002;
Double_t maxImpPro=0.0005;
Int_t nBinsSigVer=12;
Double_t minSigVer=0.;
Double_t maxSigVer=0.03;
Int_t nBinsPrat=7;
Double_t minPrat=0.;
Double_t maxPrat=1.4;
Int_t nBinsImpLc=12;
Double_t minImpLc=0.;
Double_t maxImpLc=0.03;

THnSparseF* CreateSparse();
void ConfigureSelectionsAndAxes(const char *selectionFile);
TDatime dt;

void GenerateLambdacSignalCandidates(Int_t nevents = 100000, 
				     double Eint = 160., 
				     const char *setup = "setup-10um-itssa_Eff1.txt", 
				     const char *filNamPow="/home/prino/cernbox/na60plus/POWHEG/pp20/Charm1dot5/pp0_frag-PtSpectra-Boost.root", 
				     const char *privateDecayTable = "../decaytables/USERTABLC.DEC",
				     int optPartAntiPart=3,
				     const char *selectionFile="",
				     bool writeNtuple = kFALSE, 
				     bool simulateBg=kTRUE){
  
  // Generate Lc->pKpi signals and simulate detector response for decay tracks
  // Need input Lc pt and y ditribution from POWHEG


  int refreshBg = 10000;
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
  printf("--> pt and y shape of Lc from %s\n",filNamPow);
  TFile *filPow=new TFile(filNamPow);
  TH3D* h3Dpow=(TH3D*)filPow->Get("hptyeta4122");
  TH1D *hLcpt = (TH1D*)h3Dpow->ProjectionX("hLcpt");
  TH1D *hLcy = (TH1D*)h3Dpow->ProjectionY("hLcy");
  TH3D* h3Dbarpow=(TH3D*)filPow->Get("hptyetam4122");
  if(h3Dbarpow){
    TH1D *hLcbarpt = (TH1D*)h3Dbarpow->ProjectionX("hLcbarpt");
    TH1D *hLcbary = (TH1D*)h3Dbarpow->ProjectionY("hLcbary");
    if(optPartAntiPart==3){
      hLcpt->Add(hLcbarpt);
      hLcy->Add(hLcbary);
    }else if(optPartAntiPart==2){
      hLcpt=hLcbarpt;
      hLcy=hLcbary;
    }
  }

  ConfigureSelectionsAndAxes(selectionFile);

  TH2F *hptK = new TH2F("hptK", "kaons from Lc decays", 50,0.,10.,50, 0., 10.);
  TH2F *hptPi = new TH2F("hptPi", "pions from Lc decays", 50, 0.,10.,50,0., 10.);
  TH2F *hptP = new TH2F("hptP", "protons from Lc decays", 50, 0.,10.,50,0., 10.);
  TH1D *hyK = new TH1D("hyK", "y kaons from Lc decays", 50, 0., 5.);
  TH1D *hyPi = new TH1D("hyPi", "y pions from Lc decays", 50, 0., 5.);
  TH1D *hyP = new TH1D("hyP", "y protons from Lc decays", 50, 0., 5.);
  TH2F *hyPiK = new TH2F("hyPiK", "y pions vs y Kaons from Lc decays", 50, 0., 5., 50, 0., 5.);
  
  TFile *fout = new TFile("Lc-Signal-histos.root", "recreate");
  
  //Magnetic field and detector parameters

  //int outN = nev/10;
  //if (outN<1) outN=1;

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

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
  TLorentzVector parentgen, daugen[3], parent, daurec[3]; 
  Double_t daumass[3];
  KMCProbeFwd recProbe[3];  
  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  bool privTab=kFALSE;
  if (strlen(privateDecayTable)>0){
    if(gSystem->Exec(Form("ls -l %s",privateDecayTable))==0){
      fDecayer->SetDecayTablePath((char*)privateDecayTable);
      fDecayer->ReadDecayTable();
      printf("-- Use Lc decay table from file %s\n",privateDecayTable);
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
  Int_t pdgParticle = 4122;
  
  TH2F* hYPtGen = new TH2F("hYPtGen", "Y-Pt corr match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtGen = new TH1D("hPtGen", "Pt gen", 40, ptminSG, ptmaxSG);
  TH1D* hYGen = new TH1D("hYGen", "Y full phase space", 80., 1., 5.4);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 40, ptminSG, ptmaxSG);  
  TH1D* hYRecoAll = new TH1D("hYRecoAll", "Y all match", 80., 1., 5.4);
  TH2F* hYPtRecoFake = new TH2F("hYPtRecoFake", "Y-Pt fake match", 80, 1.0, 5.4, 40, ptminSG, ptmaxSG);
  TH1D* hPtRecoFake = new TH1D("hPtRecoFake", "Pt fake match", 40, ptminSG, ptmaxSG);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 200, 1., 3.5);
  TH1D* hMassFake = new TH1D("hMassFake", "Mass fake match", 200, 1., 3.5);

  TH2D* hDauClu = new TH2D("hDauClu", "N hits kaon daughter ; N hits ; N fake hits", 11, -0.5, 10.5, 11, -0.5, 10.5);
  TH1D* hDauHitPat = new TH1D("hDauHitPat", "Hits on layer ; Layer", 11, -0.5, 10.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistZ = new TH2F("hDistZ", "", 100, 0, 0.2, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 300, -1, 1, 30, 0, 3);
  TH2F *hCospXY = new TH2F("hCospXY", "", 100, -1, 1, 30, 0, 3);
  TH2F *hSigVert = new TH2F("hSigVert", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY3 = new TH2F("hd0xy3", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hMassVsPt = new TH2F("hMassVsPt", "", 200, 1.5, 2.5, 6, 0, 3);
  TH2F *hMassVsY = new TH2F("hMassVsY", "", 200, 1.5, 2.5, 10, 0, 5);
  
  TH2F *hResVx = new TH2F("hResVx", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 200, -1000., 1000., 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResVxVsY = new TH2F("hResVxVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVyVsY = new TH2F("hResVyVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResVzVsY = new TH2F("hResVzVsY", "", 200, -1000., 1000., 50, 0, 5);
  TH2F *hResPxVsY = new TH2F("hResPxVsY", "", 100, -1, 1, 50, 0, 5); //for Kaons
  TH2F *hResPyVsY = new TH2F("hResPyVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResPzVsY = new TH2F("hResPzVsY", "", 100, -1, 1, 50, 0, 5);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntLccand = 0x0;
  if (writeNtuple){
    fnt = new TFile("Lc-Signal-ntuple.root", "recreate");
    ntLccand = new TNtuple("ntLccand", "ntLccand", "mass:pt:y:dist:cosp:d01:d02:d03:sigvert:pP:pK:pPi:d0Lc", 32000);
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
    Double_t ptGenD = hLcpt->GetRandom(); // get Lc distribution from file
    Double_t yGenD = hLcy->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t pxGenD = ptGenD * TMath::Cos(phi);
    Double_t pyGenD = ptGenD * TMath::Sin(phi);
    
    Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();
    Double_t mt = TMath::Sqrt(ptGenD * ptGenD + mass * mass);
    Double_t pzGenD = mt * TMath::SinH(yGenD);
    Double_t en = mt * TMath::CosH(yGenD);
    
    mom->SetPxPyPzE(pxGenD, pyGenD, pzGenD, en);
    Int_t np;
    do{
      fDecayer->Decay(pdgParticle, mom);
      np = fDecayer->ImportParticles(particles);
    } while (np < 0);
    
    Int_t arrpdgdau[3];
    Double_t ptK = -999.;
    Double_t ptPi = -999.;
    Double_t ptP = -999.;
    Double_t momK = -999.;
    Double_t momPi = -999.;
    Double_t momP = -999.;
    Double_t yK=-999.;
    Double_t yPi = -999.;
    Double_t yP = -999.;
    Int_t icount = 0;
    Double_t secvertgenK[3]={0.,0.,0.};
    Double_t secvertgenPi[3]={0.,0.,0.};
    Double_t secvertgenP[3]={0.,0.,0.};
    
    // loop on decay products
    for (int i = 0; i < np; i++) { 
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = TMath::Abs(iparticle1->GetPdgCode());
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
      if (kf == pdgParticle){
	// Lc particle
	hYGen->Fill(iparticle1->Y());
	hPtGen->Fill(iparticle1->Pt());
	hYPtGen->Fill(iparticle1->Y(), iparticle1->Pt());
	//printf("mother part = %d, code=%d, pt=%f, y=%f \n", i, kf, iparticle1->Pt(), iparticle1->Y());
      }	
      if (kf == 2212 || kf == 321 || kf == 211){
	// daughters that can be reconstructed: p, K and pi
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
	if (trw->GetNormChi2(kTRUE) > ChiTot) continue;
	nrec++;
	    
	nfake += trw->GetNFakeITSHits();
	if (kf == 2212){
	  // Proton daughter
	  ptP = iparticle1->Pt();
	  momP = iparticle1->P();
	  yP = iparticle1->Y();
	  hptP->Fill(ptGenD,ptP);
	  hyP->Fill(yP);
	  secvertgenP[0] = iparticle1->Vx();
	  secvertgenP[1] = iparticle1->Vy();
	  secvertgenP[2] = iparticle1->Vz();
	  daugen[0].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daumass[0] = iparticle1->GetMass();
	  recProbe[0] = *trw;
	}else if (kf == 321){
	  // Kaon daughter
	  ptK = iparticle1->Pt();
	  momK = iparticle1->P();
	  yK = iparticle1->Y();
	  hptK->Fill(ptGenD,ptK);
	  hyK->Fill(yK);
	  secvertgenK[0] = iparticle1->Vx();
	  secvertgenK[1] = iparticle1->Vy();
	  secvertgenK[2] = iparticle1->Vz();
	  daugen[1].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daumass[1] = iparticle1->GetMass();
	  recProbe[1] = *trw;
	}else if (kf == 211){
	  // Pion daughter
	  ptPi = iparticle1->Pt();
	  momPi = iparticle1->P();
	  yPi = iparticle1->Y();
	  hptPi->Fill(ptGenD,ptPi);
	  hyPi->Fill(yPi);
	  secvertgenPi[0] = iparticle1->Vx();
	  secvertgenPi[1] = iparticle1->Vy();
	  secvertgenPi[2] = iparticle1->Vz();
	  daugen[2].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daumass[2] = iparticle1->GetMass();
	  recProbe[2] = *trw;
	}
      }
    }
    if (ptK > 0 && ptPi > 0) hyPiK->Fill(yPi, yK);
    if (nrec < 3) continue;
    
    hDauClu->Fill(recProbe[0].GetNHits(),recProbe[0].GetNFakeITSHits());
    hDauClu->Fill(recProbe[1].GetNHits(),recProbe[1].GetNFakeITSHits());
    hDauClu->Fill(recProbe[2].GetNHits(),recProbe[2].GetNFakeITSHits());
    UInt_t hmap0=recProbe[0].GetHitsPatt();
    UInt_t hmap1=recProbe[1].GetHitsPatt();
    UInt_t hmap2=recProbe[2].GetHitsPatt();
    for(Int_t jl=0; jl<=10; jl++){
      if( hmap0 & (1<<jl) ) hDauHitPat->Fill(jl);
      if( hmap1 & (1<<jl) ) hDauHitPat->Fill(jl);
      if( hmap2 & (1<<jl) ) hDauHitPat->Fill(jl);
    }

    Double_t d0x1 = recProbe[0].GetX();
    Double_t d0y1 = recProbe[0].GetY();
    Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
    if (d0x1 < 0)
      d0xy1 *= -1;
    
    Double_t d0x2 = recProbe[1].GetX();
    Double_t d0y2 = recProbe[1].GetY();
    Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
    if (d0x2 < 0)
      d0xy2 *= -1;

    Double_t d0x3 = recProbe[2].GetX();
    Double_t d0y3 = recProbe[2].GetY();
    Double_t d0xy3 = TMath::Sqrt(d0x3 * d0x3 + d0y3 * d0y3);
    if (d0x3 < 0)
      d0xy3 *= -1;

    if(TMath::Abs(d0xy1)<mind0xyDau || TMath::Abs(d0xy2)<mind0xyDau || TMath::Abs(d0xy3)<mind0xyDau) continue;

    Double_t xV, yV, zV;
    Double_t sigmaVert;
    ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xV,yV,zV,sigmaVert);
    Double_t residVx=10000.*(xV - secvertgenK[0]);
    Double_t residVy=10000.*(yV - secvertgenK[1]);
    Double_t residVz=10000.*(zV - secvertgenK[2]);

   // Get daughter track momentum at decay vertex
    for(int idau=0; idau<3; idau++){
      recProbe[idau].PropagateToZBxByBz(zV);
      recProbe[idau].GetPXYZ(pxyz);
      daurec[idau].SetXYZM(pxyz[0], pxyz[1], pxyz[2], daumass[idau]);
    }

    parent = daurec[0];
    parent += daurec[1];
    parent += daurec[2];
    parentgen = daugen[0];
    parentgen += daugen[1];
    parentgen += daugen[2];

    Double_t  ptRecD=parent.Pt();
    Double_t  massRecD=parent.M();
    Double_t yRecD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
    hYPtRecoAll->Fill(yRecD, ptRecD);
    hPtRecoAll->Fill(ptRecD);
    hYRecoAll->Fill(yRecD);
    hMassAll->Fill(massRecD);
    if (nfake > 0){
      hYPtRecoFake->Fill(yRecD, ptRecD);
      hPtRecoFake->Fill(ptRecD);
      hMassFake->Fill(massRecD);
    }
    hMassVsPt->Fill(massRecD,ptRecD);
    hMassVsY->Fill(massRecD,yRecD);
    
    hResVx->Fill(residVx, ptRecD);
    hResVy->Fill(residVy, ptRecD);
    hResVz->Fill(residVz, ptRecD);
    hResVxVsY->Fill(residVx, yRecD);
    hResVyVsY->Fill(residVy, yRecD);
    hResVzVsY->Fill(residVz, yRecD);
    hSigVert->Fill(sigmaVert,ptRecD);
    
    hResPx->Fill(daurec[0].Px() - daugen[0].Px(), ptRecD);
    hResPy->Fill(daurec[0].Py() - daugen[0].Py(), ptRecD);
    hResPz->Fill(daurec[0].Pz() - daugen[0].Pz(), ptRecD);
    hResPx->Fill(daurec[1].Px() - daugen[1].Px(), ptRecD);
    hResPy->Fill(daurec[1].Py() - daugen[1].Py(), ptRecD);
    hResPz->Fill(daurec[1].Pz() - daugen[1].Pz(), ptRecD);
    hResPx->Fill(daurec[2].Px() - daugen[2].Px(), ptRecD);
    hResPy->Fill(daurec[2].Py() - daugen[2].Py(), ptRecD);
    hResPz->Fill(daurec[2].Pz() - daugen[2].Pz(), ptRecD);

    hResPxVsY->Fill(daurec[0].Px() - daugen[0].Px(), yRecD);
    hResPyVsY->Fill(daurec[0].Py() - daugen[0].Py(), yRecD);
    hResPzVsY->Fill(daurec[0].Pz() - daugen[0].Pz(), yRecD);
    hResPxVsY->Fill(daurec[1].Px() - daugen[1].Px(), yRecD);
    hResPyVsY->Fill(daurec[1].Py() - daugen[1].Py(), yRecD);
    hResPzVsY->Fill(daurec[1].Pz() - daugen[1].Pz(), yRecD);
    hResPxVsY->Fill(daurec[2].Px() - daugen[2].Px(), yRecD);
    hResPyVsY->Fill(daurec[2].Py() - daugen[2].Py(), yRecD);
    hResPzVsY->Fill(daurec[2].Pz() - daugen[2].Pz(), yRecD);
        

    Float_t dist = TMath::Sqrt(xV * xV + yV * yV + zV * zV);
    Float_t distXY = TMath::Sqrt(xV * xV + yV * yV);
    Float_t distZ = zV;
    Float_t distgen = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1] + secvertgenPi[2] * secvertgenPi[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgenPi[0] * secvertgenPi[0] + secvertgenPi[1] * secvertgenPi[1]);
    
    Double_t vsec[3] = {xV, yV, zV};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t cospxy = CosPointingAngleXY(vprim, vsec, parent);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRecD);
    hCospXY->Fill(cospxy, ptRecD);
    
    hResDist->Fill(dist - distgen, ptRecD);
    hResDistXY->Fill(distXY - distgenXY, ptRecD);
    
    hDistXY->Fill(distXY, ptRecD);
    hDistZ->Fill(vsec[2], ptRecD);
    hDist->Fill(dist, ptRecD);
    hDistgenXY->Fill(distgenXY, ptRecD);
    hDistgen->Fill(distgen, ptRecD);
      
    //AliExternalTrackParam *track1 = (AliExternalTrackParam *)recProbe[0].GetTrack();
    //AliExternalTrackParam *track2 = (AliExternalTrackParam *)recProbe[1].GetTrack();
    
    hd0XY1->Fill(d0xy1, ptRecD);
    hd0XY2->Fill(d0xy2, ptRecD);
    hd0XY3->Fill(d0xy3, ptRecD);
      
    if(cosp>cutCosPointCand && d0xy1*d0xy3<cutImpParProd){
      arrsp[0] = massRecD;
      arrsp[1] = ptRecD;
      arrsp[2] = yRecD;
      arrsp[3] = dist;
      arrsp[4] = distXY;
      arrsp[5] = cosp;
      arrsp[6] = TMath::Min(TMath::Abs(d0xy1),TMath::Min(TMath::Abs(d0xy2),TMath::Abs(d0xy3)));
      arrsp[7] = d0xy1*d0xy3; // same sign daughters
      arrsp[8] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(
      arrsp[9] = momPi/momP;
      arrsp[10] = TMath::Abs(ipD);	    
      hsp->Fill(arrsp);
      
      if (ntLccand){
	arrnt[0] = massRecD;
	arrnt[1] = ptRecD;
	arrnt[2] = yRecD;
	arrnt[3] = dist;
	arrnt[4] = cosp;
	arrnt[5] = d0xy1;
	arrnt[6] = d0xy2;
	arrnt[7] = d0xy3;
	arrnt[8] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
	arrnt[9] = momP;
	arrnt[10] = momK;
	arrnt[11] = momPi;
	arrnt[12] = TMath::Abs(ipD);	    
	ntLccand->Fill(arrnt);
      }
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
  hMassVsPt->Write();
  hMassVsY->Write();
  hYPtGen->Write();
  hPtGen->Write();
  hYGen->Write();
  hYPtRecoAll->Write();
  hYPtRecoFake->Write();
  hPtRecoAll->Write();
  hYRecoAll->Write();
  hPtRecoFake->Write();
  hDauClu->Write();
  hDauHitPat->Write();
  hDistXY->Write();
  hDistZ->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hCospXY->Write();
  hSigVert->Write();
  hResVx->Write();
  hResVy->Write();
  hResVz->Write();
  hResPx->Write();
  hResPy->Write();
  hResPz->Write();
  hResVxVsY->Write();
  hResVyVsY->Write();
  hResVzVsY->Write();
  hResPxVsY->Write();
  hResPyVsY->Write();
  hResPzVsY->Write();
  hResDist->Write();
  hResDistXY->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hd0XY3->Write();
  hNevents->Write();
  hsp->Write();
  if (ntLccand){
    fnt->cd();
    ntLccand->Write();
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
  hptP->Write();
  hyK->Write();
  hyPi->Write();
  hyP->Write();
  hYGen->Write();
  hyPiK->Write();
  fout2.Close();

  fout->Close();
}



void MakeLambdacCombinBkgCandidates(const char *setup = "setup-10um-itssa_Eff1.txt",
				    const char* trackTreeFile="treeBkgEvents.root",
				    const char *selectionFile="",
				    Int_t nevents = 999999, 
				    Int_t writeNtuple = kFALSE,
				    Bool_t usePID=kFALSE,
				    Double_t maxMomForPID=1.9){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create Lc combinatorial background candidates (= triplets of tracks)
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

  ConfigureSelectionsAndAxes(selectionFile);

  TFile *filetree = new TFile(trackTreeFile);
  TTree *tree = (TTree *)filetree->Get("tree");
  TClonesArray *arr = 0;
  tree->SetBranchAddress("tracks", &arr);
  Int_t entries = tree->GetEntries();
  printf("Number of events in tree = %d\n",entries);
  if(nevents>entries) nevents=entries;
  else printf(" --> Analysis performed on first %d events\n",nevents);
  if(usePID) printf("Rough PID cuts will be used\n");
  
  TDatime dt;
  static UInt_t seed = dt.Get();
  gRandom->SetSeed(seed);

  TFile *fout = new TFile("Lc-Bkg-histos.root", "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 40, 1., 5., 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 1., 3.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistZ = new TH2F("hDistZ", "", 100, 0, 0.2, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hSigVert = new TH2F("hSigVert", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY3 = new TH2F("hd0xy3", "", 100, -0.1, 0.1, 30, 0, 3);

  TH1D* hMomPion = new TH1D("hMomPion","",200,0.,10.);
  TH1D* hMomKaon = new TH1D("hMomKaon","",200,0.,10.);
  TH1D* hMomProton = new TH1D("hMomProton","",200,0.,10.);
  
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
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
  
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntLccand = 0x0;
  Float_t arrnt[13];
  if (writeNtuple){
    fnt = new TFile("Lc-Bkg-ntuple.root", "recreate");
    ntLccand = new TNtuple("ntLccand", "ntLccand", "mass:pt:y:dist:cosp:d01:d02:d03:sigvert:pP:pK:pPi:d0Lc", 32000);
  }

  // define mother particle
  Int_t pdgParticle = 4122;
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdgParticle)->Mass();

  KMCProbeFwd recProbe[3];
  TLorentzVector parent, daurec[3];

  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    Double_t countCandInPeak = 0;
    Double_t countCand = 0;
    tree->GetEvent(iev);
    Int_t arrentr = arr->GetEntriesFast();
    
    Double_t pxyz0[3],pxyz1[3],pxyz2[3];
    for (Int_t itr = 0; itr < arrentr; itr++){
      KMCProbeFwd *tr1 = (KMCProbeFwd *)arr->At(itr);
      // cout << "tr P=" << tr1->GetP() << endl;
      if(tr1->GetNHits()<minITShits) continue;
      Float_t ch1 = tr1->GetCharge();
      recProbe[0] = *tr1;
      recProbe[0].PropagateToZBxByBz(vprim[2]);
      recProbe[0].GetPXYZ(pxyz0);
      Double_t d0x1 = recProbe[0].GetX();
      Double_t d0y1 = recProbe[0].GetY();
      Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
      if (d0x1 < 0) d0xy1 *= -1;
      if(TMath::Abs(d0xy1)<mind0xyDau) continue;

      for (Int_t itr2 = 0; itr2 < arrentr; itr2++){
	if(itr2==itr) continue;
	KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
	if(tr2->GetNHits()<minITShits) continue;
	Float_t ch2 = tr2->GetCharge();
	// convention: charge signs are ordered as +-+ or -+-
	if (ch1 * ch2 > 0) continue;
	recProbe[1] = *tr2;
	recProbe[1].PropagateToZBxByBz(vprim[2]);
	recProbe[1].GetPXYZ(pxyz1);
	Double_t d0x2 = recProbe[1].GetX();
	Double_t d0y2 = recProbe[1].GetY();
	Double_t d0xy2 = TMath::Sqrt(d0x2 * d0x2 + d0y2 * d0y2);
	if (d0x2 < 0) d0xy2 *= -1;
	if(TMath::Abs(d0xy2)<mind0xyDau) continue;
	// recProbe[0].PropagateToDCA(&recProbe[1]);
	// Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
	// Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
	// Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
	// Float_t dca01 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
	// recProbe[0].PropagateToZBxByBz(vprim[2]);
	// recProbe[1].PropagateToZBxByBz(vprim[2]);

	for (Int_t itr3 = itr2+1; itr3 < arrentr; itr3++){
	  if(itr3==itr) continue;
	  KMCProbeFwd *tr3 = (KMCProbeFwd *)arr->At(itr3);
	  if(tr3->GetNHits()<minITShits) continue;
	  Float_t ch3 = tr3->GetCharge();
	  // convention: charge signs are ordered as +-+ or -+-
	  if (ch3 * ch2 > 0) continue;
	  recProbe[2] = *tr3;
	  recProbe[2].PropagateToZBxByBz(vprim[2]);
	  recProbe[2].GetPXYZ(pxyz2);
	  Double_t d0x3 = recProbe[2].GetX();
	  Double_t d0y3 = recProbe[2].GetY();
	  Double_t d0xy3 = TMath::Sqrt(d0x3 * d0x3 + d0y3 * d0y3);
	  if (d0x3 < 0) d0xy3 *= -1;
	  //printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
	  if(TMath::Abs(d0xy3)<mind0xyDau) continue;

	  for(Int_t iMassHyp=0; iMassHyp<2; iMassHyp++){
	    // mass hypothesis: pKpi, piKp
	    Double_t momPi=0;
	    Double_t momK=0;
	    Double_t momP=0;
	    if(iMassHyp==0){
	      daurec[0].SetXYZM(pxyz0[0], pxyz0[1], pxyz0[2], KMCDetectorFwd::kMassP);
	      daurec[1].SetXYZM(pxyz1[0], pxyz1[1], pxyz1[2], KMCDetectorFwd::kMassK);
	      daurec[2].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassPi);
	      momPi=recProbe[2].GetTrack()->P();
	      momK=recProbe[1].GetTrack()->P();
	      momP=recProbe[0].GetTrack()->P();
	      if(usePID && recProbe[0].GetMass()<0.2 && momP<maxMomForPID) continue; // rough PID (reject pions with p<thresh)
	      if(usePID && recProbe[2].GetMass()>0.7 && momPi<maxMomForPID) continue; // rough PID (reject protons with p<thresh)
	    }else{
	      daurec[0].SetXYZM(pxyz0[0], pxyz0[1], pxyz0[2], KMCDetectorFwd::kMassPi);
	      daurec[1].SetXYZM(pxyz1[0], pxyz1[1], pxyz1[2], KMCDetectorFwd::kMassK);
	      daurec[2].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassP);
	      momPi=recProbe[0].GetTrack()->P();
	      momK=recProbe[1].GetTrack()->P();
	      momP=recProbe[2].GetTrack()->P();	      
	      if(usePID && recProbe[2].GetMass()<0.2 && momP<maxMomForPID) continue; // rough PID (reject pions with p<thresh)
	      if(usePID && recProbe[0].GetMass()>0.7 && momPi<maxMomForPID) continue; // rough PID (reject protons with p<thresh)
	    }
	    parent = daurec[0];
	    parent += daurec[1];
	    parent += daurec[2];
	    countCand++;
	    Float_t ptD=parent.Pt();
	    Float_t invMassD=parent.M();
	    Float_t yD = 0.5 * TMath::Log((parent.E() + parent.Pz()) / (parent.E() - parent.Pz()));
	    hYPtRecoAll->Fill(yD, ptD);
	    hPtRecoAll->Fill(ptD);
	    hMassAll->Fill(invMassD);
	    if(invMassD>2.15  && invMassD<2.45){
	      // range to fill histos
	      if(TMath::Abs(invMassD-mass)<0.06) countCandInPeak++;
	      hMomPion->Fill(momPi);
	      hMomKaon->Fill(momK);
	      hMomProton->Fill(momP);

	      Double_t xV, yV, zV;
	      Double_t sigmaVert;
	      ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xV,yV,zV,sigmaVert);
	      hVx->Fill(xV);
	      hVy->Fill(yV);
	      hVz->Fill(zV);
	      hSigVert->Fill(sigmaVert,ptD);

	      Float_t dist = TMath::Sqrt(xV * xV + yV * yV + zV * zV);
	      Float_t distXY = TMath::Sqrt(xV * xV + yV * yV);
	      Float_t distZ = zV;
	      Double_t vsec[3] = {xV, yV, zV};
	      Double_t cosp = CosPointingAngle(vprim, vsec, parent);
	      Double_t cospxy = CosPointingAngleXY(vprim, vsec, parent);
	      Double_t ipD = ImpParXY(vprim, vsec, parent);
	      hCosp->Fill(cosp, ptD);
	      hCospXY->Fill(cospxy, ptD);
	      //printf(" ***** ***** cos point = %f \n", cosp);	    
	      hDistXY->Fill(distXY, ptD);
	      hDistZ->Fill(zV, ptD);
	      hDist->Fill(dist, ptD);
	      hd0XY1->Fill(d0xy1, ptD);
	      hd0XY2->Fill(d0xy2, ptD);
	      hd0XY3->Fill(d0xy3, ptD);
	      if(cosp>cutCosPointCand && d0xy1*d0xy3<cutImpParProd){
		arrsp[0] = invMassD;
		arrsp[1] = ptD;
		arrsp[2] = yD;
		arrsp[3] = dist;
		arrsp[4] = distXY;
		arrsp[5] = cosp;
		arrsp[6] = TMath::Min(TMath::Abs(d0xy1),TMath::Min(TMath::Abs(d0xy2),TMath::Abs(d0xy3)));
		arrsp[7] = d0xy1*d0xy3; // same sign daughters
		arrsp[8] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(
		arrsp[9] = momPi/momP;
		arrsp[10] = TMath::Abs(ipD);	    
		hsp->Fill(arrsp);
		
		if (ntLccand){
		  arrnt[0] = invMassD;
		  arrnt[1] = ptD;
		  arrnt[2] = yD;
		  arrnt[3] = dist;
		  arrnt[4] = cosp;
		  arrnt[5] = d0xy1;
		  arrnt[6] = d0xy2;
		  arrnt[7] = d0xy3;
		  arrnt[8] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
		  arrnt[9] = momP;
		  arrnt[10] = momK;
		  arrnt[11] = momPi;
		  arrnt[12] = TMath::Abs(ipD);	    
		  ntLccand->Fill(arrnt);
		}
	      }
	    } // check on inv mass
	  } // loop on mass hypothesis
	} // loop on third track
      } // loop on second track
    }// loop on first track
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot Lc candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
  }
  
  fout->cd();
  hNevents->Write();
  hcand->Write();
  hcandpeak->Write();  
  hMassAll->Write();
  hYPtRecoAll->Write();
  hPtRecoAll->Write();
  hDistXY->Write();
  hDistZ->Write();
  hDist->Write();
  hSigVert->Write();
  hVx->Write();
  hVy->Write();
  hVz->Write();
  hCosp->Write();
  hCospXY->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hd0XY3->Write();
  hMomPion->Write();
  hMomKaon->Write();
  hMomProton->Write();
  hsp->Write();
  if (ntLccand){
    fnt->cd();
    hNevents->Write();
    hcand->Write();
    hcandpeak->Write();  
    ntLccand->Write("",TObject::kOverwrite);
    fnt->Close();
  }
  fout->Close();
}

THnSparseF* CreateSparse(){
  const Int_t nAxes=11;
  TString axTit[nAxes]={"Inv. mass (GeV/c^{2})",
			"p_{T} (GeV/c)",
			"y",
			"Dec Len (cm)",
			"DecLenXY (cm)",
			"cos(#vartheta_{p})",
			"d_0^{min} (cm)",
			"d01xd03 (cm2)",
			"sigmaVert",
			"p_{#pi}/p_{p}",
			"d_0^{Lc} (cm)"};
  Int_t bins[nAxes] =   {60,  3,  16, nBinsDecLen, nBinsDecLenXY, nBinsCosPoi, nBinsImpDau, nBinsImpPro, nBinsSigVer, nBinsPrat, nBinsImpLc};
  Double_t min[nAxes] = {2.15,  0., 1., minDecLen,   minDecLenXY,   minCosPoi,   minImpDau,   minImpPro,   minSigVer,   minPrat,   minImpLc};
  Double_t max[nAxes] = {2.45,  3., 5., maxDecLen,   maxDecLenXY,   maxCosPoi,   maxImpDau,   maxImpPro,   maxSigVer,   maxPrat,   maxImpLc};
  THnSparseF *hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
  for(Int_t iax=0; iax<nAxes; iax++) hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  return hsp;
}


void ConfigureSelectionsAndAxes(const char *selectionFile){
  if(strlen(selectionFile)==0) return;
  if(gSystem->Exec(Form("ls -l %s",selectionFile)) !=0 ){
    printf("File %s with configuration of selections not found\n",selectionFile);
    return;
  }
  FILE* confFil=fopen(selectionFile,"r");
  char line[50];
  int n;
  float x;
  bool readok;
  while(!feof(confFil)){
    readok=fscanf(confFil,"%s:",line);
    if(strstr(line,"MinVThits")){
      readok=fscanf(confFil,"%d",&n);
      minITShits=n;
    }
    if(strstr(line,"MinImpParXY")){
      readok=fscanf(confFil,"%f",&x);
      mind0xyDau=x;
    }
    if(strstr(line,"CosPointCut")){
      readok=fscanf(confFil,"%f",&x);
      cutCosPointCand=x;
    }
    if(strstr(line,"ImpParProdCut")){
      readok=fscanf(confFil,"%f",&x);
      cutImpParProd=x;
    }
    else if(strstr(line,"NumOfDecLenBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsDecLen=n;
    }
    else if(strstr(line,"MinDecLenForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minDecLen=x;
    }
    else if(strstr(line,"MaxDecLenForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxDecLen=x;
    }
    else if(strstr(line,"NumOfDecLenXYBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsDecLenXY=n;
    }
    else if(strstr(line,"MinDecLenXYForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minDecLenXY=x;
    }
    else if(strstr(line,"MaxDecLenXYForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxDecLenXY=x;
    }
    else if(strstr(line,"NumOfCosPoiBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsCosPoi=n;
    }
    else if(strstr(line,"MinCosPoiForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minCosPoi=x;
    }
    else if(strstr(line,"MaxCosPoiForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxCosPoi=x;
    }
    else if(strstr(line,"NumOfImpParDauBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsImpDau=n;
    }
    else if(strstr(line,"MinImpParDauForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minImpDau=x;
    }
    else if(strstr(line,"MaxImpParDauForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxImpDau=x;
    }
    else if(strstr(line,"NumOfImpParProdBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsImpPro=n;
    }
    else if(strstr(line,"MinImpParProdForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minImpPro=x;
    }
    else if(strstr(line,"MaxImpParProdForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxImpPro=x;
    }
    else if(strstr(line,"NumOfSigVerBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsSigVer=n;
    }
    else if(strstr(line,"MinSigVerForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minSigVer=x;
    }
    else if(strstr(line,"MaxSigVerForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxSigVer=x;
    }
    else if(strstr(line,"NumOfPratioBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsPrat=n;
    }
    else if(strstr(line,"MinPratioForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minPrat=x;
    }
    else if(strstr(line,"MaxPratioForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxPrat=x;
    }
    else if(strstr(line,"NumOfImpParLcBins")){
      readok=fscanf(confFil,"%d",&n);
      nBinsImpLc=n;
    }
    else if(strstr(line,"MinImpParLcForSparse")){
      readok=fscanf(confFil,"%f",&x);
      minImpLc=x;
    }
    else if(strstr(line,"MaxImpParLcForSparse")){
      readok=fscanf(confFil,"%f",&x);
      maxImpLc=x;
    }
  }
}
