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

// Track Chi2Tot cut
double ChiTot = 1.5;

// settings for signal generation
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double ptminSG = 0.;
double ptmaxSG = 10; //Elena's change, it was 3 GeV/c

double vX = 0, vY = 0, vZ = 0; // event vertex

THnSparseF* CreateSparse();
TDatime dt;

void GenerateDplusSignalCandidates(Int_t nevents = 100000, 
				   double Eint = 160., 
				   const char *setup = "setup-10um-itssa_Eff1.txt", 
				   const char *filNamPow="/home/prino/cernbox/na60plus/POWHEG/pp20/Charm1dot5/pp0_frag-PtSpectra-Boost.root", 
				   const char *privateDecayTable = "../decaytables/USERTABLC.DEC",
				   int optPartAntiPart=3,
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
  printf("--> pt and y shape of D+ from %s\n",filNamPow);
  TFile *filPow=new TFile(filNamPow);
  TH3D* h3Dpow=(TH3D*)filPow->Get("hptyeta411");
  TH1D *hDpluspt = (TH1D*)h3Dpow->ProjectionX("hDpluspt");
  TH1D *hDplusy = (TH1D*)h3Dpow->ProjectionY("hDplusy");
  TH3D* h3Dbarpow=(TH3D*)filPow->Get("hptyetam411");
  if(h3Dbarpow){
    TH1D *hDminuspt = (TH1D*)h3Dbarpow->ProjectionX("hDminuspt");
    TH1D *hDminusy = (TH1D*)h3Dbarpow->ProjectionY("hDminusy");
    if(optPartAntiPart==3){
      hDpluspt->Add(hDminuspt);
      hDplusy->Add(hDminusy);
    }else if(optPartAntiPart==2){
      hDpluspt=hDminuspt;
      hDplusy=hDminusy;
    }
  }

  TH2F *hptK = new TH2F("hptK", "kaons from D+ decays", 50,0.,10.,50, 0., 10.);
  TH2F *hptPi = new TH2F("hptPi", "pions from D+ decays", 50, 0.,10.,50,0., 10.);
  TH1D *hyK = new TH1D("hyK", "y kaons from D+ decays", 50, 0., 5.);
  TH1D *hyPi = new TH1D("hyPi", "y pions from D+ decays", 50, 0., 5.);
  
 
  TFile *fout = new TFile("Dplus-Signal-histos.root", "recreate");
   
  //Magnetic field and detector parameters

  //int outN = nev/10;
  //if (outN<1) outN=1;

  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
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
  KMCProbeFwd recProbe[3];  
  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); //read the default decay table DECAY.DEC and particle table
  bool privTab=kFALSE;
  if (strlen(privateDecayTable)>0){
    if(gSystem->Exec(Form("ls -l %s",privateDecayTable))==0){
      fDecayer->SetDecayTablePath((char*)privateDecayTable);
      fDecayer->ReadDecayTable();
      printf("-- Use D+ decay table from file %s\n",privateDecayTable);
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
  Int_t pdgParticle = 411;
  
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

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistZ = new TH2F("hDistZ", "", 100, 0, 0.2, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 300, -1, 1, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
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
  TNtuple *ntDpcand = 0x0;
  if (writeNtuple){
    fnt = new TFile("fntSig.root", "recreate");
    ntDpcand = new TNtuple("ntDpcand", "ntDpcand", "mass:pt:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
  }
  Float_t arrnt[10];
  for (Int_t iev = 0; iev < nevents; iev++){
    hNevents->Fill(0.5);
    Double_t vprim[3] = {0, 0, 0};
    if(iev%100==0) printf(" ***************  ev = %d \n", iev);
    int nrec = 0;
    int nfake = 0;
    double pxyz[3];
    int nPions=0;
    int nKaons=0;
    
    if (simulateBg && (iev%refreshBg)==0) det->GenBgEvent(0.,0.,0.);
    Double_t ptGenD = hDpluspt->GetRandom(); // get D0 distribution from file
    Double_t yGenD = hDplusy->GetRandom();
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
    Double_t ptPi1 = -999.;
    Double_t ptPi2 = -999.;
    Double_t yK=-999.;
    Double_t yPi1 = -999.;
    Double_t yPi2 = -999.;
    Int_t icount = 0;
    Double_t secvertgenK[3]={0.,0.,0.};
    Double_t secvertgenPi1[3]={0.,0.,0.};
    Double_t secvertgenPi2[3]={0.,0.,0.};
    
    // loop on decay products
    for (int i = 0; i < np; i++) { 
      TParticle *iparticle1 = (TParticle *)particles->At(i);
      Int_t kf = TMath::Abs(iparticle1->GetPdgCode());
      vX = iparticle1->Vx();
      vY = iparticle1->Vy();
      vZ = iparticle1->Vz();
      if (kf == pdgParticle){
	// D+ particle
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
	if (trw->GetNormChi2(kTRUE) > ChiTot) continue;
	nrec++;
	    
	nfake += trw->GetNFakeITSHits();
	trw->GetPXYZ(pxyz);
	if (kf == 321){
	  // Kaon daughter
	  ptK = iparticle1->Pt();
	  yK = iparticle1->Y();
	  hptK->Fill(ptGenD,ptK);
	  hyK->Fill(iparticle1->Y());
	  secvertgenK[0] = iparticle1->Vx();
	  secvertgenK[1] = iparticle1->Vy();
	  secvertgenK[2] = iparticle1->Vz();
	  daugen[1].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	  daurec[1].SetXYZM(pxyz[0], pxyz[1], pxyz[2], iparticle1->GetMass());
	  recProbe[1] = *trw;
	  ++nKaons;
	}else if (kf == 211){
	  // Pion daughter
	  if(nPions==0){
	    ptPi1 = iparticle1->Pt();
	    yPi1 = iparticle1->Y();
	    hptPi->Fill(ptGenD,ptPi1);
	    hyPi->Fill(yPi1);
	    secvertgenPi1[0] = iparticle1->Vx();
	    secvertgenPi1[1] = iparticle1->Vy();
	    secvertgenPi1[2] = iparticle1->Vz();
	    daugen[0].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	    daurec[0].SetXYZM(pxyz[0], pxyz[1], pxyz[2],  iparticle1->GetMass());
	    recProbe[0] = *trw;
	  }else if(nPions==1){
	    ptPi2 = iparticle1->Pt();
	    yPi2 = iparticle1->Y();
	    hptPi->Fill(ptGenD,ptPi2);
	    hyPi->Fill(yPi2);
	    secvertgenPi2[0] = iparticle1->Vx();
	    secvertgenPi2[1] = iparticle1->Vy();
	    secvertgenPi2[2] = iparticle1->Vz();
	    daugen[2].SetXYZM(iparticle1->Px(), iparticle1->Py(), iparticle1->Pz(), iparticle1->GetMass());
	    daurec[2].SetXYZM(pxyz[0], pxyz[1], pxyz[2],  iparticle1->GetMass());
	    recProbe[2] = *trw;
	  }
	  ++nPions;
	}
      }
    }
    if (nrec < 3 || nPions!=2 || nKaons!=1) continue;
    
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
    
    // recProbe[0].PropagateToDCA(&recProbe[1]);
    // Float_t d1 = recProbe[1].GetX() - recProbe[0].GetX();
    // Float_t d2 = recProbe[1].GetY() - recProbe[0].GetY();
    // Float_t d3 = recProbe[1].GetZ() - recProbe[0].GetZ();
    // Float_t dca01 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // hDCA->Fill(dca01, ptRecD);
    // recProbe[0].PropagateToDCA(&recProbe[2]);
    // d1 = recProbe[2].GetX() - recProbe[0].GetX();
    // d2 = recProbe[2].GetY() - recProbe[0].GetY();
    // d3 = recProbe[2].GetZ() - recProbe[0].GetZ();
    // Float_t dca02 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // hDCA->Fill(dca02, ptRecD);
    // recProbe[1].PropagateToDCA(&recProbe[2]);
    // d1 = recProbe[2].GetX() - recProbe[1].GetX();
    // d2 = recProbe[2].GetY() - recProbe[1].GetY();
    // d3 = recProbe[2].GetZ() - recProbe[1].GetZ();
    // Float_t dca12 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    // hDCA->Fill(dca12, ptRecD);

    Double_t xV, yV, zV;
    Double_t sigmaVert;
    ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xV,yV,zV,sigmaVert);
    Double_t residVx=10000.*(xV - secvertgenK[0]);
    Double_t residVy=10000.*(yV - secvertgenK[1]);
    Double_t residVz=10000.*(zV - secvertgenK[2]);
    hResVx->Fill(residVx, ptRecD);
    hResVy->Fill(residVy, ptRecD);
    hResVz->Fill(residVz, ptRecD);
    hResVxVsY->Fill(residVx, yRecD);
    hResVyVsY->Fill(residVy, yRecD);
    hResVzVsY->Fill(residVz, yRecD);
    
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
    Float_t distgen = TMath::Sqrt(secvertgenPi1[0] * secvertgenPi1[0] + secvertgenPi1[1] * secvertgenPi1[1] + secvertgenPi1[2] * secvertgenPi1[2]);
    Float_t distgenXY = TMath::Sqrt(secvertgenPi1[0] * secvertgenPi1[0] + secvertgenPi1[1] * secvertgenPi1[1]);
    
    Double_t vsec[3] = {xV, yV, zV};
    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
    Double_t ipD = ImpParXY(vprim, vsec, parent);
    hCosp->Fill(cosp, ptRecD);
    // printf(" ***** ***** cos point = %f \n", cosp);
    //if (cosp < -0.98)
    //    printf("SMALL COSPOINT");
    
    hResDist->Fill(dist - distgen, ptRecD);
    hResDistXY->Fill(distXY - distgenXY, ptRecD);
    
    //recProbe[0].PropagateToDCA(&recProbe[1]);
    
    hDistXY->Fill(distXY, ptRecD);
    hDistZ->Fill(vsec[2], ptRecD);
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

    recProbe[2].PropagateToZBxByBz(0);
    Double_t d0x3 = recProbe[2].GetX();
    Double_t d0y3 = recProbe[2].GetY();
    Double_t d0xy3 = TMath::Sqrt(d0x3 * d0x3 + d0y3 * d0y3);
    if (d0x3 < 0)
      d0xy3 *= -1;
    
    // printf("d0xy1 = %f, d0xy2 = %f \n", d0xy1, d0xy2);
    
    hd0XY1->Fill(d0xy1, ptRecD);
    hd0XY2->Fill(d0xy2, ptRecD);
    hd0XY3->Fill(d0xy3, ptRecD);
      
    arrsp[0] = massRecD;
    arrsp[1] = ptRecD;
    arrsp[2] = dist;
    arrsp[3] = distXY;
    arrsp[4] = distZ;
    arrsp[5] = cosp;
    arrsp[6] = TMath::Min(TMath::Abs(d0xy1),TMath::Min(TMath::Abs(d0xy2),TMath::Abs(d0xy3)));
    arrsp[7] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
    arrsp[9] = TMath::Abs(ipD);	    
    hsp->Fill(arrsp);
    
    if (ntDpcand){
      arrnt[0] = massRecD;
      arrnt[1] = ptRecD;
      arrnt[2] = dist;
      arrnt[3] = cosp;
      arrnt[4] = d0xy1;
      arrnt[5] = d0xy2;
      arrnt[6] = d0xy1 * d0xy2;
      arrnt[7] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
      arrnt[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      arrnt[9] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
      ntDpcand->Fill(arrnt);
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
  hDistXY->Write();
  hDistZ->Write();
  hDist->Write();
  hDistgenXY->Write();
  hDistgen->Write();
  hCosp->Write();
  hDCA->Write();
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
  hNevents->Write();
  hsp->Write();
  if (ntDpcand){
    fnt->cd();
    ntDpcand->Write();
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
  hyK->Write();
  hyPi->Write();
  hYGen->Write();
  fout2.Close();

  fout->Close();
}



void MakeDplusCombinBkgCandidates(const char* trackTreeFile="treeBkgEvents.root",
				  Int_t nevents = 999999, 
				  Int_t writeNtuple = kFALSE,
				  Bool_t usePID=kFALSE){

  // Read the TTree of tracks produced with runBkgVT.C
  // Create D0 combinatorial background candidates (= OS pairs of tracks)
  // Store in THnSparse and (optionally) TNtuple

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

  TFile *fout = new TFile("Dplus-Bkg-histos.root", "recreate");
  TH1D* hPtRecoAll = new TH1D("hPtRecoAll", "Pt all match", 50, 0., 5.);
  TH2F* hYPtRecoAll = new TH2F("hYPtRecoAll", "Y-Pt all match", 80, 1.0, 5.4, 50, 0., 5.);
  TH1D* hMassAll = new TH1D("hMassAll", "Mass all match", 250, 1., 3.5);

  TH2F *hDistXY = new TH2F("hDistXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistZ = new TH2F("hDistZ", "", 100, 0, 0.2, 30, 0, 3);
  TH2F *hDist = new TH2F("hDist", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDistgenXY = new TH2F("hDistgenXY", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hDistgen = new TH2F("hDistgen", "", 300, 0, 10, 30, 0, 3);
  TH2F *hDCA = new TH2F("hDCA", "", 100, 0, 0.1, 30, 0, 3);
  TH2F *hd0XY1 = new TH2F("hd0xy1", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY2 = new TH2F("hd0xy2", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hd0XY3 = new TH2F("hd0xy3", "", 100, -0.1, 0.1, 30, 0, 3);

  TH1D* hMomPion = new TH1D("hMomPion","",200,0.,10.);
  TH1D* hMomKaon = new TH1D("hMomKaon","",200,0.,10.);

  TH2F *hResVx = new TH2F("hResVx", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVy = new TH2F("hResVy", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResVz = new TH2F("hResVz", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hResPx = new TH2F("hResPx", "", 100, -1, 1, 30, 0, 3); //for Kaons
  TH2F *hResPy = new TH2F("hResPy", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResPz = new TH2F("hResPz", "", 100, -1, 1, 30, 0, 3);
  TH2F *hResDist = new TH2F("hResDist", "", 100, -0.5, 0.5, 30, 0, 3);
  TH2F *hResDistXY = new TH2F("hResDistXY", "", 100, -0.1, 0.1, 30, 0, 3);
  TH2F *hCosp = new TH2F("hCosp", "", 100, -1, 1, 30, 0, 3);
  
  TH2F *hd0 = new TH2F("hd0", "", 100, 0, 0.1, 30, 0, 3);
  
  TH1D *hcand = new TH1D("hcand", "", 1000, 0, 500000000);
  TH1D *hcandpeak = new TH1D("hcandpeak", "", 500, 0, 15000000);
  TH1D *hNevents = new TH1D("hNevents", "", 1, 0, 1);
  
  THnSparseF *hsp = CreateSparse();
  const Int_t nDim=static_cast<const Int_t>(hsp->GetNdimensions());
  Double_t arrsp[nDim];

  TFile *fnt = 0x0;
  TNtuple *ntDpcand = 0x0;
  Float_t arrnt[10];
  if (writeNtuple){
    fnt = new TFile("fntBkg.root", "recreate");
    ntDpcand = new TNtuple("ntDpcand", "ntDpcand", "mass:pt:dist:cosp:d01:d02:d0prod:dca:ptMin:ptMax", 32000);
  }

  
  // define mother particle
  Int_t pdgParticle = 411;
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
      Float_t ch1 = tr1->GetCharge();
      recProbe[0] = *tr1;
      recProbe[0].PropagateToZBxByBz(vprim[2]);
      recProbe[0].GetPXYZ(pxyz0);
      Double_t d0x1 = recProbe[0].GetX();
      Double_t d0y1 = recProbe[0].GetY();
      Double_t d0xy1 = TMath::Sqrt(d0x1 * d0x1 + d0y1 * d0y1);
      if (d0x1 < 0) d0xy1 *= -1;

      for (Int_t itr2 = 0; itr2 < arrentr; itr2++){
	if(itr2==itr) continue;
	KMCProbeFwd *tr2 = (KMCProbeFwd *)arr->At(itr2);
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

	  // due to the charge sign convention (+-+ or -+-) the Kaon is the middle track
	  daurec[0].SetXYZM(pxyz0[0], pxyz0[1], pxyz0[2], KMCDetectorFwd::kMassPi);
	  daurec[1].SetXYZM(pxyz1[0], pxyz1[1], pxyz1[2], KMCDetectorFwd::kMassK);
	  daurec[2].SetXYZM(pxyz2[0], pxyz2[1], pxyz2[2], KMCDetectorFwd::kMassPi);
	  Double_t momPi1=recProbe[2].GetTrack()->P();
	  Double_t momK=recProbe[1].GetTrack()->P();
	  Double_t momPi2=recProbe[0].GetTrack()->P();
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
	  if(invMassD>1.6  && invMassD<2.1){
	    // range to fill histos
	    if(TMath::Abs(invMassD-mass)<0.06) countCandInPeak++;
	    hMomPion->Fill(momPi1);
	    hMomPion->Fill(momPi2);
	    hMomKaon->Fill(momK);

	    Double_t xV, yV, zV;
	    Double_t sigmaVert;
	    ComputeVertex(recProbe[0],recProbe[1],recProbe[2],vprim[2],xV,yV,zV,sigmaVert);
	    Float_t dist = TMath::Sqrt(xV * xV + yV * yV + zV * zV);
	    Float_t distXY = TMath::Sqrt(xV * xV + yV * yV);
	    Float_t distZ = zV;
	    Double_t vsec[3] = {xV, yV, zV};
	    Double_t cosp = CosPointingAngle(vprim, vsec, parent);
	    Double_t ipD = ImpParXY(vprim, vsec, parent);
	    hCosp->Fill(cosp, ptD);
	    //printf(" ***** ***** cos point = %f \n", cosp);	    
	    hDistXY->Fill(distXY, ptD);
	    hDistZ->Fill(zV, ptD);
	    hDist->Fill(dist, ptD);
	    hd0XY1->Fill(d0xy1, ptD);
	    hd0XY2->Fill(d0xy2, ptD);
	    hd0XY3->Fill(d0xy3, ptD);
	    arrsp[0] = invMassD;
	    arrsp[1] = ptD;
	    arrsp[2] = dist;
	    arrsp[3] = distXY;
	    arrsp[4] = distZ;
	    arrsp[5] = cosp;
	    arrsp[6] = TMath::Min(TMath::Abs(d0xy1),TMath::Min(TMath::Abs(d0xy2),TMath::Abs(d0xy3)));
	    arrsp[7] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
	    arrsp[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	    arrsp[9] = TMath::Abs(ipD);	    
	    hsp->Fill(arrsp);
	      
	    if (ntDpcand){
	      arrnt[0] = invMassD;
	      arrnt[1] = ptD;
	      arrnt[2] = dist;
	      arrnt[3] = cosp;
	      arrnt[4] = d0xy1;
	      arrnt[5] = d0xy2;
	      arrnt[6] = d0xy1 * d0xy2;
	      arrnt[7] = sigmaVert;//TMath::Max(TMath::Abs(dca01),TMath::Max(TMath::Abs(dca12),TMath::Abs(dca02)));
	      arrnt[8] = TMath::Min(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      arrnt[9] = TMath::Max(recProbe[0].GetTrack()->Pt(),recProbe[1].GetTrack()->Pt());
	      ntDpcand->Fill(arrnt);
	    }
	  } // check on inv mass
	} // loop on third track
      } // loop on second track
    }// loop on first track
    hcand->Fill(countCand);
    hcandpeak->Fill(countCandInPeak);
    printf(" --> Event %d, tot D+ candidates = %.0f  in peak = %.0f\n",iev,countCand,countCandInPeak);
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
  hDCA->Write();
  hCosp->Write();
  hd0XY1->Write();
  hd0XY2->Write();
  hd0XY3->Write();
  hMomPion->Write();
  hMomKaon->Write();
  hsp->Write();
  fout->Close();
  if (ntDpcand){
    fnt->cd();
    ntDpcand->Write();
    fnt->Close();
  }
}



THnSparseF* CreateSparse(){
  const Int_t nAxes=10;
  TString axTit[nAxes]={"Inv. mass (GeV/c^{2})","p_{T} (GeV/c)",
			"Dec Len (cm)","DecLenXY (cm)","DecLenZ (cm)",
			"cos(#vartheta_{p})",
			"d_0^{min} (cm)",
			"sigmaVert",
			"p_{T}^{min} (GeV/c)",
			"d_0^{D} (cm)"};
  Int_t bins[nAxes] =   {100,  5,  30,  10,   30,  20,   12,   12,    8,   8}; 
  Double_t min[nAxes] = {2.0,  0., 0.,  0.,   0.,  0.98, 0.,   0.0,   0.,  0.};
  Double_t max[nAxes] = {2.5,  5., 0.3, 0.05, 0.3, 1.,   0.03, 0.03,  4.,  0.04};  
  THnSparseF *hsp = new THnSparseF("hsp", "hsp", nAxes, bins, min, max);
  for(Int_t iax=0; iax<nAxes; iax++) hsp->GetAxis(iax)->SetTitle(axTit[iax].Data());
  return hsp;
}
