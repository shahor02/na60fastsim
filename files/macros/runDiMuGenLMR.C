#if !defined(__CINT__) || defined(__MAKECINT__)
//#define _NOALIROOT_
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "KMCFlukaParser.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TROOT.h"
#include "GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"
#ifdef _NOALIROOT_
#include "TLocTreeStream.h"
using TTreeSRedirector = TLocTreeSRedirector;
#else
#include "TTreeStream.h"
#endif
#endif

//TRandom *rnd;

// Track Chi2Tot cut
double ChiTot = 1.5; 


// settings for backgroung generation. Used just to create occupancy in the tracker,
// no need to this f.ph.sp...
// default values for 40 GeV - recalculated in CalcBkgPar
//double y0BG   = 1.9;   // gaussian y mean - 20 GeV 
//double y0BG   = 2.08;   // gaussian y mean - 30 GeV
double y0BG   = 2.22;   // gaussian y mean - 30 GeV
double sigyBG = 1.2;   // .. sigma
double y0BGPi  = 0.7;  
double y0BGKplus  = 0.8;  
double y0BGKminus = 0.8;  
double y0BGP   = 39.8; 
double sigyBGPi = 1.18;
double sigyBGKplus  = 0.88;
double sigyBGKminus = 0.88;
double sigyBGP = 8.07;
double yminBG = 1.5;   // min y to generate 
double ymaxBG = 4.5;   // 
double TBG    = 0.17;  // inv.slope of thermal pt distribution
double dndyBGPi = 615.;
double dndyBGK = 78.;
double dndyBGP = 150.;
double NBGPi     = 74.;
double NBGKplus  = 16.2;
double NBGKminus = 6.03;
double NBGP      = 37.5;
double Piratio = 0.91;
double TBGpi = 0.17;
double TBGK = 0.23;
double TBGP = 0.25;
double ptminBG = 0.01;
double ptmaxBG = 3;

// settings for signal generation
double MotherMass;
int generYPtPar, ProcType;
//double y0SG   = 1.9;   // gaussian y mean - 20 GeV
//double y0SG   = 2.08;   // gaussian y mean - 30 GeV
double y0SG   = 2.22;   // gaussian y mean - 40 GeV
double sigySG = 1.;   // .. sigma
double yminSG = -10.;   // min y to generate 
double ymaxSG = 10.;   // 
double TSG;            // inv.slope of thermal pt distribution
double ptminSG = 0.01;
double ptmaxSG = 3;

// Added June 2019 for J/psi generation
// parameters for parameteriazion of pt distributions for J/psi from PYTHIA6
double par_p0;
double par_n;
double par_n2;
// parameters for parameteriazion of y distributions for J/psi (Schuler)
double p1_xF, p2_xF, p3_xF, p4_xF;


TF1* dndyFunSG=0;
TF1* dndptFunSG=0;

double vX=0,vY=0,vZ=0; // event vertex
KMCDetectorFwd* det = 0;

TH2F *hYPtAll=0,*hYPtFake=0, *hYPtGen;
TH1F *hMassAll=0,*hMassFake=0, *hMassFake1=0, *hMassFake2=0, *hMassFake3=0, *hMassFake4=0, *hMassFake5=0;
TH1F *hYGen=0;
TH1F *hMassSignalZoom=0;
TH2F *hYPtAllCBMbin=0, *hYPtFakeCBMbin=0;
TH1F *hPtGen, *hPtAll, *hPtFake;
TH2F *hxySextantCH[4],*hxySextantTR[2],*hxySextantITS[5];
TH2F *hSext1VsSext0, *hSext1VsSext0Trig;
TH1F *hPhi;
TH2F *hY0Y1=0;
TH1F *hY0=0, *hY1=0;
TH1F *hnFake;

//Double_t RhoLineShapeNew(Double_t *x, Double_t *para);
Int_t GetSextant(Double_t x,Double_t y);
Int_t GetSextant1(Double_t x,Double_t y);

void SetProcessParameters(const Char_t *Process, Double_t E);
void BookHistos();
void CalcBkgPar(Double_t E);

void runDiMuGenLMR(int nev=30000,     // n events to generate
		   const Char_t *Process = "jpsi",
		   int Seed=12345,
		   double Eint=110.,   // bg particles density - 30 GeV
		   int refreshBg=10,   // generate new bg event for each refreshBg-th, 1 to refresh for every signal
		   //double dndyBG=500,   // bg particles density - 20 GeV
		   //double dndyBG=650,   // bg particles density - 20 GeV 
		   const char* setup="setup.txt", // setup to load
		   const char* flukaBGList="fluka.lst", // optional fluka background file, if empty, then use parametric background
		   const char* interactionSource="" // optional primary interaction volume in fluka files, if empty, take all
		   ){

// Process can be "etaDalitz", "eta2Body", "rho", "omega2Body", "omegaDalitz", "phi", "etaPrime", "Jpsi", "jpsisch"
//
  //  rnd = new TRandom(Seed);
  gRandom->SetSeed(Seed);

  CalcBkgPar(Eint);
  TTreeSRedirector outStream("dimuGenLMR.root"); // the output stream trees will go here
  TString flukaBG = flukaBGList;

  MagField *mag = new MagField(1);
  int BNreg = mag->GetNReg();
  const double *BzMin=mag->GetZMin();
  const double *BzMax=mag->GetZMax();
  const double *BVal;
  printf("*************************************\n");
  printf("number of magnetic field regions = %d\n",BNreg);
  for(int i = 0; i < BNreg; i++){
    BVal = mag->GetBVals(i);
    printf("*** Field region %d ***\n",i);
    if(i == 0) {
      printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n",BVal[0],BVal[1],BVal[2],BzMin[i],BzMax[i]);
    } else if(i == 1) {
      printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n",BVal[0],BVal[1],BVal[2],BzMin[i],BzMax[i]);

    }
  }

  int outN = nev/10;
  if (outN<1) outN=1;
  //
  //
  det = new KMCDetectorFwd();
  //det->SetUseRPhiErrorMS(true);
  det->ReadSetup(setup,setup);
  KMCFlukaParser* flukaParser = 0;
  //det->SetUseRPhiErrorMS(true);
  if (flukaBG.IsNull()) {
    NBGPi=0; NBGKplus=0;NBGKminus=0;NBGP=0;
    //det->InitBgGeneration(dndyBG,y0BG,sigyBG,yminBG,ymaxBG,TBG,ptminBG,ptmaxBG);
    //det->InitBgGenerationPart(NBGPi,NBGKplus,NBGKminus,NBGP,Piratio,y0BG,y0BGPi,y0BGKplus,y0BGKminus,y0BGP,sigyBGPi,sigyBGKplus,sigyBGKminus,sigyBGP,yminBG,ymaxBG,TBGpi,TBGK,TBGP,ptminBG,ptmaxBG);
    det->InitBgGenerationPart(0.,dndyBGK,0.,y0BG,sigyBG,yminBG,ymaxBG,TBGpi,TBGK,TBGP,ptminBG,ptmaxBG);

    printf("pion   multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChPi());
    printf("kaon   multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChK());
    printf("proton multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChP());
    //
  } else {
    printf("Will use Fluka background from %s\n", flukaBG.Data());
    flukaParser = new KMCFlukaParser();
    flukaParser->SetInpList(flukaBG.Data());
  }
  
  // set the min N tracker hits needed to validate the track
  printf("min number of hits in MS = %d\n",det->GetNumberOfActiveLayersMS());
  det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); //NA60+
  //det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(det->GetNumberOfActiveLayersMS()); //NA60+
  //det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(det->GetNumberOfActiveLayersTR());
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
  //
  det->SetApplyBransonPCorrection(); // Branson correction

  // create signal generation f-ns
  //TF1* fRhoLineShape = new TF1("fRhoLineShape",RhoLineShapeNew,0,2,2); 

  SetProcessParameters(Process,TMath::Abs(Eint));
  GenMUONLMR *gener = new GenMUONLMR(7,1);
//      0         1        2         3         4        5           6         7
//  "fPtPion","fPtKaon","fPtEta","fPtRho","fPtOmega","fPtPhi","fPtEtaPrime"  "fPtJPsi"
  gener->SetYParams(generYPtPar,1.,y0SG,sigySG,0.);
  gener->SetPtParams(generYPtPar,1.,TSG,MotherMass,0.);
// 2 = rho
// 3 = omega 2 body
// 4 = omega Dalitz
//
// Processes eta 2B=0  eta D=1  rho =2   omega 2B=3  omega D=4  phi=5    eta p=6       pi=7      K=8      J/psi=9
//           kEtaLMR   kEtaLMR  kRhoLMR  kOmegaLMR   kOmegaLMR  kPhiLMR  kEtaPrimeLMR  kPionLMR  kK        kJPsi
//	      
//  BR       5.8e-6    3.1e-4   4.55e-5  7.28e-5     1.3e-4     2.86e-4  1.04e-4       1         0.6344    0.05

// Added June 2019 for J/psi generation
  if(strcmp(Process,"jpsisch") == 0){
  gener->SetYParams(generYPtPar,1.,p1_xF, p2_xF, p3_xF, p4_xF);
  } else {
  gener->SetYParams(generYPtPar,1.,y0SG,sigySG,0.,0.);
  }
  
  if(strcmp(Process,"jpsi") == 0 || strcmp(Process,"jpsisch") == 0 ){
  gener->SetPtParams(generYPtPar,1.,par_p0,par_n,par_n2);
  } else {
  gener->SetPtParams(generYPtPar,1.,TSG,MotherMass,0.);
  }

  gener->GenerateSingleProcess(ProcType);
  printf("LMR generator initialization completed\n"); 

  // book some control histos, used for cuts tuning, can be skipped
  //det->BookControlHistos();
  BookHistos();
  //
  // prepare decays
  TGenPhaseSpace decay;
  TLorentzVector parent, dimuGen, muGen[2], muRec[2], muRecMS[2], dimuRecMS, dimuRec;//eli
  int npix[2]={0}, nfake[2]={0};
  float chiGlo[2] = {999.,999.};
  TParticle mupart[2];

  TTree *t1 = new TTree("t1","TLorentzVectorMuons");
  TClonesArray *arr = new TClonesArray("TParticle");
  TClonesArray &ar = *arr;
  t1->Branch("MuonsPhaseSpace",&arr,32000,1);


  Double_t prodM[2] = {KMCDetectorFwd::kMassMu, KMCDetectorFwd::kMassMu};
  const int crgMu[2] = {-1,1};
  //
  det->BookControlHistos();

  for (int iev=0;iev<nev;iev++) {
    //
    if ((iev%outN)==0) printf("Done %d out of %d\n",iev,nev);
    if (dndyBGPi>0 && (iev%refreshBg)==0) {
      if (flukaParser) {
	det->ImposeFlukaBackground(flukaParser, interactionSource, true); // allow rewind
      } else {
	det->GenBgEvent(vX,vY,vZ);
      }
    }
    //
//     double y  = dndyFunSG->GetRandom();
//     double pt = dndptFunSG->GetRandom();
//     double phi = gRandom->Rndm()*TMath::TwoPi();
//     //
//     // generate decay
//     mass=fRhoLineShape->GetRandom(); 
//     double mt = TMath::Sqrt(0.775*0.775+pt*pt);
//     parent.SetXYZM(pt*TMath::Cos(phi),pt*TMath::Sin(phi),mt*TMath::SinH(y),mass);
//     decay.SetDecay(parent, 2, prodM);
//     decay.Generate();
    //
    gener->Generate();
    
    int crgOrd = gRandom->Rndm() > 0.5 ? -1:1;
    //
    int nrec = 0;
    int nfakeHits = 0;
    double pxyz[3];
    int sext[2];
    //    if (iev==255) AliLog::SetGlobalDebugLevel(3);
    TParticle* fMu[2];
    fMu[0] = gener->GetMuon(0);
    fMu[1] = gener->GetMuon(1);
    //    printf("px0=%f py0=%f pz0=%f e0=%f m0=%f\n",fMu[0]->Px(),fMu[0]->Py(),fMu[0]->Pz(),fMu[0]->Energy(),fMu[0]->GetMass());
    //    printf("px1=%f py1=%f pz1=%f e1=%f m1=%f\n",fMu[1]->Px(),fMu[1]->Py(),fMu[1]->Pz(),fMu[1]->Energy(),fMu[1]->GetMass());
    parent.SetPxPyPzE(fMu[0]->Px()+fMu[1]->Px(),fMu[0]->Py()+fMu[1]->Py(),fMu[0]->Pz()+fMu[1]->Pz(),
		      fMu[0]->Energy()+fMu[1]->Energy());
    for (int imu=0;imu<2;imu++) {
      muGen[imu].SetXYZM(fMu[imu]->Px(),fMu[imu]->Py(),fMu[imu]->Pz(),0.105658369);//eli
    }
    dimuGen  = muGen[0];
    dimuGen += muGen[1];

    double y = parent.Rapidity();
    double pt = parent.Pt();
   
    hYGen->Fill(y);
    hPtGen->Fill(pt);
    hYPtGen->Fill(y,pt);
    //    printf("y=%f pt=%f\n",y,pt);
    for (int imu=0;imu<2;imu++) {
      //      fMu[imu] = gener->GetMuon(imu); 
      TLorentzVector* pDecMu = new TLorentzVector(0.,0.,0.,0.);
      pDecMu->SetXYZM(fMu[imu]->Px(),fMu[imu]->Py(),fMu[imu]->Pz(),fMu[imu]->GetMass());
      //      printf("iev=%d\n",iev);
      //      printf("px=%f py=%f pz=%f e=%f m=%f | phi=%e\n",pDecMu->Px(),pDecMu->Py(),pDecMu->Pz(),pDecMu->Energy(),pDecMu->M(), pDecMu->Phi());

      int crg = crgOrd*crgMu[imu];
      if (!det->SolveSingleTrack(pDecMu->Pt(),pDecMu->Rapidity(),pDecMu->Phi(),prodM[imu],crg,vX,vY,vZ, 0,1,99)) continue;  

      //     // Cluster of Muon Station
      //KMCClusterFwd* c[4];
      //for(int ch=0 ;ch<4 ;ch++) {
      //	c[ch]= det->GetLayerMS(ch)->GetMCCluster();
      //	if (!c[ch]) break; 
      //	//      sext[imu] = GetSextant(cl->GetXLab(),cl->GetYLab());
      //	sext[imu] = GetSextant(c[ch]->GetXLab(),c[ch]->GetYLab());
      //	//      printf("%d sxt %d\n",imu,sext[imu]);
      //if (ch==0 && sext[imu]<0) {
      //printf("in dead zone\n");
      //  break; // in dead zone
      //}
//       // hit map histogram of muon station qui!!!!!!!!!!!
//      hxySextantCH[ch]->Fill(c[ch]->GetXLab(),c[ch]->GetYLab());
      //    }
      
      KMCClusterFwd* cl = det->GetLayerMS(0)->GetMCCluster();
      if (!cl) break; 
      //      sext[imu] = GetSextant(cl->GetXLab(),cl->GetYLab());
      sext[imu] = GetSextant(cl->GetXLab(),cl->GetYLab());
      //      printf("%d sxt %d\n",imu,sext[imu]);
      if (sext[imu]<0) {
//	printf("in dead zone\n");
      	break; // in dead zone
      }

      // Cluster of Trigger Station
      KMCClusterFwd* t[2];
      for(int tr=0 ;tr<2 ;tr++) {
	t[tr]= det->GetLayerTR(tr)->GetMCCluster();
	if (!t[tr]) break; 
	sext[imu] = GetSextant(t[tr]->GetXLab(),t[tr]->GetYLab());
	//      printf("%d sxt %d\n",imu,sext[imu]);
	hxySextantTR[tr]->Fill(t[tr]->GetXLab(),t[tr]->GetYLab());
      }
    //   KMCClusterFwd* clITS[5];
//       for(Int_t i =0; i<5; i++){
// 	clITS[i]= det->GetLayerITS(i)->GetMCCluster();
// 	if (!clITS[i]) break; 
// 	hxySextantITS[i]->Fill(clITS[i]->GetXLab(),clITS[i]->GetYLab());
//       }



      KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw) break;
      if(trw->GetNormChi2(kTRUE)>ChiTot) continue;
      nrec++;
      nfakeHits += trw->GetNFakeITSHits();
      trw->GetPXYZ(pxyz);
      muRec[imu].SetXYZM(pxyz[0],pxyz[1],pxyz[2],prodM[imu]);
      npix[imu] = trw->GetNITSHits();
      nfake[imu] = trw->GetNFakeITSHits();
      chiGlo[imu] = trw->GetNormChi2(kTRUE);
      mupart[imu].SetMomentum(muRec[imu].Px(),muRec[imu].Py(),muRec[imu].Pz(),muRec[imu].Energy());

      KMCProbeFwd* muMS = det->GetMuBransonCorrVtx();
      if (muMS) {
	muMS->GetPXYZ(pxyz);
	muRecMS[imu].SetXYZM(pxyz[0],pxyz[1],pxyz[2],KMCDetectorFwd::kMassMu);
      }
    }

    if (nrec<2) continue;

    KMCClusterFwd* clITS[5];
    for(Int_t i =0; i<5; i++){
      clITS[i]= det->GetLayerITS(i)->GetMCCluster();
      if (!clITS[i]) break; 
      hxySextantITS[i]->Fill(clITS[i]->GetXLab(),clITS[i]->GetYLab());
    }

    hSext1VsSext0->Fill(sext[0],sext[1]);
    //if (sext[0]==sext[1]) {  
    //  printf("Dimuon killed by the same sextant exclusion\n");                                                     
    //  continue; // reject dimuon in same sextant                                                                  
    //}  
    hSext1VsSext0Trig->Fill(sext[0],sext[1]);
    dimuRec  = muRec[0];
    dimuRec += muRec[1];

    dimuRecMS = muRecMS[0];
    dimuRecMS += muRecMS[1];
    
    mupart[0].SetStatusCode(nfakeHits);
    mupart[1].SetStatusCode(nfakeHits);
    new(ar[0]) TParticle(mupart[0]);
    new(ar[1]) TParticle(mupart[1]);
    t1->Fill();

    outStream << "genrecAcc" << "gen=" << &dimuGen << "rec=" << &dimuRec << "recMS=" << &dimuRecMS
	      << "genMu0=" << &muGen[0] << "genMu1=" << &muGen[1]
	      << "recMu0=" << &muRec[0] << "recMu1=" << &muRec[1]
      	      << "recMuMS0=" << &muRecMS[0] << "recMuMS1=" << &muRecMS[1]
	      << "npix0=" << npix[0] << "npix1=" << npix[1]
	      << "nfake0=" << nfake[0] << "nfake1=" << nfake[1]
	      << "chi0=" << chiGlo[0] << "chi1=" << chiGlo[1] << "\n";
    
    
    //
    hnFake->Fill(1.*nfakeHits);
    if (nfakeHits>0) {
      hYPtFake->Fill(dimuRec.Rapidity(),dimuRec.Pt());
      hPtFake->Fill(dimuRec.Pt());
      hYPtFakeCBMbin->Fill(dimuRec.Rapidity(),dimuRec.Pt());
      hMassFake->Fill(dimuRec.M());
    }
    else if (dimuRec.M()<2.*prodM[0]) {
      printf("ev %d",iev);
      muRec[0].Print();
      muRec[1].Print();
      dimuRec.Print();
      return;
    }

    if (nfakeHits==1) {      
      hMassFake1->Fill(dimuRec.M());
    }
    else if (nfakeHits==2) {      
      hMassFake2->Fill(dimuRec.M());
    }
    else if (nfakeHits==3) {      
      hMassFake3->Fill(dimuRec.M());
    }
    else if (nfakeHits==4) {      
      hMassFake4->Fill(dimuRec.M());
    }
    else if (nfakeHits==5) {      
      hMassFake5->Fill(dimuRec.M());
    }

    hYPtAll->Fill(dimuRec.Rapidity(),dimuRec.Pt());
    hPtAll->Fill(dimuRec.Pt());
    hYPtAllCBMbin->Fill(dimuRec.Rapidity(),dimuRec.Pt());
    hMassAll->Fill(dimuRec.M());
    //hPhi->Fill(dimuRec.Phi()); 
    //printf("m=%f\n",dimuRec.M());
    if (nfakeHits==0) hMassSignalZoom->Fill(dimuRec.M());
    hY0->Fill(muRec[0].Rapidity());
    hY1->Fill(muRec[1].Rapidity());
    hY0Y1->Fill(muRec[0].Rapidity(),muRec[1].Rapidity());


  }
  //
  hMassAll->SetLineColor(kBlue);
  hMassAll->Draw();
  hMassAll->SetMinimum(0.1);
  hMassFake->SetLineColor(kRed);
  hMassFake->Draw("same");

  TFile *fout = new TFile("Matching-histos.root","recreate");
  hY0   ->Write();
  hY1   ->Write();
  hY0Y1 ->Write();
  hMassAll->Write();
  hnFake->Write();
  hMassFake->Write();
  hMassFake1->Write();
  hMassFake2->Write();
  hMassFake3->Write();
  hMassFake4->Write();
  hMassFake5->Write();
  hYPtGen->Write();
  hYPtAll->Write();
  hYPtFake->Write();
  hYGen->Write();
  hMassSignalZoom->Write();
  hYPtAllCBMbin->Write();
  hYPtFakeCBMbin->Write();
  hPtGen->Write();
  hPtAll->Write();
  hPtFake->Write();
  for(int ch=0; ch<4; ch++ )hxySextantCH[ch]->Write();
  for(int tr=0; tr<2; tr++ )hxySextantTR[tr]->Write();
  for(Int_t i=0; i<5; i++)  hxySextantITS[i]->Write();
  //hxySextantCH2->Write();
  hSext1VsSext0->Write();
  hSext1VsSext0Trig->Write();
  hPhi->Write();
  t1->Write();
  fout->Close();

}

Double_t RhoLineShapeNew(Double_t *x, Double_t *para){
  //new parameterization implemented by Hiroyuki Sako (GSI)
  Double_t mass = *x;
  double r, GammaTot;
  Double_t mRho    = 0.775; //TDatabasePDG::Instance()->GetParticle("rho0")->Mass();
  Double_t mPi     = 0.1395; //TDatabasePDG::Instance()->GetParticle("pi0")->Mass();
  Double_t mMu     = 0.1056; //TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t Gamma0  = 0.170; //TDatabasePDG::Instance()->GetParticle("rho0")->Width();

  const double Norm = 0.0744416*1.01;  

  // 0.0744416 at m = 0.72297
  // is the max number with Norm=1 (for rho)
  
  double mThreshold = 2.*mPi;

  const double T = 0.170; // Assumption of pi+ temperature [GeV/c^2]
  //const double T = 0.11; // Taken from fit to pi+ temperature [GeV/c^2]
  // with Reference: LEBC-EHS collab., Z. Phys. C 50 (1991) 405

  if (mass < mThreshold) {
    r = 0.;
    return r;
  }

  double k = sqrt(0.25*mass*mass-(mThreshold/2)*(mThreshold/2));
  double k0 = sqrt(0.25*mRho*mRho-(mThreshold/2)*(mThreshold/2));

  GammaTot = (k/k0)*(k/k0)*(k/k0)*(mRho/mass)*(mRho/mass)*Gamma0;

  double FormFactor2 = 1/((mass*mass-mRho*mRho)*(mass*mass-mRho*mRho)+
			  mass*mass*GammaTot*GammaTot);

  r = pow(mass,1.5)*pow((1-mThreshold*mThreshold/(mass*mass)),1.5)*
    ((mass*mass+2*mMu*mMu)/(mass*mass))*(pow((mass*mass-4*mMu*mMu),0.5)/mass)*FormFactor2
    *exp(-mass/T)/Norm;

  return r;
}


Int_t GetSextant(Double_t x,Double_t y)
{
  // get sextant from 1 to 6. If negative - dead zone
  const int kDeadAngle=18;
  double phi = TMath::ATan2(y,x);
  if (phi<0) phi += TMath::Pi()*2;
  int phis  = phi*TMath::RadToDeg();
  int sxt = 1 + phis/60;
  int locPhi = phis%60;
  if (locPhi<kDeadAngle/2 || locPhi>60-kDeadAngle/2) return -sxt;
  return sxt;
}

Int_t GetSextant1(Double_t x,Double_t y)
{
  
  // get sextant from 1 to 6. If negative - dead zone
  int sxt = 1;
  if((y < 20. & y > -20.) || 
     (y < (x * TMath::Tan(2.*TMath::Pi()*60./360.) + 20.) & y > (x * TMath::Tan(2.*TMath::Pi()*60./360.) -20.)) ||
     (y < (x * TMath::Tan(2.*TMath::Pi()*120./360.) +20.) & y > (x * TMath::Tan(2.*TMath::Pi()*120./360.) -20.)) ) return -100;
  
  return sxt;
}

void SetProcessParameters(const Char_t *Process, Double_t E){
  
//modified June2019 for J/psi parameterization
  if(E == 20.) {
    y0SG = 1.9;
    sigySG = 1.;
  } else if(E == 30.){
    y0SG = 2.08;
    sigySG = 1.;
  } else if(E == 40.){
    y0SG = 2.22;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.19; //PYTHIA 50 GeV
    else sigySG = 1.; 
  } else if(E == 50.){
    y0SG = 2.33;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.19; //PYTHIA 50 GeV
  }  else if(E == 60.){
    y0SG = 2.42;
    sigySG = 1.;
  } else if(E == 70.){
    y0SG = 2.50;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.28; //PYTHIA 70 GeV
  }  else if(E == 80.){
    y0SG = 2.57;
    sigySG = 1.;
  } else if(E == 90.){
    y0SG = 2.63;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.33; //PYTHIA 90 GeV
  } else if(E == 110.){
    y0SG = 2.73;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.37; //PYTHIA 110 GeV
  } else if(E == 130.){
    y0SG = 2.81;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.40; //PYTHIA 130 GeV
  } else if(E == 150.){
    y0SG = 2.88;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.42; //PYTHIA 150 GeV
  }  else if(E == 160.){
    y0SG = 2.9;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.42; //PYTHIA
    else sigySG = 1.2; //NA49
  }  else if(E == 400.){
    y0SG = 3.37;
    sigySG = 1.;
  }


  if(strcmp(Process,"etaDalitz") == 0) { // eta Dalitz
    ProcType    = 1;
    generYPtPar = 2;
    MotherMass  = 0.549;  
    if(E == 40.) {
      TSG = 0.225; 
    }
    else if(E == 160.){
      TSG = 0.24; 
    }
  } else   if(strcmp(Process,"eta2Body") == 0) { // eta Dalitz
    ProcType    = 0;
    generYPtPar = 2;
    MotherMass  = 0.549;  
    if(E == 40.) {
      TSG = 0.225; 
    }
    else if(E == 160.){
      TSG = 0.24; 
    }
  } else if(strcmp(Process,"rho") == 0) { // rho
      ProcType    = 2;
      generYPtPar = 3; 
      MotherMass  = 0.775; 
      if(E == 40.) {
	TSG = 0.25; 
      }
      else if(E == 160.){
	TSG = 0.29; 
      }
  }
  else if(strcmp(Process,"omega2Body") == 0) { // omega 2 body
    ProcType    = 3;
    generYPtPar = 4; 
    MotherMass  = 0.781;  
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
	TSG = 0.29; 
    }
  } 
  else if(strcmp(Process,"omegaDalitz") == 0) { // omega Dalitz
    ProcType    = 4;
    generYPtPar = 4; 
    MotherMass  = 0.781;
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
      TSG = 0.29; 
    }   
  } 
  else if(strcmp(Process,"phi") == 0) { // phi
    ProcType    = 5;
    generYPtPar = 5; 
    MotherMass  = 1.02; 
    if(E == 40.) {
      TSG =0.25;           
    }
    else if(E == 160.){
      TSG = 0.30; //NA49
    }  
  }   else if(strcmp(Process,"etaPrime") == 0) { // eta prime
    ProcType    = 6;
    generYPtPar = 6; 
    MotherMass  = 0.958;  
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
      TSG = 0.29; 
    } 
// modified June 2019 J/psi parameterizations    
  }  else if(strcmp(Process,"jpsi") == 0 || strcmp(Process,"jpsisch") == 0) { // J/psi
    ProcType    = 9;
    if(strcmp(Process,"jpsi") == 0 ) {
    ProcType = 9;
    generYPtPar = 7;
    } else if(strcmp(Process,"jpsisch") == 0) {
    ProcType=10;
    generYPtPar = 8;
    }
    MotherMass  = 3.0969; 
//     if(E == 40.) {
//       TSG = 0.25; 
//     }
//     else if(E == 160.){
//       TSG = 0.284; //NA50
//     }   
    if(E == 150.){   
       par_p0=5.34;
       par_n=18.8;
       par_n2=2.70;
    } else if(E==130.){
       par_p0=13.96;
       par_n=201.1;
       par_n2=2.58;
    } else if(E==110.){
       par_p0=3.73;
       par_n=8.83;
       par_n2=2.84;
    } else if(E==90.){
       par_p0=4.89;
       par_n=17.3;
       par_n2=2.80;
    } else if(E==70.){
       par_p0=3.86;
       par_n=10.6;
       par_n2=2.92;
    } else if(E==50.){
       par_p0=11.03;
       par_n=429.7;
       par_n2=3.27;
    } else if(E==40.){ // use same parameter as for E=40 GeV 
       par_p0=11.03;
       par_n=429.7;
       par_n2=3.27;
    }
    
   }
   
        if(strcmp(Process,"jpsisch") == 0) { // J/psi
      ProcType=10;
      generYPtPar = 8; 
      MotherMass  = 3.0969; 
       
       Double_t a = 13.5;
       Double_t b = 44.9;
       p4_xF = y0SG;
       p1_xF = 3.096; //mjpsi
    if(E == 150.){   
      p2_xF = 16.8;  //sqrts
    }  else if(E == 130.){
      p2_xF = 15.7;
    }  else if(E == 110.){
      p2_xF = 14.4;
    }  else if(E == 90.){
      p2_xF = 13.1;
    }  else if(E == 70.){
      p2_xF = 11.5;
    }  else if(E == 50.){
      p2_xF = 9.77;
    }  else if(E == 40.){
      p2_xF = 8.76;
    }
      p3_xF = a/(1+b/p2_xF);
      printf("p1_xF=%f, p2_xF=%f, p3_xF=%f, p4_xF=%f\n",p1_xF,p2_xF,p3_xF,p4_xF);
  }   
  

  printf("Generating %s\n",Process);
  printf("Particle mass = %f\n",MotherMass);
  printf("Inverse slope = %f\n",TSG);
  printf("rapidity sigma = %f\n",sigySG);

}

void BookHistos(){

  // Book some histos

  //hYPtGen = new TH2F("YPTGen","Y-Pt corr match",40,1.8,4.0,30,ptminSG,ptmaxSG);
  hPtGen = new TH1F("PTGen","Pt gen",30,ptminSG,ptmaxSG);
  hYGen = new TH1F("hYGen","Y full phase space",100.,yminSG,ymaxSG);
  //hYPtFake = new TH2F("YPTFake","Y-Pt fake match",40,1.8,4.0,30,ptminSG,ptmaxSG);
  //hYPtAll = new TH2F("YPTAll","Y-Pt all match",40,1.8,4.0,30,ptminSG,ptmaxSG);
  hPtFake = new TH1F("PTFake","Pt fake match",30,ptminSG,ptmaxSG);
  hPtAll = new TH1F("PTAll","Pt all match",30,ptminSG,ptmaxSG);
 
  hYPtGen = new TH2F("YPTGen","Y-Pt corr match",80,1.0,5.4,30,ptminSG,ptmaxSG);
  hYPtFake = new TH2F("YPTFake","Y-Pt fake match",80,1.0,5.4,30,ptminSG,ptmaxSG);
  hYPtAll = new TH2F("YPTAll","Y-Pt all match",80,1.0,5.4,30,ptminSG,ptmaxSG);
  
  hnFake = new TH1F("hnFake","number of fake matches",20,0.,20.);

  if(ProcType == 9) {
    hMassFake = new TH1F("MassFake","Mass fake match",100,2.,3.5);
    hMassFake1 = new TH1F("MassFake1","Mass fake match 1",100,2.,3.5);
    hMassFake2 = new TH1F("MassFake2","Mass fake match 2",100,2.,3.5);
    hMassFake3 = new TH1F("MassFake3","Mass fake match 3",100,2.,3.5);
    hMassFake4 = new TH1F("MassFake4","Mass fake match 4",100,2.,3.5);
    hMassFake5 = new TH1F("MassFake5","Mass fake match 5",100,2.,3.5);
    hMassAll = new TH1F("MassAll","Mass all match",100,2.,3.5);
  } else {
    hMassFake = new TH1F("MassFake","Mass fake match",100,0.,2.);
    hMassFake1 = new TH1F("MassFake1","Mass fake match 1",100,0.,2.);
    hMassFake2 = new TH1F("MassFake2","Mass fake match 2",100,0.,2.);
    hMassFake3 = new TH1F("MassFake3","Mass fake match 3",100,0.,2.);
    hMassFake4 = new TH1F("MassFake4","Mass fake match 4",100,0.,2.);
    hMassFake5 = new TH1F("MassFake5","Mass fake match 5",100,0.,2.);
    hMassAll = new TH1F("MassAll","Mass all match",100,0.,2.0);
  }

  hMassSignalZoom = new TH1F("hMassSignalZoom","Mass signal",50,0.7,0.85);
  hxySextantCH[0] = new TH2F("hxySextantCH0","y vs y in sextant",100.,-150.,150.,100.,-150.,150.);
  hY0   = new TH1F("Y0","Y muon 0",40,1.8,4.0);
  hY1   = new TH1F("Y1","Y muon 1",40,1.8,4.0);
  hY0Y1 = new TH2F("Y0Y1","Y muon 0 vs Y muon 1",40,1.8,4.0,40,1.8,4.0);
  hxySextantCH[1] = new TH2F("hxySextantCH1","y vs y in sextant",100.,-200.,200.,100.,-200.,200.);
  hxySextantCH[2] = new TH2F("hxySextantCH2","y vs y in sextant",100.,-300.,300.,100.,-300.,300.);
  hxySextantCH[3] = new TH2F("hxySextantCH3","y vs y in sextant",100.,-350.,350.,100.,-350.,350.);
  hxySextantTR[0] = new TH2F("hxySextantTR1","y vs y in sextant",100.,-500.,500.,100.,-500.,500.);
  hxySextantTR[1] = new TH2F("hxySextantTR2","y vs y in sextant",100.,-500.,500.,100.,-500.,500.);
  for(Int_t i=0; i<5; i++){
    hxySextantITS[i] = new TH2F(Form("hxySextantITS_%d",i),"x vs y in sextant",400,-20.,20.,400,-20.,20.);
  }
  hSext1VsSext0 = new TH2F("hSext1VsSext0","sext1 vs sext2",10,0,9,10,0,9);
  hSext1VsSext0Trig = new TH2F("hSext1VsSext0Trig","sext1 vs sext2 trig cut",10,0,9,10,0,9);
  
  hYPtAllCBMbin = new TH2F("hYPtAllCBMbin","Y-Pt all match",40,0.,4.,40,0.,2.);
  hYPtFakeCBMbin = new TH2F("hYPtFakeCBMbin","Y-Pt fake match",40,0.,4.,40,0.,2.);
  
  hPhi = new TH1F("hPhi","phi",100,0.,TMath::TwoPi());

}

void CalcBkgPar(Double_t E){

  if(E < 0.){
    printf("Beam E is negative, diable BCKG generation\n");
    y0BG   = 1.9;   // gaussian y mean - 40 GeV
    sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 0.;
    dndyBGK = 0.;
    dndyBGP = 0.;
    TBGpi = 0.17;
    TBGK = 0.22;
    TBGP = 0.25;
  }
  else if(E == 20.){
    y0BG   = 1.9;   // gaussian y mean - 40 GeV
    sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 410.;
    dndyBGK = 51.;
    dndyBGP = 148.;
    TBGpi = 0.17;
    TBGK = 0.22;
    TBGP = 0.25;
  }
  else if(E == 40.){
    y0BG   = 2.22;   // gaussian y mean - 40 GeV
    y0BGPi  = 0.666;  
    y0BGKplus   = 0.694;  
    y0BGKminus  = 0.569;  
    y0BGP   = 0.907; 
    sigyBGPi = 0.872;
    sigyBGKplus  = 0.725;
    sigyBGKminus = 0.635;
    sigyBGP = 0.798;
    sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    NBGPi     = 74.;
    NBGKplus  = 16.2;
    NBGKminus = 6.03;
    NBGP      = 37.5;
    Piratio = 0.91;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
  }
  else if(E == 60.){
    y0BG   = 2.42;   // gaussian y mean - 60 GeV
    sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
  }
 else if(E == 80.){
    y0BG   = 2.57;   // gaussian y mean - 80 GeV
    sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
 }

 else if(E == 160.){//NA49
   y0BG   = 2.9;   // gaussian y mean - 160 GeV
   y0BGPi  = 0.72;  
   y0BGKplus   = 0.839;  
   y0BGKminus   = 0.727;  
   y0BGP   = 39.8; 
   sigyBGPi = 1.18;
   sigyBGKplus = 0.88;
   sigyBGKminus = 0.81;
   sigyBGP = 8.07;
   yminBG = 1.5;   // min y to generate 
   ymaxBG = 4.5;   // 
   ptminBG = 0.01;
   ptmaxBG = 3;
   dndyBGPi = 1258.;
   dndyBGK = 155.;
   dndyBGP = 292.;
   NBGPi     = 107.6;
   NBGKplus  = 23.4;
   NBGKminus = 12.8;
   NBGP      = 2.55e+06;
   Piratio = 0.97;   
   TBGpi = 0.18; // inv.slope of thermal pt distribution
   TBGK = 0.23; // inv.slope of thermal pt distribution
   TBGP = 0.31; // inv.slope of thermal pt distribution
  }




















 else if(E == 400.){
   y0BG   = 3.37;   // gaussian y mean - 80 GeV
   sigyBG = 1.2;   // .. sigma
    yminBG = 1.5;   // min y to generate 
    ymaxBG = 4.5;   // 
    TBG    = 0.17;  // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
  }



  printf("Simulating E=%f interactions\n",E);
  printf("pions:   total multiplicity = %f ; T = %f\n",dndyBGPi,TBGpi);
  printf("kaons:   total multiplicity = %f ; T = %f\n",dndyBGK,TBGK);
  printf("protons: total multiplicity = %f ; T = %f\n",dndyBGP,TBGP);

}
