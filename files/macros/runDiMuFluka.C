#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLorentzVector.h>
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "KMCLayerFwd.h"
#include "KMCFlukaParser.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TTreeStream.h"
#endif

// Track Chi2Tot cut
double ChiTot = 1.5;

// settings for backgroung generation. Used just to create occupancy in the tracker,
// no need to this f.ph.sp...
// default values - recalculated in CalcBkgPar
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

double vX=0,vY=0,vZ=0; // event vertex
KMCDetectorFwd* det = 0;

KMCFlukaParser flukaParser;
  
//
void  Fluka2Particle(Int_t fCode, double &mass, int &charge);
void  SetupFlukaParticle(const FlukaPart& part);

Int_t GetSextant(Double_t x,Double_t y);
void CalcBkgPar(Double_t E);
//====================================================================

void runMuFluka(double Eint=40., // Elab energy 
		const char* setup="setup.txt", // setup to load
		const char* lstName="inpFluka.txt", // list of fluka files to process
		int numITS=5,         // N ITS
		int maxAcc=-1,        // max trackables to accept ( no limit if <0)
		int maxMissITS=0,   // N ITS  hits allowed to be absent in fluka record
		int maxMissMS =0,     // N MS   hits allowed to be absent in fluka record
		int maxMissTR =0,    // N Trig hits allowed to be absent in fluka record
		int refreshBg=10,    // generate new bg event for each refreshBg-th
		Int_t kNoITS = 0
		)
{
  //  gROOT->Macro("LoadLibs.C");
  TTreeSRedirector outStream("dimuFluka.root"); // the output stream trees will go here
  
  CalcBkgPar(Eint);

  double mass = KMCDetectorFwd::kMassMu;  // particle to generate
  //
  det = new KMCDetectorFwd();
  det->ReadSetup(setup,setup);
  //det->InitBgGeneration(dndyBG,y0BG,sigyBG,yminBG,ymaxBG,TBG,ptminBG,ptmaxBG);
  
  det->InitBgGenerationPart(NBGPi,NBGKplus,NBGKminus,NBGP,Piratio,y0BG,y0BGPi,y0BGKplus,y0BGKminus,y0BGP,sigyBGPi,sigyBGKplus,sigyBGKminus,sigyBGP,yminBG,ymaxBG,TBGpi,TBGK,TBGP,ptminBG,ptmaxBG);
  //det->InitBgGenerationPart(dndyBGPi,dndyBGK,dndyBGP,y0BG,sigyBG,yminBG,ymaxBG,TBGpi,TBGK,TBGP,ptminBG,ptmaxBG);
  //  det->InitBgGenerationPart(0.,dndyBGK,0.,y0BG,sigyBG,yminBG,ymaxBG,TBGpi,TBGK,TBGP,ptminBG,ptmaxBG);

  printf("pion   multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChPi());
  printf("kaon   multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChK());
  printf("proton multiplicity in %f<y<%f = %f\n",yminBG,ymaxBG,det->GetNChP());

  //
  // set the min N tracker hits needed to validate the track
  if(kNoITS == 0) {
    det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); 
    printf("running reconstruction for trigger rate estimate (neglect ITS hits)\n");
  } else if(kNoITS == 1) {
    ChiTot = 6.;
    det->SetMinITSHits(-1); // det->GetNumberOfActiveLayersITS()-1);  // for negative - ignore its hits
    printf("running standard reconstruction\n");
  }
  det->SetMinMSHits(det->GetNumberOfActiveLayersMS()); 
  det->SetMinTRHits(det->GetNumberOfActiveLayersTR()); 
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000); 
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(-20);  // fiducial cut on chi2 of convergence to vtx
  //
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1);
  det->SetApplyBransonPCorrection(0.1); 
  // book some control histos, used for cuts tuning, can be skipped
  //det->BookControlHistos();

  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  
  //
  if (!flukaParser.SetInpList(lstName)) return; // init fluka list
  det->SetExternalInput(kTRUE);
  //
  //
  int minPix = det->GetNumberOfActiveLayersITS() - maxMissITS;
  int minMS  = det->GetNumberOfActiveLayersMS()  - maxMissMS;
  int minTR  = det->GetNumberOfActiveLayersTR()  - maxMissTR;
  //
  
  const FlukaStat& flStat = flukaParser.GetStat();
  const std::vector<FlukaPart>& flukaParts = flukaParser.GetParticles();
  
  while ( (maxAcc<0 || flStat.totalRead<maxAcc) &&
	  flukaParser.GetNextGoodPair(minPix,minMS,minTR) ) {
    //
    if (dndyBGPi>0 && (flStat.totalAccepted%refreshBg)==0) det->GenBgEvent(vX,vY,vZ);
    //
    int nrec = 0, npix[2]={0}, nfake[2]={0};
    double pxyz[3];     

    TLorentzVector muRec[2],muGen[2], dimuGen, dimuRec;

    for (int imu=0;imu<2;imu++) {
      SetupFlukaParticle( flukaParts[imu] );
      
      KMCProbeFwd* genProbe = det->GetProbe(); // save generated muon kinem
      genProbe->GetPXYZ(pxyz);
      muGen[imu].SetXYZM(pxyz[0],pxyz[1],pxyz[2],KMCDetectorFwd::kMassMu);
      
      if (!det->SolveSingleTrackViaKalmanMC(999)) break;
      // access winner track
      KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack();  // save reconstructed muon kinem
      if (!trw) break;
      if(trw->GetNormChi2(kTRUE)>ChiTot) continue;
      nrec++;
      npix[imu] = trw->GetNITSHits();
      nfake[imu] = trw->GetNFakeITSHits();
      trw->GetPXYZ(pxyz);
      muRec[imu].SetXYZM(pxyz[0],pxyz[1],pxyz[2],KMCDetectorFwd::kMassMu);
    }
    if (nrec<2) continue;

    dimuGen = muGen[0];
    dimuRec = muRec[0];
    dimuGen+= muGen[1];
    dimuRec+= muRec[1];
    int ngenEv = flStat.totalRead;
    int nrecEv = flStat.totalAccepted;
    outStream << "genrecAcc" << "gen=" << &dimuGen << "rec=" << &dimuRec << "ngen=" << ngenEv << "nrec=" << nrecEv <<"\n";
    printf("GenMass: %.3f RecMass: %.3f\n",dimuGen.M(), dimuRec.M());
  }
  //
}



void Fluka2Particle(Int_t fCode, double &mass, int &charge)
{
  switch(fCode) {
  case 1:  mass = 0.94;                    charge = 1; break; // proton
  case 2:  mass = 0.94;                    charge =-1; break; // a-proton
  case 3:  mass = 0.0005;                  charge =-1; break; // e-
  case 4:  mass = 0.0005;                  charge = 1; break; // e+
  case 8:  mass = 0.94;                    charge = 0; break; // neutron 
  case 10: mass = KMCDetectorFwd::kMassMu; charge = 1; break; // mu+
  case 11: mass = KMCDetectorFwd::kMassMu; charge =-1; break; // mu-
  case 12: mass = KMCDetectorFwd::kMassK;  charge = 0; break; // KLong
  case 13: mass = KMCDetectorFwd::kMassPi; charge = 1; break; // pi+
  case 14: mass = KMCDetectorFwd::kMassPi; charge =-1; break; // pi-
  case 15: mass = KMCDetectorFwd::kMassK;  charge = 1; break; // K+
  case 16: mass = KMCDetectorFwd::kMassK;  charge =-1; break; // K-    
  default: printf("Unknown particle code %d\n",fCode); exit(1);    
  }
}

//____________________________________________________
void SetupFlukaParticle(const FlukaPart& part)
{
  // set input particle and fluka hits
  KMCLayerFwd* lr = 0;
  KMCProbeFwd* anProbe = 0;
  double r2,xyzTr[3];
  double rx,ry,mass, massImp = KMCDetectorFwd::kMassMu;
  int charge;
  anProbe = det->GetProbe();
  Fluka2Particle(part.codeOr,mass,charge);
  anProbe->ImposeKinematics(&part.recDataPrim[kX], &part.recDataPrim[kCX], part.recDataPrim[kE],mass,charge);
  //
  int maxLr = -1;
  for (int ih=0;ih<det->GetNumberOfActiveLayersITS();ih++) {
    if (part.recTypePix[ih]!=kRecDummy) {
      lr = det->GetLayerITS(ih);
      lr->GetMCCluster()->Kill();
      r2 = part.recDataPix[ih][kX]*part.recDataPix[ih][kX]+part.recDataPix[ih][kY]*part.recDataPix[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectITS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &part.recDataPix[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(),xyzTr[2]-rx*lr->GetXRes(),-1);
      anProbe = lr->GetAnProbe();
      Fluka2Particle(part.recTypePix[ih],mass,charge);
      //      printf("\nSetLr PIX "); lr->Print();
      anProbe->ImposeKinematics(&part.recDataPix[ih][kX],&part.recDataPix[ih][kCX],part.recDataPix[ih][kE],massImp,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  for (int ih=0;ih<det->GetNumberOfActiveLayersMS();ih++) {
    if (part.recTypeMS[ih]!=kRecDummy) {
      lr = det->GetLayerMS(ih);
      lr->GetMCCluster()->Kill();
      r2 = part.recDataMS[ih][kX]*part.recDataMS[ih][kX]+part.recDataMS[ih][kY]*part.recDataMS[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectMS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &part.recDataMS[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(), xyzTr[2]-rx*lr->GetXRes(),-1);
      Fluka2Particle(part.recTypeMS[ih],mass,charge);
      anProbe = lr->GetAnProbe();
      //      printf("\nSetLr MS  "); lr->Print();
      anProbe->ImposeKinematics(&part.recDataMS[ih][kX],&part.recDataMS[ih][kCX],part.recDataMS[ih][kE],massImp,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  for (int ih=0;ih<det->GetNumberOfActiveLayersTR();ih++) {
    if (part.recTypeTR[ih]!=kRecDummy) {
      lr = det->GetLayerTR(ih);
      lr->GetMCCluster()->Kill();
      r2 = part.recDataTR[ih][kX]*part.recDataTR[ih][kX]+part.recDataTR[ih][kY]*part.recDataTR[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectMS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &part.recDataTR[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(), xyzTr[2]-rx*lr->GetXRes(),-1);
      Fluka2Particle(part.recTypeTR[ih],mass,charge);
      anProbe = lr->GetAnProbe();
      //      printf("\nSetLr TR  "); lr->Print();
      anProbe->ImposeKinematics(&part.recDataTR[ih][kX],&part.recDataTR[ih][kCX],part.recDataTR[ih][kE],massImp,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  det->SetLastActiveLayerTracked(maxLr);
  //
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

void CalcBkgPar(Double_t E){

  //  if(E == 20.){
  //    y0BG   = 1.9;   // gaussian y mean - 40 GeV
  //    sigyBG = 1.2;   // .. sigma
  //    yminBG = 1.5;   // min y to generate 
  //    ymaxBG = 4.5;   // 
  //    TBG    = 0.17;  // inv.slope of thermal pt distribution
  //    ptminBG = 0.01;
  //    ptmaxBG = 3;
  //    dndyBGPi = 410.;
  //    dndyBGK = 51.;
  //    dndyBGP = 148.;
  //    TBGpi = 0.17;
  //   TBGK = 0.22;
  //    TBGP = 0.25;
  //  }
  //  else if(E == 40.){
  //    y0BG   = 2.22;   // gaussian y mean - 40 GeV
  //    sigyBG = 1.2;   // .. sigma
  //    yminBG = 1.5;   // min y to generate 
  //   ymaxBG = 4.5;   // 
  //   TBG    = 0.17;  // inv.slope of thermal pt distribution
  //    ptminBG = 0.01;
  //    ptmaxBG = 3;
  //    dndyBGPi = 615.;
  //    dndyBGK = 78.;
  //    dndyBGP = 150.;
  //    TBGpi = 0.17;
  //    TBGK = 0.23;
  //    TBGP = 0.25;
  //  }

  if(E == 40.){
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






  printf("Simulating E=%f interactions\n",E);
  printf("pions:   total multiplicity = %f ; T = %f\n",dndyBGPi,TBGpi);
  printf("kaons:   total multiplicity = %f ; T = %f\n",dndyBGK,TBGK);
  printf("protons: total multiplicity = %f ; T = %f\n",dndyBGP,TBGP);

}
