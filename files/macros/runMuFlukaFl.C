#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLorentzVector.h>
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "KMCLayerFwd.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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

TH2F *hYPtCorr=0,*hYPtFake=0, *hYPtGen=0;
TH3F *hPxPyPz=0, *hPxPyPzPos=0, *hPxPyPzNeg=0;
TH1F *hZDecReco=0; // Z of reconstructed decays
TH1F *hCharge;
TH2F *hxySextantMS0;
//TH2F *hxySextantITS[5];
//TH2F *hxySextantITS[10];
TH1F *hSextant;
TH1F *hRecType;
TH1F *hTrigType;
TH1F *hNormChi2;

TH1F *hFlMuZAll=0, *hFlMuZFake=0;

TTreeSRedirector* fDebugStreamer = 0;

//====================================================================
//
// These are settings for the fluka parser
//
enum {kRecDummy=-9999,kMaxPix=20,kMaxMS=20,kMaxTrig=10};
enum {kE,kX,kY,kZ,kCX,kCY,kCZ,kNFld};

Int_t GetNextGoodParticle(Int_t minPix=0,Int_t minMS=3,Int_t minTr=2);
Int_t SetInpList(const char* list);
//
char* readNextRecord(ifstream& strm);
Int_t readNextParticle(ifstream& strm,Bool_t verbose=kFALSE);
void  FlukaParticle(Int_t fCode, double &mass, int &charge);
void  ImposeKinematics(KMCProbeFwd* probe, double* xyzLab,double* cosinesLab, double en, double mass, int charge);
void  SetupFlukaParticle();

TObjArray inpFileList;
Int_t     nFiles=0;

int    codeOr=0; // particle fluka code
double massOr;
int    chargeOr;
double recDataPrim[kNFld];
int    recTypePix[kMaxPix];
double recDataPix[kMaxPix][kNFld];
int    recTypeMS[kMaxMS];
double recDataMS[kMaxMS][kNFld];
int    recTypeTR[kMaxTrig];
double recDataTR[kMaxTrig][kNFld];
int nPix=0,nMS=0,nTrig=0;
int totalRead = 0;
int totalAccepted = 0;
float zMuFirst = 2e6; // Z of 1st muon observation

Int_t GetSextant(Double_t x,Double_t y);
void CalcBkgPar(Double_t E);
//====================================================================

void runMuFlukaFl(double Eint=40., // Elab energy 
		  const char* setup="setup.txt", // setup to load
		  const char* lstName="inpFluka.txt", // list of fluka files to process
		  int numITS=5,         // N ITS
		  int maxAcc=-1,        // max trackables to accept ( no limit if <0)
		  int maxMissITS=999,   // N ITS  hits allowed to be absent in fluka record
		  int maxMissMS =0,     // N MS   hits allowed to be absent in fluka record
		  int maxMissTR =0,    // N Trig hits allowed to be absent in fluka record
		  int refreshBg=10,    // generate new bg event for each refreshBg-th
		  Int_t kNoITS = 0
		  )
{
  //  gROOT->Macro("LoadLibs.C");
  CalcBkgPar(Eint);

  fDebugStreamer = new TTreeSRedirector("flk_debug.root","recreate");
  
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
  if (!SetInpList(lstName)) return; // init fluka list
  det->SetExternalInput(kTRUE);
  //

  TH2F *hxySextantITS[numITS];

  // Book some histos
  hYPtGen  = new TH2F("YPTGen","Y-Pt  Gen",       20,1.5,4.5,32,0,3.2);
  hYPtCorr = new TH2F("YPTCorr","Y-Pt corr match",20,1.5,4.5,32,0,3.2);
  hYPtFake = new TH2F("YPTFake","Y-Pt fake match",20,1.5,4.5,32,0,3.2);
  hPxPyPz  = new TH3F("hPxPyPz","Px, Py, Pz",40,-2.,2.,40,-2.,2.,290,1.,30.);
  hPxPyPzPos  = new TH3F("hPxPyPzPos","Px, Py, Pz",40,-2.,2.,40,-2.,2.,290,1.,30.);
  hPxPyPzNeg  = new TH3F("hPxPyPzNeg","Px, Py, Pz",40,-2.,2.,40,-2.,2.,290,1.,30.);
  hxySextantMS0 = new TH2F("hxySextantMS0","x vs y in sextant",300,-150.,150.,300,-150.,150.);
  //for(Int_t i=0; i<5; i++){
  for(Int_t i=0; i<numITS; i++){
    hxySextantITS[i] = new TH2F(Form("hxyITS_%d",i),"x vs y ",400,-20.,20.,400,-20.,20.);
  }

  hSextant = new TH1F("hSextant","sextant",10,0.,9.);
  hZDecReco = new TH1F("zdec","zdecay of reconstruced",200,0.,400.);
  hCharge = new TH1F("hCharge","charge of reconstruced",3,-1.,2.);
  hRecType = new TH1F("hRecType","reconstructed fluka part id",20,1,20);
  hTrigType = new TH1F("hTrigType","Trig fluka part id",20,1,20);
  hNormChi2 = new TH1F("hNormChi2","norm chi2 track",100,0.,10.);

  int zMax = int( det->GetLayer( det->GetLastActiveLayer() )->GetZ() );
  
  hFlMuZAll = new TH1F("zMuAll","zMuAll",zMax,0,zMax);
  hFlMuZFake = new TH1F("zMuFake","zMuFake",zMax,0,zMax);
  
  TTree *t1 = new TTree("t1","RecSingleMuon");

  //
  int minPix = det->GetNumberOfActiveLayersITS() - maxMissITS;
  int minMS  = det->GetNumberOfActiveLayersMS()  - maxMissMS;
  int minTR  = det->GetNumberOfActiveLayersTR()  - maxMissTR;
  //
  TLorentzVector murec;
  TParticle *mu = new TParticle();

  t1->Branch("Data","TParticle",&mu,64000,0);

  det->BookControlHistos();
  while ( (maxAcc<0 || totalRead<maxAcc) && GetNextGoodParticle(minPix,minMS,minTR) ) {
    //  while ( (maxAcc<0 || totalAccepted<maxAcc) && GetNextGoodParticle(minPix,minMS,minTR) ) {
    //
    if (dndyBGPi>0 && (totalAccepted%refreshBg)==0) det->GenBgEvent(vX,vY,vZ);
    //
    KMCProbeFwd* gen = det->GetProbe();
    double pxyz[3], pxyzRec[3];     
    int sext;
    gen->GetPXYZ(pxyz);

    murec.SetXYZM(pxyz[0],pxyz[1],pxyz[2],KMCDetectorFwd::kMassMu);
    double yg = murec.Rapidity();
    double ptg = murec.Pt();
    hYPtGen->Fill(yg,ptg);
    //
    SetupFlukaParticle();
    hTrigType->Fill(recTypeTR[1]);

    TList* lstLay = det->GetLayers();
    TIter nextL(lstLay);
    KMCLayerFwd* lrq = 0;
    while ( (lrq=(KMCLayerFwd*)nextL()) ) {
      KMCClusterFwd* cl = lrq->GetMCCluster();
      if (cl->IsKilled()) continue; // no cluster
      int lrType = lrq->GetType();
      int lrID = lrq->GetActiveID();
      float xcl = cl->GetXLab(),ycl = cl->GetYLab(), zcl = cl->GetZLab();
      (*fDebugStreamer)<<"clFl"<< "lrType=" << lrType << "lrID=" << lrID
		       << "x=" << xcl << "y=" << ycl << "z=" << zcl << "\n";
    }
    
    if (!det->SolveSingleTrackViaKalmanMC(999)) continue;
    //    return;
    KMCClusterFwd* cl = det->GetLayerMS(0)->GetMCCluster();
    if (!cl) break; 
    sext = GetSextant(cl->GetXLab(),cl->GetYLab());
    //      printf("%d sxt %d\n",imu,sext[imu]);
    if (sext<0) {
      printf("in dead zone\n");
      continue; // in dead zone
    }
    hxySextantMS0->Fill(cl->GetXLab(),cl->GetYLab());
    // KMCClusterFwd* clITS[5];
//     for(Int_t i =0; i<5; i++){
//       clITS[i]= det->GetLayerITS(i)->GetMCCluster();
//       if (!clITS[i]) break; 
//       hxySextantITS[i]->Fill(clITS[i]->GetXLab(),clITS[i]->GetYLab());
//     }
    // access winner track
    KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack();
    //    det->Print("cl");
    //return;
    if(trw) hNormChi2->Fill(trw->GetNormChi2(kTRUE));
    if (trw && trw->GetNormChi2(kTRUE)<ChiTot) {
      printf("z dec = %f\n",det->GetZDecay());
      printf("ZMu1stSeen = %f\n",zMuFirst);
      hZDecReco->Fill(det->GetZDecay());
      //      trw->Print("etp");
      //      double pxyz[3];     
      //trw->GetPXYZ(pxyz);
      trw->GetPXYZ(pxyzRec);
      murec.SetXYZM(pxyzRec[0],pxyzRec[1],pxyzRec[2],KMCDetectorFwd::kMassMu);
      mu->SetMomentum(pxyzRec[0],pxyzRec[1],pxyzRec[2],murec.Energy());
      if(trw->GetCharge() < 0.) mu->SetPdgCode(13);
      if(trw->GetCharge() > 0.) mu->SetPdgCode(-13);
      //if (mu->Eta()>5.5) {
      //mu->Print();
      //return;
      //}
      double y  = murec.Rapidity();
      double pt = murec.Pt();
      hPxPyPz->Fill(pxyzRec[0],pxyzRec[1],pxyzRec[2]);
      if(trw->GetCharge() < 0.) hPxPyPzNeg->Fill(pxyzRec[0],pxyzRec[1],pxyzRec[2]);
      if(trw->GetCharge() > 0.) hPxPyPzPos->Fill(pxyzRec[0],pxyzRec[1],pxyzRec[2]);
      int nfITS = trw->GetNFakeITSHits();
      mu->SetUniqueID(nfITS);
      if (nfITS > 0) {
	hYPtFake->Fill(y,pt);
	hFlMuZFake->Fill(zMuFirst);
      }
      // KMCClusterFwd* clITS[5];
//       for(Int_t i =0; i<5; i++){
      KMCClusterFwd* clITS[numITS];
      for(Int_t i =0; i<numITS; i++){
	clITS[i]= det->GetLayerITS(i)->GetMCCluster();
	if (!clITS[i]) break; 
	hxySextantITS[i]->Fill(clITS[i]->GetXLab(),clITS[i]->GetYLab());
      }
      hYPtCorr->Fill(y,pt);
      hCharge->Fill(trw->GetCharge());
      hSextant->Fill(sext);
      hRecType->Fill(recTypeTR[1]);
      mu->SetFirstMother(recTypeTR[1]);
      hFlMuZAll->Fill(zMuFirst);
      t1->Fill();
    }
    //    else return;
  }
  //
  hYPtCorr->Draw("box");
  //
  TFile *fout;
  if(kNoITS == 0) fout = new TFile("Matching-histos-Mu.root","recreate");
  if(kNoITS == 1) fout = new TFile("Matching-histos-Mu-Trig.root","recreate");
  hYPtGen->Write();
  hYPtCorr->Write();
  hYPtFake->Write();
  hPxPyPz->Write();
  hPxPyPzPos->Write();
  hPxPyPzNeg->Write();
  hZDecReco->Write();
  hCharge->Write();
  hxySextantMS0->Write();
  //for(Int_t i=0; i<5; i++) hxySextantITS[i]->Write();
  for(Int_t i=0; i<numITS; i++) hxySextantITS[i]->Write();
  hSextant->Write();
  t1->Write();
  hRecType->Write();
  hTrigType->Write();
  hFlMuZAll->Write();
  hFlMuZFake->Write();  
  //
  det->GetHChi2LrCorr()->Write();
  det->GetHChi2VCorr()->Write();
  det->GetHChi2VFake()->Write();
  det->GetHChi2MS()->Write();
  det->GetHChi2NDFCorr()->Write();
  det->GetHChi2NDFFake()->Write();
  det->GetHNCand()->Write();
  det->GetHCandCorID()->Write();
  hNormChi2->Write();
  fout->Close();
  delete fDebugStreamer;

}

//=====================================================================
//=====================================================================
//=======                 FLUKA OUTPUT PARSER              ============
//=====================================================================
//=====================================================================


Int_t GetNextGoodParticle(Int_t minPix,Int_t minMS,Int_t minTr)
{
  static int curFile=-1;
  static ifstream* inpf = 0;
  //
  // prepare input stream if needed
  while(1) { // loop over files
    if (!inpf) {
      if (curFile>=nFiles-1) return 0;
      const char* fnm = inpFileList.At(++curFile)->GetName();
      printf("Processing %d-th file %s\n",curFile,fnm);
      inpf = new ifstream(fnm);
      if (!inpf->good()) {
	printf("Failed to open file %s\n",fnm);
	exit(1);
      }
    }
    while (readNextParticle(*inpf)) {
      totalRead++;
      if (nPix>=minPix && nMS>=minMS && nTrig>=minTr) {
	totalAccepted++;
	printf("Read %d Accepted %d | %d %d %d\n",totalRead,totalAccepted,nPix,nMS,nTrig);
	return 1;
      }
    }
    delete inpf;
    inpf = 0;
  }
  return 0;
}


Int_t readNextParticle(ifstream& strm,Bool_t verbose)
{
  //
  static TString recSav;  
  TString rec;
  double en,x,y,z,cx,cy,cz;
  int code,st;
  nPix=0;
  nMS=0;
  nTrig=0;
  //
  if (!strm.good()) return 0;
  //
  for (int i=kMaxPix;i--;)  recTypePix[i]  = kRecDummy;
  for (int i=kMaxMS ;i--;)  recTypeMS[i]   = kRecDummy;
  for (int i=kMaxTrig;i--;) recTypeTR[i] = kRecDummy;
  //
  double *datTmp=0;
  int ent = 0;
  //
  zMuFirst = 2e6; // Z of 1st muon observation
  
  while(1) {
    if (!recSav.IsNull()) { // header of next particle was already read
      rec = recSav;
      recSav = "";
    }
    else {
      rec = readNextRecord(strm);
      rec = rec.Strip(TString::kLeading,' ');
      //      printf("%s\n",rec.Data());
    }
    //
    if (rec.IsNull()) return 0;
    //    printf("rec is |%s|\n",rec.Data());
    if (rec.BeginsWith("Primary")) {
      if (ent>0) { // previous particle is over
	recSav = rec;
	return 1;
      }
      int nr = sscanf(rec.Data(),"Primary %d %lf %lf %lf %lf %lf %lf %lf",&code,&en,&x,&y,&z,&cx,&cy,&cz);
      if (nr!=8) {printf("Error in reading %s\n",rec.Data()); exit(1);}
      if (zMuFirst>1e6 && (code==10||code==11) ) zMuFirst = z; 
      codeOr = code;
      datTmp = recDataPrim;
    }
    else if (rec.BeginsWith("PixStn")) {
      int nr = sscanf(rec.Data(),"PixStn%d %d %lf %lf %lf %lf %lf %lf %lf",&st,&code,&en,&x,&y,&z,&cx,&cy,&cz);
      if (nr!=9) {printf("Error in reading %s\n",rec.Data()); exit(1);}
      //
      if (recTypePix[st]!=kRecDummy) {
	if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied PixStn%d |"
			    "Enew=%.2e Eold=%.2e\n",code,codeOr,st,en,recDataPix[st][kE]);
	if (en<recDataPix[st][kE]) continue; // ignore particle with lower energy
      }
      else nPix++;
      if (zMuFirst>1e6 && (code==10||code==11) ) zMuFirst = z; 
      datTmp = recDataPix[st];
      recTypePix[st] = code; // pixel station
    }
    else if (rec.BeginsWith("MS")) {
      int nr = sscanf(rec.Data(),"MS%d %d %lf %lf %lf %lf %lf %lf %lf",&st,&code,&en,&x,&y,&z,&cx,&cy,&cz);
      if (nr!=9) {printf("Error in reading %s\n",rec.Data()); exit(1);}
      //
      if (recTypeMS[st]!=kRecDummy) {
	if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied MS%d |"
			    "Enew=%.2e Eold=%.2e\n",code,codeOr,st,en,recDataMS[st][kE]);
	if (en<recDataMS[st][kE]) continue; // ignore particle with lower energy
      }
      else nMS++;
      if (zMuFirst>1e6 && (code==10||code==11) ) zMuFirst = z; 
      datTmp = recDataMS[st];
      recTypeMS[st] = code; // MS station
    }
    else if (rec.BeginsWith("TrigStn")) {
      int nr = sscanf(rec.Data(),"TrigStn%d %d %lf %lf %lf %lf %lf %lf %lf",&st,&code,&en,&x,&y,&z,&cx,&cy,&cz);
      if (nr!=9) {printf("Error in reading %s\n",rec.Data()); exit(1);}
      if (recTypeTR[st]!=kRecDummy) {
	if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied TrigSt%d |"
			    "Enew=%.2e Eold=%.2e\n",code,codeOr,st,en,recDataTR[st][kE]);
	if (en<recDataTR[st][kE]) continue; // ignore particle with lower energy	
      }
      else nTrig++;
      if (zMuFirst>1e6 && (code==10||code==11) ) zMuFirst = z; 
      datTmp = recDataTR[st];
      recTypeTR[st] = code; // MS station     
    }
    else {printf("Unknown record in %s\n",rec.Data()); exit(1);}
    //
    datTmp[kE] = en;
    datTmp[kX] = x;
    datTmp[kY] = y;
    datTmp[kZ] = z;    
    datTmp[kCX] = cx;
    datTmp[kCY] = cy;
    datTmp[kCZ] = cz;    
    ent++;
  }
  return 1;
}

char* readNextRecord(ifstream& strm)
{
  if (!strm.good()) return 0;
  static TString saveStr;
  static TString recStr,recFull;
  recFull = saveStr;
  //
  while ( recStr.ReadLine(strm) ) {
    //   printf("read %s\n",recStr.Data());
    recStr = recStr.Strip(TString::kLeading,' ');
    recStr = recStr.ReplaceAll("\r"," ");
    char fst = recStr.Data()[0];
    if (isalpha(fst)) { // new record keyword
      if (recFull.IsNull()) { // just starting
	recFull += recStr;
	recFull += " ";
	saveStr = "";
	continue;
      }
      // new record started, old one is over
      saveStr = recStr;
      return (char*)recFull.Data();
    }
    else {
      recFull += recStr;
    }
  }
  saveStr="";
  return (char*)recFull.Data();
}

void FlukaParticle(Int_t fCode, double &mass, int &charge)
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

Int_t SetInpList(const char* list)
{
  // files to parse
  ifstream inpf(list);
  if (!inpf.good()) {
    printf("Failed on input filename %s\n",list);
    return kFALSE;
  }
  inpFileList.Clear();
  nFiles = 0;
  TString flName;
  flName.ReadLine(inpf);
  while ( !flName.IsNull() ) {
    flName = flName.Strip(TString::kBoth,' ');
    if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
    flName = flName.Strip(TString::kBoth,',');
    flName = flName.Strip(TString::kBoth,'"');
    printf("Adding %s\n",flName.Data());
    inpFileList.AddLast(new TObjString(flName));
    nFiles++;
    flName.ReadLine(inpf);
  }
  return nFiles;
}

//____________________________________________________
void SetupFlukaParticle()
{
  // set input particle and fluka hits
  KMCLayerFwd* lr = 0;
  KMCProbeFwd* anProbe = 0;
  double r2,xyzTr[3];
  double rx,ry,mass,massImp=KMCDetectorFwd::kMassMu;
  int charge;
  anProbe = det->GetProbe();
  FlukaParticle(codeOr,mass,charge);
  ImposeKinematics(anProbe,&recDataPrim[kX],&recDataPrim[kCX],recDataPrim[kE],mass,charge);
  TList* lstLay = det->GetLayers();
  TIter nextL(lstLay);
  KMCLayerFwd* lrq = 0;
  while ( (lrq=(KMCLayerFwd*)nextL()) ) {
    lrq->GetMCCluster()->Kill();
  }
  //
  int maxLr = -1;
  for (int ih=0;ih<det->GetNumberOfActiveLayersITS();ih++) {
    if (recTypePix[ih]!=kRecDummy) {
      lr = det->GetLayerITS(ih);
      lr->GetMCCluster()->Kill();
      r2 = recDataPix[ih][kX]*recDataPix[ih][kX]+recDataPix[ih][kY]*recDataPix[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectITS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &recDataPix[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(), xyzTr[2]-rx*lr->GetXRes(),-1);
      anProbe = lr->GetAnProbe();
      FlukaParticle(recTypePix[ih],mass,charge);
      //      printf("\nSetLr PIX "); lr->Print();
      ImposeKinematics(anProbe,&recDataPix[ih][kX],&recDataPix[ih][kCX],recDataPix[ih][kE],mass,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  for (int ih=0;ih<det->GetNumberOfActiveLayersMS();ih++) {
    if (recTypeMS[ih]!=kRecDummy) {
      lr = det->GetLayerMS(ih);
      lr->GetMCCluster()->Kill();
      r2 = recDataMS[ih][kX]*recDataMS[ih][kX]+recDataMS[ih][kY]*recDataMS[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectMS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &recDataMS[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(), xyzTr[2]-rx*lr->GetXRes(),-1);
      FlukaParticle(recTypeMS[ih],mass,charge);
      anProbe = lr->GetAnProbe();
      //      printf("\nSetLr MS  "); lr->Print();
      ImposeKinematics(anProbe,&recDataMS[ih][kX],&recDataMS[ih][kCX],recDataMS[ih][kE],mass,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  for (int ih=0;ih<det->GetNumberOfActiveLayersTR();ih++) {
    if (recTypeTR[ih]!=kRecDummy) {
      lr = det->GetLayerTR(ih);
      lr->GetMCCluster()->Kill();
      r2 = recDataTR[ih][kX]*recDataTR[ih][kX]+recDataTR[ih][kY]*recDataTR[ih][kY];
      if (r2>lr->GetRMax()*lr->GetRMax() || r2<lr->GetRMin()*lr->GetRMin()) {
	printf("RejectMS r=%f z=%f theta=%f\n",TMath::Sqrt(r2),lr->GetZ(),TMath::ATan(r2/lr->GetZ()/lr->GetZ()));
	continue;
      }
      KMCProbeFwd::Lab2Trk( &recDataTR[ih][kX], xyzTr);
      gRandom->Rannor(rx,ry);
      lr->GetMCCluster()->Set(xyzTr[0], xyzTr[1]+ry*lr->GetYRes(), xyzTr[2]-rx*lr->GetXRes(),-1);
      FlukaParticle(recTypeTR[ih],mass,charge);
      anProbe = lr->GetAnProbe();
      //      printf("\nSetLr TR  "); lr->Print();
      ImposeKinematics(anProbe,&recDataTR[ih][kX],&recDataTR[ih][kCX],recDataTR[ih][kE],mass,charge);
      if (lr->GetID()>maxLr) maxLr = lr->GetID();
    }
  }
  det->SetLastActiveLayerTracked(maxLr);
  //
}

//____________________________________________
void ImposeKinematics(KMCProbeFwd* probe, double* xyzLab,double* cosinesLab, double en, double mass, int charge) 
{
  //
  double p = en*en - mass*mass;
  if (p<=0) {
    printf("Anomalous kinematics: E:%e M:%e",en,mass);
    exit(1);
  }
  p = TMath::Sqrt(p);
  double pxyz[3] = {p*cosinesLab[0],p*cosinesLab[1],p*cosinesLab[2]};
  //  printf("Imp : %f %f %f | %f\n",pxyz[0],pxyz[0],pxyz[0], en);
  probe->Init(xyzLab,pxyz,charge,1.e4);
  probe->SetMass(mass);
  probe->ResetCovariance();// reset cov.matrix
  //  printf("set at %f %f %f \n",xyzLab[0],xyzLab[1],xyzLab[2]);
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
