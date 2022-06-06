#ifndef KMCFLUKAPARSER_H
#define KMCFLUKAPARSER_H

#include <TString.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <fstream>
#include <vector>
#include "KMCLayerFwd.h"
//====================================================================
//
// These are settings for the fluka parser
//
enum {kRecDummy=-9999,kMaxPix=20,kMaxMS=20,kMaxTrig=10};
enum {kE,kX,kY,kZ,kCX,kCY,kCZ,kNFld};

struct SimHit {
  int stationType = int(KMCLayerFwd::kTypeNA);
  int stationID = -1;
  TParticle track;
  ClassDefNV(SimHit,1);
};

struct SimEvent {
  std::vector<TParticle> signal;
  std::vector<std::vector<SimHit>> signalHits;
  std::vector<SimHit> bgHits;
  void clear() {
    signal.clear();
    signalHits.clear();
    bgHits.clear();
  }
  ClassDefNV(SimEvent, 1);
};

struct FlukaPart {
  int    codeOr; // particle fluka code
  double massOr;
  int    chargeOr;
  double recDataPrim[kNFld];
  int    recTypePix[kMaxPix];
  double recDataPix[kMaxPix][kNFld];
  int    recTypeMS[kMaxMS];
  double recDataMS[kMaxMS][kNFld];
  int    recTypeTR[kMaxTrig];
  double recDataTR[kMaxTrig][kNFld];
  int nPix,nMS,nTrig;
  float zMuFirst; // Z of 1st muon observation

  FlukaPart() : codeOr(0),massOr(0.),chargeOr(0.),nPix(0),nMS(0),nTrig(0),zMuFirst(2e6) {}
  void clear() {
    nPix = nMS = nTrig = 0;
    for (int i=kMaxPix;i--;)  recTypePix[i]  = kRecDummy;
    for (int i=kMaxMS ;i--;)  recTypeMS[i]   = kRecDummy;
    for (int i=kMaxTrig;i--;) recTypeTR[i] = kRecDummy;
    zMuFirst = 2e6;
  }
  ClassDefNV(FlukaPart,1);
};

struct FlukaStat {
  int totalRead;
  int totalAccepted;
  FlukaStat() : totalRead(0), totalAccepted(0) {}
  ClassDefNV(FlukaStat,1);
};

class KMCFlukaParser {

 public:

  KMCFlukaParser() : fInpFileList(),fCurFileID(-1),fInpFile(),fStat() {}
  ~KMCFlukaParser() {fInpFileList.Delete();}
    
  Int_t SetInpList(const char* list);
  bool GetNextGoodPair(Int_t minPix=0,Int_t minMS=3,Int_t minTr=2);

  bool GetNextBackgroundEvent(const TString& interactionSource="", bool allowRewind = true);
  int readBackground(TString& interactionSource);

  bool readNextTrackGeant();
  
  const FlukaStat& GetStat() const {return fStat;}
  const std::vector<FlukaPart>& GetParticles() const {return fParts;}
  const std::vector<SimHit> &GetHits() const { return fHits; }

  int Fluka2PDG(int flcode) const;

  const SimEvent& getSimEvent() const { return simEvent; }
  SimEvent& getSimEvent() { return simEvent; }
  
 protected:

  int readNextPair(Bool_t verbose=kFALSE);
  char* readNextRecord();

  SimEvent simEvent;
  
  TObjArray fInpFileList;
  Int_t     fCurFileID;
  std::ifstream  fInpFile;
  std::vector<FlukaPart> fParts;
  std::vector<SimHit> fHits;
  FlukaStat fStat;
  static const TString fgEndEvRecord;

  ClassDefNV(KMCFlukaParser,0)
};

#endif
