#include "KMCFlukaParser.h"
#include <TObjString.h>

const TString KMCFlukaParser::fgEndEvRecord = "**************** end of event";

int KMCFlukaParser::Fluka2PDG(int flCode) const
{
  const int pdgCode[] =
    {
     0, // NA
     2212, // Proton   
     -2212, // Antiproton   
     11, // Electron   
     -11, // Positron   
     12, // Electron Neutrino  
     -12, // Electron Antineutrino  
     22, // Photon   
     2112, // Neutron   
     -2112, // Antineutron   
     13, // Positive Muon   
     -13, // Negative Muon   
     130, // Kaon-zero long  
     211, // Positive Pion   
     -211, // Negative Pion   
     321, // Positive Kaon   
     -321, // Negative Kaon   
     3122, // Lambda   
     -3122, // Antilambda   
     310, // Kaon zero short  
     -3112, // Negative Sigma  
     3222, // Positive Sigma  
     3212, // Sigma-zero   
     111, // Pion-zero   
     311, // Kaon-zero   
     -311, // Antikaon-zero
     0, // NA
     14, // Muon neutrino   
     -14, // Muon antineutrino
     0, // NA
     0, // NA
     -3222, // Antisigma-minus  
     -3212, // Antisigma-zero  
     -3112, // Antisigma-plus  
     3322, // Xi-zero   
     -3322, // Antixi-zero   
     -3312, // Negative Xi   
     3312, // Positive Xi   
     3334, // Omega-minus   
     -3334, // Antiomega
     0, // NA   
     15, // Positive Tau    
     -15, // Negative Tau    
     16, // Tau neutrino    
     -16, // Tau antineutrino  
     411, // D-plus   
     411, // D-minus   
     421, // D-zero   
     -421, // AntiD-zero   
     431, // D_s-plus   
     -431, // D_s-minus   
     4122, // Lambda_c-plus   
     4232, // Xi_c-plus   
     4112, // Xi_c-zero   
     4322, // Xi'_c-plus   
     4312, // Xi'_c-zero   
     4332, // Omega_c-zero    
     -4122, // Antilambda_c-minus  
     -4232, // AntiXi_c-minus  
     -4132, // AntiXi_c-zero   
     -4322, // AntiXi'_c-minus  
     -4312, // AntiXi'_c-zero  
     -4332  // AntiOmega_c-zero
    };
  return flCode<0 || flCode > int(sizeof(pdgCode)/sizeof(int)) ? 0 : pdgCode[flCode];
}
    
Int_t KMCFlukaParser::SetInpList(const char* list)
{
  // files to parse
  fInpFile.open(list);
  if (!fInpFile.good()) {
    printf("Failed on input filename %s\n",list);
    fInpFile.close();
    return kFALSE;
  }
  fInpFileList.Clear();
  TString flName;
  flName.ReadLine(fInpFile);
  while ( !flName.IsNull() ) {
    flName = flName.Strip(TString::kBoth,' ');
    if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(fInpFile); continue;}
    flName = flName.Strip(TString::kBoth,',');
    flName = flName.Strip(TString::kBoth,'"');
    printf("Adding %s\n",flName.Data());
    fInpFileList.AddLast(new TObjString(flName));
    flName.ReadLine(fInpFile);
  }
  fInpFile.close();
  return fInpFileList.GetEntriesFast();
}


bool KMCFlukaParser::GetNextGoodPair(Int_t minPix,Int_t minMS,Int_t minTr)
{
  // read next pair of fluka particles passing requirement of given number of hits in detetectors (each)
  // prepare input stream if needed
  if (fParts.size()<2) fParts.resize(2);
  
  while(1) { // loop over files
    if (!fInpFile.is_open()) {
      if (fCurFileID>=fInpFileList.GetEntriesFast()-1) return false;
      const char* fnm = fInpFileList.At(++fCurFileID)->GetName();
      printf("Processing %d-th file %s\n",fCurFileID,fnm);
      fInpFile.open(fnm);
      if (!fInpFile.good()) {
	printf("Failed to open file %s\n",fnm);
	exit(1);
      }
    }
    while (readNextPair()) {
      fStat.totalRead++;
      int nAcc = 0;
      for (int ip=2;ip--;) {
	if (fParts[ip].nPix<minPix || fParts[ip].nMS<minMS || fParts[ip].nTrig<minTr) break;
	nAcc++;
      }
      if (nAcc==2) {
	fStat.totalAccepted++;
	printf("\nRead %d Accepted %d | %d/%d %d/%d %d/%d\n",
	       fStat.totalRead,fStat.totalAccepted,
	       fParts[0].nPix,fParts[1].nPix, fParts[0].nMS,fParts[1].nMS,
	       fParts[0].nTrig, fParts[1].nTrig);
	return true;
      }
    }
    fInpFile.close();
  }
  return false;  
}


bool KMCFlukaParser::GetNextBackgroundEvent(const TString& interactionSource, bool allowRewind)
{
  // read next set of background hits, optionally from the interaction in the volume whose name starts with interactionSource
  
  while(1) { // loop over files
    if (!fInpFile.is_open()) {
      if (fCurFileID>=fInpFileList.GetEntriesFast()-1) {
	if (!allowRewind) return false;
	printf("Rewinding fluka backgroung files\n");
	fCurFileID = -1;
      }
      const char* fnm = fInpFileList.At(++fCurFileID)->GetName();
      printf("Processing %d-th file %s\n",fCurFileID,fnm);
      fInpFile.open(fnm);
      if (!fInpFile.good()) {
	printf("Failed to open file %s\n",fnm);
	exit(1);
      }
    }
    TString intVol = "";
    while (1) {
      int res = readBackground(intVol);
      fStat.totalRead++;
      if (res<-1) break; // no more data in this file
      if (!interactionSource.IsNull() && !intVol.BeginsWith(interactionSource)) continue;
      printf("\nRead background event with %d hits\n", res);
      return true;
    }
    fInpFile.close();
  }
  return false;  
}

int KMCFlukaParser::readBackground(TString& interactionSource)
{
  //
  TString rec;
  //
  fHits.clear();
  if (!fInpFile.good()) return -10;
  interactionSource = "";
  //
  double en, x, y, z, cx, cy, cz;
  int partCode;
  while(1) {
    SimHit hit;
    char vol0[20], volint[20];

    rec = readNextRecord();
    rec = rec.Strip(TString::kLeading,' ');
    if (rec.IsNull()) return -9;
    if (rec.BeginsWith(fgEndEvRecord)) break; // end of record reached
    //    printf("rec is |%s|\n",rec.Data());

    int nr = sscanf(rec.Data(), "%s %d %s %lf %lf %lf %lf %lf %lf %lf", vol0,
                    &partCode, volint, &en, &x, &y, &z, &cx, &cy, &cz);
    if (nr!=10) {
      printf("Expected to read %d items, got %d from\n%s\n",10,nr,rec.Data());
      exit(1);
    }
    if (interactionSource.IsNull()) interactionSource = volint;
    int pdg = Fluka2PDG(partCode);
    if (pdg == 0) {
      printf("Skipping particle with uknown code %d\n", partCode);
      continue;
    }
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
    if (particle->Charge() == 0)
      continue;
    int stType = -1, stID = -1;
    if (rec.BeginsWith("PixStn")) {
      int nr = sscanf(vol0, "PixStn%d", &stID);
      if (nr != 1) {
        printf("Failed to get VerTel plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      stType = KMCLayerFwd::kITS;
    } else if (rec.BeginsWith("MS")) {
      int nr = sscanf(vol0, "MS%d", &stID);
      if (nr != 1) {
        printf("Failed to get MS plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      stType = KMCLayerFwd::kMS;
    } else if (rec.BeginsWith("TrigStn")) {
      int nr = sscanf(vol0, "TrigStn%d", &stID);
      if (nr != 1) {
        printf("Failed to get Trigger plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      stType = KMCLayerFwd::kTRIG;
    } else {
      printf("Unknown detector keyword in %s\n", rec.Data());
      exit(1);
    }
    // register hit
    double ptot = TMath::Sqrt(
        en * (en + 2 * particle->Mass())); // file provides kinetic energy
    fHits.emplace_back(SimHit{stType, stID,
                              TParticle{pdg, 0, 0, -1, -1, -1, -1, ptot * cx,
                                        ptot * cy, ptot * cz, x, y, z, 0.}});
  }
  //
  return fHits.size();
}

// for Geant version from Maryna
bool KMCFlukaParser::readNextTrackGeant()
{
  //
  TString rec;
  double en,x,y,z,cx,cy,cz;
  int code, trID;
  char vol0[20], vol1[20];
  int signalTrackID = -1;
  while(1) {
    if (!fInpFile.is_open()) {
      if (fCurFileID>=fInpFileList.GetEntriesFast()-1) return false;
      const char* fnm = fInpFileList.At(++fCurFileID)->GetName();
      printf("Processing %d-th file %s\n",fCurFileID,fnm);
      fInpFile.open(fnm);
      if (!fInpFile.good()) {
	printf("Failed to open file %s\n",fnm);
	exit(1);
      }
    }
    
    rec = readNextRecord();
    rec = rec.Strip(TString::kLeading,' ');
    if (rec.IsNull()) {
      fInpFile.close();
      return false;
    }
    if (rec.BeginsWith(fgEndEvRecord)) break; // end of record reached
    //    printf("rec is |%s|\n",rec.Data());
    if (rec.BeginsWith("#")) continue; // comment
    
    int nr = sscanf(rec.Data(),"%s %d %s %lf %lf %lf %lf %lf %lf %lf %d",vol0,&code,vol1,&en,&x,&y,&z,&cx,&cy,&cz,&trID);
    if (nr!=11) {
      printf("Expected to read %d items, got %d from\n%s\n",11,nr,rec.Data());
      exit(1);
    }
    if (trID == 0) trID = 1; // RS: this is a hack, in the files of Maryna 1 means injected particle
    int pdg = Fluka2PDG(code);
    if (pdg==0) {
      printf("Skipping particle with uknown code %d\n", code);
      continue;
    }
    // tmp rotate by pi/2
    {
      double tmp;
      tmp = y; y = x; x = -tmp;
      tmp = cy; cy = cx; cx = -tmp;
    }
    //
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
    if (particle->Charge()==0) continue;
    double ptot = TMath::Sqrt(en*(en + 2*particle->Mass())); // file provides kinetic energy
    double etot = en + particle->Mass();
    if (!strcmp(vol0,"IP")) { // interaction point: this is a primary
      // register primary kinematics
      // register particle
      simEvent.signal.emplace_back(pdg, 0, -1, -1, -1, -1, ptot*cx, ptot*cy, ptot*cz, etot, x, y, z, 0.);
      simEvent.signalHits.emplace_back(); // add empty vector of hits for this particle
      signalTrackID = trID;
    } else { // we are seeing hits
      int stType = -1, stID = -1;
      if (rec.BeginsWith("PixStn")) {
	int nr = sscanf(vol0,"PixStn%d",&stID);
	if (nr!=1) {
	  printf("Failed to get VerTel plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	stType = KMCLayerFwd::kITS;
      } else if (rec.BeginsWith("MS")) {
	int nr = sscanf(vol0,"MS%d",&stID);
	if (nr!=1) {
	  printf("Failed to get MS plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	stType = KMCLayerFwd::kMS;
      } else if (rec.BeginsWith("TrigStn")) {
	int nr = sscanf(vol0,"TrigStn%d",&stID);
	if (nr!=1) {
	  printf("Failed to get Trigger plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	stType = KMCLayerFwd::kTRIG;
      } else {
	printf("Unknown detector keyword in %s\n",rec.Data());
	exit(1);
      }
      // register hit
      auto& dest = (trID == signalTrackID) ? simEvent.signalHits.back() : simEvent.bgHits;
      dest.emplace_back(SimHit{stType, stID, TParticle{pdg, 0, -1, -1, -1, -1, ptot*cx, ptot*cy, ptot*cz, etot, x, y, z, 0.}});
      dest.back().track.SetUniqueID(trID);      
    }
  }
  return true;
}

/*
// for old fluka version from Gianluca
int KMCFlukaParser::readNextPair(Bool_t verbose)
{
  //
  TString rec;
  double en,x,y,z,cx,cy,cz;
  int code,codeP,st;
  char key0[20];
  //
  if (!fInpFile.good()) return 0;
  //
  for (int ip=2;ip--;) fParts[ip].clear();
  //
  int nPart = 0, idPart = 0;
  double *datTmp=0;
  //  
  while(1) {
    rec = readNextRecord();
    rec = rec.Strip(TString::kLeading,' ');
    if (rec.IsNull()) return nPart;
    if (rec.BeginsWith(fgEndEvRecord)) return nPart; // end of record reached
    //    printf("rec is |%s|\n",rec.Data());

    int nr = sscanf(rec.Data(),"%s %d %d %lf %lf %lf %lf %lf %lf %lf",key0,&code,&codeP,&en,&x,&y,&z,&cx,&cy,&cz);
    if (nr!=10) {
      printf("Expected to read %d items, got %d from\n%s\n",10,nr,rec.Data());
      exit(1);
    }
    if (!strcmp(key0,"Primary")) {
      nPart++;
      if (nPart>2) {
	printf("More than 2 Primary records are found in the same event\n%s\n",rec.Data());
	exit(1);
      }
      idPart = nPart-1;
      fParts[idPart].codeOr = code;
      datTmp = fParts[idPart].recDataPrim;
    } 
    else {
      // make sure 2 primaries were read
      if (nPart<2) {
	printf("Non-Primary record is found while only %d primary is read (must be 2)\n%s\n",nPart,rec.Data());
	exit(1);
      }
      // find to which particle this record belongs to
      if (codeP==fParts[0].codeOr) {
	idPart = 0;
      }
      else if (codeP==fParts[1].codeOr) {
	idPart = 1;
      }
      else {
	printf("skip record: detector hit parent code %d does not correspond to neither of parents %d and %d\n%s\n",
	       codeP,fParts[0].codeOr,fParts[1].codeOr,rec.Data());
	//exit(1);
	continue;
      }
      if (rec.BeginsWith("PixStn")) {
	int nr = sscanf(key0,"PixStn%d",&st);
	if (nr!=1 || st>=kMaxPix) {
	  printf("Failed to get VerTel plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	if (fParts[idPart].recTypePix[st]!=kRecDummy) {
	  if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied PixStn%d |"
			      "Enew=%.2e Eold=%.2e\n",code,codeP,st,en,fParts[idPart].recDataPix[st][kE]);
	  if (en<fParts[idPart].recDataPix[st][kE]) continue; // ignore particle with lower energy
	}
	else fParts[idPart].nPix++;
	
	datTmp = fParts[idPart].recDataPix[st];
	fParts[idPart].recTypePix[st] = code; // pixel station
      }
      else if (rec.BeginsWith("MS")) {
	int nr = sscanf(key0,"MS%d",&st);
	if (nr!=1 || st>=kMaxMS) {
	  printf("Failed to get MS plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	if (fParts[idPart].recTypeMS[st]!=kRecDummy) {
	  if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied MS%d |"
			      "Enew=%.2e Eold=%.2e\n",code,codeP,st,en,fParts[idPart].recDataMS[st][kE]);
	  if (en<fParts[idPart].recDataMS[st][kE]) continue; // ignore particle with lower energy
	}
	else fParts[idPart].nMS++;

	datTmp = fParts[idPart].recDataMS[st];
	fParts[idPart].recTypeMS[st] = code; // MS station
	//
      }
      else if (rec.BeginsWith("TrigStn")) {
	int nr = sscanf(key0,"TrigStn%d",&st);
	if (nr!=1 || st>=kMaxTrig) {
	  printf("Failed to get Trigger plane number from:\n%s\n", rec.Data());
	  exit(1);
	}
	if (fParts[idPart].recTypeTR[st]!=kRecDummy) {
	  if (verbose) printf("Competing hit from part. %d (primary:%d) on occupied MS%d |"
			      "Enew=%.2e Eold=%.2e\n",code,codeP,st,en,fParts[idPart].recDataTR[st][kE]);
	  if (en<fParts[idPart].recDataTR[st][kE]) continue; // ignore particle with lower energy
	}
	else fParts[idPart].nTrig++;

	datTmp = fParts[idPart].recDataTR[st];
	fParts[idPart].recTypeTR[st] = code; // Trigger station
	//
      }
      else {
	printf("Unknown detector keyword in %s\n",rec.Data());
	exit(1);
      }
    }
    //
    if (fParts[idPart].zMuFirst>1e6 && (code==10||code==11) ) fParts[idPart].zMuFirst = z;
    datTmp[kE] = en;
    datTmp[kX] = x;
    datTmp[kY] = y;
    datTmp[kZ] = z;    
    datTmp[kCX] = cx;
    datTmp[kCY] = cy;
    datTmp[kCZ] = cz;    
    //      
  }
  return 0;
}
*/

char* KMCFlukaParser::readNextRecord()
{
  if (!fInpFile.good()) return 0;
  static TString recStr;
  //
  while ( recStr.ReadLine(fInpFile) ) {
    //   printf("read %s\n",recStr.Data());
    recStr = recStr.Strip(TString::kLeading,' ');
    recStr = recStr.ReplaceAll("\r"," ");
    if (recStr.BeginsWith("#")) continue; // commented out
    return (char*)recStr.Data();
  }
  return 0;
}
