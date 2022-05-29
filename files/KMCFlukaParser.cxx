#include "KMCFlukaParser.h"
#include <TObjString.h>

const TString KMCFlukaParser::fgEndEvRecord = "**************** end of event";


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
  while(1) {
    FlukaHit hit;
    char key0[20], volint[20];
    
    rec = readNextRecord();
    rec = rec.Strip(TString::kLeading,' ');
    if (rec.IsNull()) return -9;
    if (rec.BeginsWith(fgEndEvRecord)) break; // end of record reached
    //    printf("rec is |%s|\n",rec.Data());

    int nr = sscanf(rec.Data(),"%s %d %s %lf %lf %lf %lf %lf %lf %lf", key0, &hit.partCode, volint,
		    &hit.recData[kE], &hit.recData[kX], &hit.recData[kY], &hit.recData[kZ],
		    &hit.recData[kCX], &hit.recData[kCY],&hit.recData[kCZ]); 
    if (nr!=10) {
      printf("Expected to read %d items, got %d from\n%s\n",10,nr,rec.Data());
      exit(1);
    }
    if (interactionSource.IsNull()) interactionSource = volint;
    
    if (rec.BeginsWith("PixStn")) {
      int nr = sscanf(key0,"PixStn%d",&hit.stationID);	
      if (nr!=1 || hit.stationID>=kMaxPix) {
	printf("Failed to get VerTel plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      hit.stationType = KMCLayerFwd::kITS;
    }
    else if (rec.BeginsWith("MS")) {
      int nr = sscanf(key0,"MS%d",&hit.stationID);
      if (nr!=1 || hit.stationID>=kMaxMS) {
	printf("Failed to get MS plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      hit.stationType = KMCLayerFwd::kMS;
    }
    else if (rec.BeginsWith("TrigStn")) {
      int nr = sscanf(key0,"TrigStn%d",&hit.stationID);
      if (nr!=1 || hit.stationID>=kMaxTrig) {
	printf("Failed to get Trigger plane number from:\n%s\n", rec.Data());
	exit(1);
      }
      hit.stationType = KMCLayerFwd::kTRIG;
    }
    else {
	printf("Unknown detector keyword in %s\n",rec.Data());
	exit(1);
    }
    fHits.push_back( hit );
  }
  //
  return fHits.size();
}

// for Geant version from Maryna
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
