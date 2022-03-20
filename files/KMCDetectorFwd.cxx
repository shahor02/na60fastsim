#include <TStopwatch.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include "KMCDetectorFwd.h"
#include "KMCMagnetBuilder.h"

ClassImp(KMCDetectorFwd)

const Double_t KMCDetectorFwd::kMassP  = 0.938;
const Double_t KMCDetectorFwd::kMassK  = 0.4937;
const Double_t KMCDetectorFwd::kMassPi = 0.1396;
const Double_t KMCDetectorFwd::kMassMu = 0.1057;
const Double_t KMCDetectorFwd::kMassE  = 0.0005;

const Double_t KMCDetectorFwd::fgkFldEps = 1.e-4;

KMCDetectorFwd::KMCDetectorFwd(const char *name, const char *title) 
:  TNamed(name,title)
  ,fNLayers(0)
  ,fNActiveLayers(0)
  ,fNActiveLayersITS(0)
  ,fNActiveLayersMS(0)
  ,fNActiveLayersTR(0)
  ,fLastActiveLayerITS(0)
  ,fLastActiveLayerMS(0)
  ,fLastActiveLayer(0)
  ,fLastActiveLayerTracked(0)
  ,fIntegrateMSX2X0(kFALSE)
  ,fMSLrMinID(-1),
   fMSLrMaxID(-1)
  ,fLayers()
  ,fBeamPipe(0)
  ,fVtx(0)
  ,fMaterials()
  ,fMagFieldID(-1)
  ,fProbe()
  ,fExternalInput(kFALSE)
  ,fIncludeVertex(kTRUE)
  ,fApplyBransonPCorrection(0.1)
  ,fUseBackground(kTRUE)
  ,fMaxChi2Cl(10.)
  ,fMaxNormChi2NDF(10)
  ,fMaxChi2Vtx(20)
  ,fMinITSHits(0)
  ,fMinMSHits(0)
  ,fMinTRHits(0)
  ,fMinP2Propagate(1.)
  ,fMaxChi2ClSQ(0)
  ,fMaxSeedToPropagate(0)
   //
,fDecayZProf(0)
,fZDecay(0)
,fDecMode(kNoDecay)
,fChi2MuVtx(0)
,fZBendingMS(0)
,fZToroidStart(0)
,fZToroidEnd(0)
,fToroidB0(0)
,fFldNReg(0)
,fFldZMins(0)
,fFldZMaxs(0)
,fDefStepAir(1.0)
,fDefStepMat(1.0)
,fImposeVertexPosition(kFALSE)
,fPattITS(0)
,fNCh(-1)
,fdNdY(0)
,fdNdPt(0)
,fHChi2Branson(0)
,fHChi2LrCorr(0)
,fHChi2NDFCorr(0)
,fHChi2NDFFake(0)
,fHChi2VtxCorr(0)
,fHChi2VtxFake(0)
,fHNCand(0)
,fHChi2MS(0)
,fHCandCorID(0)
,fZSteps(0)
{
  //
  // default constructor
  //
  //  fLayers = new TObjArray();
  fRefVtx[0] = fRefVtx[1] = fRefVtx[2] = 0;
  //
  for (int i=0;i<7;i++) {
    fUseRPhiErr[i] = false;
  }
}

KMCDetectorFwd::~KMCDetectorFwd() { // 
  // virtual destructor
  //
  //  delete fLayers;
  delete fdNdY;
  delete fdNdPt;
}


void KMCDetectorFwd::Print(const Option_t *opt) const
{
  // Prints the detector layout
  printf("KMCDetectorFwd %s: \"%s\"\n",GetName(),GetTitle());
  for (Int_t i = 0; i<fLayers.GetEntries(); i++) GetLayer(i)->Print(opt);
  if (fBeamPipe) fBeamPipe->Print(opt);
}

Double_t KMCDetectorFwd::ThetaMCS ( Double_t mass, Double_t radLength, Double_t momentum ) const
{
  //
  // returns the Multiple Couloumb scattering angle (compare PDG boolet, 2010, equ. 27.14)
  //
  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( radLength > 0 ) theta  =  0.0136 * TMath::Sqrt(radLength) / ( beta * momentum ) * (1+0.038*TMath::Log(radLength)) ;
  return (theta) ;
}

//__________________________________________________________________________
void KMCDetectorFwd::ReadMaterials(const char* fnam)
{
  // Read materials description from the external file
  //
  const char formMat[]="afffff?ff|";
  const char keyMat[]="material";
  const char modMat[]="";
  const char formMix[]="dafffff?ff|";
  const char keyMix[]="mixture";
  const char modMix[]="";
  //
  const Int_t kMinArgMat=6;
  const Int_t kMinArgMix=7;
  //
  NaCardsInput* inp;
  char *name;
  Float_t arg[3];
  Int_t narg;
  inp = new NaCardsInput();
  //
  if( !(inp->OpenFile(fnam)) ) {
    delete inp;
    Error("ReadMaterials","Did not find File %s",fnam);
    return;
  }
  //
  TObjArray* arr = &fMaterials;
  //
  // loop over all materials
  while( (narg = inp->FindEntry(keyMat,modMat,formMat,0,1))>-1 ) {
    name = inp->GetArg(0,"U"); // convert material name to upper case
    if (arr->FindObject(name)) {
      Warning("ReadMaterials","%s is already defined",name); continue;
    }
    Float_t* elbuf = 0;
    if ( narg>kMinArgMat ) {
      for (int i=0;i<NaMaterial::kNELossPar;i++) arg[i]=inp->GetArgF(kMinArgMat+i); // ELoss supplied
      elbuf = arg;
    }
    arr->AddLast(new NaMaterial(name,name,inp->GetArgF(1),inp->GetArgF(2),
				inp->GetArgF(3),inp->GetArgF(4),inp->GetArgF(5),elbuf));
  }
  //
  inp->Rewind(); // rewind the file
  // loop over all mixtures
  while( (narg = inp->FindEntry(keyMix,modMix,formMix,0,1))>-1 ) {
    name = inp->GetArg(1,"U"); // convert material name to upper case
    if ( arr->FindObject(name) ) {
      Warning("ReadMaterials","%s is already defined",name); continue;
    }
    Float_t* elbuf = 0;
    if ( narg>kMinArgMix ) {
      for (int i=0;i<3;i++) arg[i]=inp->GetArgF(kMinArgMix+i); // ELoss supplied
      elbuf = arg;
    }
    NaMixture* mix = new NaMixture(name,name,inp->GetArgF(2),inp->GetArgF(3),
				   inp->GetArgF(4),inp->GetArgF(5),
				   inp->GetArgF(6),elbuf);
    Int_t nmix = inp->GetArgD(0);
    Int_t nmixa = TMath::Abs(nmix);
    char* form = new char[nmixa+2]; // create format for mixture components
    for (int i=0;i<nmixa;i++) form[i]='f';
    form[nmixa] = '|'; form[nmixa+1]='\0';
    // Now read A,Z and W arrays
    // A array
    narg = inp->NextEntry();
    if ( inp->CompareArgList(form) ) {
      Error("ReadMaterials","mixture:A %d values are expected",nmixa);
      inp->Print(); delete[] form; continue;
    }
    Float_t* arrA = new Float_t[narg];
    for (int i=0;i<narg;i++) arrA[i] = inp->GetArgF(i);
    //
    // Z array
    narg = inp->NextEntry();
    if ( inp->CompareArgList(form) ) {
      Error("ReadMaterials","mixture:Z %d values are expected",nmixa);
      inp->Print(); delete[] form; delete[] arrA; continue;
    }
    Float_t* arrZ = new Float_t[narg];
    for (int i=0;i<narg;i++) arrZ[i] = inp->GetArgF(i);
    //
    // W array
    narg = inp->NextEntry();
    if ( inp->CompareArgList(form) ) {
      Error("ReadMaterials","mixture:W %d values are expected",nmixa);
      inp->Print(); delete[] form; delete[] arrA; delete[] arrZ; continue;
    }
    Float_t* arrW = new Float_t[narg];
    for (int i=0;i<narg;i++) arrW[i] = inp->GetArgF(i);
    //
    mix->SetComponents(nmix,arrA,arrZ,arrW);
    delete[] arrA; delete[] arrZ; delete[] arrW; delete[] form;
    arr->AddLast(mix);
    //
  }
  delete inp;
  //
}

//__________________________________________________________________________
void KMCDetectorFwd::ReadSetup(const char* setup, const char* materials)
{
  // load setup from external file
  ReadMaterials(materials);
  NaMaterial* mat = 0;
  int narg;
  //
  NaCardsInput *inp = new NaCardsInput();
  if( !(inp->OpenFile(setup)) ) {Error("BuilSetup","Did not find setup File %s",setup);exit(1);}
  //
  // -------------------------------------------------------------------
  // modified (adf 07/02/2019) to read magnets geometry and field from setup, 
  // with optional data  
  //  original line commented here; new lines follow 
  //  if ( (narg=inp->FindEntry("define","magfield","d|",1,1))>0 ) fMagFieldID = inp->GetArgD(0);
  if ( (narg=inp->FindEntry("define","magfield","d?f?f?f?f?f?f?f?f",1,1))>0 ) fMagFieldID = inp->GetArgD(0);
  Double_t zminDipole  = narg > 1 ? inp->GetArgF(1) : -9999;  
  Double_t zmaxDipole  = narg > 2 ? inp->GetArgF(2) : -9999;  
  Double_t dipoleField = narg > 3 ? inp->GetArgF(3) : -9999;  
  Double_t zminToroid  = narg > 4 ? inp->GetArgF(4) : -9999;  
  Double_t zmaxToroid  = narg > 5 ? inp->GetArgF(5) : -9999;  
  Double_t toroidField = narg > 6 ? inp->GetArgF(6) : -9999;  
  Double_t toroidRmin  = narg > 7 ? inp->GetArgF(7) : -9999;  
  Double_t toroidRmax  = narg > 8 ? inp->GetArgF(8) : -9999;  
  // end modification
  // -------------------------------------------------------------------
  //====================================================================
  if ( (narg=inp->FindEntry("define","toroidBody","dffff",1,1))>0 ) {
    if (narg!=5) {
      printf("toroidBody requested with %d arguments, need 5: type,zStart, Length, radL, dens\n",narg);
      exit(1);
    }
    int torType = inp->GetArgD(0);
    if (torType == 1) {
      new KMCMagnetBuilder(this, inp->GetArgF(1), inp->GetArgF(2), inp->GetArgF(3), inp->GetArgF(4));
    }
    else {
      printf("Unknown toroidBody type %d, alowed: 1\n",torType);
      exit(1);
    }
  }

  //====================================================================
  // beampipe
  if ( (narg=inp->FindEntry("beampipe","","ffa|",1,1))<0 ) printf("No BeamPipe found in the setup %s\n",setup);
  else {
    mat = GetMaterial(inp->GetArg(2,"U"));
    if (!mat) {printf("Material %s is not defined\n",inp->GetArg(2,"U")); exit(1);}
    AddBeamPipe(inp->GetArgF(0), inp->GetArgF(1), mat->GetRadLength(), mat->GetDensity(), mat );
  }
  //
  if ( (narg=inp->FindEntry("vertex","","ffff|",1,1))<0 ) printf("No vertex found in the setup %s\n",setup);
  else {
    fVtx = AddLayer("vtx","vertex",inp->GetArgF(0),0,0,inp->GetArgF(1), inp->GetArgF(2),inp->GetArgF(3));
  }
  //
  //
  // dummy material (not the official absorber
  inp->Rewind();
  while ( (narg=inp->FindEntry("dummy","","aaff|",0,1))>0 ) {
    mat = GetMaterial(inp->GetArg(1,"U"));
    if (!mat) {printf("Material %s is not defined\n",inp->GetArg(1,"U")); exit(1);}
    KMCLayerFwd* lr = AddLayer("dummy", inp->GetArg(0,"U"),  inp->GetArgF(2), mat->GetRadLength(), mat->GetDensity(), inp->GetArgF(3));
    lr->SetMaterial(mat);
    lr->SetDead(kTRUE);
  }
  //
  // Absorber
  inp->Rewind();
  while ( (narg=inp->FindEntry("absorber","","aaff|",0,1))>0 ) {
    mat = GetMaterial(inp->GetArg(1,"U"));
    if (!mat) {printf("Material %s is not defined\n",inp->GetArg(1,"U")); exit(1);}
    KMCLayerFwd* lr = AddLayer("abs", inp->GetArg(0,"U"),  inp->GetArgF(2), mat->GetRadLength(), mat->GetDensity(), inp->GetArgF(3));
    lr->SetMaterial(mat);
    lr->SetDead(kTRUE);
  }
  //
  // active layers
  //
  inp->Rewind();
  TString fmtAct = "aaffff*fff";
  for (int i=0;i<KMCLayerFwd::kMaxAccReg-1;i++) fmtAct += "fff";
  while ( (narg=inp->FindEntry("activelayer",0,fmtAct.Data(),0,1))>0 ) {
    // expect format "activelayer:type	NAME MATERIAL Z	DZ sigmaX sigmaZ [eff RMin RMax] [sigmaX1 sigmaY1 RMax1] ...  [sigmaX4 sigmaYN RMaxN]
    // the optional triplets [sigmaXk sigmaYk RMaxk] allow to define regions with different resolutions and eff, max possible N is 
    // KMCLayerFwd::kMaxAccReg-1
    mat = GetMaterial(inp->GetArg(1,"U"));
    double eff = narg > 6 ? inp->GetArgF(6) : 1.0;
    double rmn = narg > 7 ? inp->GetArgF(7) : 0.0;
    double rmx = narg > 8 ? inp->GetArgF(8) : 1e9;    
    if (!mat) {printf("Material %s is not defined\n",inp->GetArg(1,"U")); exit(1);}
    KMCLayerFwd* lr = AddLayer(inp->GetModifier(), inp->GetArg(0,"U"),  inp->GetArgF(2),  
			       mat->GetRadLength(), mat->GetDensity(), 
			       inp->GetArgF(3), inp->GetArgF(4), inp->GetArgF(5), eff);    
    lr->SetRMin(rmn);
    lr->SetRMax(rmx);
    lr->SetMaterial(mat);
    int nExtra = narg - 9; // are there extra settings
    if (nExtra>0) {
      if ( (nExtra%3) ) {
	printf("ReadSetup: %d extra values provided for activelayer are not multiple of 3\n",nExtra);
	printf("%s\n",inp->GetLastBuffer());
	exit(1);
      }
      int nBloc = nExtra/3;
      if (nBloc>=KMCLayerFwd::kMaxAccReg) {
	printf("ReadSetup: number of extra regions in activelayer %d should not exceed %d\n",nBloc,KMCLayerFwd::kMaxAccReg-1);
	printf("%s\n",inp->GetLastBuffer());
	exit(1);
      }
      lr->SetNAccRegions(nBloc+1);
      for (int i=0;i<nBloc;i++) {
	double sigxE = inp->GetArgF(9+3*i+0);
	double sigyE = inp->GetArgF(9+3*i+1);	
	double rmaxE = inp->GetArgF(9+3*i+2);
	if (rmaxE <= lr->GetRMax(i)) {
	  printf("ReadSetup: RMax=%.3f of %d-th extra region in activelayer must exceed R=%.3f of previous region\n",rmaxE, i+1, lr->GetRMax(i));
	  printf("%s\n",inp->GetLastBuffer());
	  exit(1);
	}
	lr->SetRMin(lr->GetRMax(i), i+1);
	lr->SetRMax(rmaxE, i+1);
	lr->SetXRes(sigxE, i+1);
	lr->SetYRes(sigyE, i+1);
      }
      
    }
    
  }
  Add3x1DMSStations(inp);
  //-------------------------------------
  //
  // init mag field
  if (TGeoGlobalMagField::Instance()->GetField()) printf("Magnetic Field is already initialized\n");
  else {
    TVirtualMagField* fld = 0;
    // -------------------------------------------------------------------
    // modified (adf 07/02/2019) to read magnets geometry and field from setup, 
    // with optional data  
    fld = new MagField(TMath::Abs(fMagFieldID));
    if (zminDipole>-9999) ((MagField *) fld)->SetZMin(0,zminDipole);
    if (zmaxDipole>-9999) ((MagField *) fld)->SetZMax(0,zmaxDipole);
    if (zminToroid>-9999) ((MagField *) fld)->SetZMin(1,zminToroid);
    if (zmaxToroid>-9999) ((MagField *) fld)->SetZMax(1,zmaxToroid);
    if (dipoleField>-9999) ((MagField *) fld)->SetBVals(0,0,dipoleField);
    if (toroidField>-9999) {
      ((MagField *) fld)->SetBVals(1,0,toroidField);
      fToroidB0 = toroidField;
    }
    if (toroidRmin>-9999) ((MagField *) fld)->SetBVals(1,1,toroidRmin);
    if (toroidRmax>-9999) ((MagField *) fld)->SetBVals(1,2,toroidRmax);
    if (zminToroid>-9999 && zmaxToroid>-9999) {
      fZBendingMS = 0.5*(zminToroid + zmaxToroid);
      fZToroidStart = zminToroid;
      fZToroidEnd = zmaxToroid;
    }
    // end modification
    // -------------------------------------------------------------------
    TGeoGlobalMagField::Instance()->SetField( fld );
    TGeoGlobalMagField::Instance()->Lock();
  }
  //
  ClassifyLayers();
  //  BookControlHistos();
  //
}

//__________________________________________________________________________
void KMCDetectorFwd::Add3x1DMSStations(NaCardsInput *inp)
{
  //
  inp->Rewind();
  TString fmtAct = "aafffddf*fffffff*fffffff*fffffff*fffffff*fffffff*fffffff";
  for (int i=0;i<KMCLayerFwd::kMaxAccReg-1;i++) fmtAct += "fff";
  int narg = 0;
  while ( (narg=inp->FindEntry("active3x1D",0,fmtAct.Data(),0,1))>0 ) {
    // expect format "active3x1D:type	NAME MATERIAL Z	DZ eff NSectors NSegments RMin NSegments*[RMax sigmaR sigmaRPhi phiUV pitchUV phiW pitchW]
    // the sector description is: RMax of radial segment, its sigmaRPhi and sigmaR, phiUV of strip planes of pitchUV, phiW pitchW of wire chamber.
    // maximum KMCLayerFwd::kMaxAccReg-1 segments per sector is allowed, at least 1 should be provided 
    NaMaterial* mat = GetMaterial(inp->GetArg(1,"U"));
    int nSectors = inp->GetArgD(5);
    int nSegments = inp->GetArgD(6);
    if (nSegments<1) {
      printf("Add3x1DMSStations: at least 1 segment per sector must be defined\n");
      printf("%s\n",inp->GetLastBuffer());
      exit(1);
    }
    if (narg!=8+nSegments*7) {
      printf("Add3x1DMSStations: %d segments per sector require %d parameters\n", nSegments, 7+nSegments*8);
      printf("%s\n",inp->GetLastBuffer());
      exit(1);
    }
    double thickness = inp->GetArgF(3), eff = inp->GetArgF(4), radL = mat->GetRadLength(), density = mat->GetDensity();
    std::vector<float> rsegm, phiUV, pitchUV, phiW, pitchW, sigR, sigRPhi;
    float rMin=inp->GetArgF(7), rMax = 0.;
    rsegm.push_back(rMin);
    for (int is=0;is<nSegments;is++) {
      rMax = inp->GetArgF(8+is*7+0);
      rsegm.push_back(rMax);
      sigRPhi.push_back( inp->GetArgF(8+is*7+1) );
      sigR.push_back( inp->GetArgF(8+is*7+2) );
      phiUV.push_back( inp->GetArgF(8+is*7+3) );
      pitchUV.push_back( inp->GetArgF(8+is*7+4) );
      phiW.push_back( inp->GetArgF(8+is*7+5) );
      pitchW.push_back( inp->GetArgF(8+is*7+6) );    
    }

    KMCMSStation* lr = new KMCMSStation(inp->GetArg(0,"U"));
    lr->SetZ(inp->GetArgF(2));
    lr->SetThickness(thickness);
    lr->SetX2X0( radL>0 ? thickness*density/radL : 0);
    lr->SetXTimesRho(thickness*density);
    lr->SetXRes(sigRPhi[0]);
    lr->SetYRes(sigR[0]);
    lr->SetLayerEff(eff);
    lr->SetRMin(rMin);
    lr->SetRMax(rMax);
    lr->SetMaterial(mat);
    lr->init(nSectors, rsegm, phiUV, pitchUV, phiW, pitchW, sigR, sigRPhi);
    lr->SetRPhiError(true);
    AddLayer(lr, inp->GetModifier());
  }
}

//__________________________________________________________________________
KMCLayerFwd* KMCDetectorFwd::AddLayer(const char* type, const char *name, Float_t zPos, Float_t radL, Float_t density, 
				      Float_t thickness, Float_t xRes, Float_t yRes, Float_t eff, NaMaterial* mat) 
{
  //
  // Add additional layer to the list of layers (ordered by z position)
  // 
  KMCLayerFwd *newLayer = GetLayer(name);
  //
  if (!newLayer) {
    TString types = type;
    types.ToLower();
    newLayer = new KMCLayerFwd(name);
    newLayer->SetZ(zPos);
    newLayer->SetThickness(thickness);
    newLayer->SetX2X0( radL>0 ? thickness*density/radL : 0);
    newLayer->SetXTimesRho(thickness*density);
    newLayer->SetXRes(xRes);
    newLayer->SetYRes(yRes);
    newLayer->SetLayerEff(eff);
    if      (types=="vt")   newLayer->SetType(KMCLayerFwd::kITS);
    else if (types=="ms")   newLayer->SetType(KMCLayerFwd::kMS);
    else if (types=="tr")   newLayer->SetType(KMCLayerFwd::kTRIG);
    else if (types=="vtx")  {newLayer->SetType(KMCLayerFwd::kVTX); }
    else if (types=="abs")  {newLayer->SetType(KMCLayerFwd::kABS); newLayer->SetDead(kTRUE); }
    else if (types=="dummy")  {newLayer->SetType(KMCLayerFwd::kDUMMY); newLayer->SetDead(kTRUE); }
    if (newLayer->GetType()>=0) {
      newLayer->SetRPhiError( fUseRPhiErr[ newLayer->GetType() ] );
    }
    //
    if (!newLayer->IsDead()) newLayer->SetDead( xRes>=kVeryLarge && yRes>=kVeryLarge);
    //
    if (fLayers.GetEntries()==0) fLayers.Add(newLayer);
    else {      
      for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
	KMCLayerFwd *l = GetLayer(i);
	if (zPos<l->GetZ()) { fLayers.AddBefore(l,newLayer); break; }
	if (zPos>l->GetZ() && (i+1)==fLayers.GetEntries() ) { fLayers.Add(newLayer); } // even bigger then last one
      }      
    }
    //
  } else printf("Layer with the name %s does already exist\n",name);
  newLayer->SetMaterial(mat);
  //
  return newLayer;
}

//__________________________________________________________________________
void KMCDetectorFwd::AddLayer(KMCLayerFwd* newLayer, const char* type)
{
  //
  // Add additional layer to the list of layers (ordered by z position)
  //
  KMCLayerFwd *old = GetLayer(newLayer->GetName());
  if (old) {
    printf("Error: Layer with the name %s does already exist\n",old->GetName());
    exit(1);
  }
  TString types = type;
  types.ToLower();
  if (types=="ms")   newLayer->SetType(KMCLayerFwd::kMS);
  else if (types=="tr")   newLayer->SetType(KMCLayerFwd::kTRIG);
  else {
    printf("Error: 3x1D layer is allowed only for MS and TR layers\n");
    exit(1);
  }
  //
  if (fLayers.GetEntries()==0) fLayers.Add(newLayer);
  else {      
    for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
      KMCLayerFwd *l = GetLayer(i);
      if (newLayer->GetZ()<l->GetZ()) { fLayers.AddBefore(l,newLayer); break; }
      if (newLayer->GetZ()>l->GetZ() && (i+1)==fLayers.GetEntries() ) { fLayers.Add(newLayer); } // even bigger then last one
    }      
  }
}

//__________________________________________________________________________
void KMCDetectorFwd::SetUseRPhiError(bool v, int lrType)
{
  fUseRPhiErr[ lrType ] = v;
  for (int i=0;i<fLayers.GetEntries(); i++) {
    KMCLayerFwd *l = GetLayer(i);
    if (l->GetType() == lrType) l->SetRPhiError(v);
  }
}

//__________________________________________________________________________
KMCPolyLayer* KMCDetectorFwd::AddPolyLayer(const char* type, const char *name, Float_t zPos, Float_t radL, Float_t density, Float_t thickness)
{
  //
  // Add additional layer to the list of layers (ordered by z position)
  // 
  if (!GetLayer(name)) {
    TString types = type;
    types.ToLower();
    KMCPolyLayer *newLayer = new KMCPolyLayer(name);
    newLayer->SetZ(zPos);
    newLayer->SetThickness(thickness);
    newLayer->SetX2X0( radL>0 ? thickness*density/radL : 0);
    newLayer->SetXTimesRho(thickness*density);
    newLayer->SetXRes(999999);
    newLayer->SetYRes(999999);
    newLayer->SetLayerEff(1);
    if  (types=="mag")  newLayer->SetType(KMCLayerFwd::kMAG);
    else newLayer->SetType(KMCLayerFwd::kDUMMY);
    newLayer->SetDead(kTRUE);
    //
    if (fLayers.GetEntries()==0) fLayers.Add(newLayer);
    else {      
      for (Int_t i = 0; i<fLayers.GetEntries(); i++) {
	KMCLayerFwd *l = GetLayer(i);
	if (zPos<l->GetZ()) { fLayers.AddBefore(l,newLayer); break; }
	if (zPos>l->GetZ() && (i+1)==fLayers.GetEntries() ) { fLayers.Add(newLayer); } // even bigger then last one
      }      
    }
    //
    newLayer->SetMaterial(0);
    return newLayer;
  }
  printf("Layer with the name %s does already exist\n",name);
  return 0;
}

//______________________________________________________________________________________
void KMCDetectorFwd::AddBeamPipe(Float_t r, Float_t dr, Float_t radL, Float_t density, NaMaterial* mat) 
{
  //
  // mock-up cyl. beam pipe
  fBeamPipe = new BeamPipe((char*)"BeamPipe");
  fBeamPipe->SetRadius(r);
  fBeamPipe->SetThickness(dr);
  fBeamPipe->SetX2X0( radL>0 ? dr*density/radL : 0);
  fBeamPipe->SetXTimesRho(dr*density);  
  fBeamPipe->SetDead(kTRUE);
  fBeamPipe->SetMaterial(mat);
  //
}

//________________________________________________________________________________
void KMCDetectorFwd::ClassifyLayers()
{
  // assign active Id's, etc
  fLastActiveLayer = -1;
  fLastActiveLayerITS = -1;
  fNActiveLayers = 0;
  fNActiveLayersITS = 0;
  fNActiveLayersMS = 0;
  fNActiveLayersTR = 0;  
  //
  fNLayers = fLayers.GetEntries();
  for (int il=0;il<fNLayers;il++) {
    KMCLayerFwd* lr = GetLayer(il);
    lr->SetID(il);
    if (!lr->IsDead()) {
      fLastActiveLayer = il; 
      lr->SetActiveID(fNActiveLayers++);
      if (lr->IsITS()) {
	fLastActiveLayerITS = il;
	fNActiveLayersITS++;
	fLayersITS.AddLast(lr);
      }
      if (lr->IsMS())   {
        if (fMSLrMinID<0) fMSLrMinID = lr->GetID();
        if (fMSLrMaxID<lr->GetID()) fMSLrMaxID = lr->GetID();
	fLastActiveLayerMS = il;
	fNActiveLayersMS++;
	fLayersMS.AddLast(lr);
      }
      if (lr->IsTrig()) {
	fNActiveLayersTR++;
	fLastActiveLayerMS = il;
	fLayersTR.AddLast(lr);
      }
    }
  }
  //
  TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* fldm = (MagField*) fld;
    fFldNReg = fldm->GetNReg();
    fFldZMins = fldm->GetZMin();
    fFldZMaxs = fldm->GetZMax();
    fDefStepAir = 100;
    fDefStepMat = 5;
  }
  else {
    fDefStepAir = 1;
    fDefStepMat = 1;
  }
  fZSteps = new Double_t[fFldNReg+2];
  //
  printf("DefStep: Air %f Mat: %f | N.MagField regions: %d\n",fDefStepAir,fDefStepMat,fFldNReg);
  //
  KMCProbeFwd::SetNITSLayers(fNActiveLayersITS + ((fVtx && !fVtx->IsDead()) ? 1:0));
  printf("KMCProbeFwd is initialized with %d slots\n",KMCProbeFwd::GetNITSLayers());
}

//_________________________________________________________________
KMCProbeFwd* KMCDetectorFwd::CreateProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // create track of given kinematics
  double xyz[3] = {x,y,z};
  double pxyz[3] = {pt*TMath::Cos(phi),pt*TMath::Sin(phi),TMath::Sqrt(pt*pt+mass*mass)*TMath::SinH(yrap)};
  KMCProbeFwd* probe = new KMCProbeFwd(xyz,pxyz,charge);
  probe->SetMass(mass);
  probe->SetTrID(-1);
  return probe;
}

//_________________________________________________________________
void KMCDetectorFwd::CreateProbe(KMCProbeFwd* adr, double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // create track of given kinematics
  double xyz[3] = {x,y,z};
  double pxyz[3] = {pt*TMath::Cos(phi),pt*TMath::Sin(phi),TMath::Sqrt(pt*pt+mass*mass)*TMath::SinH(yrap)};
  KMCProbeFwd* probe = new(adr) KMCProbeFwd(xyz,pxyz,charge);
  probe->SetTrID(-1);
  probe->SetMass(mass);
}

//_________________________________________________________________
KMCProbeFwd* KMCDetectorFwd::PrepareProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // Prepare trackable Kalman track at the farthest position
  //
  fGenPnt[0]=x;
  fGenPnt[1]=y;  
  fGenPnt[2]=z;  
  //
  if (fDecayZProf) { // decay is requested
    fZDecay  = fDecayZProf->GetRandom();
    fDecMode = kDoRealDecay;
    AliDebug(2,Form("Selected %.2f as decay Z",fZDecay));
  }
  // track parameters
  // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
  fProbe.Reset();
  KMCProbeFwd* probe = CreateProbe(pt,yrap,phi,mass,charge,x,y,z);
  fProbe = *probe;     // store original track
  //
  // propagate to last layer
  fLastActiveLayerTracked = 0;
  int resp=0;
  KMCLayerFwd* lr=0,*lrP=0;
  for (Int_t j=0; j<=fLastActiveLayer; j++) {
    lrP = lr;
    lr = GetLayer(j);
    lr->Reset();
    if (!lrP) continue;
    //
    if (!(resp=PropagateToLayer(probe,lrP,lr,1))) return 0;
    double r = probe->GetR();
    //    printf("L%2d %f %f %f\n",j,r, lr->GetRMin(),lr->GetRMax());
    if (r<lr->GetRMax() && r>lr->GetRMin()) {
      if (resp>0) {
	lr->AddCluster(probe->GetX(), probe->GetY(), probe->GetZ(),probe->GetTrID(), -1);
      }
      else lr->GetCorCluster()->Kill();
    }
    else lr->GetCorCluster()->Kill();
    //
    if (!lr->IsDead()) fLastActiveLayerTracked = j;
  }
  probe->ResetCovariance();// reset cov.matrix
  //  printf("Last active layer trracked: %d (out of %d)\n",fLastActiveLayerTracked,fLastActiveLayer);
  //
  return probe;
}

//____________________________________________________________________________
Int_t KMCDetectorFwd::GetFieldReg(double z) 
{
  // field region *2
  int ir = 0;
  for (int i=0;i<fFldNReg;i++) {
    if (z<=fFldZMins[i]) return ir;
    ir++;
    if (z<=fFldZMaxs[i]) return ir;
    ir++;
  }
  return ir;
}

//____________________________________________________________________________
Bool_t KMCDetectorFwd::PropagateToZBxByBz(KMCProbeFwd* trc,double z,double maxDZ,Double_t xOverX0,Double_t xTimesRho,Bool_t modeMC) 
{
  // propagate to given Z, checking for the field boundaries
  //
  double curZ = trc->GetZ();
  double dza = curZ-z;
  if (TMath::Abs(dza)<fgkFldEps) return kTRUE;

  // even id's correspond to Z between the field regions, odd ones - inside
  int ib0 = GetFieldReg(curZ); // field region id of start point
  int ib1 = GetFieldReg(z);    // field region id of last point
  int nzst = 0;
  //  AliInfo(Form("FldRegID: %d %d (%f : %f)",ib0,ib1, curZ,z));
  if (ib1>ib0) { // fwd propagation with field boundaries crossing
    for (int ib=ib0;ib<ib1;ib++) {
      if ( ib&0x1 ) { // we are in the odd (field ON) region, go till the end of field reg.
	//	printf("Here00 | %d %f\n",ib>>1,fFldZMaxs[ib>>1]);
	fZSteps[nzst++] = fFldZMaxs[ib>>1] + fgkFldEps;
      }
      else { // we are in even (field free) region, go till the beginning of next field reg.
	//	printf("Here01 | %d %f\n",ib>>1,fFldZMins[ib>>1]);
	fZSteps[nzst++] = fFldZMins[ib>>1] + fgkFldEps;
      }
    }
  }
  else if (ib1<ib0) { // bwd propagation
    for (int ib=ib0;ib>ib1;ib--) {
      if ( ib&0x1 ) { // we are in the odd (field ON) region, go till the beginning of field reg.
	//	printf("Here10 | %d %f\n",(ib-1)>>1,fFldZMins[(ib-1)>>1]);
	fZSteps[nzst++] = fFldZMins[(ib-1)>>1] - fgkFldEps;
      }
      else { // we are in even (field free) region, go till the beginning of next field reg.
	//	printf("Here11 | %d %f\n",(ib-1)>>1,fFldZMaxs[(ib-1)>>1]);
	fZSteps[nzst++] = fFldZMaxs[(ib-1)>>1] - fgkFldEps;
      }
    } 
  }
  fZSteps[nzst++] = z; // same field region, do nothing
  //
  //  printf("ZSteps: "); for (int ist=0;ist<nzst;ist++) printf("%+.5f ",fZSteps[ist]); printf("\n");

  for (int ist=0;ist<nzst;ist++) {
    double frc = (trc->GetZ()-fZSteps[ist])/dza;
    if (!trc->PropagateToZBxByBz(fZSteps[ist], maxDZ, frc*xOverX0, frc*xTimesRho, modeMC)) return kFALSE;
  }
  return kTRUE;
  //
}

//____________________________________________________________________________
Int_t KMCDetectorFwd::PropagateToLayer(KMCProbeFwd* trc, KMCLayerFwd* lrFrom, KMCLayerFwd* lrTo, int dir, Bool_t modeMC)
{
  // bring the track to lrTo, moving in direction dir (1: forward, -1: backward)
  // if relevant, account for the materials on lr
  //
  AliDebug(2,Form("From %d to %d, dir: %d",lrFrom? lrFrom->GetUniqueID():-1, lrTo->GetUniqueID(), dir));
  if (trc->GetTrack()->GetAlpha()<0) return -1;
  if ( lrFrom && dir*(lrTo->GetZ()-lrFrom->GetZ())<0 ) AliFatal(Form("Dir:%d Zstart: %f Zend: %f\n",dir,lrFrom->GetZ(),lrTo->GetZ()));
  //
  if      (dir>0 && (lrTo->GetZ()-0.5*lrTo->GetThickness()) < trc->GetZ() ) return -1;
  else if (dir<0 && (lrTo->GetZ()+0.5*lrTo->GetThickness()) > trc->GetZ() ) return -1;
  double dstZ;
  //
  if (lrFrom) {
    //
    if (!lrFrom->IsDead()) { // active layers are thin, no need for step by step tracking. The track is always in the middle
      AliDebug(2,Form("Correcting for mat.in active layer: X/X0: %f X*rho:%f ", lrFrom->GetX2X0(), dir*lrFrom->GetXTimesRho()));
      // note: for thin layer we ignore difference between the real BB and ETP eloss params
      if (!trc->CorrectForMeanMaterial(lrFrom->GetX2X0(), -dir*lrFrom->GetXTimesRho(), modeMC)) return 0; 
    }
    else {
      //
      dstZ = lrFrom->GetZ()+0.5*dir*lrFrom->GetThickness(); // go till the end of starting layer applying corrections
      if (dir==1 && trc->GetZ()<=fZDecay && dstZ>fZDecay) { // need to perform or to apply decay
	double frac = (fZDecay-trc->GetZ())/lrFrom->GetThickness();
	// account for the difference between real BB and ETP param eloss
        
	double corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
        float x2x0=0,xrho=0;
        lrFrom->getMatBudget(trc->GetX(),trc->GetY(), x2x0, xrho);
        if (fIntegrateMSX2X0) trc->addMSMaterialsSeen(x2x0);
	if (!PropagateToZBxByBz(trc,fZDecay, fDefStepMat, frac*x2x0, -frac*xrho*corrELoss, modeMC)) return 0;
	PerformDecay(trc);
	frac = 1.-frac;
	corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
        lrFrom->getMatBudget(trc->GetX(),trc->GetY(), x2x0, xrho);
        if (fIntegrateMSX2X0) trc->addMSMaterialsSeen(x2x0);
	if (!PropagateToZBxByBz(trc,dstZ, fDefStepMat, frac*x2x0, -frac*xrho*corrELoss, modeMC)) return 0;
      }
      else {
	// account for the difference between real BB and ETP param eloss
	double corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
        float x2x0=0,xrho=0;
        lrFrom->getMatBudget(trc->GetX(),trc->GetY(), x2x0, xrho);
        if (fIntegrateMSX2X0) trc->addMSMaterialsSeen(x2x0);
	if (!PropagateToZBxByBz(trc,dstZ, fDefStepMat, x2x0, -dir*xrho*corrELoss, modeMC)) return 0;
      }
    }
  }
  //
  dstZ = lrTo->GetZ();
  if (lrTo->IsDead()) dstZ += -dir*lrTo->GetThickness()/2; // for thick dead layers go till entrance
  //
  if (dir==1 && trc->GetZ()<=fZDecay && dstZ>fZDecay) { // need to perform or to apply decay
    if (!PropagateToZBxByBz(trc,fZDecay, fDefStepAir)) return 0;
    PerformDecay(trc);
  }
  if (!PropagateToZBxByBz(trc,dstZ, fDefStepAir)) return 0;
  //
  if (AliLog::GetGlobalDebugLevel()>=2) trc->GetTrack()->Print();

  return 1;
}

//________________________________________________________________________________
Bool_t KMCDetectorFwd::SolveSingleTrackViaKalman(double pt, double yrap, double phi, 
						 double mass, int charge, double x,double y, double z)
{
  // analytical estimate of tracking resolutions
  //  fProbe.SetUseLogTermMS(kTRUE);
  //
  for (int i=3;i--;) fGenPnt[i] = 0;
  if (fMinITSHits>fNActiveLayersITS) {
    fMinITSHits = fLastActiveLayerITS; 
    printf("Redefined request of min N ITS hits to %d\n",fMinITSHits);
  }
  //
  KMCProbeFwd* probe = PrepareProbe(pt,yrap,phi,mass,charge,x,y,z);
  if (!probe) return kFALSE;
  //
  KMCLayerFwd *lr = 0;
  //
  // Start the track fitting --------------------------------------------------------
  //
  // Back-propagate the covariance matrix along the track. 
  // Kalman loop over the layers
  //
  KMCProbeFwd* currTr = 0;
  lr = GetLayer(fLastActiveLayerTracked);
  lr->SetAnProbe(*probe);
  delete probe; // rethink...
  //
  int nupd = 0;
  for (Int_t j=fLastActiveLayerTracked; j--; ) {  // Layer loop
    //
    KMCLayerFwd *lrP = lr;
    lr = GetLayer(j);
    //
    lr->SetAnProbe( *lrP->GetAnProbe() );
    currTr = lr->GetAnProbe();
    currTr->ResetHit(lrP->GetActiveID());
    //
    // if there was a measurement on prev layer, update the track
    if (!lrP->IsDead()) { // include "ideal" measurement
      //      printf("Before update on %d : ",j); currTr->Print("etp");
      KMCClusterFwd* cl = lrP->GetCorCluster();
      if (!cl->IsKilled()) {
	if (!UpdateTrack(currTr,lrP,cl))  return kFALSE;
	nupd++;
      }
      //      printf("After update on %d (%+e %+e) : ",j, lrP->GetXRes(),lrP->GetYRes()); currTr->Print("etp");

    }

    if (!PropagateToLayer(currTr,lrP,lr,-1)) return kFALSE;      // propagate to current layer
    //
  } // end loop over layers
  // is MS reco ok?
  // check trigger
  int nhMS=0,nhTR=0,nhITS=0;
  for (int ilr=fNLayers;ilr--;) {
    KMCLayerFwd *lrt = GetLayer(ilr);
    if (lrt->IsTrig()) if (!lrt->GetCorCluster()->IsKilled()) nhTR++;
    if (lrt->IsMS())   if (!lrt->GetCorCluster()->IsKilled()) nhMS++;
    if (lrt->IsITS())  if (!lrt->GetCorCluster()->IsKilled()) nhITS++;
  }
  //  printf("ITS: %d MS: %d TR: %d\n",nhITS,nhMS,nhTR);
  if (nhTR<fMinTRHits) return kFALSE;
  if (nhMS<fMinMSHits) return kFALSE;
  //
  return kTRUE;
}

//____________________________________________________________________________
Bool_t KMCDetectorFwd::UpdateTrack(KMCProbeFwd* trc, const KMCLayerFwd* lr, const KMCClusterFwd* cl) const
{
  // update track with measured cluster
  // propagate to cluster
  if (cl->IsKilled()) return kTRUE;
  // Note: we are working in the tracking frame: Lab X,Y,Z  <->  Tracking -Z,Y,X
  double meas[2] = {cl->GetYTF(), cl->GetZTF()}; // ideal cluster coordinate, tracking (AliExtTrParam frame)
  double measErr2[3] = {cl->GetSigYY(), cl->GetSigYZ(), cl->GetSigZZ()};
  //
  //  if (trc->GetZ()>620 && trc->GetZ()<1140) {
  //  printf("Update phi:%f Meas: %f %f Err: %f %f %f\n",TMath::ATan2(trc->GetY(), trc->GetX()),meas[0], meas[1], measErr2[0], measErr2[1], measErr2[2]);
  //  }
  double chi2 = trc->GetTrack()->GetPredictedChi2(meas,measErr2);
  //  printf("Update for lr:%s -> chi2=%f\n",lr->GetName(), chi2);
  //  printf("cluster at Lr:%s was [%e %e / %e %e %e]: ", lr->GetName(), meas[0],meas[1], measErr2[0], measErr2[1], measErr2[2]); cl->Print("lc");
  //  printf("track   was : %e %e\n", trc->GetY(), trc->GetZ()); trc->Print("etp");
  if (chi2>fMaxChi2Cl) return kTRUE; // chi2 is too large
    
  if (!trc->Update(meas,measErr2)) {
    AliDebug(2,Form("layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		    lr->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
    if (AliLog::GetGlobalDebugLevel()>1) trc->Print("l");
    return kFALSE;
  }
  trc->AddHit(lr, chi2, cl->GetTrID());
  //
  return kTRUE;
}

//________________________________________________________________________________
Bool_t KMCDetectorFwd::SolveSingleTrack(double pt, double yrap, double phi, 
					double mass, int charge, double x,double y, double z, 
					TObjArray* sumArr,int nMC, int offset)
{
  // analityc and fullMC (nMC trials) evaluaion of tracks with given kinematics.
  // the results are filled in KMCTrackSummary objects provided via summArr array
  //
  //
  if (!SolveSingleTrackViaKalman(pt,yrap,phi,mass,charge,x,y,z)) return kFALSE;
  //
  /*
  int nsm = sumArr ? sumArr->GetEntriesFast() : 0;
  KMCLayerFwd* vtx = GetLayer(0);
  */
  //
  /*RS
  for (int i=0;i<nsm;i++) {
    KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(i);
    if (!tsm) continue;
    tsm->SetRefProbe( GetProbeTrack() ); // attach reference track (generated)
    tsm->SetAnProbe( vtx->GetAnProbe() ); // attach analitycal solution
  }
  */
  //
  if (offset<0) offset = fNLayers;

  TStopwatch sw;
  sw.Start();
  for (int it=0;it<nMC;it++) {
    //    printf("ev: %d\n",it);
    SolveSingleTrackViaKalmanMC(offset);
    /*RS
    KMCProbeFwd* trc = vtx->GetWinnerMCTrack();
    for (int ism=nsm;ism--;) { // account the track in each of summaries
      KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(ism);
      if (!tsm) continue;
      tsm->AddUpdCalls(GetUpdCalls());
      tsm->AddTrack(trc); 
    }
    */
  }
  //
  sw.Stop();
  //  printf("Total time: "); sw.Print();
  return kTRUE;
}

//____________________________________________________________________________
void KMCDetectorFwd::ResetMCTracks(Int_t maxLr)
{
  if (maxLr<0 || maxLr>=fNLayers) maxLr = fNLayers-1;
  for (int i=maxLr+1;i--;) GetLayer(i)->ResetMCTracks();
}

//____________________________________________________________
void KMCDetectorFwd::TrackMS()
{
  // RSTMP
  KMCLayerFwd* lr0 = GetLayerMS(fNActiveLayersMS/2);
  KMCLayerFwd* lr1 = GetLayerMS(fNActiveLayersMS-1);
  for (int i0=-1;i0<lr0->GetNBgClusters();i0++) {
    KMCClusterFwd* cl0 = lr0->GetCluster(i0);
    if (cl0->IsKilled()) continue;
    for (int i1=-1;i1<lr1->GetNBgClusters();i1++) {
      KMCClusterFwd* cl1 = lr1->GetCluster(i1);
      if (cl1->IsKilled()) continue;
      //
      KMCProbeFwd* seed = CreateMSSeed(cl0, cl1);
      if (!seed) continue;
      // check if there are TR hits matching to seed
      if ( no TR ) {
	delete seed;
	continue;
      }
      
      delete seed;
    }
  }
}

//____________________________________________________________
Bool_t KMCDetectorFwd::SolveSingleTrackViaKalmanMC(int offset)
{
  // MC estimate of tracking resolutions/effiencies. Requires that the SolveSingleTrackViaKalman
  // was called before, since it uses data filled by this method
  //
  // The MC tracking will be done starting from fLastActiveLayerITS + offset (before analytical estimate will be used)
  //
  // At this point, the fProbe contains the track params generated at vertex.
  // Clone it and propagate to target layer to generate hit positions affected by MS
  //
  const float kErrScale = 100.; // RS: this is the parameter defining the initial cov.matrix error wrt sensor resolution
      
  Bool_t checkMS = kTRUE;
  fMuTrackVertex.SetUniqueID(999); // invalidate
  fMuTrackBCVertex.SetUniqueID(999); // invalidate
  fMuTrackBCLastITS.SetUniqueID(999); // invalidate
  fMuTrackLastITS.SetUniqueID(999); // invalidate
  //
  KMCProbeFwd *currTrP=0,*currTr=0;
  static KMCProbeFwd trcConstr;
  int maxLr = fLastActiveLayerITS;
  if (offset>0) maxLr += offset;
  if (maxLr>fLastActiveLayer) maxLr = fLastActiveLayer;
  if (fExternalInput) maxLr = fLastActiveLayerTracked;
  if (maxLr<0) return kFALSE;
  if (fVtx && !fImposeVertexPosition) {
    fVtx->GetMCCluster()->Set(fProbe.GetX(),fProbe.GetY(), fProbe.GetZ());
    fRefVtx[0] = fProbe.GetX();
    fRefVtx[1] = fProbe.GetY();
    fRefVtx[2] = fProbe.GetZ();
  }
  //
  //  printf("MaxLr: %d\n",maxLr);
  ResetMCTracks(maxLr);
  KMCLayerFwd* lr = GetLayer(maxLr);
  currTr = lr->AddMCTrack(&fProbe); // start with original track at vertex
  //
  //  printf("INI SEED: "); currTr->Print("etp");
  if (!fExternalInput) {if (!TransportKalmanTrackWithMS(currTr, maxLr, kFALSE)) return kFALSE;} // transport it to outermost layer where full MC is done
  else *currTr->GetTrack() = *GetLayer(maxLr)->GetAnProbe()->GetTrack();
  //printf("LastTrackedMS: "); currTr->GetTrack()->Print();

  PrepareForTracking();

  //
  if (maxLr<=fLastActiveLayerTracked && maxLr>fLastActiveLayerITS) { // prolongation from MS
    // start from correct track propagated from above till maxLr
    double *covMS = (double*)currTr->GetTrack()->GetCovariance();
    const double *covIdeal = GetLayer(maxLr)->GetAnProbe()->GetCovariance();
    for (int i=15;i--;) covMS[i] = covIdeal[i];
  }
  else { // ITS SA: randomize the starting point
    double r = currTr->GetR();
    currTr->ResetCovariance( kErrScale*TMath::Sqrt(lr->GetXRes(r)*lr->GetYRes(r)) ); // RS: this is the coeff to play with
  }
  //
  int fst = 0;
  const int fstLim = -1;

  // RSTMP
  KMCClusterFwd* cl0 = GetLayerMS(fNActiveLayersMS/2)->GetCluster(-1);
  KMCClusterFwd* cl1 = GetLayerMS(fNActiveLayersMS-1)->GetCluster(-1);
  if (!cl0->IsKilled() && !cl1->IsKilled()) {
    CreateMSSeed(-1,-1);
  }

  if (maxLr>=fLastActiveLayerMS) {
    TrackMS();
  }

  
  for (Int_t j=maxLr; j--; ) {  // Layer loop
    //
    int ncnd=0,cndCorr=-1;
    KMCLayerFwd *lrP = lr;
    lr = GetLayer(j);
    fIntegrateMSX2X0 = (j>fMSLrMinID && j<fMSLrMaxID);
    int ntPrev = lrP->GetNMCTracks();
    //
    //    printf("Lr:%d %s IsDead:%d\n",j, lrP->GetName(),lrP->IsDead());
    if (lrP->IsDead()) { // for passive layer just propagate the copy of all tracks of prev layer >>>
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); 
	if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (fst<fstLim) {
	  fst++;
	  currTr->Print("etp");
	}
	if (!PropagateToLayer(currTr,lrP,lr,-1))  {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of dead layer <<<
    //
    if (lrP->IsMS() || lrP->IsTrig()) { // we don't consider bg hits in MS, just update with MC cluster
      KMCClusterFwd* clmc = lrP->GetMCCluster();
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	//	printf("At lr %d  | ",j+1); currTrP->Print("etp");
	currTr = lr->AddMCTrack( currTrP );
	if (fst<fstLim) {
	  fst++;
	  currTr->Print("etp");
	}
	// printf("\nAt Lr:%s ",lrP->GetName()); currTr->GetTrack()->Print();
	if (!clmc->IsKilled()) {
	  //	  lrP->Print("");
	  //	  printf("BeforeMS Update: "); currTr->Print("etp");
	  if (!UpdateTrack(currTr, lrP, clmc)) {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // update with correct MC cl.
	  //	  printf("AfterMS Update: "); currTr->Print("etp");
	}
	if (!PropagateToLayer(currTr,lrP,lr,-1))            {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // propagate to current layer	
      }      
      //      currTr->Print("etp");
      continue;
    } // treatment of ideal layer <<<
    //
    // active layer under eff. study (ITS?): propagate copy of every track to MC cluster frame (to have them all in the same frame)
    // and calculate the limits of bg generation
    //    KMCClusterFwd* clMC = lrP->GetMCCluster();
    //int nseeds = 0;
    // here ITS layers
    //   
    //
    AliDebug(2,Form("From Lr: %d | %d seeds, %d bg clusters",j+1,ntPrev,lrP->GetNBgClusters()));
    for (int itrP=0;itrP<ntPrev;itrP++) { // loop over all tracks from previous layer
      currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
      //
      if (checkMS) {
	checkMS = kFALSE;
	// check if muon track is well defined
	if (currTrP->GetNTRHits()<fMinTRHits) {currTrP->Kill(); continue;}
	if (currTrP->GetNMSHits()<fMinMSHits) {currTrP->Kill(); continue;}
	//
	//	printf("Check %d of %d at lr%d Nhits:%d NTR:%d NMS:%d\n",
	//       itrP,ntPrev,j, currTrP->GetNHits(),currTrP->GetNTRHits(),currTrP->GetNMSHits());
	if (fHChi2MS) fHChi2MS->Fill(currTr->GetChi2(),currTr->GetNHits());      
      }
      //
      // Are we entering to the last ITS layer? Apply Branson plane correction if requested
      if (lrP->GetID() == fLastActiveLayerITS && fVtx && !fVtx->IsDead() && fApplyBransonPCorrection>=0) {
	//	printf("%e -> %e (%d %d) | %e\n",lrP->GetZ(),lr->GetZ(), lrP->GetID(),j, currTr->GetZ());

	trcConstr = *currTrP;
	fMuTrackLastITS = trcConstr;
	if (!PropagateToLayer(&trcConstr,lrP,fVtx,-1))  {currTrP->Kill();continue;} // propagate to vertex
	//////	trcConstr.ResetCovariance();
	// update with vertex point + eventual additional error
	float origErrX = fVtx->GetYRes(), origErrY = fVtx->GetXRes(); // !!! Ylab<->Ytracking, Xlab<->Ztracking
	KMCClusterFwd* clv = fVtx->GetMCCluster();
	double measCV[2] = {clv->GetYTF(), -clv->GetZTF()}, errCV[3] = {
          origErrY*origErrY+fApplyBransonPCorrection*fApplyBransonPCorrection,
          0.,
	  origErrX*origErrX+fApplyBransonPCorrection*fApplyBransonPCorrection
	};
	fChi2MuVtx = trcConstr.GetPredictedChi2(measCV,errCV);
	if (fHChi2Branson) fHChi2Branson->Fill(fChi2MuVtx);
	//	printf("UpdVtx: {%+e %+e}/{%e %e %e}\n",measCV[0],measCV[1],errCV[0],errCV[1],errCV[2]);
	//	printf("Muon@Vtx:  "); trcConstr.Print("etp");

	fMuTrackVertex = trcConstr;
	if (!trcConstr.Update(measCV,errCV)) {currTrP->Kill();continue;}
	//	printf("Truth@VTX: "); fProbe.Print("etp");
	//	printf("Constraint@VTX: "); trcConstr.Print("etp");
	fMuTrackBCVertex = trcConstr;
	fMuTrackBCVertex.SetUniqueID(0);
	
	
	if (!PropagateToLayer(&trcConstr,fVtx,lrP,1)) {currTrP->Kill();continue;}
	// constrain Muon Track
	//	printf("Constraint: "); trcConstr.Print("etp");

	//////	double measCM[2] = {trcConstr.GetYLoc(), trcConstr.GetZLoc()}, errCM[3] = {
	//////	  trcConstr.GetSigmaY2(),trcConstr.GetSigmaXY(),trcConstr.GetSigmaX2()
	//////	};

	//////	printf("UpdMuon: {%+e %+e}/{%e %e %e}\n",measCM[0],measCM[1],errCM[0],errCM[1],errCM[2]);
	//////	printf("Before constraint: "); currTrP->Print("etp");
	//////	if (!currTrP->Update(measCM,errCM)) {currTrP->Kill();continue;}
	fMuTrackBCLastITS = trcConstr;
	fMuTrackBCLastITS.SetUniqueID(0);
	(*currTrP->GetTrack()) = *trcConstr.GetTrack(); // override with constraint
	//	printf("After  constraint: "); currTrP->Print("etp");
	//	printf("MuTruth "); lrP->GetAnProbe()->Print("etp");
      }
      
      currTr = lr->AddMCTrack( currTrP );
      
      if (fst<fstLim) {
	fst++;
	currTr->Print("etp");
      }
      //
      AliDebug(2,Form("LastChecked before:%d",currTr->GetInnerLayerChecked()));
      CheckTrackProlongations(currTr, lrP,lr);
      AliDebug(2,Form("LastChecked after:%d",currTr->GetInnerLayerChecked()));
      ncnd++;
      if (currTr->GetNFakeITSHits()==0 && cndCorr<ncnd) cndCorr=ncnd;
      if (NeedToKill(currTr)) {currTr->Kill(); continue;}
    }
    /*
    if (ncnd>100) {
      printf("\n NCND= %d NPrev=  %d \n",ncnd,ntPrev);
      currTrP->Print("etp");
      Print("cl");
    }
    */
    if (fHNCand)     fHNCand->Fill(lrP->GetActiveID(), ncnd);
    if (fHCandCorID) fHCandCorID->Fill(lrP->GetActiveID(), cndCorr);
    //
    //  
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
    if (ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0) {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--)  lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    //
    for (int itr=ntTot;itr--;) {
      currTr = lr->GetMCTrack(itr);
      if (!PropagateToLayer(currTr,lrP,lr,-1))  {currTr->Kill();continue;} // propagate to current layer
    }
    AliDebug(1,Form("Got %d tracks on layer %s",ntTot,lr->GetName()));
    //    lr->GetMCTracks()->Print();
    //
  } // end loop over layers    
  //
  // do we use vertex constraint?
  if (fVtx && !fVtx->IsDead() && fIncludeVertex) {
    printf("Apply vertex constraint\n");
    int ntr = fVtx->GetNMCTracks();
    for (int itr=0;itr<ntr;itr++) {
      currTr = fVtx->GetMCTrack(itr);
      if (currTr->IsKilled()) continue;
      double meas[2] = {0.,0.};
      if (fImposeVertexPosition) {
	meas[0] = fRefVtx[1];
	meas[1] =-fRefVtx[0];
      }
      else {
	KMCClusterFwd* clv = fVtx->GetMCCluster();
	meas[0] = clv->GetYTF();
	meas[1] = clv->GetZTF();
      }
      double measErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(), 0, fVtx->GetXRes()*fVtx->GetXRes()};
      //double chi2v = currTr->GetPredictedChi2(meas,measErr2);
      if (!currTr->Update(meas,measErr2)) continue;
      //currTr->AddHit(fVtx->GetActiveID(), chi2v, -1);
      currTr->SetInnerLrChecked(fVtx->GetActiveID());
      //if (NeedToKill(currTr)) currTr->Kill();
      // if (fVtx->IsITS()) {if (!UpdateTrack(currTr, fVtx, fVtx->GetMCCluster(), kFALSE)) {currTr->Kill();continue;}}
    }
  }
  //  EliminateUnrelated(); //RSS
  int ntTot = lr->GetNMCTracks();
  ntTot = TMath::Min(1,ntTot);
  for (int itr=ntTot;itr--;) {
    currTr = lr->GetMCTrack(itr);
    if (currTr->IsKilled()) continue;
    if (fHChi2NDFCorr&&fHChi2NDFFake) {
      if (IsCorrect(currTr)) fHChi2NDFCorr->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
      else                   fHChi2NDFFake->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
    }
  }
  //  
  return kTRUE;
}

/*

//____________________________________________________________
Bool_t KMCDetectorFwd::SolveSingleTrackViaKalmanMC(int offset)
{
  // MC estimate of tracking resolutions/effiencies. Requires that the SolveSingleTrackViaKalman
  // was called before, since it uses data filled by this method
  //
  // The MC tracking will be done starting from fLastActiveLayerITS + offset (before analytical estimate will be used)
  //
  // At this point, the fProbe contains the track params generated at vertex.
  // Clone it and propagate to target layer to generate hit positions affected by MS
  //
  const float kErrScale = 100.; // RS: this is the parameter defining the initial cov.matrix error wrt sensor resolution
      
  Bool_t checkMS = kTRUE;
  fMuTrackVertex.SetUniqueID(999); // invalidate
  fMuTrackBCVertex.SetUniqueID(999); // invalidate
  fMuTrackBCLastITS.SetUniqueID(999); // invalidate
  fMuTrackLastITS.SetUniqueID(999); // invalidate
  //
  KMCProbeFwd *currTrP=0,*currTr=0;
  static KMCProbeFwd trcConstr;
  int maxLr = fLastActiveLayerITS;
  if (offset>0) maxLr += offset;
  if (maxLr>fLastActiveLayer) maxLr = fLastActiveLayer;
  if (fExternalInput) maxLr = fLastActiveLayerTracked;
  if (maxLr<0) return kFALSE;
  if (fVtx && !fImposeVertexPosition) {
    fVtx->GetMCCluster()->Set(fProbe.GetX(),fProbe.GetY(), fProbe.GetZ());
    fRefVtx[0] = fProbe.GetX();
    fRefVtx[1] = fProbe.GetY();
    fRefVtx[2] = fProbe.GetZ();
  }
  //
  //  printf("MaxLr: %d\n",maxLr);
  ResetMCTracks(maxLr);
  KMCLayerFwd* lr = GetLayer(maxLr);
  currTr = lr->AddMCTrack(&fProbe); // start with original track at vertex
  //
  //  printf("INI SEED: "); currTr->Print("etp");
  if (!fExternalInput) {if (!TransportKalmanTrackWithMS(currTr, maxLr, kFALSE)) return kFALSE;} // transport it to outermost layer where full MC is done
  else *currTr->GetTrack() = *GetLayer(maxLr)->GetAnProbe()->GetTrack();
  //printf("LastTrackedMS: "); currTr->GetTrack()->Print();

  PrepareForTracking();

  //
  if (maxLr<=fLastActiveLayerTracked && maxLr>fLastActiveLayerITS) { // prolongation from MS
    // start from correct track propagated from above till maxLr
    double *covMS = (double*)currTr->GetTrack()->GetCovariance();
    const double *covIdeal = GetLayer(maxLr)->GetAnProbe()->GetCovariance();
    for (int i=15;i--;) covMS[i] = covIdeal[i];
  }
  else { // ITS SA: randomize the starting point
    double r = currTr->GetR();
    currTr->ResetCovariance( kErrScale*TMath::Sqrt(lr->GetXRes(r)*lr->GetYRes(r)) ); // RS: this is the coeff to play with
  }
  //
  int fst = 0;
  const int fstLim = -1;
  
  for (Int_t j=maxLr; j--; ) {  // Layer loop
    //
    int ncnd=0,cndCorr=-1;
    KMCLayerFwd *lrP = lr;
    lr = GetLayer(j);
    fIntegrateMSX2X0 = (j>fMSLrMinID && j<fMSLrMaxID);
    int ntPrev = lrP->GetNMCTracks();
    //
    //    printf("Lr:%d %s IsDead:%d\n",j, lrP->GetName(),lrP->IsDead());
    if (lrP->IsDead()) { // for passive layer just propagate the copy of all tracks of prev layer >>>
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); 
	if (currTrP->IsKilled()) continue;
	currTr = lr->AddMCTrack( currTrP );
	if (fst<fstLim) {
	  fst++;
	  currTr->Print("etp");
	}
	if (!PropagateToLayer(currTr,lrP,lr,-1))  {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // propagate to current layer
      }
      continue;
    } // treatment of dead layer <<<
    //
    if (lrP->IsMS() || lrP->IsTrig()) { // we don't consider bg hits in MS, just update with MC cluster
      KMCClusterFwd* clmc = lrP->GetMCCluster();
      for (int itrP=ntPrev;itrP--;) { // loop over all tracks from previous layer
	currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
	//	printf("At lr %d  | ",j+1); currTrP->Print("etp");
	currTr = lr->AddMCTrack( currTrP );
	if (fst<fstLim) {
	  fst++;
	  currTr->Print("etp");
	}
	// printf("\nAt Lr:%s ",lrP->GetName()); currTr->GetTrack()->Print();
	if (!clmc->IsKilled()) {
	  //	  lrP->Print("");
	  //	  printf("BeforeMS Update: "); currTr->Print("etp");
	  if (!UpdateTrack(currTr, lrP, clmc)) {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // update with correct MC cl.
	  //	  printf("AfterMS Update: "); currTr->Print("etp");
	}
	if (!PropagateToLayer(currTr,lrP,lr,-1))            {currTr->Kill(); lr->GetMCTracks()->RemoveLast(); continue;} // propagate to current layer	
      }      
      //      currTr->Print("etp");
      continue;
    } // treatment of ideal layer <<<
    //
    // active layer under eff. study (ITS?): propagate copy of every track to MC cluster frame (to have them all in the same frame)
    // and calculate the limits of bg generation
    //    KMCClusterFwd* clMC = lrP->GetMCCluster();
    //int nseeds = 0;
    // here ITS layers
    //   
    //
    AliDebug(2,Form("From Lr: %d | %d seeds, %d bg clusters",j+1,ntPrev,lrP->GetNBgClusters()));
    for (int itrP=0;itrP<ntPrev;itrP++) { // loop over all tracks from previous layer
      currTrP = lrP->GetMCTrack(itrP); if (currTrP->IsKilled()) continue;
      //
      if (checkMS) {
	checkMS = kFALSE;
	// check if muon track is well defined
	if (currTrP->GetNTRHits()<fMinTRHits) {currTrP->Kill(); continue;}
	if (currTrP->GetNMSHits()<fMinMSHits) {currTrP->Kill(); continue;}
	//
	//	printf("Check %d of %d at lr%d Nhits:%d NTR:%d NMS:%d\n",
	//       itrP,ntPrev,j, currTrP->GetNHits(),currTrP->GetNTRHits(),currTrP->GetNMSHits());
	if (fHChi2MS) fHChi2MS->Fill(currTr->GetChi2(),currTr->GetNHits());      
      }
      //
      // Are we entering to the last ITS layer? Apply Branson plane correction if requested
      if (lrP->GetID() == fLastActiveLayerITS && fVtx && !fVtx->IsDead() && fApplyBransonPCorrection>=0) {
	//	printf("%e -> %e (%d %d) | %e\n",lrP->GetZ(),lr->GetZ(), lrP->GetID(),j, currTr->GetZ());

	trcConstr = *currTrP;
	fMuTrackLastITS = trcConstr;
	if (!PropagateToLayer(&trcConstr,lrP,fVtx,-1))  {currTrP->Kill();continue;} // propagate to vertex
	//////	trcConstr.ResetCovariance();
	// update with vertex point + eventual additional error
	float origErrX = fVtx->GetYRes(), origErrY = fVtx->GetXRes(); // !!! Ylab<->Ytracking, Xlab<->Ztracking
	KMCClusterFwd* clv = fVtx->GetMCCluster();
	double measCV[2] = {clv->GetYTF(), -clv->GetZTF()}, errCV[3] = {
          origErrY*origErrY+fApplyBransonPCorrection*fApplyBransonPCorrection,
          0.,
	  origErrX*origErrX+fApplyBransonPCorrection*fApplyBransonPCorrection
	};
	fChi2MuVtx = trcConstr.GetPredictedChi2(measCV,errCV);
	if (fHChi2Branson) fHChi2Branson->Fill(fChi2MuVtx);
	//	printf("UpdVtx: {%+e %+e}/{%e %e %e}\n",measCV[0],measCV[1],errCV[0],errCV[1],errCV[2]);
	//	printf("Muon@Vtx:  "); trcConstr.Print("etp");

	fMuTrackVertex = trcConstr;
	if (!trcConstr.Update(measCV,errCV)) {currTrP->Kill();continue;}
	//	printf("Truth@VTX: "); fProbe.Print("etp");
	//	printf("Constraint@VTX: "); trcConstr.Print("etp");
	fMuTrackBCVertex = trcConstr;
	fMuTrackBCVertex.SetUniqueID(0);
	
	
	if (!PropagateToLayer(&trcConstr,fVtx,lrP,1)) {currTrP->Kill();continue;}
	// constrain Muon Track
	//	printf("Constraint: "); trcConstr.Print("etp");

	//////	double measCM[2] = {trcConstr.GetYLoc(), trcConstr.GetZLoc()}, errCM[3] = {
	//////	  trcConstr.GetSigmaY2(),trcConstr.GetSigmaXY(),trcConstr.GetSigmaX2()
	//////	};

	//////	printf("UpdMuon: {%+e %+e}/{%e %e %e}\n",measCM[0],measCM[1],errCM[0],errCM[1],errCM[2]);
	//////	printf("Before constraint: "); currTrP->Print("etp");
	//////	if (!currTrP->Update(measCM,errCM)) {currTrP->Kill();continue;}
	fMuTrackBCLastITS = trcConstr;
	fMuTrackBCLastITS.SetUniqueID(0);
	(*currTrP->GetTrack()) = *trcConstr.GetTrack(); // override with constraint
	//	printf("After  constraint: "); currTrP->Print("etp");
	//	printf("MuTruth "); lrP->GetAnProbe()->Print("etp");
      }
      
      currTr = lr->AddMCTrack( currTrP );
      
      if (fst<fstLim) {
	fst++;
	currTr->Print("etp");
      }
      //
      AliDebug(2,Form("LastChecked before:%d",currTr->GetInnerLayerChecked()));
      CheckTrackProlongations(currTr, lrP,lr);
      AliDebug(2,Form("LastChecked after:%d",currTr->GetInnerLayerChecked()));
      ncnd++;
      if (currTr->GetNFakeITSHits()==0 && cndCorr<ncnd) cndCorr=ncnd;
      if (NeedToKill(currTr)) {currTr->Kill(); continue;}
    }
    //    if (ncnd>100) {
    //      printf("\n NCND= %d NPrev=  %d \n",ncnd,ntPrev);
    //      currTrP->Print("etp");
    //      Print("cl");
    //    }
    if (fHNCand)     fHNCand->Fill(lrP->GetActiveID(), ncnd);
    if (fHCandCorID) fHCandCorID->Fill(lrP->GetActiveID(), cndCorr);
    //
    //  
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
    if (ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0) {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--)  lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    //
    for (int itr=ntTot;itr--;) {
      currTr = lr->GetMCTrack(itr);
      if (!PropagateToLayer(currTr,lrP,lr,-1))  {currTr->Kill();continue;} // propagate to current layer
    }
    AliDebug(1,Form("Got %d tracks on layer %s",ntTot,lr->GetName()));
    //    lr->GetMCTracks()->Print();
    //
  } // end loop over layers    
  //
  // do we use vertex constraint?
  if (fVtx && !fVtx->IsDead() && fIncludeVertex) {
    printf("Apply vertex constraint\n");
    int ntr = fVtx->GetNMCTracks();
    for (int itr=0;itr<ntr;itr++) {
      currTr = fVtx->GetMCTrack(itr);
      if (currTr->IsKilled()) continue;
      double meas[2] = {0.,0.};
      if (fImposeVertexPosition) {
	meas[0] = fRefVtx[1];
	meas[1] =-fRefVtx[0];
      }
      else {
	KMCClusterFwd* clv = fVtx->GetMCCluster();
	meas[0] = clv->GetYTF();
	meas[1] = clv->GetZTF();
      }
      double measErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(), 0, fVtx->GetXRes()*fVtx->GetXRes()};
      //double chi2v = currTr->GetPredictedChi2(meas,measErr2);
      if (!currTr->Update(meas,measErr2)) continue;
      //currTr->AddHit(fVtx->GetActiveID(), chi2v, -1);
      currTr->SetInnerLrChecked(fVtx->GetActiveID());
      //if (NeedToKill(currTr)) currTr->Kill();
      // if (fVtx->IsITS()) {if (!UpdateTrack(currTr, fVtx, fVtx->GetMCCluster(), kFALSE)) {currTr->Kill();continue;}}
    }
  }
  //  EliminateUnrelated(); //RSS
  int ntTot = lr->GetNMCTracks();
  ntTot = TMath::Min(1,ntTot);
  for (int itr=ntTot;itr--;) {
    currTr = lr->GetMCTrack(itr);
    if (currTr->IsKilled()) continue;
    if (fHChi2NDFCorr&&fHChi2NDFFake) {
      if (IsCorrect(currTr)) fHChi2NDFCorr->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
      else                   fHChi2NDFFake->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
    }
  }
  //  
  return kTRUE;
}

 */

//________________________________________________________________________________
Bool_t KMCDetectorFwd::TransportKalmanTrackWithMS(KMCProbeFwd *probTr, int maxLr, Bool_t bg)
{
  // Transport track till layer maxLr, applying random MS
  //
  int resP = 0;
  for (Int_t j=0; j<maxLr; j++) {
    KMCLayerFwd* lr0 = GetLayer(j);
    KMCLayerFwd* lr  = GetLayer(j+1);
    //
    if (!bg) {
      lr->GetMCCluster()->Kill(); // will be reactivated if accepted
    }
    if (!(resP=PropagateToLayer(probTr,lr0,lr, 1, kTRUE))) return kFALSE;
    if (lr->IsDead()) continue;
    if (resP<0) continue;
    lr->AddCluster(probTr->GetX(), probTr->GetY(), probTr->GetZ(), probTr->GetTrID(), bg ? 1:0);
    //
  }
  //
  return kTRUE;
}

//____________________________________________________________________________
void KMCDetectorFwd::CheckTrackProlongations(KMCProbeFwd *probe, KMCLayerFwd* lrP, KMCLayerFwd* lr)
{
  // explore prolongation of probe from lrP to lr with all possible clusters of lrP
  // the probe is already brought to clusters frame
  // for the last ITS plane apply Branson correction
  //if (lrP->GetUniqueID()==fLastActiveLayerITS) probe->Print("etp");
  /*
  if (lrP->GetUniqueID()==fLastActiveLayerITS && fVtx) {
    printf("Before Branson: "); probe->Print("etp");
    double zP = probe->GetZ();
    if (!PropagateToZBxByBz(probe,fVtx->GetZ(),fDefStepAir)) return;
    double measVErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(),0,fVtx->GetXRes()*fVtx->GetXRes()}; // we work in tracking frame here!
    double measV[2] = {0,0};
    probe->Update(measV,measVErr2);
    if (!PropagateToZBxByBz(probe,zP,fDefStepAir)) return;
    printf("After Branson: "); probe->Print("etp");
  }
  */
  static KMCProbeFwd propVtx;
  //
  int nCl = lrP->GetNBgClusters();
  double rad = probe->GetR();
  double sgy = lrP->GetYRes(rad), sgx = lrP->GetXRes(rad); 
  double measErr2[3] = { sgy*sgy, 0, sgx*sgx}; 
  double meas[2] = {0,0};
  double tolerY = probe->GetTrack()->GetSigmaY2() + measErr2[0];
  double tolerX = probe->GetTrack()->GetSigmaZ2() + measErr2[2]; // Xlab = -Zloc
  tolerY = TMath::Sqrt(fMaxChi2Cl*tolerY);
  tolerX = TMath::Sqrt(fMaxChi2Cl*tolerX);
  double yMin = probe->GetY() - tolerY;
  double yMax = probe->GetY() + tolerY;    
  double xMin = probe->GetX() - tolerX;
  double xMax = probe->GetX() + tolerX;
  //probe->GetTrack()->Print();
  //
  probe->SetInnerLrChecked(lrP->GetActiveID());
  AliDebug(2,Form("From Lr(%d) %s to Lr(%d) %s | LastChecked %d",
		  lrP->GetActiveID(),lrP->GetName(),lr->GetActiveID(),lr->GetName(),probe->GetInnerLayerChecked()));
  for (int icl=-1;icl<nCl;icl++) {
    //
    if (gRandom->Rndm() > lrP->GetLayerEff()) continue; // generate layer eff
    //
    KMCClusterFwd *cl = icl<0 ? lrP->GetMCCluster() : lrP->GetBgCluster(icl);  // -1 is for true MC cluster
    if (cl->IsKilled()) {
      if (AliLog::GetGlobalDebugLevel()>1) {printf("Skip cluster %d ",icl); cl->Print();}
      continue;
    }
    double y = cl->GetY(); // ! tracking frame coordinates: Ylab
    double x = cl->GetX(); //                               XLab, sorted in decreasing order
    //
    //AliInfo(Form("Check against cl#%d(%d) out of %d at layer %s | y: Tr:%+8.4f Cl:%+8.4f (%+8.4f:%+8.4f) z: Tr:%+8.4f Cl: %+8.4f (%+8.4f:%+8.4f)",
    //		    icl,cl->GetTrID(),nCl,lrP->GetName(), probe->GetY(),y,yMin,yMax,probe->GetX(),x,xMin,xMax));
    //
    if (x>xMax) {if (icl==-1) continue; else break;} // all other x will be even smaller, no chance to match
    if (x<xMin) continue;
    if (y<yMin || y>yMax) continue; 
    //
    meas[0] = y; meas[1] = -x;
    double chi2 = probe->GetPredictedChi2(meas,measErr2);
    //
    //    AliInfo(Form("Seed-to-cluster chi2 = Chi2=%.2f for cl:",chi2));
    //      cl->Print("lc");
    //    AliDebug(2,Form("Seed-to-cluster chi2 = Chi2=%.2f",chi2));
    if (icl<0 && fHChi2LrCorr) fHChi2LrCorr->Fill(lrP->GetActiveID(), chi2);
    if (chi2>fMaxChi2Cl) continue;
    // 
    //    printf("Lr%d | cl%d, chi:%.3f X:%+.4f Y:%+.4f | x:%+.4f y:%+.4f |Sg: %.4f %.4f\n",
    //	   lrP->GetActiveID(),icl,chi2, (zMin+zMax)/2,(yMin+yMax)/2, z,y, tolerX/TMath::Sqrt(fMaxChi2Cl),tolerY/TMath::Sqrt(fMaxChi2Cl));
    // update track copy
    KMCProbeFwd* newTr = lr->AddMCTrack( probe );
    if (!newTr->Update(meas,measErr2)) {
      AliDebug(2,Form("Layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}",
		      lrP->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
      if (AliLog::GetGlobalDebugLevel()>1) newTr->Print("l");
      newTr->Kill();
      lr->GetMCTracks()->RemoveLast();
      continue;
    }
    if (fMinP2Propagate>0) {
      double p = newTr->GetTrack()->GetP();
      if (p<fMinP2Propagate) {
	newTr->Kill();
	lr->GetMCTracks()->RemoveLast();
	continue;
      }
    }
    newTr->AddHit(lrP, chi2, cl->GetTrID());

    //////////////////// check chi2 to vertex
    if (fVtx && !fVtx->IsDead() && fMaxChi2Vtx>0) {
      double measVErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(),0,fVtx->GetXRes()*fVtx->GetXRes()}; // we work in tracking frame here!
      propVtx = *newTr;

      if (!PropagateToZBxByBz(&propVtx,fRefVtx[2],fDefStepAir)) {/*printf("???????? Kill\n");*/ newTr->Kill(); lr->GetMCTracks()->RemoveLast();}
      double chi2V = propVtx.GetTrack()->GetPredictedChi2(fRefVtx,measVErr2);
      if (fHChi2VtxCorr && fHChi2VtxFake) {
	if (IsCorrect(newTr)) fHChi2VtxCorr->Fill(newTr->GetNITSHits(),chi2V);
	else                  fHChi2VtxFake->Fill(newTr->GetNITSHits(),chi2V);
      }
      AliDebug(2,Form("Chi2 to vertex: %f | y: Tr:%+8.4f Cl:%+8.4f  z: Tr:%+8.4f Cl: %+8.4f",chi2V,
		      propVtx.GetY(),fRefVtx[0],
		      propVtx.GetX(),fRefVtx[1]));

      if (chi2V>fMaxChi2Vtx) {
	newTr->Kill();
	lr->GetMCTracks()->RemoveLast();
	continue;
      }

      /*
      double pz = 1./TMath::Abs(propVtx.GetTrack()->Get1P());
      if (pz<1) {
	printf("LowMom: %f | chi2V=%f (%f | %f %f %f)  ",pz,chi2V, fMaxChi2Vtx, fRefVtx[0],fRefVtx[1],fRefVtx[2]);  propVtx.Print("etp");       
      }
      */

    }

    ////////////////////////////////////////

    //    if (!PropagateToLayer(newTr,lrP,lr,-1)) {newTr->Kill(); continue;} // propagate to next layer
    if (AliLog::GetGlobalDebugLevel()>1) {
      AliInfo("Cloned updated track is:");
      newTr->Print();
    }
  }
  //
}

//_________________________________________________________
//Double_t KMCDetectorFwd::HitDensity(double xLab,double ylab,double zlab ) const
Double_t KMCDetectorFwd::HitDensity(double ,double ,double  ) const
{
  // RS to do
  return 1;
}

//_________________________________________________________
Bool_t KMCDetectorFwd::NeedToKill(KMCProbeFwd* probe) const
{
  // check if the seed should be killed
  const Bool_t kModeKillMiss = kFALSE;
  Bool_t kill = kFALSE;
  int decVtx = (fVtx && !fVtx->IsDead());
  while (1) {
    int il = probe->GetInnerLayerChecked();
    int nITS = probe->GetNITSHits();
    int nITSMax = nITS + il - decVtx; // maximum it can have
    if (nITSMax<fMinITSHits) {
      kill = kTRUE; 
      break;
    }    // has no chance to collect enough ITS hits
    //
    int ngr = fPattITS.GetSize();
    if (ngr>0) { // check pattern
      UInt_t patt = probe->GetHitsPatt();
      // complete the layers not checked yet
      for (int i=il;i--;) patt |= (0x1<<i);
      for (int ig=ngr;ig--;) 
	if (!(((UInt_t)fPattITS[ig]) & patt)) {
	  kill = kTRUE; 
	  break;
	}
      //
    }
    //
    if (nITS>2) {  // check if smallest possible norm chi2/ndf is acceptable
      double chi2min = probe->GetChi2ITS();
      if (kModeKillMiss) {
	int nMiss = fNActiveLayersITS - (probe->GetInnerLayerChecked() - decVtx) - nITS; // layers already missed
	chi2min = nMiss*probe->GetMissingHitPenalty();
      }
      chi2min /= ((nITSMax<<1)-KMCProbeFwd::kNDOF);
      if (chi2min>fMaxNormChi2NDF) {
	kill = kTRUE; 
	break;
      }
    }
    //
    /*
    // loose vertex constraint
    double dst;
    if (nITS>=2) {
      probe->GetZAt(0,fBFieldG,dst);
      //printf("Zd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    if (nITS>=3) {
      probe->GetYAt(0,fBFieldG,dst);
      //printf("Dd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    */
    //
    break;
  }
  if (kill && AliLog::GetGlobalDebugLevel()>1 && probe->GetNFakeITSHits()==0) {
    printf("Killing good seed, last upd layer was %d\n",probe->GetInnerLayerChecked());
    probe->Print("l");
  }
  return kill;  
  //
}

//_________________________________________________________
void KMCDetectorFwd::PerformDecay(KMCProbeFwd* trc)
{
  // Decay track
  if (fDecMode==kNoDecay) return;
  //  printf("DecMode: %d\n",fDecMode);
  static TGenPhaseSpace decay;
  static TLorentzVector pDecParent,pDecMu,pParCM;
  //
  const double kTol = 5e-3; // 5 MeV tolerance
  double mass = trc->GetMass();  
  double pxyz[3],xyz[3]={0};
  trc->GetPXYZ(pxyz);
  AliDebug(2,Form(" >>Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],mass));
  static double ctau = 1e10;
  if (fDecMode==kDoRealDecay) {
    //
    if (TMath::Abs(mass-kMassMu)<kTol) {AliInfo(Form("Decay requested but provided mass %.4f hints to muon",mass)); exit(1);}
    if (TMath::Abs(mass-kMassPi)<kTol) {
      mass = kMassPi;
      Double_t masses[2] = {kMassMu, kMassE};
      pParCM.SetXYZM(0,0,0,mass);
      decay.SetDecay(pParCM, 2, masses);
      ctau = 780.45;
    }
    else if (TMath::Abs(mass-kMassK)<kTol) {
      mass = kMassK;
      Double_t masses[3] = {kMassMu, 0};
      pParCM.SetXYZM(0,0,0,mass);
      decay.SetDecay(pParCM, 2, masses);
      ctau = 371.2;
    }
    else {AliInfo(Form("Decay requested but provided mass %.4f is not recognized as pi or K",mass)); exit(1);}
    //
    decay.Generate();
    pDecMu = *decay.GetDecay(0); // muon kinematics in parent frame
    fDecMode = kApplyDecay;
  }
  //
  pDecParent.SetXYZM(pxyz[0],pxyz[1],pxyz[2],pParCM.M());
  pDecMu = *decay.GetDecay(0);
  pDecMu.Boost(pDecParent.BoostVector());
  //
  pxyz[0] = pDecMu.Px();
  pxyz[1] = pDecMu.Py();
  pxyz[2] = pDecMu.Pz();
  //
  static KMCProbeFwd tmpPr;
  tmpPr.Init(xyz,pxyz,trc->GetTrack()->Charge());
  double* parTr  = (double*)trc->GetTrack()->GetParameter();
  double* parNew = (double*)tmpPr.GetTrack()->GetParameter();  
  for (int i=0;i<5;i++) parTr[i] = parNew[i];
  trc->SetMass(kMassMu);
  // 
  // set decay weight
  trc->GetXYZ(xyz);
  for (int i=3;i--;) xyz[i]-=fGenPnt[i];
  double dst = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  double ctgamma = ctau*pDecParent.Gamma();
  double exparg = dst/ctgamma;
  //  double wgh =  exparg<100 ? 1.-TMath::Exp(-exparg) : 1.0;
  double wgh =  exparg<100 ? TMath::Exp(-exparg) : 0;
  // account for the losses due to the hadron interactions >>>
  double wabs = 0;
  int nb = fDecayZProf->GetNbinsX();
  for (int ib=1;ib<nb;ib++) {
    double x = fDecayZProf->GetBinCenter(ib);
    exparg = x/ctgamma;
    if (exparg>100) break;
    double wb = fDecayZProf->GetBinContent(ib);
    wabs += TMath::Exp(-exparg)*wb*fDecayZProf->GetBinWidth(ib);
    if (wb<1e-9) break;
  }
  wgh *= wabs/ctgamma;
  // account for the losses due to the hadron interactions <<<
  trc->SetWeight(wgh);
  //  printf("Decay %.3f Z:%+7.3f ctau: %f gamma %f wgh:%e\n",mass,fZDecay,ctau, pDecMu.Gamma(),wgh);
  //
  AliDebug(2,Form(" <<Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f, Wgh:%.3e",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],kMassMu,wgh));
  //AliInfo(Form(" <<Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f, Wgh:%.3e",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],kMassMu,wgh));
  //
}

//_________________________________________________________
void KMCDetectorFwd::InitDecayZHisto(double absorberLambda)
{
  // prepare a profile of Zdecay: uniform till start of the absorber, then exponentially dumped
  if (absorberLambda<1) absorberLambda = 1;
  KMCLayerFwd* labs = 0;
  for (int i=0;i<fNLayers;i++) {
    if ( (labs=GetLayer(i))->IsAbs()) break;
    labs = 0;
  }
  if (!labs) {
    AliError("Could not find beginning of the absorber, is setup loaded?");
    exit(1);
  }
  double zdmp = labs->GetZ()-labs->GetThickness()/2;
  double zmax = GetLayer(fLastActiveLayer)->GetZ();
  AliInfo(Form("Decay will be done uniformly till Z=%.1f, then dumped with Lambda=%.1f cm",zdmp,absorberLambda));
  //
  int nbn = int(zmax-zdmp+1);
  TH1F* hd = new TH1F("DecayZProf","Z decay profile",nbn,0,zmax);
  for (int i=1;i<=nbn;i++) {
    double z = hd->GetBinCenter(i);
    if (z<zdmp) hd->SetBinContent(i,1.);
    else {
      double arg = (z-zdmp)/absorberLambda;
      hd->SetBinContent(i, arg>100 ? 0 : TMath::Exp(-arg));
    }
  }
  SetDecayZProfile(hd);
  //
}

//_________________________________________________________
void KMCDetectorFwd::GenBgEvent(double x, double y, double z, int offset)
{
  if (fNChPi<0 && fNChK<0 && fNChP<0) return;
  // generate bg. events from simple thermalPt-gaussian Y parameterization
  if (!fdNdYPi || !fdNdYK || !fdNdYP || !fdNdPtPi || !fdNdPtK || !fdNdPtP) AliFatal("Background generation was not initialized");
  //
  int maxLr = fLastActiveLayerITS + offset;
  if (maxLr > fLastActiveLayer) maxLr = fLastActiveLayer;
  //
  for (int ilr=fLastActiveLayer;ilr--;) {
    KMCLayerFwd* lr = GetLayer(ilr);
    if (lr->IsDead()) continue;
    lr->ResetBgClusters();
    if (lr->InheritsFrom(KMCMSStation::Class())) {
      ((KMCMSStation*)lr)->ClearPrimaryHits();
    }
  }
  int decMode = fDecMode;
  fDecMode = kNoDecay;
  KMCProbeFwd bgtr;
  //
  //  double ntr = gRandom->Poisson( fNCh );
//   for (int itr=0;itr<ntr;itr++) {
//     double yrap  = fdNdY->GetRandom();
//     double pt = fdNdPt->GetRandom();
//     double phi = gRandom->Rndm()*TMath::Pi()*2;
//     int charge = gRandom->Rndm()>0.5 ? 1:-1;
//     CreateProbe(&bgtr, pt, yrap, phi, kMassPi, charge, x,y,z);
//     bgtr.SetTrID(itr);
//     TransportKalmanTrackWithMS(&bgtr, maxLr,kTRUE);
//   }
  int ntrTot = 0;
  // pions
  double ntrPi = gRandom->Poisson( fNChPi );
  printf("fNChPi=%f ntrPi=%f\n",fNChPi,ntrPi);
  for (int itr=0;itr<ntrPi;itr++) {
    double yrap  = fdNdYPi->GetRandom();
    double pt = fdNdPtPi->GetRandom();
    double phi = gRandom->Rndm()*TMath::Pi()*2;
    int charge = gRandom->Rndm()>0.52 ? 1:-1;
    CreateProbe(&bgtr, pt, yrap, phi, kMassPi, charge, x,y,z);
    bgtr.SetTrID(ntrTot++);
    TransportKalmanTrackWithMS(&bgtr, maxLr,kTRUE);
  }
  // kaons
  double ntrK = gRandom->Poisson( fNChK );
  printf("fNChK=%f ntrK=%f\n",fNChK,ntrK);
  for (int itr=0;itr<ntrK;itr++) {
    double yrap  = fdNdYK->GetRandom();
    double pt = fdNdPtK->GetRandom();
    double phi = gRandom->Rndm()*TMath::Pi()*2;
    int charge = gRandom->Rndm()>0.3 ? 1:-1;
    CreateProbe(&bgtr, pt, yrap, phi, kMassK, charge, x,y,z);
    bgtr.SetTrID(ntrTot++);
    TransportKalmanTrackWithMS(&bgtr, maxLr,kTRUE);
  }
  // protons
  double ntrP = gRandom->Poisson( fNChP );
  printf("fNChP=%f ntrP=%f\n",fNChP,ntrP);
  for (int itr=0;itr<ntrP;itr++) {
    double yrap  = fdNdYP->GetRandom();
    double pt = fdNdPtP->GetRandom();
    double phi = gRandom->Rndm()*TMath::Pi()*2;
    int charge = 1;
    CreateProbe(&bgtr, pt, yrap, phi, kMassP, charge, x,y,z);
    bgtr.SetTrID(ntrTot++);
    TransportKalmanTrackWithMS(&bgtr, maxLr,kTRUE);
  }
  //
  for (int ilr=maxLr;ilr--;) {
    KMCLayerFwd* lr = GetLayer(ilr);
    if (lr->IsDead()) continue;
    lr->SortBGClusters();
  }
  fDecMode = decMode;
  //  
}

//_________________________________________________________
void KMCDetectorFwd::InitBgGeneration(int dndeta, 
				      double y0, double sigy, double ymin,double ymax,
				      double T, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNCh = dndeta*0.5*(TMath::Erf((ymax-y0)/sqrt(2.)/sigy)-TMath::Erf((ymin-y0)/sqrt(2.)/sigy));
  fdNdY = new TF1("dndy","exp( -0.5*pow( (x-[0])/[1],2) )",ymin,ymax);
  fdNdY->SetParameters(y0,sigy);
  //
  fdNdPt = new TF1("dndpt","x*exp(-sqrt(x*x+1.949e-02)/[0])",ptmin,ptmax); // assume pion
  fdNdPt->SetParameter(0,T);
}

//_________________________________________________________
void KMCDetectorFwd::InitBgGenerationPart(double NPi,double NKplus,double NKminus ,double NP,double Piratio,
					  double y0, double y0Pi,double y0Kplus,double y0Kminus,double y0P,
					  double sigyPi, double sigyKplus,double sigyKminus, double sigyP,
					  double ymin,double ymax,
					  double Tpi, double TK, double TP, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNChPi = NPi*(1+Piratio)*sigyPi*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Pi)/sqrt(2.)/sigyPi)+ TMath::Erf((ymax-y0+y0Pi)/sqrt(2.)/sigyPi)-TMath::Erf((ymin-y0-y0Pi)/sqrt(2.)/sigyPi)-TMath::Erf((ymin-y0+y0Pi)/sqrt(2.)/sigyPi));
  fNChK = NKplus*sigyKplus*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Kplus)/sqrt(2.)/sigyKplus)+ TMath::Erf((ymax-y0+y0Kplus)/sqrt(2.)/sigyKplus)-TMath::Erf((ymin-y0-y0Kplus)/sqrt(2.)/sigyKplus)-TMath::Erf((ymin-y0+y0Kplus)/sqrt(2.)/sigyKplus))+
    NKminus*sigyKminus*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Kminus)/sqrt(2.)/sigyKminus)+ TMath::Erf((ymax-y0+y0Kminus)/sqrt(2.)/sigyKminus)-TMath::Erf((ymin-y0-y0Kminus)/sqrt(2.)/sigyKminus)-TMath::Erf((ymin-y0+y0Kminus)/sqrt(2.)/sigyKminus));
  fNChP = NP*sigyP*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0P)/sqrt(2.)/sigyP)+ TMath::Erf((ymax-y0+y0P)/sqrt(2.)/sigyP)-TMath::Erf((ymin-y0-y0P)/sqrt(2.)/sigyP)-TMath::Erf((ymin-y0+y0P)/sqrt(2.)/sigyP));
  //
  fdNdYPi = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYPi->SetParameters(y0,y0Pi,sigyPi);
  fdNdYK  = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYK ->SetParameters(y0,y0Kplus ,sigyKplus);
  fdNdYP  = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYP ->SetParameters(y0,y0P,sigyP);
  //
  fdNdPtPi = new TF1("dndptPi","x*exp(-sqrt(x*x+1.949e-02)/[0])",ptmin,ptmax); // pion
  fdNdPtK  = new TF1("dndptK" ,"x*exp(-sqrt(x*x+0.493*0.493)/[0])",ptmin,ptmax); // kaon
  fdNdPtP  = new TF1("dndptP" ,"x*exp(-sqrt(x*x+0.938*0.938)/[0])",ptmin,ptmax); // proton
  fdNdPtPi->SetParameter(0,Tpi);
  fdNdPtK->SetParameter(0,TK);
  fdNdPtP->SetParameter(0,TP);
}

//_____________________________________________________________________
void KMCDetectorFwd::RequirePattern(UInt_t patt)
{
  // optional pattern to satyisfy
  if (!patt) return;
  int ngr = fPattITS.GetSize();
  fPattITS.Set(ngr+1);
  fPattITS[ngr] = patt;
}

//_____________________________________________________________________
void KMCDetectorFwd::BookControlHistos()
{
  fHChi2Branson = new TH1F("chi2Branson","chi2 Mu @ vtx",      100,0,100);
  fHChi2LrCorr = new TH2F("chi2Cl","chi2 corr cluster",        fNActiveLayersITS+1,0,fNActiveLayersITS+1,100,0,fMaxChi2Cl);
  fHChi2NDFCorr = new TH2F("chi2NDFCorr","chi2/ndf corr tr.",  fNActiveLayersITS+1,0,fNActiveLayersITS+1,100,0,fMaxNormChi2NDF);
  fHChi2NDFFake = new TH2F("chi2NDFFake","chi2/ndf fake tr.",  fNActiveLayersITS+1,0,fNActiveLayersITS+1,100,0,fMaxNormChi2NDF);
  fHChi2VtxCorr = new TH2F("chi2VCorr","chi2 to VTX corr tr." ,fNActiveLayersITS+1,0,fNActiveLayersITS+1,100,0,100);
  fHChi2VtxFake = new TH2F("chi2VFake","chi2 to VTX fake tr." ,fNActiveLayersITS+1,0,fNActiveLayersITS+1,100,0,100);
  fHNCand     = new TH2F("hNCand","Ncand per layer",           fNActiveLayersITS+1,0,fNActiveLayersITS+1,200,0,-1);
  fHCandCorID = new TH2F("CandCorID","Corr.cand ID per layer", fNActiveLayersITS+1,0,fNActiveLayersITS+1,200,0,-1);
  //
  fHChi2MS = new TH2F("chi2ms","chi2ms",100,0,30,10,0,10);
  //
}
//_________________________________________________________
void KMCDetectorFwd::InitBkg(double beamenergy){

  int E=TMath::Nint(beamenergy);

  // default values (from 40 GeV)
  double y0BG = 2.22; // gaussian y mean - 40 GeV
  double y0BGPi = 0.666;
  double y0BGKplus = 0.694;
  double y0BGKminus = 0.569;
  double y0BGP = 0.907;
  double sigyBG = 1.2; // .. sigma
  double sigyBGPi = 0.872;
  double sigyBGKplus = 0.725;
  double sigyBGKminus = 0.635;
  double sigyBGP = 0.798;
  double yminBG = 1.5; // min y to generate
  double ymaxBG = 4.5; //
  double TBG = 0.17;   // inv.slope of thermal pt distribution
  double TBGpi = 0.17;
  double TBGK = 0.23;
  double TBGP = 0.26;
  double ptminBG = 0.01;
  double ptmaxBG = 5;
  double dndyBGPi = 615.;
  double dndyBGK = 78.;
  double dndyBGP = 150.;
  double NBGPi = 74.;
  double NBGKplus = 16.2;
  double NBGKminus = 6.03;
  double NBGP = 37.5;
  double Piratio = 0.91;

  if (E == 20){
    printf("--- Background parameters for E=20 GeV/nucleon ---\n");
    y0BG = 1.9;   // gaussian y mean - 40 GeV
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 410.;
    dndyBGK = 51.;
    dndyBGP = 148.;
    TBGpi = 0.17;
    TBGK = 0.22;
    TBGP = 0.26;
  }else if (E == 40){ 
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=40 GeV/nucleon ---\n");
    y0BG = 2.22; // gaussian y mean - 40 GeV
    y0BGPi = 0.666;
    y0BGKplus = 0.694;
    y0BGKminus = 0.569;
    y0BGP = 0.907;
    sigyBG = 1.2; // .. sigma
    sigyBGPi = 0.872;
    sigyBGKplus = 0.725;
    sigyBGKminus = 0.635;
    sigyBGP = 0.798;
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    NBGPi = 74.;
    NBGKplus = 16.2;
    NBGKminus = 6.03;
    NBGP = 37.5;
    Piratio = 0.91;
  }else if (E == 60){ 
    // average of values at 40 and 80
    printf("--- Background parameters for E=60 GeV/nucleon ---\n");
    y0BG = 2.42;  // gaussian y mean - 60 GeV
    y0BGPi = 0.5*(0.666+0.756);
    y0BGKplus = 0.5*(0.694+0.742);
    y0BGKminus = 0.5*(0.569+0.668);
    y0BGP = 0.5*(0.907+0.907);
    sigyBG = 1.2; // .. sigma
    sigyBGPi = 0.5*(0.872+0.974);
    sigyBGKplus = 0.5*(0.725+0.792);
    sigyBGKminus = 0.5*(0.635+0.705);
    sigyBGP = 0.5*(0.798+0.798);
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi =  0.5*(615.+920.);
    dndyBGK =  0.5*(78.+109.);
    dndyBGP =  0.5*(150.+(30.1/41.3)*150.);
    NBGPi = 0.5*(74.+97.);
    NBGKplus = 0.5*(16.2+19.3);
    NBGKminus = 0.5*(6.03+9.16);
    NBGP = 0.5*(37.5+(30.1/41.3)*37.5);
    Piratio = 0.93;
  }else if (E == 80){
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=80 GeV/nucleon ---\n");
    y0BG = 2.57;  // gaussian y mean - 80 GeV
    y0BGPi = 0.756;
    y0BGKplus = 0.742;
    y0BGKminus = 0.668;
    y0BGP = 0.907;
    sigyBGPi = 0.974;
    sigyBGKplus = 0.792;
    sigyBGKminus = 0.705;
    sigyBGP = 0.798;
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.18;   // inv.slope of thermal pt distribution
    TBGpi = 0.18;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 920.;
    dndyBGK = 109.;
    dndyBGP = (30.1/41.3)*150.;
    NBGPi = 97.;
    NBGKplus = 19.3;
    NBGKminus = 9.16;
    NBGP = (30.1/41.3)*37.5; //ratio 80/40 from PRC73, 044910 (2006)
    Piratio = 0.94;
  }else if (E == 160){
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=160 GeV/nucleon ---\n");
    y0BG = 2.9; // gaussian y mean - 160 GeV
    y0BGPi = 0.72;
    y0BGKplus = 0.839;
    y0BGKminus = 0.727;
    y0BGP = 39.8;
    sigyBGPi = 1.18;
    sigyBGKplus = 0.88;
    sigyBGKminus = 0.81;
    sigyBGP = 8.07;
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 1258.;
    dndyBGK = 155.;
    dndyBGP = 292.;
    NBGPi = 107.6;
    NBGKplus = 23.4;
    NBGKminus = 12.8;
    NBGP = 2.55e+06;
    Piratio = 0.97;
    TBGpi = 0.18; // inv.slope of thermal pt distribution
    TBGK = 0.23;  // inv.slope of thermal pt distribution
    TBGP = 0.31;  // inv.slope of thermal pt distribution
  }else if (E == 400){
    y0BG = 3.37;  // gaussian y mean - 80 GeV
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
  }else{
    printf("--- Parameters not available at this energy, use those for E=40 GeV/nucleon ---\n");
  }
  
  printf("Simulation of background at %d GeV/nucleon\n", E);
  printf("pions:   total multiplicity = %f ; T = %f\n", dndyBGPi, TBGpi);
  printf("kaons:   total multiplicity = %f ; T = %f\n", dndyBGK, TBGK);
  printf("protons: total multiplicity = %f ; T = %f\n", dndyBGP, TBGP);

  InitBgGenerationPart(NBGPi, NBGKplus, NBGKminus, NBGP, Piratio, y0BG, y0BGPi, y0BGKplus, y0BGKminus, y0BGP, sigyBGPi, sigyBGKplus, sigyBGKminus, sigyBGP, yminBG, ymaxBG, TBGpi, TBGK, TBGP, ptminBG, ptmaxBG);
  return;
}

//_____________________________________________________________________
Bool_t KMCDetectorFwd::IsCorrect(KMCProbeFwd *probTr)
{
  if (probTr->GetNFakeITSHits()) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________
void  KMCDetectorFwd::SetMinITSHits(int n)
{
  fMinITSHits = TMath::Min(n,fNActiveLayersITS);
}

//_____________________________________________________________________
void  KMCDetectorFwd::SetMinMSHits(int n)
{
  fMinMSHits = TMath::Min(n,fNActiveLayersMS);
}
//_____________________________________________________________________
void  KMCDetectorFwd::SetMinTRHits(int n)
{
  fMinTRHits = TMath::Min(n,fNActiveLayersTR);
}

//_____________________________________________________________________
void  KMCDetectorFwd::ForceLastActiveLayer(int lr)
{
  printf("Attention: overriding last active layer from %d to %d\n",fLastActiveLayer,lr);
  fLastActiveLayer = lr;
}

//_____________________________________________________________________
void KMCDetectorFwd::PrepareForTracking()
{
  for (Int_t j=fNLayers; j--;) {  // Layer loop
    KMCLayerFwd* lr = GetLayer(j);
    if (lr->IsDead()) continue;
    lr->PrepareForTracking();
  }
}

//_____________________________________________________________________
KMCLayerFwd* KMCDetectorFwd::GetActiveLayer(int i, int ltype) const
{
  if (ltype == KMCLayerFwd::kITS) return (KMCLayerFwd*)fLayersITS[i];
  if (ltype == KMCLayerFwd::kMS) return (KMCLayerFwd*)fLayersMS[i];
  if (ltype == KMCLayerFwd::kTRIG) return (KMCLayerFwd*)fLayersTR[i];
  return 0;
}

//_____________________________________________________________________
bool KMCDetectorFwd::ImposeFlukaBackground(KMCFlukaParser* fp, const TString& interactionSource, bool allowRewind)
{
  if (!fp->GetNextBackgroundEvent(interactionSource, allowRewind)) return false;
  const std::vector<FlukaHit>& hits = fp->GetHits();
  for (int ilr=fLastActiveLayer;ilr--;) {
    KMCLayerFwd* lr = GetLayer(ilr);
    lr->ResetBgClusters();
    if (lr->InheritsFrom(KMCMSStation::Class())) {
      ((KMCMSStation*)lr)->ClearPrimaryHits();
    }
  }
  int ic = 0;
  for (const auto& hit : hits) {
    KMCLayerFwd* lr = GetActiveLayer(hit.stationID, hit.stationType);
    if (!lr) {
      printf("ERROR: hit in non-existing layer %d of type %d, check setup\n",hit.stationID, hit.stationType);
      exit(1);
    }
    lr->AddCluster(hit.recData[kX],hit.recData[kY],hit.recData[kZ], ic++, 1);
  }
  printf("\nAdded Fluka background event of %d hits\n",ic);
  return true;
}

KMCProbeFwd* KMCDetectorFwd::CreateMSSeed(KMCClusterFwd* cl0, KMCClusterFwd* cl1)
{
  // try to create MS track seed from clusters ic0 and ic1 of the 1st and last layer after the magnet
  double pos0[3] = {cl0->GetXLab(), cl0->GetYLab(), cl0->GetZLab()};
  double pos1[3] = {cl1->GetXLab(), cl1->GetYLab(), cl1->GetZLab()};
  double dirc[3], dirf[3], posBend[3], dirnrm = 0., dirnrmF = 0., dzbend = fZBendingMS - pos0[2];
  for (int i=0;i<3;i++) {
    dirc[i] = pos1[i] - pos0[i];
    dirnrm += dirc[i]*dirc[i];
  }
  // check if the segment does not cross the normal from the origin ("too much" bending)
  double tcross = - (dirc[0]*pos0[0] + dirc[1]*pos0[1])/(dirnrm - dirc[2]*dirc[2]);
  if (tcross>0 && tcross<1) {
    printf("reject as crossing origin\n");
    for (int i=0;i<3;i++) printf("%+e ", pos0[i]); printf("\n");
    for (int i=0;i<3;i++) printf("%+e ", pos1[i]); printf("\n");
    return 0;
  }
  
  dirnrm = 1./TMath::Sqrt(dirnrm);
  // bending angle of line from the origin to bending point and from bending point to clusters
  for (int i=0;i<3;i++) {
    dirc[i] *= dirnrm;
    posBend[i] = pos0[i] + dirc[i]*dzbend;
    dirnrmF += posBend[i]*posBend[i];
  }
  dirnrmF = 1./TMath::Sqrt(dirnrmF);
  double theta0 = TMath::Sqrt(posBend[0]*posBend[0]+posBend[1]*posBend[1])/fZBendingMS;
  double theta1 = TMath::Sqrt(dirc[0]*dirc[0]+dirc[1]*dirc[1])/dirc[2];
  double dTheta = TMath::ATan(theta0) - TMath::ATan(theta1);  //  dTheta = 0.3e-3* B0 / pT * ln(z2/z1) == bending angle, B in kgaus
  double pTQF = TMath::Abs(dTheta)>1e-4 ? 0.3e-3*fToroidB0*TMath::Log(fZToroidEnd/fZToroidStart)/dTheta : 1e3; // estimate of q*pT at target
  double ptot = TMath::Abs(pTQF) / TMath::Sin(theta0); // full momentum
  //printf("Theta : %e %e dTheta: %e | qpt : %+e ptot : %+e\n", theta0, theta1, dTheta, pTQF, ptot);
  double pT = ptot*TMath::Sin(theta1), pz = ptot*TMath::Cos(theta1);
  double phi = TMath::ATan2(dirc[1],dirc[0]);
  double pxyz[3] = {pT*TMath::Cos(phi), pT*TMath::Sin(phi), ptot*TMath::Cos(theta1)};

  KMCProbeFwd* seed = new KMCProbeFwd(pos1, pxyz, pTQF>0 ? 1:-1);
  printf("Seed : "); seed->Print("etp");
  printf("True : "); GetLayerMS(fNActiveLayersMS-1)->GetAnProbe()->Print("etp");
  return seed;
}


//=======================================================================================


//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================


NaCardsInput::NaCardsInput() 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
  fLastComment = "";
}
//--------------------------------------------------
//
NaCardsInput::NaCardsInput(const char* fname) 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
  fCardsFile = OpenFile(fname);
  fLastComment = "";
}
//--------------------------------------------------
//
NaCardsInput::~NaCardsInput() 
{
  if (fCardsFile) fclose(fCardsFile);
  ClearArgs();
}
//--------------------------------------------------
//
FILE* NaCardsInput::OpenFile(const char *fname) 
{
  TString flstr = fname;
  gSystem->ExpandPathName(flstr);
  if (fCardsFile) fclose(fCardsFile);
  if ( !(fCardsFile = fopen(flstr.Data(),"r")) ) { 
    printf("Error in <NaCardsInput::OpenFile>: Did not find file %s\n",flstr.Data());
    return NULL;
  }  
  ClearArgs();
  return fCardsFile;
}
//--------------------------------------------------
//
void NaCardsInput::ClearArgs() 
{
  if (fArgs) {
    for (int i=0;i<fNArgs;i++) delete[] fArgs[i];
    delete[] fArgs;
    fArgs = 0;
  }
  if (fKey) {
    delete[] fKey;
    fKey = 0;
  }
  if (fModifier) {
    delete[] fModifier;
    fModifier = 0;
  }
  fNArgs = 0;
  fLastComment = "";
}
//--------------------------------------------------
//
int NaCardsInput::FindEntry(const char* ckey, const char* cmod,const char* form, int rew, int errmsg) 
{
  // Find and fill arguments list from the line starting with KEY:MODIFIER ...
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'b' - argument should be string of 0 and 1, interpreted as bit mask
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If rew != 0, the search is done from the beginning of the file
  // If errmsg is not 0, then the warning is printed in case KEY and MODIFIER were matched
  // but format - not
  // Case insensitive
  // If not found, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<NaCardsInput::FindEntry>  KEY is empty\n");
    return -1;
  }
  if (rew) rewind(fCardsFile);
  //
  while ( (narg=NextEntry()) > -1 ) {
    if ( CompareKey(ckey) ) continue; // KEY does not match
    if ( cmod && (cmod[0]=='\0' || cmod[0]==' ' || cmod[0]=='\t') 
	 && fModifier[0]!='\0' ) continue; // Empty modifier was requested
    if ( cmod && CompareModifier(cmod) ) continue; // MODIFIER does not match
    //
    if ( CompareArgList(form) ) {  // ArgList form does not match
      if (cmod && errmsg) {
	printf("<NaCardsInput::FindEntry> arguments of %s\n %s to \"%s\"\n",
	       fPoint,"do not match to format",form);
	//return -10;
      }
      continue;
    }
    return narg;
  }
  return -1;
}
//--------------------------------------------------
//
int NaCardsInput::NextEntry(const char* ckey, const char* cmod,const char* form)
{
  // Reads next entry in required format
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // Case insensitive
  // If next entry does not match, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<NaCardsInput::NextEntry(...)>  KEY is empty\n");
    return -1;
  }
  //
  narg=NextEntry();
  if (narg<0) return narg;
  if ( CompareKey(ckey) ) return -1; // KEY does not match
  if (!(cmod && (cmod[0]=='\0'||cmod[0]==' '||cmod[0]=='\t')&&fModifier[0]!='\0')) // NonEmpty modifier requested
    if ( cmod && CompareModifier(cmod) ) return -1; // MODIFIER does not match
    //
  if ( CompareArgList(form) ) {  // ArgList form does not match
    printf("<NaCardsInput::NextEntry> arguments of %s\n %s to \"%s\"\n",
	   fPoint,"do not match to format",form);
    return -1;
  }
  return narg;
  //
}
//--------------------------------------------------
//
int NaCardsInput::NextEntry() 
{
  // get next entry, skipping commented lines
  //
  int len;
  char *tmp1,*tmp2;
  char buftmp[fgkMaxLen]; 
  //
  if (!fCardsFile) {
    printf("<NaCardsInput::NextEntry> - No file was opened\n");
    return -1;
  }
  ClearArgs(); 
  while ( 1 ) {
    fLastPos = ftell(fCardsFile);
    if ( !(fPoint=fgets(fbuffer,fgkMaxLen,fCardsFile)) ) break; // end of file reached
    while(*fPoint == ' ' || *fPoint == '\t') fPoint++; // skip spaces in the beginning
    if (*fPoint == fgkComment) {fLastComment+=fbuffer ;continue;} // check if the line is commented
    //
    // Check if there is a line continuation flag
    tmp2 = fPoint;
    do {
      tmp1 = &fPoint[strlen(fPoint)-1];
      if (*tmp1=='\n') {*tmp1--='\0';}
      while(*tmp1 == ' ' || *tmp1 == '\t') *tmp1--='\0'; // skip spaces in the end
      if (*tmp1 == fgkContinuation) {
	if ( !(tmp2=fgets(tmp1,fgkMaxLen-(tmp1-fPoint),fCardsFile)) ) break; // end of file reached
      }
      else
	break;
    } while(1);
    if (!tmp2) break; // end of file reached
    //
    // Check if there is KEY:[MODIFIER]
    tmp1 = fPoint;
    if ( sscanf(tmp1,"%s",buftmp)<=0 ) continue; // Empty line
    if ( (tmp2=strchr(buftmp,fgkDelimiter)) ) { // there is a key/delimiter
      len = (tmp2-buftmp)/sizeof(char);
      fKey = new char[len+1];
      for (int i=0;i<len;i++) fKey[i] = tolower(buftmp[i]);
      fKey[len]='\0';
      //
      len = strlen(++tmp2); // skip delimiter, get modifier length
      fModifier = new char[len+1];
      for (int i=0;i<len;i++) fModifier[i] = tolower(tmp2[i]);
      fModifier[len] = '\0';
      // skip key:mod
      tmp1 += strlen(buftmp);
      while(*tmp1 == ' ' || *tmp1 == '\t') tmp1++; // skip spaces
    }
    else {
      continue; // no delimiter - not valid record
    }
    // First, throw away everything after comment (if any)
    if ( (tmp2=strchr(tmp1,fgkComment)) ) *tmp2 = '\0';
    // now, scan the string to get the number of arguments
    tmp2 = tmp1;
    while( sscanf(tmp2,"%s",buftmp)!= -1 ) {
      while( *tmp2 == ' ' || *tmp2 == '\t') tmp2++;
      tmp2 += strlen(buftmp); 
      fNArgs++;
    }
    // Fill the arguments list
    fArgs = new char*[fNArgs];
    for (int i=0;i<fNArgs;i++) {
      sscanf(tmp1,"%s",buftmp);
      while( *tmp1 == ' ' || *tmp1 == '\t') tmp1++;
      int lgt = strlen(buftmp);
      tmp1 += lgt;
      fArgs[i] = new char[lgt+1];
      strcpy(fArgs[i],buftmp);
    }
    return fNArgs;
  }
  return -1;
}
//--------------------------------------------------
//
char* NaCardsInput::GetArg(const int iarg, const char* option, int *err) 
{
  // Return argument iarg in character string format
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  //
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      carg[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      carg[i] = toupper(carg[i]);
  //
  return carg;
}
//--------------------------------------------------
//
char* NaCardsInput::GetArg(char* &dest, const int iarg, const char* option, int *err) 
{
  // Creates new string at dest and fill it with argument iarg 
  // in character string format.
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  // No check of dest == 0 is done!
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  dest = new char[lgt+1];
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      dest[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      dest[i] = toupper(carg[i]);
  else
    for (int i=0;i<lgt;i++) dest[i] = carg[i];
  //
  dest[lgt] = '\0';
  return dest;
}
//--------------------------------------------------
//
float NaCardsInput::GetArgF(const int iarg, int *err) 
{
  // Get requsted argument in float format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0.; }
  //
  char *errc = 0;
  float vald = float(strtod(fArgs[iarg],&errc));
  if (*errc) { *err = 1; return 0.; }
  return vald;
  //
}
//--------------------------------------------------
//
int NaCardsInput::GetArgD(const int iarg, int *err) 
{
  // Get requsted argument in integer format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  //
  char *errc = 0;
  int vald = int(strtol(fArgs[iarg],&errc,10));
  if (*errc) { *err = 1; return 0; }
  return vald;
  //
}
//--------------------------------------------------
//
unsigned int NaCardsInput::GetArgB(const int iarg, int *err) 
{
  // Get requsted argument in bitfield format, interpretting 0s and 1s from right side
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  //
  unsigned int valu = 0;
  char *str = fArgs[iarg];
  int lgt = strlen(str);
  for (int ibt=lgt;ibt--;) {
    if (str[ibt]=='1') valu |= 1<<(lgt-1-ibt);
    else if (str[ibt]!='0') {
      printf("Bit string must contain only 0 or 1: %s\n",str);
      *err = 1; 
      return 0;
    }
  }
  return valu;
  //
}
//--------------------------------------------------
//
int NaCardsInput::CompareKey(const char* key) 
{
  if (fKey) return strcasecmp(key,fKey);
  return 0;
}
//--------------------------------------------------
//
int NaCardsInput::CompareModifier(const char* mod) 
{
  if (fModifier) return strcasecmp(mod,fModifier);
  return 1;
}
//--------------------------------------------------
//
int NaCardsInput::CompareArgList(const char *form)
{
  // If form is not NULL or nonempty, check if the record matches following 
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'b' - argument should be string of 0 and 1, interpreted as bit mask
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If matches, return 0, otherwise -1
  //
  char cf=' ';
  const char *pf = form;
  char *err = 0;
  if ( !form ) return 0; // no for is requested, anything matches
  int iarg = 0;
  float valf=0.0;
  int vald = 0;
  int isoption = 0;
  //
  if (fNArgs == -1 ) {
    printf("<NaCardsInput::CompareArgList> There is no record in the buffer\n");
    return -1;
  }

  while ( (cf=tolower(*pf++)) ) {
    //
    if (cf=='*') { // the following arguments are optional
      isoption = 1;
      continue;
    }
    //
    if (cf=='?') { // the following argument block is optional
      if ( fNArgs > iarg+1 ) {
	isoption = 0; // number of following arguments should be respected
	continue;
      }
      else return 0; // no more arguments, but the block was optional
    }
    //
    if (cf=='|') { // no more arguments are alowed
      if ( fNArgs != iarg ) {printf("err1\n");return -1;} // number of arguments > than allowed
      else return 0; // no more arguments, OK
    }
    //
    if ( iarg >= fNArgs ) { // number of agruments is less than requested
      if ( isoption || cf=='|' ) return 0; // but this was allowed -> MATCHED
      else {printf("err2\n");return -1;} // too few arguments
    }
    //
    // Analize requsted argument type
    //
    //
    if (cf=='s') { // string argument is requested, anything matches
      iarg++;
      continue;
    }
    //
    if (cf=='d') { // Integer argument is requested
      vald = strtol(fArgs[iarg++],&err,10);
      if (*err) {printf("err3\n");return -1;} // Not integer, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='f') { // Float argument is requested
      valf = strtod(fArgs[iarg++],&err);
      if (*err) {printf("err4\n");return -1;} // Not float, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='a') { // string argument with non-numeric character is requested
      valf = strtod(fArgs[iarg++],&err);
      if (!(*err)) {printf("err5: |%s|\n",fArgs[iarg-1]);return -1;} // looks like number, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='b') { // 0,1 bitfield is requested
      char* str = fArgs[iarg++];
      int lgt = strlen(str);
      if (lgt>32) {printf("err6\n");return -1;}
      for (int i=0;i<lgt;i++) if (str[i]!='0'&&str[i]!='1') return -1; 
      continue; // this argument matches
    }    
    //
    // unknown format character
    printf("<NaCardsInput::CompareArgList> unknown format requsted \'%c\'\n",cf);
  }
  return 0; // Matched
}
//--------------------------------------------------
//
void NaCardsInput::Print() 
{
  if (fKey) printf("Key = %s\n",fKey);
  else printf("No Key\n");
  if (fModifier) printf("Modifier = %s\n",fModifier);
  else printf("No Modifiers\n");
  for (int i=0;i<fNArgs;i++) {
    printf("Arg %d: %s\n",i+1,fArgs[i]);
  }
}

const char NaCardsInput::fgkComment;   // comment identifier
const char NaCardsInput::fgkDelimiter; // delimiter between keyword and modifier
const char NaCardsInput::fgkContinuation; // delimiter between keyword and modifier
const int  NaCardsInput::fgkMaxLen;  // max. length of the entry


