#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <Riostream.h>
#include <TMaterial.h>
#include <TSystem.h>
#include <TString.h>
#include <stdio.h>
#include <TGeoGlobalMagField.h>

const double kVeryLarge = 1e16;

//====================================================================
class NaMaterial :public TMaterial {
public:
  enum {kNELossPar=3};
  NaMaterial();
  NaMaterial(const char* name, const char *title, 
	     Float_t a, Float_t z, Float_t dens, Float_t radl=0, 
	     Float_t absl=0, Float_t *elbuff=0);
  virtual ~NaMaterial();
  //
  virtual void Dump() const;
  virtual void Print(Option_t* option="") const;
  //
  virtual Float_t GetELoss(Float_t p) const;  
  virtual Float_t *GetELossPars() {return fELossPar;}
  virtual Float_t GetELossPar(int i) const 
    {return (i<kNELossPar && i>=0) ? fELossPar[i]:0.;}
  //
 protected:
  Float_t fELossPar[kNELossPar];  // Params of dEdX=c0+c1*p+c2*log(p) (per GeV/cm)
  //
  ClassDef(NaMaterial,1) // Material Object
};

//==========================================================================
class NaMixture :public NaMaterial {
public:
  NaMixture();
  NaMixture(const char* name, const char *title,
	    Float_t a, Float_t z, Float_t dens,
	    Float_t radl=0, Float_t absl=0, Float_t *elbuff=0);
  virtual void SetComponents(Int_t nmixt, Float_t* a,Float_t* z,Float_t* w);
  virtual ~NaMixture();
  //
  virtual void Print(Option_t* option="") const;
  //
 protected:
  Int_t fNMix;    // number of components to mix (<>0 a la Geant3 mixture)
  Float_t* fAMix; // [fNMix] A of components
  Float_t* fZMix; // [fNMix] Z ..
  Float_t* fWMix; // [fNMix] Weights ...
  //
  ClassDef(NaMixture,1) // Mixture Class
};

//==========================================================================

class NaCardsInput {
 public:
  NaCardsInput();
  NaCardsInput(const char* fname);
  virtual ~NaCardsInput();
  FILE* OpenFile(const char* fname);
  //
  virtual int FindEntry(const char* key, const char* mod="", 
			const char* form="", int rew=0, int errmsg=1);
  virtual int NextEntry(const char* key, const char* mod="",const char* form="");
  virtual int NextEntry();
  virtual int GetNArgs() {return fNArgs;}
  void  Rewind() {if (fCardsFile) rewind(fCardsFile);}
  void  StepBack() {if (fCardsFile) fseek(fCardsFile, fLastPos, SEEK_SET);}
  char* GetKey() {return fKey;}
  char* GetModifier() {return fModifier;}
  char* GetArg(const int iarg, const char* option="", int *err=0);
  char* GetArg(char* &dest, const int iarg, const char* option="", int *err=0);
  float GetArgF(const int iarg, int *err=0);
  int   GetArgD(const int iarg, int *err=0);
  unsigned int GetArgB(const int iarg, int *err=0);
  char  **GetArgs() {return fArgs;} 
  int   CompareKey(const char *key);
  int   CompareModifier(const char *mod);
  int   CompareArgList(const char *form);
  char* GetLastComment() const {return (char*)fLastComment.Data();}
  char* GetLastBuffer()  const {return (char*)fbuffer;}
  virtual void Print();
  //
 protected:
  virtual void ClearArgs(); 
 protected:
  static const char fgkComment='#';   // comment identifier
  static const char fgkDelimiter=':'; // delimiter between keyword and modifier
  static const char fgkContinuation='\\'; // delimiter between keyword and modifier
  static const int  fgkMaxLen = 2048;  // max. length of the entry
  //
  FILE* fCardsFile;        // pointer on the opened file
  long  fLastPos;            // position in the stream of last record read
  char  fbuffer[fgkMaxLen];// string buffer for current line
  char* fPoint;            // pointer on the beginning of data in the line
  int   fNArgs;            // number of the arguments in the entry
  char** fArgs;            // list of the arguments
  char* fKey;              // current Key
  char* fModifier;         // current Modifier
  TString fLastComment;     // last comments block read
  //
};

//==========================================================================
class MagField: public TVirtualMagField
{
 public:
  enum {kNReg=2};
  MagField(UInt_t id) {SetUniqueID(id);};
  virtual ~MagField() {}
  virtual void Field(const Double_t *xyz, Double_t *bxyz);
  //
  int           GetNReg() const {return (int)kNReg;}
  const double* GetZMin() const {return fZMin;}
  const double* GetZMax() const {return fZMax;}
  const double* GetBVals(int ir) const {return &fBVal[ir][0];}


 protected:
  static const double fZMin[kNReg]; // min z of each field region
  static const double fZMax[kNReg]; // max z of each field region
  static const double fBVal[kNReg][3]; // field values
  //
  ClassDef(MagField, 0) // custom magfield
};


#endif
