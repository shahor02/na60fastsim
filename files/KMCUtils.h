#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <TGeoGlobalMagField.h>
#include "TVirtualMagField.h"
#include "TF1.h"

//==========================================================================
class MagField: public TVirtualMagField
{
 public:
  enum {kNReg=2};
  MagField(UInt_t id);
  virtual ~MagField() {}
  virtual void Field(const Double_t *xyz, Double_t *bxyz);
  //
  virtual int GetNReg() const {return (int)kNReg;}
  virtual const double* GetZMin() const {return fZMin;}
  virtual const double* GetZMax() const {return fZMax;}
  const double* GetBVals(int ir) const {return &fBVal[ir][0];} 
  virtual void SetZMin(int nreg, double zmin) { fZMin[nreg] = zmin; }
  virtual void SetZMax(int nreg, double zmin) { fZMax[nreg] = zmin; }
  virtual void SetBVals(int nreg, int index, double val) { fBVal[nreg][index] = val; }

  void SetMagVTParam(const std::string& paramString, const std::vector<float>& parVals);
  void SetMagMSParam(const std::string& paramString, const std::vector<float>& parVals);
  
 protected:
  double fZMin[kNReg] = {0., 380.}; // min z of each field region
  double fZMax[kNReg] = {40.,680.}; // max z of each field region
  double fBVal[kNReg][3] = {{-30,0,0},{-20,0,0}}; // field values
  //
  std::array<TF1*,kNReg> fFieldFunReg{};
  
  ClassDef(MagField, 1) // custom magfield
};


#endif
