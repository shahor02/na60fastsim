#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <TGeoGlobalMagField.h>
#include "TVirtualMagField.h"

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

 protected:
  double fZMin[kNReg]; // min z of each field region
  double fZMax[kNReg]; // max z of each field region
  double fBVal[kNReg][3]; // field values
  //
  ClassDef(MagField, 1) // custom magfield
};


#endif
