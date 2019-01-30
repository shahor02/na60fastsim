#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <TGeoGlobalMagField.h>

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
