#ifndef KMCLAYERRECT_H
#define KMCLAYERRECT_H

#include "KMCLayerFwd.h"
#include <vector>
#include <utility>

// rectangular layer

class KMCLayerRect : public KMCLayerFwd
{
 public:
  KMCLayerRect(const char *name, Float_t radL, Float_t density, Float_t zPos, Float_t thickness, Float_t xRes, Float_t yRes, Float_t xHole, Float_t yHole, float xHSide, float yHSide, Float_t eff);
  KMCLayerRect();
  virtual ~KMCLayerRect() {}
  virtual void  Print(Option_t* option = "") const;
  
  int GetRegionID(float x, float y, float r=-1) const; // determine if x,y (or r if > 0) is in the acceptance  
  virtual int GetAccRegion(const KMCProbeFwd* tr) const;
  virtual bool isInAcc(float x, float y, float r=-1) const;
  virtual Float_t GetXRes(const KMCProbeFwd* tr)  const {
    int id = GetRegionID(tr->GetX(), tr->GetY());
    return id<0 ? fXRes[0] : fXRes[id];
  }
  virtual Float_t GetYRes(const KMCProbeFwd* tr)  const {
    int id = GetRegionID(tr->GetX(), tr->GetY());
    return id<0 ? fYRes[0] : fYRes[id];
  }

  float getXSideHalf(int i) { return i<fNAccReg ? fHalfX[i] : -1;}
  float getYSideHalf(int i) { return i<fNAccReg ? fHalfY[i] : -1;}  

  void setRegionData(int i, float xh, float yh, float resX, float resY)
  {
    fHalfX[i] = xh;
    fHalfY[i] = yh;
    fXRes[i] = resX;
    fYRes[i] = resY;
  }
  
  ///////////////////////////
 private :
  float fHole[2]; // central hole. hole[0]>0 && hole[1]==0 means round hole with radius hole[0], if both > 0 : rectrangular hole half-sizes
  float fHalfX[kMaxAccReg];
  float fHalfY[kMaxAccReg];
  
  ClassDef(KMCLayerRect,1);
};

#endif
