#ifndef KMCPIXELPLANE_H
#define KMCPIXELPLANE_H

#include "KMCPolyLayer.h"


struct KMCPixelPlane : public KMCPolyLayer
{
  KMCPixelPlane(const char *name, float zpos, float thickness,
		float radL, float density, NaMaterial* mat, // substrate
		float radLS, float densityS, NaMaterial* matS, // sensor
		float sizeSX, float sizeSY, // single sensor size in X, and Y
		float offX, float offsY, // offset of the inner corner of the 1st quadrant chip
		float xRes, Float_t yRes, Float_t eff);
  
  virtual bool isInAcc(float x, float y, float r=-1) const { return getPolygonID(x,y) >= 0; }
  
  NaMaterial* sensMat;
  
  ClassDef(KMCPixelPlane,1);
};

#endif
