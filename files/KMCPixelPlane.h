#ifndef KMCPIXELPLANE_H
#define KMCPIXELPLANE_H

#include "KMCPolyLayer.h"

#define _NOVTPLANEGAPS_


struct KMCPixelPlane : public KMCPolyLayer
{
  KMCPixelPlane(const char *name, float zpos, float thickness,
		float radL, float density, NaMaterial* mat, // substrate
		float radLS, float densityS, NaMaterial* matS, // sensor
		float offX, float offsY, float interChipGap, // offset of the inner corner of the 1st quadrant chip
		float xRes, Float_t yRes, Float_t eff);
  
  virtual int isInAcc(float x, float y, float r=-1) const;
  
  NaMaterial* sensMat = nullptr;
  float offsX = 0;
  float offsY = 0;
  float interChipGap = 0;
  // left/right top/bottom are with respect to the 0/0 pixel in the left bottom corner
  static const float XSizeTot;
  static const float YSizeTot; // readout side
  static const float DeadXLong;  // right end cap
  static const float DeadXShort; // left end cap (readout)
  static const float DeadYBottom; // bottom dead zone
  static const float DeadYTop; // top dead zone
  static const float DeadTopBotHalves; // dead space between top and bottom halves
  static const float DeadXTile;  // dead space between tiles in X
  static const float DeadXDataBackbone; // dead space between every segment (triplet of tiles)
  static const int NXTiles;
  static const int NXSegments;
  static const int NYSensors;
  static const float DXTile;
  static const float DXSegment;
  static const float DYSens;
  
  ClassDef(KMCPixelPlane,1);
};

#endif
