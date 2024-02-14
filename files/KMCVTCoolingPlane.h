#ifndef KMCVTCOOLINGPLANE_H
#define KMCVTCOOLINGPLANE_H

#include "KMCPolyLayer.h"


struct KMCVTCoolingPlane : KMCPolyLayer
{
  KMCVTCoolingPlane(const char *name = "", float zpos = 0, float thickness = 0,
		    float radL = 0, float density = 0, NaMaterial* mat = 0,     // subtrate (air?)
		    float radLP = 0, float densityP = 0, NaMaterial* matP = 0); // pipes

  ClassDef(KMCVTCoolingPlane,1);
};

#endif
