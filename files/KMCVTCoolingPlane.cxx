#include "KMCVTCoolingPlane.h"

KMCVTCoolingPlane::KMCVTCoolingPlane(const char *name, float zpos, float thickness,
			     float radL, float density, NaMaterial* mat,     // subtrate (air?)
			     float radLP, float densityP, NaMaterial* matP) : KMCPolyLayer(name)
{
  // subtrate settings
  SetZ(zpos);
  SetThickness(thickness);
  SetX2X0( radL>0 ? thickness*density/radL : 0);
  SetXTimesRho(thickness*density);
  SetMaterial(mat);
  //
  // sensors resolution and efficiency
  const float dpipe = 0.3f; // pipe transverse thickness in cm
  const int npoints = 16;
  // 
  const float contourX[npoints] = {
    -16.0f, // >> outer edge of the contour
    +15.5f,
    +15.5f,
     -0.5f,
     -0.5f,
    -15.5f,
    -15.5f,
    -16.0f, // << outer edge of the contour
    //
    -16.0f, // << returning part of the contour
    -15.5f-dpipe,
    -15.5f-dpipe,
    -0.5f+dpipe,
    -0.5f+dpipe,
    +15.5f-dpipe,
    +15.5f-dpipe,
    -16.0f // >> returning part of the contour
  };
  
  const float contourY[npoints] = {
    +16.0f, // >> outer edge of the contour
    +16.0f,
    -15.5f,
    -15.5f,
    +15.1f,
    +15.1f,
    -15.5f,
    -15.5f, //<< outer edge of the contour
    //
    -15.5f+dpipe, // << returning part of the contour
    -15.5f+dpipe,
    +15.1f+dpipe,
    +15.1f+dpipe,
    -15.5f+dpipe,
    -15.5f+dpipe,
    +16.0f-dpipe,
    +16.0f-dpipe
  };
  addPolygon(npoints, contourX, contourY, (radLP>0 && densityP>0) ? thickness/(radLP/densityP) : 0, thickness*densityP);
  setNSectorsPhiStart(1, 0);
  SetXRes(999999);
  SetYRes(999999);
  SetDead(kTRUE);
}

