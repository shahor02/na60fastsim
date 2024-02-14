#include "KMCPixelPlane.h"

KMCPixelPlane::KMCPixelPlane(const char *name, float zpos, float thickness,
			     float radL, float density, NaMaterial* mat, // substrate
			     float radLS, float densityS, NaMaterial* matS, // sensor
			     float sizeSX, float sizeSY, // single sensor size in X, and Y
			     float offsX, float offsY, // offset of the inner corner of the 1st quadrant chip
			     float xRes, Float_t yRes, Float_t eff
			     ) : KMCPolyLayer(name)
{
  // subtrate settings
  SetZ(zpos);
  SetThickness(thickness);
  SetX2X0( radL>0 ? thickness*density/radL : 0);
  SetXTimesRho(thickness*density);
  SetMaterial(mat);
  //
  // sensors resolution and efficiency
  SetXRes(xRes);
  SetYRes(yRes);
  SetLayerEff(eff);
  SetDead(false);
  //
  sensMat = matS;
  // define chip contour
  setXYOffsets(offsX, offsY);
  float contourX[4] = {}, contourY[4] = {};
  contourX[0] = 0; // inner corrner
  contourY[0] = 0;
  contourX[1] = sizeSX; // bottom right corrner
  contourY[1] = 0;
  contourX[2] = sizeSX; // top right corrner
  contourY[2] = sizeSY;
  contourX[3] = sizeSX; // top lefy corrner
  contourY[3] = sizeSY;  
  addPolygon(4, contourX, contourY, (radLS>0 && densityS>0) ? thickness/(radLS/densityS) : 0, thickness*densityS);
  setNSectorsPhiStart(4, 0.0); // define full set of chips as rotated copies of the 1st one
}
