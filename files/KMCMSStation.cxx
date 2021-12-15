#include "KMCMSStation.h"

MeasPlane1D::MeasPlane1D(float phiSector, float dphih, float rmin, float rmax, float phiMeas, float ptc)
  : Sector(phiSector, dphih, rmin, rmax), measAngle(phiMeas), pitch(ptc)
{
  // find offset
  float x,y;
  float chan = get1DMeasurement(rMin, rMin*TMath::Tan(dphiH));
  chan = TMath::Min(chan, get1DMeasurement(rMin, -rMin*TMath::Tan(dphiH)));
  chan = TMath::Min(chan, get1DMeasurement(rMax,  rMax*TMath::Tan(dphiH)));
  chan = TMath::Min(chan, get1DMeasurement(rMax, -rMax*TMath::Tan(dphiH)));
  offset = chan;
}

KMCMSSector::KMCMSSector(float phiSector, float dphih, float rmin, float rmax,
		     float phiUV, float pitchUV,
		     float phiW, float pitchW)
  : Sector(phiSector, dphih, rmin, rmax),
    stripPlaneU(phiSector, dphih, rmin, rmax, phiUV, pitchUV),
    stripPlaneV(phiSector, dphih, rmin, rmax, -phiUV, pitchUV),
    wirePlaneW(phiSector, dphih, rmin, rmax, phiW, pitchW)
{}
    
bool KMCMSSector::getUVW(float x, float y, float& U, float& V, float& W) const
{
  float xl,yl;
  if (!rotateToSector(x,y, xl, yl)) return false;
  U = stripPlaneU.get1DMeasurement(xl,yl);
  V = stripPlaneV.get1DMeasurement(xl,yl);
  W = wirePlaneW.get1DMeasurement(xl,yl);
  return true;
}
