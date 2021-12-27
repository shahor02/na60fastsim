#include "KMCMSStation.h"

MeasPlane1D::MeasPlane1D(float phiSector, float dphih, float rmin, float rmax, float phiMeas, float ptc)
  : Sector(phiSector, dphih, rmin, rmax), pitch(ptc), offset(0), measAngle(phiMeas)
{
  // find offset
  float chan = get1DMeasurement(rMin, rMin*TMath::Tan(dphiH));
  chan = TMath::Min(chan, get1DMeasurement(rMin, -rMin*TMath::Tan(dphiH)));
  chan = TMath::Min(chan, get1DMeasurement(rMax,  rMax*TMath::Tan(dphiH)));
  chan = TMath::Min(chan, get1DMeasurement(rMax, -rMax*TMath::Tan(dphiH)));
  offset = chan;
}

KMCMSSector::KMCMSSector(float phiSector, float dphih, float rmin, float rmax,
			 float phiUV, float pitchUV,
			 float phiW, float pitchW, float _sigR, float _sigRPhi)
  : Sector(phiSector, dphih, rmin, rmax),
    stripPlaneU(phiSector, dphih, rmin, rmax, phiUV, pitchUV),
    stripPlaneV(phiSector, dphih, rmin, rmax, -phiUV, pitchUV),
    wirePlaneW(phiSector, dphih, rmin, rmax, phiW, pitchW),
    sigR(_sigR), sigRPhi(_sigRPhi)
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


void KMCMSStation::init(int nsect, const std::vector<float>& r, const std::vector<float>& _phiUV, const std::vector<float>& _pitchUV,
	    const std::vector<float>& _phiW, const std::vector<float>& _pitchW,
	    const std::vector<float>& _sigR, const std::vector<float>& _sigRPhi)
{
  dPhi = TMath::Pi()*2./nsect;
  nRadSegments = r.size()-1;
  for (int ip=0;ip<nsect;ip++) {
    float phiSect = (ip+0.5)*dPhi;
    for (int ir=0;ir<nRadSegments;ir++) {
      sectors.emplace_back(phiSect, 0.5*dPhi, r[ir], r[ir+1], _phiUV[ir], _pitchUV[ir], _phiW[ir], _pitchW[ir], _sigR[ir], _sigRPhi[ir]);
    }
  }
}

int KMCMSStation::getSectorID(float x, float y) const
{
  float phi = TMath::ATan2(y,x);
  if (phi<0) phi += TMath::Pi()*2;
  int id = phi/dPhi;
  if (id>=nSectors) id = nSectors-1;
  id *= nRadSegments;
  const auto& s0 = sectors[id]; // lowest segment of given sector;
  float xl, yl;
  s0.rotateToSector(x,y, xl,yl); // in this frame local X axis goes along bisector of the sector
  if (xl<s0.rMin) return -1; // not in the acceptance
  for (int ir=0;ir<nRadSegments;ir++) {
    if (xl<sectors[id].rMax) {
      return id;
    }
    id++;
  }
  return -1; // missed from top
}
