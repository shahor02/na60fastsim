#include "KMCPolyLayer.h"

void KMCPolyLayer::setNSectorsPhiStart(int n, float phi)
{
  phiStart = phi;
  if (phi>TMath::Pi()) phi -= TMath::Pi()*2; // bring to -pi : pi
  nSectors = n>1 ? n : 1;
  sectorCoverage = TMath::Pi()*2./n;
  sectorCoverageInv = 1./sectorCoverage;
  for (int i=0;i<n;i++) {
    float phi = phiStart + sectorCoverage*i;
    sincosSec.push_back( std::pair<float,float>(TMath::Sin(phi), TMath::Cos(phi)) );
  }
}


int KMCPolyLayer::getPolygonID(float x,float y) const {
  //
  // determine in which sector we are
  //  int pieceID = -1;
  for (int isec=0;i<nSectors;i++) {
    float xloc = sincosSec[sector].second*x + sincosSec[sector].first*y - sectorOffsX;
    float yloc = -sincosSec[sector].first*x + sincosSec[sector].second*y - sectorOffsY;
    // did this bring us to 1st sector (0:sectorCoverage)?
    float phi = TMath::ATan2(y,x); // in -pi:pi
    if (phi<0) phi += TMath::Pi()*2;
    if (phi>sectorCoverage) continue;

    // are we inside the polygons?
    for (int i=pieces.size();i--;) {
      if (pieces[i].isInside(xloc,yloc)) {
	return i;
	//	pieceID = i;
	//	break; // inside!!! return?
      }
    }
    break;
  }
  /*
  float phi = TMath::ATan2(y,x); // in -pi:pi
  float phis = phi - phiStart;
  if (phis<0) phis += TMath::Pi()*2;
  int sector = phis*sectorCoverageInv;
  // rotate point to sector coordinates
  float xloc = sincosSec[sector].second*x + sincosSec[sector].first*y;
  float yloc = -sincosSec[sector].first*x + sincosSec[sector].second*y;
  for (int i=pieces.size();i--;) {
    if (pieces[i].isInside(xloc,yloc)) {
      return i;
    }
  }
  */
  return -1;
}
