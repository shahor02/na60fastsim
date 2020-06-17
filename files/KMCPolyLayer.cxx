#include "KMCPolyLayer.h"

KMCPolygon& KMCPolyLayer::addPolygon(int nv, const float* x, const float* y, float _x2x0, float _xrho)
{
  pieces.emplace_back(nv, x, y, _x2x0, _xrho);
  for (int i=0;i<nv;i++) {
    if (minX>x[i]) minX = x[i];
    if (minY>y[i]) minY = y[i];
    if (maxX<x[i]) maxX = x[i];
    if (maxY<y[i]) maxY = y[i];
    auto xlab = x[i]+sectorOffsX;
    auto ylab = y[i]+sectorOffsY;    
    auto r2 = xlab*xlab + ylab*ylab;
    if (maxR2<r2) maxR2 = r2;
    if (minR2>r2) minR2 = r2;
  }
  return pieces.back();
}

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
  auto r2 = x*x+y*y;
  if (r2<minR2 || r2>maxR2) return -1;
  for (int isec=0;isec<nSectors;isec++) {
    float xloc = sincosSec[isec].second*x + sincosSec[isec].first*y - sectorOffsX;
    float yloc = -sincosSec[isec].first*x + sincosSec[isec].second*y - sectorOffsY;

    if (xloc<minX || xloc>maxX || yloc<minY || yloc>maxY) {
      continue;
    }
    // are we inside the polygons?
    for (int i=pieces.size();i--;) {
      if (pieces[i].isInside(xloc,yloc)) {
	return isec*10+i;
      }
    }
    continue;
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
