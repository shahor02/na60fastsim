#ifndef KMCMSSTATION_H
#define KMCMSSTATION_H

#include <TMath.h>
#include <Rtypes.h>
#include "KMCLayerFwd.h"

struct RotAngle
{
  float phi0;
  float snphi0;
  float csphi0;
  
  RotAngle(float phi=0) : phi0(phi), snphi0(TMath::Sin(phi)), csphi0(TMath::Cos(phi)) {}
  void rotateZ(float x, float y, float& xr, float& yr) const
  {
   xr = x * csphi0 + y * snphi0;
   yr =-x * snphi0 + y * csphi0;
  }
  ClassDefNV(RotAngle, 1);
};

struct Sector
{
  Sector(float phi=0, float dphih=3.1415927, float rmin=0, float rmax=999.) : sectAngle(phi), dphiH(dphih), rMin(rmin), rMax(rmax) {}
  bool rotateToSector(float x, float y, float &xloc, float &yloc) const   // if point inside, return local coordinates
  {
    sectAngle.rotateZ(x,y, xloc, yloc);
    float phi = TMath::ATan2(yloc,xloc);
    return TMath::Abs(phi)<dphiH && x>rMin && x<rMax;
  }  
  RotAngle sectAngle; // bisector
  float dphiH; // phi half coverage
  float rMin; // min rad at bisector
  float rMax; // max rad at bisector
  ClassDefNV(Sector, 1);
};

struct MeasPlane1D : public Sector
{
  MeasPlane1D(float phiSector=0, float dphih=TMath::Pi(), float rmin=0., float rmax=999., float phiMeas=0., float pitch=0.1);
  
  float get1DMeasurement(float x, float y) const
  {
    float xloc, yloc;
    measAngle.rotateZ(x,y, xloc, yloc); // rotate to measurement plane
    return xloc/pitch - offset;
  }
  float pitch;
  float offset;
  RotAngle measAngle; // inclination of the 1d measurement direction wrt mother frame

  ClassDefNV(MeasPlane1D,1);
};

struct KMCMSSector : public Sector
{
  KMCMSSector(float phiSector=0, float dphih=TMath::Pi(), float rmin=0., float rmax=999., float phiUV=0., float pitchUV=0.1, float phiW=TMath::Pi()/2, float pitchW=0.1);

  static KMCMSSector* createSectorDef(float phiSector, float dphih, float rmin, float rmax, float stereoAngle, float pitchUV, float pitchW)
  {
    return new KMCMSSector(phiSector, dphih, rmin, rmax, stereoAngle, pitchUV, dphih+TMath::Pi(), pitchW);
  }
  
  bool getUVW(float x, float y, float& U, float& V, float& W) const;
  MeasPlane1D stripPlaneU;
  MeasPlane1D stripPlaneV;
  MeasPlane1D wirePlaneW;
  
  ClassDefNV(KMCMSSector, 1);
};

struct KMCMSStation : public KMCLayerFwd
{
  KMCMSStation(const char *name = "") : KMCLayerFwd(name), nSectors(0), nRadSegments(0), dPhi(0) {}
  void init(int nsect, const std::vector<float>& r, const std::vector<float>& _phiUV, const std::vector<float>& _pitchUV,
	    const std::vector<float>& _phiW, const std::vector<float>& _pitchW);
  int getSectorID(float x, float y) const;
  KMCMSSector* getSector(int id) { return id<int(sectors.size()) ? &sectors[id] : 0; }
  KMCMSSector* getSector(float x, float y) { return getSector( getSectorID(x,y) ); }

  int nSectors;               // number of sectors (each having 2pi/nSectors coverage)
  int nRadSegments;           // number of radial segments
  float dPhi;                 // sector coverage
  std::vector<float> radii;   // single sector may be made of a few trapezoidal segments delimited by these radii
  std::vector<KMCMSSector> sectors; // sectors ordered in phi.Each sector may have a few trapezoidal segments
  ClassDef(KMCMSStation, 1);
};

#endif
