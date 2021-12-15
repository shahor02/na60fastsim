#ifndef KMCMSSTATION_H
#define KMCMSSTATION_H

#include <TMath.h>
#include <Rtypes.h>

struct RotAngle
{
  float phi0 = 0;
  float snphi0 = 0;
  float csphi0 = 0;
  
  RotAngle(float phi) : phi0(phi), snphi0(TMath::Sin(phi)), csphi0(TMath::Cos(phi)) {}
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
  float dphiH = 0; // phi half coverage
  float rMin = 0; // min rad at bisector
  float rMax = 0; // max rad at bisector
  ClassDefNV(Sector, 1);
};

struct MeasPlane1D : public Sector
{
  MeasPlane1D(float phiSector, float dphih, float rmin, float rmax, float phiMeas, float pitch);
  
  float get1DMeasurement(float x, float y) const
  {
    float xloc, yloc;
    measAngle.rotateZ(x,y, xloc, yloc); // rotate to measurement plane
    return xloc/pitch - offset;
  }
  float pitch = 1;
  float offset = 0;
  RotAngle measAngle; // inclination of the 1d measurement direction wrt mother frame

  ClassDefNV(MeasPlane1D,1);
};

struct KMCMSSector : public Sector
{
  KMCMSSector(float phiSector, float dphih, float rmin, float rmax, float phiUV, float pitchUV, float phiW, float pitchW);

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

#endif
