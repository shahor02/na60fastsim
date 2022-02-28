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
  void rotateZToLocal(float x, float y, float& xr, float& yr) const
  {
   xr = x * csphi0 + y * snphi0;
   yr =-x * snphi0 + y * csphi0;
  }
  void rotateZFromLocal(float x, float y, float& xr, float& yr) const
  {
   xr = x * csphi0 - y * snphi0;
   yr = x * snphi0 + y * csphi0;
  }
  ClassDefNV(RotAngle, 1);
};

struct Sector
{
  Sector(float phi=0, float dphih=3.1415927, float rmin=0, float rmax=999.) : sectAngle(phi), dphiH(dphih), rMin(rmin), rMax(rmax) {}
  bool isInside(float xS, float yS) const
  {
    // check if the point with sector coordinate belongs to the sector
    float phi = TMath::ATan2(yS,xS);
    return TMath::Abs(phi)<dphiH && xS>rMin && xS<rMax;
  }
  bool rotateToSector(float xlab, float ylab, float &xloc, float &yloc) const   // if point inside, return local coordinates
  {
    sectAngle.rotateZToLocal(xlab,ylab, xloc, yloc);
    return isInside(xloc, yloc);
  }
  void sector2Lab(float xs, float ys, float &xl, float &yl) const
  {
    sectAngle.rotateZFromLocal(xs,ys, xl, yl);
  }
  
  RotAngle sectAngle; // bisector
  float dphiH; // phi half coverage
  float rMin; // min rad at bisector
  float rMax; // max rad at bisector
  ClassDefNV(Sector, 1);
};

struct FiredChannel
{
  float channel;
  short label;
  FiredChannel(float _c=0, short _l=0) : channel(_c), label(_l) {}
  ClassDefNV(FiredChannel,1);
};


struct MeasPlane1D : public Sector
{
  MeasPlane1D(float phiSector=0, float dphih=TMath::Pi(), float rmin=0., float rmax=999., float phiMeas=0., float pitch=0.1);
  
  float get1DMeasurement(float x, float y) const
  {
    // get channel corresponding to sector coordinate x,y (x axis along bisector)
    float xloc, yloc;
    measAngle.rotateZToLocal(x,y, xloc, yloc); // rotate to measurement plane
    return yloc/pitch - offset;
  }
  void local2sector(float channel, float t, float &x, float& y)
  {
    // get sector coordinates for channel and a pont on the line along the channel
    float q = (channel + offset)*pitch;
    measAngle.rotateZFromLocal(t, q, x, y); // rotate from measurement to sector frame
  }
  
  void clear() { hits.clear(); }
  
  std::vector<FiredChannel> hits;
  float pitch;
  float offset;
  RotAngle measAngle; // inclination of the 1d measurement direction wrt mother frame

  ClassDefNV(MeasPlane1D,1);
};

struct KMCMSSector : public Sector
{
  KMCMSSector(float phiSector=0, float dphih=TMath::Pi(), float rmin=0., float rmax=999., float phiUV=0., float pitchUV=0.1,
	      float phiW=TMath::Pi()/2, float pitchW=0.1, float _sigR=0.03, float _sigRPhi=0.03);

  static KMCMSSector* createSectorDef(float phiSector, float dphih, float rmin, float rmax, float stereoAngle, float pitchUV, float pitchW)
  {
    return new KMCMSSector(phiSector, dphih, rmin, rmax, stereoAngle, pitchUV, dphih+TMath::Pi(), pitchW);
  }
  
  bool getUVW(float x, float y, float& U, float& V, float& W) const;
  void clear();
  MeasPlane1D stripPlaneU;
  MeasPlane1D stripPlaneV;
  MeasPlane1D wirePlaneW;
  double cUV, sUV; // cos and sine of U,V angle differences
  double cUW, sUW; // cos and sine of U,W angle differences
  double cVW, sVW; // cos and sine of V,W angle differences
  float sigR;
  float sigRPhi;
  ClassDefNV(KMCMSSector, 1);
};

struct KMCMSStation : public KMCLayerFwd
{
  KMCMSStation(const char *name = "") : KMCLayerFwd(name), nSectors(0), nRadSegments(0), dPhi(0),
    signalU(-1.), signalV(-1.), signalW(-1.), signalSectorID(-1) {}
  
  void init(int nsect, const std::vector<float>& r, const std::vector<float>& _phiUV, const std::vector<float>& _pitchUV,
	    const std::vector<float>& _phiW, const std::vector<float>& _pitchW,
	    const std::vector<float>& _sigR, const std::vector<float>& _sigRPhi);
  int getSectorID(float x, float y) const;
  KMCMSSector* getSector(int id) { return id<int(sectors.size()) ? &sectors[id] : 0; }
  const KMCMSSector* getSector(int id) const { return id<int(sectors.size()) ? &sectors[id] : 0; }
  KMCMSSector* getSector(float x, float y) { return getSector( getSectorID(x,y) ); }
  virtual bool AddCluster(double x,double y,double z, Int_t id, int clType);
  virtual void Print(Option_t *opt) const;
  virtual void SetRPhiError(bool );
  virtual void PrepareForTracking();
  void ClearPrimaryHits();
  
  int nSectors;               // number of sectors (each having 2pi/nSectors coverage)
  int nRadSegments;           // number of radial segments
  float dPhi;                 // sector coverage
  float signalU, signalV, signalW; // signal channels cache
  int signalSectorID;
  std::vector<float> radii;   // single sector may be made of a few trapezoidal segments delimited by these radii
  std::vector<KMCMSSector> sectors; // sectors ordered in phi.Each sector may have a few trapezoidal segments
  ClassDef(KMCMSStation, 1);
};

#endif
