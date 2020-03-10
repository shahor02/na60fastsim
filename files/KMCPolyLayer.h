#ifndef KMCPOLYLAYER_H
#define KMCPOLYLAYER_H

#include "KMCLayerFwd.h"


struct KMCPolyVertex {
  float x=0,y=0;
  KMCPolyVertex(float _x=0,float _y=0) : x(_x), y(_y) {}; 
};

// poligonal projection of an obstacle on XY plane
struct KMCPolygon
{
  std::vector<KMCPolyVertex> poly;
  float x2x0 = 0.;
  float xrho = 0.;

  KMCPolygon() {}

  void addVertex(float x,float y) {
    poly.emplace_back(x,y);
  }
  
  bool isInside(float x, float y) const {
    int n = poly.size(), j = n-1;
    bool odd = false;
    for (int i=0;i<n;i++) {
      if ( ( (poly[i].y<y && poly[j].y>=y) || (poly[j].y<y && poly[i].y>=y)) && (poly[i].x<=x && poly[j].x<=x)) {
	odd ^= (poly[i].x+(y-poly[i].y)/(poly[j].y-poly[i].y)*(poly[j].x-poly[i].x)<x);
      }
      j=i;
    }
    return odd;
  }
  ClassDef(KMCPolygon,1);
};

// set of polygonal projections 
struct KMCPolyLayer : public KMCLayerFwd
{
  std::vector<KMCPolygon> pieces;
  
  KMCPolyLayer(const char *name) : KMCLayerFwd(name) {}

  // get ID of the polygon containing the point (if any)
  inline int getPolygonID(float x,float y) const
  {
    for (int i=pieces.size();i--;) {
      if (pieces[i].isInside(x,y)) {
	return i;
      }
    }
    return -1;
  }

  inline void getMatBudget(int pid, float &_x2x0, float &_xrho) const
  {
    _x2x0 = GetX2X0();
    _xrho = GetXTimesRho(); 
    if (pid<0) return;
    _x2x0 = pieces[pid].x2x0;
    _xrho = pieces[pid].xrho;
  }
  
  inline virtual void getMatBudget(float x, float y, float &_x2x0, float &_xrho) const
  {
    getMatBudget(getPolygonID(x,y), _x2x0, _xrho);
  }

  ClassDef(KMCPolyLayer,1);
};

#endif
