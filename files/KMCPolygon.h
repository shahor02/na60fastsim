#ifndef KMCPOLYGON_H
#define KMCPOLYGON_H

struct KMCPolyVertex {
  float x=0,y=0;
  KMCPolyVertex(float _x=0,float _y=0) : x(_x), y(_y) {}; 
};

struct KMCPolygon
{
  std::vector<KMCPolyVertex> poly;
  KMCPolygon() {}

  void addVertex(float x,float y) {
    poly.emplace_back(x,y);
  }
  
  bool IsInside(float x, float y) {
    int n = poly.size(), j = n-1;
    bool odd = false;
    for (int i=0;i<n;i++) {
      if ( (poly[i].y<y && poly[j].y>=y || poly[j].y<y && poly[i].y>=y) && (poly[i].x<=x && poly[j].x<=x)) {
	odd ^= (poly[i].x+(y-poly[i].y)/(poly[j].y-poly[i].y)*(poly[j].x-poly[i].x)<x);
      }
      j=i;
    }
    return odd;
  }
};

#endif
