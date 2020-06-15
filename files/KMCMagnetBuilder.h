#ifndef KMCMAGNETBUILDER_H
#define KMCMAGNETBUILDER_H

#include "Rtypes.h"

class KMCDetectorFwd;

class KMCMagnetBuilder {

 public:
  KMCMagnetBuilder(KMCDetectorFwd* det, float zIn, float dZ=430.5,
                   float matRadL = 8.897, float matDens = 2.699,
                   float defRadL = 0., float defDens = 0.);
  virtual ~KMCMagnetBuilder() {}
  ClassDef(KMCMagnetBuilder,1);
};

#endif
