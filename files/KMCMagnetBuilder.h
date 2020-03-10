#ifndef KMCMAGNETBUILDER_H
#define KMCMAGNETBUILDER_H

#include "Rtypes.h"

class KMCDetectorFwd;

class KMCMagnetBuilder {

  KMCMagnetBuilder(KMCDetectorFwd* det, float z);

  ClassDef(KMCMagnetBuilder,1);
};

#endif
