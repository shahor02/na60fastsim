#include "KMCClusterFwd.h"
#include <TString.h>

ClassImp(KMCClusterFwd)

//_________________________________________________________________________
KMCClusterFwd::KMCClusterFwd(KMCClusterFwd &src) 
:  TObject(src)
  ,fX(src.fX)
  ,fY(src.fY)
  ,fZ(src.fZ)
  ,fSigYY(src.fSigYY)
  ,fSigYZ(src.fSigYZ)
  ,fSigZZ(src.fSigZZ)
{}

//__________________________________________________________________________
KMCClusterFwd& KMCClusterFwd::operator=(const KMCClusterFwd& src) 
{
  if (this!=&src) {
    TObject::operator=(src);
    fX = src.fX;
    fY = src.fY;
    fZ = src.fZ;
    fSigYY = src.fSigYY;
    fSigYZ = src.fSigYZ;
    fSigZZ = src.fSigZZ;
  }
  return *this;
}

//_________________________________________________________________________
void KMCClusterFwd::Print(Option_t *opt) const 
{
  TString opts = opt;
  opts.ToLower();
  if (opts.Contains("lc")) 
    printf("Tr#%4d TF (%+.4e,%+.4e %+.4e / {%.4e %.4e %.4e}) %s",GetTrID(),GetXTF(), GetYTF(), GetZTF(), GetSigYY(),GetSigYZ(), GetSigZZ(), IsKilled()?"Killed":""); 
  else 
    printf("Tr#%4d Lab (%+.4e,%+.4e %+.4e / {%.4e %.4e %.4e}) %s",GetTrID(),GetX(),GetY(),GetZ(), GetSigYY(),GetSigYZ(), GetSigZZ(), IsKilled()?"Killed":""); 
  if (opts.Contains("nl")) return;
  printf("\n");
}

//_________________________________________________________________________
Bool_t KMCClusterFwd::IsEqual(const TObject* obj) const 
{
  // check if clusters are equal
  KMCClusterFwd* cl = (KMCClusterFwd*)obj;
  const double kTiny = 1e-12;
  if (TMath::Abs(GetX()-cl->GetX())>kTiny || 
      TMath::Abs(GetY()-cl->GetY())>kTiny) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________
Int_t KMCClusterFwd::Compare(const TObject* obj) const
{
  // compare 1st labx, then laby
  KMCClusterFwd* cl = (KMCClusterFwd*)obj;
  if (GetX() > cl->GetX()) return -1; // tracking -Z = labX
  if (GetX() < cl->GetX()) return  1;
  if (GetY() < cl->GetY()) return -1; // tracking Y = labY
  if (GetY() > cl->GetY()) return  1;
  return 0;
}
