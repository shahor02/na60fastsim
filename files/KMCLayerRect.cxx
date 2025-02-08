#include "KMCLayerRect.h"

KMCLayerRect::KMCLayerRect(const char *name, Float_t radL, Float_t density, Float_t zPos, Float_t thickness, Float_t xRes, Float_t yRes, Float_t xHole, Float_t yHole, float xHSide, float yHSide, Float_t eff) : KMCLayerFwd(name)
{
  SetZ(zPos);
  SetThickness(thickness);
  SetX2X0( radL>0 ? thickness*density/radL : 0);
  SetXTimesRho(thickness*density);
  SetLayerEff(eff);
  SetDead(eff > 0 ? false : true );
  fHole[0] = xHole;
  fHole[1] = yHole;
  setRegionData(0, xHSide, yHSide, xRes, yRes);
}

KMCLayerRect::KMCLayerRect()
{
  fHole[0] = fHole[1] = 0.;
  for (int i=0;i<kMaxAccReg;i++) {
    fHalfX[i] = fHalfY[i] = 0.;
  }
}

int KMCLayerRect::GetRegionID(float x, float y, float r) const // determine if x,y (or r if > 0) is in the acceptance
{
  if (isInAcc(x,y,r) == 0) return -1;
  for (int i=0;i<fNAccReg;i++) {
    if (x>-fHalfX[i] && x<fHalfX[i] && y>-fHalfY[i] && y<fHalfY[i]) return i;
  }
  return -1;
}

int KMCLayerRect::GetAccRegion(const KMCProbeFwd* tr) const
{
  return GetRegionID(tr->GetX(), tr->GetY());
}

int KMCLayerRect::isInAcc(float x, float y, float r) const
{	     
  if (fHole[0] > 0) { // there is a central hole
    if (fHole[1] <= 0.) { // it is round round 
      if (r > 0) {
	if (r < fHole[0]) return 0;
      } else {
	float r2 = (r>0.) ? r*r : x*x + y*y;
	if (r2 < fHole[0]*fHole[0]) return 0;
      } 
    } else {
      if (TMath::Abs(x)<fHole[0] && TMath::Abs(y)<fHole[1]) return 0;
    }
  }
  if (x<-fHalfX[fNAccReg-1] || x>fHalfX[fNAccReg-1] || y<-fHalfY[fNAccReg-1] || y>fHalfY[fNAccReg-1]) {
    return 0;
  }
  return 1;
}

//__________________________________________________________________________
void KMCLayerRect::Print(Option_t *opt) const
{
  printf("Lr%3d(A%3d) %15s %+7.3f<Z<%+7.3f X2X0=%.3e XRho=%.3e SigX=%.3e SigY=%.3e | Eff:%4.2f Hole:",
	 GetUniqueID(),fActiveID,GetName(), fZ-fThickness/2,fZ+fThickness/2, fx2X0,fXRho,fXRes[0],fYRes[0], fEff);
  if (fHole[0]>0) {
    if (fHole[1]>0) printf(" |X|<%4.2f |Y|<%4.2f ", fHole[0],fHole[1]);
    else printf(" R<%4.2f ", fHole[0]);
  } else {
    printf(" None. ");
  }
  printf(" AccReg: ");
  for (int ir=0;ir<fNAccReg;ir++) { // print extra regions
    printf("SigX=%.3e SigY=%.3e |X|<%.3e |Y|<%.3e ",fXRes[ir],fYRes[ir], fHalfX[ir], fHalfY[ir]);
  }
  printf("\n");
  TString opts = opt; opts.ToLower();
  //
  if (opts.Contains("cl")) {
    printf("Clusters: MC: "); fClMC.Print(opts+"nl");
    printf("  Corr: "); fClCorr.Print(opts+"nl");
    printf("  NBgCl: %3d NTrMC: %4d\n",GetNBgClusters(),GetNMCTracks());
  }
  if (opts.Contains("bcl")) fClBg.Print(opt);
}
