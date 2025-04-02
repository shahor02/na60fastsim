#include "KMCLayerFwd.h"

ClassImp(KMCLayerFwd)
ClassImp(BeamPipe)

Double_t KMCLayerFwd::fgDefEff = 1.0;

//__________________________________________________________________________
KMCLayerFwd::KMCLayerFwd(const char *name) 
  : TNamed(name,name)
  ,fZ(0)
  ,fThickness(0)
  ,fx2X0(0.)
  ,fXRho(0.)
  ,fEff(0.)
  ,fIsDead(kFALSE)
  ,fNAccReg(1)
  ,fType(-1)
  ,fActiveID(-1)
  ,fSig2EstX(999)
  ,fSig2EstY(999)
  ,fIsRPhiErr(false)
  ,fClCorr()
  ,fClMC()
  ,fClBg("KMCClusterFwd",5)
  ,fTrCorr()
  ,fTrMC("KMCProbeFwd",5)
  ,fMaterial(0)
{
  for (int i=0;i<kMaxAccReg;i++) {
    fRMin[i] = fRMax[i] = -1;
    fXRes[i] = fYRes[i] = 0;
  }
  fRMin[0] = 0;
  fRMax[0] = 999.;

  Reset();
}

//__________________________________________________________________________
void KMCLayerFwd::Reset() 
{
  fTrCorr.Reset();
  fClCorr.Reset();
  ResetMC();
  fSig2EstX = fSig2EstY = 999;
  fMaterial = 0;
  //
}

//__________________________________________________________________________
KMCProbeFwd* KMCLayerFwd::AddMCTrack(KMCProbeFwd* src) 
{
  int ntr = GetNMCTracks(); 
  KMCProbeFwd* prb = 0;
  if (src) prb = new(fTrMC[ntr]) KMCProbeFwd(*src);
  else     prb = new(fTrMC[ntr]) KMCProbeFwd();
  if (!IsDead()) prb->ResetHit(GetActiveID());
  return prb;
}

//__________________________________________________________________________
void KMCLayerFwd::Print(Option_t *opt) const
{
  printf("Lr%3d(A%3d) %15s %+7.3f<Z<%+7.3f X2X0=%.3e XRho=%.3e SigX=%.3e SigY=%.3e %s | Eff:%4.2f RMin:%.3e RMax:%.3e ",
	 GetUniqueID(),fActiveID,GetName(), fZ-fThickness/2,fZ+fThickness/2, fx2X0,fXRho,fXRes[0],fYRes[0], IsRPhiError() ? "(err in RPhi)" : "", fEff,fRMin[0],fRMax[0]);
  for (int ir=1;ir<fNAccReg;ir++) { // print extra regions
    printf("SigX=%.3e SigY=%.3e RMax:%.3e ",fXRes[ir],fYRes[ir], fRMax[ir]);
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

//__________________________________________________________________________
KMCProbeFwd* KMCLayerFwd::GetWinnerMCTrack()  
{
  if (!fTrMC.IsSorted()) fTrMC.Sort();
  KMCProbeFwd* win = fTrMC.GetEntries() ? (KMCProbeFwd*)fTrMC[0]:0;
  if (!win || win->IsKilled()) return 0;
  return win;
}

//__________________________________________________________________________
bool KMCLayerFwd::AddCluster(double x,double y,double z, Int_t id, int clType)
{
  // clTypes are: -1 : ideal MC cluster, 0: signal MC cluster, 1: bg MC clusters
  double r = TMath::Sqrt(x*x + y*y);
  int inAcc = isInAcc(x,y,r);
  double measErr2[3]={}, sgY = GetYRes(r), sgX = GetXRes(r), sgY2 = sgY*sgY, sgX2 = sgX*sgX;
  double phi = TMath::ATan2(y,x), cs = TMath::Cos(phi), sn = TMath::Sin(phi), cs2 = cs*cs, sn2 = sn*sn, cssn = cs*sn;
  if (clType==-1 || inAcc > 0) {
    if (IsRPhiError()) { // rotate error to angle phi
      measErr2[0] = sgX2*cs2+sgY2*sn2;
      measErr2[2] = sgX2*sn2+sgY2*cs2;
      measErr2[1] = (sgY2-sgX2)*cssn;
    } else { // !!! Ylab<->Ytracking, Xlab<->Ztracking
      measErr2[0] = sgY2;
      measErr2[2] = sgX2;
    }
  }
  if (clType==-1) {
    KMCClusterFwd* cl = GetCorCluster();
    cl->Kill(false);
    cl->Set(x, y, z, id);
    cl->SetErr(measErr2[0], measErr2[1], measErr2[2]);
    return true;
  }
  if (inAcc<=0) {
    return false;
  }
  // store randomized cluster local coordinates and phi
  double rx,ry;
  gRandom->Rannor(rx,ry);
  double xerr = rx*sgX, yerr = ry*sgY;
  if (IsRPhiError()) { // rotate track position to R,phi
    r += yerr;
    double rphi = xerr;
    x = r*cs - rphi*sn;
    y = r*sn + rphi*cs;
  } else {
    x += xerr;
    y += yerr;
  }
  
  if (clType==0) { // signal
    GetMCCluster()->Kill(false);
    GetMCCluster()->Set(x, y, z, id);
    GetMCCluster()->SetErr(measErr2[0], measErr2[1], measErr2[2]);    
  } else {
    int ncl = AddBgCluster(x, y, z, id);
    GetBgCluster(ncl-1)->SetErr(measErr2[0], measErr2[1], measErr2[2]);
  }
  return true;
}

//__________________________________________________________________________
void KMCLayerFwd::PrepareForTracking()
{
  //  printf("Lr:%s: Signal clusters: %d, Bg.clusters: %d\n", GetName(), fClMC.IsKilled() ? 0 : 1, fClBg.GetEntriesFast());
}
