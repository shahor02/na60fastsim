#include "KMCLayerFwd.h"

ClassImp(KMCLayerFwd)
ClassImp(BeamPipe)

Double_t KMCLayerFwd::fgDefEff = 1.0;

//__________________________________________________________________________
KMCLayerFwd::KMCLayerFwd(const char *name) 
  : TNamed(name,name)
  ,fZ(0)
  ,fRMin(0)
  ,fRMax(1e9)
  ,fThickness(0)
  ,fx2X0(0)
  ,fXRho(0)
  ,fXRes(0)
  ,fYRes(0)
  ,fEff(0)
  ,fIsDead(kFALSE)
  ,fType(-1)
  ,fActiveID(-1)
  ,fSig2EstX(999)
  ,fSig2EstY(999)
  ,fClCorr()
  ,fClMC()
  ,fClBg("KMCClusterFwd",5)
  ,fTrCorr()
  ,fTrMC("KMCProbeFwd",5)
{
  Reset();
}

//__________________________________________________________________________
void KMCLayerFwd::Reset() 
{
  fTrCorr.Reset();
  fClCorr.Reset();
  ResetMC();
  fSig2EstX = fSig2EstY = 999;
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
  printf("Lr%3d(A%3d) %15s Z=%+7.1f DZ=%7.3f X2X0=%.3e XRho=%.3e SigX=%.3e SigY=%.3e RMin:%.3e RMax:%.3e Eff:%4.2f\n",
	 GetUniqueID(),fActiveID,GetName(), fZ, fThickness, fx2X0,fXRho,fXRes,fYRes,fRMin,fRMax,fEff);
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

