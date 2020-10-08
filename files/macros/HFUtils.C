#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLorentzVector.h>
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "../KMCProbeFwd.h"
#endif

// void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, Double_t &xV, Double_t &yV, Double_t &zV);
// void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, KMCProbeFwd &t2, Double_t zStart, Double_t &xV, Double_t &yV, Double_t &zV, Double_t &sigmaVert);
// Double_t CosPointingAngle(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent);
// Double_t ImpParXY(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent);


void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, Double_t &xV, Double_t &yV, Double_t &zV){
  // Based on AliESDv0

  //Trivial estimation of the vertex parameters
  Double_t tmp[3];
  t0.GetXYZ(tmp);
  Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
  Double_t sx1=t0.GetSigmaX2()+ss, sy1=t0.GetSigmaY2()+ss;
  t1.GetXYZ(tmp);
  Double_t  x2=tmp[0],  y2=tmp[1],  z2=tmp[2];
  Double_t sx2=t1.GetSigmaX2()+ss, sy2=t1.GetSigmaY2()+ss; 
    
  Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
  Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
  Double_t wz1=0.5, wz2=1.- wz1;
  xV=wx1*x1 + wx2*x2; 
  yV=wy1*y1 + wy2*y2; 
  zV=wz1*z1 + wz2*z2;
  return;
}

void ComputeVertex(KMCProbeFwd &t0, KMCProbeFwd &t1, KMCProbeFwd &t2, Double_t zStart, Double_t &xV, Double_t &yV, Double_t &zV, Double_t &sigmaVert){

  AliStrLine **linarray = new AliStrLine* [3];
  Double_t pos[3],dir[3],sigmasq[3];
  Double_t wmat[9];
  for(Int_t iel=0; iel<9; iel++) wmat[iel]=0;
  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision

  t0.PropagateToZBxByBz(zStart);
  t1.PropagateToZBxByBz(zStart);
  t2.PropagateToZBxByBz(zStart);

  t0.GetXYZ(pos);
  t0.GetPXYZ(dir);
  sigmasq[0]=t0.GetSigmaX2()+ss;
  sigmasq[1]=t0.GetSigmaY2()+ss;
  sigmasq[2]=1;
  wmat[0]=1./sigmasq[0];
  wmat[5]=1./sigmasq[1];
  wmat[8]=1./sigmasq[2];
  linarray[0]= new AliStrLine(pos,sigmasq,wmat,dir);
  t1.GetXYZ(pos);
  t1.GetPXYZ(dir);
  sigmasq[0]=t1.GetSigmaX2()+ss;
  sigmasq[1]=t1.GetSigmaY2()+ss;
  sigmasq[2]=1;
  wmat[0]=1./sigmasq[0];
  wmat[5]=1./sigmasq[1];
  wmat[8]=1./sigmasq[2];
  linarray[1]= new AliStrLine(pos,sigmasq,wmat,dir);
  t2.GetXYZ(pos);
  t2.GetPXYZ(dir);
  sigmasq[0]=t2.GetSigmaX2()+ss;
  sigmasq[1]=t2.GetSigmaY2()+ss;
  sigmasq[2]=1;
  wmat[0]=1./sigmasq[0];
  wmat[5]=1./sigmasq[1];
  wmat[8]=1./sigmasq[2];
  linarray[2]= new AliStrLine(pos,sigmasq,wmat,dir);
  // do not use errors as weights for the moment
  AliESDVertex vv=AliVertexerTracks::TrackletVertexFinder(linarray,3,kFALSE); 
  xV=vv.GetX();
  yV=vv.GetY();
  zV=vv.GetZ();
  sigmaVert=vv.GetDispersion();
  for(Int_t k=0; k<3; k++) delete linarray[k];
  delete [] linarray;
  return;
}

Double_t ImpParXY(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent){

  Double_t k = -(vsec[0]-vprim[0])*parent.Px()-(vsec[1]-vprim[1])*parent.Py();
  k /= parent.Pt()*parent.Pt();
  Double_t dx = vsec[0]-vprim[0]+k*parent.Px();
  Double_t dy = vsec[1]-vprim[1]+k*parent.Py();
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy);
  TVector3 mom(parent.Px(), parent.Py(), parent.Pz());
  TVector3 fline(vsec[0] - vprim[0],
		 vsec[1] - vprim[1],
		 vsec[2] - vprim[2]);
  TVector3 cross = mom.Cross(fline);
  return (cross.Z()>0. ? absImpPar : -absImpPar);
}

Double_t CosPointingAngle(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent)
{

    // /// XY
    // /// Cosine of pointing angle in transverse plane assuming it is produced
    // /// at "point"

    // TVector3 momXY(parent.Px(), parent.Py(), 0.);
    // TVector3 flineXY(vsec[0] - vprim[0],
    //                  vsec[1] - vprim[1],
    //                  0.);

    // Double_t ptot2 = momXY.Mag2() * flineXY.Mag2();
    // if (ptot2 <= 0)
    // {
    //     return 0.0;
    // }
    // else
    // {
    //     Double_t cos = momXY.Dot(flineXY) / TMath::Sqrt(ptot2);
    //     if (cos > 1.0)
    //         cos = 1.0;
    //     if (cos < -1.0)
    //         cos = -1.0;
    //     return cos;
    // }
    /// Cosine of pointing angle in space assuming it is produced at "point"

    TVector3 mom(parent.Px(), parent.Py(), parent.Pz());
    TVector3 fline(vsec[0] - vprim[0],
                   vsec[1] - vprim[1],
                   vsec[2] - vprim[2]);

    Double_t ptot2 = mom.Mag2() * fline.Mag2();
    if (ptot2 <= 0)
    {
        return 0.0;
    }
    else
    {
        Double_t cos = mom.Dot(fline) / TMath::Sqrt(ptot2);
        if (cos > 1.0)
            cos = 1.0;
        if (cos < -1.0)
            cos = -1.0;
        return cos;
    }
}
