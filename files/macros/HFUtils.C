#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TLorentzVector.h>
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "KMCProbeFwd.h"
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

//--------------------------------------------------------------------------  
void GetStrLinDerivMatrix(const Double_t *p0,const Double_t *p1,Double_t (*m)[3],Double_t *d)
{
  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t kk=x12*x12+y12*y12+z12*z12;
  m[0][0]=2-2/kk*x12*x12;
  m[0][1]=-2/kk*x12*y12;
  m[0][2]=-2/kk*x12*z12;
  m[1][0]=-2/kk*x12*y12;
  m[1][1]=2-2/kk*y12*y12;
  m[1][2]=-2/kk*y12*z12;
  m[2][0]=-2/kk*x12*z12;
  m[2][1]=-2/kk*y12*z12;
  m[2][2]=2-2/kk*z12*z12;
  d[0]=2*p0[0]-2/kk*p0[0]*x12*x12-2/kk*p0[2]*x12*z12-2/kk*p0[1]*x12*y12;
  d[1]=2*p0[1]-2/kk*p0[1]*y12*y12-2/kk*p0[0]*x12*y12-2/kk*p0[2]*z12*y12;
  d[2]=2*p0[2]-2/kk*p0[2]*z12*z12-2/kk*p0[0]*x12*z12-2/kk*p0[1]*z12*y12;
}

//------------------------------------------------------------------------
Double_t GetDeterminant3X3(Double_t matr[][3])
{
  //
  Double_t det=matr[0][0]*matr[1][1]*matr[2][2]-matr[0][0]*matr[1][2]*matr[2][1]-matr[0][1]*matr[1][0]*matr[2][2]+matr[0][1]*matr[1][2]*matr[2][0]+matr[0][2]*matr[1][0]*matr[2][1]-matr[0][2]*matr[1][1]*matr[2][0];
 return det;
}
//--------------------------------------------------------------------------   
Double_t GetStrLinMinDist(const Double_t *p0,const Double_t *p1,const Double_t *x0)
{
  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t x10=p0[0]-x0[0];
  Double_t y10=p0[1]-x0[1];
  Double_t z10=p0[2]-x0[2];
  return ((y10*z12-z10*y12)*(y10*z12-z10*y12)+
	  (z10*x12-x10*z12)*(z10*x12-x10*z12)+
	  (x10*y12-y10*x12)*(x10*y12-y10*x12))
    /(x12*x12+y12*y12+z12*z12);
}
//
void ComputeVertexNew(KMCProbeFwd &t0, KMCProbeFwd &t1, KMCProbeFwd &t2, Double_t zStart, Double_t &xV, Double_t &yV, Double_t &zV, Double_t &sigmaVert){
  
  t0.PropagateToZBxByBz(zStart);
  t1.PropagateToZBxByBz(zStart);
  t2.PropagateToZBxByBz(zStart);
  int knacc = 3;
  
  Double_t initPos[3]={0.,0.,0.};

  Double_t (*vectP0)[3]=new Double_t [knacc][3];
  Double_t (*vectP1)[3]=new Double_t [knacc][3];
  
  Double_t sum[3][3];
  Double_t dsum[3]={0,0,0};
  TMatrixD sumWi(3,3);
  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<3;j++){
      sum[i][j]=0;
      sumWi(i,j)=0.;
    }
  }

  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
  KMCProbeFwd* track;
  
  for(Int_t i=0; i<knacc; i++){
    if(i==0) track=&t0;
    else if(i==1) track=&t1;
    else if(i==2) track=&t2;
    
    Double_t p0[3],cd[3],sigmasq[3];
    Double_t wmat[9];
    for(Int_t iel=0; iel<9; iel++) wmat[iel]=0;
    track->GetXYZ(p0);
    track->GetPXYZ(cd);
    sigmasq[0]=track->GetSigmaX2()+ss;
    sigmasq[1]=track->GetSigmaY2()+ss;
    wmat[0]=1./sigmasq[0];
    wmat[5]=1./sigmasq[1];
    wmat[8]=1./sigmasq[2];
    TMatrixD wWi(3,3);
    Int_t iel=0;
    for(Int_t ia=0;ia<3;ia++){
      for(Int_t ib=0;ib<3;ib++){
	wWi(ia,ib)=wmat[iel];
	iel++;
      }    
    }

    sumWi+=wWi;

    Double_t p1[3]={p0[0]+cd[0],p0[1]+cd[1],p0[2]+cd[2]};
    vectP0[i][0]=p0[0];
    vectP0[i][1]=p0[1];
    vectP0[i][2]=p0[2];
    vectP1[i][0]=p1[0];
    vectP1[i][1]=p1[1];
    vectP1[i][2]=p1[2];
    
    Double_t matr[3][3];
    Double_t dknow[3];
    GetStrLinDerivMatrix(p0,p1,matr,dknow);

    for(Int_t iii=0;iii<3;iii++){
      dsum[iii]+=dknow[iii]; 
      for(Int_t lj=0;lj<3;lj++) sum[iii][lj]+=matr[iii][lj];
    }
  }
 
  TMatrixD invsumWi(TMatrixD::kInverted,sumWi);
  Double_t covmatrix[6];
  covmatrix[0] = invsumWi(0,0);
  covmatrix[1] = invsumWi(0,1);
  covmatrix[2] = invsumWi(1,1);
  covmatrix[3] = invsumWi(0,2);
  covmatrix[4] = invsumWi(1,2);
  covmatrix[5] = invsumWi(2,2);

  Double_t vett[3][3];
  Double_t det=GetDeterminant3X3(sum);
  Double_t sigma=0;
  
  if(TMath::Abs(det) > kAlmost0){
    for(Int_t zz=0;zz<3;zz++){
      for(Int_t ww=0;ww<3;ww++){
	for(Int_t kk=0;kk<3;kk++) vett[ww][kk]=sum[ww][kk];
      }
      for(Int_t kk=0;kk<3;kk++) vett[kk][zz]=dsum[kk];
      initPos[zz]=GetDeterminant3X3(vett)/det;
    }


    for(Int_t i=0; i<knacc; i++){
      Double_t p0[3]={0,0,0},p1[3]={0,0,0};
      for(Int_t ii=0;ii<3;ii++){
	p0[ii]=vectP0[i][ii];
	p1[ii]=vectP1[i][ii];
      }
      sigma+=GetStrLinMinDist(p0,p1,initPos);
    }

    if(sigma>0.) {sigma=TMath::Sqrt(sigma);}else{sigma=999;}
  }else{
    sigma=999;
  }
  xV=initPos[0];
  yV=initPos[1];
  zV=initPos[2];
  sigmaVert=sigma;
  delete [] vectP0;
  delete [] vectP1;
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

Double_t CosPointingAngleXY(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent)
{

  /// Cosine of pointing angle in transverse plane assuming it is produced
  /// at "point"

  TVector3 momXY(parent.Px(), parent.Py(), 0.);
  TVector3 flineXY(vsec[0] - vprim[0],
		   vsec[1] - vprim[1],
		   0.);

  Double_t ptot2 = momXY.Mag2() * flineXY.Mag2();
  if (ptot2 <= 0)
    {
      return 0.0;
    }
  else
    {
      Double_t cos = momXY.Dot(flineXY) / TMath::Sqrt(ptot2);
      if (cos > 1.0)
            cos = 1.0;
      if (cos < -1.0)
	cos = -1.0;
      return cos;
    }
}

Double_t CosPointingAngle(Double_t vprim[3], Double_t vsec[3], TLorentzVector &parent)
{

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
