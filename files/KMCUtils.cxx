#include "KMCUtils.h"
#include <TMath.h>
#ifndef _NOALIROOT_
#include "AliLog.h"
#else
#include "LocLog.h"
using AliLog = LocLog;
#endif

//---------------------------------
ClassImp(MagField)


MagField::MagField(UInt_t id) {  SetUniqueID(id);  }

//__________________________________________________
void  MagField::SetMagVTParam(const std::string& paramString, const std::vector<float>& parVals, float zref)
{
  fFieldFunReg[0] = new TF1("dipVT", paramString.c_str(), fZMin[0], fZMax[0]);
  for (int i=0;i<(int)parVals.size();i++) {
    fFieldFunReg[0]->SetParameter(i, parVals[i]);
  }
  fZRef[0] = (zref<-1e5) ? (fZMin[0]+fZMax[0])/2 : zref;
  printf("Set parameterized field for VerTel dipole, reference Z = %f\n", fZRef[0]);
  fFieldFunReg[0]->Print();
}

//__________________________________________________
void  MagField::SetMagMSParam(const std::string& paramString, const std::vector<float>& parVals, float zref)
{
  fFieldFunReg[1] = new TF1("dipMS", paramString.c_str(), fZMin[1], fZMax[1]);
  for (int i=0;i<(int)parVals.size();i++) {
    fFieldFunReg[1]->SetParameter(i, parVals[i]);
  }
  fZRef[1] = (zref<-1e5) ? (fZMin[1]+fZMax[1])/2 : zref;
  printf("Set parameterized field for MuSpec dipole, reference Z = %f\n", fZRef[1]);
  fFieldFunReg[1]->Print();
}

//__________________________________________________
void MagField::Field(const Double_t *xyz, Double_t *bxyz) 
{
  bxyz[0]=bxyz[1]=bxyz[2]=0.;
  for (int ir=0;ir<kNReg;ir++) {
    if((xyz[2] >fZMin[ir]) && (xyz[2] < fZMax[ir])) {
      if (fFieldFunReg[ir]) {
	bxyz[0] = fFieldFunReg[ir]->Eval(xyz[2] - fZRef[ir]);
	return;
      }
      if(ir == 0) {
	for (int i=3;i--;) bxyz[i] = fBVal[ir][i];
      }
      else {
	if (GetUniqueID() == 0) {
	  for (int i=3;i--;) bxyz[i] = fBVal[ir][i];
	} else if (GetUniqueID() == 1) {
	  double R = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
	  if(R > fBVal[ir][1] && R < fBVal[ir][2]){
	    bxyz[0] = -fBVal[ir][0] * xyz[1]/R/R;
	    bxyz[1] =  fBVal[ir][0] * xyz[0]/R/R;
	    bxyz[2] = 0.;
	  }
	} else {
	  AliFatal(Form("Uknown field ID %d", (int)GetUniqueID()));
	}
      }
      return;
    }
  }
  //  printf("z=%f Bx=%f By=%f Bz=%f\n",xyz[2],bxyz[0],bxyz[1],bxyz[2]);
}
