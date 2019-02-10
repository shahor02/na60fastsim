#include "KMCUtils.h"
#include <TMath.h>
//---------------------------------

const double MagField::fZMin[MagField::kNReg] = {0., 380.};// cm, cm

const double MagField::fZMax[MagField::kNReg] = {40.,680.};//cm cm
const double MagField::fBVal[MagField::kNReg][3] = {{-30,0,0},{1.e2,30.,300.}}; //
//const double MagField::fZMax[MagField::kNReg] = {40.,830.};//cm cm 



//const double MagField::fBVal[MagField::kNReg][3] = {{-30,0,0},{1.e2,30.,160.}};




//__________________________________________________
void MagField::Field(const Double_t *xyz, Double_t *bxyz) 
{
  bxyz[0]=bxyz[1]=bxyz[2]=0.;
  double R = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
  for (int ir=0;ir<kNReg;ir++) {
    if((xyz[2] >fZMin[ir]) && (xyz[2] < fZMax[ir])) {
      if(ir == 0) for (int i=3;i--;) bxyz[i] = fBVal[ir][i];
      //    if(ir == 2) for (int i=3;i--;) bxyz[i] = fBVal[ir][i];
      if(ir == 1){
	//	double R = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
	if(R > fBVal[ir][1] && R < fBVal[ir][2]){
	  // Atlas/Chorus toroid
	  //	  bxyz[0] = -fBVal[ir][0] * xyz[1]/R;
	  //	  bxyz[1] =  fBVal[ir][0] * xyz[0]/R;
          //bxyz[0] = -fBVal[ir][0] * xyz[1]/R - 0.001 * xyz[1]/R/R;
          //bxyz[1] =  fBVal[ir][0] * xyz[0]/R + 0.001 * xyz[0]/R/R;

	  // ACM toroid
	  // direct polarity
	  bxyz[0] = -fBVal[ir][0] * xyz[1]/R/R;
	  bxyz[1] =  fBVal[ir][0] * xyz[0]/R/R;
	  // reverse polarity
 	  //bxyz[0] =  fBVal[ir][0] * xyz[1]/R/R;
 	  //bxyz[1] = -fBVal[ir][0] * xyz[0]/R/R;
	  bxyz[2] = 0.;
	}
      }
      //if(R<400.) 
      //      printf("ir = %d z=%f Bx=%f By=%f Bz=%f\n",ir,xyz[2],bxyz[0],bxyz[1],bxyz[2]);
      //if (xyz[2]<0) printf("ir = %d z=%f Bx=%f By=%f Bz=%f\n",ir,xyz[2],bxyz[0],bxyz[1],bxyz[2]);
      return;
    }
  }
  //  printf("z=%f Bx=%f By=%f Bz=%f\n",xyz[2],bxyz[0],bxyz[1],bxyz[2]);
}
