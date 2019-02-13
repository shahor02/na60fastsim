#define gFortran
#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

#include <stdio.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>

// For GenMUONLMR
#include "GenMUONLMR.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"


#ifndef WIN32
#define dimugenlmr dimugenlmr_
#else
#define dimugenlmr DIMUGENLMR
#endif

// settings for signal generation
double MotherMass;
int generYPtPar, ProcType;
double y0SG;   // gaussian y mean
double sigySG = 1.;   // .. sigma
double yminSG = -10.;   // min y to generate 
double ymaxSG = 10.;   // 
double TSG;            // inv.slope of thermal pt distribution
double ptminSG = 0.01;
double ptmaxSG = 3;

static TH1F *hPtGen = 0;
static TH1F *hYGen = 0;
static TH2F *hYPtGen = 0;
static TH1F *hMassGen = 0;

static TRandom *rgen = 0;
static GenMUONLMR *gener = 0;
static TParticle *fMu[2]; 
int ifirstgen=0;

extern "C" {

void SetProcessParameters(const Char_t *Process, Double_t E);
void BookHistos();

void dimugenlmr(double &px1, double &py1, double &pz1, double &px2, double &py2, double &pz2){

// Choose here the signal process to generate
const Char_t *Process = "omegaDalitz";
// Process can be "etaDalitz", "eta2Body", "rho", "omega2Body", "omegaDalitz", "phi", "etaPrime", "Jpsi"

//Choose here the collision energy in the lab
double Eint=40.; 
//

int Seed=12345;

if(ifirstgen==0){

  rgen = new TRandom();
  rgen->SetSeed(Seed);
  printf("Initialized TRandom\n");

  SetProcessParameters(Process,Eint);

  gener = new GenMUONLMR(7,1); //N.B. first argument has no meaning when second (kLowEnergy) is set!
  printf("Generator class created\n");
  
//      0         1        2         3         4        5           6         7
//  "fPtPion","fPtKaon","fPtEta","fPtRho","fPtOmega","fPtPhi","fPtEtaPrime"  "fPtJPsi"
// 2 = rho
// 3 = omega 2 body
// 4 = omega Dalitz
//
// Processes eta 2B=0  eta D=1  rho =2   omega 2B=3  omega D=4  phi=5    eta p=6       pi=7      K=8      J/psi=9
//           kEtaLMR   kEtaLMR  kRhoLMR  kOmegaLMR   kOmegaLMR  kPhiLMR  kEtaPrimeLMR  kPionLMR  kK        kJPsi
//	      
//  BR       5.8e-6    3.1e-4   4.55e-5  7.28e-5     1.3e-4     2.86e-4  1.04e-4       1         0.6344    0.05

  gener->SetYParams(generYPtPar,1.,y0SG,sigySG,0.);
  gener->SetPtParams(generYPtPar,1.,TSG,MotherMass,0.);

  gener->GenerateSingleProcess(ProcType);
  printf("LMR generator initialization completed\n"); 

  BookHistos();

  ifirstgen++;
}

    gener->Generate();
    
    fMu[0] = gener->GetMuon(0);
    fMu[1] = gener->GetMuon(1);
    px1=fMu[0]->Px();
    py1=fMu[0]->Py();
    pz1=fMu[0]->Pz();
    px2=fMu[1]->Px();
    py2=fMu[1]->Py();
    pz2=fMu[1]->Pz();
    
    printf("px0=%f py0=%f pz0=%f e0=%f m0=%f\n", fMu[0]->Px(), fMu[0]->Py(), fMu[0]->Pz(), fMu[0]->Energy(), fMu[0]->GetMass());
    printf("px1=%f py1=%f pz1=%f e1=%f m1=%f\n", fMu[1]->Px(), fMu[1]->Py(), fMu[1]->Pz(), fMu[1]->Energy(),fMu[1]->GetMass());

    TLorentzVector parentgen, mugen[2];

    for (int imu=0;imu<2;imu++) {
      mugen[imu].SetXYZM(fMu[imu]->Px(),fMu[imu]->Py(),fMu[imu]->Pz(),0.105658369);//eli
    }
    parentgen  = mugen[0];
    parentgen += mugen[1];

    double y = parentgen.Rapidity();
    double pt = parentgen.Pt();
    double m = parentgen.M();
   
    hYGen->Fill(y);
    hPtGen->Fill(pt);
    hYPtGen->Fill(y,pt);
    hMassGen->Fill(m);

}

void SetProcessParameters(const Char_t *Process, Double_t E){
  
  if(E == 20.) {
    y0SG = 1.9;
    sigySG = 1.;
  } else if(E == 30.){
    y0SG = 2.08;
    sigySG = 1.;
  } else if(E == 40.){
    y0SG = 2.22;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.19; //PYTHIA 50 GeV
    else sigySG = 1.; 
  }  else if(E == 60.){
    y0SG = 2.42;
    sigySG = 1.;
  }  else if(E == 80.){
    y0SG = 2.57;
    sigySG = 1.;
  }  else if(E == 160.){
    y0SG = 2.9;
    if(strcmp(Process,"jpsi") == 0) sigySG = 0.42; //PYTHIA
    else sigySG = 1.2; //NA49
  }  else if(E == 400.){
    y0SG = 3.37;
    sigySG = 1.;
  }


  if(strcmp(Process,"etaDalitz") == 0) { // eta Dalitz
    ProcType    = 1;
    generYPtPar = 2;
    MotherMass  = 0.549;  
    if(E == 40.) {
      TSG = 0.225; 
    }
    else if(E == 160.){
      TSG = 0.24; 
    }
  } else   if(strcmp(Process,"eta2Body") == 0) { // eta Dalitz
    ProcType    = 0;
    generYPtPar = 2;
    MotherMass  = 0.549;  
    if(E == 40.) {
      TSG = 0.225; 
    }
    else if(E == 160.){
      TSG = 0.24; 
    }
  } else if(strcmp(Process,"rho") == 0) { // rho
      ProcType    = 2;
      generYPtPar = 3; 
      MotherMass  = 0.775; 
      if(E == 40.) {
	TSG = 0.25; 
      }
      else if(E == 160.){
	TSG = 0.29; 
      }
  }
  else if(strcmp(Process,"omega2Body") == 0) { // omega 2 body
    ProcType    = 3;
    generYPtPar = 4; 
    MotherMass  = 0.781;  
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
	TSG = 0.29; 
    }
  } 
  else if(strcmp(Process,"omegaDalitz") == 0) { // omega Dalitz
    ProcType    = 4;
    generYPtPar = 4; 
    MotherMass  = 0.781;
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
      TSG = 0.29; 
    }   
  } 
  else if(strcmp(Process,"phi") == 0) { // phi
    ProcType    = 5;
    generYPtPar = 5; 
    MotherMass  = 1.02; 
    if(E == 40.) {
      TSG =0.25;           
    }
    else if(E == 160.){
      TSG = 0.30; //NA49
    }  
  }   else if(strcmp(Process,"etaPrime") == 0) { // eta prime
    ProcType    = 6;
    generYPtPar = 6; 
    MotherMass  = 0.958;  
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
      TSG = 0.29; 
    } 
  }  else if(strcmp(Process,"jpsi") == 0) { // J/psi
    ProcType    = 9;
    generYPtPar = 7; 
    MotherMass  = 3.0969; 
    if(E == 40.) {
      TSG = 0.25; 
    }
    else if(E == 160.){
      TSG = 0.284; //NA50
    }   
  }

  printf("Generating %s\n",Process);
  printf("Particle mass = %f\n",MotherMass);
  printf("Inverse slope = %f\n",TSG);
  printf("rapidity sigma = %f\n",sigySG);

}

 void BookHistos(){
 
   // Book some histos
 
   hPtGen = new TH1F("hPTGen","Pt gen",30,ptminSG,ptmaxSG);
   hYGen = new TH1F("hYGen","Y full phase space",100.,yminSG,ymaxSG);
   hYPtGen = new TH2F("hYPTGen","Y-Pt corr match",80,1.0,5.4,30,ptminSG,ptmaxSG);
   hMassGen = new TH1F("hMassGen","Mass gen",100,0.,4.);
  
 } 

}

#ifndef WIN32
#define printroothistos printroothistos_
#else
#define printroothistos PRINTROOTHISTOS
#endif

extern "C" {

void printroothistos(){

  TFile *fout = new TFile("roothistos.root","recreate");
  hYPtGen->Write();
  hYGen->Write();
  hPtGen->Write();
  hMassGen->Write();
  fout->Close();

} 

}
