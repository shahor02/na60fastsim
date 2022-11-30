#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include "KMCProbeFwd.h"
#include "KMCDetectorFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "./HFUtils.C"
#endif

constexpr int kPDGofInterest = 1010020040;
constexpr double kMassOfInterest = 4.839961;
constexpr double kC = 0.0299792458;
constexpr double kTau = 274.;
constexpr double kCTau = kTau * kC;
constexpr double kHyperMass = (3727.380 + 1115.683 - 3.102) * 1.e-3;
constexpr double kDaughterMasses[] = {3.727380, 0.938272, 0.139570};
constexpr double kDaughterCharges[] = {2., 1., -1.};
constexpr int kNdaughters = 3;
constexpr int seed = 1;


struct ParticleProperty {
  int charge;
  double mass;
  double ctau;
  std::vector<float> decayProb;
  std::vector<std::vector<int>> decayProngs;
  std::vector<std::vector<int>> decayCharges; /// easier to deal with the TGenPhaseSpace later...
  std::vector<std::vector<double>> decayMasses; /// easier to deal with the TGenPhaseSpace later...
  bool isStable() { return decayProngs.empty(); }
};

struct Particle {
  Particle(int c, int pd, double m, TLorentzVector& p) : charge{c}, pdg{pd}, mass{m}, mom{p} {}
  int charge;
  int pdg;
  double mass;
  TLorentzVector mom;
};

std::unordered_map<int, ParticleProperty> particles{
    {11, {-1, 0.000510999, 1.e12, {}, {}, {}, {}}},
    {13, {-1, 0.105658, 1.e12, {}, {}, {}, {}}},
    {211, {1, 0.139570, 1.e12, {}, {}, {}, {}}},
    {310, {0, 0.497611, 2.6843417, {0.6920f, 1.f}, {{211, -211}, {}}, {{1, -1}, {}}, {{0.139570, 0.139570}, {}}}},
    {321, {1, 0.493677, 1.e12, {}, {}, {}, {}}},
    {2212, {1, 0.938272, 1.e12, {}, {}, {}, {}}},
    {3122, {0, 1.115683, kCTau, {0.639f, 1.f}, {{2212, -211}, {}}, {{1, -1}, {}}, {{0.938272, 0.139570}, {}}}},
    {1000010020, {1, 1.87561, 1.e12, {}, {}, {}, {}}},
    {1000010030, {1, 2.80892, 1.e12, {}, {}, {}, {}}},
    {1000020030, {2, 2.80839, 1.e12, {}, {}, {}, {}}},
    {1000020040, {2, 3.72738, 1.e12, {}, {}, {}, {}}},
    {1010010030, {1, 2.991134, kCTau, {0.25f, 0.6f, 1.}, {{1000020030, -211}, {1000010020, 2212, -211}, {}}, {{2, -1}, {1, 1, -1}, {}}, {{2.80839, 0.139570}, {1.87561, 0.938272, 0.139570}, {}}}},
    {1010010040, {1, 3.922434, kCTau, {0.5f, 1.}, {{1000020040, -211}, {}}, {{2, -1}, {}}, {{3.72738, 0.139570}, {}}}},
    {1010020040, {1, 3.921728, kCTau, {0.5f, 1.}, {{1000020030, 2212, -211}, {}}, {{2, 1, -1}, {}}, {{2.80839, 0.938272, 0.139570}, {}}}}};

struct MiniH
{
  float px, py, pz;
  float x, y, z;
  float m;
  float ProngPvDCAx[kNdaughters];
  float ProngPvDCAy[kNdaughters];
  float pxMC, pyMC, pzMC;
  float xMC, yMC, zMC;
  unsigned char NClusters[kNdaughters];
  char reconstructed = 0;
  char fake = 0;
};

struct MiniBkg
{
  float m;
  float pt;
  float eta;
  float phi;
  float x, y, z;
  float ProngPvDCAx[kNdaughters];
  float ProngPvDCAy[kNdaughters];
  float cosPA;
  float ct;
};

class Decayer : public TGenPhaseSpace
{
public:
  Decayer() = default;
  double maxZdecay = 30.;
  std::vector<std::pair<std::array<float, 3>, Particle>> decayParticle(int pdg, const ParticleProperty &part, std::array<float, 3> &vtx, TLorentzVector &mom)
  {
    double prob{gRandom->Uniform()};
    // std::cout << "\t\t prob:" << prob << std::endl;
    int decay{-1};
    for (unsigned int iD{0}; iD < part.decayProb.size(); ++iD)
    {
      if (prob < part.decayProb[iD])
      {
        decay = iD;
        break;
      }
    }
    // std::cout << "\t\t decay:" << decay << std::endl;
    if (decay < 0 || part.decayProngs[decay].empty())
    {
      return {};
    }
    std::array<float, 3> xyz;
    double decayLen = gRandom->Exp(kCTau) * mom.Beta() * mom.Gamma();
    xyz[2] = vtx[2] + decayLen * mom.Pz() / mom.P();
    if (xyz[2] > maxZdecay) {
      return {};
    }
    xyz[0] = vtx[0] + decayLen * mom.Pt() * std::cos(mom.Phi()) / mom.P();
    xyz[1] = vtx[1] + decayLen * mom.Pt() * std::sin(mom.Phi()) / mom.P();
    SetDecay(mom, part.decayMasses[decay].size(), part.decayMasses[decay].data());
    Generate();
    std::vector<std::pair<std::array<float, 3>, Particle>> outVec;
    for (unsigned int d{0}; d < part.decayMasses[decay].size(); ++d)
    {
      outVec.push_back(std::make_pair(xyz, Particle((1 - 2 * std::signbit(pdg)) * part.decayCharges[decay][d], part.decayProngs[decay][d], part.decayMasses[decay][d], *GetDecay(d))));
    }
    return outVec;
  }
};

// Track cuts
double ChiTot = 1.5;
int minITShits = 4;

KMCDetectorFwd *createDetectorFwd(double Eint, const char *setup)
{
  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup, setup);
  det->InitBkg(Eint);
  det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VT

  det->SetMinITSHits(TMath::Min(minITShits, det->GetNumberOfActiveLayersITS())); // NA60+
  // det->SetMinITSHits(det->GetNumberOfActiveLayersITS()-1); //NA60
  det->SetMinMSHits(0); // NA60+
  // det->SetMinMSHits(det->GetNumberOfActiveLayersMS()-1); //NA60
  det->SetMinTRHits(0);
  //
  // max number of seeds on each layer to propagate (per muon track)
  det->SetMaxSeedToPropagate(3000);
  //
  // set chi2 cuts
  det->SetMaxChi2Cl(10.);  // max track to cluster chi2
  det->SetMaxChi2NDF(3.5); // max total chi2/ndf
  det->SetMaxChi2Vtx(20);  // fiducial cut on chi2 of convergence to vtx
  // IMPORTANT FOR NON-UNIFORM FIELDS
  det->SetDefStepAir(1);
  det->SetMinP2Propagate(1); // NA60+
  // det->SetMinP2Propagate(2); //NA60
  //
  det->SetIncludeVertex(kFALSE); // count vertex as an extra measured point
  //  det->SetApplyBransonPCorrection();
  det->ImposeVertex(0., 0., 0.);
  //

  TVirtualMagField *fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class())
  {
    MagField *mag = (MagField *)fld;
    int BNreg = mag->GetNReg();
    const double *BzMin = mag->GetZMin();
    const double *BzMax = mag->GetZMax();
    const double *BVal;
    printf("*************************************\n");
    printf("number of magnetic field regions = %d\n", BNreg);
    for (int i = 0; i < BNreg; i++)
    {
      BVal = mag->GetBVals(i);
      printf("*** Field region %d ***\n", i);
      if (i == 0)
      {
        printf("Bx = %f B = %f Bz = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
      else if (i == 1)
      {
        printf("B = %f Rmin = %f Rmax = %f zmin = %f zmax = %f\n", BVal[0], BVal[1], BVal[2], BzMin[i], BzMax[i]);
      }
    }
  }
  return det;
}

void runHyperNucleiCandidates(int nevents = 14000000,
                              double Eint = 160.,
                              const char *setup = "setup.txt",
                              const char *filNamPow = "AnalysisResults.root",
                              bool writeNtuple = true,
                              bool simulateBg = true)
{
  int refreshBg = 10000;
  gRandom->SetSeed(seed);

  // MUSIC input
  printf("--> pt and y shape of the hypernuclei from %s\n", filNamPow);
  TFile *filPow = new TFile(filNamPow);
  TH3D *ptYetaShape = (TH3D *)filPow->Get("ptYeta");
  TH1D *ptShape = (TH1D *)ptYetaShape->ProjectionX("ptShape");
  TH1D *yShape = (TH1D *)ptYetaShape->ProjectionY("yShape");

  KMCDetectorFwd *det = createDetectorFwd(Eint, setup);

  // prepare decays
  MiniH hyper;
  TGenPhaseSpace decay;
  TLorentzVector parentgen, daugen[3], parent, daurec[3];
  KMCProbeFwd recProbe[3];

  TFile fnt("Hypernuclei-Signal-ntuple.root", "recreate");
  TTree tree("hyper", "hyper");
  tree.Branch("h", &hyper);

  TH1D hDCAx("hDCAx", ";DCA_{x} (cm);", 200, -0.1, 0.1);
  TH1D hDCAy("hDCAy", ";DCA_{y} (cm);", 200, -0.1, 0.1);
  for (int iev = 0; iev < nevents; iev++)
  {
    double vprim[3] = {0, 0, 0};
    if (iev % 100 == 0)
      printf(" ***************  ev = %d \n", iev);
    double pxyz[3];

    if (simulateBg && (iev % refreshBg) == 0)
      det->GenBgEvent(0., 0., 0.);
    double ptGenD = ptShape->GetRandom(); // get Lc distribution from file
    double yGenD = yShape->GetRandom();
    double phi = gRandom->Rndm() * 2 * TMath::Pi();
    hyper.pxMC = ptGenD * std::cos(phi);
    hyper.pyMC = ptGenD * std::sin(phi);

    double mt = std::hypot(ptGenD, kHyperMass);
    hyper.pzMC = mt * TMath::SinH(yGenD);
    double totGenMom{std::hypot(ptGenD, hyper.pzMC)};
    double en = mt * TMath::CosH(yGenD);
    hyper.reconstructed = 0;
    hyper.fake = 0;

    parentgen.SetPxPyPzE(hyper.pxMC, hyper.pyMC, hyper.pzMC, en);

    double decayLenght = gRandom->Exp(kCTau) * parentgen.Beta() * parentgen.Gamma();
    hyper.zMC = decayLenght * hyper.pzMC / totGenMom;
    hyper.xMC = decayLenght * ptGenD * std::cos(phi) / totGenMom;
    hyper.yMC = decayLenght * ptGenD * std::sin(phi) / totGenMom;
    decay.SetDecay(parentgen, 3, kDaughterMasses);
    decay.Generate();

    for (int iD{0}; iD < kNdaughters; ++iD)
    {
      daugen[iD] = *(decay.GetDecay(iD));
      if (!det->SolveSingleTrack(daugen[iD].Pt(), daugen[iD].Rapidity(), daugen[iD].Phi(), kDaughterMasses[iD], kDaughterCharges[iD], hyper.xMC, hyper.yMC, hyper.zMC, 0, 1, 99))
        continue;
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw || trw->GetNormChi2(kTRUE) > ChiTot)
        continue;
      hyper.reconstructed++;

      hyper.NClusters[iD] = trw->GetNITSHits();
      hyper.fake += trw->GetNFakeITSHits();
      recProbe[iD] = *trw;
    }

    if (hyper.reconstructed < 3)
    {
      tree.Fill();
      continue;
    }

    for (int iD{0}; iD < kNdaughters; ++iD)
    {
      hyper.ProngPvDCAx[iD] = recProbe[iD].GetX();
      hyper.ProngPvDCAy[iD] = recProbe[iD].GetY();
      if (iD)
      {
        hDCAx.Fill(recProbe[iD].GetX());
        hDCAy.Fill(recProbe[iD].GetY());
      }
    }

    double xV, yV, zV;
    double sigmaVert;
    ComputeVertex(recProbe[0], recProbe[1], recProbe[2], vprim[2], xV, yV, zV, sigmaVert);
    hyper.x = xV;
    hyper.y = yV;
    hyper.z = zV;

    // Get daughter track momentum at decay vertex
    for (int idau = 0; idau < kNdaughters; idau++)
    {
      recProbe[idau].PropagateToZBxByBz(zV);
      recProbe[idau].GetPXYZ(pxyz);
      daurec[idau].SetXYZM(pxyz[0], pxyz[1], pxyz[2], kDaughterMasses[idau]);
    }

    parent = daurec[0];
    parentgen = daugen[0];
    for (int iD = 1; iD < kNdaughters; ++iD)
    {
      parent += daurec[iD];
      parentgen += daugen[iD];
    }

    hyper.px = parent.Px();
    hyper.py = parent.Py();
    hyper.pz = parent.Pz();
    hyper.m = parent.M();

    tree.Fill();
  }

  fnt.cd();
  hDCAx.Write();
  hDCAy.Write();
  tree.Write();
  fnt.Close();
}

void runHyperNucleiBackgroundCandidates(int nevents = 14000000,
                                        double Eint = 160.,
                                        const char *setup = "setup.txt",
                                        const char *filNamPow = "output.root",
                                        bool writeNtuple = true,
                                        bool simulateBg = true)
{

  float momentum{50.7 * 82.f};
  float energy{193.f + std::hypot(193.f, momentum)};
  float rapidity{0.5f * std::log((energy + momentum) / (energy - momentum))};
  float sqrtS{std::sqrt(energy * energy - momentum * momentum)};

  int refreshBg = 10000;
  gRandom->SetSeed(seed);

  // MUSIC input
  printf("--> pt and y shape of the hypernuclei from %s\n", filNamPow);
  TFile *filPow = new TFile(filNamPow);
  TTree *kineTree = (TTree *)filPow->Get("kinematics");
  int pdg{0}, event{0};
  float p0{0.f}, px{0.f}, py{0.f}, pz{0.f};
  kineTree->SetBranchAddress("event", &event);
  kineTree->SetBranchAddress("pdg", &pdg);
  kineTree->SetBranchAddress("p0", &p0);
  kineTree->SetBranchAddress("px", &px);
  kineTree->SetBranchAddress("py", &py);
  kineTree->SetBranchAddress("pz", &pz);

  KMCDetectorFwd *det = createDetectorFwd(Eint, setup);

  // prepare decays
  MiniBkg hyper;

  TFile fnt("Hypernuclei-Background-ntuple.root", "recreate");
  TTree tree("hyper", "hyper");
  tree.Branch("h", &hyper);

  int inputEntry{0};
  bool newEvent{false};
  TH1D hDCAx("hDCAx", ";DCA_{x} (cm);", 200, -0.1, 0.1);
  TH1D hDCAy("hDCAy", ";DCA_{y} (cm);", 200, -0.1, 0.1);
  TH1D hDCA("hDCA", ";DCA (cm);", 200, -0.1, 0.1);
  Decayer decayer;
  for (int iev = 0; iev < nevents; iev++)
  {
    if (iev % 10 == 0) {
      std::cout << "Processing event " << iev << "\r" << std::flush;
    }
    std::vector<std::pair<std::array<float, 3>, Particle>> kine;
    std::array<float, 3> pvtx{0.f, 0.f, 0.f};
    while (inputEntry < kineTree->GetEntries())
    {
      if (newEvent)
      {
        newEvent = event != iev;
      }
      else
      {
        kineTree->GetEntry(inputEntry++);
      }
      if (event != iev)
      {
        newEvent = true;
        break;
      }
      if (std::abs(pdg) == kPDGofInterest || particles.find(std::abs(pdg)) == particles.end()) {
        continue;
      }
      auto& prop = particles[std::abs(pdg)];
      TLorentzVector mom(px, py, pz, p0);
      if (prop.isStable()) {
        kine.push_back(std::make_pair(pvtx, Particle((1 - 2 * std::signbit(pdg)) * prop.charge, pdg, prop.mass, mom)));
      } else {
        auto out = decayer.decayParticle(pdg, prop, pvtx, mom);
        for (auto& k : out) {
          kine.push_back(k);
        }
      }
    }
    std::vector<KMCProbeFwd> nuclei[2], tracks[2];
    int nTracks{0};
    for (auto& part : kine) {
      auto& p = part.second;
      if (!det->SolveSingleTrack(p.mom.Pt(), p.mom.Rapidity() + rapidity, p.mom.Phi(), p.mass, p.charge, part.first[0], part.first[1], part.first[2], 0, 1, 99))
        continue;
      KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
      if (!trw || trw->GetNormChi2(kTRUE) > ChiTot)
        continue;

      hDCAx.Fill(trw->GetX());
      hDCAy.Fill(trw->GetY());
      hDCA.Fill(std::hypot(trw->GetY(), trw->GetX()) * (1.f - 2.f * (trw->GetX() < 0)));
      if (std::abs(trw->GetY()) < 80.e-4 || std::abs(trw->GetX()) < 60.e-4) {
        continue;
      }
      if (std::abs(p.charge) == 1) {
        tracks[std::signbit(p.charge)].push_back(*trw);
      } else {
        nuclei[std::signbit(p.charge)].push_back(*trw);
      }
      nTracks++;
    }


    TLorentzVector daurec[3], tot;
    KMCProbeFwd* recProbe[3];
    double pxyz[3];
    printf("%zu, %zu, %zu, %zu\n", nuclei[0].size(), nuclei[1].size(), tracks[0].size(), tracks[1].size());
    for (int iCharge{0}; iCharge < 2; ++iCharge)
    {
      for (auto &nucleus : nuclei[iCharge])
      {
        recProbe[0] = &nucleus;
        hyper.ProngPvDCAx[0] = nucleus.GetX();
        hyper.ProngPvDCAy[0] = nucleus.GetY();
        for (auto& track1 : tracks[iCharge])
        {
          recProbe[1] = &track1;
          hyper.ProngPvDCAx[1] = track1.GetX();
          hyper.ProngPvDCAy[1] = track1.GetY();
          for (auto &track2 : tracks[!iCharge])
          {
            recProbe[2] = &track2;
            hyper.ProngPvDCAx[2] = track1.GetX();
            hyper.ProngPvDCAy[2] = track1.GetY();
            double xV, yV, zV;
            double sigmaVert;
            ComputeVertex(nucleus, track1, track2, 0, xV, yV, zV, sigmaVert);
            hyper.x = xV;
            hyper.y = yV;
            hyper.z = zV;

            // Get daughter track momentum at decay vertex
            for (int idau = 0; idau < kNdaughters; idau++)
            {
              recProbe[idau]->PropagateToZBxByBz(zV);
              recProbe[idau]->GetPXYZ(pxyz);
              daurec[idau].SetXYZM(pxyz[0], pxyz[1], pxyz[2], kDaughterMasses[idau]);
            }
            tot = daurec[0] + daurec[1] + daurec[2];
            hyper.pt = tot.Pt();
            hyper.eta = tot.Eta();
            hyper.phi = tot.Phi();
            hyper.cosPA = (tot.Px() * xV + tot.Py() * yV + tot.Pz() * zV) / std::hypot(xV, yV ,zV) / tot.P();
            hyper.ct = std::hypot(xV, yV, zV) * kMassOfInterest / tot.P();
            hyper.m = tot.M();
            tree.Fill();
          }
        }
      }
    }
  }

  fnt.cd();
  hDCAx.Write();
  hDCAy.Write();
  hDCA.Write();
  tree.Write();
  fnt.Close();
}
