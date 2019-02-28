class MagField;
class KMCDetectorFwd;

void DisplaySetupQuickAndDirty(KMCDetectorFwd *det);
void DisplaySetupTGeo(KMCDetectorFwd *det);

//=================================================================================

void DisplaySetup(Bool_t useTGeo=kFALSE, const char* setup="setup_ams_10pixel_wall_180.txt"){
    
  //  gROOT->LoadMacro("../src/loadFastSim.C");
  //  loadFastSim("../src"); 
  KMCDetectorFwd *det = new KMCDetectorFwd();
  det->ReadSetup(setup,setup);
  det->Print();  
  if (useTGeo) DisplaySetupTGeo(det);
  else DisplaySetupQuickAndDirty(det);
}

//=================================================================================

void DisplaySetupQuickAndDirty(KMCDetectorFwd *det) { 
  Int_t nAbso = 0; 

  Int_t colorBField = kYellow - 8;
  Int_t colorAbso   = kBlue - 7;
  Int_t colorITS    = kViolet + 6;
  Int_t colorMS     = kRed;
  Int_t colorTrig   = kMagenta;
  Int_t colorVertex = kBlack;
  
  // loop over all layers (active and passive); 
  TList *layerList = det->GetLayers();
  Int_t nlayer = layerList->GetEntries(); 
  TVector r(nlayer), z(nlayer); 
  
  // check layers to calculate the dimensions of the drawing
  
  for (Int_t il=0; il < nlayer; il++) {
    KMCLayerFwd *layer = (KMCLayerFwd*) layerList->At(il); 
    printf ("layer %d ActiveID = %d rmax = %g\n",il ,layer->GetActiveID(), layer->GetRMax());
    if (layer->GetActiveID()>0) r[il] = layer->GetRMax(); 
    else r[il] = 0; 
    z[il] = layer->GetZ(); 
  }

  Double_t dz = z.Max()-z.Min(); 
  
  // create canvas and put an empty histogram on top of which objects will be drawn

  Double_t xhistmin = z.Min()-0.05*dz, xhistmax =  z.Max() + 0.2*dz, yhistmin = -1.2*r.Max(), yhistmax = 1.2*r.Max();
  Double_t dx = xhistmax - xhistmin, dy = yhistmax - yhistmin;
  TCanvas *cSetup= new TCanvas("cSetup","cSetup",1000,1000*dy/dx); 
  TH2F *hEmpty = new TH2F("hEmpty","hEmpty",1000, xhistmin, xhistmax, 1000, yhistmin, yhistmax); 
  hEmpty->SetXTitle("z (cm)"); 
  hEmpty->SetYTitle("r (cm)"); 
  hEmpty->Draw(); 

  TLegend *legend = new TLegend(0.85,0.25,1,0.88);

  // display layers

  Int_t color = kGray; 
  for (Int_t il=0; il < nlayer; il++) {
    KMCLayerFwd* layer = (KMCLayerFwd*) layerList->At(il); 
    if (layer->IsVertex()) color = kBlack; 
    else if (layer->IsITS()) color = colorITS; 
    //    else if (layer->IsAbs()) color = colorAbso + nAbso++;
    else if (layer->IsAbs()) color = colorAbso;
    else if (layer->IsMS()) color = colorMS;
    else if (layer->IsTrig()) color = colorTrig;
    else if (layer->IsDummy()) color = kCyan;
    
    Double_t z1 = layer->GetZ() - layer->GetThickness()/2., z2 = layer->GetZ() + layer->GetThickness()/2.;
    Double_t radius = layer->GetRMax();
    if (radius > 1e6) radius = 1.2*r.Max(); 
    TBox *box = new TBox(z1, -radius, z2, radius);  
    box->SetLineColor(color);
    box->SetFillColor(color);
    //    if (layer->IsAbs()) box->SetFillStyle(3002);
    box->Draw(); 
    legend->AddEntry(box,layer->GetName(),"f");
  }

  // display magnetic field as colored areas
  
  MagField *mag= (MagField*) TGeoGlobalMagField::Instance()->GetField(); 
  
  int nBreg = mag->GetNReg();
  const double *zMinB=mag->GetZMin();
  const double *zMaxB=mag->GetZMax();
  Bool_t addBtoLegend = kTRUE; 
  for (Int_t ifield = 0; ifield < nBreg; ifield++) {
    const double *bvals=mag->GetBVals(ifield); 
    printf ("Magnetic field %d: \t zmin = %g \t zmax = %g \n",ifield, zMinB[ifield], zMaxB[ifield]);
    printf ("bvals: %g %g %g\n",bvals[0], bvals[1], bvals[2]); 
    TBox *box = new TBox(zMinB[ifield], hEmpty->GetYaxis()->GetXmin(), zMaxB[ifield], hEmpty->GetYaxis()->GetXmax());  
    box->SetLineColor(colorBField); 
    box->SetFillColor(colorBField); 
    box->SetFillStyle(3003);
    box->Draw(); 
    if (addBtoLegend) {
      legend->AddEntry(box,"B Field","f"); 
      addBtoLegend = kFALSE; 
    }
  }

  //  

  legend->Draw(); 
}

//=================================================================================

void DisplaySetupTGeo(KMCDetectorFwd *det){

  new TGeoManager("geom","NA60+ simplified geometry");
  TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
  
  Int_t nAbso = 0; 
  
  Int_t colorBField = kYellow - 10;
  Int_t colorAbso   = kBlue - 7;
  Int_t colorITS    = kViolet + 6;
  Int_t colorMS     = kRed;
  Int_t colorTrig   = kMagenta;
  Int_t colorVertex = kBlack;
  
  // loop over all layers (active and passive); 
  TList *layerList = det->GetLayers();
  Int_t nlayer = layerList->GetEntries(); 
  TVector r(nlayer), z(nlayer); 

  // check layers to calculate the dimensions of the drawing
  
  for (Int_t il=0; il < nlayer; il++) {
    KMCLayerFwd *layer = (KMCLayerFwd*) layerList->At(il); 
    printf ("layer %d ActiveID = %d rmax = %g\n",il ,layer->GetActiveID(), layer->GetRMax());
    if (layer->GetActiveID()>0) r[il] = layer->GetRMax(); 
    else r[il] = 0; 
    z[il] = layer->GetZ(); 
  }

  Double_t dx = r.Max(), dy = r.Max(), dz = z.Max()-z.Min(); 
  
  TGeoVolume *top=gGeoManager->MakeBox("top",Air, 1.1*dx, 1.1*dy, 1.1*dz/2.);
  gGeoManager->SetTopVolume(top);


  //  TCanvas *cSetup= new TCanvas("cSetup","cSetup",1000,1000*dy/dx); 
  
  MagField *mag= (MagField*) TGeoGlobalMagField::Instance()->GetField(); 
  
  // TLegend *legend = new TLegend(0.85,0.25,1,0.88);
  
  int nBreg = mag->GetNReg();
  const double *zMinB=mag->GetZMin();
  const double *zMaxB=mag->GetZMax();
  Bool_t addBtoLegend = kTRUE; 
  char name[100]; 
  for (Int_t ifield = 0; ifield < nBreg; ifield++) {
    Double_t zMean  = (zMaxB[ifield]+zMinB[ifield])/2.;
    Double_t dzHalf = (zMaxB[ifield]-zMinB[ifield])/2.;
    sprintf (name,"Mag_%d",ifield); 
      
    TGeoVolume *tube=gGeoManager->MakeTubs(name,Air,0,1.2*r.Max(),dzHalf,0,360);
    tube->SetFillColor(colorBField);
    tube->SetLineColor(colorBField);
    tube->SetTransparency(1); 
    top->AddNodeOverlap(tube,1,new TGeoCombiTrans(0,0,zMean,new TGeoRotation(name,0,180,0)));
    // if (addBtoLegend) {
    //   legend->AddEntry(box,"B Field","f"); 
    //   addBtoLegend = kFALSE; 
    // }
  }

  // display layers

  Int_t color = kGray; 
  for (Int_t il=0; il < nlayer; il++) {
    KMCLayerFwd* layer = (KMCLayerFwd*) layerList->At(il); 
    if (layer->IsVertex()) color = kBlack; 
    else if (layer->IsITS()) color = colorITS; 
    else if (layer->IsAbs()) color = colorAbso;
    // else if (layer->IsAbs()) color = colorAbso + nAbso++;
    else if (layer->IsMS()) color = colorMS;
    else if (layer->IsTrig()) color = colorTrig;
    else if (layer->IsDummy()) color = kCyan;
    Double_t z1 = layer->GetZ() - layer->GetThickness()/2., z2 = layer->GetZ() + layer->GetThickness()/2.;
    Double_t zMean = (z2+z1)/2., dzHalf = (z2-z1)/2.; 
    Double_t rmin = layer->GetRMin();
    Double_t rmax = layer->GetRMax();
    if (rmax > 1e6) rmax = 1.2*r.Max(); 

    TGeoVolume *tube=gGeoManager->MakeTubs(layer->GetName(),Air,rmin,rmax,dzHalf,0,360);
    tube->SetFillColor(color);
    tube->SetLineColor(color);
    tube->SetTransparency(1); 
    top->AddNodeOverlap(tube,1,new TGeoCombiTrans(0,0,zMean,new TGeoRotation(layer->GetName(),0,180,0)));
    // legend->AddEntry(box,layer->GetName(),"f");
  }
  //  gGeoManager->SetNsegments(4); 
  top->Draw("ogl");
  //  top->Draw("x3d");
  // legend->Draw(); 
}
