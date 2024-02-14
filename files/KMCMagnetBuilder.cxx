#include "KMCMagnetBuilder.h"
#include "KMCDetectorFwd.h"
#include "KMCPolyLayer.h"

// see vf1.csv / print.pdf for details
const int nvCyan = 4;
float dataCyan[2][nvCyan] = {
  { 0.0000,  0.0000, 188.4000, 188.4000},
  {15.7200, 20.7200, 20.7200, 15.7200}};
float dataCyanO[2][nvCyan] = {
  { 0.0000,  0.0000, 340.,    340.},
  {15.7200, 20.7200, 20.7200, 15.7200}};

const int nvGreen = 9;
float dataGreen[2][nvGreen] = {
  {  5.240,  5.2400,  51.7000, 54.201,  61.3400,  60.0000, 51.7000, 50.0000, 45.7049},
  {10.4800, 15.4800,  15.4800, 15.720,  15.7200,  14.800 , 11.0000, 10.6400, 10.5629}};

const int nvBlue = 12;
float dataBlue[2][nvBlue] = {
  {10.4500, 10.4500, 45.7000,  50.0000, 51.7000,  60.0000, 61.3400, 67.5800, 60.000,  51.7000, 50.0000, 45.700 },
  { 5.2400, 10.2400, 10.2400,  10.6400, 11.0000,  14.8000, 15.7800, 15.7800,  9.5600,  5.8500,  5.6000,  5.2400}};

const int nvMagenta = 14;
float dataMagenta[2][nvMagenta] = {
  {15.7600, 15.7600, 45.7000, 50.000,  51.7000, 60.0000,   67.570, 84.000 , 78.0000, 71.3670, 63.4800, 60.0000, 51.7000,  45.700},
  {0.0000,  5.0000,  5.0000,  5.4000,   5.8500,  9.5600,  15.7500, 15.7500,  15.000, 12.2430,  5.9300, 3.7000 ,  0.7300,  0.0000}};


KMCMagnetBuilder::KMCMagnetBuilder(KMCDetectorFwd* det, float zIn, float dZ,
                                   float matRadL, float matDens, float defRadL, float defDens)
{
  const float barSect = 5.0;
  const float spacing = 0.24;

  int nSectors = 8;
  float phi0 = 0;
  
  float x2x0 = (matRadL>0 && matDens>0) ? barSect/(matRadL/matDens) : 0;
  float xrho = matDens*barSect;
  
  float zl = zIn + barSect/2;
  KMCPolyLayer* lrIn0 = det->AddPolyLayer("mag","magIn0", zl, defRadL, defDens, barSect);
  lrIn0->setXYOffsets(8.28, -15.72);
  auto& pin00 = lrIn0->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin00.poly[0].x = 52;
  pin00.poly[1].x = 52;
  pin00.poly[2].x = pin00.poly[3].x = 188.4; // fix edges
  lrIn0->setNSectorsPhiStart(nSectors, phi0);
  
  //-----------------
  zl += spacing + barSect;
  KMCPolyLayer* lrIn1 = det->AddPolyLayer("mag","magIn1", zl, defRadL, defDens, barSect);
  lrIn1->setXYOffsets(8.28, -15.72);
  auto& pin10 = lrIn1->addPolygon(nvMagenta, dataMagenta[0],dataMagenta[1], x2x0, xrho);
  auto& pin11 = lrIn1->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin11.poly[0].x = 67.570;
  pin11.poly[1].x = 78.0000;  
  pin11.poly[2].x = pin11.poly[3].x = 188.4; // fix edges
  lrIn1->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn2 = det->AddPolyLayer("mag","magIn2", zl, defRadL, defDens, barSect);
  lrIn2->setXYOffsets(8.28, -15.72);
  auto& pin20 = lrIn2->addPolygon(nvBlue, dataBlue[0],dataBlue[1], x2x0, xrho);
  auto& pin21 = lrIn2->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin21.poly[0].x = 61.3400;
  pin21.poly[1].x = 67.5800;  
  lrIn2->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn3 = det->AddPolyLayer("mag","magIn3", zl, defRadL, defDens, barSect);
  lrIn3->setXYOffsets(8.28, -15.72);
  auto& pin30 = lrIn3->addPolygon(nvBlue, dataBlue[0],dataBlue[1], x2x0, xrho);
  pin30.poly[0].x += spacing + barSect;
  pin30.poly[1].x += spacing + barSect;  
  auto& pin31 = lrIn3->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin31.poly[0].x = 61.3400;
  pin31.poly[1].x = 67.5800;  
  lrIn3->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn4 = det->AddPolyLayer("mag","magIn4", zl, defRadL, defDens, barSect);
  lrIn4->setXYOffsets(8.28, -15.72);
  auto& pin40 = lrIn4->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  auto& pin41 = lrIn4->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin41.poly[0].x = 54.201;
  pin41.poly[1].x = 61.3400;
  lrIn4->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn5 = det->AddPolyLayer("mag","magIn5", zl, defRadL, defDens, barSect);
  lrIn5->setXYOffsets(8.28, -15.72);
  auto& pin50 = lrIn5->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pin50.poly[0].x += spacing + barSect;
  pin50.poly[1].x += spacing + barSect;  
  auto& pin51 = lrIn5->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin51.poly[0].x = 54.201;
  pin51.poly[1].x = 61.3400;
  lrIn5->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn6 = det->AddPolyLayer("mag","magIn6", zl, defRadL, defDens, barSect);
  lrIn6->setXYOffsets(8.28, -15.72);
  auto& pin60 = lrIn6->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pin60.poly[0].x += 2*(spacing + barSect);
  pin60.poly[1].x += 2*(spacing + barSect);  
  auto& pin61 = lrIn6->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin61.poly[0].x = 54.201;
  pin61.poly[1].x = 61.3400;
  lrIn6->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn7 = det->AddPolyLayer("mag","magIn7", zl, defRadL, defDens, barSect);
  lrIn7->setXYOffsets(8.28, -15.72);
  auto& pin70 = lrIn7->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pin70.poly[0].x += 3*(spacing + barSect);
  pin70.poly[1].x += 3*(spacing + barSect);  
  auto& pin71 = lrIn7->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  pin71.poly[0].x = 54.201;
  pin71.poly[1].x = 61.3400;
  lrIn7->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn8 = det->AddPolyLayer("mag","magIn8", zl, defRadL, defDens, barSect);
  lrIn8->setXYOffsets(8.28, -15.72);
  auto& pin80 = lrIn8->addPolygon( nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  lrIn8->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn9 = det->AddPolyLayer("mag","magIn9", zl, defRadL, defDens, barSect);
  lrIn9->setXYOffsets(8.28, -15.72);
  auto& pin90 = lrIn9->addPolygon( nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  //  pin90.poly[0].x += (spacing + barSect);
  //  pin90.poly[1].x += (spacing + barSect);
  lrIn9->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn10 = det->AddPolyLayer("mag","magIn10", zl, defRadL, defDens, barSect);
  lrIn10->setXYOffsets(8.28, -15.72);
  auto& pin100 = lrIn10->addPolygon( nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  //  pin100.poly[0].x += 2*(spacing + barSect);
  //  pin100.poly[1].x += 2*(spacing + barSect);
  lrIn10->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn11 = det->AddPolyLayer("mag","magIn11", zl, defRadL, defDens, barSect);
  lrIn11->setXYOffsets(8.28, -15.72);
  auto& pin110 = lrIn11->addPolygon( nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  //  pin110.poly[0].x += 3*(spacing + barSect);
  //  pin110.poly[1].x += 3*(spacing + barSect);
  lrIn11->setNSectorsPhiStart(nSectors, phi0);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn12 = det->AddPolyLayer("mag","magIn12", zl, defRadL, defDens, barSect);
  lrIn12->setXYOffsets(8.28, -15.72);
  auto& pin120 = lrIn12->addPolygon( nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  //  pin120.poly[0].x += 4*(spacing + barSect);
  //  pin120.poly[1].x += 4*(spacing + barSect);
  pin120.poly[2].x = pin120.poly[3].x = 92.0; // fix edges
  lrIn12->setNSectorsPhiStart(nSectors, phi0);

  //============================================
  ///pin100.poly[2].x = pin90.poly[3].x = 188.4; // fix edges

  zIn += dZ;

  //  zl = zIn - barSect/2;
  //  KMCPolyLayer* lrOut0 = det->AddPolyLayer("mag","magOut0", zl, defRadL, defDens, barSect);
  //  auto& pout00 = lrOut0->addPolygon(nvCyan, dataCyan[0],dataCyan[1], x2x0, xrho);
  //  pout00.poly[0].x = 52;
  //  pout00.poly[2].x = pout00.poly[3].x = 188.4; // fix edges
  //  lrOut0->setNSectorsPhiStart(nSectors, phi0);

  //-----------------
  zl = zIn - barSect/2;
  // zl -= spacing + barSect;
  KMCPolyLayer* lrOut1 = det->AddPolyLayer("mag","magOut1", zl, defRadL, defDens, barSect);
  lrOut1->setXYOffsets(8.28, -15.72);
  auto& pout10 = lrOut1->addPolygon(nvMagenta, dataMagenta[0],dataMagenta[1], x2x0, xrho);
  auto& pout11 = lrOut1->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout11.poly[0].x = 67.570;
  pout11.poly[1].x = 78.0000;  
  pout11.poly[2].x = pout11.poly[3].x = 340.; // fix edges
  lrOut1->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut2 = det->AddPolyLayer("mag","magOut2", zl, defRadL, defDens, barSect); 
  lrOut2->setXYOffsets(8.28, -15.72);
  auto& pout20 = lrOut2->addPolygon(nvBlue, dataBlue[0],dataBlue[1], x2x0, xrho);
  auto& pout21 = lrOut2->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout21.poly[0].x = 61.3400;
  pout21.poly[1].x = 67.5800;  
  lrOut2->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut3 = det->AddPolyLayer("mag","magOut3", zl, defRadL, defDens, barSect);
  lrOut3->setXYOffsets(8.28, -15.72);
  auto& pout30 = lrOut3->addPolygon(nvBlue, dataBlue[0],dataBlue[1], x2x0, xrho);
  pout30.poly[0].x += spacing + barSect;
  pout30.poly[1].x += spacing + barSect;  
  auto& pout31 = lrOut3->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout31.poly[0].x = 61.3400;
  pout31.poly[1].x = 67.5800;  
  lrOut3->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut4 = det->AddPolyLayer("mag","magOut4", zl, defRadL, defDens, barSect);
  lrOut4->setXYOffsets(8.28, -15.72);
  auto& pout40 = lrOut4->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  auto& pout41 = lrOut4->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout41.poly[0].x = 54.201;
  pout41.poly[1].x = 61.3400;
  lrOut4->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut5 = det->AddPolyLayer("mag","magOut5", zl, defRadL, defDens, barSect);
  lrOut5->setXYOffsets(8.28, -15.72);
  auto& pout50 = lrOut5->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pout50.poly[0].x += spacing + barSect;
  pout50.poly[1].x += spacing + barSect;  
  auto& pout51 = lrOut5->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout51.poly[0].x = 54.201;
  pout51.poly[1].x = 61.3400;
  lrOut5->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut6 = det->AddPolyLayer("mag","magOut6", zl, defRadL, defDens, barSect);
  lrOut6->setXYOffsets(8.28, -15.72);
  auto& pout60 = lrOut6->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pout60.poly[0].x += 2*(spacing + barSect);
  pout60.poly[1].x += 2*(spacing + barSect);  
  auto& pout61 = lrOut6->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout61.poly[0].x = 54.201;
  pout61.poly[1].x = 61.3400;
  lrOut6->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut7 = det->AddPolyLayer("mag","magOut7", zl, defRadL, defDens, barSect);
  lrOut7->setXYOffsets(8.28, -15.72);
  auto& pout70 = lrOut7->addPolygon(nvGreen, dataGreen[0],dataGreen[1], x2x0, xrho);
  pout70.poly[0].x += 3*(spacing + barSect);
  pout70.poly[1].x += 3*(spacing + barSect);  
  auto& pout71 = lrOut7->addPolygon(nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  pout71.poly[0].x = 54.201;
  pout71.poly[1].x = 61.3400;
  lrOut7->setNSectorsPhiStart(nSectors, phi0);

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut8 = det->AddPolyLayer("mag","magOut8", zl, defRadL, defDens, barSect);
  lrOut8->setXYOffsets(8.28, -15.72);
  auto& pout80 = lrOut8->addPolygon( nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  lrOut8->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut9 = det->AddPolyLayer("mag","magOut9", zl, defRadL, defDens, barSect);
  lrOut9->setXYOffsets(8.28, -15.72);
  auto& pout90 = lrOut9->addPolygon( nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  //  pout90.poly[0].x += (spacing + barSect);
  //  pout90.poly[1].x += (spacing + barSect);
  lrOut9->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut10 = det->AddPolyLayer("mag","magOut10", zl, defRadL, defDens, barSect);
  lrOut10->setXYOffsets(8.28, -15.72);
  auto& pout100 = lrOut10->addPolygon( nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  //  pout100.poly[0].x += 2*(spacing + barSect);
  //  pout100.poly[1].x += 2*(spacing + barSect);
  lrOut10->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut11 = det->AddPolyLayer("mag","magOut11", zl, defRadL, defDens, barSect);
  lrOut11->setXYOffsets(8.28, -15.72);
  auto& pout110 = lrOut11->addPolygon( nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  //  pout110.poly[0].x += 3*(spacing + barSect);
  //  pout110.poly[1].x += 3*(spacing + barSect);
  lrOut11->setNSectorsPhiStart(nSectors, phi0);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut12 = det->AddPolyLayer("mag","magOut12", zl, defRadL, defDens, barSect);
  lrOut12->setXYOffsets(8.28, -15.72);
  auto& pout120 = lrOut12->addPolygon( nvCyan, dataCyanO[0],dataCyanO[1], x2x0, xrho);
  //  pout120.poly[0].x += 4*(spacing + barSect);
  //  pout120.poly[1].x += 4*(spacing + barSect);
  pout120.poly[2].x = pout120.poly[3].x = 92.0; // fix edges
  lrOut12->setNSectorsPhiStart(nSectors, phi0); 
}
