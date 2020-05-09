#include "KMCMagnetBuilder.h"
#include "KMCDetectorFwd.h"

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


KMCMagnetBuilder::KMCMagnetBuilder(KMCDetectorFwd* det, float z)
{
  const float barSect = 5.0;
  const float spacing = 0.24;
  float defRadL = 0;
  float defDens = 0;

  float zl = z+barSect/2;
  KMCPolyLayer* lrIn0 = det->AddPolyLayer("mag","magIn0", zl, defRadL, defDens, barSect);
  auto pin00 = lrIn0->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin00.poly[0].x = 52;
  pin00.poly[2].x = pin00.poly[3].x = 188.4; // fix edges

  //-----------------
  zl += spacing + barSect;
  KMCPolyLayer* lrIn1 = det->AddPolyLayer("mag","magIn1", zl, defRadL, defDens, barSect);
  auto pin10 = lrIn1->addPolygon(nvMagenta, dataMagenta[0],dataMagenta[1] );
  auto pin11 = lrIn1->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin11.poly[0].x = 67.570;
  pin11.poly[1].x = 78.0000;  
  pin11.poly[2].x = pin11.poly[3].x = 188.4; // fix edges

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn2 = det->AddPolyLayer("mag","magIn2", zl, defRadL, defDens, barSect);
  auto pin20 = lrIn2->addPolygon(nvBlue, dataBlue[0],dataBlue[1] );
  auto pin21 = lrIn2->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin21.poly[0].x = 61.3400;
  pin21.poly[1].x = 67.5800;  

  zl += spacing + barSect;
  KMCPolyLayer* lrIn3 = det->AddPolyLayer("mag","magIn3", zl, defRadL, defDens, barSect);
  auto pin30 = lrIn3->addPolygon(nvBlue, dataBlue[0],dataBlue[1] );
  pin30.poly[0].x += spacing + barSect;
  pin30.poly[1].x += spacing + barSect;  
  auto pin31 = lrIn3->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin31.poly[0].x = 61.3400;
  pin31.poly[1].x = 67.5800;  

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn4 = det->AddPolyLayer("mag","magIn4", zl, defRadL, defDens, barSect);
  auto pin40 = lrIn4->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  auto pin41 = lrIn4->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin41.poly[0].x = 54.201;
  pin41.poly[1].x = 61.3400;

  zl += spacing + barSect;
  KMCPolyLayer* lrIn5 = det->AddPolyLayer("mag","magIn5", zl, defRadL, defDens, barSect);
  auto pin50 = lrIn5->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pin50.poly[0].x += spacing + barSect;
  pin50.poly[1].x += spacing + barSect;  
  auto pin51 = lrIn5->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin51.poly[0].x = 54.201;
  pin51.poly[1].x = 61.3400;

  zl += spacing + barSect;
  KMCPolyLayer* lrIn6 = det->AddPolyLayer("mag","magIn6", zl, defRadL, defDens, barSect);
  auto pin60 = lrIn6->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pin60.poly[0].x += 2*(spacing + barSect);
  pin60.poly[1].x += 2*(spacing + barSect);  
  auto pin61 = lrIn6->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin61.poly[0].x = 54.201;
  pin61.poly[1].x = 61.3400;

  zl += spacing + barSect;
  KMCPolyLayer* lrIn7 = det->AddPolyLayer("mag","magIn7", zl, defRadL, defDens, barSect);
  auto pin70 = lrIn7->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pin70.poly[0].x += 3*(spacing + barSect);
  pin70.poly[1].x += 3*(spacing + barSect);  
  auto pin71 = lrIn7->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pin71.poly[0].x = 54.201;
  pin71.poly[1].x = 61.3400;

  //-----------------  
  zl += spacing + barSect;
  KMCPolyLayer* lrIn8 = det->AddPolyLayer("mag","magIn8", zl, defRadL, defDens, barSect);
  auto pin80 = lrIn8->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );

  zl += spacing + barSect;
  KMCPolyLayer* lrIn9 = det->AddPolyLayer("mag","magIn9", zl, defRadL, defDens, barSect);
  auto pin90 = lrIn9->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pin90.poly[0].x += (spacing + barSect);
  pin90.poly[1].x += (spacing + barSect);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn10 = det->AddPolyLayer("mag","magIn10", zl, defRadL, defDens, barSect);
  auto pin100 = lrIn10->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pin100.poly[0].x += 2*(spacing + barSect);
  pin100.poly[1].x += 2*(spacing + barSect);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn11 = det->AddPolyLayer("mag","magIn11", zl, defRadL, defDens, barSect);
  auto pin110 = lrIn11->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pin110.poly[0].x += 3*(spacing + barSect);
  pin110.poly[1].x += 3*(spacing + barSect);

  zl += spacing + barSect;
  KMCPolyLayer* lrIn12 = det->AddPolyLayer("mag","magIn12", zl, defRadL, defDens, barSect);
  auto pin120 = lrIn12->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pin120.poly[0].x += 4*(spacing + barSect);
  pin120.poly[1].x += 4*(spacing + barSect);
  pin100.poly[2].x = pin100.poly[3].x = 92.0; // fix edges
  
  //============================================
  ///pin100.poly[2].x = pin90.poly[3].x = 188.4; // fix edges

  zl = z + 430.52;
  
  zl = z-barSect/2;
  KMCPolyLayer* lrOut0 = det->AddPolyLayer("mag","magOut0", zl, defRadL, defDens, barSect);
  auto pout00 = lrOut0->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout00.poly[0].x = 52;
  pout00.poly[2].x = pout00.poly[3].x = 188.4; // fix edges

  //-----------------
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut1 = det->AddPolyLayer("mag","magOut1", zl, defRadL, defDens, barSect);
  auto pout10 = lrOut1->addPolygon(nvMagenta, dataMagenta[0],dataMagenta[1] );
  auto pout11 = lrOut1->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout11.poly[0].x = 67.570;
  pout11.poly[1].x = 78.0000;  
  pout11.poly[2].x = pout11.poly[3].x = 188.4; // fix edges

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut2 = det->AddPolyLayer("mag","magOut2", zl, defRadL, defDens, barSect);
  auto pout20 = lrOut2->addPolygon(nvBlue, dataBlue[0],dataBlue[1] );
  auto pout21 = lrOut2->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout21.poly[0].x = 61.3400;
  pout21.poly[1].x = 67.5800;  

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut3 = det->AddPolyLayer("mag","magOut3", zl, defRadL, defDens, barSect);
  auto pout30 = lrOut3->addPolygon(nvBlue, dataBlue[0],dataBlue[1] );
  pout30.poly[0].x += spacing + barSect;
  pout30.poly[1].x += spacing + barSect;  
  auto pout31 = lrOut3->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout31.poly[0].x = 61.3400;
  pout31.poly[1].x = 67.5800;  

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut4 = det->AddPolyLayer("mag","magOut4", zl, defRadL, defDens, barSect);
  auto pout40 = lrOut4->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  auto pout41 = lrOut4->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout41.poly[0].x = 54.201;
  pout41.poly[1].x = 61.3400;

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut5 = det->AddPolyLayer("mag","magOut5", zl, defRadL, defDens, barSect);
  auto pout50 = lrOut5->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pout50.poly[0].x += spacing + barSect;
  pout50.poly[1].x += spacing + barSect;  
  auto pout51 = lrOut5->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout51.poly[0].x = 54.201;
  pout51.poly[1].x = 61.3400;

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut6 = det->AddPolyLayer("mag","magOut6", zl, defRadL, defDens, barSect);
  auto pout60 = lrOut6->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pout60.poly[0].x += 2*(spacing + barSect);
  pout60.poly[1].x += 2*(spacing + barSect);  
  auto pout61 = lrOut6->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout61.poly[0].x = 54.201;
  pout61.poly[1].x = 61.3400;

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut7 = det->AddPolyLayer("mag","magOut7", zl, defRadL, defDens, barSect);
  auto pout70 = lrOut7->addPolygon(nvGreen, dataGreen[0],dataGreen[1] );
  pout70.poly[0].x += 3*(spacing + barSect);
  pout70.poly[1].x += 3*(spacing + barSect);  
  auto pout71 = lrOut7->addPolygon(nvCyan, dataCyan[0],dataCyan[1] );
  pout71.poly[0].x = 54.201;
  pout71.poly[1].x = 61.3400;

  //-----------------  
  zl -= spacing + barSect;
  KMCPolyLayer* lrOut8 = det->AddPolyLayer("mag","magOut8", zl, defRadL, defDens, barSect);
  auto pout80 = lrOut8->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut9 = det->AddPolyLayer("mag","magOut9", zl, defRadL, defDens, barSect);
  auto pout90 = lrOut9->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pout90.poly[0].x += (spacing + barSect);
  pout90.poly[1].x += (spacing + barSect);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut10 = det->AddPolyLayer("mag","magOut10", zl, defRadL, defDens, barSect);
  auto pout100 = lrOut10->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pout100.poly[0].x += 2*(spacing + barSect);
  pout100.poly[1].x += 2*(spacing + barSect);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut11 = det->AddPolyLayer("mag","magOut11", zl, defRadL, defDens, barSect);
  auto pout110 = lrOut11->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pout110.poly[0].x += 3*(spacing + barSect);
  pout110.poly[1].x += 3*(spacing + barSect);

  zl -= spacing + barSect;
  KMCPolyLayer* lrOut12 = det->AddPolyLayer("mag","magOut12", zl, defRadL, defDens, barSect);
  auto pout120 = lrOut12->addPolygon( nvCyan, dataCyan[0],dataCyan[1] );
  pout120.poly[0].x += 4*(spacing + barSect);
  pout120.poly[1].x += 4*(spacing + barSect);
  pout100.poly[2].x = pout100.poly[3].x = 92.0; // fix edges
 
  
}
