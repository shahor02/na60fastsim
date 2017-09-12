#include "KMCUtils.h"
#include <TMath.h>

//==========================================================================

ClassImp(NaMaterial)
ClassImp(MagField)

NaMaterial::NaMaterial()
{ 
// Defaul Constructor
}
//----------------------------------------------------------------
//
NaMaterial::NaMaterial(const char* name, const char* title, 
		       Float_t a, Float_t z, Float_t dens, Float_t radl, 
		       Float_t absl, Float_t *elbuff) 
  : TMaterial(name,title,a,z,dens,radl,absl) 
{
  // Constructor
  if (elbuff) // ELoss Parameters are provided
    for (int i=0;i<kNELossPar;i++) fELossPar[i] = elbuff[i];
  else
    for (int i=0;i<kNELossPar;i++) fELossPar[i] = 0.;
  //
}
//----------------------------------------------------------------
//
NaMaterial::~NaMaterial() {}

//
//----------------------------------------------------------------
//
void NaMaterial::Print(Option_t* option) const
{
//Print material information
  Float_t a=0,z=0,rho=0,rl=0,il=0;
  //  
  a = GetA(); z = GetZ(); rho = GetDensity(); 
  rl = GetRadLength(); il = GetInterLength();
  //
  if (option[0] != 'h' && option[0] != 'H' ) //header is not suppressed
    printf("%-20s %7s %7s %7s %10s %10s %10s %10s %10s\n",
	   "Material","   A   ","   Z   ","Density","  Rad.L  ",
	   "  Inter.L "," ELossC0 "," ELossC1 "," ELossC2 ");
  printf("%-20s %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 GetName(),a,z,rho,rl,il,fELossPar[0],
	 fELossPar[1],fELossPar[2]);
  //
}
//----------------------------------------------------------------
//
/*inline*/ Float_t NaMaterial::GetELoss(Float_t p) const
{
  // return energy lost (per cm) at momentum p
  return fELossPar[0]+fELossPar[1]*p+fELossPar[2]*TMath::Log(p);
}
//----------------------------------------------------------------
//
void NaMaterial::Dump() const 
{
  Print();
}
//----------------------------------------------------------------
//

//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================

ClassImp(NaMixture)

NaMixture::NaMixture()
{ 
// Defaul Constructor
  fNMix = 0;
  fAMix = 0;
  fZMix = 0;
  fWMix = 0;
  //
}
//----------------------------------------------------------------
//
NaMixture::NaMixture(const char* name, const char* title, 
		     Float_t a, Float_t z, Float_t dens,
		     Float_t radl,Float_t absl, Float_t *elbuff) 
  : NaMaterial(name,title,a,z,dens,radl,absl,elbuff) 
{
  // Constructor
  fNMix = 0;
  fAMix = fZMix = fWMix = 0;
}
//----------------------------------------------------------------
//
NaMixture::~NaMixture() 
{
  if (fAMix) delete[] fAMix;
  if (fZMix) delete[] fZMix;
  if (fWMix) delete[] fWMix;
}
//

void NaMixture::SetComponents(Int_t nmixt, Float_t* a,Float_t* z,Float_t* w)
{
  // Sets the components of the mixture (Geant3 conventions preserved)
  if (fAMix || fZMix || fWMix) { 
    Error("SetComponents","Components arrays were already initialized %s",GetName());
    return;
  }
  fNMix = TMath::Abs(nmixt);
  if ( fNMix<1 ) {
    Error("SetComponents","Number of components is 0 for %s",GetName());
    return;
  }
  fAMix = new Float_t[fNMix];
  fZMix = new Float_t[fNMix];
  fWMix = new Float_t[fNMix];
  Float_t amol = 0.;
  for (int i=0;i<fNMix;i++) {
    fAMix[i] = a[i];
    fZMix[i] = z[i];
    fWMix[i] = w[i];
    amol += fWMix[i]*fAMix[i];
  }
  if (amol<=0.) {Error("SetComponents","total weigth for %s is <=0",GetName()); return;}
  if (nmixt<0 ) for (int i=0;i<fNMix;i++) fWMix[i] *= fAMix[i]/amol; //use 'proportion by weight'
  //
}

//----------------------------------------------------------------
//
void NaMixture::Print(Option_t* option) const
{
//Print mixture information
  Float_t a=0,z=0,rho=0,rl=0,il=0;
  //  
  a = GetA(); z = GetZ(); rho = GetDensity(); 
  rl = GetRadLength(); il = GetInterLength();
  //
  if (option[0] != 'h' && option[0] != 'H' ) //header is not suppressed
    printf("%-20s %7s %7s %7s %10s %10s %10s %10s %10s\n",
	   "Material","   A   ","   Z   ","Density","  Rad.L  ",
	   "  Inter.L "," ELossC0 "," ELossC1 "," ELossC2 ");
  printf("%-20s %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %6s %6s %6s\n",
	 GetName(),a,z,rho,rl,il,fELossPar[0],
	 fELossPar[1],fELossPar[2],"  A  ","  Z  ","  W  ");
  for (int i=0;i<TMath::Abs(fNMix);i++) 
    printf("%107s %6.2f %6.2f %6.4f\n"," ",fAMix[i],fZMix[i],fWMix[i]);
  //
}
//----------------------------------------------------------------
//


//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================


NaCardsInput::NaCardsInput() 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
  fLastComment = "";
}
//--------------------------------------------------
//
NaCardsInput::NaCardsInput(const char* fname) 
{
  // Initialize only
  *fbuffer = '\0';
  fPoint = 0;
  fNArgs = 0;
  fArgs = 0;
  fKey = 0;
  fModifier = 0;
  fCardsFile = 0;
  fCardsFile = OpenFile(fname);
  fLastComment = "";
}
//--------------------------------------------------
//
NaCardsInput::~NaCardsInput() 
{
  if (fCardsFile) fclose(fCardsFile);
  ClearArgs();
}
//--------------------------------------------------
//
FILE* NaCardsInput::OpenFile(const char *fname) 
{
  TString flstr = fname;
  gSystem->ExpandPathName(flstr);
  if (fCardsFile) fclose(fCardsFile);
  if ( !(fCardsFile = fopen(flstr.Data(),"r")) ) { 
    printf("Error in <NaCardsInput::OpenFile>: Did not find file %s\n",flstr.Data());
    return NULL;
  }  
  ClearArgs();
  return fCardsFile;
}
//--------------------------------------------------
//
void NaCardsInput::ClearArgs() 
{
  if (fArgs) {
    for (int i=0;i<fNArgs;i++) delete[] fArgs[i];
    delete[] fArgs;
    fArgs = 0;
  }
  if (fKey) {
    delete[] fKey;
    fKey = 0;
  }
  if (fModifier) {
    delete[] fModifier;
    fModifier = 0;
  }
  fNArgs = 0;
  fLastComment = "";
}
//--------------------------------------------------
//
int NaCardsInput::FindEntry(const char* ckey, const char* cmod,const char* form, int rew, int errmsg) 
{
  // Find and fill arguments list from the line starting with KEY:MODIFIER ...
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'b' - argument should be string of 0 and 1, interpreted as bit mask
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If rew != 0, the search is done from the beginning of the file
  // If errmsg is not 0, then the warning is printed in case KEY and MODIFIER were matched
  // but format - not
  // Case insensitive
  // If not found, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<NaCardsInput::FindEntry>  KEY is empty\n");
    return -1;
  }
  if (rew) rewind(fCardsFile);
  //
  while ( (narg=NextEntry()) > -1 ) {
    if ( CompareKey(ckey) ) continue; // KEY does not match
    if ( cmod && (cmod[0]=='\0' || cmod[0]==' ' || cmod[0]=='\t') 
	 && fModifier[0]!='\0' ) continue; // Empty modifier was requested
    if ( cmod && CompareModifier(cmod) ) continue; // MODIFIER does not match
    //
    if ( CompareArgList(form) ) {  // ArgList form does not match
      if (cmod && errmsg) {
	printf("<NaCardsInput::FindEntry> arguments of %s\n %s to \"%s\"\n",
	       fPoint,"do not match to format",form);
	//return -10;
      }
      continue;
    }
    return narg;
  }
  return -1;
}
//--------------------------------------------------
//
int NaCardsInput::NextEntry(const char* ckey, const char* cmod,const char* form)
{
  // Reads next entry in required format
  // KEY(ckey) is obligatory, while MODIFIER(cmod) can be either 0 (any
  // modifier matches), or "" or " .." - empty modifier is requested.
  // Return number of found arguments (not counting key:modifier)
  // If form is not NULL or nonempty, then looks for the record in the following
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // Case insensitive
  // If next entry does not match, return -1;
  //
  int narg=-1;
  //
  if ( !ckey ) {
    printf("<NaCardsInput::NextEntry(...)>  KEY is empty\n");
    return -1;
  }
  //
  narg=NextEntry();
  if (narg<0) return narg;
  if ( CompareKey(ckey) ) return -1; // KEY does not match
  if (!(cmod && (cmod[0]=='\0'||cmod[0]==' '||cmod[0]=='\t')&&fModifier[0]!='\0')) // NonEmpty modifier requested
    if ( cmod && CompareModifier(cmod) ) return -1; // MODIFIER does not match
    //
  if ( CompareArgList(form) ) {  // ArgList form does not match
    printf("<NaCardsInput::NextEntry> arguments of %s\n %s to \"%s\"\n",
	   fPoint,"do not match to format",form);
    return -1;
  }
  return narg;
  //
}
//--------------------------------------------------
//
int NaCardsInput::NextEntry() 
{
  // get next entry, skipping commented lines
  //
  int len;
  char *tmp1,*tmp2;
  char buftmp[fgkMaxLen]; 
  //
  if (!fCardsFile) {
    printf("<NaCardsInput::NextEntry> - No file was opened\n");
    return -1;
  }
  ClearArgs(); 
  while ( 1 ) {
    fLastPos = ftell(fCardsFile);
    if ( !(fPoint=fgets(fbuffer,fgkMaxLen,fCardsFile)) ) break; // end of file reached
    while(*fPoint == ' ' || *fPoint == '\t') fPoint++; // skip spaces in the beginning
    if (*fPoint == fgkComment) {fLastComment+=fbuffer ;continue;} // check if the line is commented
    //
    // Check if there is a line continuation flag
    tmp2 = fPoint;
    do {
      tmp1 = &fPoint[strlen(fPoint)-1];
      if (*tmp1=='\n') {*tmp1--='\0';}
      while(*tmp1 == ' ' || *tmp1 == '\t') *tmp1--='\0'; // skip spaces in the end
      if (*tmp1 == fgkContinuation) {
	if ( !(tmp2=fgets(tmp1,fgkMaxLen-(tmp1-fPoint),fCardsFile)) ) break; // end of file reached
      }
      else
	break;
    } while(1);
    if (!tmp2) break; // end of file reached
    //
    // Check if there is KEY:[MODIFIER]
    tmp1 = fPoint;
    if ( sscanf(tmp1,"%s",buftmp)<=0 ) continue; // Empty line
    if ( (tmp2=strchr(buftmp,fgkDelimiter)) ) { // there is a key/delimiter
      len = (tmp2-buftmp)/sizeof(char);
      fKey = new char[len+1];
      for (int i=0;i<len;i++) fKey[i] = tolower(buftmp[i]);
      fKey[len]='\0';
      //
      len = strlen(++tmp2); // skip delimiter, get modifier length
      fModifier = new char[len+1];
      for (int i=0;i<len;i++) fModifier[i] = tolower(tmp2[i]);
      fModifier[len] = '\0';
      // skip key:mod
      tmp1 += strlen(buftmp);
      while(*tmp1 == ' ' || *tmp1 == '\t') tmp1++; // skip spaces
    }
    else {
      continue; // no delimiter - not valid record
    }
    // First, throw away everything after comment (if any)
    if ( (tmp2=strchr(tmp1,fgkComment)) ) *tmp2 = '\0';
    // now, scan the string to get the number of arguments
    tmp2 = tmp1;
    while( sscanf(tmp2,"%s",buftmp)!= -1 ) {
      while( *tmp2 == ' ' || *tmp2 == '\t') tmp2++;
      tmp2 += strlen(buftmp); 
      fNArgs++;
    }
    // Fill the arguments list
    fArgs = new char*[fNArgs];
    for (int i=0;i<fNArgs;i++) {
      sscanf(tmp1,"%s",buftmp);
      while( *tmp1 == ' ' || *tmp1 == '\t') tmp1++;
      int lgt = strlen(buftmp);
      tmp1 += lgt;
      fArgs[i] = new char[lgt+1];
      strcpy(fArgs[i],buftmp);
    }
    return fNArgs;
  }
  return -1;
}
//--------------------------------------------------
//
char* NaCardsInput::GetArg(const int iarg, const char* option, int *err) 
{
  // Return argument iarg in character string format
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  //
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      carg[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      carg[i] = toupper(carg[i]);
  //
  return carg;
}
//--------------------------------------------------
//
char* NaCardsInput::GetArg(char* &dest, const int iarg, const char* option, int *err) 
{
  // Creates new string at dest and fill it with argument iarg 
  // in character string format.
  // Optionally converting it to upper or lower case
  // if iarg is wrong, return 0 and set err to 1
  // No check of dest == 0 is done!
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  char *carg = fArgs[iarg];
  int lgt = strlen(carg);
  dest = new char[lgt+1];
  if (!strncasecmp(option,"l",1))  // convert to lower case
    for (int i=0;i<lgt;i++) 
      dest[i] = tolower(carg[i]);
  else if (!strncasecmp(option,"u",1)) // convert to upper case
    for (int i=0;i<lgt;i++) 
      dest[i] = toupper(carg[i]);
  else
    for (int i=0;i<lgt;i++) dest[i] = carg[i];
  //
  dest[lgt] = '\0';
  return dest;
}
//--------------------------------------------------
//
float NaCardsInput::GetArgF(const int iarg, int *err) 
{
  // Get requsted argument in float format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0.; }
  //
  char *errc = 0;
  float vald = float(strtod(fArgs[iarg],&errc));
  if (*errc) { *err = 1; return 0.; }
  return vald;
  //
}
//--------------------------------------------------
//
int NaCardsInput::GetArgD(const int iarg, int *err) 
{
  // Get requsted argument in integer format, it it is not float, set Error to 1
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  //
  char *errc = 0;
  int vald = int(strtol(fArgs[iarg],&errc,10));
  if (*errc) { *err = 1; return 0; }
  return vald;
  //
}
//--------------------------------------------------
//
unsigned int NaCardsInput::GetArgB(const int iarg, int *err) 
{
  // Get requsted argument in bitfield format, interpretting 0s and 1s from right side
  err = 0;
  if (iarg>=fNArgs || iarg<0) { *err = 1; return 0; }
  //
  unsigned int valu = 0;
  char *str = fArgs[iarg];
  int lgt = strlen(str);
  for (int ibt=lgt;ibt--;) {
    if (str[ibt]=='1') valu |= 1<<(lgt-1-ibt);
    else if (str[ibt]!='0') {
      printf("Bit string must contain only 0 or 1: %s\n",str);
      *err = 1; 
      return 0;
    }
  }
  return valu;
  //
}
//--------------------------------------------------
//
int NaCardsInput::CompareKey(const char* key) 
{
  if (fKey) return strcasecmp(key,fKey);
  return 0;
}
//--------------------------------------------------
//
int NaCardsInput::CompareModifier(const char* mod) 
{
  if (fModifier) return strcasecmp(mod,fModifier);
  return 1;
}
//--------------------------------------------------
//
int NaCardsInput::CompareArgList(const char *form)
{
  // If form is not NULL or nonempty, check if the record matches following 
  // format according to characters in the "form" string: 
  // 's' - argument should be string (any type of argument mathces)
  // 'a' - argument should be string with at least one non-numeric char
  // 'b' - argument should be string of 0 and 1, interpreted as bit mask
  // 'd' - argument should be integer
  // 'f' - argument should be float (no check for . is done so the integer too matches)
  // '?' - following arguments are optional, but their number should be not less than in 'form' 
  // '*' - following arguments are optional, and their number can be less than in 'form'
  // '|' - no more arguments is allowed
  // If matches, return 0, otherwise -1
  //
  char cf=' ';
  const char *pf = form;
  char *err = 0;
  if ( !form ) return 0; // no for is requested, anything matches
  int iarg = 0;
  float valf=0.0;
  int vald = 0;
  int isoption = 0;
  //
  if (fNArgs == -1 ) {
    printf("<NaCardsInput::CompareArgList> There is no record in the buffer\n");
    return -1;
  }

  while ( (cf=tolower(*pf++)) ) {
    //
    if (cf=='*') { // the following arguments are optional
      isoption = 1;
      continue;
    }
    //
    if (cf=='?') { // the following argument block is optional
      if ( fNArgs > iarg+1 ) {
	isoption = 0; // number of following arguments should be respected
	continue;
      }
      else return 0; // no more arguments, but the block was optional
    }
    //
    if (cf=='|') { // no more arguments are alowed
      if ( fNArgs != iarg ) {printf("err1\n");return -1;} // number of arguments > than allowed
      else return 0; // no more arguments, OK
    }
    //
    if ( iarg >= fNArgs ) { // number of agruments is less than requested
      if ( isoption || cf=='|' ) return 0; // but this was allowed -> MATCHED
      else {printf("err2\n");return -1;} // too few arguments
    }
    //
    // Analize requsted argument type
    //
    //
    if (cf=='s') { // string argument is requested, anything matches
      iarg++;
      continue;
    }
    //
    if (cf=='d') { // Integer argument is requested
      vald = strtol(fArgs[iarg++],&err,10);
      if (*err) {printf("err3\n");return -1;} // Not integer, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='f') { // Float argument is requested
      valf = strtod(fArgs[iarg++],&err);
      if (*err) {printf("err4\n");return -1;} // Not float, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='a') { // string argument with non-numeric character is requested
      valf = strtod(fArgs[iarg++],&err);
      if (!(*err)) {printf("err5: |%s|\n",fArgs[iarg-1]);return -1;} // looks like number, does not match
      else continue; // this argument matches
    }
    //
    if (cf=='b') { // 0,1 bitfield is requested
      char* str = fArgs[iarg++];
      int lgt = strlen(str);
      if (lgt>32) {printf("err6\n");return -1;}
      for (int i=0;i<lgt;i++) if (str[i]!='0'&&str[i]!='1') return -1; 
      continue; // this argument matches
    }    
    //
    // unknown format character
    printf("<NaCardsInput::CompareArgList> unknown format requsted \'%c\'\n",cf);
  }
  return 0; // Matched
}
//--------------------------------------------------
//
void NaCardsInput::Print() 
{
  if (fKey) printf("Key = %s\n",fKey);
  else printf("No Key\n");
  if (fModifier) printf("Modifier = %s\n",fModifier);
  else printf("No Modifiers\n");
  for (int i=0;i<fNArgs;i++) {
    printf("Arg %d: %s\n",i+1,fArgs[i]);
  }
}

const char NaCardsInput::fgkComment;   // comment identifier
const char NaCardsInput::fgkDelimiter; // delimiter between keyword and modifier
const char NaCardsInput::fgkContinuation; // delimiter between keyword and modifier
const int  NaCardsInput::fgkMaxLen;  // max. length of the entry


//---------------------------------

const double MagField::fZMin[MagField::kNReg] = {0., 700.};// cm, cm

const double MagField::fZMax[MagField::kNReg] = {40.,1020.};//cm cm
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
