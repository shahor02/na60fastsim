#ifndef LOCTTREESTREAM_H
#define LOCTTREESTREAM_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//  marian.ivanov@cern.ch
//
//  ------------------------------------------------------------------------------------------------
//  TLocTreeStream
//  Standard stream (cout) like input for the tree
//  Run and see TLocTreeStreamer::Test() - to see TLocTreeStreamer functionality
//  ------------------------------------------------------------------------------------------------  
//
//  -------------------------------------------------------------------------------------------------
//  TLocTreeSRedirector
//  Redirect file to  different TLocTreeStreams  
//  Run and see   TLocTreeSRedirector::Test() as an example of TLocTreeSRedirector functionality
// 

#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
class TObjArray;
class TDataType;

class TLocTreeDataElement: public TNamed {
  friend class TLocTreeStream;
 public:
  TLocTreeDataElement(Char_t type);
  TLocTreeDataElement(TDataType* type);
  TLocTreeDataElement(TClass* cl);
  void   SetPointer(void* pointer) {fPointer=pointer;} 
  Char_t GetType() const {return fType;}
 protected:

  TLocTreeDataElement(const TLocTreeDataElement & tde);
  TLocTreeDataElement & operator=(const TLocTreeDataElement & tde);

  Char_t  fType;     // type of data element
  TDataType *fDType; //data type pointer 
  TClass    *fClass; //data type pointer
  void * fPointer;  // pointer to element
  ClassDef(TLocTreeDataElement,2)
};

class TLocTreeStream: public TNamed {
  friend class TLocTreeSRedirector;
public:
  TLocTreeStream(const char *treename, TTree* externalTree=NULL);
  ~TLocTreeStream();
  void Close();
  static void Test();
  Int_t CheckIn(Char_t type, void *pointer);  
  //Int_t CheckIn(const char *type, void *pointer);
  Int_t CheckIn(TObject *o);
  void BuildTree();
  void Fill();
  Double_t GetSize(){ return fTree->GetZipBytes();}
  TLocTreeStream& Endl();
  //
  TLocTreeStream  &operator<<(Bool_t   &b){CheckIn('B',&b);return *this;}
  TLocTreeStream  &operator<<(Char_t   &c){CheckIn('B',&c);return *this;}
  TLocTreeStream  &operator<<(UChar_t  &c){CheckIn('b',&c);return *this;}
  TLocTreeStream  &operator<<(Short_t  &h){CheckIn('S',&h);return *this;}
  TLocTreeStream  &operator<<(UShort_t &h){CheckIn('s',&h);return *this;}
  TLocTreeStream  &operator<<(Int_t    &i){CheckIn('I',&i);return *this;}
  TLocTreeStream  &operator<<(UInt_t   &i){CheckIn('i',&i);return *this;}
  TLocTreeStream  &operator<<(Long_t   &l){CheckIn('L',&l);return *this;}
  TLocTreeStream  &operator<<(ULong_t  &l){CheckIn('l',&l);return *this;}
  TLocTreeStream  &operator<<(Long64_t &l){CheckIn('L',&l);return *this;}
  TLocTreeStream  &operator<<(ULong64_t &l){CheckIn('l',&l);return *this;}
  TLocTreeStream  &operator<<(Float_t   &f){CheckIn('F',&f);return *this;}
  TLocTreeStream  &operator<<(Double_t  &d){CheckIn('D',&d);return *this;}
  TLocTreeStream  &operator<<(TObject*o){CheckIn(o);return *this;} 
  TLocTreeStream  &operator<<(const Char_t *name);
  TTree * GetTree() const { return fTree;}
 protected:
  //

  TLocTreeStream(const TLocTreeStream & ts);
  TLocTreeStream & operator=(const TLocTreeStream & ts);

  TObjArray *fElements; //array of elements
  TObjArray *fBranches; //pointers to branches
  TTree *fTree;         //data storage
  Int_t fCurrentIndex;  //index of current element
  Int_t fId;            //identifier of layout
  TString fNextName;    //name for next entry
  Int_t   fNextNameCounter; //next name counter
  Int_t   fStatus;      //status of the layout
  ClassDef(TLocTreeStream,1)
};


class TLocTreeSRedirector: public TObject { 
public:
  TLocTreeSRedirector(const char *fname="", const char * option="update");
  virtual ~TLocTreeSRedirector();
  void Close();
  static void Test();
  static void Test2();
  static void UnitTestSparse(Double_t scale, Int_t testEntries);
  static void UnitTest(Int_t testEntries=5000);
  void StoreObject(TObject* object);
  TFile * GetFile() {return fDirectory->GetFile();}
  TDirectory * GetDirectory() {return fDirectory;}
  virtual   TLocTreeStream  &operator<<(Int_t id);
  virtual   TLocTreeStream  &operator<<(const char *name);
  void      SetDirectory(TDirectory *sfile); 
  void      SetFile(TFile *sfile) {SetDirectory(sfile);} 
  void SetExternalTree(const char* name, TTree* externalTree);
  static void SetDisabled(Bool_t b=kTRUE) {fgDisabled=b;}
  static Bool_t IsDisabled()        {return fgDisabled;}
    static void FixLeafNameBug(TTree* tree);
private:

  TLocTreeSRedirector(const TLocTreeSRedirector & tsr);
  TLocTreeSRedirector & operator=(const TLocTreeSRedirector & tsr);

  TDirectory* fDirectory;        //file
  Bool_t      fDirectoryOwner;   //do we own the directory?
  TObjArray *fDataLayouts;   //array of data layouts
  static Bool_t fgDisabled;  //disable - do not open any files
  ClassDef(TLocTreeSRedirector,2) 
};




#endif
