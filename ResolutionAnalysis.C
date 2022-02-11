#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TMath.h"

#include <iostream> 

using namespace std;
using namespace TMath;

void ResolutionAnalysis(){

  // INPUT
  // Declaration of the leaf types  
  Double_t	x11;
  Double_t	x12;
  Double_t	y11;
  Double_t	y12;
  Double_t	z11;
  Double_t	z12;
  
  Double_t	x21;
  Double_t	x22;
  Double_t	y21;
  Double_t	y22;
  Double_t	z21;
  Double_t	z22;
    
  // Input files
  TFile *inputFileNoGauss = new TFile("BasicNoGauss0.root", "read");
  TFile *inputFileWithGauss = new TFile("BasicWithGauss0.root", "read");
  
  // Output variables
  Double_t	d1;
  Double_t	d2;  
  
  // Create a directory
  TDirectory * dirNoGauss = (TDirectory*)inputFileNoGauss->Get("BasicNoGauss0.root:/ntuple");
  TDirectory * dirWithGauss = (TDirectory*)inputFileWithGauss->Get("BasicWithGauss0.root:/ntuple");
  
  TTree *inputTreeNoGauss;
  TTree *inputTreeWithGauss;
  
  dirNoGauss->GetObject("BasicNoGauss", inputTreeNoGauss);
  dirWithGauss->GetObject("BasicWithGauss", inputTreeWithGauss);
  
  inputTreeNoGauss->SetBranchAddress("x11", &x11);
  inputTreeNoGauss->SetBranchAddress("x12", &x12);
  inputTreeNoGauss->SetBranchAddress("y11", &y11);
  inputTreeNoGauss->SetBranchAddress("y12", &y12); 
  inputTreeNoGauss->SetBranchAddress("z11", &z11);
  inputTreeNoGauss->SetBranchAddress("z12", &z12); 
  
  inputTreeWithGauss->SetBranchAddress("x21", &x21);
  inputTreeWithGauss->SetBranchAddress("x22", &x22);
  inputTreeWithGauss->SetBranchAddress("y21", &y21);
  inputTreeWithGauss->SetBranchAddress("y22", &y22); 
  inputTreeWithGauss->SetBranchAddress("z21", &z21);
  inputTreeWithGauss->SetBranchAddress("z22", &z22);
   
  Long64_t nentries = inputTreeNoGauss->GetEntriesFast();
  
  TH1D * d = new TH1D("d",
                       "Differences between hits with and without the non-collinearity(mm)",
                       100, 0, 2000);  
  
  for(Long64_t entry = 0; entry < nentries; entry++){
  
    inputTreeNoGauss->GetEntry(entry);
    inputTreeWithGauss->GetEntry(entry);
    
    d1 = Sqrt( pow((x11-x21),2) + pow((y11-y21),2) + pow((z11-z21),2) );
    d2 = Sqrt( pow((x12-x22),2) + pow((y12-y22),2) + pow((z12-z22),2) );
    
    cout << "x-coordinate of hit 1 for no gauss deviation: " << x11 << endl;
    cout << "x-coordinate of hit 2 for no gauss deviation: " << x12 << endl;
    cout << "x-coordinate of hit 1 with gauss deviation: " << x21 << endl;
    cout << "x-coordinate of hit 2 with gauss deviation: " << x22 << endl;
    
    cout << "Distance d1 is: " << d1 << endl;
    cout << "Distance d2 is: " << d2 << endl;

    d->Fill(d1);
    d->Fill(d2);

    d->SetLineColor(kBlue);
    d->Draw();
}
