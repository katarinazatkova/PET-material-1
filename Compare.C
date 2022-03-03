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

void Compare(){

  // INPUT
  // Declaration of the leaf types  
  Double_t	x11;
  Double_t	x12;
  Double_t	y11;
  Double_t	y12;
  Double_t	z11;
  Double_t	z12;
  Double_t  E_det1;
  
  Double_t	x21;
  Double_t	x22;
  Double_t	y21;
  Double_t	y22;
  Double_t	z21;
  Double_t	z22;
  Double_t  E_det2;
    
  // List of branches
  TBranch *E_abs;
  
  // Input files
  TFile *inputFileNoGauss   = new TFile("BasicNoGaussKatarina.root", "read");
  TFile *inputFileWithGauss = new TFile("BasicWithGaussKatarina.root", "read");
  
  // Output variables
  Double_t	d1;
  Double_t	d2;  
  Double_t displacement;
  Double_t def_angle;
  
  TTree *inputTreeNoGauss;
  TTree *inputTreeWithGauss;
  
  inputFileNoGauss->GetObject("BasicWithGauss", inputTreeNoGauss);
  inputFileWithGauss->GetObject("BasicWithGauss", inputTreeWithGauss);
  
  inputTreeNoGauss->SetBranchAddress("x1", &x11);
  inputTreeNoGauss->SetBranchAddress("x2", &x12);
  inputTreeNoGauss->SetBranchAddress("y1", &y11);
  inputTreeNoGauss->SetBranchAddress("y2", &y12);
  inputTreeNoGauss->SetBranchAddress("z1", &z11);
  inputTreeNoGauss->SetBranchAddress("z2", &z12);

  inputTreeNoGauss->SetBranchAddress("E_detector", &E_det1);
  
  inputTreeWithGauss->SetBranchAddress("x1", &x21);
  inputTreeWithGauss->SetBranchAddress("x2", &x22);
  inputTreeWithGauss->SetBranchAddress("y1", &y21);
  inputTreeWithGauss->SetBranchAddress("y2", &y22); 
  inputTreeWithGauss->SetBranchAddress("z1", &z21);
  inputTreeWithGauss->SetBranchAddress("z2", &z22);

  inputTreeWithGauss->SetBranchAddress("E_detector", &E_det2);
   
  Long64_t nentries = inputTreeWithGauss->GetEntriesFast();
  
  // We want to consider only the good events for the analysis.
  // Typically, such events occurs if more than 91.8% of the original 
  // energy of 1.022 MeV is deposited in the detector.
  // energy resolution for LYSO is 8.2% 
  // good event = edep > 91.8% of 1.022MeV
  Double_t EnergyRes = 1.022*0.082;
  Double_t Threshold = 1.022 - EnergyRes;
  


  // Plots to draw:
  // - Histogram of the observed (radionuclide) displacement
  // - Number of events against the deflection angle (check of the gaussian distribution)
  // - Histogram of the hits in the X-Y plane
  // - 3D Histogram of the hits in the X-Y plane
  
  TH1D * d = new TH1D("d",
		      "Differences between hits with and without the non-collinearity(mm)",
		      100, 0, 1000);

  TH1D * gausscheck = new TH1D("gausscheck",
		      "Number of Events against deflection angle(radians)",
		      100, 0, 3);

  TH2F * h2_YZ = new TH2F("h2_YZ",
          "2D Histogram of the hits in the Y-Z plane",
          100, -40, 40, 100, -40, 40);

  TH2F * h3_YZ = new TH2F("h3_YZ",
          "3D Histogram of the hits in the Y-Z plane",
          100, -20, 20, 100, -20, 20);


  for(Long64_t entry = 0; entry < nentries; entry++){
  
    inputTreeWithGauss->GetEntry(entry);
    inputTreeNoGauss->GetEntry(entry);

    //Considering only the good events
    if ((E_det1 > Threshold) && (E_det2 > Threshold)){
    
      d1 = Sqrt( pow((x11-x21),2) + pow((y11-y21),2) + pow((z11-z21),2) ); 
      d2 = Sqrt( pow((x12-x22),2) + pow((y12-y22),2) + pow((z12-z22),2) );
      displacement = d2 - d1;

      def_angle = asin(Sqrt(pow((y22),2)) + Sqrt(pow((z22),2))/100);


      //gamma 1
      cout << endl;
      cout << " x11 = " << x11 << endl;
      cout << " x21 = " << x21 << endl;

     cout << endl;
      cout << " y11 = " << y11 << endl;
      cout << " y21 = " << y21 << endl;
    
      cout << endl;
      cout << " z12 = " << z12 << endl;
      cout << " z22 = " << z22 << endl;


      //gamma 2
      cout << endl;
      cout << " x12 = " << x12 << endl;
      cout << " x22 = " << x22 << endl;

      cout << endl;
      cout << " y12 = " << y12 << endl;
      cout << " y22 = " << y22 << endl;
    
      cout << endl;
      cout << " z12 = " << z12 << endl;
      cout << " z22 = " << z22 << endl;


      cout << "Distance d1 is: " << d1 << endl;
      cout << "Distance d2 is: " << d2 << endl;
      
      
      d->Fill(displacement);

      gausscheck->Fill(def_angle);

      h2_YZ->Fill(y22, z22, E_det2);

      h3_YZ->Fill(y22, z22, x22);

    }
  }

  TCanvas * c1 = new TCanvas();
  d->SetLineColor(kBlue);
  d->Draw();
  c1->SaveAs("d.pdf");

  TCanvas * c2 = new TCanvas();
  gausscheck->SetLineColor(kBlue);
  gausscheck->Draw();
  c2->SaveAs("gausscheck.pdf");

  TCanvas * c3 = new TCanvas();
  h2_YZ->SetContour(1000);
  h2_YZ->Draw("colz");   
  c3->SaveAs("h2_XY.pdf");

  TCanvas * c4 = new TCanvas();
  h3_YZ->Draw("surf2d");   
  c4->SaveAs("h3_YZ.pdf");

  //Creating a new file that will contain all the plots
  TFile * plots_file = new TFile("plots.root", "RECREATE");
  plots_file->cd();
  
  d->Write();
  gausscheck->Write();
  h2_YZ->Write();
  h3_YZ->Write();

  plots_file->Close();

}
