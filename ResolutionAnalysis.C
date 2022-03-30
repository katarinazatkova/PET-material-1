#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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
  TFile *inputFileNoGauss   = new TFile("BasicNoGaussC.root", "read");
  TFile *inputFileWithGauss = new TFile("BasicWithGaussC.root", "read");
  
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
  // - Histogram of the hits in the Y-Z plane
  // - 3D Histogram of the hits in the Y-Z plane
  // - Spatial resolution against the detector length
  
  
  TH1D * d = new TH1D("",
		      "Differences between hits with and without the non-collinearity (mm)",
		      100, 0, 10);
  
  TH1D * gausscheck = new TH1D("",  
		      "Number of Events against deflection angle(radians)",
		      100, 0, 10);

  TH2F * h2_YZ = new TH2F("",
          "2D Histogram of the hits in the Y-Z plane",
          100, -1000, 1000, 100, -400, 400); 

  TH2F * h3_YZ = new TH2F("",
          "3D Histogram of the hits in the Y-Z plane",
          100, -1000, 1000, 100, -400, 400);
  

  // NOTE: to obtain a plot of the spatial resolution against the detector length
  // uncomment parts of the code below
  /*
  const Int_t n = 45;//97;
  Double_t length_vals[n];
  Double_t res_vals[n];

  for (int i = 0; i <= 44; ++i) { //96
  Double_t min_range = (-1940/2) + 50 + 10*i;
  Double_t max_range = - min_range;
  length_vals[i] = 2*max_range;
  */

  // Obtaining resolution for specific detector size
  //Double_t min_range = -1940/2 + 200;
  //Double_t max_range = - min_range;


  for(Long64_t entry = 0; entry < nentries; entry++){
  
    inputTreeWithGauss->GetEntry(entry);
    inputTreeNoGauss->GetEntry(entry);

    //Considering only the good events
    if ((E_det1 > Threshold) && (E_det2 > Threshold)){
      
      // Considering shorter detector lengths
      //if (((z11 > min_range) && (z11 < max_range)) && ((z12 > min_range) && (z12 < max_range)) && ((z21 > min_range) &&  (z21 < max_range)) && ((z22 > min_range) && (z22 < max_range))){
      
      // Considering only highly oblique events at max 20 cm from the edge of the PET detector
      //if (((z11 < min_range) || (z11 > max_range)) && ((z12 < min_range) || (z12 > max_range)) && ((z21 < min_range) ||  (z21 > max_range)) && ((z22 < min_range) ||  (z22 > max_range))){
    
      d1 = Sqrt( pow((x11-x21),2) + pow((y11-y21),2) + pow((z11-z21),2) ); 
      d2 = Sqrt( pow((x12-x22),2) + pow((y12-y22),2) + pow((z12-z22),2) );
      displacement = d2 - d1;
      cout << " d = " << displacement << endl;
      def_angle = asin(Sqrt(pow((y22),2) + (pow((z22),2)))/x22);

      //gamma 1
      cout << endl;
      cout << " x11 = " << x11 << endl;
      cout << " x21 = " << x21 << endl;

      cout << endl;
      cout << " y11 = " << y11 << endl;
      cout << " y21 = " << y21 << endl;
    
      cout << endl;
      cout << " z11 = " << z11 << endl;
      cout << " z21 = " << z21 << endl;


      //gamma 2
      cout << endl;
      cout << " x21 = " << x21 << endl;
      cout << " x22 = " << x22 << endl;

      cout << endl;
      cout << " y21 = " << y21 << endl;
      cout << " y22 = " << y22 << endl;
    
      cout << endl;
      cout << " z21 = " << z21 << endl;
      cout << " z22 = " << z22 << endl;


      cout << "Distance d1 is: " << d1 << endl;
      cout << "Distance d2 is: " << d2 << endl;
      
      d->Fill(displacement);
      
      gausscheck->Fill(def_angle);

      h2_YZ->Fill(z22, y22, E_det2);

      h3_YZ->Fill(z22, y22, x22);
      
    //}
  }}
  
  // Calculating spatial resolution from the histogram
  Int_t bin = d->FindLastBinAbove(d->GetMaximum()/2);
  Double_t FWHM = 2*(d->GetBinCenter(bin));
  // Alternative way
  //Double_t stdev = d->GetStdDev();
  //Double_t FWHM = 2.3548*stdev;
  cout << endl << "FWHM(Spatial Resolution) = " << FWHM << endl;

  
  /*
  res_vals[i] = FWHM; 
  }
  
  // Plot spatial resolution against length of the detector
  TCanvas *rlc  = new TCanvas();
  TGraph *rl = new TGraph(n,length_vals,res_vals);
  rl->SetTitle("Spatial resolution as a function of detector length; Detector length (mm); FWHM (mm)");
  rl->GetXaxis()->CenterTitle(true);
  rl->GetYaxis()->CenterTitle(true);
  rl->SetLineColor(kAzure+2);
  rl->SetFillStyle(3005);
  rl->SetMarkerStyle(2);
  rl->Draw("APL");
  rlc->SaveAs("res_vs_length.pdf");
  */

  TCanvas * c1 = new TCanvas();
  d->SetTitle("Difference between the observed and the actual position of the radionuclide; Difference in position (mm); Number of events");
  d->GetXaxis()->CenterTitle(true);
  d->GetYaxis()->CenterTitle(true);
  d->SetFillColor(kAzure+1);
  d->SetLineColor(kAzure+2);
  //gStyle->SetOptStat(0);
  d->Draw();
  c1->SaveAs("d.pdf");
  
  TCanvas * c2 = new TCanvas();
  gausscheck->SetTitle("Distribution of deflection angles; Deflection angle (radians); Number of events");
  gausscheck->GetXaxis()->CenterTitle(true);
  gausscheck->GetYaxis()->CenterTitle(true);
  gausscheck->SetFillColor(kAzure+1);
  gausscheck->SetLineColor(kAzure+2);
  gausscheck->Draw();
  c2->SaveAs("gausscheck.pdf");

  TCanvas * c3 = new TCanvas();
  h2_YZ->SetTitle("Histogram of the deposited energy in the detector in the Y-Z plane; Z position (mm); Y position (mm)");
  h2_YZ->GetXaxis()->CenterTitle(true);
  h2_YZ->GetYaxis()->CenterTitle(true);
  h2_YZ->SetContour(1000);
  h2_YZ->Draw("colz");   
  c3->SaveAs("h2_YZ.pdf");
  
  TCanvas * c4 = new TCanvas();
  h3_YZ->SetTitle("Histogram of the spatial distribution of the considered events; Z position (mm); Y position (mm); X position (mm)");
  h3_YZ->GetXaxis()->CenterTitle(true);
  h3_YZ->GetYaxis()->CenterTitle(true);
  h3_YZ->GetZaxis()->CenterTitle(true);
  h3_YZ->Draw("surf2d");   
  c4->SaveAs("h3_YZ.pdf");

  /*
  //Creating a new file that will contain all the plots
  TFile * plots_file = new TFile("plots.root", "RECREATE");
  plots_file->cd();
  
  d->Write();
  gausscheck->Write();
  h2_YZ->Write();
  h3_YZ->Write();


  plots_file->Close();
  */

}
