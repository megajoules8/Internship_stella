//values from run_01_2
#include "libPerso.h"

//Energy Calibration for QDC 
Double_t myLin (Double_t *var, Double_t *par){
  // par[0]: offset
  // par[1]: slope

  return par[0] + par[1]*var[0];
}

void energy_calibration () {


// defining the calibration data
 const Int_t n = 3;
  Double_t QDC_Channel_peaks[n] = {12050 , 12830 , 13590 }; //taken from fitting procedure in event_reconstruction.C
  Double_t uncertainty_QDC[n] = {0.5 , 0.4 , 1.606 }; 
  Double_t Energy_Calibrated[n] = {5156.59 , 5485.56 , 5762.64 }; //source:nndc.bnl.gov [keV]
  Double_t uncertainty_energy[n] = {0.045 , 0.03 , 0.02 }; 

TGraphErrors* gr = new TGraphErrors (n, QDC_Channel_peaks, Energy_Calibrated, uncertainty_QDC, uncertainty_energy);  

// fit function
TF1 *f_fit = new TF1 ("f_fit", myLin, 9700., 14000., 2);

// display and fit
  
  TCanvas *c1 = new TCanvas ("c1", "Energy Calibration", 700, 500);
  gr->SetTitle("Energy Calibration");
  gr->GetXaxis()->SetTitle("Channel");
  gr->GetYaxis()->SetTitle("Energy [keV]"); 
  gr->Draw ("A*"); //ask what star stands for

   gr->Fit (f_fit, "R"); // not Q i.e. display of fit parameters in command line
  //  f_fit->Draw ("SAME"); // can check parameters when commenting out fit

  TCanvas *c_comp = new TCanvas ("c_comp", "", 700, 500);
  c_comp->Draw ();
  c_comp->cd ()->SetTickx ();
  c_comp->cd ()->SetTicky ();
  //  c_comp->cd ()->SetLogy ();
    
  gPad->SetRightMargin (.005);
  gPad->SetTopMargin (.005);
  gPad->SetLeftMargin (.13);
  gPad->SetBottomMargin (.13);

 
/*  TH1F *h_dummy = new TH1F ("h_dummy", "Energy Calibration", 11, 0., 11.);
  gr->GetXaxis()->SetTitle("Channel");
  gr->GetYaxis()->SetTitle("Energy [keV]"); 
  gr->Draw ("A*"); //ask what star stands for

  gr->Fit (f_fit, "R"); // not Q i.e. display of fit parameters in command line
  //  f_fit->Draw ("SAME"); // can check parameters when commenting out fit
  h_dummy->GetXaxis ()->SetLabelSize (.05);
  h_dummy->GetXaxis ()->SetLabelOffset (.01);
  h_dummy->GetXaxis ()->SetTitleSize (.065);
  h_dummy->GetXaxis ()->SetTitleOffset (0.9); // 65
  h_dummy->GetXaxis ()->CenterTitle ();
 
  h_dummy->GetYaxis ()->SetLabelSize (h_dummy->GetXaxis ()->GetLabelSize ());
  h_dummy->GetYaxis ()->SetLabelOffset (h_dummy->GetXaxis ()->GetLabelOffset ());
  h_dummy->GetYaxis ()->SetTitleSize (h_dummy->GetXaxis ()->GetTitleSize ());
  h_dummy->GetYaxis ()->SetTitleOffset (h_dummy->GetXaxis ()->GetTitleOffset ());
  h_dummy->GetYaxis ()->CenterTitle ();

 // h_dummy->GetXaxis ()->SetTitle ("d_{target} [cm]");
 // h_dummy->GetYaxis ()->SetTitle ("Acceptance [%]");
  h_dummy->SetTitle ("Energy Calibration");
  //  h_dummy->GetXaxis ()->SetRangeUser (1., 5.8);
  h_dummy->GetYaxis ()->SetRangeUser (0., 15.);
  h_dummy->SetLineColor (0);
  h_dummy->Draw ("SAME"); */

}


