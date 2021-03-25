
//#include "../../Libs/Root/libPerso.h"
#include "libPerso.h"

Double_t my_gaus (Double_t *var, Double_t *par){
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2));
}


void fit_charge_qdcUnsort (){

  // adapted for analyzing spectra with a triple alpha source
  
  gStyle->SetOptStat (1111110);
  const Int_t run_in = 7;
  //  const Int_t thresh[4] = {150, 150, 150, 150}; // r001
  //  const Int_t thresh[4] = {1000, 1000, 1000, 1000}; // r003
  const Int_t thresh[4] = {400, 1000, 1000, 1000}; // r007

  
  // pull data
  Long64_t entry_q;
  Long64_t entry_p[2];
  std::vector <Int_t> *det_q = 0;
  std::vector <Int_t> *det_p = 0;
  std::vector <Long64_t> *time_q = 0;
  std::vector <Long64_t> *time_p = 0;
  std::vector <Double_t> *qdc_q = 0;
  std::vector <Double_t> *qdc_p = 0;
  
  TBranch *b_entry_q, *b_entry_p, *b_det_q, *b_det_p, *b_time_q, *b_time_p, *b_qdc_q, *b_qdc_p;

  TChain *f_in = new TChain ("t_qdc_unsort");
  f_in->Add ("../Unpacker/goPixel.root");
  f_in->SetBranchAddress ("entry_q"   , &entry_q , &b_entry_q);
  f_in->SetBranchAddress ("entry_p[2]", entry_p  , &b_entry_p);
  f_in->SetBranchAddress ("det_q"     , &det_q   , &b_det_q);   // 1..4
  f_in->SetBranchAddress ("det_p"     , &det_p   , &b_det_p); 
  f_in->SetBranchAddress ("time_q"    , &time_q  , &b_time_q); 
  f_in->SetBranchAddress ("time_p"    , &time_p  , &b_time_p); 
  f_in->SetBranchAddress ("qdc_q"     , &qdc_q   , &b_qdc_q); 
  f_in->SetBranchAddress ("qdc_p"     , &qdc_p   , &b_qdc_p); 


  // data container
  const Double_t qdc_low = 0e3;
  const Double_t qdc_high = 20e3;
  const Int_t qdc_bin = (qdc_high - qdc_low)/10;
  const Double_t bin2Qdc = (qdc_high - qdc_low)/qdc_bin;

  TH1F *h_qdc_q[4];
  for (Int_t i = 0; i < 4; i++) h_qdc_q[i] = new TH1F (Form ("h_qdc_q[%i]", i), Form ("charge %i", i + 1), qdc_bin, qdc_low, qdc_high);
  

  // data loop
  for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << "  " << std::setw (3) << (Int_t) (100.*i/f_in->GetEntries () + 0.5) << "% done" << std::endl;

      if (entry_q > 0)
	{
	  if (det_q->at (0) == 5) continue; // time sync signal
	  if (det_q->at (0) < 1) continue; // ask Marc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	  //	  std::cout << det_q->at (0) - 1 << " " << qdc_q->at (0) << std::endl;
	  h_qdc_q[det_q->at (0) - 1]->Fill (qdc_q->at (0));
	}

    }


  // scan histograms for crude fit parameters
  Double_t ampl_q[4][3], mean_q[4][3], sigma_q[4][3];
  for (Int_t i = 0; i < 1; i++) // 4
    {
      Int_t peak = 0;
      
      //      for (Int_t j = 1; j < (Int_t) h_qdc_q[i]->GetNbinsX (); j++) // avoid underflow , r001
      for (Int_t j = 50; j < (Int_t) h_qdc_q[i]->GetNbinsX (); j++) // avoid base line , r003
	{
	  if (h_qdc_q[i]->GetBinContent (j) > thresh[i])
	    {	      
	      Int_t k = 4; // avoid fluctuations
	      while (h_qdc_q[i]->GetBinContent (j + k) > thresh[i]) k++;

	      h_qdc_q[i]->GetXaxis ()->SetRange (j, j + k);
	      ampl_q[i][peak] = h_qdc_q[i]->GetMaximum ();
	      mean_q[i][peak] = qdc_low + h_qdc_q[i]->GetMaximumBin ()*bin2Qdc;
	      sigma_q[i][peak] = k/2*bin2Qdc;
	      //	      std::cout << i << " " << peak << " " << ampl_q[i][peak] << " " << mean_q[i][peak] << " " << sigma_q[i][peak] << std::endl;
	      
	      peak++;
	      j += k + 4; // avoid fluctuations
	    }
	}
      h_qdc_q[i]->GetXaxis ()->SetRange (0, qdc_bin);
    }

  
  // display
  TF1 *f_fit = new TF1 ("f_fit", my_gaus, 10.e3, 20.e3, 3);
  f_fit->SetNpx (1e3);
  f_fit->SetParLimits (0, 0., 10.e6);
  f_fit->SetParLimits (2, 0., 10.e6);
  
  TCanvas *c_qdc_q = new TCanvas ("c_qdc_q", "charge side", 700, 500);
  c_qdc_q->Draw ();
  //  c_qdc_q->Divide (2, 2);
  for (Int_t i = 0; i < 1; i++) // 4
    {
      // c_qdc_q->cd ()->SetTickx ();
      // c_qdc_q->cd ()->SetTicky ();

      
      c_qdc_q->cd (i + 1)->SetTickx ();
      c_qdc_q->cd (i + 1)->SetTicky ();
      //      c_qdc_q->cd (i + 1)->SetLogy ();
      h_qdc_q[i]->Draw ();

      // fits
      for (Int_t j = 0; j < 3; j++)
	{
	  // adjustment
	  f_fit->SetRange (mean_q[i][j] - 1.*sigma_q[i][j], mean_q[i][j] + 1.*sigma_q[i][j]);
	  f_fit->SetParameters (ampl_q[i][j], mean_q[i][j], sigma_q[i][j]);
	  //	  std::cout << " " << j << " " << f_fit->GetParameter (0) << " " << f_fit->GetParameter (1) << " " << f_fit->GetParameter (2) << std::endl;

	  h_qdc_q[i]->Fit (f_fit, "RQI");
	  ampl_q[i][j] = f_fit->GetParameter (0);
	  mean_q[i][j] = f_fit->GetParameter (1);
	  sigma_q[i][j] = f_fit->GetParameter (2);

	  // fit high energy flack
	  f_fit->SetRange (mean_q[i][j] - sigma_q[i][j], mean_q[i][j] + 2.*sigma_q[i][j]);
	  h_qdc_q[i]->Fit (f_fit, "RQI");
	  ampl_q[i][j] = f_fit->GetParameter (0);
	  mean_q[i][j] = f_fit->GetParameter (1);
	  sigma_q[i][j] = f_fit->GetParameter (2);

	  // // just to make sure	  
	  // f_fit->SetRange (mean_q[i][j] - 0.5*sigma_q[i][j], mean_q[i][j] + 2.*sigma_q[i][j]);
	  // h_qdc_q[i]->Fit (f_fit, "RQI");
	  f_fit->DrawCopy ("SAME");

	  //	  std::cout << " Q" << i + 1 << " peak " << j + 1 << ": " << 100.*f_fit->GetParameter (2)/f_fit->GetParameter (1) << " % (sigma)" << std::endl;
	  std::cout << " Q" << i + 1 << " peak " << j + 1 << ": " << 235.5*f_fit->GetParameter (2)/f_fit->GetParameter (1) << " % (FWHM)" << std::endl;
	  h_qdc_q[i]->GetXaxis ()->SetRange ((mean_q[i][0] - 10*sigma_q[i][0])/bin2Qdc, (mean_q[i][2] + 10*sigma_q[i][2])/bin2Qdc); // zoom on peaks
	}
      std::cout << "\n";
    }
  c_qdc_q->Print (Form ("run%03i_fit_charge_qdcUnsort.pdf", run_in));

  

  
  
}
