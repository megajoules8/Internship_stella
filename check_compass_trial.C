 #include "libPerso.h"
//implement event recontruction
//imlement automatic fitting
Double_t my_gaus (Double_t *var, Double_t *par){
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2));
}


void trial (){

  //  gStyle->SetOptStat (1111111);
  gStyle->SetOptStat (0);
  gStyle->SetOptFit (1);
  

  std::cout << "\n Check the spectra and the event building using the CoMPASS DAQ from Caen.\n\n";
  const Int_t n_run = 9;
  const Double_t roi00_low = 6.e3;   // QDC ROI chn 00
  const Double_t roi00_high = 12.e3;
  const Double_t roi02_low = 10.e3;  // 02
  const Double_t roi02_high = 16.e3;
  const Double_t roi04_low = 7.e3;
  const Double_t roi04_high = 13.e3;
  const Int_t thresh[4] = {400, 200, 200, 200};
  
  Short_t Channel, Energy; // 16 bit int
  Long64_t Timestamp;      // 64 bit int
  
  TChain *f_in = new TChain ("Data_R");
  f_in->Add ("/home/megajoules/projects/check_compass/run01_2/RAW/SDataR_run01_2.root");
  //f_in->Add ("/home/mheine/Work/Tmp/run01_2/RAW/SDataR_run01_1.root");
  f_in->SetBranchAddress ("Energy"   , &Energy );
  f_in->SetBranchAddress ("Channel"  , &Channel );  // 0, 2, 4
  f_in->SetBranchAddress ("Timestamp", &Timestamp); // ps


  // data container
  const Double_t qdc_low = 0e3;
  const Double_t qdc_high = 20e3;
  const Int_t qdc_bin = (qdc_high - qdc_low)/10;
  const Double_t diff_low = -50;
  const Double_t diff_high = 50;
  const Int_t diff_bin = (diff_high - diff_low)/1;

  TH1F *h_qdc[6], *h_qdc_j[4], *h_qdc_q[3];
  TH2F *h_qdc_qj;
  TH1F *h_diff[3];
  for (Int_t i = 0; i < 3; i++)
    {
      h_qdc[i] = new TH1F (Form ("h_qdc[%i]", i), "", qdc_bin, qdc_low, qdc_high);
      h_qdc[i + 3] = new TH1F (Form ("h_qdc[%i]", i + 3), "", qdc_bin, qdc_low, qdc_high);
      h_qdc[i + 3]->SetLineColor (2);

     
      h_diff[i] = new TH1F (Form ("h_diff[%i]", i), "", diff_bin, diff_low, diff_high);
    }
 
  
  // data loop
  for (Long64_t i = 0; i < (Long64_t) f_in->GetEntries () - 10; i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << "  " << std::setw (3) << (Int_t) (100.*i/f_in->GetEntries () + 0.5) << "% done" << std::endl;

      h_qdc[(Int_t) Channel/2]->Fill (Energy); // QDC spectra straight
    }
 //display
 
 TCanvas *c_qdc = new TCanvas ("c_qdc", "c_qdc", 700, 500);
  c_qdc->Draw ();
  c_qdc->Divide (2, 2);

  std::cout << "\n  ------------ stats  ------------" << std::endl;
  h_qdc[0]->GetXaxis ()->SetRangeUser (roi00_low, roi00_high);
  h_qdc[3]->GetXaxis ()->SetRangeUser (roi00_low, roi00_high);
  std::cout << " 00, all  in ROI: " << h_qdc[0]->Integral () << std::endl;
  std::cout << " 00, 3coi in ROI: " << h_qdc[3]->Integral () << std::endl;
  h_qdc[1]->GetXaxis ()->SetRangeUser (roi02_low, roi02_high);
  h_qdc[4]->GetXaxis ()->SetRangeUser (roi02_low, roi02_high);
  std::cout << " 02, all  in ROI: " << h_qdc[1]->Integral () << std::endl;
 // std::cout << " 02, 3coi in ROI: " << h_qdc[4]->Integral () << std::endl;
  h_qdc[2]->GetXaxis ()->SetRangeUser (roi04_low, roi04_high);
  h_qdc[5]->GetXaxis ()->SetRangeUser (roi04_low, roi04_high);
  std::cout << " 04, all  in ROI: " << h_qdc[2]->Integral () << std::endl;
 // std::cout << " 04, 3coi in ROI: " << h_qdc[5]->Integral () << std::endl;
  std::cout << "  ------------ stats  ------------\n" << std::endl;
  
  for (Int_t i = 0; i < 3; i++)
    {
      c_qdc->cd (i + 1)->SetTickx ();
      c_qdc->cd (i + 1)->SetTicky ();
      //      c_qdc->cd (i + 1)->SetLogy ();
      h_qdc[i]->GetXaxis ()->SetRangeUser (100., 16.e3);
      h_qdc[i]->Draw ();
      h_qdc[i + 3]->Draw ("SAME");
    }
  c_qdc->Print (Form ("run%03i_check_compass_c_qdc.pdf", n_run));

  //fitting
  // scan histograms for crude fit parameters
  Double_t ampl[4][3], mean[4][3], sigma[4][3];
  for (Int_t i = 0; i < 1; i++) // 4
    {
      Int_t peak = 0;
      
      //      for (Int_t j = 1; j < (Int_t) h_qdc[i]->GetNbinsX (); j++) // avoid underflow , r001
      for (Int_t j = 50; j < (Int_t) h_qdc[i]->GetNbinsX (); j++) // avoid base line , r003
	{
	  if (h_qdc[i]->GetBinContent (j) > thresh[i])
	    {	      
	      Int_t k = 4; // avoid fluctuations
	      while (h_qdc[i]->GetBinContent (j + k) > thresh[i]) k++;

	      h_qdc[i]->GetXaxis ()->SetRange (j, j + k);
	      ampl[i][peak] = h_qdc[i]->GetMaximum ();
	      mean[i][peak] = qdc_low + h_qdc[i]->GetMaximumBin ()*diff_bin;
	      sigma[i][peak] = k/2*diff_bin;
	      
          std::cout << i << " " << peak << " " << ampl[i][peak] << " " << mean[i][peak] << " " << sigma[i][peak] << std::endl;
	      
	      peak++;
	      j += k + 4; // avoid fluctuations
	    }
	}
      h_qdc[i]->GetXaxis ()->SetRange (0, qdc_bin);
    }

    
}
