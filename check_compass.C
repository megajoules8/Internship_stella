
//#include "../../Libs/Root/libPerso.h"
#include "libPerso.h"

Double_t my_gaus (Double_t *var, Double_t *par){
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2));
}


void check_compass (){

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

      h_qdc_j[i] = new TH1F (Form ("h_qdc_j[%i]", i), "", qdc_bin, qdc_low, qdc_high);
      h_qdc_j[i]->SetLineColor (2);
      h_qdc_q[i] = new TH1F (Form ("h_qdc_q[%i]", i), "", qdc_bin, qdc_low, qdc_high);
      h_qdc_q[i]->SetLineColor (2);

      h_diff[i] = new TH1F (Form ("h_diff[%i]", i), "", diff_bin, diff_low, diff_high);
    }
  h_qdc_j[3] = new TH1F ("h_qdc_j[3]", "", 2*qdc_bin, qdc_low, 2*qdc_high);
  h_qdc_j[3]->SetLineColor (2);
  h_qdc_qj = new TH2F ("h_qdc_qj", "", 2*qdc_bin, qdc_low, 2*qdc_high, qdc_bin, qdc_low, qdc_high);

  
  // data loop
  for (Long64_t i = 0; i < (Long64_t) f_in->GetEntries () - 10; i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << "  " << std::setw (3) << (Int_t) (100.*i/f_in->GetEntries () + 0.5) << "% done" << std::endl;

      h_qdc[(Int_t) Channel/2]->Fill (Energy); // QDC spectra straight

      
      // make some events
      std::vector <Long64_t> time_coi;
      time_coi.push_back (Timestamp);
      std::vector <Short_t> qdc_coi, chn_coi;
      qdc_coi.push_back (Energy);
      chn_coi.push_back (Channel);
      do
      	{
      	  f_in->GetEntry (i + (Int_t) time_coi.size ());
      	  time_coi.push_back (Timestamp);
      	  qdc_coi.push_back (Energy);
      	  chn_coi.push_back (Channel);

      	  // std::cout << "  " << i << " " << time_coi.size () << " " << time_coi.back () << " " <<  time_coi.at (0) << " " << time_coi.back () - time_coi.at (0) << " " <<
      	  //   chn_coi.back () << " " << chn_coi.at (0) << " " <<
      	  //   qdc_coi.back () << " " << qdc_coi.at (0)  << std::endl;
      	}
      while (time_coi.back () - time_coi.at (0) < 400e3);

      
      if ((Int_t) time_coi.size () > 2) // analyze coincident events
      	{
	  // last element isn't coincident
	  time_coi.pop_back ();
      	  qdc_coi.pop_back ();
      	  chn_coi.pop_back ();
	  
      	  for (Int_t j = 1; j < (Int_t) time_coi.size (); j++) h_qdc[(Int_t) chn_coi.at (j)/2]->Fill (qdc_coi.at (j)); // QDC spectra straight (which we'd jump over)
	  
	  
	  if  ((Int_t) time_coi.size () == 2) // two coincident
	    {
	      if (chn_coi.at (0) == 2 || chn_coi.at (1) == 2) // coincidence of Ji + Q1
		{
		  for (Int_t j = 0; j < 2; j++) h_qdc_q[(Int_t) chn_coi.at (j)/2]->Fill (qdc_coi.at (j));
		}
	      else //  coincidence of J1 + J2
		{
		  for (Int_t j = 0; j < 2; j++) h_qdc_j[(Int_t) chn_coi.at (j)/2]->Fill (qdc_coi.at (j));
		}
	    }
      	  else if ((Int_t) time_coi.size () == 3) // three coincident
      	    {
      	      if (chn_coi.at (0) == chn_coi.at (1) || chn_coi.at (0) == chn_coi.at (2) || chn_coi.at (1) == chn_coi.at (2)) continue; // coincidence of J1 + J2 + Q1

	      // sort according channel number: 0, 2, 4
	      Bool_t iSwapped;
	      do
	      	{
	      	  iSwapped = false;
	      	  if (chn_coi.at (0) > chn_coi.at (1))
	      	    {
	      	      std::swap(chn_coi.at (0), chn_coi.at (1));
	      	      std::swap(qdc_coi.at (0), qdc_coi.at (1));
	      	      std::swap(time_coi.at (0), time_coi.at (1));
	      	      iSwapped = true;
	      	    }
	      	  if (chn_coi.at (1) > chn_coi.at (2))
	      	    {
	      	      std::swap(chn_coi.at (1), chn_coi.at (2));
	      	      std::swap(qdc_coi.at (1), qdc_coi.at (2));
	      	      std::swap(time_coi.at (1), time_coi.at (2));
	      	      iSwapped = true;
	      	    }
	      	}
	      while (iSwapped);

	      h_diff[0]->Fill ((time_coi.at (0) - time_coi.at (1))/1.e3); // time difference
	      h_diff[1]->Fill ((time_coi.at (1) - time_coi.at (2))/1.e3);
	      h_diff[2]->Fill ((time_coi.at (0) - time_coi.at (2))/1.e3);

	      h_qdc_j[3]->Fill (qdc_coi.at (0) + qdc_coi.at (2)); // total charge junction side
	      h_qdc_qj->Fill (qdc_coi.at (0) + qdc_coi.at (2), qdc_coi.at (1));
	      
	      
	      for (Int_t j = 0; j < 3; j++) h_qdc[(Int_t) chn_coi.at (j)/2 + 3]->Fill (qdc_coi.at (j)); // triple coincident spectra
      	    }
	  	  
      	  i += (Int_t) time_coi.size () - 1; // increment at beginning of next turn
      	}
    }
  
  
  // display
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
  std::cout << " 02, 3coi in ROI: " << h_qdc[4]->Integral () << std::endl;
  h_qdc[2]->GetXaxis ()->SetRangeUser (roi04_low, roi04_high);
  h_qdc[5]->GetXaxis ()->SetRangeUser (roi04_low, roi04_high);
  std::cout << " 04, all  in ROI: " << h_qdc[2]->Integral () << std::endl;
  std::cout << " 04, 3coi in ROI: " << h_qdc[5]->Integral () << std::endl;
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

  
  TCanvas *c_qdc_q = new TCanvas ("c_qdc_q", "c_qdc_q", 700, 500);
  c_qdc_q->Draw ();
  c_qdc_q->Divide (2, 2);
  for (Int_t i = 0; i < 3; i++)
    {
      c_qdc_q->cd (i + 1)->SetTickx ();
      c_qdc_q->cd (i + 1)->SetTicky ();
      //      c_qdc_q->cd (i + 1)->SetLogy ();
      h_qdc[i]->Draw ();
      h_qdc_q[i]->Draw ("SAME");
    }
  c_qdc_q->Print (Form ("run%03i_check_compass_c_qdc_q.pdf", n_run));

  
  TCanvas *c_qdc_j = new TCanvas ("c_qdc_j", "c_qdc_j", 700, 500);
  c_qdc_j->Draw ();
  c_qdc_j->Divide (2, 2);
  for (Int_t i = 0; i < 4; i++)
    {
      c_qdc_j->cd (i + 1)->SetTickx ();
      c_qdc_j->cd (i + 1)->SetTicky ();
      //      c_qdc_j->cd (i + 1)->SetLogy ();
      if (i < 3)
	{
	  h_qdc[i]->Draw ();
	  h_qdc_j[i]->Draw ("SAME");
	}
      else h_qdc_j[i]->Draw ();
    }
  c_qdc_j->Print (Form ("run%03i_check_compass_c_qdc_j.pdf", n_run));

  
  TCanvas *c_diff = new TCanvas ("c_diff", "c_diff", 700, 500);
  c_diff->Draw ();
  c_diff->Divide (2, 2);

  h_diff[0]->SetTitle ("[0] - [2]");
  h_diff[1]->SetTitle ("[2] - [4]");
  h_diff[2]->SetTitle ("[0] - [4]");
  for (Int_t i = 0; i < 3; i++)
    {
      c_diff->cd (i + 1)->SetTickx ();
      c_diff->cd (i + 1)->SetTicky ();
      //      c_diff->cd (i + 1)->SetLogy ();

      h_diff[i]->GetXaxis ()->SetTitle ("#Delta t [ns]");
      h_diff[i]->Draw ();

     // if (i == 2) // fit the strip resolution
      //	{
	  TF1 *f_fit = new TF1 ("f_fit", my_gaus, diff_low, diff_high, 3);
	  f_fit->SetNpx (1e3);
	  f_fit->SetParLimits (0, 0., 10.e6);
	  f_fit->SetParLimits (2, 0., 10.e6);

	  h_diff[i]->Fit (f_fit, "RQI");
	  f_fit->Draw ("SAME");
      //	}
    }
  c_diff->Print (Form ("run%03i_check_compass_c_diff.pdf", n_run));

  
  TCanvas *c_corr = new TCanvas ("c_corr", "c_corr", 700, 500);
  c_corr->Draw ();
  c_corr->cd ()->SetTickx ();
  c_corr->cd ()->SetTicky ();
  //      c_corr->cd ()->SetLogy ();

  h_qdc_qj->GetYaxis ()->SetTitle ("QDC junction");
  h_qdc_qj->GetXaxis ()->SetTitle ("QDC charge");
  h_qdc_qj->Draw ("COLZ");
  
  c_corr->Print (Form ("run%03i_check_compass_c_corr.pdf", n_run));


  
}
