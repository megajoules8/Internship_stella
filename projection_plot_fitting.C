//event reconstruction
#include "libPerso.h"
//implement event recontruction in a later iteration of this script
//imlement automatic fitting in this script on only the qdc side
Double_t my_gaus (Double_t *var, Double_t *par){ //like y(x,t)
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2)); // x,n
}



void projection_plot_fitting(){

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
  const Int_t thresh[1] = {20}; //saying fitting 1 histogram, so only need one threshold, changed from 100 on 20/04
  
  UShort_t Energy, Channel; // 16 bit int
  ULong64_t Timestamp;     // 64 bit int
  
  TChain *f_in = new TChain ("Data_R");
  f_in->Add ("/home/megajoules/projects/check_compass/run_1/RAW/SDataR_run_1.root");
  //f_in->Add ("/home/megajoules/projects/check_compass/run_2/RAW/SDataR_run_2.root");
  //f_in->Add ("/home/mheine/Work/Tmp/run01_2/RAW/SDataR_run01_1.root");
  f_in->SetBranchAddress ("Energy"   , &Energy ); //Energy in the raw data is a misnomer! it refers to the channel, i.e. the uncalibrated value
  f_in->SetBranchAddress ("Channel"  , &Channel );  // 0, 2, 4 corresponding to the junction 1, 2, and the QDC channel
  f_in->SetBranchAddress ("Timestamp", &Timestamp); // ps


  // data container
  const Double_t qdc_low = 0e3;
  const Double_t qdc_high = 25e3; //changed from 20e3 on 20/04
  const Int_t qdc_bin = (qdc_high - qdc_low)/10;
  const Double_t bin2Qdc = (qdc_high - qdc_low)/qdc_bin;

  //for calibrated energy
  const Double_t Energy_low = 0e3; //correct ranges here
  const Double_t Energy_high = 8000;
  const Int_t Energy_bin = (Energy_high - Energy_low)/10; //choose appropriate binning
  const Double_t bin2Energy = (Energy_high - Energy_low)/Energy_bin;
  std:: cout << bin2Qdc << endl;
  //const Float_t Shift = 181.063; //shift and slope privided by energy_calibration.C
  //const Float_t Slope = 0.413141;
  const Float_t Shift = 183.115; //shift and slope privided by energy_calibration.C run_01
  const Float_t Slope = 0.411245;
  const Float_t Angle = 0.9837155236; 
  const Double_t sin_Angle = 0.8325612636;
  const Double_t cos_Angle = 0.5539329764; 

 

  TH1F *h_qdc[6], *h_qdc_jj[4], *h_Energy_calib[3], *h_Energy_Ohmic_junction[1]; //declaring TH1f variables [which is an array of objects of the class with a specified length]
  TH2F *h_qdc_qj;
  //h_qdc is an array with a size 6, each of those 6 slots is a pointer to an object, those objects are in the class TH1F
 
 
  for (Int_t i = 0; i < 3; i++){
      h_qdc[i] = new TH1F (Form ("h_qdc[%i]", i), "", qdc_bin, qdc_low, qdc_high);
      h_qdc[i + 3] = new TH1F (Form ("h_qdc[%i]", i + 3), "", qdc_bin, qdc_low, qdc_high); //initialises 3, 4, 5
      h_qdc[i + 3]->SetLineColor (2);

      h_qdc_jj[i] = new TH1F (Form ("h_qdc_jj[%i]", i), "", qdc_bin, qdc_low, qdc_high);
      h_qdc_jj[i]->SetLineColor (2);
    }
  h_qdc_jj[3] = new TH1F ("h_qdc_jj[3]", "", 2*qdc_bin, qdc_low, 2*qdc_high);
  h_qdc_jj[3]->SetLineColor (2);
  h_qdc_qj = new TH2F ("h_qdc_qj", "", 2*qdc_bin, qdc_low, 2*qdc_high, qdc_bin, qdc_low, qdc_high);
  h_Energy_Ohmic_junction[1] = new TH1F ("h_Energy_Ohmic_junction", "", 2*qdc_bin, qdc_low, 2*qdc_high);
  //calibrated energy
  for (Int_t i = 0; i < 3; i++){
      h_Energy_calib[i] = new TH1F (Form ("h_Energy_calib[%i]", i), "", Energy_bin, Energy_low, Energy_high);
  }

  // data loop
  for (Long64_t i = 0; i < (Long64_t) f_in->GetEntries () - 10; i++)
    {
      f_in->GetEntry (i);
      if (i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << "  " << std::setw (3) << (Int_t) (100.*i/f_in->GetEntries () + 0.5) << "% done" << std::endl;

      h_qdc[(Int_t) Channel/2]->Fill (Energy); // QDC spectra straight

      h_Energy_calib[(Int_t) Channel/2]->Fill (Shift + Slope*Energy); //Energy spectra

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

      	 /*if ( i < 50 )
        {
         std::cout << "  " << i << " " << time_coi.size () << " " << time_coi.back () << " " <<  time_coi.at (0) << " " << time_coi.back () - time_coi.at (0) << " " <<
      	 chn_coi.back () << " " << chn_coi.at (0) << " " <<
      	 qdc_coi.back () << " " << qdc_coi.at (0)  << std::endl;
        }*/
          
        }
          while (time_coi.back () - time_coi.at (0) < 400e3);

         if ((Int_t) time_coi.size () > 2) // analyze coincident events

      {

          // last element isn't coincident
	        time_coi.pop_back ();
      	  qdc_coi.pop_back ();
      	  chn_coi.pop_back ();
	  
      	  for (Int_t j = 1; j < (Int_t) time_coi.size (); j++) { //need same number of counts in coi and non-coi hists
            h_qdc[(Int_t) chn_coi.at (j)/2]->Fill (qdc_coi.at (j)); // QDC spectra straight (which we'd jump over)
            h_Energy_calib[(Int_t) chn_coi.at (j)/2]->Fill (Shift + Slope*(qdc_coi.at (j))); //Energy spectra 
          }
          if  ((Int_t) time_coi.size () == 2) // two coincident
	        {
	        if (chn_coi.at (0) == 2 || chn_coi.at (1) == 2) // coincidence of Ji + Q1
		      {
		        for (Int_t j = 0; j < 2; j++) h_qdc_jj[(Int_t) chn_coi.at (j)/2]->Fill (qdc_coi.at (j));
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

	      h_qdc_jj[3]->Fill (qdc_coi.at (0) + qdc_coi.at (2)); // total charge junction side
	      h_qdc_qj->Fill (qdc_coi.at (1), qdc_coi.at (0) + qdc_coi.at (2)); //x,y
        h_Energy_Ohmic_junction[1]->Fill (sqrt((sin_Angle*qdc_coi.at (1)*sin_Angle*qdc_coi.at (1)) + (cos_Angle*(qdc_coi.at (0) + qdc_coi.at (2)))*(cos_Angle*(qdc_coi.at (0) + qdc_coi.at (2)))));
        
	      
	      for (Int_t j = 0; j < 3; j++) { // triple coincident spectra
          h_qdc[(Int_t) chn_coi.at (j)/2 + 3]->Fill (qdc_coi.at (j));
        } // triple coincident spectra
        //one line loop, can put a single statement but better for readability to not
	  	  
      	  i += (Int_t) time_coi.size () - 1; // increment at beginning of next turn, += refers to how much you add each increment 
            
      	}
               
      }   
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
      h_qdc[i]->GetXaxis ()->SetRangeUser (100., 16.e3);
      h_qdc[i]->Draw ();
      h_qdc[i + 3]->Draw ("SAME");
    }
  //c_qdc->Print (Form ("run%03i_check_compass_c_qdc.pdf", n_run));

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
	  h_qdc_jj[i]->Draw ("SAME");
	}
      else h_qdc_jj[i]->Draw ();
    }
  c_qdc_j->Print (Form ("run%03i_check_compass_c_qdc_j.pdf", n_run));

  //Calibrated Energy Spectra
   TCanvas *c_Energy_calib = new TCanvas ("c_Energy", "c_Energy", 700, 500);
  c_Energy_calib->Draw ();

  
      c_Energy_calib->cd ()->SetTickx ();
      c_Energy_calib->cd ()->SetTicky ();
      h_Energy_calib[1]->GetXaxis ()->SetRangeUser (100., 8000.); //ohmic side =2 then divided by 2
      h_Energy_calib[1]->Draw ();
  




 //fitting
  // scan histograms for crude fit parameters
  Double_t ampl[4][3], mean[4][3], sigma[4][3];
 for (Int_t i = 0; i < 4; i++) // 4 //change 20/04
    {
      Int_t peak = 0;
      
      //      for (Int_t j = 1; j < (Int_t) h_qdc[i]->GetNbinsX (); j++) // avoid underflow , r001
      for (Int_t j = 500; j < (Int_t) h_Energy_Ohmic_junction[1]->GetNbinsX (); j++) // avoid base line , r003 //changed j from 500 on 20/04 caused crashing
	    {
	        if (h_Energy_Ohmic_junction[1]->GetBinContent (j) > thresh[i])
	        {	      
	         Int_t k = 4; // avoid fluctuations
	         while (h_Energy_Ohmic_junction[i]->GetBinContent (j + k) > thresh[i]) k++;

	         h_Energy_Ohmic_junction[i]->GetXaxis ()->SetRange (j, j + k);
	         ampl[i][peak] = h_Energy_Ohmic_junction[1]->GetMaximum ();
	         mean[i][peak] = qdc_low + h_Energy_Ohmic_junction[1]->GetMaximumBin ()*bin2Qdc;
	         sigma[i][peak] = k/2*bin2Qdc;
	      
          //std::cout << i << " " << peak << " " << ampl[i][peak] << " " << mean[i][peak] << " " << sigma[i][peak] << ""<< j << ""<< j+k << std::endl;
	        
	        peak++;
	        j += k + 4; // avoid fluctuations
	        }
	    }
      h_Energy_Ohmic_junction[1]->GetXaxis ()->SetRange (0, qdc_bin);
    
      
    }

     TCanvas *c_Energy_Ohmic_junction = new TCanvas ("c_Energy_Ohmic_junction", "c_Energy_Ohmic_junction", 700, 500);
    c_Energy_Ohmic_junction->Draw (); 
    c_Energy_Ohmic_junction->cd ()->SetTickx ();
    c_Energy_Ohmic_junction->cd ()->SetTicky ();
    h_Energy_Ohmic_junction[1]->GetXaxis ()->SetRangeUser (0., 25000.); 
   // h_Energy_Ohmic_junction[1]->Draw ();

    h_qdc_qj->GetYaxis ()->SetTitle ("QDC junction"); //fix this
    h_qdc_qj->GetXaxis ()->SetTitle ("QDC charge");
  
   //fitting display
   TF1 *f_fit = new TF1 ("f_fit", my_gaus, 10.e3, 25.e3, 3); //instantiating new TF1 object with the parameters it needs to do that, returns the pointer f_fit //*changed 20/04
     f_fit->SetNpx (1e3); 
     f_fit->SetParLimits (0, 0., 10.e6);
     f_fit->SetParLimits (2, 0., 10.e6);

     for (Int_t i = 1; i < 2; i++) // select histogram on canvas for fitting
    {
    
      
      c_Energy_Ohmic_junction->SetTickx ();
      c_Energy_Ohmic_junction->SetTicky ();
      //h_Energy_Ohmic_junction[i]->DrawCopy ("SAME");
      h_Energy_Ohmic_junction[1]->Draw ();
      // fits
      for (Int_t j = 0; j < 1; j++){ //select peak numbers 0, 1, and 2
	        // adjustment
	        f_fit->SetRange (mean[i][j] - 1.*sigma[i][j], mean[i][j] + 1.*sigma[i][j]);
	        f_fit->SetParameters (ampl[i][j], mean[i][j], sigma[i][j]);
	        //	  std::cout << " " << j << " " << f_fit->GetParameter (0) << " " << f_fit->GetParameter (1) << " " << f_fit->GetParameter (2) << std::endl;

	        h_Energy_Ohmic_junction[1]->Fit (f_fit, "RQI");
	        ampl[i][j] = f_fit->GetParameter (0);
	        mean[i][j] = f_fit->GetParameter (1);
	        sigma[i][j] = f_fit->GetParameter (2);

	        // fit high energy flack
	        f_fit->SetRange (mean[i][j] - sigma[i][j], mean[i][j] + 2.*sigma[i][j]);
	        h_Energy_Ohmic_junction[i]->Fit (f_fit, "RQI");
	        ampl[i][j] = f_fit->GetParameter (0);
	        mean[i][j] = f_fit->GetParameter (1);
	        sigma[i][j] = f_fit->GetParameter (2);
          f_fit->DrawCopy ("SAME");
            std::cout << " " <<ampl[i][j] << " " << mean[i][j] << " " << sigma[i][j] << std::endl;
        }
    } 
    c_Energy_Ohmic_junction->Print (Form ("run%03i_check_compass_c_Energy_Ohmic_junction.pdf", n_run));

   

   
}

//uncalibrated ohmic side fitting:
//peak1:
//peak2:
//peak3: 89.92=/- 3.52, 13650+/-3.347, 57.74+/-3.20
