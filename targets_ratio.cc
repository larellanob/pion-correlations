#include "Root/RatioPlotter.cc"

bool m_simulation = false;

void targets_ratio()
{
  
  if ( m_simulation ) {
    cout << "running in simulation mode!" << endl;
  }

  std::vector <TString> modes =
    {
     //"pp",
     //"pm",
     "ptlead",
     "zhlead",
     //"zhpp",
     //"zhpm",
     "zh3",
     "zh2",
     "zh1",
     /*
     "ppp",
     "zh0002pp",
     "zh0002pm",
     "zh0203pp",
     "zh0203pm",
     "zh0305pp",
     "zh0305pm",*/
    };

  
  std::vector <TString> targets =
    {
     "C",
     "Fe",
     "Pb",
     //"D"
    };

  std::vector <TString> variables =
    {
     "phiPQboosted",
     "y",
     "thetaPQboosted",
    };
  
  TString input_dir = "/home/luciano/Physics/CLAS/pion_correlation/";
  if ( m_simulation ) {
    input_dir = input_dir+"simulations/";
  }

  TH1F hCRead;
  TH1F hPbRead;
  TH1F hFeRead;
  
  TH1F *hCSave;
  TH1F *hPbSave;
  TH1F *hFeSave;
  
  TH1F hD;
  std::vector<TH1F *>corr_histos;

  TFile *f = new TFile();
  new TCanvas();

  TString axislabel;
  for ( TString mode: modes ) {
    for ( TString var: variables ) {
      corr_histos.clear();
      if ( var == "phiPQboosted" ){
	axislabel = "#Delta#phi";
	f->Open(input_dir+"histograms/alltarg_"+mode+"__boo.root");
      } else if ( var == "y" ) {
	axislabel = "#Delta y";
	f->Open(input_dir+"histograms/alltarg_"+mode+"__rap.root");
      } else if ( var == "thetaPQboosted" ) {
	axislabel = "#Delta #theta";
	f->Open(input_dir+"histograms/alltarg_"+mode+"__boo.root");
      }
      // deuterium histogram for variable var
      // not a pointer
      hD = *(TH1F*)gDirectory->Get("nc_"+var+"_D");
      int target_counter=0;
      //hresult = (TH1F*)gDirectory->Get("nc_"+var+"_D");

      int nbins = hD.GetNbinsX();
      
      hFeRead = *(TH1F*)gDirectory->Get("nc_"+var+"_Fe");
      hPbRead = *(TH1F*)gDirectory->Get("nc_"+var+"_Pb");
      hCRead = *(TH1F*)gDirectory->Get("nc_"+var+"_C");

      hFeSave = (TH1F*)gDirectory->Get("nc_"+var+"_Fe");
      hPbSave = (TH1F*)gDirectory->Get("nc_"+var+"_Pb");
      hCSave = (TH1F*)gDirectory->Get("nc_"+var+"_C");

      hFeSave->Reset();
      hPbSave->Reset();
      hCSave->Reset();

      cout << "variable " << var << endl;
      for ( int bin = 1; bin < nbins+1; bin++ ) {
	double tvalue = 0;
	double fillvalue = 0;
	double error = 0;
	double bincenter = hD.GetXaxis()->GetBinCenter(bin);;
	double dvalue = hD.GetBinContent(bin);
	double terror, derror;
	derror = hD.GetBinError(bin);
	
	//cout << bincenter << endl;

	tvalue = hFeRead.GetBinContent(bin);
	terror = hFeRead.GetBinError(bin);
	if ( dvalue != 0 ) {
	  fillvalue = tvalue/dvalue;
	} else {
	  fillvalue = 0;
	}
	hFeSave->Fill(bincenter,fillvalue);
	if ( tvalue != 0 && dvalue != 0 ) {
	  error = (dvalue/tvalue)*sqrt( (derror*derror)/(dvalue*dvalue) + (terror*terror)/(tvalue*tvalue) );
	  hFeSave->SetBinError(bin,error);
	}
	//cout << Form("bin %i targ Fe, (%.6f/%.6f) = %.6f", bin, tvalue, dvalue,fillvalue) << endl;
	
	tvalue = hCRead.GetBinContent(bin);
	terror = hCRead.GetBinError(bin);
	if ( dvalue != 0 ) {
	  fillvalue = tvalue/dvalue;
	} else {
	  fillvalue = 0;
	}
	hCSave->Fill(bincenter,fillvalue);
	if ( tvalue != 0 && dvalue != 0 ) {
	  error = (dvalue/tvalue)*sqrt( (derror*derror)/(dvalue*dvalue) + (terror*terror)/(tvalue*tvalue) );
	  hCSave->SetBinError(bin,error);
	}
	//cout << Form("bin %i targ C, (%.6f/%.6f) = %.6f", bin, tvalue, dvalue,fillvalue) << endl;

	
	tvalue = hPbRead.GetBinContent(bin);
	terror = hPbRead.GetBinError(bin);
	if ( dvalue != 0 ) {
	  fillvalue = tvalue/dvalue;
	} else {
	  fillvalue = 0;
	}
	hPbSave->Fill(bincenter,fillvalue);
	if ( tvalue != 0 && dvalue != 0 ) {
	  error = (dvalue/tvalue)*sqrt( (derror*derror)/(dvalue*dvalue) + (terror*terror)/(tvalue*tvalue) );
	  hPbSave->SetBinError(bin,error);
	}
	//cout << Form("bin %i targ Pb, (%.6f/%.6f) = %.6f, error %.6f", bin, tvalue, dvalue,fillvalue, error) << endl;
      }
	
      //cout << "bINS : " << hC->GetNbinsX() << endl;
      //cout << hPb->GetBinContent(15) << endl;
      //hCSave->Draw("E0X0");
      corr_histos.push_back(hCSave);
      //hFeSave->Draw("E0X0");
      corr_histos.push_back(hFeSave);
      //hPbSave->Draw("E0X0");

      corr_histos.push_back(hPbSave);


      double minrange = 0;
      double maxrange = 0;
      minrange = 10000;
      for ( auto h: corr_histos ) {
	if ( h->GetMaximum() > maxrange ) {
	  maxrange = h->GetMaximum();
	}
	double minnotzero = 100000;
	for ( int bin = 1; bin < nbins+1; bin++ ) {
	  if ( h->GetBinContent(bin) < minnotzero && h->GetBinContent(bin) != 0 ) {
	    minnotzero = h->GetBinContent(bin);
	  }
	}
	if ( minnotzero < minrange ) {
	  minrange = minnotzero;
	}
      }
      cout << Form("max %.2f min %.2f",maxrange,minrange ) << endl;
      hCSave->SetMaximum(maxrange+0.2*maxrange);
      hCSave->SetMinimum(minrange-0.2*minrange);
      hCSave->SetXTitle(axislabel);
      cout << "HISTOGRAMS " << corr_histos.size() << endl;
      TString out_filename = input_dir+"plots_ratio/"+var+"_"+mode+"_Dratio.png";
      RatioPlotter(corr_histos[0],corr_histos[1],corr_histos[2],out_filename,m_simulation);
    }
    

  }
  
  
}

