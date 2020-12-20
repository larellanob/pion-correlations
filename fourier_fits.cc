#include "Root/Plotter.cc"

bool m_simulation = false;

Double_t fourier( Double_t *x, Double_t *par ) {
  Double_t xx = x[0];
  // third order fourier
  Double_t fr3 = par[0] + par[1]*cos(xx*(TMath::Pi()/180.)) + par[2]*cos(2.* xx*(TMath::Pi()/180.))+ par[3]*cos(3.* xx*(TMath::Pi()/180.));
  //Double_t fr3 = par[0] + par[1]*sin(xx) + par[2]*sin(2.* xx) + par[3]*sin(3.* xx);
  return fr3;
}


void fourier_fits()
{
  if ( m_simulation ) {
    cout << "running in simulation mode!" << endl;
  }
  std::vector <TString> modes =
    {
     //"pp",
     //"pm",
     "ptlead",
     //"zhlead",
     //"zhpp",
     //"zhpm",
     };

  // do not change the order or comment these without changing
  // Root/Plotter.cc
  std::vector <TString> targets =
    {
     //"C",
     //"Fe",
     //"Pb",
     "D"
    };


  std::vector <TString> var =
    {
     "phiPQboosted",
     //"y",
     //"thetaPQboosted",
    };
  
  TString input_dir = "/home/luciano/Physics/CLAS/pion_correlation/";
  if ( m_simulation ) {
    input_dir = input_dir+"simulations/";
  }
  
  TH1F *h1;
  TH1F *h2;
  TH1F *h3;
  std::vector<TH1F *>corr_histos;

  TFile *f = new TFile();
  new TCanvas();

  TString stylized;
  for ( TString i: modes ) {
    for ( TString j: var ) {
      corr_histos.clear();
      if ( j == "phiPQboosted" ) {
	stylized = "#phi_{PQ}";
	f->Open(input_dir+"histograms/alltarg_"+i+"__boo.root");
      } else if ( j == "y" ) {
	stylized = "y";
	f->Open(input_dir+"histograms/alltarg_"+i+"__rap.root");
      } else if ( j == "thetaPQboosted" ) {
	stylized = "#theta_{PQ}";
	f->Open(input_dir+"histograms/alltarg_"+i+"__boo.root");
      }
      int target_counter=0;
      
      double rangeh1 = 0;
      double rangeh23 = 0;
      for ( TString k: targets ) {
	h1 = (TH1F*)gDirectory->Get("nc_"+j+"_"+k);
	if ( h1->GetMaximum() > rangeh1 ) {
	  rangeh1 = h1->GetMaximum();
	}
      }

      for ( TString k: targets ) {
	h1 = (TH1F*)gDirectory->Get("nc_"+j+"_"+k);
	h1->Draw("E0X0 same");
	corr_histos.push_back(h1);
	// fit between 0.1 < Delta phi < 2.0, 4 parameters (n=3)
	TF1 * fit1 = new TF1("fit1",fourier,5.8,115.,4);


	h1->GetYaxis()->SetTitle("N_{s}/N_{m}");
	TH1F *h1clon = (TH1F*) h1->Clone();
	h1->Fit("fit1","R");
	h1->Draw();

	
	double v0 = fit1->GetParameter(0);
	double v1 = fit1->GetParameter(1);
	double v2 = fit1->GetParameter(2);
	double v3 = fit1->GetParameter(3);
	cout << v2 << endl;

	TF1 *fou0 = new TF1("fou0", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fou0->SetParameters(v0,0,0);
	TF1 *fou1 = new TF1("fou1", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fou1->SetParameters(v0,v1,1);
	TF1 *fou2 = new TF1("fou2", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fou2->SetParameters(v0,v2,2);
	TF1 *fou3 = new TF1("fou3", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fou3->SetParameters(v0,v3,3);

	new TCanvas();
	h1clon->GetXaxis()->SetRangeUser(-90.,270.);

	cout << v0 << endl;
	fou0->Draw();
	fou0->SetLineColor(kMagenta);
	h1clon->Draw("same");
	fou1->Draw("same");
        fou1->SetLineColor(kBlue);
	fou2->Draw("same");
	fou2->SetLineColor(kRed);
	fou3->Draw("same");
	fou3->SetLineColor(kGreen);

	fit1->Draw("same");
	fit1->SetRange(0,180);
	fit1->SetLineColor(kBlack);
	/*
	new TCanvas();
	h1clon->Draw();
	fou->SetParameters(v0,v1,1);
	fou->Draw("same");
	*/
      }

      
      

      //TString out_filename3 = input_dir+"plots/"+i+"_"+j+"_mult.png";
      //Plotter(mult_histos[0],mult_histos[1],mult_histos[2],mult_histos[3],out_filename3);
      
    }
    
  }
}
