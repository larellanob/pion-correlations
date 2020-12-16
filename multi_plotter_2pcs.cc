#include "Root/Plotter.cc"

bool m_simulation = true;

void multi_plotter_2pcs()
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
     "zhpp",
     "zhpm",
     /*
     "ppp",
     "zh0002pp",
     "zh0002pm",
     "zh0203pp",
     "zh0203pm",
     "zh0305pp",
     "zh0305pm",*/
    };

  // do not change the order or comment these without changing
  // Root/Plotter.cc
  std::vector <TString> targets =
    {
     "C",
     "Fe",
     "Pb",
     "D"
    };


  std::vector <TString> var =
    {
     "phiPQboosted",
     "y",
     "thetaPQboosted",
    };
  
  TString input_dir = "/home/luciano/Physics/CLAS/pion_correlation/";
  if ( m_simulation ) {
    input_dir = input_dir+"simulations/";
  }
  
  TH1F *h1;
  TH1F *h2;
  TH1F *h3;
  std::vector<TH1F *>corr_histos;
  std::vector<TH1F *>same_histos;
  std::vector<TH1F *>mult_histos;

  TFile *f = new TFile();
  new TCanvas();

  TString stylized;
  for ( TString i: modes ) {
    for ( TString j: var ) {
      corr_histos.clear();
      same_histos.clear();
      mult_histos.clear();
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
	h2 = (TH1F*)gDirectory->Get("ns_"+j+"_"+k);
	if ( h2->GetMaximum() > rangeh23 ) {
	  rangeh23 = h2->GetMaximum();
	}
	h3 = (TH1F*)gDirectory->Get("nm_"+j+"_"+k);
	if ( h3->GetMaximum() > rangeh23 ) {
	  rangeh23 = h3->GetMaximum();
	}
	if ( k == "Pb" && j == "y" ) {
	  double aux1 = h2->GetMaximum();
	  double aux2 = h3->GetMaximum();
	  cout << Form("Carbon same max %.1f",h2->GetMaximum()) << endl;
	  cout << Form("Carbon mult max %.1f",h3->GetMaximum()) << endl;
	  cout << Form("nc: %.4f" , aux1/aux2) << endl;
	}
      }

      h1->SetMaximum(2.05);
      h1->SetMinimum(0);
      h2->SetMaximum(rangeh23+0.1*rangeh23);
      h3->SetMaximum(rangeh23+0.1*rangeh23);
      cout << TString::Format("modes <%s> var <%s> max %.1f",i.Data(),j.Data(),rangeh23) << endl;
      cout << Form("maximums h2: %.1f, h3: %.1f",h2->GetMaximum(),h3->GetMaximum())<<endl;

      for ( TString k: targets ) {
	h1 = (TH1F*)gDirectory->Get("nc_"+j+"_"+k);
	h1->Draw("E0X0 same");
	corr_histos.push_back(h1);
      }
      h1->GetYaxis()->SetTitle("N_{s}/N_{m}");
      TString out_filename = input_dir+"plots/"+i+"_"+j+"_corr.png";
      Plotter(corr_histos[0],corr_histos[1],corr_histos[2],corr_histos[3],out_filename);

      for ( TString k: targets ) {
	h2 = (TH1F*)gDirectory->Get("ns_"+j+"_"+k);
	h2->Draw("E0X0 same");
	same_histos.push_back(h2);
      }
      h2->SetMinimum(0);
      h2->GetYaxis()->SetTitle("N_{s}");
      TString out_filename2 = input_dir+"plots/"+i+"_"+j+"_same.png";
      Plotter(same_histos[0],same_histos[1],same_histos[2],same_histos[3],out_filename2);

      

      for ( TString k: targets ) {
	h3 = (TH1F*)gDirectory->Get("nm_"+j+"_"+k);
	h3->Draw("E0X0 same");
	mult_histos.push_back(h3);
      }
      h3->SetMinimum(0);
      h3->GetYaxis()->SetTitle("N_{m}");
      TString out_filename3 = input_dir+"plots/"+i+"_"+j+"_mult.png";
      Plotter(mult_histos[0],mult_histos[1],mult_histos[2],mult_histos[3],out_filename3);
    }
    
  }
}
