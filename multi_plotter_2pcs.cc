#include "Root/Plotter.cc"

bool m_simulation = false;

void multi_plotter_2pcs()
{
  if ( m_simulation ) {
    cout << "running in simulation mode!" << endl;
  }
  std::vector <TString> modes =
    {
     //"pp",
     //"pm",
     //"ptlead",
     //"zhlead",
     /*
     "pt3",
     "pt2",
     "pt1",
     */
     /*
     "zh3",
     "zh2",
     "zh1",
     */
     /*
     "zh3ycut",
     "zh2ycut",
     "zh1ycut",
     */
     "pt3ycut",
     "pt2ycut",
     "pt1ycut",

       //"zhpp",
     //"zhpm",
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
     "phiPQboosted2",
     "y",
     //"thetaPQboosted",
    };
  
  TString input_dir = "/home/luciano/Physics/CLAS/pion_correlation/";
  if ( m_simulation ) {
    input_dir = input_dir+"simulations/";
  }
  
  TH1F *h1;
  TH1F *h2;
  TH1F *h3;
  TH2F *ridge;

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
      if ( j == "phiPQboosted2" ) {
	stylized = "#Delta#phi_{PQ} (deg)";
	f->Open(input_dir+"histograms/alltarg_"+i+"__rap.root");
      } else if ( j == "y" ) {
	stylized = "#Delta y";
	f->Open(input_dir+"histograms/alltarg_"+i+"__rap.root");
      } else if ( j == "thetaPQboosted" ) {
	stylized = "#Delta#theta_{PQ} (deg)";
	f->Open(input_dir+"histograms/alltarg_"+i+"__boo.root");
      }
      int target_counter=0;
      
      double rangeh1max = 0;
      double rangeh1min = 1000000;
      double rangeh23 = 0;
      for ( TString k: targets ) {
	h1 = (TH1F*)gDirectory->Get("nc_"+j+"_"+k);
	if ( h1->GetMaximum() > rangeh1max ) {
	  rangeh1max = h1->GetMaximum();
	}
	if ( h1->GetMinimum() < rangeh1min ) {
	  rangeh1min = h1->GetMinimum();
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

      //h1->SetMaximum(2.05);
      /*
      h1->SetMaximum(rangeh1max+0.05*rangeh1max);
      h1->SetMinimum(rangeh1min-0.05*rangeh1min);
      */

      h1->SetMaximum(0.40);
      h1->SetMinimum(0.02);

      h1->GetXaxis()->SetTitle(stylized);
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
      Plotter(corr_histos[0],corr_histos[1],corr_histos[2],corr_histos[3],out_filename,i);

      for ( TString k: targets ) {
	h2 = (TH1F*)gDirectory->Get("ns_"+j+"_"+k);
	h2->Draw("E0X0 same");
	same_histos.push_back(h2);
      }
      h2->SetMinimum(0);
      h2->GetYaxis()->SetTitle("N_{s}");
      TString out_filename2 = input_dir+"plots/"+i+"_"+j+"_same.png";
      Plotter(same_histos[0],same_histos[1],same_histos[2],same_histos[3],out_filename2,i);

      

      for ( TString k: targets ) {
	h3 = (TH1F*)gDirectory->Get("nm_"+j+"_"+k);
	h3->Draw("E0X0 same");
	mult_histos.push_back(h3);
      }
      h3->SetMinimum(0);
      h3->GetYaxis()->SetTitle("N_{m}");
      TString out_filename3 = input_dir+"plots/"+i+"_"+j+"_mult.png";
      Plotter(mult_histos[0],mult_histos[1],mult_histos[2],mult_histos[3],out_filename3,i);
    }
    // f->close??
    f->Open(input_dir+"histograms/alltarg_"+i+"__rap.root");


    
    for ( TString targ: targets ) {
      // same
      ridge = (TH2F*)gDirectory->Get("mirrored_n_sphiPQboosted2_y_"+targ);
      //void ridge_plot(TH2F h2, TString pre, TString post,TString target)
      TString pre = input_dir+"plots_ridge2/";
      //TString axi = (TString)ridge->GetXaxis()->GetTitle()+"_"+(TString)ridge->GetYaxis()->GetTitle();
      TString axi = "phi-y";
      TString pos = i+"_"+targ+"_"+axi+"_same.png";
      ridge_plot(*ridge,pre,pos,targ,i);

      // multi
      ridge = (TH2F*)gDirectory->Get("mirrored_n_mphiPQboosted2_y_"+targ);
      pos = i+"_"+targ+"_"+axi+"_mult.png";
      ridge_plot(*ridge,pre,pos,targ,i);

      // corr
      ridge = (TH2F*)gDirectory->Get("mirrored_n_cphiPQboosted2_y_"+targ);
      pos = i+"_"+targ+"_"+axi+"_corr.png";
      ridge_plot(*ridge,pre,pos,targ,i);
    }
    

   
  }
}
