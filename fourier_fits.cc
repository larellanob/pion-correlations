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

  /*
  gStyle->SetCanvasDefH(400);
  gStyle->SetCanvasDefW(500);
  */
  //gStyle->SetPadRightMargin(1.5);
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
     "zh3",
     "zh2",
     "zh1",
     */
     "zh3ycut",
     "zh2ycut",
     "zh1ycut",
     /*
     "pt3",
     "pt2",
     "pt1",
     */
     //"zhpp",
     //"zhpm",
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

  
  TString stylized;
  for ( TString i: modes ) {
    for ( TString j: var ) {
      corr_histos.clear();
      if ( j == "phiPQboosted2" ) {
	stylized = "#phi_{PQ}";
	f->Open(input_dir+"histograms/alltarg_"+i+"__rap.root");
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
	TCanvas *c2 = new TCanvas();
	
	
	h1 = (TH1F*)gDirectory->Get("nc_"+j+"_"+k);
	h1->Draw("E0X0 same");
	corr_histos.push_back(h1);
	// fit between 0.1 < Delta phi < 2.0, 4 parameters (n=3)
	//TF1 * fit1 = new TF1("fit1",fourier,5.8,115.,4);
	//TF1 * fit1 = new TF1("fit1",fourier,5.8,150.,4);
	TF1 * fit1 = new TF1("fit1",fourier,5.8,180.,4);
	
	

	TH1F *h1clon = (TH1F*) h1->Clone();
	h1->Fit("fit1","R");


	cout << Form("| %s | %s | %.2e \\pm %.2e |",k.Data(),i.Data(),fit1->GetParameter(2),fit1->GetParError(2)) << endl;

	
	h1->SetMarkerStyle(kFullCircle);
	h1->Draw("E0X0");
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("#Delta#phi (deg)");
	h1->GetYaxis()->SetTitle("N_{s}/N_{m}");

	TString sector;
	if  ( i  == "zhlead" ) {
	  sector = "z_{h}^{assoc.} < 0.5";
	} else if  ( i  == "zh3" ) {
	  sector = "0.25 < z_{h}^{assoc.} < 0.4";
	} else if  ( i  == "zh2" ) {
	  sector = "0.12 < z_{h}^{assoc.} < 0.25";
	} else if  ( i  == "zh1" ) {
	  sector = "z_{h}^{assoc.} < 0.12";
	}
	h1->SetTitle(k+" target, "+sector+"   ");
	gStyle->SetOptTitle(0);
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleSize(0.05);
	fit1->Draw("same");
	fit1->SetTitle("Fourier fit, n #leq 3");
	
	TLegend *myleg2;
	myleg2 = c2->BuildLegend(0.2,0.7,0.5,0.9,"","");
	h1->SetStats(0);
	c2->SaveAs(input_dir+"plots_fourier/"+i+j+k+".png");
	
	double v0 = fit1->GetParameter(0);
	double v1 = fit1->GetParameter(1);
	double v2 = fit1->GetParameter(2);
	double v3 = fit1->GetParameter(3);
	cout << Form("Saved for %s %s %s with v2 %.7f",i.Data(),j.Data(),k.Data(),v2) << endl;

	TF1 *fourier0 = new TF1("fourier0", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fourier0->SetParameters(v0,0,0);
	TF1 *fourier1 = new TF1("fourier1", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fourier1->SetParameters(v0,v1,1);
	TF1 *fourier2 = new TF1("fourier2", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fourier2->SetParameters(v0,v2,2);
	TF1 *fourier3 = new TF1("fourier3", "[0]+[1]*cos([2]*x*(TMath::Pi()/180.))",-90,270);
	fourier3->SetParameters(v0,v3,3);

	TF1 *fourierall = new TF1("fourierall", "fourier2(x)");//(fourier2+fourier3)/2.0");
	
	TCanvas *c1 = new TCanvas();
	
	h1clon->GetXaxis()->SetRangeUser(-90.,270.);


	fourier1->GetXaxis()->SetTitle("#Delta#phi (deg)");
	fourier1->GetYaxis()->SetTitle("N_{s}/N_{m}");
	//fourier1->SetTitle("F_{0}");
	cout << v0 << endl;

	//fourier0->Draw();
	//fourier0->SetLineColor(kMagenta);


	fourier1->Draw();
	fourier1->SetTitle("F_{1}");
        fourier1->SetLineColor(kBlue);
	fourier2->Draw("same");
	fourier2->SetTitle("F_{2}");
	fourier2->SetLineColor(kRed);
	fourier3->Draw("same");
	fourier3->SetTitle("F_{3}");
	fourier3->SetLineColor(kGreen);

	/*
	fourierall->Draw("same");
	fourierall->SetTitle("Fourier fit 3");
	//fourierall->SetRange(0,180);
	fourierall->SetLineColor(kBlack);
	*/
	fit1->Draw("same");
	fit1->SetTitle("Fourier fit 3");
	//fit1->SetRange(0,180);
	fit1->SetLineColor(kBlack);
	
	h1clon->Draw("E0X0 same");
	h1clon->SetMarkerStyle(kFullCircle);
	h1clon->SetTitle("Data");
	double h1max =  h1clon->GetMaximum();
	double f1max =  fourier1->GetMaximum();
	if ( h1max > f1max ) {
	  fourier1->SetMaximum(h1max+0.05*h1max);
	}
	
	cout << "f1 max" << fourier1->GetMaximum() << endl;
	cout << "h1 max" << h1clon->GetMaximum() << endl;
	
	gStyle->SetOptTitle(0);
	TLegend *myleg;
	myleg = c1->BuildLegend(0.2,0.7,0.44,0.9,"","");
	//c1->SetTitle("");


	LatexText(0.55,0.33,62,"EG2 data");
	LatexText(0.695,0.33,42,"- "+k+" target");
	if ( i != "zhlead" ) {
	  LatexText(0.55,0.27,42,"z_{h}^{trig.} > 0.4");
	} else if  ( i == "zhlead" ) {
	  LatexText(0.55,0.27,42,"z_{h}^{trig.} > 0.5");
	}
	TString zhrange;
	if  ( i  == "zhlead" ) {
	  zhrange = "z_{h}^{assoc.} < 0.5";
	} else if  ( i  == "zh3" || i == "zh3ycut" ) {
	  zhrange = "0.25 < z_{h}^{assoc.} < 0.4";
	} else if  ( i  == "zh2" || i == "zh2ycut" ) {
	  zhrange = "0.12 < z_{h}^{assoc.} < 0.25";
	} else if  ( i  == "zh1" || i == "zh1ycut" ) {
	  zhrange = "z_{h}^{assoc.} < 0.12";
	}
	LatexText(0.55,0.21,42,zhrange);
	LatexText(0.55,0.15,42,"|#Deltay-y_{bias}| > 1.0");
	
	c1->SaveAs(input_dir+"plots_fourier/full_"+i+j+k+".png");
	//myleg->SetFillStyle(0);


	
	/*
	new TCanvas();
	h1clon->Draw();
	fourier->SetParameters(v0,v1,1);
	fourier->Draw("same");
	*/
      }

      
      

      //TString out_filename3 = input_dir+"plots/"+i+"_"+j+"_mult.png";
      //Plotter(mult_histos[0],mult_histos[1],mult_histos[2],mult_histos[3],out_filename3);
      
    }
    
  }
}
