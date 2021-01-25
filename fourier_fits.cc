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
     "zh1ycut",
     "zh2ycut",
     "zh3ycut",
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
     "D",
     "C",
     "Fe",
     "Pb"
    };


  std::vector <TString> variables =
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


  double v1table[4][3];
  double v1tableerr[4][3];
  double v2table[4][3];
  double v2tableerr[4][3];

  int modcount = -1;
  int tarcount = -1;
  
  TString stylized;
  for ( TString mod: modes ) {
    modcount++;
    for ( TString var: variables ) {

      corr_histos.clear();
      if ( var == "phiPQboosted2" ) {
	stylized = "#phi_{PQ}";
	f->Open(input_dir+"histograms/alltarg_"+mod+"__rap.root");
      } else if ( var == "y" ) {
	stylized = "y";
	f->Open(input_dir+"histograms/alltarg_"+mod+"__rap.root");
      } else if ( var == "thetaPQboosted" ) {
	stylized = "#theta_{PQ}";
	f->Open(input_dir+"histograms/alltarg_"+mod+"__boo.root");
      }
      int target_counter=0;
      
      double rangeh1 = 0;
      double rangeh23 = 0;
      for ( TString tar: targets ) {
	h1 = (TH1F*)gDirectory->Get("nc_"+var+"_"+tar);
	if ( h1->GetMaximum() > rangeh1 ) {
	  rangeh1 = h1->GetMaximum();
	}
      }
      tarcount = -1;
      for ( TString tar: targets ) {
	tarcount++;
	TCanvas *c2 = new TCanvas();
	
	
	h1 = (TH1F*)gDirectory->Get("nc_"+var+"_"+tar);
	h1->Draw("E0X0 same");
	corr_histos.push_back(h1);
	// fit between 0.1 < Delta phi < 2.0, 4 parameters (n=3)
	//TF1 * fit1 = new TF1("fit1",fourier,5.8,115.,4);
	//TF1 * fit1 = new TF1("fit1",fourier,5.8,150.,4);
	TF1 * fit1 = new TF1("fit1",fourier,5.8,180.,4);
	
	

	TH1F *h1clon = (TH1F*) h1->Clone();
	h1->Fit("fit1","R");

	v1table[tarcount][modcount] = fit1->GetParameter(1);
	v1tableerr[tarcount][modcount] = fit1->GetParError(1);
	v2table[tarcount][modcount] = fit1->GetParameter(2);
	v2tableerr[tarcount][modcount] = fit1->GetParError(2);
	
	cout << Form("| %s | %s | %.2e \\pm %.2e |",tar.Data(),mod.Data(),fit1->GetParameter(2),fit1->GetParError(2)) << endl;

	
	h1->SetMarkerStyle(kFullCircle);
	h1->Draw("E0X0");
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("#Delta#phi (deg)");
	h1->GetYaxis()->SetTitle("N_{s}/N_{m}");

	TString sector;
	if  ( mod == "zhlead" ) {
	  sector = "z_{h}^{assoc.} < 0.5";
	} else if  ( mod == "zh3" ) {
	  sector = "0.25 < z_{h}^{assoc.} < 0.4";
	} else if  ( mod == "zh2" ) {
	  sector = "0.12 < z_{h}^{assoc.} < 0.25";
	} else if  ( mod == "zh1" ) {
	  sector = "z_{h}^{assoc.} < 0.12";
	}
	h1->SetTitle(tar+" target, "+sector+"   ");
	gStyle->SetOptTitle(0);
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleSize(0.05);
	fit1->Draw("same");
	fit1->SetTitle("Fourier fit, n #leq 3");
	
	TLegend *myleg2;
	myleg2 = c2->BuildLegend(0.2,0.7,0.5,0.9,"","");
	h1->SetStats(0);
	c2->SaveAs(input_dir+"plots_fourier/"+mod+var+tar+".png");
	
	double v0 = fit1->GetParameter(0);
	double v1 = fit1->GetParameter(1);
	double v2 = fit1->GetParameter(2);
	double v3 = fit1->GetParameter(3);
	cout << Form("Saved for %s %s %s with v2 %.7f",mod.Data(),var.Data(),tar.Data(),v2) << endl;

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
	LatexText(0.695,0.33,42,"- "+tar+" target");
	if ( mod!= "zhlead" ) {
	  LatexText(0.55,0.27,42,"z_{h}^{trig.} > 0.4");
	} else if  ( mod== "zhlead" ) {
	  LatexText(0.55,0.27,42,"z_{h}^{trig.} > 0.5");
	}
	TString zhrange;
	if  ( mod == "zhlead" ) {
	  zhrange = "z_{h}^{assoc.} < 0.5";
	} else if  ( mod == "zh3" || mod== "zh3ycut" ) {
	  zhrange = "0.25 < z_{h}^{assoc.} < 0.4";
	} else if  ( mod == "zh2" || mod== "zh2ycut" ) {
	  zhrange = "0.12 < z_{h}^{assoc.} < 0.25";
	} else if  ( mod == "zh1" || mod== "zh1ycut" ) {
	  zhrange = "z_{h}^{assoc.} < 0.12";
	}
	LatexText(0.55,0.21,42,zhrange);
	LatexText(0.55,0.15,42,"|#Deltay-y_{bias}| > 1.0");
	
	c1->SaveAs(input_dir+"plots_fourier/full_"+mod+var+tar+".png");
	//myleg->SetFillStyle(0);


	
	/*
	new TCanvas();
	h1clon->Draw();
	fourier->SetParameters(v0,v1,1);
	fourier->Draw("same");
	*/
      }

      
      

      //TString out_filename3 = input_dir+"plots/"mod"_"+var+"_mult.png";
      //Plotter(mult_histos[0],mult_histos[1],mult_histos[2],mult_histos[3],out_filename3);
      
    }
    
  }

  cout << "|";
  for ( auto mod: modes ) {
    cout << "|" << mod;
  }
  cout << "|" << endl;
  for ( int tar = 0; tar < 4; tar++ ) {
    cout << "|" << targets[tar];
    for ( int bin = 0; bin < 3; bin++ ) {
      cout << Form("| %.2e",v2table[tar][bin]);
    }
    cout << "|" << endl;
  }

  cout << "|";
  for ( auto mod: modes ) {
    cout << "|" << mod;
  }
  cout << "|" << endl;
  for ( int tar = 0; tar < 4; tar++ ) {
    cout << "|" << targets[tar];
    for ( int bin = 0; bin < 3; bin++ ) {
      cout << Form("| %.2e",v2tableerr[tar][bin]);
    }
    cout << "|" << endl;
  }

  TH1F htarv2z1;
  TH1F htarv2z2;
  TH1F htarv2z3;
  for ( int tar = 0; tar < 4; tar++ ) {
    htarv2z1.Fill(targets[tar],v2table[tar][0]);
    htarv2z1.SetBinError(tar+1,v2tableerr[tar][0]);

    htarv2z2.Fill(targets[tar],v2table[tar][1]);
    htarv2z2.SetBinError(tar+1,v2tableerr[tar][1]);

    htarv2z3.Fill(targets[tar],v2table[tar][2]);
    htarv2z3.SetBinError(tar+1,v2tableerr[tar][2]);
  }
  

  TH1F hbinv2d;
  TH1F hbinv2c;
  TH1F hbinv2fe;
  TH1F hbinv2pb;
  for ( int bin = 0; bin < 3; bin++ ) {
    TString binname;
    if ( bin == 0 ) {
      binname = "Z_{h}^{assoc.} < 0.12";
    } else if ( bin == 1 ) {
      binname = "0.12 < Z_{h}^{assoc.} < 0.25";
    } else if ( bin == 2 ) {
      binname = "0.25 < Z_{h}^{assoc.} < 0.4";
    }
    hbinv2d.Fill(binname,v2table[0][bin]);
    hbinv2d.SetBinError(bin+1,v2tableerr[0][bin]);
		 
    hbinv2c.Fill(binname,v2table[1][bin]);
    hbinv2c.SetBinError(bin+1,v2tableerr[1][bin]);
    
    hbinv2fe.Fill(binname,v2table[2][bin]);
    hbinv2fe.SetBinError(bin+1,v2tableerr[2][bin]);

    hbinv2pb.Fill(binname,v2table[3][bin]);
    hbinv2pb.SetBinError(bin+1,v2tableerr[3][bin]);
    cout << "ATENTO"  << modes[bin] << endl;
  }
  /*
  hbinv2d.Fill("munga",v2table[0][1]);
  hbinv2d.SetBinError(3+1,v2tableerr[0][1]);
		 
  hbinv2c.Fill("munga",v2table[1][1]);
  hbinv2c.SetBinError(3+1,v2tableerr[1][1]);
  
  hbinv2fe.Fill("munga",v2table[2][1]);
  hbinv2fe.SetBinError(3+1,v2tableerr[2][1]);
  
  hbinv2pb.Fill("munga",v2table[3][1]);
  hbinv2pb.SetBinError(3+1,v2tableerr[3][1]);
  */
  TCanvas * c1 = new TCanvas();
  c1->SetLeftMargin(0.13);
  double max;
  max = TMath::Max(htarv2z1.GetMaximum(),htarv2z2.GetMaximum());
  max = TMath::Max(max,htarv2z3.GetMaximum());


  htarv2z1.SetMarkerStyle(kFullTriangleUp);
  htarv2z1.SetMarkerColor(kRed);
  htarv2z2.SetMarkerStyle(kFullCircle);
  htarv2z2.SetMarkerColor(kGreen);
  htarv2z3.SetMarkerStyle(kFullSquare);
  htarv2z3.SetMarkerColor(kBlue);
  htarv2z1.SetTitle("Z_{h}^{assoc.} < 0.12;;v_{2}");
  htarv2z2.SetTitle("0.12 < Z_{h}^{assoc.} < 0.25;;v_{2}");
  htarv2z3.SetTitle("0.25 < Z_{h}^{assoc.} < 0.4;;v_{2}");
  
  htarv2z1.Draw("E0X0");
  htarv2z1.SetMaximum(max+0.1*max);
  htarv2z2.Draw("E0X0 same");

  htarv2z3.Draw("E0X0 same");


  gStyle->SetOptStat(0);
  TLegend *myleg;
  myleg = c1->BuildLegend(0.65,0.7,0.9,0.9,"Z_{h} bins","");
  myleg->SetFillStyle(0);

  c1->SaveAs(input_dir+"plots_fourier/v2_targ_bins.png");
  

  TCanvas * c2 = new TCanvas();
  c2->SetLeftMargin(0.13);
  
  max = 0;
  max = TMath::Max(hbinv2d.GetMaximum(),hbinv2c.GetMaximum());
  max = TMath::Max(max,hbinv2fe.GetMaximum());
  max = TMath::Max(max,hbinv2pb.GetMaximum());


  hbinv2d.SetMarkerStyle(kFullTriangleUp);
  hbinv2d.SetMarkerColor(kRed);
  hbinv2c.SetMarkerStyle(kFullCircle);
  hbinv2c.SetMarkerColor(kGreen);
  hbinv2fe.SetMarkerStyle(kFullSquare);
  hbinv2fe.SetMarkerColor(kBlue);
  hbinv2pb.SetMarkerStyle(kFullTriangleDown);
  hbinv2pb.SetMarkerColor(kMagenta);
  hbinv2d.SetTitle("D;;v_{2}");
  hbinv2c.SetTitle("C;;v_{2}");
  hbinv2fe.SetTitle("Fe;;v_{2}");
  hbinv2pb.SetTitle("Pb;;v_{2}");

  hbinv2d.Draw("E0X0");
  hbinv2d.SetMaximum(max+0.1*max);
  hbinv2d.SetMinimum(-0.004);
  hbinv2c.Draw("E0X0 same");
  hbinv2fe.Draw("E0X0 same");
  hbinv2pb.Draw("E0X0 same");


  gStyle->SetOptStat(0);
  TLegend *myleg2;
  myleg2 = c2->BuildLegend(0.65,0.7,0.9,0.9,"target","");
  myleg2->SetFillStyle(0);


  c2->SaveAs(input_dir+"plots_fourier/v2_zh_bins.png");
    
  cout <<  hbinv2d.GetNbinsX() << endl;
  
}
