// Provides plotting capabilities for the 2pc (one dimensional) analysis
void LatexText(Double_t x, Double_t y, int font, TString text);
void Plotter(TH1F * C, TH1F * Fe, TH1F * Pb, TH1F * D, TString out_filename = "" , TString gmode = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas * c1 = new TCanvas();
  c1->SetCanvasSize(700,700);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.12);

  // maybe Ill have to add setfillstyle 3353 here?
  D->SetMarkerColor(kViolet);
  D->SetMarkerStyle(kFullTriangleDown);
  D->SetMarkerSize(1.5);

  C->SetFillStyle(3353);
  C->SetMarkerStyle(kFullCircle);
  C->SetMarkerColor(kBlue);
  C->SetMarkerSize(1.5);

  Fe->SetMarkerColor(kRed);
  //Fe->SetFillColor(kRed);
  // Fe->SetFillStyle(3335);
  Fe->SetMarkerStyle(kFullSquare);
  Fe->SetMarkerSize(1.5);

  Pb->SetMarkerColor(kGreen);
  Pb->SetMarkerStyle(kFullTriangleUp);
  Pb->SetMarkerSize(1.5);

  D->Draw("E0X0");
  C->Draw("E0X0 same");
  Fe->Draw("E0X0 same");
  Pb->Draw("E0X0 same");
  //D->GetYaxis()->SetTitle("N_{s}/N_{m}(#theta_{PQ})");


  TLegend *myleg;
  // targets + entries
  /*  
  D->SetTitle(TString::Format("D  %.0f",D->GetEntries()));
  C->SetTitle(TString::Format("C  %.0f",C->GetEntries()));
  Pb->SetTitle(TString::Format("Pb  %.0f",Pb->GetEntries()));
  Fe->SetTitle(TString::Format("Fe  %.0f",Fe->GetEntries()));
  //myleg = c1->BuildLegend(0.2,0.7,0.44,0.9,"Targets / Entries","");
  */


  // targets without entries
  D->SetTitle("D");
  C->SetTitle("C");
  Pb->SetTitle("Pb");
  Fe->SetTitle("Fe");
  myleg = c1->BuildLegend(0.2,0.7,0.35,0.9,"Target","");
  
  myleg->SetFillStyle(0);

    TString zhrange;
  if  ( gmode  == "zhlead" ) {
    zhrange = "z_{h}^{assoc.} < 0.5";
  } else if  ( gmode  == "zh3" || gmode == "zh3ycut" ) {
    zhrange = "0.25 < z_{h}^{assoc.} < 0.4";
  } else if  ( gmode  == "zh2" || gmode == "zh2ycut" ) {
    zhrange = "0.12 < z_{h}^{assoc.} < 0.25";
  } else if  ( gmode  == "zh1" || gmode == "zh1ycut" ) {
    zhrange = "z_{h}^{assoc.} < 0.12";
  }

  if  ( gmode  == "ptlead" ) {
    zhrange = "p_{T}^{assoc.} < 0.5";
  } else if  ( gmode  == "pt3" || gmode == "pt3ycut" ) {
    zhrange = "0.3 < p_{T}^{assoc.} < 0.4";
  } else if  ( gmode  == "pt2" || gmode == "pt2ycut" ) {
    zhrange = "0.2 < p_{T}^{assoc.} < 0.3";
  } else if  ( gmode  == "pt1" || gmode == "pt1ycut" ) {
    zhrange = "p_{T}^{assoc.} < 0.2";
  }

  //  zhrange = "";
  LatexText(0.40,0.85,42,zhrange);    
  if ( gmode == "zh1ycut" || gmode == "zh2ycut" || gmode == "zh3ycut" ) {
     LatexText(0.40,0.79,42,"|#Deltay-y_{bias}| > 1.0");    
  }
  if ( gmode == "pt1ycut" || gmode == "pt2ycut" || gmode == "pt3ycut" ) {
    LatexText(0.40,0.79,42,"|#Deltay| > 1.0");     
  }

  
  if ( out_filename == "" ) {
    c1->SaveAs("out1.png");
  } else {
    c1->SaveAs(out_filename);
  }
}

void LatexText(Double_t x, Double_t y, int font, TString text)
{
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(font);
  l2.DrawLatex(x,y,text);
}



void ridge_plot(TH2F h2, TString pre, TString post,TString target, TString gmode, bool m_simul = false)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",1000,1000);
  //c1->SetLeftMargin(0.14);
  c1->SetLeftMargin(0.18);
  //gStyle->SetOptStat(0);
  h2.Draw("surf1");

  // aesthetics

  //soumya angles
  //double theta = -50;
  //double phi   = 70;

  // h1 angles
  double theta = 140;
  double phi   = 40;

  //gStyle->SetPadRightMargin(0.5);
  TString axisx, axisy;
  h2.SetTitle("");

  TString histname, zaxis;
  histname = (TString) h2.GetName();
  if ( histname.Contains("n_c") ) {
    zaxis = "#frac{N_{s}}{N_{m}}";
  } else if ( histname.Contains("n_m") ) {
    zaxis = "N_{m}";
  } else if ( histname.Contains("n_s") ) {
    zaxis = "N_{s}";
  }

  if ( !m_simul ) {
    LatexText(0.03,0.95,62,"EG2 data");
  } else if ( m_simul ) {
    LatexText(0.03,0.95,62,"EG2 simulation");
  }
  LatexText(0.03,0.90,42,target+" target");

  TString zhrange;
  if  ( gmode  == "zhlead" ) {
    zhrange = "z_{h}^{assoc.} < 0.5";
  } else if  ( gmode  == "zh3" || gmode == "zh3ycut" ) {
    zhrange = "0.25 < z_{h}^{assoc.} < 0.4";
  } else if  ( gmode  == "zh2" || gmode == "zh2ycut" ) {
    zhrange = "0.12 < z_{h}^{assoc.} < 0.25";
  } else if  ( gmode  == "zh1" || gmode == "zh1ycut" ) {
    zhrange = "z_{h}^{assoc.} < 0.12";
  }
  LatexText(0.03,0.85,42,zhrange);
  //LatexText(0.03,0.85,42,zhrange);

  if ( target == "D" && gmode == "zh1ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 1.578");
  } else if ( target == "D" && gmode == "zh2ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.9759");
  } else if ( target == "D" && gmode == "zh3ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.5042");
  } else if ( target == "C" && gmode == "zh1ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 1.606");
  } else if ( target == "C" && gmode == "zh2ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 1.003");
  } else if ( target == "C" && gmode == "zh3ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.4939");
  } else if ( target == "Fe" && gmode == "zh1ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 1.510");
  } else if ( target == "Fe" && gmode == "zh2ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.969");
  } else if ( target == "Fe" && gmode == "zh3ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.4806");
  } else if ( target == "Pb" && gmode == "zh1ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 1.496");
  } else if ( target == "Pb" && gmode == "zh2ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.9623");
  } else if ( target == "Pb" && gmode == "zh3ycut") {
    LatexText(0.03,0.80,42,"y_{bias} = 0.5376");
  } 


  
  
  /*
  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(42);
  l2.DrawLatex(0.01,0.9,"EG2, "+target+" target");
  */
  
  axisx = h2.GetXaxis()->GetTitle();
  axisy = h2.GetYaxis()->GetTitle();
  h2.GetZaxis()->SetTitleOffset(2.2);
  h2.GetZaxis()->SetTitle(zaxis);
  h2.GetZaxis()->CenterTitle(true);
  h2.GetXaxis()->CenterTitle(true);
  h2.GetXaxis()->SetTitleOffset(2.2);
  if ( axisx == "#Delta#phiPQboosted2" ) {
    h2.GetXaxis()->SetTitle("#Delta#phi_{CoM} [deg]");
  }
  h2.GetYaxis()->CenterTitle(true);
  h2.GetYaxis()->SetTitleOffset(2.0);
  if ( axisy == "#Delta#y" ) {
    h2.GetYaxis()->SetTitle("#Delta y_{CoM}");
    if ( gmode == "zh1ycut" || gmode == "zh2ycut" || gmode == "zh3ycut" ) {
      h2.GetYaxis()->SetTitle("#Delta y_{CoM} - y_{bias}");
    }
  }
 
  // main 2d plot 
  c1->SaveAs(pre+post);
  gPad->GetView()->RotateView(theta,phi);
  c1->SaveAs(pre+post);
  /*
  // mountain side plots
  gPad->GetView()->RotateView(0,90);
  c1->SaveAs(pre+"_hillsideY_"+post);
  gPad->GetView()->RotateView(-90,90);
  c1->SaveAs(pre+"_hillsideX_"+post);

  // projection plots
  TH1 * projx = h2.ProjectionX();
  projx->SetFillColor(kBlue);
  projx->Draw();
  c1->SaveAs(pre+"_projx_"+post);
  
  TH1 * projy = h2.ProjectionY();
  projy->SetFillColor(kBlue);
  projy->Draw();
  c1->SaveAs(pre+"_projy_"+post);
  */

  // comment next lines if not running on batch mode I guess
  delete c1;
  //delete projx;
  //delete projy;
}
