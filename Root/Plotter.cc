// Provides plotting capabilities for the 2pc (one dimensional) analysis

void Plotter(TH1F * C, TH1F * Fe, TH1F * Pb, TH1F * D, TString out_filename = "" )
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


  //TString test;
  //test = Form("D ",D->GetEntries());
  D->SetTitle(TString::Format("D  %.0f",D->GetEntries()));
  C->SetTitle(TString::Format("C  %.0f",C->GetEntries()));
  Pb->SetTitle(TString::Format("Pb  %.0f",Pb->GetEntries()));
  Fe->SetTitle(TString::Format("Fe  %.0f",Fe->GetEntries()));

  TLegend *myleg;
  myleg = c1->BuildLegend(0.2,0.7,0.44,0.9,"Targets / Entries","");
  
  myleg->SetFillStyle(0);

  if ( out_filename == "" ) {
    c1->SaveAs("out1.png");
  } else {
    c1->SaveAs(out_filename);
  }
}


void ridge_plot(TH2F h2, TString pre, TString post)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",1000,1000);
  c1->SetLeftMargin(0.14);
  //gStyle->SetOptStat(0);
  h2.Draw("surf1");

  // aesthetics

  //soumya angles
  //double theta = -50;
  //double phi   = 70;

  // h1 angles
  double theta = 140;
  double phi   = 40;

  
  TString axisx, axisy;
  axisx = h2.GetXaxis()->GetTitle();
  axisy = h2.GetYaxis()->GetTitle();
  h2.GetXaxis()->CenterTitle(true);
  h2.GetXaxis()->SetTitleOffset(2.2);
  if ( axisx == "#Delta#phiPQboosted2" ) {
    h2.GetXaxis()->SetTitle("#Delta#phi_{CoM} [deg]");
  }
  h2.GetYaxis()->CenterTitle(true);
  h2.GetYaxis()->SetTitleOffset(1.8);
  if ( axisy == "#Delta#y" ) {
    h2.GetYaxis()->SetTitle("#Delta y_{CoM}");
  }
 
  // main 2d plot 
  c1->SaveAs(pre+post);
  gPad->GetView()->RotateView(theta,phi);
  c1->SaveAs(pre+post);
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


  // comment next lines if not running on batch mode I guess
  delete c1;
  delete projx;
  delete projy;
}


