// Provides plotting capabilities for the 2pc (one dimensional) analysis

void Plotter(TH1F * C, TH1F * Fe, TH1F * Pb, TH1F * D, TString out_filename = "" )
{
  //gStyle->SetOptStat(0);
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
  D->GetYaxis()->SetTitle("N_{s}/N_{m}(#theta_{PQ})");

  D->SetTitle("D");
  C->SetTitle("C");
  Pb->SetTitle("Pb");
  Fe->SetTitle("Fe");

  TLegend *myleg;
  myleg = c1->BuildLegend(0.2,0.7,0.34,0.9,"Targets","");
  
  myleg->SetFillStyle(0);

  if ( out_filename == "" ) {
    c1->SaveAs("out1.png");
  } else {
    c1->SaveAs(out_filename);
  }
}
