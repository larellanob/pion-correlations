// Provides plotting capabilities for the 2pc (one dimensional) analysis

void RatioPlotter(TH1F * C, TH1F * Fe, TH1F * Pb, TString out_filename, bool m_simulation)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas * c1 = new TCanvas();
  c1->SetCanvasSize(700,700);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.12);

  // maybe Ill have to add setfillstyle 3353 here?
  //D->SetMarkerColor(kViolet);
  //D->SetMarkerStyle(kFullTriangleDown);
  //D->SetMarkerSize(1.5);

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

  //D->Draw("E0X0");
  C->Draw("E0X0");
  Fe->Draw("E0X0 same");
  Pb->Draw("E0X0 same");

  // shows entries
  //C->SetTitle(TString::Format("C  %.0f",C->GetEntries()));
  //Pb->SetTitle(TString::Format("Pb  %.0f",Pb->GetEntries()));
  //Fe->SetTitle(TString::Format("Fe  %.0f",Fe->GetEntries()));

  // doesnt show entries
  C->SetTitle("C/D");
  Pb->SetTitle("Pb/D");
  Fe->SetTitle("Fe/D");
  
  C->SetYTitle("N_{c}^{A}/N_{c}^{D}");
  TLegend *myleg;
  if ( m_simulation ) {
    myleg = c1->BuildLegend(0.2,0.7,0.44,0.9,"Hayk's simulation","");
  } else
    myleg = c1->BuildLegend(0.2,0.7,0.44,0.9,"eg2 data         ","");
  
  myleg->SetFillStyle(0);

  int nbins = C->GetNbinsX();
  TLine *t = new TLine(C->GetBinLowEdge(1),1,C->GetBinLowEdge(nbins+1),1);
  t->SetLineStyle(9);
  t->Draw();


  if ( out_filename == "" ) {
    c1->SaveAs("out1.png");
  } else {
    c1->SaveAs(out_filename);
  }
}
