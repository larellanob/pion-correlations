#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

TString out_path = "/home/luciano/Physics/CLAS/pion_correlation/";

void virtualframe(TLorentzVector * v1,
		  double angle1, double angle2, double boost)
{
  TLorentzVector *v2;
  v2 = v1;
  v2->RotateZ(-angle1);
  v2->RotateY(-angle2);
  v2->Boost(0,0,-boost);
}

void plot(TGraph *gr, TString title)
{
  TString corr;
  corr =  Form("%.4f",gr->GetCorrelationFactor());
  gr->SetTitle("Correlation: "+corr+" / "+title);
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(kBlue);
  new TCanvas();
  gr->Draw("AP");
  //cout <<  gr->GetCorrelationFactor() << endl;
}

void side_by_side(TH2 *hlab,TH2 *hvir, TString out_filename)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",1200,500);
  TString corr;
  c1->SetFixedAspectRatio();
  c1->Divide(2,0,0.02,0.02);
  
  c1->cd(1);
  hlab->Draw("colz");
  gStyle->SetOptStat(0);
  TPaveText pt(0.8,1.5,1.7,1.8,"TR");
  corr = Form("%.4f",hlab->GetCorrelationFactor());
  corr = "Correlation = "+corr;
  pt.AddText(corr);
  pt.Draw();
  
  c1->cd(2);
  hvir->Draw("colz");
  corr = Form("%.4f",hvir->GetCorrelationFactor());
  corr = "Correlation = "+corr;
  pt.AddText(corr);
  pt.Draw();
  
  out_filename = out_path+out_filename;
  c1->SaveAs(out_filename);
  //delete c1;
}

void pion_correlation()
{
  TDatabasePDG db;
  bool m_debug = false;
  bool m_simulation = true;

  TFile *f;
  Double_t kEbeam;
  if ( m_simulation ){
    f = new TFile("/home/luciano/Physics/CLAS/data/tree_output.root");
    kEbeam = 11.0;
  } else {
    f = new TFile("/home/luciano/Physics/CLAS/data/CFFTree_data.root");
    kEbeam = 5.014;
  }
  const Double_t kMpi   = 0.13957;
  const Double_t kMprt  = 0.938272;
  float rhomass = db.GetParticle(113)->Mass();
  
  TTreeReader th("tree_data",f);
  TTreeReader e("e_rec",f);
    
  TTreeReaderArray<Float_t> pid(th,"pid");
  TTreeReaderArray<Float_t> q2(th,"Q2");
  TTreeReaderArray<Float_t> px(th,"Px");
  TTreeReaderArray<Float_t> py(th,"Py");
  TTreeReaderArray<Float_t> pz(th,"Pz");
  TTreeReaderArray<Float_t> E(th,"E"); // this isnt on CFF

  TTreeReaderValue<Float_t> epx(e,"Pex");
  TTreeReaderValue<Float_t> epy(e,"Pey");
  TTreeReaderValue<Float_t> epz(e,"Pez");
  TTreeReaderValue<Float_t> eq2(e,"Q2");
  TTreeReaderValue<Float_t> ew(e,"W");
    
  cout << "working" << endl;

  TH1 * h1 = new TH1F("h1",
		      "Inv. mass of (#pi^{+} #pi^{-});M_{sum};counts",
		      50, 0, 1.5);
  std::vector<int> pos_pions;
  std::vector<int> neg_pions;
  std::vector<TLorentzVector> Particles;
  
  UInt_t events = 0;
  int counter=0;

  //px
  std::vector<Float_t> x_lab;
  std::vector<Float_t> y_lab;
  std::vector<Float_t> x_vir;
  std::vector<Float_t> y_vir;

  //py
  std::vector<Float_t> x_labpy;
  std::vector<Float_t> y_labpy;
  std::vector<Float_t> x_virpy;
  std::vector<Float_t> y_virpy;

  TH2 * hpy_lab = new TH2F("hpy_lab",
			   "Lab frame;p_{y}(#pi^{+});p_{y}(#pi^{-})",
			   30,-1.5,1.5, 30,-1.5,1.5);
  
  TH2 * hpy_vir = new TH2F("hpy_vir",
		       "#gamma^{*}-P frame;p_{y}(#pi^{+});p_{y}(#pi^{-})",
		       30,-1.5,1.5, 30,-1.5,1.5);

  //pz
  std::vector<Float_t> x_labpz;
  std::vector<Float_t> y_labpz;
  std::vector<Float_t> x_virpz;
  std::vector<Float_t> y_virpz;
  
  //P
  std::vector<Float_t> x_labp;
  std::vector<Float_t> y_labp;
  std::vector<Float_t> x_virp;
  std::vector<Float_t> y_virp;

  //E
  std::vector<Float_t> x_labe;
  std::vector<Float_t> y_labe;
  std::vector<Float_t> x_vire;
  std::vector<Float_t> y_vire;

  //angle
  std::vector<Float_t> x_labang;
  std::vector<Float_t> y_labang;
  std::vector<Float_t> x_virang;
  std::vector<Float_t> y_virang;
  
  while ( th.Next() ){
    e.Next();
    events++;
    
    if ( events % 1000000 == 0 )
      cout << events << " events" << endl;


    if ( *ew < 2.0 || *eq2 < 1.0 ) // DIS cut
      continue;

    // Electron
    
    pos_pions.clear();
    neg_pions.clear();
    for ( UInt_t i = 0; i < pid.GetSize(); i++ ){
      if ( pid[i] == 211 ) {
	pos_pions.push_back(i);
      }
      if ( pid[i] == -211 ){
	neg_pions.push_back(i);
      }
    }

    
    if ( pos_pions.size() != neg_pions.size() )
      continue;
    
    if ( pos_pions.size() == 0 )
      continue;
    
    for ( int i = 0; i < min(pos_pions.size(),neg_pions.size()); i++ ){
      x_lab.push_back(px[pos_pions[i]]);
      y_lab.push_back(px[neg_pions[i]]);
      x_labpy.push_back(py[pos_pions[i]]);
      y_labpy.push_back(py[neg_pions[i]]);
      hpy_lab->Fill(py[pos_pions[i]],py[neg_pions[i]]);
      x_labpz.push_back(pz[pos_pions[i]]);
      y_labpz.push_back(pz[neg_pions[i]]);
      x_labe.push_back(E[pos_pions[i]]);
      y_labe.push_back(E[neg_pions[i]]);
      double phi = 0;
      phi = atan2(py[pos_pions[i]],px[pos_pions[i]]);
      x_labang.push_back(phi);
      phi = atan2(py[neg_pions[i]],px[neg_pions[i]]);
      y_labang.push_back(phi);
    }
    
    // virtual photon frame is per event
    Double_t scattered_energy = sqrt((*epx)*(*epx)+(*epy)*(*epy)+(*epz)*(*epz));
    if ( scattered_energy < 0.005 )
      cout << "WARNING: Low e' energy " << scattered_energy << endl; 
    if ( scattered_energy > kEbeam  )
      cout << "WARNING: High e' energy = " <<  scattered_energy << endl;
    TLorentzVector virtual_photon
      (-(*epx),-(*epy),kEbeam-(*epz),kEbeam-scattered_energy);
    TLorentzVector initial_photon = virtual_photon;
    Double_t Egamma = virtual_photon.E();
    Double_t angle1 = TMath::Pi() + atan2((*epy),(*epx));
    virtual_photon.RotateZ(-angle1);


    double numer = kEbeam-(*epz);
    double denom = virtual_photon.P();

    Double_t angle2;// = acos(numer/denom);
    
    if ( numer/denom < TMath::Pi() ){
      angle2 = acos(numer/denom);
    } else if ( numer/denom >= TMath::Pi() ){
      angle2 = 2*TMath::Pi()*acos(numer/denom);
    }
    
    virtual_photon.RotateY(-angle2);
    Double_t boost  = virtual_photon.P()/(Egamma+kMprt);
    virtual_photon.Boost(0,0,-boost);

 
    // lorentz vectors
    Particles.clear();
    for ( int i = 0; i < pid.GetSize(); i++ ){
      float mass = db.GetParticle(pid[i])->Mass();
      TLorentzVector v (px[i],py[i],pz[i],E[i]);
      Particles.push_back(v);
      virtualframe(&Particles[i],angle1,angle2,boost);	
    }

    for ( int i = 0; i < min(pos_pions.size(),neg_pions.size()); i++ ){
      x_vir.push_back(Particles[pos_pions[i]].Px() );
      y_vir.push_back(Particles[neg_pions[i]].Px() );
      x_virpy.push_back(Particles[pos_pions[i]].Py() );
      y_virpy.push_back(Particles[neg_pions[i]].Py() );
      hpy_vir->Fill(Particles[pos_pions[i]].Py(),Particles[neg_pions[i]].Py());
      x_virpz.push_back(Particles[pos_pions[i]].Pz() );
      y_virpz.push_back(Particles[neg_pions[i]].Pz() );
      x_vire.push_back(Particles[pos_pions[i]].E() );
      y_vire.push_back(Particles[neg_pions[i]].E() );
      double phi = 0;
      phi = atan2(Particles[pos_pions[i]].Py(),Particles[pos_pions[i]].Px());
      x_virang.push_back(phi);
      phi = atan2(Particles[neg_pions[i]].Py(),Particles[neg_pions[i]].Px());
      y_virang.push_back(phi);
    }
    
    TLorentzVector vsum_pions;
    for ( int j = 0; j < min(pos_pions.size(),neg_pions.size()); j++ ){
      vsum_pions = Particles[neg_pions[j]]+Particles[pos_pions[j]];
      h1->Fill(vsum_pions.M() );  
    }
    
    //cout << pos_pions.size() << " " << neg_pions.size() << endl;;
    counter++;
    
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << counter << endl;

  cout << "xlab, ylab sizes " << x_lab.size() << " " << y_lab.size() << endl;
  cout << "xvir, yvir sizes " << x_vir.size() << " " << y_vir.size() << endl;
  
  TGraph *grlab = new TGraph(x_lab.size(),&x_lab[0],&y_lab[0]);
  TGraph *grvir = new TGraph(x_vir.size(),&x_vir[0],&y_vir[0]);
  TGraph *grlabpy = new TGraph(x_labpy.size(),&x_labpy[0],&y_labpy[0]);
  TGraph *grvirpy = new TGraph(x_virpy.size(),&x_virpy[0],&y_virpy[0]);
  TGraph *grlabpz = new TGraph(x_labpz.size(),&x_labpz[0],&y_labpz[0]);
  TGraph *grvirpz = new TGraph(x_virpz.size(),&x_virpz[0],&y_virpz[0]);
  TGraph *grlabe = new TGraph(x_labe.size(),&x_labe[0],&y_labe[0]);
  TGraph *grvire = new TGraph(x_vire.size(),&x_vire[0],&y_vire[0]);
  TGraph *grlabang = new TGraph(x_labang.size(),&x_labang[0],&y_labang[0]);
  TGraph *grvirang = new TGraph(x_virang.size(),&x_virang[0],&y_virang[0]);


  plot(grlab,"Lab frame;#pi^{+} P_x;#pi^{-} P_x");
  plot(grvir,"Virtual photon - Proton frame;#pi^{+} P_x;#pi^{-} P_x");
  plot(grlabpy,"Lab frame;#pi^{+} P_y;#pi^{-} P_y");
  plot(grvirpy,"Virtual photon - Proton frame;#pi^{+} P_y;#pi^{-} P_y");
  plot(grlabpz,"Lab frame;#pi^{+} P_z;#pi^{-} P_z");
  plot(grvirpz,"Virtual photon - Proton frame;#pi^{+} P_z;#pi^{-} P_z");
  plot(grlabe,"Lab frame;#pi^{+} E;#pi^{-} E");
  plot(grvire,"Virtual photon - Proton frame;#pi^{+} E;#pi^{-} E");
  plot(grlabe,"Lab frame;#pi^{+} #phi;#pi^{-} #phi");
  plot(grvire,"Virtual photon - Proton frame;#pi^{+} #phi;#pi^{-} #phi");
  new TCanvas();
  h1->Draw();
  new TCanvas();
  hpy_lab->Draw("colz");
  new TCanvas();
  hpy_vir->Draw("colz");

  side_by_side(hpy_lab,hpy_vir,"Py.png");
  
  cout << "finished correctly" << endl;
}
