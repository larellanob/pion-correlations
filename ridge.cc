#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

TString out_path = "/home/luciano/Physics/CLAS/pion_ridge/";

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
  
  hlab->Draw("surf1");
  gStyle->SetOptStat(0);
  /*
  TPaveText pt(0.8,1.5,1.7,1.8,"TR");
  corr = Form("%.4f",hlab->GetCorrelationFactor());
  corr = "Correlation = "+corr;
  pt.AddText(corr);  
  pt.Draw();*/
  double theta = -30;
  double phi   = 60;
  hlab->GetXaxis()->CenterTitle(true);
  hlab->GetXaxis()->SetTitleOffset(1.5);
  hlab->GetYaxis()->CenterTitle(true);
  gPad->GetView()->RotateView(theta,phi);  
  //gPad->Modified();
  //gPad->Update();
  

  c1->cd(2);
  hvir->Draw("surf1");
  gStyle->SetOptStat(0);
  hvir->GetXaxis()->CenterTitle(true);
  hvir->GetXaxis()->SetTitleOffset(1.5);
  hvir->GetYaxis()->CenterTitle(true);
  gPad->GetView()->RotateView(-30,60);  
  /*
  corr = Form("%.4f",hvir->GetCorrelationFactor());
  corr = "Correlation = "+corr;
  pt.AddText(corr);
  pt.Draw();*/

  
  out_filename = out_path+out_filename;
  c1->SaveAs(out_filename);
  //delete c1;
}

void ridge()
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
  //TTreeReaderArray<Float_t> E(th,"E"); // this isnt on CFF

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

  TH2 * hpy_lab = new TH2F("hpy_lab",
			   "Lab frame;p_{y}(#pi^{+});p_{y}(#pi^{-})",
			   30,-1.5,1.5, 30,-1.5,1.5);
  
  TH2 * hpy_vir = new TH2F("hpy_vir",
		       "#gamma^{*}-P frame;p_{y}(#pi^{+});p_{y}(#pi^{-})",
		       30,-1.5,1.5, 30,-1.5,1.5);

  TH2 * ns = new TH2F("Signal Distribution, lab frame",
		      "N_{s}(#Delta#phi,#Delta#eta);#Delta#phi;#Delta#theta",
		      16,-4.0,4.0,
		      16,-4.0,4.0);
  
  TH2 * ns_vir = new TH2F("Signal Distribution, #gamma^{*} frame",
			  "N_{s}(#Delta#phi,#Delta#eta), #gamma^{*} frame;#Delta#phi;#Delta#theta",
			  16,-4.0,4.0,
			  16,-4.0,4.0);
  
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

    /*
    if ( pos_pions.size() != neg_pions.size() )
      continue;
    
    */
    if ( pos_pions.size() == 0  || neg_pions.size() == 0 )
      continue;
        
    //for ( int i = 0; i < pos_pions.size(); i++ ) {
    for ( auto i: pos_pions ) {
      TVector3 piplus (px[i],py[i],pz[i]);
      Float_t phiplus    = piplus.Phi();
      Float_t thetaplus  = piplus.Theta();
      //for ( int j = 0; j < neg_pions.size(); j++ ) {
      for ( auto j: neg_pions ) {
	TVector3 piminus (px[j],py[j],pz[j]);
	Float_t phiminus    = piminus.Phi();
	Float_t thetaminus  = piminus.Theta();
	ns->Fill(phiplus-phiminus,thetaplus-thetaminus);
	if ( events < 1000 )
	  cout << phiplus << " " << phiminus << endl;
      }
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
      double E = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]-mass*mass);
      TLorentzVector v (px[i],py[i],pz[i],E);
      Particles.push_back(v);
      virtualframe(&Particles[i],angle1,angle2,boost);	
    }


    /// same as before but on gamma* frame
    for ( int i = 0; i < pos_pions.size(); i++ ) {
      TVector3 piplus ( Particles[pos_pions[i]].Px(),
		        Particles[pos_pions[i]].Py(),
		        Particles[pos_pions[i]].Pz());
      Float_t phiplus    = piplus.Phi();
      Float_t thetaplus  = piplus.Theta();
      for ( int j = 0; j < neg_pions.size(); j++ ) {
	TVector3 piminus ( Particles[neg_pions[j]].Px(),
			   Particles[neg_pions[j]].Py(),
			   Particles[neg_pions[j]].Pz());
	Float_t phiminus    = piminus.Phi();
	Float_t thetaminus  = piminus.Theta();
	ns_vir->Fill(phiplus-phiminus,thetaplus-thetaminus);
      }
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
  
  side_by_side(ns,ns_vir,"NS_lab_virtual.png");
  cout << "finished correctly" << endl;
}
