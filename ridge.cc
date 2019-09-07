#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

TString out_path = "/home/luciano/Physics/CLAS/pion_ridge/";
bool m_debug = false;
bool m_simulation = false;
TDatabasePDG db;

void virtualframe(TLorentzVector * v1,
		  double angle1, double angle2, double boost)
{
  TLorentzVector *v2;
  v2 = v1;
  v2->RotateZ(-angle1);
  v2->RotateY(-angle2);
  v2->Boost(0,0,-boost);
}

Float_t PhiPQ(TVector3 pion, TLorentzVector gamma_vir )
{
  TVector3 Vhelp(0.,0.,1.0);
  TVector3 gamma = gamma_vir.Vect();
  Double_t phi_z = TMath::Pi()-gamma.Phi();
  pion.RotateZ(phi_z);
  gamma.RotateZ(phi_z);
  Double_t phi_y = gamma.Angle(Vhelp);
  pion.RotateY(phi_y);
  //gamma.RotateY(phi_y);
  return pion.Phi();
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

void export_hist(TH2 * h2, TString out_filename, TString options = "colz") {
  auto c1 = new TCanvas();
  c1->SetCanvasSize(800,800);
  //c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  h2->Draw(options);
  c1->SaveAs(out_filename);
  delete c1;
}

void ridge_plot(TH2 *h2, TString out_name)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",1000,1000);
  c1->SetLeftMargin(0.14);
  h2->Draw("surf1");
  
  // aesthetics
  gStyle->SetOptStat(0);
  double theta = -30;
  double phi   = 60;
  h2->GetXaxis()->CenterTitle(true);
  h2->GetXaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->CenterTitle(true);

  // main 2d plot 
  TString out_filename = out_path+out_name;
  c1->SaveAs(out_filename);

  // mountain side plots
  gPad->GetView()->RotateView(0,90);
  out_filename = out_path+"mountainsideY_"+out_name;
  c1->SaveAs(out_filename);
  gPad->GetView()->RotateView(-90,90);
  out_filename = out_path+"mountainsideX_"+out_name;
  c1->SaveAs(out_filename);

  // projection plots
  TH1 * projx = h2->ProjectionX();
  projx->SetFillColor(kBlue);
  projx->Draw();
  out_filename = out_path+"projx_"+out_name;
  c1->SaveAs(out_filename);
  
  TH1 * projy = h2->ProjectionY();
  projy->SetFillColor(kBlue);
  projy->Draw();
  out_filename = out_path+"projy_"+out_name;
  c1->SaveAs(out_filename);

  // comment next line if not running on batch mode I guess
  delete c1;
  delete projx;
  delete projy;
}

// plots two th2 side by side
void side_by_side(TH2 *hlab,TH2 *hvir, TString out_filename)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",2400,1000);
  TString corr;
  c1->SetFixedAspectRatio();
  c1->Divide(2,0,0.02,0.02);
  
  c1->cd(1);
  
  hlab->Draw("surf1");
  gStyle->SetOptStat(0);
  // TODO: maybe add text with entry number
  double theta = -30;
  double phi   = 60;
  hlab->GetXaxis()->CenterTitle(true);
  hlab->GetXaxis()->SetTitleOffset(1.5);
  hlab->GetYaxis()->CenterTitle(true);
  gPad->GetView()->RotateView(theta,phi);  
  
  c1->cd(2);
  hvir->Draw("surf1");
  gStyle->SetOptStat(0);
  hvir->GetXaxis()->CenterTitle(true);
  hvir->GetXaxis()->SetTitleOffset(1.5);
  hvir->GetYaxis()->CenterTitle(true);
  gPad->GetView()->RotateView(-30,60);  

  out_filename = out_path+out_filename;
  c1->SaveAs(out_filename);
  delete c1;
}

void ridge()
{
  if ( m_simulation == true ){
    out_path+="simulation/";
  } else {
    out_path += "data/";
  }
  
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
  //TTreeReaderArray<Float_t> tpq(th,"ThetaPQ"); // not on 11GeV ascii file
  //TTreeReaderArray<Float_t> ppq(th,"PhiPQ"); // not on 11GeV ascii file

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

  TH2 * ns = new TH2F("Signal Distribution, lab frame (theta)",
		      "N_{s}(#Delta#theta, #Delta#phi);#Delta#theta;#Delta#phi",
		      16,-4.0,4.0,
		      16,-4.0,4.0);

  TH2 * ns_eta = new TH2F("Signal Distribution, lab frame",
			  "N_{s}(#Delta#eta, #Delta#phi);#Delta#eta;#Delta#phi",
			  16,-4.0,4.0,
			  16,-4.0,4.0);
  
  TH2 * ns_vir = new TH2F("Signal Distribution, #gamma^{*} frame",
			  "N_{s}(#Delta#phi, #Delta#theta), #gamma^{*} frame;#Delta#theta;#Delta#phi",
			  16,-4.0,4.0,
			  16,-4.0,4.0);

  TH2 * ns_pqf = new TH2F("Signal Distribution, PQ frame",
			  "N_{s}(#Delta#theta_{PQ}, #Delta#phi_{PQ});#Delta#theta_{PQ};#Delta#phi_{PQ}",
			  16,-4.0,4.0,
			  16,-4.0,4.0);
  
  TH2 * reco = new TH2F("Reconstructed angles",
			"Reconstructed #pi^{#pm};#phi;#theta",
			72,-TMath::Pi(),TMath::Pi(),
			36,0.0,TMath::Pi());
  TH2 * reco_PQ = new TH2F("Reconstructed angles, PQ frame",
			   "Reconstructed #pi^{#pm};#phi_{PQ};#theta_{PQ}",
			   72,-TMath::Pi(),TMath::Pi(),
			   36,0.0,TMath::Pi());

  TH2 * reco_eta = new TH2F("Reconstructed pions, eta",
			    "Reconstructed #pi^{#pm};#phi;#eta",
			    72,-TMath::Pi(),TMath::Pi(),
			    80,-2.0,4.0);
			  
  
  
  
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
    
    bool fill_neg = true;
    for ( auto i: pos_pions ) {
      TVector3 piplus (px[i],py[i],pz[i]);
      Float_t phiplus      = piplus.Phi();
      Float_t thetaplus    = piplus.Theta();
      Float_t etaplus      = piplus.Eta();
      reco->Fill(phiplus,thetaplus);
      reco_eta->Fill(phiplus,etaplus);

      TVector3 gammav_rot  = virtual_photon.Vect();
      Float_t thetaPQplus  = gammav_rot.Angle(piplus);
      Float_t phiPQplus    = PhiPQ(piplus,virtual_photon);
      reco_PQ->Fill(phiPQplus,thetaPQplus);
      
      for ( auto j: neg_pions ) {
	TVector3 piminus (px[j],py[j],pz[j]);
	Float_t phiminus     = piminus.Phi();
	Float_t thetaminus   = piminus.Theta();
	Float_t etaminus     = piminus.Eta();
	
	TVector3 gammav_rot   = virtual_photon.Vect();
	Float_t thetaPQminus  = gammav_rot.Angle(piminus);
	Float_t phiPQminus    = PhiPQ(piminus,virtual_photon);
	
	ns->Fill(thetaplus-thetaminus,phiplus-phiminus);
	ns_eta->Fill(etaplus-etaminus,phiplus-phiminus);
	ns_pqf->Fill(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	if ( fill_neg ){ // fill "acceptance" histograms once only
	  reco->Fill(phiminus,thetaminus);
	  reco_eta->Fill(phiminus,etaminus);
	  reco_PQ->Fill(phiPQminus,thetaPQminus);
	}
	fill_neg = false;
	//ns_PQ->Fill(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
      }
    }

    // this analysis doesn't use the virtual photon frame
    
    /*
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
	ns_vir->Fill(thetaplus-thetaminus,phiplus-phiminus);
      }
    }
    
    TLorentzVector vsum_pions;
    for ( int j = 0; j < min(pos_pions.size(),neg_pions.size()); j++ ){
      vsum_pions = Particles[neg_pions[j]]+Particles[pos_pions[j]];
      h1->Fill(vsum_pions.M() );  
    }
    
    //cout << pos_pions.size() << " " << neg_pions.size() << endl;;
    counter++;
    */
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << counter << endl;
  
  //side_by_side(ns,ns_vir,"NS_lab_virtual.png");

  ridge_plot(ns,"single_ridge_theta.pdf");
  ridge_plot(ns_eta,"single_ridge_eta.pdf");
  ridge_plot(ns_pqf,"single_ridge_PQ.pdf");

  ridge_plot(ns,"single_ridge_theta.png");
  ridge_plot(ns_eta,"single_ridge_eta.png");
  ridge_plot(ns_pqf,"single_ridge_PQ.png");


  
  export_hist(reco,out_path+"reco.png");
  export_hist(reco_eta,out_path+"reco_eta.png");
  export_hist(reco_PQ,out_path+"reco_PQ.png");
  
  cout << "finished correctly" << endl;
}
