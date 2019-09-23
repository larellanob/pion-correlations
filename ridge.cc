#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TView.h"
#include "TStyle.h"

TString out_path = "/home/luciano/Physics/CLAS/pion_ridge/";
bool m_debug      = false;
bool m_simulation = false;
bool old_plus     = false; // if false, use old_minus
bool gExport_pdf  = false;
bool gExport_png  = true;
TDatabasePDG db;

Float_t PhiPQ(TVector3, TLorentzVector);
void export_hist(TH2 *, TString, TString options = "colz");
void ridge_plot(TH2 *, TString);
void side_by_side(TH2 *,TH2 *, TString);
void initialize_histograms(int nsbinsx, int nsbinsy, double nsedgex, double nsedgey);

// Histograms that will be accompanying us tonight:
TH2 * ns, * ns_eta, * ns_pqf;
TH2 * reco, * reco_plus, * reco_minus, * reco_PQ, * reco_eta;
TH2 * nm, *nm_pqf;
// if you need more, add them to function 'initialize_histograms' at
// the bottom
  

// MAIN
void ridge()
{
  if ( m_simulation == true ) {
    out_path+="simulation/";
  } else {
    out_path += "data/";
  }
  
  TFile *f;
  Double_t kEbeam;
  if ( m_simulation ) {
    f = new TFile("/home/luciano/Physics/CLAS/data/tree_output.root");
    kEbeam = 11.0;
  } else {
    f = new TFile("/home/luciano/Physics/CLAS/data/CFFTree_data.root");
    kEbeam = 5.014;
  }
  const Double_t kMpi   = 0.13957;
  const Double_t kMprt  = 0.938272;
  float rhomass = db.GetParticle(113)->Mass();


  // tree reader for thrown particles and scattered electron
  TTreeReader th("tree_data",f);
  TTreeReader e("e_rec",f);
    
  TTreeReaderArray<Float_t> pid(th,"pid");
  TTreeReaderArray<Float_t> q2(th,"Q2");
  TTreeReaderArray<Float_t> px(th,"Px");
  TTreeReaderArray<Float_t> py(th,"Py");
  TTreeReaderArray<Float_t> pz(th,"Pz");

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

  std::vector<TVector3> old_piplus;
  std::vector<TVector3> old_piminus;
  
  UInt_t events = 0;
  int counter=0;


  // histograms
  int nsbinsx = 32;
  int nsbinsy = nsbinsx*2.5;
  double nsedgex = TMath::Pi()*1.0;
  double nsedgey = TMath::Pi()*2.5;
  
  initialize_histograms(nsbinsx, nsbinsy, nsedgex, nsedgey);

  TH2 * pippim = new TH2F("pippim",
			  "#pi^{+}#pi^{-} same event;N_{#pi^{+}};N_{#pi^{-}}",
			  8,-0.5,7.5,8,-0.5,7.5);
  TH2 * pippim_ppp = new TH2F("pippim_ppp",
			      "#pi^{+}#pi^{-} past #pi^{+};N_{#pi^{+}};N_{#pi^{-}}",
			      8,-0.5,7.5,8,-0.5,7.5);
  TH2 * pippim_ppm = new TH2F("pippim_ppm",
			      "#pi^{+}#pi^{-} past #pi^{-};N_{#pi^{+}};N_{#pi^{-}}",
			      8,-0.5,7.5,8,-0.5,7.5);
  
  if ( m_debug )
    cout << "DEBUG: entering main tree reading loop" << endl;
  // main file/tree reading loop
  TLorentzVector virtual_photon;
  TLorentzVector oldgamma;
  while ( th.Next() ) {
    e.Next();
    events++;
    
    Double_t scattered_energy = sqrt((*epx)*(*epx)+(*epy)*(*epy)+(*epz)*(*epz));
    if ( scattered_energy < 0.005 )
      cout << "WARNING: Low e' energy " << scattered_energy << endl; 
    if ( scattered_energy > kEbeam  )
      cout << "WARNING: High e' energy = " <<  scattered_energy << endl;


    // virtual photon frame is per event
    oldgamma = virtual_photon; // true old gamma

    virtual_photon = {-(*epx),-(*epy),kEbeam-(*epz),kEbeam-scattered_energy};
    //oldgamma = virtual_photon; // actually the same old gamma
    
    if ( events % 1000000 == 0 )
      cout << events << " events" << endl;

    if ( *ew < 2.0 || *eq2 < 1.0 ) // DIS cut
      continue;

    // Electron
    
    pos_pions.clear();
    neg_pions.clear();
    for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
      if ( pid[i] == 211 ) {
	pos_pions.push_back(i);
      }
      if ( pid[i] == -211 ) {
	neg_pions.push_back(i);
      }
    }

    /*
    if ( pos_pions.size() != neg_pions.size() )
      continue;
    */
    
    if ( pos_pions.size() == 0  || neg_pions.size() == 0 )
      continue;

    if ( m_debug )
      cout << "DEBUG: Entering mutli event filling" << endl;
    // MULTI EVENT FILLING

    // note: pos_pions and neg_pions store the integer index of the
    // pions position within the event, whereas old_piminus and
    // old_piplus store the whole TVector3 pion data. The global
    // variable old_plus defines wheter the old pions + or - are used
    
    // OLD PI MINUS
    if ( old_piminus.size() != 0 && !old_plus && events > 1 ) {
      for ( auto i: pos_pions ) {
	TVector3 piplus (px[i],py[i],pz[i]);
	Float_t phiplus      = piplus.Phi();
	Float_t thetaplus    = piplus.Theta();
	
	TVector3 gammav_rot  = virtual_photon.Vect();
	Float_t thetaPQplus  = gammav_rot.Angle(piplus);
	Float_t phiPQplus    = PhiPQ(piplus,virtual_photon);
	//Float_t etaplus      = piplus.Eta(); // just implement theta for now
	for ( auto piminus: old_piminus ) {
	  Float_t phiminus     = piminus.Phi();
	  Float_t thetaminus   = piminus.Theta();

	  TVector3 gammav_rot  = oldgamma.Vect();
	  Float_t thetaPQminus  = gammav_rot.Angle(piminus);
	  Float_t phiPQminus    = PhiPQ(piminus,oldgamma);
	  nm->Fill(thetaplus-thetaminus,phiplus-phiminus);
	  nm_pqf->Fill(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	}
      }
    }

    
    // OLD PI PLUS
    if ( old_piplus.size() != 0 && old_plus && events > 1 ) {
      for ( auto piplus: old_piplus ) {
	Float_t phiplus      = piplus.Phi();
	Float_t thetaplus    = piplus.Theta();
	
	TVector3 gammav_rot  = oldgamma.Vect();
	Float_t thetaPQplus  = gammav_rot.Angle(piplus);
	Float_t phiPQplus    = PhiPQ(piplus,oldgamma);
	//Float_t etaplus      = piplus.Eta(); // just implement theta for now
	for ( auto i: neg_pions ) {
	  TVector3 piminus (px[i],py[i],pz[i]);
	  Float_t phiminus     = piminus.Phi();
	  Float_t thetaminus   = piminus.Theta();

	  TVector3 gammav_rot  = virtual_photon.Vect();
	  Float_t thetaPQminus  = gammav_rot.Angle(piminus);
	  Float_t phiPQminus    = PhiPQ(piminus,virtual_photon);
	  nm->Fill(thetaplus-thetaminus,phiplus-phiminus); 
	  nm_pqf->Fill(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	}
      }
    }
    
    pippim_ppp->Fill(old_piplus.size(),neg_pions.size());
    pippim_ppm->Fill(pos_pions.size(),old_piminus.size());
    pippim->Fill(pos_pions.size(),neg_pions.size());
    old_piplus.clear();
    old_piminus.clear();


    if ( m_debug )
      cout << "DEBUG: Entering same event filling" << endl;
    // SAME EVENT FILLING
    bool fill_neg = true;

    if ( m_debug )
      cout << "DEBUG: entering loop of pos_pions" << endl;
    for ( auto i: pos_pions ) {
      TVector3 piplus (px[i],py[i],pz[i]);
      Float_t phiplus      = piplus.Phi();
      Float_t thetaplus    = piplus.Theta();
      Float_t etaplus      = piplus.Eta();
      if ( m_debug )
	cout << "DEBUG: Attempting to fill some histograms" << endl;
      reco->Fill(phiplus,thetaplus);
      reco_plus->Fill(phiplus,thetaplus);
      reco_eta->Fill(phiplus,etaplus);
      if ( m_debug )
	cout << "DEBUG: Successfully filled some histograms" << endl;
      old_piplus.push_back(piplus);

      TVector3 gammav_rot  = virtual_photon.Vect();
      Float_t thetaPQplus  = gammav_rot.Angle(piplus);
      Float_t phiPQplus    = PhiPQ(piplus,virtual_photon);
      reco_PQ->Fill(phiPQplus,thetaPQplus);
      
      if ( m_debug )
	cout << "DEBUG: what seems to be the problem?" << endl;
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
	if ( fill_neg ) { // fill "acceptance" histograms once only
	  old_piminus.push_back(piminus);
	  reco->Fill(phiminus,thetaminus);
	  reco_minus->Fill(phiminus,thetaminus);
	  reco_eta->Fill(phiminus,etaminus);
	  reco_PQ->Fill(phiPQminus,thetaPQminus);
	}
      }
      fill_neg = false;
    }
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << counter << endl;
  
  //side_by_side(ns,ns_vir,"NS_lab_virtual.png");

  TString multi_method;
  if ( old_plus )
    multi_method = "_pastPiPlus";
  else if ( !old_plus )
    multi_method = "_pastPiMinus";

  
  if ( gExport_pdf ) {
    ridge_plot(ns,"single_ridge_theta.pdf");
    ridge_plot(ns_eta,"single_ridge_eta.pdf");
    ridge_plot(ns_pqf,"single_ridge_PQ.pdf");
    ridge_plot(nm,"single_ridge_theta_MULTI"+multi_method+".pdf");
    ridge_plot(nm_pqf,"single_ridge_theta_MULTI_PQ"+multi_method+".pdf");
  }

  if ( gExport_png ) {
    ridge_plot(ns,"single_ridge_theta.png");
    ridge_plot(ns_eta,"single_ridge_eta.png");
    ridge_plot(ns_pqf,"single_ridge_PQ.png");
    ridge_plot(nm,"single_ridge_theta_MULTI"+multi_method+".png");
    ridge_plot(nm_pqf,"single_ridge_theta_MULTI_PQ"+multi_method+".png");
    
    export_hist(reco,out_path+"reco_theta.png");
    export_hist(reco_eta,out_path+"reco_eta.png");
    export_hist(reco_PQ,out_path+"reco_PQ.png");
    
    export_hist(reco_plus,out_path+"reco_theta_plus.png");
    export_hist(reco_minus,out_path+"reco_theta_minus.png");

    export_hist(pippim,out_path+"pions_SameEvent.png","colztext");
    export_hist(pippim_ppp,out_path+"pions_PastPiP.png","colztext");
    export_hist(pippim_ppm,out_path+"pions_PastPiM.png","colztext");
  }

  TH2 * ndiv = new TH2F("Correlations, lab frame (theta)",
			"N_{s}/N_{m}(#Delta#theta, #Delta#phi);#Delta#theta;#Delta#phi",
			nsbinsx,-nsedgex,nsedgex,
			nsbinsy,-nsedgey,nsedgey);

  TH2 * ndiv_pqf = new TH2F("Correlations, PQ frame (theta)",
			"N_{s}/N_{m}(#Delta#theta_{PQ}, #Delta#phi_{PQ});#Delta#theta_{PQ};#Delta#phi_{PQ}",
			nsbinsx,-nsedgex,nsedgex,
			nsbinsy,-nsedgey,nsedgey);

  float ratio;
  float settle_at = 1.0;
  for ( int x = 1; x < nsbinsx+1; x++ ) {
    for ( int y = 1; y < nsbinsy+1; y++ ) {
      float nsvalue = ns->GetBinContent(x,y);
      float nmvalue = nm->GetBinContent(x,y);
      float nsvalue_pqf = ns_pqf->GetBinContent(x,y);
      float nmvalue_pqf = nm_pqf->GetBinContent(x,y);

      if ( nsvalue < 500.0 || nmvalue < 500.0 ) {
	ndiv->SetBinContent(x,y,settle_at);	
      } else {
	ndiv->SetBinContent(x,y,nsvalue/nmvalue);
      }

      if ( nsvalue_pqf < 500.0 || nmvalue_pqf < 500.0 ) {
	ndiv_pqf->SetBinContent(x,y,settle_at);
      } else {
	ndiv_pqf->SetBinContent(x,y,nsvalue_pqf/nmvalue_pqf);
      }

      /*
      if ( nm->GetBinContent(x,y) != 0 ) {
	ratio = (ns->GetBinContent(x,y))/(nm->GetBinContent(x,y));
	ndiv->GetBin(x,y);
	
	if ( ratio > 1.8 || ratio < 0.2) {
	  cout << ns->GetBinContent(x,y) << " / "
	       << nm->GetBinContent(x,y) << " = "
	       << ratio << endl;
	  ndiv->SetBinContent(x,y,settle_at);
	  continue;
	}
	ndiv->SetBinContent(x,y,ratio);

	
      } else {
	ndiv->SetBinContent(x,y,settle_at);
      }
      */
    }
  }

  if ( gExport_pdf ) {
    ridge_plot(ndiv,"single_ridge_theta_NDIV"+multi_method+".pdf");
    ridge_plot(ndiv_pqf,"single_ridge_theta_NDIV_PQ"+multi_method+".pdf");
  }
  if ( gExport_png) {
    ridge_plot(ndiv,"single_ridge_theta_NDIV"+multi_method+".png");
    ridge_plot(ndiv_pqf,"single_ridge_theta_NDIV_PQ"+multi_method+".png");
  }
  cout << "finished correctly" << endl;
} // END MAIN


////////////////////////
//// FUNCTIONS
////////////////////////

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
  //gStyle->SetOptStat(0);
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

  // comment next lines if not running on batch mode I guess
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

void initialize_histograms(int nsbinsx, int nsbinsy, double nsedgex, double nsedgey)
{
  ns = new TH2F("Signal Distribution, lab frame (theta)",
		"N_{s}(#Delta#theta, #Delta#phi);#Delta#theta;#Delta#phi",
		nsbinsx,-nsedgex,nsedgex,
		nsbinsy,-nsedgey,nsedgey);
  
  ns_eta = new TH2F("Signal Distribution, lab frame",
		    "N_{s}(#Delta#eta, #Delta#phi);#Delta#eta;#Delta#phi",
		    nsbinsx,-nsedgex,nsedgex,
		    nsbinsy,-nsedgey,nsedgey);
  
  ns_pqf = new TH2F("Signal Distribution, PQ frame",
		    "N_{s}(#Delta#theta_{PQ}, #Delta#phi_{PQ});#Delta#theta_{PQ};#Delta#phi_{PQ}",
		    nsbinsx,-nsedgex,nsedgex,
		    nsbinsy,-nsedgey,nsedgey);
  
  reco = new TH2F("Reconstructed angles",
		  "Reconstructed #pi^{#pm};#phi;#theta",
		  72,-TMath::Pi(),TMath::Pi(),
		  36,0.0,TMath::Pi());
  
  reco_plus = new TH2F("Reconstructed angles, #pi^{+}",
		       "Reconstructed #pi^{+};#phi;#theta",
		       72,-TMath::Pi(),TMath::Pi(),
		       36,0.0,TMath::Pi());
  
  reco_minus = new TH2F("Reconstructed angles, #pi^{-}",
			"Reconstructed #pi^{-};#phi;#theta",
			72,-TMath::Pi(),TMath::Pi(),
			36,0.0,TMath::Pi());
  
  
  reco_PQ = new TH2F("Reconstructed angles, PQ frame",
		     "Reconstructed #pi^{#pm};#phi_{PQ};#theta_{PQ}",
		     72,-TMath::Pi(),TMath::Pi(),
		     36,0.0,TMath::Pi());
  
  reco_eta = new TH2F("Reconstructed pions, eta",
		      "Reconstructed #pi^{#pm};#phi;#eta",
		      72,-TMath::Pi(),TMath::Pi(),
		      80,-2.0,4.0);
  
  nm = new TH2F("Multi events, lab frame (theta)",
		"N_{m}(#Delta#theta, #Delta#phi);#Delta#theta;#Delta#phi",
		nsbinsx,-nsedgex,nsedgex,
		nsbinsy,-nsedgey,nsedgey);
  nm_pqf = new TH2F("Multi events, PQ frame (theta)",
		    "N_{m}(#Delta#theta_{PQ}, #Delta#phi_{PQ});#Delta#theta_{PQ};#Delta#phi_{PQ}",
		    nsbinsx,-nsedgex,nsedgex,
		    nsbinsy,-nsedgey,nsedgey);  
}
