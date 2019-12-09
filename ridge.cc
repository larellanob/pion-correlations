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
#include "TCorrelation.cxx"

TString out_path = "/home/luciano/Physics/CLAS/pion_ridge/";
bool m_debug      = false;
bool m_simulation = false;
bool old_plus     = false; // if false, use old_minus
bool gExport_pdf  = false;
bool gExport_png  = true;
Float_t gDataCap  = 0.01; // percentage of data to be used in analysis
TDatabasePDG db;

Float_t PhiPQ(TVector3, TLorentzVector);
void export_hist(TH2 *, TString, TString options = "colz");
void ridge_plot(TH2F h2, TString out_name, TString out_path_modifier = "" );
void full_ridge_plots(TCorrelation* );
void side_by_side(TH2 *,TH2 *, TString);
void initialize_histograms(int nsbinsx, int nsbinsy, double nsedgex, double nsedgey);


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
    //f = new TFile("/home/luciano/Physics/CLAS/data/CFFTree_data.root");
    f = new TFile("/home/luciano/Physics/CLAS/data/full_C_files.root");
    kEbeam = 5.014;
  }
  const Double_t kMpi   = 0.13957;
  const Double_t kMprt  = 0.938272;
  float rhomass = db.GetParticle(113)->Mass();
  TTree *auxtree = (TTree*) f->GetKey("tree_data")->ReadObj();
  Int_t Entries = 0;
  Entries = auxtree->GetEntries();
  
  // tree reader for thrown particles and scattered electron
  TTreeReader th("tree_data",f);
  TTreeReader e("e_rec",f);
    
  TTreeReaderArray<Int_t> pid(th,"pid");
  TTreeReaderValue<Float_t> q2(th,"Q2");
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
    
  if ( m_debug )
    cout << "DEBUG: entering main tree reading loop" << endl;
  // main file/tree reading loop
  TLorentzVector virtual_photon;
  TLorentzVector oldgamma;

  TCorrelation * corr_ang = new TCorrelation("phi","theta");
  TCorrelation * corr_apq = new TCorrelation("phi_{PQ}","theta_{PQ}");
  TCorrelation * corr_eta = new TCorrelation("eta","theta");
  
  while ( th.Next()  && events < (gDataCap * Entries) ) {
    e.Next();
    events++;
    
    Double_t scattered_energy = sqrt((*epx)*(*epx)+(*epy)*(*epy)+(*epz)*(*epz));
    if ( scattered_energy < 0.005 )
      cout << "WARNING: Low e' energy " << scattered_energy << endl; 
    if ( scattered_energy > kEbeam  )
      cout << "WARNING: High e' energy = " <<  scattered_energy << endl;


    // virtual photon frame is per event
    //oldgamma = virtual_photon; // true old gamma

    virtual_photon = {-(*epx),-(*epy),kEbeam-(*epz),kEbeam-scattered_energy};
    oldgamma = virtual_photon; // actually the same old gamma
    
       
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
      bool do_negativee = true;
      for ( auto i: pos_pions ) {
	TVector3 piplus (px[i],py[i],pz[i]);
	Float_t phiplus      = piplus.Phi();
	Float_t thetaplus    = piplus.Theta();
	
	TVector3 gammav_rot  = virtual_photon.Vect();
	Float_t thetaPQplus  = gammav_rot.Angle(piplus);
	Float_t phiPQplus    = PhiPQ(piplus,virtual_photon);
	//Float_t etaplus      = piplus.Eta(); // just implement theta for now
	//angle2->Fill(abs(piplus.Phi()-gammav_rot.Phi()));
	//cout << "------------------------------------" << endl;
	for ( auto piminus: old_piminus ) {
	  Float_t phiminus     = piminus.Phi();
	  Float_t thetaminus   = piminus.Theta();

	  TVector3 gammav_rot  = oldgamma.Vect();
	  Float_t thetaPQminus  = gammav_rot.Angle(piminus);
	  Float_t phiPQminus    = PhiPQ(piminus,oldgamma);
	  corr_ang->FillMulti(thetaplus-thetaminus,phiplus-phiminus);
	  corr_apq->FillMulti(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	  if ( do_negativee ) {
	    //angle1->Fill(abs(piminus.Phi()-gammav_rot.Phi()));
	    /*
	    cout << "-----" << endl;
	    cout << "gamma, piminus, phipqminus: " << endl;;
	    gammav_rot.Print();
	    piminus.Print();
	    cout << phiPQminus << " " << thetaPQminus << endl;
	    */
	  }
	}
	do_negativee = false;
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
	  corr_ang->FillMulti(thetaplus-thetaminus,phiplus-phiminus);
	  corr_apq->FillMulti(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	}
      }
    }
    /*
    pippim_ppp->Fill(old_piplus.size(),neg_pions.size());
    pippim_ppm->Fill(pos_pions.size(),old_piminus.size());
    pippim->Fill(pos_pions.size(),neg_pions.size());
    */
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
      corr_ang->FillReco(phiplus,thetaplus);
      corr_ang->FillRecoPlus(phiplus,thetaplus);
      corr_eta->FillReco(phiplus,etaplus);
      if ( m_debug )
	cout << "DEBUG: Successfully filled some histograms" << endl;
      old_piplus.push_back(piplus);
      
      TVector3 gammav_rot  = virtual_photon.Vect();
      Float_t thetaPQplus  = gammav_rot.Angle(piplus);
      Float_t phiPQplus    = PhiPQ(piplus,virtual_photon);
      corr_apq->FillReco(phiPQplus,thetaPQplus);
      
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

	corr_ang->FillSame(thetaplus-thetaminus,phiplus-phiminus);
	corr_eta->FillSame(etaplus-etaminus,phiplus-phiminus);
	corr_apq->FillSame(thetaPQplus-thetaPQminus,phiPQplus-phiPQminus);
	if ( fill_neg ) { // fill "acceptance" histograms once only
	  old_piminus.push_back(piminus);
	  corr_ang->FillReco(phiminus,thetaminus);
	  corr_apq->FillReco(phiPQminus,thetaPQminus);
	  corr_eta->FillReco(phiminus,etaminus);
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

  
  corr_ang->FillCorrelation();
  corr_apq->FillCorrelation();

  full_ridge_plots(corr_ang);
  full_ridge_plots(corr_apq);
  //full_ridge_plots(corr_eta);
  

  cout << "finished correctly" << endl;
} // END MAIN


////////////////////////
//// FUNCTIONS
////////////////////////

Float_t PhiPQ(TVector3 pion, TLorentzVector gamma_vir )
{
  TVector3 Vhelp(0.,0.,1.0);
  TVector3 gamma = gamma_vir.Vect();
  if ( m_debug ) {
    cout << "0 gamm: ";
    gamma.Print();
    cout << "0 pion: ";
    pion.Print();
    cout << "0 betw: ";
    cout << gamma.Angle(pion) << endl;
  }
  Double_t phi_z = TMath::Pi()-gamma.Phi();
  pion.RotateZ(phi_z);
  gamma.RotateZ(phi_z);
  if ( m_debug ) {
    cout << "1 gamm: ";
    gamma.Print();
    cout << "1 pion: ";
    pion.Print();
    cout << "1 betw: ";
    cout << gamma.Angle(pion) << endl;
  }
  
  Double_t phi_y = gamma.Angle(Vhelp);
  pion.RotateY(phi_y);
  if ( m_debug ) {
    gamma.RotateY(phi_y);
    cout << "2 gamm: ";
    gamma.Print();
    cout << "2 pion: ";
    pion.Print();
    cout << "2 betw: ";
    cout << gamma.Angle(pion) << endl;
  }
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


void ridge_plot(TH2F h2, TString out_name, TString out_path_modifier)
{
  TCanvas *c1 = new TCanvas("c1","mi cambas",1000,1000);
  c1->SetLeftMargin(0.14);
  h2.Draw("surf1");

  // aesthetics
  //gStyle->SetOptStat(0);
  double theta = -30;
  double phi   = 60;
  h2.GetXaxis()->CenterTitle(true);
  h2.GetXaxis()->SetTitleOffset(1.5);
  h2.GetYaxis()->CenterTitle(true);

  // main 2d plot 
  TString out_filename = out_path+out_path_modifier+out_name;
  c1->SaveAs(out_filename);

  // mountain side plots
  gPad->GetView()->RotateView(0,90);
  out_filename = out_path+out_path_modifier+"hillsideY_"+out_name;
  c1->SaveAs(out_filename);
  gPad->GetView()->RotateView(-90,90);
  out_filename = out_path+out_path_modifier+"hillsideX_"+out_name;
  c1->SaveAs(out_filename);

  // projection plots
  TH1 * projx = h2.ProjectionX();
  projx->SetFillColor(kBlue);
  projx->Draw();
  out_filename = out_path+out_path_modifier+"projx_"+out_name;
  c1->SaveAs(out_filename);
  
  TH1 * projy = h2.ProjectionY();
  projy->SetFillColor(kBlue);
  projy->Draw();
  out_filename = out_path+out_path_modifier+"projy_"+out_name;
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


void full_ridge_plots(TCorrelation * corr)
{
  TString connector = corr->GetVar(1)+"_"+corr->GetVar(2);
  TString connector2 = connector+"/";
  ridge_plot(corr->GetSE(),"ridge_same_"+connector+".png", connector2);
  ridge_plot(corr->GetME(),"ridge_multi_"+connector+".png", connector2);
  ridge_plot(corr->GetCO(),"ridge_correlation_"+connector+".png", connector2);
}
