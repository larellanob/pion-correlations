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
bool DMode        = false; // deuterium
Float_t gDataCap  = 0.1; // percentage of data to be used in analysis
TDatabasePDG db;


//Float_t PhiPQ(TVector3, TLorentzVector);
void export_hist(TH2 *, TString, TString options = "colz");
//void export_hist(TH1F, TString, TString options = "");
void ridge_plot(TH2F h2, TString out_name, TString out_path_modifier = "" );
void full_ridge_plots(TCorrelation* );
void full_1d_plots(TCorrelation*);
void side_by_side(TH2 *,TH2 *, TString);
Float_t ComputeThetaPQ(TLorentzVector, TLorentzVector);
Float_t ComputePhiPQ(TLorentzVector, TLorentzVector);
TLorentzVector VirtualFrame(TLorentzVector, double, double, double);
Float_t DeltaAngleRad(Float_t x_rad, Float_t y_rad);


// MAIN
void ridge(TString target = "")
{
  //gStyle->SetOptStat(0);
  if ( m_simulation == true ) {
    out_path+="simulation/";
  } else {
    out_path += "data/";
  }

  TChain ch("tree_data");
  TChain ech("e_rec");
  //TChain *f;
  Double_t kEbeam;
  if ( m_simulation ) {
    //f = new TChain("/home/luciano/Physics/CLAS/data/tree_output.root");
    ch.Add("/home/luciano/Physics/CLAS/data/tree_output.root");
    ech.Add("/home/luciano/Physics/CLAS/data/tree_output.root");
    kEbeam = 11.0;
  } else {
    //f = new TChain("/home/luciano/Physics/CLAS/data/CFFTree_data.root");
    if ( target == "Pb" && !DMode ) {
      //f = new TChain("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
      ch.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
      ech.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    } else if ( target == "Fe" && !DMode) {
      //f = new TChain("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
      ch.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
      ech.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    } else {
      if ( !DMode ) {
	target = "C";
	//f = new TChain("/home/luciano/Physics/CLAS/data/full_C_files.root");
	ch.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
	ech.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
      }
    }
    kEbeam = 5.014;
  }
  if ( DMode ) {
    target = "D";
    ch.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    ch.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    ch.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
  }
  
  Double_t kChargedPionMass = db.GetParticle(211)->Mass();
  Double_t kProtonMass = db.GetParticle(2212)->Mass();
  Double_t kRhoMass = db.GetParticle(113)->Mass();
  //TTree *auxtree = (TTree*) f->GetKey("tree_data")->ReadObj();
  
  Int_t Entries = 0;
  //Entries = auxtree->GetEntries();
  Entries = ch.GetEntries();
  
  // tree reader for thrown particles and scattered electron
  TTreeReader th(&ch);
  TTreeReader e(&ech);
  //TTreeReader th("tree_data",f);
  //TTreeReader e("e_rec",f);
    
  TTreeReaderArray<Int_t> pid(th,"pid");
  TTreeReaderValue<Float_t> q2(th,"Q2");
  TTreeReaderArray<Float_t> px(th,"Px");
  TTreeReaderArray<Float_t> py(th,"Py");
  TTreeReaderArray<Float_t> pz(th,"Pz");

  TTreeReaderValue<Int_t> TargType(th,"TargType");
  TTreeReaderArray<Float_t> Zh(th,"Zh");

  TTreeReaderArray<Float_t> Theta(th,"Theta");
  TTreeReaderArray<Float_t> Phi(th,"Phi");
  TTreeReaderArray<Float_t> ThetaPQ(th,"ThetaPQ");
  TTreeReaderArray<Float_t> PhiPQ(th,"PhiPQ");
  
  TTreeReaderValue<Float_t> epx(e,"Pex");
  TTreeReaderValue<Float_t> epy(e,"Pey");
  TTreeReaderValue<Float_t> epz(e,"Pez");
  TTreeReaderValue<Float_t> eq2(e,"Q2");
  TTreeReaderValue<Float_t> ew(e,"W");
    
  cout << "working" << endl;

  std::vector<int> pos_pions;
  std::vector<int> neg_pions;

  std::vector<TLorentzVector> piplus4v;
  std::vector<TLorentzVector> piminus4v;

  std::vector<TLorentzVector> piplus4v_rotated;
  std::vector<TLorentzVector> piminus4v_rotated;

  std::vector<TLorentzVector> piplus4v_boosted;
  std::vector<TLorentzVector> piminus4v_boosted;

  std::vector<TLorentzVector> old4v;
  std::vector<TLorentzVector> old4v_rotated;
  std::vector<TLorentzVector> old4v_boosted;
  
  
  UInt_t events = 0;
  int counter=0;
    
  if ( m_debug )
    cout << "DEBUG: entering main tree reading loop" << endl;
  // main file/tree reading loop
  TLorentzVector virtual_photon;
  //TLorentzVector oldgamma;

  TCorrelation * corr_ang = new TCorrelation("phi","theta",target);
  TCorrelation * corr_apq = new TCorrelation("phiPQ","thetaPQ", target);
  TCorrelation * corr_boo = new TCorrelation("phiPQboosted","thetaPQboosted",target);
  //TCorrelation * corr_eta = new TCorrelation("eta","theta", target);
  //TCorrelation * corr_vir = new TCorrelation("vphi","vtheta", target);

  ULong_t NTriggers = 0;
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
    //oldgamma = virtual_photon; // actually the same old gamma

    if ( events % 1000000 == 0 ) {
      cout << events << "/" << gDataCap*Entries << endl;
    }
    

    //////////////////////////
    // event cuts
    if ( *ew < 2.0 || *eq2 < 1.0 ) // DIS cut
      continue;
    if ( *TargType == 1 && !DMode ) { // D == 1; heavy nuclei == 2
      continue;
    } else if ( *TargType != 1 && DMode ) {
      continue;
    }
      
    

    //////////////////////////
    // pions in the event
    pos_pions.clear();
    neg_pions.clear();
    
    // pi+ pi+
    /*
      bool trigger = false;
    for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
      if ( pid[i] == 211  && !trigger ) {
	pos_pions.push_back(i);
	trigger = true;
	continue;
      }
      if ( pid[i] == 211 && trigger ) {
	neg_pions.push_back(i);
      }
    }
    */

    // pi+ pi-
    for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
      if ( pid[i] == 211 ) {
	pos_pions.push_back(i);
      }
      if ( pid[i] == -211 ) {
	neg_pions.push_back(i);
      }
    }
    
    // make sure of this, normalization might depend on it
    if ( pos_pions.size() == 0  || neg_pions.size() == 0 )
      continue;



    ///////////////////////////
    // Pions 4vectors
    piplus4v.clear();
    piminus4v.clear();
    piplus4v_rotated.clear();
    piminus4v_rotated.clear();
    piplus4v_boosted.clear();
    piminus4v_boosted.clear();
    
    for ( int i = 0; i < pid.GetSize(); i++ ) {
      if ( pid[i] == 211 ) {
	Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
	TLorentzVector pi (px[i],py[i],pz[i],e);
	piplus4v.push_back(pi);
      } else if ( pid[i] == -211 ) {
	Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
	TLorentzVector pi (px[i],py[i],pz[i],e);
	piminus4v.push_back(pi);
      }
    }

    //////////////////////////////
    // VIRTUAL PHOTON FRAME
    // event virtual photon frame

    // rotation 1
    Double_t Egamma = virtual_photon.E();
    Double_t angle1 = - atan2((*epy),(*epx));    
    virtual_photon.RotateZ(angle1);

    // rotation 2
    double numer = virtual_photon.Z();
    double denom = virtual_photon.P();
    Double_t angle2;
    angle2 = acos(numer/denom);
    virtual_photon.RotateY(angle2);
    
    // boost
    Double_t boost  = virtual_photon.P()/(Egamma+kProtonMass);
    virtual_photon.Boost(0,0,-boost);


    //////////////////////////////
    // rotation and boost for pions
    for ( int i = 0; i < piplus4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(piplus4v[i],angle1,angle2,0);
      piplus4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(piplus4v[i],angle1,angle2,-boost);
      piplus4v_boosted.push_back(aux_boost);
    }

    for ( int i = 0; i < piminus4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(piminus4v[i],angle1,angle2,0);
      piminus4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(piminus4v[i],angle1,angle2,-boost);
      piminus4v_boosted.push_back(aux_boost);
    }

    // rotation and boost for old particle
    for ( int i = 0; i < old4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(old4v[i],angle1,angle2,0);
      old4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(old4v[i],angle1,angle2,-boost);
      old4v_boosted.push_back(aux_boost);
      
    }

    // debug
    if ( m_debug ) {
      cout << "---------------------------------------------------------------" << endl;
      cout << "pospions: " << pos_pions.size() << endl;
      for ( int i = 0; i < pos_pions.size(); i++ ) {
	cout << "pz: " << pz[pos_pions[i]] << endl;
	cout << "theta: " << Theta[pos_pions[i]] << endl;
	cout << "phi: " << Phi[pos_pions[i]] << endl;
	cout << "thetapq: " << ThetaPQ[pos_pions[i]] << endl;
	cout << "phipq: " << PhiPQ[pos_pions[i]] << endl;
	cout << "-------------" << endl;
      }
      
      cout << "4v: " << piplus4v.size() << endl;
      for ( int i = 0; i < piplus4v.size(); i++ ) {
	cout << "-----------------------------" << endl;
	cout << "pz: " << piplus4v[i].Pz() << endl;
	cout << "theta: " << piplus4v[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi: " << piplus4v[i].Phi()*180./TMath::Pi() << endl;
	cout << "-------------" << endl;
	cout << "pz rot: " << piplus4v_rotated[i].Pz() << endl;
	cout << "theta rot: " << piplus4v_rotated[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi rot: " << piplus4v_rotated[i].Phi()*180./TMath::Pi() << endl;
	cout << "-------------" << endl;
	cout << "pz boo: " << piplus4v_boosted[i].Pz() << endl;
	cout << "theta boo: " << piplus4v_boosted[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi boo: " << piplus4v_boosted[i].Phi()*180./TMath::Pi()<< endl;
      }
    }
    
    
    if ( m_debug )
      cout << "DEBUG: Entering mutli event filling" << endl;


    ///////////////////////////
    // multi event filling
    if ( old4v.size() != 0 && !old_plus && events > 1 ) {
      for ( int i = 0; i < piplus4v.size(); i++ ) {
	for ( int j = 0; j < old4v.size(); j++ ) {
	  Float_t x,y;

	  x = DeltaAngleRad(piplus4v[i].Phi(),old4v[j].Phi());
	  y = DeltaAngleRad(piplus4v[i].Theta(),old4v[j].Theta());
	  corr_ang->FillMulti(x,y);

	  x = DeltaAngleRad(piplus4v_rotated[i].Phi(),old4v_rotated[j].Phi());
	  y = DeltaAngleRad(piplus4v_rotated[i].Theta(),old4v_rotated[j].Theta());
	  corr_apq->FillMulti(x,y);
	  
	  x = DeltaAngleRad(piplus4v_boosted[i].Phi(),old4v_boosted[j].Phi());
	  y = DeltaAngleRad(piplus4v_boosted[i].Theta(),old4v_boosted[j].Theta());
	  corr_boo->FillMulti(x,y);
	}
      }
    }
    

    old4v.clear();
    old4v_rotated.clear();
    old4v_boosted.clear();
    
    if ( m_debug )
      cout << "DEBUG: Entering same event filling" << endl;
    // SAME EVENT FILLING
    bool fill_neg = true;

    if ( m_debug )
      cout << "DEBUG: entering loop of pos_pions" << endl;

    ////////////////////////////////////
    /// SAME EVENT FILLING

    for ( int i = 0; i < piplus4v.size(); i++ ) {
      corr_ang->FillReco(piplus4v[i].Phi()*180./TMath::Pi(), piplus4v[i].Theta()*180./TMath::Pi());
      corr_apq->FillReco(piplus4v_rotated[i].Phi()*180./TMath::Pi(), piplus4v_rotated[i].Theta()*180./TMath::Pi() );
      corr_boo->FillReco(piplus4v_boosted[i].Phi()*180./TMath::Pi(), piplus4v_boosted[i].Theta()*180./TMath::Pi() );

      corr_ang->FillRecoPlus(piplus4v[i].Phi()*180./TMath::Pi(), piplus4v[i].Theta()*180./TMath::Pi());
      corr_apq->FillRecoPlus(piplus4v_rotated[i].Phi()*180./TMath::Pi(), piplus4v_rotated[i].Theta()*180./TMath::Pi() );
      corr_boo->FillRecoPlus(piplus4v_boosted[i].Phi()*180./TMath::Pi(), piplus4v_boosted[i].Theta()*180./TMath::Pi() );
      for ( int j = 0; j < piminus4v.size(); j++ ) {
	float x,y;

	x = DeltaAngleRad(piplus4v[i].Phi(),piminus4v[j].Phi());
	y = DeltaAngleRad(piplus4v[i].Theta(),piminus4v[j].Theta());
	corr_ang->FillSame(x,y);

	x = DeltaAngleRad(piplus4v_rotated[i].Phi(),piminus4v_rotated[j].Phi());
	y = DeltaAngleRad(piplus4v_rotated[i].Theta(),piminus4v_rotated[j].Theta());	
	corr_apq->FillSame(x,y);
	
	x = DeltaAngleRad(piplus4v_boosted[i].Phi(),piminus4v_boosted[j].Phi());
	y = DeltaAngleRad(piplus4v_boosted[i].Theta(),piminus4v_boosted[j].Theta());
	corr_boo->FillSame(x,y);

	// enters this loop only once per event
	if ( fill_neg ) {
	  corr_ang->FillReco(piplus4v[j].Phi()*180./TMath::Pi(), piplus4v[j].Theta()*180./TMath::Pi());
	  corr_apq->FillReco(piplus4v_rotated[j].Phi()*180./TMath::Pi(), piplus4v_rotated[j].Theta()*180./TMath::Pi() );
	  corr_boo->FillReco(piplus4v_boosted[j].Phi()*180./TMath::Pi(), piplus4v_boosted[j].Theta()*180./TMath::Pi() );

	  corr_ang->FillRecoMinus(piplus4v[j].Phi()*180./TMath::Pi(), piplus4v[j].Theta()*180./TMath::Pi());
	  corr_apq->FillRecoMinus(piplus4v_rotated[j].Phi()*180./TMath::Pi(), piplus4v_rotated[j].Theta()*180./TMath::Pi() );
	  corr_boo->FillRecoMinus(piplus4v_boosted[j].Phi()*180./TMath::Pi(), piplus4v_boosted[j].Theta()*180./TMath::Pi() );
	}
      }
      fill_neg = false;
    }

    //////////////////////////////////////
    // old pion saving

    // old negative pion:
    for ( int j = 0; j < piminus4v.size(); j++ ) {
      old4v.push_back(piminus4v[j]);
      //old4v_rotated.push_back(piminus4v_rotated[j]);
      //old4v_boosted.push_back(piminus4v_boosted[j]);
    }
    

    NTriggers += pos_pions.size();
    //    if ( pos_pions.size() < 1 )
    //cout << "what " << pos_pions.size() << endl;
    counter++;
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << counter << endl;
  cout << "number of trigger particles: " << NTriggers << endl;
  //side_by_side(ns,ns_vir,"NS_lab_virtual.png");

  TString multi_method;
  if ( old_plus )
    multi_method = "_pastPiPlus";
  else if ( !old_plus )
    multi_method = "_pastPiMinus";

  //corr_apq->NormalizeSame(NTriggers);
  //corr_apq->NormalizeMulti();

  
  corr_ang->FillCorrelation();
  corr_apq->FillCorrelation();
  corr_boo->FillCorrelation();
  
  full_ridge_plots(corr_ang);
  full_ridge_plots(corr_apq);
  full_ridge_plots(corr_boo);

  full_1d_plots(corr_ang);
  full_1d_plots(corr_apq);
  full_1d_plots(corr_boo);
  
  TH2F * reco1 = new TH2F("reco1","test",10,0,10,10,0,10);
  TH2F * reco2 = new TH2F("reco2","test",10,0,10,10,0,10);
  TH2F * reco3 = new TH2F("reco3","test",10,0,10,10,0,10);
  *reco1 = corr_apq->GetReco();
  *reco2 = corr_apq->GetRecoPlus();
  *reco3 = corr_apq->GetRecoMinus();

  new TCanvas();
  reco1->Draw("colz");
  new TCanvas();
  reco2->Draw("colz");
  new TCanvas();
  reco3->Draw("colz");
  /*
  new TCanvas();
  corr_apq->GetReco().Draw("colz");
  new TCanvas();
  corr_apq->GetRecoPlus().Draw("colz");
  new TCanvas();
  corr_apq->GetRecoMinus().Draw("colz");
  */
  /*
  export_hist(corr_ang->GetCorr1D(1), out_path+"1D_"+corr_ang->GetVar(1)+".png");
  export_hist(corr_ang->GetCorr1D(2), out_path+"1D_"+corr_ang->GetVar(2)+".png");
  export_hist(corr_apq->GetCorr1D(1), out_path+"1D_"+corr_apq->GetVar(1)+".png");
  export_hist(corr_apq->GetCorr1D(2), out_path+"1D_"+corr_apq->GetVar(2)+".png");
  export_hist(corr_boo->GetCorr1D(1), out_path+"1D_"+corr_boo->GetVar(1)+".png");
  export_hist(corr_boo->GetCorr1D(2), out_path+"1D_"+corr_boo->GetVar(2)+".png");
  */
  //full_ridge_plots(corr_eta);
  

  cout << "finished correctly" << endl;
} // END MAIN


////////////////////////
//// FUNCTIONS
////////////////////////

Float_t ComputePhiPQ(TLorentzVector pion, TLorentzVector gamma_vir )
{
  TVector3 Vhelp(0.,0.,1.0);
  TLorentzVector gamma = gamma_vir;
  Double_t phi_z = TMath::Pi()-gamma.Phi();
  pion.RotateZ(phi_z);
  gamma.RotateZ(phi_z);  
  Double_t phi_y = gamma.Angle(Vhelp);
  pion.RotateY(phi_y);
  return ( (pion.Phi()*180.)/TMath::Pi() );
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

void export_hist(TH1F h2, TString out_filename, TString options = "") {
  auto c1 = new TCanvas();
  c1->SetCanvasSize(800,800);
  //c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  h2.Draw(options);
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
  h2.GetYaxis()->SetTitleOffset(2.0);
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

void full_1d_plots(TCorrelation* corr)
{
  for ( int var = 1; var < 3; var++ ) {
    TString con = "1D_"+corr->GetVar(var)+"_";
    export_hist(corr->GetCorr1D(var), out_path+con+"corr.png");
    export_hist(corr->GetSame1D(var), out_path+con+"same.png");
    export_hist(corr->GetMult1D(var), out_path+con+"mult.png");
  }
}

TLorentzVector VirtualFrame(TLorentzVector v1,
			    double angle1, double angle2, double boost)
{
  /*
  TLorentzVector *v2;
  v2 = &v1;
  v2->RotateZ(-angle1);
  v2->RotateY(-angle2);
  v2->Boost(0,0,-boost);
  return *v2;
  */
  v1.RotateZ(angle1);
  v1.RotateY(angle2);
  v1.Boost(0,0,-boost);
  return v1;
}

Float_t DeltaAngleRad(Float_t x_rad, Float_t y_rad)
{
  Double_t x = x_rad*180./TMath::Pi();
  Double_t y = y_rad*180./TMath::Pi();

  Float_t result;
  
  if ( ( x < 0 && y > 0) || (x > 0 && y < 0 ) ) {
    result = abs(y+x);
  }
  if ( x > 0 && y > 0 ) {
    result = abs(x-y);
  }
  if ( x < 0 && y < 0 ) {
    result = abs(x-y);
  }

  if ( result > 180. ) {
    result = 360.-result;
  }

  return result;
}
