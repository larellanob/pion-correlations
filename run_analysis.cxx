#include "Root/TCorrelation.cxx"
#include "Root/Plotter.cc"
#include "Root/DeltaAngle.cxx"

void export_hist(TH1 * h2, TString out_filename);
void recoplots(TCorrelation *t, TString target, TString gMode, TString var);
TString gDirName = "/home/luciano/Physics/CLAS/pion_correlation/";
bool m_simulation = true;

//void run_analysis( TString filename = "/home/luciano/Physics/CLAS/pion_ridge/data/2pcPairs_D_zhpm.root" )

void run_analysis( TString centroid_s = "",TString target = "C", TString gMode = "zhpm")
{
  if ( m_simulation ) {
    cout << "running in simulation mode!" << endl;
    gDirName=gDirName+"simulations/";
  }
  TString filename = gDirName+"2pcPairs/2pcPairs_"+target+"_"+gMode+".root";
  if ( filename == "" ) {
    cout << "Please input some filename to analyze" << endl;
    return;
  }
  double centroid = centroid_s.Atof();
  TH1F h_inv_mass("h_inv_mass","Invariant mass of trigger+partner "+target+";m_{t}+m_{p} [GeV]",100,0,2);
  TFile f(filename);

  // metadata
  // mode, target, % of data, etc.
  //TString target = "D";

  
  TTreeReader triggers("triggers",&f);
  TTreeReader partners("partners",&f);

  // triggers 
  TTreeReaderArray<Int_t>  pid_t(triggers,"pid");
  TTreeReaderValue<Float_t> q2_t(triggers,"Q2");
  TTreeReaderArray<Float_t> Zh_t(triggers,"Zh");

  TTreeReaderArray<Float_t> TheBoo_t(triggers,"ThetaBoosted");
  TTreeReaderArray<Float_t> PhiBoo_t(triggers,"PhiBoosted");
  TTreeReaderArray<Float_t> RapBoo_t(triggers,"y"); // change this to eta after generating new files
  
  TTreeReaderArray<Float_t> Px_t(triggers,"Px");
  TTreeReaderArray<Float_t> Py_t(triggers,"Py");
  TTreeReaderArray<Float_t> Pz_t(triggers,"Pz");
  TTreeReaderArray<Float_t> E_t(triggers,"E");
  
  // partners 
  TTreeReaderArray<Int_t>  pid_p(partners,"pid");
  TTreeReaderValue<Float_t> q2_p(partners,"Q2");
  TTreeReaderArray<Float_t> Zh_p(partners,"Zh");

  TTreeReaderArray<Float_t> TheBoo_p(partners,"ThetaBoosted");
  TTreeReaderArray<Float_t> PhiBoo_p(partners,"PhiBoosted");
  TTreeReaderArray<Float_t> RapBoo_p(partners,"y"); // change this to eta after generating new files

  TTreeReaderArray<Float_t> Px_p(partners,"Px");
  TTreeReaderArray<Float_t> Py_p(partners,"Py");
  TTreeReaderArray<Float_t> Pz_p(partners,"Pz");
  TTreeReaderArray<Float_t> E_p(partners,"E");

  
  TCorrelation * corr_boo = new TCorrelation("phiPQboosted","thetaPQboosted",target);
  TCorrelation * corr_rap = new TCorrelation("phiPQboosted2","y",target);
  int p_size,t_size;
  UInt_t events = 0;
  Int_t entries = partners.GetEntries();
  bool rho_warning = false;
  int rhos = 0;
  //while ( triggersp.Next() && events < 1.0*entries  && events < 800000) {

  while ( triggers.Next() ) {
    events++;
    t_size = pid_t.GetSize();


    // you can get "old partners" by working here, before partners.Next()

    /// Mixed Event Pairs
    if ( centroid == 0 ) {
      rho_warning= true;
    }
    if ( events > 1  && rho_warning ) {
      for ( int t = 0; t < t_size; t++ ) {
	for ( int o = 0; o < p_size; o++ ) {
	  Float_t DPhi,DThe,DRap;
	  DPhi = DeltaAngle(PhiBoo_t[t],PhiBoo_p[o]);
	  DThe = DeltaAngle(TheBoo_t[t],TheBoo_p[o]);
	  DRap = RapBoo_t[t]-RapBoo_p[o];
	  corr_boo->FillMulti(DPhi,DThe);
	  corr_rap->FillMulti(DPhi,DRap);
	}
      }
    }

    partners.Next();
    int old_p_size = p_size;
    p_size = pid_p.GetSize();
    // rho veto

    rho_warning = false;
    if ( centroid == 0 ) {
      rho_warning= true;
    }
    for ( int t = 0; t < t_size; t++ ) {
      if ( rho_warning ) {
	break;
      }
      TLorentzVector trig4v(Px_t[t],Py_t[t],Pz_t[t],E_t[t]);
      for ( int p = 0; p < p_size; p++ ) {
	TLorentzVector part4v(Px_p[p],Py_p[p],Pz_p[p],E_p[p]);
	double invmass = (trig4v+part4v).M();
	h_inv_mass.Fill(invmass);
	if ( invmass > centroid-0.1  && invmass < centroid+0.1 ) {
	  rho_warning = true;
	  rhos++;
	}
	if ( rho_warning ) {
	  break;
	}
      }
    }
    if ( !rho_warning ) {
      continue;
    }
        
    /// Same Event Pairs

    bool fill_neg = true;
    for ( int t = 0; t < t_size; t++ ) {
      corr_boo->FillReco(PhiBoo_t[t],TheBoo_t[t]);
      corr_rap->FillReco(PhiBoo_t[t],RapBoo_t[t]);
      corr_boo->FillRecoTriggers(PhiBoo_t[t],TheBoo_t[t]);
      corr_rap->FillRecoTriggers(PhiBoo_t[t],RapBoo_t[t]);
      
      
      for ( int p = 0; p < p_size; p++ ) {
	if ( gMode == "pp" && t == p ) {
	  continue;
	}
	//if ( t == 0 ) {
	if ( fill_neg) {
	  corr_boo->FillReco(PhiBoo_p[p],TheBoo_p[p]);
	  corr_rap->FillReco(PhiBoo_p[p],RapBoo_p[p]);
	  corr_boo->FillRecoPartners(PhiBoo_p[p],TheBoo_p[p]);
	  corr_rap->FillRecoPartners(PhiBoo_p[p],RapBoo_p[p]);
	}
	
	Float_t DPhi,DThe,DRap;
	DPhi = DeltaAngle(PhiBoo_t[t],PhiBoo_p[p]);
	DThe = DeltaAngle(TheBoo_t[t],TheBoo_p[p]);
	DRap = RapBoo_t[t]-RapBoo_p[p];
	corr_boo->FillSame(DPhi,DThe);
	corr_rap->FillSame(DPhi,DRap);
      }
      fill_neg = false;
    }
    

    
  }
  cout << "rhos " << rhos << endl;
  corr_boo->FillCorrelation();
  corr_rap->FillCorrelation();

  TFile * fout = new TFile(gDirName+"histograms/hist_"+target+"_"+gMode+"_"+centroid_s+"_boo.root","recreate");
  corr_boo->GetCorr1D(1).Write();
  corr_boo->GetCorr1D(2).Write();
  corr_boo->GetSame1D(1).Write();
  corr_boo->GetSame1D(2).Write();
  corr_boo->GetMult1D(1).Write();
  corr_boo->GetMult1D(2).Write();
  
  TFile * fout2 = new TFile(gDirName+"histograms/hist_"+target+"_"+gMode+"_"+centroid_s+"_rap.root","recreate");
  corr_rap->GetCorr1D(1).Write();
  corr_rap->GetCorr1D(2).Write();
  corr_rap->GetSame1D(1).Write();
  corr_rap->GetSame1D(2).Write();
  corr_rap->GetMult1D(1).Write();
  corr_rap->GetMult1D(2).Write();

  // Reconstructed pions plots
  // initialization of arbitrary TH2Fs to get from corr_*
  recoplots(corr_boo,target, gMode, "boo");
  recoplots(corr_rap,target, gMode, "rap");



  // Ridge plots
  TString ridge_vars = corr_rap->GetVar(1)+"_"+corr_rap->GetVar(2);
  
  TString ridge_corr1 = gDirName+"plots_ridge/";
  TString ridge_corr2 = gMode+"_"+target+"_"+ridge_vars+"_corr.png";
  ridge_plot(corr_rap->GetCO(),ridge_corr1,ridge_corr2);

  TString ridge_mult1 = gDirName+"plots_ridge/";
  TString ridge_mult2 = gMode+"_"+target+"_"+ridge_vars+"_mult.png";
  ridge_plot(corr_rap->GetME(),ridge_mult1,ridge_mult2);
  
  TString ridge_same1 = gDirName+"plots_ridge/";
  TString ridge_same2 = gMode+"_"+target+"_"+ridge_vars+"_same.png";
  ridge_plot(corr_rap->GetSE(),ridge_same1,ridge_same2);
  

  TCanvas *cinv = new TCanvas();
  h_inv_mass.Draw();
  cinv->SaveAs(gDirName+"inv_mass_"+gMode+"_"+target+".png");
  /*
  TH2F * reco1 = new TH2F("reco1","test",10,0,10,10,0,10);
  TH2F * reco2 = new TH2F("reco2","test",10,0,10,10,0,10);
  TH2F * reco3 = new TH2F("reco3","test",10,0,10,10,0,10);
  *reco1 = corr_rap->GetReco();
  *reco2 = corr_rap->GetRecoTriggers();
  *reco3 = corr_rap->GetRecoPartners();

  new TCanvas();
  reco1->Draw("colz");
  new TCanvas();
  reco2->Draw("colz");
  new TCanvas();
  reco3->Draw("colz");

  TString con = "1D_reco_pm_";
  export_hist(reco1, gDirName+"/plots_reco/"+con+target+"_"+gMode+".png");
  con = "1D_reco_p_";
  export_hist(reco2, gDirName+"/plots_reco/"+con+target+"_"+gMode+".png");
  con = "1D_reco_m_";
  export_hist(reco3, gDirName+"/plots_reco/"+con+target+"_"+gMode+".png");
  */

  
  // Plotting
  /*
  full_ridge_plots(corr_boo);
  full_1d_plots(corr_boo);
  full_ridge_plots(corr_rap);
  full_1d_plots(corr_rap);
  */

}

void export_hist(TH1 * h2, TString out_filename) {
  auto c1 = new TCanvas();
  c1->SetCanvasSize(800,600);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.15);
  h2->Draw("colz");

  c1->SaveAs(out_filename);
  c1->Modified();
  c1->Update();
  delete c1;
}

void recoplots(TCorrelation * t, TString target, TString gMode, TString var) {
  TH2F * reco1 = new TH2F("reco1","test",10,0,10,10,0,10);
  TH2F * reco2 = new TH2F("reco2","test",10,0,10,10,0,10);
  TH2F * reco3 = new TH2F("reco3","test",10,0,10,10,0,10);
  *reco1 = t->GetReco();
  *reco2 = t->GetRecoTriggers();
  *reco3 = t->GetRecoPartners();

  new TCanvas();
  reco1->Draw("colz");
  new TCanvas();
  reco2->Draw("colz");
  new TCanvas();
  reco3->Draw("colz");

  TString con = "1D_reco_pm_";
  export_hist(reco1, gDirName+"/plots_reco/"+con+target+"_"+gMode+"_"+var+".png");
  con = "1D_reco_p_";
  export_hist(reco2, gDirName+"/plots_reco/"+con+target+"_"+gMode+"_"+var+".png");
  con = "1D_reco_m_";
  export_hist(reco3, gDirName+"/plots_reco/"+con+target+"_"+gMode+"_"+var+".png");
  
}
