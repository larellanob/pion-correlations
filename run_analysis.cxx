#include "Root/TCorrelation.cxx"
#include "Root/Plotter.cc"
#include "Root/DeltaAngle.cxx"

void export_hist(TH1 * h2, TString out_filename);
void recoplots(TCorrelation *t, TString target, TString gMode, TString var);
TString gDirName = "/home/luciano/Physics/CLAS/pion_correlation/";
bool m_simulation = false;
bool gMirror = true;
Float_t gDataCap  = 1.0; // fraction of data to be used in analysis (1.0 == full data)

//void run_analysis( TString filename = "/home/luciano/Physics/CLAS/pion_ridge/data/2pcPairs_D_zhpm.root" )

void run_analysis( TString centroid_s = "",TString target = "Pb", TString gMode = "pt3")
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



  std::vector<std::vector<Float_t>> old_events_phi;
  std::vector<std::vector<Float_t>> old_aux_phi;

  std::vector<std::vector<Float_t >> old_events_the;
  std::vector<std::vector<Float_t >> old_aux_the;

  std::vector<std::vector<Float_t >> old_events_rap;
  std::vector<std::vector<Float_t >> old_aux_rap;

  std::vector<int> old_sizes_vector;
  std::vector<int> old_aux_sizes;

  TH2F h2mix("h2mix","",16,0,180,8,0,2.4);
  bool high_mult = true;
  UInt_t pairs = 0;


  // zh bins biases
  double zbias = 0;
  if ( !m_simulation && gMirror ) {
    if ( target == "D" && gMode == "zh1") {
      zbias = 1.521;
    } else if ( target == "D" && gMode == "zh2") {
      zbias = 0.9432;
    } else if ( target == "D" && gMode == "zh3") {
      zbias = 0.4351;
    } else if ( target == "C" && gMode == "zh1") {
      zbias = 1.661;
    } else if ( target == "C" && gMode == "zh2") {
      zbias = 0.9911;
    } else if ( target == "C" && gMode == "zh3") {
      zbias = 0.4347;
    } else if ( target == "Fe" && gMode == "zh1") {
      zbias = 1.507;
    } else if ( target == "Fe" && gMode == "zh2") {
      zbias = 1.009;
    } else if ( target == "Fe" && gMode == "zh3") {
      zbias = 0.4779;
    } else if ( target == "Pb" && gMode == "zh1") {
      zbias = 1.415;
    } else if ( target == "Pb" && gMode == "zh2") {
      zbias = 0.9739;
    } else if ( target == "Pb" && gMode == "zh3") {
      zbias = 0.4943;
    } 
  }  else if ( m_simulation && gMirror ) {
    if ( target == "D" && gMode == "zh1") {
      zbias = 1.577;
    } else if ( target == "D" && gMode == "zh2") {
      zbias = 0.8392;
    } else if ( target == "D" && gMode == "zh3") {
      zbias = 0.4448;
    } else if ( target == "C" && gMode == "zh1") {
      zbias = 1.573;
    } else if ( target == "C" && gMode == "zh2") {
      zbias = 0.8445;
    } else if ( target == "C" && gMode == "zh3") {
      zbias = 0.4491;
    } else if ( target == "Fe" && gMode == "zh1") {
      zbias = 1.559;
    } else if ( target == "Fe" && gMode == "zh2") {
      zbias = 0.9012;
    } else if ( target == "Fe" && gMode == "zh3") {
      zbias = 0.4577;
    } else if ( target == "Pb" && gMode == "zh1") {
      zbias = 1.543;
    } else if ( target == "Pb" && gMode == "zh2") {
      zbias = 0.8243;
    } else if ( target == "Pb" && gMode == "zh3") {
      zbias = 0.4371;
    } 
  }
  
  while ( triggers.Next()  && events < (gDataCap*entries) ) {
    events++;
    t_size = pid_t.GetSize();

    
    // you can get "old partners" by working here, before partners.Next()
    
    /// Mixed Event Pairs
    if ( centroid == 0 ) {
      rho_warning= true;
    }

    if ( events > 1 && rho_warning && high_mult ) {
      std::vector<Float_t> valuephi;
      std::vector<Float_t> valuethe;
      std::vector<Float_t> valuerap;
      for ( int o = 0; o < p_size; o++ ) {
	valuephi.push_back(PhiBoo_p[o]);
	valuethe.push_back(TheBoo_p[o]);
	valuerap.push_back(RapBoo_p[o]);
      }
      
      if ( old_events_phi.size() < 10 ) {

	old_events_phi.push_back(valuephi);
	old_events_the.push_back(valuethe);
	old_events_rap.push_back(valuerap);
	old_sizes_vector.push_back(p_size);
      } else {
	old_aux_phi = old_events_phi;
	old_aux_the = old_events_the;
	old_aux_rap = old_events_rap;
	old_aux_sizes = old_sizes_vector;
	old_events_phi.clear();
	old_events_the.clear();
	old_events_rap.clear();
	old_sizes_vector.clear();
	for ( int i = 0; i < 9; i++ ) {
	  old_events_phi.push_back(old_aux_phi[i+1]);
	  old_events_the.push_back(old_aux_the[i+1]);
	  old_events_rap.push_back(old_aux_rap[i+1]);
	  old_sizes_vector.push_back(old_aux_sizes[i+1]);
	}
	old_events_phi.push_back(valuephi);
	old_events_the.push_back(valuethe);
	old_events_rap.push_back(valuerap);
	old_sizes_vector.push_back(p_size);
      }
      
    }
    
    /*    
    if ( events > 1  && rho_warning ) {
      for ( int t = 0; t < t_size; t++ ) {
	for ( int o = 0; o < p_size; o++ ) {
	  Float_t DPhi,DThe,DRap;
	  DPhi = DeltaAngle(PhiBoo_t[t],PhiBoo_p[o]);
	  DThe = DeltaAngle(TheBoo_t[t],TheBoo_p[o]);
	  DRap = abs(RapBoo_t[t]-RapBoo_p[o]);
	  corr_boo->FillMulti(DPhi,DThe);
	  corr_rap->FillMulti(DPhi,DRap);
	}
      }
    }
    */

    h2mix.Reset();
    if ( events > 1  && rho_warning && high_mult ) {
      for ( int t = 0; t < t_size; t++ ) {
	for ( int v = 0; v < old_events_phi.size(); v++ ) {
	  int old_size = old_sizes_vector[v];
	  for ( int o = 0; o < old_size; o++ ) {
	    Float_t DPhi,DThe,DRap;
	    DPhi = DeltaAngle(PhiBoo_t[t],old_events_phi[v][o]);
	    DThe = DeltaAngle(TheBoo_t[t],old_events_the[v][o]);
	    DRap = RapBoo_t[t]-old_events_rap[v][o];
	    corr_boo->FillMulti(DPhi,DThe);
	    if ( !gMirror ) {
	      corr_rap->FillMulti(DPhi,DRap-zbias);
	    } else if ( gMirror ) {
	      corr_rap->FillMulti(DPhi,abs(DRap-zbias));
	    }
	    h2mix.Fill(DPhi,DRap);
	  }
	}
      }
    }


    partners.Next();
    int old_p_size = p_size;
    p_size = pid_p.GetSize();

    // high multiplicity
    /*
    if ( p_size < 2 ) {
      high_mult = false;
      continue;
    } else {
      high_mult = true;
    }
    */
    pairs++;
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
	if ( !gMirror ) {
	  corr_rap->FillSame(DPhi,DRap-zbias);
	} else if ( gMirror ) {
	  corr_rap->FillSame(DPhi,abs(DRap-zbias));
	}
	double instant_mixed = 0;
	instant_mixed = h2mix.GetBinContent(h2mix.FindBin(DPhi,DRap));
	double instant_mixed0 = 0;
	instant_mixed = h2mix.GetBinContent(h2mix.FindBin(0,0));
	if ( instant_mixed != 0 ) {
	  //double filler = instant_mixed0/instant_mixed;
	  double filler = 1.0/instant_mixed;
	  corr_rap->FillInstantRidge(DPhi,DRap,filler);
	} else {
	  corr_rap->FillInstantRidge(DPhi,DRap);
	}
      }
      fill_neg = false;
    }
    

    
  }
  cout << "rhos " << rhos << endl;
  cout << "from total events " << events << endl;
  cout << "found this many pairs " << pairs << endl;
  corr_boo->FillCorrelation();
  corr_rap->FillCorrelation();

  TFile * fout = new TFile(gDirName+"histograms/hist_"+target+"_"+gMode+"_"+centroid_s+"_boo.root","recreate");
  corr_boo->GetCorr1D(1).Write();
  corr_boo->GetCorr1D(2).Write();
  corr_boo->GetSame1D(1).Write();
  corr_boo->GetSame1D(2).Write();
  corr_boo->GetMult1D(1).Write();
  corr_boo->GetMult1D(2).Write();
  corr_boo->GetCorrMirror().Write();
  corr_boo->GetMultMirror().Write();
  corr_boo->GetSameMirror().Write();
  corr_boo->GetCorrMirrorphi().Write();
  corr_boo->GetMultMirrorphi().Write();
  corr_boo->GetSameMirrorphi().Write();
  
  TFile * fout2 = new TFile(gDirName+"histograms/hist_"+target+"_"+gMode+"_"+centroid_s+"_rap.root","recreate");
  corr_rap->GetCorr1D(1).Write();
  corr_rap->GetCorr1D(2).Write();
  corr_rap->GetSame1D(1).Write();
  corr_rap->GetSame1D(2).Write();
  corr_rap->GetMult1D(1).Write();
  corr_rap->GetMult1D(2).Write();
  corr_rap->GetCorrMirror().Write();
  corr_rap->GetMultMirror().Write();
  corr_rap->GetSameMirror().Write();
  corr_boo->GetCorrMirrorphi().Write();
  corr_boo->GetMultMirrorphi().Write();
  corr_boo->GetSameMirrorphi().Write();
  
  // Reconstructed pions plots
  // initialization of arbitrary TH2Fs to get from corr_*
  recoplots(corr_boo,target, gMode, "boo");
  recoplots(corr_rap,target, gMode, "rap");



  // Ridge plots
  TString ridge_vars = corr_rap->GetVar(1)+"_"+corr_rap->GetVar(2);
  TString ridge_corr1 = gDirName+"plots_ridge/";
  TString ridge_corr2 = gMode+"_"+target+"_"+ridge_vars+"_corr.png";
  TString ridge_mult1 = gDirName+"plots_ridge/";
  TString ridge_mult2 = gMode+"_"+target+"_"+ridge_vars+"_mult.png";
  TString ridge_same1 = gDirName+"plots_ridge/";
  TString ridge_same2 = gMode+"_"+target+"_"+ridge_vars+"_same.png";
  // data
  if ( !m_simulation ) {
    //ridge_plot(corr_rap->GetCorr(),ridge_corr1,ridge_corr2,target,gMode);
    if ( !gMirror ) {
      ridge_plot(corr_rap->GetCO(),ridge_corr1,ridge_corr2,target,gMode);
      ridge_plot(corr_rap->GetME(),ridge_mult1,ridge_mult2,target,gMode);
      ridge_plot(corr_rap->GetSE(),ridge_same1,ridge_same2,target,gMode);
    } else if ( gMirror ) {
      ridge_plot(corr_rap->GetCorrMirror(),ridge_corr1,ridge_corr2,target,gMode);
      ridge_plot(corr_rap->GetMultMirror(),ridge_mult1,ridge_mult2,target,gMode);
      ridge_plot(corr_rap->GetSameMirror(),ridge_same1,ridge_same2,target,gMode);
    }
  } else if ( m_simulation ) {
    //ridge_plot(corr_rap->GetCorr(),ridge_corr1,ridge_corr2,target,gMode);
    if ( !gMirror ) {
      ridge_plot(corr_rap->GetCO(),ridge_corr1,ridge_corr2,target,gMode,true);
      ridge_plot(corr_rap->GetME(),ridge_mult1,ridge_mult2,target,gMode,true);
      ridge_plot(corr_rap->GetSE(),ridge_same1,ridge_same2,target,gMode,true);
    } else if ( gMirror ) {
      ridge_plot(corr_rap->GetCorrMirror(),ridge_corr1,ridge_corr2,target,gMode,true);
      ridge_plot(corr_rap->GetMultMirror(),ridge_mult1,ridge_mult2,target,gMode,true);
      ridge_plot(corr_rap->GetSameMirror(),ridge_same1,ridge_same2,target,gMode,true);
    }
  }
  corr_rap->ScaleInsta(1./pairs);

  /*
  TString ridge_inst1 = gDirName+"plots_ridge/";
  TString ridge_inst2 = gMode+"_"+target+"_"+ridge_vars+"_inst.png";
  ridge_plot(corr_rap->GetInst(),ridge_inst1,ridge_inst2,target,gMode);
  */

  TCanvas *cinv = new TCanvas();
  h_inv_mass.Draw();
  cinv->SaveAs(gDirName+"plots_invmass/inv_mass_"+gMode+"_"+target+".png");
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
