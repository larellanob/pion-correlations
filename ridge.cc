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
#include "Root/TCorrelation.cxx"

TString out_path = "/home/luciano/Physics/CLAS/pion_ridge/";
bool m_debug      = false;
bool m_simulation = false;
bool old_triggers     = false; // if false, use old_partners
bool DMode        = false; // deuterium
Float_t gDataCap  = 0.2; // fraction of data to be used in analysis (1.0 == full data)
TDatabasePDG db;
// modes
// ++: trigger pi+ no cuts, assoc pi- no cuts
// +-: trigger pi+ no cuts, assoc pi+ no cuts
// zh+-: trigger pi+ zh > 0.5, assocc pi- zh < 0.5
// zh++: trigger pi+ zh > 0.5, assocc pi+ zh < 0.5
// zhlluu+-: trigger pi+ zh>0.5, assoc pi- zh in (ll,uu)
// zhlluu++: trigger pi+ zh>0.5, assoc pi+ zh in (ll,uu)
TString gMode = "+-"; // choose +- or ++


//Float_t PhiPQ(TVector3, TLorentzVector);
void export_hist(TH2 *, TString, TString options = "colz");
void export_hist(TH1 * h2, TString out_filename);
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
  if ( target == "" ) {
    cout << "No target input, using C as default" << endl;
    target = "C";
  }
  //gStyle->SetOptStat(0);
  if ( m_simulation == true ) {
    out_path+="simulation/";
  } else {
    out_path += "data/";
  }

  // Chain for hadrons (ch) and chain for triggers electrons (ech)
  TChain ch("tree_data");
  TChain ech("e_rec");

  // Add files to the chains depending on target
  // Also sets beam energy
  Double_t kEbeam;
  if ( m_simulation ) {
    ch.Add("/home/luciano/Physics/CLAS/data/tree_output.root");
    ech.Add("/home/luciano/Physics/CLAS/data/tree_output.root");
    kEbeam = 11.0;
  } else {
    if ( target == "Pb" ) {
      ch.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
      ech.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    } else if ( target == "Fe" ) {
      ch.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
      ech.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    } else if ( target == "C" ) {
      ch.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
      ech.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
    }

    kEbeam = 5.014;
  }
  if ( target == "D" ) {
    DMode = true;
    ch.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    ch.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    ch.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_Pb_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_Fe_files.root");
    ech.Add("/home/luciano/Physics/CLAS/data/full_C_files.root");
  }

  // Known particle masses
  Double_t kChargedPionMass = db.GetParticle(211)->Mass();
  Double_t kProtonMass = db.GetParticle(2212)->Mass();
  Double_t kRhoMass = db.GetParticle(113)->Mass();

  // Full chain entries
  Int_t Entries = 0;
  Entries = ch.GetEntries();
  
  // Tree reader for hadrons and scattered electron
  TTreeReader th(&ch);
  TTreeReader e(&ech);

  // Variables to be read from trees
  // hadrons
  TTreeReaderArray<Int_t> pid(th,"pid");
  TTreeReaderValue<Float_t> q2(th,"Q2");
  TTreeReaderArray<Float_t> P(th,"P");
  TTreeReaderArray<Float_t> px(th,"Px");
  TTreeReaderArray<Float_t> py(th,"Py");
  TTreeReaderArray<Float_t> pz(th,"Pz");

  TTreeReaderValue<Int_t> TargType(th,"TargType");
  TTreeReaderArray<Float_t> Zh(th,"Zh");

  TTreeReaderArray<Float_t> Theta(th,"Theta");
  TTreeReaderArray<Float_t> Phi(th,"Phi");
  TTreeReaderArray<Float_t> ThetaPQ(th,"ThetaPQ");
  TTreeReaderArray<Float_t> PhiPQ(th,"PhiPQ");

  // electrons
  TTreeReaderValue<Float_t> epx(e,"Pex");
  TTreeReaderValue<Float_t> epy(e,"Pey");
  TTreeReaderValue<Float_t> epz(e,"Pez");
  TTreeReaderValue<Float_t> eq2(e,"Q2");
  TTreeReaderValue<Float_t> ew(e,"W");
    
  cout << "working" << endl;


  // Vectors to be filled with triggers and associated particles
  std::vector<int> triggers_index;
  std::vector<int> partners_index;

  std::vector<TLorentzVector> triggers4v;
  std::vector<TLorentzVector> partners4v;

  std::vector<TLorentzVector> triggers4v_rotated;
  std::vector<TLorentzVector> partners4v_rotated;

  std::vector<TLorentzVector> triggers4v_boosted;
  std::vector<TLorentzVector> partners4v_boosted;

  std::vector<TLorentzVector> old4v;
  std::vector<TLorentzVector> old4v_rotated;
  std::vector<TLorentzVector> old4v_boosted;

  TH1 * mom_sum = new TH1F("mom_sum","lab frame;P;entries",100,0,20);
  TH1 * mom_sum2 = new TH1F("mom_sum2","transformed;P;entries",100,0,20);

  TH1 * energy_sum = new TH1F("energy_sum","lab frame;E;entries",100,0,20);
  TH1 * energy_sum2 = new TH1F("energy_sum2","transformed;E;entries",100,0,20);
  // Setting some counters
  UInt_t events = 0;
  int processed_events=0;
    
  TLorentzVector virtual_photon;
  //TLorentzVector oldgamma;


  // TCorrelation objects initialization
  // For ridge analysis
  // angular (lab frame)  correlation
  TCorrelation * corr_ang = new TCorrelation("phi","theta",target);
  // angular (pq frame) correlation
  TCorrelation * corr_apq = new TCorrelation("phiPQ","thetaPQ", target);
  // angular (virtual photon frame) correlation
  TCorrelation * corr_boo = new TCorrelation("phiPQboosted","thetaPQboosted",target);
  //TCorrelation * corr_eta = new TCorrelation("eta","theta", target);
  //TCorrelation * corr_vir = new TCorrelation("vphi","vtheta", target);

  // Number of trigger particles (for normalization)
  ULong_t NTriggers = 0;
  
  
  // Color entanglement
  TH2 * h15
    = new TH2F("h15",
	       "3vector momentum diff #ne ("+gMode+") Lab frame;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);
  TH2 * h16
    = new TH2F("h16",
	       "3vector momentum sum ("+gMode+") Lab frame;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);

  TH2 * h151
    = new TH2F("h151",
	       "3vector momentum diff ("+gMode+") Rotated;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);
  TH2 * h161
    = new TH2F("h161",
	       "3vector momentum sum ("+gMode+") Rotated;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);

  TH2 * h152
    = new TH2F("h152",
	       "3vector momentum diff ("+gMode+") Transformed;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);
  TH2 * h162
    = new TH2F("h162",
	       "3vector momentum sum ("+gMode+") Transformed;Q^{2};P_{diff};counts",
	       35, 0.5,4 ,30,0,3);
  
  
  if ( m_debug )
    cout << "DEBUG: entering main tree reading loop" << endl;
  // main file/tree reading loop

  ///////////////////////////
  //// LOOP ////////////
  ///////////////////////////
  while ( th.Next()  && events < (gDataCap * Entries) ) {
    e.Next();
    events++;
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

    ///////////////////////////
    // Zh > 1 veto 
    double zh_sum = 0;
    for ( int i = 0; i < Zh.GetSize(); i++ ) {
      if ( pid[i] == 211 || pid[i] == -211 ) {
	zh_sum += Zh[i];
      }
    }
    if ( zh_sum > 1.0 ) {
      continue;
    }
    bool impossible_momentum = false;
    for ( int i = 0; i < P.GetSize(); i++ ) {
      if ( P[i] > 6.0 ) {
	impossible_momentum = true;
	break;
      }
    }
    if ( impossible_momentum )
      continue;
    
    
    
    // Electron energy
    Double_t scattered_energy = sqrt((*epx)*(*epx)+(*epy)*(*epy)+(*epz)*(*epz));
    if ( scattered_energy < 0.005 )
      cout << "WARNING: Low e' energy " << scattered_energy << endl; 
    if ( scattered_energy > kEbeam  )
      cout << "WARNING: High e' energy = " <<  scattered_energy << endl;

    // Virtual photon frame is per event
    // oldgamma = virtual_photon; // true old gamma

    // Virtual photon 4vector
    virtual_photon = {-(*epx),-(*epy),kEbeam-(*epz),kEbeam-scattered_energy};
    //oldgamma = virtual_photon; // actually the same old gamma

    //////////////////////////
    // Pions in the event
    triggers_index.clear();
    partners_index.clear();
    
    // modes
    // +++: trigger pi+ with highest zh, partner all other pi+
    // ++: trigger pi+ no cuts, assoc pi- no cuts
    // +-: trigger pi+ no cuts, assoc pi+ no cuts
    // zh+-: trigger pi+ zh > 0.5, assocc pi- zh < 0.5
    // zh++: trigger pi+ zh > 0.5, assocc pi+ zh < 0.5
    // zhlluu+-: trigger pi+ zh>0.5, assoc pi- zh in (ll,uu)
    // zhlluu++: trigger pi+ zh>0.5, assoc pi+ zh in (ll,uu)
    // high-low zh 
    if ( gMode == "+++" ) {
      int zh_max = 0;
      // first pi+ found is stored at triggers_index[0]
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 ) {
	  if ( Zh[i] > zh_max ) {
	    if ( triggers_index.size() > 0 ) {
	      partners_index.push_back(triggers_index[0]);
	    }
	    triggers_index.clear();
	    triggers_index.push_back(i);
	  } else {
	    partners_index.push_back(i);
	  }
	}
      }
    } else if ( gMode == "++" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 ) {
	  triggers_index.push_back(i);
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "+-" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh+-" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh++" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == +211 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0305++" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.3 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0305+-" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.3 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0203++" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.2 && Zh[i] < 0.3 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0203+-" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.2 && Zh[i] < 0.3 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0002++" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.0 && Zh[i] < 0.2 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0002+-" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.0 && Zh[i] < 0.2 ) {
	  partners_index.push_back(i);
	}
      }
    }

    
    // make sure of this, normalization might depend on it
    if ( triggers_index.size() == 0  || partners_index.size() == 0 )
      continue;


    double p_sum = 0;
    double E_sum = 0;
    for ( int i = 0; i < Zh.GetSize(); i++ ) {
      if ( pid[i] == 211 || pid[i] == -211 ) {
	p_sum += P[i];
	E_sum+= sqrt(P[i]*P[i] + kChargedPionMass*kChargedPionMass);
      }
    }
    if ( p_sum > 6.0 || E_sum > 6.0 ) {
      cout << "exiting gracefully " << endl;
      continue;
    }
    mom_sum->Fill(p_sum);
    energy_sum->Fill(E_sum);
    
    ///////////////////////////
    // Pions 4vectors
    triggers4v.clear();
    partners4v.clear();
    triggers4v_rotated.clear();
    partners4v_rotated.clear();
    triggers4v_boosted.clear();
    partners4v_boosted.clear();

    for ( int i:  triggers_index ) {
      Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
      TLorentzVector pi (px[i],py[i],pz[i],e);
      triggers4v.push_back(pi);
    }

    for ( int i: partners_index ) {
      Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
      TLorentzVector pi (px[i],py[i],pz[i],e);
      partners4v.push_back(pi);
    }
    
    /*
    for ( int i = 0; i < pid.GetSize(); i++ ) {
      if ( pid[i] == 211 ) {
	Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
	TLorentzVector pi (px[i],py[i],pz[i],e);
	triggers4v.push_back(pi);
      } else if ( pid[i] == -211 ) {
	Float_t e = px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass;
	TLorentzVector pi (px[i],py[i],pz[i],e);
	partners4v.push_back(pi);
      }
    }
    */
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
    for ( int i = 0; i < triggers4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(triggers4v[i],angle1,angle2,0);
      triggers4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(triggers4v[i],angle1,angle2,-boost);
      triggers4v_boosted.push_back(aux_boost);
    }

    for ( int i = 0; i < partners4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(partners4v[i],angle1,angle2,0);
      partners4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(partners4v[i],angle1,angle2,-boost);
      partners4v_boosted.push_back(aux_boost);
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
      cout << "pospions: " << triggers_index.size() << endl;
      for ( int i = 0; i < triggers_index.size(); i++ ) {
	cout << "pz: " << pz[triggers_index[i]] << endl;
	cout << "theta: " << Theta[triggers_index[i]] << endl;
	cout << "phi: " << Phi[triggers_index[i]] << endl;
	cout << "thetapq: " << ThetaPQ[triggers_index[i]] << endl;
	cout << "phipq: " << PhiPQ[triggers_index[i]] << endl;
	cout << "-------------" << endl;
      }
      
      cout << "4v: " << triggers4v.size() << endl;
      for ( int i = 0; i < triggers4v.size(); i++ ) {
	cout << "-----------------------------" << endl;
	cout << "pz: " << triggers4v[i].Pz() << endl;
	cout << "theta: " << triggers4v[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi: " << triggers4v[i].Phi()*180./TMath::Pi() << endl;
	cout << "-------------" << endl;
	cout << "pz rot: " << triggers4v_rotated[i].Pz() << endl;
	cout << "theta rot: " << triggers4v_rotated[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi rot: " << triggers4v_rotated[i].Phi()*180./TMath::Pi() << endl;
	cout << "-------------" << endl;
	cout << "pz boo: " << triggers4v_boosted[i].Pz() << endl;
	cout << "theta boo: " << triggers4v_boosted[i].Theta()*180./TMath::Pi() << endl;
	cout << "phi boo: " << triggers4v_boosted[i].Phi()*180./TMath::Pi()<< endl;
      }
    }
    
    
    if ( m_debug )
      cout << "DEBUG: Entering mutli event filling" << endl;


    ///////////////////////////
    // multi event filling
    // all possible combinations, i.e. all pi+ are triggers in +- mode
    if ( old4v.size() != 0 && !old_triggers && events > 1 ) {
      for ( int i = 0; i < triggers4v.size(); i++ ) {
	for ( int j = 0; j < old4v.size(); j++ ) {
	  Float_t x,y;

	  x = DeltaAngleRad(triggers4v[i].Phi(),old4v[j].Phi());
	  y = DeltaAngleRad(triggers4v[i].Theta(),old4v[j].Theta());
	  corr_ang->FillMulti(x,y);

	  x = DeltaAngleRad(triggers4v_rotated[i].Phi(),old4v_rotated[j].Phi());
	  y = DeltaAngleRad(triggers4v_rotated[i].Theta(),old4v_rotated[j].Theta());
	  corr_apq->FillMulti(x,y);
	  
	  x = DeltaAngleRad(triggers4v_boosted[i].Phi(),old4v_boosted[j].Phi());
	  y = DeltaAngleRad(triggers4v_boosted[i].Theta(),old4v_boosted[j].Theta());
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
    
    ////////////////////////////////////
    /// SAME EVENT FILLING

    for ( int i = 0; i < triggers4v.size(); i++ ) {
      // Reconstructed pions saving
      corr_ang->FillReco(triggers4v[i].Phi()*180./TMath::Pi(), triggers4v[i].Theta()*180./TMath::Pi());
      corr_apq->FillReco(triggers4v_rotated[i].Phi()*180./TMath::Pi(), triggers4v_rotated[i].Theta()*180./TMath::Pi() );
      corr_boo->FillReco(triggers4v_boosted[i].Phi()*180./TMath::Pi(), triggers4v_boosted[i].Theta()*180./TMath::Pi() );

      corr_ang->FillRecoTriggers(triggers4v[i].Phi()*180./TMath::Pi(), triggers4v[i].Theta()*180./TMath::Pi());
      corr_apq->FillRecoTriggers(triggers4v_rotated[i].Phi()*180./TMath::Pi(), triggers4v_rotated[i].Theta()*180./TMath::Pi() );
      corr_boo->FillRecoTriggers(triggers4v_boosted[i].Phi()*180./TMath::Pi(), triggers4v_boosted[i].Theta()*180./TMath::Pi() );
      for ( int j = 0; j < partners4v.size(); j++ ) {
	// removes correlation to itself
	if ( partners4v[j] == triggers4v[i] ) {
	  continue;
	}
	float x,y;

	// Correlations filling
	x = DeltaAngleRad(triggers4v[i].Phi(),partners4v[j].Phi());
	y = DeltaAngleRad(triggers4v[i].Theta(),partners4v[j].Theta());
	corr_ang->FillSame(x,y);

	x = DeltaAngleRad(triggers4v_rotated[i].Phi(),partners4v_rotated[j].Phi());
	y = DeltaAngleRad(triggers4v_rotated[i].Theta(),partners4v_rotated[j].Theta());	
	corr_apq->FillSame(x,y);
	
	x = DeltaAngleRad(triggers4v_boosted[i].Phi(),partners4v_boosted[j].Phi());
	y = DeltaAngleRad(triggers4v_boosted[i].Theta(),partners4v_boosted[j].Theta());
	corr_boo->FillSame(x,y);

	// Should enter this loop only once per event (when i == 0 and loop once over j)
	if ( fill_neg ) {
	  // Reconstructed pions saving
	  corr_ang->FillReco(triggers4v[j].Phi()*180./TMath::Pi(), triggers4v[j].Theta()*180./TMath::Pi());
	  corr_apq->FillReco(triggers4v_rotated[j].Phi()*180./TMath::Pi(), triggers4v_rotated[j].Theta()*180./TMath::Pi() );
	  corr_boo->FillReco(triggers4v_boosted[j].Phi()*180./TMath::Pi(), triggers4v_boosted[j].Theta()*180./TMath::Pi() );

	  corr_ang->FillRecoPartners(triggers4v[j].Phi()*180./TMath::Pi(), triggers4v[j].Theta()*180./TMath::Pi());
	  corr_apq->FillRecoPartners(triggers4v_rotated[j].Phi()*180./TMath::Pi(), triggers4v_rotated[j].Theta()*180./TMath::Pi() );
	  corr_boo->FillRecoPartners(triggers4v_boosted[j].Phi()*180./TMath::Pi(), triggers4v_boosted[j].Theta()*180./TMath::Pi() );
	}
      }
      fill_neg = false;
    }

    //////////////////////////////////////
    // Old pion saving

    // Old negative pion:
    for ( int j = 0; j < partners4v.size(); j++ ) {
      old4v.push_back(partners4v[j]);
      //old4v_rotated.push_back(partners4v_rotated[j]);
      //old4v_boosted.push_back(partners4v_boosted[j]);
    }


    // Color correlation filling
    for ( TLorentzVector i : triggers4v ) {
      for ( TLorentzVector j: partners4v ) {
	if ( i == j ) continue;
	h15->Fill(*eq2,(i-j).P());
	h16->Fill(*eq2,(i+j).P());
      }
    }
    for ( TLorentzVector i : triggers4v_rotated ) {
      for ( TLorentzVector j: partners4v_rotated ) {
	if ( i == j ) continue;
	h151->Fill(*eq2,(i-j).P());
	h161->Fill(*eq2,(i+j).P());
      }
    }
    double p_sum2 = 0;
    double E_sum2 = 0;
    bool parn = true;
    for ( TLorentzVector i : triggers4v_boosted ) {
      p_sum2 += i.P();
      E_sum2 += sqrt(i.P()*i.P() + kChargedPionMass*kChargedPionMass);
      
      for ( TLorentzVector j: partners4v_boosted ) {
	if ( i == j ) continue;
	if ( parn ) {
	  p_sum2 += j.P();
	  E_sum2 += sqrt(j.P()*j.P() + kChargedPionMass*kChargedPionMass);
	}
	parn = false;
	h152->Fill(*eq2,(i-j).P());
	h162->Fill(*eq2,(i+j).P());
      }
    }
    mom_sum2->Fill(p_sum2);
    energy_sum2->Fill(E_sum2);
    NTriggers += triggers_index.size();
    //    if ( triggers_index.size() < 1 )
    //cout << "what " << triggers_index.size() << endl;
    processed_events++;
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << processed_events << endl;
  cout << "number of trigger particles: " << NTriggers << endl;
  //side_by_side(ns,ns_vir,"NS_lab_virtual.png");

  TString multi_method;
  if ( old_triggers )
    multi_method = "_pastTriggers";
  else if ( !old_triggers )
    multi_method = "_pastPartners";

  //corr_apq->NormalizeSame(NTriggers);
  //corr_apq->NormalizeMulti();

  // FillCorrelation generates the "ridge like" 2D histograms
  corr_ang->FillCorrelation();
  corr_apq->FillCorrelation();
  corr_boo->FillCorrelation();

  // This draws the 2D plots (saves to disk)
  full_ridge_plots(corr_ang);
  full_ridge_plots(corr_apq);
  full_ridge_plots(corr_boo);

  // This draws the 1D plots (saves to disk)
  full_1d_plots(corr_ang);
  full_1d_plots(corr_apq);
  full_1d_plots(corr_boo);

  // Reconstructed pions plots
  // initialization of arbitrary TH2Fs to get from corr_*
  TH2F * reco1 = new TH2F("reco1","test",10,0,10,10,0,10);
  TH2F * reco2 = new TH2F("reco2","test",10,0,10,10,0,10);
  TH2F * reco3 = new TH2F("reco3","test",10,0,10,10,0,10);
  *reco1 = corr_apq->GetReco();
  *reco2 = corr_apq->GetRecoTriggers();
  *reco3 = corr_apq->GetRecoPartners();

  new TCanvas();
  reco1->Draw("colz");
  new TCanvas();
  reco2->Draw("colz");
  new TCanvas();
  reco3->Draw("colz");

  TString con = "1D_reco_pm_";
  export_hist(reco1, out_path+con+".png");
  con = "1D_reco_p_";
  export_hist(reco2, out_path+con+".png");
  con = "1D_reco_m_";
  export_hist(reco3, out_path+con+".png");
  /*
  new TCanvas();
  corr_apq->GetReco().Draw("colz");
  new TCanvas();
  corr_apq->GetRecoTriggers().Draw("colz");
  new TCanvas();
  corr_apq->GetRecoPartners().Draw("colz");
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


  TFile *fout1 = new TFile("/home/luciano/Physics/CLAS/pion_correlation/"+gMode+"_"+target+".root","recreate");
  corr_boo->GetCorr1D(1).Write();
  corr_boo->GetCorr1D(2).Write();

  // Color entanglement plotting
  new TCanvas();
  h152->Draw("colz");
  new TCanvas();
  h162->Draw("colz");
  new TCanvas();
  h151->Draw("colz");
  new TCanvas();
  h161->Draw("colz");
  new TCanvas();
  h15->Draw("colz");
  new TCanvas();
  h16->Draw("colz");
  new TCanvas();
  mom_sum->Draw();
  new TCanvas();
  mom_sum2->Draw();
  new TCanvas();
  energy_sum->Draw();
  new TCanvas();
  energy_sum2->Draw();


  TString outdir_color = "/home/luciano/Physics/CLAS/color_entanglement/";
  
  for ( int j = 0; j < 30; j+=5 ){
    double low = (3.0/30.0)*(j+0.0);
    double high = (3.0/30.0)*(j+4.0);
    TString title = Form("%.1f #leq P_{diff} [GeV] #leq %.1f",low, high);
    TH1*  hproj = h152->ProjectionX(title,j, j+4);
    hproj->SetFillColor(kBlue);
    TString export_name = Form("bins_%ito%i",j,j+4);
    TString export_path = outdir_color+gMode+"_dif_"+export_name+".png";
    export_hist(hproj,export_path);
    delete hproj;
  }
  for ( int j = 0; j < 30; j+=5 ){
    double low = (3.0/30.0)*(j+0.0);
    double high = (3.0/30.0)*(j+4.0);
    TString title = Form("%.1f #leq P_{diff} [GeV] #leq %.1f",low, high);
    TH1*  hproj = h162->ProjectionX(title,j, j+4);
    hproj->SetFillColor(kBlue);
    TString export_name = Form("bins_%ito%i",j,j+4);
    TString export_path = outdir_color+gMode+"_sum_"+export_name+".png";
    export_hist(hproj,export_path);
    delete hproj;
  }
  
  
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

void export_hist(TH1 * h2, TString out_filename) {
  auto c1 = new TCanvas();
  c1->SetCanvasSize(800,600);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.15);
  h2->Draw();

  c1->SaveAs(out_filename);
  c1->Modified();
  c1->Update();
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
  TString connector2 = connector+"/"+gMode;
  ridge_plot(corr->GetSE(),"ridge_same_"+connector+".png", connector2);
  ridge_plot(corr->GetME(),"ridge_multi_"+connector+".png", connector2);
  ridge_plot(corr->GetCO(),"ridge_correlation_"+connector+".png", connector2);
}

void full_1d_plots(TCorrelation* corr)
{
  for ( int var = 1; var < 3; var++ ) {
    TString con = "1D_"+gMode+"_"+corr->GetVar(var)+"_";
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
  v1.Boost(0,0,boost);
  return v1;
}


Float_t DeltaAngleRad(Float_t x_rad, Float_t y_rad)
{
  Double_t x = x_rad*180./TMath::Pi();
  Double_t y = y_rad*180./TMath::Pi();

  Float_t result;

  result = abs(y-x);
  
  if ( result > 360. ) {
    cout << "Angle difference " << result << ", better check your code!" << endl;
  }
  if ( result > 180. ) {
    result = 360.-result;
  }

  return result;
}
