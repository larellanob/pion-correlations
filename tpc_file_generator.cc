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
#include "Root/DeltaAngleRad.cxx"

bool m_debug      = false;
bool m_simulation = false;
bool old_triggers     = false; // if false, use old_partners
bool DMode        = false; // deuterium
Float_t gDataCap  = 1.0; // fraction of data to be used in analysis (1.0 == full data)
TDatabasePDG db;
// modes
// pp: trigger pi+ no cuts, assoc pi+ no cuts
// pm: trigger pi+ no cuts, assoc pi- no cuts
// zhpm: trigger pi+ zh > 0.5, assocc pi- zh < 0.5
// zhpp: trigger pi+ zh > 0.5, assocc pi+ zh < 0.5
// zhlluupm: trigger pi+ zh>0.5, assoc pi- zh in (ll,uu)
// zhlluupp: trigger pi+ zh>0.5, assoc pi+ zh in (ll,uu)
TString gMode = "zhpm"; // choose pm or pp

TLorentzVector VirtualFrame(TLorentzVector, double, double, double);
Float_t DeltaAngleRad(Float_t x_rad, Float_t y_rad);


// MAIN
void tpc_file_generator(TString mode = "", TString target = "")
{
  if ( mode == "" ) {
    cout << "No mode input, using zhpm as default" << endl;
    mode = gMode;
  } else {
    gMode = mode;
  }
  
  if ( target == "" ) {
    cout << "No target input, using C as default" << endl;
    target = "C";
  }
  //gStyle->SetOptStat(0);


  // Chain for hadrons (ch) and chain for triggers electrons (ech)
  TChain ch;
  if ( !m_simulation ) {
    ch.Add("tree_data");
  }
    
  TChain ech("e_rec");

  // Add files to the chains depending on target
  // Also sets beam energy
  Double_t kEbeam;
  if ( m_simulation ) {
    if ( target == "Pb" ) {
      ch.Add("/eos/user/a/arellano/CFF/full_Pb_simulations.root/tree_accept");
      ech.Add("/eos/user/a/arellano/CFF/full_Pb_simulations.root");
    } else if ( target == "Fe" ) {
      ch.Add("/eos/user/a/arellano/CFF/full_Fe_simulations.root/tree_accept");
      ech.Add("/eos/user/a/arellano/CFF/full_Fe_simulations.root");
    } else if ( target == "Pb" ) {
      ch.Add("/eos/user/a/arellano/CFF/full_Pb_simulations.root/tree_accept");
      ech.Add("/eos/user/a/arellano/CFF/full_Pb_simulations.root");
    } else if ( target == "D" ) {
      ch.Add("/eos/user/a/arellano/CFF/full_D_simulations.root/tree_accept");
      ech.Add("/eos/user/a/arellano/CFF/full_D_simulations.root");
      DMode = true;
    }
    //kEbeam = 11.0;
    kEbeam = 5.014;
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
    gDataCap *= 0.5;
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
  TTreeReaderArray<Float_t> pt(th,"Pt");

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

  // Setting some counters
  UInt_t events = 0;
  int processed_events=0;
    
  TLorentzVector virtual_photon;
  //TLorentzVector oldgamma;

  // Number of trigger particles (for normalization)
  int paired_triggers = 0;
  int total_triggers = 0;
  
  
  if ( m_debug )
    cout << "DEBUG: entering main tree reading loop" << endl;
  // main file/tree reading loop


  ///////
  // Export ttrees
  TFile out_tree;
  if ( !m_simulation ) {
    out_tree.Open("/home/luciano/Physics/CLAS/pion_correlation/2pcPairs/2pcPairs_"+target+"_"+gMode+".root", "RECREATE");
  } else {
    out_tree.Open("/eos/user/a/arellano/CLAS/pion_correlation/2pcPairs/2pcPairs_"+target+"_"+gMode+".root", "RECREATE");
  }
  
  TTree triggers_tree("triggers","trigger particles");
  TTree partners_tree("partners","partner particles");
  TTree old_evnt_tree("old_evnt","old partners for mixed events ");

  Float_t Q2_w1;
  Float_t Q2_w2;
  Float_t Q2_w3;
  TString Target_w1 = target;
  TString Target_w2 = target;
  TString Target_w3 = target;
  std::vector<Float_t> Zh_w1;
  std::vector<Float_t> Zh_w2;
  std::vector<Float_t> Zh_w3;
  std::vector<Int_t> pid_w1;
  std::vector<Int_t> pid_w2;
  std::vector<Int_t> pid_w3;
  
  std::vector<Float_t> e_w1;
  std::vector<Float_t> x_w1;
  std::vector<Float_t> y_w1;
  std::vector<Float_t> z_w1;
  std::vector<Float_t> theboo_w1;
  std::vector<Float_t> phiboo_w1;
  std::vector<Float_t> rap_w1;
  std::vector<Float_t> Pt_w1;
  
  std::vector<Float_t> e_w2;
  std::vector<Float_t> x_w2;
  std::vector<Float_t> y_w2;
  std::vector<Float_t> z_w2;
  std::vector<Float_t> theboo_w2;
  std::vector<Float_t> phiboo_w2;
  std::vector<Float_t> rap_w2;
  std::vector<Float_t> Pt_w2;
  
  std::vector<Float_t> e_w3;
  std::vector<Float_t> x_w3;
  std::vector<Float_t> y_w3;
  std::vector<Float_t> z_w3;
  
  //std::vector<TLorentzVector> pions_w1;
  //std::vector<TLorentzVector> pions_w2;
  //std::vector<TLorentzVector> pions_w3;



  
  triggers_tree.Branch("Target",&Target_w1);
  triggers_tree.Branch("Q2",&Q2_w1);
  triggers_tree.Branch("Zh",&Zh_w1);
  triggers_tree.Branch("pid",&pid_w1);
  triggers_tree.Branch("E",&e_w1);
  triggers_tree.Branch("Px",&x_w1);
  triggers_tree.Branch("Py",&y_w1);
  triggers_tree.Branch("Pz",&z_w1);
  triggers_tree.Branch("ThetaBoosted",&theboo_w1);
  triggers_tree.Branch("PhiBoosted",&phiboo_w1);
  triggers_tree.Branch("y",&rap_w1);
  triggers_tree.Branch("Pt",&Pt_w1);
  
  partners_tree.Branch("Target",&Target_w2);
  partners_tree.Branch("Q2",&Q2_w2);
  partners_tree.Branch("Zh",&Zh_w2);
  partners_tree.Branch("pid",&pid_w2);
  partners_tree.Branch("E",&e_w2);
  partners_tree.Branch("Px",&x_w2);
  partners_tree.Branch("Py",&y_w2);
  partners_tree.Branch("Pz",&z_w2);
  partners_tree.Branch("ThetaBoosted",&theboo_w2);
  partners_tree.Branch("PhiBoosted",&phiboo_w2);
  partners_tree.Branch("y",&rap_w2);
  partners_tree.Branch("Pt",&Pt_w2);

  
  old_evnt_tree.Branch("Target",&Target_w3);
  old_evnt_tree.Branch("Q2",&Q2_w3);
  old_evnt_tree.Branch("Zh",&Zh_w3);
  old_evnt_tree.Branch("pid",&pid_w3);
  old_evnt_tree.Branch("E",&e_w3);
  old_evnt_tree.Branch("Px",&x_w3);
  old_evnt_tree.Branch("Py",&y_w3);
  old_evnt_tree.Branch("Pz",&z_w3);
  

  // first event, no old 
  //old_evnt_tree.Fill();
  
  
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
    // Zh > 1 veto and pion count
    double zh_sum = 0;
    int pioncount = 0;
    for ( int i = 0; i < Zh.GetSize(); i++ ) {
      if ( pid[i] == 211 || pid[i] == -211 ) {
	zh_sum += Zh[i];
	pioncount++;
      }
    }
    if ( zh_sum > 1.1 ) { //|| pioncount < 2) {
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
    // ppp: trigger pi+ with highest zh, partner all other pi+
    // pp: trigger pi+ no cuts, assoc pi+ no cuts
    // pm: trigger pi+ no cuts, assoc pi- no cuts
    // zhpm: trigger pi+ zh > 0.5, assocc pi- zh < 0.5
    // zhpp: trigger pi+ zh > 0.5, assocc pi+ zh < 0.5
    // zhlluupm: trigger pi+ zh>0.5, assoc pi- zh in (ll,uu)
    // zhlluupp: trigger pi+ zh>0.5, assoc pi+ zh in (ll,uu)
    // high-low zh 
    if ( gMode == "zhlead" ) {
      double zh_max = 0;
      // first pi+ found is stored at triggers_index[0]
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] != 211 && pid[i] != -211 ) {
	  continue;
	}
	if ( Zh[i] > zh_max ) {
	  zh_max = Zh[i];
	  if ( triggers_index.size() > 0 ) {
	    partners_index.push_back(triggers_index[0]);
	  }
	  triggers_index.clear();
	  triggers_index.push_back(i);
	} else {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "ptlead" ) {
      double pt_lead = 0;
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] != 211 && pid[i] != -211 ) {
	  continue;
	}
	if ( pt[i] > pt_lead ) {
	  pt_lead = pt[i];
	  if ( triggers_index.size() > 0 ) {
	    partners_index.push_back(triggers_index[0]);
	  }
	  triggers_index.clear();
	  triggers_index.push_back(i);
	} else {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "pp" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 ) {
	  triggers_index.push_back(i);
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "pm" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zhpm" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zhpp" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == +211 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0305pp" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.3 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0305pm" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.3 && Zh[i] < 0.5 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0203pp" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.2 && Zh[i] < 0.3 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0203pm" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.2 && Zh[i] < 0.3 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0002pp" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == 211 && Zh[i] > 0.0 && Zh[i] < 0.2 ) {
	  partners_index.push_back(i);
	}
      }
    } else if ( gMode == "zh0002pm" ) {
      for ( UInt_t i = 0; i < pid.GetSize(); i++ ) {
	if ( pid[i] == 211 && Zh[i] > 0.5 ) {
	  triggers_index.push_back(i);
	}
	if ( pid[i] == -211 && Zh[i] > 0.0 && Zh[i] < 0.2 ) {
	  partners_index.push_back(i);
	}
      }
    }


    if ( triggers_index.size() != 0 ) {
      total_triggers++;
    }
    
    // make sure of this, normalization might depend on it
    if ( triggers_index.size() == 0  || partners_index.size() == 0 ) {
      continue;
    }

    
    ///////////////////////////
    // Pions 4vectors
    triggers4v.clear();
    partners4v.clear();
    triggers4v_rotated.clear();
    partners4v_rotated.clear();
    triggers4v_boosted.clear();
    partners4v_boosted.clear();

    for ( int i:  triggers_index ) {
      Float_t e = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass);
      TLorentzVector pi (px[i],py[i],pz[i],e);
      triggers4v.push_back(pi);
    }

    for ( int i: partners_index ) {
      Float_t e = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+kChargedPionMass*kChargedPionMass);
      TLorentzVector pi (px[i],py[i],pz[i],e);
      partners4v.push_back(pi);
    }


    Zh_w1.clear(); Zh_w2.clear();
    pid_w1.clear(); pid_w2.clear();
    e_w1.clear(); e_w2.clear();
    x_w1.clear(); x_w2.clear();
    y_w1.clear(); y_w2.clear();
    z_w1.clear(); z_w2.clear();
    theboo_w1.clear(); theboo_w2.clear();
    phiboo_w1.clear(); phiboo_w2.clear();
    rap_w1.clear(); rap_w2.clear();
    Pt_w1.clear(); Pt_w2.clear();

    for ( int i: triggers_index ) {
      Zh_w1.push_back(Zh[i]);
      pid_w1.push_back(pid[i]);
      Q2_w1 = *eq2;
    }
    for ( int j: partners_index ) {
      Zh_w2.push_back(Zh[j]);
      pid_w2.push_back(pid[j]);
      Q2_w2 = *eq2;
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
    for ( int i = 0; i < triggers4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(triggers4v[i],angle1,angle2,0);
      triggers4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(triggers4v[i],angle1,angle2,-boost);
      triggers4v_boosted.push_back(aux_boost);
      e_w1.push_back(aux_boost.E());
      x_w1.push_back(aux_boost.X());
      y_w1.push_back(aux_boost.Y());
      z_w1.push_back(aux_boost.Z());
      theboo_w1.push_back(aux_boost.Theta()*180./TMath::Pi());
      phiboo_w1.push_back(aux_boost.Phi()*180./TMath::Pi());
      rap_w1.push_back(0.5* log((aux_boost.E()+aux_boost.Z()) / (aux_boost.E()-aux_boost.Z())) );
      Pt_w1.push_back(sqrt(aux_boost.X()*aux_boost.X()+aux_boost.Y()*aux_boost.Y()));
    }

    /*
    if ( events > 10 && events < 1000 ) {
      for ( int i = 0; i < triggers_index.size(); i++ ) {
	cout << "old " << pt[i] << endl;
      }

      for ( int i = 0; i < triggers4v.size(); i++ ) {
	cout << "new " << Pt_w1[i] << endl;
      }
      cout <<  "----" << endl;
    }
    if ( events > 1000 ) {
      break;
    }
    */

    for ( int i = 0; i < partners4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(partners4v[i],angle1,angle2,0);
      partners4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(partners4v[i],angle1,angle2,-boost);
      partners4v_boosted.push_back(aux_boost);
      e_w2.push_back(aux_boost.E());
      x_w2.push_back(aux_boost.X());
      y_w2.push_back(aux_boost.Y());
      z_w2.push_back(aux_boost.Z());
      theboo_w2.push_back(aux_boost.Theta()*180./TMath::Pi());
      phiboo_w2.push_back(aux_boost.Phi()*180./TMath::Pi());
      rap_w2.push_back(0.5* log((aux_boost.E()+aux_boost.Z()) / (aux_boost.E()-aux_boost.Z())) );
      Pt_w2.push_back(sqrt(aux_boost.X()*aux_boost.X()+aux_boost.Y()*aux_boost.Y()));
    }
    triggers_tree.Fill();
    partners_tree.Fill();
    
    // rotation and boost for old particle
    //Zh_w3.clear(); pid_w3.clear(); e_w3.clear(); x_w3.clear(); y_w3.clear(); z_w3.clear();
    
    for ( int i = 0; i < old4v.size(); i++ ) {
      TLorentzVector aux_rotation;
      aux_rotation = VirtualFrame(old4v[i],angle1,angle2,0);
      old4v_rotated.push_back(aux_rotation);
      TLorentzVector aux_boost;
      aux_boost = VirtualFrame(old4v[i],angle1,angle2,-boost);
      old4v_boosted.push_back(aux_boost);
      e_w3.push_back(aux_boost.E());
      x_w3.push_back(aux_boost.X());
      y_w3.push_back(aux_boost.Y());
      z_w3.push_back(aux_boost.Z());
    }
    
    Zh_w3.clear(); pid_w3.clear(); e_w3.clear(); x_w3.clear(); y_w3.clear(); z_w3.clear();
    for ( auto k : partners_index ) {
      Zh_w3.push_back(Zh[k]);
      pid_w3.push_back(pid[k]);
    }
    Q2_w3 = *eq2;


    
    
    if ( m_debug )
      cout << "DEBUG: Entering mutli event filling" << endl;



    old4v.clear();
    old4v_rotated.clear();
    old4v_boosted.clear();
    
    if ( m_debug )
      cout << "DEBUG: Entering same event filling" << endl;
    // SAME EVENT FILLING
    bool fill_neg = true;
    

    //////////////////////////////////////
    // Old pion saving

    // Old negative pion:

    for ( int j = 0; j < partners4v.size(); j++ ) {
      old4v.push_back(partners4v[j]);
      //old4v_rotated.push_back(partners4v_rotated[j]);
      //old4v_boosted.push_back(partners4v_boosted[j]);
    }


    paired_triggers += triggers_index.size();
    //    if ( triggers_index.size() < 1 )
    //cout << "what " << triggers_index.size() << endl;
    processed_events++;
  }

  cout << "ran through this many events: " << events << endl;
  cout << "processed this many events: " << processed_events << endl;
  cout << "number of trigger particles: " << paired_triggers << endl;
  //side_by_side(ns,ns_vir,"NS_lab_virtual.png");

  TString multi_method;
  if ( old_triggers )
    multi_method = "_pastTriggers";
  else if ( !old_triggers )
    multi_method = "_pastPartners";

  triggers_tree.Write();
  partners_tree.Write();
  //old_evnt_tree.Write();
  TTree Normalization("Metadata","Miscelaneous data which might be helpful for the analysis");
  Normalization.Branch("total_triggers",&total_triggers);
  Normalization.Branch("paired_triggers",&paired_triggers);
  Normalization.Branch("data_cap",&gDataCap);
  Normalization.Branch("target",&target);
  Normalization.Branch("trigger_mode",&gMode);
  Normalization.Fill();
  Normalization.Write();

  
  out_tree.Close();

  cout << "finished correctly" << endl;
} // END MAIN


////////////////////////
//// FUNCTIONS
////////////////////////

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
