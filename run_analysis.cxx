#include "Root/TCorrelation.cxx"
#include "Root/DeltaAngle.cxx"

void run_analysis( TString filename = "/home/luciano/Physics/CLAS/pion_ridge/data/2pcPairs_C_zhpm.root" )
{
  if ( filename == "" ) {
    cout << "Please input some filename to analyze" << endl;
    return;
  }
    
  TFile f(filename);

  // metadata
  // mode, target, % of data, etc.
  TString target = "C";

  
  TTreeReader triggers("triggers",&f);
  TTreeReader partners("partners",&f);
  TTreeReader old_evnt("old_evnt",&f);

  // triggers 
  TTreeReaderArray<Int_t>  pid_t(triggers,"pid");
  TTreeReaderValue<Float_t> q2_t(triggers,"Q2");
  TTreeReaderArray<Float_t> Zh_t(triggers,"Zh");

  TTreeReaderArray<Float_t> TheBoo_t(triggers,"ThetaBoosted");
  TTreeReaderArray<Float_t> PhiBoo_t(triggers,"PhiBoosted");
  TTreeReaderArray<Float_t> EtaBoo_t(triggers,"Eta");

  // partners 
  TTreeReaderArray<Int_t>  pid_p(partners,"pid");
  TTreeReaderValue<Float_t> q2_p(partners,"Q2");
  TTreeReaderArray<Float_t> Zh_p(partners,"Zh");

  TTreeReaderArray<Float_t> TheBoo_p(partners,"ThetaBoosted");
  TTreeReaderArray<Float_t> PhiBoo_p(partners,"PhiBoosted");
  TTreeReaderArray<Float_t> EtaBoo_p(partners,"Eta");

  TCorrelation * corr_boo = new TCorrelation("phiPQboosted","thetaPQboosted",target);
  TCorrelation * corr_eta = new TCorrelation("phiPQboosted2","Eta",target);
  
  while ( triggers.Next() ) {
    // you can get "old partners" by working here, before partners.Next()
    partners.Next();


    /// Same Event Pairs
    int t_size = pid_t.GetSize();
    int p_size = pid_p.GetSize();

    for ( int t = 0; t < t_size; t++ ) {
      corr_boo->FillReco(PhiBoo_t[t],TheBoo_t[t]);
      corr_eta->FillReco(PhiBoo_t[t],EtaBoo_t[t]);
      corr_boo->FillRecoTriggers(PhiBoo_t[t],TheBoo_t[t]);
      corr_eta->FillRecoTriggers(PhiBoo_t[t],EtaBoo_t[t]);
      
      
      for ( int p = 0; p < p_size; p++ ) {
	if ( t == 0 ) {
	  corr_boo->FillReco(PhiBoo_p[p],TheBoo_p[p]);
	  corr_eta->FillReco(PhiBoo_p[p],EtaBoo_p[p]);
	  corr_boo->FillRecoPartners(PhiBoo_p[p],TheBoo_p[p]);
	  corr_eta->FillRecoPartners(PhiBoo_p[p],EtaBoo_p[p]);
	}
	
	Float_t DPhi,DThe,DEta;
	DPhi = DeltaAngle(PhiBoo_t[t],PhiBoo_p[p]);
	DThe = DeltaAngle(TheBoo_t[t],TheBoo_p[p]);
	DEta = EtaBoo_t[t]-EtaBoo_p[p];
	corr_boo->FillSame(DPhi,DThe);
	corr_eta->FillSame(DPhi,DEta);
      }
    }
    

    
  }
  corr_boo->FillCorrelation();
  corr_eta->FillCorrelation();


  // Plotting
  /*
  full_ridge_plots(corr_boo);
  full_1d_plots(corr_boo);
  full_ridge_plots(corr_eta);
  full_1d_plots(corr_eta);
  */
  
}
