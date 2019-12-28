class TCorrelation
{
private:
  TH2F fSameEvent2D;
  TH2F fMultiEvent2D;
  TH2F fCorrelation2D;
  TH2F fReco;
  TH2F fRecoPlus;
  TH2F fRecoMinus;

  TH1F fSameEvent1D1;
  TH1F fMultiEvent1D1;
  TH1F fCorrelation1D1;

  TH1F fSameEvent1D2;
  TH1F fMultiEvent1D2;
  TH1F fCorrelation1D2;

  TString fVar1;
  TString fVar2;
  TString fTarget;
  int fNbinsx = 32;
  int fNbinsy = fNbinsx;
  //double nedgex = TMath::Pi()*1.0;
  double fBinMaxX = 190.;
  double fBinMinX = -10;
  //double nedgey = TMath::Pi()*2.5;
  double fBinMaxY = 180.;
  double fBinMinY = -180.;
public:
  TCorrelation(TString, TString, TString);
  TH2F GetSE()   { return fSameEvent2D;   }
  TH2F GetME()   { return fMultiEvent2D;  }
  TH2F GetCO()   { return fCorrelation2D;  }
  TH2F GetReco()   { return fReco; }
  TH2F GetRecoPlus()   { return fRecoPlus; }
  TH2F GetRecoMinus()   { return fRecoMinus; }
  TH2F GetCorr() { return fCorrelation2D; }
  TH1F GetCorr1D(int i = 1) {
    if ( i == 2 ) {
      return fCorrelation1D2;
    } else {
      return fCorrelation1D1;
    }
  }
  TH1F GetSame1D(int i = 1) {
    if ( i == 2 ) {
      return fSameEvent1D2;
    } else {
      return fSameEvent1D1;
    }
  }
  TH1F GetMult1D(int i = 1) {
    if ( i == 2 ) {
      return fMultiEvent1D2;
    } else {
      return fMultiEvent1D1;
    }
  }
  TString GetVar( int a ) {
    if ( a == 1 ) return fVar1;
    else if ( a == 2 ) return fVar2;
    else return "use GetVar(a) with a==1 or a==2";
  }
  void FillSame(Float_t var1, Float_t var2, Float_t weight = 1) {
    fSameEvent2D.Fill(var1,var2, weight);
    fSameEvent1D1.Fill(var1,weight);
    fSameEvent1D2.Fill(var2,weight);
  }
  void FillMulti(Float_t var1, Float_t var2, Float_t weight = 1) {
    fMultiEvent2D.Fill(var1,var2, weight);
    fMultiEvent1D1.Fill(var1, weight);
    fMultiEvent1D2.Fill(var2, weight);
  }
  void FillReco(Float_t var1, Float_t var2, Float_t weight = 1) {
    fReco.Fill(var1,var2, weight);
  }
  void FillRecoPlus(Float_t var1, Float_t var2, Float_t weight = 1) {
    fRecoPlus.Fill(var1,var2, weight);
  }
  void FillRecoMinus(Float_t var1, Float_t var2, Float_t weight = 1) {
    fRecoMinus.Fill(var1,var2, weight);
  }
  void FillCorrelation();
  void SetBins2D(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  void NormalizeSame(Double_t);
  void NormalizeMulti();
  Double_t GetSameIntegral() { return fSameEvent2D.Integral();}
  Double_t GetCorrIntegral() { return fCorrelation2D.Integral();}
  Double_t GetRealCorrIntegral() {
    double retorno = 0;
    for ( int x = 0; x < fNbinsx+1; x++ ) {
      for ( int y = 0; y < fNbinsy+1; y++ ) {
	retorno += fCorrelation2D.GetBinContent(x,y);
      }
    }
    return retorno;
  }
  
};

TCorrelation::TCorrelation(TString var1, TString var2, TString target)
{
  fVar1 = var1;
  fVar2 = var2;
  fTarget = target;
  TH2F ns("ns_"+var1+"_"+var2+"_"+target,
	  "N_{s}(#Delta#"+var1+", #Delta#"+var2+");#Delta#"+var1+";#Delta#"+var2,
	  fNbinsx,fBinMinX,fBinMaxX,
	  fNbinsy,fBinMinY,fBinMaxY
	  );
  TH2F nm("nm_"+var1+"_"+var2+"_"+target,
	  "N_{m}(#Delta#"+var1+", #Delta#"+var2+");#Delta#"+var1+";#Delta#"+var2,
	  fNbinsx,fBinMinX,fBinMaxX,
	  fNbinsy,fBinMinY,fBinMaxY
	  );
  TH1F ns11("ns_"+var1+"_"+target,
	    "N_{s}(#Delta#"+var1+");#Delta#"+var1,
	    fNbinsx,fBinMinX,fBinMaxX
	    );
  TH1F ns12("ns_"+var2+"_"+target,
	    "N_{s}(#Delta#"+var2+");#Delta#"+var2,
	    fNbinsy,fBinMinY,fBinMaxY
	    );
  TH1F nm11("nm_"+var1+"_"+target,
	    "N_{m}(#Delta#"+var1+");#Delta#"+var1,
	    fNbinsx,fBinMinX,fBinMaxX
	    );
  TH1F nm12("nm_"+var2+"_"+target,
	    "N_{m}(#Delta#"+var2+");#Delta#"+var2,
	    fNbinsy,fBinMinY,fBinMaxY
	    );
  fSameEvent2D  = ns;
  fMultiEvent2D = nm;
  fSameEvent1D1 = ns11;
  fSameEvent1D2 = ns12;
  fMultiEvent1D1 = nm11;
  fMultiEvent1D2 = nm12;

  // reco

  TH2F reco1("reco_trig_partner_"+var1+"_"+var2+"_"+target,
	     "reconstructed #pi^{#pm};#Delta#"+var1+";#Delta#"+var2,
	     72,-180,180,36,0,180);
  TH2F reco2("reco_trig_"+var1+"_"+var2+"_"+target,
	     "reconstructed #pi^{+};#Delta#"+var1+";#Delta#"+var2,
	     72,-180,180,36,0,180);
  TH2F reco3("reco_partner_"+var1+"_"+var2+"_"+target,
	     "reconstructed #pi^{-};#Delta#"+var1+";#Delta#"+var2,
	     72,-180,180,36,0,180);
  fReco = reco1;
  fRecoPlus = reco2;
  fRecoMinus = reco3;
  
}

void TCorrelation::FillCorrelation()
{
  TH2F ndiv("nc_"+fVar1+"_"+fVar2+"_"+fTarget,
	    "N_{s}/N_{m}(#Delta#"+fVar1+", #Delta#"+fVar2+");#Delta#"+fVar1+";#Delta#"+fVar2,
	    fNbinsx,fBinMinX,fBinMaxX,
	    fNbinsy,fBinMinY,fBinMaxY);
  /*
  TH1F ndiv1d1("nc_"+fVar1,
	       "N_{s}/N_{m}(#Delta#"+fVar1+");#Delta#"+fVar1,
	       fNbinsy,fBinMinY,fBinMaxY);
  TH1F ndiv1d2("nc_"+fVar2,
	       "N_{s}/N_{m}(#Delta#"+fVar2+");#Delta#"+fVar2,
	       fNbinsx,fBinMinX,fBinMaxX);
  */
  float ratio;
  float settle_at = 1.0;
  for ( int x = 1; x < fNbinsx+1; x++ ) {
    for ( int y = 1; y < fNbinsy+1; y++ ) {
      float nsvalue = fSameEvent2D.GetBinContent(x,y);
      float nmvalue = fMultiEvent2D.GetBinContent(x,y);
      if ( nsvalue < 5.0 || nmvalue < 5.0 ) {
	ndiv.SetBinContent(x,y,settle_at);
      } else {
	ndiv.SetBinContent(x,y,nsvalue/nmvalue);
      }
    }
  }

  // 1 d histograms
  double num;
  double denom;
  
  TH1F h1 = fSameEvent1D1;
  for ( int i = 0; i < fNbinsx; i++ ) {
    num   = 0;
    denom = 0;
    num   = fSameEvent1D1.GetBinContent(i);
    denom = fMultiEvent1D1.GetBinContent(i);
    if ( denom != 0 ) {
      h1.Fill(fSameEvent1D1.GetBinCenter(i),
	      (num/denom));
    }
  }
  h1.Scale(1./h1.Integral());
  h1.SetNameTitle("nc_"+fVar1+"_"+fTarget,"N_{s}/N_{m}(#Delta#"+fVar1+");#Delta#"+fVar1);
  
  TH1F h2 = fSameEvent1D2;
  for ( int i = 0; i < fNbinsy; i++ ) {
    num   = 0;
    denom = 0;
    num   = fSameEvent1D2.GetBinContent(i);
    denom = fMultiEvent1D2.GetBinContent(i);
    if ( denom != 0 ) {
      h2.Fill(fSameEvent1D2.GetBinCenter(i),
	      (num/denom));
    }
  }
  h2.Scale(1./h2.Integral());
  h2.SetNameTitle("nc_"+fVar2+"_"+fTarget,"N_{s}/N_{m}(#Delta#"+fVar2+");#Delta#"+fVar2);
  
  fCorrelation1D1 = h1;
  fCorrelation1D2 = h2;
  fCorrelation2D = ndiv;
}


void TCorrelation::SetBins2D(Float_t nbinsx, Float_t binminx, Float_t binmaxx,
			     Float_t nbinsy, Float_t binminy, Float_t binmaxy)
{
  fNbinsx = nbinsx;
  fNbinsy = nbinsy;
  fBinMaxX = binmaxx;
  fBinMinX = binminx;
  fBinMaxY = binmaxy;
  fBinMinY = binminy;
}


void TCorrelation::NormalizeSame(Double_t N)
{
  for ( Int_t x = 0; x < fNbinsx; x++ ) {
    for ( Int_t y = 0; y < fNbinsy; y++ ) {
      Double_t content = fSameEvent2D.GetBinContent(x,y);
      fSameEvent2D.SetBinContent(x,y,content/N);
    }
  }
}

void TCorrelation::NormalizeMulti()
{
  Double_t N = fMultiEvent2D.GetMaximum();
  for ( Int_t x = 0; x < fNbinsx; x++ ) {
    for ( Int_t y = 0; y < fNbinsy; y++ ) {
      Double_t content = fMultiEvent2D.GetBinContent(x,y);
      fMultiEvent2D.SetBinContent(x,y,content/N);
      //fMultiEvent2D.GetBin(x,y) *= (1/N);
    }
  }
}
