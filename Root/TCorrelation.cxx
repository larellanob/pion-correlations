class TCorrelation
{
private:
  TH2F fSameEvent2D;
  TH2F fMultiEvent2D;
  TH2F fCorrelation2D;
  TH2F fReco;
  TH2F fRecoTriggers;
  TH2F fRecoPartners;

  TH1F fSameEvent1D1;
  TH1F fMultiEvent1D1;
  TH1F fCorrelation1D1;

  TH1F fSameEvent1D2;
  TH1F fMultiEvent1D2;
  TH1F fCorrelation1D2;

  TString fVar1;
  TString fVar2;
  TString fTarget;

  // binsx: phi (-180,180)
  // binsy: y (-3,4), theta (0,180)
  // binsdeltax: dphi (0,180)
  // binsdeltay: dy (), dtheta (0,180)
  int fNbinsx = 32;
  int fNbinsy = fNbinsx;
  double fBinMaxX = 180.;
  double fBinMinX = -180;
  double fBinMaxY = 180.;
  double fBinMinY = 0.;
  //int fNbinsDeltax = 32;
  int fNbinsDeltax = 16;
  int fNbinsDeltay = fNbinsDeltax;
  double fBinMaxDeltaX = 180.;
  double fBinMinDeltaX = 0;
  double fBinMaxDeltaY = 180.;
  double fBinMinDeltaY = -20;
  
public:
  TCorrelation(TString, TString, TString);
  TH2F GetSE()   { return fSameEvent2D;   }
  TH2F GetME()   { return fMultiEvent2D;  }
  TH2F GetCO()   { return fCorrelation2D;  }
  TH2F GetReco()   { return fReco; }
  TH2F GetRecoTriggers()   { return fRecoTriggers; }
  TH2F GetRecoPartners()   { return fRecoPartners; }
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
  void FillRecoTriggers(Float_t var1, Float_t var2, Float_t weight = 1) {
    fRecoTriggers.Fill(var1,var2, weight);
  }
  void FillRecoPartners(Float_t var1, Float_t var2, Float_t weight = 1) {
    fRecoPartners.Fill(var1,var2, weight);
  }
  void FillCorrelation();
  void SetBins2D(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Double_t GetSameIntegral() { return fSameEvent2D.Integral();}
  Double_t GetCorrIntegral() { return fCorrelation2D.Integral();}
  Double_t GetRealCorrIntegral() {
    double retorno = 0;
    for ( int x = 0; x < fNbinsDeltax+1; x++ ) {
      for ( int y = 0; y < fNbinsDeltay+1; y++ ) {
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
  if ( var1 == "y" ) {
    fBinMaxX = 4.;
    fBinMinX = -3;
    fBinMinDeltaX = -3;
    fBinMaxDeltaX = 3;
  }
  if ( var2 == "y" ) {
    fBinMaxY = 4.;
    fBinMinY = -3;
    fBinMinDeltaY = -3;
    fBinMaxDeltaY = 3;
  }
  
  fTarget = target;
  TH2F ns("ns_"+var1+"_"+var2+"_"+target,
	  "N_{s}(#Delta#"+var1+", #Delta#"+var2+");#Delta#"+var1+";#Delta#"+var2,
	  fNbinsDeltax,fBinMinDeltaX,fBinMaxDeltaX,
	  fNbinsDeltay,fBinMinDeltaY,fBinMaxDeltaY
	  );
  TH2F nm("nm_"+var1+"_"+var2+"_"+target,
	  "N_{m}(#Delta#"+var1+", #Delta#"+var2+");#Delta#"+var1+";#Delta#"+var2,
	  fNbinsDeltax,fBinMinDeltaX,fBinMaxDeltaX,
	  fNbinsDeltay,fBinMinDeltaY,fBinMaxDeltaY
	  );
  TH1F ns11("ns_"+var1+"_"+target,
	    "N_{s}(#Delta#"+var1+");#Delta#"+var1,
	    fNbinsDeltax,fBinMinDeltaX,fBinMaxDeltaX
	    );
  TH1F ns12("ns_"+var2+"_"+target,
	    "N_{s}(#Delta#"+var2+");#Delta#"+var2,
	    fNbinsDeltay,fBinMinDeltaY,fBinMaxDeltaY
	    );
  TH1F nm11("nm_"+var1+"_"+target,
	    "N_{m}(#Delta#"+var1+");#Delta#"+var1,
	    fNbinsDeltax,fBinMinDeltaX,fBinMaxDeltaX
	    );
  TH1F nm12("nm_"+var2+"_"+target,
	    "N_{m}(#Delta#"+var2+");#Delta#"+var2,
	    fNbinsDeltay,fBinMinDeltaY,fBinMaxDeltaY
	    );
  fSameEvent2D  = ns;
  fMultiEvent2D = nm;
  fSameEvent1D1 = ns11;
  fSameEvent1D2 = ns12;
  fMultiEvent1D1 = nm11;
  fMultiEvent1D2 = nm12;

  // reco
  TH2F reco1("reco_trig_partner_"+var1+"_"+var2+"_"+target,
	     "reconstructed triggers + partners;#"+var1+";#"+var2,
	     fNbinsx,fBinMinX,fBinMaxX,fNbinsy,fBinMinY,fBinMaxY);
  TH2F reco2("reco_trig_"+var1+"_"+var2+"_"+target,
	     "reconstructed triggers;#"+var1+";#"+var2,
	     fNbinsx,fBinMinX,fBinMaxX,fNbinsy,fBinMinY,fBinMaxY);
  TH2F reco3("reco_partner_"+var1+"_"+var2+"_"+target,
	     "reconstructed partners;#"+var1+";#"+var2,
	     fNbinsx,fBinMinX,fBinMaxX,fNbinsy,fBinMinY,fBinMaxY);
  fReco = reco1;
  fRecoTriggers = reco2;
  fRecoPartners = reco3;
  
}

void TCorrelation::FillCorrelation()
{
  TH2F ndiv("nc_"+fVar1+"_"+fVar2+"_"+fTarget,
	    "N_{s1}/N_{m}(#Delta#"+fVar1+", #Delta#"+fVar2+");#Delta#"+fVar1+";#Delta#"+fVar2,
	    fNbinsDeltax,fBinMinDeltaX,fBinMaxDeltaX,
	    fNbinsDeltay,fBinMinDeltaY,fBinMaxDeltaY);

  float ratio;
  float settle_at = 1.0;
  double error1;
  for ( int x = 1; x < fNbinsDeltax+1; x++ ) {
    double bcx = (ndiv.GetXaxis()->GetBinCenter(x));
    for ( int y = 1; y < fNbinsDeltay+1; y++ ) {
      double bcy = (ndiv.GetYaxis()->GetBinCenter(y));
      float nsvalue = fSameEvent2D.GetBinContent(x,y);
      float nmvalue = fMultiEvent2D.GetBinContent(x,y);
      if ( nsvalue < 5.0 || nmvalue < 5.0 ) {
	//ndiv.SetBinContent(x,y,settle_at);
	ndiv.Fill(bcx,bcy,settle_at);
      } else {
	//ndiv.SetBinContent(x,y,nsvalue/nmvalue);
	ndiv.Fill(bcx,bcy,nsvalue/nmvalue);
      }
      if ( nmvalue != 0 && nsvalue != 0) {
	error1 = (nsvalue/nmvalue)*sqrt( (1./nsvalue) + (1./nmvalue) );
	ndiv.SetBinError(x,y,error1);
      }

    }
  }

  // 1 d histograms
  double num;
  double denom;
  double error;
  
  TH1F h1 = fSameEvent1D1;
  h1.Reset();
  double h1norm = fMultiEvent1D1.GetBinContent(fMultiEvent1D1.GetMaximumBin());
  for ( int i = 1; i < fNbinsDeltax+1; i++ ) {
    num   = 0;
    denom = 0;
    num   = fSameEvent1D1.GetBinContent(i);
    denom = fMultiEvent1D1.GetBinContent(i);
    if ( denom != 0 && num > 50 && denom > 50 ) {
      h1.Fill(fSameEvent1D1.GetBinCenter(i),
	      (num/denom));
      error = (num/denom)*sqrt( (1./num) + (1./denom)  );
      h1.SetBinError(i,error);
    }
  }
  h1.SetNameTitle("nc_"+fVar1+"_"+fTarget,"N_{s2}/N_{m}(#Delta#"+fVar1+");#Delta#"+fVar1);
  
  TH1F h2 = fSameEvent1D2;
  h2.Reset();
  double h2norm = fMultiEvent1D2.GetBinContent(fMultiEvent1D2.GetMaximumBin());
  for ( int i = 1; i < fNbinsDeltay+1; i++ ) {
    num   = 0;
    denom = 0;
    num   = fSameEvent1D2.GetBinContent(i);
    denom = fMultiEvent1D2.GetBinContent(i);
    if ( denom != 0 && num > 50 && denom > 50 ) {
      h2.Fill(fSameEvent1D2.GetBinCenter(i),
	      (num/denom));
      error = (num/denom)*sqrt( (1./num) + (1./denom)  );
      h2.SetBinError(i,error);
    }
  }
  h2.SetNameTitle("nc_"+fVar2+"_"+fTarget,"N_{s3}/N_{m}(#Delta#"+fVar2+");#Delta#"+fVar2);
  
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

