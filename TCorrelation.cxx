class TCorrelation
{
private:
  TH2F fSameEvent2D;
  TH2F fMultiEvent2D;
  TH2F fCorrelation2D;
  TH2F fReco;
  TH2F fRecoPlus;
  TH2F fRecoMinus;
  TString fVar1;
  TString fVar2;
  int nbinsx = 32;
  int nbinsy = nbinsx*2.5;
  double nedgex = TMath::Pi()*1.0;
  double nedgey = TMath::Pi()*2.5;
public:
  TCorrelation(TString, TString);
  TH2F GetSE()   { return fSameEvent2D;   }
  TH2F GetME()   { return fMultiEvent2D;  }
  TH2F GetCO()   { return fCorrelation2D;  }
  TH2F GetReco()   { return fReco; }
  TH2F GetRecoPlus()   { return fRecoPlus; }
  TH2F GetRecoMinus()   { return fRecoMinus; }
  TH2F GetCorr() { return fCorrelation2D; }
  TString GetVar( int a ) {
    if ( a == 1 ) return fVar1;
    else if ( a == 2 ) return fVar2;
    else return "use GetVar(a) with a==1 or a==2";
  }
  void FillSame(Float_t var1, Float_t var2, Float_t weight = 1) {
    fSameEvent2D.Fill(var1,var2, weight);
  }
  void FillMulti(Float_t var1, Float_t var2, Float_t weight = 1) {
    fMultiEvent2D.Fill(var1,var2, weight);
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
};

TCorrelation::TCorrelation(TString var1, TString var2)
{
  fVar1 = var1;
  fVar2 = var2;
  TH2F ns("ns_"+var1+"_"+var2,
	  "N_{s}(#Delta#"+var2+", #Delta#"+var1+");#Delta#"+var2+";#Delta#"+var1,
	  nbinsx,-nedgex,nedgex,
	  nbinsy,-nedgey,nedgey
	  );
  TH2F nm("nm_"+var1+"_"+var2,
	  "N_{m}(#Delta#"+var2+", #Delta#"+var1+");#Delta#"+var2+";#Delta#"+var1,
	  nbinsx,-nedgex,nedgex,
	  nbinsy,-nedgey,nedgey
	  );
  fSameEvent2D  = ns;
  fMultiEvent2D = nm;
}

void TCorrelation::FillCorrelation()
{
  TH2F ndiv("nc_"+fVar1+"_"+fVar2,
	    "N_{s}/N_{m}(#Delta#"+fVar2+", #Delta#"+fVar1+");#Delta#"+fVar2+";#Delta#"+fVar1,
	    nbinsx,-nedgex,nedgex,
	    nbinsy,-nedgey,nedgey);
  float ratio;
  float settle_at = 1.0;
  for ( int x = 1; x < nbinsx+1; x++ ) {
    for ( int y = 1; y < nbinsy+1; y++ ) {
      float nsvalue = fSameEvent2D.GetBinContent(x,y);
      float nmvalue = fMultiEvent2D.GetBinContent(x,y);
      if ( nsvalue < 500.0 || nmvalue < 500.0 ) {
	ndiv.SetBinContent(x,y,settle_at);	
      } else {
	ndiv.SetBinContent(x,y,nsvalue/nmvalue);
      }   
    }
  }
  fCorrelation2D = ndiv;
}
