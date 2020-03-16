#include "Root/Plotter.cc"

void multi_plotter_2pcs()
{
  std::vector <TString> modes =
    {
     "pp",
     "pm",
     "zhpp",
     "zhpm",
     "zh0002pp",
     "zh0002pm",
     "zh0203pp",
     "zh0203pm",
     "zh0305pp",
     "zh0305pm",     
    };
	 
  std::vector <TString> targets =
    {
     "C",
     "Fe",
     "Pb",
     "D"
    };
  
  TString input_dir = "/home/luciano/Physics/CLAS/pion_correlation/";
  
  TH1F *h;
  std::vector<TH1F *>target_th1f;

  TFile *f = new TFile();
  new TCanvas();
  for ( TString i: modes ) {
    target_th1f.clear();
    f->Open(input_dir+i+"/alltarg.root");
    int target_counter=0;
    for ( TString k: targets ) {
      h = (TH1F*)gDirectory->Get("nc_phiPQboosted_"+k);
      h->GetEntries();
      h->Draw("E0X0 same");
      target_th1f.push_back(h);
    }
    TString out_filename = input_dir+i+"/"+i+".pdf";
    cout << out_filename << endl;
    Plotter(target_th1f[0],target_th1f[1],target_th1f[2],target_th1f[3],out_filename);
  }

  
}
