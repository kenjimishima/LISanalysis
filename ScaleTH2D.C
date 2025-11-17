#include "include.h"
#include "ErrorCalibration.C"

void ScaleTH2D(const char* basename="RUN45_Spatial_40Ca_Beamoff",
	       const char* calibname="RUN45_Spatial_40Ca_Beamoff",
	       const int sliceIndex = 1,
	       const char* histname = "h2_subtracted")
{

  std::string base  = strip_ext_and_dir(basename);
  std::string calib = strip_ext_and_dir(calibname);
  std::string inpath  = "./root/baseline/" + base + "_Baseline.root";
  std::string envpath = "./results/" + calib + "_slice_Y" + std::to_string(sliceIndex) + ".env";
  std::string save_envbasename = base + "_slice_Y" + std::to_string(sliceIndex) + "_Scaled";
  std::string outdir  = "./root/scaled/";
  std::string outpath = outdir + base + "_Scaled.root";  
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);

  
  //----------------------------------------------------------------
  // 1. 入力ROOTファイルを開く
  //----------------------------------------------------------------
  TFile* fin = TFile::Open(inpath.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    ::Error("ScaleTH2DwithEnv", "Cannot open input file: %s", inpath.c_str());
    return;
  }

  //----------------------------------------------------------------
  // 2. ヒストグラムを取得
  //----------------------------------------------------------------
  TH2D* h2 = dynamic_cast<TH2D*>(fin->Get(histname));
  if (!h2) {
    ::Error("ScaleTH2DwithEnv", "Histogram '%s' not found in %s", histname, inpath.c_str());
    fin->Close();
    return;
  }

  //----------------------------------------------------------------
  // 3. TEnvからスケール係数を取得
  //----------------------------------------------------------------
  TEnv env;
  Int_t status = env.ReadFile(envpath.c_str(), kEnvLocal);
  if (status < 0) {
    std::cerr << "Failed to read " << envpath.c_str() << std::endl;
    return;
  }
  
  Double_t chi2 = env.GetValue("FitResult.chi2", 0.0); 
  Double_t ndf  = env.GetValue("FitResult.ndf", 0.0); 
  Double_t chi2_over_ndf  = env.GetValue("FitResult.chi2_over_ndf", 0.0); 
  Double_t count_per_mV  = env.GetValue("FitResult.count_per_mV", 0.0); 
  Double_t base_voltage = env.GetValue("FitResult.base_voltage", -999999.);

  double scale = 1./chi2_over_ndf;
  std::cout << "Base voltage from env: " << base_voltage << std::endl;
  std::cout << "Count per mV: " << count_per_mV << std::endl;

  //----------------------------------------------------------------
  // 4. 新しいヒストグラムを作成してスケール変換
  //----------------------------------------------------------------
  TH2D* h2_scaled = (TH2D*)h2->Clone("h2_scaled");
  h2_scaled->Reset();

  int nx = h2->GetNbinsX();
  int ny = h2->GetNbinsY();

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double val = h2->GetBinContent(ix, iy);
      double newval = (val-base_voltage) * count_per_mV;
      double err = std::sqrt(std::abs(newval)); // sqrt(N)
      h2_scaled->SetBinContent(ix, iy, newval);
      h2_scaled->SetBinError(ix, iy, err);
    }
  }

  h2_scaled->SetTitle(Form("%s (scaled by %.3f)", histname, scale));

  // --- 最初のYビンを投影（ProjectionX） ---
  TH1D* h1 = h2_scaled->ProjectionX("h1_scaled", sliceIndex, sliceIndex);
  h1->Draw("EH");
  cout << save_envbasename<<endl;
  FitMultiGauss(h1, save_envbasename.c_str());

  //----------------------------------------------------------------
  // 5. 出力ROOTファイルに保存
  //----------------------------------------------------------------
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  h2_scaled->Write("h2_scaled");
  fout->Close();

  std::cout << "Saved scaled histogram to: " << outpath << std::endl;
  //  fin->Close();
}
