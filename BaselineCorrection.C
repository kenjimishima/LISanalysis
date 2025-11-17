#include "include.h"

void BaselineCorrection(const char* basename = "RUN45_Spatial_40Ca_Beamoff",
                            const char* histname = "h2"
			)
{
  std::string base = strip_ext_and_dir(basename);
  std::string indir  = "./root/rawdata/";
  std::string outdir  = "./root/baseline/";
  std::string inpath  = indir + base + ".root";
  std::string outpath = outdir +  base + "_Baseline.root";
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);

  TH2D* h2 = GetTH2Dfile(basename, indir.c_str(), histname);

  // --- 出力ヒストグラムを作成（コピー） ---
  TH2D* h2_sub = (TH2D*)h2->Clone(Form("%s_subtracted", histname));
  h2_sub->SetTitle(Form("%s (baseline subtracted)", histname));
  
  // --- X, Y ビン数と範囲を取得 ---
  Int_t nx = h2->GetNbinsX();
  Int_t ny = h2->GetNbinsY();
  Int_t bxmin = h2->GetXaxis()->FindBin(baseline_xmin);
  Int_t bxmax = h2->GetXaxis()->FindBin(baseline_xmax);
  
  std::cout << "Subtracting baseline from X = " << baseline_xmin
	    << " to " << baseline_xmax << " (" << bxmin << "–" << bxmax << " bins)" << std::endl;
  
  // --- 各Yスライスでベースライン平均を求める ---
  for (Int_t iy = 1; iy <= ny; ++iy) {
    Double_t sum = 0.0;
    Int_t ncount = 0;
    for (Int_t ix = bxmin; ix <= bxmax; ++ix) {
      sum += h2->GetBinContent(ix, iy);
      ++ncount;
    }
    Double_t baseline = 0;
    if (ncount != 0) baseline  = sum / ncount;
    
    // --- ベースラインを引く ---
    for (Int_t ix = 1; ix <= nx; ++ix) {
      Double_t val = h2->GetBinContent(ix, iy) - baseline;
      h2_sub->SetBinContent(ix, iy, val);
      h2_sub->SetBinError(ix, iy, TMath::Sqrt(TMath::Abs(val)));
    }
  }
  
  // --- 出力ファイルに保存 ---
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  h2_sub->Write();
  fout->Close();
  
  std::cout << "Saved baseline-subtracted histogram to " << outpath << std::endl;
  return;
}
