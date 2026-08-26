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
  std::string envdir  = "./results/";
  std::string envpath = envdir + base + "_baseline.env";
  std::string figdir  = "./outputs/baseline/";
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);
  if (gSystem->AccessPathName(envdir.c_str())) gSystem->mkdir(envdir.c_str(), true);
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);

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

    /*
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
      //      h2_sub->SetBinError(ix, iy, 0.1*TMath::Sqrt(TMath::Abs(val))); //適当な誤差
    }
  }
    */
  //for save
  TEnv env;
  //  env << std::setprecision(17);
  Double_t baseline = 0.;
  Double_t baseline_err = 0.;
  Double_t baseline_sigma = 0.;

  // --- 各Yスライスでベースライン平均を求める ---
  for (Int_t iy = 1; iy <= ny; ++iy) {
    // iy番目のYスライスをX方向の1次元ヒストグラムにする
    TH1D* hpx = h2->ProjectionX(Form("hpx_baseline_iy%d", iy), iy, iy);    
    // pol0 = 定数関数。ベースラインを横一直線でフィットする
    TF1* fpol0 = new TF1(Form("fpol0_iy%d", iy), "pol0", baseline_xmin,  baseline_xmax);
    TCanvas* c1 = new TCanvas("c1", "Baseline", 800, 600);
    hpx->Draw("h");
    hpx->Fit(fpol0, "WR");
    TString outpng = Form("%s/%s_Y%d.png", figdir.c_str(), basename, iy);
    c1->SaveAs(outpng);
    // fit結果からベースライン中心値と誤差を取り出す
    baseline = fpol0->GetParameter(0);
    baseline_err = fpol0->GetParError(0);
    Double_t chi2 = fpol0->GetChisquare();
    Int_t ndf = fpol0->GetNDF();
    baseline_sigma = TMath::Sqrt(chi2 / ndf);
    env.SetValue(Form("Baseline.Y%d",iy), baseline);
    env.SetValue(Form("Baseline.Err.Y%d",iy), baseline_err);
    env.SetValue(Form("Baseline.Chi2.Y%d",iy), chi2);
    env.SetValue(Form("Baseline.NDF.Y%d",iy), ndf);
    env.SetValue(Form("Baseline.Sigma.Y%d",iy), baseline_sigma);
    
    for (Int_t ix = 1; ix <= nx; ++ix) {
      Double_t val = h2->GetBinContent(ix, iy) - baseline;
      h2_sub->SetBinContent(ix, iy, val);
      h2_sub->SetBinError(ix, iy, TMath::Sqrt(TMath::Abs(val)));
      //      h2_sub->SetBinError(ix, iy, baseline_sigma); 
    }
    delete fpol0;
    delete hpx;
    delete c1;
  }

  // --- 出力ファイルに保存 ---
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  h2_sub->Write();
  fout->Close();
  env.WriteFile(envpath.c_str());
  std::cout << "Saved env file: " << envpath << std::endl;
  std::cout << "Saved baseline-subtracted histogram to " << outpath << std::endl;
  return;
}
