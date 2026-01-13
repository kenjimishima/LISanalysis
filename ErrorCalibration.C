#include "include.h"

void FitMultiGauss(TH1D* h1, const char* envbasename)
{
  if (!h1) { ::Error("FitMultiGauss", "Null histogram pointer."); return; }

  h1->Rebin(n_rebin);
  //  double xmin = 3.10;
  //  double xmax = 3.52;
  //  double xmax = 4.09;
  //  double mass[] = {24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
  //  double mass[] = {24,25,26,27,28,29};
  double xmin = 2.20;
  double xmax = 2.45;
  double mass[] = {12};
  
  int ngaus = sizeof(mass) / sizeof(mass[0]);
  
  // --- 関数式生成 ---
  TString formula = "[0]"; // pol0 baseline
  for (int i = 0; i < ngaus; ++i)
    formula += Form(" + gaus(%d)", 1 + 3*i); // 3パラ/gaus

  // 全パラ数 = 1 + 3*ngaus
  TF1* fitfunc = new TF1("fitfunc", formula, xmin, xmax);
  fitfunc->SetNpx(1000);
  // --- 初期値設定 ---
  double A0 = h1->GetMaximum();
  fitfunc->SetParameter(0, h1->GetMinimum()); // baseline
  for (int i = 0; i < ngaus; ++i) {
    double mean  = MasstoTOF(mass[i],def_tof_40Ca);
    cout << mass[i]<<" "<<mean<<endl;
    fitfunc->SetParameter(1 + 3*i + 0, A0*0.1); // amplitude
    fitfunc->SetParameter(1 + 3*i + 1, mean);   // mean
    fitfunc->SetParameter(1 + 3*i + 2, def_sigma);  // sigma

    fitfunc->SetParLimits(1 + 3*i + 0, A0*0., A0*1.);
    fitfunc->SetParLimits(1 + 3*i + 1, mean * 0.98, mean * 1.02);
    fitfunc->SetParLimits(1 + 3*i + 2, def_sigma * 0.7, def_sigma * 1.5);
  }

  // --- フィット ---
  TCanvas* c2 = new TCanvas("c2", "Multi-Gaussian Fit", 900, 700);
  // ---- 描画 ----
  h1->Print();
  h1->Draw("EH");   
  h1->Fit(fitfunc, "R");
  c2->Update();

  //refit with 3 sigma
  double mean = fitfunc->GetParameter(2);
  double sigma = fitfunc->GetParameter(3);
  TF1 *f = h1->GetFunction("fitfunc");
  if (f) {
    h1->GetListOfFunctions()->Remove(f);
    delete f;
  }

  //  fitfunc->SetRange(mean-3.*sigma, mean+3.*sigma);
  h1->Fit(fitfunc, "R+");
  fitfunc->Print();
  c2->Update();

  std::cout << mean<<" "<<sigma <<  std::endl;
  std::cout << "fitfunc ptr = " << fitfunc << std::endl;
  std::cout << "hist func ptr = "
	    << h1->GetFunction("fitfunc") << std::endl;
  return;

  double chi2 = fitfunc->GetChisquare();
  double ndf  = fitfunc->GetNDF();
  double chi2ndf = chi2 / ndf ;
  double base_voltage = fitfunc->GetParameter(0)/n_rebin;

  // --- 結果表示 ---
  std::cout << "chi2 = " << chi2 << "\n";
  std::cout << "ndf  = " << ndf  << "\n";
  std::cout << "chi2/ndf = " << chi2ndf << "\n";
  std::cout << "base voltage = " << base_voltage << std::endl;

  // --- TEnvに保存 ---
  TString outenv = Form("./results/%s.env", envbasename);
  TEnv env;
  env.SetValue("FitResult.chi2", chi2);
  env.SetValue("FitResult.ndf", ndf);
  env.SetValue("FitResult.chi2_over_ndf", chi2ndf);
  env.SetValue("FitResult.count_per_mV", 1./chi2ndf);
  env.SetValue("FitResult.base_voltage", base_voltage);
  env.WriteFile(outenv);
  std::cout << "Saved fit results to: " << outenv << std::endl;

  gPad->SetLogy();
  gPad->SetGrid();
  // --- 保存 ---
  TString outpng = Form("./outputs/ErrorCalib/%s_fit%dG.png", h1->GetName(), ngaus);
  TString outpdf = Form("./outputs/ErrorCalib/%s_fit%dG.pdf", h1->GetName(), ngaus);
  TString outroot = Form("./outputs/ErrorCalib/%s_fit%dG.root", h1->GetName(), ngaus);
  c2->SaveAs(outpng);
  c2->SaveAs(outpdf);
  c2->SaveAs(outroot);
  std::cout << "Saved fit result:\n  " << outpng << "\n  " << outpdf << std::endl;
}

void ErrorCalibration(const char* basename = "RUN52_Spatial_Beamoff_550.txt",
		      //void ErrorCalibration(const char* basename = "RUN45_Spatial_40Ca_Beamoff",
		      //void ErrorCalibration(const char* basename = "RUN51_Spatial_40Ca_Beamoff",
		      const int sliceIndex = 1,
		      const char* histname = "h2_subtracted"
		      )
{
  std::string base = strip_ext_and_dir(basename);
  std::string indir  = "./root/baseline/";
  std::string infile  = base + "_Baseline";
  std::string envbasename = base + "_slice_Y" + std::to_string(sliceIndex);
  std::string figdir  = "./outputs/ErrorCalib/";
  std::string envdir  = "./results/";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(envdir.c_str())) gSystem->mkdir(envdir.c_str(), true);
  
  TCanvas* c0 = new TCanvas("c0", "TH2 Histogram", 900, 700);
  TH2D* h2 = GetTH2Dfile(infile.c_str(), indir.c_str(), histname);
  h2->Draw("colz");
  
  // --- 最初のYビンを投影（ProjectionX） ---
  TCanvas* c1 = new TCanvas("c1", "TH2D Slice Viewer", 800, 600);
  TH1D* h1 = h2->ProjectionX(Form("%s_Y%d", h2->GetName(), sliceIndex),sliceIndex, sliceIndex);
  if (!h1) {
    Error("PlotFirstSliceFromTH2D", "Projection failed or empty slice (Y-bin=%d)", sliceIndex);
    return;
  }
  h1->SetEntries(h1->GetNbinsX()); //force to input entry
  h1->Print();
  
  // --- 描画設定 ---
  h1->SetTitle(Form("%s: X-slice at Y-bin %d", h2->GetName(), sliceIndex));
  h1->GetXaxis()->SetTitle("TOF [us]");
  h1->GetYaxis()->SetTitle("Voltage [mV]");
  std::cout << "min=" << h1->GetMinimum() << ", max=" << h1->GetMaximum() << std::endl;
  // --- フィット実行 ---
  h1->Draw("EH");
  FitMultiGauss(h1, envbasename.c_str());
  //  fin->Close();
}
