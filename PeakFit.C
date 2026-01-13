#include "include.h"

void PeakFit(const char* basename = "RUN45_Spatial_40Ca_Beamoff",
//void PeakFit(const char* basename = "RUN51_Spatial_40Ca_Beamoff",
	     const int sliceIndex = 2)
{
  const char* histname = "h2_scaled";
  std::string base = strip_ext_and_dir(basename);
  std::string indir  = "./root/scaled/";
  std::string infile  = base + "_Scaled";
  std::string outdir  = "./root/peak_fit/";
  std::string outpath = outdir +  base + "_PeakFit.root";
  std::string envbasename = base + "_slice_Y" + std::to_string(sliceIndex) + "_PeakFit";
  std::string envdir  = "./results/";
  std::string figdir  = "./outputs/peak_fit/";
  std::string figpath_pdf = figdir + base + "_slice_Y" + std::to_string(sliceIndex) + ".pdf";
  std::string figpath_png = figdir + base + "_slice_Y" + std::to_string(sliceIndex) + ".png";
  std::string figpath_root= figdir + base + "_slice_Y" + std::to_string(sliceIndex) + ".root";
  if (gSystem->AccessPathName(envdir.c_str())) gSystem->mkdir(envdir.c_str(), true);
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);
  
  TCanvas* c1 = new TCanvas("c1", "TH2D Slice Viewer", 800, 600);
  TH2D* h2 = GetTH2Dfile(infile.c_str(), indir.c_str(), histname);
  // --- IndexのYビンを投影（ProjectionX） ---
  TH1D* h1 = h2->ProjectionX(Form("%s_Y%d", h2->GetName(), sliceIndex),sliceIndex, sliceIndex);
  if (!h1) {
    Error("PlotFirstSliceFromTH2D", "Projection failed or empty slice (Y-bin=%d)", sliceIndex);
    return;
  }
  h1->SetEntries(h1->GetNbinsX()); //force to input entry
  h1->Rebin(n_rebin);
  // --- 描画設定 ---
  h1->Draw("EH");
  h1->SetTitle(Form("%s: X-slice at Y-bin %d", h2->GetName(), sliceIndex));
  h1->GetXaxis()->SetTitle("TOF [us]");
  h1->GetYaxis()->SetTitle("Count [1/bin]");
  gPad->SetLogy();
  gPad->SetGrid();
  std::cout << "min=" << h1->GetMinimum() << ", max=" << h1->GetMaximum() << std::endl;

  // --- フィット実行 ---
  double mean  = MasstoTOF(40, def_tof_40Ca, 0); //default 
  double xmin_single = mean * 0.97;
  double xmax_single = mean * 1.03;
  auto [peak_bin, peak_val] = FindMaxBinAndValueInRange(h1, xmin_single, xmax_single);
  Double_t peak_center = h1->GetXaxis()->GetBinCenter(peak_bin);
  xmin_single = peak_center - 3*def_sigma;
  xmax_single = peak_center + 3*def_sigma;
  //  cout <<  peak_bin << " "<< peak_val<<endl;
  cout << "Peak center = "<< peak_center <<" [us] "<<endl;
  cout << "Fitting region : "<<  xmin_single<<" - "<<xmax_single<< " [us]"<<endl;
  TF1* fitfunc = new TF1("fitfunc", "gaus", xmin_single, xmax_single);
  fitfunc->SetParameter(0, peak_val); // amplitude
  fitfunc->SetParameter(1, peak_center);   // mean
  fitfunc->SetParameter(2, def_sigma);  // sigma
  fitfunc->SetParLimits(0, peak_val*0.5, peak_val*2.0);
  fitfunc->SetParLimits(1, xmin_single, xmax_single);
  fitfunc->SetParLimits(2, def_sigma * 0.5, def_sigma * 1.5);
  //  h1->GetXaxis()->SetRangeUser(xmin_single,xmax_single);
  h1->Fit(fitfunc, "R");
  c1->SaveAs(figpath_pdf.c_str());
  c1->SaveAs(figpath_png.c_str());
  c1->SaveAs(figpath_root.c_str());
  double p0  = fitfunc->GetParameter(0);
  double p1  = fitfunc->GetParameter(1);
  double p2  = fitfunc->GetParameter(2);
  double ep0 = fitfunc->GetParError(0);
  double ep1 = fitfunc->GetParError(1);
  double ep2 = fitfunc->GetParError(2);
  //  cout <<p0<<" +/-"<<ep0<<" : "<<p1<<" +/-"<<ep1<<" : "<<p2<<" +/-"<<ep2<<endl;

  // --- Histogramを保存 ---
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  h1->Write();
  fout->Close();
  //  fin->Close();

  // --- TEnvに保存 ---
  TString outenv = Form("./results/%s.env", envbasename.c_str());
  TEnv env;
  env.SetValue("Ca40PeakTOF.A", p0);
  env.SetValue("Ca40PeakTOF.A.Error", ep0);
  env.SetValue("Ca40PeakTOF.Mean", p1);
  env.SetValue("Ca40PeakTOF.Mean.Error", ep1);
  env.SetValue("Ca40PeakTOF.Sigma", p2);
  env.SetValue("Ca40PeakTOF.Sigma.Error", ep2);
  env.WriteFile(outenv);
  std::cout << "Saved fit results to: " << outenv << std::endl;


  //  fin->Close();
}

