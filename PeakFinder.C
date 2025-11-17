#include "include.h"

double xmin = 3.10;
double xmax = 4.09;
//double mass[] = {24,25,26,27,28,29,36,37,38,39,40,42,43,44,48,49,50,56};
//double mass[] = {12,24,25,26,27,28,29,36,37,38,39,40,42,43,44,48,49,50,56};
//double mass[] = {12,24,26,27,28,36,37,38,39,40,42,44,48,51,52,57};
double mass[] = {12,24,26,27,28,36,37,38,39,40,42,44,48,57};

TGraphErrors* SubtractTF1FromGraph(const TGraphErrors* gr, const TF1* f)
{
    Int_t n = gr->GetN();
    TGraphErrors* gr_sub = new TGraphErrors(n);
    gr_sub->SetName(Form("%s_sub", gr->GetName()));
    gr_sub->SetTitle(Form("%s - %s", gr->GetTitle(), f->GetName()));

    for (Int_t i = 0; i < n; i++) {
        Double_t x, y;
        gr->GetPoint(i, x, y);
        Double_t ey = gr->GetErrorY(i);

        Double_t fval = f->Eval(x);
        Double_t y_sub = y - fval;

        gr_sub->SetPoint(i, x, y_sub);
        gr_sub->SetPointError(i, 0, ey); // x誤差は0と仮定
    }

    gr_sub->SetMarkerStyle(20);
    gr_sub->SetMarkerColor(kRed);
    return gr_sub;
}


//======================================================================
// gr1 - gr2 の差分グラフを作成する
//======================================================================
TGraphErrors* MakeDiffGraph(const TGraphErrors* gr1, const TGraphErrors* gr2, const char* envbasename)
{
  if (!gr1 || !gr2) {
    Error("MakeDiffGraph", "Null input graph(s)");
    return nullptr;
  }
  Int_t n = gr1->GetN();
  TGraphErrors* gr_diff = (TGraphErrors*) gr1->Clone("gr_diff");
  gr_diff->SetTitle(Form("%s - %s; mass - 40Ca; Diff [us]", gr1->GetName(), gr2->GetName()));

  for (Int_t i = 0; i < n; i++) {
    Double_t x1, y1, x2, y2;
    gr1->GetPoint(i, x1, y1);
    gr2->GetPoint(i, x2, y2);

    Double_t ey1 = gr1->GetErrorY(i);
    Double_t ey2 = gr2->GetErrorY(i);

    Double_t diff  = y1 - y2;
    Double_t ediff = TMath::Sqrt(ey1*ey1 + ey2*ey2);
    gr_diff->SetPoint(i, x1 - mass_40Ca, diff);
    gr_diff->SetPointError(i, 0., ediff);
  }

  gr_diff->SetMarkerStyle(20);

  TCanvas* c3 = new TCanvas("c3", "Diff", 800, 600);
  gr_diff->Draw("AP");
  gPad->SetGrid();
  TF1* fitfunc = new TF1("fitfunc","pol1");
  gr_diff->Fit("fitfunc");
  double p0 = fitfunc->GetParameter(0);
  double p1 = fitfunc->GetParameter(1);

  // --- TEnvに保存 ---
  TString outenv = Form("./results/%s.env", envbasename);
  TEnv env;
  env.SetValue("TOFCorrection.p0", p0);
  env.SetValue("TOFCorrection.p1", p1);
  env.WriteFile(outenv);
  std::cout << "Saved fit results to: " << outenv << std::endl;
  
  return gr_diff;
}


void PeakFinder(const char* basename = "RUN45_Spatial_40Ca_Beamoff",
		const int sliceIndex = 1,
		const char* histname = "h2_subtracted"
		)
{
  std::string base = strip_ext_and_dir(basename);
  std::string indir  = "./root/baseline/";
  std::string infile  = base + "_Baseline";
  std::string outdir  = "./root/peak_finder/";
  std::string outpath = outdir +  base + "_PeakFinder.root";
  std::string envbasename = base + "_slice_Y" + std::to_string(sliceIndex) + "_PeakFinder";
  std::string envdir  = "./results/";
  std::string figdir  = "./outputs/peak_finder/";
  std::string figpath_pdf = figdir + base + "_slice_Y" + std::to_string(sliceIndex) + ".pdf";
  std::string figpath_png = figdir + base + "_slice_Y" + std::to_string(sliceIndex) + ".png";
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
  h1->GetYaxis()->SetTitle("Voltage [mV]");
  std::cout << "min=" << h1->GetMinimum() << ", max=" << h1->GetMaximum() << std::endl;

  TGraphErrors *gr = new TGraphErrors();
  gr->SetTitle("Mass number v.s TOF; Mass number;TOF [us];");
  TGraphErrors *gr_calc = new TGraphErrors();
  gr->SetTitle("Mass number v.s TOF; Mass number;TOF [us];");

  /*  TOFtoMass(4.86,def_tof_40Ca);
  TOFtoMass(2.29,def_tof_40Ca);
  MasstoTOF(12,def_tof_40Ca,0); 
  MasstoTOF(12,def_tof_40Ca,1); 
  */
  
  // --- フィット実行 ---
  int ngaus = sizeof(mass) / sizeof(mass[0]);
  for (int i = 0; i < ngaus; i++) {
    double mean  = MasstoTOF(mass[i],def_tof_40Ca,0); //default 
    //    double corr_mean  = MasstoTOF(mass[i],def_tof_40Ca,0); //corrected
    double xmin_single = mean * 0.992;
    double xmax_single = mean * 1.008;
    auto [peak_bin, peak_val] = FindMaxBinAndValueInRange(h1, xmin_single, xmax_single);
    Double_t peak_center = h1->GetXaxis()->GetBinCenter(peak_bin);
    //   cout << mass[i]<<" "<<mean<<" "<<peak_center<<" "<< peak_val<<endl;
    TF1* fitfunc = new TF1("fitfunc", "gaus", xmin_single, xmax_single);
    fitfunc->SetParameter(0, peak_val); // amplitude
    fitfunc->SetParameter(1, peak_center);   // mean
    fitfunc->SetParameter(2, def_sigma);  // sigma
    fitfunc->SetParLimits(0, peak_val*0.5, peak_val*2.0);
    fitfunc->SetParLimits(1, peak_center*0.99, peak_center*1.01);
    fitfunc->SetParLimits(2, def_sigma * 0.5, def_sigma * 1.5);
    h1->GetXaxis()->SetRangeUser(xmin_single,xmax_single);
    h1->Fit(fitfunc, "RQ");
    c1->SaveAs(Form("%s/%s_PeakFit_Yslice%d_%0.0f.png",figdir.c_str(),base.c_str(),sliceIndex, mass[i]));
    double p0  = fitfunc->GetParameter(0);
    double p1  = fitfunc->GetParameter(1);
    double p2  = fitfunc->GetParameter(2);
    double ep0 = fitfunc->GetParError(0);
    double ep1 = fitfunc->GetParError(1);
    double ep2 = fitfunc->GetParError(2);
    cout <<p0<<" +/-"<<ep0<<" : "<<p1<<" +/-"<<ep1<<" : "<<p2<<" +/-"<<ep2<<endl;
    gr->SetPoint(i, mass[i], p1);
    gr->SetPointError(i, 0., ep1);
    //    gr_calc->SetPoint(i, mass[i], mean);
    gr_calc->SetPoint(i, mass[i], mean);
    gr_calc->SetPointError(i, 0., 0);
  }

  TCanvas* c2 = new TCanvas("c2", "Peak position", 800, 600);
  c2->SetGrid();
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.);
  gr->Draw("AP");
  gr->GetYaxis()->SetRangeUser(0,8);
  gr->GetXaxis()->SetLimits(0,60);
  //  gr_calc->Draw("Lsame");
  //  gr_calc->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.2, 0.7, 0.4, 0.85,"");
  leg->AddEntry(gr, "Peak fit", "pl");
  leg->AddEntry(gr_calc, "Calculation", "l");
  leg->Draw();

  double min_mass = *std::min_element(mass, mass + ngaus);
  double max_mass = *std::max_element(mass, mass + ngaus);
  TF1* func = new TF1("sqrtfunc","[0]*TMath::Sqrt(x)+[1]",0., max_mass+1);
  func->SetParameter(0,1.);
  func->SetParameter(1,0.1);
  // func->SetParameter(2,1.);
  gr->Fit("sqrtfunc","R+");
  func->SetLineColor(4);

  //  MakeDiffGraph(gr, gr_calc, envbasename.c_str());
  c2->Update();
  c2->SaveAs(figpath_pdf.c_str());
  c2->Update();
  c2->SaveAs(figpath_png.c_str());

  //
  TCanvas* c3 = new TCanvas("c3", "Diff", 800, 600);
  TGraphErrors* gr_subt = SubtractTF1FromGraph(gr, func);
  gr_subt->Draw("AP");

  Int_t draw_min = 10;
  Int_t draw_max = 60;
  const Int_t num = draw_max - draw_min;
  TString h_mass_name = TString(h1->GetName()) + "_mass";
  TH1D* h_outgas = new TH1D(h_mass_name,"Mass of Ions; M/e; Count [1/mass]",num, draw_min, draw_max);
  // --- vectorで動的配列を確保 ---
  std::vector<Double_t> border(num + 1);  
  // --- 計算 ---

  c1->cd();
  h1->SetStats(0);
  h1->Draw("eh");
  h1->GetXaxis()->UnZoom();
  gPad->SetLogy();
  std::vector<TLine*> lines;
  border[0] = func->Eval(draw_min - 0.5);
  Double_t ymin = 1.e-2;
  Double_t ymax = 1.e+2;
  TLine* l0 = new TLine(border[0], ymin, border[0], ymax);
  
  for (Int_t i = 0; i < num; i++) {
    Double_t err;
    border[i+1] = func->Eval(draw_min + i + 0.5);  
    Int_t minbin = h1->FindBin(border[i]);
    Int_t maxbin = h1->FindBin(border[i+1]);
    Double_t count_in_mass = h1->IntegralAndError(minbin, maxbin, err);
    h_outgas->SetBinContent(i + 1, count_in_mass);
    h_outgas->SetBinError(i + 1, err);
    //    cout << i <<" mass = "<< draw_min + i  <<", Count in " << border[i] << " - " << border[i+1]
    //       << " = " << count_in_mass << " +/- " << err << endl;

    TLine* l = new TLine(border[i+1], ymin, border[i+1], ymax);
    l->Draw("same");    
    lines.push_back(l);   // 必要なら保持
  }  
  
  TCanvas* c4 = new TCanvas("c4", "Mass counts", 800, 600);
  h_outgas -> Draw("EH");
  gPad->SetGrid();
  
  
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  gr->Write();
  gr_calc->Write();
  h1->Write();
  h_outgas->Write();
  fout->Close();

  //  fin->Close();
}

