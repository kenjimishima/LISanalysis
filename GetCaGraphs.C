#include "include.h"

Double_t draw_mass_min = 36;
Double_t draw_mass_max = 52;

/**
 * Draw nonlinear mass axis using TOFtoMass()
 *
 * @param h          TH1 (x-axis = TOF)
 * @param tof_40Ca   reference TOF for 40Ca
 */
void DrawMassAxis(TH1* h, Double_t tof_40Ca, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax)
{
  if (!h) return;
  std::vector<Double_t> masses = {38, 40, 42, 44, 46, 48};
  gPad->Update();  // ★ 必須

  // TOF range
  Double_t tof_min = h->GetXaxis()->GetXmin();
  Double_t tof_max = h->GetXaxis()->GetXmax();
  // Pad geometry
  Double_t y_top   = gPad->GetUymax();
  Double_t y_range = gPad->GetUymax() - gPad->GetUymin();
  // Mass range (axis frame only)
  Double_t mass_min = TOFtoMass(tof_min, tof_40Ca, 0);
  Double_t mass_max = TOFtoMass(tof_max, tof_40Ca, 0);

  // ---- Manual ticks ----
  Double_t tick_len     = 0.20 * ymax;
  Double_t label_offset = 0.30 * ymax;
  for (Double_t m_target : masses) {
    // --- invert TOFtoMass numerically ---
    Double_t tof_tick = -1;
    const Int_t Nscan = 5000;
    for (Int_t i = 0; i < Nscan - 1; ++i) {
      Double_t tof1 = tof_min + (tof_max - tof_min) * i / Nscan;
      Double_t tof2 = tof_min + (tof_max - tof_min) * (i + 1) / Nscan;
      Double_t m1 = TOFtoMass(tof1, tof_40Ca, 0);
      Double_t m2 = TOFtoMass(tof2, tof_40Ca, 0);
      if ((m1 - m_target) * (m2 - m_target) <= 0) {
	tof_tick = 0.5 * (tof1 + tof2);
	break;
      }
    }

    if (tof_tick < tof_min || tof_tick > tof_max) continue;
    // Tick
    TLine* l = new TLine(
			 tof_tick, ymax,
			 tof_tick, ymax + tick_len
			 );
    l->Draw();
    // Label
    TLatex* t = new TLatex(
			   tof_tick,
			   ymax + tick_len + label_offset,
			   Form("%.0f", m_target)
			   );
    t->SetTextAlign(22);
    t->SetTextSize(0.035);
    t->SetTextFont(42);
    t->Draw();
    // Title
    TLatex* titl = new TLatex(
			   0.1*xmin + 0.9*xmax,
			   ymax*2.,
			   "Mass number [e/M]"
			   );
    titl->SetTextAlign(22);
    titl->SetTextSize(0.04);
    titl->SetTextFont(42);
    titl->Draw();
  }
}

Bool_t ReadCa40PeakParam(const std::string &envpath,
			 Double_t &A,
			 Double_t &mean,
			 Double_t &sigma)
{
  TEnv env;
  Int_t status = env.ReadFile(envpath.c_str(), kEnvLocal);
  if (status < 0) {
    std::cerr << "Failed to read " << envpath << std::endl;
    return kFALSE;
  }
  A     = env.GetValue("Ca40PeakTOF.A",     -1.);
  mean  = env.GetValue("Ca40PeakTOF.Mean",  -1.);
  sigma = env.GetValue("Ca40PeakTOF.Sigma", -1.);
  return kTRUE;
}

//void GetCaGraphs(const char* basename = "RUN45_Spatial_40Ca_Beamoff")
//void GetCaGraphs(const char* basename = "RUN51_Spatial_40Ca_Beamoff")
void GetCaGraphs(const char* basename = "RUN52_Spatial_Beamoff_550")
{
  const char* histname = "h2_scaled";
  std::string base = strip_ext_and_dir(basename);
  std::string infile = base + "_Scaled.root";
  std::string indir  = "./root/scaled/";
  std::string envdir  = "./results/";
  std::string outdir  = "./root/ca_graph/";
  std::string outpath = outdir +  base + "_CaGraph.root";
  std::string textdir  = "./text/ca_graph/";
  std::string figdir  = "./outputs/ca_graph/";
  std::string figpath_pdf = figdir + base + ".pdf";
  std::string figpath_png = figdir + base + ".png";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);
  if (gSystem->AccessPathName(textdir.c_str())) gSystem->mkdir(textdir.c_str(), true);

  TH2D* h2 = GetTH2Dfile(infile.c_str(), indir.c_str(), histname);
  h2->Draw();
  // --- X, Y ビン数と範囲を取得 ---
  Int_t nx = h2->GetNbinsX();
  Int_t ny = h2->GetNbinsY();
  Double_t dsigma = 4.;  //Integral region of Ca peaks in sigma

  //--------------------------------------------------------------
  // データ格納ベクトル
  //--------------------------------------------------------------  
  std::vector<Double_t> yval;
  std::vector<Double_t> counts[num_ca]; // A, B, C
  std::vector<Double_t> errors[num_ca];

  //--------------------------------------------------------------
  // 各Yビンで処理
  //--------------------------------------------------------------  
  for (Int_t iy = 1; iy <= ny; ++iy) {
    TH1D* h1 = h2->ProjectionX(Form("h1_Y%d", iy), iy, iy);
    h1->SetTitle("");
    if (!h1) {
      cerr << "Can not read" << h1->GetName()<<endl;
      exit(0);
    }
    h1->Rebin(n_rebin);
    Double_t ycenter = h2->GetYaxis()->GetBinCenter(iy);
    yval.push_back(ycenter);
    TCanvas* canv = new TCanvas("c","c");
    canv->SetLogy();
    h1->Draw("eh");
    
    //--------------------------------------------------------------
    // Get parameters
    //--------------------------------------------------------------  
    std::string envname = envdir + base + "_slice_Y" + std::to_string(iy) + "_PeakFit.env";
    Double_t Ca40_A, Ca40_mean, Ca40_sigma;
    ReadCa40PeakParam(envname, Ca40_A, Ca40_mean, Ca40_sigma);
    Double_t peak_center  = Ca40_mean;
    Double_t sigma = Ca40_sigma;
    cout <<"peak center ="<< peak_center <<", sigma = "<< sigma<<endl;

    // --- 領域定義 ---
    Double_t ca_centers[num_ca];
    Double_t ca_xmin[num_ca];
    Double_t ca_xmax[num_ca];
    Double_t integ[num_ca];
    Double_t err[num_ca];

    std::vector<TLine*> lines;
    Double_t xmin = MasstoTOF(draw_mass_min, peak_center);
    Double_t xmax = MasstoTOF(draw_mass_max, peak_center);
    Double_t ymin = 0.1;
    Double_t ymax = h1->GetMaximum() * 5.;
    
    // --- 積分と誤差 ---
    for (int i = 0; i < num_ca; ++i) {
      ca_centers[i] = MasstoTOF(ca_mass[i], peak_center);
      ca_xmin[i] = ca_centers[i] - dsigma*sigma;
      ca_xmax[i] = ca_centers[i] + dsigma*sigma;
      Int_t bmin = h1->FindBin(ca_xmin[i]);
      Int_t bmax = h1->FindBin(ca_xmax[i]);

      integ[i] = h1->IntegralAndError(bmin, bmax, err[i]);
      counts[i].push_back(integ[i]);
      errors[i].push_back(err[i]);

      TLine* line_min = new TLine(ca_xmin[i], ymin, ca_xmin[i], ymax);
      TLine* line_max = new TLine(ca_xmax[i], ymin, ca_xmax[i], ymax);
      lines.push_back(line_min);
      lines.push_back(line_max);
    }
    h1->GetXaxis()->SetRangeUser(xmin,xmax);
    h1->GetYaxis()->SetRangeUser(ymin,ymax);
    // --- ヒストグラムを描画してから、線を重ねる ---
    for (auto line : lines) {
      line->SetLineColor(kRed);
      line->SetLineStyle(2);
      line->Draw("same");
    }

    DrawMassAxis(h1,peak_center,xmin,xmax,ymin,ymax);
    canv->SetGrid();
    canv->SaveAs(Form("%s%s_Y%d.png",figdir.c_str(),base.c_str(),iy));
    delete h1;
    delete canv;
  }

  //--------------------------------------------------------------
  // TGraphErrors配列で作成
  //--------------------------------------------------------------
  TGraphErrors* gr[num_ca];
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(Form("%s; YAG laser position [mm]; Peak count",base.c_str()));
  TLegend* leg = new TLegend(0.7, 0.8, 0.95, 0.90,"");
  Int_t n = yval.size();

  for (int i = 0; i < num_ca; ++i) {
    gr[i] = new TGraphErrors(n, &yval[0], &counts[i][0], 0, &errors[i][0]);
    gr[i]->SetName(Form("gr_ca%0.0f",ca_mass[i]));
    gr[i]->SetTitle(Form("%0.0fCa; YAG laser position [mm]; Peak count",ca_mass[i]));
    gr[i]->SetMarkerColor(i+2);
    gr[i]->SetLineColor(i+2);
    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerSize(1.);
    mg->Add(gr[i]);
    leg->AddEntry(gr[i], gr[i]->GetTitle(), "pl");
    //    gr[i]->Print();
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->SetGrid();
  c2->SetLogy();
  mg->Draw("APL");
  mg->GetYaxis()->SetRangeUser(0.1,1e6);
  leg->Draw();
  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  c2->Update();
  c2->SaveAs(figpath_png.c_str());
  c2->Update();
  c2->SaveAs(figpath_pdf.c_str());
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  for (int i = 0; i < num_ca; ++i) {
    gr[i]->Write();
    TString fname = Form("%s%s_%s.txt",textdir.c_str(), basename, gr[i]->GetName());
    SaveTGraphErrors(gr[i], fname.Data());    
  }
  fout->Close();
  //  fin->Close();
  return;
}



