#include "include.h"

void GetCaGraphs(const char* basename = "RUN45_Spatial_40Ca_Beamoff",
//void GetCaGraphs(const char* basename = "RUN45_Spatial_44Ca_48Ca_Beamon48",
		 const char* histname = "h2_scaled"
		 )
{
  std::string base = strip_ext_and_dir(basename);
  std::string infile = base + "_Scaled.root";
  std::string indir  = "./root/scaled/";
  std::string outdir  = "./root/ca_graph/";
  std::string outpath = outdir +  base + "_CaGraph.root";
  std::string figdir  = "./outputs/ca_graph/";
  std::string figpath_pdf = figdir + base + ".pdf";
  std::string figpath_png = figdir + base + ".png";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);

  TH2D* h2 = GetTH2Dfile(infile.c_str(), indir.c_str(), histname);
  h2->Draw();
  // --- X, Y ビン数と範囲を取得 ---
  Int_t nx = h2->GetNbinsX();
  Int_t ny = h2->GetNbinsY();
  Double_t xmin = def_tof_40Ca * 0.992;
  Double_t xmax = def_tof_40Ca * 1.008;
  Double_t dsigma = 3.;  //Integral region of Ca peaks in sigma

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
    
    auto [peak_bin, peak_val] = FindMaxBinAndValueInRange(h1, xmin, xmax);
    Double_t peak_center = h1->GetXaxis()->GetBinCenter(peak_bin);

    TF1 fG("fG", "gaus", xmin, xmax);
    fG.SetParameters(peak_val, peak_center, def_sigma);  // 初期値: 振幅, mean, sigma
    fG.SetParLimits(0, peak_val*0.7, peak_val*1.5);
    fG.SetParLimits(1, xmin, xmax);
    fG.SetParLimits(2, def_sigma*0.3, def_sigma*2.);
    h1->Fit(&fG, "RQ");

    // --- gaus fit ---
    h1->Fit(&fG, "QNR");
    Double_t mean  = fG.GetParameter(1);
    Double_t sigma = fG.GetParameter(2);

    // --- 領域定義 ---
    Double_t ca_centers[num_ca];
    Double_t ca_xmin[num_ca];
    Double_t ca_xmax[num_ca];
    Double_t integ[num_ca];
    Double_t err[num_ca];

    std::vector<TLine*> lines;
    Double_t xmin = MasstoTOF(ca_mass[0]-2, peak_center);
    Double_t xmax = MasstoTOF(ca_mass[num_ca-1]+2, peak_center);
    Double_t ymin = 0.1;
    Double_t ymax = h1->GetMaximum() * 100.;
    
    // --- 積分と誤差 ---
    for (int i = 0; i < num_ca; ++i) {
      Double_t eff_mass = ca_mass[i];
      eff_mass -= (ca_mass[i] - mass_40Ca)*0.014;
      ca_centers[i] = MasstoTOF(eff_mass, peak_center, 0);// without TOF correction
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
    h1->GetYaxis()->SetRangeUser(ymin,ymax*0.5);
    // --- ヒストグラムを描画してから、線を重ねる ---
    for (auto line : lines) {
      line->SetLineColor(kRed);
      line->SetLineStyle(2);
      line->Draw("same");
    }
    canv->SetGrid();
    canv->SaveAs(Form("%s%s_Y%d.png",figdir.c_str(),base.c_str(),iy));
    // --- 格納 ---

    delete h1;
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
  mg->GetYaxis()->SetRangeUser(1,1e6);
  leg->Draw();
  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  c2->Update();
  c2->SaveAs(figpath_png.c_str());
  c2->Update();
  c2->SaveAs(figpath_pdf.c_str());
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  for (int i = 0; i < num_ca; ++i) gr[i]->Write();
  fout->Close();
  //  fin->Close();
  return;
}



