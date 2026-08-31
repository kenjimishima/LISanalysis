#include "include.h"

void RatioGraph(const char* basename = "RUN72_Spatial_Beamoff_500_withlens_LD70mW")
{
  std::string base = strip_ext_and_dir(basename);
  std::string infile = base + "_CaGraphCombined.root";
  std::string indir  = "./root/ca_graph_combined/";
  std::string inpath = indir + infile;
  std::string outdir  = "./root/ratio/";
  std::string outpath = outdir +  base + "_Ratio.root";
  std::string figdir  = "./outputs/ratio/";
  std::string figpath_log_pdf = figdir + base + "_log.pdf";
  std::string figpath_log_png = figdir + base + "_log.png";
  std::string figpath_lin_pdf = figdir + base + "_lin.pdf";
  std::string figpath_lin_png = figdir + base + "_lin.png";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);
  
  //--- ファイルを開く ---
  TFile* fin = TFile::Open(inpath.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    Error("CompareGraphsByName", "Cannot open file %s", inpath.c_str());
    return;
  }


  TString gr_name[num_ca];
  TGraphErrors* gr[num_ca];
  TGraphErrors* gr_ratio[num_ca];
  //--- グラフを読み込み ---
  for (int i = 0; i < num_ca; i++) {
    gr_name[i] = Form("gr_comb_ca%0.0f",ca_mass[i]);
    gr[i] = (TGraphErrors*)fin->Get(gr_name[i]);
    if (!gr[i]) {
      Error("RatioGraph", "Cannot get %s", gr_name[i].Data());
      return;
    }
    gr_ratio[i] = (TGraphErrors*)gr[i]->Clone(Form("%s_ratio", gr_name[i].Data()));
  }
  /*
  TCanvas* c = new TCanvas("c", "Graph Ratios", 800, 600);
  gr[0]->Draw("APL");
  gr[1]->Draw("PLsame");
  gr[2]->Draw("PLsame");
  gr[3]->Draw("PLsame");
  */
  const Int_t num = gr[0]->GetN();
  std::vector<Double_t> sumY(num, 0.0);
  std::vector<Double_t> sumEsq(num, 0.0);
  //--- 点数確認 ---
  for (Int_t i = 0; i < num; i++) {
    sumY[i] = 0.0;
    sumEsq[i] = 0.0;
    for (int j = 0; j < num_ca; j++){
      //--- 総和と誤差計算 ---
      Double_t x, y;
      gr[j]->GetPoint(i, x, y);
      Double_t ey = gr[j]->GetErrorY(i);
      sumY[i] += y;
      sumEsq[i] += ey * ey;
    }
    for (int j = 0; j < num_ca; j++){
      //--- 総和と誤差計算 ---
      Double_t x, y;
      gr[j]->GetPoint(i, x, y);
      Double_t ey = gr[j]->GetErrorY(i);
      gr_ratio[j]->SetPoint(i, x, y/sumY[i]);
      Double_t rerr = ratio_err(y, ey, sumY[i], TMath::Sqrt(sumEsq[i]));
      gr_ratio[j]->SetPointError(i, 0., rerr);
    }
  }
  
  TCanvas* c1 = new TCanvas("c1", "Graph Ratios", 800, 600);
  TLegend* leg = new TLegend(0.8, 0.8, 0.95, 0.95,"");
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(Form("Ca abundance;%s;Abundance",
		    gr_ratio[0]->GetXaxis()->GetTitle()));  
  for (int i = 0; i < num_ca; i++) {
    mg->Add(gr_ratio[i]);
    leg->AddEntry(gr_ratio[i], gr_ratio[i]->GetTitle(), "pl");
  }
  mg->Draw("APL");
  gPad->SetGrid();
  gPad->SetLogy();
  mg->GetYaxis()->SetRangeUser(1.e-3,10);
  leg->Draw();
  c1->SaveAs(figpath_log_png.c_str());
  c1->SaveAs(figpath_log_pdf.c_str());
  gPad->SetLogy(0);
  mg->GetYaxis()->SetRangeUser(-0.2, 1.2);
  gPad->Update();
  c1->SaveAs(figpath_lin_png.c_str());
  c1->SaveAs(figpath_lin_pdf.c_str());
  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  for (int i = 0; i < num_ca; ++i) gr_ratio[i]->Write();
  fout->Close();
  //  fin->Close();
}
