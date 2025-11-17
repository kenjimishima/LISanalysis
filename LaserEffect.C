#include "include.h"

void LaserEffect(const char* basename1 = "RUN45_Spatial_40Ca_Beamoff",
		 const char* basename2 = "RUN45_Spatial_40Ca_Beamon48" )
{
  std::string base1 = strip_ext_and_dir(basename1);
  std::string infile1 = base1 + "_CaGraph.root";
  std::string base2 = strip_ext_and_dir(basename2);
  std::string infile2 = base2 + "_CaGraph.root";
  std::string indir  = "./root/ca_graph/";
  std::string inpath1 = indir + infile1;
  std::string inpath2 = indir + infile2;
  std::string outdir  = "./root/laser_effect/";
  std::string outpath = outdir +  base1 + " " + base2 + "_LaserEffect.root";
  std::string figdir  = "./outputs/laser_effect/";
  std::string figpath_pdf = figdir + base1 + " " + base2 + "_LaserEffect.pdf";
  std::string figpath_png = figdir + base1 + " " + base2 + "_LaserEffect.png";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);

  //--- ファイルを開く ---
  TFile* fin1 = TFile::Open(inpath1.c_str(), "READ");
  if (!fin1 || fin1->IsZombie()) {
    Error("CompareGraphsByName", "Cannot open file %s", inpath1.c_str());
    return;
  }
  TGraphErrors* gr1[num_ca];
  //--- グラフを読み込み ---
  for (int i = 0; i < num_ca; i++) {
    TString gr_name = Form("gr_ca%0.0f",ca_mass[i]);
    gr1[i] = (TGraphErrors*)fin1->Get(gr_name);
    if (!gr1[i]) {
      Error("RatioGraph", "Cannot get %s in %s", gr_name.Data(), inpath1.c_str());
      return;
    }
  }
  
  //--- ファイルを開く ---
  TFile* fin2 = TFile::Open(inpath2.c_str(), "READ");
  if (!fin2 || fin2->IsZombie()) {
    Error("CompareGraphsByName", "Cannot open file %s", inpath2.c_str());
    return;
  }
  TGraphErrors* gr2[num_ca];
  //--- グラフを読み込み ---
  for (int i = 0; i < num_ca; i++) {
    TString gr_name = Form("gr_ca%0.0f",ca_mass[i]);
    gr2[i] = (TGraphErrors*)fin2->Get(gr_name);
    if (!gr2[i]) {
      Error("RatioGraph", "Cannot get %s in %s", gr_name.Data(), inpath2.c_str());
      return;
    }
  }
  
  TCanvas* c1 = new TCanvas("c1", "Laser Effect", 800, 900);
  c1->Divide(1,3);
  TLegend* leg = new TLegend(0.15, 0.5, 0.35, 0.7,"");
  leg->AddEntry(gr1[0], "Without Laser", "pl");
  leg->AddEntry(gr2[0], "With    Laser", "pl");

  TF1* funcgaus1[num_ca];
  TF1* funcgaus2[num_ca];
  
  TMultiGraph* mg[num_ca];
  for (int i = 0; i < num_ca; i++) {
    funcgaus1[i] = new TF1(Form("func_%s",gr1[i]->GetName()),"gaus");
    funcgaus2[i] = new TF1(Form("func_%s",gr2[i]->GetName()),"gaus");
    gr1[i]->SetLineColor(kBlack);
    gr1[i]->SetMarkerColor(kBlack);
    gr2[i]->SetLineColor(kRed);
    gr2[i]->SetMarkerColor(kRed);
    funcgaus1[i]->SetLineColor(kBlack);
    funcgaus2[i]->SetLineColor(kRed);
    
    mg[i] = new TMultiGraph();
    mg[i]->SetTitle(Form("Ca%0.0f",ca_mass[i]));
    mg[i]->Add(gr1[i]);
    mg[i]->Add(gr2[i]);
    c1->cd(i+1);
    mg[i]->Draw("AP");
    gr1[i]->Fit(funcgaus1[i]);
    gr2[i]->Fit(funcgaus2[i]);

    gPad->SetGrid();
    gPad->Update();
    TPaveStats *stats1 = (TPaveStats*)gr1[i]->GetListOfFunctions()->FindObject("stats");
    TPaveStats *stats2 = (TPaveStats*)gr2[i]->GetListOfFunctions()->FindObject("stats");
    stats1->SetTextColor(1);
    stats2->SetTextColor(2);
    stats1->SetX1NDC(0.10); stats1->SetX2NDC(0.40); stats1->SetY1NDC(0.75);stats1->SetY2NDC(0.99);
    stats2->SetX1NDC(0.70); stats2->SetX2NDC(0.99); stats2->SetY1NDC(0.75);stats2->SetY2NDC(0.99);
    
    leg->Draw();
  }
  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  c1->SaveAs(figpath_png.c_str());
  c1->SaveAs(figpath_pdf.c_str());
  c1->SaveAs(outpath.c_str());
  //  fin->Close();
}
