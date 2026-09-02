#include "include.h"

void CombineRanges(
		   const char* basename_low   = "RUN72_Spatial_Beamoff_500_withlens_LD70mW_1mVscale",
		   const char* basename_high = "RUN72_Spatial_Beamoff_500_withlens_LD70mW_500mVscale",
		   const char* basename_comb = "RUN72_Spatial_Beamoff_500_withlens_LD70mW"
		   )
{
  std::string indir  = "./root/ca_graph/";
  std::string base_low = strip_ext_and_dir(basename_low);
  std::string infile_low = base_low + "_CaGraph.root";
  std::string inpath_low = indir + infile_low;
  std::string base_high = strip_ext_and_dir(basename_high);
  std::string infile_high = base_high + "_CaGraph.root";
  std::string inpath_high = indir + infile_high;
  std::string outdir  = "./root/ca_graph_combined/";
  std::string outpath = outdir +  basename_comb + "_CaGraphCombined.root";
  std::string figdir  = "./outputs/ca_graph_combined/";
  std::string figpath_pdf = figdir + basename_comb + ".pdf";
  std::string figpath_png = figdir + basename_comb + ".png";
  std::string figpath_root = figdir + basename_comb + ".root";
  if (gSystem->AccessPathName(figdir.c_str())) gSystem->mkdir(figdir.c_str(), true);
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);
  
  //--- ファイルを開く ---
  TFile* fin_low = TFile::Open(inpath_low.c_str(), "READ");
  if (!fin_low || fin_low->IsZombie()) {
    Error("CompareGraphsByName", "Cannot open file %s", inpath_low.c_str());
    return;
  }
  TFile* fin_high = TFile::Open(inpath_high.c_str(), "READ");
  if (!fin_high || fin_high->IsZombie()) {
    Error("CompareGraphsByName", "Cannot open file %s", inpath_high.c_str());
    return;
  }

  TString gr_name[num_ca];
  TGraphErrors* gr_low[num_ca];
  TGraphErrors* gr_high[num_ca];
  TGraphErrors* gr_comb[num_ca];
  //--- グラフを読み込み ---
  for (int i = 0; i < num_ca; i++) {
    gr_name[i] = Form("gr_ca%0.0f",ca_mass[i]);
    gr_low[i] = (TGraphErrors*)fin_low->Get(gr_name[i]);
    if (!gr_low[i]) {
      Error("CombinedGraph", "Cannot get %s", gr_name[i].Data());
      return;
    }
    gr_low[i]->SetName(Form("gr_low_ca%0.0f",ca_mass[i])); 

    gr_high[i] = (TGraphErrors*)fin_high->Get(gr_name[i]);
    if (!gr_high[i]) {
      Error("CombinedGraph", "Cannot get %s", gr_name[i].Data());
      return;
    }
    gr_high[i]->SetName(Form("gr_high_ca%0.0f",ca_mass[i])); 

    gr_comb[i] = (TGraphErrors*)gr_low[i]->Clone(Form("gr_comb_ca%0.0f",ca_mass[i]));
  }

  //--- 点数確認 ---
  const Int_t num = gr_comb[0]->GetN();

  // Number of Ca isotope
  for (int j = 0; j < num_ca; j++){
    // Number of x potions
    for (Int_t i = 0; i < num; i++) {
      Double_t x, y;
      gr_comb[j]->GetPoint(i, x, y);
      // Swap from low to high if the point was saturated
      if(y > PeakSaturated){
	Double_t x_high, y_high;
	gr_high[j]->GetPoint(i, x_high, y_high);
	Double_t ey = gr_high[j]->GetErrorY(i);	
	gr_comb[j]->SetPoint(i, x_high, y_high);
	gr_comb[j]->SetPointError(i, 0, ey);
	
      }
    }
  }
  
  TCanvas* c1 = new TCanvas("c1", "Ca graph combined", 800, 600);
  TLegend* leg = new TLegend(0.8, 0.8, 0.95, 0.95,"");
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(Form("Ca abundance;%s;%s",
		    gr_comb[0]->GetXaxis()->GetTitle(),gr_comb[0]->GetYaxis()->GetTitle()));  
  for (int i = 0; i < num_ca; i++) {
    mg->Add(gr_comb[i]);
    leg->AddEntry(gr_comb[i], Form("Ca%0.0f",ca_mass[i]), "pl");
  }
  mg->Draw("APL");
  gPad->SetGrid();
  gPad->SetLogy();
  mg->GetYaxis()->SetRangeUser(1.e-4,1.e4);
  leg->Draw();
  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  c1->SaveAs(figpath_png.c_str());
  c1->SaveAs(figpath_pdf.c_str());
  c1->SaveAs(figpath_root.c_str());
  TFile* fout = new TFile(outpath.c_str(), "RECREATE");
  for (int i = 0; i < num_ca; ++i) gr_comb[i]->Write();
  fout->Close();
  //  fin->Close();
}
