#include "include.h"


TH2D* ReadTH2_XYMatrix_core(const char* inpath,
			    const char* hname="h2",
			    const char* delim=" ,\t",
			    bool draw=true)
{
  std::ifstream ifs(inpath);
  if(!ifs.good()){ ::Error("ReadTH2_XYMatrix","Cannot open %s", inpath); return nullptr; }

  std::string line;
  std::vector<double> xcenters;
  while(std::getline(ifs,line)){
    if(is_blank(line) || starts_with_hash(line)) continue;
    auto toks = split_any(line, delim);
    if(toks.empty()) continue;
    for(size_t i=1;i<toks.size();++i){
      double xv; if(parse_double(toks[i],xv)) xcenters.push_back(xv);
    }
    if(!xcenters.empty()) break;
  }
  if(xcenters.empty()){
    ::Error("ReadTH2_XYMatrix","No numeric X headers found."); return nullptr;
  }

  std::vector<double> ycenters;
  std::vector<std::vector<double>> Z_rows;
  while(std::getline(ifs,line)){
    if(is_blank(line) || starts_with_hash(line)) continue;
    auto toks = split_any(line, delim); if(toks.empty()) continue;
    double yv; if(!parse_double(toks[0], yv)) continue;
    ycenters.push_back(yv);

    std::vector<double> row; row.reserve(xcenters.size());
    for(size_t i=1;i<toks.size();++i){ double zv; if(parse_double(toks[i],zv)) row.push_back(zv); }
    if(row.size()<xcenters.size()) row.resize(xcenters.size(), 0.0);
    if(row.size()>xcenters.size()) row.resize(xcenters.size());
    Z_rows.push_back(std::move(row));
  }
  if(ycenters.empty()){ ::Error("ReadTH2_XYMatrix","No Y/Z data read."); return nullptr; }

  int nx = (int)xcenters.size();
  int ny = (int)ycenters.size();
  auto Xedges = makeEdgesFromCenters(xcenters);
  auto Yedges = makeEdgesFromCenters(ycenters);

  // 固定：X = TOF, Y = X position
  TH2D* h2 = new TH2D(hname, hname, ny, &Yedges[0], nx, &Xedges[0]);
  h2->SetDirectory(nullptr);

  for(int j=0; j<ny; ++j){
    for(int i=0; i<nx; ++i){
      double z = 0.0;
      if(j < (int)Z_rows.size() && i < (int)Z_rows[j].size()) z = Z_rows[j][i];
      z *= 1000.; //V to mV
      h2->SetBinContent(j+1, i+1, z);
      h2->SetBinError(j+1, i+1, 0.);
    }
  }

  h2->GetXaxis()->SetTitle("TOF [us]");
  h2->GetYaxis()->SetTitle("X position [mm]");
  h2->GetZaxis()->SetTitle("Voltage [V]");

  ::Info("ReadTH2_XYMatrix","Axis fixed: X=TOF [us], Y=X position [mm]");

  if(draw){
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas(Form("c_%s",hname), hname, 1000, 800);
    c->SetRightMargin(0.15);
    h2->Draw("COLZ");
    c->Update();
  }

  return h2;
}

// メイン関数：basenameのみ指定
TH2D* ReadTH2_XYMatrix(const char* basename,
                       const char* hname="h2",
                       const char* delim=" ,\t",
                       bool draw=true)
{
  std::string base = strip_ext_and_dir(basename);
  std::string inpath  = "./data/" + base + ".txt";
  std::string outdir  = "./root/rawdata/";
  std::string outpath = outdir + "/" + base + ".root";
  if (gSystem->AccessPathName(outdir.c_str())) gSystem->mkdir(outdir.c_str(), true);

  TH2D* h2 = ReadTH2_XYMatrix_core(inpath.c_str(), hname, delim, draw);
  if(!h2) return nullptr;

  TFile f(outpath.c_str(), "RECREATE");
  if(!f.IsZombie()){ h2->Write(); f.Close(); ::Info("ReadTH2_XYMatrix","Saved to %s", outpath.c_str()); }
  else { ::Error("ReadTH2_XYMatrix","Failed to create %s", outpath.c_str()); }

  return h2;
}
