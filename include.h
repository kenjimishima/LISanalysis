#pragma once //read only once

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>
#include <TError.h>
#include <TKey.h>
#include <TError.h>
#include <TEnv.h>
#include <TColor.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <cmath>

gStyle->SetOptFit(1111);
gStyle->SetOptStat(1001111);

Double_t mass_40Ca = 40.;
Double_t mass_44Ca = 44.;
Double_t mass_48Ca = 48.;
const Int_t num_ca = 3;
Double_t ca_mass[num_ca] = {mass_40Ca, mass_44Ca, mass_48Ca};

Double_t baseline_xmin = 2.60;
Double_t baseline_xmax = 3.00;
Double_t def_sigma = 0.01; // [us]

//TOF parameters
Double_t trigger_delay =  0.1278; //[us]
//Double_t def_tof_40Ca = 4.09; // [us]
Double_t def_tof_40Ca = 4.124; // [us]

Int_t n_rebin = 4; //TOF rebin

namespace {
  inline bool starts_with_hash(const std::string& s){
    for(char c : s){ if(!std::isspace((unsigned char)c)) return c=='#'; }
    return false;
  }
  inline bool is_blank(const std::string& s){
    for(char c : s){ if(!std::isspace((unsigned char)c)) return false; }
    return true;
  }
  std::vector<std::string> split_any(const std::string& line, const std::string& delim){
    std::vector<std::string> out; std::string tok;
    auto is_delim = [&](char c){ return delim.find(c)!=std::string::npos; };
    for(char c : line){ if(is_delim(c)){ if(!tok.empty()){ out.push_back(tok); tok.clear(); } } else tok.push_back(c); }
    if(!tok.empty()) out.push_back(tok);
    return out;
  }
  bool parse_double(const std::string& s, double& v){
    char* endp=nullptr; v=std::strtod(s.c_str(), &endp); return endp!=s.c_str();
  }
  std::vector<double> makeEdgesFromCenters(const std::vector<double>& c){
    std::vector<double> e; size_t n=c.size(); if(n==0) return e; e.resize(n+1);
    if(n==1){ double w=1.0; e[0]=c[0]-0.5*w; e[1]=c[0]+0.5*w; return e; }
    for(size_t i=1;i<n;i++) e[i]=0.5*(c[i-1]+c[i]);
    e[0]=c[0]-0.5*(c[1]-c[0]); e[n]=c[n-1]+0.5*(c[n-1]-c[n-2]); return e;
  }
  std::string strip_ext_and_dir(const std::string& s){
    size_t p = s.find_last_of("/\\");
    std::string base = (p==std::string::npos)? s : s.substr(p+1);
    size_t dot = base.find_last_of('.');
    if(dot==std::string::npos) return base;
    return base.substr(0, dot);
  }
}

// Get TOF time from mass number
Double_t MasstoTOF(Double_t mass, Double_t tof_40Ca, const Bool_t bTOFCorr=0)
{
  Double_t tof = (tof_40Ca - trigger_delay) * TMath::Sqrt(mass/mass_40Ca) + trigger_delay;
  cout <<"Expected value for mass "<< mass <<" is "<<tof<<" [us] "<<endl;
  return tof;
}

// Get mass number from TOF time
Double_t TOFtoMass(Double_t tof, Double_t tof_40Ca, const Bool_t bTOFCorr=0)
{
  Double_t tof_nodelay = tof - trigger_delay;
  Double_t tof_40Ca_nodelay = tof_40Ca - trigger_delay;
  Double_t mass = mass_40Ca * TMath::Sq(tof_nodelay/tof_40Ca_nodelay);
  cout <<"Expected value for TOF at "<< tof <<" [us] is "<<mass<<endl;
  return mass;
}

Double_t ratio_err(Double_t a, Double_t ea, Double_t b, Double_t eb)
{
    return (a / b) * TMath::Sqrt((ea / a)*(ea / a) + (eb / b)*(eb / b));
}

Int_t FindMaxBinInRange(TH1D* h, Double_t xmin, Double_t xmax)
{
  Int_t bin_min = h->FindBin(xmin);
  Int_t bin_max = h->FindBin(xmax);

  Double_t max_val = -1;
  Int_t max_bin = bin_min;
  for (Int_t ib = bin_min; ib <= bin_max; ++ib) {
    Double_t val = h->GetBinContent(ib);
    if (val > max_val) {
      max_val = val;
      max_bin = ib;
    }
  }
  return max_bin;
}

std::pair<Int_t, Double_t> FindMaxBinAndValueInRange(TH1D* h, Double_t xmin, Double_t xmax)
{
  Int_t bin_min = h->FindBin(xmin);
  Int_t bin_max = h->FindBin(xmax);

  Double_t max_val = -1;
  Int_t max_bin = bin_min;

  for (Int_t ib = bin_min; ib <= bin_max; ++ib) {
    Double_t val = h->GetBinContent(ib);
    if (val > max_val) {
      max_val = val;
      max_bin = ib;
    }
  }
  return std::make_pair(max_bin, max_val);
}

TH2D* GetTH2Dfile(const char* basename, const char* dirname, const char* histname)
{
  std::string base = strip_ext_and_dir(basename);
  std::string inpath = std::string(dirname) + base + ".root";
  
  // --- ファイルを開く ---
  TFile* fin = new TFile(inpath.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    Error("PlotFirstSliceFromTH2D", "Cannot open file: %s", inpath.c_str());
    exit(0);
    return NULL;
  }

  // --- 指定名 or 最初のTH2Dを探す ---
  TH2D* h2 = dynamic_cast<TH2D*>(fin->Get(histname));
  if (!h2) {
    TIter nextkey(fin->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
      TObject* obj = key->ReadObj();
      if (obj->InheritsFrom(TH2::Class())) {
        h2 = (TH2D*)obj;
        std::cout << "Found TH2D: " << h2->GetName() << std::endl;
        break;
      }
    }
  }

  if (!h2) {
    Error("PlotFirstSliceFromTH2D", "No TH2D found in file.");
    fin->Close();
    return NULL;
  }

  // --- 内容確認 ---
  std::cout << "Histogram: " << h2->GetName()
            << "  Entries = " << h2->GetEntries() << std::endl;

  if (h2->GetEntries() == 0) {
    Warning("PlotFirstSliceFromTH2D", "Histogram has no entries (may appear empty).");
  }
  return h2;
}


