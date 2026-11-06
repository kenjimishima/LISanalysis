#include "include.h"
#include "CaThermalDisplacementModel.h"
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TPaveStats.h>
#include <algorithm>
#include <cmath>

// Examples:
//   root -l -q 'LaserEffect.C("beamoff","beamon",40)'
//   root -l -q 'LaserEffect.C("beamoff","beamon",42)'
//   root -l -q 'LaserEffect.C("beamoff","beamon",44)'
//   root -l -q 'LaserEffect.C("beamoff","beamon",48)'
// The selected isotope uses the thermal model only for its with-laser graph.

void LaserEffect(const char* basename1 = "RUN75_Spatial_Beamoff_500_Xs0",
		 const char* basename2 = "RUN75_Spatial_Beamon48_500_Xs0",
         int withLaserThermalFitMass = 0) // 0: Gaussian; 40, 42, 44, or 48: select thermal model
{
  if (withLaserThermalFitMass != 0 && withLaserThermalFitMass != 40 &&
      withLaserThermalFitMass != 42 && withLaserThermalFitMass != 44 &&
      withLaserThermalFitMass != 48) {
    Error("LaserEffect", "Fit option must be 0, 40, 42, 44, or 48.");
    return;
  }
  std::string indir  = "./root/ca_graph_combined/";
  std::string base1 = strip_ext_and_dir(basename1);
  std::string infile1 = base1 + "_CaGraphCombined.root";
  std::string base2 = strip_ext_and_dir(basename2);
  std::string infile2 = base2 + "_CaGraphCombined.root";
  std::string inpath1 = indir + infile1;
  std::string inpath2 = indir + infile2;
  std::string outdir  = "./root/laser_effect/";
  std::string fitSuffix = withLaserThermalFitMass == 0 ? "" :
      Form("_ThermalCa%d", withLaserThermalFitMass);
  std::string outpath = outdir + base1 + "_" + base2 + fitSuffix + "_LaserEffect.root";
  std::string figdir  = "./outputs/laser_effect/";
  std::string figpath_pdf = figdir + base1 + "_" + base2 + fitSuffix + "_LaserEffect.pdf";
  std::string figpath_png = figdir + base1 + "_" + base2 + fitSuffix + "_LaserEffect.png";
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
    TString gr_name = Form("gr_comb_ca%0.0f",ca_mass[i]);
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
    TString gr_name = Form("gr_comb_ca%0.0f",ca_mass[i]);
    gr2[i] = (TGraphErrors*)fin2->Get(gr_name);
    if (!gr2[i]) {
      Error("RatioGraph", "Cannot get %s in %s", gr_name.Data(), inpath2.c_str());
      return;
    }
  }
  
  TCanvas* c1 = new TCanvas("c1", "Laser Effect", 1200, 800);
  gr1[0]->Draw("goff");
  double gr_max = gr1[0]->GetYaxis()->GetXmax();
  cout <<"Maximum value: "<< gr_max <<" [pC]" << endl;
  c1->Divide(2,2);
  TLegend* leg = new TLegend(0.63, 0.5, 0.88, 0.7,"");
  leg->AddEntry(gr1[0], "Without Laser", "pl");
  leg->AddEntry(gr2[0], "With    Laser", "pl");

  TF1* funcgaus1[num_ca];
  TF1* funcgaus2[num_ca];
  TF1* funcThermal[num_ca] = {};
  Double_t peak1[num_ca];
  Double_t peak2[num_ca];
  Double_t mean1[num_ca];
  Double_t mean2[num_ca];
  Double_t sigma1[num_ca];
  Double_t sigma2[num_ca];

  Double_t peak1_err[num_ca];
  Double_t peak2_err[num_ca];
  Double_t mean1_err[num_ca];
  Double_t mean2_err[num_ca];
  Double_t sigma1_err[num_ca];
  Double_t sigma2_err[num_ca];

  TMultiGraph* mg[num_ca];
  for (int i = 0; i < num_ca; i++) {
    funcgaus1[i] = new TF1(Form("func_%s",gr1[i]->GetName()),"gaus+[3]");
    funcgaus2[i] = new TF1(Form("func_%s",gr2[i]->GetName()),"gaus+[3]");
    funcgaus1[i]->SetNpx(1000);
    funcgaus2[i]->SetNpx(1000);
    funcgaus1[i]->SetParNames("A","Mean","Sigma","Baseline");
    funcgaus2[i]->SetParNames("A","Mean","Sigma","Baseline");
    funcgaus1[i]->SetParameters(gr_max*ca_abundance[i], 50., 10., 0.);
    funcgaus2[i]->SetParameters(gr_max*ca_abundance[i], 50., 10., 0.);
    funcgaus1[i]->SetParLimits(0,0.,1.e7);
    funcgaus1[i]->SetParLimits(1,30,60);
    funcgaus1[i]->SetParLimits(2,1.,20.);
    funcgaus2[i]->SetParLimits(0,0.,1.e7);
    funcgaus2[i]->SetParLimits(1,30,60);
    funcgaus2[i]->SetParLimits(2,1.,20.);

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
    double ymin = mg[i]->GetYaxis()->GetXmin();
    double ymax = mg[i]->GetYaxis()->GetXmax();
    mg[i]->GetYaxis()->SetRangeUser(ymin,ymax*1.5);
    gr1[i]->Fit(funcgaus1[i]);
    peak1[i] = funcgaus1[i]->GetParameter(0);
    mean1[i] = funcgaus1[i]->GetParameter(1);
    sigma1[i] = funcgaus1[i]->GetParameter(2);
    peak1_err[i] = funcgaus1[i]->GetParError(0);
    mean1_err[i] = funcgaus1[i]->GetParError(1);
    sigma1_err[i] = funcgaus1[i]->GetParError(2);

    const int isotopeMass = static_cast<int>(std::round(ca_mass[i]));
    if (isotopeMass == withLaserThermalFitMass) {
      const double* xValues = gr2[i]->GetX();
      const auto xBounds = std::minmax_element(xValues, xValues + gr2[i]->GetN());
      const double fitXMin = *xBounds.first;
      const double fitXMax = *xBounds.second;
      const double initialMaximum = std::max(1.0, TMath::MaxElement(gr2[i]->GetN(), gr2[i]->GetY()));
      const double initialArea = initialMaximum * std::max(0.2, sigma1[i]) *
                                 std::sqrt(2.0 * std::acos(-1.0));
      constexpr bool fixTemperature =
          CaThermalDisplacement::FitConfig::fixTemperatureAtInitialValue;
      constexpr bool fitTemperature = !fixTemperature;
      constexpr double initialTemperature =
          CaThermalDisplacement::FitConfig::initialTemperature_K;
      constexpr double initialRate =
          CaThermalDisplacement::FitConfig::initialRate_photonsPerSecond;
      constexpr double initialPeakCenter =
          CaThermalDisplacement::FitConfig::initialPeakCenter_mm;
      constexpr bool fixPeakCenter =
          CaThermalDisplacement::FitConfig::fixPeakCenterAtInitialValue;
      const double fixedSigma = sigma1[i];
      const double isotopeMass_u = ca_mass[i];
      const double isotopeWavelength_nm = ca_wavelength[i];
      auto thermalModel = [=](double* x, double* p) {
        const double temperature = fitTemperature ? p[4] : initialTemperature;
        const double modelParameters[5] = {p[0], temperature, p[1],
                                            fixedSigma, p[3]};
        return CaThermalDisplacement::Spectrum(x[0], modelParameters,
            isotopeMass_u, isotopeWavelength_nm) + p[2];
      };
      funcThermal[i] = new TF1(Form("func_thermal_Ca%d", isotopeMass),
          thermalModel, fitXMin, fitXMax, 5);
      funcThermal[i]->SetParNames("Normalization", "N [photons/s]", "Baseline",
                                  "Peak center [mm]", "Temperature [K]");
      funcThermal[i]->SetParameters(initialArea, initialRate,
                                    0.0, initialPeakCenter, initialTemperature);
      funcThermal[i]->SetParLimits(0, 1e-9, std::max(1e3, initialArea * 1e4));
      const double rateMin = std::max(1.0, initialRate / 1000.0);
      const double rateMax = initialRate * 1000.0;
      funcThermal[i]->SetParLimits(1, rateMin, rateMax);
      funcThermal[i]->SetParLimits(2, -initialMaximum, initialMaximum);
      if (fixPeakCenter) {
        funcThermal[i]->FixParameter(3, initialPeakCenter);
      } else {
        funcThermal[i]->SetParLimits(3, fitXMin, fitXMax);
      }
      if (fitTemperature) {
        funcThermal[i]->SetParLimits(4,
            CaThermalDisplacement::FitConfig::minimumTemperature_K,
            CaThermalDisplacement::FitConfig::maximumTemperature_K);
      } else {
        funcThermal[i]->FixParameter(4, initialTemperature);
      }
      funcThermal[i]->SetLineColor(kRed);
      funcThermal[i]->SetLineWidth(2);
      funcThermal[i]->SetNpx(1000);
      TFitResultPtr thermalResult = gr2[i]->Fit(funcThermal[i], "RS");
      if (static_cast<int>(thermalResult) != 0)
        Warning("LaserEffect", "Thermal fit for Ca%d returned status %d", isotopeMass,
                static_cast<int>(thermalResult));

      peak2[i] = funcThermal[i]->GetParameter(0);
      peak2_err[i] = funcThermal[i]->GetParError(0);
      mean2[i] = funcThermal[i]->GetMaximumX(fitXMin, fitXMax);
      sigma2[i] = sigma1[i]; // The laser-off beam width is fixed in this model.
      sigma2_err[i] = 0.0;
      const double fittedRate = funcThermal[i]->GetParameter(1);
      const double rateError = funcThermal[i]->GetParError(1);
      const double fittedCenter = funcThermal[i]->GetParameter(3);
      const double centerError = funcThermal[i]->GetParError(3);
      const double fittedTemperature = funcThermal[i]->GetParameter(4);
      const double temperatureError = funcThermal[i]->GetParError(4);
      // Propagate the peak-position uncertainty with the fit covariance,
      // including the strong N-T correlation in this model.
      const int propagatedParameters[] = {1, 3, 4};
      double peakGradient[3] = {0.0, 0.0, 0.0};
      for (int j = 0; j < 3; ++j) {
        const int parameter = propagatedParameters[j];
        if ((parameter == 3 && fixPeakCenter) ||
            (parameter == 4 && !fitTemperature)) continue;
        const double value = funcThermal[i]->GetParameter(parameter);
        const double error = funcThermal[i]->GetParError(parameter);
        double low = value - error;
        double high = value + error;
        if (parameter == 1) { low = std::max(rateMin, low); high = std::min(rateMax, high); }
        if (parameter == 3) { low = std::max(fitXMin, low); high = std::min(fitXMax, high); }
        if (parameter == 4) {
          low = std::max(CaThermalDisplacement::FitConfig::minimumTemperature_K, low);
          high = std::min(CaThermalDisplacement::FitConfig::maximumTemperature_K, high);
        }
        if (error <= 0.0 || high <= low) continue;
        funcThermal[i]->SetParameter(parameter, low);
        const double peakLow = funcThermal[i]->GetMaximumX(fitXMin, fitXMax);
        funcThermal[i]->SetParameter(parameter, high);
        const double peakHigh = funcThermal[i]->GetMaximumX(fitXMin, fitXMax);
        peakGradient[j] = (peakHigh - peakLow) / (high - low);
        funcThermal[i]->SetParameter(parameter, value);
      }
      double peakVariance = 0.0;
      bool hasCovariance = thermalResult.Get() != nullptr;
      if (hasCovariance) {
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col)
            peakVariance += peakGradient[row] *
                thermalResult->CovMatrix(propagatedParameters[row], propagatedParameters[col]) *
                peakGradient[col];
      }
      if (hasCovariance && peakVariance >= 0.0) {
        mean2_err[i] = std::sqrt(peakVariance);
      } else {
        double diagonalVariance = 0.0;
        for (int j = 0; j < 3; ++j) {
          const double error = funcThermal[i]->GetParError(propagatedParameters[j]);
          diagonalVariance += std::pow(peakGradient[j] * error, 2);
        }
        mean2_err[i] = std::sqrt(diagonalVariance);
      }
      cout << "Ca" << isotopeMass << " with-laser thermal fit: T="
           << fittedTemperature << " +/- " << temperatureError
           << (fitTemperature ? " K, N=" : " K fixed, N=") << fittedRate << " +/- "
           << rateError << " photons/s; fitted center="
           << fittedCenter << " +/- " << centerError
           << " mm; spectrum peak x=" << mean2[i] << " +/- "
           << mean2_err[i] << " mm" << endl;
    } else {
      gr2[i]->Fit(funcgaus2[i]);
      peak2[i] = funcgaus2[i]->GetParameter(0);
      mean2[i] = funcgaus2[i]->GetParameter(1);
      sigma2[i] = funcgaus2[i]->GetParameter(2);
      peak2_err[i] = funcgaus2[i]->GetParError(0);
      mean2_err[i] = funcgaus2[i]->GetParError(1);
      sigma2_err[i] = funcgaus2[i]->GetParError(2);
    }
    gPad->SetGrid();
    gPad->Update();
    TPaveStats *stats1 = (TPaveStats*)gr1[i]->GetListOfFunctions()->FindObject("stats");
    TPaveStats *stats2 = (TPaveStats*)gr2[i]->GetListOfFunctions()->FindObject("stats");
    if (stats1) {
      stats1->SetTextColor(1);
      stats1->SetX1NDC(0.12); stats1->SetX2NDC(0.40); stats1->SetY1NDC(0.75);stats1->SetY2NDC(0.99);
    }
    if (stats2) {
      stats2->SetTextColor(2);
      stats2->SetX1NDC(0.72); stats2->SetX2NDC(0.99); stats2->SetY1NDC(0.75);stats2->SetY2NDC(0.99);
    }
    
    leg->Draw();
  }

  cout << "Peak shift : " << endl;
  for (int i = 0; i < num_ca; i++){
    Double_t error = TMath::Sqrt(mean1_err[i]*mean1_err[i] + mean2_err[i]*mean2_err[i]);
    cout << "Ca" << std::round(ca_mass[i]) <<": "<< mean1[i]  - mean2[i] << " +/- "<<error <<" [mm]"<<endl;
  }

  //--------------------------------------------------------------
  // 保存
  //--------------------------------------------------------------
  c1->SaveAs(figpath_png.c_str());
  c1->SaveAs(figpath_pdf.c_str());
  c1->SaveAs(outpath.c_str());
  //  fin->Close();
}
