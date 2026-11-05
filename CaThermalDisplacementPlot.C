// Thermal-beam laser displacement model for later ROOT fitting.
// Run from LISanalysis:
//   root -l -q 'CaThermalDisplacementPlot.C()'
// The fit parameters are count, temperature, absorption rate, beam sigma,
// and undeflected center. Laser width d and drift distance L are fixed below.

#include "include.h"
#include "CaThermalDisplacementModel.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TNamed.h>
#include <cmath>

// TF1 callback for the default Ca-48 demonstration plot.
inline double RootFunction(double* x, double* p) {
    return CaThermalDisplacement::Spectrum(x[0], p, mass_48Ca, Mu48Ca);
}

void CaThermalDisplacementPlot(
    double count_atoms = 10000.0,
    double temperature_K = 550.0,
    double absorptionRate_photonsPerSecond = 1.4e7,
    double beamSigma_mm = 3.124,
    double centerX_mm = 50.0,
    double xMin_mm = 0.0,
    double xMax_mm = 100.0)
{
    using namespace CaThermalDisplacement;
    if (count_atoms <= 0 || temperature_K <= 0 ||
        absorptionRate_photonsPerSecond <= 0 || beamSigma_mm <= 0 ||
        xMax_mm <= xMin_mm) {
        std::fprintf(stderr, "Invalid plot parameter.\n");
        return;
    }

    TF1 function("fCaThermalDisplacement", RootFunction, xMin_mm, xMax_mm, 5);
    function.SetParNames("Count [atoms]", "Temperature [K]",
                         "N [photons/s]", "Beam sigma [mm]", "Center x_{0} [mm]");
    function.SetParameters(count_atoms, temperature_K,
                           absorptionRate_photonsPerSecond, beamSigma_mm, centerX_mm);
    function.SetParLimits(0, 1.0, 1e9);
    function.SetParLimits(1, 1.0, 3000.0);
    function.SetParLimits(2, 1.0, 1e10);
    function.SetParLimits(3, 0.01, 100.0);
    function.SetNpx(1200);
    function.SetLineColor(kBlue + 1);
    function.SetLineWidth(3);

    constexpr double binWidth_mm = 1.0;
    const int nBins = static_cast<int>(std::ceil((xMax_mm - xMin_mm) / binWidth_mm));
    TH1D prediction("hCaThermalDisplacement",
        "Thermal-beam displacement model;Detector x (mm);Expected atoms / 1 mm",
        nBins, xMin_mm, xMin_mm + nBins * binWidth_mm);
    prediction.SetDirectory(nullptr);
    prediction.SetStats(false);
    prediction.SetLineColor(kBlue + 1);
    prediction.SetLineWidth(2);
    for (int bin = 1; bin <= nBins; ++bin) {
        const double lo = prediction.GetXaxis()->GetBinLowEdge(bin);
        const double hi = prediction.GetXaxis()->GetBinUpEdge(bin);
        prediction.SetBinContent(bin, function.Integral(lo, hi));
    }

    gStyle->SetOptStat(0);
    TCanvas canvas("cCaThermalDisplacement", "Thermal-beam displacement", 1000, 650);
    prediction.Draw("HIST");
    function.Draw("SAME");
    TLegend legend(0.53, 0.68, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&prediction, "1 mm-bin prediction", "l");
    legend.AddEntry(&function, "Fit-ready model F(x)", "l");
    legend.AddEntry((TObject*)nullptr, Form("T = %.0f K; N = %.2g photons/s", temperature_K,
                                            absorptionRate_photonsPerSecond), "");
    legend.AddEntry((TObject*)nullptr, Form("d = %.3f m; L = %.3f m", laserWidth_d_m,
                                            driftDistance_L_m), "");
    legend.AddEntry((TObject*)nullptr, Form("#sigma = %.2f mm; x_{0} = %.1f mm",
                                            beamSigma_mm, centerX_mm), "");
    legend.Draw();
    gPad->SetTicks(1, 1);
    gPad->Modified();
    gPad->Update();

    const char* directory = "CaDisplacement";
    if (gSystem->AccessPathName(directory) && gSystem->mkdir(directory, kTRUE) != 0) {
        std::fprintf(stderr, "Cannot create %s\n", directory);
        return;
    }
    const char* base = "CaDisplacement/CaThermalDisplacement";
    canvas.SaveAs(Form("%s.pdf", base));
    canvas.SaveAs(Form("%s.png", base));
    TFile output(Form("%s.root", base), "RECREATE");
    if (output.IsZombie()) {
        std::fprintf(stderr, "Cannot create ROOT output.\n");
        return;
    }
    canvas.Write();
    function.Write();
    prediction.Write();
    TParameter<double>("laserWidth_d_m", laserWidth_d_m).Write();
    TParameter<double>("driftDistance_L_m", driftDistance_L_m).Write();
    TParameter<double>("temperature_K", temperature_K).Write();
    TParameter<double>("absorptionRate_photonsPerSecond", absorptionRate_photonsPerSecond).Write();
    TNamed definition("modelDefinition",
        "F(x)=A*integral_0^infty f(vz;T)*G(x-(x0-N*d*vr*L/vz^2);sigma) dvz; "
        "f(vz;T)=m^2*vz^3/(2*(kB*T)^2)*exp[-m*vz^2/(2*kB*T)], "
        "vr=h/(m*lambda); d and L are fixed in CaThermalDisplacementPlot.C");
    definition.Write();
    output.Close();

    std::printf("Saved %s.{pdf,png,root}\n", base);
    std::printf("Fixed geometry: d=%.4f m, L=%.4f m; fit parameters: count, T, N, sigma, x0\n",
                laserWidth_d_m, driftDistance_L_m);
}
