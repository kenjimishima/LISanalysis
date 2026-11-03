// ROOT macro for the transverse position of a laser-deflected Ca atomic beam.
// Run from the directory containing this file:
//   root -l -q 'CaLaserDeflection.C()'
//   root -l -q 'CaLaserDeflection.C(0.20,0.25,1.0,400,700,Mu48Ca,mass_48Ca,3.0)'
// All lengths are in metres unless the argument name states another unit.

#include <TGraphErrors.h>
#include <TMath.h>
#include "include.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLine.h>
#include <TPaveStats.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

void CaLaserDeflection(
    double zLaserStart_m = 0.20,
    double zLaserEnd_m = 0.25,
    double zDetector_m = 0.80,
    int nAbsorptions = 1000, // Mean count per atom at resonance, before Doppler suppression
    double meanVz_mps = 700.0,
    double wavelength_nm = Mu48Ca,
    double mass_u = mass_48Ca,
    double sigmaDivergence_mrad = 3.0,
    double sigmaSource_mm = 2.0,
    double sigmaVz_mps = 0.0,
    bool spontaneousEmission = true,
    double directedKicksPerAbsorption = 1.0,
    int nAtoms = 10000,
    unsigned int seed = 12345,
    double sourceX_mm = 50.0,
    int nStep = 250, // Number of equal laser-region intervals
    double gamma_MHz = 60.0, // Effective absorption FWHM in MHz (not angular frequency)
    bool dopplerEffect = true,
    double initialVx_mps = +5.0) // Mean initial x velocity; laser propagates along -x
{
    gStyle->SetOptStat(1001100); // Integral, mean and standard deviation (sigma).

    // z=0 is the source/reference plane. The laser kicks in the -x direction.
    // sigmaDivergence_mrad and sigmaSource_mm are 1-sigma Gaussian widths.
    // spontaneousEmission=true adds one isotropic emission recoil per absorption.
    // For ideal counterpropagating stimulated emission, set spontaneousEmission
    // to false and directedKicksPerAbsorption to 2 (one absorption + one emission).
    if (zLaserStart_m < 0 || zLaserEnd_m <= zLaserStart_m ||
        zDetector_m < zLaserEnd_m || nAbsorptions < 0 || nAtoms <= 0 || nStep <= 0 ||
        meanVz_mps <= 0 || wavelength_nm <= 0 || mass_u <= 0 ||
        sigmaDivergence_mrad < 0 || sigmaSource_mm < 0 ||
        sigmaVz_mps < 0 || directedKicksPerAbsorption < 0 ||
        !std::isfinite(gamma_MHz) || gamma_MHz <= 0 ||
        !std::isfinite(initialVx_mps)) {
        std::fprintf(stderr, "Invalid geometry or parameter value.\n");
        return;
    }

    const double h = 6.62607015e-34;      // J s
    const double atomicMassUnit = 1.66053906660e-27; // kg
    const double recoil_mps = h / (mass_u * atomicMassUnit * wavelength_nm * 1e-9);
    const double kick_mps = -directedKicksPerAbsorption * recoil_mps;
    // Laser propagates along -x and is resonant for an atom at rest.
    // v_parallel = -vx; detuning in the atom frame is -v_parallel/lambda.
    const double wavelength_m = wavelength_nm * 1e-9;
    const double gamma_Hz = gamma_MHz * 1e6;
    auto relativeAbsorptionRate = [&](double vx) {
        if (!dopplerEffect) return 1.0;
        const double detuning_Hz = vx / wavelength_m;
        const double normalizedDetuning = 2.0 * detuning_Hz / gamma_Hz;
        return 1.0 / (1.0 + normalizedDetuning * normalizedDetuning);
    };
    const double sourceX_m = sourceX_mm * 1e-3;
    const double sigmaTheta = sigmaDivergence_mrad * 1e-3;
    const double sigmaX = sigmaSource_mm * 1e-3;
    const double dz = zLaserEnd_m - zLaserStart_m;
    const int nSteps = 250;
    const double zCenter = 0.5 * (zLaserStart_m + zLaserEnd_m);
    // Conditional Poisson draws; Doppler feedback makes the total non-Poisson.
    const double meanAbsPerStep = static_cast<double>(nAbsorptions) / nStep;
    std::vector<double> kickZ(nStep);
    double sumLeverArmSq = 0.0;
    for (int i = 0; i < nStep; ++i) {
        kickZ[i] = zLaserStart_m + dz * (i + 0.5) / nStep;
        sumLeverArmSq += std::pow(zDetector_m - kickZ[i], 2);
    }

    // Mean at the detector: each kick contributes (kick/vz)*(zDetector-zKick).
    const double expectedMean_m = sourceX_m + initialVx_mps / meanVz_mps * zDetector_m + nAbsorptions * kick_mps /
                                  meanVz_mps * (zDetector_m - zCenter);
    TRandom3 rng(seed);
    std::vector<double> hit_mm;
    hit_mm.reserve(nAtoms);
    std::vector<double> hitOff_mm;
    hitOff_mm.reserve(nAtoms);
    std::vector<double> zPlot(nSteps + 1), meanPlot(nSteps + 1), meanOffPlot(nSteps + 1);
    std::vector<double> sumKickSlope(nStep, 0.0);
    double sumX0 = 0.0, sumTheta0 = 0.0;

    double sumAbs = 0.0, sumAbsSq = 0.0;
    for (int a = 0; a < nAtoms; ++a) {
        const double vz = sigmaVz_mps == 0 ? meanVz_mps :
                          rng.Gaus(meanVz_mps, sigmaVz_mps);
        if (vz <= 0) { --a; continue; } // Reject unphysical Gaussian draws.
        const double x0 = rng.Gaus(sourceX_m, sigmaX);
        const double theta0 = rng.Gaus(0, sigmaTheta);
        double vx = initialVx_mps + theta0 * vz;
        const double initialSlope = vx / vz;
        sumX0 += x0;
        sumTheta0 += initialSlope;
        double xDetector = x0 + initialSlope * zDetector_m;
        // Paired comparison with exactly the same initial conditions.
        hitOff_mm.push_back(1000.0 * xDetector);
        double atomAbsorptions = 0.0;
        for (int step = 0; step < nStep; ++step) {
            const double zi = kickZ[step];
            // Freeze the rate within this interval; update velocity for the next.
            const double meanStep = meanAbsPerStep * relativeAbsorptionRate(vx);
            const int stepAbsorptions = rng.Poisson(meanStep);
            atomAbsorptions += stepAbsorptions;
            double dvx = stepAbsorptions * kick_mps;
            if (spontaneousEmission) {
                // One isotropic spontaneous-emission recoil per absorbed photon.
                for (int i = 0; i < stepAbsorptions; ++i)
                    dvx -= recoil_mps * rng.Uniform(-1.0, 1.0);
            }
            vx += dvx; // Include both directed and spontaneous-emission recoil.
            sumKickSlope[step] += dvx / vz;
            xDetector += (dvx / vz) * (zDetector_m - zi);
        }
        sumAbs += atomAbsorptions;
        sumAbsSq += atomAbsorptions * atomAbsorptions;
        hit_mm.push_back(1000.0 * xDetector);
    }

    // Actual ensemble mean trajectory, including velocity-dependent absorption.
    for (int j = 0; j <= nSteps; ++j) {
        zPlot[j] = zDetector_m * j / nSteps;
        meanOffPlot[j] = 1000.0 * (sumX0 + sumTheta0 * zPlot[j]) / nAtoms;
        meanPlot[j] = meanOffPlot[j];
        for (int step = 0; step < nStep; ++step)
            if (zPlot[j] >= kickZ[step])
                meanPlot[j] += 1000.0 * sumKickSlope[step] / nAtoms * (zPlot[j] - kickZ[step]);
    }

    const double sampledMeanAbs = sumAbs / nAtoms;
    const double sampledSigmaAbs = std::sqrt(std::max(0.0,
        sumAbsSq / nAtoms - sampledMeanAbs * sampledMeanAbs));
    std::printf("Mean initial vx: %.6f m/s (negative = laser direction)\n", initialVx_mps);
    std::printf("Doppler effect: %s; absorption FWHM %.3f MHz (%.3f m/s); laser resonant at rest\n",
                dopplerEffect ? "on" : "off", gamma_MHz, wavelength_m * gamma_Hz);
    std::printf("Resonant reference (before Doppler): %d steps, mean per step %.6f; total mean %d, sigma %.3f\n",
                nStep, meanAbsPerStep, nAbsorptions, std::sqrt(static_cast<double>(nAbsorptions)));
    std::printf("Sampled absorption count: mean %.3f, std.dev. %.3f\n",
                sampledMeanAbs, sampledSigmaAbs);
    double sumOff = 0, sumSqOff = 0;
    for (double x : hitOff_mm) { sumOff += x; sumSqOff += x*x; }
    const double meanOff = sumOff / nAtoms;
    const double rmsOff = std::sqrt(std::max(0.0, sumSqOff / nAtoms - meanOff*meanOff));
    std::printf("Mass: %.8f u; wavelength: %.8f nm\n", mass_u, wavelength_nm);
    std::printf("Laser off: mean %.3f mm, std.dev. %.3f mm\n", meanOff, rmsOff);
    double sum = 0, sumSq = 0;
    for (double x : hit_mm) { sum += x; sumSq += x * x; }
    const double mean_mm = sum / nAtoms;
    const double rms_mm = std::sqrt(std::max(0.0, sumSq / nAtoms - mean_mm * mean_mm));
    const double sigmaEmission_mm = spontaneousEmission ?
        std::sqrt(meanAbsPerStep * sumLeverArmSq / 3.0) * recoil_mps / meanVz_mps * 1000.0 : 0.0;
    std::printf("Recoil speed per photon: %.6f m/s\n", recoil_mps);
    std::printf("Mean directed transverse speed change: %.6f m/s\n",
                sampledMeanAbs * kick_mps);
    std::printf("No-Doppler reference detector x (fixed vz): %.3f mm\n", 1000 * expectedMean_m);
    std::printf("Simulated detector x: mean %.3f mm, std.dev. %.3f mm (%d atoms)\n",
                mean_mm, rms_mm, nAtoms);
    std::printf("No-Doppler reference emission std.dev. (fixed vz): %.3f mm\n",
                sigmaEmission_mm);

    TCanvas canvas("cDeflection", "Ca laser deflection", 1200, 500);
    canvas.Divide(2, 1);
    // Match frame widths: left pad * 0.80 == right pad * 0.60.
    const double split = 0.43;
    canvas.GetPad(1)->SetPad(0.01, 0.01, split, 0.99);
    canvas.GetPad(2)->SetPad(split, 0.01, 0.99, 0.99);
    canvas.GetPad(1)->SetLeftMargin(0.10);
    canvas.GetPad(1)->SetRightMargin(0.10);
    canvas.GetPad(2)->SetLeftMargin(0.10);

    canvas.cd(1);
    TGraph trajectory(nSteps + 1, zPlot.data(), meanPlot.data());
    trajectory.SetTitle("Simulated mean trajectory;z from source (m);x (mm)");
    trajectory.SetLineColor(kBlue + 1);
    trajectory.SetLineWidth(3);
    const double trajectoryLow = std::min(*std::min_element(meanPlot.begin(), meanPlot.end()),
                                          *std::min_element(meanOffPlot.begin(), meanOffPlot.end()));
    const double trajectoryHigh = std::max(*std::max_element(meanPlot.begin(), meanPlot.end()),
                                           *std::max_element(meanOffPlot.begin(), meanOffPlot.end()));
    const double trajectoryMargin = std::max(1.0, 0.15 * (trajectoryHigh - trajectoryLow));
    trajectory.SetMinimum(trajectoryLow - trajectoryMargin);
    trajectory.SetMaximum(trajectoryHigh + trajectoryMargin);
    trajectory.Draw("AL");
    TGraph trajectoryOff(nSteps + 1, zPlot.data(), meanOffPlot.data());
    trajectoryOff.SetLineColor(kOrange + 7);
    trajectoryOff.SetLineStyle(2);
    trajectoryOff.SetLineWidth(3);
    trajectoryOff.Draw("L SAME");
    gPad->Update();
    TLine laserStart(zLaserStart_m, gPad->GetUymin(), zLaserStart_m, gPad->GetUymax());
    laserStart.SetLineStyle(2); laserStart.SetLineColor(kRed + 1);
    laserStart.Draw("SAME");
    TLine laserEnd(zLaserEnd_m, gPad->GetUymin(), zLaserEnd_m, gPad->GetUymax());
    laserEnd.SetLineStyle(2); laserEnd.SetLineColor(kRed + 1);
    laserEnd.Draw("SAME");
    TLegend legend(0.15, 0.15, 0.54, 0.32);
    legend.AddEntry(&trajectory, "Laser recoil on", "l");
    legend.AddEntry(&trajectoryOff, "Laser recoil off", "l");
    legend.AddEntry(&laserStart, "Laser region boundaries", "l");
    legend.Draw();

    canvas.cd(2);
    gPad->SetRightMargin(0.30); // Reserve space outside the frame for the legend.
    // Use common bins covering both samples, including their tails.
    const double low = std::min(*std::min_element(hit_mm.begin(), hit_mm.end()),
                                *std::min_element(hitOff_mm.begin(), hitOff_mm.end()));
    const double high = std::max(*std::max_element(hit_mm.begin(), hit_mm.end()),
                                 *std::max_element(hitOff_mm.begin(), hitOff_mm.end()));
    const double margin = std::max(1.0, 0.05 * (high - low));
    TH1D histogram("hDetector", "Detector distribution;x at detector (mm);Atoms",
                   100, low - margin, high + margin);
    TH1D histogramOff("hDetectorOff", "Detector distribution;x at detector (mm);Atoms",
                      100, low - margin, high + margin);
    histogram.SetStats(kTRUE);
    histogramOff.SetStats(kTRUE);
    histogram.SetDirectory(nullptr);
    histogramOff.SetDirectory(nullptr);
    histogram.SetLineColor(kBlue + 1);
    histogram.SetLineWidth(2);
    histogramOff.SetLineColor(kOrange + 7);
    histogramOff.SetLineStyle(2);
    histogramOff.SetLineWidth(2);
    for (double x : hit_mm) histogram.Fill(x);
    for (double x : hitOff_mm) histogramOff.Fill(x);
    histogram.SetMaximum(1.45 * std::max(histogram.GetMaximum(), histogramOff.GetMaximum()));
    histogram.Draw("HIST");
    histogramOff.Draw("HIST SAMES"); // SAMES also creates the second stats box.
    gPad->Update();
    auto positionStats = [](TH1D& hist, double x1, double x2) {
        auto* stats = dynamic_cast<TPaveStats*>(hist.FindObject("stats"));
        if (!stats) return;
        stats->SetX1NDC(x1);
        stats->SetX2NDC(x2);
        stats->SetY1NDC(0.75);
        stats->SetY2NDC(0.89);
        stats->SetTextColor(hist.GetLineColor());
        stats->SetLineColor(hist.GetLineColor());
        stats->SetFillColor(kWhite);
        stats->SetTextSize(0.035);
    };
    positionStats(histogram, 0.12, 0.38); // Laser on: upper left.
    positionStats(histogramOff, 0.43, 0.69); // Laser off: upper right.
    TLegend detectorLegend(0.72, 0.75, 0.99, 0.89);
    detectorLegend.SetBorderSize(0);
    detectorLegend.SetFillStyle(0);
    detectorLegend.SetTextSize(0.028);
    detectorLegend.AddEntry(&histogram, "Laser recoil on", "l");
    detectorLegend.AddEntry(&histogramOff, "Laser recoil off", "l");
    detectorLegend.Draw();
    gPad->Modified();
    gPad->Update();
    canvas.SaveAs("CaLaserDeflection.pdf");
    canvas.SaveAs("CaLaserDeflection.png");
    std::printf("Saved CaLaserDeflection.pdf and CaLaserDeflection.png\n");
}
