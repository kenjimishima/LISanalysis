#ifndef CA_THERMAL_DISPLACEMENT_MODEL_H
#define CA_THERMAL_DISPLACEMENT_MODEL_H

#include <cmath>

namespace CaThermalDisplacement {
namespace FitConfig {
// Initial fit values. The fix* options hold the corresponding parameter at its initial value.
// N and T are strongly correlated when both are free.
constexpr double initialTemperature_K = 500.0 + 273.15;
constexpr double initialRate_photonsPerSecond = 1.4e7;
constexpr double initialPeakCenter_mm = 47.46;
constexpr bool fixTemperatureAtInitialValue = true;
constexpr bool fixPeakCenterAtInitialValue = false;
constexpr double minimumTemperature_K = 50.0;
constexpr double maximumTemperature_K = 3000.0;
} // namespace FitConfig

constexpr double kB = 1.380649e-23;             // J/K
constexpr double h = 6.62607015e-34;            // J s
constexpr double amu = 1.66053906660e-27;       // kg

// Fixed geometry shared by the plotting and laser-effect fit macros.
constexpr double laserWidth_d_m = 0.030;        // Ca path length through laser
constexpr double laserCenterZ_m = 0.0;
constexpr double detectorZ_m = 0.680;
constexpr double driftDistance_L_m = detectorZ_m - laserCenterZ_m;

// F(x) is the expected signal per mm. Parameters:
// p[0] = normalization, p[1] = T [K], p[2] = N [photons/s],
// p[3] = Gaussian sigma [mm], p[4] = undeflected center x0 [mm].
// The physical calculation is intentionally separated into steps:
// 1. vz follows the positive longitudinal flux distribution;
// 2. laser transit time is d/vz;
// 3. mean absorbed photons are N*d/vz;
// 4. drift time to detector is L/vz;
// 5. mean recoil displacement is -N*d*vr*L/vz^2;
// 6. the result is smeared by a Gaussian position response.
inline double Spectrum(double x_mm, const double* p, double mass_u,
                       double wavelength_nm) {
    const double count = p[0], temperature = p[1], rate = p[2];
    const double sigma = p[3], center = p[4];
    if (count <= 0 || temperature <= 0 || rate <= 0 || sigma <= 0 ||
        mass_u <= 0 || wavelength_nm <= 0) return 0.0;

    const double mass = mass_u * amu;
    const double thermalEnergy = kB * temperature;
    const double recoilSpeed = h / (mass * wavelength_nm * 1e-9);
    const double mostProbableScale = std::sqrt(2.0 * thermalEnergy / mass);

    // With u=vz/sqrt(2*kB*T/m), f(vz;T)dvz=2*u^3*exp(-u^2)du.
    constexpr int n = 600; // even, for Simpson integration
    constexpr double uMin = 0.0, uMax = 8.0;
    const double step = (uMax - uMin) / n;
    const double gaussianNorm = 1.0 / (std::sqrt(2.0 * std::acos(-1.0)) * sigma);
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double u = uMin + i * step;

        // Step 1: longitudinal speed and its probability density.
        const double vz = mostProbableScale * u;
        const double speedDensity = 2.0 * u * u * u * std::exp(-u * u);
        double value = 0.0;
        if (vz > 0.0) {
            // Steps 2-3: laser crossing time and mean absorbed photons.
            const double laserTransitTime_s = laserWidth_d_m / vz;
            const double meanAbsorbedPhotons = rate * laserTransitTime_s;

            // Step 4: flight time from the laser center to the detector.
            const double detectorFlightTime_s = driftDistance_L_m / vz;

            // Step 5: mean transverse recoil displacement, converted to mm.
            const double displacement_mm = -1000.0 * meanAbsorbedPhotons *
                                           recoilSpeed * detectorFlightTime_s;
            const double meanDetectorX_mm = center + displacement_mm;

            // Step 6: Gaussian detector-position response for this speed.
            const double residual = (x_mm - meanDetectorX_mm) / sigma;
            const double gaussian = gaussianNorm * std::exp(-0.5 * residual * residual);
            value = speedDensity * gaussian;
        }
        sum += (i == 0 || i == n) ? value : (i % 2 ? 4.0 : 2.0) * value;
    }
    return count * step * sum / 3.0;
}
} // namespace CaThermalDisplacement

#endif
