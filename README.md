# LISanalysis

This repository includes code for analyzing laser isotope separation.

> ./run.sh

makes all analysized results and files

----
How to use the run.sh

Global valiables and parameters are defined in 
include.h
Useful functions are also written in here.

Input files are written at the begining of run.sh
files=(
RUN45_Spatial_40Ca_Beamoff.txt
RUN45_Spatial_40Ca_Beamon48.txt
RUN45_Spatial_44Ca_48Ca_Beamoff.txt
RUN45_Spatial_44Ca_48Ca_Beamon48.txt
)

ReadTH2_XYMatrix.C(\"$f\")
From the .txt file to ROOT files
./data/xxx.txt ./root/rawdata/xxx.root

BaselineCorrection.C(\"$f\")
apply baseline correction 
include.h


Error calibration

Subtract the average value of the defined region from each Y (laser position) to obtain TH2D.
Double_t baseline_xmin = 2.50;
Double_t baseline_xmax = 3.00;
Yslice=1             
f0="${files[0]}"                                                                                                                   
root -l -b -q "ErrorCalibration.C(\"$f0\",$Yslice)"

Assume 1 mV = 1 count and apply error = sqrt(count)
Fit the outgas peak at M = 24-29.                                                                                              
Calculate chi2/ndf and determine count_per_mV so that chi2/ndf to be unity. In usual, 1 mV corresponds to about 10 counts/bin.
Save the env file in ./results/.

Multiply by count_per_mV and add sqrt(N) error.
Save to ./root/scaled.
root -l -b -q "ScaleTH2D.C(\"$f\",\"$f0\",$Yslice)"



Obtain the peak integrals for Ca isotopes(40, 44, and 48)
root -l -b -q "GetCaGraphs.C(\"$f\")"
Plot the ratio
root -l -b -q "RatioGraph.C(\"$f\")"

Plot the effect of the laser irradiation: file[0] and file[1]
f0="${files[0]}"
f1="${files[1]}"
root -l -b -q "LaserEffect.C(\"$f0\",\"$f1\")"

Find peaks and fit them. Get the relation of TOFs and masses, and mass species of outgas. 
root -l -b -q "PeakFinder.C(\"$f0\",$Yslice)"

Tips
In bash script, 
# is used for commenting out 

Also 
: <<'EOF'
EOF
can comment out the region between.

## Isotope-specific laser-effect fit

`LaserEffect.C` can apply the thermal displacement model to the with-laser graph of one selected calcium isotope. The other isotope graphs are fitted with the Gaussian model. The isotope selector accepts `40`, `42`, `44`, or `48`; use `0` to use Gaussian fits for every isotope.

Run the macro from the run directory (for example, `RUN75`) after `run75.sh` has created the combined graph ROOT files. Pass the beam-off and beam-on basenames without directories, file extensions, or the `_CaGraphCombined.root` suffix:

```bash
cd RUN75
beamoff=RUN75_Spatial_Beamoff_500_Xs0
beamon=RUN75_Spatial_Beamon48_500_Xs0
for isotope in 40 42 44 48; do
    root -l -b -q "../LaserEffect.C(\"${beamoff}\",\"${beamon}\",${isotope})"
done
```

The isotope mass is the third argument. Fit starting values are configured in `CaThermalDisplacementModel.h`: `initialTemperature_K`, `initialRate_photonsPerSecond`, and `initialPeakCenter_mm`. Temperature is fixed at its configured initial value when `fixTemperatureAtInitialValue` is `true`; set it to `false` to fit temperature. The peak center starts at `initialPeakCenter_mm` and is fixed there when `fixPeakCenterAtInitialValue` is `true`; set it to `false` to let the fit vary it. Results are written under `root/laser_effect/` and plots under `outputs/laser_effect/` in the run directory.
