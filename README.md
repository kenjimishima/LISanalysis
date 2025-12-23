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
