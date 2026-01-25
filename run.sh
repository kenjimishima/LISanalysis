#!/bin/bash
# List file names

#Input working dir name 
#runname=RUN45
#runname=RUN51
#runname=RUN52
#runname=RUN53
runname=RUN54

#Input files
files=(
# RUN45_Spatial_40Ca_Beamoff.txt
# RUN45_Spatial_40Ca_Beamon48.txt
# RUN45_Spatial_44Ca_48Ca_Beamoff.txt
# RUN45_Spatial_44Ca_48Ca_Beamon48.txt
# RUN51_Spatial_40Ca_Beamoff.txt
# RUN51_Spatial_40Ca_Beamon48.txt
# RUN52_Spatial_Beamoff_550.txt
# RUN52_Spatial_Beamoff_650.txt
# RUN52_Spatial_Beamon48_550.txt
# RUN52_Spatial_Beamon48_650.txt
# RUN53_Spatial_Beamoff_550.txt
# RUN53_Spatial_Beamoff_600.txt
# RUN53_Spatial_Beamon48_550.txt
# RUN53_Spatial_Beamon48_600.txt
RUN54_Spatial_Beamoff_500.txt
RUN54_Spatial_Beamoff_550.txt
RUN54_Spatial_Beamoff_600.txt
RUN54_Spatial_Beamon48_500_withlens.txt
RUN54_Spatial_Beamon48_500_wolens.txt
RUN54_Spatial_Beamon48_550_withlens.txt
RUN54_Spatial_Beamon48_550_wolens.txt
)

mkdir -p $runname/data
for f in "${files[@]}"; do
  cp "./data/$f" "$runname/data/"
done
cd $runname

: <<'EOF'

#Create ROOT files
for f in "${files[@]}"; do
    echo -e "\nStart Processing $f ..."
    #Read text
    root -l -b -q "../ReadTH2_XYMatrix.C(\"$f\")"
    #Baseline correction
    root -l -b -q "../BaselineCorrection.C(\"$f\")"
done

#Normalization to counts and errors to be chi2/ndf of unity, using Yslice of file[0]
Yslice=1
f0="${files[0]}"
echo -e "\nCalibration slice $Yslice of $f ..."
root -l -b -q "../ErrorCalibration.C(\"$f0\",$Yslice)"

#Scale TOF 
for f in "${files[@]}"; do
    echo -e "\nStart Scaling $f ..."
    root -l -b -q "../ScaleTH2D.C(\"$f\",\"$f0\",$Yslice)"
done

#Draw TOF histogram and fit 40Ca peak
for f in "${files[@]}"; do
    echo -e "\nStart making TOF histogram $f ..."
    for ((Yslice=1; Yslice<=17; Yslice++)); do
	root -l -b -q "../PeakFit.C(\"$f\",$Yslice)"
    done
done

#Identity peaks and integrate their counts
for f in "${files[@]}"; do
    f0="${files[0]}"
    for ((Yslice=1; Yslice<=17; Yslice++)); do
	root -l -b -q "../PeakFinder.C(\"$f0\",$Yslice)"
    done
done
    
#Get Ca peak ratios
for f in "${files[@]}"; do
    echo -e "\nStart Making Graphs $f ..."
    root -l -b -q "../GetCaGraphs.C(\"$f\")"
    root -l -b -q "../RatioGraph.C(\"$f\")"
done

EOF

f0="${files[0]}"
f1="${files[3]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"

f0="${files[0]}"
f1="${files[4]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"

f0="${files[1]}"
f1="${files[5]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"

f0="${files[1]}"
f1="${files[6]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"

