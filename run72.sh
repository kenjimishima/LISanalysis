#!/bin/bash
# List file names

#Input working dir name 
runname=RUN72

#Input files
files=(
RUN72_Spatial_Beamoff_500_withlens_LD70mW_1mVscale.txt
RUN72_Spatial_Beamoff_500_withlens_LD70mW_500mVscale.txt
RUN72_Spatial_Beamon48_500_withlens_LD70mW_1mVscale.txt
RUN72_Spatial_Beamon48_500_withlens_LD70mW_500mVscale.txt
)

#Combined file name
files_comb=(
RUN72_Spatial_Beamoff_500_withlens_LD70mW
RUN72_Spatial_Beamon48_500_withlens_LD70mW
)

mkdir -p $runname/data
for f in "${files[@]}"; do
  cp "./data/$f" "$runname/data/"
done
cd $runname

f="${files[0]}"
nrow=$(wc -l < "./data/$f")
ncol=$(awk 'NR==1 { print NF; exit }' "./data/$f")
echo "rows=$nrow"
echo "cols=$ncol"

#Create ROOT files
for f in "${files[@]}"; do
    echo -e "\nStart Processing $f ..."
    #Read text
    root -l -b -q "../ReadTH2_XYMatrix.C(\"$f\")"
    #Baseline correction
    root -l -b -q "../BaselineCorrection.C(\"$f\")"
done

#Normalization to counts and errors to be chi2/ndf of unity, using Yslice of each run
Yslice=1
#f0="${files[0]}"
for f in "${files[@]}"; do
echo -e "\nCalibration slice $Yslice of $f ..."
root -l -b -q "../ErrorCalibration.C(\"$f\",$Yslice)"
done

#Scale TOF 
for f in "${files[@]}"; do
    echo -e "\nStart Scaling $f ..."
    root -l -b -q "../ScaleTH2D.C(\"$f\",\"$f\",$Yslice)"
done

#Draw TOF histogram and fit 40Ca peak
for f in "${files[@]}"; do
    echo -e "\nStart making TOF histogram $f ..."
    for ((Yslice=1; Yslice<ncol; Yslice++)); do
	root -l -b -q "../PeakFit.C(\"$f\",$Yslice)"
    done
done

: <<'EOF'
#Identity peaks and integrate their counts
for f in "${files[@]}"; do
    f0="${files[0]}"
    for ((Yslice=1; Yslice<=ncol; Yslice++)); do
	root -l -b -q "../PeakFinder.C(\"$f0\",$Yslice)"
    done
done
EOF

#Get Ca peak ratios
for f in "${files[@]}"; do
    echo -e "\nStart Making Graphs $f ..."
    root -l -b -q "../GetCaGraphs.C(\"$f\")"
done

f0="${files[0]}"
f1="${files[1]}"
fc="${files_comb[0]}"
root -l -b -q "../CombineRanges.C(\"$f0\",\"$f1\",\"$fc\")"

f0="${files[2]}"
f1="${files[3]}"
fc="${files_comb[1]}"
root -l -b -q "../CombineRanges.C(\"$f0\",\"$f1\",\"$fc\")"

for f in "${files_comb[@]}"; do
    echo -e "\nStart Making Ratio $f ..."
    root -l -b -q "../RatioGraph.C(\"$f\")"
done

f0="${files_comb[0]}"
f1="${files_comb[1]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"


: <<'EOF'
EOF

