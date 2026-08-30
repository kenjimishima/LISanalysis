#!/bin/bash
# List file names

#Input working dir name 
#runname=RUN45
#runname=RUN51
#runname=RUN52
#runname=RUN53
#runname=RUN54
#runname=RUN61
#runname=RUN62
#runname=RUN69
runname=RUN72

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
# RUN54_Spatial_Beamoff_500.txt
# RUN54_Spatial_Beamoff_550.txt
# RUN54_Spatial_Beamoff_600.txt
# RUN54_Spatial_Beamon48_500_withlens.txt
# RUN54_Spatial_Beamon48_500_wolens.txt
# RUN54_Spatial_Beamon48_550_withlens.txt
# RUN54_Spatial_Beamon48_550_wolens.txt
#RUN61_Spatial_Beamoff_550_withlens_10mVscale_500mVscaleatcenter.txt
#RUN61_Spatial_Beamon48_550_withlens_10mVscale_500mVscaleatcenter.txt
#RUN62_Spatial_Beamoff_550_withlens_10mVscale.txt
#RUN62_Spatial_Beamon48_550_withlens_10mVscale.txt
#RUN62_Spatial_Beamoff_600_withlens_10mVscale.txt
#RUN62_Spatial_Beamon48_600_withlens_10mVscale.txt
#RUN62_Spatial_Beamoff_700_withlens_10mVscale.txt
#RUN62_Spatial_Beamon48_700_withlens_10mVscale.txt
#RUN69_450_P30_Beamoff40Ca.txt
#RUN69_450_P30_Beamoff44Ca.txt
#RUN69_450_P30_Beamon40Ca.txt
#RUN69_450_P30_Beamon44Ca.txt
#RUN69_450_P30_Beamoff42Ca.txt
#RUN69_450_P30_Beamoff48Ca.txt
#RUN69_450_P30_Beamon42Ca.txt
#RUN69_450_P30_Beamon48Ca.txt
RUN72_Spatial_Beamoff_500_withlens_LD70mW_1mVscale.txt
RUN72_Spatial_Beamon48_500_withlens_LD70mW_1mVscale.txt
RUN72_Spatial_Beamoff_500_withlens_LD70mW_500mVscale.txt
RUN72_Spatial_Beamon48_500_withlens_LD70mW_500mVscale.txt
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
    root -l -b -q "../RatioGraph.C(\"$f\")"
done

f0="${files[0]}"
f1="${files[1]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"

f0="${files[2]}"
f1="${files[3]}"
root -l -b -q "../LaserEffect.C(\"$f0\",\"$f1\")"



: <<'EOF'
EOF

