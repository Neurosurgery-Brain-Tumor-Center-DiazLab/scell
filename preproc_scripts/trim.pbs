#!/bin/bash

#PBS -m ae
#PBS -l nodes=1:ppn=1


export PATH=/diazlab/shared/bin:$PATH

if [ -z "${STEP+1}" ]
then
    a=($(echo $to_proc|sed 's/:/ /g'))
    cd "$indir${a[PBS_ARRAYID]}"
else
    cd $indir
    a=(Sample_*)
    idx=$((PBS_ARRAYID+STEP))
    cd "${a[idx]}"
fi
declare lf=(*R1*.gz)
declare rf=(*R2*.gz)
if [ -d tmp ]; then
    rm -fR tmp
fi
mkdir tmp

for ((i=0;i<${#lf[@]};i++))
do
    trim_galore -q 20 --nextera --length 20 -o ./tmp --paired ${lf[$i]} ${rf[$i]}
done


if [ -d trimmed ]; then
    rm -fR trimmed
fi
mkdir trimmed

lf=(./tmp/*val_1*)
rf=(./tmp/*val_2*)
for ((i=0;i<${#lf[@]};i++))
do
    trim_galore -q 20 --fastqc --length 20 -o ./trimmed --paired ${lf[$i]} ${rf[$i]}
done

cat ./tmp/*trimming_report.txt ./trimmed/*trimming_report.txt >trimming_report.txt
rm -f ./trimmed/*trimming_report.txt
mv trimming_report.txt ./trimmed/
rm -fR tmp

cd $indir

if [ -d trimmed ]; then
    rm -fR trimmed
fi