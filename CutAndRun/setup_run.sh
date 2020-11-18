#!/bin/bash


dir=/data/Seq37H/SYNC/200626_NB552201_0089_AH23TYBGXG/OUTPUT/FASTQ

mkdir -p FASTQ;

wd=`pwd`
cd FASTQ;

FILES=($(ls $dir/*gz | grep -v Undetermined))
for i in ${!FILES[@]}
do
    f=${FILES[$i]}
    echo "1: $f";
    f2=${f/_001.fastq.gz/.fastq.gz}
    outname=$(basename $f2)
    echo $outname;
    ln -s $f $outname
done

cd $wd

#python Scripts/make_config.py --prefix TEST --template Templates/template_CutAndRunConfig.json --meta pairs.tab  --exp_type CutAndRun
