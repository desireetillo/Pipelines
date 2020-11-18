#!/bin/bash


dir=/data/Seq37H/SYNC/200818_NB552201_0115_AHCWVYBGXG/OUTPUT/FASTQ

mkdir -p FASTQ;

wd=`pwd`
cd FASTQ;
i=1
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

# python Scripts/make_config.py --prefix prefix --template Templates/template_RRBS.json --meta meta.tab --exp_type RRBS
