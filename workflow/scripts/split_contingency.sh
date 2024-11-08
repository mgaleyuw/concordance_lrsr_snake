#!/bin/bash

APPEND=0

while getopts "o:i:" opt; do
  case "$opt" in
    o) OUTPUTFILE=$OPTARG ;;
    i) INPUT=$OPTARG ;;
  esac
done

if [ -z ${OUTPUTFILE+x} ]
then
    trunk=${INPUT##*/}
    filename=${trunk%*.vcf.gz}
else
    filename=$OUTPUTFILE
fi

#call
bcftools view -i 'INFO/CONC_ST="FP"' $INPUT | bcftools view -e 'INFO/CONC_ST="FN" || INFO/CONC_ST="TP"' > $filename.Unique.Call.vcf
#truth
bcftools view -i 'INFO/CONC_ST="FN"' $INPUT | bcftools view -e 'INFO/CONC_ST="TP" || INFO/CONC_ST="FP"' > $filename.Unique.Reference.vcf
#overlap
bcftools view -i 'INFO/CONC_ST="TP"' $INPUT | bcftools view -e 'INFO/CONC_ST="FP" || INFO/CONC_ST="FN"' > $filename.Overlap.Call.Reference.vcf
#multiallelic
bcftools view -i 'INFO/CONC_ST="FP" && INFO/CONC_ST="FN"' $INPUT > $filename.Multiallelic.Mismatch.vcf
#mixedZygosity
bcftools view -i 'INFO/CONC_ST="TP"' $INPUT | bcftools view -i 'INFO/CONC_ST="FP" | INFO/CONC_ST="FN"' > $filename.Mixed.Zygosity.Mismatch.vcf

