#!/bin/bash

APPEND=0

while getopts "o:i:1:2:" opt; do
  case "$opt" in
    o) OUTPUTFILE=$OPTARG ;;
    i) INPUT=$OPTARG ;;
    1) VCF1NAME=$OPTARG ;;
    2) VCF2NAME=$OPTARG ;;
  esac
done

if [ -z ${OUTPUTFILE+x} ]
then
    trunk=${INPUT##*/}
    filename=${trunk%*.vcf.gz}
else
    filename=$OUTPUTFILE
fi

if [ -z ${VCF1NAME+x} ]
then
  VCF1NAME="VCF1"
fi

if [ -z ${VCF2NAME+x} ]
then
  VCF2NAME="VCF2"
fi


#call
bcftools view -i 'INFO/CONC_ST="FP"' $INPUT | bcftools view -e 'INFO/CONC_ST="FN" || INFO/CONC_ST="TP"' > $filename.Unique.$VCF2NAME.vcf
#truth
bcftools view -i 'INFO/CONC_ST="FN"' $INPUT | bcftools view -e 'INFO/CONC_ST="TP" || INFO/CONC_ST="FP"' > $filename.Unique.$VCF1NAME.vcf
#overlap
bcftools view -i 'INFO/CONC_ST="TP"' $INPUT | bcftools view -e 'INFO/CONC_ST="FP" || INFO/CONC_ST="FN"' > $filename.Overlap.vcf
#multiallelic
bcftools view -i 'INFO/CONC_ST="FP" && INFO/CONC_ST="FN"' $INPUT > $filename.Multiallelic.Mismatch.vcf
#mixedZygosity
bcftools view -i 'INFO/CONC_ST="TP"' $INPUT | bcftools view -i 'INFO/CONC_ST="FP" | INFO/CONC_ST="FN"' > $filename.Mixed.Zygosity.Mismatch.vcf

