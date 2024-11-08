#!/bin/bash

#REFERENCE=/n/dat/hg38/hg38.fa
KEEPVCF=0

while getopts "1:2:o:m:v:r:f" option; do
  case $option in
    1) VCF1="$OPTARG";;
	  2) VCF2="$OPTARG";;
	  o) OUTPUTNAME="$OPTARG";;
    m) MASKFILE="$OPTARG";;
    v) SVBED="$OPTARG";;
    r) REFERENCE="$OPTARG";;
    f) KEEPVCF=1;;
  esac
done

VCF1TRUNK="${VCF1##*/}"
VCF1NAME="${VCF1TRUNK%*.vcf.gz}"
VCF1_OUT="${VCF1TRUNK%*.vcf.gz}.filtered.vcf"

VCF2TRUNK="${VCF2##*/}"
VCF2_OUT="${VCF2TRUNK%*.vcf.gz}.filtered.vcf"
VCF2NAME="${VCF2TRUNK%*.vcf.gz}"

TEMPOUTPUT="$OUTPUTNAME"_tmp

echo "Making temporary directory $TEMPOUTPOUT"

mkdir -p $TEMPOUTPUT


echo "Sorting"
picard SortVcf -I $VCF1 -O $TEMPOUTPUT/$VCF1NAME.sorted.vcf
picard SortVcf -I $VCF2 -O $TEMPOUTPUT/$VCF2NAME.sorted.vcf

echo "Masking to confident regions"
bedtools intersect -a $TEMPOUTPUT/$VCF1NAME.sorted.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf
bedtools intersect -a $TEMPOUTPUT/$VCF2NAME.sorted.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf

echo "Filtering non passing variants"
bcftools view -f PASS $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.temp.vcf
bcftools view -f PASS $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.temp.vcf

bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$VCF1NAME.sorted.masked.temp.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.filtered.vcf
bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$VCF2NAME.sorted.masked.temp.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.filtered.vcf

echo "running concordance"
if [ $KEEPVCF -eq 1 ]
then
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.sorted.masked.filtered.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.sorted.masked.filtered.vcf --OUTPUT $OUTPUTNAME --OUTPUT_VCF
else
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.sorted.masked.filtered.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.sorted.masked.filtered.vcf --OUTPUT $OUTPUTNAME
fi


echo "Done"
exit 0