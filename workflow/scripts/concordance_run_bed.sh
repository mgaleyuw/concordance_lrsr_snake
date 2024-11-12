#!/bin/bash

#REFERENCE=/n/dat/hg38/hg38.fa
KEEPVCF=0

while getopts "1:2:o:m:b:v:r:f" option; do
  case $option in
      1) VCF1="$OPTARG";;
	  2) VCF2="$OPTARG";;
	  o) OUTPUTNAME="$OPTARG";;
      m) MASKFILE="$OPTARG";;
      b) BED="$OPTARG";;
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

echo "Removing INFO columns and filtering to canonical autosomes"
bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 -O v -o $TEMPOUTPUT/$VCF1NAME.vcf $VCF1
bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 -O v -o $TEMPOUTPUT/$VCF2NAME.vcf $VCF2
bcftools annotate -x "INFO" -o $TEMPOUTPUT/$VCF1NAME.clean.vcf $TEMPOUTPUT/$VCF1NAME.vcf
bcftools annotate -x "INFO" -o $TEMPOUTPUT/$VCF2NAME.clean.vcf $TEMPOUTPUT/$VCF2NAME.vcf

echo "reheadering"
bcftools view -h $TEMPOUTPUT/$VCF1NAME.clean.vcf | grep -v "random" | grep -v "chrUn" > $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf
bcftools view -H $TEMPOUTPUT/$VCF1NAME.clean.vcf >> $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf
bcftools view -h $TEMPOUTPUT/$VCF2NAME.clean.vcf | grep -v "random" | grep -v "chrUn" > $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf
bcftools view -H $TEMPOUTPUT/$VCF2NAME.clean.vcf >> $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf


echo "Sorting"
picard SortVcf -I $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf -O $TEMPOUTPUT/$VCF1NAME.sorted.vcf
picard SortVcf -I $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf -O $TEMPOUTPUT/$VCF2NAME.sorted.vcf

echo "Masking to confident regions"
bedtools intersect -a $TEMPOUTPUT/$VCF1NAME.sorted.vcf -b $BED -wa -header > $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf
bedtools intersect -a $TEMPOUTPUT/$VCF2NAME.sorted.vcf -b $BED -wa -header > $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf

echo "filter to bed file"
bedtools intersect -a $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.vcf
bedtools intersect -a $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.vcf

echo "Filtering non passing variants"
bcftools view -f PASS $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.temp.vcf
bcftools view -f PASS $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.temp.vcf


if grep -q "##FORMAT=<ID=DP" <(bcftools view -h $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.temp.vcf)
then
  echo "Filtering $VCF1NAME to reads with Depth >=5"
  bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.temp.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.filtered.vcf
else
  echo "no DP tag in $TEMPOUTPUT/$VCF1NAME, skipping depth filter"
  mv $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.temp.vcf $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.filtered.vcf
fi

if grep -q "##FORMAT=<ID=DP" <(bcftools view -h $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.temp.vcf)
then
  echo "Filtering $VCF2NAME to reads with Depth >=5"
  bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.temp.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.filtered.vcf
else
  echo "no DP tag in $TEMPOUTPUT/$VCF1NAME, skipping depth filter"
  mv $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.temp.vcf $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.filtered.vcf
fi

echo "running concordance"
if [ $KEEPVCF -eq 1 ]
then
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.filtered.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.filtered.vcf --OUTPUT $OUTPUTNAME --OUTPUT_VCF
else
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.filtered.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.filtered.vcf --OUTPUT $OUTPUTNAME
fi


echo "Done"
exit 0