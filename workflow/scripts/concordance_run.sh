#!/bin/bash

set -eo pipefail

#REFERENCE=/n/dat/hg38/hg38.fa
KEEPVCF=0
AUTOSOMES=0
FILTERPASS=0
DEBIANLIKE=0

while getopts "1:2:o:m:b:v:r:fapw:" option; do
  case $option in
      1) VCF1="$OPTARG";;
	    2) VCF2="$OPTARG";;
	    o) OUTPUTNAME="$OPTARG";;
      m) MASKFILE="$OPTARG";;
      b) BED="$OPTARG";;
      v) SVBED="$OPTARG";;
      r) REFERENCE="$OPTARG";;
      f) KEEPVCF=1;;
      a) AUTOSOMES=1;;
      p) FILTERPASS=1;;
      w) WORKDIR="$OPTARG" ;;
  esac
done

OSTEXT=$( cat /etc/os-release )

if [[ $OSTEXT == *"debian"* ]]
then
  DEBIANLIKE=1
fi

VCF1TRUNK="${VCF1##*/}"
VCF1NAME="${VCF1TRUNK%*.vcf.gz}"
VCF1_OUT="${VCF1TRUNK%*.vcf.gz}.filtered.vcf"

VCF2TRUNK="${VCF2##*/}"
VCF2_OUT="${VCF2TRUNK%*.vcf.gz}.filtered.vcf"
VCF2NAME="${VCF2TRUNK%*.vcf.gz}"

if [ -z ${WORKDIR+x} ]
then
  TEMPOUTPUT="$OUTPUTNAME"_tmp
else
  TEMPOUTPUT="$WORKDIR"
fi

echo "Making temporary directory $TEMPOUTPUT"
mkdir -p $TEMPOUTPUT

if [ $AUTOSOMES -eq 1 ]
then
  chroms=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
else
  chroms=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
fi

echo "Removing INFO columns and filtering to canonical chromosomes"
bcftools view -r $chroms -O v -o $TEMPOUTPUT/$VCF1NAME.vcf $VCF1
bcftools view -r $chroms -O v -o $TEMPOUTPUT/$VCF2NAME.vcf $VCF2
bcftools annotate -x "INFO" -o $TEMPOUTPUT/$VCF1NAME.clean.vcf $TEMPOUTPUT/$VCF1NAME.vcf
bcftools annotate -x "INFO" -o $TEMPOUTPUT/$VCF2NAME.clean.vcf $TEMPOUTPUT/$VCF2NAME.vcf

echo "reheadering, removing structural variants"
bcftools view -h $TEMPOUTPUT/$VCF1NAME.clean.vcf | grep -v "random" | grep -v "chrUn" | grep -v "chrM" > $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf
bcftools view -H $TEMPOUTPUT/$VCF1NAME.clean.vcf | grep -v "HGSV" >> $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf
bcftools view -h $TEMPOUTPUT/$VCF2NAME.clean.vcf | grep -v "random" | grep -v "chrUn" | grep -v "chrM"  > $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf
bcftools view -H $TEMPOUTPUT/$VCF2NAME.clean.vcf | grep -v "HGSV" >> $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf


echo "Sorting"
picard SortVcf -I $TEMPOUTPUT/$VCF1NAME.clean.reheader.vcf -O $TEMPOUTPUT/$VCF1NAME.sorted.vcf
picard UpdateVcfSequenceDictionary -I $TEMPOUTPUT/$VCF2NAME.clean.reheader.vcf -O $TEMPOUTPUT/$VCF2NAME.clean.reheader.updated.vcf -SD $TEMPOUTPUT/$VCF1NAME.sorted.vcf
picard SortVcf -I $TEMPOUTPUT/$VCF2NAME.clean.reheader.updated.vcf -O $TEMPOUTPUT/$VCF2NAME.sorted.vcf

echo "Masking to confident regions"
bedtools intersect -a $TEMPOUTPUT/$VCF1NAME.sorted.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf
bedtools intersect -a $TEMPOUTPUT/$VCF2NAME.sorted.vcf -b $MASKFILE -wa -header > $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf

if [ -z ${BED+x} ]
then
  if [ $FILTERPASS -eq 1 ]
  then
    echo "Filtering non passing variants"
    bcftools view -f PASS $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.temp.vcf
    bcftools view -f PASS $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.temp.vcf

    PRECONCORDANCEFILE1=$VCF1NAME.sorted.masked.temp.vcf
    PRECONCORDANCEFILE2=$VCF2NAME.sorted.masked.temp.vcf
  else
    PRECONCORDANCEFILE1=$VCF1NAME.sorted.masked.vcf
    PRECONCORDANCEFILE2=$VCF2NAME.sorted.masked.vcf
  fi
else
  echo "filter to bed file"
  bedtools intersect -a $TEMPOUTPUT/$VCF1NAME.sorted.masked.vcf -b $BED -wa -header > $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.vcf
  bedtools intersect -a $TEMPOUTPUT/$VCF2NAME.sorted.masked.vcf -b $BED -wa -header > $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.vcf

  if [ $FILTERPASS -eq 1 ]
  then
    echo "Filtering non passing variants"
    bcftools view -f PASS $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.vcf > $TEMPOUTPUT/$VCF1NAME.sorted.masked.bed.temp.vcf
    bcftools view -f PASS $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.vcf > $TEMPOUTPUT/$VCF2NAME.sorted.masked.bed.temp.vcf

    PRECONCORDANCEFILE1=$VCF1NAME.sorted.masked.bed.temp.vcf
    PRECONCORDANCEFILE2=$VCF2NAME.sorted.masked.bed.temp.vcf
  else
    PRECONCORDANCEFILE1=$VCF1NAME.sorted.masked.bed.vcf
    PRECONCORDANCEFILE2=$VCF2NAME.sorted.masked.bed.vcf
  fi
fi

if grep -q "##FORMAT=<ID=DP" <(bcftools view -h $TEMPOUTPUT/$PRECONCORDANCEFILE1)
then
  echo "Filtering $VCF1NAME to reads with Depth >=5"
  bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$PRECONCORDANCEFILE1 > $TEMPOUTPUT/$VCF1NAME.preconcordance.vcf
else
  echo "no DP tag in $TEMPOUTPUT/$VCF1NAME, skipping depth filter"
  mv $TEMPOUTPUT/$PRECONCORDANCEFILE1 $TEMPOUTPUT/$VCF1NAME.preconcordance.vcf
fi

if grep -q "##FORMAT=<ID=DP" <(bcftools view -h $TEMPOUTPUT/$PRECONCORDANCEFILE2)
then
  echo "Filtering $VCF2NAME to reads with Depth >=5"
  bcftools view -e "FORMAT/DP<5" $TEMPOUTPUT/$PRECONCORDANCEFILE2 > $TEMPOUTPUT/$VCF2NAME.preconcordance.vcf
else
  echo "no DP tag in $TEMPOUTPUT/$VCF1NAME, skipping depth filter"
  mv $TEMPOUTPUT/$PRECONCORDANCEFILE2 $TEMPOUTPUT/$VCF2NAME.preconcordance.vcf
fi

echo "running concordance"
if [ $KEEPVCF -eq 1 ]
then
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.preconcordance.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.preconcordance.vcf --OUTPUT $OUTPUTNAME --OUTPUT_VCF
else
  picard GenotypeConcordance --CALL_VCF $TEMPOUTPUT/$VCF2NAME.preconcordance.vcf --TRUTH_VCF $TEMPOUTPUT/$VCF1NAME.preconcordance.vcf --OUTPUT $OUTPUTNAME
fi

if [ $DEBIANLIKE -eq 0 ]
then
  rename concordance_contingency_metrics concordance.contingency_metrics $OUTPUTNAME*
  rename concordance_detail_metrics concordance.detail_metrics $OUTPUTNAME*
  rename concordance_summary_metrics concordance.summary_metrics $OUTPUTNAME*
else
  rename s/concordance_contingency_metrics/concordance.contingency_metrics/ $OUTPUTNAME*
  rename s/concordance_detail_metrics/concordance.detail_metrics/ $OUTPUTNAME*
  rename s/concordance_summary_metrics/concordance.summary_metrics/ $OUTPUTNAME*
fi


echo "Done"
exit 0