#!/bin/bash

APPEND=0
PRINTOUTPUT=1

while getopts "o:a" opt; do
  case "$opt" in
    o) OUTPUTFILE=$OPTARG ;;
    a) APPEND=1 ;;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

inputfiles=( "$@" )

header="File,VariantType,Zygosity,Classification,Value"

if [ -z ${OUTPUTFILE+x} ]
then
    OUTPUTFILE=temp.csv
else
    PRINTOUTPUT=0
fi

if [ $APPEND -eq 0 ]
then
    echo $header > $OUTPUTFILE
fi
 
for file in ${inputfiles[@]}
do
    SNPhetTPmv=0
    SNPhetTPmv=$(grep -E "SNP.+HET.+HET.+TP$" $file | cut -f6 | paste -sd+ | bc )
    SNPhetTP=$(grep -E "SNP.+HET.+HET.+TP,TN$" $file | cut -f6 | paste -sd+ | bc)
    SNPhetTPall=$(echo "0+$SNPhetTPmv+$SNPhetTP+0" | tr -s '+' | bc)
    SNPhomTP=$(grep -E "SNP.+HOM.+HOM.+TP$" $file | cut -f6 | paste -sd+ | bc)  
    SNPhetFN=$(grep -E "SNP.+HET.+MISSING.+FN$" $file | cut -f6 | paste -sd+ | bc)
    SNPhomFN=$(grep -E "SNP.+HOM.+MISSING.+FN$" $file | cut -f6 | paste -sd+ | bc)
    SNPhetFP=$(grep -E "SNP.+MISSING.+HET.+FP$" $file | cut -f6 | paste -sd+ | bc)
    SNPhomFP=$(grep -E "SNP.+MISSING.+HOM.+FP$" $file | cut -f6 | paste -sd+ | bc)
    SNPhomMulti=$(grep -E "SNP.+HOM.+HOM.+FP,FN$" $file | cut -f6 | paste -sd+ | bc)
    SNPallMulti=$(grep -E "SNP.+FP,FN$" $file | cut -f6 | paste -sd+ | bc)
    SNPhetMulti=$(echo "$SNPallMulti-$SNPhomMulti-0" | tr -s '-' |  bc)
    SNPzygo=$(grep -E "SNP.+TP,FN$" $file | cut -f6 | paste -sd+ | bc)

    INDELhetTPmv=$(grep -E "INDEL.+HET.+HET.+TP$" $file | cut -f6 | paste -sd+ | bc)
    INDELhetTP=$(grep -E "INDEL.+HET.+HET.+TP,TN$" $file | cut -f6 | paste -sd+ | bc)
    INDELhetTPall=$(echo "0+$INDELhetTPmv+$INDELhetTP+0" | tr -s '+' | bc)
    INDELhomTP=$(grep -E "INDEL.+HOM.+HOM.+TP$" $file | cut -f6 | paste -sd+ | bc)
    INDELhetFN=$(grep -E "INDEL.+HET.+MISSING.+FN$" $file | cut -f6 | paste -sd+ | bc)
    INDELhomFN=$(grep -E "INDEL.+HOM.+MISSING.+FN$" $file | cut -f6 | paste -sd+ | bc)
    INDELhetFP=$(grep -E "INDEL.+MISSING.+HET.+FP$" $file | cut -f6 | paste -sd+ | bc)
    INDELhomFP=$(grep -E "INDEL.+MISSING.+HOM.+FP$" $file | cut -f6 | paste -sd+ | bc)
    INDELhomMulti=$(grep -E "INDEL.+HOM.+HOM.+FP,FN$" $file | cut -f6 | paste -sd+ | bc)
    INDELallMulti=$(grep -E "INDEL.+FP,FN$" $file | cut -f6 | paste -sd+ | bc)
    INDELhetMulti=$(echo "$INDELallMulti-$INDELhomMulti-0" | tr -s '-' | bc)
    INDELzygo=$(grep -E "INDEL.+TP,FN$" $file | cut -f6 | paste -sd+ | bc)

    echo "$file,SNP,hom,Overlap,$SNPhomTP" >> $OUTPUTFILE
    echo "$file,SNP,het,Overlap,$SNPhetTPall" >> $OUTPUTFILE
    echo "$file,SNP,hom,UniqueCall,$SNPhomFP" >> $OUTPUTFILE
    echo "$file,SNP,het,UniqueCall,$SNPhetFP" >> $OUTPUTFILE
    echo "$file,SNP,hom,UniqueRef,$SNPhomFN" >> $OUTPUTFILE
    echo "$file,SNP,het,UniqueRef,$SNPhetFN" >> $OUTPUTFILE
    echo "$file,SNP,hom,Multiallelic,$SNPhomMulti" >> $OUTPUTFILE
    echo "$file,SNP,het,Multiallelic,$SNPhetMulti" >> $OUTPUTFILE
    echo "$file,SNP,het,MixedZygosity,$SNPzygo" >> $OUTPUTFILE
    echo "$file,INDEL,hom,Overlap,$INDELhomTP" >> $OUTPUTFILE
    echo "$file,INDEL,het,Overlap,$INDELhetTPall" >> $OUTPUTFILE
    echo "$file,INDEL,hom,UniqueCall,$INDELhomFP" >> $OUTPUTFILE
    echo "$file,INDEL,het,UniqueCall,$INDELhetFP" >> $OUTPUTFILE
    echo "$file,INDEL,hom,UniqueRef,$INDELhomFN" >> $OUTPUTFILE
    echo "$file,INDEL,het,UniqueRef,$INDELhetFN" >> $OUTPUTFILE
    echo "$file,INDEL,hom,Multiallelic,$INDELhomMulti" >> $OUTPUTFILE
    echo "$file,INDEL,het,Multiallelic,$INDELhetMulti" >> $OUTPUTFILE
    echo "$file,INDEL,het,MixedZygosity,$INDELzygo" >> $OUTPUTFILE
    if [ $PRINTOUTPUT -eq 1 ]
    then
        cat $OUTPUTFILE
        rm $OUTPUTFILE
    fi
done