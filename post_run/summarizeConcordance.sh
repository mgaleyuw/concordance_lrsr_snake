#!/bin/bash

# expects the summary metrics files as input

help() {
    echo "summarizeConcordance.sh helper script"
    echo ""
    echo "Append the genotype and non reference concordance values for many outputs to one output csv file."
    echo ""
    echo "Usage: bash summarizeConcordance.sh -o outputPathFile.csv -a comparison_1_input_.genotype_concordance_summary_metrics comparison_2_input_.genotype_concordance_summary_metrics ..."
    echo ""
    echo "-o: name of output file"
    echo "-a: include to append to output file rather than overwrite"
    echo "-h: display this help"
    echo ""
}


APPEND=0

while getopts "o:ah" opt; do
  case "$opt" in
    o) OUTPUTFILE=$OPTARG ;;
    a) APPEND=1 ;;
    h) help
       exit 1 ;;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

inputfiles=( "$@" )

header="File,TYPE,GenotypeConcordance,NonRefConcordance"

if [ $APPEND -eq 0 ]
then
    if [ -z ${OUTPUTFILE+x} ]
    then
        echo $header
    else
        echo $header > $OUTPUTFILE
    fi
fi

for file in ${inputfiles[@]}
do
    SNPCON=$(grep ^SNP $file | cut -f13 )
    SNPNRCON=$(grep ^SNP $file | cut -f14 )
    INDELCON=$(grep ^INDEL $file|  cut -f13 )
    INDELNRCON=$(grep ^INDEL $file | cut -f14 )

    if [ -z ${OUTPUTFILE+x} ]
    then
        echo "$file,SNP,$SNPCON,$SNPNRCON"
        echo "$file,INDEL,$INDELCON,$INDELNRCON"
    else
        echo "$file,SNP,$SNPCON,$SNPNRCON" >> $OUTPUTFILE
        echo "$file,INDEL,$INDELCON,$INDELNRCON" >> $OUTPUTFILE
    fi

done