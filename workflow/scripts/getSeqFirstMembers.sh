#!/bin/bash

while getopts "d:" option; do
  case $option in
    d) INDIR="$OPTARG";;
  esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

families=( "$@" )

for family in ${families[@]}
do
    members=( $(ls -d1 /n/alignments/seqFirst-ONT-hg38/$family/M*) )
    membernames=()
    for member in ${members[@]}
    do
        prefix=${member##*/}
        membernames+=( ${prefix%%-*} )
    done
    set ${membernames[@]}
    for a; do
    shift
       for b; do
           printf "%s %s\n" "$a" "$b"
       done
    done
done