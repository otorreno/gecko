#!/bin/bash 

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# != 2 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta prefixSize(to split the work in the binary tree)"
   echo ""
   exit -1
fi

seqName=$(basename "$1")
extension="${seqName##*.}"
seqName="${seqName%.*}"

WL=32
prefixSize=$2

# create the dictionary
echo "${BINDIR}/dictionary ${seqName}.${extension} 32 ${prefixSize} ${seqName}"
${BINDIR}/dictionary ${seqName}.${extension} 32 ${prefixSize} ${seqName}
