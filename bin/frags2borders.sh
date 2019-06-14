#!/bin/bash

if [ "$#" -lt 5 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 fragsFILE.csv fastaX fastaY alignments.txt borderSize [query]"
   echo "	include "query" as ending parameter to extract also the query sequence as fasta file "
   echo ""
   exit -1
fi

FRAGS=$1
FASTAX=$2
FASTAY=$3
ALIGN=$4
BORDER=$5
ONLYQUERYSEQ=$6




BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# generate indices

$BINDIR/indexmaker $FASTAX $FASTAX.idx
$BINDIR/indexmaker $FASTAY $FASTAY.idx

$BINDIR/reverseComplement $FASTAY $FASTAY.rev
$BINDIR/indexmaker $FASTAY.rev $FASTAY.rev.idx

# extract frags

sed 's/,/ /g' $FRAGS > $FRAGS.fix

$BINDIR/csvExtractBorders $FRAGS.fix $FASTAX $FASTAY $FASTAY.rev $FASTAX.idx $FASTAY.idx $FASTAY.rev.idx $ALIGN $BORDER

if [[ $ONLYQUERYSEQ = "query" ]]
then

	grep '>\|X:' $ALIGN | awk 'BEGIN{n=0;} {if (substr($0,1,2) ~ "X:" ) { print substr($0,8);} else { printf(">S%d-%s\n", n, $0); n = n + 1; }}' > query.fasta
	sed 's/>/#/g' query.fasta | sed 's/^#/>/' > query-$BORDER.fasta
	rm query.fasta

fi


rm $FRAGS.fix
rm $FASTAX.idx
rm $FASTAY.idx
rm $FASTAY.rev
rm $FASTAY.rev.idx

