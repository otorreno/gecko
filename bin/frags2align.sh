#!/bin/bash

if [ $# != 4 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 fragsFILE.frags/.csv fastaX fastaY alignments.txt"
   echo ""
   exit -1
fi

FRAGS=$1
FASTAX=$2
FASTAY=$3
ALIGN=$4



BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# generate indices

$BINDIR/indexmaker $FASTAX $FASTAX.idx
$BINDIR/indexmaker $FASTAY $FASTAY.idx

$BINDIR/reverseComplement $FASTAY $FASTAY.rev
$BINDIR/indexmaker $FASTAY.rev $FASTAY.rev.idx

# extract frags

if [ ${FRAGS: -4} == ".csv" ]
then

	echo "....using csv"
	sed 's/,/ /g' $FRAGS > $FRAGS.fix
	$BINDIR/csvFrags2text $FRAGS.fix $FASTAX $FASTAY $FASTAY.rev $FASTAX.idx $FASTAY.idx $FASTAY.rev.idx $ALIGN
	rm $FRAGS.fix

else
	echo "....using frags"
	$BINDIR/frags2text $FRAGS $FASTAX $FASTAY $FASTAY.rev $FASTAX.idx $FASTAY.idx $FASTAY.rev.idx $ALIGN
fi


rm $FASTAX.idx
rm $FASTAY.idx
rm $FASTAY.rev
rm $FASTAY.rev.idx

