if [ $# != 4 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 ramgecko.frags ramgecko_converted.frags query_file.fasta database_file.fasta"
   echo ""
   exit -1
fi

originalfrags=$1
convertedfrags=$2
query=$3
database=$4

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Perform conversion

echo "## Conversion step 1: Creating indexes"

${BINDIR}/indexmaker $query $query.IND2
${BINDIR}/indexmaker $database $database.IND2

echo "## Conversion step 2: Converting coordinates"

${BINDIR}/ramgecko2gecko $originalfrags $convertedfrags $query.IND2 $database.IND2

#rm $query.IND2 $query.IND2.sorted $database.IND2 $database.IND2.sorted

echo "## Finished. Converted fragments file at $convertedfrags"


