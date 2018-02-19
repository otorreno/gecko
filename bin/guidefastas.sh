#!/bin/bash 

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -lt 7 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 seqXName.fasta seqYName.fasta guideFile ratio LEN SIM WL"
   echo ""
   exit -1
fi

seqNameX=$(basename "$1")
extension="${seqNameX##*.}"
seqNameX="${seqNameX%.*}"

seqNameY=$(basename "$2")
seqNameY="${seqNameY%.*}"

guided=$3
ratio=$4

LEN=$5
SIM=$6
WL=$7

lenX=$(wc -c $1 | awk '{print $1}')
lenY=$(wc -c $2 | awk '{print $1}')

counterX=0
counterY=0

turn=0

mkdir tempfastas
mkdir all-results

echo "All by-Identity Ungapped Fragments (Hits based approach)
[Abr.98/Apr.2010/Dec.2011 -- <ortrelles@uma.es>
SeqX filename        : undef
SeqY filename        : undef
SeqX name            : undef
SeqY name            : undef
SeqX length          : $lenX
SeqY length          : $lenY
Min.fragment.length  : undef
Min.Identity         : undef
Tot Hits (seeds)     : undef
Tot Hits (seeds) used: undef
Total fragments      : undef
========================================================
Total CSB: 0
========================================================
Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY" > all-results/master.csv


actualX=0
actualY=0

for i in $( tail -n +2 $guided ); do
	#echo "$i"

	if [[ turn -eq 0 ]]; then
	
		actualX=`expr $i - 1`
		actualX=$(($actualX * $ratio))
	
		echo ">nothing" > tempfastas/X_${counterX}.fasta
		(tail -c +"$actualX" "$1" | head -c "$ratio") >> tempfastas/X_${counterX}.fasta
		
	
		counterX=`expr $counterX + 1`
	fi
	
	if [[ turn -eq 1 ]]; then
		
		actualY=`expr $i - 1`
		actualY=$(($actualY * $ratio))
		
		echo ">nothing" > tempfastas/Y_${counterY}.fasta	
		(tail -c +"$actualY" "$2" | head -c "$ratio") >> tempfastas/Y_${counterY}.fasta
	
		counterY=`expr $counterY + 1`
	fi
	
	
	
	# switcher
	if [[ turn -eq 0 ]]; then
		turn=1
	else
		turn=0
		
		counterXprev=`expr $counterX - 1`
		counterYprev=`expr $counterY - 1`
		
		# run gecko on this comparison
		echo "$BINDIR/workflow.sh tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta $LEN $SIM $WL 1"
		$BINDIR/workflow.sh tempfastas/X_${counterXprev}.fasta tempfastas/Y_${counterYprev}.fasta $LEN $SIM $WL 1
		(tail -n +18 results/X_${counterXprev}-Y_${counterYprev}.csv | awk -F "," -v OFS=',' -v a="$actualX" -v b="$actualY" '{print $1,$2+a,$3+b,$4+a,$5+b,$6,$7,$8,$9,$10,$11,$12,$13,$14}') >> all-results/master.csv


		mv results/X_${counterXprev}-Y_${counterYprev}.csv all-results
		rm -rf intermediateResults/ results/
		
	fi
	
done
rm -rf intermediateResults results


