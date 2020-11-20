# GECKO

A pairwise genome comparison software for the detection of High-scoring Segment Pairs.

GECKO (GEnome Comparison with K-mers Out-of-core) is a fast, modular application designed to identify collections of High-scoring Segment Pairs in a pairwise genome comparisons. By employing novel filtering and data storing strategies, it is able to compare chromosome-sized sequences in less time.

## Requirements

GCC compiler. Please make sure that on your linux CC resolves to GCC, otherwise it might not work.

Simply download the .zip and unzip it, or clone the repository.
Then issue the following command:

```cd gecko/src && make all```

If the installation finished without errors, you are ready to go!

## Use

You can run GECKO directly by issuing:

```./gecko/bin/workflow.sh query reference length similarity wordLength 1```

Please note the ending "1". This is used for internal developing. Do not change unless actively knowing why.

## Parameters

- Query sequence: The sequence that will be compared against the reference. Use only FASTA format.
- Reference sequence: The reference sequence where to look for matches from the query. Note that the reverse strand is computed for the reference and also matched. Use only FASTA format.
- Length: This parameter is the minimum length in nucleotides for an HSP (similarity fragment) to be conserved. Any HSP below this length will be filtered out of the comparison. It is recommended to use around 40 bp for small organisms (e.g. bacterial mycoplasma or E. Coli) and around 100 bp or more for larger organisms (e.g. human chromosomes).
- Similarity: This parameter is analogous to the minimum length, however, instead of length, the similarity is used as threshold. The similarity is calculated as the score attained by an HSP divided by the maximum possible score. Use values above 50-60 to filter noise.
- Word length: This parameter is the seed size used to find HSPs. A smaller seed size will increase sensitivity and decrease performance, whereas a larger seed size will decrease sensitivity and increase performance. Recommended values are 12 or 16 for smaller organisms (bacteria) and 32 for larger organisms (chromosomes). These values must be multiples of 4.

## Output

Each time GECKO is executed, a folder structure is created containing both intermediate files and results, see below:

```
current_directory --|
                    + intermediateFiles --|
                                          + dictionaries [...]
                                          + hits         [...]

                    + results ------------|
                                          + query-reference.csv
                                          + query-reference.frags

```

- The "dictionaries" folder contains files that are used in case GECKO is run again with the same sequences, either as query or reference. 
- The "hits" folder contains the files regarding seed matches between query and reference. This folder can be heavyweight and the recommendation is to remove it after usage.- The "results" folder contains the files that represent the comparison output. In particular:
 - query-reference.csv: A CSV file that includes metadata about the sequence compared and each HSP detected. See section "Interpreting the CSV" below for more information. This file can be used to visualize the comparison in the interactive sequence visualizer GECKO-MGV (use online [here](https://pistacho.ac.uma.es/) or download and install [here](https://github.com/estebanpw/docker-geckomgv)). 
 - query-reference.frags: The same information contained in the CSV but without metadata and in binary format. This file is used for further processing (e.g. extracting the alignments). This file is made up from the `FragFile` structure found in the `structs.h` file in the `src` folder.

## Interpreting the CSV

After the metadata header, each line in the CSV represents a detected HSP (alignment). Each column is:

``` Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY ```

- Type:   currently, this field is reserved for `Frag`.
- xStart: starting coordinates of the alignment in the query sequence.
- yStart: starting coordinates of the alignment in the reference sequence.
- xEnd:   ending coordinates of the alignment in the query sequence.
- yEnd:   ending coordinates of the alignment in the reference sequence.
- strand: a character `f` or `r` encoding whether the alignment is in the forward or reverse strand.
- block:  currently reserved.
- length: the length in nucleotides of the alignment.
- score:  the raw score of the alignment calculated with +4 and -4 per match and mismatch.
- ident:  the number of identities found in the alignment (i.e. matches).
- similarity: the similarity percentage calculated as the achieved raw score divided by the maximum possible score.
- %ident:     the number of identities divided by the length of the alignment.
- SeqX:       the ID corresponding to the sequence in the query file to which the xStart and xEnd coordinates correspond (0=> first sequence, 1=> second sequence, etc).
- SeqY:       same as above but for the reference file.

Note that fragments in the reverse strand (marked with the `r` field) have their `yStart` and `yEnd` coordinates switched, i.e. `yEnd` is smaller than `yStart`.

## Extracting the alignments

Use `./gecko/bin/frags2align.sh fragsFILE.frags query reference output_alignments.txt` to extract the alignments from a comparison. The `query` and `reference` files must be the same FASTA input sequences as used to generate the `fragsFILE.frags` found in the results folder. This will create a text file `output_alignments.txt` containing the alignments.

## Further documentation

You can find a guided exercise [here](http://chirimoyo.ac.uma.es/gecko/documents/GuidedExercise-fromGENOMES2Visualization-reduced-v0.2.pdf) and the supplementary material [here](http://chirimoyo.ac.uma.es/gecko/documents/HSPWorkflow-SuppMat-submittedv2.pdf).



