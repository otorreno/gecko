CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=../bin

all: reverseComplement dictionary hits sortHits filterHits FragHits combineFrags hitsStat fragStat getInfoCSB

hits: hits.c
	$(CC) $(CFLAGS) hits.c commonFunctions.c comparisonFunctions.c dictionaryFunctions.c -o $(BIN)/hits

sortHits: sortHits.c
	$(CC) $(CFLAGS) -DBaseType=hit sortHits.c comparisonFunctions.c commonFunctions.c quicksort.c -lpthread -o $(BIN)/sortHits

filterHits: filterHits.c
	$(CC) $(CFLAGS) filterHits.c comparisonFunctions.c commonFunctions.c -o $(BIN)/filterHits

FragHits: FragHits.c
	$(CC) $(CFLAGS) FragHits.c comparisonFunctions.c commonFunctions.c -o $(BIN)/FragHits
	
combineFrags: combineFrags.c
	$(CC) $(CFLAGS) combineFrags.c comparisonFunctions.c commonFunctions.c -o $(BIN)/combineFrags
	
hitsStat: hitsStat.c
	$(CC) $(CFLAGS) hitsStat.c commonFunctions.c -o $(BIN)/hitsStat
	
fragStat: fragStat.c
	$(CC) $(CFLAGS) fragStat.c comparisonFunctions.c commonFunctions.c -o $(BIN)/fragStat
	
reverseComplement: reverseComplement.c
	$(CC) $(CFLAGS) reverseComplement.c commonFunctions.c -o $(BIN)/reverseComplement

dictionary: dictionary.c
	$(CC) $(CFLAGS) -lm commonFunctions.c dictionary.c dictionaryFunctions.c -o $(BIN)/dictionary

getInfoCSB:
	$(CC) $(CFLAGS) getInfoCSB.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/getInfo

clean:
	rm -rf $(BIN)/getInfo
	rm -rf $(BIN)/reverseComplement $(BIN)/dictionary
	rm -rf $(BIN)/hits $(BIN)/sortHits $(BIN)/filterHits $(BIN)/FragHits $(BIN)/combineFrags $(BIN)/hitsStat $(BIN)/fragStat $(BIN)/hdStat
