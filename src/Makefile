CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=../bin

all: reverseComplement words sortWords w2hd wordsStat hits sortHits filterHits FragHits combineFrags hitsStat hdStat fragStat getInfoCSB

png: af2pngrev aMf2pngrev

hits: hits.c
	$(CC) $(CFLAGS) hits.c comparisonFunctions.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/hits

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
	
hdStat: hdStat.c
	$(CC) $(CFLAGS) hdStat.c commonFunctions.c dictionaryFunctions.c -o $(BIN)/hdStat
	
fragStat: fragStat.c
	$(CC) $(CFLAGS) fragStat.c comparisonFunctions.c commonFunctions.c -o $(BIN)/fragStat
	
af2pngrev: af2pngrev.cpp
	$(CXX) $(CFLAGS) -Wno-write-strings af2pngrev.cpp commonFunctions.c comparisonFunctions.c -I/usr/include/freetype2  -L/usr/local/lib -lpngwriter -lz -lm -lfreetype -o $(BIN)/af2pngrev -lpng 

aMf2pngrev: aMf2pngrev.cpp
	$(CXX) $(CFLAGS) -Wno-write-strings aMf2pngrev.cpp commonFunctions.c comparisonFunctions.c -I/usr/include/freetype2  -L/usr/local/lib -lpngwriter -lz -lm -lfreetype -o $(BIN)/aMf2pngrev -lpng 

reverseComplement: reverseComplement.c
	$(CC) $(CFLAGS) reverseComplement.c commonFunctions.c -o $(BIN)/reverseComplement

words: words.c
	$(CC) $(CFLAGS) words.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/words

sortWords: sortWords.c
	$(CC) $(CFLAGS) -DBaseType=wentry sortWords.c dictionaryFunctions.c commonFunctions.c quicksort.c -lpthread -o $(BIN)/sortWords

w2hd: w2hd.c
	$(CC) $(CFLAGS) w2hd.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/w2hd

wordsStat: wordsStat.c
	$(CC) $(CFLAGS) wordsStat.c dictionaryFunctions.c commonFunctions.c -o $(BIN)/wordsStat

sortFragsBySeqX: sortFragsBySeqX.c
	$(CC) $(CFLAGS) -D'BaseType=struct FragFile' sortFragsBySeqX.c commonFunctions.c quicksort.c comparisonFunctions.c -lpthread -o $(BIN)/sortFragsBySeqX

getInfoCSB:
	$(CC) $(CFLAGS) getInfoCSB.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/getInfo

clean:
	rm -rf $(BIN)/getInfo
	rm -rf $(BIN)/reverseComplement $(BIN)/words $(BIN)/sortWords $(BIN)/w2hd $(BIN)/wordsStat
	rm -rf $(BIN)/hits $(BIN)/sortHits $(BIN)/filterHits $(BIN)/FragHits $(BIN)/combineFrags $(BIN)/hitsStat $(BIN)/fragStat $(BIN)/hdStat $(BIN)/af2pngrev $(BIN)/aMf2pngrev
