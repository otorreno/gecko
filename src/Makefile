CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DELAPSEDTIME
BIN=../bin

all: gecko wordsStat hitsStat ramgecko2gecko indexmaker

stats: wordsStat hitsStat

coordinates: ramgecko2gecko indexmaker

gecko: gecko.c
	$(CC) $(CFLAGS) dictionaryFunctions.c -lm commonFunctions.c -lm comparisonFunctions.c quicksortWordForward.c quicksortWordReverse.c quicksortHit.c gecko.c -lpthread -o $(BIN)/gecko

wordsStat: wordsStat.c
	$(CC) $(CFLAGS) dictionaryFunctions.c -lm commonFunctions.c comparisonFunctions.c -lm quicksortWordForward.c quicksortWordReverse.c quicksortHit.c wordsStat.c -lpthread -o $(BIN)/wordsStat

hitsStat: hitsStat.c
	$(CC) $(CFLAGS) dictionaryFunctions.c -lm commonFunctions.c comparisonFunctions.c -lm quicksortWordForward.c quicksortWordReverse.c quicksortHit.c hitsStat.c -lpthread -o $(BIN)/hitsStat

ramgecko2gecko: ramgecko2gecko.c
	$(CC) $(CFLAGS) dictionaryFunctions.c -lm commonFunctions.c comparisonFunctions.c -lm quicksortWordForward.c quicksortWordReverse.c quicksortHit.c ramgecko2gecko.c -lpthread -o $(BIN)/ramgecko2gecko

indexmaker: indexmaker.c
	$(CC) $(CFLAGS) dictionaryFunctions.c -lm commonFunctions.c comparisonFunctions.c -lm quicksortWordForward.c quicksortWordReverse.c quicksortHit.c indexmaker.c -lpthread -o $(BIN)/indexmaker



clean:
	rm -rf $(BIN)/gecko $(BIN)/wordsStat $(BIN)/hitsStat $(BIN)/ramgecko2gecko $(BIN)/indexmaker
