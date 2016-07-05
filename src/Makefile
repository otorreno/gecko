CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=../bin

all: gecko

gecko: gecko.c
	$(CC) $(CFLAGS) dictionaryFunctions.c commonFunctions.c comparisonFunctions.c quicksortWordForward.c quicksortWordReverse.c quicksortHit.c gecko.c -lpthread -o $(BIN)/gecko

clean:
	rm -rf $(BIN)/gecko