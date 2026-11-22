#CC=gcc
#CXX=g++
CC?=gcc
NVCC?=nvcc
CUDA_ARCH?=-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_89,code=sm_89
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=../bin

all: reverseComplement words sortWords w2hd wordsStat hits sortHits filterHits FragHits combineFrags hitsStat hdStat fragStat getInfoCSB filterFrags indexmaker frags2text csvFrags2text csvExtractBorders gecko_fast

gecko_fast: gecko_fast.c
	$(CC) $(CFLAGS) -Wextra -march=native gecko_fast.c comparisonFunctions.c commonFunctions.c -lpthread -lm -o $(BIN)/gecko_fast

CUDA_SRCS=cuda/gecko_cuda_main.cu cuda/cuda_io.cu cuda/cuda_index.cu cuda/cuda_query.cu cuda/cuda_extend.cu
gecko_cuda: $(CUDA_SRCS) cuda/cuda_common.cuh
	$(CC) $(CFLAGS) -c comparisonFunctions.c -o cuda/comparisonFunctions.o
	$(CC) $(CFLAGS) -c commonFunctions.c -o cuda/commonFunctions.o
	$(NVCC) -O3 -std=c++17 $(CUDA_ARCH) -D_FILE_OFFSET_BITS=64 -I. $(CUDA_SRCS) cuda/comparisonFunctions.o cuda/commonFunctions.o -o $(BIN)/gecko_cuda

gecko_wgpu:
	cargo build --release --manifest-path ../gpu-wgpu/Cargo.toml
	cp ../gpu-wgpu/target/release/gecko_wgpu $(BIN)/gecko_wgpu

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

getInfoCSB:
	$(CC) $(CFLAGS) getInfoCSB.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/getInfo

filterFrags:
	$(CC) $(CFLAGS) filterFrags.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/filterFrags	

indexmaker:
	$(CC) $(CFLAGS) indexmaker.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/indexmaker

frags2text:
	$(CC) $(CFLAGS) frags2text.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/frags2text

csvFrags2text:
	$(CC) $(CFLAGS) csvFrags2text.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/csvFrags2text

csvExtractBorders:
	$(CC) $(CFLAGS) csvExtractBorders.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/csvExtractBorders

clean:
	rm -rf $(BIN)/getInfo
	rm -rf $(BIN)/filterFrags
	rm -rf $(BIN)/indexmaker
	rm -rf $(BIN)/frags2text
	rm -rf $(BIN)/reverseComplement $(BIN)/words $(BIN)/sortWords $(BIN)/w2hd $(BIN)/wordsStat
	rm -rf $(BIN)/hits $(BIN)/sortHits $(BIN)/filterHits $(BIN)/FragHits $(BIN)/combineFrags $(BIN)/hitsStat $(BIN)/fragStat $(BIN)/hdStat
	rm -rf $(BIN)/csvExtractBorders $(BIN)/csvFrags2text.c $(BIN)/gecko_fast $(BIN)/gecko_cuda $(BIN)/gecko_wgpu cuda/*.o
