#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <ctype.h>

#include "commonFunctions.h"
#include "structs.h"

#include "dictionaryFunctions.h"

#define NODEPOOLSIZE 4000000
#define LOCATIONPOOLSIZE 4000000

#define MIN(a,b) (((a)<(b))?(a):(b))

//Internal structures to work with a memory pool
typedef struct mn {
	node *data;
	struct mn *next;
} treeMemory_node;

typedef struct ln {
	list_node *data;
	struct ln *next;
} locationMemory_node;

unsigned long time() {
	struct timeval tv;

	gettimeofday(&tv, NULL);

	return (tv.tv_usec / 1000) + (tv.tv_sec * 1000);
}

void hashValue2Word(int index, int K, char *word) {
	int resto;
	int cociente = index;
	char alph[] = { 'A', 'C', 'G', 'T' };
	int i = 0;

	for (i = 0; i < K; i++)
		word[i] = 'A';

	i = 0;
	while (cociente >= 4) { // 4 is the alphabet size
		resto = cociente % 4;
		cociente = cociente / 4;
		word[K - i - 1] = alph[resto];
		i++;
	}
	word[K - i - 1] = alph[cociente];
}

int nextPrefix(int prefixSize, int prefixIndex, char *word) {
	//Out of range?
	if (prefixIndex >= pow(4, prefixSize) || (prefixSize == 0)) {
		return -1;
	}

	//Not really..
	hashValue2Word(prefixIndex, prefixSize, word);
	return 0;
}

int main(int ac, char **av) {
	FILE *f, *fout;
	//unsigned long startt, stop;
	unsigned long i = 0, s = 0, nWords = 0;
	char *seq, fname[1024], *prefix;
	char c;
	int n = 0, Tot = 0, NoACGT = 0, NoC = 0, BAD, index, wsize, prefixSize,
			prefixIndex = 0;
	location loc;
	node *start = NULL;

	node *currentTreeRow;
	list_node *currentLocationRow;

	if (ac != 5) {
		terror(
				"USE: dictionary seq.fasta wsize prefixSize(for splitting work) prefixOutFile");
	}

	//Node Memory pool
	int actualNodes = 0;
	treeMemory_node *nodePool = (treeMemory_node *) malloc(
			sizeof(treeMemory_node));
	currentTreeRow = nodePool->data = (node *) malloc(
			NODEPOOLSIZE * sizeof(node));
	nodePool->next = NULL;
	memset(nodePool->data, 0, NODEPOOLSIZE * sizeof(node));

	//Location Memory pool
	int actualLocations = 0;
	locationMemory_node *locationPool = (locationMemory_node *) malloc(
			sizeof(locationMemory_node));
	currentLocationRow = locationPool->data = (list_node *) malloc(
			LOCATIONPOOLSIZE * sizeof(list_node));
	locationPool->next = NULL;
	memset(locationPool->data, 0, LOCATIONPOOLSIZE * sizeof(list_node));

	//fprintf(stdout, "Memory pool initialized\n");

	wsize = BAD = atoi(av[2]);
	if ((seq = (char*) malloc(wsize + 1)) == NULL)
		perror("memory for Seq");
	seq[wsize] = '\0';
	if ((f = fopen(av[1], "rt")) == NULL)
		perror("opening sequence file");
	prefixSize = atoi(av[3]);
	prefix = malloc(prefixSize + 1);
	prefix[prefixSize] = '\0';
	sprintf(fname, "%s.dict", av[4]);
	if ((fout = fopen(fname, "wb")) == NULL)
		perror("opening words file");

	//startt=time();

	while (nextPrefix(prefixSize, prefixIndex, prefix) != -1) {
		// first line (skip ID fasta Line)
		c = fgetc(f);
		while (c != '\n')
			c = fgetc(f);

		c = toupper(fgetc(f));
		while (!feof(f)) {
			if (c < 65 || c > 90) {
				NoC++;
				if (c != '>') {
					c = toupper(fgetc(f));
					continue;
				}
				// coding (a=0,c=1,g=2,t=3, others=9 )
			}
			index = MIN(Tot, wsize - 1);
			seq[index] = c;
			switch (c) {
			case 'A':
				BAD--;
				break;
			case 'C':
				BAD--;
				break;
			case 'G':
				BAD--;
				break;
			case 'T':
				BAD--;
				break;
			case '>': //skip ID fasta Line
				while (c != '\n') {
					c = fgetc(f);
				}

				s++;
				Tot = -1;
				BAD = wsize;
				break;
			default:
				NoACGT++;
				BAD = wsize;
				break;
			}
			if (Tot > (wsize - 2) && !BAD) {
				loc.pos = n - wsize + 1;
				loc.seq = s;

				if (strncmp(seq, prefix, prefixSize) == 0) {
					//Corresponds to the actual interval
					//fprintf(stdout, "inserting word: %s Address: %p, actualNodes: %d, actualLocations: %d\n",seq,nodePool[nodePoolIndex],actualNodes,actualLocations);
					nWords += insert(seq, wsize, loc, &start, currentTreeRow,
							&actualNodes, currentLocationRow, &actualLocations);
					//fprintf(stdout, "inserted word: %s, actualNodes: %d, actualLocations: %d\n",seq, actualNodes, actualLocations);

					//Check if the node pool is full
					if (actualNodes > NODEPOOLSIZE) {
						//fprintf(stdout, "Node pool full\n");
						treeMemory_node *aux = nodePool;
						while (aux->next != NULL) {
							aux = aux->next;
						}
						aux->next = (treeMemory_node *) malloc(
								sizeof(treeMemory_node));
						aux->next->data = (node *) malloc(
								NODEPOOLSIZE * sizeof(node));
						aux->next->next = NULL;
						currentTreeRow = aux->next->data;
						actualNodes = 0;
					}

					if (actualLocations > LOCATIONPOOLSIZE) {
						//fprintf(stdout, "Location pool full\n");
						locationMemory_node *aux = locationPool;
						while (aux->next != NULL) {
							aux = aux->next;
						}
						aux->next = (locationMemory_node *) malloc(
								sizeof(locationMemory_node));
						aux->next->data = (list_node *) malloc(
								NODEPOOLSIZE * sizeof(list_node));
						aux->next->next = NULL;
						currentLocationRow = aux->next->data;
						actualLocations = 0;
					}
				}

				for (i = 1; i < wsize; i++) {
					seq[i - 1] = seq[i];
				}
				BAD = 1;
			}
			Tot++;
			n++;
			c = toupper(fgetc(f));
		}
		//stop=time();

		//fprintf(stdout, "All words in the tree. Time: %lu\n",stop-startt);
		//startt=time();
		write(start, fout);
		//stop=time();
		//fprintf(stdout, "Writing the words. Time: %lu\n",stop-startt);

		//Destroy the tree and rewind the input file to start from zero
		//destroy_tree(&start);

		//fprintf(stdout, "Tree destroyed\n");

		//fprintf(stdout, "Freeing up %d rows in nodePool and %d in locationPool\n",nodePoolIndex-1,locationPoolIndex-1);

		//Free up memory from the pool
		if (nodePool->next != NULL) {
			//Free all but the first row of the pool matrix
			treeMemory_node *aux = nodePool->next;
			treeMemory_node *ant = nodePool->next;
			nodePool->next = NULL;
			do {
				aux = aux->next;
				free(ant->data);
				free(ant);
				ant = aux;
			} while (aux != NULL);
			currentTreeRow = nodePool->data;
		}

		if (locationPool->next != NULL) {
			//Free all but the first row of the pool matrix
			locationMemory_node *aux = locationPool->next;
			locationMemory_node *ant = locationPool->next;
			locationPool->next = NULL;
			do {
				aux = aux->next;
				free(ant->data);
				free(ant);
				ant = aux;
			} while (aux != NULL);
			currentLocationRow = locationPool->data;
		}
		memset(nodePool->data, 0, NODEPOOLSIZE * sizeof(node));
		memset(locationPool->data, 0, LOCATIONPOOLSIZE * sizeof(list_node));

		//fprintf(stdout, "Pool matrix free\n");

		start = NULL;
		//fprintf(stdout, "Tree destroyed\n");
		BAD = wsize;
		Tot = n = s = actualNodes = actualLocations = 0;
		prefixIndex++;
		rewind(f);
	}

	fprintf(stdout, "Number of words: %lu\n", nWords);

	free(seq);
	free(prefix);
	fclose(fout);
	fclose(f);
	return 0;
}
