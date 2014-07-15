/*
 * Convert an array of symbols to multiplicities representation efficiently.
 */

#include <stdlib.h>
#include <stdio.h>
#include "fastWords.h"
#include "mex.h"

void addWord(Node** root, wordType* word, int nAlphabet, int wordLength) {
    int k;
    Node *current = *root;

    if (!current) {
	/* printf("creating the root node.\n"); */
	current = (Node *) calloc(1, sizeof(Node));
	*root = current;
    }

    for(k = 0; k < wordLength; k++) {
	current->count++;
	/* printf("inc [depth = %d][count = %d]\n", k, current->count); */

	if(current->children == NULL) {
	    current->children = (Node**) calloc(nAlphabet, sizeof(Node));
	}
	/* printf("current word [%d]\n", word[k]); */
	if(word[k] >= nAlphabet) {
	    mexErrMsgIdAndTxt("fastWordsTree:critical", "Word contains a value greater the # of alphabets!");
	}

	if(current->children[word[k]] == NULL) {
	    /* printf("allocating a child\n"); */
	    current->children[word[k]] = (Node*) calloc(1, sizeof(Node));
	}
	current = current->children[word[k]];
    }
    current->count++;
    /* printf(">>> Leaf node count [%d]\n", current->count); */
}

Node* buildTree(wordType* words, int nAlphabet, int wordLength, int nWords) {
    int k;
    Node* root = (Node*) calloc(1, sizeof(Node*));
    /* printf("buildTree [%p]\n", words); */

    for(k = 0; k < nWords; k++) {
	/* printf("buildTree [%d] -> [%p]\n", k, &words[k*wordLength]); */
	addWord(&root, &words[k*wordLength], nAlphabet, wordLength);
    }

    return root;
}

void dfsRecursion(int *leafCount, Node** firstLeaf, Node** tail, Node* current, int depth, int nAlphabet, int wordLength) {

    int k;
    if(depth == wordLength) {
	/* printf(">>> Leaf node count [%d] at depth [%d]\n", current->count, depth); */
	if(*firstLeaf == NULL) {
	    *firstLeaf = current;
	    /* printf(">>>>> saving the first leaf node [%p]\n", current); */
	    *tail = current;
	}
	current->children = (Node**) calloc(1, sizeof(Node));
	current->children[0] = (Node*) calloc(1, sizeof(Node));
	current->children[0] = *tail; /* make the linked list */
	/* printf(">>>>> tail node [%p][%p]\n", current, *tail); */
	*tail = current;
	(*leafCount)++;
	/* printf("<< going up one [leaf count %d]\n", *leafCount); */
	return;
    }

    if(current->children) {
	for(k = 0; k < nAlphabet; k++) {
	    if(current->children[k]) {
		/* printf("<< going down one [%d]\n", k); */
		dfsRecursion(leafCount, firstLeaf, tail, current->children[k], depth+1, nAlphabet, wordLength);
	    }
	}
    } else {
	printf("ASSERTION failed: no children...?\n");
    }
    /* printf("<< going up one\n"); */
}

void getHistogram(int **counts, int *countsLength, Node* tree, int nAlphabet, int wordLength) {
    /* depth first search and collect the leaf nodes */
    int leafCount = 0;
    int k;
    Node* current = NULL;
    Node* firstLeaf = NULL;
    Node* tail = NULL;

    dfsRecursion(&leafCount, &firstLeaf, &tail, tree, 0, nAlphabet, wordLength);
    /* printf("Total leaves [%d] [%p]\n", leafCount, firstLeaf); */

    current = tail;
    *counts = (int*) calloc(leafCount, sizeof(int));
    for(k = 0; k < leafCount; k++) {
	/* printf(">>>>> current node [%p]\n", current); */
	(*counts)[k] = current->count;
	if(current->children == NULL) {
	    printf("Assertion failed. Where's the next node?\n");
	}
	current = (current->children[0]);
    }
    *countsLength = leafCount;
}

void freeTree(Node* current, int depth, int nAlphabet, int wordLength) {
    int k;
    if(depth == wordLength) {
	if(current->children) free(current->children);
	/* printf("leaf [%p]\n", current); */
	return;
    }
    for(k = 0; k < nAlphabet; k++) {
	if(current->children[k]) {
	    /* printf("children [%p]\n", current->children[k]); */
	    freeTree(current->children[k], depth+1, nAlphabet, wordLength);
	}
    }
    /* printf("children list [%p]\n", current->children); */
    if(current->children) free(current->children);
    /* printf("current [%p]\n", current); */
    if(current) free(current);
}

/*
int main() {
    Node* root = NULL;
    size_t wordLength = 4;
    size_t nAlphabet = 3;
    int *counts = NULL;
    int countsLength;
    size_t k;
    wordType words[][4] = {{0, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 2, 2}, {0, 0, 2, 2}};
    // expecting 2 1 2

    root = buildTree(&words[0][0], nAlphabet, wordLength, 5);

    getHistogram(&counts, &countsLength, root, nAlphabet, wordLength);
    printf("histogram size %d\n", countsLength);
    printf("[");
    for(k = 0; k < countsLength; k++) {
	printf("%d ", counts[k]);
    }
    printf("\b]\n"); fflush(NULL);

    freeTree(root, 0, nAlphabet, wordLength);

    return 1;
}
*/
