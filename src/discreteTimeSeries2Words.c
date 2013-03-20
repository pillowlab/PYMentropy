/*
 * Convert an array of symbols over time to a single stream of concise alphabet
 */

#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#include "fastWords.h"

static int ID; /* global variable tracking the current alphabet ID */

int getID(Node** root, wordType* word, int nAlphabet, int wordLength) {
    int k;
    Node *current = *root;

    if (!current) {
	current = (Node *) calloc(1, sizeof(Node));
	*root = current;
    }

    for(k = 0; k < wordLength; k++) {
	if(current->children == NULL) {
	    current->children = (Node**) calloc(nAlphabet, sizeof(Node));
	}
	if(current->children[word[k]] == NULL) {
	    current->children[word[k]] = (Node*) calloc(1, sizeof(Node));
	    
	    /* assign a new label to the leaf node */
	    if(k == wordLength-1) current->children[word[k]]->count = ID++;
	}
	current = current->children[word[k]];
    }

    return current->count;
}

unsigned int* discrete2words(wordType* words, int nAlphabet, int wordLength, int nWords) 
{
    int k;
    Node *root = NULL;
    unsigned int *list = (unsigned int *) calloc(nWords, sizeof(unsigned int));

    ID = 0; /* initialize the global ID */

    for(k = 0; k < nWords; k++) {
	list[k] = getID(&root, &words[k*wordLength], nAlphabet, wordLength);
	/* printf("%d ", list[k]); */
    }
    /* printf("\n"); */
    freeTree(root, 0, nAlphabet, wordLength);

    return list;
}

/*
int main() {
    size_t wordLength = 4;
    size_t nAlphabet = 3;
    unsigned int *list;
    size_t k;
    wordType words[][4] = {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 2, 2}, {0, 0, 2, 2}};
    // expecting 0 1 0 2 2

    list = discrete2words(&words[0][0], nAlphabet, wordLength, 5);

    free(list);

    return 1;
}
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* First argument is a regular array of words [m x n] where
     * m: the wordLength
     * n: number of words
     *
     * the array MUST be in uint16 (unsigned short) - to change this restriction
     *	    change the wordType definition, and recompile the mex.
     */
    mwSize wordLength, nWords;
    double *p;
    wordType *words;
    int nAlphabet;
    int k;
    unsigned int *pCounts;
    mwSignedIndex dims[2];
    unsigned int *list;
    
    if(nrhs < 2) mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "At least two arguments are required.");

    /* parse the first argument: words */
    wordLength = mxGetM(prhs[0]);
    nWords = mxGetN(prhs[0]);
    /* printf("wordLength [%d], nWords [%d]\n", wordLength, nWords); */

    if(nWords < 1) mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "at least one word needed");
    if(wordLength < 1) mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "word length has to be at least 1");

    if(mxIsSparse(prhs[0])) {
	mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "This is sparse, man...(not yet implemented!) \n");
    } else {
	if(mxGetClassID(prhs[0]) != mxUINT16_CLASS) {
	    mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "Words must be in UINT16 or sparse. Try uint16(words) [%d].", mxGetClassID(prhs[0]));
	}
    }

    /* parse the second argument: nAlphabet */
    p = mxGetPr(prhs[1]);
    nAlphabet = (int) *p;
    /* printf("nAlphabet [%d]\n", nAlphabet); */
    if(nAlphabet <= 1) mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "Number of alphabets should be at least 2 [%d]", nAlphabet);

    words = (wordType*) mxGetPr(prhs[0]);
    /* Build the tree and renumber */
    list = discrete2words(words, nAlphabet, wordLength, nWords);

    /* Generate the histogram for returning */
    dims[0] = nWords;
    dims[1] = 1;
    plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);

    if(!plhs[0]) mexErrMsgIdAndTxt("fastWords2Counts:retrunErr", "failed creating return array");

    pCounts = (unsigned int*) mxGetData(plhs[0]);
    if(!pCounts) mexErrMsgIdAndTxt("fastWords2Counts:retrunErr", "failed creating return array contents???");

    for(k = 0; k < nWords; k++) pCounts[k] = list[k];

    /* clean up */
    free(list);
}
