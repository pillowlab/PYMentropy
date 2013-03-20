#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#include "fastWords.h"

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
    int *counts = NULL;
    int countsLength;
    int k;
    Node *root = NULL;
    unsigned int *pCounts;
    mwSignedIndex dims[2];
    
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
    if(mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) {
	mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "nAlphabert should be in double [%d], sorry!", mxGetClassID(prhs[1]));
    }
    nAlphabet = (int) *p;
    /* printf("nAlphabet [%d]\n", nAlphabet); */
    if(nAlphabet <= 1) mexErrMsgIdAndTxt("fastWords2Counts:invalidInput", "Number of alphabets should be at least 2 [%d]", nAlphabet);

    words = (wordType*) mxGetPr(prhs[0]);
    /* Build the tree and count */
    root = buildTree(words, nAlphabet, wordLength, nWords);
    /* retrieve the leaves which gives the counts */
    getHistogram(&counts, &countsLength, root, nAlphabet, wordLength);

    /*
    printf("histogram size %d\n", countsLength);
    printf("[");
    for(k = 0; k < countsLength; k++) {
	printf("%d ", counts[k]);
    }
    printf("\b]\n");
    */

    /* Generate the histogram for returning */
    dims[0] = countsLength;
    dims[1] = 1;
    plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);

    if(!plhs[0]) mexErrMsgIdAndTxt("fastWords2Counts:retrunErr", "failed creating return array");

    pCounts = (unsigned int*) mxGetData(plhs[0]);
    if(!pCounts) mexErrMsgIdAndTxt("fastWords2Counts:retrunErr", "failed creating return array contents???");

    for(k = 0; k < countsLength; k++) pCounts[k] = counts[k];

    /* clean up */
    freeTree(root, 0, nAlphabet, wordLength);
}
