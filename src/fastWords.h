typedef struct treeNode Node;
typedef unsigned short wordType; /* something that allows unsigned integers */

struct treeNode {
    Node **children;
    size_t count; /* number of occurrences */
};

Node* buildTree(wordType* words, int nAlphabet, int wordLength, int nWords);
void freeTree(Node* current, int depth, int nAlphabet, int wordLength);
void getHistogram(int **counts, int *countsLength, Node* tree, int nAlphabet, int wordLength);
