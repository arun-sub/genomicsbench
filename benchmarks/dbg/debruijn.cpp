#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <errno.h>
#include <algorithm>
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include "omp.h"
#include "htslib/sam.h"
#include "common.h"
#include "htslib/faidx.h"

// #define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

inline char _getBase(uint8_t *s, int i) {
    char* baseLookup = "=ACMGRSVTWYHKDBN";
    return baseLookup[bam_seqi(s, i)];
}      

inline int Read_IsQCFail(struct alignedRead* theRead) {
    return ( (theRead->flag & BAM_FQCFAIL) != 0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int REF = 1;
int READ = 2;
int REF_AND_READ = 3;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Encapsulate information related to a node of the de-Bruijn graph,
// i.e. a sequence k-mer and related data.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward declaration of Edge struct, so it can be referred to in Node.
typedef struct Edge Edge;
typedef struct EdgeStack EdgeStack;

// Represents a node in the graph
typedef struct {
    Edge* edges[4];
    char* sequence;
    int colours;
    int position;
    int kmerSize;
    int nEdges;
    double weight;
    char dfsColour;
} Node;
    
// Simple implementation of a stack, for storing nodes.
typedef struct {
    Node** elements;
    int capacity;
    int top;
} NodeStack;

// Represents an edge in the graph
struct Edge {
    Node* startNode;
    Node* endNode;
    double weight;
};

// A stack of edges
struct EdgeStack {
    Edge** elements;
    int capacity;
    int top;
};

// A dictionary of Nodes
typedef struct {
    Node*** buckets;
    int* bucketSize;
    int nBuckets;
} NodeDict;

// Hold a path through the graph
typedef struct {
    NodeStack* nodes;
    int nNodes;
    int isBubble;
    double weight;
} Path;

// A stack of paths
typedef struct {
    Path** elements;
    int capacity;
    int top;
} PathStack;

// A graph
typedef struct {
    int kmerSize;
    NodeStack* allNodes;
    NodeDict* nodes;
} DeBruijnGraph;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NodeStack* createNodeStack(int capacity) {
    /*
    Create and return a stack for storing nodes.
    */
    NodeStack* theStack = (NodeStack*)(malloc(sizeof(NodeStack)));

    if (theStack == NULL) {
        fprintf(stderr, "Error. Could not allocate node stack with capacity %d\n", (capacity));
        exit(EXIT_FAILURE);
    }

    theStack->elements = (Node**)(malloc(sizeof(Node*)*capacity));

    if (theStack->elements == NULL) {
        fprintf(stderr, "Error. Could not allocate node stack elements with capacity %d\n", (capacity));
        exit(EXIT_FAILURE);
    }

    theStack->top = -1; // Empty;
    theStack->capacity = capacity;
    return theStack;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyNodeStack(NodeStack* theStack){
    /*
    Clears up memory in stack. Does not destroy the
    nodes.
    */
    free(theStack->elements);
    theStack->elements = NULL;
    free(theStack);
    theStack = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int NodeStack_IsEmpty(NodeStack* theStack){
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if (theStack->top == -1) {
        return 1;
    }
    else{
        return 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int NodeStack_IsFull(NodeStack* theStack){
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if (theStack->top + 1 == theStack->capacity) {
        return 1;
    }
    else{
        return 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeStack_Push(NodeStack* theStack, Node* element){
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Node** temp = NULL;

    // Need to realloc
    if (NodeStack_IsFull(theStack)) {
        temp = (Node**)(realloc(theStack->elements, 2 * sizeof(Node*) * theStack->capacity));

        if (temp == NULL) {
            fprintf(stderr, "Could not re-allocate NodeStack\n");
            exit(EXIT_FAILURE);
        }
        else{
            theStack->elements = temp;
            theStack->capacity *= 2;
        }
    }

    theStack->top += 1;
    theStack->elements[theStack->top] = element;
}
    

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Node* NodeStack_Pop(NodeStack* theStack){
    /*
    Pop and return the top element from the stack.
    */
    Node* theNode = NULL;

    if (NodeStack_IsEmpty(theStack)) {
        return NULL;
    }
    else{
        theNode = theStack->elements[theStack->top];
        theStack->top -= 1;
        return theNode;
    }
}

Path* createPath(int initialSize){
    /*
    Create and return a Path struct.
    */
    Path* thePath = (Path*)(malloc(sizeof(Path)));

    if (thePath == NULL) {
        fprintf(stderr, "Could not allocate path\n");
        exit(EXIT_FAILURE);
    }

    thePath->nodes = createNodeStack(initialSize);
    thePath->nNodes = 0;
    thePath->isBubble = 0;
    thePath->weight = 0.0;

    return thePath;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyPath(Path* thePath){
    /*
    free up memory in Path struct.
    */
    destroyNodeStack(thePath->nodes);
    thePath->nodes = NULL;
    free(thePath);
    thePath = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addNodeToPath(Path* thePath, Node* theNode, double weight){
    /*
    Add a Node to the specified path.
    */
    if (thePath == NULL) {
        fprintf(stderr, "Null path\n");
        exit(EXIT_FAILURE);
    }

    if (theNode == NULL) {
        fprintf(stderr, "Null Node\n");
        exit(EXIT_FAILURE);
    }

    NodeStack_Push(thePath->nodes, theNode);
    thePath->nNodes += 1;
    thePath->weight += weight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Node* createNode(char* sequence, int colour, int position, int kmerSize, double weight){
    /*
    Create and return a new node.
    */
    Node* theNode = (Node*)(malloc(sizeof(Node)));

    if (theNode == NULL) {
        fprintf(stderr, "Could not allocate node\n");
        exit(EXIT_FAILURE);
    }

    // Need to allocate kmerSize + 1 to store the null terminating character
    theNode->edges[0] = NULL;
    theNode->edges[1] = NULL;
    theNode->edges[2] = NULL;
    theNode->edges[3] = NULL;

    theNode->sequence = sequence;
    theNode->weight = weight;
    theNode->colours = colour;
    theNode->dfsColour = 'N';
    theNode->kmerSize = kmerSize;
    theNode->position = position;
    theNode->nEdges = 0;

    return theNode;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyEdge(Edge* theEdge) {
    /*
    free memory used by Edge struct.
    */
    free(theEdge);
    theEdge = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyNode(Node* theNode) {
    /*
    free all memory associated with a Node struct. The Node takes ownership
    of all associated string data.
    */
    // Edges belong to the node, so destroy these now.
    int i = 0;

    for (i = 0; i < theNode->nEdges; i++) {
        destroyEdge(theNode->edges[i]);
    }

    free(theNode);
    theNode = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int Node_Equal(Node* thisNode, Node* otherNode){
    /*
    Two nodes are equal if and only if their sequencese are equal.
    */
    if (thisNode->kmerSize != otherNode->kmerSize) {
        return 0;
    }
    else{
        return (strncmp(thisNode->sequence, otherNode->sequence, thisNode->kmerSize) == 0);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void Node_AddEdge(Node* theNode, Edge* theEdge){
    /*
    Add the specified edge to the specified node.
    */
    theNode->edges[theNode->nEdges] = theEdge;
    theNode->nEdges += 1;

    if (theNode->nEdges > 4) {
        fprintf(stderr, "Node struct cannot have > 4 edges. Now we have %d. This will cause problems", theNode->nEdges);
        exit(EXIT_FAILURE);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline char* Node_GetSuffix(Node* theNode){
    /*
    Return the suffix string of this node. The suffix is simply the last kmerSize-1
    elements of the node sequence.
    */
    return theNode->sequence + 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline char* Node_GetPrefix(Node* theNode){
    /*
    Return the suffix string of this node. The suffix is simply the first kmerSize-1
    elements of the node sequence.
    */
    return theNode->sequence;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline int nodePosComp(const void* x, const void* y){
    /*
    Comparison function for Node structs, for use in qsort, to sort Nodes by their
    positions.
    */
    const Node** nodeOne = (const Node**)(x);
    const Node** nodeTwo = (const Node**)(y);

    // Sort by position;
    return  nodeOne[0]->position - nodeTwo[0]->position;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline unsigned int hashKmer(char* theKmer, int size){
    /*
    Return a hash value for the specified kmer.
    */
    int i = 0;
    unsigned int hashVal = 0;

    // 1) Sum integers formed by 4-character sub-strings.
    for (i = 0; i < size - 4; i += 4) {
        hashVal += theKmer[i];
        hashVal += (theKmer[i + 1] << 8);
        hashVal += (theKmer[i + 2] << 16);
        hashVal += (theKmer[i + 3] << 24);
    }

    // 2) multiply by 101 and add new value. From Paul Larson (see StackOverflow)
    //for i in range(size):
    //    hashVal = hashVal * 101 + theKmer[i]

    return hashVal;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PathStack* createPathStack(int capacity){
    /*
    Create and return a stack for storing Paths.
    */
    PathStack* theStack = (PathStack*)(malloc(sizeof(PathStack)));

    if (theStack == NULL) {
        fprintf(stderr, "Error. Could not allocate path stack with capacity %d\n", capacity);
        exit(EXIT_FAILURE);
    }

    theStack->elements = (Path**)(malloc(sizeof(Path*)*capacity));

    if (theStack->elements == NULL) {
        fprintf(stderr, "Error. Could not allocate path stack elements with capacity %d\n", capacity);
        exit(EXIT_FAILURE);
    }
    
    theStack->top = -1; // Empty;
    theStack->capacity = capacity;
    return theStack;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyPathStack(PathStack* theStack){
    /*
    Clears up memory in stack. Does not destroy the
    Paths.
    */
    int i = 0;

    for (i = 0; i < (theStack->top + 1); i++) {
        destroyPath(theStack->elements[i]);
    }

    free(theStack->elements);
    theStack->elements = NULL;
    free(theStack);
    theStack = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int PathStack_IsEmpty(PathStack* theStack){
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if (theStack->top == -1) {
        return 1;
    }
    else{
        return 0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int PathStack_IsFull(PathStack* theStack){
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if (theStack->top + 1 == theStack->capacity) {
        return 1;
    }
    else{
        return 0;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PathStack_Push(PathStack* theStack, Path* element){
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Path** temp = NULL;

    // Need to realloc
    if (PathStack_IsFull(theStack)) {
        temp = (Path**)(realloc(theStack->elements, 2 * sizeof(Path*) * theStack->capacity));
    

        if (temp == NULL) {
            fprintf(stderr, "Could not re-allocate PathStack\n");
            exit(EXIT_FAILURE);
        }
        else{
            theStack->elements = temp;
            theStack->capacity *= 2;
        }
    }
    theStack->top += 1;
    theStack->elements[theStack->top] = element;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Path* PathStack_Pop(PathStack* theStack) {
    /*
    Pop and return the top element from the stack.
    */
    Path* thePath = NULL;

    if (PathStack_IsEmpty(theStack)) {
        return NULL;
    }
    else{
        thePath = theStack->elements[theStack->top];
        theStack->top -= 1;
        return thePath;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: Write function for printing nodes.
//    def __repr__(self):
//        /*
//        Returns a string representation of the node.
//        */
//        return "%s (%s). Pos = %s" %(self.sequence, self.colours, self.position)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Structs and functions for handling Edges.
// An Edge struct encapsulates an edge of the de-Bruijn graph.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Edge* createEdge(Node* startNode, Node* endNode, double weight){
    /*
    Create and return a new edge.
    */
    Edge* theEdge = (Edge*)(malloc(sizeof(Edge)));

    if (theEdge == NULL) {
        fprintf(stderr, "Error. Could not allocate edge\n");
        exit(EXIT_FAILURE);
    }

    theEdge->startNode = startNode;
    theEdge->endNode = endNode;
    theEdge->weight = weight;

    return theEdge;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EdgeStack* createEdgeStack(int capacity){
    /*
    Create and return a stack for storing Edges.
    */
    EdgeStack* theStack = (EdgeStack*)(malloc(sizeof(EdgeStack)));

    if (theStack == NULL) {
        fprintf(stderr, "Error. Could not allocate edge stack with capacity %d\n", (capacity));
        exit(EXIT_FAILURE);
    }

    theStack->elements = (Edge**)(malloc(sizeof(Edge*) * capacity));

    if (theStack->elements == NULL) {
        fprintf(stderr, "Error. Could not allocate edge stack elements with capacity %d\n", (capacity));
        exit(EXIT_FAILURE);
    }
        
    theStack->top = -1; // Empty
    theStack->capacity = capacity;
    return theStack;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyEdgeStack(EdgeStack* theStack){
    /*
    Clears up memory in stack. Does not destroy the
    Edges.
    */
    free(theStack->elements);
    theStack->elements = NULL;
    free(theStack);
    theStack = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int EdgeStack_IsEmpty(EdgeStack* theStack){
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if (theStack->top == -1) {
        return 1;
    }
    else{
        return 0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int EdgeStack_IsFull(EdgeStack* theStack){
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if (theStack->top == theStack->capacity - 1) {
        return 1;
    }
    else {
        return 0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EdgeStack_Push(EdgeStack* theStack, Edge* element){
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Edge** temp = NULL;

    // Need to realloc
    if (EdgeStack_IsFull(theStack)) {
        temp = (Edge**)(realloc(theStack->elements, 2 * sizeof(Edge*) * theStack->capacity));

        if (temp == NULL) {
            fprintf(stderr, "Could not re-allocate EdgeStack\n");
            exit(EXIT_FAILURE);
        }
        else {
            theStack->elements = temp;
            theStack->capacity *= 2;
        }
    }

    theStack->top += 1;
    theStack->elements[theStack->top] = element;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Edge* EdgeStack_Pop(EdgeStack* theStack){
    /*
    Pop and return the top element from the stack.
    */
    Edge* theEdge = NULL;

    if (EdgeStack_IsEmpty(theStack)) {
        return NULL;
    }
    else {
        theEdge = theStack->elements[theStack->top];
        theStack->top -= 1;
        return theEdge;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NodeDict* createNodeDict(int nBuckets){
    /*
    Create and return a dictionary of kmer/void pointers.
    */
    NodeDict* theDict = (NodeDict*)(malloc(sizeof(NodeDict)));

    if (theDict == NULL) {
        fprintf(stderr, "Error. Could not NodeDict\n");
        exit(EXIT_FAILURE);
    }

    theDict->buckets = (Node***)(malloc(nBuckets * sizeof(Node**)));

    if (theDict->buckets == NULL) {
        fprintf(stderr, "Error. Could not NodeDict. buckets of size %d\n", (nBuckets));
        exit(EXIT_FAILURE);
    }

    theDict->bucketSize = (int*)(malloc(nBuckets*sizeof(int)));

    if (theDict->bucketSize == NULL) {
        fprintf(stderr, "Error. Could not NodeDict.bucketSize of size %d\n", (nBuckets));
        exit(EXIT_FAILURE);
    }
    
    theDict->nBuckets = nBuckets;

    int i = 0;

    for (i = 0; i < (nBuckets); i++) {
        theDict->buckets[i] = NULL;
        theDict->bucketSize[i] = 0;
    }

    return theDict;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyNodeDict(NodeDict* theDict){
    /*
    free memory used by NodeDict.
    */
    int i = 0;

    for (i = 0; i < (theDict->nBuckets); i++) {
        free(theDict->buckets[i]);
    }

    free(theDict->buckets);
    free(theDict->bucketSize);
    free(theDict);
    theDict = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int NodeDict_FindOrInsert(NodeDict* theDict, Node** theNode, int keyLen, Node** nodeForUpdating){
    /*
    Either return the element which is associated with the specified key, or
    insert a new element at the relevant position. Return 0 if the element was not found,
    and 1 if it was.
    */
    int hashValue = hashKmer(theNode[0]->sequence, keyLen) % theDict->nBuckets;
    int bucketSize = 0;
    int i = 0;
    int initialBucketSize = 5;
    Node* testNode = NULL;
    Node* newNode = NULL;

    // Need to allocate new bucket
    if (theDict->buckets[hashValue] == NULL) {
        theDict->buckets[hashValue] = (Node**)(malloc(initialBucketSize * sizeof(Node*)));

        if (theDict->buckets[hashValue] == NULL) {
            fprintf(stderr, "Could not allocate hash table bucket with size %d\n", (initialBucketSize));
            exit(EXIT_FAILURE);
        }

        // Always set to NULL first
        for (i = 0; i < (initialBucketSize); i++) {
            theDict->buckets[hashValue][i] = NULL;
        }

        newNode = createNode(theNode[0]->sequence, theNode[0]->colours, theNode[0]->position, theNode[0]->kmerSize, theNode[0]->weight);
        theDict->buckets[hashValue][0] = newNode;
        theDict->bucketSize[hashValue] = initialBucketSize;
        theNode[0] = newNode;
        return 0;
    }

    // Bucket is there. Check all elements in bucket for match.
    else{
        bucketSize = theDict->bucketSize[hashValue];

        for (i = 0; i < (bucketSize); i++) {

            // Found empty slot. Insert new element
            if (theDict->buckets[hashValue][i] == NULL) {
                newNode = createNode(theNode[0]->sequence, theNode[0]->colours, theNode[0]->position, theNode[0]->kmerSize, theNode[0]->weight);
                theDict->buckets[hashValue][i] = newNode;
                theNode[0] = newNode;
                return 0;
            }

            // Check for match
            else{
                // Match. Return this element.
                testNode = theDict->buckets[hashValue][i];
                if ((theNode[0] == testNode) || strncmp(theNode[0]->sequence, testNode->sequence, keyLen) == 0) {
                    nodeForUpdating[0] = testNode;
                    return 1;
                }
            }
        }
    }

    // If we get here, then we found no empty slots, and no matches. So we need to
    // allocate more space in the relevant bucket and then insert the key/value pair
    // in the next available space..
    int oldBucketSize = theDict->bucketSize[hashValue];
    int newBucketSize = 2*oldBucketSize;
    Node** temp = (Node**)(realloc(theDict->buckets[hashValue], sizeof(Node*) * newBucketSize));

    if (temp == NULL) {
        fprintf(stderr, "Could not re-allocate bucket\n");
        exit(EXIT_FAILURE);
    }
    else{
        // Set new entries to NULL
        for (i = oldBucketSize; i < newBucketSize; i++) {
            temp[i] = NULL;
        }
        
        newNode = createNode(theNode[0]->sequence, theNode[0]->colours, theNode[0]->position, theNode[0]->kmerSize, theNode[0]->weight);
        theDict->bucketSize[hashValue] = newBucketSize;
        theDict->buckets[hashValue] = temp;
        theDict->buckets[hashValue][oldBucketSize] = newNode;
        theNode[0] = newNode;
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DeBruijnGraph* createDeBruijnGraph(int kmerSize, int nBuckets){
    /*
    Allocate memory for graph data.
    */
    DeBruijnGraph* theGraph = (DeBruijnGraph*)(malloc(sizeof(DeBruijnGraph)));

    if (theGraph == NULL) {
        fprintf(stderr, "Could not allocate memory for DeBruijnGraph\n");
        exit(EXIT_FAILURE);
    }
    theGraph->kmerSize = kmerSize;
    theGraph->allNodes = createNodeStack(nBuckets);
    theGraph->nodes = createNodeDict(nBuckets);

    return theGraph;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyDeBruijnGraph(DeBruijnGraph* theGraph) {
    /*
    Free memory used by graph.
    */
    int i = 0;

    // Destroy all nodes
    for (i = 0; i < (theGraph->allNodes->top + 1); i++) {
        destroyNode(theGraph->allNodes->elements[i]);
    }
    // These only hold pointers to memory which will be
    // freed elsewhere.
    destroyNodeStack(theGraph->allNodes);
    destroyNodeDict(theGraph->nodes);
    free(theGraph);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int DeBruijnGraph_InsertOrUpdateNode(DeBruijnGraph* theGraph, Node** theNode){
    /*
    Check if a node is already present. If it is, update it, otherwise
    insert it. Return 1 if the node was inserted and 0 if it was updated.
    */
    Node* nodeForUpdating = NULL;
    int foundNode = NodeDict_FindOrInsert(theGraph->nodes, theNode, theGraph->kmerSize, &nodeForUpdating);

    // Need to create a new node, copying values from theNode
    if (!foundNode) {
        NodeStack_Push(theGraph->allNodes, theNode[0]);
        return 1;
    }
    // Update existing node
    else{
        // Update colours of nodes already in graph.
        nodeForUpdating->colours |= theNode[0]->colours;
        nodeForUpdating->weight += theNode[0]->weight;
        theNode[0] = nodeForUpdating;
        return 0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DeBruijnGraph_AddEdge(DeBruijnGraph* theGraph, Node* startNode, Node* endNode, double weight) {
    /*
    */
    DeBruijnGraph_InsertOrUpdateNode(theGraph, &startNode);
    DeBruijnGraph_InsertOrUpdateNode(theGraph, &endNode);
    Edge* newEdge = NULL;

    int i = 0;

    // Check all outgoing edges from startNode. If it has none, then make one for this edge. Otherwise,
    // check all existing edges for a match, and update accordingly.
    for (i = 0; i < 4; i++) {
        if (startNode->edges[i] == NULL) {
            newEdge = createEdge(startNode, endNode, weight);
            Node_AddEdge(startNode, newEdge);
            break;
        }
        else if (startNode->edges[i]->endNode == endNode) {
            startNode->edges[i]->weight += weight;
            break;
        }
        else{
            continue;
        }
    }
        
    // else {
    //     pass; // This only happens when there are Ns in the sequence.
        //logger.error("Error in assembler. Could not find matching end-node for edge. Something is very wrong")
        //logger.error("Start node sequence is %s. End node sequence is %s" %(startNode[0].sequence[0:startNode[0].kmerSize], endNode[0].sequence[0:endNode[0].kmerSize]))
        //raise StandardError, "Assembly Error!!"
    // }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dfsVisit(Node* theNode, double minWeight){
    /*
    */
    int nEdges = theNode->nEdges;
    int i = 0;
    Node* nextNode = NULL;
    Edge* edge = NULL;

    theNode->dfsColour = 'g';

    for (i = 0; i < nEdges; i++) {
        edge = theNode->edges[i];

        // Ignore low-weight edges that are only in reads
        if (edge->endNode->colours == READ && edge->weight < minWeight) {
            continue;
        }
        nextNode = edge->endNode;

        if (nextNode->dfsColour == 'w') {

            // Found cycle in this path
            if (dfsVisit(nextNode, minWeight) == 1) {
                return 1;
            }
            // This path ok. Go to next edge
            else{
                continue;
            }
        }
        else if (nextNode->dfsColour == 'g') {
            // Found cycle
            //logger.debug("Found cycle. From %s to %s" %(theNode.position, nextNode.position))
            return 1;
        }

        // Black node. Already explored past this node. Go to next edge.
        else{
            continue;
        }
    }
    
    // No cycles in any path reachable from this node.
    theNode->dfsColour = 'b';
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int detectCyclesInGraph_Recursive(DeBruijnGraph* theGraph, double minWeight){
    /*
    */
    Node** allNodes = theGraph->allNodes->elements;
    Node* thisNode = NULL;
    int i = 0;
    int nNodes = theGraph->allNodes->top + 1;
    int foundCycle = 0;

    for (i = 0; i < nNodes; i++) {
        thisNode = allNodes[i];
        thisNode->dfsColour = 'w';
    }

    for (i = 0; i < nNodes; i++) {
        thisNode = allNodes[i];

        if (thisNode->dfsColour == 'w') {
            // Found cycle
            foundCycle = dfsVisit(thisNode, minWeight);

            if (foundCycle == 1) {
                return 1;
            }
        }
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int detectCyclesInGraph(DeBruijnGraph* theGraph, double minWeight){
    /*
    */
    Node* thisNode;
    Node* nextNode;
    Edge* edge;
    int i = 0;
    int j = 0;
    int nEdges = 0;
    int nNodes = theGraph->allNodes->top + 1;
    Node** allNodes = theGraph->allNodes->elements;

    int nFilledBuckets = 0;
    int nEntriesThisBucket = 0;

    for (i = 0; i < theGraph->nodes->nBuckets; i++){ 
        if (theGraph->nodes->buckets[i] != NULL) {
            nFilledBuckets += 1;
            for (j = 0; j < theGraph->nodes->bucketSize[i]; j++) {
                if (theGraph->nodes->buckets[i][j] != NULL) {
                    nEntriesThisBucket += 1;
                }
                else{
                    break;
                }
            }
        }
    }

#ifdef DEBUG
    fprintf(stderr, "nNodes = %s. nFilledBuckets = %s. mean entries/bucket = %s\n", (nNodes, nFilledBuckets, (float)(nEntriesThisBucket)/nFilledBuckets));
#endif
    qsort((void*)allNodes, nNodes, sizeof(Node*), nodePosComp);

    for (i = 0; i < nNodes; i++ ){
        thisNode = allNodes[i];
        thisNode->dfsColour = 'w';
    }

    Node* sourceNode = allNodes[0];
    Node* endNode = allNodes[nNodes-1];
    NodeStack* theStack = createNodeStack(nNodes);
    int reachedEnd = 0;

    NodeStack_Push(theStack, sourceNode);

    while (!NodeStack_IsEmpty(theStack)) {

        thisNode = NodeStack_Pop(theStack);

        if (thisNode->dfsColour == 'w') {
            thisNode->dfsColour = 'g';
        }
        else if (thisNode->dfsColour == 'g') {
            thisNode->dfsColour = 'b';
        }
        // else {
        // pass;
        // }
        
        nEdges = thisNode->nEdges;

        for (i = 0; i < nEdges; i++) {
            edge = thisNode->edges[i];
            nextNode = edge->endNode;

            // TODO{ temp hack. Replace with Nodes_Equal later
            if (Node_Equal(nextNode, endNode)) {
                reachedEnd = 1;
            }

            if (nextNode->dfsColour == 'w') {
                NodeStack_Push(theStack, nextNode);
            }
            // Found a cycle
            else if (nextNode->dfsColour == 'g') {
                destroyNodeStack(theStack);
                return 1;
            }
            // else{
            //     pass;
            // }
   
        }
        
    }
    
    destroyNodeStack(theStack);
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char* createSequenceFromPath(Path* thePath){
    /*
    Create and return a string representation of the sequence of a specific path
    through the graph.
    */
    int nNodes = thePath->nNodes;
    char* theString = (char*)(malloc( (nNodes+1)*sizeof(char)));

    if (theString == NULL) {
        fprintf(stderr, "Could not allocate memory for string of size %d\n", nNodes + 1);
        exit(EXIT_FAILURE);
    }
    int i = 0;

    for (i = 0; i < nNodes; i++) {
        theString[i] = thePath->nodes->elements[i]->sequence[0];
    }
    theString[nNodes] = 0;
    return theString;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int checkPathForCycles(Path* thePath){
    /*
    Check if this path contains a cycle.
    */
    int nNodes = thePath->nNodes;
    int i = 0;

    //logger.debug("Checking path with %s nodes for cycles" %(nNodes))

    // Set all dfs colours to white
    for (i = 0; i < nNodes; i++) {
        thePath->nodes->elements[i]->dfsColour = 'w';
    }
    // Check all nodes in order. If we see the same node twice, then
    // there is a cycle. If we get to the end without seeing any nodes twice,
    // then no cycle.
    for (i = 0; i < nNodes; i++) {
        if (thePath->nodes->elements[i]->dfsColour == 'w') {
            thePath->nodes->elements[i]->dfsColour = 'g';
        }
        else{
            //logger.debug("Found cycle")
            return 1;
        }
    }

    //logger.debug("No cycles")
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PathStack* getVariantPathsThroughGraphFromNode(DeBruijnGraph* theGraph, Path* thePath, double minWeight){
    /*
    Check all valid paths through the graph starting at the last node in "thePath". If any path
    returns to the reference sequence, then stop and add that path to the returned list. Also, if the
    path never returns to the reference, but is sufficiently long, then add this to the list.
    */
    PathStack* thePathStack = createPathStack(10);
    PathStack* finishedPaths = createPathStack(10);
    Path* pathSoFar = NULL;
    Path* newPath = NULL;
    Node* endSoFar = NULL;
    Node* newEnd = NULL;
    Edge* theEdge = NULL;
    int nEdgesThisNode = 0;
    int i = 0;
    int j = 0;
    int hasCycle = 0;

    PathStack_Push(thePathStack, thePath);

    while (!PathStack_IsEmpty(thePathStack)) {

        pathSoFar = PathStack_Pop(thePathStack);
        endSoFar = pathSoFar->nodes->elements[pathSoFar->nNodes-1];

        // TODO{ Replace with maxHaplotypes??
        if (thePathStack->top + 1 > 20 || finishedPaths->top + 1 > 20) {
            //logger.info("Too many paths %s (%s) from this node. Giving up" %(thePathStack.top, finishedPaths.top))
            destroyPath(pathSoFar);
            destroyPathStack(thePathStack);
            destroyPathStack(finishedPaths);
            return NULL;
        }
        // At the moment, don't allow any cycles.
        hasCycle = checkPathForCycles(pathSoFar);

        if (hasCycle) {
            destroyPath(pathSoFar);
        }
        // Got back to ref, this path is done with. This is a bubble.
        else if (endSoFar->colours == REF_AND_READ) {
            pathSoFar->isBubble = 1;
            PathStack_Push(finishedPaths, pathSoFar);
        }
        // No reads here. Not quite sure how this could happen. Went from ref and read to only
        // ref.
        else if (endSoFar->colours == REF) {
            destroyPath(pathSoFar);
        }
        // Keep extending path
        else{
            // Dumb check for loops
            //if pathSoFar.nNodes > 1000:

                // If we are at the end of this path, and we haven't returned to the reference yet, then still suggest a variant
                // using all the sequence so far, but only if the path is long enough. This is not a bubble.
                //if nEdgesThisNode == 0:
                //    destroyPath(pathSoFar)
                //    if pathSoFar.nNodes > 50:
                //        pathSoFar.isBubble = 0
                //        PathStack_Push(finishedPaths, pathSoFar)
                //    else:
                //        destroyPath(pathSoFar)

                //else:
            nEdgesThisNode = endSoFar->nEdges;

            for (i = 0; i < nEdgesThisNode; i++) {
                theEdge = endSoFar->edges[i];
                newEnd = theEdge->endNode;
            
                if (theEdge->weight >= minWeight || newEnd->colours == REF_AND_READ || newEnd->colours == REF) {
                    newPath = createPath(theGraph->kmerSize);

                    // Copy old path
                    for (j = 0; j < pathSoFar->nNodes; j++) {
                        addNodeToPath(newPath, pathSoFar->nodes->elements[j], 0.0);
                    }
                    // Weight for this path is weight of existing path + weight of new edge
                    newPath->weight = pathSoFar->weight;
                    addNodeToPath(newPath, newEnd, theEdge->weight);
                    PathStack_Push(thePathStack, newPath);
                }
                
            }
            destroyPath(pathSoFar);
        }
        
    }    
        
    destroyPathStack(thePathStack);
    return finishedPaths;
}  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void logPath(Path* thePath, char* refSeq, int refStart){
    /*
    Log a path through the graph.
    */
    int i = 0;
    int startPos = thePath->nodes->elements[0]->position;
    Node* theNode;
#ifdef DEBUG
    fprintf(stderr, "Logging path of %s nodes\n", (thePath->nNodes));
#endif
    for (i = 0; i < thePath->nNodes; i++) {
        theNode = thePath->nodes->elements[i];
#ifdef DEBUG
        fprintf(stderr, "Pos = %s. Seq = %s. RefSeq = %c. Colours = %s. Node weight = %s\n", (theNode.position, theNode.sequence[0: theNode.kmerSize], refSeq[startPos + i - refStart], theNode.colours, theNode.weight));
#endif
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadReferenceIntoGraph(DeBruijnGraph* theGraph, char* refSeq, int refStart, int kmerSize){
    /*
    Load k-mers from the specified reference sequence into the
    graph.
    */
    int i = 0;
    int lenRef = strlen(refSeq);
    Node tempStartNode;
    Node tempEndNode;

    for (i = 0; i < (lenRef-kmerSize) - 1; i++) {

        tempStartNode.sequence = refSeq + i;
        tempStartNode.kmerSize = kmerSize;
        tempStartNode.colours = REF;
        tempStartNode.position = refStart + i;
        tempStartNode.weight = 1;

        tempEndNode.sequence = refSeq + i + 1;
        tempEndNode.kmerSize = kmerSize;
        tempEndNode.colours = REF;
        tempEndNode.position = refStart + i + 1;
        tempEndNode.weight = 1;

        DeBruijnGraph_AddEdge(theGraph, &tempStartNode, &tempEndNode, 1);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char* createReverseComplementSequence(char* theSeq, int seqLen){
    /*
    Create and return a string containing the reverse complement of the specified
    sequence.
    */
    char* newSeq = (char*)(malloc(sizeof(char)*(seqLen+1)));
    int i = 0;

    for (i = 0; i < seqLen; i++) {
        if (theSeq[i] == 'A') {
            newSeq[seqLen - i - 1] = 'T';
        }
        else if (theSeq[i] == 'T') {
            newSeq[seqLen - i - 1] = 'A';
        }
        else if (theSeq[i] == 'C') {
            newSeq[seqLen - i - 1] = 'G';
        }
        else if (theSeq[i] == 'G') {
            newSeq[seqLen - i - 1] = 'C';
        }
        else {
            newSeq[seqLen - i - 1] = theSeq[i];
        }
    }
    newSeq[seqLen] = 0;
    return newSeq;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadReadIntoGraph(alignedRead* theRead, DeBruijnGraph* theGraph, int minQual, int kmerSize) {
    /* 
    */
    char* theSeq = theRead->seq;
    uint8_t* theQuals = theRead->qual;
    int length = theRead->rlen;
    int startPos = theRead->pos;
    int i = 0;
    int j = 0;
    int thisMinQual = 100000000;
    int NsInKmer = 0;
    Node tempStartNode;
    Node tempEndNode;
    // fprintf(stderr, "Name:%s\n", theRead->qname);

    for (i = 0; i < (length-kmerSize) - 1; i++) {
        thisMinQual = 100000000;
        NsInKmer = 0;

        for (j = i; j < i + kmerSize + 1; j++) {

            thisMinQual = std::min(thisMinQual, (int)theQuals[j]);

            if (theSeq[j] == 'N') {
                NsInKmer = 1;
            }
        }

        if (thisMinQual >= minQual && NsInKmer == 0) {

            tempStartNode.sequence = theSeq + i;
            tempStartNode.kmerSize = kmerSize;
            tempStartNode.colours = READ;
            tempStartNode.position = - 1;
            tempStartNode.weight = thisMinQual;

            tempEndNode.sequence = theSeq + i + 1;
            tempEndNode.kmerSize = kmerSize;
            tempEndNode.colours = READ;
            tempEndNode.position = -1;
            tempEndNode.weight = thisMinQual;

            DeBruijnGraph_AddEdge(theGraph, &tempStartNode, &tempEndNode, thisMinQual);
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadBAMDataIntoGraph(DeBruijnGraph* theGraph, struct alignedRead* start, struct alignedRead* end, int assembleBadReads, int assembleBrokenPairs, int minQual, int kmerSize){
    /*
    Load k-mers from the specified BAM file into the graph. K-mers containing
    Ns are ignored, as are k-mers containing low-quality bases.
    */
    alignedRead* readStart = start;
    alignedRead* readEnd = end;

    while (readStart != readEnd) {
        if (!Read_IsQCFail(readStart)) {
            loadReadIntoGraph(readStart, theGraph, minQual, kmerSize);
        }
        readStart += 1;
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void assembleReadsAndDetectVariants(int refStart, int refEnd, struct alignedRead* windowStart, struct alignedRead* windowEnd, char* refSeq) {
    /*
    Below are filled from default Platypus options
    */
    int minQual = 20;
    int minMapQual = 20;
    int largestVariant = 1500;
    int minReads = 2;
    int assembleBadReads = 1;
    int assembleBrokenPairs = 0;
    int kmerSize = 15;
    int minWeight = minReads*minQual;
    int nBuckets = 5000;

    DeBruijnGraph* theGraph = createDeBruijnGraph(kmerSize, nBuckets);

    loadReferenceIntoGraph(theGraph, refSeq, refStart, kmerSize);
    loadBAMDataIntoGraph(theGraph, windowStart, windowEnd, assembleBadReads, assembleBrokenPairs, minQual, kmerSize);

    // If this is true, then don't allow cycles in the graph.
    /*
    while (detectCyclesInGraph_Recursive(theGraph, minWeight)) {
        if (kmerSize > 50) {
            if (verbosity >= 3) {
                fprintf(stderr, "Could not assemble region %s:%d-%d without cycles. Max k-mer size tried = %d\n", chrom, assemStart, assemEnd, kmerSize);
            }
            break;
        }
        else{
            if (verbosity >= 3) {
                fprintf(stderr, "Found cycles in region %s:%d-%d with kmer size %d. Trying again with kmer size %d\n", chrom, assemStart, assemEnd, kmerSize, kmerSize+5);
            }
            kmerSize += 5;
            destroyDeBruijnGraph(theGraph);
            theGraph = createDeBruijnGraph(kmerSize, nBuckets);
            loadReferenceIntoGraph(theGraph, refSeq, refStart, kmerSize);
            loadBAMDataIntoGraph(theGraph, readBuffers, assembleBadReads, assembleBrokenPairs, minQual, kmerSize);
        }
    }
    */
    destroyDeBruijnGraph(theGraph);
    //logger.debug("Finished assembling region %s:%s-%s" %(chrom, start, end))
    //logger.debug("Found vars %s" %(theVars))
    // return sorted(theVars);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc,char** argv){
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    // check args
    if (argc != 5) {
        fprintf(stderr, "Usage %s file.bam chr:start-stop ref.fa n_threads\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    // these come from htslib/sam.h
	hts_itr_t *iter = NULL;
	hts_idx_t *idx = NULL;
	samFile *in = NULL;
	bam1_t *b = NULL;
    bam_hdr_t *header = NULL;

    // open the BAM file for reading (though called sam_open it opens bam files too :P)
    in = sam_open(argv[1], "r");
    errorCheckNULL(in);

    //get the sam header. 
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }
    
    // print the chromosome names in the header
    // see the bam_hdr_t struct in htslib/sam.h for parsing the rest of the stuff in header
    // int i;
    // for (i = 0; i < (header->n_targets); i++) {
    //     printf("Chromosome ID %d = %s\n", i, (header->target_name[i]));
    // } 
    
    faidx_t* fai = fai_load(argv[3]);

    // load the index file for BAM
	idx = sam_index_load(in, argv[1]);
	errorCheckNULL(idx);
    
    // the iterator for the BAM random access is probably initialised here. Note that we pass the region string to this function 
	iter  = sam_itr_querys(idx, header, argv[2]); 
	errorCheckNULL(iter);
    
    // this must be the initialisation for the structure that stores a read (need to verify)
	b = bam_init1();
    
    int numThreads = atoi(argv[4]);

    // my structure for a readbuffer (see common.h)
    struct bamReadBuffer readBuffer;

    for (int i = 0; i < numThreads; i++) {
        readBuffer.reads.__capacity = MAX_READS_IN_REGION;
        readBuffer.reads.__size = 0;
        readBuffer.reads.array = (struct alignedRead*) malloc(sizeof(struct alignedRead) * MAX_READS_IN_REGION);
        readBuffer.reads.windowStart = readBuffer.reads.array;
        readBuffer.reads.windowEnd = readBuffer.reads.array;
        readBuffer.reads.__longestRead = 0;
    }
    
    char* reg = argv[2];
    int beg, end;
    const char *q;
    char tmp_a[1024], *tmp = tmp_a;
    q = hts_parse_reg(reg, &beg, &end);
    if (q) {   
        if (q - reg + 1 > 1024)
            tmp = (char*)malloc(q - reg + 1);
        strncpy(tmp, reg, q - reg);
        tmp[q - reg] = 0;
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        beg = 0; end = INT_MAX;
    }

    struct timeval start_time, end_time;
    double runtime = 0;

    std::vector<Batch> batches;

    // extract reads from region    
    while (sam_itr_next(in, iter, b) >= 0) {
        getRead(readBuffer.reads.windowEnd, b); // copy the current read to the myread structure. See common.c for information
        // printRead(readBuffer.reads.windowEnd, header);  // print data in structure. See common.c for information;
        int readLength = readBuffer.reads.windowEnd->end - readBuffer.reads.windowEnd->pos;
        if (readLength > readBuffer.reads.__longestRead) {
            readBuffer.reads.__longestRead = readLength;
        }
        readBuffer.reads.__size++;
        readBuffer.reads.windowEnd++;
            
        if (readBuffer.reads.__size == readBuffer.reads.__capacity) {
            readBuffer.reads.__capacity *= 2;
            readBuffer.reads.array = (struct alignedRead*) realloc(readBuffer.reads.array, sizeof(struct alignedRead) * readBuffer.reads.__capacity);
            readBuffer.reads.windowStart = readBuffer.reads.array;
            readBuffer.reads.windowEnd = readBuffer.reads.windowStart + readBuffer.reads.__size;
            fprintf(stderr,"Buffer doubled. Old size = %d, New size = %d\n", readBuffer.reads.__size, readBuffer.reads.__capacity);
        }
    }

    // process reads
    const int assemblyRegionSize = 1500;
    int assemRegionShift = std::max(100, std::min(1000, assemblyRegionSize / 2));
    for (int k = beg; k < end; k += assemRegionShift) {
        int assemStart = k;
        int assemEnd = std::min(assemStart + assemblyRegionSize, end);
        int refStart = std::max(0, assemStart - assemblyRegionSize);
        int refEnd = assemEnd + assemblyRegionSize;
        int len;
        char* ref = faidx_fetch_seq(fai, tmp, refStart, refEnd - 1, &len);
        setWindowPointers(&readBuffer.reads, assemStart, assemEnd);
        Batch b;
        b.offset = k;
        b.ref = ref;
        b.windowStart = readBuffer.reads.windowStart;
        b.windowEnd = readBuffer.reads.windowEnd;
        batches.push_back(b);
    }

#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    if (tid == 0) {
        fprintf(stderr, "Found %d batches. Running with threads: %d\n", batches.size(), numThreads);
    }
}

    gettimeofday(&start_time, NULL);
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < batches.size(); i++) {
            int assemStart = batches[i].offset;
            int assemEnd = std::min(assemStart + assemblyRegionSize, end);
            int refStart = std::max(0, assemStart - assemblyRegionSize);
            int refEnd = assemEnd + assemblyRegionSize;
            int verbosity = 2;
            if (verbosity >= 3) {
                fprintf(stderr, "Assembling region %s:%d-%d, tid = %d\n", tmp, assemStart, assemEnd, tid);
            }
            assembleReadsAndDetectVariants(refStart, refEnd, batches[i].windowStart, batches[i].windowEnd, batches[i].ref);
        }
}
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    gettimeofday(&end_time, NULL);
    runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;

    // wrap up
    for (int i = 0; i < batches.size(); i++) {
        free(batches[i].ref);
    }
    free(readBuffer.reads.array);

    if (tmp != tmp_a) {
        free(tmp);
    }
	hts_idx_destroy(idx);
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
    fai_destroy(fai);

    fprintf(stderr, "Kernel runtime: %.2f s\n", runtime*1e-6);
    
    return 0;
}
