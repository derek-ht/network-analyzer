// Dijkstra API implementation
// COMP2521 Assignment 2

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

static ShortestPaths *initArray(Graph g);
//static PQ initPq(Graph g);

ShortestPaths dijkstra(Graph g, Vertex src) {
	// Initialize all arrays
	ShortestPaths *sps = initArray(g);
	
	sps->src = src;
	
	sps->numNodes = GraphNumVertices(g);
	
	// Initialize priority queue of vertices adjacten to source

	PQ vSet = PQNew();
	
	sps->dist[src] = 0;

	PQInsert(vSet, src, 0);
	
	while (!PQIsEmpty(vSet)) {
		Vertex ver = PQDequeue(vSet);
				
		AdjList curr = GraphOutIncident(g, ver);
		
		while (curr != NULL) { 
			// EDGE RELAXING
			// If distance between the currrent 'ver' and its neighbour 'curr->v',
			// plus the distance from the source to the 'ver' is less than the
			// current distance from the source to the neighbouring vertex. Then
			// reinitialize the dist[] and pred[] values of 'curr->v'
			if ((curr->weight + sps->dist[ver]) < (sps->dist[curr->v])) {
				// Replace dist[] index with n
				sps->dist[curr->v] = curr->weight + sps->dist[ver];

				PredNode *currPred = sps->pred[curr->v];
				
				// Clear old predecessors
				while (currPred != NULL) {

					PredNode *tmp = currPred;

					currPred = currPred->next;

					free(tmp);
				}

				// Initialize new pred node, adding to the head and adding the
				// current 'ver' to its neighbour's predecessor array.
				sps->pred[curr->v] = malloc(sizeof(struct PredNode));

				sps->pred[curr->v]->next = NULL;
				
				sps->pred[curr->v]->v = ver;
				
				// Add the vertex to the queue
				PQInsert(vSet, curr->v, sps->dist[curr->v]);
			
			// If the current 'ver' has the same distance to its neighbour as the
			// currently saved distance in the array then add it as another predecessor
			} else if ((curr->weight + sps->dist[ver]) == sps->dist[curr->v]) {
				
				// Initialize PredNode and add it to the head of the list
				PredNode *temp = malloc(sizeof(struct PredNode));

                temp->v = ver;

                temp->next = sps->pred[curr->v];
				
                sps->pred[curr->v] = temp;
			}

			curr = curr->next;
		}
	}
	PQFree(vSet);

	return *sps;

}

// Initialize ShortestPath's arrays
static ShortestPaths *initArray(Graph g) {
	int numVertices = GraphNumVertices(g);

	ShortestPaths *newPath = malloc(sizeof(struct ShortestPaths));

	newPath->dist = malloc(numVertices * sizeof(int));

	newPath->pred = malloc(numVertices * sizeof(struct PredNode));
	
	for (int i = 0; i < numVertices; i++) {
	
		newPath->dist[i] = INFINITY;
	
		newPath->pred[i] = NULL;
	
	}
	
	return newPath;
}

void showShortestPaths(ShortestPaths sps) {

}

// Free the ShortestPaths struct
void freeShortestPaths(ShortestPaths sps) {
    free(sps.dist);
    for (int i = 0; i < sps.numNodes; i++) {   
        PredNode *curr = sps.pred[i];
        while (curr != NULL) {
            PredNode *tmp = curr;
            curr = curr->next;
            free(tmp);
        }
    }
    free(sps.pred);
}
