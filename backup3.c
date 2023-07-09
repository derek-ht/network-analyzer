// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX
#define X 0
#define Y 1

// static bool incomplete(Dendrogram *dendA, int numV);
static Dendrogram createDnode(int v);
static void mergeDend(double **dist, Dendrogram *dendA, int numV, int i[2], 
					  int method);
static void findMin(double **dist, int numV, int index[2]);
static double lanceFormula(double ci, double cj, int method);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */

Dendrogram LanceWilliamsHAC(Graph g, int method) {
	int numVerts = GraphNumVertices(g);
	double **dist = malloc(numVerts * sizeof(*dist));
	
	for (int i = 0; i < numVerts; i++) {
		dist[i] = malloc(numVerts * sizeof(double));
		for (int j = 0; j < numVerts; j++) {
			dist[i][j] = INFINITY;
		}
	}

	for (int v = 0; v < numVerts; v++) {
		AdjList edgeOut = GraphOutIncident(g, v);
		while (edgeOut != NULL) {
			double distOut = (1.0 / edgeOut->weight);
			if (distOut < dist[v][edgeOut->v]) {
				dist[v][edgeOut->v] = distOut;
			}
			edgeOut = edgeOut->next;
		}
		AdjList edgeIn = GraphOutIncident(g, v);
		while (edgeIn != NULL) {
			double distIn = (1.0 / edgeIn->weight);
			if (distIn < dist[v][edgeIn->v]) {
				dist[v][edgeIn->v] = distIn;
			}
			edgeIn = edgeIn->next;
		}
	}
	
	Dendrogram *dendA = malloc(numVerts * sizeof(struct DNode));
	for (int i = 0; i < numVerts; i++) {
		dendA[i] = createDnode(i);
	}
	
	for (int i = 1; i < numVerts; i++) {
		while (dendA[i] != NULL) {
			int index[2];
			findMin(dist, numVerts, index);
			mergeDend(dist, dendA, numVerts, index, method);
		}
	}

	for (int i = 0; i < numVerts; i++) {
		free(dist[i]);
	}
	free(dist);
	Dendrogram root = dendA[0];
	free(dendA);
	// Return root Dendrogram
	return root;
}

// static int incomplete(Dendrogram *dendA, int numV) {

// 	// If there are dendrograms found in the array (besides index 0 of array),
// 	// return true.
// 	for (int i = 1; i < numV; i++) {
// 		if (dendA[i] != NULL) {
// 		return TRUE;
// 		}
// 	}

// 	return false;
// }

static Dendrogram createDnode(int v) {
	Dendrogram newNode = malloc(sizeof (struct DNode));
	newNode->vertex = v;
	newNode->left = NULL;
	newNode->right = NULL;
	return newNode;
}

// Merge 2 dendrograms.
static void mergeDend(double **dist, Dendrogram *dendA, int numV, int *i, 
					  int method) {
	// If x is greater than y,
	// swap x and y.

	int x, y;

	if (i[Y] < i[X]) {
		x = i[Y];
		y = i[X];
	} else {
		x = i[X];
		y = i[Y];
	}

	// Update distances for all vertices adjacent to x and y
	for (int i = 0; i < numV; i++) {
		if ((i != x) && (i != y) && (dendA[i] != NULL)) {
			dist[i][x] = lanceFormula(dist[i][x], dist[i][y], method);
			dist[x][i] = lanceFormula(dist[x][i], dist[y][i], method);
		}
		dist[i][y] = INFINITY;
		dist[y][i] = INFINITY;
	}

	// Create new dendrogram containing x and y.
	Dendrogram newNode = createDnode(0);
	newNode->left = dendA[x];
	newNode->right = dendA[y];
	dendA[x] = newNode;
	dendA[y] = NULL;

}

static void findMin(double **dist, int numV, int index[2]) {
	double min = INFINITY;
	for (int i = 0; i < numV; i++) {
		for (int j = 0; j < numV; j++) {
			if ((dist[i][j] <= min) && (dist[i][j] > 0)) {
				min = dist[i][j];
				index[X] = i;
				index[Y] = j;
			}
		}
	}
	return;
}

static double lanceFormula(double ci, double cj, int method) {
	// If one of the distances is infinity,
	// the other distance is the distance.
	if (ci == INFINITY) {
		return cj;
	} else if (cj == INFINITY) {
		return ci;
	}

	double alpha = 0.5;
	double gamma;
	if (method == SINGLE_LINKAGE) {
		gamma = -0.5;
	} else if (method == COMPLETE_LINKAGE) {
		gamma = 0.5;
	}

	double modSubtract = (ci - cj);
	if ((modSubtract) < 0) {
		modSubtract *= -1;
	}

	double answer = ((alpha * ci) + (alpha * cj) + (gamma * modSubtract));
	return answer;
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram node) {
	if (node != NULL) {
        freeDendrogram(node->left);
        freeDendrogram(node->right);
        free(node);
    }
}

