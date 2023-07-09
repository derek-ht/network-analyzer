// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2
// Derek Tran z5359557
// August 2021

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

static Dendrogram createDnode(int v);
static void mergeDend(double **dist, Dendrogram *dendA, int numV, int *i, 
					  int method);
static double lanceFormula(double ci, double cj, int method);
static void findMin(double **dist, int numV, int *index);

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
	
	// Initialize distance array and make all distances infinity
	for (int i = 0; i < numVerts; i++) {
		dist[i] = malloc(numVerts * sizeof(double));
		for (int j = 0; j < numVerts; j++) {
			dist[i][j] = INFINITY;
		}
	}

	// Fill distance array with the distance '1/wt', where wt is
	// weight of the direct edge between 2 vertices.
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
	
    // Initialize Dendrogram array
	Dendrogram *dendA = malloc(numVerts * sizeof(struct DNode));
	for (int i = 0; i < numVerts; i++) {
		dendA[i] = createDnode(i);
	}
	
    // Merge dendrograms until the root has been reached, denda[0]
	for (int i = 1; i < numVerts; i++) {
		while (dendA[i] != NULL) {
			int index[2];
			findMin(dist, numVerts, index);
			mergeDend(dist, dendA, numVerts, index, method);
		}
	}

	// Free dynamically allocated memory
	for (int i = 0; i < numVerts; i++) {
		free(dist[i]);
	}
	free(dist);

	Dendrogram root = dendA[0];
	free(dendA);

	// Return root Dendrogram
	return root;
}

// Creates a new Dendrogram with the given vertice.
static Dendrogram createDnode(int v) {
	Dendrogram newNode = malloc(sizeof (struct DNode));
	newNode->vertex = v;
	newNode->left = NULL;
	newNode->right = NULL;
	return newNode;
}

// Find the index of the closest clusters clusters
// and store the index into the array given
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

// Merges dendrogram function
static void mergeDend(double **dist, Dendrogram *dendA, int numV, int *i, 
					  int method) {
	int x, y;
	// Ensure x is the larger value
	if (i[Y] < i[X]) {
		x = i[Y];
		y = i[X];
	} else {
		x = i[X];
		y = i[Y];
	}
	
	// Update dist array using the Lance William formula
	for (int i = 0; i < numV; i++) {
		if ((i != x) && (i != y)) {
			dist[i][x] = lanceFormula(dist[i][x], dist[i][y], method);
			dist[x][i] = lanceFormula(dist[x][i], dist[y][i], method);
		}
		dist[i][y] = INFINITY;
		dist[y][i] = INFINITY;
	}

	// Merge 2 clusters and add new merged dendrogram to the array dendA
	Dendrogram newNode = createDnode(0);
	newNode->left = dendA[x];
	newNode->right = dendA[y];
	dendA[x] = newNode;
	dendA[y] = NULL;
}

// Uses the general Lance Williams Method to determine which distance
// should be used.
static double lanceFormula(double ci, double cj, int method) {
	// If one of the distances is infinity change the distance
	// to the one that's not.
	if (ci == INFINITY) {
		return cj;
	} else if (cj == INFINITY) {
		return ci;
	}

	// Single Linkage
	// Choose the smaller (min) of the 2 distances
	if (method == SINGLE_LINKAGE) {
		if (ci < cj) {
			return ci;
		} else {
			return cj;
		}

	// Complete Linkage
	// Choose the larger (max) of the 2 distances
	} else {
		if (ci > cj) {
			return ci;
		} else {
			return cj;
		}
	}
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

