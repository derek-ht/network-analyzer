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

static Dendrogram createDnode(int v);
static void mergeDend(double **dist, Dendrogram *dendA, int numV, int i[2], 
					  int method);
static void findMin(double **dist, int numV, int index[2]);

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
	
	// Initialize distance array make all 
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

// Merges 2 Dendrograms
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

	// Update dist array using the Lance William algorithm
	for (int i = 0; i < numV; i++) {
		if ((i != x) && (i != y)) {
			float distX = dist[i][x];
			float distX2 = dist[x][i];
			float distY = dist[i][y];
			 distY2 = dist[y][i];
			
			// If one of the distances is infinity change the distance
			// to the one that's not.
			if (distX == INFINITY) {
				dist[i][x] = distY;
				dist[x][i] = distY2;
			} else if (distY == INFINITY) {
				dist[i][x] = distX;
				dist[x][i] = distX2;
			// Complete Linkage
			// Choose the larger (max) of the 2 distances
			} else if (method == COMPLETE_LINKAGE) {
				if (distX < distY) {
					dist[i][x] = distY;
					dist[x][i] = distY2;
				} else {
					dist[i][x] = distX;
					dist[x][i] = distX2;
				}

			// Single Linkage
			// Choose the smaller (min) of the 2 distances
			} else if (method == SINGLE_LINKAGE) {
				if (distX < distY) {
					dist[i][x] = distX;
					dist[x][i] = distX2;
				} else {
					dist[i][x] = distY;
					dist[x][i] = distY2;
				}
			}
		}
		dist[i][y] = INFINITY;
		dist[y][i] = INFINITY;
	}

	// Merge the Dendrograms into one
	Dendrogram newNode = createDnode(0);
	newNode->left = dendA[x];
	newNode->right = dendA[y];
	dendA[x] = newNode;
	dendA[y] = NULL;
}

// Find the index of the minimum value and store the index into the
// array given
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

