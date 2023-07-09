// Centrality Measures API implementation
// COMP2521 Assignment 2
// Derek Tran z5359557
// August 2021

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

static double wafFormula (int n, int N, double pathSum);
static double numPaths(int start, int end, PredNode **pred);

/**
 * Finds the closeness centrality for each vertex in the given graph and
 * returns the results in a NodeValues structure.
 */
NodeValues closenessCentrality(Graph g) {
	NodeValues nvs = {0};

	nvs.numNodes = GraphNumVertices(g);
	
	nvs.values = calloc(nvs.numNodes, sizeof(double));
	
	for (int i = 0; i < nvs.numNodes; i++) {
		// Get Shortestpaths from vertex 'i' to every other vertex
		ShortestPaths sps = dijkstra(g, i);

		int pathSum = 0;

		int reachableNodes = 0;

		// Count number of reachable nodes from vertex 'k', the node itself
		// is not counted
		for (int k = 0; k < sps.numNodes; k++) {
			if (sps.dist[k] != 0 && sps.dist[k] != INFINITY) {
				reachableNodes++;
				pathSum = pathSum + sps.dist[k];
			}
		}

		// Closeness is 0 if the path sum is 0
		if (pathSum <= 0.0) {
			nvs.values[i] = 0.0;
		} else {
			nvs.values[i] = wafFormula(reachableNodes, nvs.numNodes, pathSum);
		}
		freeShortestPaths(sps);
	}
	return nvs;
}

// Wasserman and Faust formula
static double wafFormula (int n, int N, double pathSum) {
	// Don't need (n - 1) in numerator because when increment the reachable nodes
	// the node itself is not included.
	return ((n * n) / ((N - 1) * pathSum));
}

/**
 * Finds  the  betweenness centrality for each vertex in the given graph
 * and returns the results in a NodeValues structure.
 */
NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs;
	
	nvs.numNodes = GraphNumVertices(g);
	
	nvs.values = calloc(nvs.numNodes, sizeof(double));

	// Iterate through all possible intermediate nodes 'v'
	for (int v = 0; v < nvs.numNodes; v++) {
		// Iterate through all possible starting nodes 's'
		for (int s = 0; s < nvs.numNodes; s++) {
			// Ensures that node to pass through isn't the starting node
			if (v == s) {
				continue;
			}
			ShortestPaths sps = dijkstra(g, s);
			// Iterate through all possible ending nodes 't'
			for (int t = 0; t < sps.numNodes; t++) {
				// If s != v != t, and t is reachable from s add to sum
				if (s != t && v != t) {
					// Multiplying the number of paths from s to intermediate node 
					// 'v' by the number of paths from 'v' to 't' gives us the 
					// number of paths from 's' to 't' that pass through 'v'
					// Denominator is  is the total number of shortest paths from
					// node s to t 
					double sv = numPaths(s, v, sps.pred);
					double vt = numPaths(v, t, sps.pred);
					double totalPaths = numPaths(s, t, sps.pred);
					if (totalPaths) {
						nvs.values[v] += ((sv * vt)/(totalPaths));
					}
				}		
			}
			freeShortestPaths(sps);
		}
	}
	return nvs;
}

// Recursively count the number of paths between two nodes
static double numPaths(int start, int end, PredNode **pred) {
	if (pred[end] == NULL) {
		return 0;
	}

	double pathCount = 0;

	// From the end node unwind the predecessor array until there are no
	// more predecessors.
	PredNode *curr = pred[end];
	if (curr->v == start) {
		pathCount = 1;
		curr = curr->next;
	}

	// If there are multiple predecessors, explore each pred node's 'path'
	// until the last pred in the list.
	while (curr != NULL) {
		if (curr->v == start) {
			pathCount++;
		}

		pathCount = (numPaths(start, curr->v , pred) + pathCount);
		curr = curr->next;
	}
	return pathCount;
}

/**
 * Finds  the  normalised  betweenness centrality for each vertex in the
 * given graph and returns the results in a NodeValues structure.
 */
NodeValues betweennessCentralityNormalised(Graph g) {
	NodeValues nvs = betweennessCentrality(g);
	double numV = GraphNumVertices(g);
	for (int i = 0; i < nvs.numNodes; i++) {
		nvs.values[i] = nvs.values[i] / ((numV - 1) * (numV - 2));
	}
	return nvs;
}

void showNodeValues(NodeValues nvs) {

}

void freeNodeValues(NodeValues nvs) {
	free(nvs.values);
}

