// Centrality Measures API implementation
// COMP2521 Assignment 2

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

//static double wafFormula (int n, int N, double pathSum);

// n is the number of nodes reachable from 'ver'
// N is the number of vertices in the graph
// if a node is not connected to any other node
NodeValues closenessCentrality(Graph g) {
	NodeValues nvs = {0};

	nvs.numNodes = GraphNumVertices(g);
	
	nvs.values = malloc(nvs.numNodes * sizeof (double));

	for (int i = 0; i < nvs.numNodes; i++) {
		// Get Shortestpaths from vertex 'i' to every other vertex
		ShortestPaths sps = dijkstra(g, i);

		double pathSum = 0.0;


		for (int i = 0; i < nvs.numNodes; i++) {
            if (curr != NULL) {
                reachableNodes++;
            }
        }

		for (int i = 0; i < sps.numNodes; i++) {
			
			pathSum += sps.dist[i];
		}

		if (pathSum == 0.0) {
			// printf("BBBBBBBBB\n");
			nvs.values[i] = 0.0;
		} else {
			// printf("CCCCCCCCCC\n");
			// printf("%f\n%d\n%f\n", reachableNodes, nvs.numNodes, pathSum);
			nvs.values[i] = //wafFormula(reachableNodes, nvs.numNodes, pathSum);
			((reachableNodes * reachableNodes)/((nvs.numNodes - 1) * pathSum));
		}
	}
	return nvs;
}

// static double wafFormula (int n, int N, double pathSum) {
// 	return ((n * n) / ((N - 1) * pathSum));
// }

NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs = {0};
	
	return nvs;
}

NodeValues betweennessCentralityNormalised(Graph g) {
	NodeValues nvs = {0};
	return nvs;
}

void showNodeValues(NodeValues nvs) {

}

void freeNodeValues(NodeValues nvs) {

}
