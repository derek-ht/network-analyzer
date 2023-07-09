/*Lance-Williams Algorithm for Agglomerative Clustering
  Written by 
  COMP2521 2019 T1
*/
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "LanceWilliamsHAC.h"
#include "Graph.h"

#define INFINITY DBL_MAX

// #define numVertices GraphNumVertices

typedef struct coordinates {
  int x;
  int y;
} coordinates;

static bool incomplete(Dendrogram *array, int num_verticies);
static coordinates find_min(double **distances, int num_verticies);
static void combine_dendrograms (double **distances, Dendrogram *array, int num_verticies, int x, int y, int method);
static double lance_williams(double ik, double jk, int method);


/* 
 * Finds Dendrogram using Lance-Williams algorithm (as discussed in the specs) 
   for the given graph g and the specified method for agglomerative clustering.
   Method value 1 represents 'Single linkage', and 2 represents 'Complete linkage'.
   For this assignment, you only need to implement the above two methods. 
   The function returns 'Dendrogram' structure (binary tree) with the required information.
 * 
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {

  int num_verticies = GraphNumVertices(g);
  double **distances = malloc(num_verticies * sizeof(*distances));

  // Set all distances in 2d array to inf.
  // Set all distances from vertice to itself to 0.
  for (int i = 0; i < num_verticies; i++) {
    distances[i] = malloc(num_verticies * sizeof(double));
    for (int j = 0; j < num_verticies; j++) {
      
      if (i == j) {
        distances[i][j] = 0;
      } else {
        distances[i][j] = INFINITY;
      }
    }
  }

  // Check if vertices are adjacent to other vertices
  // and update the distances.
  for (int vertex = 0; vertex < num_verticies; vertex++) {
    
    for (AdjList out = GraphOutIncident(g, vertex); out != NULL; out = out->next) {
      
      if ((1.0 / out->weight) < distances[vertex][out->v]) {

        distances[vertex][out->v] = (1.0 / out->weight);
      }
    }

    for (AdjList in = GraphInIncident(g, vertex); in != NULL; in = in->next) {

      if ((1.0 / in->weight) < distances[vertex][in->v]) {

        distances[vertex][in->v] = (1.0 / in->weight);
      }
    }
  }

  // Array of dendrograms.
  Dendrogram *array = malloc(num_verticies * sizeof(Dendrogram));

  // Initialise all dendrograms.
  for (int i = 0; i < num_verticies; i++) {
    array[i] = malloc(sizeof(DNode));
    array[i]->vertex = i;
    array[i]->left = NULL;
    array[i]->right = NULL;
  }

  // While all vertices are not in the same dendrogram,
  // continue mergine dendrograms.
  while (incomplete(array, num_verticies)) {
    coordinates min = find_min(distances, num_verticies);
    combine_dendrograms(distances, array, num_verticies, min.x, min.y, method);
  }

  // Final dendrogram will be at index 0 of array.
  Dendrogram ans = array[0];
  
  // Free 2d array used to calculate distances.
  for (int i = 0; i < num_verticies; i++) {
    free(distances[i]);
  }
  free(distances);

  return ans;
}


void freeDendrogram(Dendrogram d) {

  if (d->left != NULL) {
    freeDendrogram(d->left);
  }

  if (d->right != NULL) {
    freeDendrogram(d->right);
  }

  free(d);
}

// Check if all vertices are in the same dendrogram.
static bool incomplete(Dendrogram *array, int num_verticies) {

  // If there are dendrograms found in the array (besides index 0 of array),
  // return true.
  for (int i = 1; i < num_verticies; i++) {
    if (array[i] != NULL) {
      return true;
    }
  }

  return false;
}

// Find the minimum distance in a 2d array.
static coordinates find_min(double **distances, int num_verticies){
  coordinates answer;
  double curr_min = INFINITY;

  for (int i = 0; i < num_verticies; i++) {
    for (int j = 0; j < num_verticies; j++) {
      if ((distances[i][j] <= curr_min) && (distances[i][j] > 0)) {
        curr_min = distances[i][j];
        answer.x = i;
        answer.y = j;
      }
    }
  }
  return answer;
}

// Merge 2 dendrograms.
static void combine_dendrograms (double **distances, Dendrogram *array, int num_verticies, int x, int y, int method) {
  
  // If x is greater than y,
  // swap x and y.
  if (y < x) {
    int tmp = x;
    x = y;
    y = tmp;
  }

  // Update distances for all vertices adjacent to x and y
  for (int i = 0; i < num_verticies; i++) {
    if ((i != x) && (i != y) && (array[i] != NULL)) {
      distances[i][x] = lance_williams(distances[i][x], distances[i][y], method);
      distances[x][i] = lance_williams(distances[x][i], distances[y][i], method);
    }
    distances[i][y] = INFINITY;
    distances[y][i] = INFINITY;
  }

  // Create new dendrogram containing x and y.
  Dendrogram new = malloc(sizeof(DNode));
  new->left = array[x];
  new->right = array[y];
  array[x] = new;
  array[y] = NULL;
}

static double lance_williams(double ik, double jk, int method) {

  // If one of the distances is infinity,
  // the other distance is the distance.
  if (ik == INFINITY) {
    return jk;
  } else if (jk == INFINITY) {
    return ik;
  }

  double var = 0;
  if (method == 1) {
    var = (-0.5);
  } else if (method == 2) {
    var = (0.5);
  }

  double abs = (ik - jk);
  if (abs < 0) {
    abs = -abs;
  }

  double answer = (((0.5) * ik) + ((0.5) * jk) + (var * abs));
  return answer;
}