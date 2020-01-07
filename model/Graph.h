//
// Created by Baoxing song on 2019-01-08.
//

#ifndef SONG_CNS_GRAPH_H
#define SONG_CNS_GRAPH_H
#include <list>
#include "Node.h"
#include <stack>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <limits.h>

class Graph {
    private:
        int V;
        // Pointer to an array containing adjacency lists
        std::list<Node> *adj;
        std::list<Node> *revAdj;

        // A function used by longestPath
        void topologicalSortUtil(int v, bool visited[], std::stack<int> &Stack);
    public:
        Graph(int V);   // Constructor
        // function to add an edge to graph
        void addEdge(int u, int v, int weight);
        // Finds longest distances from given source vertex
        void longestPath(int s, int **matrix, std::string *new_mems, int count, int no_seq, std::vector<int *> &LP_Matrix, std::vector<std::string> &LP_MEM);
};

#endif //SONG_CNS_GRAPH_H
