#include <iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <cstdint>
#include <queue>
#ifndef MINCOST
#define MINCOST


using namespace std;

int INT_MAX = std::numeric_limits<std::int32_t>::max();

// An edge is represented as a struct
// Fields:
//      destination   -  denotes the ending node of an edge. For example, 'v' in u-->v
//      capacity      -  the maximum capacity of an edge
//      residualFlow  -  the residual amount of flow that can flow through the edge
//      counterEdge   -  a pointer to the counter edge in residual graph for performance optimization
struct Edge{
    int destination;
    int capacity;
    int residualFlow;
    double cost;
    double realCost; 
    Edge* counterEdge;
};

// A graph is represented as a struct
// Fields:
//      numVertices - denotes the number of vertices in the graph
//      adj         - Adjacency list : Collection of unordered lists one for each vertex
struct Graph{
    int numVertices;
    vector<Edge*> *adj;
};

// Graph g;
// Graph resGraph;

// Generates a new edge (allocating space dynamically) and returns a pointed to the edge
Edge* genEdge(int destination, int capacity, int residualFlow, double cost, double realCost){
    Edge* e1 = new Edge;
    e1->destination = destination;
    e1->capacity = capacity;
    e1->residualFlow = residualFlow;
    e1->cost = cost;
    e1->realCost = realCost; 
    return e1;
}

// Prints all the edges in the graph
// Output:
//      List of edges where each edge is represented by
//          u(start node)   v(end node)   flow   capacity
int printGraph(Graph g){
    for(int i=0; i<g.numVertices; i++){
        for(int j=0; j<g.adj[i].size(); j++){
            cout << i+1 << " " << g.adj[i][j]->destination+1 << " " << (g.adj[i][j]->capacity - g.adj[i][j]->residualFlow) << " " << g.adj[i][j]->cost << endl;
        }
    }
    return 0;
}

void printParams(int numVertices, int parentVertex[], int distance[]){
    cout << "Parents Vertex:\n";
    for(int i=0; i<numVertices; i++){
        cout << "\tparentVertex[" << i << "]: " << parentVertex[i] << " ";
    }
    cout << endl;
    cout << "Distance Vertex:\n";
    for(int i=0; i<numVertices; i++){
        cout << "\tdistance[" << i << "]: " << distance[i] << " ";
    }
    cout << endl;
}

// Detects the presence of negative cycles in graph
// Output:
//      -1          if no negative cycles present
//      node_num    index of a node in negative cycle
// Lists parentVertex and parentEdge are updated and can be used to reconstruct the negative cycle
int BFCycleDetection(Graph resGraph, int source, int parentVertex[], Edge* parentEdge[]){
    // Initialize variables that will be needed
    int cycle_node = -1;
    int numVertices = resGraph.numVertices;
    vector<Edge*> *adj = resGraph.adj;
    int distance[numVertices];
    // Initialize visited, parentVertex and distance
    for(int i=0; i<numVertices; i++){
        parentVertex[i] = -1;
        distance[i] = INT_MAX;
    }
    // BF - Relax edges repeatedly
    distance[source] = 0;
    for(int i=0; i<numVertices-1; i++){
        // loop on all edges
        for(int u=0; u<numVertices; u++){
            for(int e=0; e<adj[u].size(); e++){
                if(adj[u][e]->residualFlow > 0){
                    int v = adj[u][e]->destination;
                    int w = adj[u][e]->cost;
                    if(distance[v]>distance[u]+w){
                        distance[v] = distance[u]+w;
                        parentVertex[v] = u;
                        parentEdge[v] = adj[u][e];
                    }
                }
            }
        }
    }
    // Check for negative weight cycle - loop on all edges
    for(int u=0; u<numVertices; u++){
        for(int e=0; e<adj[u].size(); e++){
            if(adj[u][e]->residualFlow > 0){
                int v = adj[u][e]->destination;
                int w = adj[u][e]->cost;
                if(distance[v]>distance[u]+w){
                    return v;       // Negative cycle detected!
                }
            }
        }
    }
    return cycle_node;
}

// Cancels all negative cycles
// Output:
//      costSaved       amount of cost saved by cycle detection and cancellation
int cancelNegativeCycles(Graph resGraph){
    int costSaved=0, cyclePossible=1, u, v;
    Edge *te1, *te2;
    int numVertices = resGraph.numVertices;
    while(cyclePossible){
        cyclePossible=0;
        for(int i=0; i<numVertices; i++){
            int parent[resGraph.numVertices];
            Edge* parentEdge[resGraph.numVertices];
            int node_num = BFCycleDetection(resGraph, i, parent, parentEdge);
            if(node_num!=-1){               // A cycle is detected
                cyclePossible=1;
                // Calculate path flow
                int path_flow = INT_MAX;
                v=node_num; u = parent[v]; te1 = parentEdge[v];
                path_flow = min(path_flow, te1->residualFlow);
                for (v=u; v!=node_num; v=parent[v]){
                    u = parent[v];
                    te1 = parentEdge[v];
                    path_flow = min(path_flow, te1->residualFlow);
                }
                // Update graph by removing the cycle
                v=node_num; u = parent[v];
                te1 = parentEdge[v];
                te2 = te1->counterEdge;
                te1->residualFlow -= path_flow;
                te2->residualFlow += path_flow;
                costSaved += path_flow*(te1->cost);
                for (v=u; v != node_num; v=parent[v]){
                    u = parent[v];
                    te1 = parentEdge[v];
                    te2 = te1->counterEdge;
                    te1->residualFlow -= path_flow;
                    te2->residualFlow += path_flow;
                    costSaved += path_flow*(te1->cost);
                }
            }
        }
    }
    return -1*costSaved;
}


// Finds the shortest path from source to sink
// Output:
//      0           if no path exists from source to sink
//      1           if there is a path from source to sink
// Lists parentVertex and parentEdge are updated and can be used to reconstruct the shortest path
bool BF(Graph resGraph, int source, int sink, int parentVertex[], Edge* parentEdge[]){
    // Initialize variables that will be needed
    int numVertices = resGraph.numVertices;
    vector<Edge*> *adj = resGraph.adj;
    int distance[numVertices];
    // Initialize visited, parentVertex and distance
    for(int i=0; i<numVertices; i++){
        parentVertex[i] = -1;
        distance[i] = INT_MAX;
    }
    // printParams(numVertices, parentVertex, distance);
    // BF
    distance[source] = 0;
    for(int i=0; i<numVertices-1; i++){
        // loop on all edges
        for(int u=0; u<numVertices; u++){
            for(int e=0; e<adj[u].size(); e++){
                if(adj[u][e]->residualFlow > 0){
                    int v = adj[u][e]->destination;
                    int w = adj[u][e]->cost;
                    // std::cout << "find an edge (" << u << "," << v << ")" << std::endl;
                    // std::cout << "distance[v]: " << distance[v] << ", distance[u]+w: " << (distance[u]+w) << std::endl;
                    if(distance[u] != INT_MAX && distance[v]>distance[u]+w){
                        distance[v] = distance[u]+w;
                        parentVertex[v] = u;
                        parentEdge[v] = adj[u][e];
                    }
                    // printParams(numVertices, parentVertex, distance);
                }
            }
        }
    }
    // std::cout << "parentVertex[sink] == -1: " << (parentVertex[sink] == -1) << std::endl;
    if(parentVertex[sink] == -1){
        return false;
    }else{
        return true;
    }
}

// Calculates the cost of flow 'requiredFlow' from 's' to 't'
// Returns 'INT_MAX' if such a flow is not possible
double calcMinCostFlow(Graph resGraph, int s, int t, int requiredFlow){
    int u, v, currFlow=0;
    double runningCost=0;
    double currCost=0;
    Edge *te1, *te2;
    // Detect negative cycles and remove
    int benifit = cancelNegativeCycles(resGraph);
    if(benifit){
        cout << "Negative cycle detected and removed. Resulting cost benifit is " << benifit << endl;
    }
    // Run shortest path augmentation
    int parent[resGraph.numVertices];
    Edge* parentEdge[resGraph.numVertices];
    while (BF(resGraph, s, t, parent, parentEdge)){
        currCost=0;
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v]){
            u = parent[v];
            te1 = parentEdge[v];
            path_flow = min(path_flow, te1->residualFlow);
        }
        path_flow = min(path_flow, requiredFlow-currFlow);
        for (v=t; v != s; v=parent[v]){
            u = parent[v];
            te1 = parentEdge[v];
            te2 = te1->counterEdge;
            te1->residualFlow -= path_flow;
            te2->residualFlow += path_flow;
            runningCost += double(path_flow)*(te1->realCost);
            currCost += double(path_flow)*(te1->realCost);
        }
        currFlow += path_flow;
        // std::cout << "path: "; 
        // for (v=t; v != s; v=parent[v]){
        //     std::cout << v;
        //     if(v != s)
        //         std::cout << " <- ";
        // }
        // std::cout << "0 , flow: " << path_flow << ", cost: " << currCost << std::endl;
        if(currFlow == requiredFlow){
            break;
        }
    }
    if(currFlow == requiredFlow){
        return runningCost;
    }else{
        return INT_MAX;
    }
}

#endif

