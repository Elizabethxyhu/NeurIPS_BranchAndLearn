#include <iostream>
#include <fstream>
#include <set>
// #include <limits>
#include "mincost.h"
#include "paramincost.h" 

// global variables 
// graph data: nodeNum, edgeNum, source, dest
int nodeNum, edgeNum, source, sink; 
// data for trainning
const int maxFeatureNum = 100;

int main(int argc, char *argv[]){
	if(argc != 5){
		std::cerr << "usage: ./test [graph file] [test data] [result file] [flow]" << std::endl;
		exit(0);
	}

 	// read graph information 
	std::ifstream infile;
  infile.open(argv[1]);
  infile >> nodeNum; 
  infile >> edgeNum; 
  infile >> source; 
  infile >> sink; 
  int graph[edgeNum][3]; // 2d-array for edges: (source, destination, capacity)
  for (int i = 0; i < edgeNum; i++) {
    for (int j = 0; j < 3; j++)
        infile >> graph[i][j];
  }
  infile.close(); 

  // read test data 
  int benchmarkNum;  
  infile.open(argv[2]);
  //infile >> benchmarkNum;
  benchmarkNum = 30;
  double predictedCost[benchmarkNum][edgeNum];
  double realCost[benchmarkNum][edgeNum];
  
  for(int i = 0; i < benchmarkNum; i++) {
  	double benchmarkId; 
  	for(int j = 0; j < edgeNum; j++) {
  		infile >> benchmarkId; 
      infile >> realCost[i][j]; 
      infile >> predictedCost[i][j]; 
  	}
  }

  int requireFlow = atoi(argv[4]); 
  // int iteration = atoi(argv[4]);
  
  // compute real objective
  double totalRealObj;
  double realObj[benchmarkNum];
  for(int i = 0; i < benchmarkNum; i++) {
    Graph resGraph; 
    resGraph.numVertices = nodeNum;
    resGraph.adj = new vector<Edge*>[nodeNum];
    for(int j = 0; j < edgeNum; j++ ) {
      Edge* tmpEdge1 = genEdge(graph[j][1], graph[j][2], graph[j][2], realCost[i][j], realCost[i][j]);
      Edge* tmpEdge2 = genEdge(graph[j][0], graph[j][2], 0, -1*realCost[i][j], -1*realCost[i][j]);
      tmpEdge1->counterEdge = tmpEdge2;
      tmpEdge2->counterEdge = tmpEdge1;
      resGraph.adj[graph[j][0]].push_back(tmpEdge1);
      resGraph.adj[graph[j][1]].push_back(tmpEdge2);
    }
    realObj[i] = calcMinCostFlow(resGraph, source, sink, requireFlow);
    // std::cout << "cost for instance " << i << ": " << realObj[i] << std::endl;
    if(realObj[i] == INT_MAX)
      std::cout << "instance " << i << " cannot send flow " << requireFlow << std::endl;
    totalRealObj += realObj[i]; 
  }

  double tempObj = 0.0; 
  for(int i = 0; i < benchmarkNum; i++) {
    Graph resGraph; 
    resGraph.numVertices = nodeNum;
    resGraph.adj = new vector<Edge*>[nodeNum];
    for(int j = 0; j < edgeNum; j++ ) {
      Edge* tmpEdge1 = genEdge(graph[j][1], graph[j][2], graph[j][2], predictedCost[i][j], realCost[i][j]);
      Edge* tmpEdge2 = genEdge(graph[j][0], graph[j][2], 0, -1*predictedCost[i][j], -1*realCost[i][j]);
      tmpEdge1->counterEdge = tmpEdge2;
      tmpEdge2->counterEdge = tmpEdge1;
      resGraph.adj[graph[j][0]].push_back(tmpEdge1);
      resGraph.adj[graph[j][1]].push_back(tmpEdge2);
    }
    double instanceCost = calcMinCostFlow(resGraph, source, sink, requireFlow);
    // std::cout << "predicted cost for instance " << i << ": " << instanceCost << ", real cost: " << realObj[i] << ", regret: " << instanceCost - realObj[i] << std::endl;
    tempObj += instanceCost; 
  }

  ofstream ofile;
  ofile.open(argv[3]);
  ofile << benchmarkNum << " " << argv[2] << ", Average regret error: " << (tempObj - totalRealObj)/benchmarkNum << std::endl;
  //ofile << benchmarkNum << " " << argv[2] << ", Average real Obj: " << totalRealObj/benchmarkNum << std::endl;
  ofile.close();

}
