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
	if(argc != 7){
		std::cerr << "usage: ./train [graph file] [train data] [test data] [result file] [flow] [iteration]" << std::endl;
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

  // verify the graph 
  for (int i = 0; i < edgeNum; i++) {
    for (int j = 0; j < 3; j++) 
        std::cout << graph[i][j] << " ";
    std::cout << std::endl;
  }

  // read train data 
  int benchmarkNum, featureNum; 
  infile.open(argv[2]);
  // infile >> benchmarkNum; 
  // infile >> featureNum; 
  benchmarkNum = 70;
  featureNum = 8;
  double features[benchmarkNum][edgeNum][featureNum];
  double realCost[benchmarkNum][edgeNum];
  double realObj[benchmarkNum];
  for(int i = 0; i < benchmarkNum; i++) {
  	double benchmarkId; 
  	for(int j = 0; j < edgeNum; j++) {
  		infile >> benchmarkId; 
  		for(int k = 0; k < featureNum; k++)
  			infile >> features[i][j][k]; 
  		infile >> realCost[i][j]; 
  	}
  }
  infile.close();

  int requireFlow = atoi(argv[5]); 
  int iteration = atoi(argv[6]);
  
  // compute real objective
  double totalRealObj = 0; 
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
  //   if(realObj[i] == std::numeric_limits<int>::infinity())
  //     std::cout << "instance " << i << " cannot send flow " << requireFlow << std::endl;
  //   else
  //     std::cout << "realObj[" << i << "]=" << realObj[i] << std::endl;
    totalRealObj = totalRealObj + realObj[i];
  }

  // initialize linear prediction function 
  // double alpha[8] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  double alpha[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  // double alpha[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  
  // double alpha[4] = {0.1, 0.1, 0.1, 0.1};
  // double alpha[4] = {1, 1, 1, 1};
  // double alpha[4] = {-1, -1, -1, -1};


  bool converge = false;
  while(iteration > 0 && !converge) {
    converge=true;
    // for each feature, compute a piecewise function mapping alpha to regret value 
    // for(int feaCount = 6; feaCount < featureNum; feaCount++) {
    for(int feaCount = 0; feaCount < featureNum; feaCount++) {
      // assumption: we assume that the graph is an acylic graph 
      double tempObj = 0.0; 
      for(int i = 0; i < benchmarkNum; i++) {
        Graph resGraph; 
        resGraph.numVertices = nodeNum;
        resGraph.adj = new vector<Edge*>[nodeNum];
        for(int j = 0; j < edgeNum; j++ ) {
          // compute the cost with predict function
          double tempCost = 0.0; 
          for(int k = 0; k < featureNum; k++)
            tempCost += features[i][j][k] * alpha[k];
          // std::cout << "predicted cost for edge " << j << " is " << tempCost << std::endl;
          Edge* tmpEdge1 = genEdge(graph[j][1], graph[j][2], graph[j][2], tempCost, realCost[i][j]);
          Edge* tmpEdge2 = genEdge(graph[j][0], graph[j][2], 0, -1*tempCost, -1*realCost[i][j]);
          tmpEdge1->counterEdge = tmpEdge2;
          tmpEdge2->counterEdge = tmpEdge1;
          resGraph.adj[graph[j][0]].push_back(tmpEdge1);
          resGraph.adj[graph[j][1]].push_back(tmpEdge2);
        }
        double instanceCost = calcMinCostFlow(resGraph, source, sink, requireFlow);
        // std::cout << "predicted cost for instance " << i << "with alpha " << feaCount << " being " << alpha[feaCount] << ": " << instanceCost << std::endl;
        tempObj += instanceCost; 
      }
      // std::cout << "optimize alpha " << feaCount << ", original value: " << alpha[feaCount] << ", original total cost: " << std::fixed << tempObj << std::endl; 
      std::cout << "optimize alpha " << feaCount << ", original value: " << alpha[feaCount] << ", loss: " << std::fixed << tempObj - totalRealObj << std::endl; 

      piecewiseConst paraCost; 

      // std::cout << "total number of benchmarkNum: " << benchmarkNum << std::endl;
      // for(int i = 338; i < benchmarkNum; i++) {
      for(int i = 0; i < benchmarkNum; i++) {
        // construct the parameterized graph 
        // construct the parameterized table 
        // std::cout << "number of benchmark: " << i << "/" << benchmarkNum << std::endl;
        double lowerBound = -std::numeric_limits<double>::infinity(); 
        double upperBound = std::numeric_limits<double>::infinity();
        ParaGraph resGraph; 
        resGraph.numVertices = nodeNum; 
        resGraph.adj = new vector<ParaEdge*>[nodeNum]; 
        for(int j = 0; j < edgeNum; j++) {
          // computet the linear function of the (feaCount)^th 
          double w = features[i][j][feaCount];
          double b = 0.0; 
          for(int k = 0; k < featureNum; k++) {
            if(k != feaCount)
              b += features[i][j][k] * alpha[k];
          }
          // check linear function 
          // std::cout << "predict function for edge " << j << " is: " << w << "* alpha + " << b << ", value at alpha = " << alpha[feaCount] << " is " << (b + w*alpha[feaCount]) << std::endl;
          ParaEdge* tmpEdge1 = genEdge(graph[j][1], graph[j][2], graph[j][2], realCost[i][j], w, b);
          ParaEdge* tmpEdge2 = genEdge(graph[j][0], graph[j][2], 0, -realCost[i][j], -w, -b);
          tmpEdge1->counterEdge = tmpEdge2;
          tmpEdge2->counterEdge = tmpEdge1;
          resGraph.adj[graph[j][0]].push_back(tmpEdge1);
          resGraph.adj[graph[j][1]].push_back(tmpEdge2);
          
          // restrict the range of alpha to enforce the cost to be positive 
          if(w > 0) 
            lowerBound = (-b)/w;
          else if(w < 0)
            upperBound = (-b)/w;
        }

        // std::cout << "lowerBound: " << lowerBound << ", upperBound: " << upperBound << std::endl;
        piecewiseConst returnCost = paraMinCost(resGraph, source, sink, requireFlow, lowerBound, upperBound); 

        // adding two piecewise linear function 
        piecewiseConst::iterator iter1 = paraCost.begin(); 
        piecewiseConst::iterator iter2 = returnCost.begin();
        piecewiseConst resultCost; 
        double last1 = 0;
        double last2 = 0; 
        while( iter1 != paraCost.end() || iter2 != returnCost.end() ){
          if( iter1 == paraCost.end() ) {
            resultCost[(*iter2).first] = (*iter2).second + last1; 
            iter2++; 
          }
          else if( iter2 == returnCost.end() ) {
            resultCost[(*iter1).first] = (*iter1).second+ last2; 
            iter1++; 
          }
          else {
            double nextPt = min((*iter1).first, (*iter2).first); 
            resultCost[nextPt] = (*iter1).second + (*iter2).second;
            if((*iter1).first < (*iter2).first || (*iter1).first == (*iter2).first || abs((*iter1).first - (*iter2).first) < 1e-10 ){
              last1=(*iter1).second;
              iter1++;
            }
            if((*iter2).first < (*iter1).first || (*iter2).first == (*iter1).first || abs((*iter2).first - (*iter1).first) < 1e-10 ){
              last2=(*iter2).second;
              iter2++;
            }
          }
        }
        paraCost = resultCost; 

        double currentValue = 0; 
        for(piecewiseConst::iterator iter = returnCost.begin(); iter != returnCost.end(); iter++) {
          if(iter->first > alpha[feaCount] ){
            currentValue = iter->second; 
            break;
          }
        }
        // std::cout << "predicted cost for instance " << i << ": " << currentValue << std::endl;

        // double tempPt = -std::numeric_limits<double>::infinity(); 
        // std::cout << "piecewise constant function of instance " << i << " for alpha " << feaCount << std::endl;
        // for(piecewiseConst::iterator iter = returnCost.begin(); iter != returnCost.end(); iter++) {
        //     std::cout << "\t(" << tempPt << "," << iter->first << "): " << iter->second << std::endl;
        //     tempPt = iter->first; 
        // }

        // tempPt = -std::numeric_limits<double>::infinity(); 
        // std::cout << "piecewise constant function after instance " << i << " for alpha " << feaCount << std::endl;
        // for(piecewiseConst::iterator iter = paraCost.begin(); iter != paraCost.end(); iter++) {
        //     std::cout << "\t(" << tempPt << "," << iter->first << "): " << iter->second << std::endl;
        //     tempPt = iter->first; 
        // }
        // cin.ignore();
      }

      // double tempPt = -std::numeric_limits<double>::infinity(); 
      // std::cout << "piecewise constant function for alpha " << feaCount << std::endl;
      // for(piecewiseConst::iterator iter = paraCost.begin(); iter != paraCost.end(); iter++) {
      //     std::cout << "\t(" << tempPt << "," << iter->first << "): " << std::fixed << iter->second << std::endl;
      //     tempPt = iter->first; 
      // }

      double lowerBound=-std::numeric_limits<double>::infinity();
      double minValue=std::numeric_limits<double>::infinity();
      std::pair<double, double> minRange;
      for(piecewiseConst::iterator iter = paraCost.begin(); iter != paraCost.end(); iter++) {
        if(minValue > (*iter).second) {
          minValue=(*iter).second;
          minRange=std::make_pair(lowerBound, (*iter).first); 
        }
        lowerBound=(*iter).first;
      }
      
      // std::cout << "alpha " << feaCount << " minRange: (" << minRange.first << ", " << minRange.second << "), minValue: " << std::fixed << minValue << std::endl;
      
      if( !(alpha[feaCount] > minRange.first && alpha[feaCount] < minRange.second) ) {
        // update alpha to minimize the cost 
        double backup = alpha[feaCount];
        converge=false;
        if(minRange.first != -std::numeric_limits<double>::infinity() && minRange.second != std::numeric_limits<double>::infinity()) 
          alpha[feaCount] = (minRange.first + minRange.second) / 2;
        else if (minRange.first != -std::numeric_limits<double>::infinity()) 
          alpha[feaCount] = minRange.first + DBL_PRECISION; 
        else if (minRange.second != -std::numeric_limits<double>::infinity())
          alpha[feaCount] = minRange.second + DBL_PRECISION; 

        // std::cout << "change alpha " << feaCount << " from " << backup << " to " << alpha[feaCount] << std::endl;
        // std::cout << "current alpha = [";
        // for(int j=0; j != featureNum; j++) {
        //   std::cout << alpha[j] << ","; 
        // }
        // std::cout << "]" << std::endl;

        // std::cin.ignore();
      }
    }
    iteration--;
  }

  std::cout << "alpha = ["; 
  for(int i=0; i < featureNum; i++) {
    std::cout << alpha[i];
    if(i+1 != featureNum)
      std::cout << ",";
  }
  std::cout << "]" << std::endl;

  // read test data 
  // int testbenchmarkNum, testfeatureNum; 
  infile.open(argv[3]);
  // infile >> benchmarkNum; 
  // infile >> featureNum; 
  benchmarkNum = 30;
  featureNum = 8;
  double testFeatures[benchmarkNum][edgeNum][featureNum];
  double testRealCost[benchmarkNum][edgeNum];
  ofstream ofile;
  ofile.open(argv[4]); 
  // std::cout << "test data" << std::endl;
  //ofile << benchmarkNum << std::endl;
  for(int i = 0; i < benchmarkNum; i++) {
    for(int j = 0; j < edgeNum; j++) {
      double benchmarkId; 
      infile >> benchmarkId; 
      double predictedCost = 0.0; 
      for(int k = 0; k < featureNum; k++) {
        infile >> features[i][j][k]; 
        // std::cout << features[i][j][k] << " "; 
        predictedCost += features[i][j][k] * alpha[k]; 
      }
      // std::cout << std::endl;
      infile >> realCost[i][j]; 
      ofile << i  << " " << realCost[i][j] << " " << predictedCost << std::endl;
    }
  }
  ofile.close();
  infile.close();

}

